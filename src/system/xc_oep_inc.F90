!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

! ---------------------------------------------------------
!> This file handles the evaluation of the OEP potential, in the KLI or full OEP
!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)
!!
!! This file has to be outside the module xc, for it requires the Hpsi.
!! This is why it needs the xc_functl module. I prefer to put it here since
!! the rest of the Hamiltonian module does not know about the gory details
!! of how xc is defined and calculated.
subroutine X(xc_oep_calc)(oep, xcs, apply_sic_pz, gr, hm, st, ex, ec, vxc)
  type(xc_oep_t),      intent(inout) :: oep
  type(xc_t),          intent(in)    :: xcs
  logical,             intent(in)    :: apply_sic_pz
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(inout) :: st
  FLOAT,               intent(inout) :: ex, ec
  FLOAT, optional,     intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  FLOAT :: eig
  integer :: is, ist, ixc
  logical, save :: first = .true.

  if(oep%level == XC_OEP_NONE) return

  call profiling_in(C_PROFILING_XC_OEP)
  PUSH_SUB(X(xc_oep_calc))

  ! initialize oep structure
  SAFE_ALLOCATE(oep%eigen_type (1:st%nst))
  SAFE_ALLOCATE(oep%eigen_index(1:st%nst))
  SAFE_ALLOCATE(oep%X(lxc)(1:gr%mesh%np, st%st_start:st%st_end))
  SAFE_ALLOCATE(oep%uxc_bar    (1:st%nst))

  ! this part handles the (pure) orbital functionals
  spin: do is = 1, min(st%d%nspin, 2)
    oep%X(lxc) = M_ZERO

    ! get lxc
    functl_loop: do ixc = 1, 2
      if(xcs%functl(ixc, 1)%family .ne. XC_FAMILY_OEP) cycle

      eig = M_ZERO
      select case(xcs%functl(ixc,1)%id)
      case(XC_OEP_X)
        call X(oep_x) (gr, st, is, oep, eig, xcs%exx_coef)
        ex = ex + eig
      end select
    end do functl_loop

    ! SIC a la PZ is handled here
    if(apply_sic_pz) then
      call X(oep_sic) (xcs, gr, st, is, oep, ex, ec)
    end if

    ! get the HOMO state
    call xc_oep_AnalyzeEigen(oep, st, is)

    ! calculate uxc_bar for the occupied states
    do ist = st%st_start, st%st_end
      oep%uxc_bar(ist) = X(mf_dotp)(gr%mesh, R_CONJ(st%X(psi)(1:gr%mesh%np, 1, ist, is)), oep%X(lxc)(1:gr%mesh%np, ist))
    end do

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
      do ist = 1, st%nst
        call MPI_Bcast(oep%uxc_bar(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    if(present(vxc)) then
      ! solve the KLI equation
      if(oep%level .ne. XC_OEP_FULL .or. first) then
        first = .false.

        oep%vxc = M_ZERO
        call X(xc_KLI_solve) (gr%mesh, st, is, oep)
      end if

      ! if asked, solve the full OEP equation
      if(oep%level == XC_OEP_FULL) then
        call X(xc_oep_solve)(gr, hm, st, is, vxc(:,is), oep)
      end if

      vxc(1:gr%mesh%np, is) = vxc(1:gr%mesh%np, is) + oep%vxc(1:gr%mesh%np)
    end if
  end do spin

  SAFE_DEALLOCATE_P(oep%eigen_type)
  SAFE_DEALLOCATE_P(oep%eigen_index)
  SAFE_DEALLOCATE_P(oep%X(lxc))
  SAFE_DEALLOCATE_P(oep%uxc_bar)

  POP_SUB(X(xc_oep_calc))
  call profiling_out(C_PROFILING_XC_OEP)
end subroutine X(xc_OEP_calc)


! ---------------------------------------------------------
subroutine X(xc_oep_solve) (gr, hm, st, is, vxc, oep)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: is
  FLOAT,               intent(inout) :: vxc(:) !< (gr%mesh%np)
  type(xc_oep_t),      intent(inout) :: oep

  integer :: iter, ist
  FLOAT :: vxc_bar, ff
  FLOAT, allocatable :: ss(:), vxc_old(:)
  R_TYPE, allocatable :: bb(:,:)

  call profiling_in(C_PROFILING_XC_OEP_FULL)
  PUSH_SUB(X(xc_oep_solve))

  SAFE_ALLOCATE(     bb(1:gr%mesh%np, 1:1))
  SAFE_ALLOCATE(     ss(1:gr%mesh%np))
  SAFE_ALLOCATE(vxc_old(1:gr%mesh%np))

  vxc_old(1:gr%mesh%np) = vxc(1:gr%mesh%np)

  if(.not. lr_is_allocated(oep%lr)) then
    call lr_allocate(oep%lr, st, gr%mesh)
    ! initialize to something non-zero
    oep%lr%X(dl_psi)(:,:, :, :) = M_ONE
  end if

  ! fix xc potential (needed for Hpsi)
  vxc(1:gr%mesh%np) = vxc_old(1:gr%mesh%np) + oep%vxc(1:gr%mesh%np)

  do iter = 1, oep%scftol%max_iter
    ! iteration over all states
    ss = M_ZERO
    do ist = 1, st%nst
      ! evaluate right-hand side
      vxc_bar = dmf_dotp(gr%mesh, (R_ABS(st%X(psi)(1:gr%mesh%np, 1, ist, is)))**2, oep%vxc(1:gr%mesh%np))
      bb(1:gr%mesh%np, 1) = -(oep%vxc(1:gr%mesh%np) - (vxc_bar - oep%uxc_bar(ist))) * &
        R_CONJ(st%X(psi)(1:gr%mesh%np, 1, ist, is)) + oep%X(lxc)(1:gr%mesh%np, ist) 

      call X(lr_orth_vector) (gr%mesh, st, bb, ist, is, R_TOTYPE(M_ZERO))

      call X(solve_HXeY) (oep%solver, hm, gr, st, ist, is, oep%lr%X(dl_psi)(:,:, ist, is), bb, &
           R_TOTYPE(-st%eigenval(ist, is)), CNST(1e-6))
      
      call X(lr_orth_vector) (gr%mesh, st, oep%lr%X(dl_psi)(:,:, ist, is), ist, is, R_TOTYPE(M_ZERO))

      ! calculate this funny function ss
      ss(1:gr%mesh%np) = ss(1:gr%mesh%np) + &
        M_TWO*R_REAL(oep%lr%X(dl_psi)(1:gr%mesh%np, 1, ist, is)*st%X(psi)(1:gr%mesh%np, 1, ist, is))
    end do

    oep%vxc(1:gr%mesh%np) = oep%vxc(1:gr%mesh%np) + oep%mixing*ss(1:gr%mesh%np)

    do ist = 1, st%nst
      if(oep%eigen_type(ist) == 2) then
        vxc_bar = dmf_dotp(gr%mesh, (R_ABS(st%X(psi)(1:gr%mesh%np, 1, ist, is)))**2, oep%vxc(1:gr%mesh%np))
        oep%vxc(1:gr%mesh%np) = oep%vxc(1:gr%mesh%np) - (vxc_bar - oep%uxc_bar(ist))
      end if
    end do

    ff = dmf_nrm2(gr%mesh, ss)
    if(ff < oep%scftol%conv_abs_dens) exit
  end do

  write(message(1), '(a,i4,a,es14.6)') "Info: After ", iter, " iterations, the OEP converged to ", ff
  message(2) = ''
  call messages_info(2)

  vxc(1:gr%mesh%np) = vxc_old(1:gr%mesh%np)
  SAFE_DEALLOCATE_A(bb)
  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(vxc_old)

  POP_SUB(X(xc_oep_solve))
  call profiling_out(C_PROFILING_XC_OEP_FULL)
end subroutine X(xc_oep_solve)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
