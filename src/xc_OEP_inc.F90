!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!!! This file handles the evaluation of the OEP potential, in the KLI or full OEP
!!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)

!!! This file has to be outside the module xc, for it requires the Hpsi
!!! This is why it needs the xc_functl module. I prefer to put it here since
!!! the rest of the hamiltonian module does not know about the gory details
!!! of how xc is defined and calculated.

subroutine X(xc_oep_calc)(oep, xcs, m, f_der, h, st, vxc, ex, ec)
  use xc_functl

  type(xc_oep_type),      intent(inout) :: oep
  type(xc_type),          intent(in)    :: xcs
  type(mesh_type),        intent(in)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(hamiltonian_type), intent(inout) :: h
  type(states_type),      intent(inout) :: st
  FLOAT,                  intent(inout) :: vxc(m%np, st%d%nspin)
  FLOAT,                  intent(inout) :: ex, ec
  
  FLOAT :: e
  integer :: is, ist, ixc, ifunc, ierr
  logical, save :: first = .true.

  if(oep%level == XC_OEP_NONE) return

  call push_sub('h_xc_oep')

  ! initialize oep structure
  allocate(oep%eigen_type(st%nst))
  allocate(oep%eigen_index(st%nst))
  allocate(oep%X(lxc)(m%np, st%st_start:st%st_end))
  allocate(oep%uxc_bar(st%nst))

  ! this part handles the (pure) orbital functionals
  spin: do is = 1, min(st%d%nspin, 2)
    oep%X(lxc) = M_ZERO

    ! get lxc
    functl_loop: do ixc = 1, 2
      if(xcs%functl(ixc, 1)%family.ne.XC_FAMILY_OEP) cycle

      e = M_ZERO
      select case(xcs%functl(ixc,1)%id)
      case(XC_OEP_X)
        call X(oep_x) (m, f_der, st, is, oep, e)
        ex = ex + e
      end select
    end do functl_loop
  
    ! SIC a la PZ is handled here
    if(oep%sic_correction.ne.0) then
      call X(oep_sic) (xcs, m, f_der, st, is, oep, vxc, ex, ec)
    end if
    
    ! get the HOMO state
    call xc_oep_AnalizeEigen(oep, st, is)

    ! calculate uxc_bar for the occupied states
    do ist = st%st_start, st%st_end
      oep%uxc_bar(ist) = sum(R_REAL(st%X(psi)(:, 1, ist, is) * oep%X(lxc)(:, ist))*m%vol_pp(:))
    end do

#if defined(HAVE_MPI)
    if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      call mpi_barrier(mpi_comm_world, ierr)
      do ist = 1, st%nst
         call mpi_bcast(oep%uxc_bar(ist), 1, MPI_FLOAT, st%node(ist), MPI_COMM_WORLD, ierr)
      enddo
    endif
#endif
    
    ! solve the KLI equation
    if(oep%level.ne.XC_OEP_FULL.or.first) then
      first = .false.

      oep%vxc = M_ZERO
      call X(xc_KLI_solve) (m, st, is, oep)
    end if

    ! if asked, solve the full OEP equation
    if(oep%level == XC_OEP_FULL) then
      call X(xc_oep_solve)(m, f_der, h, st, is, vxc(:,is), oep)
    end if

    vxc(:, is) = vxc(:, is) + oep%vxc(:)
  end do spin

  deallocate(oep%eigen_type, oep%eigen_index, oep%X(lxc), oep%uxc_bar)

  call pop_sub()
end subroutine X(xc_OEP_calc)


subroutine X(xc_oep_solve) (m, f_der, h, st, is, vxc, oep)
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(hamiltonian_type), intent(inout) :: h
  type(states_type),      intent(IN)    :: st
  integer,                intent(in)    :: is
  FLOAT,                  intent(inout) :: vxc(m%np)
  type(xc_oep_type),      intent(inout) :: oep

  integer :: iter, ist, ierr
  FLOAT :: vxc_bar, f
  FLOAT, allocatable :: s(:), vxc_old(:)
  R_TYPE, allocatable :: b(:,:)

  allocate(b(m%np, 1), s(m%np), vxc_old(m%np))

  vxc_old = vxc(:)

  ierr = X(lr_alloc_psi) (st, m, oep%lr)
  do ist = 1, st%nst
    call X(lr_orth_vector) (m, st, oep%lr%X(dl_psi)(:,:, ist, is), is)
  end do

  ! fix xc potential (needed for Hpsi)
  vxc(:) = vxc_old(:) + oep%vxc(:)

  do iter = 1, 10
    ! iteration over all states
    s = M_ZERO
    do ist = 1, st%nst
      ! evaluate right-hand side
      vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:) * m%vol_pp(:))
      b(:,1) =  -(oep%vxc(:) - (vxc_bar - oep%uxc_bar(ist)))*R_CONJ(st%X(psi)(:, 1, ist, is)) &
         + oep%X(lxc)(:, ist)

      call X(lr_orth_vector) (m, st, b, is)
      
      ! and we now solve the equation [h-eps_i] psi_i = b_i
      call X(lr_solve_HXeY) (oep%lr, h, m, f_der, st%d, is, oep%lr%X(dl_psi)(:,:, ist, is), b, &
         R_TOTYPE(-st%eigenval(ist, is)))

      call X(lr_orth_vector) (m, st, oep%lr%X(dl_psi)(:,:, ist, is), is)

      ! calculate this funny function s
      s(:) = s(:) + M_TWO*R_REAL(oep%lr%X(dl_psi)(:, 1, ist, is)*st%X(psi)(:, 1, ist, is))
    end do

    oep%vxc(:) = oep%vxc(:) + oep%mixing*s(:)

    do ist = 1, st%nst
      if(oep%eigen_type(ist) == 2) then
        vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:) * m%vol_pp(:))
        oep%vxc(:) = oep%vxc(:) - (vxc_bar - oep%uxc_bar(ist))
      end if
    end do

    f = dmf_nrm2(m, s)
    if(f < oep%lr%conv_abs_dens) exit
  end do

  write(message(1), '(a,i4,a,es14.6)') "Info: After ", iter, " iterations, the OEP converged to ", f
  message(2) = ''
  call write_info(2)

  vxc(:) = vxc_old(:)
  deallocate(b, s, vxc_old)

end subroutine X(xc_oep_solve)
