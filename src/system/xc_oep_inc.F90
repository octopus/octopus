!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!!
! ---------------------------------------------------------
!> This file handles the evaluation of the OEP potential, in the KLI or full OEP
!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)
!!
!! This file has to be outside the module xc, for it requires the Hpsi.
!! This is why it needs the xc_functl module. I prefer to put it here since
!! the rest of the Hamiltonian module does not know about the gory details
!! of how xc is defined and calculated.
subroutine X(xc_oep_calc)(oep, namespace, xcs, apply_sic_pz, gr, hm, st, ex, ec, vxc)
  type(xc_oep_t),           intent(inout) :: oep
  type(namespace_t),        intent(in)    :: namespace
  type(xc_t),               intent(inout) :: xcs
  logical,                  intent(in)    :: apply_sic_pz
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(inout) :: st
  FLOAT,                    intent(inout) :: ex, ec
  FLOAT, optional,          intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  FLOAT :: eig
  integer :: is, ist, ixc, nspin_, isp, idm, jdm
  logical, save :: first = .true.
  R_TYPE, allocatable :: psi(:)

  if(oep%level == XC_OEP_NONE) return

  call profiling_in(C_PROFILING_XC_OEP, 'XC_OEP')
  PUSH_SUB(X(xc_oep_calc))

  ! initialize oep structure
  nspin_ = min(st%d%nspin, 2)
  SAFE_ALLOCATE(oep%eigen_type (1:st%nst))
  SAFE_ALLOCATE(oep%eigen_index(1:st%nst))
  SAFE_ALLOCATE(oep%X(lxc)(1:gr%mesh%np, st%st_start:st%st_end, 1:nspin_))
  SAFE_ALLOCATE(oep%uxc_bar(1:st%nst, 1:nspin_))

  ! this part handles the (pure) orbital functionals
  oep%X(lxc) = M_ZERO
  spin: do is = 1, nspin_
    !
    ! distinguish between 'is' being the spin_channel index (collinear)
    ! and being the spinor (noncollinear)
    if (st%d%ispin==SPINORS) then
      isp = 1
      idm = is
    else
      isp = is
      idm = 1
    end if
    ! get lxc
    functl_loop: do ixc = 1, 2
      if(xcs%functional(ixc, 1)%family /= XC_FAMILY_OEP) cycle

      eig = M_ZERO
      select case(xcs%functional(ixc,1)%id)
      case(XC_OEP_X)
        sum_comp: do jdm = 1, st%d%dim
          call X(oep_x) (namespace, gr%der, hm%psolver, st, is, jdm, oep%X(lxc), eig, xcs%cam_alpha)
        end do sum_comp
        ex = ex + eig
      end select
    end do functl_loop

    ! SIC a la PZ is handled here
    if(apply_sic_pz) then
      call X(oep_sic) (xcs, gr, hm%psolver, namespace, st, is, oep, ex, ec)
    end if
    ! calculate uxc_bar for the occupied states

    SAFE_ALLOCATE(psi(1:gr%mesh%np))

    oep%uxc_bar(:, is) = M_ZERO
    do ist = st%st_start, st%st_end
      call states_elec_get_state(st, gr%mesh, idm, ist, isp, psi)
      oep%uxc_bar(ist, is) = R_REAL(X(mf_dotp)(gr%mesh, R_CONJ(psi), oep%X(lxc)(1:gr%mesh%np, ist, is), reduce = .false.))
    end do
    if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, oep%uxc_bar(1:st%st_end, is), dim = st%st_end)

    SAFE_DEALLOCATE_A(psi)

  end do spin

#if defined(HAVE_MPI) 
  if(st%parallel_in_states) then
    call MPI_Barrier(st%mpi_grp%comm, mpi_err)
    do ist = 1, st%nst
      call MPI_Bcast(oep%uxc_bar(ist, isp), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
    end do
  end if
#endif

  if (st%d%ispin==SPINORS) then
    call xc_KLI_Pauli_solve(gr%mesh, namespace, st, oep)
    if(present(vxc)) then
      vxc(1:gr%mesh%np,:) = oep%vxc(1:gr%mesh%np,:)
    end if
    ! full OEP not implemented!
  else
    spin2: do is = 1, nspin_
      ! get the HOMO state
      call xc_oep_AnalyzeEigen(oep, st, is)
      !
      if(present(vxc)) then
        ! solve the KLI equation
        if(oep%level /= XC_OEP_FULL .or. first) then
          oep%vxc = M_ZERO
          call X(xc_KLI_solve) (namespace, gr%mesh, gr, hm, st, is, oep, first)
          if(present(vxc)) then
            vxc(1:gr%mesh%np, is) = vxc(1:gr%mesh%np, is) + oep%vxc(1:gr%mesh%np, 1)
          end if
        end if
        ! if asked, solve the full OEP equation
        if(oep%level == XC_OEP_FULL .and. (.not. first)) then
          call X(xc_oep_solve)(namespace, gr, hm, st, is, vxc(:,is), oep)
          if(present(vxc)) then
            vxc(1:gr%mesh%np, is) = vxc(1:gr%mesh%np, is) + oep%vxc(1:gr%mesh%np, is)
          end if
        end if
        if (is == nspin_) first = .false.
      end if
    end do spin2
  end if
  SAFE_DEALLOCATE_P(oep%eigen_type)
  SAFE_DEALLOCATE_P(oep%eigen_index)
  SAFE_DEALLOCATE_P(oep%X(lxc))
  SAFE_DEALLOCATE_P(oep%uxc_bar)

  POP_SUB(X(xc_oep_calc))
  call profiling_out(C_PROFILING_XC_OEP)
end subroutine X(xc_OEP_calc)


! ---------------------------------------------------------
subroutine X(xc_oep_solve) (namespace, gr, hm, st, is, vxc, oep)
  type(namespace_t),        intent(in)    :: namespace
  type(grid_t),             intent(in)    :: gr
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  FLOAT,                    intent(inout) :: vxc(:) !< (gr%mesh%np, given for the spin is)
  type(xc_oep_t),           intent(inout) :: oep

  integer :: iter, ist, iter_used
  FLOAT :: vxc_bar, ff, residue
  FLOAT, allocatable :: ss(:), vxc_old(:)
  R_TYPE, allocatable :: bb(:,:), psi(:, :), psi2(:,:)
  R_TYPE, allocatable :: phi1(:,:,:)
  logical, allocatable :: orthogonal(:)
  
  call profiling_in(C_PROFILING_XC_OEP_FULL, 'XC_OEP_FULL')
  PUSH_SUB(X(xc_oep_solve))

  if(st%parallel_in_states) &
    call messages_not_implemented("Full OEP parallel in states", namespace=namespace)

  SAFE_ALLOCATE(     bb(1:gr%mesh%np, 1:1))
  SAFE_ALLOCATE(     ss(1:gr%mesh%np))
  SAFE_ALLOCATE(vxc_old(1:gr%mesh%np))
  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi2(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(orthogonal(1:oep%noccst))

  if (oep%has_photons) then
    SAFE_ALLOCATE(phi1(1:gr%mesh%np, 1:st%d%dim, 1:oep%noccst))
  end if

  call lalg_copy(gr%mesh%np, vxc, vxc_old)

  if(.not. lr_is_allocated(oep%lr)) then
    call lr_allocate(oep%lr, st, gr%mesh)
    oep%lr%X(dl_psi)(:,:, :, :) = M_ZERO
  end if

  if (oep%has_photons) then
    if(.not. lr_is_allocated(oep%photon_lr)) then
      call lr_allocate(oep%photon_lr, st, gr%mesh)
      oep%photon_lr%X(dl_psi)(:, :, :, :) = M_ZERO
    end if
    call X(xc_oep_pt_phi)(namespace, gr, hm, st, is, oep, phi1)
  end if

  ! fix xc potential (needed for Hpsi)
  call lalg_axpy(gr%mesh%np, M_ONE, vxc_old(:), oep%vxc(:, is))

  do iter = 1, oep%scftol%max_iter
    ! iteration over all states
    ss = M_ZERO
    do ist = 1, st%nst !only over occupied states

      if(abs(st%occ(ist,1))<= M_EPSILON) cycle
      call states_elec_get_state(st, gr%mesh, ist, is, psi)
      psi2(:, 1) = R_CONJ(psi(:, 1))*psi(:,1)

      ! evaluate right-hand side
      vxc_bar = X(mf_integrate)(gr%mesh, psi2(:, 1)*oep%vxc(1:gr%mesh%np, is))

      bb(1:gr%mesh%np, 1) = -(oep%vxc(1:gr%mesh%np, is) - (vxc_bar - oep%uxc_bar(ist, is)))* &
        R_CONJ(psi(:, 1)) + oep%X(lxc)(1:gr%mesh%np, ist, is)

      if (oep%has_photons) call X(xc_oep_pt_rhs)(gr, st, is, oep, phi1, ist, bb)

      if (oep%has_photons) then
        orthogonal = .true.
        orthogonal(ist) = .false.
        call X(states_elec_orthogonalize_single)(st, gr%mesh, st%nst, is, bb, normalize = .false., mask = orthogonal)
      else
        call X(lr_orth_vector) (gr%mesh, st, bb, ist, is, R_TOTYPE(M_ZERO))
      end if

      call X(linear_solver_solve_HXeY)(oep%solver, namespace, hm, gr, st, ist, is, oep%lr%X(dl_psi)(:,:, ist, is), bb, &
           R_TOTYPE(-st%eigenval(ist, is)), oep%scftol%final_tol, residue, iter_used)

      if (oep%has_photons) then
        orthogonal = .true.
        orthogonal(ist) = .false.
        call X(states_elec_orthogonalize_single)(st, gr%mesh, st%nst, is, &
        oep%lr%X(dl_psi)(:,:, ist, is), normalize = .false., mask = orthogonal)
      else
        call X(lr_orth_vector) (gr%mesh, st, oep%lr%X(dl_psi)(:,:, ist, is), ist, is, R_TOTYPE(M_ZERO))
      end if

      ! calculate this funny function ss
      ! ss = ss + 2*dl_psi*psi
      call lalg_axpy(gr%mesh%np, M_TWO, R_REAL(oep%lr%X(dl_psi)(1:gr%mesh%np, 1, ist, is)*psi(:, 1)), ss(:))
      if (oep%has_photons) then
        call X(xc_oep_pt_inhomog)(gr, st, is, phi1, ist, ss)
      end if
    end do

    select case (oep%mixing_scheme)
    case (OEP_MIXING_SCHEME_CONST)
      call lalg_axpy(gr%mesh%np, oep%mixing, ss(:), oep%vxc(:, is))
    case (OEP_MIXING_SCHEME_DENS)
      call lalg_axpy(gr%mesh%np, oep%mixing, ss(:)/st%rho(:,is), oep%vxc(:, is))
    case (OEP_MIXING_SCHEME_BB)
      if (dmf_nrm2(gr%mesh, oep%vxc_old(1:gr%mesh%np,is)) > M_EPSILON ) then ! do not do it for the first run
        oep%mixing = -dmf_dotp(gr%mesh, oep%vxc(1:gr%mesh%np,is) - oep%vxc_old(1:gr%mesh%np,is), ss - oep%ss_old(:, is)) &
          / dmf_dotp(gr%mesh, ss - oep%ss_old(:, is), ss - oep%ss_old(:, is))
      end if

      if(debug%info) then
        write(message(1), '(a,es14.6,a,es14.8)') "Info: oep%mixing:", oep%mixing, " norm2ss: ", dmf_nrm2(gr%mesh, ss)
       call messages_info(1)
      end if

      call lalg_copy(gr%mesh%np, is, oep%vxc, oep%vxc_old)
      call lalg_copy(gr%mesh%np, ss, oep%ss_old(:, is))
      call lalg_axpy(gr%mesh%np, oep%mixing, ss(:), oep%vxc(:, is))
    end select

    do ist = 1, st%nst
      if(oep%eigen_type(ist) == 2) then
        call states_elec_get_state(st, gr%mesh, ist, is, psi)
        psi2(:, 1) = R_CONJ(psi(:, 1))*psi(:,1)
        vxc_bar = X(mf_integrate)(gr%mesh, psi2(:, 1)*oep%vxc(1:gr%mesh%np, is))
        if (oep%has_photons) then
          call X(xc_oep_pt_uxcbar)(gr, st, is, oep, phi1, ist, vxc_bar)
	end if
        oep%vxc(1:gr%mesh%np,is) = oep%vxc(1:gr%mesh%np,is) - (vxc_bar - oep%uxc_bar(ist,is))
      end if
    end do

    ff = dmf_nrm2(gr%mesh, ss)
    if(ff < oep%scftol%conv_abs_dens) exit
  end do

  if (is == 1) then
    oep%norm2ss = ff
  else
    oep%norm2ss = oep%norm2ss + ff !adding up spin up and spin down component
  end if

  if(ff > oep%scftol%conv_abs_dens) then
    write(message(1), '(a)') "OEP did not converge."
    call messages_warning(1, namespace=namespace)

    ! otherwise the number below will be one too high
    iter = iter - 1
  end if

  write(message(1), '(a,i4,a,es14.6)') "Info: After ", iter, " iterations, the OEP residual = ", ff
  message(2) = ''
  call messages_info(2)

  call lalg_copy(gr%mesh%np, vxc_old, vxc)

  SAFE_DEALLOCATE_A(bb)
  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(vxc_old)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(orthogonal)
  if (oep%has_photons) then
    SAFE_DEALLOCATE_A(phi1)
  end if

  POP_SUB(X(xc_oep_solve))
  call profiling_out(C_PROFILING_XC_OEP_FULL)
end subroutine X(xc_oep_solve)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
