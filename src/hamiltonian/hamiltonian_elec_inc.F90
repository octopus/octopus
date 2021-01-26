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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
subroutine X(hamiltonian_elec_apply) (hm, namespace, mesh, psib, hpsib, terms, set_bc)
  class(hamiltonian_elec_t),   intent(in)    :: hm
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  class(batch_t),      target, intent(inout) :: psib
  class(batch_t),      target, intent(inout) :: hpsib
  integer,           optional, intent(in)    :: terms
  logical,           optional, intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.

  PUSH_SUB(X(hamiltonian_elec_apply))

  select type (psib)
  class is (wfs_elec_t)
    select type (hpsib)
    class is (wfs_elec_t)
      call X(hamiltonian_elec_apply_batch) (hm, namespace, mesh, psib, hpsib, terms, set_bc)
    class default
      message(1) = "Internal error: imcompatible batch_t in hamiltonian_elec_apply for argument hpsib."
      call messages_fatal(1)
    end select
  class default
    message(1) = "Internal error: imcompatible batch_t in hamiltonian_elec_apply for argument psib."
    call messages_fatal(1)
  end select

  POP_SUB(X(hamiltonian_elec_apply))
end subroutine X(hamiltonian_elec_apply)

! ---------------------------------------------------------
subroutine X(hamiltonian_elec_magnus_apply) (hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
  class(hamiltonian_elec_t),   intent(in)    :: hm
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  class(batch_t),              intent(inout) :: psib
  class(batch_t),              intent(inout) :: hpsib
  FLOAT,                       intent(in)    :: vmagnus(:, :, :)
  logical,           optional, intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  PUSH_SUB(X(hamiltonian_elec_magnus_apply))

  select type (psib)
  class is (wfs_elec_t)
    select type (hpsib)
    class is (wfs_elec_t)
      call X(hamiltonian_elec_magnus_apply_batch) (hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
    class default
      message(1) = "Internal error: imcompatible batch_t in hamiltonian_elec_magnus_apply for argument hpsib."
      call messages_fatal(1)
    end select
  class default
    message(1) = "Internal error: imcompatible batch_t in hamiltonian_elec_magnus_apply for argument psib."
    call messages_fatal(1)
  end select

  POP_SUB(X(hamiltonian_elec_magnus_apply))
end subroutine X(hamiltonian_elec_magnus_apply)

! ---------------------------------------------------------
subroutine X(hamiltonian_elec_apply_batch) (hm, namespace, mesh, psib, hpsib, terms, set_bc)
  type(hamiltonian_elec_t),    intent(in)    :: hm
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(wfs_elec_t),    target, intent(inout) :: psib
  type(wfs_elec_t),    target, intent(inout) :: hpsib
  integer,           optional, intent(in)    :: terms
  logical,           optional, intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.

  logical :: apply_phase, pack, set_phase
  type(wfs_elec_t), pointer :: epsib
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection
  
  call profiling_in(prof_hamiltonian, TOSTRING(X(HAMILTONIAN)))
  PUSH_SUB(X(hamiltonian_elec_apply_batch))

  ASSERT(psib%status() == hpsib%status())

  ! all terms are enabled by default
  terms_ = optional_default(terms, TERM_ALL)

  !By default the phase is dictated by apply_phase
  !and for electronic wavefunctions, we get set_phase from the has_phase flag
  set_phase = .not. psib%has_phase

  ASSERT(psib%nst == hpsib%nst)
  ASSERT(psib%ik >= hm%d%kpt%start .and. psib%ik <= hm%d%kpt%end)

  apply_phase = associated(hm%hm_base%phase)

  pack = hamiltonian_elec_apply_packed(hm) &
    .and. (accel_is_enabled() .or. psib%nst_linear > 1) &
    .and. terms_ == TERM_ALL

  if(pack) then
    call psib%do_pack()
    call hpsib%do_pack(copy = .false.)
  end if

  if(optional_default(set_bc, .true.)) then
    if(apply_phase .and. .not.set_phase) then
      ! apply phase correction while setting boundary -> memory needs to be
      ! accessed only once
      ASSERT(psib%has_phase)
      call boundaries_set(hm%der%boundaries, psib, phase_correction = hm%hm_base%phase_corr(:, psib%ik))
    else
      call boundaries_set(hm%der%boundaries, psib)
    end if
  else
    ! This should never happen: the phase correction should be always only
    ! applied when also setting the boundary conditions.
    ASSERT(.not.(apply_phase .and. .not. set_phase))
  end if


  if(apply_phase .and. set_phase) then
    SAFE_ALLOCATE(epsib)
    call psib%copy_to(epsib)
  else
    epsib => psib
  end if

  if(apply_phase .and. set_phase) then ! we copy psi to epsi applying the exp(i k.r) phase
    call hamiltonian_elec_base_phase(hm%hm_base, mesh, mesh%np_part, .false., epsib, src = psib)
    hpsib%has_phase = .true.
  end if

  !Apply the spiral BC if needed
  if(hm%der%boundaries%spiral .and. apply_phase) then
    call hamiltonian_elec_base_phase_spiral(hm%hm_base, hm%der, epsib)
  end if

  if(bitand(TERM_KINETIC, terms_) /= 0) then
    ASSERT(associated(hm%hm_base%kinetic))
    call profiling_in(prof_kinetic_start, TOSTRING(X(KINETIC_START)))
    call X(derivatives_batch_start)(hm%hm_base%kinetic, hm%der, epsib, hpsib, handle, set_bc = .false., factor = -M_HALF/hm%mass)
    call profiling_out(prof_kinetic_start)
  end if

  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_elec_base_nlocal_start)(hm%hm_base, mesh, hm%d, hm%der%boundaries, epsib, projection)
    end if
  end if

  if(bitand(TERM_KINETIC, terms_) /= 0) then
    call profiling_in(prof_kinetic_finish, TOSTRING(X(KINETIC_FINISH)))
    call X(derivatives_batch_finish)(handle)
    call profiling_out(prof_kinetic_finish)
  else
    call batch_set_zero(hpsib)
  end if

  ! apply the local potential
  if (bitand(TERM_LOCAL_POTENTIAL, terms_) /= 0) then
    call X(hamiltonian_elec_base_local)(hm%hm_base, mesh, hm%d, states_elec_dim_get_spin_index(hm%d, psib%ik), epsib, hpsib)
  else if(bitand(TERM_LOCAL_EXTERNAL, terms_) /= 0) then
    call X(hamiltonian_elec_external)(hm, mesh, epsib, hpsib)
  end if

  ! and the non-local one
  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_elec_base_nlocal_finish)(hm%hm_base, mesh, hm%der%boundaries, hm%d, projection, hpsib)
    else
      call X(project_psi_batch)(mesh, hm%der%boundaries, hm%ep%proj, hm%ep%natoms, hm%d%dim, epsib, hpsib)
    end if
  end if
  
  if (bitand(TERM_OTHERS, terms_) /= 0 .and. hamiltonian_elec_base_has_magnetic(hm%hm_base)) then
    call X(hamiltonian_elec_base_magnetic)(hm%hm_base, mesh, hm%der, hm%d, hm%ep, &
             states_elec_dim_get_spin_index(hm%d, psib%ik), epsib, hpsib)
  end if
  
  if (bitand(TERM_OTHERS, terms_) /= 0 ) then
    call hamiltonian_elec_base_rashba(hm%hm_base, mesh, hm%der, hm%d, epsib, hpsib)
  end if

  ! multiply with occupation number
  if (hm%theory_level == RDMFT .and. bitand(TERM_RDMFT_OCC, terms_) /= 0) then
    call exchange_operator_rdmft_occ_apply(hm%exxop, mesh, hpsib)
  endif
  
  if (bitand(TERM_OTHERS, terms_) /= 0) then

    call profiling_in(prof_exx, TOSTRING(X(EXCHANGE_OPERATOR)))
    select case(hm%theory_level)

    case(HARTREE)
      call X(exchange_operator_hartree_apply)(hm%exxop, namespace, hm%der, hm%d, hm%exxop%cam_alpha, epsib, hpsib)

    case(HARTREE_FOCK)
      if(hm%scdm_EXX)  then
        call X(exchange_operator_scdm_apply)(hm%exxop, namespace, hm%scdm, hm%der, hm%d, epsib, hpsib, hm%exxop%cam_alpha, &
                          hm%theory_level == HARTREE)
      else
        ! standard HF 
        call X(exchange_operator_apply)(hm%exxop, namespace, hm%der, hm%d, epsib, hpsib, .false.)
      end if

    case(RDMFT)
      call X(exchange_operator_apply)(hm%exxop, namespace, hm%der, hm%d, epsib, hpsib, .true.)
    end select
    call profiling_out(prof_exx)
    
  end if

  if (bitand(TERM_MGGA, terms_) /= 0 .and. family_is_mgga_with_exc(hm%xc)) then
    call X(h_mgga_terms)(hm, mesh, epsib, hpsib)
  end if

  if(bitand(TERM_OTHERS, terms_) /= 0 .and. hm%scissor%apply) then
    call X(scissor_apply)(hm%scissor, mesh, epsib, hpsib)
  end if

  if(iand(TERM_DFT_U, terms_) /= 0 .and. hm%lda_u_level /= DFT_U_NONE) then
    call X(lda_u_apply)(hm%lda_u, hm%d, mesh, epsib, hpsib)
  end if  

  if(apply_phase .and. set_phase) then
    call hamiltonian_elec_base_phase(hm%hm_base, mesh, mesh%np, .true., hpsib)
    call epsib%end(copy = .false.)
    SAFE_DEALLOCATE_P(epsib)
  end if

  if(pack) then
    call psib%do_unpack(copy = .false.)
    call hpsib%do_unpack()
  end if

  POP_SUB(X(hamiltonian_elec_apply_batch))
  call profiling_out(prof_hamiltonian)
end subroutine X(hamiltonian_elec_apply_batch)

! ---------------------------------------------------------

subroutine X(hamiltonian_elec_external)(this, mesh, psib, vpsib)
  type(hamiltonian_elec_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(wfs_elec_t),            intent(in)    :: psib
  type(wfs_elec_t),            intent(inout) :: vpsib

  FLOAT, allocatable :: vpsl_spin(:,:)
  integer :: pnp, offset, ispin
  type(accel_mem_t) :: vpsl_buff

  PUSH_SUB(X(hamiltonian_elec_external))

  SAFE_ALLOCATE(vpsl_spin(1:mesh%np, 1:this%d%nspin))

  vpsl_spin(1:mesh%np, 1) = this%ep%vpsl(1:mesh%np)
  if(this%d%ispin == SPINORS) then
    ! yes this means a little unnecessary computation in the later call,
    ! but with the great benefit of being able to reuse an existing routine
    vpsl_spin(1:mesh%np, 2) = this%ep%vpsl(1:mesh%np)
    vpsl_spin(1:mesh%np, 3) = M_ZERO
    vpsl_spin(1:mesh%np, 4) = M_ZERO
  end if

  if(psib%status() == BATCH_DEVICE_PACKED) then
    pnp = accel_padded_size(mesh%np)
    call accel_create_buffer(vpsl_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, pnp * this%d%nspin)
    call accel_write_buffer(vpsl_buff, mesh%np, this%ep%vpsl)

    offset = 0
    do ispin = 1, this%d%nspin
       call accel_write_buffer(vpsl_buff, mesh%np, vpsl_spin(:, ispin), offset = offset)
       offset = offset + pnp
    end do

    call X(hamiltonian_elec_base_local_sub)(vpsl_spin, mesh, this%d, 1, &
      psib, vpsib, potential_opencl = vpsl_buff)

    call accel_release_buffer(vpsl_buff)
  else
    call X(hamiltonian_elec_base_local_sub)(vpsl_spin, mesh, this%d, 1, psib, vpsib)
  end if

  SAFE_DEALLOCATE_A(vpsl_spin)

  POP_SUB(X(hamiltonian_elec_external))
end subroutine X(hamiltonian_elec_external)

! ---------------------------------------------------------

subroutine X(hamiltonian_elec_apply_single) (hm, namespace, mesh, psi, hpsi, ist, ik, terms, set_bc, set_phase)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ist       !< the index of the state
  integer,                  intent(in)    :: ik        !< the index of the k-point
  R_TYPE, contiguous, target, intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
  R_TYPE, contiguous, target, intent(inout) :: hpsi(:,:) !< (gr%mesh%np, hm%d%dim)
  integer, optional,        intent(in)    :: terms
  logical, optional,        intent(in)    :: set_bc
  logical, optional,        intent(in)    :: set_phase

  type(wfs_elec_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_elec_apply_single))
  
  call wfs_elec_init(psib, hm%d%dim, ist, ist, psi, ik)
  call wfs_elec_init(hpsib, hm%d%dim, ist, ist, hpsi, ik)

  if(present(set_phase)) then
    psib%has_phase = .not. set_phase
    hpsib%has_phase = .not. set_phase
  end if 

  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib, hpsib, terms = terms, set_bc = set_bc)

  call psib%end()
  call hpsib%end()
  

  POP_SUB(X(hamiltonian_elec_apply_single))
end subroutine X(hamiltonian_elec_apply_single)


subroutine X(hamiltonian_elec_magnus_apply_batch) (hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(wfs_elec_t),         intent(inout) :: psib
  type(wfs_elec_t),         intent(inout) :: hpsib
  FLOAT,                    intent(in)    :: vmagnus(:, :, :)
  logical, optional,        intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  integer :: ispin
  type(wfs_elec_t) :: auxpsib, aux2psib

  PUSH_SUB(X(hamiltonian_elec_magnus_apply_batch))

  ! We will assume, for the moment, no spinors.
  if (hm%d%dim /= 1) call messages_not_implemented("Magnus with spinors", namespace=namespace)

  ASSERT(psib%nst == hpsib%nst)

  ispin = states_elec_dim_get_spin_index(hm%d, psib%ik)

  call hpsib%copy_to(auxpsib, copy_data=.false.)
  call hpsib%copy_to(aux2psib, copy_data=.false.)
  
  ! Compute (T + Vnl)|psi> and store it
  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib, auxpsib, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL)

  ! H|psi>  =  (T + Vnl)|psi> + Vpsl|psi> + Vmagnus(t2)|psi> + Vborders|psi>
  call auxpsib%copy_data_to(mesh%np, hpsib)
  call X(hamiltonian_elec_external)(hm, mesh, psib, hpsib)
  if (hm%bc%abtype == IMAGINARY_ABSORBING) then
    call batch_mul(mesh%np, hm%bc%mf(1:mesh%np), psib, aux2psib)
    call batch_axpy(mesh%np, M_zI, aux2psib, hpsib)
  end if
  call batch_mul(mesh%np, vmagnus(1:mesh%np, ispin, 2), psib, aux2psib)
  call batch_axpy(mesh%np, M_ONE, aux2psib, hpsib)

  ! Add first term of the commutator:  - i Vmagnus(t1) (T + Vnl) |psi>
  call batch_mul(mesh%np, vmagnus(1:mesh%np, ispin, 1), auxpsib, aux2psib)
  call batch_axpy(mesh%np, -M_zI, aux2psib, hpsib)

  ! Add second term of commutator:  i (T + Vnl) Vmagnus(t1) |psi>
  call batch_mul(mesh%np, vmagnus(1:mesh%np, ispin, 1), psib, auxpsib)
  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, auxpsib, aux2psib, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL)
  call batch_axpy(mesh%np, M_zI, aux2psib, hpsib)

  call auxpsib%end()
  call aux2psib%end()

  POP_SUB(X(hamiltonian_elec_magnus_apply_batch))
end subroutine X(hamiltonian_elec_magnus_apply_batch)

! ---------------------------------------------------------
subroutine X(h_mgga_terms) (hm, mesh, psib, hpsib)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(mesh_t),             intent(in)    :: mesh
  type(wfs_elec_t),         intent(inout) :: psib
  type(wfs_elec_t),         intent(inout) :: hpsib

  integer :: ispin, ii, idir
  R_TYPE, allocatable :: grad(:,:), diverg(:)
  type(wfs_elec_t) :: divb
  class(wfs_elec_t), allocatable :: gradb(:)
  
  PUSH_SUB(X(h_mgga_terms))

  ispin = states_elec_dim_get_spin_index(hm%d, psib%ik)

  SAFE_ALLOCATE(grad(1:mesh%np_part, 1:mesh%sb%dim))
  SAFE_ALLOCATE(diverg(1:mesh%np))

  SAFE_ALLOCATE_TYPE_ARRAY(wfs_elec_t, gradb, (1:mesh%sb%dim))

  call hpsib%copy_to(divb)
  
  do idir = 1, mesh%sb%dim
    call hpsib%copy_to(gradb(idir))
  end do
  call X(derivatives_batch_grad)(hm%der, psib, gradb, ghost_update = .false., set_bc = .false.)
  
  do ii = 1, psib%nst_linear

    do idir = 1, mesh%sb%dim
      call batch_get_state(gradb(idir), ii, mesh%np, grad(:, idir))
    end do

    do idir = 1, mesh%sb%dim
      grad(1:mesh%np, idir) = grad(1:mesh%np, idir)*hm%vtau(1:mesh%np, ispin)
    end do

    call X(derivatives_div)(hm%der, grad, diverg)

    call batch_set_state(divb, ii, mesh%np, diverg)

  end do

  call batch_axpy(mesh%np, CNST(-1.0), divb, hpsib)

  do idir = 1, mesh%sb%dim
    call gradb(idir)%end()
  end do

  call divb%end()
  
  SAFE_DEALLOCATE_A(gradb)
  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(diverg)

  POP_SUB(X(h_mgga_terms))
end subroutine X(h_mgga_terms)


! ---------------------------------------------------------
subroutine X(vmask) (mesh, hm, st)
  type(mesh_t),        intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t), intent(inout) :: st

  integer :: ik, ist, idim
  R_TYPE, allocatable :: psi(:)

  PUSH_SUB(X(vmask))

  SAFE_ALLOCATE(psi(1:mesh%np))

  if(hm%bc%abtype == MASK_ABSORBING) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_elec_get_state(st, mesh, idim, ist, ik, psi)
          psi(1:mesh%np) = psi(1:mesh%np)*hm%bc%mf(1:mesh%np)
          call states_elec_set_state(st, mesh, idim, ist, ik, psi)
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(vmask))
end subroutine X(vmask)


! ---------------------------------------------------------
subroutine X(hamiltonian_elec_diagonal) (hm, mesh, diag, ik)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(mesh_t),             intent(in)    :: mesh
  R_TYPE,              intent(out)   :: diag(:,:) !< hpsi(gr%mesh%np, hm%d%dim)
  integer,             intent(in)    :: ik

  integer :: idim, ip, ispin

  FLOAT, allocatable  :: ldiag(:)

  PUSH_SUB(X(hamiltonian_elec_diagonal))

  SAFE_ALLOCATE(ldiag(1:mesh%np))

  diag = M_ZERO

  call derivatives_lapl_diag(hm%der, ldiag)

  do idim = 1, hm%d%dim
    diag(1:mesh%np, idim) = -M_HALF/hm%mass*ldiag(1:mesh%np)
  end do

  select case(hm%d%ispin)

  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_elec_dim_get_spin_index(hm%d, ik)
    diag(1:mesh%np, 1) = diag(1:mesh%np, 1) + hm%vhxc(1:mesh%np, ispin) + hm%ep%vpsl(1:mesh%np)

  case(SPINORS)
    do ip = 1, mesh%np
      diag(ip, 1) = diag(ip, 1) + hm%vhxc(ip, 1) + hm%ep%vpsl(ip)
      diag(ip, 2) = diag(ip, 2) + hm%vhxc(ip, 2) + hm%ep%vpsl(ip)
    end do

  end select

  POP_SUB(X(hamiltonian_elec_diagonal))
end subroutine X(hamiltonian_elec_diagonal)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
