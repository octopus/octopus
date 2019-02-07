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
subroutine X(hamiltonian_apply_batch) (hm, der, psib, hpsib, ik, terms, set_bc, set_phase)
  type(hamiltonian_t),   intent(in)    :: hm
  type(derivatives_t),   intent(in)    :: der
  type(batch_t), target, intent(inout) :: psib
  type(batch_t), target, intent(inout) :: hpsib
  integer,               intent(in)    :: ik
  integer, optional,     intent(in)    :: terms
  logical, optional,     intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.
  logical, optional,     intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  logical :: apply_phase, pack, set_phase_
  type(batch_t), pointer :: epsib
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection
  logical :: copy_at_end
  
  call profiling_in(prof_hamiltonian, "HAMILTONIAN")
  PUSH_SUB(X(hamiltonian_apply_batch))

  ASSERT(batch_status(psib) == batch_status(hpsib))

  ! all terms are enabled by default
  terms_ = optional_default(terms, TERM_ALL)
  set_phase_ = optional_default(set_phase, .true.)
  ! OpenCL is not supported for the phase correction at the moment
  if(.not.set_phase_) then
    if(batch_status(psib) == BATCH_DEVICE_PACKED) set_phase_ = .true.
  end if

  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)
  ASSERT(ik >= hm%d%kpt%start .and. ik <= hm%d%kpt%end)

  apply_phase = associated(hm%hm_base%phase)

  pack = hamiltonian_apply_packed(hm, der%mesh) &
    .and. (accel_is_enabled() .or. psib%nst_linear > 1) &
    .and. terms_ == TERM_ALL

  copy_at_end = .false.
  if(pack) then
    ! unpack at end with copying only if the status on entry is unpacked
    copy_at_end = .not. batch_is_packed(psib)
    call batch_pack(psib)
    call batch_pack(hpsib, copy = .false.)
  end if

  if(optional_default(set_bc, .true.)) then
    if(apply_phase .and. .not.set_phase_) then
      ! apply phase correction while setting boundary -> memory needs to be
      ! accessed only once
      call boundaries_set(der%boundaries, psib, phase_correction=hm%hm_base%phase_corr(:, ik))
    else
      call boundaries_set(der%boundaries, psib)
    end if
  else
    ! This should never happen: the phase correction should be always only
    ! applied when also setting the boundary conditions.
    ASSERT(.not.(apply_phase .and. .not.set_phase_))
  end if


  if(apply_phase .and. set_phase_) then
    SAFE_ALLOCATE(epsib)
    call batch_copy(psib, epsib, fill_zeros = .false.)
  else
    epsib => psib
  end if

  if(apply_phase .and. set_phase_) then ! we copy psi to epsi applying the exp(i k.r) phase
    call X(hamiltonian_base_phase)(hm%hm_base, der, der%mesh%np_part, ik, .false., epsib, src = psib)
  end if

  if(bitand(TERM_KINETIC, terms_) /= 0) then
    ASSERT(associated(hm%hm_base%kinetic))
    call profiling_in(prof_kinetic_start, "KINETIC_START")
    call X(derivatives_batch_start)(hm%hm_base%kinetic, der, epsib, hpsib, handle, set_bc = .false., factor = -M_HALF)
    call profiling_out(prof_kinetic_start)
  end if

  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_start)(hm%hm_base, der%mesh, hm%d, ik, epsib, projection)
    end if
  end if

  if(bitand(TERM_KINETIC, terms_) /= 0) then
    call profiling_in(prof_kinetic_finish, "KINETIC_FINISH")
    call X(derivatives_batch_finish)(handle)
    call profiling_out(prof_kinetic_finish)
  else
    call batch_set_zero(hpsib)
  end if

  ! apply the local potential
  if (bitand(TERM_LOCAL_POTENTIAL, terms_) /= 0) then
    call X(hamiltonian_base_local)(hm%hm_base, der%mesh, hm%d, states_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  else if(bitand(TERM_LOCAL_EXTERNAL, terms_) /= 0) then
    call X(hamiltonian_external)(hm, der%mesh, epsib, hpsib)
  end if

  ! and the non-local one
  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_finish)(hm%hm_base, der%mesh, hm%d, ik, projection, hpsib)
    else
      call X(project_psi_batch)(der%mesh, hm%ep%proj, hm%ep%natoms, hm%d%dim, epsib, hpsib, ik)
    end if
  end if

  if (bitand(TERM_OTHERS, terms_) /= 0) then

    call profiling_in(prof_exx, "EXCHANGE_OPERATOR")
    select case(hm%theory_level)

    case(HARTREE_FOCK)
      call X(exchange_operator)(hm, der, ik, hm%exx_coef, epsib, hpsib)

    end select
    call profiling_out(prof_exx)
    
  end if

  if(apply_phase .and. set_phase_) then
    call X(hamiltonian_base_phase)(hm%hm_base, der, der%mesh%np, ik, .true., hpsib)
    call batch_end(epsib, copy = .false.)
    SAFE_DEALLOCATE_P(epsib)
  end if

  if(pack) then
    call batch_unpack(psib, copy = .false.)
    call batch_unpack(hpsib, copy=copy_at_end)
  end if

  POP_SUB(X(hamiltonian_apply_batch))
  call profiling_out(prof_hamiltonian)

end subroutine X(hamiltonian_apply_batch)

! ---------------------------------------------------------

subroutine X(hamiltonian_external)(this, mesh, psib, vpsib)
  type(hamiltonian_t),         intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  FLOAT, dimension(:), pointer :: vpsl
  FLOAT, allocatable :: vpsl_spin(:,:), Imvpsl_spin(:,:)
  integer :: pnp, offset, ispin
  type(accel_mem_t) :: vpsl_buff

  PUSH_SUB(X(hamiltonian_external))

  SAFE_ALLOCATE(vpsl_spin(1:mesh%np, 1:this%d%nspin))

  nullify(vpsl)
  ! Sets the vpsl pointer to the total potential.
  vpsl => this%ep%vpsl

  vpsl_spin(1:mesh%np, 1) = vpsl(1:mesh%np)
  if(this%d%ispin == SPINORS) then
    ! yes this means a little unnecessary computation in the later call,
    ! but with the great benefit of being able to reuse an existing routine
    vpsl_spin(1:mesh%np, 2) = vpsl(1:mesh%np)
    vpsl_spin(1:mesh%np, 3) = M_ZERO
    vpsl_spin(1:mesh%np, 4) = M_ZERO
  end if

  if(batch_status(psib) == BATCH_DEVICE_PACKED) then
    pnp = accel_padded_size(mesh%np)
    call accel_create_buffer(vpsl_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, pnp * this%d%nspin)
    call accel_write_buffer(vpsl_buff, mesh%np, vpsl)

    offset = 0
    do ispin = 1, this%d%nspin
       call accel_write_buffer(vpsl_buff, mesh%np, vpsl_spin(:, ispin), offset = offset)
       offset = offset + pnp
    end do

    call X(hamiltonian_base_local_sub)(vpsl_spin, mesh, this%d, 1, &
      psib, vpsib, potential_opencl = vpsl_buff)

    call accel_release_buffer(vpsl_buff)
  else
    call X(hamiltonian_base_local_sub)(vpsl_spin, mesh, this%d, 1, psib, vpsib)
  end if

  SAFE_DEALLOCATE_A(vpsl_spin)

  POP_SUB(X(hamiltonian_external))
end subroutine X(hamiltonian_external)

! ---------------------------------------------------------

subroutine X(hamiltonian_apply) (hm, der, psi, hpsi, ist, ik, terms, set_bc, set_phase)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ist       !< the index of the state
  integer,             intent(in)    :: ik        !< the index of the k-point
  R_TYPE,   target,    intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
  R_TYPE,   target,    intent(inout) :: hpsi(:,:) !< (gr%mesh%np, hm%d%dim)
  integer,  optional,  intent(in)    :: terms
  logical,  optional,  intent(in)    :: set_bc
  logical, optional,     intent(in)    :: set_phase

  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_apply))

  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, terms = terms, set_bc = set_bc, &
                                     set_phase = set_phase)

  call batch_end(psib)
  call batch_end(hpsib)

  POP_SUB(X(hamiltonian_apply))
end subroutine X(hamiltonian_apply)

! ---------------------------------------------------------

subroutine X(exchange_operator)(hm, der, ik, exx_coef, psib, hpsib)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ik
  FLOAT,               intent(in)    :: exx_coef
  type(batch_t),       intent(inout) :: psib
  type(batch_t),       intent(inout) :: hpsib

  integer :: ibatch, jst, ip, idim, ik2, ib, ii, ist
  type(batch_t), pointer :: psi2b
  FLOAT                              :: ff
  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :), psi(:, :), hpsi(:, :)

  PUSH_SUB(X(exchange_operator))

  ASSERT(associated(hm%hf_st))

  if(der%mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points")

  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:hm%d%dim))

  do ibatch = 1, psib%nst
    ist = psib%states(ibatch)%ist
    call batch_get_state(psib, ibatch, der%mesh%np, psi)
    call batch_get_state(hpsib, ibatch, der%mesh%np, hpsi)

    do ik2 = 1, hm%d%nik
      if(states_dim_get_spin_index(hm%d, ik2) /= states_dim_get_spin_index(hm%d, ik)) cycle

      do ib = 1, hm%hf_st%group%nblocks

        call states_parallel_get_block(hm%hf_st, der%mesh, ib, ik2, psi2b)

        do ii = 1, psi2b%nst

          jst = psi2b%states(ii)%ist

          if(hm%hf_st%occ(jst, ik2) < M_EPSILON) cycle

          pot = R_TOTYPE(M_ZERO)
          rho = R_TOTYPE(M_ZERO)

          call batch_get_state(psi2b, ii, der%mesh%np, psi2)

          do idim = 1, hm%hf_st%d%dim
            forall(ip = 1:der%mesh%np)
              rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
            end forall
          end do

          call X(poisson_solve)(psolver, pot, rho, all_nodes = .false.)

          ff = hm%hf_st%occ(jst, ik2)
          if(hm%d%ispin == UNPOLARIZED) ff = M_HALF*ff

          do idim = 1, hm%hf_st%d%dim
            forall(ip = 1:der%mesh%np)
              hpsi(ip, idim) = hpsi(ip, idim) - exx_coef*ff*psi2(ip, idim)*pot(ip)
            end forall
          end do

        end do

        call states_parallel_release_block(hm%hf_st, ib, ik2, psi2b)

      end do
    end do
    call batch_set_state(hpsib, ibatch, der%mesh%np, hpsi)

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator))
end subroutine X(exchange_operator)

! ---------------------------------------------------------
subroutine X(vmask) (gr, hm, st)
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim
  R_TYPE, allocatable :: psi(:)

  PUSH_SUB(X(vmask))

  SAFE_ALLOCATE(psi(1:gr%mesh%np))

  if(hm%bc%abtype == MASK_ABSORBING) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_get_state(st, gr%mesh, idim, ist, ik, psi)
          psi(1:gr%mesh%np) = psi(1:gr%mesh%np)*hm%bc%mf(1:gr%mesh%np)
          call states_set_state(st, gr%mesh, idim, ist, ik, psi)
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(vmask))
end subroutine X(vmask)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
