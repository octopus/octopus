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
subroutine X(hamiltonian_elec_apply_batch) (hm, namespace, mesh, psib, hpsib, ik, terms, set_bc, set_phase)
  type(hamiltonian_elec_t),   intent(in)    :: hm
  type(namespace_t),          intent(in)    :: namespace
  type(mesh_t),               intent(in)    :: mesh
  type(batch_t),      target, intent(inout) :: psib
  type(batch_t),      target, intent(inout) :: hpsib
  integer,                    intent(in)    :: ik
  integer, optional,          intent(in)    :: terms
  logical, optional,          intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.
  logical, optional,          intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  logical :: apply_phase, pack, set_phase_
  type(batch_t), pointer :: epsib
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection
  
  call profiling_in(prof_hamiltonian, "HAMILTONIAN")
  PUSH_SUB(X(hamiltonian_elec_apply_batch))

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

  pack = hamiltonian_elec_apply_packed(hm) &
    .and. (accel_is_enabled() .or. psib%nst_linear > 1) &
    .and. terms_ == TERM_ALL

  if(pack) then
    call batch_pack(psib)
    call batch_pack(hpsib, copy = .false.)
  end if

  if(optional_default(set_bc, .true.)) then
    if(apply_phase .and. .not.set_phase_) then
      ! apply phase correction while setting boundary -> memory needs to be
      ! accessed only once
      call boundaries_set(hm%der%boundaries, psib, phase_correction=hm%hm_base%phase_corr(:, ik))
    else
      call boundaries_set(hm%der%boundaries, psib)
    end if
  else
    ! This should never happen: the phase correction should be always only
    ! applied when also setting the boundary conditions.
    ASSERT(.not.(apply_phase .and. .not.set_phase_))
  end if


  if(apply_phase .and. set_phase_) then
    SAFE_ALLOCATE(epsib)
    call batch_copy(psib, epsib)
  else
    epsib => psib
  end if

  if(apply_phase .and. set_phase_) then ! we copy psi to epsi applying the exp(i k.r) phase
    call X(hamiltonian_elec_base_phase)(hm%hm_base, mesh, mesh%np_part, ik, .false., epsib, src = psib)
  end if

  !Apply the spiral BC if needed
  if(hm%der%boundaries%spiral .and. apply_phase) then
    call X(hamiltonian_elec_base_phase_spiral)(hm%hm_base, hm%der, epsib, ik)
  end if

  if(bitand(TERM_KINETIC, terms_) /= 0) then
    ASSERT(associated(hm%hm_base%kinetic))
    call profiling_in(prof_kinetic_start, "KINETIC_START")
    call X(derivatives_batch_start)(hm%hm_base%kinetic, hm%der, epsib, hpsib, handle, set_bc = .false., factor = -M_HALF/hm%mass)
    call profiling_out(prof_kinetic_start)
  end if

  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_elec_base_nlocal_start)(hm%hm_base, mesh, hm%d, hm%der%boundaries, ik, epsib, projection)
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
    call X(hamiltonian_elec_base_local)(hm%hm_base, mesh, hm%d, states_elec_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  else if(bitand(TERM_LOCAL_EXTERNAL, terms_) /= 0) then
    call X(hamiltonian_elec_external)(hm, mesh, epsib, hpsib)
  end if

  ! and the non-local one
  if (hm%ep%non_local .and. bitand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_elec_base_nlocal_finish)(hm%hm_base, mesh, hm%der%boundaries, hm%d, ik, projection, hpsib)
    else
      call X(project_psi_batch)(mesh, hm%der%boundaries, hm%ep%proj, hm%ep%natoms, hm%d%dim, epsib, hpsib, ik)
    end if
  end if
  
  if (bitand(TERM_OTHERS, terms_) /= 0 .and. hamiltonian_elec_base_has_magnetic(hm%hm_base)) then
    call X(hamiltonian_elec_base_magnetic)(hm%hm_base, mesh, hm%der, hm%d, hm%ep, &
             states_elec_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  end if
  
  if (bitand(TERM_OTHERS, terms_) /= 0 ) then
    call X(hamiltonian_elec_base_rashba)(hm%hm_base, mesh, hm%der, hm%d, epsib, hpsib)
  end if

  ! multiply with occupation number
  if (hm%theory_level == RDMFT .and. bitand(TERM_RDMFT_OCC, terms_) /= 0) then
    call X(hamiltonian_elec_rdmft_occ_apply) (hm, mesh, ik, hpsib)
  endif
  
  if (bitand(TERM_OTHERS, terms_) /= 0) then

    call profiling_in(prof_exx, "EXCHANGE_OPERATOR")
    select case(hm%theory_level)

    case(HARTREE)
      call X(exchange_operator_hartree)(hm, namespace, mesh, ik, epsib, hpsib)

    case(HARTREE_FOCK)
      if(hm%scdm_EXX)  then
        call X(scdm_exchange_operator)(hm, namespace, mesh, epsib, hpsib, ik)
      else
        ! standard HF 
        call X(exchange_operator)(hm, namespace, mesh, ik, epsib, hpsib)
      end if

    case(RDMFT)
      call X(exchange_operator)(hm, namespace, mesh, ik, epsib, hpsib)
    end select
    call profiling_out(prof_exx)
    
  end if

  if (bitand(TERM_MGGA, terms_) /= 0 .and. family_is_mgga_with_exc(hm%xc)) then
    call X(h_mgga_terms)(hm, mesh, ik, epsib, hpsib)
  end if

  if(bitand(TERM_OTHERS, terms_) /= 0 .and. hm%scissor%apply) then
    call X(scissor_apply)(hm%scissor, mesh, ik, epsib, hpsib)
  end if

  if(iand(TERM_DFT_U, terms_) /= 0 .and. hm%lda_u_level /= DFT_U_NONE) then
    call X(lda_u_apply)(hm%lda_u, hm%d, mesh, mesh%sb, ik, epsib, hpsib, apply_phase)
  end if  

  if(apply_phase .and. set_phase_) then
    call X(hamiltonian_elec_base_phase)(hm%hm_base, mesh, mesh%np, ik, .true., hpsib)
    call batch_end(epsib, copy = .false.)
    SAFE_DEALLOCATE_P(epsib)
  end if

  if(pack) then
    call batch_unpack(psib, copy = .false.)
    call batch_unpack(hpsib)
  end if

  POP_SUB(X(hamiltonian_elec_apply_batch))
  call profiling_out(prof_hamiltonian)
end subroutine X(hamiltonian_elec_apply_batch)

! ---------------------------------------------------------

subroutine X(hamiltonian_elec_external)(this, mesh, psib, vpsib)
  type(hamiltonian_elec_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  FLOAT, dimension(:), pointer :: vpsl
  FLOAT, allocatable :: vpsl_spin(:,:)
  integer :: pnp, offset, ispin
  type(accel_mem_t) :: vpsl_buff

  PUSH_SUB(X(hamiltonian_elec_external))

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

subroutine X(hamiltonian_elec_apply) (hm, namespace, mesh, psi, hpsi, ist, ik, terms, set_bc, set_phase)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ist       !< the index of the state
  integer,                  intent(in)    :: ik        !< the index of the k-point
  R_TYPE,           target, intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
  R_TYPE,           target, intent(inout) :: hpsi(:,:) !< (gr%mesh%np, hm%d%dim)
  integer, optional,        intent(in)    :: terms
  logical, optional,        intent(in)    :: set_bc
  logical, optional,        intent(in)    :: set_phase

  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_elec_apply))
  
  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib, hpsib, ik, terms = terms, set_bc = set_bc, &
                                       set_phase = set_phase)

  call batch_end(psib)
  call batch_end(hpsib)
  

  POP_SUB(X(hamiltonian_elec_apply))
end subroutine X(hamiltonian_elec_apply)


subroutine X(hamiltonian_elec_rdmft_occ_apply) (hm, mesh, ik, hpsib)
  type(hamiltonian_elec_t), intent(in) :: hm
  type(mesh_t),             intent(in) :: mesh
  integer,             intent(in)    :: ik
  type(batch_t),       intent(inout) :: hpsib
  
  PUSH_SUB(X(hamiltonian_elec_rdmft_occ_apply))
  
  ! multiply linear terms in hamiltonian with occupation number
  ! nonlinear occupation number dependency occurs only in the exchange, which is treated there
  call batch_scal(mesh%np, hm%hf_st%occ(:, ik), hpsib)
  
  POP_SUB(X(hamiltonian_elec_rdmft_occ_apply))
end subroutine X(hamiltonian_elec_rdmft_occ_apply)


! ---------------------------------------------------------
subroutine X(hamiltonian_elec_apply_all) (hm, namespace, mesh, st, hst)
  type(hamiltonian_elec_t), intent(inout) :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  type(states_elec_t),      intent(inout) :: hst

  integer :: ik, ib, ist
  R_TYPE, allocatable :: psi(:, :)
  CMPLX,  allocatable :: psiall(:, :, :, :)
  
  PUSH_SUB(X(hamiltonian_elec_apply_all))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%group%block_start, st%group%block_end
      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), hst%group%psib(ib, ik), ik)
    end do
  end do

  if(oct_exchange_enabled(hm%oct_exchange)) then

    SAFE_ALLOCATE(psiall(mesh%np_part, 1:hst%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

    call states_elec_get_state(st, mesh, psiall)
    
    call oct_exchange_prepare(hm%oct_exchange, mesh, psiall, hm%xc, hm%psolver, namespace)

    SAFE_DEALLOCATE_A(psiall)
    
    SAFE_ALLOCATE(psi(mesh%np_part, 1:hst%d%dim))
    
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        call states_elec_get_state(hst, mesh, ist, ik, psi)
        call X(oct_exchange_operator)(hm%oct_exchange, namespace, mesh, psi, ist, ik)
        call states_elec_set_state(hst, mesh, ist, ik, psi)
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    
  end if

  POP_SUB(X(hamiltonian_elec_apply_all))
end subroutine X(hamiltonian_elec_apply_all)


! ---------------------------------------------------------

subroutine X(exchange_operator_single)(hm, namespace, mesh, ist, ik, psi, hpsi, exx_coef)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ist
  integer,                  intent(in)    :: ik
  R_TYPE,                   intent(inout) :: psi(:, :)
  R_TYPE,                   intent(inout) :: hpsi(:, :)
  FLOAT,          optional, intent(in)    :: exx_coef

  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(exchange_operator_single))

  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  call X(exchange_operator)(hm, namespace, mesh, ik, psib, hpsib, exx_coef)

  call batch_end(psib)
  call batch_end(hpsib)

  POP_SUB(X(exchange_operator_single))
end subroutine X(exchange_operator_single)

! ---------------------------------------------------------

subroutine X(exchange_operator)(hm, namespace, mesh, ik, psib, hpsib, exx_coef)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ik
  type(batch_t),            intent(inout) :: psib
  type(batch_t),            intent(inout) :: hpsib
  FLOAT,          optional, intent(in)    :: exx_coef

  integer :: ibatch, jst, ip, idim, ik2, ib, ii, ist
  type(batch_t), pointer :: psi2b
  FLOAT :: ff, exx_coef_
  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :), psi(:, :), hpsi(:, :)

  PUSH_SUB(X(exchange_operator))

  ASSERT(associated(hm%hf_st))

  if(mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)

  exx_coef_ = optional_default(exx_coef, hm%exx_coef)

  SAFE_ALLOCATE(psi(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(rho(1:mesh%np))
  SAFE_ALLOCATE(pot(1:mesh%np))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:hm%d%dim))

  do ibatch = 1, psib%nst
    ist = psib%states(ibatch)%ist
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib, ibatch, mesh%np, hpsi)

    do ik2 = 1, hm%d%nik
      if(states_elec_dim_get_spin_index(hm%d, ik2) /= states_elec_dim_get_spin_index(hm%d, ik)) cycle

      do ib = 1, hm%hf_st%group%nblocks

        call states_elec_parallel_get_block(hm%hf_st, mesh, ib, ik2, psi2b)

        do ii = 1, psi2b%nst

          jst = psi2b%states(ii)%ist

          if(hm%hf_st%occ(jst, ik2) < M_EPSILON) cycle

          pot = R_TOTYPE(M_ZERO)
          rho = R_TOTYPE(M_ZERO)

          call batch_get_state(psi2b, ii, mesh%np, psi2)

          do idim = 1, hm%hf_st%d%dim
            forall(ip = 1:mesh%np)
              rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
            end forall
          end do

          call X(poisson_solve)(hm%psolver, pot, rho, all_nodes = .false.)
          
          if (hm%theory_level /= RDMFT ) then
            ff = hm%hf_st%occ(jst, ik2)
            if(hm%d%ispin == UNPOLARIZED) ff = M_HALF*ff
          else ! RDMFT
            ff = sqrt(hm%hf_st%occ(ist, ik)*hm%hf_st%occ(jst, ik2)) ! Mueller functional
          end if

          do idim = 1, hm%hf_st%d%dim
            forall(ip = 1:mesh%np)
              hpsi(ip, idim) = hpsi(ip, idim) - exx_coef_*ff*psi2(ip, idim)*pot(ip)
            end forall
          end do

        end do

        call states_elec_parallel_release_block(hm%hf_st, ib, ik2, psi2b)

      end do
    end do
    call batch_set_state(hpsib, ibatch, mesh%np, hpsi)

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator))
end subroutine X(exchange_operator)

! ---------------------------------------------------------

subroutine X(exchange_operator_hartree) (hm, namespace, mesh, ik, psib, hpsib)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ik
  type(batch_t),            intent(inout) :: psib
  type(batch_t),            intent(inout) :: hpsib

  integer :: ibatch, ip, idim, ik2, ist
  FLOAT   :: ff
  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :), psi(:, :), hpsi(:, :)

  PUSH_SUB(X(exchange_operator))

  if(mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)
  if(hm%hf_st%parallel_in_states) call messages_not_implemented("exchange operator parallel in states", namespace=namespace)

  SAFE_ALLOCATE(psi(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(rho(1:mesh%np))
  SAFE_ALLOCATE(pot(1:mesh%np))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:hm%d%dim))

  do ibatch = 1, psib%nst
    ist = psib%states(ibatch)%ist
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib, ibatch, mesh%np, hpsi)
    
    do ik2 = 1, hm%d%nik
      if(states_elec_dim_get_spin_index(hm%d, ik2) /= states_elec_dim_get_spin_index(hm%d, ik)) cycle

      if(hm%hf_st%occ(ist, ik2) < M_EPSILON) cycle

      pot = R_TOTYPE(M_ZERO)
      rho = R_TOTYPE(M_ZERO)

      call states_elec_get_state(hm%hf_st, mesh, ist, ik2, psi2)

      do idim = 1, hm%hf_st%d%dim
        forall(ip = 1:mesh%np)
          rho(ip) = rho(ip) + R_CONJ(psi2(ip, idim))*psi(ip, idim)
        end forall
      end do

      call X(poisson_solve)(hm%psolver, pot, rho, all_nodes = .false.)

      ff = hm%hf_st%occ(ist, ik2)
      if(hm%d%ispin == UNPOLARIZED) ff = M_HALF*ff

      do idim = 1, hm%hf_st%d%dim
        forall(ip = 1:mesh%np)
          hpsi(ip, idim) = hpsi(ip, idim) - hm%exx_coef*ff*psi2(ip, idim)*pot(ip)
        end forall
      end do

    end do

    call batch_set_state(hpsib, ibatch, mesh%np, hpsi)
    
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator))
end subroutine X(exchange_operator_hartree)

! scdm_EXX
! ---------------------------------------------------------
subroutine X(scdm_exchange_operator) (hm, namespace, mesh, psib, hpsib, ik)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(batch_t),            intent(inout) :: psib
  type(batch_t),            intent(inout) :: hpsib
  integer,                  intent(in)    :: ik

  integer :: ist, jst, ip, idim, ik2, ibatch
  integer :: ii, jj, kk, ll, count
  FLOAT :: ff, rr(3), dist
  R_TYPE, allocatable :: rho_l(:), pot_l(:), psil(:, :), hpsil(:, :), psi(:, :), hpsi(:, :), temp_state_global(:, :)

  PUSH_SUB(X(scdm_exchange_operator))
  
  call profiling_in(prof_exx_scdm, 'SCDM_EXX_OPERATOR')

  if(mesh%sb%kpoints%full%npoints > 1) call messages_not_implemented("exchange operator with k-points", namespace=namespace)
  
  ! make sure scdm is localized
  call X(scdm_localize)(hm%scdm, namespace, hm%hf_st, mesh)
  
  SAFE_ALLOCATE(psil(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsil(1:mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(psi(1:mesh%np_global, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np_global, 1:hm%d%dim))
  SAFE_ALLOCATE(temp_state_global(mesh%np_global, hm%hf_st%d%dim))
  SAFE_ALLOCATE(rho_l(1:hm%scdm%full_box))
  SAFE_ALLOCATE(pot_l(1:hm%scdm%full_box))
  
  do ibatch = 1, psib%nst
    ist = psib%states(ibatch)%ist
    
    call batch_get_state(psib, ibatch, mesh%np, psil)
    call batch_get_state(hpsib, ibatch, mesh%np, hpsil)

    if(mesh%parallel_in_domains) then
#ifdef HAVE_MPI
      ! the gathering is done for the domain distribution, the states are still local to the st%mpi_grp
      call vec_allgather(mesh%vp, psi(:, 1), psil(:, 1))
      call vec_allgather(mesh%vp, hpsi(:, 1), hpsil(:, 1))
#endif
    else
      psi(1:mesh%np, 1:hm%d%dim) = psil(1:mesh%np, 1:hm%d%dim)
      hpsi(1:mesh%np, 1:hm%d%dim) = hpsil(1:mesh%np, 1:hm%d%dim)
    end if
    
    ! accumulate exchange contribution to Hpsi in a temp array and add to Hpsi at the end
    temp_state_global(:,:) = M_ZERO

    do ik2 = 1, hm%d%nik
      if(states_elec_dim_get_spin_index(hm%d, ik2) /= states_elec_dim_get_spin_index(hm%d, ik)) cycle
      count = 0
      do jst = hm%scdm%st_exx_start, hm%scdm%st_exx_end

        if(hm%hf_st%occ(jst, ik2) < M_EPSILON) cycle
        ! for psi in scdm representation check if it overlaps with the box of jst
        ! NOTE: this can be faster by building an array with overlapping index pairs
        !       within the radius of scdm%box_size
        if(hm%scdm%psi_scdm) then
          do ii = 1, 3
            rr(1:3) = hm%scdm%center(ii,ist) - hm%scdm%center(ii,jst)
          end do
          dist = sqrt(dot_product(rr, rr))
          if(dist .gt. hm%scdm%box_size) cycle
        end if

        ! in Hartree we just remove the self-interaction
        if(hm%theory_level == HARTREE .and. jst /= ist) cycle

        ! for scdm do product only in the local box
        rho_l(:) = M_ZERO

        ! copy density to local box
        do jj = 1, hm%scdm%box_size*2 + 1
          do kk = 1, hm%scdm%box_size*2 + 1
            do ll = 1, hm%scdm%box_size*2 + 1
              ip = (jj - 1)*((hm%scdm%box_size*2 + 1))**2+(kk - 1)*((hm%scdm%box_size*2 + 1)) + ll
              rho_l(ip) = R_CONJ(hm%scdm%X(psi)(ip, jst))*psi(hm%scdm%box(jj, kk, ll, jst), 1)
            end do
          end do
        end do

        call X(poisson_solve)(hm%scdm%poisson, pot_l, rho_l, all_nodes=.false.)

        ff = hm%hf_st%occ(jst, ik2)
        if(hm%d%ispin == UNPOLARIZED) ff = M_HALF*ff
 
        do idim = 1, hm%hf_st%d%dim
          ! potential in local box to full H*psi 
          do jj =1, hm%scdm%box_size*2 + 1
            do kk =1, hm%scdm%box_size*2 + 1
              do ll =1, hm%scdm%box_size*2 + 1
                ip = (jj - 1)*((hm%scdm%box_size*2 + 1))**2 + (kk - 1)*((hm%scdm%box_size*2 + 1)) + ll
                temp_state_global(hm%scdm%box(jj, kk, ll, jst), idim) = &
                  temp_state_global(hm%scdm%box(jj, kk, ll, jst), idim) - hm%exx_coef*ff*hm%scdm%X(psi)(ip, jst)*pot_l(ip)
              end do
            end do
          end do

        end do

      end do
    end do

    ! sum contributions to hpsi from all processes in the st_exx_grp group
    call comm_allreduce(hm%scdm%st_exx_grp%comm, temp_state_global)
    
    ! add exchange contribution to the input state
    hpsi(1:mesh%np_global, 1) =  hpsi(1:mesh%np_global, 1) + temp_state_global(1:mesh%np_global, 1)

    if(mesh%parallel_in_domains) then
#ifdef HAVE_MPI
      call vec_scatter(mesh%vp, 0, hpsil(:, 1), hpsi(:, 1))
#endif
    else
      hpsil(1:mesh%np, 1:hm%d%dim) = hpsi(1:mesh%np, 1:hm%d%dim)
    end if

    call batch_set_state(hpsib, ibatch, mesh%np, hpsil)
  end do
  
  SAFE_DEALLOCATE_A(psil)
  SAFE_DEALLOCATE_A(hpsil)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(temp_state_global)
  SAFE_DEALLOCATE_A(rho_l)
  SAFE_DEALLOCATE_A(pot_l)

  call profiling_out(prof_exx_scdm)
  
  POP_SUB(X(scdm_exchange_operator))
end subroutine X(scdm_exchange_operator)

! ---------------------------------------------------------

subroutine X(magnus) (hm, namespace, mesh, psi, hpsi, ik, vmagnus, set_phase)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ik
  R_TYPE,                   intent(inout) :: psi(:,:)
  R_TYPE,                   intent(out)   :: hpsi(:,:)
  FLOAT,                    intent(in)    :: vmagnus(:, :, :)
  logical, optional,        intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim, ispin

  PUSH_SUB(X(magnus))

  ! We will assume, for the moment, no spinors.
  if(hm%d%dim /= 1) &
    call messages_not_implemented("Magnus with spinors", namespace=namespace)

  SAFE_ALLOCATE( auxpsi(1:mesh%np_part, 1:hm%d%dim))
  SAFE_ALLOCATE(aux2psi(1:mesh%np,      1:hm%d%dim))

  ispin = states_elec_dim_get_spin_index(hm%d, ik)

  ! Compute (T + Vnl)|psi> and store it
  call X(hamiltonian_elec_apply)(hm, namespace, mesh, psi, auxpsi, 1, ik, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, &
    set_phase = set_phase)

  ! H|psi>  =  (T + Vnl)|psi> + Vpsl|psi> + Vmagnus(t2)|psi> + Vborders
  do idim = 1, hm%d%dim
    call lalg_copy(mesh%np, auxpsi(:, idim), hpsi(:, idim))
    hpsi(1:mesh%np, idim) = hpsi(1:mesh%np, idim) + hm%ep%Vpsl(1:mesh%np)*psi(1:mesh%np,idim)
    call X(vborders)(mesh, hm, psi(:, idim), hpsi(:, idim))
  end do
  hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) + vmagnus(1:mesh%np, ispin, 2)*psi(1:mesh%np, 1)

  ! Add first term of the commutator:  - i Vmagnus(t1) (T + Vnl) |psi>
  hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) - M_zI*vmagnus(1:mesh%np, ispin, 1)*auxpsi(1:mesh%np, 1)

  ! Add second term of commutator:  i (T + Vnl) Vmagnus(t1) |psi>
  auxpsi(1:mesh%np, 1) = vmagnus(1:mesh%np, ispin, 1)*psi(1:mesh%np, 1)
  call X(hamiltonian_elec_apply)(hm, namespace, mesh, auxpsi, aux2psi, 1, ik, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, &
    set_phase = set_phase)
  hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) + M_zI*aux2psi(1:mesh%np, 1)

  SAFE_DEALLOCATE_A(auxpsi)
  SAFE_DEALLOCATE_A(aux2psi)
  POP_SUB(X(magnus))
end subroutine X(magnus)

subroutine X(hamiltonian_elec_apply_magnus) (hm, namespace, mesh, psib, hpsib, ik, vmagnus, set_phase)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  integer,                  intent(in)    :: ik
  type(batch_t),            intent(inout) :: psib
  type(batch_t),            intent(inout) :: hpsib
  FLOAT,                    intent(in)    :: vmagnus(:, :, :)
  logical, optional,        intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

  integer :: ispin
  type(batch_t) :: auxpsib, aux2psib

  PUSH_SUB(X(hamiltonian_elec_apply_magnus))

  ! We will assume, for the moment, no spinors.
  if (hm%d%dim /= 1) call messages_not_implemented("Magnus with spinors", namespace=namespace)

  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)

  ispin = states_elec_dim_get_spin_index(hm%d, ik)

  call batch_copy(hpsib, auxpsib, copy_data=.false.)
  call batch_copy(hpsib, aux2psib, copy_data=.false.)
  
  ! Compute (T + Vnl)|psi> and store it
  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib, auxpsib, ik, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, &
    set_phase = set_phase)

  ! H|psi>  =  (T + Vnl)|psi> + Vpsl|psi> + Vmagnus(t2)|psi> + Vborders|psi>
  call batch_copy_data(mesh%np, auxpsib, hpsib)
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
  call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, auxpsib, aux2psib, ik, terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, &
    set_phase = set_phase)
  call batch_axpy(mesh%np, M_zI, aux2psib, hpsib)

  call batch_end(auxpsib)
  call batch_end(aux2psib)

  POP_SUB(X(hamiltonian_elec_apply_magnus))
end subroutine X(hamiltonian_elec_apply_magnus)


! ---------------------------------------------------------
subroutine X(vborders) (mesh, hm, psi, hpsi)
  type(mesh_t),             intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  R_TYPE,              intent(in)    :: psi(:)
  R_TYPE,              intent(inout) :: hpsi(:)

  integer :: ip

  PUSH_SUB(X(vborders))

  if(hm%bc%abtype == IMAGINARY_ABSORBING) then
    forall(ip = 1:mesh%np) hpsi(ip) = hpsi(ip) + M_zI*hm%bc%mf(ip)*psi(ip)
  end if
   
  POP_SUB(X(vborders))
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(h_mgga_terms) (hm, mesh, ik, psib, hpsib)
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(mesh_t),             intent(in)    :: mesh
  integer,             intent(in)    :: ik
  type(batch_t),       intent(inout) :: psib
  type(batch_t),       intent(inout) :: hpsib

  integer :: ispin, ii, idir, ip
  R_TYPE, allocatable :: grad(:,:), diverg(:)
  type(batch_t) :: divb
  type(batch_t), allocatable :: gradb(:)
  
  PUSH_SUB(X(h_mgga_terms))

  ASSERT(.not. batch_is_packed(psib))
  
  ispin = states_elec_dim_get_spin_index(hm%d, ik)

  SAFE_ALLOCATE(grad(1:mesh%np_part, 1:mesh%sb%dim))
  SAFE_ALLOCATE(diverg(1:mesh%np))

  SAFE_ALLOCATE(gradb(1:mesh%sb%dim))

  call batch_copy(hpsib, divb)
  
  do idir = 1, mesh%sb%dim
    call batch_copy(hpsib, gradb(idir))
    call X(derivatives_batch_perform)(hm%der%grad(idir), hm%der, psib, gradb(idir), ghost_update = .false., set_bc = .false.)
  end do
  
  do ii = 1, psib%nst_linear

    do idir = 1, mesh%sb%dim
      call batch_get_state(gradb(idir), ii, mesh%np, grad(:, idir))
    end do
    
    ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
    if (simul_box_is_periodic(mesh%sb) .and. mesh%sb%nonorthogonal ) then
      forall (ip = 1:mesh%np)
        grad(ip, 1:hm%der%dim) = matmul(mesh%sb%klattice_primitive(1:hm%der%dim, 1:hm%der%dim),grad(ip, 1:hm%der%dim))
      end forall
    end if

    do idir = 1, mesh%sb%dim
      grad(1:mesh%np, idir) = grad(1:mesh%np, idir)*hm%vtau(1:mesh%np, ispin)
    end do
    
    call X(derivatives_div)(hm%der, grad, diverg)

    call batch_set_state(divb, ii, mesh%np, diverg)

  end do

  call batch_axpy(mesh%np, CNST(-1.0), divb, hpsib)

  do idir = 1, mesh%sb%dim
    call batch_end(gradb(idir))
  end do

  call batch_end(divb)
  
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
