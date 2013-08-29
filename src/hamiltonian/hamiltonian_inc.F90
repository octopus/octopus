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
!! $Id$

! ---------------------------------------------------------
subroutine X(hamiltonian_apply_batch) (hm, der, psib, hpsib, ik, time, Imtime, terms)
  type(hamiltonian_t),   intent(in)    :: hm
  type(derivatives_t),   intent(in)    :: der
  type(batch_t), target, intent(inout) :: psib
  type(batch_t), target, intent(inout) :: hpsib
  integer,               intent(in)    :: ik
  FLOAT, optional,       intent(in)    :: time
  FLOAT, optional,       intent(in)    :: Imtime
  integer, optional,     intent(in)    :: terms

  integer :: nst, bs, sp
  R_TYPE, pointer     :: epsi(:)
  R_TYPE, allocatable :: psi_copy(:, :, :)

  logical :: kinetic_only_, apply_phase, pack
  integer :: ii, ist, idim, ip
  R_TYPE, pointer :: psi(:), hpsi(:)
  type(batch_t), pointer :: epsib
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection
    
  call profiling_in(prof_hamiltonian, "HAMILTONIAN")
  PUSH_SUB(X(hamiltonian_apply_batch))

  ASSERT(batch_status(psib) == batch_status(hpsib))

  if(present(time)) then
    if(abs(time - hm%current_time) > CNST(1e-10)) then
      write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
      write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', hm%current_time
      call messages_fatal(2)
    endif
  end if

  if(present(Imtime)) then
    ASSERT(abs(Imtime - hm%Imcurrent_time) < CNST(1e-10))
  end if

  ! all terms are enabled by default
  terms_ = optional_default(terms, TERM_ALL)

  kinetic_only_ = (ieor(terms_, TERM_KINETIC) == 0)

  if(present(time) .and. hm%d%cdft) call messages_not_implemented('TDCDFT')

  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)
  ASSERT(ik >= hm%d%kpt%start .and. ik <= hm%d%kpt%end)
  nst = psib%nst_linear

  apply_phase = associated(hm%phase)

  pack = hamiltonian_apply_packed(hm, der%mesh) &
    .and. (opencl_is_enabled() .or. nst > 1) &
    .and. terms_ == TERM_ALL

  if(pack) then
    call batch_pack(psib)
    call batch_pack(hpsib, copy = .false.)
  end if

  call X(derivatives_batch_set_bc)(der, psib)

  if(apply_phase) then
    SAFE_ALLOCATE(epsib)
    SAFE_ALLOCATE(psi_copy(1:der%mesh%np_part, 1:hm%d%dim, 1:psib%nst))
    call batch_init(epsib, hm%d%dim, psib%nst)
    do ii = 1, psib%nst
      call batch_add_state(epsib, psib%states(ii)%ist, psi_copy(:, :, ii))
    end do
    ASSERT(batch_is_ok(epsib))
    if(batch_is_packed(psib)) call batch_pack(epsib, copy = .false.)
  else
    epsib => psib
  end if

  bs = hardware%X(block_size)
 
  if(apply_phase) then ! we copy psi to epsi applying the exp(i k.r) phase
    call X(hamiltonian_phase)(hm, der, der%mesh%np_part, ik, .false., epsib, src = psib)
  end if


  if(iand(TERM_KINETIC, terms_) /= 0) then
    ASSERT(associated(hm%hm_base%kinetic))
    call profiling_in(prof_kinetic_start, "KINETIC_START")
    call X(derivatives_batch_start)(hm%hm_base%kinetic, der, epsib, hpsib, handle, set_bc = .false., factor = -M_HALF/hm%mass)
    call profiling_out(prof_kinetic_start)
  end if

  if (hm%ep%non_local .and. iand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_start)(hm%hm_base, der%mesh, hm%d, ik, epsib, projection)
    end if
  end if

  if(iand(TERM_KINETIC, terms_) /= 0) then
    call profiling_in(prof_kinetic_finish, "KINETIC_FINISH")
    call X(derivatives_batch_finish)(handle)
    call profiling_out(prof_kinetic_finish)

    if(hm%cmplxscl%space) then !cmplxscl
    !complex scale the laplacian 
      do sp = 1, der%mesh%np, bs   
        if(batch_is_packed(hpsib)) then
          forall (ist = 1:hpsib%nst_linear, ip = sp:min(sp + bs - 1, der%mesh%np))
            hpsib%pack%X(psi)(ist, ip) = exp(-M_TWO*M_zI*hm%cmplxscl%theta)*hpsib%pack%X(psi)(ist, ip)
          end forall
        else
          do ii = 1, nst
            call set_pointers()
            forall(ip = sp:min(sp + bs - 1, der%mesh%np)) hpsi(ip) = exp(-M_TWO*M_zI*hm%cmplxscl%theta)*hpsi(ip)
          end do
        end if
      end do    
    end if
  else
    call batch_set_zero(hpsib)
  end if

  ! apply the local potential
  if (iand(TERM_LOCAL_POTENTIAL, terms_) /= 0) then
    call X(hamiltonian_base_local)(hm%hm_base, der%mesh, hm%d, states_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  else if(iand(TERM_LOCAL_EXTERNAL, terms_) /= 0) then
    call X(hamiltonian_external)(hm, der%mesh, epsib, hpsib)
  end if

  ! and the non-local one 
  if (hm%ep%non_local .and. iand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_finish)(hm%hm_base, der%mesh, hm%d, ik, projection, hpsib)
    else
      ASSERT(.not. batch_is_packed(hpsib))
      call X(project_psi_batch)(der%mesh, hm%ep%proj, hm%ep%natoms, hm%d%dim, epsib, hpsib, ik)
    end if
  end if

  if (iand(TERM_OTHERS, terms_) /= 0 .and. hamiltonian_base_has_magnetic(hm%hm_base)) then
    call X(hamiltonian_base_magnetic)(hm%hm_base, der, hm%d, hm%ep, states_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  end if

  if (iand(TERM_OTHERS, terms_) /= 0 ) then 
    call X(hamiltonian_base_rashba)(hm%hm_base, der, hm%d, epsib, hpsib)
  end if

  if (iand(TERM_OTHERS, terms_) /= 0) then

    if(hm%theory_level == HARTREE .or. hm%theory_level == HARTREE_FOCK) then
      ASSERT(.not. batch_is_packed(hpsib))
      do ii = 1, psib%nst
        call X(exchange_operator)(hm, der, epsib%states(ii)%X(psi), hpsib%states(ii)%X(psi), psib%states(ii)%ist, ik, hm%exx_coef)
      end do
    end if

    if(hm%ab == IMAGINARY_ABSORBING) then
      ASSERT(.not. batch_is_packed(hpsib))
      do ii = 1, nst
        call set_pointers()
        if(present(time)) call X(vborders)(der, hm, epsi, hpsi)      
      end do
    end if

  end if

  if (iand(TERM_MGGA, terms_) /= 0 .and. iand(hm%xc_family, XC_FAMILY_MGGA) /= 0) then
    ASSERT(.not. batch_is_packed(hpsib))
    do ii = 1, nst
      call set_pointers()
      call X(h_mgga_terms)(hm, der, epsi, hpsi, ik)
    end do
  end if

  if(apply_phase) then
    call X(hamiltonian_phase)(hm, der, der%mesh%np, ik, .true., hpsib)
    if(batch_is_packed(epsib)) call batch_unpack(epsib, copy = .false.)
    call batch_end(epsib)
    SAFE_DEALLOCATE_P(epsib)
  end if

  if(pack) then
    call batch_unpack(psib, copy = .false.)
    call batch_unpack(hpsib)
  end if

  POP_SUB(X(hamiltonian_apply_batch))
  call profiling_out(prof_hamiltonian)

contains

  subroutine set_pointers
    ! no push_sub since called very frequently

    psi  => psib%states_linear(ii)%X(psi)
    epsi => epsib%states_linear(ii)%X(psi)
    hpsi => hpsib%states_linear(ii)%X(psi)
    ist  =  batch_linear_to_ist(psib, ii)
    idim =  batch_linear_to_idim(psib, ii)
  end subroutine set_pointers

end subroutine X(hamiltonian_apply_batch)

! ---------------------------------------------------------

subroutine X(hamiltonian_external)(this, mesh, psib, vpsib)
  type(hamiltonian_t),         intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  FLOAT, allocatable :: vpsl_spin(:,:), Imvpsl_spin(:,:)
#ifdef HAVE_OPENCL
  integer :: pnp, offset, ispin
  type(opencl_mem_t) :: vpsl_buff
#endif

  PUSH_SUB(X(hamiltonian_external))

  SAFE_ALLOCATE(vpsl_spin(mesh%np, this%d%nspin))
  vpsl_spin(1:mesh%np, 1) = this%ep%vpsl(1:mesh%np)
  if(this%d%ispin == SPINORS) then
    ! yes this means a little unnecessary computation in the later call,
    ! but with the great benefit of being able to reuse an existing routine
    vpsl_spin(1:mesh%np, 2) = this%ep%vpsl(1:mesh%np)
    vpsl_spin(1:mesh%np, 3) = M_ZERO
    vpsl_spin(1:mesh%np, 4) = M_ZERO
  endif

  if(this%cmplxscl%space) then
    SAFE_ALLOCATE(Imvpsl_spin(mesh%np, this%d%nspin))
    Imvpsl_spin(1:mesh%np, 1) = this%ep%Imvpsl(1:mesh%np)
    if(this%d%ispin == SPINORS) then
      Imvpsl_spin(1:mesh%np, 2) = this%ep%Imvpsl(1:mesh%np)
      Imvpsl_spin(1:mesh%np, 3) = M_ZERO
      Imvpsl_spin(1:mesh%np, 4) = M_ZERO
    endif
  endif

  if(batch_status(psib) == BATCH_CL_PACKED) then
#ifdef HAVE_OPENCL
    pnp = opencl_padded_size(mesh%np)
    call opencl_create_buffer(vpsl_buff, CL_MEM_READ_ONLY, TYPE_FLOAT, pnp * this%d%nspin)
    call opencl_write_buffer(vpsl_buff, mesh%np, this%ep%vpsl)

    offset = 0
    do ispin = 1, this%d%nspin
       call opencl_write_buffer(vpsl_buff, mesh%np, vpsl_spin(:, ispin), offset = offset)
       offset = offset + pnp
    end do

    call X(hamiltonian_base_local_sub)(vpsl_spin, mesh, this%d, 1, &
      psib, vpsib, potential_opencl = vpsl_buff)

    call opencl_release_buffer(vpsl_buff)
#endif
  else
    if(this%cmplxscl%space) then
      call X(hamiltonian_base_local_sub)(vpsl_spin, mesh, this%d, 1, &
        psib, vpsib, Impotential = Imvpsl_spin)
    else
      call X(hamiltonian_base_local_sub)(vpsl_spin, mesh, this%d, 1, psib, vpsib)
    endif
  endif

  SAFE_DEALLOCATE_A(vpsl_spin)
  if(this%cmplxscl%space) then
    SAFE_DEALLOCATE_A(Imvpsl_spin)
  endif
  
  POP_SUB(X(hamiltonian_external))
end subroutine X(hamiltonian_external)

! ---------------------------------------------------------

subroutine X(hamiltonian_apply) (hm, der, psi, hpsi, ist, ik, time, terms, Imtime)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ist       !< the index of the state
  integer,             intent(in)    :: ik        !< the index of the k-point
  R_TYPE,   target,    intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
  R_TYPE,   target,    intent(out)   :: hpsi(:,:) !< (gr%mesh%np, hm%d%dim)
  FLOAT,    optional,  intent(in)    :: time
  integer,  optional,  intent(in)    :: terms
  FLOAT,    optional,  intent(in)    :: Imtime

  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_apply))

  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  if(present(time)) then
    if(present(terms)) then
      if(present(Imtime)) then
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time, terms = terms, Imtime = Imtime)
      else       
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time, terms = terms)
      end if
    else
      if(present(Imtime)) then
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time, Imtime = Imtime)
      else 
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time)
      end if
    endif
  else
    if(present(terms)) then
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, terms = terms)
    else
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)
    endif
  endif

  call batch_end(psib)
  call batch_end(hpsib)

  POP_SUB(X(hamiltonian_apply))
end subroutine X(hamiltonian_apply)


! ---------------------------------------------------------
subroutine X(hamiltonian_apply_all) (hm, der, st, hst, time, Imtime)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  type(states_t),      intent(inout) :: hst
  FLOAT, optional,     intent(in)    :: time
  FLOAT, optional,     intent(in)    :: Imtime

  integer :: ik, ib

  PUSH_SUB(X(hamiltonian_apply_all))

  if(present(Imtime)) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call X(hamiltonian_apply_batch)(hm, der, st%group%psib(ib, ik), hst%group%psib(ib, ik), ik, time, Imtime)
      end do
    end do
  else 
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call X(hamiltonian_apply_batch)(hm, der, st%group%psib(ib, ik), hst%group%psib(ib, ik), ik, time)
      end do
    end do
  end if
  
  if(hamiltonian_oct_exchange(hm)) call X(oct_exchange_operator_all)(hm, der, st, hst)

  POP_SUB(X(hamiltonian_apply_all))
end subroutine X(hamiltonian_apply_all)


! ---------------------------------------------------------
subroutine X(exchange_operator) (hm, der, psi, hpsi, ist, ik, exx_coef)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: psi(:,:)
  R_TYPE,              intent(inout) :: hpsi(:,:)
  integer,             intent(in)    :: ist
  integer,             intent(in)    :: ik
  FLOAT,               intent(in)    :: exx_coef

  R_TYPE, allocatable :: rho(:), pot(:), psi2(:, :)
  integer :: jst, ip, idim, ik2

  FLOAT :: ff

  PUSH_SUB(X(exchange_operator))

  if(der%mesh%sb%kpoints%reduced%npoints > 1) call messages_not_implemented("exchange operator with k-points")
  if(hm%hf_st%parallel_in_states) call messages_not_implemented("exchange operator parallel in states")

  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:hm%d%dim))

  do ik2 = 1, hm%d%nik
    if(states_dim_get_spin_index(hm%d, ik2) /= states_dim_get_spin_index(hm%d, ik)) cycle

    do jst = 1, hm%hf_st%nst
      if(hm%hf_st%occ(jst, ik2) < M_EPSILON) cycle

      ! in Hartree we just remove the self-interaction
      if(hm%theory_level == HARTREE .and. jst /= ist) cycle

      pot = R_TOTYPE(M_ZERO)
      rho = R_TOTYPE(M_ZERO)

      call states_get_state(hm%hf_st, der%mesh, jst, ik2, psi2)
    
      if(hm%cmplxscl%space) psi2 = R_CONJ(psi2)

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
  end do

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  SAFE_DEALLOCATE_A(psi2)

  POP_SUB(X(exchange_operator))
end subroutine X(exchange_operator)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator_all) (hm, der, st, hst)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(in)    :: st
  type(states_t),      intent(inout) :: hst

  integer :: ik, ist
  integer :: ik2
  FLOAT,  allocatable :: rho(:, :), pot(:, :)
  R_TYPE, allocatable :: psi(:, :), psi2(:, :), hpsi(:, :)
  integer :: jst, ip

  PUSH_SUB(X(oct_exchange_operator_all))


  SAFE_ALLOCATE(rho(1:der%mesh%np, 1:st%d%nspin))
  SAFE_ALLOCATE(pot(1:der%mesh%np, 1:st%d%nspin))
  SAFE_ALLOCATE(psi(1:der%mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(psi2(1:der%mesh%np, 1:hm%d%dim))
  SAFE_ALLOCATE(hpsi(1:der%mesh%np, 1:hm%d%dim))

  select case(hm%d%ispin)
  case(UNPOLARIZED)
    ASSERT(st%d%nik  ==  1)

    pot = M_ZERO
    rho = M_ZERO
    do jst = 1, hm%oct_st%nst
      call states_get_state(st, der%mesh, jst, 1, psi)
      call states_get_state(hm%oct_st, der%mesh, jst, 1, psi2)
      forall (ip = 1:der%mesh%np)
        rho(ip, 1) = rho(ip, 1) + hm%oct_st%occ(jst, 1)*R_AIMAG(R_CONJ(psi2(ip, 1))*psi(ip, 1))
      end forall
    end do
    call dpoisson_solve(psolver, pot(:, 1), rho(:, 1), all_nodes = .false.)
    do ist = st%st_start, st%st_end
      call states_get_state(hst, der%mesh, ist, 1, hpsi)
      call states_get_state(hm%oct_st, der%mesh, ist, 1, psi2)
      forall(ip = 1:der%mesh%np)
        hpsi(ip, 1) = hpsi(ip, 1) + M_TWO*M_zI*psi2(ip, 1)*(pot(ip, 1) + hm%oct_fxc(ip, 1, 1)*rho(ip, 1))
      end forall
      call states_set_state(hst, der%mesh, ist, 1, hpsi)
    end do

  case(SPIN_POLARIZED)
    ASSERT(st%d%nik  ==  2)

    pot = M_ZERO
    rho = M_ZERO
    do ik = 1, 2
      do jst = 1, hm%oct_st%nst
        call states_get_state(st, der%mesh, jst, ik, psi)
        call states_get_state(hm%oct_st, der%mesh, jst, ik, psi2)
        forall (ip = 1:der%mesh%np)
          rho(ip, ik) = rho(ip, ik) + hm%oct_st%occ(jst, ik) * R_AIMAG(R_CONJ(psi2(ip, 1))*psi(ip, 1))
        end forall
      end do
    end do

    do ik = 1, 2
      call dpoisson_solve(psolver, pot(:, ik), rho(:, ik), all_nodes = .false.)
    end do

    do ik = 1, 2
      do ist = st%st_start, st%st_end
        call states_get_state(hst, der%mesh, ist, ik, hpsi)
        call states_get_state(hm%oct_st, der%mesh, ist, ik, psi2)

        do ik2 = 1, 2
          forall(ip = 1:der%mesh%np)
            hpsi(ip, 1) = hpsi(ip, 1) + M_TWO * M_zI * hm%oct_st%occ(ist, ik) * &
              psi2(ip, 1) * (pot(ip, ik2) + hm%oct_fxc(ip, ik, ik2)*rho(ip, ik2))
          end forall
        end do

        call states_set_state(hst, der%mesh, ist, ik, hpsi)
      end do
    end do

  case(SPINORS)
    call messages_not_implemented("Function oct_exchange_operator_all for spin_polarized or spinors")
  end select

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi2)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  POP_SUB(X(oct_exchange_operator_all))
end subroutine X(oct_exchange_operator_all)



! ---------------------------------------------------------
subroutine X(magnus) (hm, der, psi, hpsi, ik, vmagnus)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: psi(:,:)
  R_TYPE,              intent(out)   :: hpsi(:,:)
  FLOAT,               intent(in)    :: vmagnus(:, :, :)

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim, ispin

  ! We will assume, for the moment, no spinors.

  PUSH_SUB(X(magnus))

  SAFE_ALLOCATE( auxpsi(1:der%mesh%np_part, 1:hm%d%dim))
  SAFE_ALLOCATE(aux2psi(1:der%mesh%np,      1:hm%d%dim))

  ispin = states_dim_get_spin_index(hm%d, ik)

  call X(hamiltonian_apply)(hm, der, psi, hpsi, ist = 1, ik = ik, terms = TERM_KINETIC)

  do idim = 1, hm%d%dim
    call lalg_copy(der%mesh%np, hpsi(:, idim), auxpsi(:, idim))
  end do

  if (hm%ep%non_local) call X(hamiltonian_apply)(hm, der, psi, hpsi, ist = 1, ik = ik, terms = TERM_NON_LOCAL_POTENTIAL)

  hpsi(1:der%mesh%np, 1) = hpsi(1:der%mesh%np, 1) -  M_zI*vmagnus(1:der%mesh%np, ispin, 1)*auxpsi(1:der%mesh%np, 1)
  auxpsi(1:der%mesh%np, 1) = vmagnus(1:der%mesh%np, ispin, 1)*psi(1:der%mesh%np, 1)

  call X(hamiltonian_apply)(hm, der, auxpsi, aux2psi, ist = 1, ik = ik, terms = TERM_KINETIC)

  if (hm%ep%non_local) call X(hamiltonian_apply)(hm, der, psi, hpsi, ist = 1, ik = ik, terms = TERM_NON_LOCAL_POTENTIAL)

  hpsi(1:der%mesh%np, 1) = hpsi(1:der%mesh%np, 1) + M_zI*aux2psi(1:der%mesh%np, 1)

  do idim = 1, hm%d%dim
    hpsi(1:der%mesh%np, idim) = hpsi(1:der%mesh%np, idim) + hm%ep%Vpsl(1:der%mesh%np)*psi(1:der%mesh%np,idim)
  end do

  hpsi(1:der%mesh%np, 1) = hpsi(1:der%mesh%np, 1) + vmagnus(1:der%mesh%np, ispin, 2)*psi(1:der%mesh%np, 1)

  if (hm%ep%non_local) call X(hamiltonian_apply)(hm, der, psi, hpsi, ist = 1, ik = ik, terms = TERM_NON_LOCAL_POTENTIAL)

  do idim = 1, hm%d%dim
    call X(vborders)(der, hm, psi(:, idim), hpsi(:, idim))
  end do

  SAFE_DEALLOCATE_A(auxpsi)
  SAFE_DEALLOCATE_A(aux2psi)
  POP_SUB(X(magnus))
end subroutine X(magnus)


! ---------------------------------------------------------
subroutine X(vborders) (der, hm, psi, hpsi)
  type(derivatives_t), intent(in)    :: der
  type(hamiltonian_t), intent(in)    :: hm
  R_TYPE,              intent(in)    :: psi(:)
  R_TYPE,              intent(inout) :: hpsi(:)

  integer :: ip

  PUSH_SUB(X(vborders))

  if(hm%ab == IMAGINARY_ABSORBING) then
    forall(ip = 1:der%mesh%np) hpsi(ip) = hpsi(ip) + M_zI*hm%ab_pot(ip)*psi(ip)
  end if
  
  POP_SUB(X(vborders))
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(h_mgga_terms) (hm, der, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: psi(:)
  R_TYPE,              intent(inout) :: hpsi(:)
  integer,             intent(in)    :: ik

  integer :: ispace, ispin
  R_TYPE, allocatable :: cgrad(:,:), diverg(:), grad(:, :)

  PUSH_SUB(X(h_mgga_terms))

  ispin = states_dim_get_spin_index(hm%d, ik)

  SAFE_ALLOCATE(grad(1:der%mesh%np_part, 1:MAX_DIM))
  SAFE_ALLOCATE(cgrad(1:der%mesh%np_part, 1:MAX_DIM))
  SAFE_ALLOCATE(diverg(1:der%mesh%np))

  call X(derivatives_grad)(der, psi, grad, ghost_update = .false., set_bc = .false.)
  
  do ispace = 1, der%mesh%sb%dim
    cgrad(1:der%mesh%np, ispace) = grad(1:der%mesh%np, ispace)*hm%vtau(1:der%mesh%np, ispin)
  end do
  
  diverg(1:der%mesh%np) = M_ZERO
  call X(derivatives_div)(der, cgrad, diverg)
  
  hpsi(1:der%mesh%np) = hpsi(1:der%mesh%np) - diverg(1:der%mesh%np)

  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(cgrad)
  SAFE_DEALLOCATE_A(diverg)

  POP_SUB(X(h_mgga_terms))
end subroutine X(h_mgga_terms)


! ---------------------------------------------------------
subroutine X(vmask) (gr, hm, st)
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim
  R_TYPE, allocatable :: psi(:)

  PUSH_SUB(X(vmask))

  SAFE_ALLOCATE(psi(1:gr%mesh%np))

  if(hm%ab == MASK_ABSORBING) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_get_state(st, gr%mesh, idim, ist, ik, psi)
          psi(1:gr%mesh%np) = psi(1:gr%mesh%np)*(M_ONE - hm%ab_pot(1:gr%mesh%np))
          call states_set_state(st, gr%mesh, idim, ist, ik, psi)
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(vmask))
end subroutine X(vmask)


! ---------------------------------------------------------
subroutine X(hamiltonian_diagonal) (hm, der, diag, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(out)   :: diag(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  integer,             intent(in)    :: ik

  integer :: idim, ip, ispin

  FLOAT, allocatable  :: ldiag(:)

  PUSH_SUB(X(hamiltonian_diagonal))
  
  SAFE_ALLOCATE(ldiag(1:der%mesh%np))

  diag = M_ZERO

  call derivatives_lapl_diag(der, ldiag)

  do idim = 1, hm%d%dim
    diag(1:der%mesh%np, idim) = -M_HALF/hm%mass*ldiag(1:der%mesh%np)
  end do

  select case(hm%d%ispin)

  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_dim_get_spin_index(hm%d, ik)
    diag(1:der%mesh%np, 1) = diag(1:der%mesh%np, 1) + hm%vhxc(1:der%mesh%np, ispin) + hm%ep%vpsl(1:der%mesh%np)
    
  case(SPINORS)
    do ip = 1, der%mesh%np
      diag(ip, 1) = diag(ip, 1) + hm%vhxc(ip, 1) + hm%ep%vpsl(ip)
      diag(ip, 2) = diag(ip, 2) + hm%vhxc(ip, 2) + hm%ep%vpsl(ip)
    end do
    
  end select
    
  POP_SUB(X(hamiltonian_diagonal))
end subroutine X(hamiltonian_diagonal)


! ---------------------------------------------------------
subroutine X(hamiltonian_phase)(this, der, np, iqn, conjugate, psib, src)
  type(hamiltonian_t),                   intent(in)    :: this
  type(derivatives_t),                   intent(in)    :: der
  integer,                               intent(in)    :: np
  integer,                               intent(in)    :: iqn
  logical,                               intent(in)    :: conjugate
  type(batch_t),                 target, intent(inout) :: psib
  type(batch_t),       optional, target, intent(in)    :: src

  integer :: ip, ii
  type(batch_t), pointer :: src_
  type(profile_t), save :: phase_prof
#ifdef HAVE_OPENCL
  integer :: wgsize
  type(octcl_kernel_t), save :: ker_phase
  type(cl_kernel) :: kernel
#endif

  PUSH_SUB(X(hamiltonian_phase))
  call profiling_in(phase_prof, "PBC_PHASE_APPLY")

  ASSERT(np <= der%mesh%np_part)

  src_ => psib
  if(present(src)) src_ => src

  select case(batch_status(psib))
  case(BATCH_PACKED)
    
    if(conjugate) then

      do ip = 1, np
        forall (ii = 1:psib%nst_linear)
          psib%pack%X(psi)(ii, ip) = conjg(this%phase(ip, iqn))*src_%pack%X(psi)(ii, ip)
        end forall
      end do
      
    else

      do ip = 1, np
        forall (ii = 1:psib%nst_linear)
          psib%pack%X(psi)(ii, ip) = this%phase(ip, iqn)*src_%pack%X(psi)(ii, ip)
        end forall
      end do

    end if
    
  case(BATCH_NOT_PACKED)

    if(conjugate) then

      do ii = 1, psib%nst_linear
        forall(ip = 1:np)
          psib%states_linear(ii)%X(psi)(ip) = conjg(this%phase(ip, iqn))*src_%states_linear(ii)%X(psi)(ip)
        end forall
      end do
      
    else

      do ii = 1, psib%nst_linear
        forall(ip = 1:np)
          psib%states_linear(ii)%X(psi)(ip) = this%phase(ip, iqn)*src_%states_linear(ii)%X(psi)(ip)
        end forall
      end do

    end if

  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(ker_phase, 'phase.cl', 'phase_hamiltonian')
    kernel = octcl_kernel_get_ref(ker_phase)

    if(conjugate) then
      call opencl_set_kernel_arg(kernel, 0, 1_4)
    else
      call opencl_set_kernel_arg(kernel, 0, 0_4)
    end if
      
    call opencl_set_kernel_arg(kernel, 1, (iqn - this%d%kpt%start)*der%mesh%np_part)
    call opencl_set_kernel_arg(kernel, 2, np)
    call opencl_set_kernel_arg(kernel, 3, this%buff_phase)
    call opencl_set_kernel_arg(kernel, 4, src_%pack%buffer)
    call opencl_set_kernel_arg(kernel, 5, log2(src_%pack%size(1)))
    call opencl_set_kernel_arg(kernel, 6, psib%pack%buffer)
    call opencl_set_kernel_arg(kernel, 7, log2(psib%pack%size(1)))
    
    wgsize = opencl_kernel_workgroup_size(kernel)/psib%pack%size(1)

    call opencl_kernel_run(kernel, (/psib%pack%size(1), pad(np, wgsize)/), (/psib%pack%size(1), wgsize/))

    call opencl_finish()
#endif
  end select

  call batch_pack_was_modified(psib)


  call profiling_out(phase_prof)
  POP_SUB(X(hamiltonian_phase))
end subroutine X(hamiltonian_phase)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
