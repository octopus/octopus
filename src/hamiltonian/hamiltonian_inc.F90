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
subroutine X(hamiltonian_apply_batch) (hm, der, psib, hpsib, ik, time, terms)
  type(hamiltonian_t),   intent(in)    :: hm
  type(derivatives_t),   intent(in)    :: der
  type(batch_t), target, intent(inout) :: psib
  type(batch_t),         intent(inout) :: hpsib
  integer,               intent(in)    :: ik
  FLOAT, optional,       intent(in)    :: time
  integer, optional,     intent(in)    :: terms

  integer :: nst, bs, sp
  R_TYPE, pointer     :: epsi(:,:)
  R_TYPE, pointer     :: grad(:, :, :)
  R_TYPE, allocatable :: psi_copy(:, :, :)

  type(profile_t), save :: phase_prof
  logical :: kinetic_only_, apply_phase, pack
  integer :: ii, ist, idim, ip
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)
  type(batch_t), pointer :: epsib
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection

  call profiling_in(prof_hamiltonian, "HAMILTONIAN")
  PUSH_SUB(X(hamiltonian_apply_batch))

  if(present(time)) then
    ASSERT(abs(time - hm%current_time) < CNST(1e-10))
  end if

  ! all terms are enabled by default
  terms_ = optional_default(terms, TERM_ALL)

  kinetic_only_ = (ieor(terms_, TERM_KINETIC) == 0)

  if(present(time) .and. hm%d%cdft) then
    message(1) = "TDCDFT not yet implemented."
    call messages_fatal(1)
  end if

  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)
  ASSERT(ik >= hm%d%kpt%start .and. ik <= hm%d%kpt%end)
  nst = psib%nst

  apply_phase = associated(hm%phase)

  call X(derivatives_batch_set_bc)(der, psib)

  pack = hamiltonian_apply_packed(hm, der%mesh) &
    .and. (opencl_is_enabled() .or. nst > 1) &
    .and. terms_ == TERM_ALL

  if(pack) then
    call batch_pack(psib)
    call batch_pack(hpsib, copy = .false.)
  end if

  if(apply_phase) then
    SAFE_ALLOCATE(epsib)
    SAFE_ALLOCATE(psi_copy(1:der%mesh%np_part, 1:hm%d%dim, 1:nst))
    call batch_init(epsib, hm%d%dim, nst)
    do ii = 1, nst
      call batch_add_state(epsib, psib%states(ii)%ist, psi_copy(:, :, ii))
    end do
    ASSERT(batch_is_ok(epsib))
    if(batch_is_packed(psib)) call batch_pack(epsib, copy = .false.)
  else
    epsib => psib
  end if

  bs = hardware%X(block_size)

  if(apply_phase) then ! we copy psi to epsi applying the exp(i k.r) phase

    call profiling_in(phase_prof, "PBC_PHASE_APPLY")
    do sp = 1, der%mesh%np_part, bs
      if(batch_is_packed(psib)) then
        forall (ist = 1:psib%nst_linear, ip = sp:min(sp + bs - 1, der%mesh%np_part))
          epsib%pack%X(psi)(ist, ip) = hm%phase(ip, ik)*psib%pack%X(psi)(ist, ip)
        end forall
      else
        do ii = 1, nst
          call set_pointers()
          forall (idim = 1:hm%d%dim, ip = sp:min(sp + bs - 1, der%mesh%np_part)) 
            psi_copy(ip, idim, ii) = hm%phase(ip, ik)*psi(ip, idim)
          end forall
        end do
      end if
    end do
    call profiling_out(phase_prof)
  end if


  if(iand(TERM_KINETIC, terms_) /= 0) then
    ASSERT(associated(hm%hm_base%kinetic))
    call profiling_in(prof_kinetic, "KINETIC")
    call X(derivatives_batch_start)(hm%hm_base%kinetic, der, epsib, hpsib, handle, set_bc = .false., factor = -M_HALF/hm%mass)
    call profiling_out(prof_kinetic)
  end if

  if (hm%ep%non_local .and. iand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_start)(hm%hm_base, der%mesh, hm%d, ik, epsib, projection)
    end if
  end if

  if(iand(TERM_KINETIC, terms_) /= 0) then
    call profiling_in(prof_kinetic, "KINETIC")
    call X(derivatives_batch_finish)(handle)
    call profiling_out(prof_kinetic)
  else
    call batch_set_zero(hpsib)
  end if

  ! apply the local potential
  if (iand(TERM_LOCAL_POTENTIAL, terms_) /= 0) then
    call X(hamiltonian_base_local)(hm%hm_base, der%mesh, hm%d, states_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  else if(iand(TERM_LOCAL_EXTERNAL, terms_) /= 0) then
    do ii = 1, nst
      call set_pointers()
      do idim = 1, hm%d%dim
        hpsi(1:der%mesh%np, idim) = hpsi(1:der%mesh%np, idim) + &
          hm%ep%vpsl(1:der%mesh%np)*epsi(1:der%mesh%np, idim)
      end do
    end do
  end if

  ! and the non-local one 
  if (hm%ep%non_local .and. iand(TERM_NON_LOCAL_POTENTIAL, terms_) /= 0) then
    if(hm%hm_base%apply_projector_matrices) then
      call X(hamiltonian_base_nlocal_finish)(hm%hm_base, der%mesh, hm%d, ik, projection, hpsib)
    else
      call X(project_psi_batch)(der%mesh, hm%ep%proj, hm%ep%natoms, hm%d%dim, epsib, hpsib, ik)
    end if
  end if

  if (iand(TERM_OTHERS, terms_) /= 0 .and. hamiltonian_base_has_magnetic(hm%hm_base)) then
    call X(hamiltonian_base_magnetic)(hm%hm_base, der, hm%d, hm%ep, states_dim_get_spin_index(hm%d, ik), epsib, hpsib)
  end if

  do ii = 1, nst
    call set_pointers()

    if (iand(TERM_OTHERS, terms_) /= 0) then

      ! all functions that require the gradient or other derivatives of
      ! epsi should go after this point and calculate it using Xget_grad

      ! FIX ME: here not all the possible couplings between vector
      ! potentials are considered

      nullify(grad)

      if(hm%theory_level == HARTREE .or. hm%theory_level == HARTREE_FOCK) call X(exchange_operator)(hm, der, epsi, hpsi, ist, ik)

      if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) &
        call X(h_mgga_terms) (hm, der, epsi, hpsi, ik, grad)

      if(present(time)) call X(vborders) (der, hm, epsi, hpsi)

      SAFE_DEALLOCATE_P(grad)

    end if
  end do

  if(apply_phase) then
    call profiling_in(phase_prof)
    ! now we need to remove the exp(-i k.r) factor
    do sp = 1, der%mesh%np, bs
      
      if(batch_is_packed(hpsib)) then
        forall (ist = 1:hpsib%nst_linear, ip = sp:min(sp + bs - 1, der%mesh%np))
          hpsib%pack%X(psi)(ist, ip) = conjg(hm%phase(ip, ik))*hpsib%pack%X(psi)(ist, ip)
        end forall
      else
        do ii = 1, nst
          call set_pointers()
          
          forall (idim = 1:hm%d%dim, ip = sp:min(sp + bs - 1, der%mesh%np))
            hpsi(ip, idim) = conjg(hm%phase(ip, ik))*hpsi(ip, idim)
          end forall
        end do
      end if
    end do

    call profiling_out(phase_prof)

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

    psi  => psib%states(ii)%X(psi)
    epsi => epsib%states(ii)%X(psi)
    hpsi => hpsib%states(ii)%X(psi)
    ist  =  psib%states(ii)%ist
  end subroutine set_pointers

end subroutine X(hamiltonian_apply_batch)


! ---------------------------------------------------------
subroutine X(get_grad)(hm, der, psi, grad)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: psi(:, :)
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: idim

  PUSH_SUB(X(get_grad))

  if( .not. associated(grad)) then
    SAFE_ALLOCATE(grad(1:der%mesh%np, 1:MAX_DIM, 1:hm%d%dim))
    do idim = 1, hm%d%dim 
      ! boundary points were already set by the Laplacian
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
  end if
  
  POP_SUB(X(get_grad))
end subroutine X(get_grad)       


! ---------------------------------------------------------
subroutine X(hamiltonian_apply) (hm, der, psi, hpsi, ist, ik, time, terms)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ist       ! the index of the state
  integer,             intent(in)    :: ik        ! the index of the k-point
  R_TYPE,   target,    intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  FLOAT,    optional,  intent(in)    :: time
  integer,  optional,  intent(in)    :: terms

  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_apply))

  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  if(present(time)) then
    if(present(terms)) then
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time, terms = terms)
    else
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time)
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
subroutine X(hamiltonian_apply_all) (hm, der, psi, hpsi, time)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: psi
  type(states_t),      intent(inout) :: hpsi
  FLOAT, optional,     intent(in)    :: time

  integer :: ik
  type(batch_t) :: psib, hpsib

  PUSH_SUB(X(hamiltonian_apply_all))

  do ik = psi%d%kpt%start, psi%d%kpt%end
    call batch_init(psib, hm%d%dim, psi%st_start, psi%st_end, psi%X(psi)(:, :, :, ik))
    call batch_init(hpsib, hm%d%dim, hpsi%st_start, hpsi%st_end, hpsi%X(psi)(:, :, :, ik))
    if(present(time)) then
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik, time = time)
    else
      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)
    endif
    call batch_end(psib)
    call batch_end(hpsib)
  end do
  
  if(hamiltonian_oct_exchange(hm)) call X(oct_exchange_operator_all)(hm, der, psi, hpsi)

  POP_SUB(X(hamiltonian_apply_all))
end subroutine X(hamiltonian_apply_all)


! ---------------------------------------------------------
subroutine X(exchange_operator) (hm, der, psi, hpsi, ist, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: psi(:,:)
  R_TYPE,              intent(inout) :: hpsi(:,:)
  integer,             intent(in)    :: ist
  integer,             intent(in)    :: ik

  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: jst, ip, idim

  FLOAT :: ff

  PUSH_SUB(X(exchange_operator))

  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))

  do jst = 1, hm%hf_st%nst
    if(hm%hf_st%occ(jst, ik) <= M_ZERO) cycle

    ! in Hartree we just remove the self-interaction
    if(hm%theory_level == HARTREE .and. jst .ne. ist) cycle

    pot = M_ZERO
    rho = M_ZERO

    do idim = 1, hm%hf_st%d%dim
      forall(ip = 1:der%mesh%np)
        rho(ip) = rho(ip) + R_CONJ(hm%hf_st%X(psi)(ip, idim, jst, ik))*psi(ip, idim)
      end forall
    end do
    
    call X(poisson_solve)(psolver, pot, rho)

    ff = hm%hf_st%occ(jst, ik)
    if(hm%d%ispin == UNPOLARIZED) ff = M_HALF*ff

    do idim = 1, hm%hf_st%d%dim
      forall(ip = 1:der%mesh%np)
        hpsi(ip, idim) = hpsi(ip, idim) - hm%exx_coef*ff*hm%hf_st%X(psi)(ip, idim, jst, ik)*pot(ip)
      end forall
    end do

  end do

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  POP_SUB(X(exchange_operator))
end subroutine X(exchange_operator)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator_all) (hm, der, psi, hpsi)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(in)    :: psi
  type(states_t),      intent(inout) :: hpsi

  integer :: ik, ist
  FLOAT, allocatable :: rho(:), pot(:)
  integer :: jst, ip

  PUSH_SUB(X(oct_exchange_operator_all))

  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))

  select case(hm%d%ispin)
  case(UNPOLARIZED)

    do ik = psi%d%kpt%start, psi%d%kpt%end

      pot = M_ZERO
      rho = M_ZERO
      do jst = 1, hm%oct_st%nst
        forall (ip = 1:der%mesh%np)
          rho(ip) = rho(ip) + hm%oct_st%occ(jst, 1) * &
            R_AIMAG(R_CONJ(hm%oct_st%X(psi)(ip, 1, jst, ik))*psi%X(psi)(ip, 1, jst, ik))
        end forall
      end do 
      call dpoisson_solve(psolver, pot, rho)
      do ist = psi%st_start, psi%st_end
        forall(ip = 1:der%mesh%np)
          hpsi%X(psi)(ip, 1, ist, ik) = hpsi%X(psi)(ip, 1, ist, ik) + M_TWO*M_zI* &
            hm%oct_st%X(psi)(ip, 1, ist, ik)*(pot(ip) + hm%oct_fxc(ip, 1, 1)*rho(ip))
        end forall
      end do

    end do

  case(SPIN_POLARIZED, SPINORS)
    call messages_not_implemented("Function oct_exchange_operator_all for spin_polarized or spinors")
  end select

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  POP_SUB(X(oct_exchange_operator_all))
end subroutine X(oct_exchange_operator_all)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator) (hm, der, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(in)    :: psi(:, :)
  R_TYPE,              intent(inout) :: hpsi(:, :)
  integer,             intent(in)    :: ik

  FLOAT, allocatable :: rho(:), pot(:)
  integer :: jst, ip

  PUSH_SUB(X(oct_exchange_operator))

  SAFE_ALLOCATE(rho(1:der%mesh%np))
  SAFE_ALLOCATE(pot(1:der%mesh%np))

  select case(hm%d%ispin)
  case(UNPOLARIZED)
    do jst = 1, hm%oct_st%nst
      pot = M_ZERO
      forall (ip = 1:der%mesh%np)
        rho(ip) = hm%oct_st%occ(jst, 1) * R_AIMAG(R_CONJ(hm%oct_st%X(psi)(ip, 1, jst, ik))*psi(ip, 1))
      end forall
      call dpoisson_solve(psolver, pot, rho)
      forall(ip = 1:der%mesh%np)
        hpsi(ip, 1) = hpsi(ip, 1) + M_TWO*M_zI * &
          hm%oct_st%X(psi)(ip, 1, jst, ik)*(pot(ip) + hm%oct_fxc(ip, 1, 1)*rho(ip))
      end forall
    end do 

  case(SPIN_POLARIZED)
    do jst = 1, hm%oct_st%nst
      if(hm%oct_st%occ(jst, ik) <= M_ZERO) cycle
      pot = M_ZERO
      do ip = 1, der%mesh%np
        rho(ip) = R_AIMAG(R_CONJ(hm%oct_st%X(psi)(ip, 1, jst, ik))*psi(ip, 1))
      end do
      call dpoisson_solve(psolver, pot, rho)
      do ip = 1, der%mesh%np
        hpsi(ip, 1) = hpsi(ip, 1) +  M_TWO*M_zI*hm%oct_st%occ(ip, ik)* &
          hm%oct_st%X(psi)(ip, 1, jst, ik)*pot(ip)
      end do
    end do 

  case(SPINORS)
    call messages_not_implemented("Function oct_exchange_operator for spinors")
  end select

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  POP_SUB(X(oct_exchange_operator))
end subroutine X(oct_exchange_operator)


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

  call X(vborders)(der, hm, psi, hpsi)

  SAFE_DEALLOCATE_A(auxpsi)
  SAFE_DEALLOCATE_A(aux2psi)
  POP_SUB(X(magnus))
end subroutine X(magnus)


! ---------------------------------------------------------
subroutine X(vborders) (der, hm, psi, hpsi)
  type(derivatives_t), intent(in)    :: der
  type(hamiltonian_t), intent(in)    :: hm
  R_TYPE,              intent(in)    :: psi(:,:)
  R_TYPE,              intent(inout) :: hpsi(:,:)

  integer :: idim

  PUSH_SUB(X(vborders))

  if(hm%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, hm%d%dim
      hpsi(1:der%mesh%np, idim) = hpsi(1:der%mesh%np, idim) + &
        M_zI*hm%ab_pot(1:der%mesh%np)*psi(1:der%mesh%np, idim)
    end do
  end if

  POP_SUB(X(vborders))
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(h_mgga_terms) (hm, der, psi, hpsi, ik, grad)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: psi(:,:)
  R_TYPE,              intent(inout) :: hpsi(:,:)
  integer,             intent(in)    :: ik
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: ispace, idim, ispin
  R_TYPE, allocatable :: cgrad(:,:), diverg(:)

  PUSH_SUB(X(h_mgga_terms))

  call X(get_grad)(hm, der, psi, grad)
  ispin = states_dim_get_spin_index(hm%d, ik)

  SAFE_ALLOCATE(cgrad(1:der%mesh%np_part, 1:MAX_DIM))
  SAFE_ALLOCATE(diverg(1:der%mesh%np))

  do idim = 1, hm%d%dim
    do ispace = 1, der%mesh%sb%dim
      cgrad(1:der%mesh%np, ispace) = M_TWO*grad(1:der%mesh%np, ispace, idim)*hm%vtau(1:der%mesh%np, ispin)
    end do

    diverg(1:der%mesh%np) = M_ZERO
    call X(derivatives_div)(der, cgrad, diverg)

    hpsi(1:der%mesh%np, idim) = hpsi(1:der%mesh%np, idim) - diverg(1:der%mesh%np)
  end do

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

  PUSH_SUB(X(vmask))

  if(hm%ab == MASK_ABSORBING) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
           st%X(psi)(1:gr%mesh%np, idim, ist, ik) = st%X(psi)(1:gr%mesh%np, idim, ist, ik)*(M_ONE - hm%ab_pot(1:gr%mesh%np))
        end do
      end do
    end do
  end if

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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
