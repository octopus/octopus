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
subroutine X(hamiltonian_apply_batch) (hm, gr, psib, hpsib, ik, t, kinetic_only)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(batch_t),       intent(inout) :: psib
  type(batch_t),       intent(inout) :: hpsib
  integer,             intent(in)    :: ik
  FLOAT, optional,     intent(in)    :: t
  logical, optional,   intent(in)    :: kinetic_only
  
  integer :: nst, bs, sp
  R_TYPE, pointer     :: epsi(:,:)
  R_TYPE, allocatable :: lapl(:, :, :)
  R_TYPE, pointer     :: grad(:, :, :)
  R_TYPE, allocatable :: psi_copy(:, :, :)

  type(profile_t), save :: phase_prof
  logical :: kinetic_only_, apply_kpoint
  integer :: ii, ist, idim, ip
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)
  type(batch_t) :: epsib, laplb
  type(der_handle_t), allocatable :: handles(:, :)
  type(profile_t), save :: prof

  call profiling_in(prof, "HAMILTONIAN")
  call push_sub('hamiltonian_inc.Xhpsi_batch')

  kinetic_only_ = .false.
  if(present(kinetic_only)) kinetic_only_ = kinetic_only

  if(present(t).and.hm%d%cdft) then
    message(1) = "TDCDFT not yet implemented"
    call write_fatal(1)
  end if
  
  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)

  nst = psib%nst
  
  ALLOCATE(lapl(1:gr%mesh%np, 1:hm%d%dim, 1:nst), gr%mesh%np*hm%d%dim*nst)
  call batch_init(laplb, hm%d%dim, psib%states(1)%ist, psib%states(nst)%ist, lapl)

  apply_kpoint = simul_box_is_periodic(gr%sb) .and. .not. kpoint_is_gamma(hm%d, ik)

  if(apply_kpoint) then
    ALLOCATE(psi_copy(1:gr%mesh%np_part, 1:hm%d%dim, 1:nst), gr%mesh%np_part*hm%d%dim*nst)
    call batch_init(epsib, hm%d%dim, psib%states(1)%ist, psib%states(nst)%ist, psi_copy)
  else
    call batch_copy(psib, epsib)
  end if

  call X(set_bc_batch)(gr%der, psib)


  bs = hardware%X(block_size)

  do sp = 1, gr%mesh%np_part, bs
    do ii = 1, nst
      call set_pointers()
      
      if(apply_kpoint) then ! we copy psi to epsi applying the exp(i k.r) phase
        call profiling_in(phase_prof, "PBC_PHASE_APPLY")
      
        forall (idim = 1:hm%d%dim, ip = sp:min(sp + bs - 1, gr%mesh%np_part)) 
          psi_copy(ip, idim, ii) = hm%phase(ip, ik)*psi(ip, idim)
        end forall
        
        call profiling_out(phase_prof)
      end if
      
    end do
  end do

  ! start the calculation of the laplacian
  ALLOCATE(handles(hm%d%dim, nst), hm%d%dim*nst)

  do ii = 1, nst
    do idim = 1, hm%d%dim
      call der_handle_init(handles(idim, ii), gr%der)
    end do
  end do

  call X(derivatives_lapl_batch_start)(gr%der, handles, epsib, laplb, set_bc = .false.)
    
  if (.not. kinetic_only_) then
    ! apply the potential
    call X(vlpsi_batch)(hm, gr%mesh, epsib, hpsib, ik)
    if(hm%ep%non_local) call X(vnlpsi_batch)(hm, gr%mesh, epsib, hpsib, ik)
  end if

  call X(derivatives_lapl_batch_finish)(gr%der, handles, epsib, laplb)

  do ii = 1, nst
    call set_pointers()

    if (kinetic_only_) hpsi(1:gr%mesh%np, 1:hm%d%dim) = M_ZERO

    ! finish the calculation of the laplacian
    call profiling_in(C_PROFILING_KINETIC)
    do idim = 1, hm%d%dim
      call lalg_axpy(gr%mesh%np, -M_HALF/hm%mass, lapl(:, idim, ii), hpsi(:, idim))
      call der_handle_end(handles(idim, ii))
    end do
    call profiling_out(C_PROFILING_KINETIC)
    
    if (.not. kinetic_only_) then
      
      ! all functions that require the gradient or other derivatives of
      ! epsi should go after this point and calculate it using Xget_grad

      ! FIX ME: here not all the possible couplings between vector
      ! potentials are considered

      nullify(grad)
      
      if (present(t)) then
        ! lasers
        if(hm%ep%no_lasers > 0) then
          if(any(laser_requires_gradient(hm%ep%lasers(1:hm%ep%no_lasers)))) then
            call X(get_grad)(hm, gr, epsi, grad)
          end if
          if(associated(hm%ep%a_static)) then
            call X(vlasers)(hm%ep%lasers, hm%ep%no_lasers, gr, hm%d, epsi, hpsi, grad, &
              ik, hm%ep%gyromagnetic_ratio, hm%ep%a_static, t)
          else
            call X(vlasers)(hm%ep%lasers, hm%ep%no_lasers, gr, hm%d, epsi, hpsi, grad, &
              ik, hm%ep%gyromagnetic_ratio, t = t)
          end if
        end if
#ifdef R_TCOMPLEX
        if (gauge_field_is_applied(hm%ep%gfield)) then
          call X(get_grad)(hm, gr, psi, grad)
          call gauge_field_apply(hm%ep%gfield, gr, hm%d%dim, epsi, grad, hpsi)
        end if
#endif
      end if
      
      if(hm%theory_level == HARTREE .or. hm%theory_level == HARTREE_FOCK) &
        call X(exchange_operator)(hm, gr, epsi, hpsi, ist, ik)
      
      if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) &
        call X(h_mgga_terms) (hm, gr, epsi, hpsi, ik, grad)

      if(hamiltonian_oct_exchange(hm)) call X(oct_exchange_operator)(hm, gr, epsi, hpsi, ik)
      
      call X(magnetic_terms) (gr, hm, epsi, hpsi, grad, ik)
      
      if(present(t)) call X(vborders) (gr, hm, epsi, hpsi)

      SAFE_DEALLOCATE_P(grad)
   
      if(hm%scissor%apply) call X(scissor_apply)(hm%scissor, gr%mesh, ik, epsi, hpsi)

    end if
  end do

  SAFE_DEALLOCATE_A(handles)

  if(apply_kpoint) then
    ! now we need to remove the exp(-i k.r) factor
    do sp = 1, gr%mesh%np, bs
      do ii = 1, nst
        call set_pointers()
        
        call profiling_in(phase_prof)
        
        forall (idim = 1:hm%d%dim, ip = sp:min(sp + bs - 1, gr%mesh%np))
          hpsi(ip, idim) = conjg(hm%phase(ip, ik))*hpsi(ip, idim)
        end forall

        call profiling_out(phase_prof)
        
      end do
    end do
  end if
  
  call batch_end(epsib)
  call batch_end(laplb)

  SAFE_DEALLOCATE_A(lapl)

  call pop_sub()
  call profiling_out(prof)

contains

  subroutine set_pointers
    psi  => psib%states(ii)%X(psi)
    epsi => epsib%states(ii)%X(psi)
    hpsi => hpsib%states(ii)%X(psi)
    ist  =  psib%states(ii)%ist
  end subroutine set_pointers

end subroutine X(hamiltonian_apply_batch)

! ---------------------------------------------------------

subroutine X(get_grad)(hm, gr, psi, grad)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:, :)
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: idim

  if( .not. associated(grad)) then
    ALLOCATE(grad(1:gr%mesh%np, 1:MAX_DIM, 1:hm%d%dim), gr%mesh%np*MAX_DIM*hm%d%dim)
    do idim = 1, hm%d%dim 
      ! boundary points were already set by the Laplacian
      call X(derivatives_grad)(gr%der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
  end if
  
end subroutine X(get_grad)       

! ---------------------------------------------------------
subroutine X(hamiltonian_apply) (hm, gr, psi, hpsi, ist, ik, t, kinetic_only)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ist       ! the index of the state
  integer,             intent(in)    :: ik        ! the index of the k-point
  R_TYPE, target,      intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  FLOAT, optional,     intent(in)    :: t
  logical, optional,   intent(in)    :: kinetic_only

  logical :: kinetic_only_
  type(batch_t) :: psib, hpsib

  kinetic_only_ = .false.
  if(present(kinetic_only)) kinetic_only_ = kinetic_only

  call batch_init(psib, hm%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, hm%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  if(present(t)) then
    call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik, t, kinetic_only = kinetic_only_)
  else
    call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik, kinetic_only = kinetic_only_)
  end if

  call batch_end(psib)
  call batch_end(hpsib)

end subroutine X(hamiltonian_apply)

! ---------------------------------------------------------
subroutine X(exchange_operator) (hm, gr, psi, hpsi, ist, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  integer,             intent(in)    :: ist       ! the index of the state
  integer,             intent(in)    :: ik        ! the index of the k-point
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('hamiltonian_inc.Xexchange_operator')

  ALLOCATE(rho(gr%mesh%np), gr%mesh%np)
  ALLOCATE(pot(gr%mesh%np), gr%mesh%np)

  ! WARNING: this can be very condensed
  select case(hm%d%ispin)
  case(UNPOLARIZED)
    do j = 1, hm%st%nst
      if(hm%st%occ(j, ik) <= M_ZERO) cycle

      ! in Hartree we just remove the self-interaction
      if(hm%theory_level == HARTREE .and. j .ne. ist) cycle

      pot = M_ZERO
      do k = 1, gr%mesh%np
        rho(k) = R_CONJ(hm%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%mesh%np
        hpsi(k, 1) = hpsi(k, 1) - hm%exx_coef * (hm%st%occ(j, ik)/M_TWO) * hm%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPIN_POLARIZED)
    do j = 1, hm%st%nst
      if(hm%st%occ(j, ik) <= M_ZERO) cycle

      ! in Hartree we just remove the self-interaction
      if(hm%theory_level == HARTREE .and. j .ne. ist) cycle

      pot = M_ZERO
      do k = 1, gr%mesh%np
        rho(k) = R_CONJ(hm%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%mesh%np
        hpsi(k, 1) = hpsi(k, 1) - hm%exx_coef * hm%st%occ(j, ik) * hm%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPINORS)

  end select

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  call pop_sub()
end subroutine X(exchange_operator)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator) (hm, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('hamiltonian_inc.Xexchange_operator')

  ALLOCATE(rho(gr%mesh%np), gr%mesh%np)
  ALLOCATE(pot(gr%mesh%np), gr%mesh%np)

  select case(hm%d%ispin)
  case(UNPOLARIZED)
    do j = 1, hm%oct_st%nst
      pot = M_ZERO
      do k = 1, gr%mesh%np
        rho(k) = R_CONJ(hm%oct_st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%mesh%np
        hpsi(k, 1) = hpsi(k, 1) +  M_TWO * M_z1 * (hm%oct_st%occ(j, ik)/M_TWO) * &
                     hm%oct_st%X(psi)(k, 1, j, ik) * R_AIMAG(pot(k))
      end do
    end do 

  case(SPIN_POLARIZED)
    do j = 1, hm%oct_st%nst
      if(hm%oct_st%occ(j, ik) <= M_ZERO) cycle
      pot = M_ZERO
      do k = 1, gr%mesh%np
        rho(k) = R_CONJ(hm%oct_st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%mesh%np
        hpsi(k, 1) = hpsi(k, 1) +  M_TWO * M_z1 * hm%oct_st%occ(j, ik) * &
          hm%oct_st%X(psi)(k, 1, j, ik)*R_AIMAG(pot(k))
      end do
    end do 

  case(SPINORS)

  end select

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(pot)
  call pop_sub()
end subroutine X(oct_exchange_operator)


! ---------------------------------------------------------
subroutine X(magnus) (hm, gr, psi, hpsi, ik, vmagnus)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(gr%mesh%np, hm%d%dim)
  FLOAT,               intent(in)    :: vmagnus(gr%mesh%np, hm%d%nspin, 2)

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim, ispin

  ! We will assume, for the moment, no spinors.

  call push_sub('hamiltonian_inc.Xmagnus')

  ALLOCATE( auxpsi(gr%mesh%np_part, hm%d%dim), gr%mesh%np_part*hm%d%dim)
  ALLOCATE(aux2psi(gr%mesh%np,      hm%d%dim), gr%mesh%np*hm%d%dim)

  ispin = states_dim_get_spin_index(hm%d, ik)

  call X(hamiltonian_apply)(hm, gr, psi, hpsi, ist = 1, ik = ik, kinetic_only = .true.)

  do idim = 1, hm%d%dim
    call lalg_copy(gr%mesh%np, hpsi(:, idim), auxpsi(:, idim))
  end do

  if (hm%ep%non_local) call X(vnlpsi)(hm, gr%mesh, psi, auxpsi, ik)

  hpsi(1:gr%mesh%np, 1) = hpsi(1:gr%mesh%np, 1) -  M_zI*vmagnus(1:gr%mesh%np, ispin, 1)*auxpsi(1:gr%mesh%np, 1)
  auxpsi(1:gr%mesh%np, 1) = vmagnus(1:gr%mesh%np, ispin, 1)*psi(1:gr%mesh%np, 1)

  call X(hamiltonian_apply)(hm, gr, auxpsi, aux2psi, ist = 1, ik = ik, kinetic_only = .true.)

  if (hm%ep%non_local) call X(vnlpsi)(hm, gr%mesh, auxpsi, aux2psi, ik)

  hpsi(1:gr%mesh%np, 1) = hpsi(1:gr%mesh%np, 1) + M_zI*aux2psi(1:gr%mesh%np, 1)

  do idim = 1, hm%d%dim
    hpsi(1:gr%mesh%np, idim) = hpsi(1:gr%mesh%np, idim) + hm%ep%Vpsl(1:gr%mesh%np)*psi(1:gr%mesh%np,idim)
  end do

  hpsi(1:gr%mesh%np, 1) = hpsi(1:gr%mesh%np, 1) + vmagnus(1:gr%mesh%np, ispin, 2)*psi(1:gr%mesh%np, 1)

  if (hm%ep%non_local) call X(vnlpsi)(hm, gr%mesh, psi, Hpsi, ik)

  call X(vborders) (gr, hm, psi, hpsi)

  SAFE_DEALLOCATE_A(auxpsi)
  SAFE_DEALLOCATE_A(aux2psi)
  call pop_sub()
end subroutine X(magnus)

! ---------------------------------------------------------
! Here we the terms arising from the presence of a possible static external
! magnetic field, and the terms that come from CDFT.
subroutine X(magnetic_terms) (gr, hm, psi, hpsi, grad, ik)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: hm
  R_TYPE,              intent(inout) :: psi(:, :)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:, :) ! hpsi(gr%mesh%np, hm%d%dim)
  R_TYPE,              pointer       :: grad(:, :, :)
  integer,             intent(in)    :: ik

  integer :: k, idim
  FLOAT,  allocatable :: div(:), tmp(:,:)
  R_TYPE, allocatable :: lhpsi(:, :)

  call push_sub('hamiltonian_inc.Xmagnetic_terms')

  if(hm%d%cdft .or. associated(hm%ep%A_static)) then
    call X(get_grad)(hm, gr, psi, grad)
  else
    call pop_sub(); return
  endif

  ! If we are using CDFT:
  if(hm%d%cdft) then

    ALLOCATE(div(gr%mesh%np), gr%mesh%np)
    ALLOCATE(tmp(gr%mesh%np_part, gr%mesh%sb%dim), gr%mesh%np_part*gr%mesh%sb%dim)
    select case (hm%d%ispin)
    case(UNPOLARIZED)
      tmp(1:gr%mesh%np, :) = hm%axc(1:gr%mesh%np, :, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        tmp(1:gr%mesh%np, :) = hm%axc(1:gr%mesh%np, :, 1)
      else
        tmp(1:gr%mesh%np, :) = hm%axc(1:gr%mesh%np, :, 2)
      end if
    case(SPINORS)
      write(message(1),'(a)') 'Current DFT not yet functional in spinors mode, sorry.'
      call write_fatal(2)
    end select
    call dderivatives_div(gr%der, tmp, div)
    hpsi(1:gr%mesh%np, 1) = hpsi(1:gr%mesh%np, 1) - M_HALF*M_zI*div*psi(1:gr%mesh%np, 1)
    SAFE_DEALLOCATE_A(div)
    SAFE_DEALLOCATE_A(tmp)

    select case (hm%d%ispin)
    case(UNPOLARIZED)
      do k = 1, gr%mesh%np
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(hm%axc(k, 1:gr%mesh%sb%dim, 1), grad(k, 1:gr%mesh%sb%dim, 1))
      end do
    case(SPIN_POLARIZED)
      do k = 1, gr%mesh%np
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) - &
               M_zI*dot_product(hm%axc(k, 1:gr%mesh%sb%dim, 1), grad(k, 1:gr%mesh%sb%dim, 1))
        else
          hpsi(k, 1) = hpsi(k, 1) - &
               M_zI*dot_product(hm%axc(k, 1:gr%mesh%sb%dim, 2), grad(k, 1:gr%mesh%sb%dim, 1))
        end if
      end do
    case(SPINORS)
      ! Not yet implemented
    end select

  endif !CDT

  ! If we have an external magnetic field
  if (associated(hm%ep%A_static)) then
    do k = 1, gr%mesh%np
      hpsi(k, :) = hpsi(k, :) + &
           M_HALF*dot_product(hm%ep%A_static(k, 1:gr%mesh%sb%dim), hm%ep%A_static(k, 1:gr%mesh%sb%dim))*psi(k, :)
      select case(hm%d%ispin)
      case(UNPOLARIZED, SPIN_POLARIZED)
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(hm%ep%A_static(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, 1))
      case (SPINORS)
        do idim = 1, hm%d%dim
          hpsi(k, idim) = hpsi(k, idim) - &
               M_zI*dot_product(hm%ep%A_static(k, 1:gr%mesh%sb%dim), grad(k, 1:gr%mesh%sb%dim, idim))
        end do
      end select
    end do
  end if

  ! Zeeman term
  if (associated(hm%ep%B_field) .and. hm%d%ispin /= UNPOLARIZED) then
    ALLOCATE(lhpsi(gr%mesh%np, hm%d%dim), gr%mesh%np*hm%d%dim)
    select case (hm%d%ispin)
    case (SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:gr%mesh%np, 1) = - M_HALF/P_C*sqrt(dot_product(hm%ep%B_field, hm%ep%B_field))*psi(1:gr%mesh%np, 1)
      else
        lhpsi(1:gr%mesh%np, 1) = + M_HALF/P_C*sqrt(dot_product(hm%ep%B_field, hm%ep%B_field))*psi(1:gr%mesh%np, 1)
      end if
    case (SPINORS)
      lhpsi(1:gr%mesh%np, 1) = M_HALF/P_C*( hm%ep%B_field(3)*psi(1:gr%mesh%np, 1) &
                                 + (hm%ep%B_field(1) - M_zI*hm%ep%B_field(2))*psi(1:gr%mesh%np, 2))
      lhpsi(1:gr%mesh%np, 2) = M_HALF/P_C*(-hm%ep%B_field(3)*psi(1:gr%mesh%np, 2) &
                                 + (hm%ep%B_field(1) + M_zI*hm%ep%B_field(2))*psi(1:gr%mesh%np, 1))
    end select
    hpsi(1:gr%mesh%np, :) = hpsi(1:gr%mesh%np, :) + (hm%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:gr%mesh%np, :)
    SAFE_DEALLOCATE_A(lhpsi)
  end if

  call pop_sub()
end subroutine X(magnetic_terms)


! ---------------------------------------------------------
subroutine X(vnlpsi) (hm, mesh, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(mesh_t),        intent(inout) :: mesh
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np_part, hm%d%dim)
  integer,             intent(in)    :: ik

  call profiling_in(C_PROFILING_VNLPSI)
  call push_sub('hamiltonian_inc.Xvnlpsi')

  call X(project_psi)(mesh, hm%ep%proj, hm%ep%natoms, hm%d%dim, psi, hpsi, ik)
    
  call pop_sub()
  call profiling_out(C_PROFILING_VNLPSI)
end subroutine X(vnlpsi)


! ---------------------------------------------------------
subroutine X(vnlpsi_batch) (hm, mesh, psib, hpsib, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(mesh_t),        intent(in)    :: mesh
  type(batch_t),       intent(in)    :: psib
  type(batch_t),       intent(inout) :: hpsib
  integer,             intent(in)    :: ik

  integer :: ii
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)

  call profiling_in(C_PROFILING_VNLPSI)
  call push_sub('hamiltonian_inc.Xvnlpsi_batch')

  do ii = 1, psib%nst
    psi  => psib%states(ii)%X(psi)
    hpsi => hpsib%states(ii)%X(psi)

    call X(project_psi)(mesh, hm%ep%proj, hm%ep%natoms, hm%d%dim, psi, hpsi, ik)
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_VNLPSI)
end subroutine X(vnlpsi_batch)


! ---------------------------------------------------------
subroutine X(vlpsi_batch) (hm, m, psib, hpsib, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(mesh_t),        intent(in)    :: m
  integer,             intent(in)    :: ik
  type(batch_t),       intent(in)    :: psib
  type(batch_t),       intent(inout) :: hpsib

  integer :: ip, ii, ispin
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('hamiltonian_inc.Xvlpsi_batch')

  select case(hm%d%ispin)
  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_dim_get_spin_index(hm%d, ik)
    call X(em_field_apply_batch)(hm%total(ispin), m, psib, hpsib)

  case(SPINORS)
    !the spinor case is more complicated since it mixes the two components.
    do ii = 1, psib%nst
      psi  => psib%states(ii)%X(psi)
      hpsi => hpsib%states(ii)%X(psi)

      do ip = 1, m%np
        hpsi(ip, 1) = (hm%vhxc(ip, 1) + hm%ep%vpsl(ip))*psi(ip, 1) + (hm%vhxc(ip, 3) + M_zI*hm%vhxc(ip, 4))*psi(ip, 2)
        hpsi(ip, 2) = (hm%vhxc(ip, 2) + hm%ep%vpsl(ip))*psi(ip, 2) + (hm%vhxc(ip, 3) - M_zI*hm%vhxc(ip, 4))*psi(ip, 1)
      end do

    end do
    call profiling_count_operations((6*R_ADD + 2*R_MUL)*m%np*psib%nst)

  end select

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vlpsi_batch)

! ---------------------------------------------------------
subroutine X(vexternal) (hm, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np_part, hm%d%dim)
  integer,             intent(in)    :: ik

  integer :: idim

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('hamiltonian_inc.Xvlpsi')

  do idim = 1, hm%d%dim
    hpsi(1:gr%mesh%np, idim) = hpsi(1:gr%mesh%np, idim) + hm%ep%vpsl(1:gr%mesh%np)*psi(1:gr%mesh%np, idim)
  end do

  if(hm%ep%non_local) call X(vnlpsi)(hm, gr%mesh, psi, hpsi, ik)

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vexternal)

! ---------------------------------------------------------
subroutine X(vborders) (gr, hm, psi, hpsi)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: hm
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np_part, hm%d%dim)

  integer :: idim

  call push_sub('hamiltonian_inc.Xvborders')

  if(hm%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, hm%d%dim
      hpsi(1:gr%mesh%np, idim) = hpsi(1:gr%mesh%np, idim) + M_zI*hm%ab_pot(1:gr%mesh%np)*psi(1:gr%mesh%np, idim)
    end do
  end if

  call pop_sub()
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(h_mgga_terms) (hm, gr, psi, hpsi, ik, grad)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(gr%mesh%np_part, hm%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(gr%mesh%np_part, hm%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: ispace, idim, ispin
  R_TYPE, allocatable :: cgrad(:,:), diverg(:)

  call push_sub('hamiltonian_inc.Xh_mgga_terms')

  call X(get_grad)(hm, gr, psi, grad)
  ispin = states_dim_get_spin_index(hm%d, ik)

  ALLOCATE(cgrad(1:gr%mesh%np_part, 1:MAX_DIM), gr%mesh%np_part*MAX_DIM)
  ALLOCATE(diverg(1:gr%mesh%np), gr%mesh%np)

  do idim = 1, hm%d%dim
    do ispace = 1, gr%sb%dim
      cgrad(1:gr%mesh%np, ispace) = M_TWO*grad(1:gr%mesh%np, ispace, idim)*hm%vtau(1:gr%mesh%np, ispin)
      call X(set_bc)(gr%der, cgrad(:, ispace))
    end do

    diverg(1:gr%mesh%np) = M_ZERO
    call X(derivatives_div)(gr%der, cgrad, diverg)

    hpsi(1:gr%mesh%np, idim) = hpsi(1:gr%mesh%np, idim) - diverg(1:gr%mesh%np)
  end do

  SAFE_DEALLOCATE_A(cgrad)
  SAFE_DEALLOCATE_A(diverg)

  call pop_sub()
end subroutine X(h_mgga_terms)


! ---------------------------------------------------------
subroutine X(vmask) (gr, hm, st)
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim

  call push_sub('hamiltonian_inc.Xvmask')

  if(hm%ab == MASK_ABSORBING) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
           st%X(psi)(1:gr%mesh%np, idim, ist, ik) = st%X(psi)(1:gr%mesh%np, idim, ist, ik)*(M_ONE - hm%ab_pot(1:gr%mesh%np))
        end do
      end do
    end do
  end if

  call pop_sub()
end subroutine X(vmask)

! ---------------------------------------------------------
subroutine X(hamiltonian_diagonal) (hm, gr, diag, ik)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(out)   :: diag(:,:) ! hpsi(gr%mesh%np, hm%d%dim)

  integer :: idim, ip, ispin

  FLOAT, allocatable  :: ldiag(:)

  call push_sub('hamiltonian_inc.Xhpsi_diag')
  
  ALLOCATE(ldiag(gr%mesh%np), gr%mesh%np)

  diag = M_ZERO

  call derivatives_lapl_diag(gr%der, ldiag)

! MJV: this is where the MASS enters the kinetic energy.
  do idim = 1, hm%d%dim
    diag(1:gr%mesh%np, idim) = -M_HALF/hm%mass*ldiag(1:gr%mesh%np)
  end do

  select case(hm%d%ispin)

  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_dim_get_spin_index(hm%d, ik)
    diag(1:gr%mesh%np, 1) = diag(1:gr%mesh%np, 1) + hm%vhxc(1:gr%mesh%np, ispin) + hm%ep%vpsl(1:gr%mesh%np)
    
  case(SPINORS)
    do ip = 1, gr%mesh%np
      diag(ip, 1) = diag(ip, 1) + hm%vhxc(ip, 1) + hm%ep%vpsl(ip)
      diag(ip, 2) = diag(ip, 2) + hm%vhxc(ip, 2) + hm%ep%vpsl(ip)
    end do
    
  end select
    
  call pop_sub()
end subroutine X(hamiltonian_diagonal)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
