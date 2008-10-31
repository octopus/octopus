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
subroutine X(hpsi_batch) (h, gr, psib, hpsib, ik, t, kinetic_only)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  type(batch_t),       intent(inout) :: psib
  type(batch_t),       intent(inout) :: hpsib
  integer,             intent(in)    :: ik
  FLOAT, optional,     intent(in)    :: t
  logical, optional,   intent(in)    :: kinetic_only
  
  integer :: nst
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


  call profiling_in(C_PROFILING_HPSI)
  call push_sub('hamiltonian_inc.Xhpsi_batch')

  kinetic_only_ = .false.
  if(present(kinetic_only)) kinetic_only_ = kinetic_only

  if(present(t).and.h%d%cdft) then
    message(1) = "TDCDFT not yet implemented"
    call write_fatal(1)
  end if
  
  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)

  nst = psib%nst
  
  ALLOCATE(lapl(1:NP, 1:h%d%dim, 1:nst), NP*h%d%dim*nst)
  call batch_init(laplb, h%d%dim, psib%states(1)%ist, psib%states(nst)%ist, lapl)

  apply_kpoint = simul_box_is_periodic(gr%sb) .and. .not. kpoint_is_gamma(h%d, ik)

  if(apply_kpoint) then
    ALLOCATE(psi_copy(1:NP_PART, 1:h%d%dim, 1:nst), NP_PART*h%d%dim*nst)
    call batch_init(epsib, h%d%dim, psib%states(1)%ist, psib%states(nst)%ist, psi_copy)
  else
    call batch_copy(psib, epsib)
  end if

  call X(set_bc_batch)(gr%der, psib)

  do ii = 1, nst
    call set_pointers()

    if(apply_kpoint) then ! we copy psi to epsi applying the exp(i k.r) phase
      call profiling_in(phase_prof, "PBC_PHASE_APPLY")
      
      forall (idim = 1:h%d%dim, ip = 1:NP_PART)
        psi_copy(ip, idim, ii) = h%phase(ip, ik)*psi(ip, idim)
      end forall
      
      call profiling_out(phase_prof)
    end if

  end do

  ! start the calculation of the laplacian
  ALLOCATE(handles(h%d%dim, nst), h%d%dim*nst)

  do ii = 1, nst
    do idim = 1, h%d%dim
      call der_handle_init(handles(idim, ii), gr%der)
    end do
  end do

  call X(derivatives_lapl_batch_start)(gr%der, handles, epsib, laplb, set_bc = .false.)
    
  if (.not. kinetic_only_) then
    ! apply the potential
    call X(vlpsi_batch)(h, gr%m, epsib, hpsib, ik)
    if(h%ep%non_local) call X(vnlpsi_batch)(h, gr%m, epsib, hpsib, ik)
  end if

  call X(derivatives_lapl_batch_finish)(gr%der, handles, epsib, laplb)

  do ii = 1, nst
    call set_pointers()

    if (kinetic_only_) hpsi(1:gr%m%np, 1:h%d%dim) = M_ZERO

    ! finish the calculation of the laplacian
    call profiling_in(C_PROFILING_KINETIC)
    do idim = 1, h%d%dim
      call lalg_axpy(NP, -M_HALF/h%mass, lapl(:, idim, ii), hpsi(:, idim))
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
        if(h%ep%no_lasers > 0) then
          if(any(laser_requires_gradient(h%ep%lasers(1:h%ep%no_lasers)))) call X(get_grad)(h, gr, epsi, grad)
          call X(vlasers)(h%ep%lasers, h%ep%no_lasers, gr, h%d, epsi, hpsi, grad, ik, t, h%ep%gyromagnetic_ratio, h%ep%a_static)
        end if

        if (gauge_field_is_applied(h%ep%gfield)) call X(vgauge)(gr, h, epsi, hpsi, grad)
      end if
      
      if(h%theory_level == HARTREE .or. h%theory_level == HARTREE_FOCK) &
        call X(exchange_operator)(h, gr, epsi, hpsi, ist, ik)
      
      if(iand(h%xc_family, XC_FAMILY_MGGA).ne.0) &
        call X(h_mgga_terms) (h, gr, epsi, hpsi, ik, grad)

      if(hamiltonian_oct_exchange(h)) call X(oct_exchange_operator)(h, gr, epsi, hpsi, ik)
      
      call X(magnetic_terms) (gr, h, epsi, hpsi, grad, ik)
      
      if(present(t)) call X(vborders) (gr, h, epsi, hpsi)

      if(associated(grad)) deallocate(grad)
   
    end if
    
    if(apply_kpoint) then
      ! now we need to remove the exp(-i k.r) factor
      call profiling_in(phase_prof)

      forall (idim = 1:h%d%dim, ip = 1:NP)
        hpsi(ip, idim) = conjg(h%phase(ip, ik))*hpsi(ip, idim)
      end forall
  
      call profiling_out(phase_prof)
    end if

  end do

  call batch_end(epsib)
  call batch_end(laplb)

  call pop_sub()
  call profiling_out(C_PROFILING_HPSI)

contains

  subroutine set_pointers
    psi  => psib%states(ii)%X(psi)
    epsi => epsib%states(ii)%X(psi)
    hpsi => hpsib%states(ii)%X(psi)
    ist  =  psib%states(ii)%ist
  end subroutine set_pointers

end subroutine X(hpsi_batch)

! ---------------------------------------------------------

subroutine X(get_grad)(h, gr, psi, grad)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:, :)
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: idim

  if( .not. associated(grad)) then
    ALLOCATE(grad(1:NP, 1:MAX_DIM, 1:h%d%dim), NP*MAX_DIM*h%d%dim)
    do idim = 1, h%d%dim 
      ! boundary points were already set by the Laplacian
      call X(derivatives_grad)(gr%der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
  end if
  
end subroutine X(get_grad)       

! ---------------------------------------------------------
subroutine X(hpsi) (h, gr, psi, hpsi, ist, ik, t, kinetic_only)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ist       ! the index of the state
  integer,             intent(in)    :: ik        ! the index of the k-point
  R_TYPE, target,      intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  FLOAT, optional,     intent(in)    :: t
  logical, optional,   intent(in)    :: kinetic_only

  logical :: kinetic_only_
  type(batch_t) :: psib, hpsib

  kinetic_only_ = .false.
  if(present(kinetic_only)) kinetic_only_ = kinetic_only

  call batch_init(psib, h%d%dim, 1)
  call batch_add_state(psib, ist, psi)
  call batch_init(hpsib, h%d%dim, 1)
  call batch_add_state(hpsib, ist, hpsi)

  if(present(t)) then
    call X(hpsi_batch)(h, gr, psib, hpsib, ik, t, kinetic_only = kinetic_only_)
  else
    call X(hpsi_batch)(h, gr, psib, hpsib, ik, kinetic_only = kinetic_only_)
  end if

  call batch_end(psib)
  call batch_end(hpsib)

end subroutine X(hpsi)

! ---------------------------------------------------------
subroutine X(exchange_operator) (h, gr, psi, hpsi, ist, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  integer,             intent(in)    :: ist       ! the index of the state
  integer,             intent(in)    :: ik        ! the index of the k-point
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('hamiltonian_inc.Xexchange_operator')

  ALLOCATE(rho(gr%m%np), gr%m%np)
  ALLOCATE(pot(gr%m%np), gr%m%np)

  ! WARNING: this can be very condensed
  select case(h%d%ispin)
  case(UNPOLARIZED)
    do j = 1, h%st%nst
      if(h%st%occ(j, ik) <= M_ZERO) cycle

      ! in Hartree we just remove the self-interaction
      if(h%theory_level == HARTREE .and. j .ne. ist) cycle

      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) - h%exx_coef * (h%st%occ(j, ik)/M_TWO) * h%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPIN_POLARIZED)
    do j = 1, h%st%nst
      if(h%st%occ(j, ik) <= M_ZERO) cycle

      ! in Hartree we just remove the self-interaction
      if(h%theory_level == HARTREE .and. j .ne. ist) cycle

      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) - h%exx_coef * h%st%occ(j, ik) * h%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPINORS)

  end select

  deallocate(rho, pot)
  call pop_sub()
end subroutine X(exchange_operator)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('hamiltonian_inc.Xexchange_operator')

  ALLOCATE(rho(gr%m%np), gr%m%np)
  ALLOCATE(pot(gr%m%np), gr%m%np)

  select case(h%d%ispin)
  case(UNPOLARIZED)
    do j = 1, h%oct_st%nst
      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%oct_st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) +  M_TWO * M_z1 * (h%oct_st%occ(j, ik)/M_TWO) * &
                     h%oct_st%X(psi)(k, 1, j, ik) * R_AIMAG(pot(k))
      end do
    end do 

  case(SPIN_POLARIZED)
    do j = 1, h%oct_st%nst
      if(h%oct_st%occ(j, ik) <= M_ZERO) cycle
      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%oct_st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) +  M_TWO * M_z1 * h%oct_st%occ(j, ik) * &
          h%oct_st%X(psi)(k, 1, j, ik)*R_AIMAG(pot(k))
      end do
    end do 

  case(SPINORS)

  end select

  deallocate(rho, pot)
  call pop_sub()
end subroutine X(oct_exchange_operator)


! ---------------------------------------------------------
subroutine X(magnus) (h, gr, psi, hpsi, ik, vmagnus)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  FLOAT,               intent(in)    :: vmagnus(NP, h%d%nspin, 2)

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim, ispin

  ! We will assume, for the moment, no spinors.

  call push_sub('hamiltonian_inc.Xmagnus')

  ALLOCATE( auxpsi(NP_PART, h%d%dim), NP_PART*h%d%dim)
  ALLOCATE(aux2psi(NP,      h%d%dim), NP*h%d%dim)

  ispin = states_dim_get_spin_index(h%d, ik)

  call X(hpsi)(h, gr, psi, hpsi, ist = 1, ik = ik, kinetic_only = .true.)

  do idim = 1, h%d%dim
    call lalg_copy(NP, hpsi(:, idim), auxpsi(:, idim))
  end do

  if (h%ep%non_local) call X(vnlpsi)(h, gr%m, psi, auxpsi, ik)

  hpsi(1:NP, 1) = hpsi(1:NP, 1) -  M_zI*vmagnus(1:NP, ispin, 1)*auxpsi(1:NP, 1)
  auxpsi(1:NP, 1) = vmagnus(1:NP, ispin, 1)*psi(1:NP, 1)

  call X(hpsi)(h, gr, auxpsi, aux2psi, ist = 1, ik = ik, kinetic_only = .true.)

  if (h%ep%non_local) call X(vnlpsi)(h, gr%m, auxpsi, aux2psi, ik)

  hpsi(1:NP, 1) = hpsi(1:NP, 1) + M_zI*aux2psi(1:NP, 1)

  do idim = 1, h%d%dim
    hpsi(1:NP, idim) = hpsi(1:NP, idim) + h%ep%Vpsl(1:NP)*psi(1:NP,idim)
  end do

  hpsi(1:NP, 1) = hpsi(1:NP, 1) + vmagnus(1:NP, ispin, 2)*psi(1:NP, 1)

  if (h%ep%non_local) call X(vnlpsi)(h, gr%m, psi, Hpsi, ik)

  call X(vborders) (gr, h, psi, hpsi)

  deallocate(auxpsi, aux2psi)
  call pop_sub()
end subroutine X(magnus)

! ---------------------------------------------------------
! Here we the terms arising from the presence of a possible static external
! magnetic field, and the terms that come from CDFT.
subroutine X(magnetic_terms) (gr, h, psi, hpsi, grad, ik)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:, :)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:, :) ! hpsi(NP, h%d%dim)
  R_TYPE,              pointer       :: grad(:, :, :)
  integer,             intent(in)    :: ik

  integer :: k, idim
  FLOAT,  allocatable :: div(:), tmp(:,:)
  R_TYPE, allocatable :: lhpsi(:, :)

  call push_sub('hamiltonian_inc.Xmagnetic_terms')

  if(h%d%cdft .or. associated(h%ep%A_static)) then
    call X(get_grad)(h, gr, psi, grad)
  else
    call pop_sub()
    return
  endif

  ! If we are using CDFT:
  if(h%d%cdft) then

    ALLOCATE(div(NP), NP)
    ALLOCATE(tmp(NP_PART, NDIM), NP_PART*NDIM)
    select case (h%d%ispin)
    case(UNPOLARIZED)
      tmp(1:NP, :) = h%axc(1:NP, :, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        tmp(1:NP, :) = h%axc(1:NP, :, 1)
      else
        tmp(1:NP, :) = h%axc(1:NP, :, 2)
      end if
    case(SPINORS)
      write(message(1),'(a)') 'Current DFT not yet functional in spinors mode, sorry.'
      call write_fatal(2)
    end select
    call dderivatives_div(gr%der, tmp, div)
    hpsi(1:NP, 1) = hpsi(1:NP, 1) - M_HALF*M_zI*div*psi(1:NP, 1)
    deallocate(div, tmp)

    select case (h%d%ispin)
    case(UNPOLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%axc(k, 1:NDIM, 1), grad(k, 1:NDIM, 1))
      end do
    case(SPIN_POLARIZED)
      do k = 1, NP
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) -  M_zI*dot_product(h%axc(k, 1:NDIM, 1), grad(k, 1:NDIM, 1))
        else
          hpsi(k, 1) = hpsi(k, 1) -  M_zI*dot_product(h%axc(k, 1:NDIM, 2), grad(k, 1:NDIM, 1))
        end if
      end do
    case(SPINORS)
      ! Not yet implemented
    end select

  endif !CDT

  ! If we have an external magnetic field
  if (associated(h%ep%A_static)) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(h%ep%A_static(k, 1:NDIM), h%ep%A_static(k, 1:NDIM))*psi(k, :)
      select case(h%d%ispin)
      case(UNPOLARIZED, SPIN_POLARIZED)
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%ep%A_static(k, 1:NDIM), grad(k, 1:NDIM, 1))
      case (SPINORS)
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(h%ep%A_static(k, 1:NDIM), grad(k, 1:NDIM, idim))
        end do
      end select
    end do
  end if

  ! Zeeman term
  if (associated(h%ep%B_field) .and. h%d%ispin /= UNPOLARIZED) then
    ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
    select case (h%d%ispin)
    case (SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(h%ep%B_field, h%ep%B_field))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(h%ep%B_field, h%ep%B_field))*psi(1:NP, 1)
      end if
    case (SPINORS)
      lhpsi(1:NP, 1) = M_HALF/P_C*( h%ep%B_field(3)*psi(1:NP, 1) &
                                 + (h%ep%B_field(1) - M_zI*h%ep%B_field(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-h%ep%B_field(3)*psi(1:NP, 2) &
                                 + (h%ep%B_field(1) + M_zI*h%ep%B_field(2))*psi(1:NP, 1))
    end select
    hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
    deallocate(lhpsi)
  end if

  call pop_sub()
end subroutine X(magnetic_terms)


! ---------------------------------------------------------
subroutine X(vnlpsi) (h, mesh, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(mesh_t),        intent(inout) :: mesh
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  call profiling_in(C_PROFILING_VNLPSI)
  call push_sub('hamiltonian_inc.Xvnlpsi')

  call X(project_psi)(mesh, h%ep%proj, h%ep%natoms, h%d%dim, psi, hpsi, ik)
    
  call pop_sub()
  call profiling_out(C_PROFILING_VNLPSI)
end subroutine X(vnlpsi)


! ---------------------------------------------------------
subroutine X(vnlpsi_batch) (h, mesh, psib, hpsib, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(mesh_t),        intent(in)    :: mesh
  type(batch_t),       intent(in)    :: psib
  type(batch_t),       intent(inout) :: hpsib
  integer,             intent(in)    :: ik

  integer :: ii
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)

  call profiling_in(C_PROFILING_VNLPSI)
  call push_sub('hamiltonian_inc.Xvnlpsi')

  do ii = 1, psib%nst
    psi  => psib%states(ii)%X(psi)
    hpsi => hpsib%states(ii)%X(psi)

    call X(project_psi)(mesh, h%ep%proj, h%ep%natoms, h%d%dim, psi, hpsi, ik)
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_VNLPSI)
end subroutine X(vnlpsi_batch)


! ---------------------------------------------------------
subroutine X(vlpsi_batch) (h, m, psib, hpsib, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(mesh_t),        intent(in)    :: m
  integer,             intent(in)    :: ik
  type(batch_t),       intent(in)    :: psib
  type(batch_t),       intent(inout) :: hpsib

  integer :: ip, ip2, ii, ispin, bs, ipmax
  R_TYPE, pointer :: psi(:, :), hpsi(:, :)
  FLOAT, allocatable  :: vv(:)

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('hamiltonian_inc.Xvlpsi')

  select case(h%d%ispin)
  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_dim_get_spin_index(h%d, ik)

    bs = hardware%dblock_size

    ALLOCATE(vv(1:bs), bs) 

    do ip = 1, m%np, bs
      ipmax = min(m%np, ip + bs - 1)

      forall (ip2 = ip:ipmax) vv(ip2 - ip + 1) = h%vhxc(ip2, ispin) + h%ep%vpsl(ip2)

      forall (ii = 1:psib%nst, ip2 = ip:ipmax)
        hpsib%states(ii)%X(psi)(ip2, 1) = vv(ip2 - ip + 1)*psib%states(ii)%X(psi)(ip2, 1)
      end forall
      
    end do

    deallocate(vv)

    call profiling_count_operations((R_ADD + R_MUL*psib%nst)*m%np)
    call profiling_count_transfers(m%np, M_ONE)
    call profiling_count_transfers(m%np*psib%nst, R_TOTYPE(M_ONE))

  case(SPINORS)
    do ii = 1, psib%nst
      psi  => psib%states(ii)%X(psi)
      hpsi => hpsib%states(ii)%X(psi)

      do ip = 1, m%np
        hpsi(ip, 1) = (h%vhxc(ip, 1) + h%ep%vpsl(ip))*psi(ip, 1) + (h%vhxc(ip, 3) + M_zI*h%vhxc(ip, 4))*psi(ip, 2)
        hpsi(ip, 2) = (h%vhxc(ip, 2) + h%ep%vpsl(ip))*psi(ip, 2) + (h%vhxc(ip, 3) - M_zI*h%vhxc(ip, 4))*psi(ip, 1)
      end do

    end do
    call profiling_count_operations((6*R_ADD + 2*R_MUL)*m%np*psib%nst)

  end select

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vlpsi_batch)

! ---------------------------------------------------------
subroutine X(vexternal) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  integer :: idim

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('hamiltonian_inc.Xvlpsi')

  do idim = 1, h%d%dim
    hpsi(1:NP, idim) = hpsi(1:NP, idim) + h%ep%vpsl(1:NP)*psi(1:NP, idim)
  end do

  if(h%ep%non_local) call X(vnlpsi)(h, gr%m, psi, hpsi, ik)

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vexternal)

! ---------------------------------------------------------
subroutine X(vgauge) (gr, h, psi, hpsi, grad)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: ip, idim, a2
  FLOAT :: vecpot(1:MAX_DIM)

  call push_sub('hamiltonian_inc.Xvgauge')

  ASSERT(gauge_field_is_applied(h%ep%gfield))

  call X(get_grad)(h, gr, psi, grad)

  vecpot = gauge_field_get_vec_pot(h%ep%gfield)/P_c
  a2 = sum(vecpot(1:MAX_DIM)**2)
  
  do idim = 1, h%d%dim
    do ip = 1, NP
      hpsi(ip, idim) = hpsi(ip, idim) + &
           M_HALF*a2*psi(ip, idim) + M_zI*dot_product(vecpot(1:MAX_DIM), grad(ip, 1:MAX_DIM, idim))
    end do
  end do

  call pop_sub()
end subroutine X(vgauge)


! ---------------------------------------------------------
subroutine X(vborders) (gr, h, psi, hpsi)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: idim

  call push_sub('hamiltonian_inc.Xvborders')

  if(h%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, h%d%dim
      hpsi(1:NP, idim) = hpsi(1:NP, idim) + M_zI*h%ab_pot(1:NP)*psi(1:NP, idim)
    end do
  end if

  call pop_sub()
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(h_mgga_terms) (h, gr, psi, hpsi, ik, grad)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE,              pointer       :: grad(:, :, :)

  integer :: ispace, idim, ispin
  R_TYPE, allocatable :: cgrad(:,:), diverg(:)

  call push_sub('hamiltonian_inc.Xh_mgga_terms')

  call X(get_grad)(h, gr, psi, grad)
  ispin = states_dim_get_spin_index(h%d, ik)

  ALLOCATE(cgrad(1:NP_PART, 1:MAX_DIM), NP_PART*MAX_DIM)
  ALLOCATE(diverg(1:NP), NP)

  do idim = 1, h%d%dim
    do ispace = 1, gr%sb%dim
      cgrad(1:NP, ispace) = grad(1:NP, ispace, idim)*h%vtau(1:NP, ispin)
      call X(set_bc)(gr%der, cgrad(:, ispace))
    end do

    call X(derivatives_div)(gr%der, cgrad, diverg)

    !hpsi(1:NP, idim) = hpsi(1:NP, idim) -  diverg(1:NP)
  end do

  deallocate(cgrad, diverg)

  call pop_sub()
end subroutine X(h_mgga_terms)


! ---------------------------------------------------------
subroutine X(vmask) (gr, h, st)
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_t), intent(in)    :: h
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim

  call push_sub('hamiltonian_inc.Xvmask')

  if(h%ab == MASK_ABSORBING) then
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
           st%X(psi)(1:NP, idim, ist, ik) = st%X(psi)(1:NP, idim, ist, ik)*(M_ONE - h%ab_pot(1:NP))
        end do
      end do
    end do
  end if

  call pop_sub()
end subroutine X(vmask)

! ---------------------------------------------------------
subroutine X(hpsi_diag) (h, gr, diag, ik, t)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(out)   :: diag(:,:) ! hpsi(NP, h%d%dim)
  FLOAT, optional,     intent(in)    :: t

  integer :: idim, ip, ispin

  FLOAT, allocatable  :: ldiag(:)

  call push_sub('hamiltonian_inc.Xhpsi_diag')
  
  ALLOCATE(ldiag(NP), NP)

  diag = M_ZERO

  call derivatives_lapl_diag(gr%der, ldiag)

  do idim = 1, h%d%dim
    diag(1:NP, idim) = -M_HALF/h%mass*ldiag(1:NP)
  end do

  select case(h%d%ispin)

  case(UNPOLARIZED, SPIN_POLARIZED)
    ispin = states_dim_get_spin_index(h%d, ik)
    diag(1:NP, 1) = diag(1:NP, 1) + h%vhxc(1:NP, ispin) + h%ep%vpsl(1:NP)
    
  case(SPINORS)
    do ip = 1, NP
      diag(ip, 1) = diag(ip, 1) + h%vhxc(ip, 1) + h%ep%vpsl(ip)
      diag(ip, 2) = diag(ip, 2) + h%vhxc(ip, 2) + h%ep%vpsl(ip)
    end do
    
  end select
    
  call pop_sub()
end subroutine X(hpsi_diag)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
