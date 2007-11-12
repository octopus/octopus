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
! calculates the eigenvalues of the real orbitals
subroutine X(hamiltonian_eigenval)(h, gr, st, t)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t) ,       intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  FLOAT, intent(in), optional :: t

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('h_inc.Xhamiltonian_eigenval')
  ALLOCATE(Hpsi(NP, st%d%dim), NP*st%d%dim)

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      if(present(t)) then
        call X(hpsi) (h, gr, st%X(psi)(:, :, ist, ik), hpsi, ik, t)
      else
        call X(hpsi) (h, gr, st%X(psi)(:, :, ist, ik), hpsi, ik)
      end if
      e = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
end subroutine X(hamiltonian_eigenval)


! ---------------------------------------------------------
subroutine X(hpsi) (h, gr, psi, hpsi, ik, t, E)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE, target,      intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  FLOAT, optional,     intent(in)    :: t
  FLOAT, optional,     intent(in)    :: E

  integer :: idim, ip

  R_TYPE, pointer :: epsi(:,:)
  FLOAT :: kpoint(MAX_DIM)
  type(profile_t), save :: phase_prof

  call profiling_in(C_PROFILING_HPSI)
  call push_sub('h_inc.Xhpsi')

  ASSERT(ubound(psi, DIM=1) == NP_PART)

  if(present(t).and.h%d%cdft) then
    message(1) = "TDCDFT not yet implemented"
    call write_fatal(1)
  end if

  ! first of all, set boundary conditions
  do idim = 1, h%d%dim
    call X(set_bc)(gr%f_der%der_discr, psi(:, idim))
  end do

  if(simul_box_is_periodic(gr%sb)) then ! we multiply psi by exp(i k.r)
    call profiling_in(phase_prof, "PBC_PHASE_APPLY")
    kpoint = M_ZERO
    kpoint(1:gr%sb%periodic_dim) = h%d%kpoints(1:gr%sb%periodic_dim, ik)

    ALLOCATE(epsi(1:NP_PART, 1:h%d%dim), NP_PART*h%d%dim)

    do idim = 1, h%d%dim
      do ip = 1, NP_PART
        epsi(ip, idim) = exp(-M_zI * sum(gr%m%x(ip, 1:MAX_DIM) * kpoint(1:MAX_DIM))) * psi(ip, idim)
      end do
    end do
    call profiling_out(phase_prof)
  else
    ! for finite systems we do nothing
    epsi => psi
  end if

  do idim = 1, h%d%dim
    !$omp parallel workshare
    hpsi(:, idim) = M_ZERO
    !$omp end parallel workshare
  end do
  
#if defined(HAVE_LIBNBC)
  if(gr%m%parallel_in_domains) then
    call X(kinetic_prepare)(h, gr, epsi)
    
    call X(vlpsi)(h, gr%m, epsi, hpsi, ik)
    call X(kinetic_keep_going)(h)
    if(h%ep%nvnl > 0) then
      call X(vnlpsi)(h, gr, epsi, hpsi, ik)
      call X(kinetic_keep_going)(h)
    end if
    if(present(t)) call X(vlasers)(gr, h, epsi, hpsi, ik, t)
    
    call X(kinetic_wait)(h)
    call X(kinetic_calculate)(h, gr, epsi, hpsi, ik)
  else
#endif
    call X(kinetic)(h, gr, epsi, hpsi, ik)
    call X(vlpsi)(h, gr%m, epsi, hpsi, ik)
    if(h%ep%nvnl > 0) call X(vnlpsi)(h, gr, epsi, hpsi, ik)
    if(present(t)) call X(vlasers)(gr, h, epsi, hpsi, ik, t)
#if defined(HAVE_LIBNBC)
  end if
#endif
  
  if(h%theory_level .eq. HARTREE_FOCK) then
    call X(exchange_operator)(h, gr, epsi, hpsi, ik)
  end if

  if(hamiltonian_oct_exchange(h)) then
    call X(oct_exchange_operator)(h, gr, epsi, hpsi, ik)
  end if
  
  call X(magnetic_terms) (gr, h, epsi, hpsi, ik)
  if (h%ep%with_gauge_field) call X(vgauge) (gr, h, epsi, hpsi)
  if(present(t)) call X(vborders) (gr, h, epsi, hpsi)
  
  if(present(E)) then
    ! compute (H-E) epsi = hpsi
    do idim = 1, h%d%dim
      call lalg_axpy(NP, E, epsi(:, idim), hpsi(:, idim))
    end do
  end if
  
  if(simul_box_is_periodic(gr%sb)) then
    ! now we need to remove the exp(-i k.r) factor
    call profiling_in(phase_prof)
    do idim = 1, h%d%dim
      do ip = 1, NP
        hpsi(ip, idim) = exp(M_zI * sum(gr%m%x(ip, 1:MAX_DIM) * kpoint(1:MAX_DIM))) * hpsi(ip, idim)
      end do
    end do
    
    deallocate(epsi)
    call profiling_out(phase_prof)
  end if
  
  call pop_sub()
  call profiling_out(C_PROFILING_HPSI)
end subroutine X(hpsi)


! ---------------------------------------------------------
subroutine X(exchange_operator) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('h_inc.Xexchange_operator')

  ALLOCATE(rho(gr%m%np), gr%m%np)
  ALLOCATE(pot(gr%m%np), gr%m%np)

  select case(h%d%ispin)
  case(UNPOLARIZED)
    do j = 1, h%st%nst
      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) - (h%st%occ(j, ik)/M_TWO) * h%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPIN_POLARIZED)
    do j = 1, h%st%nst
      if(h%st%occ(j, ik) <= M_ZERO) cycle
      pot = M_ZERO
      do k = 1, gr%m%np
        rho(k) = R_CONJ(h%st%X(psi)(k, 1, j, ik)) * psi(k, 1)
      end do
      call X(poisson_solve)(gr, pot, rho)
      do k = 1, gr%m%np
        hpsi(k, 1) = hpsi(k, 1) - h%st%occ(j, ik) * h%st%X(psi)(k, 1, j, ik)*pot(k)
      end do
    end do 

  case(SPINORS)

  end select

  deallocate(rho, pot)
  call pop_sub()
end subroutine X(exchange_operator)


! ---------------------------------------------------------
subroutine X(oct_exchange_operator) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  integer,             intent(in)    :: ik
  R_TYPE, allocatable :: rho(:), pot(:)
  integer :: j, k

  call push_sub('h_inc.Xexchange_operator')

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
        hpsi(k, 1) = hpsi(k, 1) + M_TWO * M_z1 * (h%oct_st%occ(j, ik)/M_TWO) * &
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
        hpsi(k, 1) = hpsi(k, 1) + M_TWO * M_z1 * h%oct_st%occ(j, ik) * h%oct_st%X(psi)(k, 1, j, ik)*R_AIMAG(pot(k))
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
  integer :: idim

  ! We will assume, for the moment, no spinors.

  call push_sub('h_inc.Xmagnus')

  ALLOCATE( auxpsi(NP_PART, h%d%dim), NP_PART*h%d%dim)
  ALLOCATE(aux2psi(NP,      h%d%dim), NP*h%d%dim)

  call X(kinetic) (h, gr, psi, hpsi, ik)

  auxpsi = hpsi
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, gr, psi, auxpsi, ik)
  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(1:NP, 1) = hpsi(1:NP, 1) -  M_zI*vmagnus(1:NP, 1, 1)*auxpsi(1:NP, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      hpsi(1:NP, 1) = hpsi(1:NP, 1) - M_zI*vmagnus(1:NP, 1, 1)*auxpsi(1:NP, 1)
    else
      hpsi(1:NP, 1) = hpsi(1:NP, 1) - M_zI*vmagnus(1:NP, 2, 1)*auxpsi(1:NP, 1)
    end if
  end select

  select case(h%d%ispin)
  case(UNPOLARIZED)
    auxpsi(1:NP, 1) = vmagnus(1:NP, 1, 1)*psi(1:NP, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      auxpsi(1:NP, 1) = vmagnus(1:NP, 1, 1)*psi(1:NP, 1)
    else
      auxpsi(1:NP, 1) = vmagnus(1:NP, 2, 1)*psi(1:NP, 1)
    end if
  end select
  call X(kinetic) (h, gr, auxpsi, aux2psi, ik)
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, gr, auxpsi, aux2psi, ik)
  hpsi(1:NP, 1) = hpsi(1:NP, 1) + M_zI*aux2psi(1:NP, 1)

  do idim = 1, h%d%dim
    hpsi(1:NP, idim) = hpsi(1:NP, idim) + h%ep%Vpsl(1:NP)*psi(1:NP,idim)
  end do

  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(1:NP, 1) = hpsi(1:NP, 1) + vmagnus(1:NP, 1, 2)*psi(1:NP, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      hpsi(1:NP, 1) = hpsi(1:NP, 1) + vmagnus(1:NP, 1, 2)*psi(1:NP, 1)
    else
      hpsi(1:NP, 1) = hpsi(1:NP, 1) + vmagnus(1:NP, 2, 2)*psi(1:NP, 1)
    end if
  end select
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, gr, psi, Hpsi, ik)
  call X(vborders) (gr, h, psi, hpsi)

  deallocate(auxpsi, aux2psi)
  call pop_sub()
end subroutine X(magnus)


! ---------------------------------------------------------
! Exchange the ghost points and write the boundary points.
subroutine X(kinetic_prepare)(h, gr, psi)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)

  integer :: idim

  call push_sub('h_inc.Xkinetic_prepare')

  do idim = 1, h%d%dim
    if(gr%m%parallel_in_domains) then
#if defined(HAVE_MPI)
#if defined(HAVE_LIBNBC)
      call X(vec_ighost_update)(gr%m%vp, psi(:, idim), h%handles(idim))
#else
      call X(vec_ghost_update)(gr%m%vp, psi(:, idim))
#endif
#endif
    end if
  end do

  call pop_sub()
end subroutine X(kinetic_prepare)

! ---------------------------------------------------------
! Call NBC_Test to give MPI a chance to process pending messages.
! This improves the obverlap but is only necessary due to deficiencies
! in MPI implementations.
#if defined(HAVE_LIBNBC)
subroutine X(kinetic_keep_going)(h)
  type(hamiltonian_t), intent(in) :: h

  integer :: idim

  call push_sub('h_inc.Xkinetic_wait')

  do idim = 1, h%d%dim
    call NBCF_Test(h%handles(idim), mpi_err)
  end do

  call pop_sub()
end subroutine X(kinetic_keep_going)
#endif


! ---------------------------------------------------------
! Wait for ghost point exchange to finish.
#if defined(HAVE_LIBNBC)
subroutine X(kinetic_wait)(h)
  type(hamiltonian_t), intent(in) :: h

  integer :: idim

  call push_sub('h_inc.Xkinetic_wait')

  do idim = 1, h%d%dim
    call NBCF_Wait(h%handles(idim), mpi_err)
  end do

  call pop_sub()
end subroutine X(kinetic_wait)
#endif


! ---------------------------------------------------------
subroutine X(kinetic) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  call push_sub('h_inc.Xkinetic')

  call X(kinetic_prepare)(h, gr, psi)
#if defined(HAVE_LIBNBC)
  if(gr%m%parallel_in_domains) then
    call X(kinetic_wait)(h)
  end if
#endif
  call X(kinetic_calculate) (h, gr, psi, hpsi, ik)

  call pop_sub()
end subroutine X(kinetic)


! ---------------------------------------------------------
subroutine X(kinetic_calculate) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  integer             :: idim
  R_TYPE, allocatable :: lapl(:, :)

  call profiling_in(C_PROFILING_KINETIC)
  call push_sub('h_inc.Xkinetic_calculate')

  ALLOCATE(lapl(gr%m%np, h%d%dim), gr%m%np*h%d%dim)

  do idim = 1, h%d%dim
    call X(f_laplacian) (gr%sb, gr%f_der, psi(:, idim), lapl(:, idim), &
         cutoff_ = M_TWO*h%cutoff, have_bndry=.true., ghost_update=.false.)
    call lalg_axpy(NP, -M_HALF/h%mass, lapl(:, idim), hpsi(:, idim))
  end do
  
  deallocate(lapl)

  call pop_sub()
  call profiling_out(C_PROFILING_KINETIC)
end subroutine X(kinetic_calculate)

! ---------------------------------------------------------
! Here we the terms arising from the presence of a possible static external
! magnetic field, and the terms that come from CDFT.
subroutine X(magnetic_terms) (gr, h, psi, hpsi, ik)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(inout) :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP, h%d%dim)
  integer,             intent(in)    :: ik

  integer :: k, idim
  FLOAT,  allocatable :: div(:), tmp(:,:)
  R_TYPE, allocatable :: grad(:, :, :), lhpsi(:, :)

  call push_sub('h_inc.Xmagnetic_terms')

  if(h%d%cdft .or. associated(h%ep%A_static)) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)
    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do
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
    call df_divergence(gr%sb, gr%f_der, tmp, div)
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

  endif

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

  deallocate(grad)
  call pop_sub()
end subroutine X(magnetic_terms)


! ---------------------------------------------------------
subroutine X(vnlpsi) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  integer :: ipj

  call profiling_in(C_PROFILING_VNLPSI)
  call push_sub('h_inc.Xvnlpsi')

  do ipj = 1, h%ep%nvnl
    call X(project_psi)(gr%m, h%ep%p(ipj), h%d%dim, psi, hpsi, &
       h%reltype, periodic = simul_box_is_periodic(gr%m%sb), ik = ik)
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_VNLPSI)
end subroutine X(vnlpsi)

! ---------------------------------------------------------
subroutine X(vlpsi) (h, m, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(mesh_t),        intent(in)    :: m
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: idim

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('h_inc.Xvlpsi')

  select case(h%d%ispin)
  case(UNPOLARIZED)
    !$omp parallel workshare
    hpsi(1:m%np, 1) = hpsi(1:m%np, 1) + (h%vhxc(1:m%np, 1) + h%ep%vpsl(1:m%np))*psi(1:m%np, 1)
    !$omp end parallel workshare
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then ! we have a spin down
      hpsi(1:m%np, 1) = hpsi(1:m%np, 1) + (h%vhxc(1:m%np, 1) + h%ep%vpsl(1:m%np))*psi(1:m%np, 1)
    else
      hpsi(1:m%np, 1) = hpsi(1:m%np, 1) + (h%vhxc(1:m%np, 2) + h%ep%vpsl(1:m%np))*psi(1:m%np, 1)
    end if
  case(SPINORS)
    hpsi(1:m%np, 1) = hpsi(1:m%np, 1) + (h%vhxc(1:m%np, 1) + &
      h%ep%vpsl(1:m%np))*psi(1:m%np, 1) +                    &
      (h%vhxc(1:m%np, 3) + M_zI*h%vhxc(1:m%np, 4))*psi(1:m%np, 2)
    hpsi(1:m%np, 2) = hpsi(1:m%np, 2) + (h%vhxc(1:m%np, 2) + &
      h%ep%vpsl(1:m%np))*psi(1:m%np, 2) +                    &
      (h%vhxc(1:m%np, 3) - M_zI*h%vhxc(1:m%np, 4))*psi(1:m%np, 1)
  end select

  if (associated(h%ep%E_field)) then
    do idim = 1, h%d%dim
      hpsi(1:m%np, idim) = hpsi(1:m%np, idim) + h%ep%v_static*psi(1:m%np, idim)
    end do
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vlpsi)


! ---------------------------------------------------------
subroutine X(vexternal) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: idim

  call profiling_in(C_PROFILING_VLPSI)
  call push_sub('h_inc.Xvlpsi')

  select case(h%d%ispin)
  case(UNPOLARIZED, SPIN_POLARIZED)
    hpsi(1:NP, 1) = hpsi(1:NP, 1) + (h%ep%vpsl(1:NP))*psi(1:NP, 1)
  case(SPINORS)
    hpsi(1:NP, 1) = hpsi(1:NP, 1) + (h%ep%vpsl(1:NP))*psi(1:NP, 1)
    hpsi(1:NP, 2) = hpsi(1:NP, 2) + (h%ep%vpsl(1:NP))*psi(1:NP, 2)
  end select

  if (associated(h%ep%E_field)) then
    do idim = 1, h%d%dim
      hpsi(1:NP, idim) = hpsi(1:NP, idim) + h%ep%v_static*psi(1:NP, idim)
    end do
  end if

  if(h%ep%nvnl > 0) call X(vnlpsi)  (h, gr, psi, hpsi, ik)

  call pop_sub()
  call profiling_out(C_PROFILING_VLPSI)
end subroutine X(vexternal)

! ---------------------------------------------------------
subroutine X(vlaser_operator_quadratic) (gr, h, psi, hpsi, ik, laser_number)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik
  integer,             intent(in)    :: laser_number

  integer :: i, k, idim
  logical :: electric_field, vector_potential, magnetic_field

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)
  FLOAT, allocatable :: a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlasers')

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.
  i = laser_number

  select case(h%ep%lasers(i)%field)
  case(E_FIELD_ELECTRIC) ! do nothing
  case(E_FIELD_MAGNETIC)
    if(.not. allocated(a)) then 
      ALLOCATE(a(NP_PART, NDIM), NP_PART*NDIM)
      a = M_ZERO
      ALLOCATE(a_prime(NP_PART, NDIM), NP_PART*NDIM)
      a_prime = M_ZERO
    end if
    call laser_vector_potential(h%ep%lasers(i), gr%m, a_prime)
    a = a + a_prime
    call laser_field(gr%sb, h%ep%lasers(i), b_prime)
    b = b + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    call laser_field(gr%sb, h%ep%lasers(i), a_field_prime)
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(magnetic_field) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a(k, 1:NDIM), a(k, 1:NDIM))*psi(k, :) / P_c**2
    end do
    deallocate(a, a_prime)
  end if
  if(vector_potential) then
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a_field(1:NDIM), a_field(1:NDIM))*psi(k, :) / P_c**2
    end do
  end if

  call pop_sub()
end subroutine X(vlaser_operator_quadratic)


! ---------------------------------------------------------
subroutine X(vlaser_operator_linear) (gr, h, psi, hpsi, ik, laser_number)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik
  integer,             intent(in)    :: laser_number

  integer :: i, k, idim
  logical :: electric_field, vector_potential, magnetic_field
  R_TYPE, allocatable :: grad(:, :, :), lhpsi(:, :)

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)

  FLOAT, allocatable :: v(:), pot(:), a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlasers')

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.
  i = laser_number

  select case(h%ep%lasers(i)%field)
  case(E_FIELD_SCALAR_POTENTIAL)
    if(.not. allocated(v)) then 
      ALLOCATE(v(NP), NP)
      v = M_ZERO
    end if
    call laser_potential(gr%sb, h%ep%lasers(i), gr%m, v)
    electric_field = .true.
  case(E_FIELD_ELECTRIC)
    if(.not. allocated(v)) then 
      ALLOCATE(v(NP), NP)
      v = M_ZERO
      ALLOCATE(pot(NP), NP)
      pot = M_ZERO
    end if
    call laser_potential(gr%sb, h%ep%lasers(i), gr%m, pot)
    v = v + pot
    electric_field = .true.
    deallocate(pot)
  case(E_FIELD_MAGNETIC)
    if(.not. allocated(a)) then 
      ALLOCATE(a(NP_PART, NDIM), NP_PART*NDIM)
      a = M_ZERO
      ALLOCATE(a_prime(NP_PART, NDIM), NP_PART*NDIM)
      a_prime = M_ZERO
    end if
    call laser_vector_potential(h%ep%lasers(i), gr%m, a_prime)
    a = a + a_prime
    call laser_field(gr%sb, h%ep%lasers(i), b_prime)
    b = b + b_prime
    magnetic_field = .true.
  case(E_FIELD_VECTOR_POTENTIAL)
    call laser_field(gr%sb, h%ep%lasers(i), a_field_prime)
    a_field = a_field + a_field_prime
    vector_potential = .true.
  end select

  if(electric_field) then
    do k = 1, h%d%dim
      hpsi(1:NP, k)= hpsi(1:NP, k) + v(1:NP) * psi(1:NP, k)
    end do
    deallocate(v)
  end if

  if(magnetic_field) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)

    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do

    ! If there is a static magnetic field, its associated vector potential is coupled with
    ! the time-dependent one defined as a "laser" (ideally one should just add them all and
    ! do the calculation only once...). Note that h%ep%a_static already has been divided
    ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
    if(associated(h%ep%a_static)) then
      do k = 1, NP
        hpsi(k, :) = hpsi(k, :) + dot_product(a(k, 1:NDIM), h%ep%a_static(k, 1:NDIM))*psi(k, :) / P_c
      end do
    end if

    select case(h%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a(k, 1:NDIM), grad(k, 1:NDIM, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a(k, 1:NDIM), grad(k, 1:NDIM, idim)) / P_c
        end do
      end do
    end select


    select case (h%d%ispin)
    case (SPIN_POLARIZED)
      ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      end if
      hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    case (SPINORS)
      ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
      lhpsi(1:NP, 1) = M_HALF/P_C*( b(3)*psi(1:NP, 1) &
           + (b(1) - M_zI*b(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-b(3)*psi(1:NP, 2) &
           + (b(1) + M_zI*b(2))*psi(1:NP, 1))
      hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    end select

    deallocate(grad)
    deallocate(a, a_prime)
  end if

  if(vector_potential) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)

    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do

    select case(h%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a_field(1:NDIM), grad(k, 1:NDIM, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a_field(1:NDIM), grad(k, 1:NDIM, idim)) / P_c
        end do
      end do
    end select
    deallocate(grad)
  end if

  call pop_sub()
end subroutine X(vlaser_operator_linear)


! ---------------------------------------------------------
subroutine X(vlasers) (gr, h, psi, hpsi, ik, t)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik
  FLOAT, optional                    :: t

  integer :: i, k, idim
  logical :: electric_field, vector_potential, magnetic_field
  R_TYPE, allocatable :: grad(:, :, :), lhpsi(:, :)

  FLOAT :: a_field(1:MAX_DIM), a_field_prime(1:MAX_DIM), b(1:MAX_DIM), b_prime(1:MAX_DIM)

  FLOAT, allocatable :: v(:), pot(:), a(:, :), a_prime(:, :)

  call push_sub('h_inc.Xvlasers')

  a_field = M_ZERO

  electric_field = .false.
  vector_potential = .false.
  magnetic_field = .false.

  do i = 1, h%ep%no_lasers

    select case(h%ep%lasers(i)%field)
    case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
      if(.not. allocated(v)) then 
        ALLOCATE(v(NP), NP)
        v = M_ZERO
      end if
      ALLOCATE(pot(NP), NP)
      pot = M_ZERO
      call laser_potential(gr%sb, h%ep%lasers(i), gr%m, pot, t)
      v = v + pot
      electric_field = .true.
      deallocate(pot)

    case(E_FIELD_MAGNETIC)
      if(.not. allocated(a)) then 
        ALLOCATE(a(NP_PART, NDIM), NP_PART*NDIM)
        a = M_ZERO
        ALLOCATE(a_prime(NP_PART, NDIM), NP_PART*NDIM)
        a_prime = M_ZERO
      end if

      call laser_vector_potential(h%ep%lasers(i), gr%m, a_prime, t)
      a = a + a_prime
      if(present(t)) then
        call laser_field(gr%sb, h%ep%lasers(i), b_prime, t)
      else 
        call laser_field(gr%sb, h%ep%lasers(i), b_prime)
      end if
      b = b + b_prime
      magnetic_field = .true.

    case(E_FIELD_VECTOR_POTENTIAL)
      call laser_field(gr%sb, h%ep%lasers(i), a_field_prime, t)
      a_field = a_field + a_field_prime
      vector_potential = .true.

    end select

  end do

  if(electric_field) then
    do k = 1, h%d%dim
      hpsi(1:NP, k)= hpsi(1:NP, k) + v(1:NP) * psi(1:NP, k)
    end do
    deallocate(v)
  end if

  if(magnetic_field) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)

    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do

    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a(k, 1:NDIM), a(k, 1:NDIM))*psi(k, :) / P_c**2
    end do

    ! If there is a static magnetic field, its associated vector potential is coupled with
    ! the time-dependent one defined as a "laser" (ideally one should just add them all and
    ! do the calculation only once...). Note that h%ep%a_static already has been divided
    ! by P_c, and therefore here we only divide by P_c, and not P_c**2.
    if(associated(h%ep%a_static)) then
      do k = 1, NP
        hpsi(k, :) = hpsi(k, :) + dot_product(a(k, 1:NDIM), h%ep%a_static(k, 1:NDIM))*psi(k, :) / P_c
      end do
    end if

    select case(h%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a(k, 1:NDIM), grad(k, 1:NDIM, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a(k, 1:NDIM), grad(k, 1:NDIM, idim)) / P_c
        end do
      end do
    end select


    select case (h%d%ispin)
    case (SPIN_POLARIZED)
      ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      end if
      hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    case (SPINORS)
      ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
      lhpsi(1:NP, 1) = M_HALF/P_C*( b(3)*psi(1:NP, 1) &
           + (b(1) - M_zI*b(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-b(3)*psi(1:NP, 2) &
           + (b(1) + M_zI*b(2))*psi(1:NP, 1))
      hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
      deallocate(lhpsi)
    end select

    deallocate(grad)
    deallocate(a, a_prime)
  end if

  if(vector_potential) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)

    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do

    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a_field(1:NDIM), a_field(1:NDIM))*psi(k, :) / P_c**2
    end do

    select case(h%d%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do k = 1, NP
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a_field(1:NDIM), grad(k, 1:NDIM, 1)) / P_c
      end do
    case (SPINORS)
      do k = 1, NP
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(a_field(1:NDIM), grad(k, 1:NDIM, idim)) / P_c
        end do
      end do
    end select
    deallocate(grad)
  end if

  call pop_sub()
end subroutine X(vlasers)


! ---------------------------------------------------------
subroutine X(vgauge) (gr, h, psi, hpsi)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: k, idim
  R_TYPE, allocatable :: grad(:,:,:)

  call push_sub('h_inc.Xvgauge')

  if (associated(h%ep%A_gauge)) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)
    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do
    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*sum(h%ep%A_gauge(:)**2)*psi(k, :)
      select case(h%d%ispin)
      case(UNPOLARIZED, SPIN_POLARIZED)
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%ep%A_gauge(1:NDIM), grad(k, 1:NDIM, 1))
      case (SPINORS)
        do idim = 1, h%d%dim
          hpsi(k, idim) = hpsi(k, idim) - M_zI*dot_product(h%ep%A_gauge(1:NDIM), grad(k, 1:NDIM, idim))
        end do
      end select
    end do
  end if

  call pop_sub()
end subroutine X(vgauge)


! ---------------------------------------------------------
subroutine X(vborders) (gr, h, psi, hpsi)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)

  integer :: idim

  call push_sub('h_inc.Xvborders')

  if(h%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, h%d%dim
      hpsi(1:NP, idim) = hpsi(1:NP, idim) + M_zI*h%ab_pot(1:NP)*psi(1:NP, idim)
    end do
  end if

  call pop_sub()
end subroutine X(vborders)


! ---------------------------------------------------------
subroutine X(vmask) (gr, h, st)
  type(grid_t),        intent(in)    :: gr
  type(hamiltonian_t), intent(in)    :: h
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim

  call push_sub('h_inc.Xvmask')

  if(h%ab == MASK_ABSORBING) then
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
           st%X(psi)(1:NP, idim, ist, ik) = st%X(psi)(1:NP, idim, ist, ik) * &
             (M_ONE - h%ab_pot(1:NP))
        end do
      end do
    end do
  end if

  call pop_sub()
end subroutine X(vmask)


! ---------------------------------------------------------
FLOAT function X(electronic_kinetic_energy)(h, gr, st) result(t0)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: tpsi(:, :)
  FLOAT, allocatable :: t(:, :)

  call push_sub('h.electronic_kinetic_energy')

  ALLOCATE(tpsi(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(t(st%st_start:st%st_end, st%d%nik), st%nst*st%d%nik)
  t = M_ZERO

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      tpsi = R_TOTYPE(M_ZERO)
      call X(kinetic) (h, gr, st%X(psi)(:, :, ist, ik), tpsi, ik)
      t(ist, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), tpsi)
    end do
  end do

  t0 = states_eigenvalues_sum(st, t)

  deallocate(tpsi, t)
  call pop_sub()
end function X(electronic_kinetic_energy)


! ---------------------------------------------------------
FLOAT function X(electronic_external_energy)(h, gr, st) result(v)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: vpsi(:, :)
  FLOAT, allocatable :: t(:, :)

  call push_sub('h.electronic_external_energy')

  ALLOCATE(vpsi(NP_PART, st%d%dim), NP_PART*st%d%dim)
  ALLOCATE(t(st%st_start:st%st_end, st%d%nik), st%nst*st%d%nik)
  t = M_ZERO

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      vpsi = R_TOTYPE(M_ZERO)
      call X(vexternal) (h, gr, st%X(psi)(:, :, ist, ik), vpsi, ik)
      t(ist, ik) = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), vpsi)
    end do
  end do

  v = states_eigenvalues_sum(st, t)

  deallocate(vpsi, t)
  call pop_sub()
end function X(electronic_external_energy)

! ---------------------------------------------------------
subroutine X(hpsi_diag) (h, gr, diag, ik, t, E)
  type(hamiltonian_t), intent(inout) :: h
  type(grid_t),        intent(inout) :: gr
  integer,             intent(in)    :: ik
  R_TYPE,              intent(out)   :: diag(:,:) ! hpsi(NP, h%d%dim)
  FLOAT, optional,     intent(in)    :: t
  FLOAT, optional,     intent(in)    :: E

  integer :: idim

  R_TYPE, allocatable :: psi(:,:)
  FLOAT, allocatable  :: ldiag(:)

  call push_sub('h_inc.Xhpsi_diag')
  
  ALLOCATE(psi(NP, h%d%dim), NP*h%d%dim)
  ALLOCATE(ldiag(NP), NP)

  psi = M_ONE
  diag = M_ZERO

  call f_laplacian_diag(gr%sb, gr%f_der, ldiag)

  do idim = 1, h%d%dim
    diag(1:NP, idim) = -M_HALF/h%mass * ldiag(1:NP)
  end do

  call X(vlpsi)(h, gr%m, psi, diag, ik)
  call X(vnlpsi_diag)(h, gr, psi, diag, ik)

  call pop_sub()
end subroutine X(hpsi_diag)

! ---------------------------------------------------------
subroutine X(vnlpsi_diag) (h, gr, psi, hpsi, ik)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(in)    :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik

  integer :: ipj

  call push_sub('h_inc.Xvnlpsi')

  do ipj = 1, h%ep%nvnl
    if( h%ep%p(ipj)%type == M_LOCAL) then 
      call X(project_psi)(gr%m, h%ep%p(ipj), h%d%dim, psi, hpsi, &
        h%reltype, periodic = simul_box_is_periodic(gr%m%sb), ik = ik)
    end if
  end do

  call pop_sub()
end subroutine X(vnlpsi_diag)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
