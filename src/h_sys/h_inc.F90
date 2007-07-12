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
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(out)   :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  FLOAT, optional,     intent(in)    :: t
  FLOAT, optional,     intent(in)    :: E

  call profiling_in(C_PROFILING_HPSI)
  call push_sub('h_inc.Xhpsi')

  ASSERT(ubound(psi, DIM=1) == NP_PART)

  if(present(t).and.h%d%cdft) then
    message(1) = "TDCDFT not yet implemented"
    call write_fatal(1)
  end if

  hpsi = R_TOTYPE(M_ZERO)

#if defined(HAVE_LIBNBC)
  if(gr%m%parallel_in_domains) then
    call X(kinetic_prepare)(h, gr, psi)

    call X(vlpsi)(h, gr%m, psi, hpsi, ik)
    if(h%ep%nvnl > 0) call X(vnlpsi)(h, gr, psi, hpsi, ik)
    if(present(t)) call X(vlasers)(gr, h, psi, hpsi, ik, t)

    call X(kinetic_wait)(h)
    call X(kinetic_calculate)(h, gr, psi, hpsi, ik)
  else
#endif
    call X(kinetic)(h, gr, psi, hpsi, ik)
    call X(vlpsi)(h, gr%m, psi, hpsi, ik)
    if(h%ep%nvnl > 0) call X(vnlpsi)(h, gr, psi, hpsi, ik)
    if(present(t)) call X(vlasers)(gr, h, psi, hpsi, ik, t)
#if defined(HAVE_LIBNBC)
  end if
#endif

  call X(magnetic_terms) (gr, h, psi, hpsi, ik)
  if (h%ep%with_gauge_field) call X(vgauge) (gr, h, psi, hpsi)
  if(present(t)) call X(vborders) (gr, h, psi, hpsi)

  if(present(E)) then
    ! compute (H-E) psi = hpsi
    hpsi = hpsi - E*psi
  end if
  
  call pop_sub()
  call profiling_out(C_PROFILING_HPSI)
end subroutine X(Hpsi)


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
#if defined(HAVE_MPI)
subroutine X(kinetic_prepare) (h, gr, psi)
  type(hamiltonian_t), intent(in)    :: h
  type(grid_t),        intent(inout) :: gr
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)

  integer :: idim

  call push_sub('h_inc.Xkinetic_prepare')

  do idim = 1, h%d%dim
#if defined(HAVE_LIBNBC)
    call X(vec_ighost_update)(gr%m%vp, psi(:, idim), h%handles(idim))
#else
    call X(vec_ghost_update)(gr%m%vp, psi(:, idim))
#endif
  end do

  call pop_sub()
end subroutine X(kinetic_prepare)
#endif


! ---------------------------------------------------------
#if defined(HAVE_LIBNBC)
subroutine X(kinetic_wait) (h)
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

#if defined(HAVE_MPI)
  if(gr%m%parallel_in_domains) then
    call X(kinetic_prepare)(h, gr, psi)
#if defined(HAVE_LIBNBC)
    call X(kinetic_wait)(h)
#endif
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

#if defined(R_TCOMPLEX)
  integer :: i
  R_TYPE, allocatable :: grad(:,:)
  FLOAT :: k2
#endif

  call profiling_in(C_PROFILING_KINETIC)
  call push_sub('h_inc.Xkinetic_calculate')

  ALLOCATE(lapl(gr%m%np, h%d%dim), gr%m%np*h%d%dim)

  if(simul_box_is_periodic(gr%sb)) then
#if defined(R_TCOMPLEX)
    ALLOCATE(grad(NP, NDIM), NP*NDIM)
    k2 = sum(h%d%kpoints(:, ik)**2)
    do idim = 1, h%d%dim
      call X(f_laplacian) (gr%sb, gr%f_der, psi(:, idim), lapl(:, idim), &
        cutoff_ = M_TWO*h%cutoff, ghost_update=.false.)
      call X(f_gradient)  (gr%sb, gr%f_der, psi(:, idim), grad(:, :), ghost_update=.false.)
      do i = 1, NP
        hpsi(i, idim) = hpsi(i, idim) - M_HALF*(lapl(i, idim)       &
          + M_TWO*M_zI*sum(h%d%kpoints(1:NDIM, ik)*grad(i, 1:NDIM)) &
          - k2*psi(i, idim))
      end do
    end do
    deallocate(grad)
#else
    idim = ik ! This lines sole purpose is to avoid unsed variable messages
              ! if R_TCOMPLEX is undefined.
    message(1) = "Real wave-function for ground state not yet implemented for polymers:"
    message(2) = "use complex wave-functions instead."
    call write_fatal(2)
#endif
  else
    do idim = 1, h%d%dim
      call X(f_laplacian) (gr%sb, gr%f_der, psi(:, idim), lapl(:, idim), &
        cutoff_ = M_TWO*h%cutoff, ghost_update=.false.)
      call lalg_axpy(NP, -M_HALF, lapl(:, idim), hpsi(:, idim))
    end do
  end if

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
    hpsi(1:m%np, 1) = hpsi(1:m%np, 1) + (h%vhxc(1:m%np, 1) + h%ep%vpsl(1:m%np))*psi(1:m%np, 1)
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
subroutine X(vlasers) (gr, h, psi, hpsi, ik, t)
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  R_TYPE,              intent(inout) :: psi(:,:)  ! psi(NP_PART, h%d%dim)
  R_TYPE,              intent(inout) :: hpsi(:,:) ! hpsi(NP_PART, h%d%dim)
  integer,             intent(in)    :: ik
  FLOAT,               intent(in)    :: t

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
    case(E_FIELD_ELECTRIC)
      if(.not. allocated(v)) then 
        ALLOCATE(v(NP), NP)
        v = M_ZERO
        ALLOCATE(pot(NP), NP)
        pot = M_ZERO
      end if

      call laser_potential(gr%sb, h%ep%lasers(i), t, gr%m, pot)
      v = v + pot
      electric_field = .true.

    case(E_FIELD_MAGNETIC)
      if(.not. allocated(a)) then 
        ALLOCATE(a(NP_PART, NDIM), NP_PART*NDIM)
        a = M_ZERO
        ALLOCATE(a_prime(NP_PART, NDIM), NP_PART*NDIM)
        a_prime = M_ZERO
      end if

      call laser_vector_potential(h%ep%lasers(i), t, gr%m, a_prime)
      a = a + a_prime
      call laser_field(gr%sb, h%ep%lasers(i), t, b_prime)
      b = b + b_prime
      magnetic_field = .true.

    case(E_FIELD_VECTOR_POTENTIAL)
      call laser_field(gr%sb, h%ep%lasers(i), t, a_field_prime)
      a_field = a_field + a_field_prime
      vector_potential = .true.

    end select

  end do

  if(electric_field) then
    do k = 1, h%d%dim
      hpsi(1:NP, k)= hpsi(1:NP, k) + v(1:NP) * psi(1:NP, k)
    end do
    deallocate(v, pot)
  end if

  if(magnetic_field) then
    ALLOCATE(grad(NP_PART, NDIM, h%d%dim), NP_PART*h%d%dim*NDIM)

    do idim = 1, h%d%dim
      call X(f_gradient)(gr%sb, gr%f_der, psi(:, idim), grad(:, :, idim))
    end do

    do k = 1, NP
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a(k, 1:NDIM), a(k, 1:NDIM))*psi(k, :) / P_c**2
    end do

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

    ALLOCATE(lhpsi(NP, h%d%dim), NP*h%d%dim)
    select case (h%d%ispin)
    case (SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(1:NP, 1) = - M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      else
        lhpsi(1:NP, 1) = + M_HALF/P_C*sqrt(dot_product(b, b))*psi(1:NP, 1)
      end if
    case (SPINORS)
      lhpsi(1:NP, 1) = M_HALF/P_C*( b(3)*psi(1:NP, 1) &
           + (b(1) - M_zI*b(2))*psi(1:NP, 2))
      lhpsi(1:NP, 2) = M_HALF/P_C*(-b(3)*psi(1:NP, 2) &
           + (b(1) + M_zI*b(2))*psi(1:NP, 1))
    end select
    hpsi(1:NP, :) = hpsi(1:NP, :) + (h%ep%gyromagnetic_ratio * M_HALF) * lhpsi(1:NP, :)
    deallocate(lhpsi)
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
  ALLOCATE(t(st%nst, st%d%nik), st%nst*st%d%nik)
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
  ALLOCATE(t(st%nst, st%d%nik), st%nst*st%d%nik)
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
