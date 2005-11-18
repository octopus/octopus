!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

subroutine X(v_ks_calc)(gr, ks, h, st, calc_eigenval)
  type(grid_type),        intent(inout) :: gr
  type(v_ks_type),        intent(inout) :: ks
  type(hamiltonian_type), intent(inout) :: h
  type(states_type),      intent(inout) :: st
  logical,      optional, intent(in)    :: calc_eigenval

  FLOAT :: amaldi_factor

  call push_sub('v_ks_inc.Xv_ks_calc')

  h%epot = M_ZERO
  h%vhxc = M_ZERO
  if(h%d%cdft) h%ahxc = M_ZERO

  ! check if we should introduce the amaldi SIC correction
  amaldi_factor = M_ONE
  if(ks%sic_type == sic_amaldi) amaldi_factor = (st%qtot-1)/st%qtot

  ! No Hartree or xc if independent electrons
  if((.not.ks%ip_app).and.(amaldi_factor>M_ZERO)) then
    call v_hartree()
!    if(h%d%cdft) call a_hartree()
    call v_a_xc()
  end if

  if(present(calc_eigenval)) call X(hamiltonian_eigenval) (h, gr, st)

  call pop_sub()

contains

  ! ---------------------------------------------------------
  ! Hartree contribution to the xc potential
  subroutine v_hartree()
    FLOAT, allocatable :: rho(:)
    integer :: is

    allocate(rho(NP))

    ! calculate the total density
    rho(:) = st%rho(:, 1)
    do is = 2, h%d%spin_channels
      rho(:) = rho(:) + st%rho(:, is)
    end do

    ! Amaldi correction
    if(ks%sic_type == sic_amaldi) rho(:) = amaldi_factor*rho(:)

    ! solve the poisson equation
    call dpoisson_solve(gr, h%vhartree, rho)

    h%vhxc(:, 1) = h%vhxc(:, 1) + h%vhartree(:)
    if(h%d%ispin > UNPOLARIZED) h%vhxc(:, 2) = h%vhxc(:, 2) + h%vhartree(:)

    ! We first add 1/2 int n vH, to then subtract int n (vxc + vH)
    ! this yields the correct formula epot = - int n (vxc + vH/2)
    do is = 1, h%d%spin_channels
      h%epot = h%epot + M_HALF*dmf_dotp(gr%m, st%rho(:, is), h%vhartree)
    end do

    deallocate(rho)
  end subroutine v_hartree


  ! ---------------------------------------------------------
  ! Hartree contribution from the current to Axc
  subroutine a_hartree()
    FLOAT, allocatable :: j(:,:), ahartree(:,:)
    integer :: is, idim

    allocate(j(NP, NDIM), ahartree(NP, NDIM))
    ahartree = M_ZERO

    ! calculate the total current
    j(:,:) = st%j(:,:, 1)
    do is = 2, h%d%spin_channels
      j(:,:) = j(:,:) + st%j(:, :, is)
    end do

    ! solve poisson equation
    do idim = 1, NDIM
      call dpoisson_solve(gr, ahartree(:,idim), j(:, idim)) ! solve the poisson equation
    end do

    h%ahxc(:,:, 1) = h%ahxc(:,:, 1) + ahartree(:,:)
    if(h%d%ispin > UNPOLARIZED) h%ahxc(:,:, 2) = h%ahxc(:,:, 2) + ahartree(:,:)

    ! We first add 1/2 int j.aH, to then subtract int j.(axc + aH)
    ! this yields the correct formula epot = - int j.(axc + aH/2)
    ! WARNING 1: the axc we store is in fact axc/c. So, to get the energy rigth, we have to multiply by c.
    do is = 1, h%d%spin_channels
      do idim = 1, NDIM
        h%epot = h%epot + M_HALF*P_c*dmf_dotp(gr%m, st%j(idim, :, is), ahartree(:, idim))
      end do
    end do

    deallocate(j, ahartree)
  end subroutine a_hartree


  ! ---------------------------------------------------------
  subroutine v_a_xc()
    FLOAT, allocatable :: rho(:, :)
    CMPLX, allocatable :: ztmp1(:), ztmp2(:)
    integer :: is

    ! next 5 lines are for an RPA calculation
    !
    !FLOAT, save, pointer ::  RPA_Vhxc(:, :)
    !logical, save :: RPA_first = .true.
    !
    !if(RPA_first) then
    !  allocate(RPA_Vhxc(NP, h%d%nspin))
    !  RPA_Vhxc = h%vhxc

    call profiling_in(C_PROFILING_XC)

    h%ex = M_ZERO
    h%ec = M_ZERO
    h%exc_j = M_ZERO

    ! get density taking into account non-linear core corrections, and the Amaldi SIC correction
    allocate(rho(NP, st%d%nspin))
    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
      end do
    else
      rho = st%rho
    end if

    ! Amaldi correction
    if(ks%sic_type == sic_amaldi) rho(:,:) = amaldi_factor*rho(:,:)

    if(h%d%cdft) then
      call xc_get_vxc_and_axc(gr, ks%xc, rho, st%j, st%d%ispin, h%vhxc, h%ahxc, &
         h%ex, h%ec, h%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
    else
      call xc_get_vxc(gr, ks%xc, rho, st%d%ispin, h%vhxc, h%ex, h%ec, &
         -minval(st%eigenval(st%nst, :)), st%qtot)
    end if
    deallocate(rho)

    ! The OEP family has to handle specially
    call X(xc_oep_calc)(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
       gr, h, st, h%vhxc, h%ex, h%ec)

    ! next 5 lines are for an RPA calculation
    !  RPA_Vhxc = h%vhxc - RPA_Vhxc  ! RPA_Vhxc now includes the xc potential
    !  RPA_first = .false.
    !else
    !  h%vhxc = h%vhxc + RPA_Vhxc
    !end if

    select case(h%d%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot - dmf_dotp(gr%m, st%rho(:, 1), h%vhxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot - dmf_dotp(gr%m, st%rho(:, 1), h%vhxc(:, 1)) &
         - dmf_dotp(gr%m, st%rho(:, 2), h%vhxc(:, 2))
    case(SPINORS)
      h%epot = h%epot - dmf_dotp(gr%m, st%rho(:, 1), h%vhxc(:, 1)) &
         - dmf_dotp(gr%m, st%rho(:, 2), h%vhxc(:, 2))

      allocate(ztmp1(NP), ztmp2(NP))
      ztmp1 = st%rho(:, 3) + M_zI*st%rho(:, 4)
      ztmp2 = h%vhxc(:, 3) - M_zI*h%vhxc(:, 4)
      ! WARNING missing real() ???
      h%epot = h%epot - M_TWO*zmf_dotp(gr%m, ztmp1, ztmp2)
      deallocate(ztmp1, ztmp2)
    end select

    call profiling_out(C_PROFILING_XC)
  end subroutine v_a_xc
end subroutine X(v_ks_calc)
