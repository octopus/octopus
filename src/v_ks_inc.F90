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
  type(grid_t),        intent(inout) :: gr
  type(v_ks_t),        intent(inout) :: ks
  type(hamiltonian_t), intent(inout) :: h
  type(states_t),      intent(inout) :: st
  logical,      optional, intent(in) :: calc_eigenval

  FLOAT :: amaldi_factor

  call push_sub('v_ks_inc.Xv_ks_calc')

  h%epot     = M_ZERO
  h%ehartree = M_ZERO
  h%vhxc     = M_ZERO
  if(h%d%cdft) h%axc = M_ZERO

  ! check if we should introduce the amaldi SIC correction
  amaldi_factor = M_ONE
  if(ks%sic_type == sic_amaldi) amaldi_factor = (st%qtot-1)/st%qtot

  ! No Hartree or xc if independent electrons
  if((.not.ks%ip_app).and.(amaldi_factor>M_ZERO)) then
    call v_hartree()
    h%vhxc(1:NP, 1) = h%vhxc(1:NP, 1) + h%vhartree(1:NP)
    if(h%d%ispin > UNPOLARIZED) h%vhxc(1:NP, 2) = h%vhxc(1:NP, 2) + h%vhartree(1:NP)
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

    ALLOCATE(rho(NP), NP)

    ! calculate the total density
    rho(1:NP) = st%rho(1:NP, 1)
    do is = 2, h%d%spin_channels
      rho(1:NP) = rho(1:NP) + st%rho(1:NP, is)
    end do

    ! Amaldi correction
    if(ks%sic_type == sic_amaldi) rho = amaldi_factor*rho

    ! solve the poisson equation
    call dpoisson_solve(gr, h%vhartree, rho)

    ! Get the Hartree energy
    h%ehartree = M_HALF*dmf_dotp(gr%m, rho, h%vhartree)

    deallocate(rho)
  end subroutine v_hartree

  ! ---------------------------------------------------------
  subroutine v_a_xc()
    FLOAT, allocatable :: rho(:, :)
    integer :: is
    call profiling_in(C_PROFILING_XC)

    h%ex = M_ZERO
    h%ec = M_ZERO
    h%exc_j = M_ZERO

    ! get density taking into account non-linear core corrections, and the Amaldi SIC correction
    ALLOCATE(rho(NP, st%d%nspin), NP*st%d%nspin)
    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(1:NP, is) = st%rho(1:NP, is) + st%rho_core(1:NP)/st%d%spin_channels
      end do
    else
      rho(1:NP, :) = st%rho(1:NP, :)
    end if

    ! Amaldi correction
    if(ks%sic_type == sic_amaldi) rho(1:NP,:) = amaldi_factor*rho(1:NP,:)

    ! Get the *local* xc term, which is added in h%vhxc to the Hartree term.
    if(h%d%cdft) then
      call xc_get_vxc_and_axc(gr, ks%xc, rho, st%j, st%d%ispin, h%vhxc, h%axc, &
         h%ex, h%ec, h%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
    else
      call xc_get_vxc(gr, ks%xc, rho, st%d%ispin, h%vhxc, h%ex, h%ec, &
         -minval(st%eigenval(st%nst, :)), st%qtot)
    end if
    deallocate(rho)

    ! The OEP family has to handle specially
    call X(xc_oep_calc)(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
       gr, h, st, h%vhxc, h%ex, h%ec)

    ! Get vxc, by substracting the Hartree term.
    h%vxc = h%vhxc
    h%vxc(1:NP, 1) = h%vxc(1:NP, 1) - h%vhartree(1:NP)
    if(h%d%ispin > UNPOLARIZED) h%vxc(1:NP, 2) = h%vxc(1:NP, 2) - h%vhartree(1:NP)

    ! Now we calculate Int[n vxc] = h%epot
    select case(h%d%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1)) &
         + dmf_dotp(gr%m, st%rho(:, 2), h%vxc(:, 2))
    case(SPINORS)
      h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1)) &
         + dmf_dotp(gr%m, st%rho(:, 2), h%vxc(:, 2)) &
         + M_TWO*dmf_dotp(gr%m, st%rho(:, 3), h%vxc(:, 3)) &
         + M_TWO*dmf_dotp(gr%m, st%rho(:, 4), h%vxc(:, 4))

    end select

    call profiling_out(C_PROFILING_XC)
  end subroutine v_a_xc
end subroutine X(v_ks_calc)
