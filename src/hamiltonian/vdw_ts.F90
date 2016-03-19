!! Copyright (C) 2015 X. Andrade, R. A. DiStasio Jr., B. Santra
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

! This code is based on the Quantum Espresso implementation of TS.

#include "global.h"

module vdw_ts_oct_m
  use derivatives_oct_m
  use geometry_oct_m
  use global_oct_m
  use hirshfeld_oct_m
  use io_function_oct_m
  use messages_oct_m
  use mesh_function_oct_m
  use profiling_oct_m
  use ps_oct_m
  use species_oct_m
  use states_oct_m
  use unit_oct_m
  use unit_system_oct_m
  
  implicit none

  private

  public ::                               &
    vdw_ts_t,                             &
    vdw_ts_init,                          &
    vdw_ts_end,                           &
    vdw_ts_calculate
  
  type vdw_ts_t
    private
    FLOAT, allocatable :: c6free(:)    !> Free atomic volumes for each atomic species.
    FLOAT, allocatable :: dpfree(:)      !> Free atomic static dipole polarizability for each atomic species.
    FLOAT, allocatable :: r0free(:)      !> Free atomic vdW radius for each atomic species.
    FLOAT, allocatable :: c6abfree(:, :) !> Free atomic heteronuclear C6 coefficient for each atom pair.
    FLOAT, allocatable :: volfree(:)
    type(hirshfeld_t) :: hirshfeld
  end type vdw_ts_t

contains

  subroutine vdw_ts_init(this, geo, der, st)
    type(vdw_ts_t),      intent(out)   :: this
    type(geometry_t),    intent(in)    :: geo
    type(derivatives_t), intent(in)    :: der
    type(states_t),      intent(in)    :: st
    
    integer :: ispecies, jspecies
    FLOAT :: num, den

    PUSH_SUB(vdw_ts_init)
    
    SAFE_ALLOCATE(this%c6free(1:geo%nspecies))
    SAFE_ALLOCATE(this%dpfree(1:geo%nspecies))
    SAFE_ALLOCATE(this%r0free(1:geo%nspecies))
    SAFE_ALLOCATE(this%volfree(1:geo%nspecies))
    SAFE_ALLOCATE(this%c6abfree(1:geo%nspecies, 1:geo%nspecies))

    do ispecies = 1, geo%nspecies
      call get_vdw_param(species_label(geo%species(ispecies)), &
        this%c6free(ispecies), this%dpfree(ispecies), this%r0free(ispecies))
      this%volfree(ispecies) = ps_density_volume(species_ps(geo%species(ispecies)))
    end do

    do ispecies = 1, geo%nspecies
      do jspecies = 1, geo%nspecies
        num = CNST(2.0)*this%c6free(ispecies)*this%c6free(jspecies)
        den = (this%dpfree(jspecies)/this%dpfree(ispecies))*this%c6free(ispecies) &
          + (this%dpfree(ispecies)/this%dpfree(jspecies))*this%c6free(jspecies)
        this%c6abfree(ispecies, jspecies) = num/den
      end do
    end do

    call hirshfeld_init(this%hirshfeld, der%mesh, geo, st)

    POP_SUB(vdw_ts_init)
  end subroutine vdw_ts_init

  !------------------------------------------
  
  subroutine vdw_ts_end(this)
    type(vdw_ts_t), intent(inout) :: this

    PUSH_SUB(vdw_ts_end)
    
    call hirshfeld_end(this%hirshfeld)
    
    SAFE_DEALLOCATE_A(this%c6free)
    SAFE_DEALLOCATE_A(this%dpfree)
    SAFE_DEALLOCATE_A(this%r0free)
    SAFE_DEALLOCATE_A(this%volfree)
    SAFE_DEALLOCATE_A(this%c6abfree)

    POP_SUB(vdw_ts_end)
  end subroutine vdw_ts_end

  !------------------------------------------

  subroutine vdw_ts_calculate(this, geo, der, density, energy, potential, force)
    type(vdw_ts_t),      intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(derivatives_t), intent(in)    :: der
    FLOAT,               intent(in)    :: density(:, :)
    FLOAT,               intent(out)   :: energy
    FLOAT,               intent(out)   :: potential(:)
    FLOAT,               intent(out)   :: force(:, :)

    interface 
     subroutine f90_vdw_calculate(natoms, zatom, coordinates, volume_ratio, &
       energy, force, derivative_coeff)
        integer, intent(in)  :: natoms
        integer, intent(in)  :: zatom
        real(8), intent(in)  :: coordinates
        real(8), intent(in)  :: volume_ratio
        real(8), intent(out) :: energy
        real(8), intent(out) :: force
        real(8), intent(out) :: derivative_coeff
      end subroutine f90_vdw_calculate
    end interface

    integer :: iatom, jatom, ispecies, jspecies, ip, idir
    FLOAT :: rr, c6ab, c6abfree, ff, dffdrr, dffdr0
    FLOAT, allocatable :: c6(:), r0(:), volume_ratio(:), dvadens(:), dvadrr(:), coordinates(:,:), derivative_coeff(:)

    integer, allocatable :: zatom(:)

    PUSH_SUB(vdw_ts_calculate)

    SAFE_ALLOCATE(c6(1:geo%natoms))
    SAFE_ALLOCATE(r0(1:geo%natoms))
    SAFE_ALLOCATE(volume_ratio(1:geo%natoms))
    SAFE_ALLOCATE(coordinates(1:3, 1:geo%natoms))
    SAFE_ALLOCATE(zatom(1:geo%natoms))
    SAFE_ALLOCATE(derivative_coeff(1:geo%natoms))
    SAFE_ALLOCATE(dvadens(1:der%mesh%np))
    SAFE_ALLOCATE(dvadrr(1:3))
    
    do iatom = 1, geo%natoms
      ispecies = species_index(geo%atom(iatom)%species)
      call hirshfeld_volume_ratio(this%hirshfeld, iatom, density, volume_ratio(iatom))
      
      c6(iatom) = volume_ratio(iatom)**2*this%c6free(ispecies)
      r0(iatom) = volume_ratio(iatom)**(CNST(1.0)/CNST(3.0))*this%r0free(ispecies)
      
      !print*, species_label(geo%atom(iatom)%species), "c6 ", c6(iatom), this%c6free(ispecies)
      !print*, species_label(geo%atom(iatom)%species), "r0 ", r0(iatom), this%r0free(ispecies)

      coordinates(1:3, iatom) = geo%atom(iatom)%x(1:3)
      !print*, "iatom ", iatom, ": ", coordinates(1:3, iatom)
      zatom(iatom) = species_z(geo%atom(iatom)%species)
      !print*, species_label(geo%atom(iatom)%species),zatom(iatom)

    end do
 
    call f90_vdw_calculate(geo%natoms, zatom(1), coordinates(1, 1), volume_ratio(1), &
      energy, force(1, 1), derivative_coeff(1))

    ! add the extra term to the force
!    force = CNST(0.0)
    do jatom = 1, geo%natoms
      do iatom = 1, geo%natoms
        call hirshfeld_position_derivative(this%hirshfeld, der, iatom, jatom, density, dvadrr)
        force(1:3, jatom) = force(1:3, jatom) - derivative_coeff(iatom)*dvadrr(1:3)
      end do
    end do

    ! and calculate the potential
    potential = CNST(0.0)
    do iatom = 1, geo%natoms
      call hirshfeld_density_derivative(this%hirshfeld, iatom, dvadens)
      potential(1:der%mesh%np) = potential(1:der%mesh%np) + derivative_coeff(iatom)*dvadens(1:der%mesh%np)
    end do

    print*, "fxx", energy, force(1, 1), volume_ratio(1),  derivative_coeff(1)

#if 0
    call dio_function_output(1, "./", "vvdw", der%mesh, potential, unit_one, ip)
    
    print*, dmf_integrate(der%mesh, potential(1:der%mesh%np)*density(1:der%mesh%np, 1))
#endif    

    SAFE_DEALLOCATE_A(c6)
    SAFE_DEALLOCATE_A(r0)
    SAFE_DEALLOCATE_A(volume_ratio)
    SAFE_DEALLOCATE_A(dvadens)
    
    POP_SUB(vdw_ts_calculate)
  end subroutine vdw_ts_calculate

  !------------------------------------------
  
  subroutine get_vdw_param(atom, c6, alpha, r0)
    character(len=*), intent(in)  :: atom
    FLOAT,            intent(out) :: c6
    FLOAT,            intent(out) :: alpha
    FLOAT,            intent(out) :: r0

    PUSH_SUB(get_vdw_param)
    
    select case(trim(atom))

    case('H')
      alpha = CNST(4.500000)
      c6 = CNST(6.500000)
      r0 = CNST(3.100000)

    case('He')
      alpha = CNST(1.380000)
      c6 = CNST(1.460000)
      r0 = CNST(2.650000)
      
    case('Li')
      alpha = CNST(164.200000)
      c6 = CNST(1387.000000)
      r0 = CNST(4.160000)
      
    case('Be')
      alpha = CNST(38.000000)
      c6 = CNST(214.000000)
      r0 = CNST(4.170000)
      
    case('B')
      alpha = CNST(21.000000)
      c6 = CNST(99.500000)
      r0 = CNST(3.890000)
      
    case('C')
      alpha = CNST(12.000000)
      c6 = CNST(46.600000)
      r0 = CNST(3.590000)
      
    case('N')
      alpha = CNST(7.400000)
      c6 = CNST(24.200000)
      r0 = CNST(3.340000)
      
    case('O')
      alpha = CNST(5.400000)
      c6 = CNST(15.600000)
      r0 = CNST(3.190000)
      
    case('F')
      alpha = CNST(3.800000)
      c6 = CNST(9.520000)
      r0 = CNST(3.040000)
      
    case('Ne')
      alpha = CNST(2.670000)
      c6 = CNST(6.380000)
      r0 = CNST(2.910000)
      
    case('Na')
      alpha = CNST(162.700000)
      c6 = CNST(1556.000000)
      r0 = CNST(3.730000)
      
    case('Mg')
      alpha = CNST(71.000000)
      c6 = CNST(627.000000)
      r0 = CNST(4.270000)
      
    case('Al')
      alpha = CNST(60.000000)
      c6 = CNST(528.000000)
      r0 = CNST(4.330000)
      
    case('Si')
      alpha = CNST(37.000000)
      c6 = CNST(305.000000)
      r0 = CNST(4.200000)
      
    case('P')
      alpha = CNST(25.000000)
      c6 = CNST(185.000000)
      r0 = CNST(4.010000)
      
    case('S')
      alpha = CNST(19.600000)
      c6 = CNST(134.000000)
      r0 = CNST(3.860000)
      
    case('Cl')
      alpha = CNST(15.000000)
      c6 = CNST(94.600000)
      r0 = CNST(3.710000)
      
    case('Ar')
      alpha = CNST(11.100000)
      c6 = CNST(64.300000)
      r0 = CNST(3.550000)
      
    case('K')
      alpha = CNST(292.900000)
      c6 = CNST(3897.000000)
      r0 = CNST(3.710000)
      
    case('Ca')
      alpha = CNST(160.000000)
      c6 = CNST(2221.000000)
      r0 = CNST(4.650000)
      
    case('Sc')
      alpha = CNST(120.000000)
      c6 = CNST(1383.000000)
      r0 = CNST(4.590000)
      
    case('Ti')
      alpha = CNST(98.000000)
      c6 = CNST(1044.000000)
      r0 = CNST(4.510000)
      
    case('V')
      alpha = CNST(84.000000)
      c6 = CNST(832.000000)
      r0 = CNST(4.440000)
      
    case('Cr')
      alpha = CNST(78.000000)
      c6 = CNST(602.000000)
      r0 = CNST(3.990000)
      
    case('Mn')
      alpha = CNST(63.000000)
      c6 = CNST(552.000000)
      r0 = CNST(3.970000)
      
    case('Fe')
      alpha = CNST(56.000000)
      c6 = CNST(482.000000)
      r0 = CNST(4.230000)
      
    case('Co')
      alpha = CNST(50.000000)
      c6 = CNST(408.000000)
      r0 = CNST(4.180000)
      
    case('Ni')
      alpha = CNST(48.000000)
      c6 = CNST(373.000000)
      r0 = CNST(3.820000)
      
    case('Cu')
      alpha = CNST(42.000000)
      c6 = CNST(253.000000)
      r0 = CNST(3.760000)
      
    case('Zn')
      alpha = CNST(40.000000)
      c6 = CNST(284.000000)
      r0 = CNST(4.020000)
      
    case('Ga')
      alpha = CNST(60.000000)
      c6 = CNST(498.000000)
      r0 = CNST(4.190000)
      
    case('Ge')
      alpha = CNST(41.000000)
      c6 = CNST(354.000000)
      r0 = CNST(4.200000)
      
    case('As')
      alpha = CNST(29.000000)
      c6 = CNST(246.000000)
      r0 = CNST(4.110000)
      
    case('Se')
      alpha = CNST(25.000000)
      c6 = CNST(210.000000)
      r0 = CNST(4.040000)
      
    case('Br')
      alpha = CNST(20.000000)
      c6 = CNST(162.000000)
      r0 = CNST(3.930000)
      
    case('Kr')
      alpha = CNST(16.800000)
      c6 = CNST(129.600000)
      r0 = CNST(3.820000)
      
    case('Rb')
      alpha = CNST(319.200000)
      c6 = CNST(4691.000000)
      r0 = CNST(3.720000)
      
    case('Sr')
      alpha = CNST(199.000000)
      c6 = CNST(3170.000000)
      r0 = CNST(4.540000)
      
    case('Rh')
      alpha = CNST(56.1)
      c6 = CNST(469.0)
      r0 = CNST(3.95)
      
    case('Pd')
      alpha = CNST(23.680000)
      c6 = CNST(157.500000)
      r0 = CNST(3.66000)
      
    case('Ag')
      alpha = CNST(50.600000)
      c6 = CNST(339.000000)
      r0 = CNST(3.820000)
      
    case('Cd')
      alpha = CNST(39.7)
      c6 = CNST(452.0)
      r0 = CNST(3.99)
      
    case('Te')
      alpha = CNST(37.65)
      c6 = CNST(396.0)
      r0 = CNST(4.22)
      
    case('I')
      alpha = CNST(35.000000)
      c6 = CNST(385.000000)
      r0 = CNST(4.170000)
      
    case('Xe')
      alpha = CNST(27.300000)
      c6 = CNST(285.900000)
      r0 = CNST(4.080000)
      
    case('Ba')
      alpha = CNST(275.0)
      c6 = CNST(5727.0)
      r0 = CNST(4.77)
      
    case('Ir')
      alpha = CNST(42.51)
      c6 = CNST(359.1)
      r0 = CNST(4.00)
      
    case('Pt')
      alpha = CNST(39.68)
      c6 = CNST(347.1)
      r0 = CNST(3.92)
      
    case('Au')
      alpha = CNST(36.5)
      c6 = CNST(298.0)
      r0 = CNST(3.86)
      
    case('Hg')
      alpha = CNST(33.9)
      c6 = CNST(392.0)
      r0 = CNST(3.98)
      
    case('Pb')
      alpha = CNST(61.8)
      c6 = CNST(697.0)
      r0 = CNST(4.31)
      
    case('Bi')
      alpha = CNST(49.02)
      c6 = CNST(571.0)
      r0 = CNST(4.32)
      
    case default
      
      call messages_write('vdw ts: reference free atom parameters not available for species '//trim(atom))
      call messages_fatal()
      
    end select

    POP_SUB(get_vdw_param)
  end subroutine get_vdw_param
  
end module vdw_ts_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
