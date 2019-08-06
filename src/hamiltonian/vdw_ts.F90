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

! This code is based on the Quantum Espresso implementation of TS.

#include "global.h"

module vdw_ts_oct_m
  use derivatives_oct_m
  use geometry_oct_m
  use global_oct_m
  use hirshfeld_oct_m
  use io_oct_m
  use io_function_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use ps_oct_m
  use simul_box_oct_m
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
    vdw_ts_write_c6ab,                    &
    vdw_ts_force_calculate,               &
    vdw_ts_calculate
  
  type vdw_ts_t
    private
    FLOAT, allocatable :: c6free(:)        !> Free atomic volumes for each atomic species.
    FLOAT, allocatable :: dpfree(:)        !> Free atomic static dipole polarizability for each atomic species.
    FLOAT, allocatable :: r0free(:)        !> Free atomic vdW radius for each atomic species.
    FLOAT, allocatable :: c6abfree(:, :)   !> Free atomic heteronuclear C6 coefficient for each atom pair.
    FLOAT, allocatable :: volfree(:)
    FLOAT, allocatable :: c6ab(:, :)       !> Effective atomic heteronuclear C6 coefficient for each atom pair.
    FLOAT              :: cutoff           !> Cutoff value for the calculation of the VdW TS correction in periodic system.
    FLOAT              :: damping          !> Parameter for the damping function steepness.
    FLOAT              :: sr               !> Parameter for the damping function. Can depend on the XC correction used.

    FLOAT, allocatable :: derivative_coeff(:) !> A pre-calculated coefficient for fast derivative evaluation
  end type vdw_ts_t

contains

  subroutine vdw_ts_init(this, namespace, geo, der)
    type(vdw_ts_t),      intent(out)   :: this
    type(namespace_t),   intent(in)    :: namespace
    type(geometry_t),    intent(in)    :: geo
    type(derivatives_t), intent(in)    :: der
    
    integer :: ispecies, jspecies
    FLOAT :: num, den

    PUSH_SUB(vdw_ts_init)
    
    !%Variable VDW_TS_cutoff
    !%Type float
    !%Default 10.0
    !%Section Hamiltonian::XC
    !%Description
    !% Set the value of the cutoff (unit of length) for the VDW correction in periodic system 
    !% in the Tkatchenko and Scheffler (vdw_ts) scheme only. 
    !%End
    call parse_variable(namespace, 'VDW_TS_cutoff', CNST(10.0), this%cutoff, units_inp%length)


    !%Variable VDW_TS_damping
    !%Type float
    !%Default 20.0
    !%Section Hamiltonian::XC
    !%Description
    !% Set the value of the damping function (in unit of 1/length) steepness for the VDW correction in the 
    !% Tkatchenko-Scheffler scheme. See Equation (12) of Phys. Rev. Lett. 102 073005 (2009). 
    !%End
    call parse_variable(namespace, 'VDW_TS_damping', CNST(20.0), this%damping, units_inp%length**(-1))

    !%Variable VDW_TS_sr
    !%Type float
    !%Default 0.94
    !%Section Hamiltonian::XC
    !%Description
    !% Set the value of the sr parameter in the damping function of the VDW correction in the 
    !% Tkatchenko-Scheffler scheme. See Equation (12) of Phys. Rev. Lett. 102 073005 (2009).
    !% This parameter depends on the xc functional used. 
    !% The default value is 0.94, which holds for PBE. For PBE0, a value of 0.96 should be used.
    !%End
    call parse_variable(namespace, 'VDW_TS_sr', CNST(0.94), this%sr)


    SAFE_ALLOCATE(this%c6free(1:geo%nspecies))
    SAFE_ALLOCATE(this%dpfree(1:geo%nspecies))
    SAFE_ALLOCATE(this%r0free(1:geo%nspecies))
    SAFE_ALLOCATE(this%volfree(1:geo%nspecies))
    SAFE_ALLOCATE(this%c6abfree(1:geo%nspecies, 1:geo%nspecies))
    SAFE_ALLOCATE(this%c6ab(1:geo%natoms, 1:geo%natoms))
    SAFE_ALLOCATE(this%derivative_coeff(1:geo%natoms))

    do ispecies = 1, geo%nspecies
      call get_vdw_param(species_label(geo%species(ispecies)), &
        this%c6free(ispecies), this%dpfree(ispecies), this%r0free(ispecies))
      this%volfree(ispecies) = ps_density_volume(species_ps(geo%species(ispecies)))
    end do

    do ispecies = 1, geo%nspecies
      do jspecies = 1, geo%nspecies
        num = M_TWO*this%c6free(ispecies)*this%c6free(jspecies)
        den = (this%dpfree(jspecies)/this%dpfree(ispecies))*this%c6free(ispecies) &
          + (this%dpfree(ispecies)/this%dpfree(jspecies))*this%c6free(jspecies)
        this%c6abfree(ispecies, jspecies) = num/den
      end do
    end do
    POP_SUB(vdw_ts_init)
  end subroutine vdw_ts_init

  !------------------------------------------
  
  subroutine vdw_ts_end(this)
    type(vdw_ts_t), intent(inout) :: this

    PUSH_SUB(vdw_ts_end)
    
    SAFE_DEALLOCATE_A(this%c6free)
    SAFE_DEALLOCATE_A(this%dpfree)
    SAFE_DEALLOCATE_A(this%r0free)
    SAFE_DEALLOCATE_A(this%volfree)
    SAFE_DEALLOCATE_A(this%c6abfree)
    SAFE_DEALLOCATE_A(this%c6ab)
    SAFE_DEALLOCATE_A(this%derivative_coeff)

    POP_SUB(vdw_ts_end)
  end subroutine vdw_ts_end

  !------------------------------------------

  subroutine vdw_ts_calculate(this, namespace, geo, der, sb, st, density, energy, potential, force)
    type(vdw_ts_t),      intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(geometry_t),    intent(in)    :: geo
    type(derivatives_t), intent(in)    :: der
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    FLOAT,               intent(in)    :: density(:, :)
    FLOAT,               intent(out)   :: energy
    FLOAT,               intent(out)   :: potential(:)
    FLOAT,               intent(out)   :: force(:, :)

    interface
      subroutine f90_vdw_calculate(natoms, dd, sr, zatom, coordinates, vol_ratio, &
        energy, force, derivative_coeff)
        integer, intent(in)  :: natoms
        real(8), intent(in)  :: dd
        real(8), intent(in)  :: sr
        integer, intent(in)  :: zatom
        real(8), intent(in)  :: coordinates
        real(8), intent(in)  :: vol_ratio
        real(8), intent(out) :: energy
        real(8), intent(out) :: force
        real(8), intent(out) :: derivative_coeff
      end subroutine f90_vdw_calculate
    end interface

    type(periodic_copy_t) :: pc
    integer :: iatom, jatom, ispecies, jspecies, jcopy, ip 
    FLOAT :: rr, rr2, rr6, dffdr0, ee, ff, dee, dffdrab, dffdvra, deabdvra
    FLOAT, allocatable :: coordinates(:,:), vol_ratio(:), dvadens(:), dvadrr(:), & 
                          dr0dvra(:), r0ab(:,:)
    type(hirshfeld_t) :: hirshfeld
    integer, allocatable :: zatom(:)
    FLOAT :: x_j(1:MAX_DIM)

    PUSH_SUB(vdw_ts_calculate)

    SAFE_ALLOCATE(vol_ratio(1:geo%natoms))
    SAFE_ALLOCATE(dvadens(1:der%mesh%np))
    SAFE_ALLOCATE(dvadrr(1:3))
    SAFE_ALLOCATE(dr0dvra(1:geo%natoms))

    energy=M_ZERO
    force(1:sb%dim, 1:geo%natoms) = M_ZERO
    this%derivative_coeff(1:geo%natoms) = M_ZERO
    call hirshfeld_init(hirshfeld, namespace, der%mesh, geo, st)

    do iatom = 1, geo%natoms
      call hirshfeld_volume_ratio(hirshfeld, iatom, density, vol_ratio(iatom))
    end do

    do iatom = 1, geo%natoms
      ispecies = species_index(geo%atom(iatom)%species)
      dr0dvra(iatom) = this%r0free(ispecies)/(CNST(3.0)*(vol_ratio(iatom)**(M_TWO/CNST(3.0))))
      do jatom = 1, geo%natoms
        jspecies = species_index(geo%atom(jatom)%species)
        this%c6ab(iatom,jatom) = vol_ratio(iatom)*vol_ratio(jatom)*this%c6abfree(ispecies,jspecies) !this operation is done again inside the .c part for the non periodic case
      end do
    end do
  
    if(sb%periodic_dim > 0) then ! periodic case
      SAFE_ALLOCATE(r0ab(1:geo%natoms,1:geo%natoms))

      !Precomputing some quantities
      do iatom = 1, geo%natoms
        ispecies = species_index(geo%atom(iatom)%species)
        do jatom = 1, geo%natoms
         jspecies = species_index(geo%atom(jatom)%species)

         r0ab(iatom,jatom) = (vol_ratio(iatom)**(M_ONE/CNST(3.0)))*this%r0free(ispecies) &
                           + (vol_ratio(jatom)**(M_ONE/CNST(3.0)))*this%r0free(jspecies)
        end do
      end do


      do jatom = 1, geo%natoms
        jspecies = species_index(geo%atom(jatom)%species)
                
        call periodic_copy_init(pc, sb, geo%atom(jatom)%x, this%cutoff)
        do jcopy = 1, periodic_copy_num(pc) ! one of the periodic copy is the initial atom  
          x_j(1:sb%dim) = periodic_copy_position(pc, sb, jcopy)
          do iatom = 1, geo%natoms
            ispecies = species_index(geo%atom(iatom)%species) 
            rr2 =  sum( (x_j(1:sb%dim) - geo%atom(iatom)%x(1:sb%dim))**2 )
            rr =  sqrt(rr2)
            rr6 = rr2**3

            if(rr < CNST(1.0e-10)) cycle !To avoid self interaction

            ee = exp( -this%damping*((rr/( this%sr*r0ab(iatom, jatom))) - M_ONE))
            ff = M_ONE/(M_ONE + ee)
            dee = ee*ff**2

            !Calculate the derivative of the damping function with respect to the distance between atoms A and B.
            dffdrab = (this%damping/(this%sr*r0ab(iatom, jatom)))*dee
            !Calculate the derivative of the damping function with respect to the distance between the van der Waals radius.
            dffdr0 =  -this%damping*rr/(this%sr*r0ab(iatom, jatom)**2)*dee

            energy = energy - M_HALF*ff*this%c6ab(iatom, jatom)/rr6

            ! Derivative of the damping function with respecto to the volume ratio of atom A.
            dffdvra = dffdr0*dr0dvra(iatom); ! Ces termes sont bon

            ! Calculation of the pair-wise partial energy derivative with respect to the volume ratio of atom A.
            deabdvra = (dffdvra*this%c6ab(iatom, jatom) + ff*vol_ratio(jatom)*this%c6abfree(ispecies, jspecies))/rr6 
               
            this%derivative_coeff(iatom) = this%derivative_coeff(iatom) + deabdvra;

          end do
        end do
        call periodic_copy_end(pc)
      end do

      SAFE_DEALLOCATE_A(r0ab)
    else ! Non periodic case 
      SAFE_ALLOCATE(coordinates(1:sb%dim, 1:geo%natoms))
      SAFE_ALLOCATE(zatom(1:geo%natoms))

      do iatom = 1, geo%natoms
        coordinates(1:sb%dim, iatom) = geo%atom(iatom)%x(1:sb%dim)
        zatom(iatom) = species_z(geo%atom(iatom)%species)

      end do
      
      call f90_vdw_calculate(geo%natoms,  this%damping, this%sr, zatom(1), coordinates(1, 1), &
                             vol_ratio(1), energy, force(1, 1), this%derivative_coeff(1))


      SAFE_DEALLOCATE_A(coordinates)
      SAFE_DEALLOCATE_A(zatom)
    end if

    ! Calculate the potential
    potential = M_ZERO
    do iatom = 1, geo%natoms
      call hirshfeld_density_derivative(hirshfeld, iatom, dvadens)
      potential(1:der%mesh%np) = potential(1:der%mesh%np) - this%derivative_coeff(iatom)*dvadens(1:der%mesh%np) 
    end do

    if(debug%info) then
      call dio_function_output(1_8, "./", "vvdw", namespace, der%mesh, potential, unit_one, ip)
    end if

    call hirshfeld_end(hirshfeld)

    SAFE_DEALLOCATE_A(vol_ratio)
    SAFE_DEALLOCATE_A(dvadens)
    SAFE_DEALLOCATE_A(dvadrr)
    SAFE_DEALLOCATE_A(dr0dvra)

    POP_SUB(vdw_ts_calculate)
    end subroutine vdw_ts_calculate



  !------------------------------------------
  subroutine vdw_ts_force_calculate(this, namespace, force_vdw, geo, der, sb, st, density)
    type(vdw_ts_t),      intent(in)    :: this
    type(namespace_t),   intent(in)    :: namespace
    FLOAT,               intent(inout) :: force_vdw(:,:)
    type(geometry_t),    intent(in)    :: geo
    type(derivatives_t), intent(in)    :: der
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    FLOAT,               intent(in)    :: density(:, :)

    type(hirshfeld_t) :: hirshfeld
    type(periodic_copy_t) :: pc

    integer :: iatom, jatom, ispecies, jspecies, jcopy
    FLOAT :: rr, rr2, rr6,  dffdr0, ee, ff, dee, dffdvra, deabdvra, deabdrab, x_j(1:MAX_DIM) 
    FLOAT, allocatable ::  vol_ratio(:), dvadrr(:), dr0dvra(:), r0ab(:,:), derivative_coeff(:), c6ab(:,:)

    PUSH_SUB(vdw_ts_force_calculate)


    SAFE_ALLOCATE(vol_ratio(1:geo%natoms))
    SAFE_ALLOCATE(dvadrr(1:3))
    SAFE_ALLOCATE(dr0dvra(1:geo%natoms))
    SAFE_ALLOCATE(r0ab(1:geo%natoms,1:geo%natoms))
    SAFE_ALLOCATE(derivative_coeff(1:geo%natoms))
    SAFE_ALLOCATE(c6ab(1:geo%natoms,1:geo%natoms))


    force_vdw(1:sb%dim, 1:geo%natoms) = M_ZERO
    derivative_coeff(1:geo%natoms) = M_ZERO
    c6ab(1:geo%natoms,1:geo%natoms) = M_ZERO
    r0ab(1:geo%natoms,1:geo%natoms) = M_ZERO
    dr0dvra(1:geo%natoms) = M_ZERO
    dvadrr(1:3) = M_ZERO
    vol_ratio(1:geo%natoms) = M_ZERO


    call hirshfeld_init(hirshfeld, namespace, der%mesh, geo, st)


    do iatom = 1, geo%natoms
      call hirshfeld_volume_ratio(hirshfeld, iatom, density, vol_ratio(iatom))
    end do

    do iatom = 1, geo%natoms
      ispecies = species_index(geo%atom(iatom)%species)
      dr0dvra(iatom) = this%r0free(ispecies)/(CNST(3.0)*(vol_ratio(iatom)**(M_TWO/CNST(3.0))))
      do jatom = 1, geo%natoms
        jspecies = species_index(geo%atom(jatom)%species)
        c6ab(iatom, jatom) = vol_ratio(iatom)*vol_ratio(jatom)*this%c6abfree(ispecies, jspecies) 
      end do
    end do

    !Precomputing some quantities
    do iatom = 1, geo%natoms
      ispecies = species_index(geo%atom(iatom)%species)
      do jatom = iatom, geo%natoms
       jspecies = species_index(geo%atom(jatom)%species)

       r0ab(iatom, jatom) = (vol_ratio(iatom)**(M_ONE/CNST(3.0)))*this%r0free(ispecies) &
                          + (vol_ratio(jatom)**(M_ONE/CNST(3.0)))*this%r0free(jspecies)
       if(iatom /= jatom) r0ab(jatom, iatom) = r0ab(iatom, jatom)
      end do
    end do


    do jatom = 1, geo%natoms
      jspecies = species_index(geo%atom(jatom)%species)
      
      call periodic_copy_init(pc, sb, geo%atom(jatom)%x, this%cutoff)
      do jcopy = 1, periodic_copy_num(pc) ! one of the periodic copy is the initial atom  
        x_j(1:sb%dim) = periodic_copy_position(pc, sb, jcopy)
        do iatom = 1, geo%natoms
          ispecies = species_index(geo%atom(iatom)%species)
          rr2 =  sum((x_j(1:sb%dim) - geo%atom(iatom)%x(1:sb%dim))**2)
          rr  =  sqrt(rr2)
          rr6 = rr2**3

          if(rr < TOL_HIRSHFELD) cycle !To avoid self interaction

          ee = exp(-this%damping*(rr/(this%sr*r0ab(iatom, jatom)) - M_ONE))
          ff = M_ONE/(M_ONE + ee)
          dee = ee*ff**2
          !Calculate the derivative of the damping function with respect to the van der Waals radius.
          dffdr0 =  -this%damping*rr/( this%sr*r0ab(iatom, jatom)**2)*dee
          ! Calculation of the pair-wise partial energy derivative with respect to the distance between atoms A and B.
          deabdrab = c6ab(iatom,jatom)*(this%damping/(this%sr*r0ab(iatom, jatom))*dee - CNST(6.0)*ff/rr)/rr6;
          ! Derivative of the damping function with respecto to the volume ratio of atom A.
          dffdvra = dffdr0*dr0dvra(iatom);
          ! Calculation of the pair-wise partial energy derivative with respect to the volume ratio of atom A.
          deabdvra = (dffdvra*c6ab(iatom, jatom) + ff*vol_ratio(jatom)*this%c6abfree(ispecies, jspecies))/rr6
          !Summing for using later
          derivative_coeff(iatom) = derivative_coeff(iatom) + deabdvra;

          ! Calculation of the pair-wise partial energy derivative with respect to the distance between atoms A and B.
          deabdrab = c6ab(iatom, jatom)*(this%damping/(this%sr*r0ab(iatom, jatom))*dee - CNST(6.0)*ff/rr)/rr6;
          force_vdw(1:sb%dim, iatom)= force_vdw(1:sb%dim, iatom) + M_HALF*deabdrab*(geo%atom(iatom)%x(1:sb%dim) - x_j(1:sb%dim))/rr;
        end do
      end do
      call periodic_copy_end(pc)
    end do

    do iatom = 1, geo%natoms
      do jatom = 1, geo%natoms
        call hirshfeld_position_derivative(hirshfeld, der, iatom, jatom, density, dvadrr) !dvadrr_ij = \frac{\delta V_i}{\delta \vec{x_j}}
        force_vdw(1:sb%dim, jatom)= force_vdw(1:sb%dim, jatom) + derivative_coeff(iatom)*dvadrr(1:sb%dim)  ! geo%atom(jatom)%f_vdw(1:sb%dim) = sum_i coeff_i * dvadrr_ij
      end do
    end do

    call hirshfeld_end(hirshfeld)

    SAFE_DEALLOCATE_A(vol_ratio)
    SAFE_DEALLOCATE_A(dvadrr)
    SAFE_DEALLOCATE_A(dr0dvra)
    SAFE_DEALLOCATE_A(r0ab)
    SAFE_DEALLOCATE_A(derivative_coeff)
    SAFE_DEALLOCATE_A(c6ab)

    POP_SUB(vdw_ts_force_calculate)
  end subroutine vdw_ts_force_calculate

  !------------------------------------------

  subroutine vdw_ts_write_c6ab(this, geo, dir, fname)
     type(vdw_ts_t)  , intent(inout) :: this
     type(geometry_t),    intent(in) :: geo
     character(len=*), intent(in)    :: dir, fname
 
     integer :: iunit, iatom, jatom

     PUSH_SUB(vdw_ts_write_c6ab)

     if(mpi_grp_is_root(mpi_world)) then  
       call io_mkdir_old(dir)
       iunit = io_open_old(trim(dir) // "/" // trim(fname), action='write')  
        write(iunit, '(a)') ' # Atom1 Atom2 C6_{12}^{eff}'


       do iatom = 1, geo%natoms
         do jatom = 1, geo%natoms
           write(iunit, '(3x, i5, i5, e15.6)') iatom, jatom, this%c6ab(iatom, jatom)
         end do
       end do
       call io_close(iunit)
       end if
     POP_SUB(vdw_ts_write_c6ab)

  end subroutine vdw_ts_write_c6ab




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
