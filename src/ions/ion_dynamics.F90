!! Copyright (C) 2008 X. Andrade
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

#include "global.h"

module ion_dynamics_oct_m
  use iso_c_binding
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use read_coords_oct_m
  use simul_box_oct_m
  use species_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                                &
    ion_dynamics_t,                        &
    ion_state_t,                           &
    ion_dynamics_init,                     &
    ion_dynamics_end,                      &
    ion_dynamics_propagate,                &
    ion_dynamics_propagate_vel,            &
    ion_dynamics_save_state,               &
    ion_dynamics_restore_state,            &
    ion_dynamics_ions_move,                &
    ion_dynamics_temperature,              &
    ion_dynamics_kinetic_energy,           &
    ion_dynamics_freeze,                   &
    ion_dynamics_unfreeze,                 &
    ion_dynamics_verlet_step1,             &
    ion_dynamics_verlet_step2

  type ion_dynamics_t
    private
    logical          :: move_ions
    logical          :: constant_velocity
    integer          :: thermostat
    FLOAT            :: dt
    FLOAT            :: current_temperature
    FLOAT, pointer   :: oldforce(:, :)
    FLOAT, pointer   :: old_pos(:, :)    
  end type ion_dynamics_t

  type ion_state_t
    private
    FLOAT, pointer :: pos(:, :)
    FLOAT, pointer :: vel(:, :)
    FLOAT, pointer :: old_pos(:, :)
  end type ion_state_t  

contains

  ! ---------------------------------------------------------
  subroutine ion_dynamics_init(this, geo)
    type(ion_dynamics_t), intent(out)   :: this
    type(geometry_t),     intent(inout) :: geo

    integer :: i, j
    FLOAT   :: x(MAX_DIM), temperature, sigma, kin1, kin2
    type(c_ptr) :: random_gen_pointer
    type(read_coords_info) :: xyz
    logical :: have_velocities

    PUSH_SUB(ion_dynamics_init)

    nullify(this%oldforce)

    have_velocities = .false.

    !%Variable IonsConstantVelocity
    !%Type logical
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% (Experimental) If this variable is set to yes, the ions will
    !% move with a constant velocity given by the initial
    !% conditions. They will not be affected by any forces.
    !%End
    call parse_variable('IonsConstantVelocity', .false., this%constant_velocity)
    call messages_print_var_value(stdout, 'IonsConstantVelocity', this%constant_velocity)

    if(this%constant_velocity) then
      call messages_experimental('IonsConstantVelocity')
      have_velocities = .true.
    end if
        
    !%Variable RandomVelocityTemp
    !%Type float
    !%Default 0.0
    !%Section System::Velocities
    !%Description
    !% If this variable is present, <tt>Octopus</tt> will assign random
    !% velocities to the atoms following a Boltzmann distribution with
    !% temperature given by <tt>RandomVelocityTemp</tt> (in degrees Kelvin).
    !%End
    if(parse_is_defined('RandomVelocityTemp')) then

      have_velocities = .true.

      if( mpi_grp_is_root(mpi_world)) then
        call loct_ran_init(random_gen_pointer)
        call parse_variable('RandomVelocityTemp', M_ZERO, temperature, unit = unit_kelvin)
      end if

      do i = 1, geo%natoms
        !generate the velocities in the root node
        if( mpi_grp_is_root(mpi_world)) then
          sigma = sqrt(temperature / species_mass(geo%atom(i)%species) )
          do j = 1, 3
             geo%atom(i)%v(j) = loct_ran_gaussian(random_gen_pointer, sigma)
          end do
        end if
#ifdef HAVE_MPI
        !and send them to the others
        call MPI_Bcast(geo%atom(i)%v, geo%space%dim, MPI_FLOAT, 0, mpi_world%comm, mpi_err)
#endif
      end do

      if( mpi_grp_is_root(mpi_world)) then
        call loct_ran_end(random_gen_pointer)
      end if

      kin1 = ion_dynamics_kinetic_energy(geo)

      call cm_vel(geo, x)
      do i = 1, geo%natoms
        geo%atom(i)%v = geo%atom(i)%v - x
      end do

      kin2 = ion_dynamics_kinetic_energy(geo)

      do i = 1, geo%natoms
        geo%atom(i)%v(:) =  sqrt(kin1/kin2)*geo%atom(i)%v(:)
      end do

      write(message(1),'(a,f10.4,1x,a)') 'Info: Initial velocities randomly distributed with T =', &
        units_from_atomic(unit_kelvin, temperature), units_abbrev(unit_kelvin)
      write(message(2),'(2x,a,f8.4,1x,a)') '<K>       =', &
        units_from_atomic(units_out%energy, ion_dynamics_kinetic_energy(geo)/geo%natoms), &
        units_abbrev(units_out%energy)
      write(message(3),'(2x,a,f8.4,1x,a)') '3/2 k_B T =', &
        units_from_atomic(units_out%energy, (M_THREE/M_TWO)*temperature), &
        units_abbrev(units_out%energy)
      call messages_info(3)

    else
      !%Variable XYZVelocities
      !%Type string
      !%Section System::Velocities
      !%Description
      !% <tt>Octopus</tt> will try to read the starting velocities of the atoms from the XYZ file 
      !% specified by the variable <tt>XYZVelocities</tt>.
      !% Note that you do not need to specify initial velocities if you are not going
      !% to perform ion dynamics; if you are going to allow the ions to move but the velocities
      !% are not specified, they are considered to be null.
      !% Note: It is important for the velocities to maintain the ordering 
      !% in which the atoms were defined in the coordinates specifications.
      !%End

      !%Variable XSFVelocities
      !%Type string
      !%Section System::Velocities
      !%Description
      !% Like <tt>XYZVelocities</tt> but in XCrySDen format, as in <tt>XSFCoordinates</tt>.
      !%End

      !%Variable PDBVelocities
      !%Type string
      !%Section System::Velocities
      !%Description
      !% Like <tt>XYZVelocities</tt> but in PDB format, as in <tt>PDBCoordinates</tt>.
      !%End

      !%Variable Velocities
      !%Type block
      !%Section System::Velocities
      !%Description
      !% If <tt>XYZVelocities</tt>, <tt>PDBVelocities</tt>, and <tt>XSFVelocities</tt>
      !% are not present, <tt>Octopus</tt> will try to fetch the initial 
      !% atomic velocities from this block. If this block is not present, <tt>Octopus</tt>
      !% will set the initial velocities to zero. The format of this block can be
      !% illustrated by this example:
      !%
      !% <tt>%Velocities
      !% <br>&nbsp;&nbsp;'C'  |      -1.7 | 0.0 | 0.0
      !% <br>&nbsp;&nbsp;'O'  | &nbsp;1.7 | 0.0 | 0.0
      !% <br>%</tt>
      !%
      !% It describes one carbon and one oxygen moving at the relative
      !% velocity of 3.4 velocity units.
      !%
      !% Note: It is important for the velocities to maintain the ordering 
      !% in which the atoms were defined in the coordinates specifications.
      !%End

      call read_coords_init(xyz)
      call read_coords_read('Velocities', xyz, geo%space)
      if(xyz%source /= READ_COORDS_ERR) then
        
        have_velocities = .true.

        if(geo%natoms /= xyz%n) then
          write(message(1), '(a,i4,a,i4)') 'I need exactly ', geo%natoms, ' velocities, but I found ', xyz%n
          call messages_fatal(1)
        end if

        ! copy information and adjust units
        do i = 1, geo%natoms
          geo%atom(i)%v = units_to_atomic(units_inp%velocity/units_inp%length, xyz%atom(i)%x)
        end do

        call read_coords_end(xyz)

      else
        do i = 1, geo%natoms
          geo%atom(i)%v = M_ZERO
        end do
      end if
    end if

    geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)

    !%Variable MoveIons
    !%Type logical
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable controls whether atoms are moved during a time
    !% propagation run. The default is yes when the ion velocity is
    !% set explicitly or implicitly, otherwise is no.
    !%End
    call parse_variable('MoveIons', have_velocities, this%move_ions)
    call messages_print_var_value(stdout, 'MoveIons', this%move_ions)

    if(ion_dynamics_ions_move(this)) then 
      SAFE_ALLOCATE(this%oldforce(1:geo%space%dim, 1:geo%natoms))
    end if

    POP_SUB(ion_dynamics_init)
  end subroutine ion_dynamics_init

  ! ---------------------------------------------------------
  subroutine ion_dynamics_end(this)
    type(ion_dynamics_t), intent(inout) :: this

    PUSH_SUB(ion_dynamics_end)
    SAFE_DEALLOCATE_P(this%oldforce)

    POP_SUB(ion_dynamics_end)
  end subroutine ion_dynamics_end

  ! ---------------------------------------------------------
  subroutine ion_dynamics_propagate(this, sb, geo, time, dt)
    type(ion_dynamics_t), intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(geometry_t),     intent(inout) :: geo
    FLOAT,                intent(in)    :: time
    FLOAT,                intent(in)    :: dt

    integer :: iatom
    FLOAT   :: DR(1:3)

    if(.not. ion_dynamics_ions_move(this)) return

    PUSH_SUB(ion_dynamics_propagate)
    
    DR = M_ZERO

    this%dt = dt

    this%current_temperature = CNST(0.0)
    
    ! integrate using verlet
    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle
      
      geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) &
        + dt*geo%atom(iatom)%v(1:geo%space%dim) + &
        M_HALF*dt**2 / species_mass(geo%atom(iatom)%species) * geo%atom(iatom)%f(1:geo%space%dim)
      this%oldforce(1:geo%space%dim, iatom) = geo%atom(iatom)%f(1:geo%space%dim)
    end do

    call simul_box_atoms_in_box(sb, geo, .false.)

    POP_SUB(ion_dynamics_propagate)
  end subroutine ion_dynamics_propagate

  ! ---------------------------------------------------------
  subroutine ion_dynamics_propagate_vel(this, geo)
    type(ion_dynamics_t), intent(inout) :: this
    type(geometry_t),     intent(inout) :: geo

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    PUSH_SUB(ion_dynamics_propagate_vel)

    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle
      
      geo%atom(iatom)%v(1:geo%space%dim) = geo%atom(iatom)%v(1:geo%space%dim) &
        + this%dt/species_mass(geo%atom(iatom)%species) * M_HALF * (this%oldforce(1:geo%space%dim, iatom) + &
        geo%atom(iatom)%f(1:geo%space%dim))
    end do

    POP_SUB(ion_dynamics_propagate_vel)
  end subroutine ion_dynamics_propagate_vel

  ! ---------------------------------------------------------
  !> A bare verlet integrator.
  subroutine ion_dynamics_verlet_step1(geo, q, v, fold, dt)
    type(geometry_t),     intent(in)    :: geo
    FLOAT,                intent(inout) :: q(:, :)
    FLOAT,                intent(inout) :: v(:, :)
    FLOAT,                intent(in)    :: fold(:, :)
    FLOAT,                intent(in)    :: dt

    integer :: iatom

    PUSH_SUB(ion_dynamics_verlet_step1)

    ! First transform momenta to velocities
    do iatom = 1, geo%natoms
      v(iatom, 1:geo%space%dim) = v(iatom, 1:geo%space%dim) / species_mass(geo%atom(iatom)%species)
    end do

    ! integrate using verlet
    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle
      q(iatom, 1:geo%space%dim) = q(iatom, 1:geo%space%dim) &
        + dt * v(iatom, 1:geo%space%dim) + &
        M_HALF*dt**2 / species_mass(geo%atom(iatom)%species) * fold(iatom, 1:geo%space%dim)
    end do

    ! And back to momenta.
    do iatom = 1, geo%natoms
      v(iatom, 1:geo%space%dim) = species_mass(geo%atom(iatom)%species) * v(iatom, 1:geo%space%dim)
    end do

    POP_SUB(ion_dynamics_verlet_step1)
  end subroutine ion_dynamics_verlet_step1

  ! ---------------------------------------------------------
  !> A bare verlet integrator.
  subroutine ion_dynamics_verlet_step2(geo, v, fold, fnew, dt)
    type(geometry_t),     intent(in)    :: geo
    FLOAT,                intent(inout) :: v(:, :)
    FLOAT,                intent(in)    :: fold(:, :)
    FLOAT,                intent(in)    :: fnew(:, :)
    FLOAT,                intent(in)    :: dt

    integer :: iatom

    PUSH_SUB(ion_dynamics_verlet_step2)

    ! First transform momenta to velocities
    do iatom = 1, geo%natoms
      v(iatom, 1:geo%space%dim) = v(iatom, 1:geo%space%dim) / species_mass(geo%atom(iatom)%species)
    end do

    ! velocity verlet
    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle
      v(iatom, 1:geo%space%dim) = v(iatom, 1:geo%space%dim) &
        + dt / species_mass(geo%atom(iatom)%species) * M_HALF * (fold(iatom, 1:geo%space%dim) + &
        fnew(iatom, 1:geo%space%dim))
    end do

    ! And back to momenta.
    do iatom = 1, geo%natoms
      v(iatom, 1:geo%space%dim) = species_mass(geo%atom(iatom)%species) * v(iatom, 1:geo%space%dim)
    end do

    POP_SUB(ion_dynamics_verlet_step2)
  end subroutine ion_dynamics_verlet_step2

  ! ---------------------------------------------------------
  subroutine ion_dynamics_save_state(this, geo, state)
    type(ion_dynamics_t), intent(in)    :: this
    type(geometry_t),     intent(in)    :: geo
    type(ion_state_t),    intent(out)   :: state

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    PUSH_SUB(ion_dynamics_save_state)

    SAFE_ALLOCATE(state%pos(1:geo%space%dim, 1:geo%natoms))
    SAFE_ALLOCATE(state%vel(1:geo%space%dim, 1:geo%natoms))

    do iatom = 1, geo%natoms
      state%pos(1:geo%space%dim, iatom) = geo%atom(iatom)%x(1:geo%space%dim)
      state%vel(1:geo%space%dim, iatom) = geo%atom(iatom)%v(1:geo%space%dim)
    end do
    
    POP_SUB(ion_dynamics_save_state)
  end subroutine ion_dynamics_save_state

  ! ---------------------------------------------------------
  subroutine ion_dynamics_restore_state(this, geo, state)
    type(ion_dynamics_t), intent(inout) :: this
    type(geometry_t),     intent(inout) :: geo
    type(ion_state_t),    intent(inout) :: state

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    PUSH_SUB(ion_dynamics_restore_state)

    do iatom = 1, geo%natoms
      geo%atom(iatom)%x(1:geo%space%dim) = state%pos(1:geo%space%dim, iatom)
      geo%atom(iatom)%v(1:geo%space%dim) = state%vel(1:geo%space%dim, iatom)
    end do

    SAFE_DEALLOCATE_P(state%pos)
    SAFE_DEALLOCATE_P(state%vel)
    
    POP_SUB(ion_dynamics_restore_state)
  end subroutine ion_dynamics_restore_state

  ! ---------------------------------------------------------
  logical pure function ion_dynamics_ions_move(this) result(ions_move)
    type(ion_dynamics_t), intent(in)    :: this
    
    ions_move = this%move_ions
  end function ion_dynamics_ions_move

  ! ---------------------------------------------------------
  FLOAT pure function ion_dynamics_kinetic_energy(geo) result(kinetic_energy)
    type(geometry_t),      intent(in) :: geo

    integer :: iatom

    kinetic_energy = M_ZERO
    do iatom = 1, geo%natoms
      kinetic_energy = kinetic_energy + &
        M_HALF * species_mass(geo%atom(iatom)%species) * sum(geo%atom(iatom)%v(1:geo%space%dim)**2)
    end do
  end function ion_dynamics_kinetic_energy

  ! ---------------------------------------------------------
  !> This function returns the ionic temperature in energy units.
  FLOAT pure function ion_dynamics_temperature(geo) result(temperature)
    type(geometry_t),      intent(in) :: geo

    temperature = CNST(2.0)/CNST(3.0)*ion_dynamics_kinetic_energy(geo)/geo%natoms
  end function ion_dynamics_temperature

  ! ---------------------------------------------------------
  !> Freezes the ionic movement.
  logical function ion_dynamics_freeze(this) result(freeze)
    type(ion_dynamics_t), intent(inout)   :: this
    if(this%move_ions) then
      this%move_ions = .false.
      freeze = .true.
    else
      freeze = .false.
    end if
  end function ion_dynamics_freeze

  ! ---------------------------------------------------------
  !> Unfreezes the ionic movement.
  subroutine ion_dynamics_unfreeze(this)
    type(ion_dynamics_t), intent(inout)   :: this
    this%move_ions = .true.
  end subroutine ion_dynamics_unfreeze

end module ion_dynamics_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
