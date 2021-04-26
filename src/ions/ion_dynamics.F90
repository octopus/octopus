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
  use loct_math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use read_coords_oct_m
  use simul_box_oct_m
  use species_oct_m
  use tdfunction_oct_m
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

  integer, parameter ::   &
    THERMO_NONE     = 0,  &
    THERMO_SCAL     = 1,  &
    THERMO_NH       = 2

  type nose_hoover_t
    private
    FLOAT :: mass
    FLOAT :: pos
    FLOAT :: vel
  end type nose_hoover_t

  type ion_td_displacement_t
    private
    logical     :: move
    type(tdf_t) :: fx 
    type(tdf_t) :: fy 
    type(tdf_t) :: fz       
  end type ion_td_displacement_t

  type ion_dynamics_t
    private
    logical          :: move_ions
    logical          :: constant_velocity
    integer          :: thermostat
    FLOAT            :: dt
    FLOAT            :: current_temperature

    FLOAT, allocatable :: oldforce(:, :)

    !> the old positions for Verlet (used for the Nose-Hoover)
    FLOAT, allocatable :: old_pos(:, :)    

    !> variables for the Nose-Hoover thermostat
    type(nose_hoover_t) :: nh(1:2)
    type(tdf_t) :: temperature_function
      
    logical :: drive_ions  
    type(ion_td_displacement_t), allocatable ::  td_displacements(:) !> Time-dependent displacements driving the ions
    type(geometry_t), pointer :: geo_t0
  end type ion_dynamics_t

  type ion_state_t
    private
    FLOAT, allocatable :: pos(:, :)
    FLOAT, allocatable :: vel(:, :)
    FLOAT, allocatable :: old_pos(:, :)
    type(nose_hoover_t) :: nh(1:2)
  end type ion_state_t  

contains

  ! ---------------------------------------------------------
  subroutine ion_dynamics_init(this, namespace, geo)
    type(ion_dynamics_t), intent(out)   :: this
    type(namespace_t),    intent(in)    :: namespace
    type(geometry_t),     intent(inout) :: geo

    integer :: i, j, iatom, ierr
    FLOAT   :: x(MAX_DIM), temperature, sigma, kin1, kin2
    type(c_ptr) :: random_gen_pointer
    type(read_coords_info) :: xyz
    character(len=100)  :: temp_function_name
    logical :: have_velocities

    type(block_t)      :: blk
    integer            :: ndisp
    character(len=200) :: expression

    PUSH_SUB(ion_dynamics_init)

    have_velocities = .false.
    this%drive_ions = .false.

    !%Variable IonsConstantVelocity
    !%Type logical
    !%Default no
    !%Section Time-Dependent::Propagation
    !%Description
    !% (Experimental) If this variable is set to yes, the ions will
    !% move with a constant velocity given by the initial
    !% conditions. They will not be affected by any forces.
    !%End
    call parse_variable(namespace, 'IonsConstantVelocity', .false., this%constant_velocity)
    call messages_print_var_value(stdout, 'IonsConstantVelocity', this%constant_velocity)

    if(this%constant_velocity) then
      call messages_experimental('IonsConstantVelocity')
      have_velocities = .true.
      this%drive_ions = .true.
    end if
    
    !%Variable IonsTimeDependentDisplacements
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% (Experimental) This variable allows you to specify a
    !% time-dependent function describing the displacement of the ions
    !% from their equilibrium position: <math>r(t) = r_0 + \Delta
    !% r(t)</math>.  Specify the displacements dx(t), dy(t), dz(t) as
    !% follows, for some or all of the atoms:
    !% 
    !% <tt>%IonsTimeDependentDisplacements
    !% <br>&nbsp;&nbsp; atom_index | "dx(t)" | "dy(t)" | "dz(t)"
    !% <br>%</tt>
    !%
    !% The displacement functions are time-dependent functions and should match one
    !% of the function names given in the first column of the <tt>TDFunctions</tt> block.
    !% If this block is set, the ions will not be affected by any forces.
    !%End

    
    ndisp = 0
    if(parse_block(namespace, 'IonsTimeDependentDisplacements', blk) == 0) then
      call messages_experimental("IonsTimeDependentDisplacements")
      ndisp= parse_block_n(blk)
      SAFE_ALLOCATE(this%td_displacements(1:geo%natoms))
      this%td_displacements(1:geo%natoms)%move = .false.
      if (ndisp > 0) this%drive_ions =.true.
      
      do i = 1, ndisp
        call parse_block_integer(blk, i-1, 0, iatom)
        this%td_displacements(iatom)%move = .true.
        
        call parse_block_string(blk, i-1, 1, expression)
        call tdf_read(this%td_displacements(iatom)%fx, namespace, trim(expression), ierr)
        if (ierr /= 0) then            
          write(message(1),'(3A)') 'Could not find "', trim(expression), '" in the TDFunctions block:'
          call messages_warning(1, namespace=namespace)
        end if
        
        
        call parse_block_string(blk, i-1, 2, expression)
        call tdf_read(this%td_displacements(iatom)%fy, namespace, trim(expression), ierr)
        if (ierr /= 0) then            
          write(message(1),'(3A)') 'Could not find "', trim(expression), '" in the TDFunctions block:'
          call messages_warning(1, namespace=namespace)
        end if
        
        call parse_block_string(blk, i-1, 3, expression)
        call tdf_read(this%td_displacements(iatom)%fz, namespace, trim(expression), ierr)
        if (ierr /= 0) then            
          write(message(1),'(3A)') 'Could not find "', trim(expression), '" in the TDFunctions block:'
          call messages_warning(1, namespace=namespace)
        end if
        
      end do
      
      SAFE_ALLOCATE(this%geo_t0)
      call geometry_copy(this%geo_t0, geo)
      
    end if
    
    


    !%Variable Thermostat
    !%Type integer
    !%Default none
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable selects the type of thermostat applied to
    !% control the ionic temperature. 
    !%Option none 0
    !% No thermostat is applied. This is the default.
    !%Option velocity_scaling 1
    !% Velocities are scaled to control the temperature.
    !%Option nose_hoover 2
    !% Nose-Hoover thermostat.
    !%End
    
    call parse_variable(namespace, 'Thermostat', THERMO_NONE, this%thermostat)
    if(.not.varinfo_valid_option('Thermostat', this%thermostat)) call messages_input_error(namespace, 'Thermostat')
    call messages_print_var_option(stdout, 'Thermostat', this%thermostat)
    
    if(this%thermostat /= THERMO_NONE) then
      
      have_velocities = .true.

      if(this%drive_ions) then
        call messages_write('You cannot use a Thermostat and IonsConstantVelocity or IonsTimeDependentDisplacements')
        call messages_write('at the same time.')
        call messages_fatal(namespace=namespace)
      end if

      call messages_experimental('Thermostat')

      !%Variable TemperatureFunction
      !%Type integer
      !%Default "temperature"
      !%Section Time-Dependent::Propagation
      !%Description
      !% If a thermostat is used, this variable indicates the name of the
      !% function in the <tt>TDFunctions</tt> block that will be used to control the
      !% temperature. The values of the temperature are given in
      !% degrees Kelvin.
      !%End
      call parse_variable(namespace, 'TemperatureFunction', 'temperature', temp_function_name)

      call tdf_read(this%temperature_function, namespace, temp_function_name, ierr)

      if(ierr /= 0) then
        message(1) = "You have enabled a thermostat but Octopus could not find"
        message(2) = "the '"//trim(temp_function_name)//"' function in the TDFunctions block."
        call messages_fatal(2, namespace=namespace)
      end if

      if(this%thermostat == THERMO_NH) then
        !%Variable ThermostatMass
        !%Type float
        !%Default 1.0
        !%Section Time-Dependent::Propagation
        !%Description
        !% This variable sets the fictitious mass for the Nose-Hoover
        !% thermostat.
        !%End
        call messages_obsolete_variable(namespace, 'NHMass', 'ThermostatMass')

        call parse_variable(namespace, 'ThermostatMass', CNST(1.0), this%nh(1)%mass)
        this%nh(2)%mass = this%nh(1)%mass

        this%nh(1:2)%pos = M_ZERO
        this%nh(1:2)%vel = M_ZERO
        
        SAFE_ALLOCATE(this%old_pos(1:geo%space%dim, 1:geo%natoms))
        
        do iatom = 1, geo%natoms
          this%old_pos(1:geo%space%dim, iatom) = geo%atom(iatom)%x(1:geo%space%dim)
        end do
      end if

    end if

    !now initialize velocities

    !%Variable RandomVelocityTemp
    !%Type float
    !%Default 0.0
    !%Section System::Velocities
    !%Description
    !% If this variable is present, <tt>Octopus</tt> will assign random
    !% velocities to the atoms following a Boltzmann distribution with
    !% temperature given by <tt>RandomVelocityTemp</tt> (in degrees Kelvin).
    !%End

    ! we now load the velocities, either from the temperature, from the input, or from a file
    if(parse_is_defined(namespace, 'RandomVelocityTemp')) then

      have_velocities = .true.

      if( mpi_grp_is_root(mpi_world)) then
        call loct_ran_init(random_gen_pointer)
        call parse_variable(namespace, 'RandomVelocityTemp', M_ZERO, temperature, unit = unit_kelvin)
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

      x = geometry_center_of_mass_vel(geo)
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
      call read_coords_read('Velocities', xyz, geo%space, namespace)
      if(xyz%source /= READ_COORDS_ERR) then
        
        have_velocities = .true.

        if(geo%natoms /= xyz%n) then
          write(message(1), '(a,i4,a,i4)') 'I need exactly ', geo%natoms, ' velocities, but I found ', xyz%n
          call messages_fatal(1, namespace=namespace)
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
    call parse_variable(namespace, 'MoveIons', have_velocities, this%move_ions)
    call messages_print_var_value(stdout, 'MoveIons', this%move_ions)

    if (this%move_ions .and. geo%space%periodic_dim == 1) then
      call messages_input_error(namespace, 'MoveIons', &
        'Moving ions for a 1D periodic system is not allowed, as forces are incorrect.')
    end if

    if(ion_dynamics_ions_move(this)) then 
      SAFE_ALLOCATE(this%oldforce(1:geo%space%dim, 1:geo%natoms))
    end if

    POP_SUB(ion_dynamics_init)
  end subroutine ion_dynamics_init


  ! ---------------------------------------------------------
  subroutine ion_dynamics_end(this)
    type(ion_dynamics_t), intent(inout) :: this

    PUSH_SUB(ion_dynamics_end)
    SAFE_DEALLOCATE_A(this%oldforce)

    if(this%thermostat /= THERMO_NONE) then
      call tdf_end(this%temperature_function)
    end if

    if (this%drive_ions .and. allocated(this%td_displacements) ) then
      if (any (this%td_displacements(1:this%geo_t0%natoms)%move)) then
        ! geometry end cannot be called here, otherwise the species are destroyed twice
        ! call geometry_end(this%geo_t0)
      end if
      SAFE_DEALLOCATE_A(this%td_displacements)
      if (any (this%td_displacements(:)%move)) then
        SAFE_DEALLOCATE_P(this%geo_t0)
      end if
    end if

    POP_SUB(ion_dynamics_end)
  end subroutine ion_dynamics_end


  ! ---------------------------------------------------------
  subroutine ion_dynamics_propagate(this, sb, geo, time, dt, namespace)
    type(ion_dynamics_t), intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(geometry_t),     intent(inout) :: geo
    FLOAT,                intent(in)    :: time
    FLOAT,                intent(in)    :: dt
    type(namespace_t),    intent(in)    :: namespace

    integer :: iatom
    FLOAT   :: DR(1:3)

    if(.not. ion_dynamics_ions_move(this)) return

    PUSH_SUB(ion_dynamics_propagate)
    
    DR = M_ZERO

    this%dt = dt
    

    ! get the temperature from the tdfunction for the current time
    if(this%thermostat /= THERMO_NONE) then
      this%current_temperature = units_to_atomic(unit_kelvin, tdf(this%temperature_function, time))

      if(this%current_temperature < M_ZERO) then 
        write(message(1), '(a, f10.3, 3a, f10.3, 3a)') &
          "Negative temperature (", &
          units_from_atomic(unit_kelvin, this%current_temperature), " ", units_abbrev(unit_kelvin), &
          ") at time ", &
          units_from_atomic(units_out%time, time), " ", trim(units_abbrev(units_out%time)), "."
        call messages_fatal(1, namespace=namespace)
      end if
    else
      this%current_temperature = CNST(0.0)
    end if

    if(this%thermostat /= THERMO_NH) then
      ! integrate using verlet
      do iatom = 1, geo%natoms
        if(.not. geo%atom(iatom)%move) cycle

        if(.not. this%drive_ions) then

          geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) &
            + dt*geo%atom(iatom)%v(1:geo%space%dim) + &
            M_HALF*dt**2 / species_mass(geo%atom(iatom)%species) * geo%atom(iatom)%f(1:geo%space%dim)
          
          this%oldforce(1:geo%space%dim, iatom) = geo%atom(iatom)%f(1:geo%space%dim)
          
        else
          if(this%constant_velocity) then
            geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) &
                                                + dt*geo%atom(iatom)%v(1:geo%space%dim)
          end if


          if (this%td_displacements(iatom)%move) then
            
            DR(1:3)=(/TOFLOAT(tdf(this%td_displacements(iatom)%fx,time)), &
                      TOFLOAT(tdf(this%td_displacements(iatom)%fy,time)), &
                      TOFLOAT(tdf(this%td_displacements(iatom)%fz,time)) /)

            geo%atom(iatom)%x(1:geo%space%dim) = this%geo_t0%atom(iatom)%x(1:geo%space%dim) + DR(1:geo%space%dim)
          end if
            
        end if

      end do

    else
      ! for the Nose-Hoover thermostat we use a special integrator

      ! The implementation of the Nose-Hoover thermostat is based on
      ! Understanding Molecular Simulations by Frenkel and Smit,
      ! Appendix E, page 540-542.

      call nh_chain(this, geo)

      do iatom = 1, geo%natoms
        geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) + M_HALF*dt*geo%atom(iatom)%v(1:geo%space%dim)
      end do

    end if

    ! When the system is periodic in some directions, the atoms might have moved to a an adjacent cell, so we need to move them back to the original cell
    call geometry_fold_atoms_into_cell(geo)

    POP_SUB(ion_dynamics_propagate)
  end subroutine ion_dynamics_propagate
  

  ! ---------------------------------------------------------
  subroutine nh_chain(this, geo)
    type(ion_dynamics_t), intent(inout) :: this
    type(geometry_t),     intent(inout) :: geo

    FLOAT :: g1, g2, ss, uk, dt, temp
    integer :: iatom

    PUSH_SUB(nh_chain)

    dt = this%dt

    uk = ion_dynamics_kinetic_energy(geo)

    temp = this%current_temperature
    
    g2 = (this%nh(1)%mass*this%nh(1)%vel**2 - temp)/this%nh(2)%mass
    this%nh(2)%vel = this%nh(2)%vel + g2*dt/CNST(4.0)
    this%nh(1)%vel = this%nh(1)%vel*exp(-this%nh(2)%vel*dt/CNST(8.0))

    g1 = (CNST(2.0)*uk - M_THREE*geo%natoms*temp)/this%nh(1)%mass
    this%nh(1)%vel = this%nh(1)%vel + g1*dt/CNST(4.0)
    this%nh(1)%vel = this%nh(1)%vel*exp(-this%nh(2)%vel*dt/CNST(8.0))
    this%nh(1)%pos = this%nh(1)%pos + this%nh(1)%vel*dt/CNST(2.0)
    this%nh(2)%pos = this%nh(2)%pos + this%nh(2)%vel*dt/CNST(2.0)

    ss = exp(-this%nh(1)%vel*dt/CNST(2.0))
    
    do iatom = 1, geo%natoms
      geo%atom(iatom)%v(1:geo%space%dim) = ss*geo%atom(iatom)%v(1:geo%space%dim)
    end do
    
    uk = uk*ss**2

    this%nh(1)%vel = this%nh(1)%vel*exp(-this%nh(2)%vel*dt/CNST(8.0))
    g1 = (CNST(2.0)*uk - M_THREE*geo%natoms*temp)/this%nh(1)%mass
    this%nh(1)%vel = this%nh(1)%vel + g1*dt/CNST(4.0)
    this%nh(1)%vel = this%nh(1)%vel*exp(-this%nh(2)%vel*dt/CNST(8.0))

    g2 = (this%nh(1)%mass*this%nh(1)%vel**2 - temp)/this%nh(2)%mass
    this%nh(2)%vel = this%nh(2)%vel + g2*dt/CNST(4.0)
    
    POP_SUB(nh_chain)
  end subroutine nh_chain
  

  ! ---------------------------------------------------------
  subroutine ion_dynamics_propagate_vel(this, geo, atoms_moved)
    type(ion_dynamics_t), intent(inout) :: this
    type(geometry_t),     intent(inout) :: geo
    logical, optional,    intent(out)   :: atoms_moved !< Returns true if the atoms were moved by this function.

    integer :: iatom
    FLOAT   :: scal

    if(.not. ion_dynamics_ions_move(this)) return
    if(this%drive_ions) return

    PUSH_SUB(ion_dynamics_propagate_vel)
    
    if(present(atoms_moved)) atoms_moved = this%thermostat == THERMO_NH

    if(this%thermostat /= THERMO_NH) then
      ! velocity verlet
      
      do iatom = 1, geo%natoms
        if(.not. geo%atom(iatom)%move) cycle
        
        geo%atom(iatom)%v(1:geo%space%dim) = geo%atom(iatom)%v(1:geo%space%dim) &
          + this%dt/species_mass(geo%atom(iatom)%species) * M_HALF * (this%oldforce(1:geo%space%dim, iatom) + &
          geo%atom(iatom)%f(1:geo%space%dim))
        
      end do
      
    else
      ! the nose-hoover integration
      do iatom = 1, geo%natoms
        geo%atom(iatom)%v(1:geo%space%dim) = geo%atom(iatom)%v(1:geo%space%dim) + &
          this%dt*geo%atom(iatom)%f(1:geo%space%dim) / species_mass(geo%atom(iatom)%species)
        geo%atom(iatom)%x(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim) + M_HALF*this%dt*geo%atom(iatom)%v(1:geo%space%dim)
      end do
      
      call nh_chain(this, geo)

    end if

    if(this%thermostat == THERMO_SCAL) then
      scal = sqrt(this%current_temperature/ion_dynamics_temperature(geo))

      do iatom = 1, geo%natoms
        geo%atom(iatom)%v(1:geo%space%dim) = scal*geo%atom(iatom)%v(1:geo%space%dim)
      end do
    end if

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

    if(this%thermostat == THERMO_NH) then
      SAFE_ALLOCATE(state%old_pos(1:geo%space%dim, 1:geo%natoms))
      state%old_pos(1:geo%space%dim, 1:geo%natoms) = this%old_pos(1:geo%space%dim, 1:geo%natoms)
      state%nh(1:2)%pos = this%nh(1:2)%pos
      state%nh(1:2)%vel = this%nh(1:2)%vel
    end if
    
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

    if(this%thermostat == THERMO_NH) then
      this%old_pos(1:geo%space%dim, 1:geo%natoms) = state%old_pos(1:geo%space%dim, 1:geo%natoms)
      this%nh(1:2)%pos = state%nh(1:2)%pos
      this%nh(1:2)%vel = state%nh(1:2)%vel
      SAFE_DEALLOCATE_A(state%old_pos)
    end if

    SAFE_DEALLOCATE_A(state%pos)
    SAFE_DEALLOCATE_A(state%vel)
    
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
