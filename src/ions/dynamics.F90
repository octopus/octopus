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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: dynamics.F90 3813 2008-03-04 20:40:01Z xavier $

#include "global.h"

module ion_dynamics_m
  use c_pointer_m
  use global_m
  use io_m
  use datasets_m
  use loct_math_m
  use loct_parser_m
  use units_m
  use messages_m
  use mpi_m
  use geometry_m
  use loct_m
  use simul_box_m
  use profiling_m
  use varinfo_m
  use math_m
  use xyz_file_m

  implicit none

  private

  public ::                                &
    ion_dynamics_t,                        &
    ion_state_t,                           &
    ion_dynamics_init,                     &
    ion_dynamics_end,                      &
    ion_dynamics_propagate,                &
    ion_dynamics_propagate_vel,     &
    ion_dynamics_save_state,               &
    ion_dynamics_restore_state,            &
    ion_dynamics_ions_move,                &
    ion_dynamics_kinetic_energy

  integer, parameter ::   &
    STATIC_IONS     = 0,  &
    VELOCITY_VERLET = 1

  type ion_dynamics_t
    private
    integer          :: method
    FLOAT, pointer   :: oldforce(:, :)
    FLOAT            :: dt
  end type ion_dynamics_t

  type ion_state_t
    private
    FLOAT, pointer :: pos(:, :)
    FLOAT, pointer :: vel(:, :)
  end type ion_state_t

contains

  subroutine ion_dynamics_init(this, geo)
    type(ion_dynamics_t), intent(out)   :: this
    type(geometry_t),     intent(inout) :: geo

    integer :: i, j
    FLOAT   :: x(MAX_DIM), temperature, sigma, kin1, kin2
    type(c_pointer_t) :: random_gen_pointer
    type(xyz_file_info) :: xyz

    !%Variable MoveIons
    !%Type integer
    !%Default static_ions
    !%Section Time Dependent::Propagation
    !%Description
    !% This variable specifies how to treat the dynamic of the ions
    !% during a time-dependent run. By default they will remain static.
    !%Option static_ions 0
    !% Do not move the ions.
    !%Option vel_verlet 1
    !% Newtonian dynamics using the velocity Verlet integrator.
    !%End
    
    call loct_parse_int(check_inp('MoveIons'), STATIC_IONS, this%method)
    if(.not.varinfo_valid_option('MoveIons', this%method)) call input_error('MoveIons')
    call messages_print_var_option(stdout, 'MoveIons', this%method)

    nullify(this%oldforce)
    
    if(ion_dynamics_ions_move(this)) then 

      ALLOCATE(this%oldforce(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)

    end if

    !now initialize velocities

    !%Variable RandomVelocityTemp
    !%Type float
    !%Section System::Velocities
    !%Description
    !% If this variable is present, octopus will assign random
    !% velocities to the atoms following a Bolzmann distribution with
    !% temperature given by RandomVelocityTemp (in degrees Kelvin).
    !%End

    ! we now load the velocities, either from the temperature, from the input, or from a file
    if(loct_parse_isdef(check_inp('RandomVelocityTemp')).ne.0) then

      if( mpi_grp_is_root(mpi_world)) then
        call loct_ran_init(random_gen_pointer)
        call loct_parse_float(check_inp('RandomVelocityTemp'), M_ZERO, temperature)
      end if

      do i = 1, geo%natoms
        !generate the velocities in the root node
        if( mpi_grp_is_root(mpi_world)) then
          sigma = sqrt( P_Kb*temperature / geo%atom(i)%spec%weight )
          do j = 1, 3
             geo%atom(i)%v(j) = loct_ran_gaussian(random_gen_pointer, sigma)
          end do
        end if
#ifdef HAVE_MPI
        !and send them to the others
        call MPI_Bcast(geo%atom(i)%v, MAX_DIM, MPI_FLOAT, 0, mpi_world%comm, mpi_err)
#endif
      end do

      if( mpi_grp_is_root(mpi_world)) then
        call loct_ran_end(random_gen_pointer)
      end if

      kin1 = ion_dynamics_kinetic_energy(this, geo)

      call cm_vel(geo, x)
      do i = 1, geo%natoms
        geo%atom(i)%v = geo%atom(i)%v - x
      end do

      kin2 = ion_dynamics_kinetic_energy(this, geo)

      do i = 1, geo%natoms
        geo%atom(i)%v(:) =  sqrt(kin1/kin2)*geo%atom(i)%v(:)
      end do

      write(message(1),'(a,f10.4,1x,a)') 'Info: Initial velocities randomly distributed with T =', &
        temperature, 'K'
      write(message(2),'(2x,a,f8.4,1x,a)') '<K>       =', &
        (ion_dynamics_kinetic_energy(this, geo)/geo%natoms)/units_out%energy%factor, &
        units_out%energy%abbrev
      write(message(3),'(2x,a,f8.4,1x,a)') '3/2 k_B T =', &
        (M_THREE/M_TWO)*P_Kb*temperature/units_out%energy%factor, &
        units_out%energy%abbrev
      call write_info(3)

    else
      !%Variable XYZVelocities
      !%Type string
      !%Section System::Velocities
      !%Description
      !% octopus will try to read the starting velocities of the atoms from the XYZ file 
      !% specified by the variable XYZVelocities.
      !% Note that you do not need to specify initial velocities if you are not going
      !% to perform ion dynamics; if you are going to allow the ions to move but the velocities
      !% are not specified, they are considered to be null.
      !%End

      !%Variable Velocities
      !%Type block
      !%Section System::Velocities
      !%Description
      !% If XYZVelocities is not present, octopus will try to fetch the initial 
      !% atomic velocities from this block. If this block is not present, octopus
      !% will reset the initial velocities to zero. The format of this block can be
      !% illustrated by this example:
      !%
      !% <tt>%Velocities
      !% <br>&nbsp;&nbsp;'C'  |      -1.7 | 0.0 | 0.0
      !% <br>&nbsp;&nbsp;'O'  | &nbsp;1.7 | 0.0 | 0.0
      !% <br>%</tt>
      !%
      !% It describes one Carbon and one Oxygen moving at the relative
      !% velocity of 3.4, velocity units.
      !%
      !% Note: It is important for the velocities to maintain the ordering 
      !% in which the species were defined in the coordinates specifications.
      !%End

      call xyz_file_init(xyz)
      call xyz_file_read('Velocities', xyz)
      if(xyz%file_type.ne.XYZ_FILE_ERR) then
        if(geo%natoms.ne.xyz%n) then
          write(message(1), '(a,i4,a,i4)') 'I need exactly ', geo%natoms, ' velocities, but I found ', xyz%n
          call write_fatal(1)
        end if

        ! copy information and adjust units
        do i = 1, geo%natoms
          geo%atom(i)%v = xyz%atom(i)%x * (units_inp%velocity%factor / units_inp%length%factor)
        end do
        call xyz_file_end(xyz)

      else
        do i = 1, geo%natoms
          geo%atom(i)%v = M_ZERO
        end do
      end if
    end if

    geo%kinetic_energy = ion_dynamics_kinetic_energy(this, geo)

  end subroutine ion_dynamics_init

  subroutine ion_dynamics_end(this)
    type(ion_dynamics_t), intent(inout) :: this

    if(associated(this%oldforce)) deallocate(this%oldforce)

  end subroutine ion_dynamics_end

  subroutine ion_dynamics_propagate(this, sb, geo, dt)
    type(ion_dynamics_t), intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(geometry_t),     intent(inout) :: geo
    FLOAT,                intent(in)    :: dt
    
    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    this%dt = dt

    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle

      geo%atom(iatom)%x(1:MAX_DIM) = geo%atom(iatom)%x(1:MAX_DIM) &
        + dt*geo%atom(iatom)%v(1:MAX_DIM) + M_HALF*dt**2/geo%atom(iatom)%spec%weight*geo%atom(iatom)%f(1:MAX_DIM)

      this%oldforce(1:MAX_DIM, iatom) = geo%atom(iatom)%f(1:MAX_DIM)

    end do

    call simul_box_atoms_in_box(sb, geo)
    
  end subroutine ion_dynamics_propagate

  subroutine ion_dynamics_propagate_vel(this, geo)
    type(ion_dynamics_t), intent(in)    :: this
    type(geometry_t),     intent(inout) :: geo

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle

      geo%atom(iatom)%v(1:MAX_DIM) = geo%atom(iatom)%v(1:MAX_DIM) &
        + this%dt/geo%atom(iatom)%spec%weight*M_HALF*(this%oldforce(1:MAX_DIM, iatom) + geo%atom(iatom)%f(1:MAX_DIM))

    end do

  end subroutine ion_dynamics_propagate_vel

  subroutine ion_dynamics_save_state(this, geo, state)
    type(ion_dynamics_t), intent(in)    :: this
    type(geometry_t),     intent(in)    :: geo
    type(ion_state_t),    intent(out)   :: state

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    ALLOCATE(state%pos(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
    ALLOCATE(state%vel(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)

    do iatom = 1, geo%natoms
      state%pos(1:MAX_DIM, iatom) = geo%atom(iatom)%x(1:MAX_DIM)
      state%vel(1:MAX_DIM, iatom) = geo%atom(iatom)%v(1:MAX_DIM)
    end do

  end subroutine ion_dynamics_save_state

  subroutine ion_dynamics_restore_state(this, geo, state)
    type(ion_dynamics_t), intent(in)    :: this
    type(geometry_t),     intent(inout) :: geo
    type(ion_state_t),    intent(inout) :: state

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    do iatom = 1, geo%natoms
      geo%atom(iatom)%x(1:MAX_DIM) = state%pos(1:MAX_DIM, iatom)
      geo%atom(iatom)%v(1:MAX_DIM) = state%vel(1:MAX_DIM, iatom)
    end do

    deallocate(state%pos, state%vel)
    
  end subroutine ion_dynamics_restore_state

  logical pure function ion_dynamics_ions_move(this) result(ions_move)
    type(ion_dynamics_t), intent(in)    :: this
    
    ions_move = (this%method /= STATIC_IONS)
    
  end function ion_dynamics_ions_move
  
  ! ---------------------------------------------------------
  FLOAT pure function ion_dynamics_kinetic_energy(this, geo) result(kinetic_energy)
    type(ion_dynamics_t),  intent(in) :: this
    type(geometry_t),      intent(in) :: geo

    integer :: iatom

    kinetic_energy = M_ZERO
    do iatom = 1, geo%natoms
      kinetic_energy = kinetic_energy + M_HALF*geo%atom(iatom)%spec%weight*sum(geo%atom(iatom)%v(1:MAX_DIM)**2)
    end do

  end function ion_dynamics_kinetic_energy

end module ion_dynamics_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
