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
  use global_m
  use io_m
  use datasets_m
  use loct_math_m
  use loct_parser_m
  use units_m
  use messages_m
  use geometry_m
  use loct_m
  use profiling_m
  use varinfo_m
  use math_m

  implicit none

  private

  public ::                                &
    ion_dynamics_t,                        &
    ion_state_t,                           &
    ion_dynamics_init,                     &
    ion_dynamics_end,                      &
    ion_dynamics_propagate,                &
    ion_dynamics_propagate_velocities,     &
    ion_dynamics_save_state,               &
    ion_dynamics_restore_state,            &
    ion_dynamics_ions_move

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
    type(geometry_t),     intent(in)    :: geo

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

  end subroutine ion_dynamics_init

  subroutine ion_dynamics_end(this)
    type(ion_dynamics_t), intent(inout) :: this

    if(associated(this%oldforce)) deallocate(this%oldforce)

  end subroutine ion_dynamics_end

  subroutine ion_dynamics_propagate(this, geo, dt)
    type(ion_dynamics_t), intent(inout) :: this
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

  end subroutine ion_dynamics_propagate

  subroutine ion_dynamics_propagate_velocities(this, geo)
    type(ion_dynamics_t), intent(in)    :: this
    type(geometry_t),     intent(inout) :: geo

    integer :: iatom

    if(.not. ion_dynamics_ions_move(this)) return

    do iatom = 1, geo%natoms
      if(.not. geo%atom(iatom)%move) cycle

      geo%atom(iatom)%v(1:MAX_DIM) = geo%atom(iatom)%v(1:MAX_DIM) &
        + this%dt/geo%atom(iatom)%spec%weight*M_HALF*(this%oldforce(1:MAX_DIM, iatom) + geo%atom(iatom)%f(1:MAX_DIM))

    end do

  end subroutine ion_dynamics_propagate_velocities

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
  
end module ion_dynamics_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
