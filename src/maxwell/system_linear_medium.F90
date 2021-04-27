!! Copyright (C) 2020 F. Bonafé, H. Appel
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

module system_linear_medium_oct_m
  use algorithm_oct_m
  use clock_oct_m
#ifdef HAVE_DFTBPLUS
  use dftbplus
#endif
  use geometry_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use iso_c_binding
  use lasers_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use species_oct_m
  use system_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    system_linear_medium_t,    &
    system_linear_medium_init

   type, extends(system_t) :: system_linear_medium_t
  contains
    procedure :: init_interaction => system_linear_medium_init_interaction
    procedure :: initial_conditions => system_linear_medium_initial_conditions
    procedure :: do_td_operation => system_linear_medium_do_td
    procedure :: iteration_info => system_linear_medium_iteration_info
    procedure :: output_start => system_linear_medium_output_start
    procedure :: output_write => system_linear_medium_output_write
    procedure :: output_finish => system_linear_medium_output_finish
    procedure :: is_tolerance_reached => system_linear_medium_is_tolerance_reached
    procedure :: update_quantity => system_linear_medium_update_quantity
    procedure :: update_exposed_quantity => system_linear_medium_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => system_linear_medium_copy_quantities_to_interaction
    procedure :: update_interactions_start => system_linear_medium_update_interactions_start
    procedure :: update_interactions_finish => system_linear_medium_update_interactions_finish
    final :: system_linear_medium_finalize
  end type system_linear_medium_t

  interface system_linear_medium_t
    procedure system_linear_medium_constructor
  end interface system_linear_medium_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function system_linear_medium_constructor(namespace) result(sys)
    class(system_linear_medium_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(system_linear_medium_constructor)

    SAFE_ALLOCATE(sys)

    call system_linear_medium_init(sys, namespace)

    POP_SUB(system_linear_medium_constructor)
  end function system_linear_medium_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine system_linear_medium_init(this, namespace)
    class(system_linear_medium_t), target, intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace


    PUSH_SUB(system_linear_medium_init)

    this%namespace = namespace

    POP_SUB(system_linear_medium_init)
  end subroutine system_linear_medium_init

  ! ---------------------------------------------------------
  subroutine system_linear_medium_init_interaction(this, interaction)
    class(system_linear_medium_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(system_linear_medium_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_init_interaction)
  end subroutine system_linear_medium_init_interaction

  ! ---------------------------------------------------------
  subroutine system_linear_medium_initial_conditions(this, from_scratch)
    class(system_linear_medium_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    PUSH_SUB(system_linear_medium_initial_conditions)

    POP_SUB(system_linear_medium_initial_conditions)
  end subroutine system_linear_medium_initial_conditions

  ! ---------------------------------------------------------
  subroutine system_linear_medium_do_td(this, operation)
    class(system_linear_medium_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(system_linear_medium_do_td)

   POP_SUB(system_linear_medium_do_td)
  end subroutine system_linear_medium_do_td

  ! ---------------------------------------------------------
  logical function system_linear_medium_is_tolerance_reached(this, tol) result(converged)
    class(system_linear_medium_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(system_linear_medium_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(system_linear_medium_is_tolerance_reached)
  end function system_linear_medium_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine system_linear_medium_iteration_info(this)
    class(system_linear_medium_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(system_linear_medium_iteration_info)

    write(message(1),'(2X,A,1X,A)') "Linar medium system:", trim(this%namespace%get())

    write(message(5),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(6),'(4x,A,e14.6)') 'Simulation time: ', this%clock%time()
    call messages_info(6)

    POP_SUB(system_linear_medium_iteration_info)
  end subroutine system_linear_medium_iteration_info

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_start(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_output_start)

    ! Output info for first iteration
    call this%output_write()

    POP_SUB(system_linear_medium_output_start)
  end subroutine system_linear_medium_output_start

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_finish(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_output_finish)

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle(1))
      call write_iter_end(this%output_handle(2))
    end if

    POP_SUB(system_linear_medium_output_finish)
  end subroutine system_linear_medium_output_finish

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_write(this)
    class(system_linear_medium_t), intent(inout) :: this

    integer :: idir, iat, iout
    character(len=50) :: aux
    character(1) :: out_label(2)
    FLOAT :: tmp(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(system_linear_medium_output_write)

    POP_SUB(system_linear_medium_output_write)
  end subroutine system_linear_medium_output_write

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_quantity(this, iq)
    class(system_linear_medium_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_linear_medium_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_update_quantity)
  end subroutine system_linear_medium_update_quantity

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_exposed_quantity(partner, iq)
    class(system_linear_medium_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_linear_medium_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_update_exposed_quantity)
  end subroutine system_linear_medium_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine system_linear_medium_copy_quantities_to_interaction(partner, interaction)
    class(system_linear_medium_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(system_linear_medium_copy_quantities_to_interaction)

    select type (interaction)
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_copy_quantities_to_interaction)
  end subroutine system_linear_medium_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_interactions_start(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_update_interactions_start)

    POP_SUB(system_linear_medium_update_interactions_start)
  end subroutine system_linear_medium_update_interactions_start

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_interactions_finish(this)
    class(system_linear_medium_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_linear_medium_update_interactions_finish)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class default
        message(1) = "Interactions not implemented for linear medium systems."
        call messages_fatal(1)
      end select
    end do

    POP_SUB(system_linear_medium_update_interactions_finish)
  end subroutine system_linear_medium_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine system_linear_medium_finalize(this)
    type(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_finalize)

    call system_end(this)

    POP_SUB(system_linear_medium_finalize)
  end subroutine system_linear_medium_finalize

end module system_linear_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
