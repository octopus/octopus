#include "global.h"

module ssys_handle_oct_m

  use refcount_oct_m

  use base_functional_oct_m
  use base_hamiltonian_oct_m
  use base_handle_oct_m
  use frozen_handle_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use live_handle_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

  implicit none

  private

  public ::             &
    ssys_handle_init,   &
    ssys_handle_start,  &
    ssys_handle_update, &
    ssys_handle_stop,   &
    ssys_handle_copy,   &
    ssys_handle_end

  integer, public, parameter :: HNDL_TYPE_SSYS = 9

contains
  
  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian__init__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_hamiltonian_t), pointer :: tnadd, hmlt
    type(base_functional_t),  pointer :: fnct

    PUSH_SUB(ssys_hamiltonian__init__)

    nullify(tnadd, hmlt, fnct)
    call base_hamiltonian_get(this, "tnadd", tnadd)
    ASSERT(associated(tnadd))
    call base_hamiltonian_get(this, "kinetic", fnct)
    ASSERT(associated(fnct))
    call base_hamiltonian_set(tnadd, "total", fnct)
    nullify(fnct)
    call base_hamiltonian_gets(this, "live", hmlt)
    ASSERT(associated(hmlt))
    call base_hamiltonian_get(hmlt, "kinetic", fnct)
    ASSERT(associated(fnct))
    call base_hamiltonian_set(tnadd, "live", fnct)
    nullify(tnadd, hmlt, fnct)

    POP_SUB(ssys_hamiltonian__init__)
  end subroutine ssys_hamiltonian__init__

  ! ---------------------------------------------------------
  subroutine ssys_handle__init__(this)
    type(base_handle_t), intent(inout) :: this

    type(base_hamiltonian_t), pointer :: hmlt

    PUSH_SUB(ssys_handle__init__)

    nullify(hmlt)
    call base_handle_get(this, hmlt)
    ASSERT(associated(hmlt))
    call ssys_hamiltonian__init__(hmlt)

    POP_SUB(ssys_handle__init__)
  end subroutine ssys_handle__init__

  ! ---------------------------------------------------------
  subroutine ssys_handle_init(this, geo, config)
    type(base_handle_t), intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(ssys_handle_init)

    call base_handle__init__(this, config)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_SSYS)
    call base_handle__init__(this, init)
    call ssys_handle__init__(this)

    POP_SUB(ssys_handle_init)
    
  contains

    subroutine init(this, config)
      type(base_handle_t), intent(out) :: this
      type(json_object_t), intent(in)  :: config

      integer :: type, ierr

      PUSH_SUB(ssys_handle_init.init)

      call json_get(config, "type", type, ierr)
      ASSERT(ierr==JSON_OK)
      select case(type)
      case(HNDL_TYPE_FRZN)
        call frozen_handle_init(this, config)
      case(HNDL_TYPE_LIVE)
        call live_handle_init(this, geo, config)
      case default
        ASSERT(.false.)
      end select

      POP_SUB(ssys_handle_init.init)
    end subroutine init
    
  end subroutine ssys_handle_init

  ! ---------------------------------------------------------
  subroutine ssys_handle__start__(this, sim)
    type(base_handle_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim

    integer :: type

    PUSH_SUB(ssys_handle__start__)

    call base_handle_get(this, type)
    select case(type)
    case(HNDL_TYPE_FRZN)
      call frozen_handle_start(this, sim)
    case(HNDL_TYPE_LIVE)
      call live_handle_start(this, sim)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(ssys_handle__start__)
  end subroutine ssys_handle__start__

  ! ---------------------------------------------------------
  subroutine ssys_handle_start(this, grid)
    type(base_handle_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid

    PUSH_SUB(ssys_handle_start)

    call base_handle__start__(this, simstr)
    call base_handle_start(this, ssys_handle__start__)
    call base_handle_sets(this, "frozen", lock=.true.)

    POP_SUB(ssys_handle_start)
    
  contains

    subroutine simstr(this)
      type(simulation_t), intent(inout) :: this

      PUSH_SUB(ssys_handle_start.simstr)

      call simulation_start(this, grid)

      POP_SUB(ssys_handle_start.simstr)
    end subroutine simstr

  end subroutine ssys_handle_start

  ! ---------------------------------------------------------
  subroutine ssys_handle__update__(this)
    type(base_handle_t), intent(inout) :: this

    integer :: type

    PUSH_SUB(ssys_handle__update__)

    call base_handle_get(this, type=type)
    select case(type)
    case(HNDL_TYPE_FRZN)
    case(HNDL_TYPE_LIVE)
    case(HNDL_TYPE_SSYS)
      call base_handle__update__(this)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(ssys_handle__update__)
  end subroutine ssys_handle__update__

  ! ---------------------------------------------------------
  subroutine ssys_handle_update(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(ssys_handle_update)

    call base_handle__apply__(this, ssys_handle__update__)

    POP_SUB(ssys_handle_update)
  end subroutine ssys_handle_update

  ! ---------------------------------------------------------
  subroutine ssys_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(ssys_handle_stop)

    call base_handle_sets(this, "frozen", lock=.false.)
    call base_handle_stop(this)

    POP_SUB(ssys_handle_stop)
  end subroutine ssys_handle_stop

  ! ---------------------------------------------------------
  subroutine ssys_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(ssys_handle_copy)

    call base_handle_copy(this, that)

    POP_SUB(ssys_handle_copy)
  end subroutine ssys_handle_copy

  ! ---------------------------------------------------------
  subroutine ssys_handle__end__(this)
    type(base_handle_t), intent(inout) :: this

    integer :: type

    PUSH_SUB(ssys_handle__end__)

    call base_handle_get(this, type=type)
    select case(type)
    case(HNDL_TYPE_FRZN)
      call frozen_handle_end(this)
    case(HNDL_TYPE_LIVE)
      call live_handle_end(this)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(ssys_handle__end__)
  end subroutine ssys_handle__end__

  ! ---------------------------------------------------------
  subroutine ssys_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    integer :: type

    PUSH_SUB(ssys_handle_end)

    call base_handle_get(this, type=type)
    ASSERT(type==HNDL_TYPE_SSYS)
    call base_handle_end(this, ssys_handle__end__)

    POP_SUB(ssys_handle_end)
  end subroutine ssys_handle_end

end module ssys_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
