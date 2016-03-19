#include "global.h"

module ssys_handle_oct_m

  use base_handle_oct_m
  use frozen_handle_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use live_handle_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m

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
  subroutine ssys_handle_init(this, geo, config)
    type(base_handle_t), intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(json_object_t), intent(in)  :: config

    type(json_array_iterator_t)  :: iter
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_init)

    nullify(cnfg, list, hndl)
    call base_handle__init__(this, config)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_SSYS)
    call json_get(config, "systems", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
      nullify(cnfg, hndl)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call json_get(cnfg, "type", type, ierr)
      ASSERT(ierr==JSON_OK)
      call base_handle_new(this, hndl)
      select case(type)
      case(HNDL_TYPE_FRZN)
        call frozen_handle_init(hndl, cnfg)
      case(HNDL_TYPE_LIVE)
        call live_handle_init(hndl, geo, cnfg)
      case default
        message(1) = "Unknown subsystems type."
        call messages_fatal(1)
      end select
      call base_handle_sets(this, hndl, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list, hndl)
    call base_handle__init__(this)

    POP_SUB(ssys_handle_init)
  end subroutine ssys_handle_init

  ! ---------------------------------------------------------
  subroutine ssys_handle_start(this, grid)
    type(base_handle_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_start)

    nullify(hndl)
    call base_handle__start__(this, grid)
    call base_handle_init(iter, this)
    do
      nullify(hndl)
      call base_handle_next(iter, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_get(hndl, type)
      select case(type)
      case(HNDL_TYPE_FRZN)
        call frozen_handle_start(hndl, grid, mpi_world)
      case(HNDL_TYPE_LIVE)
        call live_handle_start(hndl, grid)
      case default
        message(1)="Unknown subsystems type."
        call messages_fatal(1)
      end select
    end do
    call base_handle_end(iter)
    nullify(hndl)

    POP_SUB(ssys_handle_start)
  end subroutine ssys_handle_start

  ! ---------------------------------------------------------
  subroutine ssys_handle_update(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_update)

    nullify(hndl)
    call base_handle_init(iter, this)
    do
      nullify(hndl)
      call base_handle_next(iter, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_get(hndl, type)
      select case(type)
      case(HNDL_TYPE_FRZN)
      case(HNDL_TYPE_LIVE)
      case default
        message(1)="Unknown subsystems type."
        call messages_fatal(1)
      end select
    end do
    call base_handle_end(iter)
    nullify(hndl)
    call base_handle__update__(this)

    POP_SUB(ssys_handle_update)
  end subroutine ssys_handle_update

  ! ---------------------------------------------------------
  subroutine ssys_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(ssys_handle_stop)

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
  subroutine ssys_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_end)

    call base_handle_init(iter, this)
    do
      nullify(hndl)
      call base_handle_next(iter, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call base_handle_get(hndl, type)
      select case(type)
      case(HNDL_TYPE_FRZN)
        call frozen_handle_end(hndl)
      case(HNDL_TYPE_LIVE)
        call live_handle_end(hndl)
      case default
        message(1)="Unknown subsystems type."
        call messages_fatal(1)
      end select
    end do
    call base_handle_end(iter)
    nullify(hndl)
    call base_handle_end(this)

    POP_SUB(ssys_handle_end)
  end subroutine ssys_handle_end

end module ssys_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
