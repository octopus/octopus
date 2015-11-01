#include "global.h"

module ssys_handle_m

  use global_m
  use messages_m
  use mpi_m
  use profiling_m

  use geometry_m
  use grid_m, only: grid_t
  use json_m, only: JSON_OK, json_object_t, json_array_t, json_get
  use json_m, only: json_array_iterator_t, json_init, json_next, json_end

  use frozen_handle_m, only: &
    frozen_handle_init,      &
    frozen_handle_start

  use frozen_handle_m, only: &
    HNDL_TYPE_FRZN

  use live_handle_m, only: &
    live_handle_init,      &
    live_handle_start

  use live_handle_m, only: &
    HNDL_TYPE_LIVE

  use base_handle_m, only: &
    base_handle__init__,   &
    base_handle__start__,  &
    base_handle__update__

  use base_handle_m, only: &
    base_handle_t

  use base_handle_m, only: &
    base_handle_init,      &
    base_handle_get

  use base_handle_m, only:          &
    ssys_handle_t => base_handle_t

  use base_handle_m, only:                &
    ssys_handle_stop => base_handle_stop, &
    ssys_handle_next => base_handle_next, &
    ssys_handle_get  => base_handle_get,  &
    ssys_handle_copy => base_handle_copy, &
    ssys_handle_end  => base_handle_end

  use base_handle_m, only:                            &
    ssys_handle_iterator_t => base_handle_iterator_t

  use base_handle_m, only:                              &
    SSYS_HANDLE_OK          => BASE_HANDLE_OK,          &
    SSYS_HANDLE_KEY_ERROR   => BASE_HANDLE_KEY_ERROR,   &
    SSYS_HANDLE_EMPTY_ERROR => BASE_HANDLE_EMPTY_ERROR

  implicit none

  private

  public ::        &
    ssys_handle_t

  public ::             &
    ssys_handle_init,   &
    ssys_handle_start,  &
    ssys_handle_update, &
    ssys_handle_stop,   &
    ssys_handle_next,   &
    ssys_handle_get,    &
    ssys_handle_copy,   &
    ssys_handle_end

  public ::                 &
    ssys_handle_iterator_t

  public ::                  &
    SSYS_HANDLE_OK,          &
    SSYS_HANDLE_KEY_ERROR,   &
    SSYS_HANDLE_EMPTY_ERROR

  interface ssys_handle_init
    module procedure ssys_handle_init_handle
    module procedure ssys_handle_iterator_init
  end interface ssys_handle_init

  integer, public, parameter :: HNDL_TYPE_SSYS = 9

contains
  
  ! ---------------------------------------------------------
  subroutine ssys_handle_init_handle(this, geo, config)
    type(ssys_handle_t), intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(json_object_t), intent(in)  :: config

    type(json_array_iterator_t)  :: iter
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_init_handle)

    nullify(cnfg, list, hndl)
    call base_handle__init__(this, config)
    call ssys_handle_get(this, type)
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
      SAFE_ALLOCATE(hndl)
      select case(type)
      case(HNDL_TYPE_FRZN)
        call frozen_handle_init(hndl, cnfg)
      case(HNDL_TYPE_LIVE)
        call live_handle_init(hndl, geo, cnfg)
      case default
        message(1)="Unknown subsystems type."
        call messages_fatal(1)
      end select
      !call base_handle_sets(this, hndl, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list, hndl)
    call base_handle__init__(this)

    POP_SUB(ssys_handle_init_handle)
  end subroutine ssys_handle_init_handle

  ! ---------------------------------------------------------
  recursive subroutine ssys_handle_start(this, grid)
    type(ssys_handle_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid

    type(ssys_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_start)

    nullify(hndl)
    call base_handle__start__(this, grid)
    call ssys_handle_init(iter, this)
    do
      nullify(hndl)
      call ssys_handle_next(iter, hndl, ierr)
      if(ierr/=SSYS_HANDLE_OK)exit
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
    call ssys_handle_end(iter)
    nullify(hndl)

    POP_SUB(ssys_handle_start)
  end subroutine ssys_handle_start

  ! ---------------------------------------------------------
  subroutine ssys_handle_update(this)
    type(ssys_handle_t), intent(inout) :: this

    type(ssys_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: type, ierr

    PUSH_SUB(ssys_handle_update)

    nullify(hndl)
    call ssys_handle_init(iter, this)
    do
      nullify(hndl)
      call ssys_handle_next(iter, hndl, ierr)
      if(ierr/=SSYS_HANDLE_OK)exit
      call base_handle_get(hndl, type)
      select case(type)
      case(HNDL_TYPE_FRZN)
      case(HNDL_TYPE_LIVE)
      case default
        message(1)="Unknown subsystems type."
        call messages_fatal(1)
      end select
    end do
    call ssys_handle_end(iter)
    nullify(hndl)
    call base_handle__update__(this)

    POP_SUB(ssys_handle_update)
  end subroutine ssys_handle_update

  ! ---------------------------------------------------------
  subroutine ssys_handle_iterator_init(this, that)
    type(ssys_handle_iterator_t), intent(out) :: this
    type(ssys_handle_t),          intent(in)  :: that

    PUSH_SUB(ssys_handle_iterator_init)

    call base_handle_init(this, that)

    POP_SUB(ssys_handle_iterator_init)
  end subroutine ssys_handle_iterator_init

end module ssys_handle_m

!! Local Variables:
!! mode: f90
!! End:
