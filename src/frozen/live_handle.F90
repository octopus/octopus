#include "global.h"

module live_handle_m

  use base_handle_m
  use global_m
  use grid_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::         &
    HNDL_TYPE_LIVE

  public ::             &
    live_handle_init,   &
    live_handle_start,  &
    live_handle_stop,   &
    live_handle_copy,   &
    live_handle_end

  integer, parameter :: HNDL_TYPE_LIVE = 8

contains

  ! ---------------------------------------------------------
  subroutine live_handle_init(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(live_handle_init)

    call base_handle__init__(this, config)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_LIVE)
    call base_handle__init__(this)

    POP_SUB(live_handle_init)
  end subroutine live_handle_init

  ! ---------------------------------------------------------
  subroutine live_handle_start(this, grid)
    type(base_handle_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid

    PUSH_SUB(live_handle_start)

    call base_handle__start__(this, grid)
    call base_handle__update__(this)

    POP_SUB(live_handle_start)
  end subroutine live_handle_start

  ! ---------------------------------------------------------
  subroutine live_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(live_handle_start)

    call base_handle__stop__(this)

    POP_SUB(live_handle_stop)
  end subroutine live_handle_stop

  ! ---------------------------------------------------------
  subroutine live_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(live_handle_copy)

    call base_handle__copy__(this, that)

    POP_SUB(live_handle_copy)
  end subroutine live_handle_copy

  ! ---------------------------------------------------------
  subroutine live_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(live_handle_end)

    call base_handle__end__(this)

    POP_SUB(live_handle_end)
  end subroutine live_handle_end

end module live_handle_m

!! Local Variables:
!! mode: f90
!! End:
