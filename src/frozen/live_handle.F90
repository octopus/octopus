#include "global.h"

module live_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use base_handle_m, only: &
    base_handle__init__

  use base_handle_m, only:                       &
    live_handle_start  => base_handle__start__,  &
    live_handle_update => base_handle__update__, &
    live_handle_stop   => base_handle__stop__,   &
    live_handle_copy   => base_handle__copy__,   &
    live_handle_end    => base_handle__end__

  use base_handle_m, only:              &
    live_handle_t   => base_handle_t,   &
    live_handle_get => base_handle_get

  implicit none

  private

  public ::             &
    live_handle_t,      &
    live_handle_init,   &
    live_handle_start,  &
    live_handle_update, &
    live_handle_stop,   &
    live_handle_get,    &
    live_handle_copy,   &
    live_handle_end

  integer, public, parameter :: HNDL_TYPE_LIVE = 8

contains

  ! ---------------------------------------------------------
  subroutine live_handle_init(this, config)
    type(live_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    integer :: type
    !
    PUSH_SUB(live_handle_init)
    call base_handle__init__(this, config)
    call live_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_LIVE)
    call base_handle__init__(this)
    POP_SUB(live_handle_init)
    return
  end subroutine live_handle_init

end module live_handle_m

!! Local Variables:
!! mode: f90
!! End:
