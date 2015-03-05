#include "global.h"

module fio_system_m

  use global_m
  use messages_m
  use profiling_m

  use fio_states_m, only: &
    fio_states_t,         &
    fio_states_update

  use base_system_m, only:                    &
    fio_system_init  => base_system__init__,  &
    fio_system_start => base_system__start__, &
    fio_system_stop  => base_system__stop__,  &
    fio_system_copy  => base_system__copy__,  &
    fio_system_end   => base_system__end__

  use base_system_m, only:             &
    fio_system_t   => base_system_t,   &
    fio_system_get => base_system_get

  implicit none

  private
  public ::            &
    fio_system_t,      &
    fio_system_init,   &
    fio_system_start,  &
    fio_system_update, &
    fio_system_stop,   &
    fio_system_get,    &
    fio_system_copy,   &
    fio_system_end

contains

  ! ---------------------------------------------------------
  subroutine fio_system_update(this)
    type(fio_system_t), intent(inout) :: this
    !
    type(fio_states_t), pointer :: st
    !
    PUSH_SUB(fio_system_update)
    nullify(st)
    call fio_system_get(this, st)
    ASSERT(associated(st))
    call fio_states_update(st)
    nullify(st)
    POP_SUB(fio_system_update)
    return
  end subroutine fio_system_update
    
end module fio_system_m

!! Local Variables:
!! mode: f90
!! End:
