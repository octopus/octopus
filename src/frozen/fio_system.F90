#include "global.h"

module fio_system_m

  use global_m
  use messages_m
  use profiling_m

  use space_m, only: space_t, space_init, space_copy, space_end
  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use fio_states_m, only: &
    fio_states_t,         &
    fio_states_update

  use bsyst_m, only:                 &
    fio_system_t     => bsyst_t,     &
    fio_system_init  => bsyst_init,  &
    fio_system_start => bsyst_start, &
    fio_system_stop  => bsyst_stop,  &
    fio_system_get   => bsyst_get,   &
    fio_system_copy  => bsyst_copy,  &
    fio_system_end   => bsyst_end

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
