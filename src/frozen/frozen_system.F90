#include "global.h"

module frozen_system_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use fio_states_m, only: &
    fio_states_t

  use fio_system_m, only: &
    fio_system_t,         &
    fio_system_get

  use frozen_states_m, only: &
    frozen_states__update__

  use frozen_states_m, only: &
    frozen_states_t

  use base_system_m, only:                         &
    frozen_system_init   => base_system__init__,   &
    frozen_system_start  => base_system__start__,  &
    frozen_system_update => base_system__update__, &
    frozen_system_stop   => base_system__stop__,   &
    frozen_system_copy   => base_system__copy__,   &
    frozen_system_end    => base_system__end__

  use base_system_m, only:                &
    frozen_system_t   => base_system_t,   &
    frozen_system_get => base_system_get

  implicit none

  private
  public ::                  &
    frozen_system__update__

  public ::               &
    frozen_system_t,      &
    frozen_system_init,   &
    frozen_system_start,  &
    frozen_system_update, &
    frozen_system_stop,   &
    frozen_system_get,    &
    frozen_system_copy,   &
    frozen_system_end

contains

  ! ---------------------------------------------------------
  subroutine frozen_system__update__(this, that, config)
    type(frozen_system_t), intent(inout) :: this
    type(fio_system_t),    intent(in)    :: that
    type(json_object_t),   intent(in)    :: config
    !
    type(frozen_states_t), pointer :: mst
    type(fio_states_t),    pointer :: sst
    !
    PUSH_SUB(frozen_system__update__)
    nullify(mst, sst)
    call frozen_system_get(this, mst)
    ASSERT(associated(mst))
    call fio_system_get(that, sst)
    ASSERT(associated(sst))
    call frozen_states__update__(mst, sst, config)
    nullify(mst, sst)
    POP_SUB(frozen_system__update__)
    return
  end subroutine frozen_system__update__

end module frozen_system_m

!! Local Variables:
!! mode: f90
!! End:
