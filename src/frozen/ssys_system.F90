#include "global.h"

module ssys_system_m

  use global_m
  use messages_m
  use profiling_m

  use ssys_states_m, only: &
    ssys_states_t,         &
    ssys_states_update

  use base_system_m, only:                     &
    ssys_system_start => base_system__start__, &
    ssys_system_stop  => base_system__stop__

  use base_system_m, only:                &
    ssys_system_t    => base_system_t,    &
    ssys_system_new  => base_system_new,  &
    ssys_system_del  => base_system_del,  &
    ssys_system_init => base_system_init, &
    ssys_system_next => base_system_next, &
    ssys_system_get  => base_system_get,  &
    ssys_system_copy => base_system_copy, &
    ssys_system_end  => base_system_end

  use base_system_m, only:                            &
    ssys_system_iterator_t => base_system_iterator_t

  use base_system_m, only:                              &
    SSYS_SYSTEM_OK          => BASE_SYSTEM_OK,          &
    SSYS_SYSTEM_KEY_ERROR   => BASE_SYSTEM_KEY_ERROR,   &
    SSYS_SYSTEM_EMPTY_ERROR => BASE_SYSTEM_EMPTY_ERROR

  implicit none

  private
  public ::             &
    ssys_system_t,      &
    ssys_system_new,    &
    ssys_system_del,    &
    ssys_system_init,   &
    ssys_system_start,  &
    ssys_system_update, &
    ssys_system_stop,   &
    ssys_system_next,   &
    ssys_system_get,    &
    ssys_system_copy,   &
    ssys_system_end

  public ::                 &
    ssys_system_iterator_t

  public ::                  &
    SSYS_SYSTEM_OK,          &
    SSYS_SYSTEM_KEY_ERROR,   &
    SSYS_SYSTEM_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_system_update(this)
    type(ssys_system_t), intent(inout) :: this
    !
    type(ssys_states_t), pointer :: subs
    !
    PUSH_SUB(ssys_system_update)
    nullify(subs)
    call ssys_system_get(this, subs)
    ASSERT(associated(subs))
    call ssys_states_update(subs)
    POP_SUB(ssys_system_update)
    return
  end subroutine ssys_system_update

end module ssys_system_m

!! Local Variables:
!! mode: f90
!! End:
