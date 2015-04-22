#include "global.h"

module ssys_system_m

  use global_m
  use messages_m
  use profiling_m

  use ssys_states_m, only: &
    ssys_states_t

  use base_system_m, only:          &
    ssys_system_t => base_system_t

  use base_system_m, only: &
    base_system__update__, &
    base_system__reset__,  &
    base_system__acc__

  use base_system_m, only:                    &
    ssys_system_new    => base_system_new,    &
    ssys_system_del    => base_system_del,    &
    ssys_system_init   => base_system_init,   &
    ssys_system_start  => base_system_start,  &
    ssys_system_update => base_system_update, &
    ssys_system_stop   => base_system_stop,   &
    ssys_system_next   => base_system_next,   &
    ssys_system_get    => base_system_get,    &
    ssys_system_copy   => base_system_copy,   &
    ssys_system_end    => base_system_end

  use base_system_m, only:                            &
    ssys_system_iterator_t => base_system_iterator_t

  use base_system_m, only:                              &
    SSYS_SYSTEM_OK          => BASE_SYSTEM_OK,          &
    SSYS_SYSTEM_KEY_ERROR   => BASE_SYSTEM_KEY_ERROR,   &
    SSYS_SYSTEM_EMPTY_ERROR => BASE_SYSTEM_EMPTY_ERROR

  implicit none

  private
  public ::        &
    ssys_system_t

  public ::          &
    ssys_system_acc

  public ::             &
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
  subroutine ssys_system_acc(this)
    type(ssys_system_t), intent(inout) :: this

    type(ssys_system_iterator_t)  :: iter
    type(ssys_system_t),  pointer :: subs
    integer                       :: ierr

    PUSH_SUB(ssys_system_acc)

    nullify(subs)
    call base_system__reset__(this)
    call ssys_system_init(iter, this)
    do
      nullify(subs)
      call ssys_system_next(iter, subs, ierr)
      if(ierr/=SSYS_SYSTEM_OK)exit
      call base_system__acc__(this, subs)
    end do
    call ssys_system_end(iter)
    nullify(subs)
    call base_system__update__(this)

    POP_SUB(ssys_system_acc)
  end subroutine ssys_system_acc

end module ssys_system_m

!! Local Variables:
!! mode: f90
!! End:
