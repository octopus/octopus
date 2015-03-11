#include "global.h"

module ssys_external_m

  use global_m
  use messages_m
  use profiling_m

  use base_external_m, only: &
    base_external__update__, &
    base_external__reset__,  &
    base_external__acc__

  use base_external_m, only: &
    base_external_t

  use base_external_m, only:                       &
    ssys_external_start => base_external__start__, &
    ssys_external_stop  => base_external__stop__

  use base_external_m, only:                  &
    ssys_external_t    => base_external_t,    &
    ssys_external_init => base_external_init, &
    ssys_external_next => base_external_next, &
    ssys_external_eval => base_external_eval, &
    ssys_external_get  => base_external_get,  &
    ssys_external_copy => base_external_copy, &
    ssys_external_end  => base_external_end

  use base_external_m, only:                             &
    ssys_external_iterator_t => base_external_iterator_t

  use base_external_m, only:                         &
    ssys_external_intrpl_t => base_external_intrpl_t

  use base_external_m, only:                                &
    SSYS_EXTERNAL_OK          => BASE_EXTERNAL_OK,          &
    SSYS_EXTERNAL_KEY_ERROR   => BASE_EXTERNAL_KEY_ERROR,   &
    SSYS_EXTERNAL_EMPTY_ERROR => BASE_EXTERNAL_EMPTY_ERROR

  implicit none

  private
  public ::               &
    ssys_external_t,      &
    ssys_external_init,   &
    ssys_external_start,  &
    ssys_external_update, &
    ssys_external_stop,   &
    ssys_external_next,   &
    ssys_external_eval,   &
    ssys_external_get,    &
    ssys_external_copy,   &
    ssys_external_end

  public ::                   &
    ssys_external_iterator_t

  public ::                 &
    ssys_external_intrpl_t

  public ::                    &
    SSYS_EXTERNAL_OK,          &
    SSYS_EXTERNAL_KEY_ERROR,   &
    SSYS_EXTERNAL_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_external_update(this)
    type(ssys_external_t), intent(inout) :: this
    !
    type(ssys_external_iterator_t) :: iter
    type(base_external_t), pointer :: subs
    integer                        :: ierr
    !
    PUSH_SUB(ssys_external_update)
    nullify(subs)
    call base_external__reset__(this)
    call ssys_external_init(iter, this)
    do
      nullify(subs)
      call ssys_external_next(iter, subs, ierr)
      if(ierr/=SSYS_EXTERNAL_OK)exit
      call base_external__acc__(this, subs)
    end do
    call ssys_external_end(iter)
    nullify(subs)
    call base_external__update__(this)
    POP_SUB(ssys_external_update)
    return
  end subroutine ssys_external_update

end module ssys_external_m

!! Local Variables:
!! mode: f90
!! End:
