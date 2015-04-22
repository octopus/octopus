#include "global.h"

module ssys_density_m

  use global_m
  use messages_m
  use profiling_m

  use base_density_m, only:           &
    ssys_density_t => base_density_t

  use base_density_m, only: &
    base_density__update__, &
    base_density__reset__,  &
    base_density__acc__

  use base_density_m, only:                     &
    ssys_density_new    => base_density_new,    &
    ssys_density_del    => base_density_del,    &
    ssys_density_init   => base_density_init,   &
    ssys_density_start  => base_density_start,  &
    ssys_density_update => base_density_update, &
    ssys_density_stop   => base_density_stop,   &
    ssys_density_next   => base_density_next,   &
    ssys_density_get    => base_density_get,    &
    ssys_density_copy   => base_density_copy,   &
    ssys_density_end    => base_density_end

  use base_density_m, only:                             &
    ssys_density_iterator_t => base_density_iterator_t

  use base_density_m, only:                         &
    SSYS_DENSITY_NAME_LEN => BASE_DENSITY_NAME_LEN

  use base_density_m, only:                               &
    SSYS_DENSITY_OK          => BASE_DENSITY_OK,          &
    SSYS_DENSITY_KEY_ERROR   => BASE_DENSITY_KEY_ERROR,   &
    SSYS_DENSITY_EMPTY_ERROR => BASE_DENSITY_EMPTY_ERROR

  implicit none

  private
  public ::         &
    ssys_density_t

  public ::           &
    ssys_density_acc

  public ::              &
    ssys_density_new,    &
    ssys_density_del,    &
    ssys_density_init,   &
    ssys_density_start,  &
    ssys_density_update, &
    ssys_density_stop,   &
    ssys_density_next,   &
    ssys_density_get,    &
    ssys_density_copy,   &
    ssys_density_end

  public ::                  &
    ssys_density_iterator_t

  public ::                &
    SSYS_DENSITY_NAME_LEN

  public ::                   &
    SSYS_DENSITY_OK,          &
    SSYS_DENSITY_KEY_ERROR,   &
    SSYS_DENSITY_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_density_acc(this)
    type(ssys_density_t), intent(inout) :: this

    type(ssys_density_iterator_t) :: iter
    type(ssys_density_t), pointer :: subs
    integer                       :: ierr

    PUSH_SUB(ssys_density_acc)

    nullify(subs)
    call base_density__reset__(this)
    call ssys_density_init(iter, this)
    do
      nullify(subs)
      call ssys_density_next(iter, subs, ierr)
      if(ierr/=SSYS_DENSITY_OK)exit
      call base_density__acc__(this, subs)
    end do
    call ssys_density_end(iter)
    nullify(subs)
    call base_density__update__(this)

    POP_SUB(ssys_density_acc)
  end subroutine ssys_density_acc

end module ssys_density_m

!! Local Variables:
!! mode: f90
!! End:
