#include "global.h"

module ssys_ionic_m

  use global_m
  use messages_m
  use profiling_m

  use base_ionic_m, only:                   &
    ssys_ionic_t      => base_ionic_t,      &
    ssys_ionic_new    => base_ionic_new,    &
    ssys_ionic_del    => base_ionic_del,    &
    ssys_ionic_init   => base_ionic_init,   &
    ssys_ionic_update => base_ionic_update, &
    ssys_ionic_calc   => base_ionic_calc,   &
    ssys_ionic_next   => base_ionic_next,   &
    ssys_ionic_get    => base_ionic_get,    &
    ssys_ionic_copy   => base_ionic_copy,   &
    ssys_ionic_end    => base_ionic_end

  use base_ionic_m, only:                           &
    ssys_ionic_iterator_t => base_ionic_iterator_t

  use base_ionic_m, only:                       &
    SSYS_IONIC_NAME_LEN => BASE_IONIC_NAME_LEN

  use base_ionic_m, only:                             &
    SSYS_IONIC_OK          => BASE_IONIC_OK,          &
    SSYS_IONIC_KEY_ERROR   => BASE_IONIC_KEY_ERROR,   &
    SSYS_IONIC_EMPTY_ERROR => BASE_IONIC_EMPTY_ERROR

  implicit none

  private
  public ::            &
    ssys_ionic_t,      &
    ssys_ionic_new,    &
    ssys_ionic_del,    &
    ssys_ionic_init,   &
    ssys_ionic_update, &
    ssys_ionic_calc,   &
    ssys_ionic_next,   &
    ssys_ionic_get,    &
    ssys_ionic_copy,   &
    ssys_ionic_end

  public ::                &
    ssys_ionic_iterator_t

  public ::               &
    SSYS_IONIC_NAME_LEN

  public ::                 &
    SSYS_IONIC_OK,          &
    SSYS_IONIC_KEY_ERROR,   &
    SSYS_IONIC_EMPTY_ERROR

end module ssys_ionic_m

!! Local Variables:
!! mode: f90
!! End:
