#include "global.h"

module live_external_m

  use global_m
  use messages_m
  use profiling_m

  use base_external_m, only:                         &
    live_external_start  => base_external__start__,  &
    live_external_update => base_external__update__, &
    live_external_stop   => base_external__stop__

  use base_external_m, only:                  &
    live_external_t    => base_external_t,    &
    live_external_init => base_external_init, &
    live_external_eval => base_external_eval, &
    live_external_get  => base_external_get,  &
    live_external_copy => base_external_copy, &
    live_external_end  => base_external_end

  use base_external_m, only:                         &
    live_external_intrpl_t => base_external_intrpl_t

  implicit none

  private
  public ::               &
    live_external_t,      &
    live_external_init,   &
    live_external_start,  &
    live_external_update, &
    live_external_stop,   &
    live_external_eval,   &
    live_external_get,    &
    live_external_copy,   &
    live_external_end

end module live_external_m

!! Local Variables:
!! mode: f90
!! End:
