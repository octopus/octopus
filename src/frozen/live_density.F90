#include "global.h"

module live_density_m

  use global_m
  use messages_m
  use profiling_m

  use base_density_m, only:                        &
    live_density_init   => base_density__init__,   &
    live_density_start  => base_density__start__,  &
    live_density_update => base_density__update__, &
    live_density_stop   => base_density__stop__,   &
    live_density_copy   => base_density__copy__,   &
    live_density_end    => base_density__end__

  use base_density_m, only:                 &
    live_density_t    => base_density_t,    &
    live_density_eval => base_density_eval, &
    live_density_get  => base_density_get

  use base_density_m, only:                         &
    live_density_intrpl_t => base_density_intrpl_t

  implicit none

  private
  public ::              &
    live_density_t,      &
    live_density_init,   &
    live_density_start,  &
    live_density_update, &
    live_density_stop,   &
    live_density_eval,   &
    live_density_get,    &
    live_density_copy,   &
    live_density_end

end module live_density_m

!! Local Variables:
!! mode: f90
!! End:
