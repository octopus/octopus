#include "global.h"

module live_functional_m

  use global_m
  use messages_m
  use profiling_m

  use base_functional_m, only:                           &
    live_functional_init   => base_functional__init__,   &
    live_functional_start  => base_functional__start__,  &
    live_functional_update => base_functional__update__, &
    live_functional_stop   => base_functional__stop__,   &
    live_functional_copy   => base_functional__copy__,   &
    live_functional_end    => base_functional__end__

  use base_functional_m, only:                    &
    live_functional_t    => base_functional_t,    &
    live_functional_eval => base_functional_eval, &
    live_functional_get  => base_functional_get

  use base_functional_m, only:                            &
    live_functional_intrpl_t => base_functional_intrpl_t

  implicit none

  private
  public ::                 &
    live_functional_t,      &
    live_functional_init,   &
    live_functional_start,  &
    live_functional_update, &
    live_functional_stop,   &
    live_functional_eval,   &
    live_functional_get,    &
    live_functional_copy,   &
    live_functional_end

end module live_functional_m

!! Local Variables:
!! mode: f90
!! End:
