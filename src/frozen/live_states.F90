#include "global.h"

module live_states_m

  use global_m
  use messages_m
  use profiling_m

  use base_states_m, only:                       &
    live_states_start  => base_states__start__,  &
    live_states_update => base_states__update__, &
    live_states_stop   => base_states__stop__

  use base_states_m, only:                &
    live_states_t    => base_states_t,    &
    live_states_init => base_states_init, &
    live_states_get  => base_states_get,  &
    live_states_copy => base_states_copy, &
    live_states_end  => base_states_end

  implicit none

  private
  public ::             &
    live_states_t,      &
    live_states_init,   &
    live_states_start,  &
    live_states_update, &
    live_states_stop,   &
    live_states_get,    &
    live_states_copy,   &
    live_states_end

end module live_states_m

!! Local Variables:
!! mode: f90
!! End:
