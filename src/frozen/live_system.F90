#include "global.h"

module live_system_m

  use global_m
  use messages_m
  use profiling_m

  use base_system_m, only:                       &
    live_system_start  => base_system__start__,  &
    live_system_update => base_system__update__, &
    live_system_stop   => base_system__stop__

  use base_system_m, only:                &
    live_system_t    => base_system_t,    &
    live_system_init => base_system_init, &
    live_system_get  => base_system_get,  &
    live_system_copy => base_system_copy, &
    live_system_end  => base_system_end

  implicit none

  private
  public ::             &
    live_system_t,      &
    live_system_init,   &
    live_system_start,  &
    live_system_update, &
    live_system_stop,   &
    live_system_get,    &
    live_system_copy,   &
    live_system_end

end module live_system_m

!! Local Variables:
!! mode: f90
!! End:
