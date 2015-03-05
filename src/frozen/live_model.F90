#include "global.h"

module live_model_m

  use global_m
  use messages_m
  use profiling_m

  use base_model_m, only:                      &
    live_model_init   => base_model__init__,   &
    live_model_start  => base_model__start__,  &
    live_model_update => base_model__update__, &
    live_model_stop   => base_model__stop__,   &
    live_model_copy   => base_model__copy__,   &
    live_model_end    => base_model__end__

  use base_model_m, only:             &
    live_model_t   => base_model_t,   &
    live_model_get => base_model_get

  implicit none

  private
  public ::            &
    live_model_t,      &
    live_model_init,   &
    live_model_start,  &
    live_model_update, &
    live_model_stop,   &
    live_model_get,    &
    live_model_copy,   &
    live_model_end

end module live_model_m

!! Local Variables:
!! mode: f90
!! End:
