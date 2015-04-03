#include "global.h"

module ssys_functional_m

  use global_m
  use messages_m
  use profiling_m

  use base_functional_m, only:                           &
    ssys_functional_init   => base_functional__init__,   &
    ssys_functional_start  => base_functional__start__,  &
    ssys_functional_update => base_functional__update__, &
    ssys_functional_stop   => base_functional__stop__,   &
    ssys_functional_copy   => base_functional__copy__,   &
    ssys_functional_end    => base_functional__end__

  use base_functional_m, only:                    &
    ssys_functional_t    => base_functional_t,    &
    ssys_functional_eval => base_functional_eval, &
    ssys_functional_get  => base_functional_get

  use base_functional_m, only:                            &
    ssys_functional_intrpl_t => base_functional_intrpl_t

  implicit none

  private
  public ::                 &
    ssys_functional_t,      &
    ssys_functional_init,   &
    ssys_functional_start,  &
    ssys_functional_update, &
    ssys_functional_stop,   &
    ssys_functional_eval,   &
    ssys_functional_get,    &
    ssys_functional_copy,   &
    ssys_functional_end

end module ssys_functional_m

!! Local Variables:
!! mode: f90
!! End:
