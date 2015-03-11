#include "global.h"

module frozen_ionic_m

  use global_m
  use messages_m
  use profiling_m

  use base_ionic_m, only:                        &
    frozen_ionic_init   => base_ionic__init__,   &
    frozen_ionic_update => base_ionic__update__, &
    frozen_ionic_copy   => base_ionic__copy__,   &
    frozen_ionic_end    => base_ionic__end__

  use base_ionic_m, only:                 &
    frozen_ionic_t    => base_ionic_t,    &
    frozen_ionic_calc => base_ionic_calc, &
    frozen_ionic_get  => base_ionic_get

  implicit none

  private
  public ::              &
    frozen_ionic_t,      &
    frozen_ionic_init,   &
    frozen_ionic_update, &
    frozen_ionic_calc,   &
    frozen_ionic_get,    &
    frozen_ionic_copy,   &
    frozen_ionic_end

end module frozen_ionic_m

!! Local Variables:
!! mode: f90
!! End:
