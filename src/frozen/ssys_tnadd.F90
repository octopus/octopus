#include "global.h"

module ssys_tnadd_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  use base_functional_m, only:         &
    ssys_tnadd_t => base_functional_t

  use base_functional_m, only: &
    base_functional__sub__

  use base_functional_m, only: &
    base_functional_calc,      &
    base_functional_set

  use base_functional_m, only:               &
    ssys_tnadd_init => base_functional_init, &
    ssys_tnadd_get  => base_functional_get,  &
    ssys_tnadd_copy => base_functional_copy, &
    ssys_tnadd_end  => base_functional_end

  implicit none

  private
  public ::       &
    ssys_tnadd_t

  public ::          &
    ssys_tnadd_calc

  public ::            &
    ssys_tnadd_init,   &
    ssys_tnadd_get,    &
    ssys_tnadd_copy,   &
    ssys_tnadd_end

contains

  ! ---------------------------------------------------------
  subroutine ssys_tnadd_calc(this)
    type(ssys_tnadd_t), intent(inout) :: this

    type(ssys_tnadd_t), pointer :: live

    PUSH_SUB(ssys_tnadd_calc)

    nullify(live)
    call base_functional_calc(this)
    !call base_functional_gets(this, "live", live)
    ASSERT(associated(live))
    call base_functional__sub__(this, live)
    nullify(live)

    POP_SUB(ssys_tnadd_calc)
  end subroutine ssys_tnadd_calc

end module ssys_tnadd_m

!! Local Variables:
!! mode: f90
!! End:
