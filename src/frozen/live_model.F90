#include "global.h"

module live_model_oct_m

  use base_model_oct_m
  use base_system_oct_m
  use geometry_oct_m
  use global_oct_m
  use live_system_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::             &
    live_model__init__

contains

  ! ---------------------------------------------------------
  subroutine live_model__init__(this, geo)
    type(base_model_t), intent(inout) :: this
    type(geometry_t),   intent(in)    :: geo

    type(base_system_t), pointer :: that

    PUSH_SUB(live_model__init__)

    nullify(that)
    call base_model_get(this, that)
    ASSERT(associated(that))
    call live_system__init__(that, geo)
    nullify(that)

    POP_SUB(live_model__init__)
  end subroutine live_model__init__

end module live_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
