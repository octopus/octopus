#include "global.h"

module live_geometry_m

  use base_geometry_m
  use geo_intrf_m
  use geometry_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::                &
    live_geometry__init__

contains

  ! ---------------------------------------------------------
  subroutine live_geometry__init__(this, geo)
    type(base_geometry_t), intent(inout) :: this
    type(geometry_t),      intent(in)    :: geo

    type(geo_intrf_t), pointer :: igeo

    PUSH_SUB(live_geometry__init__)

    nullify(igeo)
    call base_geometry_get(this, igeo)
    ASSERT(associated(igeo))
    call geo_intrf_set(igeo, geo)
    nullify(igeo)

    POP_SUB(live_geometry__init__)
  end subroutine live_geometry__init__

end module live_geometry_m

!! Local Variables:
!! mode: f90
!! End:
