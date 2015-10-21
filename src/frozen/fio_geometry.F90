#include "global.h"

module fio_geometry_m

  use base_geometry_m
  use geo_intrf_m
  use geometry_m
  use global_m
  use json_m
  use messages_m
  use profiling_m
  use space_m

  implicit none

  private

  public ::               &
    fio_geometry__init__

contains

  ! ---------------------------------------------------------
  subroutine fio_geometry__init__(this)
    type(base_geometry_t), intent(inout) :: this

    type(geo_intrf_t), pointer :: igeo
    type(geometry_t),  pointer :: pgeo

    PUSH_SUB(fio_geometry__init__)

    nullify(igeo, pgeo)
    call base_geometry_get(this, igeo)
    ASSERT(associated(igeo))
    call geo_intrf_new(igeo, pgeo, geometry__init__)
    ASSERT(associated(pgeo))
    nullify(igeo, pgeo)

    POP_SUB(fio_geometry__init__)
  end subroutine fio_geometry__init__

  subroutine geometry__init__(this, space, config)
    type(geometry_t),    intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(geometry__init__)

    call geometry_init_from_data_object(this, space, config)

    POP_SUB(geometry__init__)

  end subroutine geometry__init__

end module fio_geometry_m

