#include "global.h"

module ssys_geometry_oct_m

  use base_geometry_oct_m
  use geo_build_oct_m
  use geo_intrf_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m

  implicit none

  private

  public ::                   &
    ssys_geometry__build__

contains

  ! ---------------------------------------------------------
  subroutine ssys_geometry__build__(this)
    type(base_geometry_t), intent(inout) :: this

    type(geo_intrf_t), pointer :: igeo
    type(geo_build_t)          :: bgeo

    PUSH_SUB(ssys_geometry__build__)

    nullify(igeo)
    call base_geometry_get(this, igeo)
    ASSERT(associated(igeo))
    call geo_intrf_new(igeo, init)
    nullify(igeo)

    POP_SUB(ssys_geometry__build__)

  contains
    
    ! ---------------------------------------------------------
    subroutine init(geom, space, gcfg)
      type(geometry_t),    intent(out) :: geom
      type(space_t),       intent(in)  :: space
      type(json_object_t), intent(in)  :: gcfg


      PUSH_SUB(ssys_geometry__build__.init)

      ASSERT(json_len(gcfg)==0)
      call geo_build_init(bgeo, space)
      call base_geometry__build__(this, build)
      call geo_build_export(bgeo, geom)
      call geo_build_end(bgeo)

      POP_SUB(ssys_geometry__build__.init)
    end subroutine init

    ! ---------------------------------------------------------
    subroutine build(this, that, name)
      type(base_geometry_t), intent(inout) :: this
      type(base_geometry_t), intent(in)    :: that
      character(len=*),      intent(in)    :: name

      type(geometry_t), pointer :: pgeo

      PUSH_SUB(ssys_geometry__build__.build)

      nullify(pgeo)
      call base_geometry_get(that, pgeo)
      ASSERT(associated(pgeo))
      call geo_build_extend(bgeo, pgeo)
      nullify(pgeo)

      POP_SUB(ssys_geometry__build__.build)
    end subroutine build

  end subroutine ssys_geometry__build__

end module ssys_geometry_oct_m

!! Local Variables:
!! mode: f90
!! End:

