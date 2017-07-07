#include "global.h"

module frozen_geometry_oct_m

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
    frozen_geometry__build__

contains

  ! ---------------------------------------------------------
  subroutine frozen_geometry__build__(this, config)
    type(base_geometry_t), intent(inout) :: this
    type(json_object_t),   intent(in)    :: config

    type(geo_intrf_t), pointer :: igeo
    type(geometry_t),  pointer :: pgeo

    PUSH_SUB(frozen_geometry__build__)

    nullify(igeo, pgeo)
    call base_geometry_get(this, igeo)
    ASSERT(associated(igeo))
    call geo_intrf_new(igeo, pgeo, geometry__init__)
    ASSERT(associated(pgeo))
    nullify(igeo, pgeo)

    POP_SUB(frozen_geometry__build__)

  contains
    
    ! ---------------------------------------------------------
    subroutine geometry__init__(geom, space, gcfg)
      type(geometry_t),    intent(out) :: geom
      type(space_t),       intent(in)  :: space
      type(json_object_t), intent(in)  :: gcfg

      type(json_object_iterator_t)          :: iter
      character(len=BASE_GEOMETRY_NAME_LEN) :: name
      type(json_object_t),          pointer :: cnfg
      type(json_array_t),           pointer :: list
      type(geometry_t),             pointer :: pgeo
      type(geo_build_t)                     :: bgeo
      integer                               :: ierr

      PUSH_SUB(frozen_geometry__build__.geometry__init__)

      call geo_build_init(bgeo, space)
      call json_init(iter, config)
      do
        nullify(cnfg, list, pgeo)
        call json_next(iter, name, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call json_get(cnfg, "positions", list, ierr)
        ASSERT(ierr==JSON_OK)
        ASSERT(json_len(list)>0)
        call base_geometry_gets(this, name, pgeo)
        ASSERT(associated(pgeo))
        call geo_build_extend(bgeo, pgeo, list)
      end do
      call json_end(iter)
      nullify(cnfg, list, pgeo)
      call geo_build_export(bgeo, geom)
      call geo_build_end(bgeo)

      POP_SUB(frozen_geometry__build__.geometry__init__)

    end subroutine geometry__init__

  end subroutine frozen_geometry__build__

end module frozen_geometry_oct_m

!! Local Variables:
!! mode: f90
!! End:

