#include "global.h"

module frozen_geometry_oct_m

  use base_geometry_oct_m
  use geo_build_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m

  implicit none

  private

  public ::                  &
    frozen_geometry__init__

contains

  ! ---------------------------------------------------------
  subroutine frozen_geometry__init__(this, config)
    type(base_geometry_t), intent(inout) :: this
    type(json_object_t),   intent(in)    :: config

    type(json_object_t), pointer :: dict
    type(space_t),       pointer :: space
    type(geo_build_t)            :: bgeo
    integer                      :: ierr

    PUSH_SUB(frozen_geometry__init__)

    nullify(dict, space)
    call base_geometry_get(this, space)
    ASSERT(associated(space))
    call geo_build_init(bgeo, space)
    nullify(space)
    call json_get(config, "subsystems", dict, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(dict)>0)
    call base_geometry__reset__(this)
    call base_geometry__build__(this, build)
    call base_geometry__update__(this)
    call base_geometry__init__(this, init)
    call geo_build_end(bgeo)
    nullify(dict)

    POP_SUB(frozen_geometry__init__)

  contains
    
    ! ---------------------------------------------------------
    subroutine init(this)
      type(geometry_t), intent(out) :: this

      PUSH_SUB(frozen_geometry__init__.init)

      call geo_build_export(bgeo, this)

      POP_SUB(frozen_geometry__init__.init)
    end subroutine init

    ! ---------------------------------------------------------
    subroutine build(this, name, that)
      type(base_geometry_t), intent(inout) :: this
      character(len=*),      intent(in)    :: name
      type(base_geometry_t), intent(in)    :: that

      type(json_object_t), pointer :: cnfg
      type(json_array_t),  pointer :: list
      type(geometry_t),    pointer :: pgeo
      integer                      :: ierr

      PUSH_SUB(frozen_geometry__init__.build)

      nullify(cnfg, list, pgeo)
      call base_geometry_get(this, pgeo)
      ASSERT(.not.associated(pgeo))
      nullify(pgeo)
      call json_get(dict, trim(adjustl(name)), cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "positions", list, ierr)
      ASSERT(ierr==JSON_OK)
      ASSERT(associated(list))
      ASSERT(json_len(list)>0)
      nullify(cnfg)
      call base_geometry_get(that, pgeo)
      ASSERT(associated(pgeo))
      call geo_build_extend(bgeo, pgeo, list)
      nullify(list, pgeo)

      POP_SUB(frozen_geometry__init__.build)
    end subroutine build

  end subroutine frozen_geometry__init__

end module frozen_geometry_oct_m

!! Local Variables:
!! mode: f90
!! End:

