#include "global.h"

module frozen_system_oct_m

  use base_geometry_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use frozen_geometry_oct_m
  use frozen_states_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::                 &
    frozen_system__build__, &
    frozen_system__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_system__build__(this, config)
    type(base_system_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: config

    type(base_geometry_t), pointer :: geo

    PUSH_SUB(frozen_system__build__)

    nullify(geo)
    call base_system_get(this, geo)
    ASSERT(associated(geo))
    call frozen_geometry__build__(geo, config)
    nullify(geo)

    POP_SUB(frozen_system__build__)
  end subroutine frozen_system__build__
    
  ! ---------------------------------------------------------
  subroutine frozen_system__acc__(this, that, config)
    type(base_system_t), intent(inout) :: this !> frozen
    type(base_system_t), intent(in)    :: that !> fio
    type(json_object_t), intent(in)    :: config

    type(base_states_t),   pointer :: mst !> frozen
    type(base_states_t),   pointer :: sst !> fio

    PUSH_SUB(frozen_system__acc__)

    nullify(mst, sst)
    call base_system_get(this, mst)
    ASSERT(associated(mst))
    call base_system_get(that, sst)
    ASSERT(associated(sst))
    call frozen_states__acc__(mst, sst, config)
    nullify(mst, sst)

    POP_SUB(frozen_system__acc__)
  end subroutine frozen_system__acc__

end module frozen_system_oct_m

!! Local Variables:
!! mode: f90
!! End:
