#include "global.h"

module fio_system_oct_m

  use base_geometry_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use fio_geometry_oct_m
  use fio_states_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::             &
    fio_system__init__, &
    fio_system__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_system__init__(this)
    type(base_system_t), intent(inout) :: this

    type(base_geometry_t), pointer :: geom
    type(base_states_t),   pointer :: stat

    PUSH_SUB(fio_system__init__)

    nullify(geom, stat)
    call base_system_get(this, geom)
    ASSERT(associated(geom))
    call fio_geometry__init__(geom)
    nullify(geom)
    call base_system_get(this, stat)
    ASSERT(associated(stat))
    call fio_states__init__(stat)
    nullify(stat)

    POP_SUB(fio_system__init__)
  end subroutine fio_system__init__
    
  ! ---------------------------------------------------------
  subroutine fio_system__load__(this)
    type(base_system_t), intent(inout) :: this

    type(base_states_t), pointer :: st

    PUSH_SUB(fio_system__load__)

    nullify(st)
    call base_system_get(this, st)
    ASSERT(associated(st))
    call fio_states__load__(st)
    nullify(st)

    POP_SUB(fio_system__load__)
  end subroutine fio_system__load__
    
end module fio_system_oct_m

!! Local Variables:
!! mode: f90
!! End:
