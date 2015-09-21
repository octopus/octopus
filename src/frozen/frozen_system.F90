#include "global.h"

module frozen_system_m

  use base_states_m
  use base_system_m
  use frozen_states_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::               &
    frozen_system__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_system__acc__(this, that, config)
    type(base_system_t), intent(inout) :: this !> frozen
    type(base_system_t), intent(in)    :: that !> fio
    type(json_object_t), intent(in)    :: config

    type(base_states_t), pointer :: mst !> frozen
    type(base_states_t), pointer :: sst !> fio

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

end module frozen_system_m

!! Local Variables:
!! mode: f90
!! End:
