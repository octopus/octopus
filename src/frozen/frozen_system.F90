#include "global.h"

module frozen_system_m

  use base_states_m
  use base_system_m
  use frozen_states_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  use base_system_m, only:            &
    frozen_system_t => base_system_t

  use base_system_m, only:                &
    frozen_system_get => base_system_get

  implicit none

  private

  public ::          &
    frozen_system_t

  public ::               &
    frozen_system__acc__

  public ::            &
    frozen_system_get

contains

  ! ---------------------------------------------------------
  subroutine frozen_system__acc__(this, that, config)
    type(frozen_system_t), intent(inout) :: this
    type(frozen_system_t), intent(in)    :: that !> fio
    type(json_object_t),   intent(in)    :: config

    type(frozen_states_t), pointer :: mst
    type(base_states_t),   pointer :: sst !> fio

    PUSH_SUB(frozen_system__update__)

    nullify(mst, sst)
    call frozen_system_get(this, mst)
    ASSERT(associated(mst))
    call frozen_system_get(that, sst)
    ASSERT(associated(sst))
    call frozen_states__acc__(mst, sst, config)
    nullify(mst, sst)

    POP_SUB(frozen_system__acc__)
  end subroutine frozen_system__acc__

end module frozen_system_m

!! Local Variables:
!! mode: f90
!! End:
