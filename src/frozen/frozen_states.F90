#include "global.h"

module frozen_states_m

  use base_density_m
  use base_states_m
  use frozen_density_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  use base_states_m, only:            & 
    frozen_states_t => base_states_t

  use base_states_m, only:                & 
    frozen_states_get => base_states_get

  implicit none

  private

  public ::          &
    frozen_states_t

  public ::               &
    frozen_states__acc__

  public ::            &
    frozen_states_get

contains
    
  ! ---------------------------------------------------------
  subroutine frozen_states__acc__(this, that, config)
    type(frozen_states_t), intent(inout) :: this
    type(frozen_states_t), intent(in)    :: that !> fio
    type(json_object_t),   intent(in)    :: config

    type(frozen_density_t), pointer :: mrho
    type(base_density_t),   pointer :: srho !> fio

    PUSH_SUB(frozen_states__acc__)

    nullify(mrho, srho)
    call frozen_states_get(this, mrho)
    ASSERT(associated(mrho))
    call frozen_states_get(that, srho)
    ASSERT(associated(srho))
    call frozen_density__acc__(mrho, srho, config)
    nullify(mrho, srho)

    POP_SUB(frozen_states__acc__)
  end subroutine frozen_states__acc__

end module frozen_states_m

!! Local Variables:
!! mode: f90
!! End:
