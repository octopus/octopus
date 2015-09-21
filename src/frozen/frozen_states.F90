#include "global.h"

module frozen_states_m

  use base_density_m
  use base_states_m
  use frozen_density_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::               &
    frozen_states__acc__

contains
    
  ! ---------------------------------------------------------
  subroutine frozen_states__acc__(this, that, config)
    type(base_states_t), intent(inout) :: this !> frozen
    type(base_states_t), intent(in)    :: that !> fio
    type(json_object_t), intent(in)    :: config

    type(base_density_t), pointer :: mrho !> frozen
    type(base_density_t), pointer :: srho !> fio

    PUSH_SUB(frozen_states__acc__)

    nullify(mrho, srho)
    call base_states_get(this, mrho)
    ASSERT(associated(mrho))
    call base_states_get(that, srho)
    ASSERT(associated(srho))
    call frozen_density__acc__(mrho, srho, config)
    nullify(mrho, srho)

    POP_SUB(frozen_states__acc__)
  end subroutine frozen_states__acc__

end module frozen_states_m

!! Local Variables:
!! mode: f90
!! End:
