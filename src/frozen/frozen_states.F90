#include "global.h"

module frozen_states_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_get

  use fio_density_m, only: &
    fio_density_t

  use fio_states_m, only: &
    fio_states_t,         &
    fio_states_get

  use frozen_density_m, only: &
    frozen_density__update__

  use frozen_density_m, only: &
    frozen_density_t

  use base_states_m, only:                         & 
    frozen_states_start  => base_states__start__,  &
    frozen_states_update => base_states__update__, &
    frozen_states_stop   => base_states__stop__

  use base_states_m, only:                  & 
    frozen_states_t    => base_states_t,    &
    frozen_states_init => base_states_init, &
    frozen_states_get  => base_states_get,  &
    frozen_states_copy => base_states_copy, &
    frozen_states_end  => base_states_end

  implicit none

  private
  public ::                  &
    frozen_states__update__

  public ::               &
    frozen_states_t,      &
    frozen_states_init,   &
    frozen_states_start,  &
    frozen_states_update, &
    frozen_states_stop,   &
    frozen_states_get,    &
    frozen_states_copy,   &
    frozen_states_end

contains
    
  ! ---------------------------------------------------------
  subroutine frozen_states__update__(this, that, config)
    type(frozen_states_t), intent(inout) :: this
    type(fio_states_t),    intent(in)    :: that
    type(json_object_t),   intent(in)    :: config
    !
    type(frozen_density_t), pointer :: mrho
    type(fio_density_t),    pointer :: srho
    !
    PUSH_SUB(frozen_states__update__)
    nullify(mrho, srho)
    call frozen_states_get(this, mrho)
    ASSERT(associated(mrho))
    call fio_states_get(that, srho)
    ASSERT(associated(srho))
    call frozen_density__update__(mrho, srho, config)
    nullify(mrho, srho)
    POP_SUB(frozen_states__update__)
    return
  end subroutine frozen_states__update__

end module frozen_states_m

!! Local Variables:
!! mode: f90
!! End:
