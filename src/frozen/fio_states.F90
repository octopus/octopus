#include "global.h"

module fio_states_m

  use global_m
  use messages_m
  use profiling_m

  use fio_density_m, only: &
    fio_density_t,         &
    fio_density_update

  use bstts_m, only:                           &
    fio_states_t          => bstts_t,          &
    fio_states_init       => bstts_init,       &
    fio_states_start      => bstts_start,      &
    fio_states_stop       => bstts_stop,       &
    fio_states_get        => bstts_get,        &
    fio_states_get_charge => bstts_get_charge, &
    fio_states_copy       => bstts_copy,       &
    fio_states_end        => bstts_end

  implicit none

  private
  public ::                &
    fio_states_t,          &
    fio_states_init,       &
    fio_states_start,      &
    fio_states_update,     &
    fio_states_stop,       &
    fio_states_get,        &
    fio_states_get_charge, &
    fio_states_copy,       &
    fio_states_end

contains

  ! ---------------------------------------------------------
  subroutine fio_states_update(this)
    type(fio_states_t), intent(inout) :: this
    !
    type(fio_density_t), pointer :: density
    !
    PUSH_SUB(fio_states_update)
    nullify(density)
    call fio_states_get(this, density)
    ASSERT(associated(density))
    call fio_density_update(density)
    POP_SUB(fio_states_update)
    return
  end subroutine fio_states_update

end module fio_states_m

!! Local Variables:
!! mode: f90
!! End:
