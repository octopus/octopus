#include "global.h"

module fio_states_m

  use global_m
  use messages_m
  use profiling_m

  use fio_density_m, only: &
    fio_density_t,         &
    fio_density_update

  use base_states_m, only:                    &
    fio_states_init  => base_states__init__,  &
    fio_states_start => base_states__start__, &
    fio_states_stop  => base_states__stop__,  &
    fio_states_copy  => base_states__copy__,  &
    fio_states_end   => base_states__end__

  use base_states_m, only:             &
    fio_states_t   => base_states_t,   &
    fio_states_get => base_states_get

  implicit none

  private
  public ::            &
    fio_states_t,      &
    fio_states_init,   &
    fio_states_start,  &
    fio_states_update, &
    fio_states_stop,   &
    fio_states_get,    &
    fio_states_copy,   &
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
