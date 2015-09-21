#include "global.h"

module fio_states_m

  use base_density_m
  use base_states_m
  use fio_density_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::             &
    fio_states__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_states__load__(this)
    type(base_states_t), intent(inout) :: this

    type(base_density_t), pointer :: density

    PUSH_SUB(fio_states__load__)

    nullify(density)
    call base_states_get(this, density)
    ASSERT(associated(density))
    call fio_density__load__(density)

    POP_SUB(fio_states__load__)
  end subroutine fio_states__load__

end module fio_states_m

!! Local Variables:
!! mode: f90
!! End:
