#include "global.h"

module fio_states_m

  use global_m
  use messages_m
  use profiling_m

  use fio_density_m, only: &
    fio_density_t

  use fio_density_m, only: &
    fio_density__load__

  use base_states_m, only:         &
    fio_states_t => base_states_t

  use base_states_m, only:             &
    fio_states_get => base_states_get

  implicit none

  private
  public ::       &
    fio_states_t

  public ::             &
    fio_states__load__

   public ::         &
     fio_states_get

contains

  ! ---------------------------------------------------------
  subroutine fio_states__load__(this)
    type(fio_states_t), intent(inout) :: this

    type(fio_density_t), pointer :: density

    PUSH_SUB(fio_states__load__)

    nullify(density)
    call fio_states_get(this, density)
    ASSERT(associated(density))
    call fio_density__load__(density)

    POP_SUB(fio_states__load__)
  end subroutine fio_states__load__

end module fio_states_m

!! Local Variables:
!! mode: f90
!! End:
