#include "global.h"

module fio_system_m

  use base_system_m
  use base_states_m
  use fio_states_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::             &
    fio_system__load__

contains

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
    
end module fio_system_m

!! Local Variables:
!! mode: f90
!! End:
