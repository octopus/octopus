#include "global.h"

module fio_system_m

  use global_m
  use messages_m
  use profiling_m

  use fio_states_m, only: &
    fio_states_t

  use fio_states_m, only: &
    fio_states__load__

  use base_system_m, only:         &
    fio_system_t => base_system_t

  use base_system_m, only:             &
    fio_system_get => base_system_get

  implicit none

  private
  public ::       &
    fio_system_t

  public ::             &
    fio_system__load__

  public ::         &
    fio_system_get

contains

  ! ---------------------------------------------------------
  subroutine fio_system__load__(this)
    type(fio_system_t), intent(inout) :: this

    type(fio_states_t), pointer :: st

    PUSH_SUB(fio_system__load__)

    nullify(st)
    call fio_system_get(this, st)
    ASSERT(associated(st))
    call fio_states__load__(st)
    nullify(st)

    POP_SUB(fio_system__load__)
  end subroutine fio_system__load__
    
end module fio_system_m

!! Local Variables:
!! mode: f90
!! End:
