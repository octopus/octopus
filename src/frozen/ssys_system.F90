#include "global.h"

module ssys_system_oct_m

  use base_system_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::          &
    ssys_system_acc

contains

  ! ---------------------------------------------------------
  subroutine ssys_system_acc(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_iterator_t)  :: iter
    type(base_system_t),  pointer :: subs
    integer                       :: ierr

    PUSH_SUB(ssys_system_acc)

    call base_system__reset__(this)
    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, subs, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call base_system__acc__(this, subs)
    end do
    call base_system_end(iter)
    nullify(subs)
    call base_system__update__(this)

    POP_SUB(ssys_system_acc)
  end subroutine ssys_system_acc

end module ssys_system_oct_m

!! Local Variables:
!! mode: f90
!! End:
