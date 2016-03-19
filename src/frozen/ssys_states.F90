#include "global.h"

module ssys_states_oct_m

  use base_states_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::          &
    ssys_states_acc

contains

  ! ---------------------------------------------------------
  subroutine ssys_states_acc(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_iterator_t)  :: iter
    type(base_states_t),  pointer :: subs
    integer                       :: ierr

    PUSH_SUB(ssys_states_acc)

    call base_states__reset__(this)
    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states__acc__(this, subs)
    end do
    call base_states_end(iter)
    nullify(subs)
    call base_states__update__(this)

    POP_SUB(ssys_states_acc)
  end subroutine ssys_states_acc

end module ssys_states_oct_m

!! Local Variables:
!! mode: f90
!! End:
