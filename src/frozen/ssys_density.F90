#include "global.h"

module ssys_density_oct_m

  use base_density_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::           &
    ssys_density_acc

contains

  ! ---------------------------------------------------------
  subroutine ssys_density_acc(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs
    integer                       :: ierr

    PUSH_SUB(ssys_density_acc)

    call base_density__reset__(this)
    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, subs, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density__acc__(this, subs)
    end do
    call base_density_end(iter)
    nullify(subs)
    call base_density__update__(this)

    POP_SUB(ssys_density_acc)
  end subroutine ssys_density_acc

end module ssys_density_oct_m

!! Local Variables:
!! mode: f90
!! End:
