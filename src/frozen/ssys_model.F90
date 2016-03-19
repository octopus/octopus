#include "global.h"

module ssys_model_oct_m

  use base_model_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::         &
    ssys_model_acc

contains

  ! ---------------------------------------------------------
  subroutine ssys_model_acc(this)
    type(base_model_t), intent(inout) :: this

    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(ssys_model_acc)

    nullify(subs)
    call base_model__reset__(this)
    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model__acc__(this, subs)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__update__(this)

    POP_SUB(ssys_model_acc)
  end subroutine ssys_model_acc

end module ssys_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
