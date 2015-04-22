#include "global.h"

module ssys_model_m

  use global_m
  use messages_m
  use profiling_m

  use ssys_system_m, only: &
    ssys_system_t

  use ssys_hamiltonian_m, only: &
    ssys_hamiltonian_t

  use base_model_m, only:         &
    ssys_model_t => base_model_t

  use base_model_m, only:  &
     base_model__update__, &
     base_model__reset__,  &
     base_model__acc__

 use base_model_m, only:                    &
    ssys_model_init   => base_model_init,   &
    ssys_model_start  => base_model_start,  &
    ssys_model_update => base_model_update, &
    ssys_model_stop   => base_model_stop,   &
    ssys_model_next   => base_model_next,   &
    ssys_model_get    => base_model_get,    &
    ssys_model_copy   => base_model_copy,   &
    ssys_model_end    => base_model_end

  use base_model_m, only:                           &
    ssys_model_iterator_t => base_model_iterator_t

  use base_model_m, only:                             &
    SSYS_MODEL_OK          => BASE_MODEL_OK,          &
    SSYS_MODEL_KEY_ERROR   => BASE_MODEL_KEY_ERROR,   &
    SSYS_MODEL_EMPTY_ERROR => BASE_MODEL_EMPTY_ERROR

  implicit none

  private
  public ::       &
    ssys_model_t

  public ::         &
    ssys_model_acc

  public ::            &
    ssys_model_init,   &
    ssys_model_start,  &
    ssys_model_update, &
    ssys_model_stop,   &
    ssys_model_next,   &
    ssys_model_get,    &
    ssys_model_copy,   &
    ssys_model_end

  public ::                &
    ssys_model_iterator_t

  public ::                 &
    SSYS_MODEL_OK,          &
    SSYS_MODEL_KEY_ERROR,   &
    SSYS_MODEL_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_model_acc(this)
    type(ssys_model_t), intent(inout) :: this

    type(ssys_model_iterator_t) :: iter
    type(ssys_model_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(ssys_model_acc)

    nullify(subs)
    call base_model__reset__(this)
    call ssys_model_init(iter, this)
    do
      nullify(subs)
      call ssys_model_next(iter, subs, ierr)
      if(ierr/=SSYS_MODEL_OK)exit
      call base_model__acc__(this, subs)
    end do
    call ssys_model_end(iter)
    nullify(subs)
    call base_model__update__(this)

    POP_SUB(ssys_model_update)
  end subroutine ssys_model_acc

end module ssys_model_m

!! Local Variables:
!! mode: f90
!! End:
