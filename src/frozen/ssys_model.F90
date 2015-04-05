#include "global.h"

module ssys_model_m

  use global_m
  use messages_m
  use profiling_m

  use ssys_system_m, only: &
    ssys_system_t,         &
    ssys_system_update

  use ssys_hamiltonian_m, only: &
    ssys_hamiltonian_t,         &
    ssys_hamiltonian_update

  use base_model_m, only:                    &
    ssys_model_start => base_model__start__, &
    ssys_model_stop  => base_model__stop__

  use base_model_m, only:               &
    ssys_model_t    => base_model_t,    &
    ssys_model_init => base_model_init, &
    ssys_model_next => base_model_next, &
    ssys_model_get  => base_model_get,  &
    ssys_model_copy => base_model_copy, &
    ssys_model_end  => base_model_end

  use base_model_m, only:                           &
    ssys_model_iterator_t => base_model_iterator_t

  use base_model_m, only:                             &
    SSYS_MODEL_OK          => BASE_MODEL_OK,          &
    SSYS_MODEL_KEY_ERROR   => BASE_MODEL_KEY_ERROR,   &
    SSYS_MODEL_EMPTY_ERROR => BASE_MODEL_EMPTY_ERROR

  implicit none

  private
  public ::            &
    ssys_model_t,      &
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
  subroutine ssys_model_update(this)
    type(ssys_model_t), intent(inout) :: this

    type(ssys_hamiltonian_t), pointer :: hml
    type(ssys_system_t),      pointer :: sys

    PUSH_SUB(ssys_model_update)

    call ssys_model_get(this, sys)
    ASSERT(associated(sys))
    call ssys_system_update(sys)
    nullify(sys)
    call ssys_model_get(this, hml)
    ASSERT(associated(hml))
    call ssys_hamiltonian_update(hml)
    nullify(hml)

    POP_SUB(ssys_model_update)
  end subroutine ssys_model_update

end module ssys_model_m

!! Local Variables:
!! mode: f90
!! End:
