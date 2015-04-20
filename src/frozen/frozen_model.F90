#include "global.h"

module frozen_model_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get

  use fio_system_m, only: &
    fio_system_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t

  use fio_model_m, only: &
    fio_model_t

  use fio_model_m, only: &
    fio_model_get

  use frozen_system_m, only: &
    frozen_system_t

  use frozen_system_m, only: &
    frozen_system__acc__

  use frozen_hamiltonian_m, only: &
    frozen_hamiltonian_t

  use frozen_hamiltonian_m, only: &
    frozen_hamiltonian__acc__

  use base_model_m, only:           &
    frozen_model_t => base_model_t

  use base_model_m, only:               &
    frozen_model_get => base_model_get

  implicit none

  private
  public ::         &
    frozen_model_t

  public ::              &
    frozen_model__acc__

  public ::           &
    frozen_model_get

contains

  ! ---------------------------------------------------------
  subroutine frozen_model__acc__(this, that, config)
    type(frozen_model_t), intent(inout) :: this
    type(fio_model_t),    intent(in)    :: that
    type(json_object_t),  intent(in)    :: config

    type(frozen_system_t),      pointer :: msys
    type(frozen_hamiltonian_t), pointer :: mhml
    type(fio_system_t),         pointer :: ssys
    type(fio_hamiltonian_t),    pointer :: shml

    PUSH_SUB(frozen_model__acc__)

    nullify(msys, mhml, ssys, shml)
    call frozen_model_get(this, msys)
    ASSERT(associated(msys))
    call fio_model_get(that, ssys)
    ASSERT(associated(ssys))
    call frozen_system__acc__(msys, ssys, config)
    nullify(msys, ssys)
    call frozen_model_get(this, mhml)
    ASSERT(associated(mhml))
    call fio_model_get(that, shml)
    ASSERT(associated(shml))
    call frozen_hamiltonian__acc__(mhml, shml, config)
    nullify(mhml, shml)

    POP_SUB(frozen_model__acc__)
  end subroutine frozen_model__acc__

end module frozen_model_m

!! Local Variables:
!! mode: f90
!! End:
