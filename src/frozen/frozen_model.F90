#include "global.h"

module frozen_model_oct_m

  use base_hamiltonian_oct_m
  use base_model_oct_m
  use base_system_oct_m
  use frozen_hamiltonian_oct_m
  use frozen_system_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::              &
    frozen_model__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_model__acc__(this, that, config)
    type(base_model_t),  intent(inout) :: this !> frozen
    type(base_model_t),  intent(in)    :: that !> fio
    type(json_object_t), intent(in)    :: config

    type(base_system_t),      pointer :: msys !> frozen
    type(base_hamiltonian_t), pointer :: mhml !> frozen
    type(base_system_t),      pointer :: ssys !> fio
    type(base_hamiltonian_t), pointer :: shml !> fio

    PUSH_SUB(frozen_model__acc__)

    nullify(msys, mhml, ssys, shml)
    call base_model_get(this, msys)
    ASSERT(associated(msys))
    call base_model_get(that, ssys)
    ASSERT(associated(ssys))
    call frozen_system__acc__(msys, ssys, config)
    nullify(msys, ssys)
    call base_model_get(this, mhml)
    ASSERT(associated(mhml))
    call base_model_get(that, shml)
    ASSERT(associated(shml))
    call frozen_hamiltonian__acc__(mhml, shml, config)
    nullify(mhml, shml)

    POP_SUB(frozen_model__acc__)
  end subroutine frozen_model__acc__

end module frozen_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
