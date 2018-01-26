#include "global.h"

module frozen_model_oct_m

  use base_density_oct_m
  use base_hamiltonian_oct_m
  use base_model_oct_m
  use frozen_hamiltonian_oct_m
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

    type(base_density_t),     pointer :: mdns !> frozen
    type(base_hamiltonian_t), pointer :: mhml !> frozen
    type(base_density_t),     pointer :: sdns !> fio
    type(base_hamiltonian_t), pointer :: shml !> fio
    logical                           :: accu


    PUSH_SUB(frozen_model__acc__)

    nullify(mdns, mhml, sdns, shml)
    call base_model_get(this, mdns)
    ASSERT(associated(mdns))
    call base_model_get(that, sdns)
    ASSERT(associated(sdns))
    call base_density_get(sdns, use=accu)
    if(accu)then
      call base_density__acc__(mdns, sdns, config)
      call base_density_notify(mdns)
    end if
    nullify(mdns, sdns)
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
