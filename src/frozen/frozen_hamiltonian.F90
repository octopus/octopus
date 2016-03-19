#include "global.h"

module frozen_hamiltonian_oct_m

  use base_hamiltonian_oct_m
  use base_potential_oct_m
  use base_system_oct_m
  use frozen_external_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m

  implicit none

  private

  public ::                    &
    frozen_hamiltonian__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian__acc__(this, that, config)
    type(base_hamiltonian_t), intent(inout) :: this !> frozen
    type(base_hamiltonian_t), intent(in)    :: that !> fio
    type(json_object_t),      intent(in)    :: config

    type(base_potential_t), pointer :: mept !> frozen
    type(base_potential_t), pointer :: sept !> fio
    real(kind=wp)                   :: energy, enrg

    PUSH_SUB(frozen_hamiltonian__acc__)

    call base_hamiltonian_get(this, energy=energy)
    call base_hamiltonian_get(that, energy=enrg)
    call base_hamiltonian_set(this, energy=(energy+enrg))
    nullify(mept, sept)
    call base_hamiltonian_get(this, "external", mept)
    ASSERT(associated(mept))
    call base_hamiltonian_get(that, "external", sept)
    ASSERT(associated(sept))
    call frozen_external__acc__(mept, sept, config)
    nullify(mept, sept)

    POP_SUB(frozen_hamiltonian__acc__)
  end subroutine frozen_hamiltonian__acc__

end module frozen_hamiltonian_oct_m

!! Local Variables:
!! mode: f90
!! End:
