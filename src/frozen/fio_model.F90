#include "global.h"

module fio_model_oct_m

  use base_density_oct_m
  use base_hamiltonian_oct_m
  use base_model_oct_m
  use fio_density_oct_m
  use fio_hamiltonian_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::            &
    fio_model__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_model__load__(this)
    type(base_model_t), intent(inout) :: this

    type(base_density_t),     pointer :: pdns
    type(base_hamiltonian_t), pointer :: hml

    PUSH_SUB(fio_model__load__)

    nullify(pdns, hml)
    call base_model_get(this, pdns)
    ASSERT(associated(pdns))
    call fio_density__load__(pdns)
    nullify(pdns)
    call base_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian__load__(hml)
    nullify(hml)

    POP_SUB(fio_model__load__)
  end subroutine fio_model__load__

end module fio_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
