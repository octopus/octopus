#include "global.h"

module fio_model_m

  use base_hamiltonian_m
  use base_model_m
  use base_system_m
  use fio_hamiltonian_m
  use fio_system_m
  use global_m
  use messages_m
  use mpi_m
  use profiling_m

  implicit none

  private

  public ::             &
    fio_model__init__,  &
    fio_model__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_model__init__(this)
    type(base_model_t), intent(inout) :: this

    type(base_system_t), pointer :: sys

    PUSH_SUB(fio_model__init__)

    nullify(sys)
    call base_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system__init__(sys)
    nullify(sys)

    POP_SUB(fio_model__init__)
  end subroutine fio_model__init__

  ! ---------------------------------------------------------
  subroutine fio_model__load__(this)
    type(base_model_t), intent(inout) :: this

    type(base_system_t),      pointer :: sys
    type(base_hamiltonian_t), pointer :: hml

    PUSH_SUB(fio_model__load__)

    nullify(sys, hml)
    call base_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system__load__(sys)
    nullify(sys)
    call base_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian__load__(hml)
    nullify(hml)

    POP_SUB(fio_model__load__)
  end subroutine fio_model__load__

end module fio_model_m

!! Local Variables:
!! mode: f90
!! End:
