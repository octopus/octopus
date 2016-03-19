#include "global.h"

module fio_model_oct_m

  use base_hamiltonian_oct_m
  use base_model_oct_m
  use base_system_oct_m
  use fio_hamiltonian_oct_m
  use fio_system_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m

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

end module fio_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
