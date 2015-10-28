#include "global.h"

module fio_model_m

  use base_hamiltonian_m
  use base_model_m
  use base_system_m
  use fio_hamiltonian_m
  use fio_simulation_m
  use fio_system_m
  use global_m
  use grid_m
  use messages_m
  use mpi_m
  use profiling_m
  use simulation_m

  implicit none

  private

  public ::             &
    fio_model__init__,  &
    fio_model__build__, &
    fio_model__start__, &
    fio_model__load__,  &
    fio_model__stop__,  &
    fio_model__copy__,  &
    fio_model__end__

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
  subroutine fio_model__build__(this)
    type(base_model_t), intent(inout) :: this

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model__build__)

    nullify(sim)
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__build__(sim)
    nullify(sim)

    POP_SUB(fio_model__build__)
  end subroutine fio_model__build__

  ! ---------------------------------------------------------
  subroutine fio_model__start__(this, mpi_grp)
    type(base_model_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: mpi_grp

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model__start__)

    nullify(sim)
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__start__(sim, mpi_grp)
    nullify(sim)

    POP_SUB(fio_model__start__)
  end subroutine fio_model__start__

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

  ! ---------------------------------------------------------
  subroutine fio_model__stop__(this)
    type(base_model_t), intent(inout) :: this

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model__stop__)

    nullify(sim)
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__stop__(sim)
    nullify(sim)

    POP_SUB(fio_model__stop__)
  end subroutine fio_model__stop__

  ! ---------------------------------------------------------
  subroutine fio_model__copy__(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    type(simulation_t), pointer :: osim, isim

    PUSH_SUB(fio_model__copy__)

    nullify(osim, isim)
    call base_model_get(that, isim)
    ASSERT(associated(isim))
    call base_model_get(this, osim)
    ASSERT(associated(osim))
    call fio_simulation__copy__(osim, isim)
    nullify(osim, isim)

    POP_SUB(fio_model__copy__)
  end subroutine fio_model__copy__

  ! ---------------------------------------------------------
  subroutine fio_model__end__(this)
    type(base_model_t), intent(inout) :: this

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model__end__)

    nullify(sim)
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__end__(sim)
    nullify(sim)

    POP_SUB(fio_model__end__)
  end subroutine fio_model__end__

end module fio_model_m

!! Local Variables:
!! mode: f90
!! End:
