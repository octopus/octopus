#include "global.h"

module fio_model_m

  use global_m
  use messages_m
  use profiling_m

  use mpi_m, only: mpi_grp_t

  use base_geom_m, only: &
    base_geom_t

  use fio_grid_m, only: &
    fio_grid_t

  use fio_simulation_m, only: &
    fio_simulation_t

  use fio_simulation_m, only: &
    fio_simulation__new__,    &
    fio_simulation__del__,    &
    fio_simulation__start__,  &
    fio_simulation__stop__

  use fio_system_m, only: &
    fio_system_t

  use fio_system_m, only: &
    fio_system__load__

  use fio_system_m, only: &
    fio_system_get

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian__load__

  use base_model_m, only:        &
    fio_model_t => base_model_t

  use base_model_m, only: &
    base_model__start__,  &
    base_model__stop__,   &
    base_model__end__

  use base_model_m, only:                 &
    fio_model_init => base_model__init__, &
    fio_model_copy => base_model__copy__

  use base_model_m, only:            &
    fio_model_get => base_model_get

  implicit none

  private
  public ::      &
    fio_model_t

  public ::            &
    fio_model__load__

  public ::           &
    fio_model_init,   &
    fio_model_start,  &
    fio_model_stop,   &
    fio_model_get,    &
    fio_model_copy,   &
    fio_model_end

contains

  ! ---------------------------------------------------------
  subroutine fio_model_start(this, mpi_grp)
    type(fio_model_t), intent(inout) :: this
    type(mpi_grp_t),   intent(in)    :: mpi_grp

    type(fio_system_t),     pointer :: sys
    type(base_geom_t),      pointer :: geom
    type(fio_simulation_t), pointer :: sim
    type(fio_grid_t),       pointer :: grid

    PUSH_SUB(fio_model_start)

    nullify(sys, geom, sim, grid)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_get(sys, geom)
    ASSERT(associated(geom))
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__new__(sim, grid, geom)
    call fio_simulation__start__(sim, grid, mpi_grp)
    nullify(sys, geom, sim)
    call base_model__start__(this, grid)

    POP_SUB(fio_model_start)
  end subroutine fio_model_start

  ! ---------------------------------------------------------
  subroutine fio_model__load__(this)
    type(fio_model_t), intent(inout) :: this

    type(fio_system_t),      pointer :: sys
    type(fio_hamiltonian_t), pointer :: hml

    PUSH_SUB(fio_model__load__)

    nullify(sys, hml)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system__load__(sys)
    nullify(sys)
    call fio_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian__load__(hml)
    nullify(hml)

    POP_SUB(fio_model__load__)
  end subroutine fio_model__load__

  ! ---------------------------------------------------------
  subroutine fio_model_stop(this)
    type(fio_model_t), intent(inout) :: this

    type(fio_simulation_t), pointer :: sim

    PUSH_SUB(fio_model_stop)

    nullify(sim)
    call base_model__stop__(this)
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__stop__(sim)
    nullify(sim)

    POP_SUB(fio_model_stop)
  end subroutine fio_model_stop

  ! ---------------------------------------------------------
  subroutine fio_model_end(this)
    type(fio_model_t), intent(inout) :: this

    type(fio_simulation_t), pointer :: sim

    PUSH_SUB(fio_model_end)

    nullify(sim)
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__del__(sim)
    nullify(sim)
    call base_model__end__(this)

    POP_SUB(fio_model_end)
  end subroutine fio_model_end

end module fio_model_m

!! Local Variables:
!! mode: f90
!! End:
