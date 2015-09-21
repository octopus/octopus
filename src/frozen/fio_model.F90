#include "global.h"

module fio_model_m

  use base_geom_m
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

  public ::            &
    fio_model__load__

  public ::           &
    !fio_model_init,   &
    fio_model_start,  &
    fio_model_stop,   &
    !fio_model_get,    &
    !fio_model_copy,   &
    fio_model_end

contains

  ! ---------------------------------------------------------
  subroutine fio_model_start(this, mpi_grp)
    type(base_model_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: mpi_grp

    type(base_system_t), pointer :: sys
    type(base_geom_t),   pointer :: geom
    type(simulation_t),  pointer :: sim
    type(grid_t),        pointer :: grid

    PUSH_SUB(fio_model_start)

    nullify(sys, geom, sim, grid)
    call base_model_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__new__(sim, grid, geom)
    call fio_simulation__start__(sim, grid, mpi_grp)
    nullify(sys, geom, sim)
    call base_model__start__(this, grid)

    POP_SUB(fio_model_start)
  end subroutine fio_model_start

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
  subroutine fio_model_stop(this)
    type(base_model_t), intent(inout) :: this

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model_stop)

    nullify(sim)
    call base_model__stop__(this)
    call base_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation__stop__(sim)
    nullify(sim)

    POP_SUB(fio_model_stop)
  end subroutine fio_model_stop

  ! ---------------------------------------------------------
  subroutine fio_model_end(this)
    type(base_model_t), intent(inout) :: this

    type(simulation_t), pointer :: sim

    PUSH_SUB(fio_model_end)

    nullify(sim)
    call base_model_get(this, sim)
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
