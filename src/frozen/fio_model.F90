#include "global.h"

module fio_model_m

  use global_m
  use messages_m
  use profiling_m

  use geometry_m, only: geometry_t
  use mpi_m,      only: mpi_grp_t
  use space_m,    only: space_t

  use base_geom_m, only:           &
    fio_geom_t   => base_geom_t,   &
    fio_geom_get => base_geom_get

  use fio_simulation_m, only: &
    fio_simulation_t,         &
    fio_simulation_start,     &
    fio_simulation_stop

  use fio_system_m, only: &
    fio_system_t,         &
    fio_system_start,     &
    fio_system_update,    &
    fio_system_stop,      &
    fio_system_get

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t,         &
    fio_hamiltonian_start,     &
    fio_hamiltonian_update,    &
    fio_hamiltonian_stop,      &
    fio_hamiltonian_end

  use base_model_m, only:              &
    fio_model_t    => base_model_t,    &
    fio_model_init => base_model_init, &
    fio_model_get  => base_model_get,  &
    fio_model_copy => base_model_copy, &
    fio_model_end  => base_model_end

  implicit none

  private
  public ::           &
    fio_model_t,      &
    fio_model_init,   &
    fio_model_start,  &
    fio_model_update, &
    fio_model_stop,   &
    fio_model_get,    &
    fio_model_copy,   &
    fio_model_end

contains

  ! ---------------------------------------------------------
  subroutine fio_model_start(this, mpi_grp)
    type(fio_model_t), intent(inout) :: this
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    !
    type(fio_system_t),      pointer :: sys
    type(fio_geom_t),        pointer :: geom
    type(geometry_t),        pointer :: geo
    type(fio_simulation_t),  pointer :: sim
    type(fio_hamiltonian_t), pointer :: hml
    !
    PUSH_SUB(fio_model_start)
    nullify(sys, geom, geo, sim, hml)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_get(sys, geom)
    ASSERT(associated(geom))
    call fio_geom_get(geom, geo)
    ASSERT(associated(geo))
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation_start(sim, geo, mpi_grp)
    nullify(geom, geo)
    call fio_system_start(sys, sim)
    nullify(sys)
    call fio_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian_start(hml, sim)
    nullify(sim, hml)
    POP_SUB(fio_model_start)
    return
  end subroutine fio_model_start

  ! ---------------------------------------------------------
  subroutine fio_model_update(this)
    type(fio_model_t), intent(inout) :: this
    !
    type(fio_system_t),      pointer :: sys
    type(fio_hamiltonian_t), pointer :: hml
    !
    PUSH_SUB(fio_model_update)
    nullify(sys, hml)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_update(sys)
    nullify(sys)
    call fio_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian_update(hml)
    nullify(hml)
    POP_SUB(fio_model_update)
    return
  end subroutine fio_model_update

  ! ---------------------------------------------------------
  subroutine fio_model_stop(this)
    type(fio_model_t), intent(inout) :: this
    !
    type(fio_simulation_t),  pointer :: sim
    type(fio_system_t),      pointer :: sys
    type(fio_hamiltonian_t), pointer :: hml
    !
    PUSH_SUB(fio_model_stop)
    nullify(sim, sys, hml)
    call fio_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_hamiltonian_stop(hml)
    nullify(hml)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_stop(sys)
    nullify(sys)
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation_stop(sim)
    nullify(sim)
    POP_SUB(fio_model_stop)
    return
  end subroutine fio_model_stop

end module fio_model_m

!! Local Variables:
!! mode: f90
!! End:
