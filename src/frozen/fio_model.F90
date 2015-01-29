#include "global.h"

module fio_model_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use mpi_m,   only: mpi_grp_t
  use space_m, only: space_t

  use fio_grid_m, only: &
    fio_grid_t

  use bgeom_m, only:           &
    fio_geom_t   => bgeom_t,   &
    fio_geom_get => bgeom_get

  use fio_simulation_m, only: &
    fio_simulation_t,         &
    fio_simulation_start,     &
    fio_simulation_stop,      &
    fio_simulation_copy,      &
    fio_simulation_end

  use fio_system_m, only: &
    fio_system_t,         &
    fio_system_get,       &
    fio_system_start,     &
    fio_system_update,    &
    fio_system_stop,      &
    fio_system_copy,      &
    fio_system_end

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t,         &
    fio_hamiltonian_init,      &
    fio_hamiltonian_start,     &
    fio_hamiltonian_update,    &
    fio_hamiltonian_stop,      &
    fio_hamiltonian_copy,      &
    fio_hamiltonian_end

  use bmodl_m, only:              &
    fio_model_t    => bmodl_t,    &
    fio_model_get  => bmodl_get

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
  subroutine fio_model_init(this, config)
    type(fio_model_t),   intent(inout) :: this
    type(json_object_t), intent(in)    :: config
    !
    type(json_object_t),     pointer :: cnfg
    type(fio_hamiltonian_t), pointer :: hml
    type(fio_system_t),      pointer :: sys
    integer                          :: ierr
    !
    PUSH_SUB(fio_model_init)
    print *, "enter fio_model_init"
    nullify(cnfg, hml, sys)
    call fio_model_get(this, hml)
    ASSERT(associated(hml))
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call json_get(config, "hamiltonian", cnfg, ierr)
    if(ierr==JSON_OK)call fio_hamiltonian_init(hml, sys, cnfg)
    nullify(cnfg, hml, sys)
    print *, "exit fio_model_init"
    POP_SUB(fio_model_init)
    return
  end subroutine fio_model_init

  ! ---------------------------------------------------------
  subroutine fio_model_start(this, grid, mpi_grp)
    type(fio_model_t), intent(inout) :: this
    type(fio_grid_t), pointer        :: grid
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    !
    type(fio_simulation_t), pointer :: sim
    type(fio_system_t),     pointer :: sys
    type(fio_geom_t),       pointer :: geom
    type(space_t),          pointer :: space
    !
    PUSH_SUB(fio_model_start)
    nullify(sim, sys, geom, space)
    call fio_model_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_get(sys, space)
    ASSERT(associated(space))
    call fio_system_get(sys, geom)
    ASSERT(associated(geom))
    nullify(sys)
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation_start(sim, grid, geom, space, mpi_grp)
    nullify(geom, space)
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
    nullify(sim)
    call fio_model_get(this, sim)
    ASSERT(associated(sim))
    call fio_simulation_stop(sim)
    nullify(sim)
    POP_SUB(fio_model_stop)
    return
  end subroutine fio_model_stop

  ! ---------------------------------------------------------
  subroutine fio_model_copy(this, that)
    type(fio_model_t), intent(inout) :: this
    type(fio_model_t), intent(in)    :: that
    !
    type(fio_simulation_t),  pointer :: osim, isim
    type(fio_hamiltonian_t), pointer :: ohml, ihml
    !
    PUSH_SUB(fio_model_copy)
    nullify(osim, isim, ohml, ihml)
    call fio_model_get(this, osim)
    ASSERT(associated(osim))
    call fio_model_get(that, isim)
    ASSERT(associated(isim))
    call fio_simulation_copy(osim, isim)
    nullify(osim, isim)
    call fio_model_get(this, ohml)
    ASSERT(associated(ohml))
    call fio_model_get(that, ihml)
    ASSERT(associated(ihml))
    call fio_hamiltonian_copy(ohml, ihml)
    nullify(ohml, ihml)
    POP_SUB(fio_model_copy)
    return
  end subroutine fio_model_copy

  ! ---------------------------------------------------------
  subroutine fio_model_end(this)
    type(fio_model_t), intent(inout) :: this
    !
    type(fio_simulation_t),  pointer :: sim
    type(fio_hamiltonian_t), pointer :: hml
    !
    PUSH_SUB(fio_model_end)
    nullify(sim, hml)
    call fio_model_get(this, sim)
    if(associated(sim))then
      call fio_simulation_end(sim)
      nullify(sim)
    end if
    call fio_model_get(this, hml)
    if(associated(hml))then
      call fio_hamiltonian_end(hml)
      nullify(hml)
    end if
    POP_SUB(fio_model_end)
    return
  end subroutine fio_model_end

end module fio_model_m

!! Local Variables:
!! mode: f90
!! End:
