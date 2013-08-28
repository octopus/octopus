#include "global.h"

module fio_simulation_m

  use global_m
  use messages_m
  use profiling_m

  use basis_m, only: basis_t
  use json_m,  only: JSON_OK, json_object_t, json_get

  use fio_geometry_m, only: fio_geometry_t

  use fio_grid_m, only: fio_grid_t, fio_grid_init

  use simulation_m, only: &
    simulation_start

  use simulation_m, only:                   &
    fio_simulation_t    => simulation_t,    &
    fio_simulation_init => simulation_init, &
    fio_simulation_get  => simulation_get,  &
    fio_simulation_copy => simulation_copy, &
    fio_simulation_end  => simulation_end

  implicit none

  private
  public ::               &
    fio_simulation_t,     &
    fio_simulation_init,  &
    fio_simulation_start, &
    fio_simulation_get,   &
    fio_simulation_copy,  &
    fio_simulation_end
  
contains

  ! ---------------------------------------------------------
  subroutine fio_simulation_start(this)
    type(fio_simulation_t), intent(inout) :: this
    !
    type(json_object_t),  pointer :: config, cnfg
    type(fio_grid_t),     pointer :: gr
    type(fio_geometry_t), pointer :: geo
    integer                       :: ierr
    !
    PUSH_SUB(fio_simulation_start)
    nullify(cnfg, gr, geo)
    call simulation_start(this)
    call fio_simulation_get(this, config)
    ASSERT(associated(config))
    call json_get(config, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_simulation_get(this, gr)
    ASSERT(associated(gr))
    call fio_simulation_get(this, geo)
    ASSERT(associated(geo))
    call fio_grid_init(gr, geo, cnfg)
    nullify(cnfg, gr, geo)
    POP_SUB(fio_simulation_start)
    return
  end subroutine fio_simulation_start

end module fio_simulation_m

!! Local Variables:
!! mode: f90
!! End:
