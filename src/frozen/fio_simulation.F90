#include "global.h"

module fio_simulation_m

  use global_m
  use messages_m
  use profiling_m

  use geometry_m, only: geometry_t
  use json_m,     only: JSON_OK, json_object_t, json_get
  use mpi_m,      only: mpi_grp_t

  use base_geom_m, only: &
    base_geom_t

  use base_geom_m, only: &
    base_geom_get

  use fio_grid_m, only: &
    fio_grid_t

  use fio_grid_m, only: &
    fio_grid_init,      &
    fio_grid_start,     &
    fio_grid_stop,      &
    fio_grid_end

  use simulation_m, only:             &
    fio_simulation_t => simulation_t

  use simulation_m, only: &
    simulation_get

  implicit none

  private
  public ::           &
    fio_simulation_t

  public ::                  &
    fio_simulation__new__,   &
    fio_simulation__del__,   &
    fio_simulation__start__, &
    fio_simulation__stop__

contains

  ! ---------------------------------------------------------
  subroutine fio_simulation__new__(this, grid, geom)
    type(fio_simulation_t), intent(in) :: this
    type(fio_grid_t),      pointer     :: grid
    type(base_geom_t),      intent(in) :: geom

    type(json_object_t), pointer :: scfg, gcfg
    type(geometry_t),    pointer :: pgeo
    integer                      :: ierr

    PUSH_SUB(fio_simulation__new__)

    nullify(grid, scfg, gcfg, pgeo)
    call simulation_get(this, scfg)
    ASSERT(associated(scfg))
    call base_geom_get(geom, pgeo)
    ASSERT(associated(pgeo))
    call json_get(scfg, "grid", gcfg, ierr)
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(grid)
    call fio_grid_init(grid, pgeo, gcfg)
    nullify(scfg, gcfg, pgeo)

    POP_SUB(fio_simulation__new__)
  end subroutine fio_simulation__new__
  
  ! ---------------------------------------------------------
  subroutine fio_simulation__del__(this)
    type(fio_simulation_t), intent(in) :: this

    type(fio_grid_t), pointer :: grid

    PUSH_SUB(fio_simulation__del__)

    nullify(grid)
    call simulation_get(this, grid)
    ASSERT(associated(grid))
    call fio_grid_end(grid)
    SAFE_DEALLOCATE_P(grid)
    nullify(grid)

    POP_SUB(fio_simulation__del__)
  end subroutine fio_simulation__del__
  
  ! ---------------------------------------------------------
  subroutine fio_simulation__start__(this, grid, mpi_grp)
    type(fio_simulation_t), intent(in)    :: this
    type(fio_grid_t),       intent(inout) :: grid
    type(mpi_grp_t),        intent(in)    :: mpi_grp

    type(json_object_t), pointer :: scfg, gcfg
    integer                      :: ierr

    PUSH_SUB(fio_simulation__start__)


    nullify(scfg, gcfg)
    call simulation_get(this, scfg)
    ASSERT(associated(scfg))
    call json_get(scfg, "grid", gcfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_grid_start(grid, mpi_grp, gcfg)
    nullify(scfg, gcfg)

    POP_SUB(fio_simulation__start__)
  end subroutine fio_simulation__start__

  ! ---------------------------------------------------------
  subroutine fio_simulation__stop__(this)
    type(fio_simulation_t), intent(in) :: this

    type(fio_grid_t), pointer :: grid

    PUSH_SUB(fio_simulation__stop__)

    nullify(grid)
    call simulation_get(this, grid)
    ASSERT(associated(grid))
    call fio_grid_stop(grid)
    nullify(grid)

    POP_SUB(fio_simulation__stop__)
  end subroutine fio_simulation__stop__

end module fio_simulation_m

!! Local Variables:
!! mode: f90
!! End:
