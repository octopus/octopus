#include "global.h"

module fio_simulation_m

  use base_geometry_m
  use fio_grid_m
  use geometry_m
  use global_m
  use grid_m
  use json_m
  use messages_m
  use mpi_m
  use profiling_m
  use simulation_m

  implicit none

  private

  public ::                  &
    fio_simulation__new__,   &
    fio_simulation__del__,   &
    fio_simulation__start__, &
    fio_simulation__stop__

contains

  ! ---------------------------------------------------------
  subroutine fio_simulation__new__(this, grid, geom)
    type(simulation_t),    intent(in) :: this
    type(grid_t),         pointer     :: grid
    type(base_geometry_t), intent(in) :: geom

    type(json_object_t), pointer :: scfg, gcfg
    type(geometry_t),    pointer :: pgeo
    integer                      :: ierr

    PUSH_SUB(fio_simulation__new__)

    nullify(grid, scfg, gcfg, pgeo)
    call simulation_get(this, scfg)
    ASSERT(associated(scfg))
    call base_geometry_get(geom, pgeo)
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
    type(simulation_t), intent(in) :: this

    type(grid_t), pointer :: grid

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
    type(simulation_t), intent(in)    :: this
    type(grid_t),       intent(inout) :: grid
    type(mpi_grp_t),    intent(in)    :: mpi_grp

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
    type(simulation_t), intent(in) :: this

    type(grid_t), pointer :: grid

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
