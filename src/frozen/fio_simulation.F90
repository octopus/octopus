#include "global.h"

module fio_simulation_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use mpi_m,   only: mpi_grp_t
  use space_m, only: space_t

  use bgeom_m, only:       &
    fio_geom_t => bgeom_t

  use fio_grid_m, only: &
    fio_grid_t,         &
    fio_grid_init,      &
    fio_grid_start,     &
    fio_grid_stop,      &
    fio_grid_copy,      &
    fio_grid_end

  use simulation_m, only: &
    simulation_set

  use simulation_m, only:                   &
    fio_simulation_t    => simulation_t,    &
    fio_simulation_init => simulation_init, &
    fio_simulation_get  => simulation_get

  implicit none

  private
  public ::               &
    fio_simulation_t,     &
    fio_simulation_init,  &
    fio_simulation_start, &
    fio_simulation_stop,  &
    fio_simulation_get,   &
    fio_simulation_copy,  &
    fio_simulation_end

contains

  ! ---------------------------------------------------------
  subroutine fio_simulation_start(this, grid, geom, space, mpi_grp)
    type(fio_simulation_t), intent(inout) :: this
    type(fio_grid_t),      pointer        :: grid
    type(fio_geom_t),       intent(in)    :: geom
    type(space_t),          intent(in)    :: space
    type(mpi_grp_t),        intent(in)    :: mpi_grp
    !
    type(json_object_t), pointer :: scfg, gcfg
    integer                      :: ierr
    !
    PUSH_SUB(fio_simulation_start)
    nullify(grid, scfg, gcfg)
    call fio_simulation_get(this, scfg)
    ASSERT(associated(scfg))
    call json_get(scfg, "grid", gcfg, ierr)
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(grid)
    call fio_grid_init(grid, geom, gcfg)
    call fio_grid_start(grid, mpi_grp, gcfg)
    nullify(scfg, gcfg)
    POP_SUB(fio_simulation_start)
    return
  end subroutine fio_simulation_start
  
  ! ---------------------------------------------------------
  subroutine fio_simulation_stop(this)
    type(fio_simulation_t), intent(inout) :: this
    !
    type(fio_grid_t), pointer :: grid
    integer                   :: ierr
    !
    PUSH_SUB(fio_simulation_stop)
    call fio_simulation_get(this, grid)
    ASSERT(associated(grid))
    call fio_grid_stop(grid)
    POP_SUB(fio_simulation_stop)
    return
  end subroutine fio_simulation_stop
  
  ! ---------------------------------------------------------
  subroutine fio_simulation_copy(this, that)
    type(fio_simulation_t), intent(inout) :: this
    type(fio_simulation_t), intent(in)    :: that
    !
    type(fio_grid_t), pointer :: ogrd, igrd
    !
    PUSH_SUB(fio_simulation_copy)
    nullify(ogrd, igrd)
    call fio_simulation_end(this)
    call fio_simulation_get(that, igrd)
    if(associated(igrd))then
      SAFE_ALLOCATE(ogrd)
      call fio_grid_copy(ogrd, igrd)
    end if
    nullify(igrd)
    call simulation_set(this, ogrd)
    nullify(ogrd)
    POP_SUB(fio_simulation_copy)
    return
  end subroutine fio_simulation_copy

  ! ---------------------------------------------------------
  subroutine fio_simulation_end(this)
    type(fio_simulation_t), intent(inout) :: this
    !
    type(fio_grid_t), pointer :: grid
    !
    PUSH_SUB(fio_simulation_end)
    nullify(grid)
    call fio_simulation_get(this, grid)
    if(associated(grid))then
      call fio_grid_end(grid)
      SAFE_DEALLOCATE_P(grid)
    end if
    nullify(grid)
    POP_SUB(fio_simulation_end)
    return
  end subroutine fio_simulation_end

end module fio_simulation_m

!! Local Variables:
!! mode: f90
!! End:
