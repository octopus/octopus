#include "global.h"

module fio_simulation_m

  use fio_grid_m
  use geometry_m
  use global_m
  use grid_m
  use grid_intrf_m
  use json_m
  use messages_m
  use mpi_m
  use profiling_m
  use simulation_m
  use space_m

  implicit none

  private

  public ::                  &
    fio_simulation__build__, &
    fio_simulation__start__, &
    fio_simulation__stop__,  &
    fio_simulation__copy__,  &
    fio_simulation__end__

contains

  ! ---------------------------------------------------------
  subroutine grid__init__(this, geo, space, config)
    type(grid_t),        intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(grid__init__)

    call fio_grid_init(this, geo, space, config)

    POP_SUB(grid__init__)
  end subroutine grid__init__

  ! ---------------------------------------------------------
  subroutine fio_simulation__build__(this)
    type(simulation_t), intent(inout) :: this

    type(grid_intrf_t), pointer :: igrd
    type(grid_t),       pointer :: grid

    PUSH_SUB(fio_simulation__build__)

    nullify(igrd, grid)
    call simulation_get(this, igrd)
    ASSERT(associated(igrd))
    call grid_intrf_new(igrd, grid, grid__init__)
    ASSERT(associated(grid))
    nullify(igrd, grid)

    POP_SUB(fio_simulation__build__)
  end subroutine fio_simulation__build__

  ! ---------------------------------------------------------
  subroutine fio_simulation__start__(this, mpi_grp)
    type(simulation_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: mpi_grp

    type(grid_intrf_t),  pointer :: igrd
    type(json_object_t), pointer :: cnfg
    type(grid_t),        pointer :: grid

    PUSH_SUB(fio_simulation__start__)

    nullify(igrd, cnfg, grid)
    call simulation_get(this, igrd)
    ASSERT(associated(igrd))
    call grid_intrf_get(igrd, cnfg)
    ASSERT(associated(cnfg))
    call grid_intrf_get(igrd, grid)
    ASSERT(associated(grid))
    call fio_grid_start(grid, mpi_grp, cnfg)
    nullify(igrd, cnfg, grid)

    POP_SUB(fio_simulation__start__)
  end subroutine fio_simulation__start__

  ! ---------------------------------------------------------
  subroutine fio_simulation__stop__(this)
    type(simulation_t), intent(inout) :: this

    type(grid_t), pointer :: grid

    PUSH_SUB(fio_simulation__stop__)

    nullify(grid)
    call simulation_get(this, grid)
    ASSERT(associated(grid))
    call fio_grid_stop(grid)
    nullify(grid)

    POP_SUB(fio_simulation__stop__)
  end subroutine fio_simulation__stop__

  ! ---------------------------------------------------------
  subroutine fio_simulation__copy__(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that

    type(grid_intrf_t), pointer :: oigrid, iigrid
    type(grid_t),       pointer :: ogrid, igrid

    PUSH_SUB(fio_simulation__copy__)

    nullify(oigrid, iigrid, ogrid, igrid)
    call simulation_get(that, iigrid)
    ASSERT(associated(iigrid))
    call grid_intrf_get(iigrid, igrid)
    ASSERT(associated(igrid))
    nullify(iigrid)
    call simulation_get(this, oigrid)
    ASSERT(associated(oigrid))
    call grid_intrf_get(oigrid, ogrid)
    ASSERT(associated(ogrid))
    nullify(oigrid)
    call fio_grid_copy(ogrid, igrid)
    nullify(ogrid, igrid)

    POP_SUB(fio_simulation__copy__)
  end subroutine fio_simulation__copy__

  ! ---------------------------------------------------------
  subroutine fio_simulation__end__(this)
    type(simulation_t), intent(inout) :: this

    type(grid_intrf_t), pointer :: igrd

    PUSH_SUB(fio_simulation__end__)

    nullify(igrd)
    call simulation_get(this, igrd)
    ASSERT(associated(igrd))
    call grid_intrf_del(igrd, fio_grid_end)
    nullify(igrd)

    POP_SUB(fio_simulation__end__)
  end subroutine fio_simulation__end__

end module fio_simulation_m

!! Local Variables:
!! mode: f90
!! End:
