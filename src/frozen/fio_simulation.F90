#include "global.h"

module fio_simulation_oct_m

  use fio_grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

  implicit none

  private

  public ::                  &
    fio_simulation__start__, &
    fio_simulation__stop__,  &
    fio_simulation__copy__,  &
    fio_simulation__end__

contains

  ! ---------------------------------------------------------
  subroutine fio_simulation__start__(this, group)
    type(simulation_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: group

    PUSH_SUB(fio_simulation__start__)

    call simulation_start(this, start)

    POP_SUB(fio_simulation__start__)

  contains

    subroutine start(this, geo, space, config)
      type(grid_t),        intent(out) :: this
      type(geometry_t),    intent(in)  :: geo
      type(space_t),       intent(in)  :: space
      type(json_object_t), intent(in)  :: config

      PUSH_SUB(fio_simulation__start__.start)

      call fio_grid_init(this, geo, space, config)
      call fio_grid_start(this, group, config)

      POP_SUB(fio_simulation__start__.start)
    end subroutine start

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

    PUSH_SUB(fio_simulation__copy__)

    call simulation_copy(this, that, copy)

    POP_SUB(fio_simulation__copy__)

  contains

    subroutine copy(this, that)
      type(grid_t), intent(inout) :: this
      type(grid_t), intent(in)    :: that

      PUSH_SUB(fio_simulation__copy__.copy)

      call fio_grid_copy(this, that)

      POP_SUB(fio_simulation__copy__.copy)
    end subroutine copy

  end subroutine fio_simulation__copy__

  ! ---------------------------------------------------------
  subroutine fio_simulation__end__(this)
    type(simulation_t), intent(inout) :: this

    PUSH_SUB(fio_simulation__end__)

    call simulation_end(this, fio_grid_end)

    POP_SUB(fio_simulation__end__)
  end subroutine fio_simulation__end__

end module fio_simulation_oct_m

!! Local Variables:
!! mode: f90
!! End:
