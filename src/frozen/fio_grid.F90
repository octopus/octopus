#include "global.h"

module fio_grid_oct_m

  use boundaries_oct_m
  use curvilinear_oct_m
  use derivatives_oct_m
  use double_grid_oct_m
  use fio_curvilinear_oct_m
  use fio_mesh_oct_m
  use fio_simul_box_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use multigrid_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use stencil_oct_m
  use transfer_table_oct_m

  implicit none

  private

  public ::         &
    fio_grid_init,  &
    fio_grid_start, &
    fio_grid_stop,  &
    fio_grid_copy,  &
    fio_grid_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_grid_init(this, geo, space, config)
    type(grid_t), target, intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(space_t),        intent(in)  :: space
    type(json_object_t),  intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_grid_init)

    nullify(cnfg)
    call json_get(config, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_simul_box_init(this%sb, geo, space, cnfg)
    nullify(cnfg)
    call json_get(config, "curvilinear", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_curvilinear_init(this%cv, this%sb, geo, cnfg)
    nullify(cnfg)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    this%fine%mesh => this%mesh
    this%fine%der => this%der
    nullify(this%mgrid)

    POP_SUB(fio_grid_init)
  end subroutine fio_grid_init

  ! ---------------------------------------------------------
  subroutine fio_grid_start(this, mpi_grp, config)
    type(grid_t),        intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: mpi_grp
    type(json_object_t), intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_grid_start)

    nullify(cnfg)
    call json_get(config, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_mesh_init(this%mesh, this%sb, this%cv, mpi_grp, cnfg)
    nullify(cnfg)

    POP_SUB(fio_grid_start)
  end subroutine fio_grid_start

  ! ---------------------------------------------------------
  subroutine fio_grid_stop(this)
    type(grid_t),    intent(inout) :: this

    PUSH_SUB(fio_grid_stop)

    call fio_mesh_end(this%mesh)

    POP_SUB(fio_grid_stop)
  end subroutine fio_grid_stop

  ! ---------------------------------------------------------
  subroutine fio_grid_copy(this, that)
    type(grid_t), target, intent(inout) :: this
    type(grid_t),         intent(in)    :: that

    PUSH_SUB(fio_grid_copy)

    ASSERT(.not.that%have_fine_mesh)
    call fio_simul_box_copy(this%sb, that%sb)
    call fio_curvilinear_copy(this%cv, that%cv)
    call fio_mesh_copy(this%mesh, that%mesh)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    this%fine%mesh => this%mesh
    this%fine%der => this%der
    nullify(this%mgrid)

    POP_SUB(fio_grid_copy)
  end subroutine fio_grid_copy

  !-------------------------------------------------------------------
  subroutine fio_grid_end(this)
    type(grid_t), intent(inout) :: this

    PUSH_SUB(fio_grid_end)

    call fio_mesh_end(this%mesh)
    call fio_curvilinear_end(this%cv)
    call fio_simul_box_end(this%sb)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    nullify(this%fine%mesh, this%fine%der, this%mgrid)

    POP_SUB(fio_grid_end)
  end subroutine fio_grid_end

end module fio_grid_oct_m

!! Local Variables:
!! mode: f90
!! End:
