#include "global.h"

module fio_grid_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp
  use mpi_m,   only: mpi_grp_t

  use igrid_m, only:          &
    fio_grid_t   => grid_t,   &
    fio_grid_get => grid_get

  use base_geom_m, only:       &
    fio_geom_t => base_geom_t

  use fio_simul_box_m, only: &
    fio_simul_box_t,         &
    fio_simul_box_init,      &
    fio_simul_box_copy,      &
    fio_simul_box_end

  use fio_mesh_m, only: &
    fio_mesh_t,         &
    fio_mesh_init,      &
    fio_mesh_start,     &
    fio_mesh_stop,      &
    fio_mesh_copy,      &
    fio_mesh_end

  use fio_curvilinear_m, only: &
    fio_curvilinear_t,         &
    fio_curvilinear_init,      &
    fio_curvilinear_copy,      &
    fio_curvilinear_end

  implicit none

  private
  public ::         &
    fio_grid_t,     &
    fio_grid_init,  &
    fio_grid_start, &
    fio_grid_stop,  &
    fio_grid_get,   &
    fio_grid_copy,  &
    fio_grid_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_grid_init(this, geom, config)
    type(fio_grid_t), target, intent(out) :: this
    type(fio_geom_t),         intent(in)  :: geom
    type(json_object_t),      intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(fio_grid_init)
    nullify(cnfg)
    call json_get(config, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_simul_box_init(this%sb, geom, cnfg)
    nullify(cnfg)
    call json_get(config, "curvilinear", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_curvilinear_init(this%cv, this%sb, geom, cnfg)
    nullify(cnfg)
    call json_get(config, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_mesh_init(this%mesh, this%sb, this%cv, cnfg)
    nullify(cnfg)
    this%have_fine_mesh=.false.
    this%fine%mesh=>this%mesh
    this%fine%der=>this%der
    nullify(this%mgrid)
    POP_SUB(fio_grid_init)
    return
  end subroutine fio_grid_init

  ! ---------------------------------------------------------
  subroutine fio_grid_start(this, mpi_grp, config)
    type(fio_grid_t),    intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: mpi_grp
    type(json_object_t), intent(in)    :: config
    !
    type(json_object_t),  pointer :: cnfg
    integer                       :: ierr
    !
    PUSH_SUB(fio_grid_start)
    nullify(cnfg)
    call json_get(config, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_mesh_start(this%mesh, mpi_grp, cnfg)
    nullify(cnfg)
    POP_SUB(fio_grid_start)
    return
  end subroutine fio_grid_start

  ! ---------------------------------------------------------
  subroutine fio_grid_stop(this)
    type(fio_grid_t), intent(inout) :: this
    !
    PUSH_SUB(fio_grid_stop)
    call fio_mesh_stop(this%mesh)
    POP_SUB(fio_grid_stop)
    return
  end subroutine fio_grid_stop

  ! ---------------------------------------------------------
  subroutine fio_grid_copy(this, that)
    type(fio_grid_t), target, intent(out) :: this
    type(fio_grid_t),         intent(in)  :: that
    !
    PUSH_SUB(fio_grid_copy)
    call fio_simul_box_copy(this%sb, that%sb)
    call fio_mesh_copy(this%mesh, that%mesh)
    call fio_curvilinear_copy(this%cv, that%cv)
    this%have_fine_mesh=.false.
    this%fine%mesh=>this%mesh
    this%fine%der=>this%der
    nullify(this%mgrid)
    POP_SUB(fio_grid_copy)
    return
  end subroutine fio_grid_copy

  ! ---------------------------------------------------------
  subroutine fio_grid_end(this)
    type(fio_grid_t), intent(inout) :: this
    !
    PUSH_SUB(fio_grid_end)
    nullify(this%mgrid, this%fine%der, this%fine%mesh)
    this%have_fine_mesh=.false.
    call fio_mesh_end(this%mesh)
    call fio_curvilinear_end(this%cv)
    call fio_simul_box_end(this%sb)
    POP_SUB(fio_grid_end)
    return
  end subroutine fio_grid_end

end module fio_grid_m

!! Local Variables:
!! mode: f90
!! End:
