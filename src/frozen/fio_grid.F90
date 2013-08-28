#include "global.h"

module fio_grid_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m, only: curvilinear_init, curvilinear_end
  use json_m,        only: JSON_OK, json_object_t, json_get
  use kinds_m,       only: wp

  use fio_geometry_m, only: &
    fio_geometry_t

  use fio_mesh_m, only: &
    fio_mesh_init,      &
    fio_mesh_end

  use fio_simul_box_m, only: &
    fio_simul_box_init,      &
    fio_simul_box_end

  use igrid_m, only: &
    grid_init,       &
    grid_end

  use igrid_m, only:            &
    fio_grid_t    => grid_t,    &
    fio_grid_get  => grid_get,  &
    fio_grid_copy => grid_copy

  implicit none

  private
  public ::        &
    fio_grid_t,    &
    fio_grid_init, &
    fio_grid_get,  &
    fio_grid_copy, &
    fio_grid_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_grid_init(this, geo, config)
    type(fio_grid_t), target, intent(out) :: this
    type(fio_geometry_t),     intent(in)  :: geo
    type(json_object_t),      intent(in)  :: config
    !
    real(kind=wp), dimension(MAX_DIM) :: spacing
    type(json_object_t),      pointer :: cnfg
    integer                           :: i, ierr
    !
    PUSH_SUB(fio_grid_init)
    call grid_init(this, geo, config)
    call json_get(config, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_simul_box_init(this%sb, geo, cnfg)
    nullify(cnfg)
    call json_get(config, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(cnfg, "spacing", spacing(1:this%sb%dim), ierr)
    ASSERT(ierr==JSON_OK)
    spacing=(/spacing(1:this%sb%dim),(0.0_wp,i=this%sb%dim+1,MAX_DIM)/)
    call curvilinear_init(this%cv, this%sb, geo, spacing)
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
  subroutine fio_grid_end(this)
    type(fio_grid_t), intent(inout) :: this
    !
    PUSH_SUB(fio_grid_end)
    call fio_mesh_end(this%mesh)
    call curvilinear_end(this%cv)
    call fio_simul_box_end(this%sb)
    call grid_end(this)
    POP_SUB(fio_grid_end)
    return
  end subroutine fio_grid_end

end module fio_grid_m

!! Local Variables:
!! mode: f90
!! End:
