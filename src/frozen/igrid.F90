#include "global.h"

module igrid_m

  use global_m
  use messages_m
  use profiling_m

  use geometry_m,  only: geometry_t
  use json_m,      only: JSON_OK, json_object_t, json_get
  use kinds_m,     only: wp
  use mesh_m,      only: mesh_t
  use ob_grid_m,   only: ob_grid_init, ob_grid_end
  use simul_box_m, only: simul_box_t

  use grid_m, only: grid_t

  implicit none

  private
  public ::    &
    grid_t,    &
    grid_init, &
    grid_get,  &
    grid_copy, &
    grid_end

  interface grid_get
    module procedure grid_get_simul_box
    module procedure grid_get_mesh
  end interface grid_get

contains
  
  ! ---------------------------------------------------------
  subroutine grid_init(this, geo, config)
    type(grid_t), target, intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(json_object_t),  intent(in)  :: config
    !
    call ob_grid_init(this%ob_grid)
    ASSERT(.not.this%ob_grid%open_boundaries)
    this%have_fine_mesh=.false.
    this%fine%mesh=>this%mesh
    this%fine%der=>this%der
    nullify(this%mgrid)
    return
  end subroutine grid_init

  ! ---------------------------------------------------------
  subroutine grid_get_simul_box(this, that)
    type(grid_t),       target, intent(in) :: this
    type(simul_box_t), pointer             :: that
    !
    that=>this%sb
    return
  end subroutine grid_get_simul_box

  ! ---------------------------------------------------------
  subroutine grid_get_mesh(this, mesh)
    type(grid_t),  target, intent(in) :: this
    type(mesh_t), pointer             :: mesh
    !
    mesh=>this%mesh
    return
  end subroutine grid_get_mesh

  ! ---------------------------------------------------------
  subroutine grid_copy(this, that)
    type(grid_t), intent(out) :: this
    type(grid_t), intent(in)  :: that
    !
    ASSERT(.false.)
    return
  end subroutine grid_copy

  ! ---------------------------------------------------------
  subroutine grid_end(this)
    type(grid_t), intent(inout) :: this
    !
    nullify(this%mgrid, this%fine%der, this%fine%mesh)
    this%have_fine_mesh=.false.
    call ob_grid_end(this%ob_grid)
    return
  end subroutine grid_end

end module igrid_m

!! Local Variables:
!! mode: f90
!! End:
