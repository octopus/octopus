#include "global.h"

module igrid_m

  use global_m
  use messages_m
  use profiling_m

  use derivatives_m, only: derivatives_t
  use json_m,        only: json_object_t
  use mesh_m,        only: mesh_t
  use simul_box_m,   only: simul_box_t

  use grid_m, only: &
    grid_t,         &
    grid_end

  implicit none

  private
  public ::    &
    grid_t,    &
    grid_init, &
    grid_get,  &
    grid_copy, &
    grid_end

  interface grid_get
    module procedure igrid_get_simul_box
    module procedure igrid_get_mesh
    module procedure igrid_get_der
  end interface grid_get

contains
  
  ! ---------------------------------------------------------
  subroutine grid_init(this, config)
    type(grid_t),                  intent(out) :: this
    type(json_object_t), optional, intent(in)  :: config
    !
    PUSH_SUB(grid_init)
    ASSERT(.false.)
    POP_SUB(grid_init)
    return
  end subroutine grid_init

  ! ---------------------------------------------------------
  subroutine igrid_get_simul_box(this, that)
    type(grid_t),       target, intent(in) :: this
    type(simul_box_t), pointer             :: that
    !
    PUSH_SUB(igrid_get_simul_box)
    that=>this%sb
    POP_SUB(igrid_get_simul_box)
    return
  end subroutine igrid_get_simul_box

  ! ---------------------------------------------------------
  subroutine igrid_get_mesh(this, that)
    type(grid_t),  target, intent(in) :: this
    type(mesh_t), pointer             :: that
    !
    PUSH_SUB(igrid_get_mesh)
    that=>this%mesh
    POP_SUB(igrid_get_mesh)
    return
  end subroutine igrid_get_mesh

  ! ---------------------------------------------------------
  subroutine igrid_get_der(this, that)
    type(grid_t),         target, intent(in) :: this
    type(derivatives_t), pointer             :: that
    !
    PUSH_SUB(igrid_get_der)
    that=>this%der
    POP_SUB(igrid_get_der)
    return
  end subroutine igrid_get_der

  ! ---------------------------------------------------------
  subroutine grid_copy(this, that)
    type(grid_t), intent(inout) :: this
    type(grid_t), intent(in)    :: that
    !
    PUSH_SUB(grid_copy)
    ASSERT(.false.)
    POP_SUB(grid_copy)
    return
  end subroutine grid_copy

end module igrid_m

!! Local Variables:
!! mode: f90
!! End:
