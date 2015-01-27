#include "global.h"

module imesh_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m, only: curvilinear_t
  use json_m,        only: json_object_t
  use simul_box_m,   only: simul_box_t

  use mesh_m, only: &
    mesh_t,         &
    mesh_end

  implicit none

  private
  public ::    &
    mesh_t,    &
    mesh_init, &
    mesh_copy, &
    mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine mesh_init(this, sb, cv, config)
    type(mesh_t),                  intent(out) :: this
    type(simul_box_t),             intent(in)  :: sb
    type(curvilinear_t),           intent(in)  :: cv
    type(json_object_t), optional, intent(in)  :: config
    !
    PUSH_SUB(mesh_init)
    ASSERT(.false.)
    POP_SUB(mesh_init)
    return
  end subroutine mesh_init

  ! ---------------------------------------------------------
  subroutine mesh_copy(this, that)
    type(mesh_t), intent(out) :: this
    type(mesh_t), intent(in)  :: that
    !
    PUSH_SUB(mesh_copy)
    ASSERT(.false.)
    POP_SUB(mesh_copy)
    return
  end subroutine mesh_copy

end module imesh_m

!! Local Variables:
!! mode: f90
!! End:
