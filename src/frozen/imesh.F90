#include "global.h"

module imesh_m

  use global_m
  use messages_m
  use profiling_m

  use mesh_m, only: &
    mesh_t,         &
    mesh_end

  implicit none

  private
  public ::    &
    mesh_t,    &
    mesh_copy, &
    mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine mesh_copy(this, that)
    type(mesh_t), intent(out) :: this
    type(mesh_t), intent(in)  :: that
    !
    PUSH_SUB(mesh_copy)
    ASSERT(.false.)
    this=that
    POP_SUB(mesh_copy)
    return
  end subroutine mesh_copy

end module imesh_m

!! Local Variables:
!! mode: f90
!! End:
