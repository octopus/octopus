#include "global.h"

module live_density_m
  use global_m
  use kinds_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::           &
    live_density_t,   &
    live_density_get

  type :: live_density_t
    private
  end type live_density_t

  interface live_density_get
    module procedure live_density_get_density_1d
    module procedure live_density_get_density_2d
  end interface live_density_get

contains

  ! ---------------------------------------------------------
  subroutine live_density_get_density_1d(this, that, total)
    type(live_density_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(live_density_get_density_1d)

    ASSERT(.false.)

    POP_SUB(live_density_get_density_1d)
  end subroutine live_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine live_density_get_density_2d(this, that)
    type(live_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_live_density_2d)

    ASSERT(.false.)

    POP_SUB(live_density_get_density_2d)
  end subroutine live_density_get_density_2d

end module live_density_m

!! Local Variables:
!! mode: f90
!! End:
