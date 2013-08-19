#include "global.h"

module frozen_grid_m

  use grid_m, only: &
    grid_t

  use grid_m, only:              &
    frozen_grid_end => grid_end

  private
  public :: &
    grid_t

  public ::           &
    frozen_grid_init, &
    frozen_grid_copy, &
    frozen_grid_end

contains

  subroutine frozen_grid_init
    ASSERT(.false.)
    return
  end subroutine frozen_grid_init

  subroutine frozen_grid_copy(this_out, this_in)
    type(grid_t), intent(out) :: this_out
    type(grid_t), intent(in)  :: this_in
    !
    ASSERT(.false.)
    return
  end subroutine frozen_grid_copy

end module frozen_grid_m

!! Local Variables:
!! mode: f90
!! End:
