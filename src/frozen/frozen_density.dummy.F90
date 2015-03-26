#include "global.h"

module frozen_density_m
  use global_m
  use kinds_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::           &
    frozen_density_t,   &
    frozen_density_get

  type :: frozen_density_t
    private
  end type frozen_density_t

  interface frozen_density_get
    module procedure frozen_density_get_density_1d
    module procedure frozen_density_get_density_2d
  end interface frozen_density_get

contains

  ! ---------------------------------------------------------
  subroutine frozen_density_get_density_1d(this, that, total)
    type(frozen_density_t),                 intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(frozen_density_get_density_1d)

    ASSERT(.false.)

    POP_SUB(frozen_density_get_density_1d)
  end subroutine frozen_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine frozen_density_get_density_2d(this, that)
    type(frozen_density_t),         intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_frozen_density_2d)

    ASSERT(.false.)

    POP_SUB(frozen_density_get_density_2d)
  end subroutine frozen_density_get_density_2d

end module frozen_density_m
 
!! Local Variables:
!! mode: f90
!! End:
