#include "global.h"

module live_external_m
  use global_m
  use kinds_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::            &
    live_external_t,   &
    live_external_get

  type :: live_external_t
    private
  end type live_external_t

  interface live_external_get
    module procedure live_external_get_potential_1d
    module procedure live_external_get_potential_md
  end interface live_external_get

contains

  ! ---------------------------------------------------------
  subroutine live_external_get_potential_1d(this, that)
    type(live_external_t),        intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(live_external_get_potential_1d)

    ASSERT(.false.)

    POP_SUB(live_external_get_potential_1d)
  end subroutine live_external_get_potential_1d

  ! ---------------------------------------------------------
  subroutine live_external_get_potential_md(this, that)
    type(live_external_t),          intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(live_external_get_potential_md)

    ASSERT(.false.)

    POP_SUB(live_external_get_potential_md)
  end subroutine live_external_get_potential_md

end module live_external_m

!! Local Variables:
!! mode: f90
!! End:
