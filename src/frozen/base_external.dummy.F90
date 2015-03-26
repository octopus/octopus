#include "global.h"

module base_external_m
  use global_m
  use kinds_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::                 &
    BASE_EXTERNAL_NAME_LEN

  public ::            &
    base_external_t,   &
    base_external_get

  integer, parameter :: BASE_EXTERNAL_NAME_LEN = 63

  type :: base_external_t
    private
  end type base_external_t

  interface base_external_get
    module procedure base_external_get_potential_1d
    module procedure base_external_get_potential_md
  end interface base_external_get

contains

  ! ---------------------------------------------------------
  subroutine base_external_get_potential_1d(this, that)
    type(base_external_t),        intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_external_get_potential_1d)
    
    ASSERT(.false.)

    POP_SUB(base_external_get_potential_1d)
  end subroutine base_external_get_potential_1d

  ! ---------------------------------------------------------
  subroutine base_external_get_potential_md(this, that)
    type(base_external_t),          intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(base_external_get_potential_md)

    ASSERT(.false.)

    POP_SUB(base_external_get_potential_md)
  end subroutine base_external_get_potential_md

end module base_external_m

!! Local Variables:
!! mode: f90
!! End:
