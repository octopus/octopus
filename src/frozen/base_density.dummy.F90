#include "global.h"

module base_density_m
  use global_m
  use kinds_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::          &
    BASE_DENSITY_OK

  public ::                &
    BASE_DENSITY_NAME_LEN

  public ::           &
    base_density_t,   &
    base_density_get

  integer, parameter :: BASE_DENSITY_OK = 0

  integer, parameter :: BASE_DENSITY_NAME_LEN = 63

  type :: base_density_t
    private
  end type base_density_t

  interface base_density_get
    module procedure base_density_get_info
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
  end interface base_density_get

contains

  ! ---------------------------------------------------------
  subroutine base_density_get_info(this, fine, size, nspin)
    type(base_density_t), intent(in)  :: this
    logical,    optional, intent(out) :: fine
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin

    PUSH_SUB(base_density_get_info)

    ASSERT(.false.)

    POP_SUB(base_density_get_info)
  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, total)
    type(base_density_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(base_density_get_density_1d)

    ASSERT(.false.)

    POP_SUB(base_density_get_density_1d)
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that)
    type(base_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_base_density_2d)

    ASSERT(.false.)

    POP_SUB(base_density_get_density_2d)
  end subroutine base_density_get_density_2d

end module base_density_m

!! Local Variables:
!! mode: f90
!! End:

