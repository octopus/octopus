#include "global.h"

module ssys_density_m
  use base_density_m
  use frozen_density_m
  use global_m
  use kinds_m
  use live_density_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::          &
    SSYS_DENSITY_OK

  public ::              &
    ssys_density_t,      &
    ssys_density_init,   &
    ssys_density_next,   &
    ssys_density_get,    &
    ssys_density_end

  public ::                  &
    ssys_density_iterator_t

  integer, parameter :: SSYS_DENSITY_OK = 0

  type :: ssys_density_t
    private
  end type ssys_density_t

  type :: ssys_density_iterator_t
    private
  end type ssys_density_iterator_t

  interface ssys_density_init
    module procedure ssys_density_iterator_init
  end interface ssys_density_init

  interface ssys_density_next
    module procedure ssys_density_iterator_next_name_that
    module procedure ssys_density_iterator_next_that
  end interface ssys_density_next

  interface ssys_density_get
    module procedure ssys_density_get_live_density
    module procedure ssys_density_get_frozen_density
    module procedure ssys_density_get_density_1d
    module procedure ssys_density_get_density_2d
  end interface ssys_density_get

  interface ssys_density_end
    module procedure ssys_density_iterator_end
  end interface ssys_density_end

contains

  ! ---------------------------------------------------------
  subroutine ssys_density_get_live_density(this, name, that)
    type(ssys_density_t),  intent(in) :: this
    character(len=*),      intent(in) :: name
    type(live_density_t), pointer     :: that

    PUSH_SUB(ssys_density_get_live_density)

    ASSERT(.false.)

    POP_SUB(ssys_density_get_live_density)
  end subroutine ssys_density_get_live_density

  ! ---------------------------------------------------------
  subroutine ssys_density_get_frozen_density(this, name, that)
    type(ssys_density_t),    intent(in) :: this
    character(len=*),        intent(in) :: name
    type(frozen_density_t), pointer     :: that

    PUSH_SUB(ssys_density_get_frozen_density)

    ASSERT(.false.)

    POP_SUB(ssys_density_get_frozen_density)
  end subroutine ssys_density_get_frozen_density

  ! ---------------------------------------------------------
  subroutine ssys_density_get_density_1d(this, that, total)
    type(ssys_density_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(ssys_density_get_density_1d)

    ASSERT(.false.)

    POP_SUB(ssys_density_get_density_1d)
  end subroutine ssys_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine ssys_density_get_density_2d(this, that)
    type(ssys_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_ssys_density_2d)

    ASSERT(.false.)

    POP_SUB(ssys_density_get_density_2d)
  end subroutine ssys_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine ssys_density_iterator_init(this, that)
    type(ssys_density_iterator_t), intent(out) :: this
    type(ssys_density_t),          intent(in)  :: that

    PUSH_SUB(ssys_density_iterator_init)
    
    ASSERT(.false.)

    POP_SUB(ssys_density_iterator_init)
  end subroutine ssys_density_iterator_init

  ! ---------------------------------------------------------
  subroutine ssys_density_iterator_next_that(this, that, ierr)
    type(ssys_density_iterator_t), intent(inout) :: this
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)   :: ierr

    PUSH_SUB(ssys_density_iterator_next_that)
    
    ASSERT(.false.)

    POP_SUB(ssys_density_iterator_next_that)
  end subroutine ssys_density_iterator_next_that

  ! ---------------------------------------------------------
  subroutine ssys_density_iterator_next_name_that(this, name, that, ierr)
    type(ssys_density_iterator_t), intent(inout) :: this
    character(len=*),              intent(out)   :: name
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)   :: ierr

    PUSH_SUB(ssys_density_iterator_next_name_that)
    
    ASSERT(.false.)

    POP_SUB(ssys_density_iterator_next_name_that)
  end subroutine ssys_density_iterator_next_name_that

  ! ---------------------------------------------------------
  subroutine ssys_density_iterator_end(this)
    type(ssys_density_iterator_t), intent(inout) :: this

    PUSH_SUB(ssys_density_iterator_end)

    ASSERT(.false.)

    POP_SUB(ssys_density_iterator_end)
  end subroutine ssys_density_iterator_end

end module ssys_density_m

!! Local Variables:
!! mode: f90
!! End:
