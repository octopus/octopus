#include "global.h"

module ssys_external_m
  use base_external_m
  use global_m
  use kinds_m
  use live_external_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::           &
    SSYS_EXTERNAL_OK

  public ::               &
    ssys_external_t,      &
    ssys_external_init,   &
    ssys_external_update, &
    ssys_external_next,   &
    ssys_external_get,    &
    ssys_external_end

  public ::                   &
    ssys_external_iterator_t

  integer, parameter :: SSYS_EXTERNAL_OK = 0

  type :: ssys_external_t
    private
  end type ssys_external_t

  type :: ssys_external_iterator_t
    private
  end type ssys_external_iterator_t

  interface ssys_external_init
    module procedure ssys_external_iterator_init
  end interface ssys_external_init

  interface ssys_external_next
    module procedure ssys_external_iterator_next_name_that
  end interface ssys_external_next

  interface ssys_external_get
    module procedure ssys_external_get_external
    module procedure ssys_external_get_potential_1d
    module procedure ssys_external_get_potential_md
  end interface ssys_external_get

  interface ssys_external_end
    module procedure ssys_external_iterator_end
  end interface ssys_external_end

contains

  ! ---------------------------------------------------------
  subroutine ssys_external_update(this)
    type(ssys_external_t), intent(inout) :: this

    PUSH_SUB(ssys_external_update)
   
    ASSERT(.false.)

    POP_SUB(ssys_external_update)
  end subroutine ssys_external_update

  ! ---------------------------------------------------------
  subroutine ssys_external_get_external(this, name, that)
    type(ssys_external_t),  intent(in) :: this
    character(len=*),       intent(in) :: name
    type(live_external_t), pointer     :: that

    PUSH_SUB(ssys_external_get_external)

    ASSERT(.false.)

    POP_SUB(ssys_external_get_external)
  end subroutine ssys_external_get_external

  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_1d(this, that)
    type(ssys_external_t),        intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(ssys_external_get_potential_1d)

    ASSERT(.false.)

    POP_SUB(ssys_external_get_potential_1d)
  end subroutine ssys_external_get_potential_1d

  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_md(this, that)
    type(ssys_external_t),          intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(ssys_external_get_potential_md)

    ASSERT(.false.)

    POP_SUB(ssys_external_get_potential_md)
  end subroutine ssys_external_get_potential_md

  ! ---------------------------------------------------------
  subroutine ssys_external_iterator_init(this, that)
    type(ssys_external_iterator_t), intent(out) :: this
    type(ssys_external_t),          intent(in)  :: that

    PUSH_SUB(ssys_external_iterator_init)
    
    ASSERT(.false.)

    POP_SUB(ssys_external_iterator_init)
  end subroutine ssys_external_iterator_init

  ! ---------------------------------------------------------
  subroutine ssys_external_iterator_next_name_that(this, name, that, ierr)
    type(ssys_external_iterator_t), intent(inout) :: this
    character(len=*),               intent(out)   :: name
    type(base_external_t),         pointer        :: that
    integer,              optional, intent(out)   :: ierr

    PUSH_SUB(ssys_external_iterator_next_name_that)
    
    ASSERT(.false.)

    POP_SUB(ssys_external_iterator_next_name_that)
  end subroutine ssys_external_iterator_next_name_that

  ! ---------------------------------------------------------
  subroutine ssys_external_iterator_end(this)
    type(ssys_external_iterator_t), intent(inout) :: this

    PUSH_SUB(ssys_external_iterator_end)

    ASSERT(.false.)

    POP_SUB(ssys_external_iterator_end)
  end subroutine ssys_external_iterator_end

end module ssys_external_m

!! Local Variables:
!! mode: f90
!! End:
