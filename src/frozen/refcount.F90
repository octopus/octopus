#include "global.h"

module refcount_oct_m

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::     &
    refcount_t

  public ::        &
    refcount_new,  &
    refcount_del,  &
    refcount_init, &
    refcount_set,  &
    refcount_get,  &
    refcount_inc,  &
    refcount_dec,  &
    refcount_copy, & 
    refcount_end

  integer, parameter :: REFCOUNT_NULL = 0
  integer, parameter :: REFCOUNT_STAT = 1
  integer, parameter :: REFCOUNT_ALLC = 2

  type :: refcount_t
    private
    integer :: type  = REFCOUNT_NULL
    integer :: count = 0
  end type refcount_t

contains

  ! ---------------------------------------------------------
  function refcount_new() result(this)
    type(refcount_t), pointer :: this

    PUSH_SUB(refcount_new)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call refcount_init(this)

    POP_SUB(refcount_new)
  end function refcount_new

  ! ---------------------------------------------------------
  subroutine refcount_del(this)
    type(refcount_t), pointer :: this

    PUSH_SUB(refcount_del)

    if(associated(this))then
      call refcount_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)
    
    POP_SUB(refcount_del)
  end subroutine refcount_del

  ! ---------------------------------------------------------
  subroutine refcount_init(this)
    type(refcount_t), intent(out) :: this

    PUSH_SUB(refcount_init)

    this%type  = REFCOUNT_STAT
    this%count = 0

    POP_SUB(refcount_init)
  end subroutine refcount_init

  ! ---------------------------------------------------------
  subroutine refcount_set(this, static, dynamic)
    type(refcount_t),  intent(inout) :: this
    logical, optional, intent(in)    :: static
    logical, optional, intent(in)    :: dynamic

    logical :: stat, allc

    PUSH_SUB(refcount_set)

    ASSERT(.not.this%count<0)
    stat = .false.
    if(present(static)) stat = static
    allc = .false.
    if(present(dynamic)) allc = dynamic
    if(stat)then
      if(allc)then
        ASSERT(.false.)
      else
        this%type = REFCOUNT_STAT
      end if
    else
      if(allc)then
        this%type = REFCOUNT_ALLC
      else
        ASSERT(.false.)
      end if
    end if

    POP_SUB(refcount_set)
  end subroutine refcount_set

  ! ---------------------------------------------------------
  subroutine refcount_get(this, free)
    type(refcount_t), intent(in)  :: this
    logical,          intent(out) :: free

    PUSH_SUB(refcount_get)

    ASSERT(this%type/=REFCOUNT_NULL)
    ASSERT(.not.this%count<0)
    select case(this%type)
    case(REFCOUNT_STAT)
      free = .false.
    case(REFCOUNT_ALLC)
      free = (this%count==0)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(refcount_get)
  end subroutine refcount_get

  ! ---------------------------------------------------------
  subroutine refcount_inc(this)
    type(refcount_t), intent(inout) :: this

    PUSH_SUB(refcount_inc)

    ASSERT(this%type/=REFCOUNT_NULL)
    this%count = this%count + 1
    ASSERT(this%count>0)

    POP_SUB(refcount_inc)
  end subroutine refcount_inc

  ! ---------------------------------------------------------
  subroutine refcount_dec(this)
    type(refcount_t), intent(inout) :: this

    PUSH_SUB(refcount_dec)

    ASSERT(this%type/=REFCOUNT_NULL)
    ASSERT(this%count>0)
    this%count = this%count - 1

    POP_SUB(refcount_dec)
  end subroutine refcount_dec

  ! ---------------------------------------------------------
  subroutine refcount_copy(this, that)
    type(refcount_t), intent(inout) :: this
    type(refcount_t), intent(in)    :: that

    PUSH_SUB(refcount_copy)

    call refcount_end(this)
    this%type  = that%type
    this%count = that%count

    POP_SUB(refcount_copy)
  end subroutine refcount_copy

  ! ---------------------------------------------------------
  subroutine refcount_end(this)
    type(refcount_t), intent(inout) :: this

    PUSH_SUB(refcount_end)

    ASSERT(this%count==0)
    this%count = 0
    this%type  = REFCOUNT_NULL

    POP_SUB(refcount_end)
  end subroutine refcount_end

end module refcount_oct_m

!! Local Variables:
!! mode: f90
!! End:
