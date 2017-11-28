#include "global.h"

module envelope_oct_m

  use global_oct_m
  use message_oct_m
  use messages_oct_m
  use profiling_oct_m
  use refcount_oct_m

  implicit none

  private

  public ::     &
    envelope_t

  public ::        &
    envelope_new,  &
    envelope_del,  &
    envelope_init, &
    envelope_get,  &
    envelope_copy, & 
    envelope_end

  integer, parameter :: ENVL_STAT_DISA = 0
  integer, parameter :: ENVL_STAT_NULL = 1
  integer, parameter :: ENVL_STAT_ASSC = 2

  type :: envelope_t
    private
    type(message_t),  pointer :: self =>null()
    type(refcount_t), pointer :: rcnt =>null()
    integer                   :: stat = ENVL_STAT_DISA
  end type envelope_t

  interface envelope_new
    module procedure envelope_new_type
    module procedure envelope_new_copy
  end interface envelope_new

  interface envelope_init
    module procedure envelope_init_type
    module procedure envelope_init_copy
  end interface envelope_init

contains

  ! ---------------------------------------------------------
  function envelope_new_type(that) result(this)
    type(message_t), optional, intent(in) :: that

    type(envelope_t), pointer :: this
    
    PUSH_SUB(envelope_new_type)

    nullify(this)
    SAFE_ALLOCATE(this)
    call envelope_init(this, that)

    POP_SUB(envelope_new_type)
  end function envelope_new_type

  ! ---------------------------------------------------------
  function envelope_new_copy(that) result(this)
    type(envelope_t), intent(in) :: that

    type(envelope_t), pointer :: this
    
    PUSH_SUB(envelope_new_copy)

    nullify(this)
    SAFE_ALLOCATE(this)
    call envelope_init(this, that)

    POP_SUB(envelope_new_copy)
  end function envelope_new_copy

  ! ---------------------------------------------------------
  subroutine envelope_del(this)
    type(envelope_t), pointer, intent(inout) :: this

    PUSH_SUB(envelope_del)

    if(associated(this))then
      call envelope_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(envelope_del)
  end subroutine envelope_del

  ! ---------------------------------------------------------
  function envelope_assoc(this) result(that)
    type(envelope_t), intent(in) :: this

    logical :: that

    PUSH_SUB(envelope_assoc)

    select case(this%stat)
    case(ENVL_STAT_DISA)
      that = .false.
    case(ENVL_STAT_NULL)
      ASSERT(.not.associated(this%self))
      ASSERT(.not.associated(this%rcnt))
      that = .false.
    case(ENVL_STAT_ASSC)
      ASSERT(associated(this%self))
      ASSERT(associated(this%rcnt))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(envelope_assoc)
  end function envelope_assoc

  ! ---------------------------------------------------------
  subroutine envelope_init_type(this, that)
    type(envelope_t),          intent(out) :: this
    type(message_t), optional, intent(in)  :: that

    PUSH_SUB(envelope_init_type)

    nullify(this%self, this%rcnt)
    this%stat = ENVL_STAT_NULL
    if(present(that))then
      call envelope_set(this, that)
    else
      call envelope_set(this, message_new())
    end if

    POP_SUB(envelope_init_type)
  end subroutine envelope_init_type

  ! ---------------------------------------------------------
  subroutine envelope_init_copy(this, that)
    type(envelope_t), intent(out) :: this
    type(envelope_t), intent(in)  :: that

    type(message_t), pointer :: mssg

    PUSH_SUB(envelope_init_copy)

    nullify(mssg)
    call envelope_get(that, mssg)
    ASSERT(associated(mssg))
    call envelope_init(this, mssg)
    nullify(mssg)

    POP_SUB(envelope_init_copy)
  end subroutine envelope_init_copy

  ! ---------------------------------------------------------
  subroutine envelope_set(this, that)
    type(envelope_t),        intent(inout) :: this
    type(message_t), target, intent(in)    :: that

    PUSH_SUB(envelope_set)

    ASSERT(this%stat==ENVL_STAT_NULL)
    this%self => that
    call message_reg(that, this%rcnt)
    ASSERT(associated(this%rcnt))
    call refcount_inc(this%rcnt)
    this%stat = ENVL_STAT_ASSC

    POP_SUB(envelope_set)
  end subroutine envelope_set

  ! ---------------------------------------------------------
  subroutine envelope_get(this, that)
    type(envelope_t),         intent(in)  :: this
    type(message_t), pointer, intent(out) :: that

    PUSH_SUB(envelope_get)

    ASSERT(this%stat>ENVL_STAT_DISA)
    nullify(that)
    if(envelope_assoc(this)) that => this%self

    POP_SUB(envelope_get)
  end subroutine envelope_get

  ! ---------------------------------------------------------
  subroutine envelope_copy(this, that)
    type(envelope_t), intent(inout) :: this
    type(envelope_t), intent(in)    :: that

    PUSH_SUB(envelope_copy)

    call envelope_end(this)
    this%stat = that%stat
    if(envelope_assoc(that)) call envelope_set(this, that%self)

    POP_SUB(envelope_copy)
  end subroutine envelope_copy

  ! ---------------------------------------------------------
  subroutine envelope_end(this)
    type(envelope_t), intent(inout) :: this

    PUSH_SUB(envelope_end)

    if(envelope_assoc(this)) call refcount_dec(this%rcnt)
    call message_del(this%self)
    nullify(this%self, this%rcnt)
    this%stat = ENVL_STAT_DISA

    POP_SUB(envelope_end)
  end subroutine envelope_end

end module envelope_oct_m

!! Local Variables:
!! mode: f90
!! End:

