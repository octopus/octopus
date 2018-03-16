#include "global.h"

module message_oct_m

  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use refcount_oct_m
  use uuid_oct_m

  implicit none

  private
  
  public ::    &
    message_t

  public ::         &
    message_new,    &
    message_del,    &
    message_init,   &
    message_attach, &
    message_detach, &
    message_get,    &
    message_copy,   & 
    message_end

  integer, parameter :: MSSG_STAT_NULL = 0
  integer, parameter :: MSSG_STAT_LIVE = 1

  type :: message_t
    private
    integer             :: stat = MSSG_STAT_NULL
    type(refcount_t)    :: rcnt
    type(json_object_t) :: data
  end type message_t

contains

  ! ---------------------------------------------------------
  function message_new() result(this)
    type(message_t), pointer :: this
    
    PUSH_SUB(message_new)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call message_init(this)
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(message_new)
  end function message_new

  ! ---------------------------------------------------------
  subroutine message_del(this)
    type(message_t), pointer, intent(inout) :: this

    logical :: free
    
    PUSH_SUB(message_del)
    
    if(associated(this))then
      ASSERT(this%stat>MSSG_STAT_NULL)
      call refcount_get(this%rcnt, free=free)
      if(free)then
        call message_end(this)
        SAFE_DEALLOCATE_P(this)
        nullify(this)
      end if
    end if

    POP_SUB(message_del)
  end subroutine message_del

  ! ---------------------------------------------------------
  subroutine message_init(this)
    type(message_t), intent(out) :: this

    PUSH_SUB(message_init)

    this%stat = MSSG_STAT_LIVE
    call refcount_init(this%rcnt)
    call json_init(this%data)
    
    POP_SUB(message_init)
  end subroutine message_init

  ! ---------------------------------------------------------
  subroutine message_attach(this, that)
    type(message_t), intent(inout) :: this
    type(uuid_t),    intent(in)    :: that

    PUSH_SUB(message_attach)

    ASSERT(this%stat>MSSG_STAT_NULL)
    call refcount_attach(this%rcnt, that)

    POP_SUB(message_attach)
  end subroutine message_attach

  ! ---------------------------------------------------------
  subroutine message_detach(this, that)
    type(message_t), intent(inout) :: this
    type(uuid_t),    intent(in)    :: that

    PUSH_SUB(message_detach)

    ASSERT(this%stat>MSSG_STAT_NULL)
    call refcount_detach(this%rcnt, that)

    POP_SUB(message_detach)
  end subroutine message_detach

  ! ---------------------------------------------------------
  subroutine message_get(this, data)
    type(message_t),      target, intent(in)  :: this
    type(json_object_t), pointer, intent(out) :: data

    PUSH_SUB(message_get)

    nullify(data)
    if(this%stat==MSSG_STAT_LIVE) data => this%data

    POP_SUB(message_get)
  end subroutine message_get

  ! ---------------------------------------------------------
  subroutine message_copy(this, that)
    type(message_t), intent(inout) :: this
    type(message_t), intent(in)    :: that

    type(json_object_t), pointer :: data

    PUSH_SUB(message_copy)

    call json_end(this%data)
    if(that%stat==MSSG_STAT_LIVE)then
      call message_get(that, data)
      ASSERT(associated(data))
      call json_copy(this%data, data)
    else
      call json_init(this%data)
    end if
    
    POP_SUB(message_copy)
  end subroutine message_copy
  
  ! ---------------------------------------------------------
  subroutine message_end(this)
    type(message_t), intent(inout) :: this

    PUSH_SUB(message_end)

    this%stat = MSSG_STAT_NULL
    call refcount_end(this%rcnt)
    call json_end(this%data)
    
    POP_SUB(message_end)
  end subroutine message_end
  
end module message_oct_m

!! Local Variables:
!! mode: f90
!! End:

