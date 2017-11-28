#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME envl
#define LIST_TYPE_NAME envelope_t
#define LIST_TYPE_MODULE_NAME envelope_oct_m
#include "tlist_inc.F90"
#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

module channel_oct_m

  use envelope_oct_m
  use envl_list_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::              &
    CHANNEL_OK,          &
    CHANNEL_EMPTY_ERROR

  public ::    &
    channel_t

  public ::             &
    channel_iterator_t

  public ::         &
    channel_new,    &
    channel_del,    &
    channel_init,   &
    channel_send,   &
    channel_recv,   &
    channel_set,    &
    channel_get,    &
    channel_remove, &
    channel_next,   &
    channel_copy,   & 
    channel_end

  integer, parameter :: CHANNEL_OK          = ENVL_LIST_OK
  integer, parameter :: CHANNEL_EMPTY_ERROR = ENVL_LIST_EMPTY_ERROR

  type :: channel_t
    private
    integer           :: cid = 0
    type(envl_list_t) :: list
  end type channel_t

  type :: channel_iterator_t
    private
    type(channel_t),   pointer :: self => null()
    type(envl_list_iterator_t) :: iter
  end type channel_iterator_t

  interface channel_new
    module procedure channel_new_type
    module procedure channel_new_copy
  end interface channel_new

  interface channel_init
    module procedure channel_init_type
    module procedure channel_iterator_init
  end interface channel_init

  interface channel_remove
    module procedure channel_iterator_remove
  end interface channel_remove

  interface channel_next
    module procedure channel_iterator_next
  end interface channel_next

  interface channel_copy
    module procedure channel_copy_type
    module procedure channel_iterator_copy
  end interface channel_copy

  interface channel_end
    module procedure channel_end_type
    module procedure channel_iterator_end
  end interface channel_end

contains

  ! ---------------------------------------------------------
  function channel_new_type(id) result(this)
    integer, optional, intent(in) :: id
    
    type(channel_t), pointer :: this
    
    PUSH_SUB(channel_new_type)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call channel_init(this, id)

    POP_SUB(channel_new_type)
  end function channel_new_type

  ! ---------------------------------------------------------
  function channel_new_copy(that, id) result(this)
    type(channel_t),   intent(in) :: that
    integer, optional, intent(in) :: id
    
    type(channel_t), pointer :: this
    
    PUSH_SUB(channel_new_type)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call channel_copy(this, that)
    if(present(id)) call channel_set(this, id=id)

    POP_SUB(channel_new_copy)
  end function channel_new_copy

  ! ---------------------------------------------------------
  subroutine channel_del(this)
    type(channel_t), pointer, intent(inout) :: this

    PUSH_SUB(channel_del)
    
    if(associated(this))then
      call channel_end(this)
      SAFE_DEALLOCATE_P(this)
      nullify(this)
    end if

    POP_SUB(channel_del)
  end subroutine channel_del

  ! ---------------------------------------------------------
  subroutine channel_init_type(this, id)
    type(channel_t),   intent(out) :: this
    integer, optional, intent(in)  :: id

    PUSH_SUB(channel_init_type)

    this%cid = 0
    if(present(id)) this%cid = id
    call envl_list_init(this%list)
    
    POP_SUB(channel_init_type)
  end subroutine channel_init_type

  ! ---------------------------------------------------------
  subroutine channel_send(this, envelope)
    type(channel_t),  intent(inout) :: this
    type(envelope_t), intent(in)    :: envelope

    PUSH_SUB(channel_send)

    call envl_list_append(this%list, envelope)

    POP_SUB(channel_send)
  end subroutine channel_send

  ! ---------------------------------------------------------
  subroutine channel_recv(this, envelope)
    type(channel_t),           intent(inout) :: this
    type(envelope_t), pointer, intent(out)   :: envelope

    PUSH_SUB(channel_recv)

    nullify(envelope)
    call envl_list_pop(this%list, envelope)

    POP_SUB(channel_recv)
  end subroutine channel_recv

  ! ---------------------------------------------------------
  subroutine channel_set(this, id)
    type(channel_t), intent(inout) :: this
    integer,         intent(in)    :: id

    PUSH_SUB(channel_set)

    this%cid = id
    
    POP_SUB(channel_set)
  end subroutine channel_set

  ! ---------------------------------------------------------
  subroutine channel_get(this, id, number)
    type(channel_t),   intent(in)  :: this
    integer, optional, intent(out) :: id
    integer, optional, intent(out) :: number

    PUSH_SUB(channel_get)

    if(present(id)) id = this%cid
    if(present(number)) number = envl_list_len(this%list)
    
    POP_SUB(channel_get)
  end subroutine channel_get

  ! ---------------------------------------------------------
  subroutine channel_copy_type(this, that)
    type(channel_t), intent(inout) :: this
    type(channel_t), intent(in)    :: that

    type(envl_list_iterator_t) :: iter
    type(envelope_t),  pointer :: envl
    integer                    :: chid

    PUSH_SUB(channel_copy_type)

    call channel_end(this)
    call channel_get(that, id=chid)
    call channel_init(this, chid)
    call envl_list_init(iter, that%list)
    do
      nullify(envl)
      call envl_list_next(iter, envl)
      if(.not.associated(envl)) exit
      call channel_send(this, envelope_new(envl))
    end do
    nullify(envl)
    
    POP_SUB(channel_copy_type)
  end subroutine channel_copy_type
  
  ! ---------------------------------------------------------
  subroutine channel_end_type(this)
    type(channel_t), intent(inout) :: this

    type(envelope_t), pointer :: envl

    PUSH_SUB(channel_end_type)

    do
      nullify(envl)
      call channel_recv(this, envl)
      if(.not.associated(envl)) exit
      call envelope_del(envl)
    end do
    nullify(envl)
    ASSERT(envl_list_len(this%list)==0)
    call envl_list_end(this%list)
    
    POP_SUB(channel_end_type)
  end subroutine channel_end_type

  ! ---------------------------------------------------------
  subroutine channel_iterator_init(this, that)
    type(channel_iterator_t), intent(out) :: this
    type(channel_t),  target, intent(in)  :: that

    PUSH_SUB(channel_iterator_init)

    this%self => that
    call envl_list_init(this%iter, that%list)
    
    POP_SUB(channel_iterator_init)
  end subroutine channel_iterator_init

  ! ---------------------------------------------------------
  subroutine channel_iterator_remove(this, ierr)
    type(channel_iterator_t), intent(inout) :: this
    integer,        optional, intent(out)   :: ierr

    type(envelope_t), pointer :: envl
    
    PUSH_SUB(channel_iterator_remove)

    nullify(envl)
    call envl_list_remove(this%iter, envl, ierr)
    ASSERT(associated(envl))
    call envelope_del(envl)

    POP_SUB(channel_iterator_remove)
  end subroutine channel_iterator_remove

  ! ---------------------------------------------------------
  subroutine channel_iterator_next(this, that, ierr)
    type(channel_iterator_t),  intent(inout) :: this
    type(envelope_t), pointer, intent(out)   :: that
    integer,         optional, intent(out)   :: ierr

    PUSH_SUB(channel_iterator_next)

    nullify(that)
    call envl_list_next(this%iter, that, ierr)

    POP_SUB(channel_iterator_next)
  end subroutine channel_iterator_next

  ! ---------------------------------------------------------
  subroutine channel_iterator_copy(this, that)
    type(channel_iterator_t), intent(inout) :: this
    type(channel_iterator_t), intent(in)    :: that

    PUSH_SUB(channel_iterator_copy)

    this%self => that%self
    call envl_list_copy(this%iter, that%iter)

    POP_SUB(channel_iterator_copy)
  end subroutine channel_iterator_copy

  ! ---------------------------------------------------------
  subroutine channel_iterator_end(this)
    type(channel_iterator_t), intent(inout) :: this

    PUSH_SUB(channel_iterator_end)

    nullify(this%self)
    call envl_list_end(this%iter)

    POP_SUB(channel_iterator_end)
  end subroutine channel_iterator_end

end module channel_oct_m

!! Local Variables:
!! mode: f90
!! End:

