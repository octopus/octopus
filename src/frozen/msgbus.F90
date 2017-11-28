#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME chnl
#define LIST_TYPE_NAME channel_t
#define LIST_TYPE_MODULE_NAME channel_oct_m
#include "tlist_inc.F90"
#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#define LIST_TEMPLATE_NAME msgb
#define LIST_TYPE_NAME msgbus_t

module msgbus_oct_m

  use channel_oct_m
  use chnl_list_oct_m
  use envelope_oct_m
  use global_oct_m
  use message_oct_m
  use messages_oct_m
  use profiling_oct_m

#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX

  implicit none

  private

#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER

  public ::             &
    MSGBUS_OK,          &
    MSGBUS_EMPTY_ERROR
    
  public ::   &
    msgbus_t
    
  public ::            &
    msgbus_iterator_t
  
  public ::         &
    msgbus_new,     &
    msgbus_del,     &
    msgbus_init,    &
    msgbus_attach,  &
    msgbus_detach,  &
    msgbus_notify,  &
    msgbus_publish, &
    msgbus_remove,  & 
    msgbus_next,    & 
    msgbus_copy,    & 
    msgbus_end

  integer, parameter :: MSGBUS_OK          = CHANNEL_OK
  integer, parameter :: MSGBUS_EMPTY_ERROR = CHANNEL_EMPTY_ERROR

  type :: msgbus_t
    private
    type(msgb_list_t) :: list
    type(chnl_list_t) :: plst
    type(chnl_list_t) :: slst
  end type msgbus_t

  type :: msgbus_iterator_t
    private
    type(msgbus_t),  pointer :: self => null()
    type(channel_iterator_t) :: iter
  end type msgbus_iterator_t

  interface msgbus_init
    module procedure msgbus_init_type
    module procedure msgbus_iterator_init
  end interface msgbus_init

  interface msgbus_remove
    module procedure msgbus_iterator_remove
  end interface msgbus_remove

  interface msgbus_next
    module procedure msgbus_iterator_next
  end interface msgbus_next

  interface msgbus_copy
    module procedure msgbus_copy_type
    module procedure msgbus_iterator_copy
  end interface msgbus_copy

  interface msgbus_end
    module procedure msgbus_end_type
    module procedure msgbus_iterator_end
  end interface msgbus_end

contains

#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY

  ! ---------------------------------------------------------
  function msgbus_new(number) result(this)
    integer, optional, intent(in) :: number
    
    type(msgbus_t), pointer :: this
    
    PUSH_SUB(msgbus_new)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call msgbus_init(this, number)

    POP_SUB(msgbus_new)
  end function msgbus_new

  ! ---------------------------------------------------------
  subroutine msgbus_del(this)
    type(msgbus_t), pointer, intent(inout) :: this

    PUSH_SUB(msgbus_del)

    if(associated(this))then
      call msgbus_end(this)
      SAFE_DEALLOCATE_P(this)
      nullify(this)
    end if
    nullify(this)
    
    POP_SUB(msgbus_del)
  end subroutine msgbus_del

  ! ---------------------------------------------------------
  subroutine msgbus_init_type(this, number)
    type(msgbus_t),    intent(out) :: this
    integer, optional, intent(in)  :: number

    integer :: ichn, cnum

    PUSH_SUB(msgbus_init_type)

    call msgb_list_init(this%list)
    call chnl_list_init(this%plst)
    call chnl_list_init(this%slst)
    cnum = 1
    if(present(number)) cnum = number
    ASSERT(cnum>0)
    do ichn = cnum, 1, -1
      call chnl_list_push(this%slst, channel_new(ichn))
    end do
    
    POP_SUB(msgbus_init_type)
  end subroutine msgbus_init_type

  ! ---------------------------------------------------------
  recursive subroutine msgbus__attach__(this, that)
    type(msgbus_t), target, intent(in) :: this
    type(channel_t),        intent(in) :: that

    type(msgb_list_iterator_t) :: iter
    type(chnl_list_t), pointer :: plst
    type(msgbus_t),    pointer :: msgb

    PUSH_SUB(msgbus__attach__)

    nullify(plst)
    plst => this%plst
    call chnl_list_push(plst, that)
    nullify(plst)
    call msgb_list_init(iter, this%list)
    do
      nullify(msgb)
      call msgb_list_next(iter, msgb)
      if(.not.associated(msgb)) exit
      call msgbus__attach__(msgb, that)
    end do
    call msgb_list_end(iter)
    nullify(msgb)
    
    POP_SUB(msgbus__attach__)
  end subroutine msgbus__attach__

  ! ---------------------------------------------------------
  subroutine msgbus_attach(this, that, id)
    type(msgbus_t),    intent(inout) :: this
    type(msgbus_t),    intent(in)    :: that
    integer, optional, intent(in)    :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid

    PUSH_SUB(msgbus_attach)

    nullify(chnl)
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(this%slst))
    call chnl_list_get(this%slst, chid, chnl)
    ASSERT(associated(chnl))
    call msgb_list_push(this%list, that)
    call msgbus__attach__(that, chnl)
    nullify(chnl)
    
    POP_SUB(msgbus_attach)
  end subroutine msgbus_attach

  ! ---------------------------------------------------------
  recursive subroutine msgbus__detach__(this, that)
    type(msgbus_t), target, intent(in) :: this
    type(channel_t),        intent(in) :: that

    type(msgb_list_iterator_t) :: iter
    type(chnl_list_t), pointer :: plst
    type(msgbus_t),    pointer :: msgb

    PUSH_SUB(msgbus__detach__)

    nullify(plst)
    plst => this%plst
    call chnl_list_del(plst, that)
    nullify(plst)
    call msgb_list_init(iter, this%list)
    do
      nullify(msgb)
      call msgb_list_next(iter, msgb)
      if(.not.associated(msgb)) exit
      call msgbus__detach__(msgb, that)
    end do
    call msgb_list_end(iter)
    nullify(msgb)

    POP_SUB(msgbus__detach__)
  end subroutine msgbus__detach__

  ! ---------------------------------------------------------
  subroutine msgbus_detach(this, that, id)
    type(msgbus_t),    intent(inout) :: this
    type(msgbus_t),    intent(in)    :: that
    integer, optional, intent(in)    :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid

    PUSH_SUB(msgbus_detach)

    nullify(chnl)
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(this%slst))
    call chnl_list_get(this%slst, chid, chnl)
    ASSERT(associated(chnl))
    call msgb_list_del(this%list, that)
    call msgbus__detach__(that, chnl)
    nullify(chnl)

    POP_SUB(msgbus_detach)
  end subroutine msgbus_detach

  ! ---------------------------------------------------------
  subroutine msgbus_notify(this, message, id)
    type(msgbus_t),    intent(inout) :: this
    type(message_t),   intent(in)    :: message
    integer, optional, intent(in)    :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid

    PUSH_SUB(msgbus_notify)

    nullify(chnl)
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(this%slst))
    call chnl_list_get(this%slst, chid, chnl)
    ASSERT(associated(chnl))
    call channel_send(chnl, envelope_new(message))
    nullify(chnl)

    POP_SUB(msgbus_notify)
  end subroutine msgbus_notify

  ! ---------------------------------------------------------
  subroutine msgbus_publish(this, message)
    type(msgbus_t),  intent(inout) :: this
    type(message_t), intent(in)    :: message

    type(chnl_list_iterator_t) :: iter
    type(channel_t),   pointer :: chnl

    PUSH_SUB(msgbus_publish)

    call chnl_list_init(iter, this%plst)
    do
      nullify(chnl)
      call chnl_list_next(iter, chnl)
      if(.not.associated(chnl)) exit
      call channel_send(chnl, envelope_new(message))
    end do
    call chnl_list_end(iter)
    nullify(chnl)

    POP_SUB(msgbus_publish)
  end subroutine msgbus_publish

  ! ---------------------------------------------------------
  subroutine msgbus_copy_type(this, that)
    type(msgbus_t), intent(inout) :: this
    type(msgbus_t), intent(in)    :: that

    type(msgb_list_iterator_t) :: itrb
    type(chnl_list_iterator_t) :: itrc
    type(msgbus_t),    pointer :: msgb
    type(channel_t),   pointer :: chnl

    PUSH_SUB(msgbus_copy_type)

    call msgbus_end(this)
    call msgb_list_init(this%list)
    call msgb_list_init(itrb, that%list)
    do
      nullify(msgb)
      call msgb_list_next(itrb, msgb)
      if(.not.associated(msgb)) exit
      call msgb_list_push(this%list, msgb)
    end do
    call msgb_list_end(itrb)
    nullify(msgb)
    call chnl_list_init(this%plst)
    call chnl_list_init(itrc, that%plst)
    do
      nullify(chnl)
      call chnl_list_next(itrc, chnl)
      if(.not.associated(chnl)) exit
      call chnl_list_push(this%plst, chnl)
    end do
    call chnl_list_end(itrc)
    nullify(chnl)
    call chnl_list_init(this%slst)
    call chnl_list_init(itrc, that%slst)
    do
      nullify(chnl)
      call chnl_list_next(itrc, chnl)
      if(.not.associated(chnl)) exit
      call chnl_list_push(this%slst, channel_new(chnl))
    end do
    call chnl_list_end(itrc)
    nullify(chnl)
    
    POP_SUB(msgbus_copy_type)
  end subroutine msgbus_copy_type

  ! ---------------------------------------------------------
  subroutine msgbus_end_type(this)
    type(msgbus_t), intent(inout) :: this

    type(channel_t), pointer :: chnl

    PUSH_SUB(msgbus_end_type)

    call msgb_list_end(this%list)
    call chnl_list_end(this%plst)
    do
      nullify(chnl)
      call chnl_list_pop(this%slst, chnl)
      if(.not.associated(chnl)) exit
      call channel_del(chnl)
    end do
    nullify(chnl)
    ASSERT(chnl_list_len(this%slst)==0)
    call chnl_list_end(this%slst)
    
    POP_SUB(msgbus_end_type)
  end subroutine msgbus_end_type

  ! ---------------------------------------------------------
  subroutine msgbus_iterator_init(this, that, id)
    type(msgbus_iterator_t), intent(out) :: this
    type(msgbus_t),  target, intent(in)  :: that
    integer,       optional, intent(in)  :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid

    PUSH_SUB(msgbus_iterator_init)

    nullify(chnl)
    this%self => that
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(that%slst))
    call chnl_list_get(that%slst, chid, chnl)
    ASSERT(associated(chnl))
    call channel_init(this%iter, chnl)
    
    POP_SUB(msgbus_iterator_init)
  end subroutine msgbus_iterator_init

  ! ---------------------------------------------------------
  subroutine msgbus_iterator_remove(this, ierr)
    type(msgbus_iterator_t), intent(inout) :: this
    integer,       optional, intent(out)   :: ierr

    PUSH_SUB(msgbus_iterator_remove)

    call channel_remove(this%iter, ierr)

    POP_SUB(msgbus_iterator_remove)
  end subroutine msgbus_iterator_remove

  ! ---------------------------------------------------------
  subroutine msgbus_iterator_next(this, that, ierr)
    type(msgbus_iterator_t),  intent(inout) :: this
    type(message_t), pointer, intent(out)   :: that
    integer,        optional, intent(out)   :: ierr

    type(envelope_t), pointer :: envl

    PUSH_SUB(msgbus_iterator_next)

    nullify(that, envl)
    call channel_next(this%iter, envl, ierr)
    if(associated(envl))then
      call envelope_get(envl, that)
      ASSERT(associated(that))
    end if
    nullify(envl)

    POP_SUB(msgbus_iterator_next)
  end subroutine msgbus_iterator_next

  ! ---------------------------------------------------------
  subroutine msgbus_iterator_copy(this, that)
    type(msgbus_iterator_t), intent(inout) :: this
    type(msgbus_iterator_t), intent(in)    :: that

    PUSH_SUB(msgbus_iterator_copy)

    this%self => that%self
    call channel_copy(this%iter, that%iter)

    POP_SUB(msgbus_iterator_copy)
  end subroutine msgbus_iterator_copy

  ! ---------------------------------------------------------
  subroutine msgbus_iterator_end(this)
    type(msgbus_iterator_t), intent(inout) :: this

    PUSH_SUB(msgbus_iterator_end)

    nullify(this%self)
    call channel_end(this%iter)

    POP_SUB(msgbus_iterator_end)
  end subroutine msgbus_iterator_end

end module msgbus_oct_m

!! Local Variables:
!! mode: f90
!! End:

