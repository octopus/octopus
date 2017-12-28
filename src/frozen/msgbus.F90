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
#undef LIST_TYPE_MODULE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

module msgbus_oct_m

  use envelope_oct_m
  use envl_list_oct_m
  use global_oct_m
  use message_oct_m
  use messages_oct_m
  use profiling_oct_m

#define LIST_TEMPLATE_NAME chnl
#define LIST_TYPE_NAME channel_t
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

#define LIST_TEMPLATE_NAME msgb
#define LIST_TYPE_NAME msgbus_t
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

  implicit none

  private

#define LIST_TEMPLATE_NAME chnl
#define LIST_TYPE_NAME channel_t
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

#define LIST_TEMPLATE_NAME msgb
#define LIST_TYPE_NAME msgbus_t
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

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

  integer, parameter :: CHANNEL_OK          = ENVL_LIST_OK
  integer, parameter :: CHANNEL_EMPTY_ERROR = ENVL_LIST_EMPTY_ERROR

  integer, parameter :: MSGBUS_OK          = CHANNEL_OK
  integer, parameter :: MSGBUS_EMPTY_ERROR = CHANNEL_EMPTY_ERROR

  integer, parameter :: CHNL_STAT_NULL = 0
  integer, parameter :: CHNL_STAT_INIT = 1

  type :: channel_t
    private
    integer           :: stat = CHNL_STAT_NULL
    integer           :: chid = 0
    type(msgb_list_t) :: plst
    type(envl_list_t) :: list
  end type channel_t

  type :: msgbus_t
    private
    type(chnl_list_t) :: plst
    type(chnl_list_t) :: slst
  end type msgbus_t

  type :: channel_iterator_t
    private
    type(channel_t),   pointer :: self => null()
    type(envl_list_iterator_t) :: iter
  end type channel_iterator_t

  type :: msgbus_iterator_t
    private
    type(msgbus_t),  pointer :: self => null()
    type(channel_iterator_t) :: iter
  end type msgbus_iterator_t

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

#define LIST_TEMPLATE_NAME chnl
#define LIST_TYPE_NAME channel_t
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

#define LIST_TEMPLATE_NAME msgb
#define LIST_TYPE_NAME msgbus_t
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TYPE_NAME
#undef LIST_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function channel_new_type(id) result(this)
    integer, optional, intent(in) :: id
    
    type(channel_t), pointer :: this
    
    PUSH_SUB(channel_new_type)
    
    nullify(this)
    SAFE_ALLOCATE(this)
    call channel_init(this, id=id)

    POP_SUB(channel_new_type)
  end function channel_new_type

  ! ---------------------------------------------------------
  function channel_new_copy(that, id) result(this)
    type(channel_t),   intent(in) :: that
    integer, optional, intent(in) :: id
    
    type(channel_t), pointer :: this
    
    PUSH_SUB(channel_new_type)
    
    ASSERT(that%stat==CHNL_STAT_INIT)
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

    integer :: chid

    PUSH_SUB(channel_init_type)

    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    this%stat = CHNL_STAT_INIT
    this%chid = 0
    call msgb_list_init(this%plst)
    call envl_list_init(this%list)
    call channel_set(this, id=chid)
    
    POP_SUB(channel_init_type)
  end subroutine channel_init_type

  ! ---------------------------------------------------------
  subroutine channel_send(this, envelope)
    type(channel_t),  intent(inout) :: this
    type(envelope_t), intent(in)    :: envelope

    PUSH_SUB(channel_send)

    ASSERT(this%stat==CHNL_STAT_INIT)
    call envl_list_append(this%list, envelope)

    POP_SUB(channel_send)
  end subroutine channel_send

  ! ---------------------------------------------------------
  subroutine channel_recv(this, envelope)
    type(channel_t),           intent(inout) :: this
    type(envelope_t), pointer, intent(out)   :: envelope

    integer :: ierr

    PUSH_SUB(channel_recv)

    ASSERT(this%stat==CHNL_STAT_INIT)
    nullify(envelope)
    call envl_list_pop(this%list, envelope, ierr)
    if(ierr/=ENVL_LIST_OK) nullify(envelope)

    POP_SUB(channel_recv)
  end subroutine channel_recv

  ! ---------------------------------------------------------
  subroutine channel_set(this, id)
    type(channel_t), intent(inout) :: this
    integer,         intent(in)    :: id

    PUSH_SUB(channel_set)

    ASSERT(this%stat==CHNL_STAT_INIT)
    ASSERT(id>0)
    this%chid = id
    
    POP_SUB(channel_set)
  end subroutine channel_set

  ! ---------------------------------------------------------
  subroutine channel_get(this, id, len)
    type(channel_t),   intent(in)  :: this
    integer, optional, intent(out) :: id
    integer, optional, intent(out) :: len

    PUSH_SUB(channel_get)

    ASSERT(this%stat==CHNL_STAT_INIT)
    if(present(id))  id  = this%chid
    if(present(len)) len = envl_list_len(this%list)
    
    POP_SUB(channel_get)
  end subroutine channel_get

  ! ---------------------------------------------------------
  subroutine channel_copy_type(this, that)
    type(channel_t), intent(inout) :: this
    type(channel_t), intent(in)    :: that

    type(msgb_list_iterator_t) :: itrb
    type(envl_list_iterator_t) :: itrm
    type(msgbus_t),    pointer :: msgb
    type(envelope_t),  pointer :: envl
    integer                    :: chid, ierr

    PUSH_SUB(channel_copy_type)

    nullify(msgb, envl)
    call channel_end(this)
    select case(that%stat)
    case(CHNL_STAT_NULL)
    case(CHNL_STAT_INIT)
      call channel_get(that, id=chid)
      call channel_init(this, chid)
      call envl_list_init(itrm, that%list)
      do
        nullify(envl)
        call envl_list_next(itrm, envl, ierr)
        if(ierr/=ENVL_LIST_OK)exit
        ASSERT(associated(envl))
        call channel_send(this, envelope_new(envl))
      end do
      call envl_list_end(itrm)
      nullify(envl)
      call msgb_list_init(itrb, that%plst)
      do
        nullify(msgb)
        call msgb_list_next(itrb, msgb, ierr)
        if(ierr/=MSGB_LIST_OK)exit
        ASSERT(associated(msgb))
        call msgbus__attach__(msgb, this)
      end do
      call msgb_list_end(itrb)
      nullify(msgb)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(channel_copy_type)
  end subroutine channel_copy_type
  
  ! ---------------------------------------------------------
  subroutine channel_end_type(this)
    type(channel_t), intent(inout) :: this

    type(msgbus_t),   pointer :: msgb
    type(envelope_t), pointer :: envl
    integer                   :: ierr

    PUSH_SUB(channel_end_type)

    do
      nullify(msgb)
      call msgb_list_pop(this%plst, msgb, ierr)
      if(ierr/=MSGB_LIST_OK)exit
      ASSERT(associated(msgb))
      call chnl_list_del(msgb%plst, this, ierr)
      ASSERT(ierr==CHNL_LIST_OK)
    end do
    nullify(msgb)
    do
      nullify(envl)
      call channel_recv(this, envl)
      if(.not.associated(envl)) exit
      call envelope_del(envl)
    end do
    nullify(envl)
    this%stat = CHNL_STAT_NULL
    this%chid = 0
    ASSERT(msgb_list_len(this%plst)==0)
    call msgb_list_end(this%plst)
    ASSERT(envl_list_len(this%list)==0)
    call envl_list_end(this%list)
    
    POP_SUB(channel_end_type)
  end subroutine channel_end_type

  ! ---------------------------------------------------------
  subroutine channel_iterator_init(this, that)
    type(channel_iterator_t), intent(out) :: this
    type(channel_t),  target, intent(in)  :: that

    PUSH_SUB(channel_iterator_init)

    ASSERT(that%stat==CHNL_STAT_INIT)
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

    ASSERT(associated(this%self))
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

    ASSERT(associated(this%self))
    nullify(that)
    call envl_list_next(this%iter, that, ierr)

    POP_SUB(channel_iterator_next)
  end subroutine channel_iterator_next

  ! ---------------------------------------------------------
  subroutine channel_iterator_copy(this, that)
    type(channel_iterator_t), intent(inout) :: this
    type(channel_iterator_t), intent(in)    :: that

    PUSH_SUB(channel_iterator_copy)

    call channel_iterator_end(this)
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
  subroutine msgbus__attach__(this, that)
    type(msgbus_t),  intent(inout) :: this
    type(channel_t), intent(inout) :: that

    PUSH_SUB(msgbus__attach__)

    call chnl_list_push(this%plst, that)
    call msgb_list_push(that%plst, this)
    
    POP_SUB(msgbus__attach__)
  end subroutine msgbus__attach__

  ! ---------------------------------------------------------
  subroutine msgbus_attach(this, that, id)
    type(msgbus_t),    intent(in)    :: this
    type(msgbus_t),    intent(inout) :: that
    integer, optional, intent(in)    :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid, icid

    PUSH_SUB(msgbus_attach)

    nullify(chnl)
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(this%slst))
    call chnl_list_get(this%slst, chid, chnl)
    ASSERT(associated(chnl))
    call channel_get(chnl, id=icid)
    ASSERT(chid==icid)
    call msgbus__attach__(that, chnl)
    nullify(chnl)
    
    POP_SUB(msgbus_attach)
  end subroutine msgbus_attach

  ! ---------------------------------------------------------
  subroutine msgbus__detach__(this, that)
    type(msgbus_t),  intent(inout) :: this
    type(channel_t), intent(inout) :: that

    integer :: ierr

    PUSH_SUB(msgbus__detach__)

    call chnl_list_del(this%plst, that, ierr)
    ASSERT(ierr==CHNL_LIST_OK)
    call msgb_list_del(that%plst, this, ierr)
    ASSERT(ierr==MSGB_LIST_OK)

    POP_SUB(msgbus__detach__)
  end subroutine msgbus__detach__

  ! ---------------------------------------------------------
  subroutine msgbus_detach(this, that, id)
    type(msgbus_t),    intent(in)    :: this
    type(msgbus_t),    intent(inout) :: that
    integer, optional, intent(in)    :: id

    type(channel_t), pointer :: chnl
    integer                  :: chid, icid

    PUSH_SUB(msgbus_detach)

    nullify(chnl)
    chid = 1
    if(present(id)) chid = id
    ASSERT(chid>0)
    ASSERT(chid<=chnl_list_len(this%slst))
    call chnl_list_get(this%slst, chid, chnl)
    ASSERT(associated(chnl))
    call channel_get(chnl, id=icid)
    ASSERT(chid==icid)
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

    type(chnl_list_iterator_t) :: iter
    type(channel_t),   pointer :: chnl
    integer                    :: ierr

    PUSH_SUB(msgbus_copy_type)

    call msgbus_end(this)
    call chnl_list_init(this%plst)
    call chnl_list_init(iter, that%plst)
    do
      nullify(chnl)
      call chnl_list_next(iter, chnl, ierr)
      if(ierr/=CHNL_LIST_OK)exit
      ASSERT(associated(chnl))
      call msgbus__attach__(this, chnl)
    end do
    call chnl_list_end(iter)
    nullify(chnl)
    call chnl_list_init(this%slst)
    call chnl_list_init(iter, that%slst)
    do
      nullify(chnl)
      call chnl_list_next(iter, chnl, ierr)
      if(ierr/=CHNL_LIST_OK)exit
      ASSERT(associated(chnl))
      call chnl_list_append(this%slst, channel_new(chnl))
    end do
    call chnl_list_end(iter)
    nullify(chnl)
    
    POP_SUB(msgbus_copy_type)
  end subroutine msgbus_copy_type

  ! ---------------------------------------------------------
  subroutine msgbus_end_type(this)
    type(msgbus_t), intent(inout) :: this

    type(channel_t), pointer :: chnl
    integer                  :: ierr

    PUSH_SUB(msgbus_end_type)

    do
      nullify(chnl)
      call chnl_list_pop(this%plst, chnl, ierr)
      if(ierr/=CHNL_LIST_OK)exit
      ASSERT(associated(chnl))
      call msgb_list_del(chnl%plst, this, ierr)
      ASSERT(ierr==MSGB_LIST_OK)
    end do
    nullify(chnl)
    ASSERT(chnl_list_len(this%plst)==0)
    call chnl_list_end(this%plst)
    do
      nullify(chnl)
      call chnl_list_pop(this%slst, chnl, ierr)
      if(ierr/=CHNL_LIST_OK)exit
      ASSERT(associated(chnl))
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

    ASSERT(associated(this%self))
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

    ASSERT(associated(this%self))
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

