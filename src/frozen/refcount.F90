#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#define LIST_TEMPLATE_NAME uuid

module refcount_oct_m

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use uuid_oct_m

#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX

  implicit none

  private

#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER

  public ::     &
    refcount_t

  public ::          &
    refcount_new,    &
    refcount_del,    &
    refcount_init,   &
    refcount_set,    &
    refcount_get,    &
    refcount_attach, &
    refcount_detach, &
    refcount_copy,   & 
    refcount_end

  integer, parameter :: REFCOUNT_NULL = 0
  integer, parameter :: REFCOUNT_STAT = 1
  integer, parameter :: REFCOUNT_ALLC = 2

  type :: refcount_t
    private
    integer           :: type = REFCOUNT_NULL
    type(uuid_list_t) :: list
  end type refcount_t

contains

#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY

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

    this%type = REFCOUNT_STAT
    call uuid_list_init(this%list)

    POP_SUB(refcount_init)
  end subroutine refcount_init

  ! ---------------------------------------------------------
  subroutine refcount_set(this, static, dynamic)
    type(refcount_t),  intent(inout) :: this
    logical, optional, intent(in)    :: static
    logical, optional, intent(in)    :: dynamic

    logical :: stat, allc

    PUSH_SUB(refcount_set)

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(.not.uuid_list_len(this%list)<0)
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

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(.not.uuid_list_len(this%list)<0)
    select case(this%type)
    case(REFCOUNT_STAT)
      free = .false.
    case(REFCOUNT_ALLC)
      free = (uuid_list_len(this%list)==0)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(refcount_get)
  end subroutine refcount_get

  ! ---------------------------------------------------------
  function refcount_find(this, that) result(find)
    type(refcount_t), intent(in) :: this
    type(uuid_t),     intent(in) :: that

    logical :: find

    type(uuid_list_iterator_t) :: iter
    type(uuid_t),      pointer :: uuid
    integer                    :: ierr

    PUSH_SUB(refcount_find)

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(that/=uuid_nil)
    find = .false.
    call uuid_list_init(iter, this%list)
    do
      nullify(uuid)
      call uuid_list_next(iter, uuid, ierr)
      if(ierr/=UUID_LIST_OK)exit
      ASSERT(associated(uuid))
      ASSERT(uuid/=uuid_nil)
      find = (that==uuid)
      if(find)exit
    end do
    call uuid_list_end(iter)
    nullify(uuid)

    POP_SUB(refcount_find)
  end function refcount_find

  ! ---------------------------------------------------------
  subroutine refcount_rem(this, that, ierr)
    type(refcount_t),  intent(inout) :: this
    type(uuid_t),      intent(in)    :: that
    integer, optional, intent(out)   :: ierr

    type(uuid_list_iterator_t) :: iter
    type(uuid_t),      pointer :: uuid
    integer                    :: jerr

    PUSH_SUB(refcount_rem)

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(that/=uuid_nil)
    call uuid_list_init(iter, this%list)
    do
      nullify(uuid)
      call uuid_list_next(iter, uuid, jerr)
      if(jerr/=UUID_LIST_OK)exit
      ASSERT(associated(uuid))
      ASSERT(uuid/=uuid_nil)
      if(that==uuid)exit
    end do
    nullify(uuid)
    if(jerr==UUID_LIST_OK)then
      call uuid_list_remove(iter, uuid, jerr)
      if(jerr==UUID_LIST_OK)then
        ASSERT(associated(uuid))
        ASSERT(uuid/=uuid_nil)
        call uuid_del(uuid)
      end if
      nullify(uuid)
    end if
    if(present(ierr)) ierr = jerr
    call uuid_list_end(iter)

    POP_SUB(refcount_rem)
  end subroutine refcount_rem

  ! ---------------------------------------------------------
  subroutine refcount_attach(this, that)
    type(refcount_t), intent(inout) :: this
    type(uuid_t),     intent(in)    :: that

    PUSH_SUB(refcount_attach)

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(that/=uuid_nil)
    ASSERT(.not.refcount_find(this, that))
    call uuid_list_push(this%list, uuid_new(that))
    ASSERT(uuid_list_len(this%list)>0)

    POP_SUB(refcount_attach)
  end subroutine refcount_attach

  ! ---------------------------------------------------------
  subroutine refcount_detach(this, that)
    type(refcount_t), intent(inout) :: this
    type(uuid_t),     intent(in)    :: that

    integer :: ierr
    
    PUSH_SUB(refcount_detach)

    ASSERT(this%type>REFCOUNT_NULL)
    ASSERT(that/=uuid_nil)
    call refcount_rem(this, that, ierr)
    ASSERT(ierr==UUID_LIST_OK)
    ASSERT(.not.uuid_list_len(this%list)<0)

    POP_SUB(refcount_detach)
  end subroutine refcount_detach

  ! ---------------------------------------------------------
  subroutine refcount_copy(this, that)
    type(refcount_t), intent(inout) :: this
    type(refcount_t), intent(in)    :: that

    type(uuid_list_iterator_t) :: iter
    type(uuid_t),      pointer :: uuid
    integer                    :: ierr

    PUSH_SUB(refcount_copy)

    nullify(uuid)
    call refcount_end(this)
    if(that%type>REFCOUNT_NULL)then
      call uuid_list_init(iter, that%list)
      do
        nullify(uuid)
        call uuid_list_next(iter, uuid, ierr)
        if(ierr/=UUID_LIST_OK)exit
        ASSERT(associated(uuid))
        ASSERT(uuid/=uuid_nil)
        call uuid_list_push(this%list, uuid_new(uuid))
      end do
      call uuid_list_end(iter)
      this%type = that%type
      nullify(uuid)
    end if

    POP_SUB(refcount_copy)
  end subroutine refcount_copy

  ! ---------------------------------------------------------
  subroutine refcount_end(this)
    type(refcount_t), intent(inout) :: this

    PUSH_SUB(refcount_end)

    ASSERT(uuid_list_len(this%list)==0)
    this%type = REFCOUNT_NULL
    call uuid_list_end(this%list)

    POP_SUB(refcount_end)
  end subroutine refcount_end

end module refcount_oct_m

#undef LIST_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
