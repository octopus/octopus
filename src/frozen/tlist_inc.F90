#include "global.h"
#include "util.h"

!LIST: LIST_TEMPLATE_NAME
!LIST: LIST_TYPE_NAME
!LIST: LIST_TYPE_MODULE_NAME

#if defined(LIST_TEMPLATE_NAME)
#if !defined(LIST_TYPE_NAME)
#define LIST_TYPE_NAME DECORATE(LIST_TEMPLATE_NAME,t)
#endif
#if !defined(LIST_TYPE_MODULE_NAME)
#define LIST_TYPE_MODULE_NAME DECORATE(LIST_TEMPLATE_NAME,m)
#endif
#else
#error "'LIST_TEMPLATE_NAME' must be defined"
#endif

#undef LIST_INCLUDE_MODULE
#if !defined(LIST_INCLUDE_PREFIX) && !defined(LIST_INCLUDE_HEADER) && !defined(LIST_INCLUDE_BODY)
#define LIST_INCLUDE_MODULE
#endif

#if defined(LIST_INCLUDE_PREFIX) && defined(LIST_INCLUDE_HEADER)
#error "Only one off 'LIST_INCLUDE_PREFIX' or 'LIST_INCLUDE_HEADER' can be defined."
#endif

#if defined(LIST_INCLUDE_PREFIX) && defined(LIST_INCLUDE_BODY)
#error "Only one off 'LIST_INCLUDE_PREFIX' or 'LIST_INCLUDE_BODY' can be defined."
#endif

#if defined(LIST_INCLUDE_HEADER) && defined(LIST_INCLUDE_BODY)
#error "Only one off 'LIST_INCLUDE_HEADER' or 'LIST_INCLUDE_BODY' can be defined."
#endif

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX LIST_TEMPLATE_NAME
#include "template.h"

#if defined(LIST_INCLUDE_MODULE)

module TEMPLATE(list_m)

  use global_m
  use messages_m
  use profiling_m

  use LIST_TYPE_MODULE_NAME

  implicit none

  private

  public ::                     &
    TEMPLATE(LIST_OK),          &
    TEMPLATE(LIST_EMPTY_ERROR), &
    TEMPLATE(LIST_INDEX_ERROR)

  public ::           &
    TEMPLATE(list_t)

  public ::                    &
    TEMPLATE(list_iterator_t)

  public ::                &
    TEMPLATE(list_len),    &
    TEMPLATE(list_init),   &
    TEMPLATE(list_next),   &
    TEMPLATE(list_index),  &
    TEMPLATE(list_set),    &
    TEMPLATE(list_get),    &
    TEMPLATE(list_del),    &
    TEMPLATE(list_push),   &
    TEMPLATE(list_pop),    &
    TEMPLATE(list_append), &
    TEMPLATE(list_extend), &
    TEMPLATE(list_copy),   &
    TEMPLATE(list_end)

#endif
#if defined(LIST_INCLUDE_HEADER) || defined(LIST_INCLUDE_MODULE)

  integer, parameter :: TEMPLATE(LIST_OK)          = 0
  integer, parameter :: TEMPLATE(LIST_EMPTY_ERROR) =-1
  integer, parameter :: TEMPLATE(LIST_INDEX_ERROR) =-2

  type :: INTERNAL(node_t)
    private
    type(LIST_TYPE_NAME),   pointer :: data =>null()
    type(INTERNAL(node_t)), pointer :: next =>null()
  end type INTERNAL(node_t)

  type :: TEMPLATE(list_t)
    private
    integer                         :: size = 0
    type(INTERNAL(node_t)), pointer :: head =>null()
    type(INTERNAL(node_t)), pointer :: tail =>null()
  end type TEMPLATE(list_t)

  type :: TEMPLATE(list_iterator_t)
    private
    type(INTERNAL(node_t)), pointer :: prev =>null()
    type(INTERNAL(node_t)), pointer :: node =>null()
    type(INTERNAL(node_t)), pointer :: next =>null()
  end type TEMPLATE(list_iterator_t)

  interface TEMPLATE(list_init)
    module procedure INTERNAL(list_init)
    module procedure INTERNAL(list_iterator_init_list)
    module procedure INTERNAL(list_iterator_init_iterator)
  end interface TEMPLATE(list_init)

  interface TEMPLATE(list_next)
    module procedure INTERNAL(list_iterator_next)
  end interface TEMPLATE(list_next)

  interface TEMPLATE(list_del)
    module procedure INTERNAL(list_del_data)
    module procedure INTERNAL(list_del_index)
  end interface TEMPLATE(list_del)

  interface TEMPLATE(list_copy)
    module procedure INTERNAL(list_copy)
    module procedure INTERNAL(list_iterator_copy)
  end interface TEMPLATE(list_copy)

  interface TEMPLATE(list_end)
    module procedure INTERNAL(list_end)
    module procedure INTERNAL(list_iterator_end)
  end interface TEMPLATE(list_end)

#endif
#if defined(LIST_INCLUDE_MODULE)

contains
  
#endif
#if defined(LIST_INCLUDE_BODY) || defined(LIST_INCLUDE_MODULE)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_new)(this, that)
    type(INTERNAL(node_t)), pointer     :: this
    type(LIST_TYPE_NAME),    intent(in) :: that

    PUSH_SUB(INTERNAL(node_new))

    nullify(this)
    SAFE_ALLOCATE(this)
    call INTERNAL(node_init)(this, that)

    POP_SUB(INTERNAL(node_new))
  end subroutine INTERNAL(node_new)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_del)(this)
    type(INTERNAL(node_t)), pointer :: this

    PUSH_SUB(INTERNAL(node_del))

    call INTERNAL(node_end)(this)
    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(INTERNAL(node_del))
  end subroutine INTERNAL(node_del)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_init)(this, that)
    type(INTERNAL(node_t)), intent(out) :: this
    type(LIST_TYPE_NAME),   intent(in)  :: that

    PUSH_SUB(INTERNAL(node_init))

    call INTERNAL(node_set)(this, that)
    nullify(this%next)

    POP_SUB(INTERNAL(node_init))
  end subroutine INTERNAL(node_init)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_set)(this, that)
    type(INTERNAL(node_t)),       intent(inout) :: this
    type(LIST_TYPE_NAME), target, intent(in)    :: that

    PUSH_SUB(INTERNAL(node_set))

    this%data => that

    POP_SUB(INTERNAL(node_set))
  end subroutine INTERNAL(node_set)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_get)(this, that)
    type(INTERNAL(node_t)), intent(in) :: this
    type(LIST_TYPE_NAME),  pointer     :: that

    PUSH_SUB(INTERNAL(node_get))

    nullify(that)
    if(associated(this%data)) that => this%data

    POP_SUB(INTERNAL(node_get))
  end subroutine INTERNAL(node_get)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_copy)(this, that)
    type(INTERNAL(node_t)),         intent(inout) :: this
    type(INTERNAL(node_t)), target, intent(in)    :: that

    PUSH_SUB(INTERNAL(node_copy))

    this%data => that%data
    nullify(this%next)

    POP_SUB(INTERNAL(node_copy))
  end subroutine INTERNAL(node_copy)

  ! -----------------------------------------------------
  elemental subroutine INTERNAL(node_end)(this)
    type(INTERNAL(node_t)), intent(inout) :: this

    nullify(this%data, this%next)

  end subroutine INTERNAL(node_end)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(list_len)(this) result(len)
    type(TEMPLATE(list_t)), intent(in) :: this

    integer :: len

    len = this%size

  end function TEMPLATE(list_len)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_init)(this)
    type(TEMPLATE(list_t)), intent(out) :: this

    this%size = 0
    nullify(this%head, this%tail)

  end subroutine INTERNAL(list_init)

  ! -----------------------------------------------------
  function TEMPLATE(list_index)(this, that) result(index)
    type(TEMPLATE(list_t)),       intent(in) :: this
    type(LIST_TYPE_NAME), target, intent(in) :: that

    integer :: index

    type(TEMPLATE(list_iterator_t)) :: iter
    type(LIST_TYPE_NAME),   pointer :: pval
    integer                         :: icnt, ierr

    PUSH_SUB(TEMPLATE(list_index))

    icnt = 0
    index = 0
    call INTERNAL(list_iterator_init_list)(iter, this)
    do
      nullify(pval)
      call INTERNAL(list_iterator_next)(iter, pval, ierr)
      if(ierr/=TEMPLATE(LIST_OK))exit
      icnt = icnt + 1
      if(associated(pval, that))then
        index = icnt
        exit
      end if
    end do
    call INTERNAL(list_iterator_end)(iter)

    POP_SUB(TEMPLATE(list_index))
  end function TEMPLATE(list_index)
  
  ! -----------------------------------------------------
  subroutine INTERNAL(list_walk_node)(this, that, prev, node, next, ierr)
    type(TEMPLATE(list_t)),       intent(in)  :: this
    type(LIST_TYPE_NAME), target, intent(in)  :: that
    type(INTERNAL(node_t)),      pointer      :: prev
    type(INTERNAL(node_t)),      pointer      :: node
    type(INTERNAL(node_t)),      pointer      :: next
    integer,           optional,  intent(out) :: ierr

    type(TEMPLATE(list_iterator_t)) :: iter
    integer                         :: jerr

    PUSH_SUB(INTERNAL(list_walk_node))

    nullify(prev, node, next)
    call INTERNAL(list_iterator_init_list)(iter, this)
    do
      nullify(prev, node, next)
      call INTERNAL(list_iterator_next_node)(iter, prev, node, next, jerr)
      if(jerr/=TEMPLATE(LIST_OK))then
        nullify(prev, node, next)
        exit
      end if
      ASSERT(associated(node))
      if(associated(node%data, that))exit
    end do
    call INTERNAL(list_iterator_end)(iter)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_walk_node))
  end subroutine INTERNAL(list_walk_node)
  
  ! -----------------------------------------------------
  subroutine INTERNAL(list_walk_index)(this, index, prev, node, next, ierr)
    type(TEMPLATE(list_t)),  intent(in)  :: this
    integer,                 intent(in)  :: index
    type(INTERNAL(node_t)), pointer      :: prev
    type(INTERNAL(node_t)), pointer      :: node
    type(INTERNAL(node_t)), pointer      :: next
    integer,      optional,  intent(out) :: ierr

    type(TEMPLATE(list_iterator_t)) :: iter
    integer                         :: jerr, indx

    PUSH_SUB(INTERNAL(list_walk_index))

    nullify(prev, node, next)
    jerr = TEMPLATE(LIST_INDEX_ERROR)
    if((0<index).and.(index<=this%size))then
      indx = 1
      call INTERNAL(list_iterator_init_list)(iter, this)
      do
        nullify(prev, node, next)
        call INTERNAL(list_iterator_next_node)(iter, prev, node, next, jerr)
        if(jerr/=TEMPLATE(LIST_OK))then
          nullify(prev, node, next)
          exit
        end if
        ASSERT(associated(node))
        if(indx==index)exit
        indx = indx + 1
      end do
      call INTERNAL(list_iterator_end)(iter)
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_walk_index))
  end subroutine INTERNAL(list_walk_index)
  
  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_set)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,      optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(TEMPLATE(list_set))

    nullify(prev, node, next)
    call INTERNAL(list_walk_index)(this, index, prev, node, next, ierr)
    if(associated(node))call INTERNAL(node_set)(node, that)
    nullify(prev, node, next)

    POP_SUB(TEMPLATE(list_set))
  end subroutine TEMPLATE(list_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_get)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(in)  :: this
    integer,                intent(in)  :: index
    type(LIST_TYPE_NAME),  pointer      :: that
    integer,     optional,  intent(out) :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(TEMPLATE(list_get))

    nullify(that, prev, node, next)
    call INTERNAL(list_walk_index)(this, index, prev, node, next, ierr)
    if(associated(node))call INTERNAL(node_get)(node, that)
    nullify(prev, node, next)

    POP_SUB(TEMPLATE(list_get))
  end subroutine TEMPLATE(list_get)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_node)(this, prev, node, next)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: prev
    type(INTERNAL(node_t)), pointer        :: node
    type(INTERNAL(node_t)), pointer        :: next

    PUSH_SUB(INTERNAL(list_del_node))

    if(associated(node))then
      if(associated(prev))then
        ASSERT(associated(prev%next,node))
        prev%next => next
      else
        ASSERT(associated(node,this%head))
        this%head => next
      end if
      if(associated(next))then
        ASSERT(associated(node%next,next))
      else
        ASSERT(.not.associated(node%next))
        this%tail => prev
      end if
      nullify(node%next)
    end if

    POP_SUB(INTERNAL(list_del_node))
  end subroutine INTERNAL(list_del_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_node_data)(this, that, node, ierr)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(LIST_TYPE_NAME),    intent(in)    :: that
    type(INTERNAL(node_t)), pointer        :: node
    integer,       optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, next

    PUSH_SUB(INTERNAL(list_del_node_data))

    nullify(prev, node, next)
    call INTERNAL(list_walk_node)(this, that, prev, node, next, ierr)
    if(associated(node))call INTERNAL(list_del_node)(this, prev, node, next)

    POP_SUB(INTERNAL(list_del_node_data))
  end subroutine INTERNAL(list_del_node_data)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_node_index)(this, index, node, ierr)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    integer,                 intent(in)    :: index
    type(INTERNAL(node_t)), pointer        :: node
    integer,       optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, next

    PUSH_SUB(INTERNAL(list_del_node_index))

    nullify(prev, node, next)
    call INTERNAL(list_walk_index)(this, index, prev, node, next, ierr)
    if(associated(node))call INTERNAL(list_del_node)(this, prev, node, next)

    POP_SUB(INTERNAL(list_del_node_index))
  end subroutine INTERNAL(list_del_node_index)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_data)(this, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: that
    integer,     optional,  intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(INTERNAL(list_del_data))

    nullify(node)
    call INTERNAL(list_del_node_data)(this, that, node, ierr)
    if(associated(node)) call INTERNAL(node_del)(node)

    POP_SUB(INTERNAL(list_del_data))
  end subroutine INTERNAL(list_del_data)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_index)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,     optional,  intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(INTERNAL(list_del_index))

    nullify(that, node)
    call INTERNAL(list_del_node_index)(this, index, node, ierr)
    if(associated(node))then
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_del)(node)
    end if

    POP_SUB(INTERNAL(list_del_index))
  end subroutine INTERNAL(list_del_index)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_push_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that

    PUSH_SUB(INTERNAL(list_push_node))

    if(associated(this%head))then
      that%next => this%head
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%tail))
      this%tail => that
    end if
    this%head => that
    this%size = this%size + 1
    ASSERT(this%size>0)

    POP_SUB(INTERNAL(list_push_node))
  end subroutine INTERNAL(list_push_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_push)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: that

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(TEMPLATE(list_push))

    nullify(node)
    call INTERNAL(node_new)(node, that)
    call INTERNAL(list_push_node)(this, node)
    nullify(node)

    POP_SUB(TEMPLATE(list_push))
  end subroutine TEMPLATE(list_push)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_pop_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that

    PUSH_SUB(INTERNAL(list_pop_node))

    nullify(that)
    if(associated(this%head))then
      that => this%head
      this%head => this%head%next
      if(.not.associated(this%head))&
        nullify(this%tail)
      this%size = this%size - 1
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%tail))
    end if
    ASSERT(this%size>=0)

    POP_SUB(INTERNAL(list_pop_node))
  end subroutine INTERNAL(list_pop_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_pop)(this, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,      optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(TEMPLATE(list_pop))

    nullify(that, node)
    if(present(ierr)) ierr = TEMPLATE(LIST_EMPTY_ERROR)
    call INTERNAL(list_pop_node)(this, node)
    if(associated(node))then
      if(present(ierr)) ierr = TEMPLATE(LIST_OK)
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_del)(node)
    end if

    POP_SUB(TEMPLATE(list_pop))
  end subroutine TEMPLATE(list_pop)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_append_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that

    PUSH_SUB(INTERNAL(list_append_node))

    nullify(that%next)
    if(associated(this%tail))then
      this%tail%next => that
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%head))
      this%head => that
    end if
    this%tail => that
    this%size = this%size + 1
    ASSERT(this%size>0)

    POP_SUB(INTERNAL(list_append_node))
  end subroutine INTERNAL(list_append_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_append)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(TEMPLATE(list_append))

    nullify(node)
    call INTERNAL(node_new)(node, that)
    call INTERNAL(list_append_node)(this, node)
    nullify(node)

    POP_SUB(TEMPLATE(list_append))
  end subroutine TEMPLATE(list_append)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_extend)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(TEMPLATE(list_t)), intent(in)    :: that

    type(TEMPLATE(list_iterator_t)) :: iter
    type(LIST_TYPE_NAME),   pointer :: pval
    integer                         :: ierr

    PUSH_SUB(TEMPLATE(list_extend))

    call INTERNAL(list_iterator_init_list)(iter, that)
    do
      nullify(pval)
      call INTERNAL(list_iterator_next)(iter, pval, ierr)
      if(ierr/=TEMPLATE(LIST_OK))exit
      call TEMPLATE(list_append)(this, pval)
    end do
    call INTERNAL(list_iterator_end)(iter)

    POP_SUB(TEMPLATE(list_extend))
  end subroutine TEMPLATE(list_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_copy)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(TEMPLATE(list_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(list_copy))

    call INTERNAL(list_end)(this)
    call TEMPLATE(list_extend)(this, that)

    POP_SUB(INTERNAL(list_copy))
  end subroutine INTERNAL(list_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_end)(this)
    type(TEMPLATE(list_t)), intent(inout) :: this

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(INTERNAL(list_end))

    do
      nullify(node)
      call INTERNAL(list_pop_node)(this, node)
      if(.not.associated(node))exit
      call INTERNAL(node_del)(node)
    end do
    ASSERT(this%size==0)
    ASSERT(.not.associated(this%head))
    ASSERT(.not.associated(this%tail))
    nullify(node)

    POP_SUB(INTERNAL(list_end))
  end subroutine INTERNAL(list_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_list)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_t)),          intent(in)  :: that

    PUSH_SUB(INTERNAL(list_iterator_init_list))

    call INTERNAL(list_iterator_end)(this)
    if(associated(that%head))then
       this%node => that%head
       this%next => this%node%next
    end if

    POP_SUB(INTERNAL(list_iterator_init_list))
  end subroutine INTERNAL(list_iterator_init_list)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_iterator)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)  :: that

    PUSH_SUB(INTERNAL(list_iterator_init_iterator))

    call INTERNAL(list_iterator_copy)(this, that)

    POP_SUB(INTERNAL(list_iterator_init_iterator))
  end subroutine INTERNAL(list_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next_node)(this, prev, node, next, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(INTERNAL(node_t)),         pointer        :: prev
    type(INTERNAL(node_t)),         pointer        :: node
    type(INTERNAL(node_t)),         pointer        :: next
    integer,               optional, intent(out)   :: ierr

    PUSH_SUB(INTERNAL(list_iterator_next_node))

    nullify(prev, node, next)
    if(present(ierr)) ierr = TEMPLATE(LIST_EMPTY_ERROR)
    if(associated(this%node))then
      if(present(ierr)) ierr = TEMPLATE(LIST_OK)
      prev => this%prev
      node => this%node
      next => this%next
      if(associated(this%next))then
        this%prev => this%node
        this%node => this%next
        this%next => this%node%next
      else
        call INTERNAL(list_iterator_end)(this)
      end if
    end if

    POP_SUB(INTERNAL(list_iterator_next_node))
  end subroutine INTERNAL(list_iterator_next_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(INTERNAL(list_iterator_next))

    nullify(that, prev, node, next)
    call INTERNAL(list_iterator_next_node)(this, prev, node, next, ierr)
    if(associated(node)) call INTERNAL(node_get)(node, that)
    nullify(node)

    POP_SUB(INTERNAL(list_iterator_next))
  end subroutine INTERNAL(list_iterator_next)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_copy)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(list_iterator_copy))

    call INTERNAL(list_iterator_end)(this)
    if(associated(that%node))then
      this%prev => that%prev
      this%node => that%node
      this%next => that%next
    end if

    POP_SUB(INTERNAL(list_iterator_copy))
  end subroutine INTERNAL(list_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(list_iterator_end)(this)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this

    nullify(this%prev, this%node, this%next)

  end subroutine INTERNAL(list_iterator_end)

#endif
#if defined(LIST_INCLUDE_MODULE)

end module TEMPLATE(list_m)

#endif

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:

