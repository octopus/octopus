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
#define LIST_TYPE_MODULE_NAME DECORATE(LIST_TEMPLATE_NAME,oct_m)
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

#undef TEMPLATE_NAME
#define TEMPLATE_NAME LIST_TEMPLATE_NAME
#include "template.h"

#if defined(LIST_INCLUDE_MODULE)

module TEMPLATE(list_oct_m)

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

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
    TEMPLATE(list_insert), &
    TEMPLATE(list_remove), &
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
    integer                         :: size = -1
    type(INTERNAL(node_t)), pointer :: head =>null()
    type(INTERNAL(node_t)), pointer :: tail =>null()
  end type TEMPLATE(list_t)

  type :: TEMPLATE(list_iterator_t)
    private
    type(TEMPLATE(list_t)), pointer :: self =>null()
    type(INTERNAL(node_t)), pointer :: prev =>null()
    type(INTERNAL(node_t)), pointer :: curr =>null()
    type(INTERNAL(node_t)), pointer :: next =>null()
  end type TEMPLATE(list_iterator_t)

  interface TEMPLATE(list_init)
    module procedure INTERNAL(list_init)
    module procedure INTERNAL(list_iterator_init_list)
    module procedure INTERNAL(list_iterator_init_iterator)
  end interface TEMPLATE(list_init)

  interface TEMPLATE(list_insert)
    module procedure INTERNAL(list_ins_data)
    module procedure INTERNAL(list_ins_index)
    module procedure INTERNAL(list_iterator_insert)
  end interface TEMPLATE(list_insert)

  interface TEMPLATE(list_remove)
    module procedure INTERNAL(list_del_data)
    module procedure INTERNAL(list_del_index)
    module procedure INTERNAL(list_iterator_remove)
  end interface TEMPLATE(list_remove)

  interface TEMPLATE(list_del)
    module procedure INTERNAL(list_del_data)
    module procedure INTERNAL(list_del_index)
  end interface TEMPLATE(list_del)

  interface TEMPLATE(list_push)
    module procedure INTERNAL(list_push_type)
  end interface TEMPLATE(list_push)

  interface TEMPLATE(list_pop)
    module procedure INTERNAL(list_pop_type)
  end interface TEMPLATE(list_pop)

  interface TEMPLATE(list_next)
    module procedure INTERNAL(list_iterator_next)
  end interface TEMPLATE(list_next)

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

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
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
  subroutine INTERNAL(list_ins_node)(this, that, prev, node)
    type(TEMPLATE(list_t)),          intent(inout) :: this
    type(INTERNAL(node_t)),  target, intent(inout) :: that
    type(INTERNAL(node_t)), pointer                :: prev
    type(INTERNAL(node_t)), pointer                :: node

    PUSH_SUB(INTERNAL(list_ins_node))

    ASSERT(associated(node))
    if(associated(prev))then
      ASSERT(associated(prev%next,node))
      prev%next => that
    else
      ASSERT(associated(this%head,node))
      this%head => that
    end if
    this%size = this%size + 1
    that%next => node

    POP_SUB(INTERNAL(list_ins_node))
  end subroutine INTERNAL(list_ins_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_ins_data)(this, match, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: match
    type(LIST_TYPE_NAME),   intent(in)    :: that
    integer,     optional,  intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node, prev, curr, next
    integer                         :: jerr

    PUSH_SUB(INTERNAL(list_ins_data))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(node, prev, curr, next)
    call INTERNAL(list_walk_node)(this, match, prev, curr, next, jerr)
    if(jerr==TEMPLATE(LIST_OK))then
      ASSERT(associated(curr))
      call INTERNAL(node_new)(node, that)
      call INTERNAL(list_ins_node)(this, node, prev, curr)
    end if
    nullify(node, prev, curr, next)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_ins_data))
  end subroutine INTERNAL(list_ins_data)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_ins_index)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),   intent(in)    :: that
    integer,      optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node, prev, curr, next
    integer                         :: jerr

    PUSH_SUB(INTERNAL(list_ins_index))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(node, prev, curr, next)
    call INTERNAL(list_walk_index)(this, index, prev, curr, next, jerr)
    if(jerr==TEMPLATE(LIST_OK))then
      ASSERT(associated(curr))
      call INTERNAL(node_new)(node, that)
      call INTERNAL(list_ins_node)(this, node, prev, curr)
    end if
    nullify(node, prev, curr, next)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_ins_index))
  end subroutine INTERNAL(list_ins_index)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_set)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),   intent(in)    :: that
    integer,      optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(TEMPLATE(list_set))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
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

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
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

    ASSERT(associated(node))
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
    this%size = this%size - 1
    nullify(node%next)

    POP_SUB(INTERNAL(list_del_node))
  end subroutine INTERNAL(list_del_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_data)(this, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: that
    integer,     optional,  intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(INTERNAL(list_del_data))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(prev, node, next)
    call INTERNAL(list_walk_node)(this, that, prev, node, next, ierr)
    if(associated(node))then
      call INTERNAL(list_del_node)(this, prev, node, next)
      call INTERNAL(node_del)(node)
    end if
    nullify(prev, node, next)

    POP_SUB(INTERNAL(list_del_data))
  end subroutine INTERNAL(list_del_data)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_index)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,     optional,  intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, node, next

    PUSH_SUB(INTERNAL(list_del_index))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(that, prev, node, next)
    call INTERNAL(list_walk_index)(this, index, prev, node, next, ierr)
    if(associated(node))then
      call INTERNAL(list_del_node)(this, prev, node, next)
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_del)(node)
    end if
    nullify(prev, node, next)

    POP_SUB(INTERNAL(list_del_index))
  end subroutine INTERNAL(list_del_index)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_push_node)(this, that)
    type(TEMPLATE(list_t)),         intent(inout) :: this
    type(INTERNAL(node_t)), target, intent(inout) :: that

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
  subroutine INTERNAL(list_push_type)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: that

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(INTERNAL(list_push_type))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(node)
    call INTERNAL(node_new)(node, that)
    call INTERNAL(list_push_node)(this, node)
    nullify(node)

    POP_SUB(INTERNAL(list_push_type))
  end subroutine INTERNAL(list_push_type)

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
  subroutine INTERNAL(list_pop_type)(this, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,      optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(INTERNAL(list_pop_type))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    nullify(that, node)
    if(present(ierr)) ierr = TEMPLATE(LIST_EMPTY_ERROR)
    call INTERNAL(list_pop_node)(this, node)
    if(associated(node))then
      if(present(ierr)) ierr = TEMPLATE(LIST_OK)
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_del)(node)
    end if

    POP_SUB(INTERNAL(list_pop_type))
  end subroutine INTERNAL(list_pop_type)

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
    type(LIST_TYPE_NAME),   intent(in)    :: that

    type(INTERNAL(node_t)), pointer :: node

    PUSH_SUB(TEMPLATE(list_append))

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
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

    ASSERT(.not.TEMPLATE(list_len)(this)<0)
    ASSERT(.not.TEMPLATE(list_len)(that)<0)
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
    this%size = 0
    if(TEMPLATE(list_len)(that)>0)&
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
    this%size = -1
    nullify(node)

    POP_SUB(INTERNAL(list_end))
  end subroutine INTERNAL(list_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_list)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_t)),  target, intent(in)  :: that

    PUSH_SUB(INTERNAL(list_iterator_init_list))

    ASSERT(.not.TEMPLATE(list_len)(that)<0)
    this%self => that
    nullify(this%prev, this%curr, this%next)
    if(associated(that%head)) this%next => that%head

    POP_SUB(INTERNAL(list_iterator_init_list))
  end subroutine INTERNAL(list_iterator_init_list)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_iterator)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)  :: that

    PUSH_SUB(INTERNAL(list_iterator_init_iterator))

    ASSERT(associated(that%self))
    call INTERNAL(list_iterator_copy)(this, that)

    POP_SUB(INTERNAL(list_iterator_init_iterator))
  end subroutine INTERNAL(list_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_insert)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),            intent(in)    :: that
    integer,               optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: node
    integer                         :: jerr

    PUSH_SUB(INTERNAL(list_iterator_insert))

    ASSERT(associated(this%self))
    nullify(node)
    jerr = TEMPLATE(LIST_EMPTY_ERROR)
    if(associated(this%curr))then
      jerr = TEMPLATE(LIST_OK)
      call INTERNAL(node_new)(node, that)
      call INTERNAL(list_ins_node)(this%self, node, this%prev, this%curr)
      this%prev => node
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_iterator_insert))
  end subroutine INTERNAL(list_iterator_insert)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_remove)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(INTERNAL(list_iterator_remove))

    ASSERT(associated(this%self))
    nullify(that)
    jerr = TEMPLATE(LIST_EMPTY_ERROR)
    if(associated(this%curr))then
      jerr = TEMPLATE(LIST_OK)
      call INTERNAL(node_get)(this%curr, that)
      call INTERNAL(list_del_node)(this%self, this%prev, this%curr, this%next)
      call INTERNAL(node_del)(this%curr)
      nullify(this%curr)
      if(associated(this%next))then
        this%curr => this%next
        this%next => this%curr%next
      end if
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_iterator_remove))
  end subroutine INTERNAL(list_iterator_remove)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next_node)(this, prev, curr, next, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(INTERNAL(node_t)),         pointer        :: prev
    type(INTERNAL(node_t)),         pointer        :: curr
    type(INTERNAL(node_t)),         pointer        :: next
    integer,               optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(INTERNAL(list_iterator_next_node))

    ASSERT(associated(this%self))
    nullify(prev, curr, next)
    jerr = TEMPLATE(LIST_EMPTY_ERROR)
    if(associated(this%next))then
      jerr = TEMPLATE(LIST_OK)
      this%prev => this%curr
      this%curr => this%next
      this%next => this%curr%next
      prev => this%prev
      curr => this%curr
      next => this%next
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(list_iterator_next_node))
  end subroutine INTERNAL(list_iterator_next_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr

    type(INTERNAL(node_t)), pointer :: prev, curr, next

    PUSH_SUB(INTERNAL(list_iterator_next))

    ASSERT(associated(this%self))
    nullify(that, prev, curr, next)
    call INTERNAL(list_iterator_next_node)(this, prev, curr, next, ierr)
    if(associated(curr)) call INTERNAL(node_get)(curr, that)
    nullify(prev, curr, next)

    POP_SUB(INTERNAL(list_iterator_next))
  end subroutine INTERNAL(list_iterator_next)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_copy)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(list_iterator_copy))

    call INTERNAL(list_iterator_end)(this)
    this%self => that%self
    if(associated(that%self))then
      this%prev => that%prev
      this%curr => that%curr
      this%next => that%next
    end if

    POP_SUB(INTERNAL(list_iterator_copy))
  end subroutine INTERNAL(list_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(list_iterator_end)(this)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this

    nullify(this%self, this%prev, this%curr, this%next)

  end subroutine INTERNAL(list_iterator_end)

#endif
#if defined(LIST_INCLUDE_MODULE)

end module TEMPLATE(list_oct_m)

#endif

#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

