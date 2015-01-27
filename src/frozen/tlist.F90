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

#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_TYPE_MODULE_NAME
#undef SINGLE_INCLUDE_PREFIX
#undef SINGLE_INCLUDE_HEADER
#undef SINGLE_INCLUDE_BODY

#define SINGLE_TEMPLATE_NAME LIST_TEMPLATE_NAME
#define SINGLE_TYPE_NAME LIST_TYPE_NAME
#define SINGLE_INCLUDE_PREFIX
#include "tsingle.F90"
#undef SINGLE_INCLUDE_PREFIX

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX LIST_TEMPLATE_NAME
#include "template.h"

#if !defined(LIST_INCLUDE_PREFIX)
#if !defined(LIST_INCLUDE_HEADER) && !defined(LIST_INCLUDE_BODY)

module TEMPLATE(list_m)

  use global_m
  use messages_m
  use profiling_m

  use LIST_TYPE_MODULE_NAME, only: &
    LIST_TYPE_NAME

  implicit none

  private

  public ::                     &
    TEMPLATE(LIST_OK),          &
    TEMPLATE(LIST_EMPTY_ERROR), &
    TEMPLATE(LIST_INDEX_ERROR)

  public ::                &
    TEMPLATE(list_t),      &
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

  public ::                    &
    TEMPLATE(list_iterator_t)

#endif
#if !defined(LIST_INCLUDE_BODY)
#define SINGLE_INCLUDE_HEADER
#include "tsingle.F90"
#undef SINGLE_INCLUDE_HEADER
#define TEMPLATE_PREFIX LIST_TEMPLATE_NAME
#include "template.h"

  integer, parameter :: TEMPLATE(LIST_OK)          = 0
  integer, parameter :: TEMPLATE(LIST_EMPTY_ERROR) =-1
  integer, parameter :: TEMPLATE(LIST_INDEX_ERROR) =-2

  type :: INTERNAL(node_t)
    private
    type(EXTERNAL(single_t))        :: data
    type(INTERNAL(node_t)), pointer :: next
  end type INTERNAL(node_t)

  type :: TEMPLATE(list_t)
    private
    integer                         :: size = 0
    type(INTERNAL(node_t)), pointer :: head =>null()
    type(INTERNAL(node_t)), pointer :: tail =>null()
  end type TEMPLATE(list_t)

  type :: TEMPLATE(list_iterator_t)
    private
    type(INTERNAL(node_t)), pointer :: node =>null()
  end type TEMPLATE(list_iterator_t)

  interface INTERNAL(node_associated)
    module procedure INTERNAL(node_associated)
    module procedure INTERNAL(node_node_associated)
  end interface INTERNAL(node_associated)

  interface INTERNAL(node_set)
    module procedure INTERNAL(node_set_data)
    module procedure INTERNAL(node_set_next)
  end interface INTERNAL(node_set)

  interface INTERNAL(node_get)
    module procedure INTERNAL(node_get_data)
    module procedure INTERNAL(node_get_next)
  end interface INTERNAL(node_get)

  interface TEMPLATE(list_init)
    module procedure INTERNAL(list_init)
    module procedure INTERNAL(list_iterator_init_list)
    module procedure INTERNAL(list_iterator_init_iterator)
  end interface TEMPLATE(list_init)

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
#if !defined(LIST_INCLUDE_HEADER) && !defined(LIST_INCLUDE_BODY)

contains
  
#endif
#if !defined(LIST_INCLUDE_HEADER)
#define SINGLE_INCLUDE_BODY
#include "tsingle.F90"
#undef SINGLE_INCLUDE_BODY
#define TEMPLATE_PREFIX LIST_TEMPLATE_NAME
#include "template.h"

  ! -----------------------------------------------------
  subroutine INTERNAL(node_init)(this, that)
    type(INTERNAL(node_t)), intent(out) :: this
    type(LIST_TYPE_NAME),   intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(node_init))
    call EXTERNAL(single_init)(this%data, that)
    nullify(this%next)
    POP_SUB(INTERNAL(node_init))
    return
  end subroutine INTERNAL(node_init)

  ! -----------------------------------------------------
  elemental function INTERNAL(node_associated)(this) result(eqv)
    type(INTERNAL(node_t)), intent(in) :: this
    !
    logical :: eqv
    !
    eqv=associated(this%next)
    return
  end function INTERNAL(node_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(node_node_associated)(this, that) result(eqv)
    type(INTERNAL(node_t)),         intent(in) :: this
    type(INTERNAL(node_t)), target, intent(in) :: that
    !
    logical :: eqv
    !
    eqv=associated(this%next, that)
    return
  end function INTERNAL(node_node_associated)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_set_data)(this, that)
    type(INTERNAL(node_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),   intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(node_set_data))
    call EXTERNAL(single_set)(this%data, that)
    POP_SUB(INTERNAL(node_set_data))
    return
  end subroutine INTERNAL(node_set_data)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_set_next)(this, that)
    type(INTERNAL(node_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that
    !
    PUSH_SUB(INTERNAL(node_set_next))
    this%next=>that
    POP_SUB(INTERNAL(node_set_next))
    return
  end subroutine INTERNAL(node_set_next)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_get_data)(this, that)
    type(INTERNAL(node_t)), intent(in) :: this
    type(LIST_TYPE_NAME),  pointer     :: that
    !
    PUSH_SUB(INTERNAL(node_get_data))
    call EXTERNAL(single_get)(this%data, that)
    POP_SUB(INTERNAL(node_get_data))
    return
  end subroutine INTERNAL(node_get_data)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_get_next)(this, that)
    type(INTERNAL(node_t)),  intent(in) :: this
    type(INTERNAL(node_t)), pointer     :: that
    !
    PUSH_SUB(INTERNAL(node_get_next))
    nullify(that)
    if(associated(this%next))&
      that=>this%next
    POP_SUB(INTERNAL(node_get_next))
    return
  end subroutine INTERNAL(node_get_next)

  ! -----------------------------------------------------
  subroutine INTERNAL(node_copy)(this, that)
    type(INTERNAL(node_t)), intent(inout) :: this
    type(INTERNAL(node_t)), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(node_copy))
    call EXTERNAL(single_copy)(this%data, that%data)
    nullify(this%next)
    POP_SUB(INTERNAL(node_copy))
    return
  end subroutine INTERNAL(node_copy)

  ! -----------------------------------------------------
  elemental subroutine INTERNAL(node_end)(this)
    type(INTERNAL(node_t)), intent(inout) :: this
    !
    call EXTERNAL(single_end)(this%data)
    nullify(this%next)
    return
  end subroutine INTERNAL(node_end)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(list_len)(this) result(len)
    type(TEMPLATE(list_t)), intent(in) :: this
    !
    integer :: len
    !
    len=this%size
    return
  end function TEMPLATE(list_len)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_init)(this)
    type(TEMPLATE(list_t)), intent(out) :: this
    !
    this%size=0
    nullify(this%head, this%tail)
    return
  end subroutine INTERNAL(list_init)

  ! -----------------------------------------------------
  function TEMPLATE(list_index)(this, value) result(index)
    type(TEMPLATE(list_t)),       intent(in) :: this
    type(LIST_TYPE_NAME), target, intent(in) :: value
    !
    integer :: index
    !
    type(TEMPLATE(list_iterator_t)) :: iter
    type(LIST_TYPE_NAME),   pointer :: pval
    integer                         :: icnt, ierr
    !
    PUSH_SUB(TEMPLATE(list_index))
    icnt=0
    index=0
    call INTERNAL(list_iterator_init_list)(iter, this)
    do
      nullify(pval)
      call INTERNAL(list_iterator_next)(iter, pval, ierr)
      if(ierr/=TEMPLATE(LIST_OK))exit
      icnt=icnt+1
      if(associated(pval, value))then
        index=icnt
        exit
      end if
    end do
    call INTERNAL(list_iterator_end)(iter)
    POP_SUB(TEMPLATE(list_index))
    return
  end function TEMPLATE(list_index)
  
  ! -----------------------------------------------------
  subroutine INTERNAL(list_walk)(this, index, node, ierr)
    type(TEMPLATE(list_t)),  intent(in)  :: this
    integer,                 intent(in)  :: index
    type(INTERNAL(node_t)), pointer      :: node
    integer,      optional,  intent(out) :: ierr
    !
    type(TEMPLATE(list_iterator_t)) :: iter
    integer                         :: indx
    !
    PUSH_SUB(INTERNAL(list_walk))
    nullify(node)
    if(present(ierr))ierr=TEMPLATE(LIST_INDEX_ERROR)
    if((0<index).and.(index<=this%size))then
      indx=1
      call INTERNAL(list_iterator_init_list)(iter, this)
      do
        nullify(node)
        call INTERNAL(list_iterator_next_node)(iter, node)
        if((.not.associated(node)).or.(indx==index))exit
        indx=indx+1
      end do
      call INTERNAL(list_iterator_end)(iter)
      ASSERT(associated(node))
      if(present(ierr))ierr=TEMPLATE(LIST_OK)
    end if
    POP_SUB(INTERNAL(list_walk))
    return
  end subroutine INTERNAL(list_walk)
  
  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_set)(this, index, value, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),  pointer        :: value
    integer,      optional, intent(out)   :: ierr
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_set))
    nullify(node)
    call INTERNAL(list_walk)(this, index, node, ierr)
    if(associated(node))call INTERNAL(node_set)(node, value)
    nullify(node)
    POP_SUB(TEMPLATE(list_set))
    return
  end subroutine TEMPLATE(list_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_get)(this, index, value, ierr)
    type(TEMPLATE(list_t)), intent(in)  :: this
    integer,                intent(in)  :: index
    type(LIST_TYPE_NAME),  pointer      :: value
    integer,     optional,  intent(out) :: ierr
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_get))
    nullify(value, node)
    call INTERNAL(list_walk)(this, index, node, ierr)
    if(associated(node))call INTERNAL(node_get)(node, value)
    nullify(node)
    POP_SUB(TEMPLATE(list_get))
    return
  end subroutine TEMPLATE(list_get)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_del_node)(this, index, node)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    integer,                 intent(in)    :: index
    type(INTERNAL(node_t)), pointer        :: node
     !
    type(INTERNAL(node_t)), pointer :: prev, next
    integer                         :: indx
    !
    PUSH_SUB(INTERNAL(list_del_node))
    nullify(prev, node, next)
    node=>this%head
    next=>node%next
    do indx=2, index
      prev=>node
      node=>next
      next=>node%next
    end do
    if(associated(prev))then
      prev%next=>next
    else
      this%head=>next
    end if
    if(.not.associated(next))this%tail=>prev
    nullify(node%next)
    POP_SUB(INTERNAL(list_del_node))
    return
  end subroutine INTERNAL(list_del_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_del)(this, index, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    integer,                intent(in)    :: index
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,     optional,  intent(out)   :: ierr
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_del))
    nullify(that, node)
    if(present(ierr))ierr=TEMPLATE(LIST_INDEX_ERROR)
    if((0<index).and.(index<=this%size))then
      if(present(ierr))ierr=TEMPLATE(LIST_OK)
      call INTERNAL(list_del_node)(this, index, node)
      ASSERT(associated(node))
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_end)(node)
      SAFE_DEALLOCATE_P(node)
      nullify(node)
    end if
    POP_SUB(TEMPLATE(list_del))
    return
  end subroutine TEMPLATE(list_del)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_push_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that
    !
    PUSH_SUB(INTERNAL(list_push_node))
    if(associated(this%head))then
      that%next=>this%head
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%tail))
      this%tail=>that
    end if
    this%head=>that
    this%size=this%size+1
    ASSERT(this%size>0)
    POP_SUB(INTERNAL(list_push_node))
    return
  end subroutine INTERNAL(list_push_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_push)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_push))
    nullify(node)
    SAFE_ALLOCATE(node)
    call INTERNAL(node_init)(node, that)
    call INTERNAL(list_push_node)(this, node)
    nullify(node)
    POP_SUB(TEMPLATE(list_push))
    return
  end subroutine TEMPLATE(list_push)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_pop_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that
    !
    PUSH_SUB(INTERNAL(list_pop_node))
    nullify(that)
    if(associated(this%head))then
      that=>this%head
      this%head=>this%head%next
      if(.not.associated(this%head))&
        nullify(this%tail)
      this%size=this%size-1
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%tail))
    end if
    ASSERT(this%size>=0)
    POP_SUB(INTERNAL(list_pop_node))
    return
  end subroutine INTERNAL(list_pop_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_pop)(this, that, ierr)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_pop))
    nullify(that, node)
    if(present(ierr))ierr=TEMPLATE(LIST_EMPTY_ERROR)
    call INTERNAL(list_pop_node)(this, node)
    if(associated(node))then
      if(present(ierr))ierr=TEMPLATE(LIST_OK)
      call INTERNAL(node_get)(node, that)
      call INTERNAL(node_end)(node)
      SAFE_DEALLOCATE_P(node)
      nullify(node)
    end if
    POP_SUB(TEMPLATE(list_pop))
    return
  end subroutine TEMPLATE(list_pop)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_append_node)(this, that)
    type(TEMPLATE(list_t)),  intent(inout) :: this
    type(INTERNAL(node_t)), pointer        :: that
    !
    PUSH_SUB(INTERNAL(list_append_node))
    nullify(that%next)
    if(associated(this%tail))then
      this%tail%next=>that
    else
      ASSERT(this%size==0)
      ASSERT(.not.associated(this%head))
      this%head=>that
    end if
    this%tail=>that
    this%size=this%size+1
    ASSERT(this%size>0)
    POP_SUB(INTERNAL(list_append_node))
    return
  end subroutine INTERNAL(list_append_node)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_append)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),  pointer        :: that
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_append))
    SAFE_ALLOCATE(node)
    call INTERNAL(node_init)(node, that)
    call INTERNAL(list_append_node)(this, node)
    POP_SUB(TEMPLATE(list_append))
    return
  end subroutine TEMPLATE(list_append)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_extend)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(TEMPLATE(list_t)), intent(in)    :: that
    !
    type(TEMPLATE(list_iterator_t)) :: iter
    type(LIST_TYPE_NAME),   pointer :: pval
    integer                         :: ierr
    !
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
    return
  end subroutine TEMPLATE(list_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_copy)(this, that)
    type(TEMPLATE(list_t)), intent(inout) :: this
    type(TEMPLATE(list_t)), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(list_copy))
    call TEMPLATE(list_end)(this)
    call TEMPLATE(list_extend)(this, that)
    POP_SUB(INTERNAL(list_copy))
    return
  end subroutine INTERNAL(list_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_end)(this)
    type(TEMPLATE(list_t)), intent(inout) :: this
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(INTERNAL(list_end))
    do
      nullify(node)
      call INTERNAL(list_pop_node)(this, node)
      if(.not.associated(node))exit
      call INTERNAL(node_end)(node)
      SAFE_DEALLOCATE_P(node)
    end do
    ASSERT(this%size==0)
    ASSERT(.not.associated(this%head))
    ASSERT(.not.associated(this%tail))
    nullify(node)
    POP_SUB(INTERNAL(list_end))
    return
  end subroutine INTERNAL(list_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_node)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(INTERNAL(node_t)),  target, intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(list_iterator_init_node))
    this%node=>that
    POP_SUB(INTERNAL(list_iterator_init_node))
    return
  end subroutine INTERNAL(list_iterator_init_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_iterator)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(list_iterator_init_iterator))
    call INTERNAL(list_iterator_init_node)(this, that%node)
    POP_SUB(INTERNAL(list_iterator_init_iterator))
    return
  end subroutine INTERNAL(list_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_init_list)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_t)),          intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(list_iterator_init_list))
    call INTERNAL(list_iterator_init_node)(this, that%head)
    POP_SUB(INTERNAL(list_iterator_init_list))
    return
  end subroutine INTERNAL(list_iterator_init_list)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next_node)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(INTERNAL(node_t)),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(INTERNAL(node_t)), pointer :: next
    !
    PUSH_SUB(INTERNAL(list_iterator_next_node))
    nullify(that, next)
    if(present(ierr))ierr=TEMPLATE(LIST_EMPTY_ERROR)
    if(associated(this%node))then
      if(present(ierr))ierr=TEMPLATE(LIST_OK)
      call INTERNAL(node_get)(this%node, next)
      that=>this%node
      this%node=>next
      nullify(next)
    end if
    POP_SUB(INTERNAL(list_iterator_next_node))
    return
  end subroutine INTERNAL(list_iterator_next_node)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_next)(this, that, ierr)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(LIST_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(INTERNAL(list_iterator_next))
    nullify(that, node)
    call INTERNAL(list_iterator_next_node)(this, node, ierr)
    if(associated(node))&
      call INTERNAL(node_get)(node, that)
    nullify(node)
    POP_SUB(INTERNAL(list_iterator_next))
    return
  end subroutine INTERNAL(list_iterator_next)

#if 0
  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_map_1)(this, yield)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    interface
      function yield(this) result(cont)
        use TEMPLATE(iterator_node_m), only: INTERNAL(node_t) => TEMPLATE(INTERNAL(node_t))
        implicit none
        type(INTERNAL(node_t)), intent(in) :: this
        !
        logical :: cont
      end function yield
    end interface
    !
    type(INTERNAL(node_t)), pointer :: node
    !
    PUSH_SUB(INTERNAL(list_iterator_map_1))
    do
      nullify(node)
      call INTERNAL(list_iterator_next_node)(this, node)
      if(.not.associated(node))exit
      if(.not.yield(node))exit
    end do
    POP_SUB(INTERNAL(list_iterator_map_1))
    return
  end subroutine INTERNAL(list_iterator_map_1)

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_map_2)(this, that, yield)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(TEMPLATE(iterator_t)), intent(inout) :: that
    interface
      function yield(this, that) result(cont)
        use INTERNAL(node_m), only: INTERNAL(node_t) => INTERNAL(node_t)
        implicit none
        type(t), intent(in) :: this
        type(INTERNAL(node_t)), intent(in) :: that
        !
        logical :: cont
      end function yield
    end interface
    !
    type(INTERNAL(node_t)), pointer :: nod1, nod2
    !
    PUSH_SUB(INTERNAL(list_iterator_map_2))
    do
      nullify(nod1, nod2)
      call INTERNAL(list_iterator_next_node)(this, nod1)
      call INTERNAL(list_iterator_next_node)(that, nod2)
      if((.not.associated(nod1)).and.(.not.associated(nod2)))exit
      if(.not.yield(nod1, nod2))exit
    end do
    POP_SUB(INTERNAL(list_iterator_map_2))
    return
  end subroutine INTERNAL(list_iterator_map_2)
#endif

  ! ---------------------------------------------------------
  subroutine INTERNAL(list_iterator_copy)(this, that)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(list_iterator_copy))
    nullify(this%node)
    if(associated(that%node))&
      this%node=>that%node
    POP_SUB(INTERNAL(list_iterator_copy))
    return
  end subroutine INTERNAL(list_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(list_iterator_end)(this)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    !
    nullify(this%node)
    return
  end subroutine INTERNAL(list_iterator_end)

#endif
#if !defined(LIST_INCLUDE_HEADER) && !defined(LIST_INCLUDE_BODY)

end module TEMPLATE(list_m)

#endif
#endif

#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:

