#include "global.h"
#include "template.h"

#ifdef MODULE_TYPE
#ifdef TYPE
#define MODULE_INVOCATION use MODULE_TYPE, only: value_t=>TYPE
#define ITYPE type(value_t)
#else
#define MODULE_INVOCATION use MODULE_TYPE, only: value_t=>TEMPLATE(t)
#define ITYPE type(value_t)
#endif
#else
#ifdef TYPE
#define MODULE_INVOCATION
#define ITYPE TYPE
#else
#define MODULE_INVOCATION use TEMPLATE(m), only: value_t=>TEMPLATE(t)
#define ITYPE type(value_t)
#endif
#endif

module TEMPLATE(list_m)

  use global_m
  use messages_m
  use profiling_m

  MODULE_INVOCATION

  implicit none

  private

  public ::       &
    operator(==), &
    operator(/=)

  public ::                &
    TEMPLATE(list_init),   &
    TEMPLATE(list_len),    &
    TEMPLATE(list_index),  &
    TEMPLATE(list_push),   &
    TEMPLATE(list_pop),    &
    TEMPLATE(list_append), &
    TEMPLATE(list_copy),   &
    TEMPLATE(list_end)

  public ::                       &
    TEMPLATE(list_iterator_init), &
    TEMPLATE(list_iterator_next), &
    TEMPLATE(list_iterator_copy), &
    TEMPLATE(list_iterator_end)

  interface operator(==)
    module procedure TEMPLATE(list_equal)
  end interface operator(==)

  interface operator(/=)
    module procedure TEMPLATE(list_not_equal)
  end interface operator(/=)

  type, private :: node_t
    private
    ITYPE,        pointer :: value =>null()
    type(node_t), pointer :: next  =>null()
  end type node_t

  type, public :: TEMPLATE(list_t)
    private
    integer               :: size = 0
    type(node_t), pointer :: head =>null()
    type(node_t), pointer :: tail =>null()
  end type TEMPLATE(list_t)

  type, public :: TEMPLATE(list_iterator_t)
    private
    type(node_t), pointer :: node =>null()
  end type TEMPLATE(list_iterator_t)

  interface TEMPLATE(list_iterator_init)
    module procedure TEMPLATE(list_iterator_init_list)
    module procedure TEMPLATE(list_iterator_init_iterator)
  end interface TEMPLATE(list_iterator_init)

contains

  ! ---------------------------------------------------------
  elemental subroutine TEMPLATE(list_init)(this)
    type(TEMPLATE(list_t)), intent(out) :: this
    !
    this%size=0
    this%head=>null()
    this%tail=>null()
    return
  end subroutine TEMPLATE(list_init)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(list_len)(this) result(len)
    type(TEMPLATE(list_t)), intent(in) :: this
    !
    integer :: len
    !
    len=this%size
    return
  end function TEMPLATE(list_len)

  ! -----------------------------------------------------
  function TEMPLATE(list_equal)(this, that) result(eqv)
    type(TEMPLATE(list_t)), intent(in) :: this
    type(TEMPLATE(list_t)), intent(in) :: that
    !
    logical :: eqv
    !
    type(node_t), pointer :: n1, n2
    !
    PUSH_SUB(TEMPLATE(list_equal))
    eqv=.false.
    nullify(n1, n2)
    if(this%size==that%size)then
      eqv=.true.
      n1=>this%head
      n2=>that%head
      do
        if((.not.associated(n1)).and.(.not.associated(n2)))exit
        if(.not.associated(n1%value,n2%value))then
          eqv=.false.
          exit
        end if
        n1=>n1%next
        n2=>n2%next
      end do
    end if
    POP_SUB(TEMPLATE(list_equal))
    return
  end function TEMPLATE(list_equal)
  
  ! -----------------------------------------------------
  function TEMPLATE(list_not_equal)(this, that) result(neqv)
    type(TEMPLATE(list_t)), intent(in) :: this
    type(TEMPLATE(list_t)), intent(in) :: that
    !
    logical :: neqv
    !
    PUSH_SUB(TEMPLATE(list_not_equal))
    neqv=.not.TEMPLATE(list_equal)(this, that)
    POP_SUB(TEMPLATE(list_not_equal))
    return
  end function TEMPLATE(list_not_equal)

  ! -----------------------------------------------------
  function TEMPLATE(list_index)(this, value) result(index)
    type(TEMPLATE(list_t)), intent(in) :: this
    ITYPE,          target, intent(in) :: value
    !
    integer :: index
    !
    type(node_t), pointer :: node
    integer               :: cnt
    !
    PUSH_SUB(TEMPLATE(list_index))
    cnt=0
    index=0
    node=>this%head
    do
      if(.not.associated(node))exit
      cnt=cnt+1
      if(associated(node%value,value))then
        index=cnt
        exit
      end if
      node=>node%next
    end do
    POP_SUB(TEMPLATE(list_index))
    return
  end function TEMPLATE(list_index)
  
  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_push)(this, value)
    type(TEMPLATE(list_t)), intent(inout) :: this
    ITYPE,          target, intent(in)    :: value
    !
    type(node_t), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_push))
    SAFE_ALLOCATE(node)
    node%value=>value
    node%next=>this%head
    if(.not.associated(this%head))this%tail=>node
    this%head=>node
    this%size=this%size+1
    ASSERT(this%size>0)
    POP_SUB(TEMPLATE(list_push))
    return
  end subroutine TEMPLATE(list_push)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_pop)(this, value)
    type(TEMPLATE(list_t)), intent(inout) :: this
    ITYPE,                 pointer        :: value
    !
    type(node_t), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_pop))
    value=>null()
    node=>this%head
    if(associated(node))then
       value=>node%value
       this%head=>node%next
       this%size=this%size-1
       if(.not.associated(this%head))this%tail=>null()
       SAFE_DEALLOCATE_P(node)
       node=>null()
    end if
    ASSERT(this%size>=0)
    POP_SUB(TEMPLATE(list_pop))
    return
  end subroutine TEMPLATE(list_pop)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_append)(this, value)
    type(TEMPLATE(list_t)), intent(inout) :: this
    ITYPE,          target, intent(in)    :: value
    !
    type(node_t), pointer :: node
    !
    PUSH_SUB(TEMPLATE(list_append))
    SAFE_ALLOCATE(node)
    node%value=>value
    node%next=>null()
    if(associated(this%head))then
      this%tail%next=>node
    else
      this%head=>node
    end if
    this%tail=>node
    this%size=this%size+1
    POP_SUB(TEMPLATE(list_append))
    return
  end subroutine TEMPLATE(list_append)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_extend)(this_out, this_in)
    type(TEMPLATE(list_t)), intent(inout) :: this_out
    type(TEMPLATE(list_t)), intent(in)    :: this_in
    !
    type(TEMPLATE(list_iterator_t)) :: iter
    ITYPE,                  pointer :: value
    !
    PUSH_SUB(TEMPLATE(list_extend))
    call TEMPLATE(list_iterator_init_list)(iter, this_in)
    do
      value=>null()
      call TEMPLATE(list_iterator_next)(iter, value)
      if(.not.associated(value))exit
      call TEMPLATE(list_append)(this_out, value)
    end do
    POP_SUB(TEMPLATE(list_extend))
    return
  end subroutine TEMPLATE(list_extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_copy)(this_out, this_in)
    type(TEMPLATE(list_t)), intent(out) :: this_out
    type(TEMPLATE(list_t)), intent(in)  :: this_in
    !
    type(TEMPLATE(list_iterator_t)) :: iter
    ITYPE,                  pointer :: value
    !
    PUSH_SUB(TEMPLATE(list_copy))
    call TEMPLATE(list_init)(this_out)
    call TEMPLATE(list_iterator_init_list)(iter, this_in)
    do
      value=>null()
      call TEMPLATE(list_iterator_next)(iter, value)
      if(.not.associated(value))exit
      call TEMPLATE(list_append)(this_out, value)
    end do
    POP_SUB(TEMPLATE(list_copy))
    return
  end subroutine TEMPLATE(list_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_end)(this)
    type(TEMPLATE(list_t)), intent(inout) :: this
    !
    ITYPE, pointer :: value
    !
    PUSH_SUB(TEMPLATE(list_end))
    value=>null()
    do
      call TEMPLATE(list_pop)(this, value)
      if(.not.associated(value))exit
    end do
    ASSERT(this%size==0)
    ASSERT(.not.associated(this%head))
    ASSERT(.not.associated(this%tail))
    POP_SUB(TEMPLATE(list_end))
    return
  end subroutine TEMPLATE(list_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_iterator_init_list)(this, list)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_t)),          intent(in)  :: list
    !
    PUSH_SUB(TEMPLATE(list_iterator_init_list))
    this%node=>list%head
    POP_SUB(TEMPLATE(list_iterator_init_list))
    return
  end subroutine TEMPLATE(list_iterator_init_list)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_iterator_init_iterator)(this, iter)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this
    type(TEMPLATE(list_iterator_t)), intent(in)  :: iter
    !
    PUSH_SUB(TEMPLATE(list_iterator_init_iterator))
    this%node=>iter%node
    POP_SUB(TEMPLATE(list_iterator_init_iterator))
    return
  end subroutine TEMPLATE(list_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_iterator_next)(this, value)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    ITYPE,                  pointer, intent(out)   :: value
    !
    PUSH_SUB(TEMPLATE(list_iterator_next))
    value=>null()
    if(associated(this%node))then
      value=>this%node%value
      this%node=>this%node%next
    end if
    POP_SUB(TEMPLATE(list_iterator_next))
    return
  end subroutine TEMPLATE(list_iterator_next)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(list_iterator_copy)(this_out, this_in)
    type(TEMPLATE(list_iterator_t)), intent(out) :: this_out
    type(TEMPLATE(list_iterator_t)), intent(in)  :: this_in
    !
    PUSH_SUB(TEMPLATE(list_iterator_copy))
    this_out%node=>this_in%node
    POP_SUB(TEMPLATE(list_iterator_copy))
    return
  end subroutine TEMPLATE(list_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine TEMPLATE(list_iterator_end)(this)
    type(TEMPLATE(list_iterator_t)), intent(inout) :: this
    !
    this%node=>null()
    return
  end subroutine TEMPLATE(list_iterator_end)

end module TEMPLATE(list_m)

#undef ITYPE
#undef MODULE_INVOCATION

!! Local Variables:
!! mode: f90
!! End:
