#include "global.h"
#include "template.h"

#ifdef MODULE_TYPE_1
#ifdef MODULE_IMPLEMENT_OPS
#define MODULE_INVOCATION_1 use MODULE_TYPE_1, only: operator(==), operator(/=), first_t=>TYPE_1
#else
#define MODULE_INVOCATION_1 use MODULE_TYPE_1, only: first_t=>TYPE_1
#endif
#define ITYPE_1 type(first_t)
#else
#define MODULE_INVOCATION_1
#define ITYPE_1 TYPE_1
#endif

#ifdef MODULE_TYPE_2
#ifdef MODULE_IMPLEMENT_OPS
#define MODULE_INVOCATION_2 use MODULE_TYPE_2, only: operator(==), operator(/=), second_t=>TYPE_2
#else
#define MODULE_INVOCATION_2 use MODULE_TYPE_2, only: second_t=>TYPE_2
#endif
#define ITYPE_2 type(second_t)
#else
#define MODULE_INVOCATION_2
#define ITYPE_2 TYPE_2
#endif

module TEMPLATE(pair_m)

  use global_m
  use messages_m
  use profiling_m

  MODULE_INVOCATION_1
  MODULE_INVOCATION_2

  implicit none

  private

#ifdef MODULE_IMPLEMENT_OPS
  public ::       &
    operator(==), &
    operator(/=)
#endif

  public ::                    &
    TEMPLATE(pair_init),       &
    TEMPLATE(pair_get_first),  &
    TEMPLATE(pair_get_second), &
    TEMPLATE(pair_set_first),  &
    TEMPLATE(pair_set_second), &
    TEMPLATE(pair_copy),       &
    TEMPLATE(pair_end)

#ifdef MODULE_IMPLEMENT_OPS
  interface operator(==)
    module procedure TEMPLATE(pair_equal)
  end interface operator(==)

  interface operator(/=)
    module procedure TEMPLATE(pair_not_equal)
  end interface operator(/=)
#endif

  type, public :: TEMPLATE(pair_t)
    private
    ITYPE_1, pointer :: a =>null()
    ITYPE_2, pointer :: b =>null()
  end type TEMPLATE(pair_t)

contains

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_init)(this, a, b)
    type(TEMPLATE(pair_t)), intent(out) :: this
    ITYPE_1,        target, intent(in)  :: a
    ITYPE_2,        target, intent(in)  :: b
    !
    PUSH_SUB(TEMPLATE(pair_init))
    this%a=>a
    this%b=>b
    POP_SUB(TEMPLATE(pair_init))
    return
  end subroutine TEMPLATE(pair_init)

#ifdef MODULE_IMPLEMENT_OPS
  ! -----------------------------------------------------
  elemental function TEMPLATE(pair_equal)(this, that) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    type(TEMPLATE(pair_t)), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=((this%a==that%a).and.(this%b==that%b))
    return
  end function TEMPLATE(pair_equal)
  
  ! -----------------------------------------------------
  elemental function TEMPLATE(pair_not_equal)(this, that) result(neqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    type(TEMPLATE(pair_t)), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=((this%a/=that%a).or.(this%b/=that%b))
    return
  end function TEMPLATE(pair_not_equal)
#endif

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_get_first)(this, that)
    type(TEMPLATE(pair_t)), intent(in) :: this
    ITYPE_1,               pointer     :: that
    !
    PUSH_SUB(TEMPLATE(pair_get_first))
    that=>null()
    if(associated(this%a))&
      that=>this%a
    POP_SUB(TEMPLATE(pair_get_first))
    return
  end subroutine TEMPLATE(pair_get_first)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_get_second)(this, that)
    type(TEMPLATE(pair_t)), intent(in) :: this
    ITYPE_2,               pointer     :: that
    !
    PUSH_SUB(TEMPLATE(pair_get_second))
    that=>null()
    if(associated(this%b))&
      that=>this%b
    POP_SUB(TEMPLATE(pair_get_second))
    return
  end subroutine TEMPLATE(pair_get_second)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_set_first)(this, that)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    ITYPE_1,        target, intent(in)    :: that
    !
    PUSH_SUB(TEMPLATE(pair_set_first))
    this%a=>that
    POP_SUB(TEMPLATE(pair_set_first))
    return
  end subroutine TEMPLATE(pair_set_first)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_set_second)(this, that)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    ITYPE_2,        target, intent(in)    :: that
    !
    PUSH_SUB(TEMPLATE(pair_set_second))
    this%b=>that
    POP_SUB(TEMPLATE(pair_set_second))
    return
  end subroutine TEMPLATE(pair_set_second)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_copy)(this_out, this_in)
    type(TEMPLATE(pair_t)), intent(out) :: this_out
    type(TEMPLATE(pair_t)), intent(in)  :: this_in
    !
    PUSH_SUB(TEMPLATE(pair_copy))
    this_out%a=>this_in%a
    this_out%b=>this_in%b
    POP_SUB(TEMPLATE(pair_copy))
    return
  end subroutine TEMPLATE(pair_copy)

  ! -----------------------------------------------------
  elemental subroutine TEMPLATE(pair_end)(this)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    !
    nullify(this%a, this%b)
    return
  end subroutine TEMPLATE(pair_end)

end module TEMPLATE(pair_m)

#undef ITYPE_2
#undef MODULE_INVOCATION_2
#undef ITYPE_1
#undef MODULE_INVOCATION_1

!! Local Variables:
!! mode: f90
!! End:
