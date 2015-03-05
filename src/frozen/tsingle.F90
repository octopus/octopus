#include "global.h"
#include "util.h"

!SINGLE: SINGLE_TEMPLATE_NAME
!SINGLE: SINGLE_TYPE_NAME
!SINGLE: SINGLE_TYPE_MODULE_NAME

#if defined(SINGLE_TEMPLATE_NAME)
#if !defined(SINGLE_TYPE_NAME)
#define SINGLE_TYPE_NAME DECORATE(SINGLE_TEMPLATE_NAME,t)
#endif
#if !defined(SINGLE_TYPE_MODULE_NAME)
#define SINGLE_TYPE_MODULE_NAME DECORATE(SINGLE_TEMPLATE_NAME,m)
#endif
#else
#error "'SINGLE_TEMPLATE_NAME' must be defined"
#endif

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX SINGLE_TEMPLATE_NAME
#include "template.h"

#if !defined(SINGLE_INCLUDE_PREFIX)
#if !defined(SINGLE_INCLUDE_HEADER) && !defined(SINGLE_INCLUDE_BODY)

module TEMPLATE(single_m)

  use global_m
  use messages_m
  use profiling_m

  use SINGLE_TYPE_MODULE_NAME, only: &
    SINGLE_TYPE_NAME

  implicit none

  private

  public ::                      &
    TEMPLATE(single_t),          &
    TEMPLATE(single_init),       &
    TEMPLATE(single_associated), &
    TEMPLATE(single_set),        &
    TEMPLATE(single_get),        &
    TEMPLATE(single_copy),       &
    TEMPLATE(single_end)

#endif
#if !defined(SINGLE_INCLUDE_BODY)

  type :: TEMPLATE(single_t)
    private
    type(SINGLE_TYPE_NAME), pointer :: value =>null()
  end type TEMPLATE(single_t)

  interface TEMPLATE(single_init)
    module procedure INTERNAL(single_init_value)
    module procedure INTERNAL(single_init_single)
  end interface TEMPLATE(single_init)

  interface TEMPLATE(single_associated)
    module procedure INTERNAL(single_associated)
    module procedure INTERNAL(single_value_associated)
    module procedure INTERNAL(single_single_associated)
  end interface TEMPLATE(single_associated)

  interface TEMPLATE(single_set)
    module procedure INTERNAL(single_set_null)
    module procedure INTERNAL(single_set_value)
  end interface TEMPLATE(single_set)

#endif
#if !defined(SINGLE_INCLUDE_HEADER) && !defined(SINGLE_INCLUDE_BODY)

contains

#endif
#if !defined(SINGLE_INCLUDE_HEADER)

  ! -----------------------------------------------------	
  subroutine INTERNAL(single_init_value)(this, that)
    type(TEMPLATE(single_t)),         intent(out) :: this
    type(SINGLE_TYPE_NAME), optional, intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(single_init_value))
    call TEMPLATE(single_end)(this)
    if(present(that))&
      call TEMPLATE(single_set)(this, that)
    POP_SUB(INTERNAL(single_init_value))
    return
  end subroutine INTERNAL(single_init_value)

  ! -----------------------------------------------------	
  subroutine INTERNAL(single_init_single)(this, that)
    type(TEMPLATE(single_t)), intent(out) :: this
    type(TEMPLATE(single_t)), intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(single_init_single))
    call TEMPLATE(single_end)(this)
    call TEMPLATE(single_copy)(this, that)
    POP_SUB(INTERNAL(single_init_single))
    return
  end subroutine INTERNAL(single_init_single)

  ! -----------------------------------------------------
  elemental function INTERNAL(single_associated)(this) result(eqv)
    type(TEMPLATE(single_t)), intent(in) :: this
    !
    logical :: eqv
    !
    eqv=(associated(this%value))
    return
  end function INTERNAL(single_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(single_value_associated)(this, that) result(eqv)
    type(TEMPLATE(single_t)),       intent(in) :: this
    type(SINGLE_TYPE_NAME), target, intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(INTERNAL(single_associated)(this))&
      eqv=associated(this%value, that)
    return
  end function INTERNAL(single_value_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(single_single_associated)(this, that) result(eqv)
    type(TEMPLATE(single_t)), intent(in) :: this
    type(TEMPLATE(single_t)), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(INTERNAL(single_associated)(that))&
      eqv=INTERNAL(single_value_associated)(this, that%value)
    return
  end function INTERNAL(single_single_associated)

  ! -----------------------------------------------------
  subroutine INTERNAL(single_set_null)(this)
    type(TEMPLATE(single_t)), intent(inout) :: this
    !
    PUSH_SUB(INTERNAL(single_set_null))
    nullify(this%value)
    POP_SUB(INTERNAL(single_set_null))
    return
  end subroutine INTERNAL(single_set_null)

  ! -----------------------------------------------------
  subroutine INTERNAL(single_set_value)(this, that)
    type(TEMPLATE(single_t)),       intent(inout) :: this
    type(SINGLE_TYPE_NAME), target, intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(single_set_value))
    this%value=>that
    POP_SUB(INTERNAL(single_set_value))
    return
  end subroutine INTERNAL(single_set_value)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_get)(this, that)
    type(TEMPLATE(single_t)), intent(in) :: this
    type(SINGLE_TYPE_NAME),  pointer     :: that
    !
    PUSH_SUB(TEMPLATE(single_get))
    nullify(that)
    if(INTERNAL(single_associated)(this))&
      that=>this%value
    POP_SUB(TEMPLATE(single_get))
    return
  end subroutine TEMPLATE(single_get)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_copy)(this, that)
    type(TEMPLATE(single_t)), intent(inout) :: this
    type(TEMPLATE(single_t)), intent(in)    :: that
    !
    PUSH_SUB(TEMPLATE(single_copy))
    call TEMPLATE(single_set)(this, that%value)
    POP_SUB(TEMPLATE(single_copy))
    return
  end subroutine TEMPLATE(single_copy)

  ! -----------------------------------------------------
  elemental subroutine TEMPLATE(single_end)(this)
    type(TEMPLATE(single_t)), intent(inout) :: this
    !
    nullify(this%value)
    return
  end subroutine TEMPLATE(single_end)

#endif
#if !defined(SINGLE_INCLUDE_HEADER) && !defined(SINGLE_INCLUDE_BODY)

end module TEMPLATE(single_m)

#endif
#endif

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:


