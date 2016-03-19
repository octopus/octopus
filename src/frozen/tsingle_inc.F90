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
#define SINGLE_TYPE_MODULE_NAME DECORATE(SINGLE_TEMPLATE_NAME, oct_m)
#endif
#else
#error "'SINGLE_TEMPLATE_NAME' must be defined"
#endif

#undef SINGLE_INCLUDE_MODULE
#if !defined(SINGLE_INCLUDE_PREFIX) && !defined(SINGLE_INCLUDE_HEADER) && !defined(SINGLE_INCLUDE_BODY)
#define SINGLE_INCLUDE_MODULE
#endif

#if defined(SINGLE_INCLUDE_PREFIX) && defined(SINGLE_INCLUDE_HEADER)
#error "Only one off 'SINGLE_INCLUDE_PREFIX' or 'SINGLE_INCLUDE_HEADER' can be defined."
#endif

#if defined(SINGLE_INCLUDE_PREFIX) && defined(SINGLE_INCLUDE_BODY)
#error "Only one off 'SINGLE_INCLUDE_PREFIX' or 'SINGLE_INCLUDE_BODY' can be defined."
#endif

#if defined(SINGLE_INCLUDE_HEADER) && defined(SINGLE_INCLUDE_BODY)
#error "Only one off 'SINGLE_INCLUDE_HEADER' or 'SINGLE_INCLUDE_BODY' can be defined."
#endif

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX SINGLE_TEMPLATE_NAME
#include "template.h"

#if defined(SINGLE_INCLUDE_MODULE)

module TEMPLATE(single_m)

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  use SINGLE_TYPE_MODULE_NAME

  implicit none

  private

  public ::             &
    TEMPLATE(single_t)

  public ::                      &
    TEMPLATE(single_new),        &
    TEMPLATE(single_del),        &
    TEMPLATE(single_init),       &
    TEMPLATE(single_associated), &
    TEMPLATE(single_set),        &
    TEMPLATE(single_get),        &
    TEMPLATE(single_copy),       &
    TEMPLATE(single_end)

#endif
#if defined(SINGLE_INCLUDE_HEADER) || defined(SINGLE_INCLUDE_MODULE)

  type :: TEMPLATE(single_t)
    private
    type(SINGLE_TYPE_NAME), pointer :: value =>null()
  end type TEMPLATE(single_t)

  interface TEMPLATE(single_new)
    module procedure INTERNAL(single_new_single)
    module procedure INTERNAL(single_new_init)
    module procedure INTERNAL(single_new_copy)
  end interface TEMPLATE(single_new)

  interface TEMPLATE(single_associated)
    module procedure INTERNAL(single_associated)
    module procedure INTERNAL(single_value_associated)
    module procedure INTERNAL(single_single_associated)
  end interface TEMPLATE(single_associated)

#endif
#if defined(SINGLE_INCLUDE_MODULE)

contains

#endif
#if defined(SINGLE_INCLUDE_BODY) || defined(SINGLE_INCLUDE_MODULE)

  ! -----------------------------------------------------
  subroutine INTERNAL(single_new_single)(this)
    type(TEMPLATE(single_t)), pointer :: this

    PUSH_SUB(INTERNAL(single_new_single))

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(INTERNAL(single_new_single))
  end subroutine INTERNAL(single_new_single)

  ! -----------------------------------------------------
  subroutine INTERNAL(single_new_init)(this, that)
    type(TEMPLATE(single_t)), pointer     :: this
    type(SINGLE_TYPE_NAME),    intent(in) :: that

    PUSH_SUB(INTERNAL(single_new_init))

    call INTERNAL(single_new_single)(this)
    call TEMPLATE(single_init)(this, that)

    POP_SUB(INTERNAL(single_new_init))
  end subroutine INTERNAL(single_new_init)

  ! -----------------------------------------------------
  subroutine INTERNAL(single_new_copy)(this, that)
    type(TEMPLATE(single_t)), pointer     :: this
    type(TEMPLATE(single_t)),  intent(in) :: that

    PUSH_SUB(INTERNAL(single_new_copy))

    call INTERNAL(single_new_single)(this)
    call TEMPLATE(single_copy)(this, that)

    POP_SUB(INTERNAL(single_new_copy))
  end subroutine INTERNAL(single_new_copy)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_del)(this)
    type(TEMPLATE(single_t)), pointer :: this

    PUSH_SUB(TEMPLATE(single_del))

    if(associated(this))then
      call TEMPLATE(single_end)(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(TEMPLATE(single_del))
  end subroutine TEMPLATE(single_del)

  ! -----------------------------------------------------	
  subroutine TEMPLATE(single_init)(this, that)
    type(TEMPLATE(single_t)), intent(out) :: this
    type(SINGLE_TYPE_NAME),   intent(in)  :: that

    PUSH_SUB(TEMPLATE(single_init))

    call TEMPLATE(single_set)(this, that)

    POP_SUB(TEMPLATE(single_init))
  end subroutine TEMPLATE(single_init)

  ! -----------------------------------------------------
  elemental function INTERNAL(single_associated)(this) result(eqv)
    type(TEMPLATE(single_t)), intent(in) :: this

    logical :: eqv

    eqv = (associated(this%value))

  end function INTERNAL(single_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(single_value_associated)(this, that) result(eqv)
    type(TEMPLATE(single_t)),       intent(in) :: this
    type(SINGLE_TYPE_NAME), target, intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(associated(this%value))&
      eqv = associated(this%value, that)

  end function INTERNAL(single_value_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(single_single_associated)(this, that) result(eqv)
    type(TEMPLATE(single_t)), intent(in) :: this
    type(TEMPLATE(single_t)), intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(associated(that%value))&
      eqv = INTERNAL(single_value_associated)(this, that%value)

  end function INTERNAL(single_single_associated)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_set)(this, that)
    type(TEMPLATE(single_t)),       intent(inout) :: this
    type(SINGLE_TYPE_NAME), target, intent(in)    :: that

    PUSH_SUB(TEMPLATE(single_set))

    this%value => that

    POP_SUB(TEMPLATE(single_set))
  end subroutine TEMPLATE(single_set)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_get)(this, that)
    type(TEMPLATE(single_t)), intent(in) :: this
    type(SINGLE_TYPE_NAME),  pointer     :: that

    PUSH_SUB(TEMPLATE(single_get))

    nullify(that)
    if(associated(this%value)) that => this%value

    POP_SUB(TEMPLATE(single_get))
  end subroutine TEMPLATE(single_get)

  ! -----------------------------------------------------
  subroutine TEMPLATE(single_copy)(this, that)
    type(TEMPLATE(single_t)), intent(inout) :: this
    type(TEMPLATE(single_t)), intent(in)    :: that

    PUSH_SUB(TEMPLATE(single_copy))

    call TEMPLATE(single_end)(this)
    if(associated(that%value))&
      call TEMPLATE(single_set)(this, that%value)

    POP_SUB(TEMPLATE(single_copy))
  end subroutine TEMPLATE(single_copy)

  ! -----------------------------------------------------
  elemental subroutine TEMPLATE(single_end)(this)
    type(TEMPLATE(single_t)), intent(inout) :: this

    nullify(this%value)

  end subroutine TEMPLATE(single_end)

#endif
#if defined(SINGLE_INCLUDE_MODULE)

end module TEMPLATE(single_m)

#endif

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:


