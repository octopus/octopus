#include "global.h"
#include "util.h"

!PAIR: PAIR_TEMPLATE_NAME

!PAIR: PAIR_FRST_TEMPLATE_NAME
!PAIR: PAIR_FRST_TYPE_NAME
!PAIR: PAIR_FRST_TYPE_MODULE_NAME

!PAIR: PAIR_SCND_TEMPLATE_NAME
!PAIR: PAIR_SCND_TYPE_NAME
!PAIR: PAIR_SCND_TYPE_MODULE_NAME

!PAIR: PAIR_INCLUDE_PREFIX
!PAIR: PAIR_INCLUDE_HEADER
!PAIR: PAIR_INCLUDE_BODY

#if defined(PAIR_FRST_TEMPLATE_NAME)
#if !defined(PAIR_FRST_TYPE_NAME)
#define PAIR_FRST_TYPE_NAME DECORATE(PAIR_FRST_TEMPLATE_NAME,t)
#endif
#if !defined(PAIR_FRST_TYPE_MODULE_NAME)
#define PAIR_FRST_TYPE_MODULE_NAME DECORATE(PAIR_FRST_TEMPLATE_NAME,m)
#endif
#else
#error "'PAIR_FRST_TEMPLATE_NAME' must be defined"
#endif

#if defined(PAIR_SCND_TEMPLATE_NAME)
#if !defined(PAIR_SCND_TYPE_NAME)
#define PAIR_SCND_TYPE_NAME DECORATE(PAIR_SCND_TEMPLATE_NAME,t)
#endif
#if !defined(PAIR_SCND_TYPE_MODULE_NAME)
#define PAIR_SCND_TYPE_MODULE_NAME DECORATE(PAIR_SCND_TEMPLATE_NAME,m)
#endif
#else
#error "'PAIR_SCND_TEMPLATE_NAME' must be defined"
#endif

#undef IPAIR_TMPL_PREFIX
#undef IPAIR_TMPL_NAME
#define IPAIR_TMPL_PREFIX DECORATE(PAIR_FRST_TEMPLATE_NAME,PAIR_SCND_TEMPLATE_NAME)
#if defined(PAIR_TEMPLATE_NAME)
#define IPAIR_TMPL_NAME PAIR_TEMPLATE_NAME
#else
#define IPAIR_TMPL_NAME IPAIR_TMPL_PREFIX
#endif

#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_TYPE_MODULE_NAME
#undef SINGLE_INCLUDE_PREFIX
#undef SINGLE_INCLUDE_HEADER
#undef SINGLE_INCLUDE_BODY

#define SINGLE_TEMPLATE_NAME PAIR_FRST_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_FRST_TYPE_NAME
#define SINGLE_INCLUDE_PREFIX
#include "tsingle.F90"
#undef SINGLE_INCLUDE_PREFIX
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME

#define SINGLE_TEMPLATE_NAME PAIR_SCND_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_SCND_TYPE_NAME
#define SINGLE_INCLUDE_PREFIX
#include "tsingle.F90"
#undef SINGLE_INCLUDE_PREFIX
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME

#undef FRST
#undef SCND
#define FRST PAIR_FRST_TEMPLATE_NAME
#define SCND PAIR_SCND_TEMPLATE_NAME

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX IPAIR_TMPL_NAME
#include "template.h"

#if !defined(PAIR_INCLUDE_PREFIX)
#if !defined(PAIR_INCLUDE_HEADER) && !defined(PAIR_INCLUDE_BODY)

module TEMPLATE(pair_m)

  use global_m
  use messages_m
  use profiling_m

  use PAIR_FRST_TYPE_MODULE_NAME, only: &
    PAIR_FRST_TYPE_NAME

  use PAIR_SCND_TYPE_MODULE_NAME, only: &
    PAIR_SCND_TYPE_NAME

  implicit none

  private

  public ::                    &
    TEMPLATE(pair_t),          &
    TEMPLATE(pair_init),       &
    TEMPLATE(pair_associated), &
    TEMPLATE(pair_get),        &
    TEMPLATE(pair_set),        &
    TEMPLATE(pair_copy),       &
    TEMPLATE(pair_end)

#endif
#if !defined(PAIR_INCLUDE_BODY)
#define SINGLE_INCLUDE_HEADER
#define SINGLE_TEMPLATE_NAME PAIR_FRST_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_FRST_TYPE_NAME
#include "tsingle.F90"
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#define SINGLE_TEMPLATE_NAME PAIR_SCND_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_SCND_TYPE_NAME
#include "tsingle.F90"
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_INCLUDE_HEADER
#define TEMPLATE_PREFIX IPAIR_TMPL_NAME
#include "template.h"

  interface TEMPLATE(pair_get)
    module procedure INTERNAL(pair_get_frst)
    module procedure INTERNAL(pair_get_scnd)
    module procedure INTERNAL(pair_get_both)
  end interface TEMPLATE(pair_get)

  interface TEMPLATE(pair_set)
    module procedure INTERNAL(pair_set_null)
    module procedure INTERNAL(pair_set_frst)
    module procedure INTERNAL(pair_set_scnd)
    module procedure INTERNAL(pair_set_both)
  end interface TEMPLATE(pair_set)

  interface TEMPLATE(pair_associated)
    module procedure INTERNAL(pair_associated)
    module procedure INTERNAL(pair_pair_frst_associated)
    module procedure INTERNAL(pair_pair_scnd_associated)
    module procedure INTERNAL(pair_pair_associated)
  end interface TEMPLATE(pair_associated)

  type :: TEMPLATE(pair_t)
    private
    type(EXTERNAL(FRST,single_t)) :: frst
    type(EXTERNAL(SCND,single_t)) :: scnd
  end type TEMPLATE(pair_t)

#endif
#if !defined(PAIR_INCLUDE_HEADER) && !defined(PAIR_INCLUDE_BODY)

contains

#endif
#if !defined(PAIR_INCLUDE_HEADER)
#define SINGLE_INCLUDE_BODY
#define SINGLE_TEMPLATE_NAME PAIR_FRST_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_FRST_TYPE_NAME
#include "tsingle.F90"
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#define SINGLE_TEMPLATE_NAME PAIR_SCND_TEMPLATE_NAME
#define SINGLE_TYPE_NAME PAIR_SCND_TYPE_NAME
#include "tsingle.F90"
#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_INCLUDE_BODY
#define TEMPLATE_PREFIX IPAIR_TMPL_NAME
#include "template.h"

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_init)(this, first, second)
    type(TEMPLATE(pair_t)),              intent(out) :: this
    type(PAIR_FRST_TYPE_NAME), optional, intent(in)  :: first
    type(PAIR_SCND_TYPE_NAME), optional, intent(in)  :: second
    !
    PUSH_SUB(TEMPLATE(pair_init))
    call EXTERNAL(FRST,single_init)(this%frst, first)
    call EXTERNAL(SCND,single_init)(this%scnd, second)
    POP_SUB(TEMPLATE(pair_init))
    return
  end subroutine TEMPLATE(pair_init)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_frst_associated)(this) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    !
    logical :: eqv
    !
    eqv=EXTERNAL(FRST,single_associated)(this%frst)
    return
  end function INTERNAL(pair_frst_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(pair_scnd_associated)(this) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    !
    logical :: eqv
    !
    eqv=EXTERNAL(SCND,single_associated)(this%scnd)
    return
  end function INTERNAL(pair_scnd_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(pair_associated)(this) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    !
    logical :: eqv
    !
    eqv=(INTERNAL(pair_frst_associated)(this).and.&
      INTERNAL(pair_scnd_associated)(this))
    return
  end function INTERNAL(pair_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_frst_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)),        intent(in) :: this
    type(EXTERNAL(FRST,single_t)), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(INTERNAL(pair_frst_associated)(this))&
      eqv=EXTERNAL(FRST,single_associated)(this%frst, that)
    return
  end function INTERNAL(pair_pair_frst_associated)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_scnd_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)),        intent(in) :: this
    type(EXTERNAL(SCND,single_t)), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(INTERNAL(pair_scnd_associated)(this))&
      eqv=EXTERNAL(SCND,single_associated)(this%scnd, that)
    return
  end function INTERNAL(pair_pair_scnd_associated)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    type(TEMPLATE(pair_t)), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(INTERNAL(pair_associated)(this).and.INTERNAL(pair_associated)(that))&
      eqv=(EXTERNAL(FRST,single_associated)(this%frst, that%frst).and.&
      EXTERNAL(SCND,single_associated)(this%scnd, that%scnd))
    return
  end function INTERNAL(pair_pair_associated)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_frst)(this, that)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_FRST_TYPE_NAME), pointer     :: that
    !
    PUSH_SUB(INTERNAL(pair_get_frst))
    call EXTERNAL(FRST,single_get)(this%frst, that)
    POP_SUB(INTERNAL(pair_get_frst))
    return
  end subroutine INTERNAL(pair_get_frst)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_scnd)(this, that)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_SCND_TYPE_NAME), pointer     :: that
    !
    PUSH_SUB(INTERNAL(pair_get_scnd))
    ASSERT(INTERNAL(pair_frst_associated)(this))
    call EXTERNAL(SCND,single_get)(this%scnd, that)
    POP_SUB(INTERNAL(pair_get_scnd))
    return
  end subroutine INTERNAL(pair_get_scnd)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_both)(this, frst, scnd)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_FRST_TYPE_NAME), pointer     :: frst
    type(PAIR_SCND_TYPE_NAME), pointer     :: scnd
    !
    PUSH_SUB(INTERNAL(pair_get_both))
    call EXTERNAL(FRST,single_get)(this%frst, frst)
    call EXTERNAL(SCND,single_get)(this%scnd, scnd)
    POP_SUB(INTERNAL(pair_get_both))
    return
  end subroutine INTERNAL(pair_get_both)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_null)(this)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    !
    PUSH_SUB(INTERNAL(pair_set_null))
    call EXTERNAL(FRST,single_set)(this%frst)
    call EXTERNAL(SCND,single_set)(this%scnd)
    POP_SUB(INTERNAL(pair_set_null))
    return
  end subroutine INTERNAL(pair_set_null)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_frst)(this, that)
    type(TEMPLATE(pair_t)),    intent(inout) :: this
    type(PAIR_FRST_TYPE_NAME), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(pair_set_frst))
    call EXTERNAL(FRST,single_set)(this%frst, that)
    POP_SUB(INTERNAL(pair_set_frst))
    return
  end subroutine INTERNAL(pair_set_frst)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_scnd)(this, that)
    type(TEMPLATE(pair_t)),    intent(inout) :: this
    type(PAIR_SCND_TYPE_NAME), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(pair_set_scnd))
    ASSERT(INTERNAL(pair_frst_associated)(this))
    call EXTERNAL(SCND,single_set)(this%scnd, that)
    POP_SUB(INTERNAL(pair_set_scnd))
    return
  end subroutine INTERNAL(pair_set_scnd)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_both)(this, frst, scnd)
    type(TEMPLATE(pair_t)),    intent(inout) :: this
    type(PAIR_FRST_TYPE_NAME), intent(in)    :: frst
    type(PAIR_SCND_TYPE_NAME), intent(in)    :: scnd
    !
    PUSH_SUB(INTERNAL(pair_set_both))
    call EXTERNAL(FRST,single_set)(this%frst, frst)
    call EXTERNAL(SCND,single_set)(this%scnd, scnd)
    POP_SUB(INTERNAL(pair_set_both))
    return
  end subroutine INTERNAL(pair_set_both)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_copy)(this, that)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    type(TEMPLATE(pair_t)), intent(in)    :: that
    !
    PUSH_SUB(TEMPLATE(pair_copy))
    call EXTERNAL(FRST,single_copy)(this%frst, that%frst)
    call EXTERNAL(SCND,single_copy)(this%scnd, that%scnd)
    POP_SUB(TEMPLATE(pair_copy))
    return
  end subroutine TEMPLATE(pair_copy)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_end)(this)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    !
    PUSH_SUB(TEMPLATE(pair_end))
    call EXTERNAL(FRST,single_end)(this%frst)
    call EXTERNAL(SCND,single_end)(this%scnd)
    POP_SUB(TEMPLATE(pair_end))
    return
  end subroutine TEMPLATE(pair_end)

#endif

#if !defined(PAIR_INCLUDE_HEADER) && !defined(PAIR_INCLUDE_BODY)

end module TEMPLATE(pair_m)

#endif
#endif

#undef FRST
#undef SCND

#undef IPAIR_TMPL_PREFIX
#undef IPAIR_TMPL_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:



