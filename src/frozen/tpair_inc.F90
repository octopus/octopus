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
#define PAIR_FRST_TYPE_MODULE_NAME DECORATE(PAIR_FRST_TEMPLATE_NAME, oct_m)
#endif
#else
#if !defined(PAIR_FRST_TYPE_NAME)
#error "'PAIR_FRST_TEMPLATE_NAME' or 'PAIR_FRST_TYPE_NAME' must be defined"
#endif
#endif

#if defined(PAIR_SCND_TEMPLATE_NAME)
#if !defined(PAIR_SCND_TYPE_NAME)
#define PAIR_SCND_TYPE_NAME DECORATE(PAIR_SCND_TEMPLATE_NAME,t)
#endif
#if !defined(PAIR_SCND_TYPE_MODULE_NAME)
#define PAIR_SCND_TYPE_MODULE_NAME DECORATE(PAIR_SCND_TEMPLATE_NAME, oct_m)
#endif
#else
#if !defined(PAIR_SCND_TYPE_NAME)
#error "'PAIR_SCND_TEMPLATE_NAME' or 'PAIR_SCND_TYPE_NAME' must be defined"
#endif
#endif

#undef IPAIR_TMPL_PREFIX
#undef IPAIR_TMPL_NAME
#define IPAIR_TMPL_PREFIX DECORATE(PAIR_FRST_TEMPLATE_NAME,PAIR_SCND_TEMPLATE_NAME)
#if defined(PAIR_TEMPLATE_NAME)
#define IPAIR_TMPL_NAME PAIR_TEMPLATE_NAME
#else
#define IPAIR_TMPL_NAME IPAIR_TMPL_PREFIX
#endif

#undef PAIR_INCLUDE_MODULE
#if !defined(PAIR_INCLUDE_PREFIX) && !defined(PAIR_INCLUDE_HEADER) && !defined(PAIR_INCLUDE_BODY)
#define PAIR_INCLUDE_MODULE
#endif

#if defined(PAIR_INCLUDE_PREFIX) && defined(PAIR_INCLUDE_HEADER)
#error "Only one off 'PAIR_INCLUDE_PREFIX' or 'PAIR_INCLUDE_HEADER' can be defined."
#endif

#if defined(PAIR_INCLUDE_PREFIX) && defined(PAIR_INCLUDE_BODY)
#error "Only one off 'PAIR_INCLUDE_PREFIX' or 'PAIR_INCLUDE_BODY' can be defined."
#endif

#if defined(PAIR_INCLUDE_HEADER) && defined(PAIR_INCLUDE_BODY)
#error "Only one off 'PAIR_INCLUDE_HEADER' or 'PAIR_INCLUDE_BODY' can be defined."
#endif

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX IPAIR_TMPL_NAME
#include "template.h"

#if defined(PAIR_INCLUDE_MODULE)

module TEMPLATE(pair_m)

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  use PAIR_FRST_TYPE_MODULE_NAME

  use PAIR_SCND_TYPE_MODULE_NAME

  implicit none

  private

  public ::           &
    TEMPLATE(pair_t)

  public ::                    &
    TEMPLATE(pair_new),        &
    TEMPLATE(pair_del),        &
    TEMPLATE(pair_init),       &
    TEMPLATE(pair_associated), &
    TEMPLATE(pair_get),        &
    TEMPLATE(pair_set),        &
    TEMPLATE(pair_copy),       &
    TEMPLATE(pair_end)

#endif
#if defined(PAIR_INCLUDE_HEADER) || defined(PAIR_INCLUDE_MODULE)

  interface TEMPLATE(pair_new)
    module procedure INTERNAL(pair_new_pair)
    module procedure INTERNAL(pair_new_init)
    module procedure INTERNAL(pair_new_copy)
  end interface TEMPLATE(pair_new)

  interface TEMPLATE(pair_associated)
    module procedure INTERNAL(pair_associated)
    module procedure INTERNAL(pair_pair_frst_associated)
    module procedure INTERNAL(pair_pair_scnd_associated)
    module procedure INTERNAL(pair_pair_associated)
  end interface TEMPLATE(pair_associated)

  interface TEMPLATE(pair_set)
    module procedure INTERNAL(pair_set_frst)
    module procedure INTERNAL(pair_set_scnd)
    module procedure INTERNAL(pair_set_both)
  end interface TEMPLATE(pair_set)

  interface TEMPLATE(pair_get)
    module procedure INTERNAL(pair_get_frst)
    module procedure INTERNAL(pair_get_scnd)
    module procedure INTERNAL(pair_get_both)
  end interface TEMPLATE(pair_get)

  type :: TEMPLATE(pair_t)
    private
    type(PAIR_FRST_TYPE_NAME), pointer :: frst => null()
    type(PAIR_SCND_TYPE_NAME), pointer :: scnd => null()
  end type TEMPLATE(pair_t)

#endif
#if defined(PAIR_INCLUDE_MODULE)

contains

#endif
#if defined(PAIR_INCLUDE_BODY) || defined(PAIR_INCLUDE_MODULE)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_new_pair)(this)
    type(TEMPLATE(pair_t)), pointer :: this

    PUSH_SUB(INTERNAL(pair_new_pair))

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(INTERNAL(pair_new_pair))
  end subroutine INTERNAL(pair_new_pair)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_new_init)(this, frst, scnd)
    type(TEMPLATE(pair_t)),   pointer     :: this
    type(PAIR_FRST_TYPE_NAME), intent(in) :: frst
    type(PAIR_SCND_TYPE_NAME), intent(in) :: scnd

    PUSH_SUB(INTERNAL(pair_new_init))

    call INTERNAL(pair_new_pair)(this)
    call TEMPLATE(pair_init)(this, frst, scnd)

    POP_SUB(INTERNAL(pair_new_init))
  end subroutine INTERNAL(pair_new_init)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_new_copy)(this, that)
    type(TEMPLATE(pair_t)), pointer     :: this
    type(TEMPLATE(pair_t)),  intent(in) :: that

    PUSH_SUB(INTERNAL(pair_new_copy))

    call INTERNAL(pair_new_pair)(this)
    call TEMPLATE(pair_copy)(this, that)

    POP_SUB(INTERNAL(pair_new_copy))
  end subroutine INTERNAL(pair_new_copy)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_del)(this)
    type(TEMPLATE(pair_t)), pointer :: this

    PUSH_SUB(TEMPLATE(pair_del))

    if(associated(this))then
      call TEMPLATE(pair_end)(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(TEMPLATE(pair_del))
  end subroutine TEMPLATE(pair_del)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_init)(this, frst, scnd)
    type(TEMPLATE(pair_t)),    intent(out) :: this
    type(PAIR_FRST_TYPE_NAME), intent(in)  :: frst
    type(PAIR_SCND_TYPE_NAME), intent(in)  :: scnd

    PUSH_SUB(TEMPLATE(pair_init))

    call INTERNAL(pair_set_both)(this, frst, scnd)

    POP_SUB(TEMPLATE(pair_init))
  end subroutine TEMPLATE(pair_init)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_associated)(this) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this

    logical :: eqv

    eqv = (associated(this%frst).and.associated(this%scnd))

  end function INTERNAL(pair_associated)
  
  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_frst_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)),            intent(in) :: this
    type(PAIR_FRST_TYPE_NAME), target, intent(in) :: that

    logical :: eqv

    eqv=.false.
    if(associated(this%frst)) eqv = associated(this%frst, that)

  end function INTERNAL(pair_pair_frst_associated)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_scnd_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)),            intent(in) :: this
    type(PAIR_SCND_TYPE_NAME), target, intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(associated(this%scnd)) eqv = associated(this%scnd, that)

  end function INTERNAL(pair_pair_scnd_associated)

  ! -----------------------------------------------------
  elemental function INTERNAL(pair_pair_associated)(this, that) result(eqv)
    type(TEMPLATE(pair_t)), intent(in) :: this
    type(TEMPLATE(pair_t)), intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(INTERNAL(pair_associated)(this).and.INTERNAL(pair_associated)(that))&
      eqv = (associated(this%frst, that%frst).and.associated(this%scnd, that%scnd))

  end function INTERNAL(pair_pair_associated)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_frst)(this, that)
    type(TEMPLATE(pair_t)),            intent(inout) :: this
    type(PAIR_FRST_TYPE_NAME), target, intent(in)    :: that

    PUSH_SUB(INTERNAL(pair_set_frst))

    nullify(this%frst)
    this%frst => that

    POP_SUB(INTERNAL(pair_set_frst))
  end subroutine INTERNAL(pair_set_frst)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_scnd)(this, that)
    type(TEMPLATE(pair_t)),            intent(inout) :: this
    type(PAIR_SCND_TYPE_NAME), target, intent(in)    :: that

    PUSH_SUB(INTERNAL(pair_set_scnd))

    nullify(this%scnd)
    ASSERT(associated(this%frst))
    this%scnd => that

    POP_SUB(INTERNAL(pair_set_scnd))
  end subroutine INTERNAL(pair_set_scnd)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_both)(this, frst, scnd)
    type(TEMPLATE(pair_t)),            intent(inout) :: this
    type(PAIR_FRST_TYPE_NAME), target, intent(in)    :: frst
    type(PAIR_SCND_TYPE_NAME), target, intent(in)    :: scnd

    PUSH_SUB(INTERNAL(pair_set_both))

    nullify(this%frst, this%scnd)
    this%frst => frst
    this%scnd => scnd

    POP_SUB(INTERNAL(pair_set_both))
  end subroutine INTERNAL(pair_set_both)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_frst)(this, that)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_FRST_TYPE_NAME), pointer     :: that

    PUSH_SUB(INTERNAL(pair_get_frst))

    nullify(that)
    if(associated(this%frst)) that => this%frst

    POP_SUB(INTERNAL(pair_get_frst))
  end subroutine INTERNAL(pair_get_frst)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_scnd)(this, that)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_SCND_TYPE_NAME), pointer     :: that

    PUSH_SUB(INTERNAL(pair_get_scnd))

    nullify(that)
    ASSERT(associated(this%frst))
    if(associated(this%scnd)) that => this%scnd

    POP_SUB(INTERNAL(pair_get_scnd))
  end subroutine INTERNAL(pair_get_scnd)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_both)(this, frst, scnd)
    type(TEMPLATE(pair_t)),     intent(in) :: this
    type(PAIR_FRST_TYPE_NAME), pointer     :: frst
    type(PAIR_SCND_TYPE_NAME), pointer     :: scnd

    PUSH_SUB(INTERNAL(pair_get_both))

    nullify(frst, scnd)
    if(associated(this%frst))then
      frst => this%frst
      if(associated(this%scnd)) scnd => this%scnd
    end if

    POP_SUB(INTERNAL(pair_get_both))
  end subroutine INTERNAL(pair_get_both)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_copy)(this, that)
    type(TEMPLATE(pair_t)), intent(inout) :: this
    type(TEMPLATE(pair_t)), intent(in)    :: that

    PUSH_SUB(TEMPLATE(pair_copy))

    call TEMPLATE(pair_end)(this)
    if(associated(that%frst))then
      if(associated(that%scnd))then
        call TEMPLATE(pair_set)(this, that%frst, that%scnd)
      else
        call TEMPLATE(pair_set)(this, that%frst)
      end if
    end if

    POP_SUB(TEMPLATE(pair_copy))
  end subroutine TEMPLATE(pair_copy)

  ! -----------------------------------------------------
  subroutine TEMPLATE(pair_end)(this)
    type(TEMPLATE(pair_t)), intent(inout) :: this

    PUSH_SUB(TEMPLATE(pair_end))

    nullify(this%frst, this%scnd)

    POP_SUB(TEMPLATE(pair_end))
  end subroutine TEMPLATE(pair_end)

#endif
#if defined(PAIR_INCLUDE_MODULE)

end module TEMPLATE(pair_m)

#endif

#undef IPAIR_TMPL_PREFIX
#undef IPAIR_TMPL_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:



