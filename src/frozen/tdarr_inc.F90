#include "global.h"
#include "util.h"

!DARR: DARR_TEMPLATE_NAME
!DARR: DARR_TYPE_NAME
!DARR: DARR_TYPE_MODULE_NAME
!DARR: DARR_INITIAL_SIZE
!DARR: DARR_GROWTH_FACTOR

#if defined(DARR_TEMPLATE_NAME)
#if !defined(DARR_TYPE_NAME)
#define DARR_TYPE_NAME DECORATE(DARR_TEMPLATE_NAME,t)
#endif
#if !defined(DARR_TYPE_MODULE_NAME)
#define DARR_TYPE_MODULE_NAME DECORATE(DARR_TEMPLATE_NAME,m)
#endif
#else
#error "'DARR_TEMPLATE_NAME' must be defined"
#endif

#if !defined(DARR_INITIAL_SIZE)
#define DARR_INITIAL_SIZE 63
#endif
#if !defined(DARR_GROWTH_FACTOR)
#define DARR_GROWTH_FACTOR 1.1
#endif

#undef DARR_INCLUDE_MODULE
#if !defined(DARR_INCLUDE_PREFIX) && !defined(DARR_INCLUDE_HEADER) && !defined(DARR_INCLUDE_BODY)
#define DARR_INCLUDE_MODULE
#endif

#if defined(DARR_INCLUDE_PREFIX) && defined(DARR_INCLUDE_HEADER)
#error "Only one off 'DARR_INCLUDE_PREFIX' or 'DARR_INCLUDE_HEADER' can be defined."
#endif

#if defined(DARR_INCLUDE_PREFIX) && defined(DARR_INCLUDE_BODY)
#error "Only one off 'DARR_INCLUDE_PREFIX' or 'DARR_INCLUDE_BODY' can be defined."
#endif

#if defined(DARR_INCLUDE_HEADER) && defined(DARR_INCLUDE_BODY)
#error "Only one off 'DARR_INCLUDE_HEADER' or 'DARR_INCLUDE_BODY' can be defined."
#endif

#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_TYPE_MODULE_NAME
#undef SINGLE_INCLUDE_PREFIX
#undef SINGLE_INCLUDE_HEADER
#undef SINGLE_INCLUDE_BODY

#define SINGLE_TEMPLATE_NAME DARR_TEMPLATE_NAME
#define SINGLE_TYPE_NAME DARR_TYPE_NAME

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX DARR_TEMPLATE_NAME
#include "template.h"

#if defined(DARR_INCLUDE_MODULE)

module TEMPLATE(darr_m)

  use global_m
  use messages_m
  use profiling_m

#endif
#if defined(DARR_INCLUDE_PREFIX) || defined(DARR_INCLUDE_MODULE)

#define SINGLE_INCLUDE_PREFIX
#include "tsingle_inc.F90"
#undef SINGLE_INCLUDE_PREFIX

  use kinds_m

#endif
#if defined(DARR_INCLUDE_MODULE)

  use DARR_TYPE_MODULE_NAME

  implicit none

  private

  public ::                     &
    TEMPLATE(DARR_OK),          &
    TEMPLATE(DARR_INDEX_ERROR), &
    TEMPLATE(DARR_EMPTY_ERROR)

  public ::           &
    TEMPLATE(darr_t)

  public ::                    &
    TEMPLATE(darr_iterator_t)

  public ::                 &
    TEMPLATE(darr_size),    &
    TEMPLATE(darr_init),    &
    TEMPLATE(darr_realloc), &
    TEMPLATE(darr_next),    &
    TEMPLATE(darr_pop),     &
    TEMPLATE(darr_get),     &
    TEMPLATE(darr_set),     &
    TEMPLATE(darr_append),  &
    TEMPLATE(darr_extend),  &
    TEMPLATE(darr_copy),    &
    TEMPLATE(darr_end)

#endif
#if defined(DARR_INCLUDE_HEADER) || defined(DARR_INCLUDE_MODULE)

#define SINGLE_INCLUDE_HEADER
#include "tsingle_inc.F90"
#undef SINGLE_INCLUDE_HEADER

#define TEMPLATE_PREFIX DARR_TEMPLATE_NAME
#include "template.h"

  real(kind=wp), parameter :: INTERNAL(DARR_FACTOR) = DECORATE(DARR_GROWTH_FACTOR,wp)
  integer,       parameter :: INTERNAL(DARR_SIZE)   = DARR_INITIAL_SIZE

  integer, parameter :: TEMPLATE(DARR_OK)          = 0
  integer, parameter :: TEMPLATE(DARR_INDEX_ERROR) =-1
  integer, parameter :: TEMPLATE(DARR_EMPTY_ERROR) =-2

  type :: TEMPLATE(darr_t)
    private
    integer                                         :: size = 0
    integer                                         :: used = 0
    type(TEMPLATE(single_t)), dimension(:), pointer :: data =>null()
  end type TEMPLATE(darr_t)

  type :: TEMPLATE(darr_iterator_t)
    private
    integer                         :: ipos = 0
    type(TEMPLATE(darr_t)), pointer :: darr =>null()
  end type TEMPLATE(darr_iterator_t)

  interface TEMPLATE(darr_init)
    module procedure INTERNAL(darr_init)
    module procedure INTERNAL(darr_iterator_init_darr)
    module procedure INTERNAL(darr_iterator_init_iterator)
  end interface TEMPLATE(darr_init)

  interface TEMPLATE(darr_next)
    module procedure INTERNAL(darr_iterator_next)
  end interface TEMPLATE(darr_next)

  interface TEMPLATE(darr_copy)
    module procedure INTERNAL(darr_copy)
    module procedure INTERNAL(darr_iterator_copy)
  end interface TEMPLATE(darr_copy)

  interface TEMPLATE(darr_end)
    module procedure INTERNAL(darr_end)
    module procedure INTERNAL(darr_iterator_end)
  end interface TEMPLATE(darr_end)

#endif
#if defined(DARR_INCLUDE_MODULE)

contains
  
#endif
#if defined(DARR_INCLUDE_BODY) || defined(DARR_INCLUDE_MODULE)

#define SINGLE_INCLUDE_BODY
#include "tsingle_inc.F90"
#undef SINGLE_INCLUDE_BODY

#define TEMPLATE_PREFIX DARR_TEMPLATE_NAME
#include "template.h"

  ! ---------------------------------------------------------
  elemental function TEMPLATE(darr_size)(this) result(that)
    type(TEMPLATE(darr_t)), intent(in) :: this

    integer :: that

    that = this%used

  end function TEMPLATE(darr_size)

  ! ---------------------------------------------------------
  elemental function INTERNAL(darr_calc_size)(this) result(that)
    integer, intent(in) :: this

    integer :: that

    real(kind=wp) :: tmp
    integer       :: pwr

    that=INTERNAL(DARR_SIZE)
    if(INTERNAL(DARR_SIZE)<this+1)then
      tmp = log(real(this,kind=wp)) - log(real(INTERNAL(DARR_SIZE),kind=wp))
      tmp = tmp / log(INTERNAL(DARR_FACTOR))
      pwr = max(ceiling(tmp), 1)
      that = ceiling((INTERNAL(DARR_FACTOR)**pwr)*real(INTERNAL(DARR_SIZE),kind=wp))
    end if

  end function INTERNAL(darr_calc_size)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_allocate)(this, that)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    integer,                  intent(in)    :: that

    PUSH_SUB(INTERNAL(darr_allocate))

    this%size = that
    SAFE_ALLOCATE(this%data(this%size))

    POP_SUB(INTERNAL(darr_allocate))
  end subroutine INTERNAL(darr_allocate)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_realloc)(this, that)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    integer,                intent(in)    :: that

    type(TEMPLATE(single_t)), dimension(:), pointer :: buff
    integer                                         :: indx, size, used

    PUSH_SUB(TEMPLATE(darr_realloc))

    nullify(buff)
    size = INTERNAL(darr_calc_size)(that)
    if(associated(this%data))then
      if(this%size/=size)then
        used = this%used
        buff => this%data
        nullify(this%data)
        call INTERNAL(darr_allocate)(this, size)
        this%used = min(used, size)
        do indx=1, this%used
          call TEMPLATE(single_copy)(this%data(indx), buff(indx))
        end do
        do indx = this%used+1, used
          call TEMPLATE(single_end)(buff(indx))
        end do
        nullify(buff)
      end if
    else
      call INTERNAL(darr_allocate)(this, size)
    end if

    POP_SUB(TEMPLATE(darr_realloc))
  end subroutine TEMPLATE(darr_realloc)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_init)(this, that)
    type(TEMPLATE(darr_t)), intent(out) :: this
    integer,      optional, intent(in)  :: that

    PUSH_SUB(INTERNAL(darr_init))

    call INTERNAL(darr_end)(this)
    ASSERT(INTERNAL(DARR_SIZE)>1)
    ASSERT(INTERNAL(DARR_FACTOR)>1.0_wp)
    this%used = that
    call TEMPLATE(darr_realloc)(this, that)

    POP_SUB(INTERNAL(darr_init))
  end subroutine INTERNAL(darr_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_pop)(this, value, ierr)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    type(DARR_TYPE_NAME),  pointer        :: value
    integer,      optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(darr_pop))

    nullify(value)
    if(present(ierr)) ierr = TEMPLATE(DARR_EMPTY_ERROR)
    if(this%used>0)then
      call TEMPLATE(single_get)(this%data(this%used), value)
      call TEMPLATE(single_end)(this%data(this%used))
      this%used = this%used - 1
      if(present(ierr)) ierr = TEMPLATE(DARR_OK)
    end if

    POP_SUB(TEMPLATE(darr_pop))
  end subroutine TEMPLATE(darr_pop)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_get)(this, index, value, ierr)
    type(TEMPLATE(darr_t)), intent(in)  :: this
    integer,                intent(in)  :: index
    type(DARR_TYPE_NAME),  pointer      :: value
    integer,      optional, intent(out) :: ierr

    PUSH_SUB(TEMPLATE(darr_get))

    nullify(value)
    if(present(ierr)) ierr = TEMPLATE(DARR_INDEX_ERROR)
    if((0<index).and.(index<=this%used))then
      call TEMPLATE(single_get)(this%data(index), value)
      if(present(ierr)) ierr = TEMPLATE(DARR_OK)
    end if

    POP_SUB(TEMPLATE(darr_get))
  end subroutine TEMPLATE(darr_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_set)(this, index, value, ierr)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    integer,                  intent(in)    :: index
    type(DARR_TYPE_NAME),  pointer        :: value
    integer,        optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(darr_set))

    if(present(ierr)) ierr = TEMPLATE(DARR_INDEX_ERROR)
    if((0<index).and.(index<=this%used))then
      call TEMPLATE(single_set)(this%data(index), value)
      if(present(ierr)) ierr = TEMPLATE(DARR_OK)
    end if

    POP_SUB(TEMPLATE(darr_set))
  end subroutine TEMPLATE(darr_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_append)(this, that)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    type(DARR_TYPE_NAME),  pointer        :: that

    PUSH_SUB(TEMPLATE(darr_append))

    call TEMPLATE(darr_realloc)(this, this%used+1)
    this%used = this%used + 1
    call TEMPLATE(single_init)(this%data(this%used), that)

    POP_SUB(TEMPLATE(darr_append))
  end subroutine TEMPLATE(darr_append)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(darr_extend)(this, that)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    type(TEMPLATE(darr_t)), intent(in)    :: that

    integer :: indx

    PUSH_SUB(TEMPLATE(darr_extend))

    call TEMPLATE(darr_realloc)(this, this%used+that%used)
    do indx = 1, that%used
      call TEMPLATE(single_copy)(this%data(this%used+indx), that%data(indx))
    end do
    this%used = this%used + that%used

    POP_SUB(TEMPLATE(darr_extend))
  end subroutine TEMPLATE(darr_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_copy)(this, that)
    type(TEMPLATE(darr_t)), intent(inout) :: this
    type(TEMPLATE(darr_t)), intent(in)    :: that

    integer :: indx

    PUSH_SUB(INTERNAL(darr_copy))

    call TEMPLATE(darr_realloc)(this, that%used)
    do indx = 1, that%used
      call TEMPLATE(single_copy)(this%data(indx), that%data(indx))
    end do
    do indx = that%used+1, this%used
      call TEMPLATE(single_end)(this%data(indx))
    end do
    this%used = that%used

    POP_SUB(INTERNAL(darr_copy))
  end subroutine INTERNAL(darr_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_end)(this)
    type(TEMPLATE(darr_t)), intent(inout) :: this

    type(DARR_TYPE_NAME), pointer :: value
    integer                       :: ierr

    PUSH_SUB(INTERNAL(darr_end))

    do
      nullify(value)
      call TEMPLATE(darr_pop)(this, value, ierr)
      if(ierr/=TEMPLATE(DARR_OK))exit
    end do
    nullify(value)
    ASSERT(this%used==0)
    SAFE_DEALLOCATE_P(this%data)
    nullify(this%data)
    this%size = 0

    POP_SUB(INTERNAL(darr_end))
  end subroutine INTERNAL(darr_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_iterator_init_darr)(this, that)
    type(TEMPLATE(darr_iterator_t)), intent(out) :: this
    type(TEMPLATE(darr_t)),  target, intent(in)  :: that

    PUSH_SUB(INTERNAL(darr_iterator_init_darr))

    this%ipos = 0
    this%darr => that

    POP_SUB(INTERNAL(darr_iterator_init_darr))
  end subroutine INTERNAL(darr_iterator_init_darr)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_iterator_init_iterator)(this, that)
    type(TEMPLATE(darr_iterator_t)), intent(out) :: this
    type(TEMPLATE(darr_iterator_t)), intent(in)  :: that

    PUSH_SUB(INTERNAL(darr_iterator_init_iterator))

    call INTERNAL(darr_iterator_copy)(this, that)

    POP_SUB(INTERNAL(darr_iterator_init_iterator))
  end subroutine INTERNAL(darr_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_iterator_next)(this, that, ierr)
    type(TEMPLATE(darr_iterator_t)), intent(inout) :: this
    type(DARR_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr

    PUSH_SUB(INTERNAL(darr_iterator_next))

    nullify(that)
    ASSERT(associated(this%darr))
    if(present(ierr)) ierr = TEMPLATE(DARR_EMPTY_ERROR)
    if(this%ipos<this%darr%used)then
      this%ipos = this%ipos + 1
      call TEMPLATE(single_get)(this%darr%data(this%ipos), that)
      if(present(ierr)) ierr = TEMPLATE(DARR_OK)
    end if

    POP_SUB(INTERNAL(darr_iterator_next))
  end subroutine INTERNAL(darr_iterator_next)

  ! ---------------------------------------------------------
  subroutine INTERNAL(darr_iterator_copy)(this, that)
    type(TEMPLATE(darr_iterator_t)), intent(inout) :: this
    type(TEMPLATE(darr_iterator_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(darr_iterator_copy))

    call INTERNAL(darr_iterator_end)(this)
    this%ipos = that%ipos
    this%darr => that%darr

    POP_SUB(INTERNAL(darr_iterator_copy))
  end subroutine INTERNAL(darr_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(darr_iterator_end)(this)
    type(TEMPLATE(darr_iterator_t)), intent(inout) :: this

    this%ipos = 0
    nullify(this%darr)

  end subroutine INTERNAL(darr_iterator_end)

#endif
#if defined(DARR_INCLUDE_MODULE)

end module TEMPLATE(darr_m)

#endif

#undef SINGLE_TEMPLATE_NAME
#undef SINGLE_TYPE_NAME
#undef SINGLE_TYPE_MODULE_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:

