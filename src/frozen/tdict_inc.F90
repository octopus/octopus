#include "global.h"

!DICT: DICT_TEMPLATE_NAME
!DICT: DICT_TYPE_NAME
!DICT: DICT_TYPE_MODULE_NAME

#if defined(DICT_TEMPLATE_NAME)
#if !defined(DICT_TYPE_NAME)
#define DICT_TYPE_NAME DECORATE(DICT_TEMPLATE_NAME,t)
#endif
#if !defined(DICT_TYPE_MODULE_NAME)
#define DICT_TYPE_MODULE_NAME DECORATE(DICT_TEMPLATE_NAME, oct_m)
#endif
#else
#error "'DICT_TEMPLATE_NAME' must be defined"
#endif

#undef HASH_TEMPLATE_NAME
#undef HASH_INITIAL_SIZE
#undef HASH_GROWTH_FACTOR
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_KEY_TYPE_MODULE_NAME
#undef HASH_KEY_INCLUDE
#undef HASH_KEY_FUNCTION_NAME
#undef HASH_KEY_FUNCTION_MODULE_NAME
#undef HASH_VAL_TEMPLATE_NAME
#undef HASH_VAL_TYPE_NAME
#undef HASH_VAL_TYPE_MODULE_NAME
#undef HASH_VAL_INCLUDE
#undef HASH_INCLUDE_PREFIX
#undef HASH_INCLUDE_HEADER
#undef HASH_INCLUDE_BODY

#define hash_t             dict_t 
#define hash_iterator_t    dict_iterator_t

#define HASH_TEMPLATE_NAME DICT_TEMPLATE_NAME
#if defined(DICT_INITIAL_SIZE)
#define HASH_INITIAL_SIZE DICT_INITIAL_SIZE
#endif
#if defined(DICT_GROWTH_FACTOR)
#define HASH_GROWTH_FACTOR DICT_GROWTH_FACTOR
#endif
#define HASH_KEY_TEMPLATE_NAME strng
#define HASH_VAL_TEMPLATE_NAME DICT_TEMPLATE_NAME
#define HASH_VAL_TYPE_NAME DICT_TYPE_NAME

#undef DICT_INCLUDE_MODULE
#if !defined(DICT_INCLUDE_PREFIX) && !defined(DICT_INCLUDE_HEADER) && !defined(DICT_INCLUDE_BODY)
#define DICT_INCLUDE_MODULE
#endif

#if defined(DICT_INCLUDE_PREFIX) && defined(DICT_INCLUDE_HEADER)
#error "Only one off 'DICT_INCLUDE_PREFIX' or 'DICT_INCLUDE_HEADER' can be defined."
#endif

#if defined(DICT_INCLUDE_PREFIX) && defined(DICT_INCLUDE_BODY)
#error "Only one off 'DICT_INCLUDE_PREFIX' or 'DICT_INCLUDE_BODY' can be defined."
#endif

#if defined(DICT_INCLUDE_HEADER) && defined(DICT_INCLUDE_BODY)
#error "Only one off 'DICT_INCLUDE_HEADER' or 'DICT_INCLUDE_BODY' can be defined."
#endif

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"
 
#if defined(DICT_INCLUDE_MODULE)

module TEMPLATE(dict_oct_m)

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

#endif
#if defined(DICT_INCLUDE_PREFIX) || defined(DICT_INCLUDE_MODULE)

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"
 
  use strng_oct_m

#if defined(DICT_TYPE_EXTERNAL) || defined(DICT_INCLUDE_MODULE)

  use DICT_TYPE_MODULE_NAME

#endif

#endif
#if defined(DICT_INCLUDE_MODULE)

  implicit none

  private

  public ::                  &
    TEMPLATE(DICT_NAME_LEN)

  public ::                     &
    TEMPLATE(DICT_OK),          &
    TEMPLATE(DICT_KEY_ERROR),   &
    TEMPLATE(DICT_EMPTY_ERROR)

  public ::           &
    TEMPLATE(dict_t)

  public ::                    &
    TEMPLATE(dict_iterator_t)

  public ::                 &
    TEMPLATE(dict_len),     &
    TEMPLATE(dict_has_key), &
    TEMPLATE(dict_init),    &
    TEMPLATE(dict_next),    &
    TEMPLATE(dict_pop),     &
    TEMPLATE(dict_set),     &
    TEMPLATE(dict_get),     &
    TEMPLATE(dict_del),     &
    TEMPLATE(dict_extend),  &
    TEMPLATE(dict_copy),    &
    TEMPLATE(dict_end)

#endif
#if defined(DICT_INCLUDE_HEADER) || defined(DICT_INCLUDE_MODULE)

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"

  integer, parameter :: TEMPLATE(DICT_NAME_LEN) = 63

  integer, parameter :: TEMPLATE(DICT_OK)          = TEMPLATE(HASH_OK)
  integer, parameter :: TEMPLATE(DICT_KEY_ERROR)   = TEMPLATE(HASH_KEY_ERROR)
  integer, parameter :: TEMPLATE(DICT_EMPTY_ERROR) = TEMPLATE(HASH_EMPTY_ERROR)

  interface TEMPLATE(dict_len)
    module procedure TEMPLATE(hash_len)
  end interface TEMPLATE(dict_len)

  interface TEMPLATE(dict_init)
    module procedure INTERNAL(hash_init)
    module procedure INTERNAL(hash_iterator_init_hash)
    module procedure INTERNAL(hash_iterator_init_iterator)
  end interface TEMPLATE(dict_init)

  interface TEMPLATE(dict_pop)
    module procedure INTERNAL(dict_pop_item)
    module procedure INTERNAL(dict_pop_key)
    module procedure INTERNAL(dict_pop_val)
  end interface TEMPLATE(dict_pop)

  interface TEMPLATE(dict_del)
    module procedure INTERNAL(dict_del_empty)
    module procedure INTERNAL(dict_del_val)
  end interface TEMPLATE(dict_del)

  interface TEMPLATE(dict_next)
    module procedure INTERNAL(dict_iterator_next_item)
    module procedure INTERNAL(dict_iterator_next_key)
    module procedure INTERNAL(dict_iterator_next_val)
  end interface TEMPLATE(dict_next)

  interface TEMPLATE(dict_copy)
    module procedure INTERNAL(dict_copy)
    module procedure INTERNAL(hash_iterator_copy)
  end interface TEMPLATE(dict_copy)

  interface TEMPLATE(dict_end)
    module procedure INTERNAL(dict_end)
    module procedure INTERNAL(hash_iterator_end)
  end interface TEMPLATE(dict_end)

#endif
#if defined(DICT_INCLUDE_MODULE)

contains

#endif
#if defined(DICT_INCLUDE_BODY) || defined(DICT_INCLUDE_MODULE)

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_item)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(out)   :: key
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr

    type(strng_t), pointer :: str
    integer                :: jerr

    PUSH_SUB(INTERNAL(dict_pop_item))

    nullify(val, str)
    call TEMPLATE(hash_pop)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK))then
      call strng_get(str, key)
      call strng_del(str)
    end if
    nullify(str)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_pop_item))
  end subroutine INTERNAL(dict_pop_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_key)(this, key, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(out)   :: key
    integer,      optional, intent(out)   :: ierr

    type(strng_t), pointer :: str
    integer                :: jerr

    PUSH_SUB(INTERNAL(dict_pop_key))

    nullify(str)
    call TEMPLATE(hash_pop)(this, str, jerr)
    if(jerr==TEMPLATE(DICT_OK))then
      call strng_get(str, key)
      call strng_del(str)
    end if
    nullify(str)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_pop_key))
  end subroutine INTERNAL(dict_pop_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_val)(this, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr
 
    type(strng_t), pointer :: str
    integer                :: jerr
 
    PUSH_SUB(INTERNAL(dict_pop_val))
 
    nullify(val, str)
    call TEMPLATE(hash_pop)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK)) call strng_del(str)
    nullify(str)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_pop_val))
  end subroutine INTERNAL(dict_pop_val)

  ! ---------------------------------------------------------
  function TEMPLATE(dict_has_key)(this, key) result(has)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key

    logical :: has

    type(strng_t) :: str

    PUSH_SUB(TEMPLATE(dict_has_key))

    call strng_init(str, key)
    call strng_tolower(str)
    has = TEMPLATE(hash_has_key)(this, str)
    call strng_end(str)

    POP_SUB(TEMPLATE(dict_has_key))
  end function TEMPLATE(dict_has_key)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_set)(this, key, val)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key
    type(DICT_TYPE_NAME),   intent(in)    :: val

    type(strng_t), pointer :: str
    logical                :: has

    PUSH_SUB(TEMPLATE(dict_set))

    nullify(str)
    call strng_new(str, key)
    call strng_tolower(str)
    has = TEMPLATE(hash_has_key)(this, str)
    call TEMPLATE(hash_set)(this, str, val)
    if(has) call strng_del(str)
    nullify(str)

    POP_SUB(TEMPLATE(dict_set))
  end subroutine TEMPLATE(dict_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_get)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(in)  :: this
    character(len=*),       intent(in)  :: key
    type(DICT_TYPE_NAME),  pointer      :: val
    integer,      optional, intent(out) :: ierr

    type(strng_t) :: str

    PUSH_SUB(TEMPLATE(dict_get))

    nullify(val)
    call strng_init(str, key)
    call strng_tolower(str)
    call TEMPLATE(hash_get)(this, str, val, ierr)
    call strng_end(str)

    POP_SUB(TEMPLATE(dict_get))
  end subroutine TEMPLATE(dict_get)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_del_val)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr

    type(strng_t), pointer :: pkey
    type(strng_t)          :: str
    integer                :: jerr

    PUSH_SUB(INTERNAL(dict_del_val))

    nullify(pkey, val)
    call strng_init(str, key)
    call strng_tolower(str)
    call TEMPLATE(hash_del)(this, str, pkey, val, jerr)
    if(jerr==TEMPLATE(DICT_OK))call strng_del(pkey)
    call strng_end(str)
    nullify(pkey)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_del_val))
  end subroutine INTERNAL(dict_del_val)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_del_empty)(this, key, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key
    integer,      optional, intent(out)   :: ierr

    type(DICT_TYPE_NAME), pointer :: pval

    PUSH_SUB(INTERNAL(dict_del_empty))

    nullify(pval)
    call INTERNAL(dict_del_val)(this, key, pval, ierr)
    nullify(pval)

    POP_SUB(INTERNAL(dict_del_empty))
  end subroutine INTERNAL(dict_del_empty)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_extend)(this, that)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(TEMPLATE(dict_t)), intent(in)    :: that

    type(TEMPLATE(dict_iterator_t)) :: iter
    type(strng_t),          pointer :: ikey, okey
    type(DICT_TYPE_NAME),   pointer :: val
    integer                         :: ierr

    PUSH_SUB(TEMPLATE(dict_extend))

    call TEMPLATE(dict_init)(iter, that)
    do
      nullify(ikey, okey, val)
      call TEMPLATE(hash_next)(iter, ikey, val, ierr)
      if(ierr/=TEMPLATE(DICT_OK))exit
      call strng_new(okey, ikey)
      call TEMPLATE(hash_set)(this, okey, val)
    end do
    nullify(ikey, val)
    call TEMPLATE(dict_end)(iter)

    POP_SUB(TEMPLATE(dict_extend))
  end subroutine TEMPLATE(dict_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_copy)(this, that)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(TEMPLATE(dict_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(dict_copy))

    call INTERNAL(dict_purge)(this)
    call TEMPLATE(dict_extend)(this, that)

    POP_SUB(INTERNAL(dict_copy))
  end subroutine INTERNAL(dict_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_purge)(this)
    type(TEMPLATE(dict_t)), intent(inout) :: this

    type(DICT_TYPE_NAME), pointer :: pval
    integer                       :: ierr

    PUSH_SUB(INTERNAL(dict_purge))

    do
      nullify(pval)
      call INTERNAL(dict_pop_val)(this, pval, ierr)
      if(ierr/=TEMPLATE(DICT_OK))exit
    end do
    nullify(pval)

    POP_SUB(INTERNAL(dict_purge))
  end subroutine INTERNAL(dict_purge)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_end)(this)
    type(TEMPLATE(dict_t)), intent(inout) :: this

    PUSH_SUB(INTERNAL(dict_end))

    call INTERNAL(dict_purge)(this)
    call TEMPLATE(hash_end)(this)

    POP_SUB(INTERNAL(dict_end))
  end subroutine INTERNAL(dict_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_item)(this, key, val, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: key
    type(DICT_TYPE_NAME),           pointer        :: val
    integer,               optional, intent(out)   :: ierr

    type(strng_t), pointer :: str
    integer                :: jerr

    PUSH_SUB(INTERNAL(dict_iterator_next_item))

    nullify(str)
    call TEMPLATE(hash_next)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK)) call strng_get(str, key)
    nullify(str)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_iterator_next_item))
  end subroutine INTERNAL(dict_iterator_next_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_key)(this, that, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: that
    integer,               optional, intent(out)   :: ierr

    type(strng_t), pointer :: str
    integer                :: jerr

    PUSH_SUB(INTERNAL(dict_iterator_next_key))

    nullify(str)
    call TEMPLATE(hash_next)(this, str, jerr)
    if(jerr==TEMPLATE(DICT_OK)) call strng_get(str, that)
    nullify(str)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(dict_iterator_next_key))
  end subroutine INTERNAL(dict_iterator_next_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_val)(this, that, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    type(DICT_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr

    PUSH_SUB(INTERNAL(dict_iterator_next_val))

    nullify(that)
    call TEMPLATE(hash_next)(this, that, ierr)

    POP_SUB(INTERNAL(dict_iterator_next_val))
  end subroutine INTERNAL(dict_iterator_next_val)

#endif
#if defined(DICT_INCLUDE_MODULE)

end module TEMPLATE(dict_oct_m)

#endif

#undef HASH_OK
#undef HASH_KEY_ERROR

#undef HASH_EMPTY_ERROR
#undef hash_t
#undef hash_iterator_t

#undef HASH_TEMPLATE_NAME
#undef HASH_INITIAL_SIZE
#undef HASH_GROWTH_FACTOR
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_KEY_TYPE_MODULE_NAME
#undef HASH_KEY_INCLUDE
#undef HASH_KEY_FUNCTION_NAME
#undef HASH_KEY_FUNCTION_MODULE_NAME
#undef HASH_VAL_TEMPLATE_NAME
#undef HASH_VAL_TYPE_NAME
#undef HASH_VAL_TYPE_MODULE_NAME
#undef HASH_VAL_INCLUDE

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:

