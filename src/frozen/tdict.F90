#include "global.h"

!DICT: DICT_TEMPLATE_NAME
!DICT: DICT_TYPE_NAME
!DICT: DICT_TYPE_MODULE_NAME

#if defined(DICT_TEMPLATE_NAME)
#if !defined(DICT_TYPE_NAME)
#define DICT_TYPE_NAME DECORATE(DICT_TEMPLATE_NAME,t)
#endif
#if !defined(DICT_TYPE_MODULE_NAME)
#define DICT_TYPE_MODULE_NAME DECORATE(DICT_TEMPLATE_NAME,m)
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

#define HASH_OK            DICT_OK
#define HASH_KEY_ERROR     DICT_KEY_ERROR
#define HASH_EMPTY_ERROR   DICT_EMPTY_ERROR 
#define hash_t             dict_t 
#define hash_iterator_t    dict_iterator_t

#define HASH_TEMPLATE_NAME DICT_TEMPLATE_NAME
#if defined(DICT_INITIAL_SIZE)
#define HASH_INITIAL_SIZE DICT_INITIAL_SIZE
#endif
#if defined(DICT_GROWTH_FACTOR)
#define HASH_GROWTH_FACTOR DICT_GROWTH_FACTOR
#endif
#define HASH_KEY_TEMPLATE_NAME string
#define HASH_VAL_TEMPLATE_NAME DICT_TEMPLATE_NAME
#define HASH_VAL_TYPE_NAME DICT_TYPE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"
 
#if !defined(DICT_INCLUDE_PREFIX)
#if !defined(DICT_INCLUDE_HEADER) && !defined(DICT_INCLUDE_BODY)

module TEMPLATE(dict_m)

  use global_m
  use messages_m
  use profiling_m

  use kinds_m,  only: wp
  use strng_m, only: operator(==), string_hash=>strng_hash

  use DICT_TYPE_MODULE_NAME, only: &
    DICT_TYPE_NAME

  use strng_m, only:                   &
    string_t       => strng_t,         &
    string_init    => strng_init,      &
    string_get     => strng_get,       &
    string_tolower => strng_tolower,   &
    string_copy    => strng_copy,      &
    string_end     => strng_end

  implicit none

  private
    
  public ::                     &
    TEMPLATE(DICT_OK),          &
    TEMPLATE(DICT_KEY_ERROR),   &
    TEMPLATE(DICT_EMPTY_ERROR)

  public ::                &
    TEMPLATE(dict_t),      &
    TEMPLATE(dict_len),    &
    TEMPLATE(dict_init),   &
    TEMPLATE(dict_next),   &
    TEMPLATE(dict_pop),    &
    TEMPLATE(dict_set),    &
    TEMPLATE(dict_get),    &
    TEMPLATE(dict_del),    &
    TEMPLATE(dict_extend), &
    TEMPLATE(dict_copy),   &
    TEMPLATE(dict_end)

  public ::                    &
    TEMPLATE(dict_iterator_t)

#endif
#if !defined(DICT_INCLUDE_BODY)
#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER
#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"

  interface TEMPLATE(dict_len)
    module procedure EXTERNAL(hash_len)
  end interface TEMPLATE(dict_len)

  interface TEMPLATE(dict_init)
    module procedure INTERNAL(hash_init)
    module procedure INTERNAL(hash_iterator_init)
  end interface TEMPLATE(dict_init)

  interface TEMPLATE(dict_pop)
    module procedure INTERNAL(dict_pop_item)
    module procedure INTERNAL(dict_pop_key)
    module procedure INTERNAL(dict_pop_val)
  end interface TEMPLATE(dict_pop)

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
#if !defined(DICT_INCLUDE_HEADER) && !defined(DICT_INCLUDE_BODY)

contains

#endif
#if !defined(DICT_INCLUDE_HEADER)
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY
#define TEMPLATE_PREFIX DICT_TEMPLATE_NAME
#include "template.h"

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_item)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(out)   :: key
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr
    !
    type(string_t), pointer :: str
    integer                 :: jerr
    !
    PUSH_SUB(INTERNAL(dict_pop_item))
    nullify(val, str)
    call EXTERNAL(hash_pop)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK))then
      call string_get(str, key)
      call string_end(str)
      SAFE_DEALLOCATE_P(str)
    end if
    nullify(str)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(dict_pop_item))
    return
  end subroutine INTERNAL(dict_pop_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_key)(this, key, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(out)   :: key
    integer,      optional, intent(out)   :: ierr
    !
    type(string_t), pointer :: str
    integer                 :: jerr
    !
    PUSH_SUB(INTERNAL(dict_pop_key))
    nullify(str)
    call EXTERNAL(hash_pop)(this, str, jerr)
    if(jerr==TEMPLATE(DICT_OK))then
      call string_get(str, key)
      call string_end(str)
      SAFE_DEALLOCATE_P(str)
    end if
    nullify(str)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(dict_pop_key))
    return
  end subroutine INTERNAL(dict_pop_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_pop_val)(this, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr
    !
    type(string_t), pointer :: str
    integer                 :: jerr
    !
    PUSH_SUB(INTERNAL(dict_pop_val))
    nullify(val, str)
    call EXTERNAL(hash_pop)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK))then
      call string_end(str)
      SAFE_DEALLOCATE_P(str)
    end if
    nullify(str)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(dict_pop_val))
    return
  end subroutine INTERNAL(dict_pop_val)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_set)(this, key, val)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key
    type(DICT_TYPE_NAME),   intent(in)    :: val
    !
    type(string_t), pointer :: str
    !
    PUSH_SUB(TEMPLATE(dict_set))
    nullify(str)
    SAFE_ALLOCATE(str)
    call string_init(str, key)
    call string_tolower(str)
    call EXTERNAL(hash_set)(this, str, val)
    nullify(str)
    POP_SUB(TEMPLATE(dict_set))
    return
  end subroutine TEMPLATE(dict_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_get)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(in)  :: this
    character(len=*),       intent(in)  :: key
    type(DICT_TYPE_NAME),  pointer      :: val
    integer,      optional, intent(out) :: ierr
    !
    type(string_t) :: str
    !
    PUSH_SUB(TEMPLATE(dict_get))
    nullify(val)
    call string_init(str, key)
    call string_tolower(str)
    call EXTERNAL(hash_get)(this, str, val, ierr)
    call string_end(str)
    POP_SUB(TEMPLATE(dict_get))
    return
  end subroutine TEMPLATE(dict_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_del)(this, key, val, ierr)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    character(len=*),       intent(in)    :: key
    type(DICT_TYPE_NAME),  pointer        :: val
    integer,      optional, intent(out)   :: ierr
    !
    type(string_t) :: str
    !
    PUSH_SUB(TEMPLATE(dict_del))
    nullify(val)
    call string_init(str, key)
    call string_tolower(str)
    call EXTERNAL(hash_del)(this, str, val, ierr)
    call string_end(str)
    POP_SUB(TEMPLATE(dict_del))
    return
  end subroutine TEMPLATE(dict_del)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(dict_extend)(this, that)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(TEMPLATE(dict_t)), intent(in)    :: that
    !
    type(TEMPLATE(dict_iterator_t)) :: iter
    type(string_t),         pointer :: ikey, okey
    type(DICT_TYPE_NAME),   pointer :: val
    integer                         :: ierr
    !
    PUSH_SUB(TEMPLATE(dict_extend))
    call TEMPLATE(dict_init)(iter, that)
    do
      nullify(ikey, okey, val)
      call EXTERNAL(hash_next)(iter, ikey, val, ierr)
      if(ierr/=TEMPLATE(DICT_OK))exit
      SAFE_ALLOCATE(okey)
      call string_init(okey, ikey)
      call EXTERNAL(hash_set)(this, okey, val)
    end do
    nullify(ikey, val)
    call TEMPLATE(dict_end)(iter)
    POP_SUB(TEMPLATE(dict_extend))
    return
  end subroutine TEMPLATE(dict_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_copy)(this, that)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    type(TEMPLATE(dict_t)), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(dict_copy))
    call INTERNAL(dict_purge)(this)
    call TEMPLATE(dict_extend)(this, that)
    POP_SUB(INTERNAL(dict_copy))
    return
  end subroutine INTERNAL(dict_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_purge)(this)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    !
    type(DICT_TYPE_NAME), pointer :: pval
    integer                       :: ierr
    !
    PUSH_SUB(INTERNAL(dict_purge))
    do
      nullify(pval)
      call INTERNAL(dict_pop_val)(this, pval, ierr)
      if(ierr/=TEMPLATE(DICT_OK))exit
    end do
    nullify(pval)
    POP_SUB(INTERNAL(dict_purge))
    return
  end subroutine INTERNAL(dict_purge)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_end)(this)
    type(TEMPLATE(dict_t)), intent(inout) :: this
    !
    PUSH_SUB(INTERNAL(dict_end))
    call INTERNAL(dict_purge)(this)
    call EXTERNAL(hash_end)(this)
    POP_SUB(INTERNAL(dict_end))
    return
  end subroutine INTERNAL(dict_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_item)(this, key, val, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: key
    type(DICT_TYPE_NAME),           pointer        :: val
    integer,               optional, intent(out)   :: ierr
    !
    type(string_t), pointer :: str
    integer                 :: jerr
    !
    PUSH_SUB(INTERNAL(dict_iterator_next_item))
    nullify(str)
    call EXTERNAL(hash_next)(this, str, val, jerr)
    if(jerr==TEMPLATE(DICT_OK))call string_get(str, key)
    nullify(str)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(dict_iterator_next_item))
    return
  end subroutine INTERNAL(dict_iterator_next_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_key)(this, that, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(string_t), pointer :: str
    integer                 :: jerr
    !
    PUSH_SUB(INTERNAL(dict_iterator_next_key))
    nullify(str)
    call EXTERNAL(hash_next)(this, str, jerr)
    if(jerr==TEMPLATE(DICT_OK))call string_get(str, that)
    nullify(str)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(dict_iterator_next_key))
    return
  end subroutine INTERNAL(dict_iterator_next_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(dict_iterator_next_val)(this, that, ierr)
    type(TEMPLATE(dict_iterator_t)), intent(inout) :: this
    type(DICT_TYPE_NAME),           pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(INTERNAL(dict_iterator_next_val))
    nullify(that)
    call EXTERNAL(hash_next)(this, that, ierr)
    POP_SUB(INTERNAL(dict_iterator_next_val))
    return
  end subroutine INTERNAL(dict_iterator_next_val)

#endif
#if !defined(DICT_INCLUDE_HEADER) && !defined(DICT_INCLUDE_BODY)

end module TEMPLATE(dict_m)

#endif
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

