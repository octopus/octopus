#include "global.h"
#include "util.h"

!HASH: HASH_TEMPLATE_NAME
!HASH: HASH_INITIAL_SIZE
!HASH: HASH_GROWTH_FACTOR
!HASH: HASH_KEY_TEMPLATE_NAME
!HASH: HASH_KEY_TYPE_NAME
!HASH: HASH_KEY_TYPE_MODULE_NAME
!HASH: HASH_KEY_FUNCTION_NAME
!HASH: HASH_KEY_FUNCTION_MODULE_NAME
!HASH: HASH_VAL_TEMPLATE_NAME
!HASH: HASH_VAL_TYPE_NAME
!HASH: HASH_VAL_TYPE_MODULE_NAME

#if defined(HASH_KEY_TEMPLATE_NAME)
#if !defined(HASH_KEY_TYPE_NAME)
#define HASH_KEY_TYPE_NAME DECORATE(HASH_KEY_TEMPLATE_NAME,t)
#endif
#if !defined(HASH_KEY_TYPE_MODULE_NAME)
#define HASH_KEY_TYPE_MODULE_NAME DECORATE(HASH_KEY_TEMPLATE_NAME,m)
#endif
#if !defined(HASH_KEY_FUNCTION_NAME)
#define HASH_KEY_FUNCTION_NAME DECORATE(HASH_KEY_TEMPLATE_NAME,hash)
#endif
#if !defined(HASH_KEY_FUNCTION_MODULE_NAME)
#define HASH_KEY_FUNCTION_MODULE_NAME DECORATE(HASH_KEY_TEMPLATE_NAME,m)
#endif
#else
#error "'HASH_KEY_TEMPLATE_NAME' must be defined"
#endif

#if defined(HASH_VAL_TEMPLATE_NAME)
#if !defined(HASH_VAL_TYPE_NAME)
#define HASH_VAL_TYPE_NAME DECORATE(HASH_VAL_TEMPLATE_NAME,t)
#endif
#if !defined(HASH_VAL_TYPE_MODULE_NAME)
#define HASH_VAL_TYPE_MODULE_NAME DECORATE(HASH_VAL_TEMPLATE_NAME,m)
#endif
#else
#error "'HASH_VAL_TEMPLATE_NAME' must be defined"
#endif

#undef IHASH_TMPL_NAME
#if defined(HASH_TEMPLATE_NAME)
#define IHASH_TMPL_NAME HASH_TEMPLATE_NAME
#else
#define IHASH_TMPL_NAME DECORATE(HASH_KEY_TEMPLATE_NAME,HASH_VAL_TEMPLATE_NAME)
#endif

#if !defined(HASH_INITIAL_SIZE)
#define HASH_INITIAL_SIZE 7
#endif
#if !defined(HASH_GROWTH_FACTOR)
#define HASH_GROWTH_FACTOR 1.5
#endif

#undef HASH_INCLUDE_MODULE
#if !defined(HASH_INCLUDE_PREFIX) && !defined(HASH_INCLUDE_HEADER) && !defined(HASH_INCLUDE_BODY)
#define HASH_INCLUDE_MODULE
#endif

#if defined(HASH_INCLUDE_PREFIX) && defined(HASH_INCLUDE_HEADER)
#error "Only one off 'HASH_INCLUDE_PREFIX' or 'HASH_INCLUDE_HEADER' can be defined."
#endif

#if defined(HASH_INCLUDE_PREFIX) && defined(HASH_INCLUDE_BODY)
#error "Only one off 'HASH_INCLUDE_PREFIX' or 'HASH_INCLUDE_BODY' can be defined."
#endif

#if defined(HASH_INCLUDE_HEADER) && defined(HASH_INCLUDE_BODY)
#error "Only one off 'HASH_INCLUDE_HEADER' or 'HASH_INCLUDE_BODY' can be defined."
#endif

#undef PAIR_TEMPLATE_NAME
#undef PAIR_FRST_TEMPLATE_NAME
#undef PAIR_FRST_TYPE_NAME
#undef PAIR_FRST_TYPE_MODULE_NAME
#undef PAIR_SCND_TEMPLATE_NAME
#undef PAIR_SCND_TYPE_NAME
#undef PAIR_SCND_TYPE_MODULE_NAME
#undef PAIR_INCLUDE_PREFIX
#undef PAIR_INCLUDE_HEADER
#undef PAIR_INCLUDE_BODY

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME
#undef LIST_INCLUDE_PREFIX
#undef LIST_INCLUDE_HEADER
#undef LIST_INCLUDE_BODY

#undef DARR_TEMPLATE_NAME
#undef DARR_TYPE_NAME
#undef DARR_TYPE_MODULE_NAME
#undef DARR_INCLUDE_PREFIX
#undef DARR_INCLUDE_HEADER
#undef DARR_INCLUDE_BODY

#define PAIR_TEMPLATE_NAME DECORATE(IHASH_TMPL_NAME,hash)
#define PAIR_FRST_TYPE_NAME HASH_KEY_TYPE_NAME
#define PAIR_SCND_TYPE_NAME HASH_VAL_TYPE_NAME

#define LIST_TEMPLATE_NAME DECORATE(IHASH_TMPL_NAME,hash_pair)

#define DARR_TEMPLATE_NAME DECORATE(IHASH_TMPL_NAME,hash_pair_list)

#undef TEMPLATE_PREFIX
#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

#if defined(HASH_INCLUDE_MODULE)

module TEMPLATE(hash_m)

  use global_m
  use messages_m
  use profiling_m

#endif
#if defined(HASH_INCLUDE_PREFIX) || defined(HASH_INCLUDE_MODULE)

#define PAIR_INCLUDE_PREFIX
#include "tpair_inc.F90"
#undef PAIR_INCLUDE_PREFIX

#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX

#define DARR_INCLUDE_PREFIX
#include "tdarr_inc.F90"
#undef DARR_INCLUDE_PREFIX

#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

  use kinds_m

  use HASH_KEY_TYPE_MODULE_NAME

  use HASH_KEY_FUNCTION_MODULE_NAME

#endif
#if defined(HASH_INCLUDE_MODULE)

  use HASH_VAL_TYPE_MODULE_NAME

  implicit none

  private

  public ::                     &
    TEMPLATE(HASH_OK),          &
    TEMPLATE(HASH_KEY_ERROR),   &
    TEMPLATE(HASH_EMPTY_ERROR)

  public ::           &
    TEMPLATE(hash_t)

  public ::                    &
    TEMPLATE(hash_iterator_t)

  public ::                 &
    TEMPLATE(hash_len),     &
    TEMPLATE(hash_has_key), &
    TEMPLATE(hash_init),    &
    TEMPLATE(hash_next),    &
    TEMPLATE(hash_pop),     &
    TEMPLATE(hash_set),     &
    TEMPLATE(hash_get),     &
    TEMPLATE(hash_del),     &
    TEMPLATE(hash_extend),  &
    TEMPLATE(hash_copy),    &
    TEMPLATE(hash_end)
    
#endif
#if defined(HASH_INCLUDE_HEADER) || defined(HASH_INCLUDE_MODULE)

#define PAIR_INCLUDE_HEADER
#include "tpair_inc.F90"
#undef PAIR_INCLUDE_HEADER

#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER

#define DARR_INCLUDE_HEADER
#include "tdarr_inc.F90"
#undef DARR_INCLUDE_HEADER

#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

  real(kind=wp), parameter :: INTERNAL(HASH_FACTOR) = DECORATE(HASH_GROWTH_FACTOR,wp)
  integer,       parameter :: INTERNAL(HASH_SIZE)   = HASH_INITIAL_SIZE

  integer, parameter :: TEMPLATE(HASH_OK)          = 0
  integer, parameter :: TEMPLATE(HASH_KEY_ERROR)   =-1
  integer, parameter :: TEMPLATE(HASH_EMPTY_ERROR) =-2

  type :: TEMPLATE(hash_t)
    private
    integer                               :: size = 0
    integer                               :: used = 0
    type(TEMPLATE(hash_pair_list_darr_t)) :: data
  end type TEMPLATE(hash_t)

  type :: TEMPLATE(hash_iterator_t)
    private
    type(TEMPLATE(hash_t)),                pointer :: hash =>null()
    type(TEMPLATE(hash_pair_list_darr_iterator_t)) :: itrd
    type(TEMPLATE(hash_pair_list_iterator_t))      :: itrl
  end type TEMPLATE(hash_iterator_t)

  interface TEMPLATE(hash_init)
    module procedure INTERNAL(hash_init)
    module procedure INTERNAL(hash_iterator_init_hash)
    module procedure INTERNAL(hash_iterator_init_iterator)
  end interface TEMPLATE(hash_init)

  interface TEMPLATE(hash_next)
    module procedure INTERNAL(hash_iterator_next_item)
    module procedure INTERNAL(hash_iterator_next_key)
    module procedure INTERNAL(hash_iterator_next_val)
  end interface TEMPLATE(hash_next)

  interface TEMPLATE(hash_pop)
    module procedure INTERNAL(hash_pop_item)
    module procedure INTERNAL(hash_pop_key)
    module procedure INTERNAL(hash_pop_val)
  end interface TEMPLATE(hash_pop)

  interface TEMPLATE(hash_del)
    module procedure INTERNAL(hash_del_item)
    module procedure INTERNAL(hash_del_key)
    module procedure INTERNAL(hash_del_val)
  end interface TEMPLATE(hash_del)

  interface TEMPLATE(hash_copy)
    module procedure INTERNAL(hash_copy)
    module procedure INTERNAL(hash_iterator_copy)
  end interface TEMPLATE(hash_copy)

  interface TEMPLATE(hash_end)
    module procedure INTERNAL(hash_end)
    module procedure INTERNAL(hash_iterator_end)
  end interface TEMPLATE(hash_end)

#endif
#if defined(HASH_INCLUDE_MODULE)

contains

#endif
#if defined(HASH_INCLUDE_BODY) || defined(HASH_INCLUDE_MODULE)

#define PAIR_INCLUDE_BODY
#include "tpair_inc.F90"
#undef PAIR_INCLUDE_BODY

#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY

#define DARR_INCLUDE_BODY
#include "tdarr_inc.F90"
#undef DARR_INCLUDE_BODY

#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

  ! ---------------------------------------------------------
  elemental function TEMPLATE(hash_len)(this) result(len)
    type(TEMPLATE(hash_t)), intent(in) :: this

    integer :: len

    len = this%used

  end function TEMPLATE(hash_len)

  ! ---------------------------------------------------------
  elemental function INTERNAL(hash_calc_size)(this) result(that)
    integer, intent(in) :: this

    integer :: that

    that = int(ceiling(INTERNAL(HASH_FACTOR)*real(this,kind=wp)))

  end function INTERNAL(hash_calc_size)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_rehash)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this

    type(HASH_KEY_TYPE_NAME),         pointer :: pkey
    type(TEMPLATE(hash_pair_list_t))          :: list
    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    type(TEMPLATE(hash_pair_t)),      pointer :: pair

    PUSH_SUB(INTERNAL(hash_rehash))

    if(this%size<INTERNAL(hash_calc_size)(this%used))then
      call TEMPLATE(hash_pair_list_init)(list)
      do
        nullify(pair)
        call INTERNAL(hash_pop_pair)(this, pair)
        if(.not.associated(pair)) exit
        call TEMPLATE(hash_pair_list_push)(list, pair)
      end do
      this%size = INTERNAL(hash_calc_size)(this%size)
      call TEMPLATE(hash_pair_list_darr_realloc)(this%data, this%size)
      do
        nullify(pkey, plst, pair)
        call TEMPLATE(hash_pair_list_pop)(list, pair)
        if(.not.associated(pair)) exit
        call TEMPLATE(hash_pair_get)(pair, pkey)
        ASSERT(associated(pkey))
        call INTERNAL(hash_get_list)(this, pkey, plst)
        call TEMPLATE(hash_pair_list_push)(plst, pair)
        this%used = this%used + 1
      end do
      call TEMPLATE(hash_pair_list_end)(list)
      nullify(pkey, plst, pair)
    end if

    POP_SUB(INTERNAL(hash_rehash))
  end subroutine INTERNAL(hash_rehash)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_init)(this, size)
    type(TEMPLATE(hash_t)), intent(out) :: this
    integer,      optional, intent(in)  :: size

    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    integer                                   :: indx

    PUSH_SUB(INTERNAL(hash_init))

    ASSERT(INTERNAL(HASH_SIZE)>1)
    ASSERT(INTERNAL(HASH_FACTOR)>1.0_wp)
    this%size=INTERNAL(HASH_SIZE)
    if(present(size))this%size=max(size,this%size)
    call TEMPLATE(hash_pair_list_darr_init)(this%data, this%size)
    do indx=1, this%size
      nullify(plst)
      SAFE_ALLOCATE(plst)
      call TEMPLATE(hash_pair_list_init)(plst)
      call TEMPLATE(hash_pair_list_darr_set)(this%data, indx, plst)
    end do

    POP_SUB(INTERNAL(hash_init))
  end subroutine INTERNAL(hash_init)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_pair)(this, pair)
    type(TEMPLATE(hash_t)),       intent(inout) :: this
    type(TEMPLATE(hash_pair_t)), pointer        :: pair

    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    integer                                   :: indx

    PUSH_SUB(INTERNAL(hash_pop_pair))

    nullify(pair)
    if(this%used>0)then
      do indx = 1, this%size
        nullify(pair, plst)
        call TEMPLATE(hash_pair_list_darr_get)(this%data, indx, plst)
        call TEMPLATE(hash_pair_list_pop)(plst, pair)
        if(associated(pair))then
          this%used = this%used - 1
          exit
        end if
      end do
      ASSERT(associated(pair))
      ASSERT(this%used>=0)
    end if

    POP_SUB(INTERNAL(hash_pop_pair))
  end subroutine INTERNAL(hash_pop_pair)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_item)(this, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr

    type(TEMPLATE(hash_pair_t)), pointer :: pair
    integer                              :: jerr

    PUSH_SUB(INTERNAL(hash_pop_item))

    jerr = TEMPLATE(HASH_EMPTY_ERROR)
    nullify(key, val, pair)
    call INTERNAL(hash_pop_pair)(this, pair)
    if(associated(pair))then
      call TEMPLATE(hash_pair_get)(pair, key, val)
      call TEMPLATE(hash_pair_del)(pair)
      jerr = TEMPLATE(HASH_OK)
    end if
    nullify(pair)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(hash_pop_item))
  end subroutine INTERNAL(hash_pop_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_key)(this, key, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    integer,         optional, intent(out)   :: ierr

    type(HASH_VAL_TYPE_NAME), pointer :: val

    PUSH_SUB(INTERNAL(hash_pop_key))

    nullify(key, val)
    call INTERNAL(hash_pop_item)(this, key, val, ierr)
    nullify(val)

    POP_SUB(INTERNAL(hash_pop_key))
  end subroutine INTERNAL(hash_pop_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_val)(this, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr

    type(HASH_KEY_TYPE_NAME), pointer :: key

    PUSH_SUB(INTERNAL(hash_pop_val))

    nullify(key, val)
    call INTERNAL(hash_pop_item)(this, key, val, ierr)
    nullify(key)

    POP_SUB(INTERNAL(hash_pop_val))
  end subroutine INTERNAL(hash_pop_val)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_get_list)(this, key, list)
    type(TEMPLATE(hash_t)),   intent(in) :: this
    type(HASH_KEY_TYPE_NAME), intent(in) :: key
    type(TEMPLATE(hash_pair_list_t)),  pointer     :: list

    integer :: ipos

    PUSH_SUB(INTERNAL(hash_get_list))

    nullify(list)
    ipos = modulo(HASH_KEY_FUNCTION_NAME(key), this%size) + 1
    call TEMPLATE(hash_pair_list_darr_get)(this%data, ipos, list)
    ASSERT(associated(list))

    POP_SUB(INTERNAL(hash_get_list))
  end subroutine INTERNAL(hash_get_list)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_get_list_pair)(this, key, list, pair)
    type(TEMPLATE(hash_t)),            intent(in) :: this
    type(HASH_KEY_TYPE_NAME),          intent(in) :: key
    type(TEMPLATE(hash_pair_list_t)), pointer     :: list
    type(TEMPLATE(hash_pair_t)),      pointer     :: pair

    type(HASH_KEY_TYPE_NAME),         pointer :: pkey
    type(TEMPLATE(hash_pair_list_iterator_t)) :: iter

    PUSH_SUB(INTERNAL(hash_get_list_pair))

    call INTERNAL(hash_get_list)(this, key, list)
    ASSERT(associated(list))
    call TEMPLATE(hash_pair_list_init)(iter, list)
    do
      nullify(pair, pkey)
      call TEMPLATE(hash_pair_list_next)(iter, pair)
      if(.not.associated(pair)) exit
      call TEMPLATE(hash_pair_get)(pair, pkey)
      ASSERT(associated(pkey))
      if(key==pkey)exit
    end do
    call TEMPLATE(hash_pair_list_end)(iter)
    nullify(pkey)

    POP_SUB(INTERNAL(hash_get_list_pair))
  end subroutine INTERNAL(hash_get_list_pair)

  ! ---------------------------------------------------------
  function TEMPLATE(hash_has_key)(this, key) result(has)
    type(TEMPLATE(hash_t)),   intent(in) :: this
    type(HASH_KEY_TYPE_NAME), intent(in) :: key

    logical :: has

    type(TEMPLATE(hash_pair_list_t)), pointer :: list
    type(TEMPLATE(hash_pair_t)),      pointer :: pair

    PUSH_SUB(TEMPLATE(hash_has_key))

    nullify(list, pair)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    has = associated(pair)
    nullify(list, pair)

    POP_SUB(TEMPLATE(hash_has_key))
  end function TEMPLATE(hash_has_key)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_set)(this, key, val)
    type(TEMPLATE(hash_t)),   intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), intent(in)    :: key
    type(HASH_VAL_TYPE_NAME), intent(in)    :: val

    type(TEMPLATE(hash_pair_list_t)), pointer :: list
    type(TEMPLATE(hash_pair_t)),      pointer :: pair

    PUSH_SUB(TEMPLATE(hash_set))

    nullify(list, pair)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(pair))then
      call TEMPLATE(hash_pair_set)(pair, val)
    else
      call TEMPLATE(hash_pair_new)(pair, key, val)
      call TEMPLATE(hash_pair_list_push)(list, pair)
      this%used = this%used + 1
      call INTERNAL(hash_rehash)(this)
    end if
    nullify(list, pair)

    POP_SUB(TEMPLATE(hash_set))
  end subroutine TEMPLATE(hash_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_get)(this, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(in)  :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)  :: key
    type(HASH_VAL_TYPE_NAME), pointer      :: val
    integer,         optional, intent(out) :: ierr

    type(TEMPLATE(hash_pair_list_t)), pointer :: list
    type(TEMPLATE(hash_pair_t)),      pointer :: pair
    integer                                   :: jerr

    PUSH_SUB(TEMPLATE(hash_get))

    jerr = TEMPLATE(HASH_KEY_ERROR)
    nullify(val, list, pair)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(pair))then
      call TEMPLATE(hash_pair_get)(pair, val)
      jerr = TEMPLATE(HASH_OK)
    end if
    nullify(list, pair)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(hash_get))
  end subroutine TEMPLATE(hash_get)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_del_pair)(this, key, pair, ierr)
    type(TEMPLATE(hash_t)),       intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),     intent(in)    :: key
    type(TEMPLATE(hash_pair_t)), pointer        :: pair
    integer,            optional, intent(out)   :: ierr

    type(TEMPLATE(hash_pair_list_t)), pointer :: list
    integer                                   :: jerr

    PUSH_SUB(INTERNAL(hash_del_pair))

    nullify(pair, list)
    jerr = TEMPLATE(HASH_KEY_ERROR)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(list).and.associated(pair))then
      call TEMPLATE(hash_pair_list_del)(list, pair, jerr)
      ASSERT(jerr==TEMPLATE(HASH_PAIR_LIST_OK))
      this%used = this%used - 1
      jerr = TEMPLATE(HASH_OK)
    end if
    nullify(list)
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(hash_del_pair))
  end subroutine INTERNAL(hash_del_pair)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_del_item)(this, that, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)    :: that
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr

    type(TEMPLATE(hash_pair_t)), pointer :: pair

    PUSH_SUB(INTERNAL(hash_del_item))

    nullify(key, val, pair)
    call INTERNAL(hash_del_pair)(this, that, pair, ierr)
    if(ierr==TEMPLATE(HASH_OK))then
      call TEMPLATE(hash_pair_get)(pair, key, val)
      call TEMPLATE(hash_pair_del)(pair)
    end if
    nullify(pair)

    POP_SUB(INTERNAL(hash_del_item))
  end subroutine INTERNAL(hash_del_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_del_key)(this, that, key, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)    :: that
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    integer,         optional, intent(out)   :: ierr

    type(HASH_VAL_TYPE_NAME), pointer :: val

    PUSH_SUB(INTERNAL(hash_del_key))

    nullify(key, val)
    call INTERNAL(hash_del_item)(this, that, key, val, ierr)
    nullify(val)

    POP_SUB(INTERNAL(hash_del_key))
  end subroutine INTERNAL(hash_del_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_del_val)(this, that, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)    :: that
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr

    type(HASH_KEY_TYPE_NAME), pointer :: key

    PUSH_SUB(INTERNAL(hash_del_val))

    nullify(val, key)
    call INTERNAL(hash_del_item)(this, that, key, val, ierr)
    nullify(key)

    POP_SUB(INTERNAL(hash_del_val))
  end subroutine INTERNAL(hash_del_val)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_extend)(this, that)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    type(TEMPLATE(hash_t)), intent(in)    :: that

    type(TEMPLATE(hash_iterator_t))   :: iter
    type(HASH_KEY_TYPE_NAME), pointer :: key
    type(HASH_VAL_TYPE_NAME), pointer :: val
    integer                           :: ierr

    PUSH_SUB(TEMPLATE(hash_extend))

    call INTERNAL(hash_iterator_init_hash)(iter, that)
    do
      nullify(key, val)
      call INTERNAL(hash_iterator_next_item)(iter, key, val, ierr)
      if(ierr/=TEMPLATE(HASH_OK)) exit
      call TEMPLATE(hash_set)(this, key, val)
    end do
    nullify(key, val)
    call INTERNAL(hash_iterator_end)(iter)

    POP_SUB(TEMPLATE(hash_extend))
  end subroutine TEMPLATE(hash_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_copy)(this, that)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    type(TEMPLATE(hash_t)), intent(in)    :: that

    PUSH_SUB(INTERNAL(hash_copy))

    call INTERNAL(hash_purge)(this)
    call TEMPLATE(hash_extend)(this, that)

    POP_SUB(INTERNAL(hash_copy))
  end subroutine INTERNAL(hash_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_purge)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this

    type(HASH_KEY_TYPE_NAME), pointer :: pkey
    type(HASH_VAL_TYPE_NAME), pointer :: pval
    integer                           :: ierr

    PUSH_SUB(INTERNAL(hash_purge))

    do
      nullify(pkey, pval)
      call INTERNAL(hash_pop_item)(this, pkey, pval, ierr)
      if(ierr/=TEMPLATE(HASH_OK))exit
    end do
    nullify(pkey, pval)
    ASSERT(this%used==0)

    POP_SUB(INTERNAL(hash_purge))
  end subroutine INTERNAL(hash_purge)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_end)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this

    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    integer                                   :: indx

    PUSH_SUB(INTERNAL(hash_end))

    call INTERNAL(hash_purge)(this)
    do indx=1, this%size
      nullify(plst)
      call TEMPLATE(hash_pair_list_darr_get)(this%data, indx, plst)
      call TEMPLATE(hash_pair_list_end)(plst)
      SAFE_DEALLOCATE_P(plst)
    end do
    call TEMPLATE(hash_pair_list_darr_end)(this%data)
    this%size = 0

    POP_SUB(INTERNAL(hash_end))
  end subroutine INTERNAL(hash_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_init_hash)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_t)),  target, intent(in)  :: that

    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    integer                                   :: ierr

    PUSH_SUB(INTERNAL(hash_iterator_init_hash))

    nullify(plst)
    call INTERNAL(hash_iterator_end)(this)
    if(that%used>0)then
      this%hash => that
      call TEMPLATE(hash_pair_list_darr_init)(this%itrd, that%data)
      call TEMPLATE(hash_pair_list_darr_next)(this%itrd, plst, ierr)
      ASSERT(ierr==TEMPLATE(HASH_PAIR_LIST_DARR_OK))
      ASSERT(associated(plst))
      call TEMPLATE(hash_pair_list_init)(this%itrl, plst)
      nullify(plst)
    end if

    POP_SUB(INTERNAL(hash_iterator_init_hash))
  end subroutine INTERNAL(hash_iterator_init_hash)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_init_iterator)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_iterator_t)), intent(in)  :: that

    PUSH_SUB(INTERNAL(hash_iterator_init_iterator))

    call INTERNAL(hash_iterator_copy)(this, that)

    POP_SUB(INTERNAL(hash_iterator_init_iterator))
  end subroutine INTERNAL(hash_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_pair)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(TEMPLATE(hash_pair_t)),    pointer        :: that

    type(TEMPLATE(hash_pair_list_t)), pointer :: plst
    integer                                   :: ierr

    PUSH_SUB(INTERNAL(hash_iterator_next_pair))

    nullify(that, plst)
    if(associated(this%hash))then
      do
        nullify(that)
        call TEMPLATE(hash_pair_list_next)(this%itrl, that)
        if(associated(that))exit
        nullify(that)
        call TEMPLATE(hash_pair_list_end)(this%itrl)
        call TEMPLATE(hash_pair_list_darr_next)(this%itrd, plst, ierr)
        if(ierr/=TEMPLATE(HASH_PAIR_LIST_DARR_OK))exit
        ASSERT(associated(plst))
        call TEMPLATE(hash_pair_list_init)(this%itrl, plst)
        nullify(plst)
      end do
    end if

    POP_SUB(INTERNAL(hash_iterator_next_pair))
  end subroutine INTERNAL(hash_iterator_next_pair)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_item)(this, key, val, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),       pointer        :: key
    type(HASH_VAL_TYPE_NAME),       pointer        :: val
    integer,               optional, intent(out)   :: ierr

    type(TEMPLATE(hash_pair_t)), pointer :: pair
    integer                              :: jerr

    PUSH_SUB(INTERNAL(hash_iterator_next_item))
    nullify(key, val, pair)
    jerr = TEMPLATE(HASH_EMPTY_ERROR)
    call INTERNAL(hash_iterator_next_pair)(this, pair)
    if(associated(pair))then
      call TEMPLATE(hash_pair_get)(pair, key, val)
      jerr = TEMPLATE(HASH_OK)
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(INTERNAL(hash_iterator_next_item))
  end subroutine INTERNAL(hash_iterator_next_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_key)(this, that, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),       pointer        :: that
    integer,               optional, intent(out)   :: ierr

    type(HASH_VAL_TYPE_NAME), pointer :: pval

    PUSH_SUB(INTERNAL(hash_iterator_next_key))

    nullify(that, pval)
    call INTERNAL(hash_iterator_next_item)(this, that, pval, ierr)
    nullify(pval)

    POP_SUB(INTERNAL(hash_iterator_next_key))
  end subroutine INTERNAL(hash_iterator_next_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_val)(this, that, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_VAL_TYPE_NAME),       pointer        :: that
    integer,               optional, intent(out)   :: ierr

    type(HASH_KEY_TYPE_NAME), pointer :: pkey

    PUSH_SUB(INTERNAL(hash_iterator_next_val))

    nullify(that, pkey)
    call INTERNAL(hash_iterator_next_item)(this, pkey, that, ierr)
    nullify(pkey)

    POP_SUB(INTERNAL(hash_iterator_next_val))
  end subroutine INTERNAL(hash_iterator_next_val)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_copy)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_iterator_t)), intent(in)  :: that

    PUSH_SUB(INTERNAL(hash_iterator_copy))

    call INTERNAL(hash_iterator_end)(this)
    this%hash => that%hash
    call TEMPLATE(hash_pair_list_darr_copy)(this%itrd, that%itrd)
    call TEMPLATE(hash_pair_list_copy)(this%itrl, that%itrl)

    POP_SUB(INTERNAL(hash_iterator_copy))
  end subroutine INTERNAL(hash_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(hash_iterator_end)(this)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this

    nullify(this%hash)
    call TEMPLATE(hash_pair_list_darr_end)(this%itrd)
    call TEMPLATE(hash_pair_list_end)(this%itrl)

  end subroutine INTERNAL(hash_iterator_end)

#endif
#if defined(HASH_INCLUDE_MODULE)

end module TEMPLATE(hash_m)

#endif

#undef PAIR_TEMPLATE_NAME
#undef PAIR_FRST_TYPE_NAME
#undef PAIR_SCND_TEMPLATE_NAME

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME

#undef DARR_TEMPLATE_NAME

#undef IHASH_TMPL_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:
