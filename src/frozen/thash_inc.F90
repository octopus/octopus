#include "global.h"
#include "util.h"

!HASH: HASH_TEMPLATE_NAME
!
!HASH: HASH_INITIAL_SIZE
!HASH: HASH_GROWTH_FACTOR
!
!HASH: HASH_KEY_TEMPLATE_NAME
!HASH: HASH_KEY_TYPE_NAME
!HASH: HASH_KEY_TYPE_MODULE_NAME
!
!HASH: HASH_KEY_FUNCTION_NAME
!HASH: HASH_KEY_FUNCTION_MODULE_NAME
!
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

#define LIST_TEMPLATE_NAME DECORATE(IHASH_TMPL_NAME,pair)
#define LIST_TYPE_NAME DECORATE(IHASH_TMPL_NAME,i_pair_t)

#define DARR_TEMPLATE_NAME DECORATE(IHASH_TMPL_NAME,pair_list)

#undef LIST
#undef DARR
#define LIST DECORATE(IHASH_TMPL_NAME,pair_list)
#define DARR DECORATE(IHASH_TMPL_NAME,pair_list_darr)

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

#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX

#define DARR_INCLUDE_PREFIX
#include "tdarr_inc.F90"
#undef DARR_INCLUDE_PREFIX

#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

  use kinds_m, only: wp

  use HASH_KEY_TYPE_MODULE_NAME
  use HASH_KEY_FUNCTION_MODULE_NAME

#endif
#if defined(HASH_INCLUDE_MODULE)

  use HASH_VAL_TYPE_MODULE_NAME

  implicit none

  private
  public ::                &
    TEMPLATE(hash_t),      &
    TEMPLATE(hash_len),    &
    TEMPLATE(hash_init),   &
    TEMPLATE(hash_next),   &
    TEMPLATE(hash_pop),    &
    TEMPLATE(hash_set),    &
    TEMPLATE(hash_get),    &
    TEMPLATE(hash_del),    &
    TEMPLATE(hash_extend), &
    TEMPLATE(hash_copy),   &
    TEMPLATE(hash_end)
    
  public ::                     &
    TEMPLATE(HASH_OK),          &
    TEMPLATE(HASH_KEY_ERROR),   &
    TEMPLATE(HASH_EMPTY_ERROR)

  public ::                    &
    TEMPLATE(hash_iterator_t)

#endif
#if defined(HASH_INCLUDE_HEADER) || defined(HASH_INCLUDE_MODULE)

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

  interface INTERNAL(pair_set)
    module procedure INTERNAL(pair_set_key)
    module procedure INTERNAL(pair_set_val)
  end interface INTERNAL(pair_set)

  interface INTERNAL(pair_get)
    module procedure INTERNAL(pair_get_key)
    module procedure INTERNAL(pair_get_val)
  end interface INTERNAL(pair_get)

  type :: INTERNAL(pair_t)
    private
    type(HASH_KEY_TYPE_NAME), pointer :: key =>null()
    type(HASH_VAL_TYPE_NAME), pointer :: val =>null()
  end type INTERNAL(pair_t)

  type :: TEMPLATE(hash_t)
    private
    integer                :: size = 0
    integer                :: used = 0
    type(EXTERNAL(DARR,t)) :: data
  end type TEMPLATE(hash_t)

  type :: TEMPLATE(hash_iterator_t)
    private
    type(TEMPLATE(hash_t)), pointer :: hash =>null()
    type(EXTERNAL(DARR,iterator_t)) :: itrd
    type(EXTERNAL(LIST,iterator_t)) :: itrl
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

#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY

#define DARR_INCLUDE_BODY
#include "tdarr_inc.F90"
#undef DARR_INCLUDE_BODY

#define TEMPLATE_PREFIX IHASH_TMPL_NAME
#include "template.h"

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_init)(this, key, val)
    type(INTERNAL(pair_t)),           intent(out) :: this
    type(HASH_KEY_TYPE_NAME), target, intent(in)  :: key
    type(HASH_VAL_TYPE_NAME), target, intent(in)  :: val
    !
    PUSH_SUB(INTERNAL(pair_init))
    this%key=>key
    this%val=>val
    POP_SUB(INTERNAL(pair_init))
    return
  end subroutine INTERNAL(pair_init)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_key)(this, that)
    type(INTERNAL(pair_t)),           intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), target, intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(pair_set_key))
    this%key=>that
    POP_SUB(INTERNAL(pair_set_key))
    return
  end subroutine INTERNAL(pair_set_key)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_set_val)(this, that)
    type(INTERNAL(pair_t)),           intent(inout) :: this
    type(HASH_VAL_TYPE_NAME), target, intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(pair_set_val))
    ASSERT(associated(this%key))
    this%val=>that
    POP_SUB(INTERNAL(pair_set_val))
    return
  end subroutine INTERNAL(pair_set_val)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_key)(this, that)
    type(INTERNAL(pair_t)),    intent(in) :: this
    type(HASH_KEY_TYPE_NAME), pointer     :: that
    !
    PUSH_SUB(INTERNAL(pair_get_key))
    nullify(that)
    if(associated(this%key))&
      that=>this%key
    POP_SUB(INTERNAL(pair_get_key))
    return
  end subroutine INTERNAL(pair_get_key)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_get_val)(this, that)
    type(INTERNAL(pair_t)),    intent(in) :: this
    type(HASH_VAL_TYPE_NAME), pointer     :: that
    !
    PUSH_SUB(INTERNAL(pair_get_val))
    ASSERT(associated(this%key))
    nullify(that)
    if(associated(this%val))&
      that=>this%val
    POP_SUB(INTERNAL(pair_get_val))
    return
  end subroutine INTERNAL(pair_get_val)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_copy)(this, that)
    type(INTERNAL(pair_t)),         intent(inout) :: this
    type(INTERNAL(pair_t)), target, intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(pair_copy))
    this%key=>that%key
    this%val=>that%val
    POP_SUB(INTERNAL(pair_copy))
    return
  end subroutine INTERNAL(pair_copy)

  ! -----------------------------------------------------
  subroutine INTERNAL(pair_end)(this)
    type(INTERNAL(pair_t)), intent(inout) :: this
    !
    PUSH_SUB(INTERNAL(pair_end))
    nullify(this%key, this%val)
    POP_SUB(INTERNAL(pair_end))
    return
  end subroutine INTERNAL(pair_end)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(hash_len)(this) result(len)
    type(TEMPLATE(hash_t)), intent(in) :: this
    !
    integer :: len
    !
    len=this%used
    return
  end function TEMPLATE(hash_len)

  ! ---------------------------------------------------------
  elemental function INTERNAL(hash_calc_size)(this) result(that)
    integer, intent(in) :: this
    !
    integer :: that
    !
    that=int(ceiling(INTERNAL(HASH_FACTOR)*real(this,kind=wp)))
    return
  end function INTERNAL(hash_calc_size)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_rehash)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    !
    type(HASH_KEY_TYPE_NAME), pointer :: pkey
    type(EXTERNAL(LIST,t))            :: list
    type(EXTERNAL(LIST,t)),   pointer :: plst
    type(INTERNAL(pair_t)),   pointer :: pair
    !
    PUSH_SUB(INTERNAL(hash_rehash))
    if(this%size<INTERNAL(hash_calc_size)(this%used))then
      call EXTERNAL(LIST,init)(list)
      do
        nullify(pair)
        call INTERNAL(hash_pop_pair)(this, pair)
        if(.not.associated(pair))exit
        call EXTERNAL(LIST,push)(list, pair)
      end do
      this%size=INTERNAL(hash_calc_size)(this%size)
      call EXTERNAL(DARR,realloc)(this%data, this%size)
      do
        nullify(pkey, plst, pair)
        call EXTERNAL(LIST,pop)(list, pair)
        if(.not.associated(pair))exit
        call INTERNAL(pair_get)(pair, pkey)
        ASSERT(associated(pkey))
        call INTERNAL(hash_get_list)(this, pkey, plst)
        call EXTERNAL(LIST,push)(plst, pair)
        this%used=this%used+1
      end do
      call EXTERNAL(LIST,end)(list)
      nullify(pkey, plst, pair)
    end if
    POP_SUB(INTERNAL(hash_rehash))
    return
  end subroutine INTERNAL(hash_rehash)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_init)(this, size)
    type(TEMPLATE(hash_t)), intent(out) :: this
    integer,      optional, intent(in)  :: size
    !
    type(EXTERNAL(LIST,t)), pointer :: plst
    integer                         :: indx
    !
    PUSH_SUB(INTERNAL(hash_init))
    ASSERT(INTERNAL(HASH_SIZE)>1)
    ASSERT(INTERNAL(HASH_FACTOR)>1.0_wp)
    this%size=INTERNAL(HASH_SIZE)
    if(present(size))this%size=max(size,this%size)
    call EXTERNAL(DARR,init)(this%data, this%size)
    do indx=1, this%size
      nullify(plst)
      SAFE_ALLOCATE(plst)
      call EXTERNAL(LIST,init)(plst)
      call EXTERNAL(DARR,set)(this%data, indx, plst)
    end do
    POP_SUB(INTERNAL(hash_init))
    return
  end subroutine INTERNAL(hash_init)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_pair)(this, pair)
    type(TEMPLATE(hash_t)),  intent(inout) :: this
    type(INTERNAL(pair_t)), pointer        :: pair
    !
    type(EXTERNAL(LIST,t)), pointer :: plst
    integer                         :: indx
    !
    PUSH_SUB(INTERNAL(hash_pop_pair))
    nullify(pair)
    if(this%used>0)then
      do indx = 1, this%size
        nullify(pair, plst)
        call EXTERNAL(DARR,get)(this%data, indx, plst)
        call EXTERNAL(LIST,pop)(plst, pair)
        if(associated(pair))then
          this%used=this%used-1
          exit
        end if
      end do
      ASSERT(associated(pair))
      ASSERT(this%used>=0)
    end if
    POP_SUB(INTERNAL(hash_pop_pair))
    return
  end subroutine INTERNAL(hash_pop_pair)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_item)(this, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr
    !
    type(INTERNAL(pair_t)), pointer :: pair
    integer                         :: jerr
    !
    PUSH_SUB(INTERNAL(hash_pop_item))
    jerr=TEMPLATE(HASH_EMPTY_ERROR)
    nullify(key, val, pair)
    call INTERNAL(hash_pop_pair)(this, pair)
    if(associated(pair))then
      call INTERNAL(pair_get)(pair, key)
      call INTERNAL(pair_get)(pair, val)
      call INTERNAL(pair_end)(pair)
      SAFE_DEALLOCATE_P(pair)
      jerr=TEMPLATE(HASH_OK)
    end if
    nullify(pair)
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(hash_pop_item))
    return
  end subroutine INTERNAL(hash_pop_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_key)(this, key, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), pointer        :: key
    integer,         optional, intent(out)   :: ierr
    !
    type(HASH_VAL_TYPE_NAME), pointer :: val
    !
    PUSH_SUB(INTERNAL(hash_pop_key))
    nullify(key, val)
    call INTERNAL(hash_pop_item)(this, key, val, ierr)
    nullify(val)
    POP_SUB(INTERNAL(hash_pop_key))
    return
  end subroutine INTERNAL(hash_pop_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_pop_val)(this, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr
    !
    type(HASH_KEY_TYPE_NAME), pointer :: key
    !
    PUSH_SUB(INTERNAL(hash_pop_val))
    nullify(key, val)
    call INTERNAL(hash_pop_item)(this, key, val, ierr)
    nullify(key)
    POP_SUB(INTERNAL(hash_pop_val))
    return
  end subroutine INTERNAL(hash_pop_val)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_get_list)(this, key, list)
    type(TEMPLATE(hash_t)),   intent(in) :: this
    type(HASH_KEY_TYPE_NAME), intent(in) :: key
    type(EXTERNAL(LIST,t)),  pointer     :: list
    !
    integer :: ipos
    !
    PUSH_SUB(INTERNAL(hash_get_list))
    nullify(list)
    ipos=modulo(HASH_KEY_FUNCTION_NAME(key), this%size)+1
    call EXTERNAL(DARR,get)(this%data, ipos, list)
    ASSERT(associated(list))
    POP_SUB(INTERNAL(hash_get_list))
    return
  end subroutine INTERNAL(hash_get_list)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_get_list_pair)(this, key, list, pair)
    type(TEMPLATE(hash_t)),   intent(in) :: this
    type(HASH_KEY_TYPE_NAME), intent(in) :: key
    type(EXTERNAL(LIST,t)),  pointer     :: list
    type(INTERNAL(pair_t)),  pointer     :: pair
    !
    type(HASH_KEY_TYPE_NAME), pointer :: pkey
    type(EXTERNAL(LIST,iterator_t))   :: iter
    !
    PUSH_SUB(INTERNAL(hash_get_list_pair))
    call INTERNAL(hash_get_list)(this, key, list)
    ASSERT(associated(list))
    call EXTERNAL(LIST,init)(iter, list)
    do
      nullify(pair, pkey)
      call EXTERNAL(LIST,next)(iter, pair)
      if(.not.associated(pair))exit
      call INTERNAL(pair_get)(pair, pkey)
      ASSERT(associated(pkey))
      if(key==pkey)exit
    end do
    call EXTERNAL(LIST,end)(iter)
    nullify(pkey)
    POP_SUB(INTERNAL(hash_get_list_pair))
    return
  end subroutine INTERNAL(hash_get_list_pair)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_set)(this, key, val)
    type(TEMPLATE(hash_t)),   intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), intent(in)    :: key
    type(HASH_VAL_TYPE_NAME), intent(in)    :: val
    !
    type(EXTERNAL(LIST,t)), pointer :: list
    type(INTERNAL(pair_t)), pointer :: pair
    !
    PUSH_SUB(TEMPLATE(hash_set))
    nullify(list, pair)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(pair))then
      call INTERNAL(pair_set)(pair, val)
    else
      SAFE_ALLOCATE(pair)
      call INTERNAL(pair_init)(pair, key, val)
      call EXTERNAL(LIST,push)(list, pair)
      this%used=this%used+1
      call INTERNAL(hash_rehash)(this)
    end if
    nullify(list, pair)
    POP_SUB(TEMPLATE(hash_set))
    return
  end subroutine TEMPLATE(hash_set)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_get)(this, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(in)  :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)  :: key
    type(HASH_VAL_TYPE_NAME), pointer      :: val
    integer,         optional, intent(out) :: ierr
    !
    type(EXTERNAL(LIST,t)), pointer :: list
    type(INTERNAL(pair_t)), pointer :: pair
    integer                         :: jerr
    !
    PUSH_SUB(TEMPLATE(hash_get))
    jerr=TEMPLATE(HASH_KEY_ERROR)
    nullify(val, list, pair)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(pair))then
      call INTERNAL(pair_get)(pair, val)
      jerr=TEMPLATE(HASH_OK)
    end if
    nullify(list, pair)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(hash_get))
    return
  end subroutine TEMPLATE(hash_get)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_del_pair)(this, key, pair)
    type(TEMPLATE(hash_t)),   intent(inout) :: this
    type(HASH_KEY_TYPE_NAME), intent(in)    :: key
    type(INTERNAL(pair_t)),  pointer        :: pair
    !
    type(EXTERNAL(LIST,t)), pointer :: list
    integer                         :: indx
    !
    PUSH_SUB(INTERNAL(hash_del_pair))
    nullify(pair, list)
    call INTERNAL(hash_get_list_pair)(this, key, list, pair)
    if(associated(pair))then
      indx=EXTERNAL(LIST,index)(list, pair)
      if(indx>0)then
        nullify(pair)
        call EXTERNAL(LIST,del)(list, indx, pair)
      end if
    end if
    nullify(list)
    POP_SUB(INTERNAL(hash_del_pair))
    return
  end subroutine INTERNAL(hash_del_pair)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_del)(this, key, val, ierr)
    type(TEMPLATE(hash_t)),    intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),  intent(in)    :: key
    type(HASH_VAL_TYPE_NAME), pointer        :: val
    integer,         optional, intent(out)   :: ierr
    !
    type(INTERNAL(pair_t)), pointer :: pair
    integer                         :: jerr
    !
    PUSH_SUB(TEMPLATE(hash_del))
    jerr=TEMPLATE(HASH_KEY_ERROR)
    nullify(val, pair)
    call INTERNAL(hash_del_pair)(this, key, pair)
    if(associated(pair))then
      call INTERNAL(pair_get)(pair, val)
      call INTERNAL(pair_end)(pair)
      SAFE_DEALLOCATE_P(pair)
      this%used=this%used-1
      jerr=TEMPLATE(HASH_OK)
    end if
    nullify(pair)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(hash_del))
    return
  end subroutine TEMPLATE(hash_del)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(hash_extend)(this, that)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    type(TEMPLATE(hash_t)), intent(in)    :: that
    !
    type(TEMPLATE(hash_iterator_t))   :: iter
    type(HASH_KEY_TYPE_NAME), pointer :: key
    type(HASH_VAL_TYPE_NAME), pointer :: val
    integer                           :: ierr
    !
    PUSH_SUB(TEMPLATE(hash_extend))
    call INTERNAL(hash_iterator_init_hash)(iter, that)
    do
      nullify(key, val)
      call INTERNAL(hash_iterator_next_item)(iter, key, val, ierr)
      if(ierr/=TEMPLATE(HASH_OK))exit
      call TEMPLATE(hash_set)(this, key, val)
    end do
    nullify(key, val)
    call INTERNAL(hash_iterator_end)(iter)
    POP_SUB(TEMPLATE(hash_extend))
    return
  end subroutine TEMPLATE(hash_extend)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_copy)(this, that)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    type(TEMPLATE(hash_t)), intent(in)    :: that
    !
    PUSH_SUB(INTERNAL(hash_copy))
    call INTERNAL(hash_purge)(this)
    call TEMPLATE(hash_extend)(this, that)
    POP_SUB(INTERNAL(hash_copy))
    return
  end subroutine INTERNAL(hash_copy)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_purge)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    !
    type(HASH_KEY_TYPE_NAME), pointer :: pkey
    type(HASH_VAL_TYPE_NAME), pointer :: pval
    integer                           :: ierr
    !
    PUSH_SUB(INTERNAL(hash_purge))
    do
      nullify(pkey, pval)
      call INTERNAL(hash_pop_item)(this, pkey, pval, ierr)
      if(ierr/=TEMPLATE(HASH_OK))exit
    end do
    nullify(pkey, pval)
    ASSERT(this%used==0)
    POP_SUB(INTERNAL(hash_purge))
    return
  end subroutine INTERNAL(hash_purge)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_end)(this)
    type(TEMPLATE(hash_t)), intent(inout) :: this
    !
    PUSH_SUB(INTERNAL(hash_end))
    call INTERNAL(hash_purge)(this)
    call EXTERNAL(DARR,end)(this%data)
    this%size=0
    POP_SUB(INTERNAL(hash_end))
    return
  end subroutine INTERNAL(hash_end)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_init_hash)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_t)),  target, intent(in)  :: that
    !
    type(EXTERNAL(LIST,t)), pointer :: plst
    integer                         :: ierr
    !
    PUSH_SUB(INTERNAL(hash_iterator_init_hash))
    nullify(plst)
    call INTERNAL(hash_iterator_end)(this)
    if(that%used>0)then
      this%hash=>that
      call EXTERNAL(DARR,init)(this%itrd, that%data)
      call EXTERNAL(DARR,next)(this%itrd, plst, ierr)
      ASSERT(ierr==EXTERNAL(DARR,OK))
      ASSERT(associated(plst))
      call EXTERNAL(LIST,init)(this%itrl, plst)
      nullify(plst)
    end if
    POP_SUB(INTERNAL(hash_iterator_init_hash))
    return
  end subroutine INTERNAL(hash_iterator_init_hash)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_init_iterator)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_iterator_t)), intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(hash_iterator_init_iterator))
    call INTERNAL(hash_iterator_copy)(this, that)
    POP_SUB(INTERNAL(hash_iterator_init_iterator))
    return
  end subroutine INTERNAL(hash_iterator_init_iterator)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_pair)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(INTERNAL(pair_t)),         pointer        :: that
    !
    type(EXTERNAL(LIST,t)), pointer :: plst
    integer                         :: ierr
    !
    PUSH_SUB(INTERNAL(hash_iterator_next_pair))
    nullify(that, plst)
    if(associated(this%hash))then
      do
        nullify(that)
        call EXTERNAL(LIST,next)(this%itrl, that)
        if(associated(that))exit
        nullify(that)
        call EXTERNAL(LIST,end)(this%itrl)
        call EXTERNAL(DARR,next)(this%itrd, plst, ierr)
        if(ierr/=EXTERNAL(DARR,OK))exit
        ASSERT(associated(plst))
        call EXTERNAL(LIST,init)(this%itrl, plst)
        nullify(plst)
      end do
    end if
    POP_SUB(INTERNAL(hash_iterator_next_pair))
    return
  end subroutine INTERNAL(hash_iterator_next_pair)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_item)(this, key, val, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),       pointer        :: key
    type(HASH_VAL_TYPE_NAME),       pointer        :: val
    integer,               optional, intent(out)   :: ierr
    !
    type(INTERNAL(pair_t)), pointer :: pair
    integer                         :: jerr
    !
    PUSH_SUB(INTERNAL(hash_iterator_next_item))
    nullify(key, val, pair)
    jerr=TEMPLATE(HASH_EMPTY_ERROR)
    call INTERNAL(hash_iterator_next_pair)(this, pair)
    if(associated(pair))then
      call INTERNAL(pair_get)(pair, key)
      call INTERNAL(pair_get)(pair, val)
      jerr=TEMPLATE(HASH_OK)
    end if
    if(present(ierr))ierr=jerr
    POP_SUB(INTERNAL(hash_iterator_next_item))
    return
  end subroutine INTERNAL(hash_iterator_next_item)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_key)(this, that, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_KEY_TYPE_NAME),       pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(HASH_VAL_TYPE_NAME), pointer :: pval
    !
    PUSH_SUB(INTERNAL(hash_iterator_next_key))
    nullify(that, pval)
    call INTERNAL(hash_iterator_next_item)(this, that, pval, ierr)
    nullify(pval)
    POP_SUB(INTERNAL(hash_iterator_next_key))
    return
  end subroutine INTERNAL(hash_iterator_next_key)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_next_val)(this, that, ierr)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    type(HASH_VAL_TYPE_NAME),       pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(HASH_KEY_TYPE_NAME), pointer :: pkey
    !
    PUSH_SUB(INTERNAL(hash_iterator_next_val))
    nullify(that, pkey)
    call INTERNAL(hash_iterator_next_item)(this, pkey, that, ierr)
    nullify(pkey)
    POP_SUB(INTERNAL(hash_iterator_next_val))
    return
  end subroutine INTERNAL(hash_iterator_next_val)

  ! ---------------------------------------------------------
  subroutine INTERNAL(hash_iterator_copy)(this, that)
    type(TEMPLATE(hash_iterator_t)), intent(out) :: this
    type(TEMPLATE(hash_iterator_t)), intent(in)  :: that
    !
    PUSH_SUB(INTERNAL(hash_iterator_copy))
    call INTERNAL(hash_iterator_end)(this)
    this%hash=>that%hash
    call EXTERNAL(DARR,copy)(this%itrd, that%itrd)
    call EXTERNAL(LIST,copy)(this%itrl, that%itrl)
    POP_SUB(INTERNAL(hash_iterator_copy))
    return
  end subroutine INTERNAL(hash_iterator_copy)

  ! ---------------------------------------------------------
  elemental subroutine INTERNAL(hash_iterator_end)(this)
    type(TEMPLATE(hash_iterator_t)), intent(inout) :: this
    !
    nullify(this%hash)
    call EXTERNAL(DARR,end)(this%itrd)
    call EXTERNAL(LIST,end)(this%itrl)
    return
  end subroutine INTERNAL(hash_iterator_end)

#endif
#if defined(HASH_INCLUDE_MODULE)

end module TEMPLATE(hash_m)

#endif

#undef LIST
#undef DARR

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef DARR_TEMPLATE_NAME

#undef IHASH_TMPL_NAME

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:
