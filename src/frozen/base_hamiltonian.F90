#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_KEY_TYPE_MODULE_NAME
#undef HASH_KEY_FUNCTION_NAME
#undef HASH_KEY_FUNCTION_MODULE_NAME
#undef HASH_VAL_TEMPLATE_NAME
#undef HASH_VAL_TYPE_NAME
#undef HASH_VAL_TYPE_MODULE_NAME
#undef HASH_INCLUDE_PREFIX
#undef HASH_INCLUDE_HEADER
#undef HASH_INCLUDE_BODY

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_PREFIX
#include "tdict.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

module base_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use strng_m, only: &
    operator(==)

  use strng_m, only:                   &
    string_t       => strng_t,         &
    string_init    => strng_init,      &
    string_hash    => strng_hash,      &
    string_tolower => strng_tolower,   &
    string_get     => strng_get,       &
    string_copy    => strng_copy,      &
    string_end     => strng_end

  use json_m, only: JSON_OK, JSON_TYPE_ERROR
  use json_m, only: json_object_t, json_object_iterator_t
  use json_m, only: json_init, json_get, json_next, json_end

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_set,       &
    config_dict_get

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t

  use base_term_m, only: &
    base_term__init__,   &
    base_term__update__, &
    base_term__add__,    &
    base_term__copy__,   &
    base_term__end__

  use base_term_m, only: &
    base_term_t

  use base_potential_m, only: &
    base_potential__init__,   &
    base_potential__start__,  &
    base_potential__update__, &
    base_potential__stop__,   &
    base_potential__add__,    &
    base_potential__copy__,   &
    base_potential__end__

  use base_potential_m, only: &
    base_potential_t

#define TEMPLATE_NAME base_hamiltonian
#define INCLUDE_PREFIX
#include "iterator_code.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_NAME

  implicit none

  private
  public ::                     &
    base_hamiltonian__init__,   &
    base_hamiltonian__start__,  &
    base_hamiltonian__update__, &
    base_hamiltonian__stop__,   &
    base_hamiltonian__add__,    &
    base_hamiltonian__copy__,   &
    base_hamiltonian__end__

  public ::                  &
    base_hamiltonian_new,    &
    base_hamiltonian_del,    &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_next,   &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_HEADER
#include "tdict.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

  integer, public, parameter :: HMLT_TYPE_NONE = 0 
  integer, public, parameter :: HMLT_TYPE_TERM = 1
  integer, public, parameter :: HMLT_TYPE_POTN = 2
  integer, public, parameter :: HMLT_TYPE_HMLT = 3
  
  type, private :: hterm_t
    private
    type(json_object_t),      pointer :: config =>null()
    type(base_system_t),      pointer :: sys    =>null()
    type(simulation_t),       pointer :: sim    =>null()
    type(base_term_t),        pointer :: term   =>null()
    type(base_potential_t),   pointer :: potn   =>null()
    type(base_hamiltonian_t), pointer :: hmlt   =>null()
    integer                           :: type   = HMLT_TYPE_NONE
  end type hterm_t

  type, public :: base_hamiltonian_t
    private
    type(json_object_t),      pointer :: config =>null()
    type(base_system_t),      pointer :: sys    =>null()
    type(simulation_t),       pointer :: sim    =>null()
    type(base_hamiltonian_t), pointer :: prnt   =>null()
    type(hterm_dict_t)                :: hdct
    type(config_dict_t)               :: dict
    type(base_hamiltonian_hash_t)     :: hash
    type(base_hamiltonian_list_t)     :: list
  end type base_hamiltonian_t

  interface hterm__init__
    module procedure hterm__init__hterm
    module procedure hterm__init__copy
  end interface hterm__init__

  interface hterm_get
    module procedure hterm_get_term
    module procedure hterm_get_potn
    module procedure hterm_get_hmlt
  end interface hterm_get

  interface base_hamiltonian__init__
    module procedure base_hamiltonian__init__hamiltonian
    module procedure base_hamiltonian__init__copy
  end interface base_hamiltonian__init__

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_hamiltonian
    module procedure base_hamiltonian_init_copy
  end interface base_hamiltonian_init

  interface base_hamiltonian_get
    module procedure base_hamiltonian_get_config
    module procedure base_hamiltonian_get_system
    module procedure base_hamiltonian_get_simulation
    module procedure base_hamiltonian_get_term
    module procedure base_hamiltonian_get_potn
    module procedure base_hamiltonian_get_hmlt
  end interface base_hamiltonian_get

  interface base_hamiltonian_copy
    module procedure base_hamiltonian_copy_hamiltonian
  end interface base_hamiltonian_copy

  interface base_hamiltonian_end
    module procedure base_hamiltonian_end_hamiltonian
  end interface base_hamiltonian_end

  integer, public, parameter :: BASE_HAMILTONIAN_OK          = BASE_HAMILTONIAN_HASH_OK
  integer, public, parameter :: BASE_HAMILTONIAN_KEY_ERROR   = BASE_HAMILTONIAN_HASH_KEY_ERROR
  integer, public, parameter :: BASE_HAMILTONIAN_EMPTY_ERROR = BASE_HAMILTONIAN_HASH_EMPTY_ERROR

#define TEMPLATE_NAME base_hamiltonian
#define INCLUDE_HEADER
#include "iterator_code.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_NAME

contains

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_BODY
#include "tdict.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine hterm__inew__(this)
    type(hterm_t), pointer :: this
    !
    PUSH_SUB(hterm__inew__)
    nullify(this)
    SAFE_ALLOCATE(this)
    call hterm__inull__(this)
    POP_SUB(hterm__inew__)
    return
  end subroutine hterm__inew__

  ! ---------------------------------------------------------
  subroutine hterm__idel__(this)
    type(hterm_t), pointer :: this
    !
    PUSH_SUB(hterm__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(hterm__idel__)
    return
  end subroutine hterm__idel__

  ! ---------------------------------------------------------
  subroutine hterm__inull__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__inull__)
    nullify(this%config, this%sys, this%sim, this%term, this%potn, this%hmlt)
    this%type = HMLT_TYPE_NONE
    POP_SUB(hterm__inull__)
    return
  end subroutine hterm__inull__

  ! ---------------------------------------------------------
  subroutine hterm__iinit__(this, sys, config)
    type(hterm_t),               intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: ierr
    !
    PUSH_SUB(hterm__iinit__)
    call hterm__inull__(this)
    this%config=>config
    this%sys=>sys
    call json_get(config, "type", this%type, ierr)
    ASSERT(ierr==JSON_OK)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      SAFE_ALLOCATE(this%term)
    case(HMLT_TYPE_POTN)
      SAFE_ALLOCATE(this%potn)
    case(HMLT_TYPE_HMLT)
      SAFE_ALLOCATE(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__iinit__)
    return
  end subroutine hterm__iinit__

  ! ---------------------------------------------------------
  recursive subroutine hterm__init__hterm(this, sys, config)
    type(hterm_t),       intent(out) :: this
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(hterm__init__hterm)
    call hterm__iinit__(this, sys, config)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__init__(this%term, sys, config)
    case(HMLT_TYPE_POTN)
      call base_potential__init__(this%potn, sys, config)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__init__(this%hmlt, sys, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__init__hterm)
    return
  end subroutine hterm__init__hterm

  ! ---------------------------------------------------------
  recursive subroutine hterm__init__copy(this, that)
    type(hterm_t), intent(out) :: this
    type(hterm_t), intent(in)  :: that
    !
    PUSH_SUB(hterm__init__copy)
    if(associated(that%config).and.associated(that%sys))&
      call hterm__init__(this, that%sys, that%config)
    POP_SUB(hterm__init__copy)
    return
  end subroutine hterm__init__copy

  ! ---------------------------------------------------------
  recursive subroutine hterm__start__(this, sim)
    type(hterm_t),                intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim
    !
    PUSH_SUB(hterm__start__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call base_potential__start__(this%potn, sim)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__start__(this%hmlt, sim)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__start__)
    return
  end subroutine hterm__start__

  ! ---------------------------------------------------------
  recursive subroutine hterm__update__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__update__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__update__(this%term)
    case(HMLT_TYPE_POTN)
      call base_potential__update__(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__update__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__update__)
    return
  end subroutine hterm__update__

  ! ---------------------------------------------------------
  recursive subroutine hterm__stop__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__stop__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call base_potential__stop__(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__stop__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__stop__)
    return
  end subroutine hterm__stop__

  ! ---------------------------------------------------------
  recursive subroutine hterm__add__(this, that, config)
    type(hterm_t),       intent(inout) :: this
    type(hterm_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(hterm__add__)
    ASSERT(this%type==that%type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__add__(this%term, that%term, config)
    case(HMLT_TYPE_POTN)
      call base_potential__add__(this%potn, that%potn, config)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__add__(this%hmlt, that%hmlt, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__add__)
    return
  end subroutine hterm__add__

  ! ---------------------------------------------------------
  subroutine hterm_get_term(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_term_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_term)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_TERM)
      ASSERT(associated(this%term))
      that=>this%term
    end if
    POP_SUB(hterm_get_term)
    return
  end subroutine hterm_get_term

  ! ---------------------------------------------------------
  subroutine hterm_get_potn(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_potential_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_potn)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_POTN)
      ASSERT(associated(this%potn))
      that=>this%potn
    end if
    POP_SUB(hterm_get_potn)
    return
  end subroutine hterm_get_potn

  ! ---------------------------------------------------------
  subroutine hterm_get_hmlt(this, that)
    type(hterm_t),             intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_hmlt)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_HMLT)
      ASSERT(associated(this%hmlt))
      that=>this%hmlt
    end if
    POP_SUB(hterm_get_hmlt)
    return
  end subroutine hterm_get_hmlt

  ! ---------------------------------------------------------
  subroutine hterm__icopy__(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm__icopy__)
    call hterm__iend__(this)
    if(associated(that%config).and.associated(that%sys))&
      call hterm__iinit__(this, that%sys, that%config)
    POP_SUB(hterm__icopy__)
    return
  end subroutine hterm__icopy__

  ! ---------------------------------------------------------
  recursive subroutine hterm__copy__(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm__copy__)
    call hterm__icopy__(this, that)
    select case(this%type)
    case(HMLT_TYPE_TERM)
      call base_term__copy__(this%term, that%term)
    case(HMLT_TYPE_POTN)
      call base_potential__copy__(this%potn, that%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__copy__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__copy__)
    return
  end subroutine hterm__copy__

  ! ---------------------------------------------------------
  subroutine hterm__iend__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__iend__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      SAFE_DEALLOCATE_P(this%term)
    case(HMLT_TYPE_POTN)
      SAFE_DEALLOCATE_P(this%potn)
    case(HMLT_TYPE_HMLT)
      SAFE_DEALLOCATE_P(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    call hterm__inull__(this)
    POP_SUB(hterm__iend__)
    return
  end subroutine hterm__iend__

  ! ---------------------------------------------------------
  recursive subroutine hterm__end__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__end__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__end__(this%term)
    case(HMLT_TYPE_POTN)
      call base_potential__end__(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__end__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    call hterm__iend__(this)
    POP_SUB(hterm__end__)
    return
  end subroutine hterm__end__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_new(this, that)
    type(base_hamiltonian_t),  target, intent(inout) :: this
    type(base_hamiltonian_t), pointer                :: that
    !
    PUSH_SUB(base_hamiltonian_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call base_hamiltonian_list_push(this%list, that)
    POP_SUB(base_hamiltonian_new)
    return
  end subroutine base_hamiltonian_new

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__idel__(this)
    type(base_hamiltonian_t), pointer :: this
    !
    PUSH_SUB(base_hamiltonian__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(base_hamiltonian__idel__)
    return
  end subroutine base_hamiltonian__idel__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del(this)
    type(base_hamiltonian_t), pointer :: this
    !
    PUSH_SUB(base_hamiltonian_del)
    if(associated(this))then
      if(associated(this%prnt))then
        call base_hamiltonian_list_del(this%prnt%list, this)
        call base_hamiltonian_end(this)
        call base_hamiltonian__idel__(this)
      end if
    end if
    POP_SUB(base_hamiltonian_del)
    return
  end subroutine base_hamiltonian_del

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__inull__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian__inull__)
    nullify(this%config, this%sys, this%sim, this%prnt)
    POP_SUB(base_hamiltonian__inull__)
    return
  end subroutine base_hamiltonian__inull__
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iinit__(this, sys, config)
    type(base_hamiltonian_t),    intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: type, ierr
    !
    PUSH_SUB(base_hamiltonian__iinit__)
    call base_hamiltonian__inull__(this)
    this%config=>config
    this%sys=>sys
    call json_get(this%config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(type==HMLT_TYPE_HMLT)
    call hterm_dict_init(this%hdct)
    call config_dict_init(this%dict)
    call base_hamiltonian_hash_init(this%hash)
    call base_hamiltonian_list_init(this%list)
    POP_SUB(base_hamiltonian__iinit__)
    return
  end subroutine base_hamiltonian__iinit__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__hamiltonian(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config
    !
    type(json_object_iterator_t)             :: iter
    type(json_object_t),             pointer :: cnfg
    type(hterm_t),                   pointer :: htrm
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(base_hamiltonian__init__hamiltonian)
    nullify(cnfg, htrm)
    call base_hamiltonian__iinit__(this, sys, config)
    call json_init(iter, config)
    do
      nullify(cnfg, htrm)
      call json_next(iter, name, cnfg, ierr)
      if(ierr==JSON_TYPE_ERROR)cycle
      if(ierr/=JSON_OK)exit
      call hterm__inew__(htrm)
      ASSERT(associated(htrm))
      call hterm__init__(htrm, sys, cnfg)
      call hterm_dict_set(this%hdct, trim(adjustl(name)), htrm)
    end do
    call json_end(iter)
    nullify(cnfg, htrm)
    POP_SUB(base_hamiltonian__init__hamiltonian)
    return
  end subroutine base_hamiltonian__init__hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that
    !
    type(hterm_dict_iterator_t)              :: iter
    type(hterm_t),                   pointer :: isub, osub
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(base_hamiltonian__init__copy)
    nullify(osub, isub)
    if(associated(that%config).and.associated(that%sys))then
      call base_hamiltonian__iinit__(this, that%sys, that%config)
      call hterm_dict_init(iter, that%hdct)
      do
        nullify(osub, isub)
        call hterm_dict_next(iter, name, isub, ierr)
        call hterm_dict_next(iter, name, isub, ierr)
        if(ierr/=HTERM_DICT_OK)exit
        call hterm__inew__(osub)
        call hterm__init__(osub, isub)
        call hterm_dict_set(this%hdct, trim(adjustl(name)), osub)
      end do
      call hterm_dict_end(iter)
      nullify(osub, isub)
    end if
    POP_SUB(base_hamiltonian__init__copy)
    return
  end subroutine base_hamiltonian__init__copy

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_hamiltonian(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config
    !
    PUSH_SUB(base_hamiltonian_init_hamiltonian)
    call base_hamiltonian__init__(this, sys, config)
    POP_SUB(base_hamiltonian_init_hamiltonian)
    return
  end subroutine base_hamiltonian_init_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_init_copy)
    nullify(cnfg, osub, isub)
    call base_hamiltonian__init__(this, that)
    call base_hamiltonian_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_hamiltonian_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_new(this, osub)
      call base_hamiltonian_init(osub, isub)
      call base_hamiltonian__add__(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_hamiltonian_init_copy)
    return
  end subroutine base_hamiltonian_init_copy

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__istart__(this, sim)
    type(base_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_hamiltonian__istart__)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    POP_SUB(base_hamiltonian__istart__)
    return
  end subroutine base_hamiltonian__istart__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__start__(this, sim)
    type(base_hamiltonian_t),     intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian__start__)
    nullify(htrm)
    if(present(sim))then
      call base_hamiltonian__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sim))
        call base_hamiltonian__istart__(this, this%prnt%sim)
      end if
    end if
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__start__(htrm, sim)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(base_hamiltonian__start__)
    return
  end subroutine base_hamiltonian__start__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_start(this, sim)
    type(base_hamiltonian_t), intent(inout) :: this
    type(simulation_t),       intent(in)    :: sim
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_start)
    nullify(subs)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_start(subs, sim)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call base_hamiltonian__start__(this, sim)
    POP_SUB(base_hamiltonian_start)
    return
  end subroutine base_hamiltonian_start

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__update__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian__update__)
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__update__(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(base_hamiltonian__update__)
    return
  end subroutine base_hamiltonian__update__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_update(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_update)
    nullify(subs)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_update(subs)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call base_hamiltonian__update__(this)
    POP_SUB(base_hamiltonian_update)
    return
  end subroutine base_hamiltonian_update

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__stop__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian__stop__)
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__stop__(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(base_hamiltonian__stop__)
    return
  end subroutine base_hamiltonian__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_stop(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_stop)
    nullify(subs)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_stop(subs)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call base_hamiltonian__stop__(this)
    POP_SUB(base_hamiltonian_stop)
    return
  end subroutine base_hamiltonian_stop

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__add__(this, that, config)
    type(base_hamiltonian_t), intent(inout) :: this
    type(json_object_t),      intent(in)    :: config
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)              :: iter
    type(hterm_t),                   pointer :: mhtr, shtr
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(base_hamiltonian__add__)
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(mhtr, shtr)
      call hterm_dict_next(iter, name, shtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_dict_get(this%hdct, name, mhtr, ierr)
      if(ierr/=HTERM_DICT_OK)cycle
      call hterm__add__(mhtr, shtr, config)
    end do
    call hterm_dict_end(iter)
    nullify(mhtr, shtr)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_hamiltonian_hash_set(this%hash, config, that)
    POP_SUB(base_hamiltonian__add__)
    return
  end subroutine base_hamiltonian__add__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__get__(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_hamiltonian__get__)
    nullify(config, that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_hamiltonian_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)nullify(that)
    end if
    POP_SUB(base_hamiltonian__get__)
    return
  end subroutine base_hamiltonian__get__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hterm(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(hterm_t),           pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(base_hamiltonian_get_hterm)
    nullify(that)
    call hterm_dict_get(this%hdct, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(base_hamiltonian_get_hterm)
    return
  end subroutine base_hamiltonian_get_hterm

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_term(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_term_t),       pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_term)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_get_term)
    return
  end subroutine base_hamiltonian_get_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_potn(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_potential_t),  pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_potn)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_get_potn)
    return
  end subroutine base_hamiltonian_get_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_hmlt)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    if(.not.associated(that))&
      call base_hamiltonian__get__(this, name, that)
    POP_SUB(base_hamiltonian_get_hmlt)
    return
  end subroutine base_hamiltonian_get_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_config(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(json_object_t),     pointer     :: that
    !
    PUSH_SUB(base_hamiltonian_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_hamiltonian_get_config)
    return
  end subroutine base_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_system(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(base_system_t),     pointer     :: that
    !
    PUSH_SUB(base_hamiltonian_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(base_hamiltonian_get_system)
    return
  end subroutine base_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_simulation(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that
    !
    PUSH_SUB(base_hamiltonian_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_hamiltonian_get_simulation)
    return
  end subroutine base_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__icopy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    PUSH_SUB(base_hamiltonian__icopy__)
    call base_hamiltonian__iend__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_hamiltonian__iinit__(this, that%sys, that%config)
      if(associated(that%sim))&
        call base_hamiltonian__istart__(this, that%sim)
    end if
    POP_SUB(base_hamiltonian__icopy__)
    return
  end subroutine base_hamiltonian__icopy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__copy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)              :: iter
    type(hterm_t),                   pointer :: isub, osub
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(base_hamiltonian__copy__)
    call base_hamiltonian__icopy__(this, that)
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(osub, isub)
      call hterm_dict_next(iter, name, isub, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__inew__(osub)
      call hterm__copy__(osub, isub)
      call hterm_dict_set(this%hdct, trim(adjustl(name)), osub)
    end do
    call hterm_dict_end(iter)
    nullify(osub, isub)
    POP_SUB(base_hamiltonian__copy__)
    return
  end subroutine base_hamiltonian__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_copy_hamiltonian(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr
    !
    PUSH_SUB(base_hamiltonian_copy_hamiltonian)
    nullify(cnfg, osub, isub)
    call base_hamiltonian_end(this)
    call base_hamiltonian__copy__(this, that)
    call base_hamiltonian_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_hamiltonian_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian_new(this, osub)
      call base_hamiltonian_copy(osub, isub)
      call base_hamiltonian__add__(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_hamiltonian_copy_hamiltonian)
    return
  end subroutine base_hamiltonian_copy_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iend__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian__iend__)
    call base_hamiltonian__inull__(this)
    call hterm_dict_end(this%hdct)
    call config_dict_end(this%dict)
    call base_hamiltonian_hash_end(this%hash)
    call base_hamiltonian_list_end(this%list)
    POP_SUB(base_hamiltonian__iend__)
    return
  end subroutine base_hamiltonian__iend__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__end__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_t), pointer :: htrm
    integer                :: ierr
    !
    PUSH_SUB(base_hamiltonian__end__)
    do
      nullify(htrm)
      call hterm_dict_pop(this%hdct, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__end__(htrm)
      call hterm__idel__(htrm)
    end do
    nullify(htrm)
    call base_hamiltonian__iend__(this)
    POP_SUB(base_hamiltonian__end__)
    return
  end subroutine base_hamiltonian__end__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_end_hamiltonian(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(base_hamiltonian_t), pointer :: subs
    !
    PUSH_SUB(base_hamiltonian_end_hamiltonian)
    do
      nullify(subs)
      call base_hamiltonian_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_hamiltonian_end(subs)
      call base_hamiltonian__idel__(subs)
    end do
    nullify(subs)
    call base_hamiltonian__end__(this)
    POP_SUB(base_hamiltonian_end_hamiltonian)
    return
  end subroutine base_hamiltonian_end_hamiltonian

#define TEMPLATE_NAME base_hamiltonian
#define INCLUDE_BODY
#include "iterator_code.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_NAME

end module base_hamiltonian_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
