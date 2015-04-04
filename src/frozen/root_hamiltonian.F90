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

#define LIST_TEMPLATE_NAME root_hamiltonian
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME root_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME root_hamiltonian
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

module root_hamiltonian_m

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
    config_dict_init,      &
    config_dict_set,       &
    config_dict_get,       &
    config_dict_end

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use storage_m, only:  &
    storage_t,          &
    storage_init,       &
    storage_start,      &
    storage_update,     &
    storage_stop,       &
    storage_reset,      &
    storage_accumulate, &
    storage_get,        &
    storage_copy,       &
    storage_end

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t,         &
    base_system_get

  use base_term_m, only: &
    base_term__init__,   &
    base_term__update__, &
    base_term__reset__,  &
    base_term__acc__,    &
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
    base_potential__reset__,  &
    base_potential__acc__,    &
    base_potential__add__,    &
    base_potential__copy__,   &
    base_potential__end__

  use base_potential_m, only: &
    base_potential_t

  use base_functional_m, only: &
    base_functional__init__,   &
    base_functional__start__,  &
    base_functional__update__, &
    base_functional__stop__,   &
    base_functional__reset__,  &
    base_functional__acc__,    &
    base_functional__add__,    &
    base_functional__copy__,   &
    base_functional__end__
  
  use base_functional_m, only: &
    base_functional_t

#define TEMPLATE_NAME root_hamiltonian
#define INCLUDE_PREFIX
!#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_NAME

  implicit none

  private
  public ::                   &
    root_hamiltonian__aget__

  public ::                    &
    root_hamiltonian__rpush__, &
    root_hamiltonian__rpop__,  &
    root_hamiltonian__rdel__,  &
    root_hamiltonian__rset__,  &
    root_hamiltonian__rget__

  public ::                     &
    root_hamiltonian__new__,    &
    root_hamiltonian__del__,    &
    root_hamiltonian__init__,   &
    root_hamiltonian__start__,  &
    root_hamiltonian__update__, &
    root_hamiltonian__stop__,   &
    root_hamiltonian__reset__,  &
    root_hamiltonian__acc__,    &
    root_hamiltonian__add__,    &
    root_hamiltonian__get__,    &
    root_hamiltonian__copy__,   &
    root_hamiltonian__end__

  public ::               &
    root_hamiltonian_t,   &
    root_hamiltonian_set, &
    root_hamiltonian_get

  public ::         &
    HMLT_TYPE_NONE, &
    HMLT_TYPE_TERM, &
    HMLT_TYPE_POTN, &
    HMLT_TYPE_FNCT, &
    HMLT_TYPE_HMLT

  public ::                       &
    ROOT_HAMILTONIAN_OK,          &
    ROOT_HAMILTONIAN_KEY_ERROR,   &
    ROOT_HAMILTONIAN_EMPTY_ERROR

#define LIST_TEMPLATE_NAME root_hamiltonian
#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME root_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME root_hamiltonian
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

  integer, parameter :: HMLT_TYPE_NONE = 0 
  integer, parameter :: HMLT_TYPE_TERM = 1
  integer, parameter :: HMLT_TYPE_POTN = 2
  integer, parameter :: HMLT_TYPE_FNCT = 3
  integer, parameter :: HMLT_TYPE_HMLT = 4

  integer, parameter :: ROOT_HAMILTONIAN_OK          = ROOT_HAMILTONIAN_HASH_OK
  integer, parameter :: ROOT_HAMILTONIAN_KEY_ERROR   = ROOT_HAMILTONIAN_HASH_KEY_ERROR
  integer, parameter :: ROOT_HAMILTONIAN_EMPTY_ERROR = ROOT_HAMILTONIAN_HASH_EMPTY_ERROR

  type, public :: root_hamiltonian_raii_t
    private
    type(root_hamiltonian_t), pointer :: prnt =>null()
    type(root_hamiltonian_list_t)     :: list
  end type root_hamiltonian_raii_t

  type, private :: hterm_t
    private
    type(base_term_t),        pointer :: term =>null()
    type(base_potential_t),   pointer :: potn =>null()
    type(base_functional_t),  pointer :: fnct =>null()
    type(root_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: type = HMLT_TYPE_NONE
  end type hterm_t

  type :: root_hamiltonian_t
    private
    type(json_object_t),      pointer :: config =>null()
    type(base_system_t),      pointer :: sys    =>null()
    type(simulation_t),       pointer :: sim    =>null()
    real(kind=wp)                     :: energy = 0.0_wp
    type(storage_t)                   :: data
    type(hterm_dict_t)                :: hdct
    type(config_dict_t)               :: dict
    type(root_hamiltonian_hash_t)     :: hash
    type(root_hamiltonian_raii_t)     :: raii
  end type root_hamiltonian_t

  interface hterm__init__
    module procedure hterm__init__hterm
    module procedure hterm__init__copy
  end interface hterm__init__

  interface hterm__get__
    module procedure hterm__get__term
    module procedure hterm__get__potn
    module procedure hterm__get__fnct
    module procedure hterm__get__hmlt
  end interface hterm__get__

  interface root_hamiltonian__init__
    module procedure root_hamiltonian__init__hamiltonian
    module procedure root_hamiltonian__init__copy
  end interface root_hamiltonian__init__

  interface root_hamiltonian__new__
    module procedure root_hamiltonian__new__term
    module procedure root_hamiltonian__new__potn
    module procedure root_hamiltonian__new__fnct
    module procedure root_hamiltonian__new__hmlt
  end interface root_hamiltonian__new__

  interface root_hamiltonian__get__
    module procedure root_hamiltonian__get__term
    module procedure root_hamiltonian__get__potn
    module procedure root_hamiltonian__get__fnct
    module procedure root_hamiltonian__get__hmlt
  end interface root_hamiltonian__get__

  interface root_hamiltonian__end__
    module procedure root_hamiltonian__end__hamiltonian
  end interface root_hamiltonian__end__

  interface root_hamiltonian_set
    module procedure root_hamiltonian_set_info
  end interface root_hamiltonian_set

  interface root_hamiltonian_get
    module procedure root_hamiltonian_get_hamiltonian_by_config
    module procedure root_hamiltonian_get_hamiltonian_by_name
    module procedure root_hamiltonian_get_info
    module procedure root_hamiltonian_get_config
    module procedure root_hamiltonian_get_system
    module procedure root_hamiltonian_get_simulation
    module procedure root_hamiltonian_get_hamiltonian_1d
    module procedure root_hamiltonian_get_hamiltonian_md
  end interface root_hamiltonian_get

#define TEMPLATE_NAME root_hamiltonian
#define INCLUDE_HEADER
!#include "intrpl_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_NAME

contains

#define LIST_TEMPLATE_NAME root_hamiltonian
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME root_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME root_hamiltonian
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
  subroutine root_hamiltonian__aget__(this, that)
    type(root_hamiltonian_t), target, intent(in) :: this
    type(config_dict_t),     pointer             :: that
    !
    PUSH_SUB(root_hamiltonian__aget__)
    that=>this%dict
    POP_SUB(root_hamiltonian__aget__)
    return
  end subroutine root_hamiltonian__aget__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__raii_init__(this, that)
    type(root_hamiltonian_raii_t),              intent(out) :: this
    type(root_hamiltonian_t), optional, target, intent(in)  :: that
    !
    PUSH_SUB(root_hamiltonian__raii_init__)
    nullify(this%prnt)
    if(present(that))this%prnt=>that
    call root_hamiltonian_list_init(this%list)
    POP_SUB(root_hamiltonian__raii_init__)
    return
  end subroutine root_hamiltonian__raii_init__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__rpush__(this, that)
    type(root_hamiltonian_t), intent(inout) :: this
    type(root_hamiltonian_t), intent(in)    :: that
    !
    PUSH_SUB(root_hamiltonian__rpush__)
    call root_hamiltonian_list_push(this%raii%list, that)
    POP_SUB(root_hamiltonian__rpush__)
    return
  end subroutine root_hamiltonian__rpush__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__rpop__(this, that)
    type(root_hamiltonian_t),  intent(inout) :: this
    type(root_hamiltonian_t), pointer        :: that
    !
    PUSH_SUB(root_hamiltonian__rpop__)
    call root_hamiltonian_list_pop(this%raii%list, that)
    POP_SUB(root_hamiltonian__rpop__)
    return
  end subroutine root_hamiltonian__rpop__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__rdel__(this, that)
    type(root_hamiltonian_t),  intent(inout) :: this
    type(root_hamiltonian_t), pointer        :: that
    !
    PUSH_SUB(root_hamiltonian__rdel__)
    call root_hamiltonian_list_del(this%raii%list, that)
    POP_SUB(root_hamiltonian__rdel__)
    return
  end subroutine root_hamiltonian__rdel__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__rset__(this, that)
    type(root_hamiltonian_t),         intent(inout) :: this
    type(root_hamiltonian_t), target, intent(in)    :: that
    !
    PUSH_SUB(root_hamiltonian__rset__)
    this%raii%prnt=>that
    POP_SUB(root_hamiltonian__rset__)
    return
  end subroutine root_hamiltonian__rset__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__rget__(this, that)
    type(root_hamiltonian_t),  intent(in) :: this
    type(root_hamiltonian_t), pointer     :: that
    !
    PUSH_SUB(root_hamiltonian__rget__)
    that=>this%raii%prnt
    POP_SUB(root_hamiltonian__rget__)
    return
  end subroutine root_hamiltonian__rget__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__raii_end__(this)
    type(root_hamiltonian_raii_t), intent(inout) :: this
    !
    type(root_hamiltonian_t), pointer :: subs
    !
    PUSH_SUB(root_hamiltonian__raii_end__)
    nullify(this%prnt)
    do
      nullify(subs)
      call root_hamiltonian_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call root_hamiltonian__raii_end__(subs%raii)
      !call root_hamiltonian__raii_del__(subs)
      SAFE_DEALLOCATE_P(subs)
      nullify(subs)
    end do
    nullify(subs)
    call root_hamiltonian_list_end(this%list)
    POP_SUB(root_hamiltonian__raii_end__)
    return
  end subroutine root_hamiltonian__raii_end__

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
    nullify(this%term, this%potn, this%fnct, this%hmlt)
    this%type = HMLT_TYPE_NONE
    POP_SUB(hterm__inull__)
    return
  end subroutine hterm__inull__

  ! ---------------------------------------------------------
  subroutine hterm__iinit__(this, type)
    type(hterm_t), intent(out) :: this
    integer,       intent(in)  :: type
    !
    PUSH_SUB(hterm__iinit__)
    call hterm__inull__(this)
    this%type=type
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      SAFE_ALLOCATE(this%term)
    case(HMLT_TYPE_POTN)
      SAFE_ALLOCATE(this%potn)
    case(HMLT_TYPE_FNCT)
      SAFE_ALLOCATE(this%fnct)
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
    type(hterm_t),               intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t),         intent(in)  :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    type(base_system_t),        pointer :: psys
    integer                             :: type, ierr
    !
    PUSH_SUB(hterm__init__hterm)
    nullify(psys)
    call json_get(config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    call hterm__iinit__(this, type)
    call json_get(config, "name", name, ierr)
    if(ierr==JSON_OK)then
      call base_system_get(sys, name, psys)
      ASSERT(associated(psys))
    else
      psys=>sys
    end if
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__init__(this%term, psys, config)
    case(HMLT_TYPE_POTN)
      call base_potential__init__(this%potn, psys, config)
    case(HMLT_TYPE_FNCT)
      call base_functional__init__(this%fnct, psys, config)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__init__(this%hmlt, psys, config)
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
    call hterm__iinit__(this, that%type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__init__(this%term, that%term)
    case(HMLT_TYPE_POTN)
      call base_potential__init__(this%potn, that%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__init__(this%fnct, that%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__init__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
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
    case(HMLT_TYPE_FNCT)
      call base_functional__start__(this%fnct, sim)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__start__(this%hmlt, sim)
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
    case(HMLT_TYPE_FNCT)
      call base_functional__update__(this%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__update__(this%hmlt)
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
    case(HMLT_TYPE_FNCT)
      call base_functional__stop__(this%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__stop__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__stop__)
    return
  end subroutine hterm__stop__

  ! ---------------------------------------------------------
  recursive subroutine hterm__reset__(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm__reset__)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__reset__(this%term)
    case(HMLT_TYPE_POTN)
      call base_potential__reset__(this%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__reset__(this%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__reset__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__reset__)
    return
  end subroutine hterm__reset__

  ! ---------------------------------------------------------
  recursive subroutine hterm__acc__(this, that)
    type(hterm_t),       intent(inout) :: this
    type(hterm_t),       intent(in)    :: that
    !
    PUSH_SUB(hterm__acc__)
    ASSERT(this%type==that%type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__acc__(this%term, that%term)
    case(HMLT_TYPE_POTN)
      call base_potential__acc__(this%potn, that%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__acc__(this%fnct, that%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__acc__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__acc__)
    return
  end subroutine hterm__acc__

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
    case(HMLT_TYPE_FNCT)
      call base_functional__add__(this%fnct, that%fnct, config)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__add__(this%hmlt, that%hmlt, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm__add__)
    return
  end subroutine hterm__add__

  ! ---------------------------------------------------------
  subroutine hterm__get__term(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_term_t), pointer     :: that
    !
    PUSH_SUB(hterm__get__term)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_TERM)
      ASSERT(associated(this%term))
      that=>this%term
    end if
    POP_SUB(hterm__get__term)
    return
  end subroutine hterm__get__term

  ! ---------------------------------------------------------
  subroutine hterm__get__potn(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_potential_t), pointer     :: that
    !
    PUSH_SUB(hterm__get__potn)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_POTN)
      ASSERT(associated(this%potn))
      that=>this%potn
    end if
    POP_SUB(hterm__get__potn)
    return
  end subroutine hterm__get__potn

  ! ---------------------------------------------------------
  subroutine hterm__get__fnct(this, that)
    type(hterm_t),            intent(in) :: this
    type(base_functional_t), pointer     :: that
    !
    PUSH_SUB(hterm__get__fnct)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_FNCT)
      ASSERT(associated(this%fnct))
      that=>this%fnct
    end if
    POP_SUB(hterm__get__fnct)
    return
  end subroutine hterm__get__fnct

  ! ---------------------------------------------------------
  subroutine hterm__get__hmlt(this, that)
    type(hterm_t),             intent(in) :: this
    type(root_hamiltonian_t), pointer     :: that
    !
    PUSH_SUB(hterm__get__hmlt)
    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_HMLT)
      ASSERT(associated(this%hmlt))
      that=>this%hmlt
    end if
    POP_SUB(hterm__get__hmlt)
    return
  end subroutine hterm__get__hmlt

  ! ---------------------------------------------------------
  subroutine hterm__icopy__(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm__icopy__)
    call hterm__iend__(this)
    call hterm__iinit__(this, that%type)
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
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__copy__(this%term, that%term)
    case(HMLT_TYPE_POTN)
      call base_potential__copy__(this%potn, that%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__copy__(this%fnct, that%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__copy__(this%hmlt, that%hmlt)
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
    case(HMLT_TYPE_FNCT)
      SAFE_DEALLOCATE_P(this%fnct)
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
    case(HMLT_TYPE_FNCT)
      call base_functional__end__(this%fnct)
    case(HMLT_TYPE_HMLT)
      call root_hamiltonian__end__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    call hterm__iend__(this)
    POP_SUB(hterm__end__)
    return
  end subroutine hterm__end__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__inew__(this, name, that)
    type(root_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(hterm_t),           pointer        :: that
    !
    PUSH_SUB(root_hamiltonian__inew__)
    call hterm__inew__(that)
    ASSERT(associated(that))
    call hterm_dict_set(this%hdct, trim(adjustl(name)), that)
    POP_SUB(root_hamiltonian__inew__)
    return
  end subroutine root_hamiltonian__inew__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__new__term(this, name, that)
    type(root_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_term_t),       pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__new__term)
    nullify(that, htrm)
    call root_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_TERM)
    call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__new__term)
    return
  end subroutine root_hamiltonian__new__term

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__new__potn(this, name, that)
    type(root_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_potential_t),  pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__new__potn)
    nullify(that, htrm)
    call root_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_POTN)
    call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__new__potn)
    return
  end subroutine root_hamiltonian__new__potn

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__new__fnct(this, name, that)
    type(root_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_functional_t), pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__new__fnct)
    nullify(that, htrm)
    call root_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_FNCT)
    call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__new__fnct)
    return
  end subroutine root_hamiltonian__new__fnct

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__new__hmlt(this, name, that)
    type(root_hamiltonian_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: name
    type(root_hamiltonian_t), pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__new__hmlt)
    nullify(that, htrm)
    call root_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_HMLT)
    call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__new__hmlt)
    return
  end subroutine root_hamiltonian__new__hmlt

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__del__(this, name)
    type(root_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__del__)
    nullify(htrm)
    call hterm_dict_del(this%hdct, trim(adjustl(name)), htrm)
    call hterm__end__(htrm)
    call hterm__idel__(htrm)
    nullify(htrm)
    POP_SUB(root_hamiltonian__del__)
    return
  end subroutine root_hamiltonian__del__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__inull__(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(root_hamiltonian__inull__)
    nullify(this%config, this%sys, this%sim)
    this%energy=0.0_wp
    POP_SUB(root_hamiltonian__inull__)
    return
  end subroutine root_hamiltonian__inull__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__iinit__(this, sys, config)
    type(root_hamiltonian_t),    intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: type, nspin, ierr
    logical :: alloc
    !
    PUSH_SUB(root_hamiltonian__iinit__)
    call root_hamiltonian__inull__(this)
    this%config=>config
    this%sys=>sys
    call json_get(this%config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(type==HMLT_TYPE_HMLT)
    call base_system_get(this%sys, nspin=nspin)
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK)alloc=.false.
    call storage_init(this%data, nspin, full=.false., allocate=alloc)
    call hterm_dict_init(this%hdct)
    call config_dict_init(this%dict)
    call root_hamiltonian_hash_init(this%hash)
    call root_hamiltonian__raii_init__(this%raii)
    POP_SUB(root_hamiltonian__iinit__)
    return
  end subroutine root_hamiltonian__iinit__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__init__hamiltonian(this, sys, config)
    type(root_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config
    !
    type(json_object_iterator_t)        :: iter
    type(json_object_t),        pointer :: cnfg
    type(hterm_t),              pointer :: htrm
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(root_hamiltonian__init__hamiltonian)
    nullify(cnfg, htrm)
    call root_hamiltonian__iinit__(this, sys, config)
    call json_init(iter, config)
    do
      nullify(cnfg, htrm)
      call json_next(iter, name, cnfg, ierr)
      if(ierr==JSON_TYPE_ERROR)cycle
      if(ierr/=JSON_OK)exit
      call root_hamiltonian__inew__(this, name, htrm)
      call hterm__init__(htrm, sys, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, htrm)
    POP_SUB(root_hamiltonian__init__hamiltonian)
    return
  end subroutine root_hamiltonian__init__hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__init__copy(this, that)
    type(root_hamiltonian_t), intent(out) :: this
    type(root_hamiltonian_t), intent(in)  :: that
    !
    type(hterm_dict_iterator_t)              :: iter
    type(hterm_t),                   pointer :: isub, osub
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(root_hamiltonian__init__copy)
    nullify(osub, isub)
    if(associated(that%config).and.associated(that%sys))then
      call root_hamiltonian__iinit__(this, that%sys, that%config)
      call hterm_dict_init(iter, that%hdct)
      do
        nullify(osub, isub)
        call hterm_dict_next(iter, name, isub, ierr)
        if(ierr/=HTERM_DICT_OK)exit
        call root_hamiltonian__inew__(this, name, osub)
        call hterm__init__(osub, isub)
      end do
      call hterm_dict_end(iter)
      nullify(osub, isub)
    end if
    POP_SUB(root_hamiltonian__init__copy)
    return
  end subroutine root_hamiltonian__init__copy

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__istart__(this, sim)
    type(root_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(root_hamiltonian__istart__)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call storage_start(this%data, sim)
    POP_SUB(root_hamiltonian__istart__)
    return
  end subroutine root_hamiltonian__istart__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__start__(this, sim)
    type(root_hamiltonian_t),     intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim
    !
    type(hterm_dict_iterator_t)       :: iter
    type(root_hamiltonian_t), pointer :: prnt
    type(hterm_t),            pointer :: htrm
    integer                           :: ierr
    !
    PUSH_SUB(root_hamiltonian__start__)
    nullify(prnt, htrm)
    if(present(sim))then
      call root_hamiltonian__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        call root_hamiltonian__rget__(this, prnt)
        ASSERT(associated(prnt))
        ASSERT(associated(prnt%sim))
        call root_hamiltonian__istart__(this, prnt%sim)
        nullify(prnt)
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
    POP_SUB(root_hamiltonian__start__)
    return
  end subroutine root_hamiltonian__start__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__update__(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(root_hamiltonian__update__)
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__update__(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    call storage_update(this%data)
    POP_SUB(root_hamiltonian__update__)
    return
  end subroutine root_hamiltonian__update__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__stop__(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(root_hamiltonian__stop__)
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__stop__(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    call storage_stop(this%data)
    POP_SUB(root_hamiltonian__stop__)
    return
  end subroutine root_hamiltonian__stop__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__reset__(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(root_hamiltonian__reset__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call hterm_dict_init(iter, this%hdct)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__reset__(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    this%energy=0.0_wp
    call storage_reset(this%data)
    POP_SUB(root_hamiltonian__reset__)
    return
  end subroutine root_hamiltonian__reset__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__acc__(this, that)
    type(root_hamiltonian_t), intent(inout) :: this
    type(root_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: mhtr, shtr
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(root_hamiltonian__acc__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(mhtr, shtr)
      call hterm_dict_next(iter, name, shtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_dict_get(this%hdct, name, mhtr, ierr)
      if(ierr/=HTERM_DICT_OK)cycle
      call hterm__acc__(mhtr, shtr)
    end do
    call hterm_dict_end(iter)
    nullify(mhtr, shtr)
    this%energy=this%energy+that%energy
    call storage_accumulate(this%data, that%data)
    POP_SUB(root_hamiltonian__acc__)
    return
  end subroutine root_hamiltonian__acc__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__add__(this, that, config)
    type(root_hamiltonian_t), intent(inout) :: this
    type(json_object_t),      intent(in)    :: config
    type(root_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)              :: iter
    type(hterm_t),                   pointer :: mhtr, shtr
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                                  :: ierr
    !
    PUSH_SUB(root_hamiltonian__add__)
    ASSERT(associated(this%config))
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
    call root_hamiltonian_hash_set(this%hash, config, that)
    POP_SUB(root_hamiltonian__add__)
    return
  end subroutine root_hamiltonian__add__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_hamiltonian_by_config(this, config, that)
    type(root_hamiltonian_t),  intent(in) :: this
    type(json_object_t),       intent(in) :: config
    type(root_hamiltonian_t), pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(root_hamiltonian_get_hamiltonian_by_config)
    nullify(that)
    ASSERT(associated(this%config))
    call root_hamiltonian_hash_get(this%hash, config, that, ierr)
    if(ierr/=ROOT_HAMILTONIAN_OK)nullify(that)
    POP_SUB(root_hamiltonian_get_hamiltonian_by_config)
    return
  end subroutine root_hamiltonian_get_hamiltonian_by_config

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_hamiltonian_by_name(this, name, that)
    type(root_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(root_hamiltonian_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(root_hamiltonian_get_hamiltonian_by_name)
    nullify(config, that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)&
      call root_hamiltonian_get(this, config, that)
    POP_SUB(root_hamiltonian_get_hamiltonian_by_name)
    return
  end subroutine root_hamiltonian_get_hamiltonian_by_name

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_set_info(this, energy)
    type(root_hamiltonian_t), intent(inout) :: this
    real(kind=wp),  optional, intent(in)    :: energy
    !
    PUSH_SUB(root_hamiltonian_set_info)
    if(present(energy))&
      this%energy=energy
    POP_SUB(root_hamiltonian_set_info)
    return
  end subroutine root_hamiltonian_set_info

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__iget__(this, name, that)
    type(root_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(hterm_t),           pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(root_hamiltonian__iget__)
    nullify(that)
    call hterm_dict_get(this%hdct, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(root_hamiltonian__iget__)
    return
  end subroutine root_hamiltonian__iget__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__get__term(this, name, that)
    type(root_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_term_t),       pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__get__term)
    nullify(that, htrm)
    call root_hamiltonian__iget__(this, name, htrm)
    if(associated(htrm))&
      call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__get__term)
    return
  end subroutine root_hamiltonian__get__term

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__get__potn(this, name, that)
    type(root_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_potential_t),  pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__get__potn)
    nullify(that, htrm)
    call root_hamiltonian__iget__(this, name, htrm)
    if(associated(htrm))&
      call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__get__potn)
    return
  end subroutine root_hamiltonian__get__potn

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__get__fnct(this, name, that)
    type(root_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_functional_t), pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__get__fnct)
    nullify(that, htrm)
    call root_hamiltonian__iget__(this, name, htrm)
    if(associated(htrm))&
      call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__get__fnct)
    return
  end subroutine root_hamiltonian__get__fnct

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__get__hmlt(this, name, that)
    type(root_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(root_hamiltonian_t), pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(root_hamiltonian__get__hmlt)
    nullify(that, htrm)
    call root_hamiltonian__iget__(this, name, htrm)
    if(associated(htrm))&
      call hterm__get__(htrm, that)
    nullify(htrm)
    POP_SUB(root_hamiltonian__get__hmlt)
    return
  end subroutine root_hamiltonian__get__hmlt

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_info(this, size, nspin, energy)
    type(root_hamiltonian_t), intent(in)  :: this
    integer,        optional, intent(out) :: size
    integer,        optional, intent(out) :: nspin
    real(kind=wp),  optional, intent(out) :: energy
    !
    PUSH_SUB(root_hamiltonian_get_info)
    if(present(energy))&
      energy=this%energy
    call storage_get(this%data, size=size, dim=nspin)
    POP_SUB(root_hamiltonian_get_info)
    return
  end subroutine root_hamiltonian_get_info

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_config(this, that)
    type(root_hamiltonian_t), intent(in) :: this
    type(json_object_t),     pointer     :: that
    !
    PUSH_SUB(root_hamiltonian_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(root_hamiltonian_get_config)
    return
  end subroutine root_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_system(this, that)
    type(root_hamiltonian_t), intent(in) :: this
    type(base_system_t),     pointer     :: that
    !
    PUSH_SUB(root_hamiltonian_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(root_hamiltonian_get_system)
    return
  end subroutine root_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_simulation(this, that)
    type(root_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that
    !
    PUSH_SUB(root_hamiltonian_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(root_hamiltonian_get_simulation)
    return
  end subroutine root_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_hamiltonian_1d(this, that)
    type(root_hamiltonian_t),     intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(root_hamiltonian_get_hamiltonian_1d)
    call storage_get(this%data, that)
    POP_SUB(root_hamiltonian_get_hamiltonian_1d)
    return
  end subroutine root_hamiltonian_get_hamiltonian_1d

  ! ---------------------------------------------------------
  subroutine root_hamiltonian_get_hamiltonian_md(this, that)
    type(root_hamiltonian_t),       intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(root_hamiltonian_get_hamiltonian_md)
    call storage_get(this%data, that)
    POP_SUB(root_hamiltonian_get_hamiltonian_md)
    return
  end subroutine root_hamiltonian_get_hamiltonian_md

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__icopy__(this, that)
    type(root_hamiltonian_t), intent(inout) :: this
    type(root_hamiltonian_t), intent(in)    :: that
    !
    PUSH_SUB(root_hamiltonian__icopy__)
    call root_hamiltonian__iend__(this)
    if(associated(that%config).and.associated(that%sys))then
      call root_hamiltonian__iinit__(this, that%sys, that%config)
      this%energy=that%energy
      if(associated(that%sim))then
        call root_hamiltonian__istart__(this, that%sim)
        call storage_copy(this%data, that%data)
      end if
    end if
    POP_SUB(root_hamiltonian__icopy__)
    return
  end subroutine root_hamiltonian__icopy__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__copy__(this, that)
    type(root_hamiltonian_t), intent(inout) :: this
    type(root_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: isub, osub
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(root_hamiltonian__copy__)
    call root_hamiltonian__icopy__(this, that)
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(osub, isub)
      call hterm_dict_next(iter, name, isub, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call root_hamiltonian__inew__(this, name, osub)
      call hterm__copy__(osub, isub)
    end do
    call hterm_dict_end(iter)
    nullify(osub, isub)
    POP_SUB(root_hamiltonian__copy__)
    return
  end subroutine root_hamiltonian__copy__

  ! ---------------------------------------------------------
  subroutine root_hamiltonian__iend__(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(root_hamiltonian__iend__)
    call root_hamiltonian__inull__(this)
    call storage_end(this%data)
    call hterm_dict_end(this%hdct)
    call config_dict_end(this%dict)
    call root_hamiltonian_hash_end(this%hash)
    call root_hamiltonian__raii_end__(this%raii)
    POP_SUB(root_hamiltonian__iend__)
    return
  end subroutine root_hamiltonian__iend__

  ! ---------------------------------------------------------
  recursive subroutine root_hamiltonian__end__hamiltonian(this)
    type(root_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_t), pointer :: htrm
    integer                :: ierr
    !
    PUSH_SUB(root_hamiltonian__end__hamiltonian)
    do
      nullify(htrm)
      call hterm_dict_pop(this%hdct, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__end__(htrm)
      call hterm__idel__(htrm)
    end do
    nullify(htrm)
    call root_hamiltonian__iend__(this)
    POP_SUB(root_hamiltonian__end__hamiltonian)
    return
  end subroutine root_hamiltonian__end__hamiltonian

#define TEMPLATE_NAME root_hamiltonian
#define INCLUDE_BODY
!#include "intrpl_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_NAME

end module root_hamiltonian_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
