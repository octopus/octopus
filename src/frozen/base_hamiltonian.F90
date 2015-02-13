#include "global.h"

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
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_init,      &
    config_dict_set,       &
    config_dict_get,       &
    config_dict_copy,      &
    config_dict_end

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t

  use base_term_m, only: &
    base_term__update__, &
    base_term__add__

  use base_term_m, only: &
    base_term_t,         &
    base_term_init,      &
    base_term_copy,      &
    base_term_end

  use base_potential_m, only: &
    base_potential__start__,  &
    base_potential__update__, &
    base_potential__stop__,   &
    base_potential__add__

  use base_potential_m, only: &
    base_potential_t,         &
    base_potential_init,      &
    base_potential_copy,      &
    base_potential_end

  implicit none

  private
  public ::                     &
    base_hamiltonian__start__,  &
    base_hamiltonian__update__, &
    base_hamiltonian__stop__,   &
    base_hamiltonian__add__,    &
    base_hamiltonian__get__

  public ::                  &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_next,   &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

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
    type(base_term_t),        pointer :: term =>null()
    type(base_potential_t),   pointer :: potn =>null()
    type(base_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: type = HMLT_TYPE_NONE
  end type hterm_t

  type, public :: base_hamiltonian_t
    private
    type(json_object_t),  pointer :: config => null()
    type(base_system_t),  pointer :: sys    => null()
    type(simulation_t),   pointer :: sim    => null()
    type(hterm_dict_t)            :: dict
    type(config_dict_t)           :: cdct
    type(base_hamiltonian_hash_t) :: hash
  end type base_hamiltonian_t

  type, public :: base_hamiltonian_iterator_t
    private
    type(base_hamiltonian_t),      pointer :: self =>null()
    type(hterm_dict_iterator_t)            :: ditr
    type(base_hamiltonian_hash_iterator_t) :: hitr
  end type base_hamiltonian_iterator_t

  interface hterm_get
    module procedure hterm_get_term
    module procedure hterm_get_potn
    module procedure hterm_get_hmlt
  end interface hterm_get

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_hamiltonian
    module procedure base_hamiltonian_iterator_init
  end interface base_hamiltonian_init

  interface base_hamiltonian_next
    module procedure base_hamiltonian_iterator_next_config_hamiltonian
    module procedure base_hamiltonian_iterator_next_config
    module procedure base_hamiltonian_iterator_next_hamiltonian
  end interface base_hamiltonian_next

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
    module procedure base_hamiltonian_iterator_copy
  end interface base_hamiltonian_copy

  interface base_hamiltonian_end
    module procedure base_hamiltonian_end_hamiltonian
    module procedure base_hamiltonian_iterator_end
  end interface base_hamiltonian_end

  integer, parameter :: HMLT_NAME_LEN = CONFIG_DICT_NAME_LEN

  integer, public, parameter :: BASE_HAMILTONIAN_OK          = BASE_HAMILTONIAN_HASH_OK
  integer, public, parameter :: BASE_HAMILTONIAN_KEY_ERROR   = BASE_HAMILTONIAN_HASH_KEY_ERROR
  integer, public, parameter :: BASE_HAMILTONIAN_EMPTY_ERROR = BASE_HAMILTONIAN_HASH_EMPTY_ERROR

contains

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
  recursive subroutine hterm_init(this, sys, config)
    type(hterm_t),       intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t), intent(in)  :: config
    !
    integer :: ierr
    !
    PUSH_SUB(hterm_init)
    ASSERT(this%type==HMLT_TYPE_NONE)
    call json_get(config, "type", this%type, ierr)
    ASSERT(ierr==JSON_OK)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      SAFE_ALLOCATE(this%term)
      call base_term_init(this%term, sys, config)
    case(HMLT_TYPE_POTN)
      SAFE_ALLOCATE(this%potn)
      call base_potential_init(this%potn, sys, config)
    case(HMLT_TYPE_HMLT)
      SAFE_ALLOCATE(this%hmlt)
      call base_hamiltonian_init(this%hmlt, sys, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm_init)
    return
  end subroutine hterm_init

  ! ---------------------------------------------------------
  recursive subroutine hterm__start__(this, sim)
    type(hterm_t),      intent(inout) :: this
    type(simulation_t), intent(in)    :: sim
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
    ASSERT(this%type==HMLT_TYPE_TERM)
    ASSERT(associated(this%term))
    that=>this%term
    POP_SUB(hterm_get_term)
    return
  end subroutine hterm_get_term

  ! ---------------------------------------------------------
  subroutine hterm_get_potn(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_potential_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_potn)
    ASSERT(this%type==HMLT_TYPE_POTN)
    ASSERT(associated(this%potn))
    that=>this%potn
    POP_SUB(hterm_get_potn)
    return
  end subroutine hterm_get_potn

  ! ---------------------------------------------------------
  subroutine hterm_get_hmlt(this, that)
    type(hterm_t),             intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_hmlt)
    ASSERT(this%type==HMLT_TYPE_HMLT)
    ASSERT(associated(this%hmlt))
    that=>this%hmlt
    POP_SUB(hterm_get_hmlt)
    return
  end subroutine hterm_get_hmlt

  ! ---------------------------------------------------------
  recursive subroutine hterm_copy(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm_copy)
    call hterm_end(this)
    this%type=that%type
    select case(this%type)
    case(HMLT_TYPE_TERM)
      SAFE_ALLOCATE(this%term)
      call base_term_copy(this%term, that%term)
    case(HMLT_TYPE_POTN)
      SAFE_ALLOCATE(this%potn)
      call base_potential_copy(this%potn, that%potn)
    case(HMLT_TYPE_HMLT)
      SAFE_ALLOCATE(this%hmlt)
      call base_hamiltonian_copy(this%hmlt, that%hmlt)
    end select
    POP_SUB(hterm_copy)
    return
  end subroutine hterm_copy

  ! ---------------------------------------------------------
  recursive subroutine hterm_end(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm_end)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term_end(this%term)
      SAFE_DEALLOCATE_P(this%term)
    case(HMLT_TYPE_POTN)
      call base_potential_end(this%potn)
      SAFE_DEALLOCATE_P(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian_end(this%hmlt)
      SAFE_DEALLOCATE_P(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    nullify(this%term,this%potn,this%hmlt)
    this%type=HMLT_TYPE_NONE
    POP_SUB(hterm_end)
    return
  end subroutine hterm_end

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_hamiltonian(this, sys, config)
    type(base_hamiltonian_t),    intent(out) :: this
    type(base_system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_iterator_t) :: iter
    type(json_object_t), pointer :: cnfg
    type(hterm_t),       pointer :: htrm
    character(len=HMLT_NAME_LEN) :: name
    integer                      :: type, ierr
    !
    PUSH_SUB(base_hamiltonian_init_hamiltonian)
    nullify(cnfg, htrm)
    this%config=>config
    this%sys=>sys
    nullify(this%sim)
    call json_get(this%config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(type==HMLT_TYPE_HMLT)
    call hterm_dict_init(this%dict)
    call json_init(iter, config)
    do
      nullify(cnfg, htrm)
      call json_next(iter, name, cnfg, ierr)
      if(ierr==JSON_TYPE_ERROR)cycle
      if(ierr/=JSON_OK)exit
      SAFE_ALLOCATE(htrm)
      call hterm_init(htrm, sys, cnfg)
      call hterm_dict_set(this%dict, trim(adjustl(name)), htrm)
    end do
    call json_end(iter)
    nullify(cnfg, htrm)
    call config_dict_init(this%cdct)
    call base_hamiltonian_hash_init(this%hash)
    POP_SUB(base_hamiltonian_init_hamiltonian)
    return
  end subroutine base_hamiltonian_init_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__start__(this, sim)
    type(base_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian__start__)
    nullify(htrm)
    this%sim=>sim
    call hterm_dict_init(iter, this%dict)
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
    call hterm_dict_init(iter, this%dict)
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
    call hterm_dict_init(iter, this%dict)
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
    type(hterm_dict_iterator_t)  :: iter
    type(hterm_t),       pointer :: mhtr, shtr
    character(len=HMLT_NAME_LEN) :: tnam, cnam
    integer                      :: ierr
    !
    PUSH_SUB(base_hamiltonian__add__)
    call hterm_dict_init(iter, that%dict)
    do
      nullify(mhtr, shtr)
      call hterm_dict_next(iter, tnam, shtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_dict_get(this%dict, tnam, mhtr, ierr)
      if(ierr/=HTERM_DICT_OK)cycle
      call hterm__add__(mhtr, shtr, config)
    end do
    call hterm_dict_end(iter)
    nullify(mhtr, shtr)
    call json_get(config, "name", cnam, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%cdct, trim(adjustl(cnam)), config)
    call base_hamiltonian_hash_set(this%hash, config, that)
    POP_SUB(base_hamiltonian__add__)
    return
  end subroutine base_hamiltonian__add__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__get__(this, name, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: name
    type(base_hamiltonian_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_hamiltonian__get__)
    nullify(that)
    call config_dict_get(this%cdct, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_hamiltonian_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)nullify(that)
    end if
    POP_SUB(base_hamiltonian__get__)
    return
  end subroutine base_hamiltonian__get__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hterm(this, name, that, ierr)
    type(base_hamiltonian_t), intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(hterm_t),           pointer      :: that
    integer,        optional, intent(out) :: ierr
    !
    integer :: jerr
    !
    PUSH_SUB(base_hamiltonian_get_hterm)
    nullify(that)
    call hterm_dict_get(this%dict, name, that, jerr)
    if(jerr/=HTERM_DICT_OK)nullify(that)
    if(present(ierr))ierr=jerr
    POP_SUB(base_hamiltonian_get_hterm)
    return
  end subroutine base_hamiltonian_get_hterm

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_term(this, name, that, ierr)
    type(base_hamiltonian_t), intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(base_term_t),       pointer      :: that
    integer,        optional, intent(out) :: ierr
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_term)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm, ierr)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_get_term)
    return
  end subroutine base_hamiltonian_get_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_potn(this, name, that, ierr)
    type(base_hamiltonian_t), intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(base_potential_t),  pointer      :: that
    integer,        optional, intent(out) :: ierr
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_potn)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm, ierr)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_get_potn)
    return
  end subroutine base_hamiltonian_get_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hmlt(this, name, that, ierr)
    type(base_hamiltonian_t),  intent(in)  :: this
    character(len=*),          intent(in)  :: name
    type(base_hamiltonian_t), pointer      :: that
    integer,         optional, intent(out) :: ierr
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_get_hmlt)
    nullify(that, htrm)
    call base_hamiltonian_get_hterm(this, name, htrm, ierr)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
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
  subroutine base_hamiltonian__copy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    PUSH_SUB(base_hamiltonian__copy__)
    this%config=>that%config
    this%sys=>that%sys
    this%sim=>that%sim
    call hterm_dict_init(this%dict, hterm_dict_len(that%dict))
    call config_dict_copy(this%cdct, that%cdct)
    call base_hamiltonian_hash_copy(this%hash, that%hash)
    POP_SUB(base_hamiltonian__copy__)
    return
  end subroutine base_hamiltonian__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_copy_hamiltonian(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)  :: iter
    type(hterm_t),       pointer :: ihtr, ohtr
    character(len=HMLT_NAME_LEN) :: name
    integer                      :: ierr
    !
    PUSH_SUB(base_hamiltonian_copy_hamiltonian)
    call base_hamiltonian__copy__(this, that)
    call hterm_dict_init(iter, that%dict)
    do
      nullify(ohtr, ihtr)
      call hterm_dict_next(iter, name, ihtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      SAFE_ALLOCATE(ohtr)
      call hterm_copy(ohtr, ihtr)
      call hterm_dict_set(this%dict, name, ohtr)
    end do
    call hterm_dict_end(iter)
    nullify(ohtr, ihtr)
    POP_SUB(base_hamiltonian_copy_hamiltonian)
    return
  end subroutine base_hamiltonian_copy_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__end__(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian__end__)
    nullify(this%config, this%sys, this%sim)
    call hterm_dict_end(this%dict)
    call config_dict_end(this%cdct)
    call base_hamiltonian_hash_end(this%hash)
    POP_SUB(base_hamiltonian__end__)
    return
  end subroutine base_hamiltonian__end__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_end_hamiltonian(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_t), pointer :: htrm
    integer                :: ierr
    !
    PUSH_SUB(base_hamiltonian_end_hamiltonian)
    do
      nullify(htrm)
      call hterm_dict_pop(this%dict, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
    end do
    nullify(htrm)
    call base_hamiltonian__end__(this)
    POP_SUB(base_hamiltonian_end_hamiltonian)
    return
  end subroutine base_hamiltonian_end_hamiltonian

 ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_init(this, that)
    type(base_hamiltonian_iterator_t), intent(out) :: this
    type(base_hamiltonian_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_hamiltonian_iterator_init)
    this%self=>that
    call hterm_dict_init(this%ditr, that%dict)
    call base_hamiltonian_hash_init(this%hitr, that%hash)
    POP_SUB(base_hamiltonian_iterator_init)
    return
  end subroutine base_hamiltonian_iterator_init

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_next_config_hamiltonian(this, config, that, ierr)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    type(json_object_t),              pointer        :: config
    type(base_hamiltonian_t),         pointer        :: that
    integer,                 optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_hamiltonian_iterator_next_config_hamiltonian)
    call base_hamiltonian_hash_next(this%hitr, config, that, ierr)
    POP_SUB(base_hamiltonian_iterator_next_config_hamiltonian)
    return
  end subroutine base_hamiltonian_iterator_next_config_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_next_config(this, that, ierr)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    type(json_object_t),              pointer        :: that
    integer,                 optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_hamiltonian_iterator_next_config)
    call base_hamiltonian_hash_next(this%hitr, that, ierr)
    POP_SUB(base_hamiltonian_iterator_next_config)
    return
  end subroutine base_hamiltonian_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_next_hamiltonian(this, that, ierr)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    type(base_hamiltonian_t),         pointer        :: that
    integer,                 optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_hamiltonian_iterator_next_hamiltonian)
    call base_hamiltonian_hash_next(this%hitr, that, ierr)
    POP_SUB(base_hamiltonian_iterator_next_hamiltonian)
    return
  end subroutine base_hamiltonian_iterator_next_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_copy(this, that)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    type(base_hamiltonian_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_hamiltonian_iterator_copy)
    this%self=>that%self
    call hterm_dict_copy(this%ditr, that%ditr)
    call base_hamiltonian_hash_copy(this%hitr, that%hitr)
    POP_SUB(base_hamiltonian_iterator_copy)
    return
  end subroutine base_hamiltonian_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_end(this)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian_iterator_end)
    nullify(this%self)
    call hterm_dict_end(this%ditr)
    call base_hamiltonian_hash_end(this%hitr)
    POP_SUB(base_hamiltonian_iterator_end)
    return
  end subroutine base_hamiltonian_iterator_end

end module base_hamiltonian_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
