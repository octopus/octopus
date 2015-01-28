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

#define HASH_TEMPLATE_NAME bhmlt
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bhmlt
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

module bhmlt_m

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

  use json_m,   only: JSON_OK, json_object_t, json_object_iterator_t
  use json_m,   only: json_init, json_get, json_next, json_end

  use simulation_m, only: &
    simulation_t

  use bsyst_m, only:     &
    system_t => bsyst_t

  use bterm_m, only:         &
    term_t    => bterm_t,    &
    term_init => bterm_init, &
    term_copy => bterm_copy, &
    term_end  => bterm_end

  use bpotn_m, only:                  &
    potential_t      => bpotn_t,      &
    potential_init   => bpotn_init,   &
    potential_start  => bpotn_start,  &
    potential_update => bpotn_update, &
    potential_stop   => bpotn_stop,   &
    potential_copy   => bpotn_copy,   &
    potential_end    => bpotn_end

  implicit none

  private
  public ::                &
    bhmlt_init,            &
    bhmlt_start,           &
    bhmlt_update,          &
    bhmlt_stop,            &
    bhmlt_get,             &
    bhmlt_setn,            &
    bhmlt_getn,            &
    bhmlt_deln,            &
    bhmlt_copy,            &
    bhmlt_end

#define HASH_TEMPLATE_NAME bhmlt
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bhmlt
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

  integer, parameter :: NAME_LEN = 63

  integer, public, parameter :: NONE_TYPE = 0
  integer, public, parameter :: TERM_TYPE = 1
  integer, public, parameter :: POTN_TYPE = 2
  integer, public, parameter :: HMLT_TYPE = 3
  
  type, private :: hterm_t
    private
    type(term_t),      pointer :: term =>null()
    type(potential_t), pointer :: potn =>null()
    type(bhmlt_t),     pointer :: hmlt =>null()
    integer                    :: type = NONE_TYPE
  end type hterm_t

  type, public :: bhmlt_t
    private
    type(json_object_t), pointer :: config => null()
    type(system_t),      pointer :: sys    => null()
    type(simulation_t),  pointer :: sim    => null()
    type(hterm_dict_t)           :: dict
    type(bhmlt_hash_t)           :: hash
  end type bhmlt_t

  interface bhmlt_init
    module procedure bhmlt_init_bhmlt
    module procedure bhmlt_init_build
  end interface bhmlt_init

  interface bhmlt_get
    module procedure bhmlt_get_config
    module procedure bhmlt_get_system
    module procedure bhmlt_get_simulation
  end interface bhmlt_get

  interface bhmlt_setn
    module procedure bhmlt_setn_term
    module procedure bhmlt_setn_potn
    module procedure bhmlt_setn_hmlt
  end interface bhmlt_setn

  interface bhmlt_getn
    module procedure bhmlt_getn_term
    module procedure bhmlt_getn_potn
    module procedure bhmlt_getn_hmlt
  end interface bhmlt_getn

  interface bhmlt_deln
    module procedure bhmlt_deln_term
    module procedure bhmlt_deln_potn
    module procedure bhmlt_deln_hmlt
  end interface bhmlt_deln

  interface hterm_init
    module procedure hterm_init_term
    module procedure hterm_init_potn
    module procedure hterm_init_hmlt
  end interface hterm_init

  interface hterm_get
    module procedure hterm_get_term
    module procedure hterm_get_potn
    module procedure hterm_get_hmlt
  end interface hterm_get

contains

#define HASH_TEMPLATE_NAME bhmlt
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bhmlt
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
  subroutine hterm_init_term(this, that)
    type(hterm_t),        intent(inout) :: this
    type(term_t), target, intent(in)    :: that
    !
    PUSH_SUB(hterm_init_term)
    ASSERT(this%type==NONE_TYPE)
    this%type=TERM_TYPE
    this%term=>that
    POP_SUB(hterm_init_term)
    return
  end subroutine hterm_init_term

  ! ---------------------------------------------------------
  subroutine hterm_init_potn(this, that)
    type(hterm_t),             intent(inout) :: this
    type(potential_t), target, intent(in)    :: that
    !
    PUSH_SUB(hterm_init_potn)
    ASSERT(this%type==NONE_TYPE)
    this%type=POTN_TYPE
    this%potn=>that
    POP_SUB(hterm_init_potn)
    return
  end subroutine hterm_init_potn

  ! ---------------------------------------------------------
  subroutine hterm_init_hmlt(this, that)
    type(hterm_t),          intent(inout) :: this
    type(bhmlt_t), target, intent(in)    :: that
    !
    PUSH_SUB(hterm_init_hmlt)
    ASSERT(this%type==NONE_TYPE)
    this%type=HMLT_TYPE
    this%hmlt=>that
    POP_SUB(hterm_init_hmlt)
    return
  end subroutine hterm_init_hmlt

  ! ---------------------------------------------------------
  recursive subroutine hterm_start(this, sim)
    type(hterm_t),      intent(inout) :: this
    type(simulation_t), intent(in)    :: sim
    !
    PUSH_SUB(hterm_start)
    select case(this%type)
    case(NONE_TYPE)
    case(TERM_TYPE)
    case(POTN_TYPE)
      call potential_start(this%potn, sim)
    case(HMLT_TYPE)
      call bhmlt_start(this%hmlt, sim)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm_start)
    return
  end subroutine hterm_start

  ! ---------------------------------------------------------
  recursive subroutine hterm_update(this)
    type(hterm_t), intent(inout) :: this
    !
    integer :: ierr
    !
    PUSH_SUB(hterm_update)
    select case(this%type)
    case(NONE_TYPE)
    case(TERM_TYPE)
    case(POTN_TYPE)
      call potential_update(this%potn)
    case(HMLT_TYPE)
      call bhmlt_update(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm_update)
    return
  end subroutine hterm_update

  ! ---------------------------------------------------------
  recursive subroutine hterm_stop(this)
    type(hterm_t), intent(inout) :: this
    !
    PUSH_SUB(hterm_stop)
    select case(this%type)
    case(NONE_TYPE)
    case(TERM_TYPE)
    case(POTN_TYPE)
      call potential_stop(this%potn)
    case(HMLT_TYPE)
      call bhmlt_stop(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    POP_SUB(hterm_stop)
    return
  end subroutine hterm_stop

  ! ---------------------------------------------------------
  subroutine hterm_get_term(this, that)
    type(hterm_t), intent(in) :: this
    type(term_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_term)
    ASSERT(this%type==TERM_TYPE)
    ASSERT(associated(this%term))
    that=>this%term
    POP_SUB(hterm_get_term)
    return
  end subroutine hterm_get_term

  ! ---------------------------------------------------------
  subroutine hterm_get_potn(this, that)
    type(hterm_t),      intent(in) :: this
    type(potential_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_potn)
    ASSERT(this%type==POTN_TYPE)
    ASSERT(associated(this%potn))
    that=>this%potn
    POP_SUB(hterm_get_potn)
    return
  end subroutine hterm_get_potn

  ! ---------------------------------------------------------
  subroutine hterm_get_hmlt(this, that)
    type(hterm_t),  intent(in) :: this
    type(bhmlt_t), pointer     :: that
    !
    PUSH_SUB(hterm_get_hmlt)
    ASSERT(this%type==HMLT_TYPE)
    ASSERT(associated(this%hmlt))
    that=>this%hmlt
    POP_SUB(hterm_get_hmlt)
    return
  end subroutine hterm_get_hmlt

  ! ---------------------------------------------------------
  subroutine hterm_copy(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm_copy)
    call hterm_end(this)
    this%type=that%type
    select case(this%type)
    case(TERM_TYPE)
      this%term=>that%term
    case(POTN_TYPE)
      this%potn=>that%potn
    case(HMLT_TYPE)
      this%hmlt=>that%hmlt
    end select
    POP_SUB(hterm_copy)
    return
  end subroutine hterm_copy

  ! ---------------------------------------------------------
  subroutine hterm_end(this)
    type(hterm_t), intent(inout) :: this
    !
    nullify(this%term,this%potn,this%hmlt)
    this%type=NONE_TYPE
    return
  end subroutine hterm_end

  ! ---------------------------------------------------------
  subroutine bhmlt_init_bhmlt(this, sys, config)
    type(bhmlt_t),               intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(bhmlt_init_bhmlt)
    this%config=>config
    this%sys=>sys
    nullify(this%sim)
    call hterm_dict_init(this%dict)
    call bhmlt_hash_init(this%hash)
    POP_SUB(bhmlt_init_bhmlt)
    return
  end subroutine bhmlt_init_bhmlt

  ! ---------------------------------------------------------
  recursive subroutine bhmlt_init_build(this, that, config)
    type(bhmlt_t),       intent(inout) :: this
    type(bhmlt_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bhmlt_init_build)
    call bhmlt_hash_set(this%hash, config, that)
    POP_SUB(bhmlt_init_build)
    return
  end subroutine bhmlt_init_build

  ! ---------------------------------------------------------
  recursive subroutine bhmlt_start(this, sim)
    type(bhmlt_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(bhmlt_start)
    this%sim=>sim
    call hterm_dict_init(iter, this%dict)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_start(htrm, sim)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(bhmlt_start)
    return
  end subroutine bhmlt_start

  ! ---------------------------------------------------------
  recursive subroutine bhmlt_update(this)
    type(bhmlt_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(bhmlt_update)
    call hterm_dict_init(iter, this%dict)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_update(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(bhmlt_update)
    return
  end subroutine bhmlt_update

  ! ---------------------------------------------------------
  recursive subroutine bhmlt_stop(this)
    type(bhmlt_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(bhmlt_stop)
    call hterm_dict_init(iter, this%dict)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_stop(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(bhmlt_stop)
    return
  end subroutine bhmlt_stop

  ! ---------------------------------------------------------
  subroutine bhmlt_get_config(this, that)
    type(bhmlt_t),        intent(in) :: this
    type(json_object_t), pointer     :: that
    !
    PUSH_SUB(bhmlt_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bhmlt_get_config)
    return
  end subroutine bhmlt_get_config

  ! ---------------------------------------------------------
  subroutine bhmlt_get_system(this, that)
    type(bhmlt_t),   intent(in) :: this
    type(system_t), pointer     :: that
    !
    PUSH_SUB(bhmlt_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(bhmlt_get_system)
    return
  end subroutine bhmlt_get_system

  ! ---------------------------------------------------------
  subroutine bhmlt_get_simulation(this, that)
    type(bhmlt_t),       intent(in) :: this
    type(simulation_t), pointer     :: that
    !
    PUSH_SUB(bhmlt_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(bhmlt_get_simulation)
    return
  end subroutine bhmlt_get_simulation

  ! ---------------------------------------------------------
  subroutine bhmlt_setn_hterm(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(hterm_t),    intent(in)    :: that
    !
    PUSH_SUB(bhmlt_setn_hterm)
    call hterm_dict_set(this%dict, name, that)
    POP_SUB(bhmlt_setn_hterm)
    return
  end subroutine bhmlt_setn_hterm

  ! ---------------------------------------------------------
  subroutine bhmlt_setn_term(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(term_t),     intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_setn_term)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call bhmlt_setn_hterm(this, name, htrm)
    POP_SUB(bhmlt_setn_term)
    return
  end subroutine bhmlt_setn_term

  ! ---------------------------------------------------------
  subroutine bhmlt_setn_potn(this, name, that)
    type(bhmlt_t),     intent(inout) :: this
    character(len=*),  intent(in)    :: name
    type(potential_t), intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_setn_potn)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call bhmlt_setn_hterm(this, name, htrm)
    nullify(htrm)
    POP_SUB(bhmlt_setn_potn)
    return
  end subroutine bhmlt_setn_potn

  ! ---------------------------------------------------------
  subroutine bhmlt_setn_hmlt(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(bhmlt_t),    intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_setn_hmlt)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call bhmlt_setn_hterm(this, name, htrm)
    nullify(htrm)
    POP_SUB(bhmlt_setn_hmlt)
    return
  end subroutine bhmlt_setn_hmlt

  ! ---------------------------------------------------------
  subroutine bhmlt_getn_hterm(this, name, that)
    type(bhmlt_t),    intent(in) :: this
    character(len=*), intent(in) :: name
    type(hterm_t),   pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(bhmlt_getn_hterm)
    nullify(that)
    call hterm_dict_get(this%dict, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(bhmlt_getn_hterm)
    return
  end subroutine bhmlt_getn_hterm

  ! ---------------------------------------------------------
  subroutine bhmlt_getn_term(this, name, that)
    type(bhmlt_t),    intent(in) :: this
    character(len=*), intent(in) :: name
    type(term_t),    pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_getn_term)
    nullify(that, htrm)
    call bhmlt_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(bhmlt_getn_term)
    return
  end subroutine bhmlt_getn_term

  ! ---------------------------------------------------------
  subroutine bhmlt_getn_potn(this, name, that)
    type(bhmlt_t),      intent(in) :: this
    character(len=*),   intent(in) :: name
    type(potential_t), pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_getn_potn)
    nullify(that, htrm)
    call bhmlt_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get_potn(htrm, that)
    nullify(htrm)
    POP_SUB(bhmlt_getn_potn)
    return
  end subroutine bhmlt_getn_potn

  ! ---------------------------------------------------------
  subroutine bhmlt_getn_hmlt(this, name, that)
    type(bhmlt_t),    intent(in) :: this
    character(len=*), intent(in) :: name
    type(bhmlt_t),   pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_getn_hmlt)
    nullify(that, htrm)
    call bhmlt_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(bhmlt_getn_hmlt)
    return
  end subroutine bhmlt_getn_hmlt

  ! ---------------------------------------------------------
  subroutine bhmlt_deln_hterm(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(hterm_t),   pointer        :: that
    !
    integer :: ierr
    !
    PUSH_SUB(bhmlt_deln_hterm)
    nullify(that)
    call hterm_dict_del(this%dict, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(bhmlt_deln_hterm)
    return
  end subroutine bhmlt_deln_hterm

  ! ---------------------------------------------------------
  subroutine bhmlt_deln_term(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(term_t),    pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_deln_term)
    nullify(that, htrm)
    call bhmlt_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(bhmlt_deln_term)
    return
  end subroutine bhmlt_deln_term

  ! ---------------------------------------------------------
  subroutine bhmlt_deln_potn(this, name, that)
    type(bhmlt_t),      intent(inout) :: this
    character(len=*),   intent(in)    :: name
    type(potential_t), pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_deln_potn)
    nullify(that, htrm)
    call bhmlt_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(bhmlt_deln_potn)
    return
  end subroutine bhmlt_deln_potn

  ! ---------------------------------------------------------
  subroutine bhmlt_deln_hmlt(this, name, that)
    type(bhmlt_t),    intent(inout) :: this
    character(len=*), intent(in)    :: name
    type(bhmlt_t),   pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(bhmlt_deln_hmlt)
    nullify(that, htrm)
    call bhmlt_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(bhmlt_deln_hmlt)
    return
  end subroutine bhmlt_deln_hmlt

  ! ---------------------------------------------------------
  subroutine bhmlt_copy(this, that)
    type(bhmlt_t), intent(inout) :: this
    type(bhmlt_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: ihtr, ohtr
    character(len=NAME_LEN)     :: name
    integer                     :: ierr
    !
    PUSH_SUB(bhmlt_copy)
    call bhmlt_end(this)
    this%config=>that%config
    this%sys=>that%sys
    this%sim=>that%sim
    call hterm_dict_init(this%dict, hterm_dict_len(that%dict))
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
    call bhmlt_hash_copy(this%hash, that%hash)
    POP_SUB(bhmlt_copy)
    return
  end subroutine bhmlt_copy

  ! ---------------------------------------------------------
  recursive subroutine bhmlt_end(this)
    type(bhmlt_t), intent(inout) :: this
    !
    type(hterm_t), pointer :: htrm
    integer                :: ierr
    !
    PUSH_SUB(bhmlt_end)
    nullify(this%config, this%sys, this%sim)
    do
      nullify(htrm)
      call hterm_dict_pop(this%dict, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
    end do
    nullify(htrm)
    call hterm_dict_end(this%dict)
    call bhmlt_hash_end(this%hash)
    POP_SUB(bhmlt_end)
    return
  end subroutine bhmlt_end

end module bhmlt_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
