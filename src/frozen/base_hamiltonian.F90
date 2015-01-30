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

  use json_m,   only: JSON_OK, json_object_t, json_object_iterator_t
  use json_m,   only: json_init, json_get, json_next, json_end

  use simulation_m, only: &
    simulation_t

  use base_system_m, only:     &
    system_t => base_system_t

  use base_term_m, only:         &
    term_t    => base_term_t,    &
    term_init => base_term_init, &
    term_copy => base_term_copy, &
    term_end  => base_term_end

  use base_potential_m, only:                  &
    potential_t      => base_potential_t,      &
    potential_init   => base_potential_init,   &
    potential_start  => base_potential_start,  &
    potential_update => base_potential_update, &
    potential_stop   => base_potential_stop,   &
    potential_copy   => base_potential_copy,   &
    potential_end    => base_potential_end

  implicit none

  private
  public ::                  &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_get,    &
    base_hamiltonian_setn,   &
    base_hamiltonian_getn,   &
    base_hamiltonian_deln,   &
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

  integer, parameter :: HMLT_NAME_LEN = 63

  integer, public, parameter :: HMLT_TYPE_NONE = 0
  integer, public, parameter :: HMLT_TYPE_TERM = 1
  integer, public, parameter :: HMLT_TYPE_POTN = 2
  integer, public, parameter :: HMLT_TYPE_HMLT = 3
  
  type, private :: hterm_t
    private
    type(term_t),             pointer :: term =>null()
    type(potential_t),        pointer :: potn =>null()
    type(base_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: type = HMLT_TYPE_NONE
  end type hterm_t

  type, public :: base_hamiltonian_t
    private
    type(json_object_t),  pointer :: config => null()
    type(system_t),       pointer :: sys    => null()
    type(simulation_t),   pointer :: sim    => null()
    type(hterm_dict_t)            :: dict
    type(base_hamiltonian_hash_t) :: hash
  end type base_hamiltonian_t

  type, public :: base_hamiltonian_iterator_t
    private
    type(base_hamiltonian_t),      pointer :: self =>null()
    type(base_hamiltonian_hash_iterator_t) :: iter
  end type base_hamiltonian_iterator_t

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

  interface base_hamiltonian_copy
    module procedure base_hamiltonian_copy_hamiltonian
    module procedure base_hamiltonian_iterator_copy
  end interface base_hamiltonian_copy

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_hamiltonian
    module procedure base_hamiltonian_init_build
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
  end interface base_hamiltonian_get

  interface base_hamiltonian_setn
    module procedure base_hamiltonian_setn_term
    module procedure base_hamiltonian_setn_potn
    module procedure base_hamiltonian_setn_hmlt
  end interface base_hamiltonian_setn

  interface base_hamiltonian_getn
    module procedure base_hamiltonian_getn_term
    module procedure base_hamiltonian_getn_potn
    module procedure base_hamiltonian_getn_hmlt
  end interface base_hamiltonian_getn

  interface base_hamiltonian_deln
    module procedure base_hamiltonian_deln_term
    module procedure base_hamiltonian_deln_potn
    module procedure base_hamiltonian_deln_hmlt
  end interface base_hamiltonian_deln

  interface base_hamiltonian_end
    module procedure base_hamiltonian_end_hamiltonian
    module procedure base_hamiltonian_iterator_end
  end interface base_hamiltonian_end

  integer, parameter :: BASE_HAMILTONIAN_OK          = BASE_HAMILTONIAN_HASH_OK
  integer, parameter :: BASE_HAMILTONIAN_KEY_ERROR   = BASE_HAMILTONIAN_HASH_KEY_ERROR
  integer, parameter :: BASE_HAMILTONIAN_EMPTY_ERROR = BASE_HAMILTONIAN_HASH_EMPTY_ERROR

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
  subroutine hterm_init_term(this, that)
    type(hterm_t),        intent(inout) :: this
    type(term_t), target, intent(in)    :: that
    !
    PUSH_SUB(hterm_init_term)
    ASSERT(this%type==HMLT_TYPE_NONE)
    this%type=HMLT_TYPE_TERM
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
    ASSERT(this%type==HMLT_TYPE_NONE)
    this%type=HMLT_TYPE_POTN
    this%potn=>that
    POP_SUB(hterm_init_potn)
    return
  end subroutine hterm_init_potn

  ! ---------------------------------------------------------
  subroutine hterm_init_hmlt(this, that)
    type(hterm_t),                    intent(inout) :: this
    type(base_hamiltonian_t), target, intent(in)    :: that
    !
    PUSH_SUB(hterm_init_hmlt)
    ASSERT(this%type==HMLT_TYPE_NONE)
    this%type=HMLT_TYPE_HMLT
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
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call potential_start(this%potn, sim)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian_start(this%hmlt, sim)
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
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call potential_update(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian_update(this%hmlt)
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
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call potential_stop(this%potn)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian_stop(this%hmlt)
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
    ASSERT(this%type==HMLT_TYPE_TERM)
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
  subroutine hterm_copy(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that
    !
    PUSH_SUB(hterm_copy)
    call hterm_end(this)
    this%type=that%type
    select case(this%type)
    case(HMLT_TYPE_TERM)
      this%term=>that%term
    case(HMLT_TYPE_POTN)
      this%potn=>that%potn
    case(HMLT_TYPE_HMLT)
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
    this%type=HMLT_TYPE_NONE
    return
  end subroutine hterm_end

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_init_hamiltonian(this, sys, config)
    type(base_hamiltonian_t),    intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(base_hamiltonian_init_hamiltonian)
    this%config=>config
    this%sys=>sys
    nullify(this%sim)
    call hterm_dict_init(this%dict)
    call base_hamiltonian_hash_init(this%hash)
    POP_SUB(base_hamiltonian_init_hamiltonian)
    return
  end subroutine base_hamiltonian_init_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_init_build(this, that, config)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    type(json_object_t),      intent(in)    :: config
    !
    PUSH_SUB(base_hamiltonian_init_build)
    call base_hamiltonian_hash_set(this%hash, config, that)
    POP_SUB(base_hamiltonian_init_build)
    return
  end subroutine base_hamiltonian_init_build

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_start(this, sim)
    type(base_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian_start)
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
    POP_SUB(base_hamiltonian_start)
    return
  end subroutine base_hamiltonian_start

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_update(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian_update)
    call hterm_dict_init(iter, this%dict)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_update(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(base_hamiltonian_update)
    return
  end subroutine base_hamiltonian_update

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_stop(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr
    !
    PUSH_SUB(base_hamiltonian_stop)
    call hterm_dict_init(iter, this%dict)
    do
      nullify(htrm)
      call hterm_dict_next(iter, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_stop(htrm)
    end do
    call hterm_dict_end(iter)
    nullify(htrm)
    POP_SUB(base_hamiltonian_stop)
    return
  end subroutine base_hamiltonian_stop

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
    type(system_t),          pointer     :: that
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
  subroutine base_hamiltonian_setn_hterm(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(hterm_t),            intent(in)    :: that
    !
    PUSH_SUB(base_hamiltonian_setn_hterm)
    call hterm_dict_set(this%dict, name, that)
    POP_SUB(base_hamiltonian_setn_hterm)
    return
  end subroutine base_hamiltonian_setn_hterm

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_setn_term(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(term_t),             intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_setn_term)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call base_hamiltonian_setn_hterm(this, name, htrm)
    POP_SUB(base_hamiltonian_setn_term)
    return
  end subroutine base_hamiltonian_setn_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_setn_potn(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(potential_t),        intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_setn_potn)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call base_hamiltonian_setn_hterm(this, name, htrm)
    nullify(htrm)
    POP_SUB(base_hamiltonian_setn_potn)
    return
  end subroutine base_hamiltonian_setn_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_setn_hmlt(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_setn_hmlt)
    nullify(htrm)
    SAFE_ALLOCATE(htrm)
    call hterm_init(htrm, that)
    call base_hamiltonian_setn_hterm(this, name, htrm)
    nullify(htrm)
    POP_SUB(base_hamiltonian_setn_hmlt)
    return
  end subroutine base_hamiltonian_setn_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_getn_hterm(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(hterm_t),           pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(base_hamiltonian_getn_hterm)
    nullify(that)
    call hterm_dict_get(this%dict, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(base_hamiltonian_getn_hterm)
    return
  end subroutine base_hamiltonian_getn_hterm

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_getn_term(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(term_t),            pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_getn_term)
    nullify(that, htrm)
    call base_hamiltonian_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_getn_term)
    return
  end subroutine base_hamiltonian_getn_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_getn_potn(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(potential_t),       pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_getn_potn)
    nullify(that, htrm)
    call base_hamiltonian_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get_potn(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_getn_potn)
    return
  end subroutine base_hamiltonian_getn_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_getn_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_getn_hmlt)
    nullify(that, htrm)
    call base_hamiltonian_getn_hterm(this, name, htrm)
    if(associated(htrm))&
      call hterm_get(htrm, that)
    nullify(htrm)
    POP_SUB(base_hamiltonian_getn_hmlt)
    return
  end subroutine base_hamiltonian_getn_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_deln_hterm(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(hterm_t),           pointer        :: that
    !
    integer :: ierr
    !
    PUSH_SUB(base_hamiltonian_deln_hterm)
    nullify(that)
    call hterm_dict_del(this%dict, name, that, ierr)
    if(ierr/=HTERM_DICT_OK)nullify(that)
    POP_SUB(base_hamiltonian_deln_hterm)
    return
  end subroutine base_hamiltonian_deln_hterm

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_deln_term(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(term_t),            pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_deln_term)
    nullify(that, htrm)
    call base_hamiltonian_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(base_hamiltonian_deln_term)
    return
  end subroutine base_hamiltonian_deln_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_deln_potn(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(potential_t),       pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_deln_potn)
    nullify(that, htrm)
    call base_hamiltonian_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(base_hamiltonian_deln_potn)
    return
  end subroutine base_hamiltonian_deln_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_deln_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: name
    type(base_hamiltonian_t), pointer        :: that
    !
    type(hterm_t), pointer :: htrm
    !
    PUSH_SUB(base_hamiltonian_deln_hmlt)
    nullify(that, htrm)
    call base_hamiltonian_deln_hterm(this, name, htrm)
    if(associated(htrm))then
      call hterm_get(htrm, that)
      call hterm_end(htrm)
      SAFE_DEALLOCATE_P(htrm)
      nullify(htrm)
    end if
    POP_SUB(base_hamiltonian_deln_hmlt)
    return
  end subroutine base_hamiltonian_deln_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_copy_hamiltonian(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that
    !
    type(hterm_dict_iterator_t)  :: iter
    type(hterm_t),       pointer :: ihtr, ohtr
    character(len=HMLT_NAME_LEN) :: name
    integer                      :: ierr
    !
    PUSH_SUB(base_hamiltonian_copy_hamiltonian)
    call base_hamiltonian_end(this)
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
    call base_hamiltonian_hash_copy(this%hash, that%hash)
    POP_SUB(base_hamiltonian_copy_hamiltonian)
    return
  end subroutine base_hamiltonian_copy_hamiltonian

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_end_hamiltonian(this)
    type(base_hamiltonian_t), intent(inout) :: this
    !
    type(hterm_t), pointer :: htrm
    integer                :: ierr
    !
    PUSH_SUB(base_hamiltonian_end_hamiltonian)
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
    call base_hamiltonian_hash_end(this%hash)
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
    call base_hamiltonian_hash_init(this%iter, that%hash)
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
    call base_hamiltonian_hash_next(this%iter, config, that, ierr)
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
    call base_hamiltonian_hash_next(this%iter, that, ierr)
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
    call base_hamiltonian_hash_next(this%iter, that, ierr)
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
    call base_hamiltonian_hash_copy(this%iter, that%iter)
    POP_SUB(base_hamiltonian_iterator_copy)
    return
  end subroutine base_hamiltonian_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_iterator_end(this)
    type(base_hamiltonian_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_hamiltonian_iterator_end)
    nullify(this%self)
    call base_hamiltonian_hash_end(this%iter)
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
