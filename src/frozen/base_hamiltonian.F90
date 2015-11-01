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

module base_hamiltonian_m

  use base_functional_m
  use base_potential_m
  use base_system_m
  use base_term_m
  use config_dict_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use storage_m

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

#define TEMPLATE_PREFIX base_hamiltonian
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::         &
    HMLT_TYPE_NONE, &
    HMLT_TYPE_TERM, &
    HMLT_TYPE_POTN, &
    HMLT_TYPE_FNCT, &
    HMLT_TYPE_HMLT

  public ::                       &
    BASE_HAMILTONIAN_OK,          &
    BASE_HAMILTONIAN_KEY_ERROR,   &
    BASE_HAMILTONIAN_EMPTY_ERROR

  public ::             &
    base_hamiltonian_t

  public ::                     &
    base_hamiltonian__new__,    &
    base_hamiltonian__del__,    &
    base_hamiltonian__init__,   &
    base_hamiltonian__start__,  &
    base_hamiltonian__update__, &
    base_hamiltonian__stop__,   &
    base_hamiltonian__reset__,  &
    base_hamiltonian__acc__,    &
    base_hamiltonian__sub__,    &
    base_hamiltonian__copy__,   &
    base_hamiltonian__end__

  public ::                  &
    base_hamiltonian_new,    &
    base_hamiltonian_del,    &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_sets,   &
    base_hamiltonian_gets,   &
    base_hamiltonian_set,    &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

  integer, parameter :: HMLT_TYPE_NONE = 0 
  integer, parameter :: HMLT_TYPE_TERM = 1
  integer, parameter :: HMLT_TYPE_POTN = 2
  integer, parameter :: HMLT_TYPE_FNCT = 3
  integer, parameter :: HMLT_TYPE_HMLT = 4

  integer, parameter :: BASE_HAMILTONIAN_OK          = BASE_HAMILTONIAN_HASH_OK
  integer, parameter :: BASE_HAMILTONIAN_KEY_ERROR   = BASE_HAMILTONIAN_HASH_KEY_ERROR
  integer, parameter :: BASE_HAMILTONIAN_EMPTY_ERROR = BASE_HAMILTONIAN_HASH_EMPTY_ERROR

  type :: base_hamiltonian_raii_t
    private
    type(base_hamiltonian_t), pointer :: prnt =>null()
    type(base_hamiltonian_list_t)     :: list
  end type base_hamiltonian_raii_t

  type :: hterm_t
    private
    type(base_term_t),        pointer :: term =>null()
    type(base_potential_t),   pointer :: potn =>null()
    type(base_functional_t),  pointer :: fnct =>null()
    type(base_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: type = HMLT_TYPE_NONE
  end type hterm_t

  type :: base_hamiltonian_t
    private
    type(json_object_t),      pointer :: config =>null()
    type(base_system_t),      pointer :: sys    =>null()
    type(simulation_t),       pointer :: sim    =>null()
    real(kind=wp)                     :: energy = 0.0_wp
    type(storage_t)                   :: data
    type(hterm_dict_t)                :: hdct
    type(config_dict_t)               :: dict
    type(base_hamiltonian_hash_t)     :: hash
    type(base_hamiltonian_raii_t)     :: raii
  end type base_hamiltonian_t

  interface hterm__init__
    module procedure hterm__init__type
    module procedure hterm__init__copy
  end interface hterm__init__

  interface hterm__get__
    module procedure hterm__get__term
    module procedure hterm__get__potn
    module procedure hterm__get__fnct
    module procedure hterm__get__hmlt
  end interface hterm__get__

  interface base_hamiltonian__init__
    module procedure base_hamiltonian__init__type
    module procedure base_hamiltonian__init__copy
  end interface base_hamiltonian__init__

  interface base_hamiltonian__new__
    module procedure base_hamiltonian__new__term
    module procedure base_hamiltonian__new__potn
    module procedure base_hamiltonian__new__fnct
    module procedure base_hamiltonian__new__hmlt
  end interface base_hamiltonian__new__

  interface base_hamiltonian__acc__
    module procedure base_hamiltonian__acc__term
    module procedure base_hamiltonian__acc__potn
    module procedure base_hamiltonian__acc__fnct
    module procedure base_hamiltonian__acc__hmlt
  end interface base_hamiltonian__acc__

  interface base_hamiltonian__sub__
    module procedure base_hamiltonian__sub__term
    module procedure base_hamiltonian__sub__potn
    module procedure base_hamiltonian__sub__fnct
    module procedure base_hamiltonian__sub__hmlt
  end interface base_hamiltonian__sub__

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_type
    module procedure base_hamiltonian_init_copy
  end interface base_hamiltonian_init

  interface base_hamiltonian_set
    module procedure base_hamiltonian_set_info
  end interface base_hamiltonian_set

  interface base_hamiltonian_gets
    module procedure base_hamiltonian_gets_config
    module procedure base_hamiltonian_gets_name
  end interface base_hamiltonian_gets

  interface base_hamiltonian_get
    module procedure base_hamiltonian_get_info
    module procedure base_hamiltonian_get_config
    module procedure base_hamiltonian_get_system
    module procedure base_hamiltonian_get_simulation
    module procedure base_hamiltonian_get_term
    module procedure base_hamiltonian_get_potn
    module procedure base_hamiltonian_get_fnct
    module procedure base_hamiltonian_get_hmlt
    module procedure base_hamiltonian_get_hamiltonian_1d
    module procedure base_hamiltonian_get_hamiltonian_md
  end interface base_hamiltonian_get

  interface base_hamiltonian_copy
    module procedure base_hamiltonian_copy_type
  end interface base_hamiltonian_copy

  interface base_hamiltonian_end
    module procedure base_hamiltonian_end_type
  end interface base_hamiltonian_end

#define TEMPLATE_PREFIX base_hamiltonian
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_hamiltonian
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_hamiltonian
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_hamiltonian
#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY
#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME hterm
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__raii_init__(this, that)
    type(base_hamiltonian_raii_t),              intent(out) :: this
    type(base_hamiltonian_t), optional, target, intent(in)  :: that

    PUSH_SUB(base_hamiltonian__raii_init__)

    nullify(this%prnt)
    if(present(that)) this%prnt => that
    call base_hamiltonian_list_init(this%list)

    POP_SUB(base_hamiltonian__raii_init__)
  end subroutine base_hamiltonian__raii_init__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__rpush__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    PUSH_SUB(base_hamiltonian__rpush__)

    call base_hamiltonian_list_push(this%raii%list, that)

    POP_SUB(base_hamiltonian__rpush__)
  end subroutine base_hamiltonian__rpush__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__rpop__(this, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    type(base_hamiltonian_t), pointer        :: that

    PUSH_SUB(base_hamiltonian__rpop__)

    call base_hamiltonian_list_pop(this%raii%list, that)

    POP_SUB(base_hamiltonian__rpop__)
  end subroutine base_hamiltonian__rpop__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__rdel__(this, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    type(base_hamiltonian_t), pointer        :: that

    PUSH_SUB(base_hamiltonian__rdel__)

    call base_hamiltonian_list_del(this%raii%list, that)

    POP_SUB(base_hamiltonian__rdel__)
  end subroutine base_hamiltonian__rdel__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__rset__(this, that)
    type(base_hamiltonian_t),         intent(inout) :: this
    type(base_hamiltonian_t), target, intent(in)    :: that

    PUSH_SUB(base_hamiltonian__rset__)

    this%raii%prnt => that

    POP_SUB(base_hamiltonian__rset__)
  end subroutine base_hamiltonian__rset__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__rget__(this, that)
    type(base_hamiltonian_t),  intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that

    PUSH_SUB(base_hamiltonian__rget__)

    that => this%raii%prnt

    POP_SUB(base_hamiltonian__rget__)
  end subroutine base_hamiltonian__rget__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__raii_end__(this)
    type(base_hamiltonian_raii_t), intent(inout) :: this

    type(base_hamiltonian_t), pointer :: subs

    PUSH_SUB(base_hamiltonian__raii_end__)

    nullify(this%prnt)
    do
      nullify(subs)
      call base_hamiltonian_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_hamiltonian__raii_end__(subs%raii)
      !call base_hamiltonian__raii_del__(subs)
      SAFE_DEALLOCATE_P(subs)
      nullify(subs)
    end do
    nullify(subs)
    call base_hamiltonian_list_end(this%list)

    POP_SUB(base_hamiltonian__raii_end__)
  end subroutine base_hamiltonian__raii_end__

  ! ---------------------------------------------------------
  subroutine hterm__new__(this)
    type(hterm_t), pointer :: this

    PUSH_SUB(hterm__new__)

    nullify(this)
    SAFE_ALLOCATE(this)
    call hterm__inull__(this)

    POP_SUB(hterm__new__)
  end subroutine hterm__new__

  ! ---------------------------------------------------------
  subroutine hterm__del__(this)
    type(hterm_t), pointer :: this

    PUSH_SUB(hterm__del__)

    if(associated(this))then
      call hterm__end__(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(hterm__del__)
  end subroutine hterm__del__

  ! ---------------------------------------------------------
  subroutine hterm__inull__(this)
    type(hterm_t), intent(inout) :: this

    PUSH_SUB(hterm__inull__)

    nullify(this%term, this%potn, this%fnct, this%hmlt)
    this%type = HMLT_TYPE_NONE

    POP_SUB(hterm__inull__)
  end subroutine hterm__inull__

  ! ---------------------------------------------------------
  subroutine hterm__iinit__(this, type)
    type(hterm_t), intent(out) :: this
    integer,       intent(in)  :: type

    PUSH_SUB(hterm__iinit__)

    call hterm__inull__(this)
    this%type = type
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
  end subroutine hterm__iinit__

  ! ---------------------------------------------------------
  recursive subroutine hterm__init__type(this, sys, config)
    type(hterm_t),               intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t),         intent(in)  :: config

    integer :: type, ierr

    PUSH_SUB(hterm__init__type)

    call json_get(config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    call hterm__iinit__(this, type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__init__(this%term, sys, config)
    case(HMLT_TYPE_POTN)
      call base_potential__init__(this%potn, sys, config)
    case(HMLT_TYPE_FNCT)
      call base_functional__init__(this%fnct, sys, config)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__init__(this%hmlt, sys, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__init__type)
  end subroutine hterm__init__type

  ! ---------------------------------------------------------
  recursive subroutine hterm__init__copy(this, that)
    type(hterm_t), intent(out) :: this
    type(hterm_t), intent(in)  :: that

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
      call base_hamiltonian__init__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__init__copy)
  end subroutine hterm__init__copy

  ! ---------------------------------------------------------
  recursive subroutine hterm__start__(this, sim)
    type(hterm_t),                intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(hterm__start__)

    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call base_potential__start__(this%potn, sim)
    case(HMLT_TYPE_FNCT)
      call base_functional__start__(this%fnct, sim)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__start__(this%hmlt, sim)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__start__)
  end subroutine hterm__start__

  ! ---------------------------------------------------------
  recursive subroutine hterm__update__(this)
    type(hterm_t), intent(inout) :: this

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
      call base_hamiltonian__update__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__update__)
  end subroutine hterm__update__

  ! ---------------------------------------------------------
  recursive subroutine hterm__stop__(this)
    type(hterm_t), intent(inout) :: this

    PUSH_SUB(hterm__stop__)

    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
    case(HMLT_TYPE_POTN)
      call base_potential__stop__(this%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__stop__(this%fnct)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__stop__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__stop__)
  end subroutine hterm__stop__

  ! ---------------------------------------------------------
  recursive subroutine hterm__reset__(this)
    type(hterm_t), intent(inout) :: this

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
      call base_hamiltonian__reset__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__reset__)
  end subroutine hterm__reset__

  ! ---------------------------------------------------------
  recursive subroutine hterm__acc__(this, that)
    type(hterm_t),       intent(inout) :: this
    type(hterm_t),       intent(in)    :: that

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
      call base_hamiltonian__acc__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__acc__)
  end subroutine hterm__acc__

  ! ---------------------------------------------------------
  recursive subroutine hterm__sub__(this, that)
    type(hterm_t),       intent(inout) :: this
    type(hterm_t),       intent(in)    :: that

    PUSH_SUB(hterm__sub__)

    ASSERT(this%type==that%type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term__sub__(this%term, that%term)
    case(HMLT_TYPE_POTN)
      call base_potential__sub__(this%potn, that%potn)
    case(HMLT_TYPE_FNCT)
      call base_functional__sub__(this%fnct, that%fnct)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian__sub__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__sub__)
  end subroutine hterm__sub__

  ! ---------------------------------------------------------
  recursive subroutine hterm__sets__(this, that, config)
    type(hterm_t),       intent(inout) :: this
    type(hterm_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(hterm__sets__)

    ASSERT(this%type==that%type)
    select case(this%type)
    case(HMLT_TYPE_NONE)
    case(HMLT_TYPE_TERM)
      call base_term_sets(this%term, that%term, config)
    case(HMLT_TYPE_POTN)
      call base_potential_sets(this%potn, that%potn, config)
    case(HMLT_TYPE_FNCT)
      call base_functional_sets(this%fnct, that%fnct, config)
    case(HMLT_TYPE_HMLT)
      call base_hamiltonian_sets(this%hmlt, that%hmlt, config)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__sets__)
  end subroutine hterm__sets__

  ! ---------------------------------------------------------
  subroutine hterm__get__term(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_term_t), pointer     :: that

    PUSH_SUB(hterm__get__term)

    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_TERM)
      ASSERT(associated(this%term))
      that => this%term
    end if

    POP_SUB(hterm__get__term)
  end subroutine hterm__get__term

  ! ---------------------------------------------------------
  subroutine hterm__get__potn(this, that)
    type(hterm_t),      intent(in) :: this
    type(base_potential_t), pointer     :: that

    PUSH_SUB(hterm__get__potn)

    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_POTN)
      ASSERT(associated(this%potn))
      that => this%potn
    end if

    POP_SUB(hterm__get__potn)
  end subroutine hterm__get__potn

  ! ---------------------------------------------------------
  subroutine hterm__get__fnct(this, that)
    type(hterm_t),            intent(in) :: this
    type(base_functional_t), pointer     :: that

    PUSH_SUB(hterm__get__fnct)

    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_FNCT)
      ASSERT(associated(this%fnct))
      that => this%fnct
    end if

    POP_SUB(hterm__get__fnct)
  end subroutine hterm__get__fnct

  ! ---------------------------------------------------------
  subroutine hterm__get__hmlt(this, that)
    type(hterm_t),             intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that

    PUSH_SUB(hterm__get__hmlt)

    nullify(that)
    if(this%type/=HMLT_TYPE_NONE)then
      ASSERT(this%type==HMLT_TYPE_HMLT)
      ASSERT(associated(this%hmlt))
      that => this%hmlt
    end if

    POP_SUB(hterm__get__hmlt)
  end subroutine hterm__get__hmlt

  ! ---------------------------------------------------------
  subroutine hterm__icopy__(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that

    PUSH_SUB(hterm__icopy__)

    call hterm__iend__(this)
    call hterm__iinit__(this, that%type)

    POP_SUB(hterm__icopy__)
  end subroutine hterm__icopy__

  ! ---------------------------------------------------------
  recursive subroutine hterm__copy__(this, that)
    type(hterm_t), intent(inout) :: this
    type(hterm_t), intent(in)    :: that

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
      call base_hamiltonian__copy__(this%hmlt, that%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select

    POP_SUB(hterm__copy__)
  end subroutine hterm__copy__

  ! ---------------------------------------------------------
  subroutine hterm__iend__(this)
    type(hterm_t), intent(inout) :: this

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
  end subroutine hterm__iend__

  ! ---------------------------------------------------------
  recursive subroutine hterm__end__(this)
    type(hterm_t), intent(inout) :: this

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
      call base_hamiltonian__end__(this%hmlt)
    case default
      message(1)="Unknown Hamiltonian term type."
      call messages_fatal(1)
    end select
    call hterm__iend__(this)

    POP_SUB(hterm__end__)
  end subroutine hterm__end__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__inew__(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(hterm_t),           pointer        :: that

    PUSH_SUB(base_hamiltonian__inew__)

    call hterm__new__(that)
    ASSERT(associated(that))
    call hterm_dict_set(this%hdct, trim(adjustl(name)), that)

    POP_SUB(base_hamiltonian__inew__)
  end subroutine base_hamiltonian__inew__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__new__term(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_term_t),       pointer        :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian__new__term)

    nullify(that, htrm)
    call base_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_TERM)
    call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian__new__term)
  end subroutine base_hamiltonian__new__term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__new__potn(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_potential_t),  pointer        :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian__new__potn)

    nullify(that, htrm)
    call base_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_POTN)
    call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian__new__potn)
  end subroutine base_hamiltonian__new__potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__new__fnct(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_functional_t), pointer        :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian__new__fnct)

    nullify(that, htrm)
    call base_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_FNCT)
    call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian__new__fnct)
  end subroutine base_hamiltonian__new__fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__new__hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: name
    type(base_hamiltonian_t), pointer        :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian__new__hmlt)

    nullify(that, htrm)
    call base_hamiltonian__inew__(this, name, htrm)
    ASSERT(associated(htrm))
    call hterm__iinit__(htrm, HMLT_TYPE_HMLT)
    call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian__new__hmlt)
  end subroutine base_hamiltonian__new__hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__del__(this, name)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian__del__)

    nullify(htrm)
    call hterm_dict_del(this%hdct, trim(adjustl(name)), htrm)
    call hterm__del__(htrm)
    nullify(htrm)

    POP_SUB(base_hamiltonian__del__)
  end subroutine base_hamiltonian__del__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_new(this, that)
    type(base_hamiltonian_t),  target, intent(inout) :: this
    type(base_hamiltonian_t), pointer                :: that

    PUSH_SUB(base_hamiltonian_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    call base_hamiltonian__rset__(that, this)
    call base_hamiltonian__rpush__(this, that)

    POP_SUB(base_hamiltonian_new)
  end subroutine base_hamiltonian_new

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__idel__(this)
    type(base_hamiltonian_t), pointer :: this

    PUSH_SUB(base_hamiltonian__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_hamiltonian__idel__)
  end subroutine base_hamiltonian__idel__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del(this)
    type(base_hamiltonian_t), pointer :: this

    type(base_hamiltonian_t), pointer :: prnt

    PUSH_SUB(base_hamiltonian_del)

    nullify(prnt)
    if(associated(this))then
      call base_hamiltonian__rget__(this, prnt)
      if(associated(prnt))then
        call base_hamiltonian__rdel__(prnt, this)
        call base_hamiltonian_end(this)
        call base_hamiltonian__idel__(this)
      end if
      nullify(prnt)
    end if

    POP_SUB(base_hamiltonian_del)
  end subroutine base_hamiltonian_del

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__inull__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__inull__)

    nullify(this%config, this%sys, this%sim)
    this%energy = 0.0_wp

    POP_SUB(base_hamiltonian__inull__)
  end subroutine base_hamiltonian__inull__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iinit__(this, sys, config)
    type(base_hamiltonian_t),    intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    integer :: type, nspin, ierr
    logical :: alloc

    PUSH_SUB(base_hamiltonian__iinit__)

    call base_hamiltonian__inull__(this)
    this%config => config
    this%sys => sys
    call json_get(this%config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(type==HMLT_TYPE_HMLT)
    call base_system_get(this%sys, nspin=nspin)
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK) alloc = .false.
    call storage_init(this%data, nspin, full=.false., allocate=alloc)
    call hterm_dict_init(this%hdct)
    call config_dict_init(this%dict)
    call base_hamiltonian_hash_init(this%hash)
    call base_hamiltonian__raii_init__(this%raii)

    POP_SUB(base_hamiltonian__iinit__)
  end subroutine base_hamiltonian__iinit__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__type(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config

    type(json_object_iterator_t)        :: iter
    type(json_object_t),        pointer :: cnfg
    type(base_system_t),        pointer :: psys
    type(hterm_t),              pointer :: htrm
    character(len=CONFIG_DICT_NAME_LEN) :: attr, sysn
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__init__type)

    nullify(cnfg, psys, htrm)
    call base_hamiltonian__iinit__(this, sys, config)
    call json_init(iter, config)
    do
      nullify(cnfg, psys, htrm)
      call json_next(iter, attr, cnfg, ierr)
      if(ierr==JSON_TYPE_ERROR)cycle
      if(ierr/=JSON_OK)exit
      psys => this%sys
      call json_get(cnfg, "system", sysn, ierr)
      if(ierr==JSON_OK) call base_system_gets(sys, sysn, psys)
      ASSERT(associated(psys))
      call base_hamiltonian__inew__(this, attr, htrm)
      call hterm__init__(htrm, psys, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, psys, htrm)

    POP_SUB(base_hamiltonian__init__type)
  end subroutine base_hamiltonian__init__type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: isub, osub
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__init__copy)

    nullify(osub, isub)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_hamiltonian__iinit__(this, that%sys, that%config)
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(osub, isub)
      call hterm_dict_next(iter, name, isub, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call base_hamiltonian__inew__(this, name, osub)
      call hterm__init__(osub, isub)
    end do
    call hterm_dict_end(iter)
    nullify(osub, isub)

    POP_SUB(base_hamiltonian__init__copy)
  end subroutine base_hamiltonian__init__copy

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_type(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config

    PUSH_SUB(base_hamiltonian_init_type)

    call base_hamiltonian__init__(this, sys, config)

    POP_SUB(base_hamiltonian_init_type)
  end subroutine base_hamiltonian_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_init_copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr

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
      call base_hamiltonian_sets(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_hamiltonian_init_copy)
  end subroutine base_hamiltonian_init_copy

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__istart__(this, sim)
    type(base_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_hamiltonian__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call storage_start(this%data, sim)

    POP_SUB(base_hamiltonian__istart__)
  end subroutine base_hamiltonian__istart__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__start__(this, sim)
    type(base_hamiltonian_t),     intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    type(hterm_dict_iterator_t)       :: iter
    type(base_hamiltonian_t), pointer :: prnt
    type(hterm_t),            pointer :: htrm
    integer                           :: ierr

    PUSH_SUB(base_hamiltonian__start__)

    nullify(prnt, htrm)
    if(present(sim))then
      call base_hamiltonian__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        call base_hamiltonian__rget__(this, prnt)
        ASSERT(associated(prnt))
        ASSERT(associated(prnt%sim))
        call base_hamiltonian__istart__(this, prnt%sim)
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

    POP_SUB(base_hamiltonian__start__)
  end subroutine base_hamiltonian__start__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_start(this, sim)
    type(base_hamiltonian_t), intent(inout) :: this
    type(simulation_t),       intent(in)    :: sim

    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr

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
  end subroutine base_hamiltonian_start

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__update__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
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

    POP_SUB(base_hamiltonian__update__)
  end subroutine base_hamiltonian__update__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_update(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr

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
  end subroutine base_hamiltonian_update

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__stop__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
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

    POP_SUB(base_hamiltonian__stop__)
  end subroutine base_hamiltonian__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_stop(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: subs
    integer                           :: ierr

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
  end subroutine base_hamiltonian_stop

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__reset__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(hterm_dict_iterator_t) :: iter
    type(hterm_t),      pointer :: htrm
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__reset__)

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
    this%energy = 0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_hamiltonian__reset__)
  end subroutine base_hamiltonian__reset__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__acc__term(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_term_t),        intent(in)    :: that

    real(kind=wp) :: energy

    PUSH_SUB(base_hamiltonian__acc__term)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_term_get(that, energy=energy)
    this%energy = this%energy + energy

    POP_SUB(base_hamiltonian__acc__term)
  end subroutine base_hamiltonian__acc__term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__acc__potn(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_potential_t),   intent(in)    :: that

    type(storage_t), pointer :: data
    real(kind=wp)            :: energy

    PUSH_SUB(base_hamiltonian__acc__potn)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential_get(that, energy=energy)
    this%energy = this%energy + energy
    call base_potential_get(that, data)
    ASSERT(associated(data))
    call storage_add(this%data, data)

    POP_SUB(base_hamiltonian__acc__potn)
  end subroutine base_hamiltonian__acc__potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__acc__fnct(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_functional_t),  intent(in)    :: that

    type(storage_t), pointer :: data
    real(kind=wp)            :: energy

    PUSH_SUB(base_hamiltonian__acc__fnct)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_functional_get(that, energy=energy)
    this%energy = this%energy + energy
    call base_functional_get(that, data)
    ASSERT(associated(data))
    call storage_add(this%data, data)

    POP_SUB(base_hamiltonian__acc__fnct)
  end subroutine base_hamiltonian__acc__fnct

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__acc__hmlt(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: mhtr, shtr
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__acc__hmlt)

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
    this%energy = this%energy + that%energy
    call storage_add(this%data, that%data)

    POP_SUB(base_hamiltonian__acc__hmlt)
  end subroutine base_hamiltonian__acc__hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__sub__term(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_term_t),        intent(in)    :: that

    real(kind=wp) :: energy

    PUSH_SUB(base_hamiltonian__sub__term)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_term_get(that, energy=energy)
    this%energy = this%energy - energy

    POP_SUB(base_hamiltonian__sub__term)
  end subroutine base_hamiltonian__sub__term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__sub__potn(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_potential_t),   intent(in)    :: that

    type(storage_t), pointer :: data
    real(kind=wp)            :: energy

    PUSH_SUB(base_hamiltonian__sub__potn)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential_get(that, energy=energy)
    this%energy = this%energy - energy
    call base_potential_get(that, data)
    ASSERT(associated(data))
    call storage_sub(this%data, data)

    POP_SUB(base_hamiltonian__sub__potn)
  end subroutine base_hamiltonian__sub__potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__sub__fnct(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_functional_t),  intent(in)    :: that

    type(storage_t), pointer :: data
    real(kind=wp)            :: energy

    PUSH_SUB(base_hamiltonian__sub__fnct)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_functional_get(that, energy=energy)
    this%energy = this%energy - energy
    call base_functional_get(that, data)
    ASSERT(associated(data))
    call storage_sub(this%data, data)

    POP_SUB(base_hamiltonian__sub__fnct)
  end subroutine base_hamiltonian__sub__fnct

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__sub__hmlt(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: mhtr, shtr
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__sub__hmlt)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(mhtr, shtr)
      call hterm_dict_next(iter, name, shtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_dict_get(this%hdct, name, mhtr, ierr)
      if(ierr/=HTERM_DICT_OK)cycle
      call hterm__sub__(mhtr, shtr)
    end do
    call hterm_dict_end(iter)
    nullify(mhtr, shtr)
    this%energy = this%energy - that%energy
    call storage_sub(this%data, that%data)

    POP_SUB(base_hamiltonian__sub__hmlt)
  end subroutine base_hamiltonian__sub__hmlt

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__sets__(this, that, config)
    type(base_hamiltonian_t), intent(inout) :: this
    type(json_object_t),      intent(in)    :: config
    type(base_hamiltonian_t), intent(in)    :: that

    type(hterm_dict_iterator_t)         :: iter
    character(len=CONFIG_DICT_NAME_LEN) :: name
    type(hterm_t),              pointer :: mhtr, shtr
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__sets__)

    call hterm_dict_init(iter, that%hdct)
    do
      nullify(mhtr, shtr)
      call hterm_dict_next(iter, name, shtr, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm_dict_get(this%hdct, name, mhtr, ierr)
      if(ierr/=HTERM_DICT_OK)cycle
      call hterm__sets__(mhtr, shtr, config)
    end do
    call hterm_dict_end(iter)
    nullify(mhtr, shtr)

    POP_SUB(base_hamiltonian__sets__)
  end subroutine base_hamiltonian__sets__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_sets(this, that, config)
    type(base_hamiltonian_t), intent(inout) :: this
    type(json_object_t),      intent(in)    :: config
    type(base_hamiltonian_t), intent(in)    :: that

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_hamiltonian_hash_set(this%hash, config, that)
    call base_hamiltonian__sets__(this, that, config)

    POP_SUB(base_hamiltonian_sets)
  end subroutine base_hamiltonian_sets

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_gets_config(this, config, that)
    type(base_hamiltonian_t),  intent(in) :: this
    type(json_object_t),       intent(in) :: config
    type(base_hamiltonian_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_hamiltonian_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_hamiltonian_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_HAMILTONIAN_OK) nullify(that)

    POP_SUB(base_hamiltonian_gets_config)
  end subroutine base_hamiltonian_gets_config

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_gets_name(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_hamiltonian_gets_name)

    nullify(config, that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_hamiltonian_gets(this, config, that)

    POP_SUB(base_hamiltonian_gets_name)
  end subroutine base_hamiltonian_gets_name

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_info(this, energy)
    type(base_hamiltonian_t), intent(inout) :: this
    real(kind=wp),  optional, intent(in)    :: energy

    PUSH_SUB(base_hamiltonian_set_info)

    if(present(energy)) this%energy = energy

    POP_SUB(base_hamiltonian_set_info)
  end subroutine base_hamiltonian_set_info

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_info(this, size, nspin, energy)
    type(base_hamiltonian_t), intent(in)  :: this
    integer,        optional, intent(out) :: size
    integer,        optional, intent(out) :: nspin
    real(kind=wp),  optional, intent(out) :: energy

    PUSH_SUB(base_hamiltonian_get_info)

    call storage_get(this%data, size=size, dim=nspin)
    if(present(energy)) energy = this%energy

    POP_SUB(base_hamiltonian_get_info)
  end subroutine base_hamiltonian_get_info

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_config(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(json_object_t),     pointer     :: that

    PUSH_SUB(base_hamiltonian_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_hamiltonian_get_config)
  end subroutine base_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_system(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(base_system_t),     pointer     :: that

    PUSH_SUB(base_hamiltonian_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_hamiltonian_get_system)
  end subroutine base_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_simulation(this, that)
    type(base_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that

    PUSH_SUB(base_hamiltonian_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_hamiltonian_get_simulation)
  end subroutine base_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__get__(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(hterm_t),           pointer     :: that

    integer :: ierr

    PUSH_SUB(base_hamiltonian__get__)

    nullify(that)
    call hterm_dict_get(this%hdct, name, that, ierr)
    if(ierr/=HTERM_DICT_OK) nullify(that)

    POP_SUB(base_hamiltonian__get__)
  end subroutine base_hamiltonian__get__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_term(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_term_t),       pointer     :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_term)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, name, htrm)
    if(associated(htrm)) call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_term)
  end subroutine base_hamiltonian_get_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_potn(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_potential_t),  pointer     :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_potn)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, name, htrm)
    if(associated(htrm)) call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_potn)
  end subroutine base_hamiltonian_get_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_fnct(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_functional_t), pointer     :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_fnct)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, name, htrm)
    if(associated(htrm)) call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_fnct)
  end subroutine base_hamiltonian_get_fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that

    type(hterm_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_hmlt)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, name, htrm)
    if(associated(htrm)) call hterm__get__(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_hmlt)
  end subroutine base_hamiltonian_get_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hamiltonian_1d(this, that)
    type(base_hamiltonian_t),     intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_hamiltonian_get_hamiltonian_1d)

    call storage_get(this%data, that)

    POP_SUB(base_hamiltonian_get_hamiltonian_1d)
  end subroutine base_hamiltonian_get_hamiltonian_1d

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hamiltonian_md(this, that)
    type(base_hamiltonian_t),       intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(base_hamiltonian_get_hamiltonian_md)

    call storage_get(this%data, that)

    POP_SUB(base_hamiltonian_get_hamiltonian_md)
  end subroutine base_hamiltonian_get_hamiltonian_md

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__icopy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    PUSH_SUB(base_hamiltonian__icopy__)

    call base_hamiltonian__iend__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_hamiltonian__iinit__(this, that%sys, that%config)
      this%energy = that%energy
      if(associated(that%sim))then
        call base_hamiltonian__istart__(this, that%sim)
        call storage_copy(this%data, that%data)
      end if
    end if

    POP_SUB(base_hamiltonian__icopy__)
  end subroutine base_hamiltonian__icopy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__copy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(hterm_dict_iterator_t)         :: iter
    type(hterm_t),              pointer :: isub, osub
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_hamiltonian__copy__)

    call base_hamiltonian__icopy__(this, that)
    call hterm_dict_init(iter, that%hdct)
    do
      nullify(osub, isub)
      call hterm_dict_next(iter, name, isub, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call base_hamiltonian__inew__(this, name, osub)
      call hterm__copy__(osub, isub)
    end do
    call hterm_dict_end(iter)
    nullify(osub, isub)

    POP_SUB(base_hamiltonian__copy__)
  end subroutine base_hamiltonian__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_copy_type(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(base_hamiltonian_iterator_t) :: iter
    type(base_hamiltonian_t), pointer :: osub, isub
    type(json_object_t),      pointer :: cnfg
    integer                           :: ierr

    PUSH_SUB(base_hamiltonian_copy_type)

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
      call base_hamiltonian_sets(this, osub, cnfg)
    end do
    call base_hamiltonian_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_hamiltonian_copy_type)
  end subroutine base_hamiltonian_copy_type

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iend__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__iend__)

    call base_hamiltonian__inull__(this)
    call storage_end(this%data)
    call hterm_dict_end(this%hdct)
    call config_dict_end(this%dict)
    call base_hamiltonian_hash_end(this%hash)
    call base_hamiltonian__raii_end__(this%raii)

    POP_SUB(base_hamiltonian__iend__)
  end subroutine base_hamiltonian__iend__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__end__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(hterm_t), pointer :: htrm
    integer                :: ierr

    PUSH_SUB(base_hamiltonian__end__)

    do
      nullify(htrm)
      call hterm_dict_pop(this%hdct, htrm, ierr)
      if(ierr/=HTERM_DICT_OK)exit
      call hterm__del__(htrm)
    end do
    nullify(htrm)
    call base_hamiltonian__iend__(this)

    POP_SUB(base_hamiltonian__end__)
  end subroutine base_hamiltonian__end__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_end_type(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_hamiltonian_t), pointer :: subs

    PUSH_SUB(base_hamiltonian_end_type)

    do
      nullify(subs)
      call base_hamiltonian__rpop__(this, subs)
      if(.not.associated(subs))exit
      call base_hamiltonian_end(subs)
      call base_hamiltonian__idel__(subs)
    end do
    nullify(subs)
    call base_hamiltonian__end__(this)

    POP_SUB(base_hamiltonian_end_type)
  end subroutine base_hamiltonian_end_type

#define TEMPLATE_PREFIX base_hamiltonian
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_hamiltonian_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
