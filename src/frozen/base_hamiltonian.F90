#include "global.h"

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_hamiltonian_oct_m

  use base_density_oct_m
  use base_functional_oct_m
  use base_potential_oct_m
  use base_system_oct_m
  use base_term_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use memo_oct_m
  use message_oct_m
  use messages_oct_m
  use msgbus_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m
  use uuid_oct_m

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_husk_t
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TYPE_NAME
#undef DICT_TEMPLATE_NAME

#define BASE_TEMPLATE_NAME base_hamiltonian
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::         &
    TERM_TYPE_TERM, &
    TERM_TYPE_POTN, &
    TERM_TYPE_FNCT, &
    TERM_TYPE_HMLT

  public ::             &
    base_hamiltonian_t

  public ::                     &
    base_hamiltonian__init__,   &
    base_hamiltonian__start__,  &
    base_hamiltonian__update__, &
    base_hamiltonian__reset__,  &
    base_hamiltonian__stop__,   &
    base_hamiltonian__copy__,   &
    base_hamiltonian__end__

  public ::                  &
    base_hamiltonian_new,    &
    base_hamiltonian_del,    &
    base_hamiltonian_init,   &
    base_hamiltonian_start,  &
    base_hamiltonian_calc,   &
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_set,    &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_husk_t
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TYPE_NAME
#undef DICT_TEMPLATE_NAME

#define BASE_TEMPLATE_NAME base_hamiltonian
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  integer, parameter :: TERM_TYPE_NONE = 0
  integer, parameter :: TERM_TYPE_TERM = 1
  integer, parameter :: TERM_TYPE_POTN = 2
  integer, parameter :: TERM_TYPE_FNCT = 3
  integer, parameter :: TERM_TYPE_HMLT = 4

  integer, parameter :: TERM_STAT_DISA = 0
  integer, parameter :: TERM_STAT_NULL = 1
  integer, parameter :: TERM_STAT_ASSC = 2
  
  integer, parameter :: default_nspin = 1

  type :: term_intrf_t
    private
    type(base_term_t),        pointer :: term =>null()
    type(base_potential_t),   pointer :: potn =>null()
    type(base_functional_t),  pointer :: fnct =>null()
    type(base_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: type = TERM_TYPE_NONE
    integer                           :: stat = TERM_STAT_DISA
    type(uuid_t)                      :: uuid
  end type term_intrf_t

  type :: term_husk_t
    private
    type(term_intrf_t), pointer :: self =>null()
    integer                     :: stat = TERM_STAT_DISA
    logical                     :: lock = .false.
    logical                     :: actv = .true.
    real(kind=wp)               :: fctr = 1.0_wp
    type(uuid_t)                :: uuid
  end type term_husk_t

  type :: base_hamiltonian_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(base_system_t),  pointer :: sys    =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(refcount_t),     pointer :: rcnt   =>null()
    integer                       :: nspin  = default_nspin
    type(memo_t)                  :: memo
    type(storage_t)               :: data
    type(term_dict_t)             :: hdct
    type(msgbus_t)                :: msgb
    type(base_hamiltonian_dict_t) :: dict
  end type base_hamiltonian_t

  interface term_intrf_new
    module procedure term_intrf_new_type
    module procedure term_intrf_new_copy
    module procedure term_intrf_new_term
    module procedure term_intrf_new_potn
    module procedure term_intrf_new_fnct
    module procedure term_intrf_new_hmlt
  end interface term_intrf_new

  interface term_intrf_init
    module procedure term_intrf_init_type
    module procedure term_intrf_init_copy
    module procedure term_intrf_init_term
    module procedure term_intrf_init_potn
    module procedure term_intrf_init_fnct
    module procedure term_intrf_init_hmlt
  end interface term_intrf_init

  interface term_intrf_sets
    module procedure term_intrf_sets_info
    module procedure term_intrf_sets_type
  end interface term_intrf_sets

  interface term_intrf_set
    module procedure term_intrf_set_term
    module procedure term_intrf_set_potn
    module procedure term_intrf_set_fnct
    module procedure term_intrf_set_hmlt
  end interface term_intrf_set

  interface term_intrf_get
    module procedure term_intrf_get_term
    module procedure term_intrf_get_potn
    module procedure term_intrf_get_fnct
    module procedure term_intrf_get_hmlt
    module procedure term_intrf_get_energy
    module procedure term_intrf_get_msgbus
    module procedure term_intrf_get_storage
  end interface term_intrf_get

  interface term_husk_set
    module procedure term_husk_set_info
    module procedure term_husk_set_type
  end interface term_husk_set

  interface term_husk_get
    module procedure term_husk_get_info
    module procedure term_husk_get_type
  end interface term_husk_get

  interface base_hamiltonian__init__
    module procedure base_hamiltonian__init__type
    module procedure base_hamiltonian__init__copy
  end interface base_hamiltonian__init__

  interface base_hamiltonian__sets__
    module procedure base_hamiltonian__sets__info
    module procedure base_hamiltonian__sets__type
  end interface base_hamiltonian__sets__

  interface base_hamiltonian__set__
    module procedure base_hamiltonian__set__info
    module procedure base_hamiltonian__set__type
  end interface base_hamiltonian__set__

  interface base_hamiltonian_new
    module procedure base_hamiltonian_new_type
    module procedure base_hamiltonian_new_pass
  end interface base_hamiltonian_new

  interface base_hamiltonian_del
    module procedure base_hamiltonian_del_none
    module procedure base_hamiltonian_del_term
    module procedure base_hamiltonian_del_potn
    module procedure base_hamiltonian_del_fnct
    module procedure base_hamiltonian_del_hmlt
  end interface base_hamiltonian_del

  interface base_hamiltonian_init
    module procedure base_hamiltonian_init_type
  end interface base_hamiltonian_init

  interface base_hamiltonian_notify
    module procedure base_hamiltonian_notify_type
    module procedure base_hamiltonian_notify_subs
  end interface base_hamiltonian_notify

  interface base_hamiltonian_set
    module procedure base_hamiltonian_set_term
    module procedure base_hamiltonian_set_potn
    module procedure base_hamiltonian_set_fnct
    module procedure base_hamiltonian_set_hmlt
  end interface base_hamiltonian_set
  
  interface base_hamiltonian_get
    module procedure base_hamiltonian_get_info
    module procedure base_hamiltonian_get_energy
    module procedure base_hamiltonian_get_config
    module procedure base_hamiltonian_get_system
    module procedure base_hamiltonian_get_density
    module procedure base_hamiltonian_get_simulation
    module procedure base_hamiltonian_get_msgbus
    module procedure base_hamiltonian_get_term
    module procedure base_hamiltonian_get_potn
    module procedure base_hamiltonian_get_fnct
    module procedure base_hamiltonian_get_hmlt
    module procedure base_hamiltonian_get_storage
    module procedure base_hamiltonian_get_data_r1
    module procedure base_hamiltonian_get_data_r2
  end interface base_hamiltonian_get

contains

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_husk_t
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TYPE_NAME
#undef DICT_TEMPLATE_NAME

#define BASE_TEMPLATE_NAME base_hamiltonian
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function term_intrf_new_type(sys, config) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(term_intrf_t), pointer :: this
    
    PUSH_SUB(term_intrf_new_type)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_intrf_init(this, sys, config)
    
    POP_SUB(term_intrf_new_type)
  end function term_intrf_new_type

  ! ---------------------------------------------------------
  recursive function term_intrf_new_copy(source, mold) result(this)
    type(term_intrf_t), optional, intent(in) :: source
    type(term_intrf_t), optional, intent(in) :: mold

    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_new_copy)

    ASSERT(.not.(present(source).eqv.present(mold)))
    nullify(this)
    SAFE_ALLOCATE(this)
    if(present(source))then
      call term_intrf_copy(this, source)
    elseif(present(mold))then
      call term_intrf_init(this, mold)
    else
      ASSERT(.FALSE.)
    end if

    POP_SUB(term_intrf_new_copy)
  end function term_intrf_new_copy

  ! ---------------------------------------------------------
  function term_intrf_new_term(that) result(this)
    type(base_term_t), intent(in) :: that

    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_new_term)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_intrf_init(this, that)

    POP_SUB(term_intrf_new_term)
  end function term_intrf_new_term

  ! ---------------------------------------------------------
  function term_intrf_new_potn(that) result(this)
    type(base_potential_t), intent(in) :: that

    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_new_potn)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_intrf_init(this, that)

    POP_SUB(term_intrf_new_potn)
  end function term_intrf_new_potn

  ! ---------------------------------------------------------
  function term_intrf_new_fnct(that) result(this)
    type(base_functional_t), intent(in) :: that

    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_new_fnct)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_intrf_init(this, that)

    POP_SUB(term_intrf_new_fnct)
  end function term_intrf_new_fnct

  ! ---------------------------------------------------------
  function term_intrf_new_hmlt(that) result(this)
    type(base_hamiltonian_t), intent(in) :: that

    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_new_hmlt)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_intrf_init(this, that)

    POP_SUB(term_intrf_new_hmlt)
  end function term_intrf_new_hmlt

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_del(this)
    type(term_intrf_t), pointer, intent(inout) :: this

    PUSH_SUB(term_intrf_del)

    if(associated(this))then
      ASSERT(term_intrf_assoc(this))
      call term_intrf_end(this)
      SAFE_DEALLOCATE_P(this)
      nullify(this)
    end if

    POP_SUB(term_intrf_del)
  end subroutine term_intrf_del
  
  ! ---------------------------------------------------------
  function term_intrf_assoc(this) result(that)
    type(term_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(term_intrf_assoc)

    that = .false.
    select case(this%stat)
    case(TERM_STAT_DISA, TERM_STAT_NULL)
      select case(this%type)
      case(TERM_TYPE_NONE)
        ASSERT(.not.associated(this%term))
        ASSERT(.not.associated(this%potn))
        ASSERT(.not.associated(this%fnct))
        ASSERT(.not.associated(this%hmlt))
        that = .false.
      case default
        ASSERT(.false.)
      end select
    case(TERM_STAT_ASSC)
      select case(this%type)
      case(TERM_TYPE_TERM)
        ASSERT(associated(this%term))
        ASSERT(.not.associated(this%potn))
        ASSERT(.not.associated(this%fnct))
        ASSERT(.not.associated(this%hmlt))
        that = .true.
      case(TERM_TYPE_POTN)
        ASSERT(.not.associated(this%term))
        ASSERT(associated(this%potn))
        ASSERT(.not.associated(this%fnct))
        ASSERT(.not.associated(this%hmlt))
        that = .true.
      case(TERM_TYPE_FNCT)
        ASSERT(.not.associated(this%term))
        ASSERT(.not.associated(this%potn))
        ASSERT(associated(this%fnct))
        ASSERT(.not.associated(this%hmlt))
        that = .true.
      case(TERM_TYPE_HMLT)
        ASSERT(.not.associated(this%term))
        ASSERT(.not.associated(this%potn))
        ASSERT(.not.associated(this%fnct))
        ASSERT(associated(this%hmlt))
        that = .true.
      case default
        ASSERT(.false.)
      end select
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_assoc)
  end function term_intrf_assoc

  ! ---------------------------------------------------------
  subroutine term_intrf_init_type(this, sys, config)
    type(term_intrf_t),  intent(out) :: this
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    type(base_term_t),        pointer :: term
    type(base_potential_t),   pointer :: potn
    type(base_functional_t),  pointer :: fnct
    type(base_hamiltonian_t), pointer :: hmlt
    integer                           :: type, ierr

    PUSH_SUB(term_intrf_init_type)

    nullify(term, potn, fnct, hmlt)
    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    call json_get(config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    select case(type)
    case(TERM_TYPE_TERM)
      term => base_term_new(sys, config)
      call term_intrf_set(this, term)
      nullify(term)
    case(TERM_TYPE_POTN)
      potn => base_potential_new(sys, config)
      call term_intrf_set(this, potn)
      nullify(potn)
    case(TERM_TYPE_FNCT)
      fnct => base_functional_new(sys, config)
      call term_intrf_set(this, fnct)
      nullify(fnct)
    case(TERM_TYPE_HMLT)
      hmlt => base_hamiltonian_new(sys, config)
      call term_intrf_set(this, hmlt)
      nullify(hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_init_type)
  end subroutine term_intrf_init_type

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_init_copy(this, that)
    type(term_intrf_t), intent(out) :: this
    type(term_intrf_t), intent(in)  :: that

    type(base_term_t),        pointer :: term
    type(base_potential_t),   pointer :: potn
    type(base_functional_t),  pointer :: fnct
    type(base_hamiltonian_t), pointer :: hmlt

    PUSH_SUB(term_intrf_init_copy)

    ASSERT(term_intrf_assoc(that))
    nullify(term, potn, fnct, hmlt)
    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    select case(that%type)
    case(TERM_TYPE_TERM)
      call term_intrf_get(that, term)
      ASSERT(associated(term))
      call term_intrf_set(this, base_term_new(mold=term))
      nullify(term)
    case(TERM_TYPE_POTN)
      call term_intrf_get(that, potn)
      ASSERT(associated(potn))
      call term_intrf_set(this, base_potential_new(mold=potn))
      nullify(potn)
    case(TERM_TYPE_FNCT)
      call term_intrf_get(that, fnct)
      ASSERT(associated(fnct))
      call term_intrf_set(this,  base_functional_new(mold=fnct))
      nullify(fnct)
    case(TERM_TYPE_HMLT)
      call term_intrf_get(that, hmlt)
      ASSERT(associated(hmlt))
      call term_intrf_set(this,  base_hamiltonian_new(mold=hmlt))
      nullify(hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_init_copy)
  end subroutine term_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine term_intrf_init_term(this, that)
    type(term_intrf_t), intent(out) :: this
    type(base_term_t),  intent(in)  :: that

    PUSH_SUB(term_intrf_init_term)

    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_term)
  end subroutine term_intrf_init_term

  ! ---------------------------------------------------------
  subroutine term_intrf_init_potn(this, that)
    type(term_intrf_t),     intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_potn)

    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_potn)
  end subroutine term_intrf_init_potn

  ! ---------------------------------------------------------
  subroutine term_intrf_init_fnct(this, that)
    type(term_intrf_t),      intent(out) :: this
    type(base_functional_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_fnct)

    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_fnct)
  end subroutine term_intrf_init_fnct

  ! ---------------------------------------------------------
  subroutine term_intrf_init_hmlt(this, that)
    type(term_intrf_t),       intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_hmlt)

    this%type = TERM_TYPE_NONE
    this%stat = TERM_STAT_NULL
    call uuid_init(this%uuid)
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_hmlt)
  end subroutine term_intrf_init_hmlt

  ! ---------------------------------------------------------
  subroutine term_intrf__attach__(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(uuid_t),       intent(in)    :: that

    PUSH_SUB(term_intrf__attach__)

    ASSERT(term_intrf_assoc(this))  
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__attach__(this%term, that)
    case(TERM_TYPE_POTN)
      call base_potential__attach__(this%potn, that)
    case(TERM_TYPE_FNCT)
      call base_functional__attach__(this%fnct, that)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__attach__(this%hmlt, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__attach__)
  end subroutine term_intrf__attach__

  ! ---------------------------------------------------------
  subroutine term_intrf__detach__(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(uuid_t),       intent(in)    :: that

    PUSH_SUB(term_intrf__detach__)

    ASSERT(term_intrf_assoc(this))  
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__detach__(this%term, that)
    case(TERM_TYPE_POTN)
      call base_potential__detach__(this%potn, that)
    case(TERM_TYPE_FNCT)
      call base_functional__detach__(this%fnct, that)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__detach__(this%hmlt, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__detach__)
  end subroutine term_intrf__detach__

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__start__(this, sim)
    type(term_intrf_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: sim

    PUSH_SUB(term_intrf__start__)

    ASSERT(term_intrf_assoc(this))  
    select case(this%type)
    case(TERM_TYPE_TERM)
    case(TERM_TYPE_POTN)
      call base_potential__start__(this%potn, sim)
    case(TERM_TYPE_FNCT)
      call base_functional__start__(this%fnct, sim)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__start__(this%hmlt, sim)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__start__)
  end subroutine term_intrf__start__

  ! ---------------------------------------------------------
  recursive function term_intrf_calc(this, that) result(energy)
    type(term_intrf_t),             intent(inout) :: this
    type(base_density_t), optional, intent(in)    :: that

    real(kind=wp) :: energy

    PUSH_SUB(term_intrf_calc)

    ASSERT(term_intrf_assoc(this))
    energy = 0.0_wp
    select case(this%type)
    case(TERM_TYPE_TERM)
    case(TERM_TYPE_POTN)
      energy = base_potential_calc(this%potn, that)
    case(TERM_TYPE_FNCT)
      energy = base_functional_calc(this%fnct, that)
    case(TERM_TYPE_HMLT)
      energy = base_hamiltonian_calc(this%hmlt, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_calc)
  end function term_intrf_calc

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__update__(this, energy)
    type(term_intrf_t), intent(inout) :: this
    logical,  optional, intent(in)    :: energy

    PUSH_SUB(term_intrf__update__)

    ASSERT(term_intrf_assoc(this))
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__update__(this%term)
    case(TERM_TYPE_POTN)
      call base_potential__update__(this%potn, energy)
    case(TERM_TYPE_FNCT)
      call base_functional__update__(this%fnct, energy)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__update__(this%hmlt, energy)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__update__)
  end subroutine term_intrf__update__

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__reset__(this)
    type(term_intrf_t), intent(inout) :: this

    PUSH_SUB(term_intrf__reset__)

    ASSERT(term_intrf_assoc(this))
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__reset__(this%term)
    case(TERM_TYPE_POTN)
      call base_potential__reset__(this%potn)
    case(TERM_TYPE_FNCT)
      call base_functional__reset__(this%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__reset__(this%hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__reset__)
  end subroutine term_intrf__reset__

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__stop__(this)
    type(term_intrf_t), intent(inout) :: this

    PUSH_SUB(term_intrf__stop__)

    ASSERT(term_intrf_assoc(this))
    select case(this%type)
    case(TERM_TYPE_TERM)
    case(TERM_TYPE_POTN)
      call base_potential__stop__(this%potn)
    case(TERM_TYPE_FNCT)
      call base_functional__stop__(this%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__stop__(this%hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__stop__)
  end subroutine term_intrf__stop__

  ! ---------------------------------------------------------
  recursive function term_intrf_exists(this, name) result(that)
    type(term_intrf_t), intent(inout) :: this
    character(len=*),   intent(in)    :: name
    
    logical :: that

    PUSH_SUB(term_intrf_exists)

    ASSERT(term_intrf_assoc(this))
    that = .false.
    select case(this%type)
    case(TERM_TYPE_TERM)
      that = base_term_exists(this%term, trim(adjustl(name)))
    case(TERM_TYPE_POTN)
      that = base_potential_exists(this%potn, trim(adjustl(name)))
    case(TERM_TYPE_FNCT)
      that = base_functional_exists(this%fnct, trim(adjustl(name)))
    case(TERM_TYPE_HMLT)
      that = base_hamiltonian_exists(this%hmlt, trim(adjustl(name)))
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_exists)
  end function term_intrf_exists

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_sets_info(this, name, lock, active)
    type(term_intrf_t), intent(inout) :: this
    character(len=*),   intent(in)    :: name
    logical,  optional, intent(in)    :: lock
    logical,  optional, intent(in)    :: active

    PUSH_SUB(term_intrf_sets_info)

    ASSERT(term_intrf_assoc(this))
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_sets(this%term, trim(adjustl(name)), lock=lock, active=active)
    case(TERM_TYPE_POTN)
      call base_potential_sets(this%potn, trim(adjustl(name)), lock=lock, active=active)
    case(TERM_TYPE_FNCT)
      call base_functional_sets(this%fnct, trim(adjustl(name)), lock=lock, active=active)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_sets(this%hmlt, trim(adjustl(name)), lock=lock, active=active)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_sets_info)
  end subroutine term_intrf_sets_info

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_sets_type(this, name, that, config, lock, active)
    type(term_intrf_t),  intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(term_intrf_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    logical,   optional, intent(in)    :: lock
    logical,   optional, intent(in)    :: active

    PUSH_SUB(term_intrf_sets_type)

    ASSERT(term_intrf_assoc(this))
    ASSERT(term_intrf_assoc(that))
    ASSERT(this%type==that%type)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_sets(this%term, trim(adjustl(name)), that%term, config, lock=lock, active=active)
    case(TERM_TYPE_POTN)
      call base_potential_sets(this%potn, trim(adjustl(name)), that%potn, config, lock=lock, active=active)
    case(TERM_TYPE_FNCT)
      call base_functional_sets(this%fnct, trim(adjustl(name)), that%fnct, config, lock=lock, active=active)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_sets(this%hmlt, trim(adjustl(name)), that%hmlt, config, lock=lock, active=active)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_sets_type)
  end subroutine term_intrf_sets_type

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_dels(this, name, that)
    type(term_intrf_t), intent(inout) :: this
    character(len=*),   intent(in)    :: name
    type(term_intrf_t), intent(in)    :: that

    PUSH_SUB(term_intrf_dels)

    ASSERT(term_intrf_assoc(this))
    ASSERT(term_intrf_assoc(that))
    ASSERT(this%type==that%type)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_dels(this%term, trim(adjustl(name)), that%term)
    case(TERM_TYPE_POTN)
      call base_potential_dels(this%potn, trim(adjustl(name)), that%potn)
    case(TERM_TYPE_FNCT)
      call base_functional_dels(this%fnct, trim(adjustl(name)), that%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_dels(this%hmlt, trim(adjustl(name)), that%hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_dels)
  end subroutine term_intrf_dels
  
  ! ---------------------------------------------------------
  subroutine term_intrf_set_term(this, that)
    type(term_intrf_t),        intent(inout) :: this
    type(base_term_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_term)

    ASSERT(this%type==TERM_TYPE_NONE)
    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(.not.term_intrf_assoc(this))
    this%term => that
    call base_term__attach__(this%term, this%uuid)
    this%type = TERM_TYPE_TERM
    this%stat = TERM_STAT_ASSC

    POP_SUB(term_intrf_set_term)
  end subroutine term_intrf_set_term

  ! ---------------------------------------------------------
  subroutine term_intrf_set_potn(this, that)
    type(term_intrf_t),             intent(inout) :: this
    type(base_potential_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_potn)

    ASSERT(this%type==TERM_TYPE_NONE)
    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(.not.term_intrf_assoc(this))
    this%potn => that
    call base_potential__attach__(this%potn, this%uuid)
    this%type = TERM_TYPE_POTN
    this%stat = TERM_STAT_ASSC

    POP_SUB(term_intrf_set_potn)
  end subroutine term_intrf_set_potn

  ! ---------------------------------------------------------
  subroutine term_intrf_set_fnct(this, that)
    type(term_intrf_t),              intent(inout) :: this
    type(base_functional_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_fnct)

    ASSERT(this%type==TERM_TYPE_NONE)
    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(.not.term_intrf_assoc(this))
    this%fnct => that
    call base_functional__attach__(this%fnct, this%uuid)
    this%type = TERM_TYPE_FNCT
    this%stat = TERM_STAT_ASSC

    POP_SUB(term_intrf_set_fnct)
  end subroutine term_intrf_set_fnct

  ! ---------------------------------------------------------
  subroutine term_intrf_set_hmlt(this, that)
    type(term_intrf_t),               intent(inout) :: this
    type(base_hamiltonian_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_hmlt)

    ASSERT(this%type==TERM_TYPE_NONE)
    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(.not.term_intrf_assoc(this))
    this%hmlt => that
    call base_hamiltonian__attach__(this%hmlt, this%uuid)
    this%type = TERM_TYPE_HMLT
    this%stat = TERM_STAT_ASSC

    POP_SUB(term_intrf_set_hmlt)
  end subroutine term_intrf_set_hmlt

  ! ---------------------------------------------------------
  subroutine term_intrf_get_term(this, that)
    type(term_intrf_t),         intent(in)  :: this
    type(base_term_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_term)

    nullify(that)
    if(term_intrf_assoc(this))then
      ASSERT(this%type==TERM_TYPE_TERM)
      ASSERT(associated(this%term))
      that => this%term
    end if

    POP_SUB(term_intrf_get_term)
  end subroutine term_intrf_get_term

  ! ---------------------------------------------------------
  subroutine term_intrf_get_potn(this, that)
    type(term_intrf_t),              intent(in)  :: this
    type(base_potential_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_potn)

    nullify(that)
    if(term_intrf_assoc(this))then
      ASSERT(this%type==TERM_TYPE_POTN)
      ASSERT(associated(this%potn))
      that => this%potn
    end if

    POP_SUB(term_intrf_get_potn)
  end subroutine term_intrf_get_potn

  ! ---------------------------------------------------------
  subroutine term_intrf_get_fnct(this, that)
    type(term_intrf_t),               intent(in)  :: this
    type(base_functional_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_fnct)

    nullify(that)
    if(term_intrf_assoc(this))then
      ASSERT(this%type==TERM_TYPE_FNCT)
      ASSERT(associated(this%fnct))
      that => this%fnct
    end if

    POP_SUB(term_intrf_get_fnct)
  end subroutine term_intrf_get_fnct

  ! ---------------------------------------------------------
  subroutine term_intrf_get_hmlt(this, that)
    type(term_intrf_t),                intent(in)  :: this
    type(base_hamiltonian_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_hmlt)

    nullify(that)
    if(term_intrf_assoc(this))then
      ASSERT(this%type==TERM_TYPE_HMLT)
      ASSERT(associated(this%hmlt))
      that => this%hmlt
    end if

    POP_SUB(term_intrf_get_hmlt)
  end subroutine term_intrf_get_hmlt

  ! ---------------------------------------------------------
  subroutine term_intrf_get_energy(this, energy)
    type(term_intrf_t), intent(in)  :: this
    real(kind=wp),      intent(out) :: energy

    PUSH_SUB(term_intrf_get_energy)

    ASSERT(term_intrf_assoc(this))
    energy = 0.0_wp
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_get(this%term, energy)
    case(TERM_TYPE_POTN)
      call base_potential_get(this%potn, energy)
    case(TERM_TYPE_FNCT)
      call base_functional_get(this%fnct, energy)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_get(this%hmlt, energy)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_get_energy)
  end subroutine term_intrf_get_energy

  ! ---------------------------------------------------------
  subroutine term_intrf_get_msgbus(this, that)
    type(term_intrf_t),      intent(in)  :: this
    type(msgbus_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_msgbus)

    ASSERT(term_intrf_assoc(this))
    nullify(that)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_get(this%term, that)
    case(TERM_TYPE_POTN)
      call base_potential_get(this%potn, that)
    case(TERM_TYPE_FNCT)
      call base_functional_get(this%fnct, that)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_get(this%hmlt, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_get_msgbus)
  end subroutine term_intrf_get_msgbus
  
  ! ---------------------------------------------------------
  subroutine term_intrf_get_storage(this, that)
    type(term_intrf_t),       intent(in)  :: this
    type(storage_t), pointer, intent(out) :: that

    PUSH_SUB(term_intrf_get_storage)

    ASSERT(term_intrf_assoc(this))
    nullify(that)
    select case(this%type)
    case(TERM_TYPE_TERM)
    case(TERM_TYPE_POTN)
      call base_potential_get(this%potn, that)
    case(TERM_TYPE_FNCT)
      call base_functional_get(this%fnct, that)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_get(this%hmlt, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_get_storage)
  end subroutine term_intrf_get_storage

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_copy(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(term_intrf_t), intent(in)    :: that

    type(base_term_t),        pointer :: term
    type(base_potential_t),   pointer :: potn
    type(base_functional_t),  pointer :: fnct
    type(base_hamiltonian_t), pointer :: hmlt

    PUSH_SUB(term_intrf_copy)

    nullify(term, potn, fnct, hmlt)
    call term_intrf_end(this)
    call uuid_init(this%uuid)
    select case(that%stat)
    case(TERM_STAT_DISA)
    case(TERM_STAT_ASSC)
      select case(that%type)
      case(TERM_TYPE_TERM)
        call term_intrf_get(that, term)
        ASSERT(associated(term))
        call term_intrf_set(this, base_term_new(source=term))
        nullify(term)
      case(TERM_TYPE_POTN)
        call term_intrf_get(that, potn)
        ASSERT(associated(potn))
        call term_intrf_set(this, base_potential_new(source=potn))
        nullify(potn)
      case(TERM_TYPE_FNCT)
        call term_intrf_get(that, fnct)
        ASSERT(associated(fnct))
        call term_intrf_set(this, base_functional_new(source=fnct))
        nullify(fnct)
      case(TERM_TYPE_HMLT)
        call term_intrf_get(that, hmlt)
        ASSERT(associated(hmlt))
        call term_intrf_set(this, base_hamiltonian_new(source=hmlt))
        nullify(hmlt)
      case default
        ASSERT(.false.)
      end select
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_copy)
  end subroutine term_intrf_copy

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_end(this)
    type(term_intrf_t), intent(inout) :: this

    PUSH_SUB(term_intrf_end)

    select case(this%stat)
    case(TERM_STAT_DISA)
    case(TERM_STAT_ASSC)
      ASSERT(term_intrf_assoc(this))
      select case(this%type)
      case(TERM_TYPE_TERM)
        call base_term__detach__(this%term, this%uuid)
        call base_term_del(this%term)
      case(TERM_TYPE_POTN)
        call base_potential__detach__(this%potn, this%uuid)
        call base_potential_del(this%potn)
      case(TERM_TYPE_FNCT)
        call base_functional__detach__(this%fnct, this%uuid)
        call base_functional_del(this%fnct)
      case(TERM_TYPE_HMLT)
        call base_hamiltonian__detach__(this%hmlt, this%uuid)
        call base_hamiltonian_del(this%hmlt)
      case default
        ASSERT(.false.)
      end select
    case default
      ASSERT(.false.)
    end select
    nullify(this%term, this%potn, this%fnct, this%hmlt)
    this%type = TERM_STAT_DISA
    call uuid_end(this%uuid)

    POP_SUB(term_intrf_end)
  end subroutine term_intrf_end
  
  ! ---------------------------------------------------------
  function term_husk_new(that, factor, lock, active) result(this)
    type(term_intrf_t),      intent(in) :: that
    real(kind=wp), optional, intent(in) :: factor
    logical,       optional, intent(in) :: lock
    logical,       optional, intent(in) :: active

    type(term_husk_t), pointer :: this
    
    PUSH_SUB(term_husk_new)

    nullify(this)
    SAFE_ALLOCATE(this)
    call term_husk_init(this, that, factor=factor, lock=lock, active=active)

    POP_SUB(term_husk_new)
  end function term_husk_new

  ! ---------------------------------------------------------
  subroutine term_husk_del(this)
    type(term_husk_t), pointer, intent(inout) :: this

    PUSH_SUB(term_husk_del)

    if(associated(this))then
      call term_husk_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(term_husk_del)
  end subroutine term_husk_del

  ! ---------------------------------------------------------
  function term_husk_assoc(this) result(that)
    type(term_husk_t), intent(in) :: this

    logical :: that

    PUSH_SUB(term_husk_assoc)

    select case(this%stat)
    case(TERM_STAT_DISA)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(TERM_STAT_ASSC)
      ASSERT(associated(this%self))
      ASSERT(term_intrf_assoc(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_husk_assoc)
  end function term_husk_assoc

  ! ---------------------------------------------------------
  subroutine term_husk_init(this, that, factor, lock, active)
    type(term_husk_t),       intent(out) :: this
    type(term_intrf_t),      intent(in)  :: that
    real(kind=wp), optional, intent(in)  :: factor
    logical,       optional, intent(in)  :: lock
    logical,       optional, intent(in)  :: active

    PUSH_SUB(term_husk_init)

    nullify(this%self)
    this%stat = TERM_STAT_NULL
    this%lock = .false.
    this%actv = .true.
    this%fctr = 1.0_wp
    call uuid_init(this%uuid)
    call term_husk_set(this, that, factor=factor, lock=lock, active=active)

    POP_SUB(term_husk_init)
  end subroutine term_husk_init

  ! ---------------------------------------------------------
  subroutine term_husk_set_info(this, factor, lock, active)
    type(term_husk_t),       intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: factor
    logical,       optional, intent(in)    :: lock
    logical,       optional, intent(in)    :: active

    PUSH_SUB(term_husk_set_info)

    ASSERT(term_husk_assoc(this))
    if(present(factor)) this%fctr = factor
    if(present(lock))   this%lock = lock
    if(present(active)) this%actv = active

    POP_SUB(term_husk_set_info)
  end subroutine term_husk_set_info

  ! ---------------------------------------------------------
  subroutine term_husk_set_type(this, that, factor, lock, active)
    type(term_husk_t),          intent(inout) :: this
    type(term_intrf_t), target, intent(in)    :: that
    real(kind=wp),    optional, intent(in)    :: factor
    logical,          optional, intent(in)    :: lock
    logical,          optional, intent(in)    :: active

    PUSH_SUB(term_husk_set_type)

    ASSERT(this%stat==TERM_STAT_NULL)
    this%self => that
    call term_intrf__attach__(this%self, this%uuid)
    this%stat = TERM_STAT_ASSC
    if(present(factor).or.present(lock).or.present(active))&
      call term_husk_set(this, factor=factor, lock=lock, active=active)

    POP_SUB(term_husk_set_type)
  end subroutine term_husk_set_type

  ! ---------------------------------------------------------
  subroutine term_husk_get_info(this, factor, lock, active)
    type(term_husk_t),       intent(in)  :: this
    real(kind=wp), optional, intent(out) :: factor
    logical,       optional, intent(out) :: lock
    logical,       optional, intent(out) :: active

    PUSH_SUB(term_husk_get_info)

    ASSERT(term_husk_assoc(this))
    if(present(factor)) factor = this%fctr
    if(present(lock))   lock   = this%lock
    if(present(active)) active = this%actv

    POP_SUB(term_husk_get_info)
  end subroutine term_husk_get_info

  ! ---------------------------------------------------------
  subroutine term_husk_get_type(this, that, factor, lock, active)
    type(term_husk_t),           intent(in)  :: this
    type(term_intrf_t), pointer, intent(out) :: that
    real(kind=wp),     optional, intent(out) :: factor
    logical,           optional, intent(out) :: lock
    logical,           optional, intent(out) :: active

    PUSH_SUB(term_husk_get_type)

    nullify(that)
    if(term_husk_assoc(this)) that => this%self
    if(present(factor).or.present(lock).or.present(active))&
      call term_husk_get(this, factor=factor, lock=lock, active=active)

    POP_SUB(term_husk_get_type)
  end subroutine term_husk_get_type

  ! ---------------------------------------------------------
  subroutine term_husk_copy(this, that)
    type(term_husk_t), intent(inout) :: this
    type(term_husk_t), intent(in)    :: that

    PUSH_SUB(term_husk_copy)

    call term_husk_end(this)
    call uuid_init(this%uuid)
    select case(that%stat)
    case(TERM_STAT_DISA)
    case(TERM_STAT_ASSC)
      ASSERT(term_husk_assoc(that))
      this%stat = TERM_STAT_NULL
      call term_husk_set(this, that%self, factor=that%fctr, lock=that%lock, active=that%actv)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_husk_copy)
  end subroutine term_husk_copy

  ! ---------------------------------------------------------
  subroutine term_husk_end(this)
    type(term_husk_t), intent(inout) :: this

    PUSH_SUB(term_husk_end)

    if(term_husk_assoc(this)) call term_intrf__detach__(this%self, this%uuid)
    nullify(this%self)
    this%stat = TERM_STAT_DISA
    this%lock = .false.
    this%actv = .true.
    this%fctr = 1.0_wp
    call uuid_end(this%uuid)
    
    POP_SUB(term_husk_end)
  end subroutine term_husk_end

  ! ---------------------------------------------------------
  recursive subroutine term__apply__(this, operation, enforce_lock, enforce_active)
    type(base_hamiltonian_t), intent(inout) :: this
    logical,        optional, intent(in)    :: enforce_lock
    logical,        optional, intent(in)    :: enforce_active

    interface
      subroutine operation(this)
        import :: term_intrf_t
        type(term_intrf_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(term_dict_iterator_t)  :: iter
    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    logical                     :: lock, actv, elck, eact
    integer                     :: ierr

    PUSH_SUB(term__apply__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    eact = .false.
    if(present(enforce_active)) eact = enforce_active
    elck = .not.eact
    if(present(enforce_lock)) elck = enforce_lock
    call term_dict_init(iter, this%hdct)
    do
      nullify(husk, term)
      call term_dict_next(iter, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, term, lock=lock, active=actv)
      ASSERT(associated(term))
      if((.not.lock.and.elck).or.(actv.and.eact)) call operation(term)
    end do
    call term_dict_end(iter)
    nullify(husk, term)

    POP_SUB(term__apply__)
  end subroutine term__apply__

  ! ---------------------------------------------------------
  function base_hamiltonian_new_type(sys, config) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_hamiltonian_t), pointer :: this

    PUSH_SUB(base_hamiltonian_new_type)

    this => base_hamiltonian_new(sys, config, base_hamiltonian_init_type)
    
    POP_SUB(base_hamiltonian_new_type)
  end function base_hamiltonian_new_type

  ! ---------------------------------------------------------
  function base_hamiltonian_new_pass(sys, config, init) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_hamiltonian_t), pointer :: this
    
    interface
      subroutine init(this, sys, config)
        use json_oct_m
        use base_system_oct_m
        import :: base_hamiltonian_t
        type(base_hamiltonian_t), intent(out) :: this
        type(base_system_t),      intent(in)  :: sys
        type(json_object_t),      intent(in)  :: config
      end subroutine init
    end interface
    
    PUSH_SUB(base_hamiltonian_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, sys, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_hamiltonian_new_pass)
  end function base_hamiltonian_new_pass

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iinit__(this, sys, config)
    type(base_hamiltonian_t),    intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr

    PUSH_SUB(base_hamiltonian__iinit__)

    nullify(cnfg)
    this%config => config
    this%sys => sys
    this%rcnt => refcount_new()
    call json_get(this%config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(type==TERM_TYPE_HMLT)
    call memo_init(this%memo)
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "full", .false.)
    call storage_init(this%data, cnfg)
    nullify(cnfg)
    call term_dict_init(this%hdct)
    call msgbus_init(this%msgb, number=2)
    call base_hamiltonian_dict_init(this%dict)

    POP_SUB(base_hamiltonian__iinit__)
  end subroutine base_hamiltonian__iinit__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__type(this, sys, config)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_system_t),      intent(in)  :: sys
    type(json_object_t),      intent(in)  :: config

    type(json_object_iterator_t)             :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    type(json_object_t),             pointer :: trms, cnfg, defn
    real(kind=wp)                            :: fctr
    logical                                  :: lock, actv
    integer                                  :: ierr

    PUSH_SUB(base_hamiltonian__init__type)

    nullify(trms, cnfg, defn)
    call base_hamiltonian__iinit__(this, sys, config)
    call json_get(this%config, "terms", trms, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(associated(trms))
    call json_init(iter, trms)
    do
      nullify(cnfg, defn)
      call json_next(iter, name, cnfg, ierr)
      if(ierr/=JSON_OK) exit
      ASSERT(associated(cnfg))
      call json_get(cnfg, "factor", fctr, ierr)
      if(ierr/=JSON_OK) fctr = 1.0_wp
      call json_get(cnfg, "lock", lock, ierr)
      if(ierr/=JSON_OK) lock = .false.
      call json_get(cnfg, "active", actv, ierr)
      if(ierr/=JSON_OK) actv = .true.
      call json_get(cnfg, "definition", defn, ierr)
      ASSERT(ierr==JSON_OK)
      ASSERT(associated(defn))
      call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(sys, defn), factor=fctr, lock=lock, active=actv)
    end do
    call json_end(iter)
    nullify(trms, cnfg, defn)

    POP_SUB(base_hamiltonian__init__type)
  end subroutine base_hamiltonian__init__type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    type(term_husk_t),               pointer :: husk
    type(term_intrf_t),              pointer :: term
    real(kind=wp)                            :: fctr
    logical                                  :: lock, actv
    integer                                  :: ierr

    PUSH_SUB(base_hamiltonian__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_hamiltonian__iinit__(this, that%sys, that%config)
    call term_dict_init(iter, that%hdct)
    do
      nullify(husk, term)
      call term_dict_next(iter, name, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, term, factor=fctr, lock=lock, active=actv)
      ASSERT(associated(term))
      call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(mold=term), factor=fctr, lock=lock, active=actv)
    end do
    call term_dict_end(iter)
    nullify(husk, term)
    if(associated(that%sim)) call base_hamiltonian__start__(this, that%sim)

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
  recursive subroutine base_hamiltonian__start__(this, sim)
    type(base_hamiltonian_t),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    logical :: uspn
    integer :: ierr

    PUSH_SUB(base_hamiltonian__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    call term__apply__(this, start)
    this%sim => sim
    this%nspin = default_nspin
    call json_get(this%config, "spin", uspn, ierr)
    if(ierr/=JSON_OK) uspn = .true.
    if(uspn) call base_system_get(this%sys, nspin=this%nspin)
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    call storage_start(this%data, sim, ndim=this%nspin)

    POP_SUB(base_hamiltonian__start__)

  contains

    recursive subroutine start(this)
      type(term_intrf_t), intent(inout) :: this

      PUSH_SUB(base_hamiltonian__start__.start)

      call term_intrf__start__(this, sim)
      
      POP_SUB(base_hamiltonian__start__.start)
    end subroutine start
  
  end subroutine base_hamiltonian__start__

  ! ---------------------------------------------------------
  function base_hamiltonian__calc__(this, that) result(energy)
    type(base_hamiltonian_t), intent(in) :: this
    type(base_density_t),     intent(in) :: that

    real(kind=wp) :: energy

    type(term_dict_iterator_t)  :: iter
    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    type(storage_t),    pointer :: data
    real(kind=wp)               :: fctr
    logical                     :: fuse, actv
    integer                     :: nspn, ierr

    PUSH_SUB(base_hamiltonian__calc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(husk, term, data)
    energy = 0.0_wp
    call base_hamiltonian_get(this, use=fuse, nspin=nspn)
    ASSERT(this%nspin<=nspn)
    if(fuse)then
      call base_density_get(that, data, total=(this%nspin<nspn))
      ASSERT(associated(data))
      call storage_integrate(this%data, data, energy)
      nullify(data)
    else
      call term_dict_init(iter, this%hdct)
      do
        nullify(husk, term)
        call term_dict_next(iter, husk, ierr)
        if(ierr/=TERM_DICT_OK)exit
        ASSERT(associated(husk))
        ASSERT(term_husk_assoc(husk))
        call term_husk_get(husk, term, factor=fctr, active=actv)
        ASSERT(associated(term))
        if(actv.and.(abs(fctr)>tiny(fctr)))&
          energy = energy + fctr * term_intrf_calc(term, that)
      end do
      call term_dict_end(iter)
      nullify(husk, term)
    end if

    POP_SUB(base_hamiltonian__calc__)
  end function base_hamiltonian__calc__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__update__(this, energy)
    type(base_hamiltonian_t), intent(inout) :: this
    logical,        optional, intent(in)    :: energy

    type(term_dict_iterator_t)  :: iter
    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    real(kind=wp)               :: fctr, enrg
    logical                     :: fuse, updt, calc, actv
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    calc = .false.
    if(present(energy)) calc = energy
    call term__apply__(this, update, enforce_active=.true.)
    call base_hamiltonian__recv__(this, updt, channel=2)
    if(updt.or.calc)then
      if(updt) call base_hamiltonian__reset__(this)
      call base_hamiltonian_get(this, use=fuse)
      updt = (fuse.and.updt)
      calc = (.not.memo_in(this%memo, "energy").and.calc)
      if(calc) call memo_set(this%memo, "energy", 0.0_wp)
      call term_dict_init(iter, this%hdct)
      do
        nullify(husk, term)
        call term_dict_next(iter, husk, ierr)
        if(ierr/=TERM_DICT_OK)exit
        ASSERT(associated(husk))
        ASSERT(term_husk_assoc(husk))
        call term_husk_get(husk, term, factor=fctr, active=actv)
        ASSERT(associated(term))
        if(actv.and.(abs(fctr)>tiny(fctr)))then
          if(updt) call base_hamiltonian__madd__(this, fctr, term)
          if(calc)then
            call term_intrf_get(term, energy=enrg)
            call memo_acc(this%memo, "energy", fctr*enrg, ierr)
            ASSERT(ierr==MEMO_OK)
          end if
        end if
      end do
      call term_dict_end(iter)
      nullify(husk, term)
      if(updt)then
        call storage_update(this%data)
        call base_hamiltonian__publish__(this)
      end if
    end if

    POP_SUB(base_hamiltonian__update__)

  contains

    recursive subroutine update(this)
      type(term_intrf_t), intent(inout) :: this

      PUSH_SUB(base_hamiltonian__update__.update)

      call term_intrf__update__(this, energy)

      POP_SUB(base_hamiltonian__update__.update)
    end subroutine update

  end subroutine base_hamiltonian__update__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__reset__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call memo_del(this%memo, "energy")
    call storage_reset(this%data)

    POP_SUB(base_hamiltonian__reset__)
  end subroutine base_hamiltonian__reset__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__stop__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call term__apply__(this, term_intrf__stop__)
    call storage_stop(this%data)

    POP_SUB(base_hamiltonian__stop__)
  end subroutine base_hamiltonian__stop__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__recv__(this, update, channel)
    type(base_hamiltonian_t), intent(in)  :: this
    logical,                intent(out) :: update
    integer,      optional, intent(in)  :: channel

    type(msgbus_iterator_t)      :: iter
    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg
    logical                      :: updt
    integer                      :: chid, ierr

    PUSH_SUB(base_hamiltonian__recv__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    chid = 1
    if(present(channel)) chid = channel
    update = .false.
    call msgbus_init(iter, this%msgb, id=chid)
    do
      nullify(data, mssg)
      call msgbus_next(iter, mssg)
      if(.not.associated(mssg))exit
      call message_get(mssg, data)
      ASSERT(associated(data))
      call json_get(data, "update", updt, ierr)
      if(ierr==JSON_OK)then
        update = (update .or. updt)
        call msgbus_remove(iter, ierr=ierr)
        ASSERT(ierr==MSGBUS_OK)
      end if
    end do
    nullify(data, mssg)

    POP_SUB(base_hamiltonian__recv__)
  end subroutine base_hamiltonian__recv__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__publish__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg

    PUSH_SUB(base_hamiltonian__publish__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(data, mssg)
    mssg => message_new()
    call message_get(mssg, data)
    ASSERT(associated(data))
    call json_set(data, "update", .true.)
    nullify(data)
    call msgbus_publish(this%msgb, mssg)
    call message_del(mssg)
    nullify(mssg)

    POP_SUB(base_hamiltonian__publish__)
  end subroutine base_hamiltonian__publish__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_notify_type(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg

    PUSH_SUB(base_hamiltonian_notify_type)

    ASSERT(associated(this%config))
    nullify(data, mssg)
    mssg => message_new()
    call message_get(mssg, data)
    ASSERT(associated(data))
    call json_set(data, "update", .true.)
    nullify(data)
    call msgbus_notify(this%msgb, mssg)
    call message_del(mssg)
    nullify(mssg)

    POP_SUB(base_hamiltonian_notify_type)
  end subroutine base_hamiltonian_notify_type

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_notify_subs(this, name)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),       intent(in)    :: name

    type(base_hamiltonian_t), pointer :: subs

    PUSH_SUB(base_hamiltonian_notify_subs)

    ASSERT(associated(this%config))
    nullify(subs)
    call base_hamiltonian_gets(this, trim(adjustl(name)), subs)
    ASSERT(associated(subs))
    call base_hamiltonian_notify(subs)
    nullify(subs)

    POP_SUB(base_hamiltonian_notify_subs)
  end subroutine base_hamiltonian_notify_subs

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_start(this, sim)
    type(base_hamiltonian_t), intent(inout) :: this
    type(simulation_t),       intent(in)    :: sim
    PUSH_SUB(base_hamiltonian_start)

    call base_hamiltonian__apply__(this, start)

    POP_SUB(base_hamiltonian_start)

  contains
    
    subroutine start(this)
      type(base_hamiltonian_t), intent(inout) :: this

      PUSH_SUB(base_hamiltonian_start.start)

      call base_hamiltonian__start__(this, sim)
      
      POP_SUB(base_hamiltonian_start.start)
    end subroutine start
    
  end subroutine base_hamiltonian_start

  ! ---------------------------------------------------------
  function base_hamiltonian_calc(this, that) result(energy)
    type(base_hamiltonian_t),       intent(inout) :: this
    type(base_density_t), optional, intent(in)    :: that

    real(kind=wp) :: energy

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_hamiltonian_calc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(dnst)
    energy = 0.0_wp
    call base_hamiltonian_update(this)
    if(present(that))then
      energy = base_hamiltonian__calc__(this, that)
    else
      call base_hamiltonian_get(this, dnst)
      ASSERT(associated(dnst))
      energy = base_hamiltonian__calc__(this, dnst)
      nullify(dnst)
    end if

    POP_SUB(base_hamiltonian_calc)
  end function base_hamiltonian_calc

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_update(this, energy)
    type(base_hamiltonian_t), intent(inout) :: this
    logical,        optional, intent(in)    :: energy

    PUSH_SUB(base_hamiltonian_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_hamiltonian__apply__(this, update, enforce_active=.true.)

    POP_SUB(base_hamiltonian_update)

  contains

    recursive subroutine update(this)
      type(base_hamiltonian_t), intent(inout) :: this

      PUSH_SUB(base_hamiltonian_update.update)

      call base_hamiltonian__update__(this, energy)

      POP_SUB(base_hamiltonian_update.update)
    end subroutine update

  end subroutine base_hamiltonian_update

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_reset(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian_reset)

    call base_hamiltonian__apply__(this, base_hamiltonian__reset__)
    
    POP_SUB(base_hamiltonian_reset)
  end subroutine base_hamiltonian_reset

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian_stop(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian_stop)

    call base_hamiltonian__apply__(this, base_hamiltonian__stop__)
    
    POP_SUB(base_hamiltonian_stop)
  end subroutine base_hamiltonian_stop

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__madd__(this, factor, that)
    type(base_hamiltonian_t), intent(inout) :: this
    real(kind=wp),            intent(in)    :: factor
    type(term_intrf_t),       intent(in)    :: that

    type(storage_t), pointer :: data

    PUSH_SUB(base_hamiltonian__madd__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(data)
    if(abs(factor)>tiny(factor))then
      call term_intrf_get(that, data)
      if(associated(data))then
        if(abs(1.0_wp-abs(factor))>epsilon(factor))then
          call storage_madd(this%data, factor, data)
        else
          if(factor>0.0_wp)then
            call storage_add(this%data, data)
          else
            call storage_sub(this%data, data)
          end if
        end if
        nullify(data)
      end if
    end if

    POP_SUB(base_hamiltonian__madd__)
  end subroutine base_hamiltonian__madd__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__sets__info(this, name, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    type(term_dict_iterator_t)  :: iter
    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__sets__info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    call term_dict_init(iter, this%hdct)
    do
      nullify(husk, term)
      call term_dict_next(iter, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, term)
      ASSERT(associated(term))
      ASSERT(term_intrf_assoc(term))
      nullify(husk)
      if(term_intrf_exists(term, trim(adjustl(name))))&
        call term_intrf_sets(term, trim(adjustl(name)), lock=lock, active=active)
    end do
    call term_dict_end(iter)
    nullify(husk, term)

    POP_SUB(base_hamiltonian__sets__info)
  end subroutine base_hamiltonian__sets__info
  
  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__sets__type(this, name, that, config, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_hamiltonian_t), intent(in)    :: that
    type(json_object_t),      intent(in)    :: config
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: snam
    type(term_husk_t),               pointer :: husk
    type(term_intrf_t),              pointer :: mtrm, strm
    logical                                  :: ilck, actv
    integer                                  :: ierr

    PUSH_SUB(base_hamiltonian__sets__type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call term_dict_init(iter, that%hdct)
    do
      nullify(husk, mtrm, strm)
      call term_dict_next(iter, snam, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, strm)
      ASSERT(associated(strm))
      ASSERT(term_intrf_assoc(strm))
      nullify(husk)
      call term_dict_get(this%hdct, trim(adjustl(snam)), husk, ierr)
      if(ierr/=TERM_DICT_OK)cycle
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, mtrm, lock=ilck, active=actv)
      ASSERT(associated(mtrm))
      ASSERT(term_intrf_assoc(mtrm))
      if(.not.actv.or.ilck)cycle
      call term_intrf_sets(mtrm, trim(adjustl(name)), strm, config, lock=lock, active=active)
    end do
    call term_dict_end(iter)
    nullify(husk, mtrm, strm)

    POP_SUB(base_hamiltonian__sets__type)
  end subroutine base_hamiltonian__sets__type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__dels__(this, name, that, config, lock, active)
    type(base_hamiltonian_t),      intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(base_hamiltonian_t),      intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: snam
    type(term_husk_t),               pointer :: husk
    type(term_intrf_t),              pointer :: mtrm, strm
    integer                                  :: ierr

    PUSH_SUB(base_hamiltonian__dels__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call term_dict_init(iter, that%hdct)
    do
      nullify(husk, mtrm, strm)
      call term_dict_next(iter, snam, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, strm)
      ASSERT(associated(strm))
      ASSERT(term_intrf_assoc(strm))
      nullify(husk)
      call term_dict_get(this%hdct, trim(adjustl(snam)), husk, ierr)
      if(ierr/=TERM_DICT_OK)cycle
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, mtrm)
      ASSERT(associated(mtrm))
      ASSERT(term_intrf_assoc(mtrm))
      call term_intrf_dels(mtrm, trim(adjustl(name)), strm)
    end do
    call term_dict_end(iter)
    nullify(husk, mtrm, strm)
    if(present(config)) continue
    if(present(lock))   continue
    if(present(active)) continue
    
    POP_SUB(base_hamiltonian__dels__)
  end subroutine base_hamiltonian__dels__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__set__info(this, name, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    type(msgbus_t),     pointer :: msgb
    logical                     :: actv
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__set__info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    nullify(husk, term, msgb)
     call term_dict_get(this%hdct, trim(adjustl(name)), husk, ierr=ierr)
    ASSERT(ierr==TERM_DICT_OK)
    ASSERT(associated(husk))
    ASSERT(term_husk_assoc(husk))
    if(present(active))then
      call term_husk_get(husk, term, active=actv)
      ASSERT(associated(term))
      if(actv.neqv.active)then
        call base_hamiltonian_notify(this)
        call term_intrf_get(term, msgb)
        ASSERT(associated(msgb))
        if(active)then
          call msgbus_attach(this%msgb, msgb, id=2)
        else
          call msgbus_detach(this%msgb, msgb, id=2)
        end if
        nullify(term, msgb)
      end if
    end if
    call term_husk_set(husk, factor=factor, lock=lock, active=active)
    nullify(husk)

    POP_SUB(base_hamiltonian__set__info)
  end subroutine base_hamiltonian__set__info
  
  ! ---------------------------------------------------------
  subroutine base_hamiltonian__set__type(this, name, that, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(term_intrf_t),       intent(in)    :: that
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    type(term_intrf_t), pointer :: term
    type(msgbus_t),     pointer :: msgb
    logical                     :: actv

    PUSH_SUB(base_hamiltonian__set__type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    nullify(term, msgb)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term)) call term_intrf_del(term)
    nullify(term)
    ASSERT(term_intrf_assoc(that))
    call base_hamiltonian_notify(this)
    call term_dict_set(this%hdct, trim(adjustl(name)), term_husk_new(that, factor=factor, lock=lock, active=active))
    actv = .true.
    if(present(active)) actv = active
    if(actv)then
      call term_intrf_get(that, msgb)
      ASSERT(associated(msgb))
      call msgbus_attach(this%msgb, msgb, id=2)
      nullify(msgb)
    end if

    POP_SUB(base_hamiltonian__set__type)
  end subroutine base_hamiltonian__set__type
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian__iget__(this, name, that, ierr)
    type(base_hamiltonian_t),   intent(in)  :: this
    character(len=*),           intent(in)  :: name
    type(term_husk_t), pointer, intent(out) :: that
    integer,          optional, intent(out) :: ierr

    type(base_hamiltonian_t), pointer :: hmlt
    integer                           :: ipos, jerr

    PUSH_SUB(base_hamiltonian__iget__)

    nullify(that, hmlt)
    ipos = index(name, "/", back=.true.)
    if(ipos>0)then
      call base_hamiltonian_gets(this, trim(adjustl(name(1:ipos-1))), hmlt, ierr=jerr)
      if(jerr==BASE_HAMILTONIAN_OK)then
        ASSERT(associated(hmlt))
        call term_dict_get(hmlt%hdct, trim(adjustl(name(ipos+1:))), that, ierr=jerr)
      end if
      nullify(hmlt)
    else
      call term_dict_get(this%hdct, trim(adjustl(name)), that, ierr=jerr)
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(base_hamiltonian__iget__)
  end subroutine base_hamiltonian__iget__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__get__(this, name, that, factor, lock, active, ierr)
    type(base_hamiltonian_t),    intent(in)  :: this
    character(len=*),            intent(in)  :: name
    type(term_intrf_t), pointer, intent(out) :: that
    real(kind=wp),     optional, intent(out) :: factor
    logical,           optional, intent(out) :: lock
    logical,           optional, intent(out) :: active
    integer,           optional, intent(out) :: ierr

    type(term_husk_t), pointer :: husk
    integer                    :: jerr

    PUSH_SUB(base_hamiltonian__get__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(that, husk)
    call base_hamiltonian__iget__(this, trim(adjustl(name)), husk, ierr=jerr)
    if(jerr==TERM_DICT_OK)then
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, that, factor=factor, lock=lock, active=active)
      ASSERT(associated(that))
      ASSERT(term_intrf_assoc(that))
    end if
    nullify(husk)
    if(present(ierr)) ierr = jerr

    POP_SUB(base_hamiltonian__get__)
  end subroutine base_hamiltonian__get__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__del__(this, name, that, factor, lock, active)
    type(base_hamiltonian_t),    intent(inout) :: this
    character(len=*),            intent(in)    :: name
    type(term_intrf_t), pointer, intent(out)   :: that
    real(kind=wp),     optional, intent(out)   :: factor
    logical,           optional, intent(out)   :: lock
    logical,           optional, intent(out)   :: active

    type(term_husk_t), pointer :: husk
    type(msgbus_t),    pointer :: msgb
    logical                    :: actv
    integer                    :: ierr

    PUSH_SUB(base_hamiltonian__del__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(that, husk, msgb)
    call term_dict_del(this%hdct, trim(adjustl(name)), husk, ierr)
    if(ierr==TERM_DICT_OK)then
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call base_hamiltonian_notify(this)
      call term_husk_get(husk, that, factor=factor, lock=lock, active=actv)
      call term_husk_del(husk)
      ASSERT(associated(that))
      ASSERT(term_intrf_assoc(that))
      if(actv)then
        call term_intrf_get(that, msgb)
        ASSERT(associated(msgb))
        call msgbus_detach(this%msgb, msgb, id=2)
        nullify(msgb)
      end if
      if(present(active)) active = actv
    end if
    nullify(husk)

    POP_SUB(base_hamiltonian__del__)
  end subroutine base_hamiltonian__del__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_info(this, name, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian_set_info)

    call base_hamiltonian__set__(this, name, factor=factor, lock=lock, active=active)

    POP_SUB(base_hamiltonian_set_info)
  end subroutine base_hamiltonian_set_info
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_term(this, name, that, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_term_t),        intent(in)    :: that
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian_set_term)

    call base_hamiltonian__set__(this, name, term_intrf_new(that), factor=factor, lock=lock, active=active)

    POP_SUB(base_hamiltonian_set_term)
  end subroutine base_hamiltonian_set_term
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_potn(this, name, that, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_potential_t),   intent(in)    :: that
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian_set_potn)

    call base_hamiltonian__set__(this, name, term_intrf_new(that), factor=factor, lock=lock, active=active)

    POP_SUB(base_hamiltonian_set_potn)
  end subroutine base_hamiltonian_set_potn
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_fnct(this, name, that, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_functional_t),  intent(in)    :: that
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian_set_fnct)

    call base_hamiltonian__set__(this, name, term_intrf_new(that), factor=factor, lock=lock, active=active)

    POP_SUB(base_hamiltonian_set_fnct)
  end subroutine base_hamiltonian_set_fnct
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_hmlt(this, name, that, factor, lock, active)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_hamiltonian_t), intent(in)    :: that
    real(kind=wp),  optional, intent(in)    :: factor
    logical,        optional, intent(in)    :: lock
    logical,        optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian_set_hmlt)

    call base_hamiltonian__set__(this, name, term_intrf_new(that), factor=factor, lock=lock, active=active)

    POP_SUB(base_hamiltonian_set_hmlt)
  end subroutine base_hamiltonian_set_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_info(this, size, nspin, use)
    type(base_hamiltonian_t), intent(in)  :: this
    integer,        optional, intent(out) :: size
    integer,        optional, intent(out) :: nspin
    logical,        optional, intent(out) :: use

    PUSH_SUB(base_hamiltonian_get_info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    call storage_get(this%data, dim=nspin, size=size, alloc=use)

    POP_SUB(base_hamiltonian_get_info)
  end subroutine base_hamiltonian_get_info

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_energy(this, energy)
    type(base_hamiltonian_t), intent(inout) :: this
    real(kind=wp),            intent(out)   :: energy

    integer :: ierr
    
    PUSH_SUB(base_hamiltonian_get_energy)
    
    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    call base_hamiltonian_update(this, energy=.true.)
    call memo_get(this%memo, "energy", energy, ierr=ierr)
    ASSERT(ierr==MEMO_OK)
    
    POP_SUB(base_hamiltonian_get_energy)
  end subroutine base_hamiltonian_get_energy

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_config(this, that)
    type(base_hamiltonian_t),     intent(in)  :: this
    type(json_object_t), pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_hamiltonian_get_config)
  end subroutine base_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_system(this, that)
    type(base_hamiltonian_t),     intent(in)  :: this
    type(base_system_t), pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_hamiltonian_get_system)
  end subroutine base_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_density(this, that)
    type(base_hamiltonian_t), target, intent(in)  :: this
    type(base_density_t),    pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_density)

    nullify(that)
    if(associated(this%sys))&
      call base_system_get(this%sys, that)

    POP_SUB(base_hamiltonian_get_density)
  end subroutine base_hamiltonian_get_density

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_simulation(this, that)
    type(base_hamiltonian_t),    intent(in)  :: this
    type(simulation_t), pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_hamiltonian_get_simulation)
  end subroutine base_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_storage(this, that)
    type(base_hamiltonian_t), target, intent(in)  :: this
    type(storage_t),         pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_storage)

    nullify(that)
    if(associated(this%config)) that => this%data

    POP_SUB(base_hamiltonian_get_storage)
  end subroutine base_hamiltonian_get_storage

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_msgbus(this, that)
    type(base_hamiltonian_t), target, intent(in)  :: this
    type(msgbus_t),          pointer, intent(out) :: that

    PUSH_SUB(base_hamiltonian_get_msgbus)

    nullify(that)
    if(associated(this%config)) that => this%msgb

    POP_SUB(base_hamiltonian_get_msgbus)
  end subroutine base_hamiltonian_get_msgbus

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_term(this, name, that, factor, lock, active)
    type(base_hamiltonian_t),   intent(in)  :: this
    character(len=*),           intent(in)  :: name
    type(base_term_t), pointer, intent(out) :: that
    real(kind=wp),    optional, intent(out) :: factor
    logical,          optional, intent(out) :: lock
    logical,          optional, intent(out) :: active

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_get_term)

    nullify(that, term)
    call base_hamiltonian__get__(this, trim(adjustl(name)), term, factor=factor, lock=lock, active=active)
    if(associated(term)) call term_intrf_get(term, that)
    nullify(term)

    POP_SUB(base_hamiltonian_get_term)
  end subroutine base_hamiltonian_get_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_potn(this, name, that, factor, lock, active)
    type(base_hamiltonian_t),        intent(in)  :: this
    character(len=*),                intent(in)  :: name
    type(base_potential_t), pointer, intent(out) :: that
    real(kind=wp),         optional, intent(out) :: factor
    logical,               optional, intent(out) :: lock
    logical,               optional, intent(out) :: active

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_get_potn)

    nullify(that, term)
    call base_hamiltonian__get__(this, trim(adjustl(name)), term, factor=factor, lock=lock, active=active)
    if(associated(term)) call term_intrf_get(term, that)
    nullify(term)

    POP_SUB(base_hamiltonian_get_potn)
  end subroutine base_hamiltonian_get_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_fnct(this, name, that, factor, lock, active)
    type(base_hamiltonian_t),         intent(in)  :: this
    character(len=*),                 intent(in)  :: name
    type(base_functional_t), pointer, intent(out) :: that
    real(kind=wp),          optional, intent(out) :: factor
    logical,                optional, intent(out) :: lock
    logical,                optional, intent(out) :: active

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_get_fnct)

    nullify(that, term)
    call base_hamiltonian__get__(this, trim(adjustl(name)), term, factor=factor, lock=lock, active=active)
    if(associated(term)) call term_intrf_get(term, that)
    nullify(term)

    POP_SUB(base_hamiltonian_get_fnct)
  end subroutine base_hamiltonian_get_fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hmlt(this, name, that, factor, lock, active)
    type(base_hamiltonian_t),          intent(in)  :: this
    character(len=*),                  intent(in)  :: name
    type(base_hamiltonian_t), pointer, intent(out) :: that
    real(kind=wp),           optional, intent(out) :: factor
    logical,                 optional, intent(out) :: lock
    logical,                 optional, intent(out) :: active

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_get_hmlt)

    nullify(that, term)
    call base_hamiltonian__get__(this, trim(adjustl(name)), term, factor=factor, lock=lock, active=active)
    if(associated(term)) call term_intrf_get(term, that)
    nullify(term)

    POP_SUB(base_hamiltonian_get_hmlt)
  end subroutine base_hamiltonian_get_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_data_r1(this, that)
    type(base_hamiltonian_t),             intent(inout) :: this
    real(kind=wp), dimension(:), pointer, intent(out)   :: that

    logical :: fuse

    PUSH_SUB(base_hamiltonian_get_data_r1)

    nullify(that)
    call base_hamiltonian_get(this, use=fuse)
    if(fuse)then
      call base_hamiltonian_update(this)
      call storage_get(this%data, that)
    end if

    POP_SUB(base_hamiltonian_get_data_r1)
  end subroutine base_hamiltonian_get_data_r1

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_data_r2(this, that)
    type(base_hamiltonian_t),               intent(inout) :: this
    real(kind=wp), dimension(:,:), pointer, intent(out)   :: that

    logical :: fuse

    PUSH_SUB(base_hamiltonian_get_data_r2)

    nullify(that)
    call base_hamiltonian_get(this, use=fuse)
    if(fuse)then
      call base_hamiltonian_update(this)
      call storage_get(this%data, that)
    end if

    POP_SUB(base_hamiltonian_get_data_r2)
  end subroutine base_hamiltonian_get_data_r2

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_none(this, name)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_del_none)

    nullify(term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term)) call term_intrf_del(term)
    nullify(term)

    POP_SUB(base_hamiltonian_del_none)
  end subroutine base_hamiltonian_del_none

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_term(this, name, that)
    type(base_hamiltonian_t),   intent(inout) :: this
    character(len=*),           intent(in)    :: name
    type(base_term_t), pointer, intent(out)   :: that

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_del_term)

    nullify(that, term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term))then
      call term_intrf_get(term, that)
      ASSERT(associated(that))
      call term_intrf_del(term)
    end if
    nullify(term)

    POP_SUB(base_hamiltonian_del_term)
  end subroutine base_hamiltonian_del_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_potn(this, name, that)
    type(base_hamiltonian_t),        intent(inout) :: this
    character(len=*),                intent(in)    :: name
    type(base_potential_t), pointer, intent(out)   :: that

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_del_potn)

    nullify(that, term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term))then
      call term_intrf_get(term, that)
      ASSERT(associated(that))
      call term_intrf_del(term)
    end if
    nullify(term)

    POP_SUB(base_hamiltonian_del_potn)
  end subroutine base_hamiltonian_del_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_fnct(this, name, that)
    type(base_hamiltonian_t),         intent(inout) :: this
    character(len=*),                 intent(in)    :: name
    type(base_functional_t), pointer, intent(out)   :: that

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_del_fnct)

    nullify(that, term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term))then
      call term_intrf_get(term, that)
      ASSERT(associated(that))
      call term_intrf_del(term)
    end if
    nullify(term)

    POP_SUB(base_hamiltonian_del_fnct)
  end subroutine base_hamiltonian_del_fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_hmlt(this, name, that)
    type(base_hamiltonian_t),          intent(inout) :: this
    character(len=*),                  intent(in)    :: name
    type(base_hamiltonian_t), pointer, intent(out)   :: that

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian_del_hmlt)

    nullify(that, term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term))then
      call term_intrf_get(term, that)
      ASSERT(associated(that))
      call term_intrf_del(term)
    end if
    nullify(term)

    POP_SUB(base_hamiltonian_del_hmlt)
  end subroutine base_hamiltonian_del_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__copy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    type(term_husk_t),               pointer :: husk
    type(term_intrf_t),              pointer :: term
    type(refcount_t),                pointer :: rcnt
    real(kind=wp)                            :: fctr
    logical                                  :: lock, actv
    integer                                  :: ierr

    PUSH_SUB(base_hamiltonian__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt, husk, term)
    call base_hamiltonian__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_hamiltonian__init__(this, that)
      call refcount_del(this%rcnt)
      call term_dict_init(iter, that%hdct)
      do
        nullify(husk, term)
        call term_dict_next(iter, name, husk, ierr)
        if(ierr/=TERM_DICT_OK)exit
        ASSERT(associated(husk))
        ASSERT(term_husk_assoc(husk))
        call term_husk_get(husk, term, factor=fctr, lock=lock, active=actv)
        ASSERT(associated(term))
        call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(source=term), factor=fctr, lock=lock, active=actv)
      end do
      call term_dict_end(iter)
      nullify(husk, term)
      call memo_copy(this%memo, that%memo)
      if(associated(that%sim)) call storage_copy(this%data, that%data)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_hamiltonian__copy__)
  end subroutine base_hamiltonian__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__end__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(term_husk_t),  pointer :: husk
    type(term_intrf_t), pointer :: term
    integer                     :: ierr

    PUSH_SUB(base_hamiltonian__end__)

    nullify(this%config, this%sys, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    this%nspin = default_nspin
    do
      nullify(husk, term)
      call term_dict_pop(this%hdct, husk, ierr)
      if(ierr/=TERM_DICT_OK)exit
      ASSERT(associated(husk))
      ASSERT(term_husk_assoc(husk))
      call term_husk_get(husk, term)
      ASSERT(associated(term))
      call term_husk_del(husk)
      call term_intrf_del(term)
    end do
    nullify(husk, term)
    call memo_end(this%memo)
    call storage_end(this%data)
    ASSERT(term_dict_len(this%hdct)==0)
    call term_dict_end(this%hdct)
    ASSERT(base_hamiltonian_dict_len(this%dict)==0)
    call msgbus_end(this%msgb)
    call base_hamiltonian_dict_end(this%dict)

    POP_SUB(base_hamiltonian__end__)
  end subroutine base_hamiltonian__end__

end module base_hamiltonian_oct_m

!! Local Variables:
!! mode: f90
!! End:
