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

  use base_functional_oct_m
  use base_potential_oct_m
  use base_system_oct_m
  use base_term_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_intrf_t
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
    base_hamiltonian__acc__,    &
    base_hamiltonian__sub__,    &
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
    base_hamiltonian_update, &
    base_hamiltonian_stop,   &
    base_hamiltonian_set,    &
    base_hamiltonian_get,    &
    base_hamiltonian_copy,   &
    base_hamiltonian_end

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_intrf_t
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

  integer, parameter :: TERM_STAT_DISA = 0
  integer, parameter :: TERM_STAT_NULL = 1
  integer, parameter :: TERM_STAT_ASSC = 2
  integer, parameter :: TERM_STAT_ALLC = 3
  
  integer, parameter :: TERM_TYPE_NONE = 0
  integer, parameter :: TERM_TYPE_TERM = 1
  integer, parameter :: TERM_TYPE_POTN = 2
  integer, parameter :: TERM_TYPE_FNCT = 3
  integer, parameter :: TERM_TYPE_HMLT = 4

  integer, parameter :: default_nspin = 1

  type :: term_intrf_t
    private
    type(base_term_t),        pointer :: term =>null()
    type(base_potential_t),   pointer :: potn =>null()
    type(base_functional_t),  pointer :: fnct =>null()
    type(base_hamiltonian_t), pointer :: hmlt =>null()
    integer                           :: stat = TERM_STAT_DISA
    integer                           :: type = TERM_TYPE_NONE
  end type term_intrf_t

  type :: base_hamiltonian_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(base_system_t),  pointer :: sys    =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(refcount_t),     pointer :: rcnt   =>null()
    integer                       :: nspin  = default_nspin
    real(kind=wp)                 :: energy = 0.0_wp
    type(storage_t)               :: data
    type(term_dict_t)             :: hdct
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
  end interface term_intrf_get

  interface base_hamiltonian__init__
    module procedure base_hamiltonian__init__type
    module procedure base_hamiltonian__init__copy
  end interface base_hamiltonian__init__

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

  interface base_hamiltonian__sets__
    module procedure base_hamiltonian__sets__info
    module procedure base_hamiltonian__sets__type
  end interface base_hamiltonian__sets__

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

  interface base_hamiltonian_set
    module procedure base_hamiltonian_set_info
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
    module procedure base_hamiltonian_get_simulation
    module procedure base_hamiltonian_get_term
    module procedure base_hamiltonian_get_potn
    module procedure base_hamiltonian_get_fnct
    module procedure base_hamiltonian_get_hmlt
    module procedure base_hamiltonian_get_storage
    module procedure base_hamiltonian_get_hamiltonian_1d
    module procedure base_hamiltonian_get_hamiltonian_md
  end interface base_hamiltonian_get

contains

#define DICT_TEMPLATE_NAME term
#define DICT_TYPE_NAME term_intrf_t
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
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

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

    ASSERT(present(source).or.present(mold))
    ASSERT(.not.(present(source).and.present(mold)))
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
    type(term_intrf_t), pointer :: this

    PUSH_SUB(term_intrf_del)

    if(associated(this))then
      call term_intrf_end(this)
      SAFE_DEALLOCATE_P(this)
    end if

    POP_SUB(term_intrf_del)
  end subroutine term_intrf_del

  ! ---------------------------------------------------------
  function term_intrf_assoc(this) result(that)
    type(term_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(term_intrf_assoc)

    select case(this%stat)
    case(TERM_STAT_DISA)
      that = .false.
    case(TERM_STAT_NULL)
      ASSERT(.not.associated(this%term))
      ASSERT(.not.associated(this%potn))
      ASSERT(.not.associated(this%fnct))
      ASSERT(.not.associated(this%hmlt))
      that = .false.
    case(TERM_STAT_ASSC, TERM_STAT_ALLC)
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
  function term_intrf_alloc(this) result(that)
    type(term_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(term_intrf_alloc)

    if(term_intrf_assoc(this))then
      select case(this%stat)
      case(TERM_STAT_ASSC)
        that = .false.
      case(TERM_STAT_ALLC)
        that = .true.
      case default
        ASSERT(.false.)
      end select
    end if

    POP_SUB(term_intrf_alloc)
  end function term_intrf_alloc

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_init_type(this, sys, config)
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
    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    call json_get(config, "type", type, ierr)
    ASSERT(ierr==JSON_OK)
    select case(type)
    case(TERM_TYPE_TERM)
      SAFE_ALLOCATE(term)
      call base_term_init(term, sys, config)
      call term_intrf_set(this, term)
      nullify(term)
    case(TERM_TYPE_POTN)
      SAFE_ALLOCATE(potn)
      call base_potential_init(potn, sys, config)
      call term_intrf_set(this, potn)
      nullify(potn)
    case(TERM_TYPE_FNCT)
      SAFE_ALLOCATE(fnct)
      call base_functional_init(fnct, sys, config)
      call term_intrf_set(this, fnct)
      nullify(fnct)
    case(TERM_TYPE_HMLT)
      SAFE_ALLOCATE(hmlt)
      call base_hamiltonian_init(hmlt, sys, config)
      call term_intrf_set(this, hmlt)
      nullify(hmlt)
    case default
      ASSERT(.false.)
    end select
    this%stat = TERM_STAT_ALLC

    POP_SUB(term_intrf_init_type)
  end subroutine term_intrf_init_type

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_init_copy(this, that)
    type(term_intrf_t), intent(out) :: this
    type(term_intrf_t), intent(in)  :: that

    type(base_term_t),        pointer :: itrm, otrm
    type(base_potential_t),   pointer :: iptn, optn
    type(base_functional_t),  pointer :: ifnc, ofnc
    type(base_hamiltonian_t), pointer :: ihml, ohml

    PUSH_SUB(term_intrf_init_copy)

    nullify(itrm, iptn, ifnc, ihml)
    nullify(otrm, optn, ofnc, ohml)
    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    select case(that%stat)
    case(TERM_STAT_NULL)
    case(TERM_STAT_ASSC, TERM_STAT_ALLC)
      select case(that%type)
      case(TERM_TYPE_TERM)
        call term_intrf_get(that, itrm)
        ASSERT(associated(itrm))
        SAFE_ALLOCATE(otrm)
        call base_term_init(otrm, itrm)
        call term_intrf_set(this, otrm)
        nullify(itrm, otrm)
      case(TERM_TYPE_POTN)
        call term_intrf_get(that, iptn)
        ASSERT(associated(iptn))
        SAFE_ALLOCATE(optn)
        call base_potential_init(optn, iptn)
        call term_intrf_set(this, optn)
        nullify(iptn, optn)
      case(TERM_TYPE_FNCT)
        call term_intrf_get(that, ifnc)
        ASSERT(associated(ifnc))
        SAFE_ALLOCATE(ofnc)
        call base_functional_copy(ofnc, ifnc)
        call term_intrf_set(this, ofnc)
        nullify(ifnc, ofnc)
      case(TERM_TYPE_HMLT)
        call term_intrf_get(that, ihml)
        ASSERT(associated(ihml))
        SAFE_ALLOCATE(ohml)
        call base_hamiltonian_init(ohml, ihml)
        call term_intrf_set(this, ohml)
        nullify(ihml, ohml)
      case default
        ASSERT(.false.)
      end select
      this%stat = TERM_STAT_ALLC
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(term_intrf_init_copy)
  end subroutine term_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine term_intrf_init_term(this, that)
    type(term_intrf_t), intent(out) :: this
    type(base_term_t),  intent(in)  :: that

    PUSH_SUB(term_intrf_init_term)

    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_term)
  end subroutine term_intrf_init_term

  ! ---------------------------------------------------------
  subroutine term_intrf_init_potn(this, that)
    type(term_intrf_t),     intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_potn)

    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_potn)
  end subroutine term_intrf_init_potn

  ! ---------------------------------------------------------
  subroutine term_intrf_init_fnct(this, that)
    type(term_intrf_t),      intent(out) :: this
    type(base_functional_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_fnct)

    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_fnct)
  end subroutine term_intrf_init_fnct

  ! ---------------------------------------------------------
  subroutine term_intrf_init_hmlt(this, that)
    type(term_intrf_t),       intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    PUSH_SUB(term_intrf_init_hmlt)

    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    call term_intrf_set(this, that)

    POP_SUB(term_intrf_init_hmlt)
  end subroutine term_intrf_init_hmlt

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
  recursive subroutine term_intrf__acc__(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(term_intrf_t), intent(in)    :: that

    PUSH_SUB(term_intrf__acc__)

    ASSERT(term_intrf_assoc(this))
    ASSERT(term_intrf_assoc(that))
    ASSERT(this%type==that%type)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__acc__(this%term, that%term)
    case(TERM_TYPE_POTN)
      call base_potential__acc__(this%potn, that%potn)
    case(TERM_TYPE_FNCT)
      call base_functional__acc__(this%fnct, that%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__acc__(this%hmlt, that%hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__acc__)
  end subroutine term_intrf__acc__

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__sub__(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(term_intrf_t), intent(in)    :: that

    PUSH_SUB(term_intrf__sub__)

    ASSERT(term_intrf_assoc(this))
    ASSERT(term_intrf_assoc(that))
    ASSERT(this%type==that%type)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__sub__(this%term, that%term)
    case(TERM_TYPE_POTN)
      call base_potential__sub__(this%potn, that%potn)
    case(TERM_TYPE_FNCT)
      call base_functional__sub__(this%fnct, that%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__sub__(this%hmlt, that%hmlt)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf__sub__)
  end subroutine term_intrf__sub__

  ! ---------------------------------------------------------
  recursive subroutine term_intrf__update__(this)
    type(term_intrf_t), intent(inout) :: this

    PUSH_SUB(term_intrf__update__)

    ASSERT(term_intrf_assoc(this))
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term__update__(this%term)
    case(TERM_TYPE_POTN)
      call base_potential__update__(this%potn)
    case(TERM_TYPE_FNCT)
      call base_functional__update__(this%fnct)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian__update__(this%hmlt)
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
  recursive subroutine term_intrf_sets(this, name, that, config, lock, active)
    type(term_intrf_t),            intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(term_intrf_t),            intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    PUSH_SUB(term_intrf_sets)

    ASSERT(term_intrf_assoc(this))
    ASSERT(term_intrf_assoc(that))
    ASSERT(this%type==that%type)
    select case(this%type)
    case(TERM_TYPE_TERM)
      call base_term_sets(this%term, trim(adjustl(name)), that%term, config=config, lock=lock, active=active)
    case(TERM_TYPE_POTN)
      call base_potential_sets(this%potn, trim(adjustl(name)), that%potn, config=config, lock=lock, active=active)
    case(TERM_TYPE_FNCT)
      call base_functional_sets(this%fnct, trim(adjustl(name)), that%fnct, config=config, lock=lock, active=active)
    case(TERM_TYPE_HMLT)
      call base_hamiltonian_sets(this%hmlt, trim(adjustl(name)), that%hmlt, config=config, lock=lock, active=active)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(term_intrf_sets)
  end subroutine term_intrf_sets

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

    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(this%type==TERM_TYPE_NONE)
    this%term => that
    this%stat = TERM_STAT_ASSC
    this%type = TERM_TYPE_TERM

    POP_SUB(term_intrf_set_term)
  end subroutine term_intrf_set_term

  ! ---------------------------------------------------------
  subroutine term_intrf_set_potn(this, that)
    type(term_intrf_t),             intent(inout) :: this
    type(base_potential_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_potn)

    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(this%type==TERM_TYPE_NONE)
    this%potn => that
    this%stat = TERM_STAT_ASSC
    this%type = TERM_TYPE_POTN

    POP_SUB(term_intrf_set_potn)
  end subroutine term_intrf_set_potn

  ! ---------------------------------------------------------
  subroutine term_intrf_set_fnct(this, that)
    type(term_intrf_t),              intent(inout) :: this
    type(base_functional_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_fnct)

    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(this%type==TERM_TYPE_NONE)
    this%fnct => that
    this%stat = TERM_STAT_ASSC
    this%type = TERM_TYPE_FNCT

    POP_SUB(term_intrf_set_fnct)
  end subroutine term_intrf_set_fnct

  ! ---------------------------------------------------------
  subroutine term_intrf_set_hmlt(this, that)
    type(term_intrf_t),               intent(inout) :: this
    type(base_hamiltonian_t), target, intent(in)    :: that

    PUSH_SUB(term_intrf_set_hmlt)

    ASSERT(this%stat==TERM_STAT_NULL)
    ASSERT(this%type==TERM_TYPE_NONE)
    this%hmlt => that
    this%stat = TERM_STAT_ASSC
    this%type = TERM_TYPE_HMLT

    POP_SUB(term_intrf_set_hmlt)
  end subroutine term_intrf_set_hmlt

  ! ---------------------------------------------------------
  subroutine term_intrf_get_term(this, that)
    type(term_intrf_t), intent(in) :: this
    type(base_term_t), pointer     :: that

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
    type(term_intrf_t),      intent(in) :: this
    type(base_potential_t), pointer     :: that

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
    type(term_intrf_t),       intent(in) :: this
    type(base_functional_t), pointer     :: that

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
    type(term_intrf_t),        intent(in) :: this
    type(base_hamiltonian_t), pointer     :: that

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
  recursive subroutine term_intrf_copy(this, that)
    type(term_intrf_t), intent(inout) :: this
    type(term_intrf_t), intent(in)    :: that

    type(base_term_t),        pointer :: itrm, otrm
    type(base_potential_t),   pointer :: iptn, optn
    type(base_functional_t),  pointer :: ifnc, ofnc
    type(base_hamiltonian_t), pointer :: ihml, ohml

    PUSH_SUB(term_intrf_copy)

    nullify(itrm, iptn, ifnc, ihml)
    nullify(otrm, optn, ofnc, ohml)
    call term_intrf_end(this)
    this%stat = TERM_STAT_NULL
    this%type = TERM_TYPE_NONE
    select case(that%stat)
    case(TERM_STAT_DISA)
      select case(that%type)
      case(TERM_TYPE_NONE)
      case default
        ASSERT(.false.)
      end select
    case(TERM_STAT_NULL)
      select case(that%type)
      case(TERM_TYPE_NONE)
      case default
        ASSERT(.false.)
      end select
    case(TERM_STAT_ASSC)
      select case(that%type)
      case(TERM_TYPE_TERM)
        call term_intrf_get(that, itrm)
        ASSERT(associated(itrm))
        call term_intrf_set(this, itrm)
        nullify(itrm)
      case(TERM_TYPE_POTN)
        call term_intrf_get(that, iptn)
        ASSERT(associated(iptn))
        call term_intrf_set(this, iptn)
        nullify(iptn)
      case(TERM_TYPE_FNCT)
        call term_intrf_get(that, ifnc)
        ASSERT(associated(ifnc))
        call term_intrf_set(this, ifnc)
        nullify(ifnc)
      case(TERM_TYPE_HMLT)
        call term_intrf_get(that, ihml)
        ASSERT(associated(ihml))
        call term_intrf_set(this, ihml)
        nullify(ihml)
      case default
        ASSERT(.false.)
      end select
    case(TERM_STAT_ALLC)
      select case(that%type)
      case(TERM_TYPE_TERM)
        call term_intrf_get(that, itrm)
        ASSERT(associated(itrm))
        SAFE_ALLOCATE(otrm)
        call base_term_copy(otrm, itrm)
        call term_intrf_set(this, otrm)
        nullify(itrm, otrm)
      case(TERM_TYPE_POTN)
        call term_intrf_get(that, iptn)
        ASSERT(associated(iptn))
        SAFE_ALLOCATE(optn)
        call base_potential_copy(optn, iptn)
        call term_intrf_set(this, optn)
        nullify(iptn, optn)
      case(TERM_TYPE_FNCT)
        call term_intrf_get(that, ifnc)
        ASSERT(associated(ifnc))
        SAFE_ALLOCATE(ofnc)
        call base_functional_copy(ofnc, ifnc)
        call term_intrf_set(this, ofnc)
        nullify(ofnc)
      case(TERM_TYPE_HMLT)
        call term_intrf_get(that, ihml)
        ASSERT(associated(ihml))
        SAFE_ALLOCATE(ohml)
        call base_hamiltonian_init(ohml, ihml)
        call term_intrf_set(this, ihml)
        nullify(ihml, ohml)
      case default
        ASSERT(.false.)
      end select
    case default
      ASSERT(.false.)
    end select
    this%stat = that%stat
    this%type = that%type

    POP_SUB(term_intrf_copy)
  end subroutine term_intrf_copy

  ! ---------------------------------------------------------
  recursive subroutine term_intrf_end(this)
    type(term_intrf_t), intent(inout) :: this

    PUSH_SUB(term_intrf_end)

    select case(this%stat)
    case(TERM_STAT_ALLC)
      select case(this%type)
      case(TERM_TYPE_TERM)
        call base_term_end(this%term)
        SAFE_DEALLOCATE_P(this%term)
      case(TERM_TYPE_POTN)
        call base_potential_end(this%potn)
        SAFE_DEALLOCATE_P(this%potn)
      case(TERM_TYPE_FNCT)
        call base_functional_end(this%fnct)
        SAFE_DEALLOCATE_P(this%fnct)
      case(TERM_TYPE_HMLT)
        call base_hamiltonian_end(this%hmlt)
        SAFE_DEALLOCATE_P(this%hmlt)
      end select
    end select
    nullify(this%term, this%potn, this%fnct, this%hmlt)
    this%stat = TERM_STAT_DISA
    this%type = TERM_TYPE_NONE

    POP_SUB(term_intrf_end)
  end subroutine term_intrf_end

  ! ---------------------------------------------------------
  recursive subroutine term__apply__(this, operation)
    type(base_hamiltonian_t), intent(inout) :: this

    interface
      subroutine operation(this)
        import :: term_intrf_t
        type(term_intrf_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(term_dict_iterator_t)  :: iter
    type(term_intrf_t), pointer :: term

    PUSH_SUB(term__apply__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    call term_dict_init(iter, this%hdct)
    do
      nullify(term)
      call term_dict_next(iter, term)
      if(.not.associated(term)) exit
      ASSERT(term_intrf_assoc(term))
      if(.not.term_intrf_alloc(term)) cycle
      call operation(term)
    end do
    call term_dict_end(iter)
    nullify(term)

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
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "full", .false.)
    call storage_init(this%data, cnfg)
    nullify(cnfg)
    call term_dict_init(this%hdct)
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
    type(json_object_t),             pointer :: cnfg
    integer                                  :: type, ierr

    PUSH_SUB(base_hamiltonian__init__type)

    call base_hamiltonian__iinit__(this, sys, config)
    call json_init(iter, config)
    do
      nullify(cnfg)
      call json_next(iter, name, cnfg, ierr)
      if(ierr==JSON_TYPE_ERROR) cycle
      if(ierr/=JSON_OK) exit
      call json_get(cnfg, "type", type, ierr)
      if(ierr/=JSON_OK) cycle
      call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(sys, cnfg))
    end do
    call json_end(iter)
    nullify(cnfg)

    POP_SUB(base_hamiltonian__init__type)
  end subroutine base_hamiltonian__init__type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__init__copy(this, that)
    type(base_hamiltonian_t), intent(out) :: this
    type(base_hamiltonian_t), intent(in)  :: that

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    type(term_intrf_t),              pointer :: term

    PUSH_SUB(base_hamiltonian__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_hamiltonian__iinit__(this, that%sys, that%config)
    call term_dict_init(iter, that%hdct)
    do
      nullify(term)
      call term_dict_next(iter, name, term)
      if(.not.associated(term))exit
      call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(mold=term))
    end do
    call term_dict_end(iter)
    nullify(term)
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

    integer :: ierr
    logical :: uspn

    PUSH_SUB(base_hamiltonian__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    call term__apply__(this, start)
    this%sim => sim
    call base_system_get(this%sys, nspin=this%nspin)
    call json_get(this%config, "spin", uspn, ierr)
    if(ierr/=JSON_OK) uspn = .true.
    if(.not.uspn) this%nspin = default_nspin
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    call storage_start(this%data, sim, ndim=this%nspin)

    POP_SUB(base_hamiltonian__start__)

  contains

    subroutine start(this)
      type(term_intrf_t), intent(inout) :: this

      PUSH_SUB(base_hamiltonian__start__.start)

      call term_intrf__start__(this, sim)
      
      POP_SUB(base_hamiltonian__start__.start)
    end subroutine start
  
  end subroutine base_hamiltonian__start__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__update__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call term__apply__(this, term_intrf__update__)
    call storage_update(this%data)

    POP_SUB(base_hamiltonian__update__)
  end subroutine base_hamiltonian__update__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__reset__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call term__apply__(this, term_intrf__reset__)
    this%energy = 0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_hamiltonian__reset__)
  end subroutine base_hamiltonian__reset__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__stop__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call term__apply__(this, term_intrf__stop__)
    nullify(this%sim)
    call storage_stop(this%data)

    POP_SUB(base_hamiltonian__stop__)
  end subroutine base_hamiltonian__stop__

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
  recursive subroutine base_hamiltonian_update(this)
    type(base_hamiltonian_t), intent(inout) :: this

    PUSH_SUB(base_hamiltonian_update)

    call base_hamiltonian__apply__(this, base_hamiltonian__update__)
    
    POP_SUB(base_hamiltonian_update)
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

    type(term_dict_iterator_t)               :: iter
    type(term_intrf_t),              pointer :: mtrm, strm
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name

    PUSH_SUB(base_hamiltonian__acc__hmlt)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sim))
    this%energy = this%energy + that%energy
    call term_dict_init(iter, that%hdct)
    do
      nullify(mtrm, strm)
      call term_dict_next(iter, name, strm)
      if(.not.associated(strm)) exit
      call base_hamiltonian__get__(this, trim(adjustl(name)), mtrm)
      if(.not.associated(mtrm)) cycle
      call term_intrf__acc__(mtrm, strm)
    end do
    call term_dict_end(iter)
    nullify(mtrm, strm)
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

    type(term_dict_iterator_t)               :: iter
    type(term_intrf_t),              pointer :: mtrm, strm
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name

    PUSH_SUB(base_hamiltonian__sub__hmlt)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sim))
    this%energy = this%energy - that%energy
    call term_dict_init(iter, that%hdct)
    do
      nullify(mtrm, strm)
      call term_dict_next(iter, name, strm)
      if(.not.associated(strm)) exit
      call base_hamiltonian__get__(this, trim(adjustl(name)), mtrm)
      if(.not.associated(mtrm)) cycle
      call term_intrf__sub__(mtrm, strm)
    end do
    call term_dict_end(iter)
    nullify(mtrm, strm)
    call storage_sub(this%data, that%data)

    POP_SUB(base_hamiltonian__sub__hmlt)
  end subroutine base_hamiltonian__sub__hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__sets__info(this, name, config, lock, active)
    type(base_hamiltonian_t),      intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    PUSH_SUB(base_hamiltonian__sets__info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    if(present(config)) continue
    if(present(lock)) continue
    if(present(active)) continue

    POP_SUB(base_hamiltonian__sets__info)
  end subroutine base_hamiltonian__sets__info
  
  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__sets__type(this, name, that, config, lock, active)
    type(base_hamiltonian_t),      intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(base_hamiltonian_t),      intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: snam
    type(term_intrf_t),              pointer :: mtrm, strm

    PUSH_SUB(base_hamiltonian__sets__type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call term_dict_init(iter, that%hdct)
    do
      nullify(mtrm, strm)
      call term_dict_next(iter, snam, strm)
      if(.not.associated(strm)) exit
      call term_dict_get(this%hdct, trim(adjustl(snam)), mtrm)
      if(.not.associated(mtrm)) cycle
      call term_intrf_sets(mtrm, trim(adjustl(name)), strm, config=config, lock=lock, active=active)
    end do
    call term_dict_end(iter)
    nullify(mtrm, strm)

    POP_SUB(base_hamiltonian__sets__type)
  end subroutine base_hamiltonian__sets__type

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__dels__(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_hamiltonian_t), intent(in)    :: that

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: snam
    type(term_intrf_t),              pointer :: mtrm, strm

    PUSH_SUB(base_hamiltonian__dels__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call term_dict_init(iter, that%hdct)
    do
      nullify(mtrm, strm)
      call term_dict_next(iter, snam, strm)
      if(.not.associated(strm)) exit
      call term_dict_get(this%hdct, trim(adjustl(snam)), mtrm)
      if(.not.associated(mtrm)) cycle
      call term_intrf_dels(mtrm, trim(adjustl(name)), strm)
    end do
    call term_dict_end(iter)
    nullify(mtrm, strm)
    
    POP_SUB(base_hamiltonian__dels__)
  end subroutine base_hamiltonian__dels__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__set__(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(term_intrf_t),       intent(in)    :: that

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian__set__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(term)
    call base_hamiltonian__del__(this, trim(adjustl(name)), term)
    if(associated(term)) call term_intrf_del(term)
    nullify(term)
    call term_dict_set(this%hdct, trim(adjustl(name)), that)

    POP_SUB(base_hamiltonian__set__)
  end subroutine base_hamiltonian__set__
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian__get__(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(term_intrf_t),      pointer     :: that

    integer :: ierr

    PUSH_SUB(base_hamiltonian__get__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(that)
    call term_dict_get(this%hdct, name, that, ierr)
    if(ierr/=TERM_DICT_OK) nullify(that)

    POP_SUB(base_hamiltonian__get__)
  end subroutine base_hamiltonian__get__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__del__(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(term_intrf_t),      pointer        :: that

    integer :: ierr

    PUSH_SUB(base_hamiltonian__del__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    nullify(that)
    call term_dict_del(this%hdct, trim(adjustl(name)), that, ierr)
    if(ierr/=TERM_DICT_OK) nullify(that)

    POP_SUB(base_hamiltonian__del__)
  end subroutine base_hamiltonian__del__

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_info(this, energy)
    type(base_hamiltonian_t), intent(inout) :: this
    real(kind=wp),  optional, intent(in)    :: energy

    PUSH_SUB(base_hamiltonian_set_info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    if(present(energy)) this%energy = energy

    POP_SUB(base_hamiltonian_set_info)
  end subroutine base_hamiltonian_set_info

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_term(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_term_t),        intent(in)    :: that

    PUSH_SUB(base_hamiltonian_set_term)

    call base_hamiltonian__set__(this, name, term_intrf_new(that))

    POP_SUB(base_hamiltonian_set_term)
  end subroutine base_hamiltonian_set_term
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_potn(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_potential_t),   intent(in)    :: that

    PUSH_SUB(base_hamiltonian_set_potn)

    call base_hamiltonian__set__(this, name, term_intrf_new(that))

    POP_SUB(base_hamiltonian_set_potn)
  end subroutine base_hamiltonian_set_potn
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_fnct(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_functional_t),  intent(in)    :: that

    PUSH_SUB(base_hamiltonian_set_fnct)

    call base_hamiltonian__set__(this, name, term_intrf_new(that))

    POP_SUB(base_hamiltonian_set_fnct)
  end subroutine base_hamiltonian_set_fnct
    
  ! ---------------------------------------------------------
  subroutine base_hamiltonian_set_hmlt(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_hamiltonian_t), intent(in)    :: that

    PUSH_SUB(base_hamiltonian_set_hmlt)

    call base_hamiltonian__set__(this, name, term_intrf_new(that))

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
    type(base_hamiltonian_t), intent(in)  :: this
    real(kind=wp),            intent(out) :: energy
    
    PUSH_SUB(base_hamiltonian_get_energy)
    
    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    energy = this%energy
    
    POP_SUB(base_hamiltonian_get_energy)
  end subroutine base_hamiltonian_get_energy

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
  subroutine base_hamiltonian_get_storage(this, that)
    type(base_hamiltonian_t), target, intent(in) :: this
    type(storage_t),         pointer             :: that

    PUSH_SUB(base_hamiltonian_get_storage)

    nullify(that)
    if(associated(this%config)) that => this%data

    POP_SUB(base_hamiltonian_get_storage)
  end subroutine base_hamiltonian_get_storage

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hamiltonian_1d(this, that)
    type(base_hamiltonian_t),     intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_hamiltonian_get_hamiltonian_1d)

    nullify(that)
    call storage_get(this%data, that)

    POP_SUB(base_hamiltonian_get_hamiltonian_1d)
  end subroutine base_hamiltonian_get_hamiltonian_1d

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hamiltonian_md(this, that)
    type(base_hamiltonian_t),       intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(base_hamiltonian_get_hamiltonian_md)

    nullify(that)
    call storage_get(this%data, that)

    POP_SUB(base_hamiltonian_get_hamiltonian_md)
  end subroutine base_hamiltonian_get_hamiltonian_md

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_term(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_term_t),       pointer     :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_term)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, trim(adjustl(name)), htrm)
    if(associated(htrm)) call term_intrf_get(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_term)
  end subroutine base_hamiltonian_get_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_potn(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_potential_t),  pointer     :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_potn)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, trim(adjustl(name)), htrm)
    if(associated(htrm)) call term_intrf_get(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_potn)
  end subroutine base_hamiltonian_get_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_fnct(this, name, that)
    type(base_hamiltonian_t), intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_functional_t), pointer     :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_fnct)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, trim(adjustl(name)), htrm)
    if(associated(htrm)) call term_intrf_get(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_fnct)
  end subroutine base_hamiltonian_get_fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_get_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(in) :: this
    character(len=*),          intent(in) :: name
    type(base_hamiltonian_t), pointer     :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_get_hmlt)

    nullify(that, htrm)
    call base_hamiltonian__get__(this, trim(adjustl(name)), htrm)
    if(associated(htrm)) call term_intrf_get(htrm, that)
    nullify(htrm)

    POP_SUB(base_hamiltonian_get_hmlt)
  end subroutine base_hamiltonian_get_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_none(this, name)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_del_none)

    nullify(htrm)
    call base_hamiltonian__del__(this, trim(adjustl(name)), htrm)
    if(associated(htrm)) call term_intrf_del(htrm)
    nullify(htrm)

    POP_SUB(base_hamiltonian_del_none)
  end subroutine base_hamiltonian_del_none

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_term(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_term_t),       pointer        :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_del_term)

    nullify(that, htrm)
    call base_hamiltonian__del__(this, trim(adjustl(name)), htrm)
    if(associated(htrm))then
      call term_intrf_get(htrm, that)
      call term_intrf_del(htrm)
    end if
    nullify(htrm)

    POP_SUB(base_hamiltonian_del_term)
  end subroutine base_hamiltonian_del_term

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_potn(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_potential_t),  pointer        :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_del_potn)

    nullify(that, htrm)
    call base_hamiltonian__del__(this, trim(adjustl(name)), htrm)
    if(associated(htrm))then
      call term_intrf_get(htrm, that)
      call term_intrf_del(htrm)
    end if
    nullify(htrm)

    POP_SUB(base_hamiltonian_del_potn)
  end subroutine base_hamiltonian_del_potn

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_fnct(this, name, that)
    type(base_hamiltonian_t), intent(inout) :: this
    character(len=*),         intent(in)    :: name
    type(base_functional_t), pointer        :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_del_fnct)

    nullify(that, htrm)
    call base_hamiltonian__del__(this, trim(adjustl(name)), htrm)
    if(associated(htrm))then
      call term_intrf_get(htrm, that)
      call term_intrf_del(htrm)
    end if
    nullify(htrm)

    POP_SUB(base_hamiltonian_del_fnct)
  end subroutine base_hamiltonian_del_fnct

  ! ---------------------------------------------------------
  subroutine base_hamiltonian_del_hmlt(this, name, that)
    type(base_hamiltonian_t),  intent(inout) :: this
    character(len=*),          intent(in)    :: name
    type(base_hamiltonian_t), pointer        :: that

    type(term_intrf_t), pointer :: htrm

    PUSH_SUB(base_hamiltonian_del_hmlt)

    nullify(that, htrm)
    call base_hamiltonian__del__(this, trim(adjustl(name)), htrm)
    if(associated(htrm))then
      call term_intrf_get(htrm, that)
      call term_intrf_del(htrm)
    end if
    nullify(htrm)

    POP_SUB(base_hamiltonian_del_hmlt)
  end subroutine base_hamiltonian_del_hmlt

  ! ---------------------------------------------------------
  subroutine base_hamiltonian__copy__(this, that)
    type(base_hamiltonian_t), intent(inout) :: this
    type(base_hamiltonian_t), intent(in)    :: that

    type(term_dict_iterator_t)               :: iter
    character(len=BASE_HAMILTONIAN_NAME_LEN) :: name
    type(term_intrf_t),              pointer :: term
    type(refcount_t),                pointer :: rcnt

    PUSH_SUB(base_hamiltonian__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt, term)
    call base_hamiltonian__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_hamiltonian__init__(this, that)
      call refcount_del(this%rcnt)
      this%energy = that%energy
      call term_dict_init(iter, that%hdct)
      do
        nullify(term)
        call term_dict_next(iter, name, term)
        if(.not.associated(term)) exit
        call base_hamiltonian__set__(this, trim(adjustl(name)), term_intrf_new(source=term))
      end do
      call term_dict_end(iter)
      nullify(term)
      if(associated(that%sim)) call storage_copy(this%data, that%data)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_hamiltonian__copy__)
  end subroutine base_hamiltonian__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_hamiltonian__end__(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(term_intrf_t), pointer :: term

    PUSH_SUB(base_hamiltonian__end__)

    nullify(this%config, this%sys, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    this%nspin = default_nspin
    this%energy = 0.0_wp
    do
      nullify(term)
      call term_dict_pop(this%hdct, term)
      if(.not.associated(term))exit
      call term_intrf_del(term)
    end do
    nullify(term)
    call storage_end(this%data)
    ASSERT(term_dict_len(this%hdct)==0)
    call term_dict_end(this%hdct)
    ASSERT(base_hamiltonian_dict_len(this%dict)==0)
    call base_hamiltonian_dict_end(this%dict)

    POP_SUB(base_hamiltonian__end__)
  end subroutine base_hamiltonian__end__

end module base_hamiltonian_oct_m

!! Local Variables:
!! mode: f90
!! End:
