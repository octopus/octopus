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

#define HASH_TEMPLATE_NAME base_functional
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_functional

module base_functional_m

  use base_density_m
  use base_states_m
  use base_system_m
  use config_dict_m
  use functional_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use storage_m

#define LIST_TEMPLATE_NAME base_functional
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_functional
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                      &
    BASE_FUNCTIONAL_OK,          &
    BASE_FUNCTIONAL_KEY_ERROR,   &
    BASE_FUNCTIONAL_EMPTY_ERROR

  public ::            &
    base_functional_t

  public ::                    &
    base_functional__init__,   &
    base_functional__start__,  &
    base_functional__update__, &
    base_functional__stop__,   &
    base_functional__calc__,   &
    base_functional__reset__,  &
    base_functional__acc__,    &
    base_functional__sub__,    &
    base_functional__copy__,   &
    base_functional__end__

  public ::                 &
    base_functional_new,    &
    base_functional_del,    &
    base_functional_init,   &
    base_functional_start,  &
    base_functional_update, &
    base_functional_stop,   &
    base_functional_calc,   &
    base_functional_sets,   &
    base_functional_gets,   &
    base_functional_set,    &
    base_functional_get,    &
    base_functional_copy,   &
    base_functional_end

#define LIST_TEMPLATE_NAME base_functional
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_FUNCTIONAL_OK          = BASE_FUNCTIONAL_HASH_OK
  integer, parameter :: BASE_FUNCTIONAL_KEY_ERROR   = BASE_FUNCTIONAL_HASH_KEY_ERROR
  integer, parameter :: BASE_FUNCTIONAL_EMPTY_ERROR = BASE_FUNCTIONAL_HASH_EMPTY_ERROR

  type :: base_functional_raii_t
    private
    type(base_functional_t), pointer :: prnt =>null()
    type(base_functional_list_t)     :: list
  end type base_functional_raii_t

  type :: base_functional_t
    private
    type(json_object_t),  pointer :: config  =>null()
    type(base_system_t),  pointer :: sys     =>null()
    type(simulation_t),   pointer :: sim     =>null()
    real(kind=wp)                 :: factor  = 1.0_wp
    real(kind=wp)                 :: energy  = 0.0_wp
    type(functional_t)            :: funct
    type(storage_t)               :: data
    type(config_dict_t)           :: dict
    type(base_functional_hash_t)  :: hash
    type(base_functional_raii_t)  :: raii
  end type base_functional_t

  interface base_functional__init__
    module procedure base_functional__init__type
    module procedure base_functional__init__copy
  end interface base_functional__init__

  interface base_functional_init
    module procedure base_functional_init_type
    module procedure base_functional_init_copy
  end interface base_functional_init

  interface base_functional_set
    module procedure base_functional_set_info
  end interface base_functional_set

  interface base_functional_gets
    module procedure base_functional_gets_config
    module procedure base_functional_gets_name
  end interface base_functional_gets

  interface base_functional_get
    module procedure base_functional_get_info
    module procedure base_functional_get_config
    module procedure base_functional_get_system
    module procedure base_functional_get_density
    module procedure base_functional_get_simulation
    module procedure base_functional_get_storage
    module procedure base_functional_get_functional_1d
    module procedure base_functional_get_functional_md
  end interface base_functional_get

  interface base_functional_copy
    module procedure base_functional_copy_type
  end interface base_functional_copy

  interface base_functional_end
    module procedure base_functional_end_type
  end interface base_functional_end

#define TEMPLATE_PREFIX base_functional
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_functional
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_functional_new(this, that)
    type(base_functional_t),  target, intent(inout) :: this
    type(base_functional_t), pointer                :: that

    PUSH_SUB(base_functional_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%raii%prnt => this
    call base_functional_list_push(this%raii%list, that)

    POP_SUB(base_functional_new)
  end subroutine base_functional_new

  ! ---------------------------------------------------------
  subroutine base_functional__idel__(this)
    type(base_functional_t), pointer :: this

    PUSH_SUB(base_functional__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_functional__idel__)
  end subroutine base_functional__idel__

  ! ---------------------------------------------------------
  subroutine base_functional_del(this)
    type(base_functional_t), pointer :: this

    PUSH_SUB(base_functional_del)

    if(associated(this))then
      if(associated(this%raii%prnt))then
        call base_functional_list_del(this%raii%prnt%raii%list, this)
        call base_functional_end(this)
        call base_functional__idel__(this)
      end if
    end if

    POP_SUB(base_functional_del)
  end subroutine base_functional_del

  ! ---------------------------------------------------------
  subroutine base_functional__init__type(this, sys, config)
    type(base_functional_t),     intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    integer :: id, nspin, ierr

    PUSH_SUB(base_functional__init__type)

    this%config => config
    this%sys => sys
    call json_get(config, "factor", this%factor, ierr)
    if(ierr/=JSON_OK) this%factor = 1.0_wp
    call json_get(config, "functional", id, ierr)
    if(ierr/=JSON_OK) id = FUNCT_XC_NONE
    call base_system_get(this%sys, nspin=nspin)
    call functional_init(this%funct, id, nspin)
    call storage_init(this%data, nspin, full=.false.)
    call config_dict_init(this%dict)
    call base_functional_hash_init(this%hash)
    call base_functional_list_init(this%raii%list)

    POP_SUB(base_functional__init__type)
  end subroutine base_functional__init__type

  ! ---------------------------------------------------------
  subroutine base_functional__init__copy(this, that)
    type(base_functional_t), intent(out) :: this
    type(base_functional_t), intent(in)  :: that

    PUSH_SUB(base_functional__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_functional__init__(this, that%sys, that%config)

    POP_SUB(base_functional__init__copy)
  end subroutine base_functional__init__copy

  ! ---------------------------------------------------------
  subroutine base_functional_init_type(this, sys, config)
    type(base_functional_t), intent(out) :: this
    type(base_system_t),     intent(in)  :: sys
    type(json_object_t),     intent(in)  :: config

    PUSH_SUB(base_functional_init_type)

    call base_functional__init__(this, sys, config)

    POP_SUB(base_functional_init_type)
  end subroutine base_functional_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_functional_init_copy(this, that)
    type(base_functional_t), intent(out) :: this
    type(base_functional_t), intent(in)  :: that

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: osub, isub
    type(json_object_t),     pointer :: cnfg
    integer                          :: ierr

    PUSH_SUB(base_functional_init_copy)

    nullify(cnfg, osub, isub)
    call base_functional__init__(this, that)
    call base_functional_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_functional_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_new(this, osub)
      call base_functional_init(osub, isub)
      call base_functional_sets(this, osub, cnfg)
    end do
    call base_functional_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_functional_init_copy)
  end subroutine base_functional_init_copy

  ! ---------------------------------------------------------
  subroutine base_functional__istart__(this, sim)
    type(base_functional_t),    intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_functional__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call functional_start(this%funct, sim, fine=.true.)
    call storage_start(this%data, sim)

    POP_SUB(base_functional__istart__)
  end subroutine base_functional__istart__

  ! ---------------------------------------------------------
  subroutine base_functional__start__(this, sim)
    type(base_functional_t),      intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(base_functional__start__)

    if(present(sim))then
      call base_functional__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%raii%prnt))
        ASSERT(associated(this%raii%prnt%sim))
        call base_functional__istart__(this, this%raii%prnt%sim)
      end if
    end if

    POP_SUB(base_functional__start__)
  end subroutine base_functional__start__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_start(this, sim)
    type(base_functional_t),      intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: subs
    integer                          :: ierr

    PUSH_SUB(base_functional_start)

    nullify(subs)
    call base_functional_init(iter, this)
    do
      nullify(subs)
      call base_functional_next(iter, subs, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_start(subs, sim)
    end do
    call base_functional_end(iter)
    nullify(subs)
    call base_functional__start__(this, sim)

    POP_SUB(base_functional_start)
  end subroutine base_functional_start

  ! ---------------------------------------------------------
  subroutine base_functional__update__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)

    POP_SUB(base_functional__update__)
  end subroutine base_functional__update__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_update(this)
    type(base_functional_t), intent(inout) :: this

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: subs
    integer                          :: ierr

    PUSH_SUB(base_functional_update)

    nullify(subs)
    call base_functional_init(iter, this)
    do
      nullify(subs)
      call base_functional_next(iter, subs, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_update(subs)
    end do
    call base_functional_end(iter)
    nullify(subs)
    call base_functional__update__(this)

    POP_SUB(base_functional_update)
  end subroutine base_functional_update

  ! ---------------------------------------------------------
  subroutine base_functional__stop__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)

    POP_SUB(base_functional__stop__)
  end subroutine base_functional__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_stop(this)
    type(base_functional_t), intent(inout) :: this

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: subs
    integer                         :: ierr

    PUSH_SUB(base_functional_stop)

    nullify(subs)
    call base_functional_init(iter, this)
    do
      nullify(subs)
      call base_functional_next(iter, subs, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_stop(subs)
    end do
    call base_functional_end(iter)
    nullify(subs)
    call base_functional__stop__(this)

    POP_SUB(base_functional_stop)
  end subroutine base_functional_stop

  ! ---------------------------------------------------------
  subroutine base_functional__calc__(this)
    type(base_functional_t), intent(inout) :: this

    real(kind=wp), dimension(:,:), pointer :: fptn, potn, dnst
    type(base_density_t),          pointer :: density
    type(storage_t)                        :: data
    integer                                :: kind
    logical                                :: fine

    PUSH_SUB(base_functional__calc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(fptn, potn, dnst, density)
    call base_functional__reset__(this)
    call base_functional_get(this, kind=kind)
    if(kind>FUNCT_XC_NONE)then
      call storage_get(this%data, potn)
      ASSERT(associated(potn))
      fptn=>potn
      call base_functional_get(this, density)
      ASSERT(associated(density))
      call base_density_get(density, fine=fine)
      if(fine)then
        call storage_init(data, this%data)
        call storage_start(data, this%sim, fine)
        call storage_get(data, fptn)
        ASSERT(associated(fptn))
      end if
      call base_density_get(density, dnst)
      ASSERT(associated(dnst))
      nullify(density)
      call functional_calc(this%funct, dnst, this%energy, fptn)
      if(fine)then
        call storage_transfer(this%data, data)
        call storage_end(data)
        nullify(fptn)
        fptn=>potn
      end if
      if(abs(this%factor-1.0_wp)>epsilon(this%factor))then
        this%energy=this%factor*this%energy
        fptn=this%factor*fptn
      end if
      call base_functional__update__(this)
    end if

    POP_SUB(base_functional__calc__)
  end subroutine base_functional__calc__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_calc(this)
    type(base_functional_t), intent(inout) :: this

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: subs
    integer                          :: ierr

    PUSH_SUB(base_functional_calc)

    nullify(subs)
    call base_functional_init(iter, this)
    do
      nullify(subs)
      call base_functional_next(iter, subs, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_calc(subs)
    end do
    call base_functional_end(iter)
    nullify(subs)
    call base_functional__calc__(this)

    POP_SUB(base_functional_calc)
  end subroutine base_functional_calc

  ! ---------------------------------------------------------
  subroutine base_functional__reset__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = 0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_functional__reset__)
  end subroutine base_functional__reset__

  ! ---------------------------------------------------------
  subroutine base_functional__acc__(this, that)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that

    PUSH_SUB(base_functional__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = this%energy + that%energy
    call storage_add(this%data, that%data)

    POP_SUB(base_functional__acc__)
  end subroutine base_functional__acc__

  ! ---------------------------------------------------------
  subroutine base_functional__sub__(this, that)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that

    PUSH_SUB(base_functional__sub__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = this%energy - that%energy
    call storage_sub(this%data, that%data)

    POP_SUB(base_functional__sub__)
  end subroutine base_functional__sub__

  ! ---------------------------------------------------------
  subroutine base_functional_sets(this, that, config)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_functional_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_functional_hash_set(this%hash, config, that)

    POP_SUB(base_functional_sets)
  end subroutine base_functional_sets

  ! ---------------------------------------------------------
  subroutine base_functional_gets_config(this, config, that)
    type(base_functional_t),  intent(in) :: this
    type(json_object_t),      intent(in) :: config
    type(base_functional_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_functional_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_functional_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_FUNCTIONAL_OK) nullify(that)

    POP_SUB(base_functional_gets_config)
  end subroutine base_functional_gets_config

  ! ---------------------------------------------------------
  subroutine base_functional_gets_name(this, name, that)
    type(base_functional_t),  intent(in) :: this
    character(len=*),         intent(in) :: name
    type(base_functional_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_functional_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_functional_gets(this, config, that)

    POP_SUB(base_functional_gets_name)
  end subroutine base_functional_gets_name

  ! ---------------------------------------------------------
  subroutine base_functional_set_info(this, energy)
    type(base_functional_t), intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_functional_set_info)

    if(present(energy)) this%energy = energy

    POP_SUB(base_functional_set_info)
  end subroutine base_functional_set_info

  ! ---------------------------------------------------------
  subroutine base_functional_get_info(this, id, family, kind, size, nspin, energy)
    type(base_functional_t), intent(in)  :: this
    integer,       optional, intent(out) :: id
    integer,       optional, intent(out) :: family
    integer,       optional, intent(out) :: kind
    integer,       optional, intent(out) :: size
    integer,       optional, intent(out) :: nspin
    real(kind=wp), optional, intent(out) :: energy

    PUSH_SUB(base_functional_get_info)

    call functional_get(this%funct, id=id, family=family, kind=kind)
    call storage_get(this%data, size=size, dim=nspin)
    if(present(energy)) energy = this%energy

    POP_SUB(base_functional_get_info)
  end subroutine base_functional_get_info

  ! ---------------------------------------------------------
  subroutine base_functional_get_config(this, that)
    type(base_functional_t), target, intent(in) :: this
    type(json_object_t),    pointer             :: that

    PUSH_SUB(base_functional_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_functional_get_config)
  end subroutine base_functional_get_config

  ! ---------------------------------------------------------
  subroutine base_functional_get_simulation(this, that)
    type(base_functional_t), target, intent(in) :: this
    type(simulation_t),     pointer             :: that

    PUSH_SUB(base_functional_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_functional_get_simulation)
  end subroutine base_functional_get_simulation

  ! ---------------------------------------------------------
  subroutine base_functional_get_storage(this, that)
    type(base_functional_t), target, intent(in) :: this
    type(storage_t),        pointer             :: that

    PUSH_SUB(base_functional_get_storage)

    that => this%data

    POP_SUB(base_functional_get_storage)
  end subroutine base_functional_get_storage

  ! ---------------------------------------------------------
  subroutine base_functional_get_system(this, that)
    type(base_functional_t), target, intent(in) :: this
    type(base_system_t),    pointer             :: that

    PUSH_SUB(base_functional_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_functional_get_system)
  end subroutine base_functional_get_system

  ! ---------------------------------------------------------
  subroutine base_functional_get_density(this, that)
    type(base_functional_t), intent(in) :: this
    type(base_density_t),   pointer     :: that

    type(base_states_t), pointer :: st

    PUSH_SUB(base_functional_get_density)

    nullify(that)
    if(associated(this%sys))then
      call base_system_get(this%sys, st)
      if(associated(st)) call base_states_get(st, that)
    end if

    POP_SUB(base_functional_get_density)
  end subroutine base_functional_get_density

  ! ---------------------------------------------------------
  subroutine base_functional_get_functional_1d(this, that)
    type(base_functional_t),      intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_functional_get_functional_1d)

    call storage_get(this%data, that)

    POP_SUB(base_functional_get_functional_1d)
  end subroutine base_functional_get_functional_1d

  ! ---------------------------------------------------------
  subroutine base_functional_get_functional_md(this, that)
    type(base_functional_t),        intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(base_functional_get_functional_md)

    call storage_get(this%data, that)

    POP_SUB(base_functional_get_functional_md)
  end subroutine base_functional_get_functional_md

  ! ---------------------------------------------------------
  subroutine base_functional__copy__(this, that)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that

    PUSH_SUB(base_functional__copy__)

    call base_functional__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_functional__init__(this, that%sys, that%config)
      this%energy = that%energy
      if(associated(that%sim))then
        call base_functional__start__(this, that%sim)
        call storage_copy(this%data, that%data)
      end if
    end if

    POP_SUB(base_functional__copy__)
  end subroutine base_functional__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_copy_type(this, that)
    type(base_functional_t), intent(inout) :: this
    type(base_functional_t), intent(in)    :: that

    type(base_functional_iterator_t) :: iter
    type(base_functional_t), pointer :: osub, isub
    type(json_object_t),    pointer :: cnfg
    integer                         :: ierr

    PUSH_SUB(base_functional_copy_type)

    nullify(cnfg, osub, isub)
    call base_functional_end(this)
    call base_functional__copy__(this, that)
    call base_functional_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_functional_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_FUNCTIONAL_OK)exit
      call base_functional_new(this, osub)
      call base_functional_copy(osub, isub)
      call base_functional_sets(this, osub, cnfg)
    end do
    call base_functional_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_functional_copy_type)
  end subroutine base_functional_copy_type

  ! ---------------------------------------------------------
  subroutine base_functional__end__(this)
    type(base_functional_t), intent(inout) :: this

    PUSH_SUB(base_functional__end__)

    nullify(this%config, this%sys, this%sim, this%raii%prnt)
    this%factor = 1.0_wp
    this%energy = 0.0_wp
    call functional_end(this%funct)
    call storage_end(this%data)
    call config_dict_end(this%dict)
    call base_functional_hash_end(this%hash)
    call base_functional_list_end(this%raii%list)

    POP_SUB(base_functional__end__)
  end subroutine base_functional__end__

  ! ---------------------------------------------------------
  recursive subroutine base_functional_end_type(this)
    type(base_functional_t), intent(inout) :: this

    type(base_functional_t), pointer :: subs

    PUSH_SUB(base_functional_end_type)

    do
      nullify(subs)
      call base_functional_list_pop(this%raii%list, subs)
      if(.not.associated(subs))exit
      call base_functional_end(subs)
      call base_functional__idel__(subs)
    end do
    nullify(subs)
    call base_functional__end__(this)

    POP_SUB(base_functional_end_type)
  end subroutine base_functional_end_type

#define TEMPLATE_PREFIX base_functional
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_functional_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
