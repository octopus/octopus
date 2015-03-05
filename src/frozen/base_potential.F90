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

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_potential
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_potential

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_potential_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,   only: JSON_OK, json_object_t, json_get

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

  use storage_m, only:     &
    storage_t,             &
    storage_init,          &
    storage_start,         &
    storage_update,        &
    storage_stop,          &
    storage_reset,         &
    storage_accumulate,    &
    storage_eval,          &
    storage_get,           &
    storage_get_size,      &
    storage_get_dimension, &
    storage_copy,          &
    storage_end

  use storage_m, only: &
    storage_intrpl_t

  use storage_m, only:                             &
    BASE_POTENTIAL_INTRPL_OK => STORAGE_INTRPL_OK, &
    BASE_POTENTIAL_INTRPL_OD => STORAGE_INTRPL_OD, &
    BASE_POTENTIAL_INTRPL_NI => STORAGE_INTRPL_NI

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t

  implicit none

  private
  public ::                   &
    base_potential__init__,   &
    base_potential__start__,  &
    base_potential__update__, &
    base_potential__stop__,   &
    base_potential__reset__,  &
    base_potential__acc__,    &
    base_potential__add__,    &
    base_potential__copy__,   &
    base_potential__end__

  public ::                &
    base_potential_new,    &
    base_potential_del,    &
    base_potential_init,   &
    base_potential_start,  &
    base_potential_update, &
    base_potential_stop,   &
    base_potential_next,   &
    base_potential_eval,   &
    base_potential_set,    &
    base_potential_get,    &
    base_potential_copy,   &
    base_potential_end

  public ::                   &
    BASE_POTENTIAL_INTRPL_OK, &
    BASE_POTENTIAL_INTRPL_OD, &
    BASE_POTENTIAL_INTRPL_NI

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_potential_t
    private
    type(json_object_t),    pointer :: config =>null()
    type(base_system_t),    pointer :: sys    =>null()
    type(simulation_t),     pointer :: sim    =>null()
    type(base_potential_t), pointer :: prnt   =>null()
    real(kind=wp)                   :: energy = 0.0_wp
    type(storage_t)                 :: data
    type(config_dict_t)             :: dict
    type(base_potential_hash_t)     :: hash
    type(base_potential_list_t)     :: list
  end type base_potential_t

  type, public :: base_potential_iterator_t
    private
    type(base_potential_t),      pointer :: self =>null()
    type(base_potential_hash_iterator_t) :: iter
  end type base_potential_iterator_t

  type, public :: base_potential_intrpl_t
    private
    type(base_potential_t), pointer :: self =>null()
    type(storage_intrpl_t)          :: intrp
  end type base_potential_intrpl_t

  interface base_potential_init
    module procedure base_potential_init_potential
    module procedure base_potential_init_copy
    module procedure base_potential_iterator_init
    module procedure base_potential_intrpl_init
  end interface base_potential_init

  interface base_potential_set
    module procedure base_potential_set_energy
  end interface base_potential_set

  interface base_potential_get
    module procedure base_potential_get_potential
    module procedure base_potential_get_info
    module procedure base_potential_get_config
    module procedure base_potential_get_system
    module procedure base_potential_get_simulation
    module procedure base_potential_get_potential_1d
    module procedure base_potential_get_potential_md
    module procedure base_potential_intrpl_get
  end interface base_potential_get

  interface base_potential_next
    module procedure base_potential_iterator_next_config_potential
    module procedure base_potential_iterator_next_config
    module procedure base_potential_iterator_next_potential
  end interface base_potential_next

  interface base_potential_eval
    module procedure base_potential_intrpl_eval_1d
    module procedure base_potential_intrpl_eval_md
  end interface base_potential_eval

  interface base_potential_copy
    module procedure base_potential_copy_potential
    module procedure base_potential_iterator_copy
    module procedure base_potential_intrpl_copy
  end interface base_potential_copy

  interface base_potential_end
    module procedure base_potential_end_potential
    module procedure base_potential_iterator_end
    module procedure base_potential_intrpl_end
  end interface base_potential_end

  integer, public, parameter :: BASE_POTENTIAL_OK          = BASE_POTENTIAL_HASH_OK
  integer, public, parameter :: BASE_POTENTIAL_KEY_ERROR   = BASE_POTENTIAL_HASH_KEY_ERROR
  integer, public, parameter :: BASE_POTENTIAL_EMPTY_ERROR = BASE_POTENTIAL_HASH_EMPTY_ERROR

contains

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_potential_new(this, that)
    type(base_potential_t),  target, intent(inout) :: this
    type(base_potential_t), pointer                :: that
    !
    PUSH_SUB(base_potential_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call base_potential_list_push(this%list, that)
    POP_SUB(base_potential_new)
    return
  end subroutine base_potential_new

  ! ---------------------------------------------------------
  subroutine base_potential__idel__(this)
    type(base_potential_t), pointer :: this
    !
    PUSH_SUB(base_potential__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(base_potential__idel__)
    return
  end subroutine base_potential__idel__

  ! ---------------------------------------------------------
  subroutine base_potential_del(this)
    type(base_potential_t), pointer :: this
    !
    PUSH_SUB(base_potential_del)
    if(associated(this))then
      if(associated(this%prnt))then
        call base_potential_list_del(this%prnt%list, this)
        call base_potential_end(this)
        call base_potential__idel__(this)
      end if
    end if
    POP_SUB(base_potential_del)
    return
  end subroutine base_potential_del

  ! ---------------------------------------------------------
  subroutine base_potential__inull__(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential__inull__)
    nullify(this%config, this%sys, this%sim, this%prnt)
    this%energy=0.0_wp
    POP_SUB(base_potential__inull__)
    return
  end subroutine base_potential__inull__

  ! ---------------------------------------------------------
  subroutine base_potential__init__(this, sys, config)
    type(base_potential_t),      intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: nspin, ierr
    logical :: alloc
    !
    PUSH_SUB(base_potential__init__)
    call base_potential__inull__(this)
    this%config=>config
    this%sys=>sys
    call json_get(this%config, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=1
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK)alloc=.true.
    call storage_init(this%data, nspin, full=.false., allocate=alloc)
    call config_dict_init(this%dict)
    call base_potential_hash_init(this%hash)
    call base_potential_list_init(this%list)
    POP_SUB(base_potential__init__)
    return
  end subroutine base_potential__init__

  ! ---------------------------------------------------------
  subroutine base_potential_init_potential(this, sys, config)
    type(base_potential_t), intent(out) :: this
    type(base_system_t),    intent(in)  :: sys
    type(json_object_t),    intent(in)  :: config
    !
    PUSH_SUB(base_potential_init_potential)
    call base_potential__init__(this, sys, config)
    POP_SUB(base_potential_init_potential)
    return
  end subroutine base_potential_init_potential

  ! ---------------------------------------------------------
  subroutine base_potential_init_copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that
    !
    PUSH_SUB(base_potential_init_copy)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    ASSERT(associated(that%sim))
    call base_potential__init__(this, that%sys, that%config)
    call base_potential__start__(this, that%sim)
    POP_SUB(base_potential_init_copy)
    return
  end subroutine base_potential_init_copy

  ! ---------------------------------------------------------
  subroutine base_potential__start__(this, sim)
    type(base_potential_t),     intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_potential__start__)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call storage_start(this%data, sim)
    POP_SUB(base_potential__start__)
    return
  end subroutine base_potential__start__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_start(this, sim)
    type(base_potential_t), intent(inout) :: this
    type(simulation_t),     intent(in)    :: sim
    !
    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr
    !
    PUSH_SUB(base_potential_start)
    nullify(subs)
    call base_potential_init(iter, this)
    do
      nullify(subs)
      call base_potential_next(iter, subs, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_start(subs, sim)
    end do
    call base_potential_end(iter)
    nullify(subs)
    call base_potential__start__(this, sim)
    POP_SUB(base_potential_start)
    return
  end subroutine base_potential_start

  ! ---------------------------------------------------------
  subroutine base_potential__update__(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential__update__)
    ASSERT(associated(this%sim))
    call storage_update(this%data)
    POP_SUB(base_potential__update__)
    return
  end subroutine base_potential__update__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_update(this)
    type(base_potential_t), intent(inout) :: this
    !
    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr
    !
    PUSH_SUB(base_potential_update)
    nullify(subs)
    this%energy=0.0_wp
    call base_potential__reset__(this)
    call base_potential_init(iter, this)
    do
      nullify(subs)
      call base_potential_next(iter, subs, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_update(subs)
      call base_potential__acc__(this, subs)
      this%energy=this%energy+subs%energy
    end do
    call base_potential_end(iter)
    nullify(subs)
    call base_potential__update__(this)
    POP_SUB(base_potential_update)
    return
  end subroutine base_potential_update

  ! ---------------------------------------------------------
  subroutine base_potential__stop__(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential__stop__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)
    POP_SUB(base_potential__stop__)
    return
  end subroutine base_potential__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_stop(this)
    type(base_potential_t), intent(inout) :: this
    !
    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr
    !
    PUSH_SUB(base_potential_stop)
    nullify(subs)
    call base_potential_init(iter, this)
    do
      nullify(subs)
      call base_potential_next(iter, subs, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_stop(subs)
    end do
    call base_potential_end(iter)
    nullify(subs)
    call base_potential__stop__(this)
    POP_SUB(base_potential_stop)
    return
  end subroutine base_potential_stop

  ! ---------------------------------------------------------
  subroutine base_potential__reset__(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential__reset__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_reset(this%data)
    POP_SUB(base_potential__reset__)
    return
  end subroutine base_potential__reset__

  ! ---------------------------------------------------------
  subroutine base_potential__acc__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    !
    PUSH_SUB(base_potential__acc__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_accumulate(this%data, that%data)
    POP_SUB(base_potential__acc__)
    return
  end subroutine base_potential__acc__

  ! ---------------------------------------------------------
  subroutine base_potential__add__(this, that, config)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(base_potential__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_potential_hash_set(this%hash, config, that)
    POP_SUB(base_potential__add__)
    return
  end subroutine base_potential__add__

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential(this, name, that)
    type(base_potential_t),  intent(in) :: this
    character(len=*),        intent(in) :: name
    type(base_potential_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_potential_get_potential)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_potential_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_POTENTIAL_OK)nullify(that)
    end if
    POP_SUB(base_potential_get_potential)
    return
  end subroutine base_potential_get_potential

  ! ---------------------------------------------------------
  elemental subroutine base_potential_set_energy(this, that)
    type(base_potential_t), intent(inout) :: this
    real(kind=wp),          intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine base_potential_set_energy

  ! ---------------------------------------------------------
  elemental function base_potential_get_size(this) result(np)
    type(base_potential_t), intent(in) :: this
    !
    integer :: np
    !
    np=storage_get_size(this%data)
    return
  end function base_potential_get_size

  ! ---------------------------------------------------------
  elemental function base_potential_get_nspin(this) result(that)
    type(base_potential_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%data)
    return
  end function base_potential_get_nspin

  ! ---------------------------------------------------------
  elemental function base_potential_get_energy(this) result(that)
    type(base_potential_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function base_potential_get_energy

  ! ---------------------------------------------------------
  !elemental 
  subroutine base_potential_get_info(this, size, nspin, energy)
    type(base_potential_t),  intent(in)  :: this
    integer,       optional, intent(out) :: size
    integer,       optional, intent(out) :: nspin
    real(kind=wp), optional, intent(out) :: energy
    !
    if(present(size))&
      size=base_potential_get_size(this)
    if(present(nspin))&
      nspin=base_potential_get_nspin(this)
    if(present(energy))&
      energy=base_potential_get_energy(this)
    return
  end subroutine base_potential_get_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_config(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(json_object_t),   pointer             :: that
    !
    PUSH_SUB(base_potential_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_potential_get_config)
    return
  end subroutine base_potential_get_config

  ! ---------------------------------------------------------
  subroutine base_potential_get_system(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(base_system_t),   pointer             :: that
    !
    PUSH_SUB(base_potential_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(base_potential_get_system)
    return
  end subroutine base_potential_get_system

  ! ---------------------------------------------------------
  subroutine base_potential_get_simulation(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(simulation_t),    pointer             :: that
    !
    PUSH_SUB(base_potential_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_potential_get_simulation)
    return
  end subroutine base_potential_get_simulation

  ! ---------------------------------------------------------
  subroutine base_potential_get_storage(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(storage_t),       pointer             :: that
    !
    PUSH_SUB(base_potential_get_storage)
    that=>this%data
    POP_SUB(base_potential_get_storage)
    return
  end subroutine base_potential_get_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_1d(this, that)
    type(base_potential_t),       intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(base_potential_get_potential_1d)
    call storage_get(this%data, that)
    POP_SUB(base_potential_get_potential_1d)
    return
  end subroutine base_potential_get_potential_1d

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_md(this, that)
    type(base_potential_t),         intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(base_potential_get_potential_md)
    call storage_get(this%data, that)
    POP_SUB(base_potential_get_potential_md)
    return
  end subroutine base_potential_get_potential_md

  ! ---------------------------------------------------------
  subroutine base_potential__copy__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    !
    PUSH_SUB(base_potential__copy__)
    call base_potential__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_potential_init(this, that%sys, that%config)
      this%energy=that%energy
      if(associated(that%sim))then
        call base_potential_start(this, that%sim)
        call storage_copy(this%data, that%data)
      end if
    end if
    POP_SUB(base_potential__copy__)
    return
  end subroutine base_potential__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_copy_potential(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    !
    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: osub, isub
    type(json_object_t),    pointer :: cnfg
    integer                         :: ierr
    !
    PUSH_SUB(base_potential_copy_potential)
    nullify(cnfg, osub, isub)
    call base_potential_end(this)
    call base_potential__copy__(this, that)
    call base_potential_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_potential_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_new(this, osub)
      call base_potential_copy(osub, isub)
      call base_potential__add__(this, osub, cnfg)
    end do
    call base_potential_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(base_potential_copy_potential)
    return
  end subroutine base_potential_copy_potential

  ! ---------------------------------------------------------
  subroutine base_potential__end__(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential__end__)
    call base_potential__inull__(this)
    call storage_end(this%data)
    call config_dict_end(this%dict)
    call base_potential_hash_end(this%hash)
    call base_potential_list_end(this%list)
    POP_SUB(base_potential__end__)
    return
  end subroutine base_potential__end__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_end_potential(this)
    type(base_potential_t), intent(inout) :: this
    !
    type(base_potential_t), pointer :: subs
    !
    PUSH_SUB(base_potential_end_potential)
    do
      nullify(subs)
      call base_potential_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_potential_end(subs)
      call base_potential__idel__(subs)
    end do
    nullify(subs)
    call base_potential__end__(this)
    POP_SUB(base_potential_end_potential)
    return
  end subroutine base_potential_end_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_init(this, that)
    type(base_potential_iterator_t), intent(out) :: this
    type(base_potential_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_potential_iterator_init)
    This%self=>that
    call base_potential_hash_init(this%iter, that%hash)
    POP_SUB(base_potential_iterator_init)
    return
  end subroutine base_potential_iterator_init

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_config_potential(this, config, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(json_object_t),            pointer        :: config
    type(base_potential_t),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_config_potential)
    call base_potential_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_potential_iterator_next_config_potential)
    return
  end subroutine base_potential_iterator_next_config_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_config(this, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(json_object_t),            pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_config)
    call base_potential_hash_next(this%iter, that, ierr)
    POP_SUB(base_potential_iterator_next_config)
    return
  end subroutine base_potential_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_potential(this, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(base_potential_t),         pointer        :: that
    integer,               optional, intent(out)    :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_potential)
    call base_potential_hash_next(this%iter, that, ierr)
    POP_SUB(base_potential_iterator_next_potential)
    return
  end subroutine base_potential_iterator_next_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_copy(this, that)
    type(base_potential_iterator_t), intent(inout) :: this
    type(base_potential_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_potential_iterator_copy)
    this%self=>that%self
    call base_potential_hash_copy(this%iter, that%iter)
    POP_SUB(base_potential_iterator_copy)
    return
  end subroutine base_potential_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_end(this)
    type(base_potential_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_iterator_end)
    nullify(this%self)
    call base_potential_hash_end(this%iter)
    POP_SUB(base_potential_iterator_end)
    return
  end subroutine base_potential_iterator_end

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_init(this, that, type)
    type(base_potential_intrpl_t),  intent(out) :: this
    type(base_potential_t), target, intent(in)  :: that
    integer,              optional, intent(in)  :: type
    !
    PUSH_SUB(base_potential_intrpl_init)
    this%self=>that
    call storage_init(this%intrp, that%data, type)
    POP_SUB(base_potential_intrpl_init)
    return
  end subroutine base_potential_intrpl_init

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_eval_1d(this, x, v, ierr)
    type(base_potential_intrpl_t), intent(in)  :: this
    real(kind=wp),   dimension(:), intent(in)  :: x
    real(kind=wp),                 intent(out) :: v
    integer,                       intent(out) :: ierr
    !
    PUSH_SUB(base_potential_intrpl_eval_1d)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(base_potential_intrpl_eval_1d)
    return
  end subroutine base_potential_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_eval_md(this, x, v, ierr)
    type(base_potential_intrpl_t), intent(in)  :: this
    real(kind=wp),   dimension(:), intent(in)  :: x
    real(kind=wp),   dimension(:), intent(out) :: v
    integer,                       intent(out) :: ierr
    !
    PUSH_SUB(base_potential_intrpl_eval_md)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(base_potential_intrpl_eval_md)
    return
  end subroutine base_potential_intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_get(this, that)
    type(base_potential_intrpl_t), intent(in) :: this
    type(base_potential_t),       pointer     :: that
    !
    PUSH_SUB(base_potential_intrpl_get)
    nullify(that)
    if(associated(this%self))&
      that=>this%self
    POP_SUB(base_potential_intrpl_get)
    return
  end subroutine base_potential_intrpl_get

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_copy(this, that)
    type(base_potential_intrpl_t), intent(out) :: this
    type(base_potential_intrpl_t), intent(in)  :: that
    !
    PUSH_SUB(base_potential_intrpl_copy)
    call base_potential_intrpl_end(this)
    if(associated(that%self))then
      this%self=>that%self
      call storage_copy(this%intrp, that%intrp)
    end if
    POP_SUB(base_potential_intrpl_copy)
    return
  end subroutine base_potential_intrpl_copy

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_end(this)
    type(base_potential_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_intrpl_end)
    if(associated(this%self))&
      call storage_end(this%intrp)
    nullify(this%self)
    POP_SUB(base_potential_intrpl_end)
    return
  end subroutine base_potential_intrpl_end

end module base_potential_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
