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

#define HASH_TEMPLATE_NAME base_potential
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_potential

module base_potential_m

  use base_density_m
  use base_states_m
  use base_system_m
  use config_dict_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use storage_m

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_potential
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                     &
    BASE_POTENTIAL_OK,          &
    BASE_POTENTIAL_KEY_ERROR,   &
    BASE_POTENTIAL_EMPTY_ERROR

  public ::           &
    base_potential_t

  public ::                   &
    base_potential__init__,   &
    base_potential__start__,  &
    base_potential__update__, &
    base_potential__stop__,   &
    base_potential__reset__,  &
    base_potential__acc__,    &
    base_potential__sub__,    &
    base_potential__copy__,   &
    base_potential__end__

  public ::                &
    base_potential_new,    &
    base_potential_del,    &
    base_potential_init,   &
    base_potential_start,  &
    base_potential_update, &
    base_potential_stop,   &
    base_potential_sets,   &
    base_potential_gets,   &
    base_potential_set,    &
    base_potential_get,    &
    base_potential_copy,   &
    base_potential_end

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_POTENTIAL_OK          = BASE_POTENTIAL_HASH_OK
  integer, parameter :: BASE_POTENTIAL_KEY_ERROR   = BASE_POTENTIAL_HASH_KEY_ERROR
  integer, parameter :: BASE_POTENTIAL_EMPTY_ERROR = BASE_POTENTIAL_HASH_EMPTY_ERROR

  type :: base_potential_t
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

  interface base_potential__init__
    module procedure base_potential__init__type
    module procedure base_potential__init__copy
  end interface base_potential__init__

  interface base_potential_init
    module procedure base_potential_init_type
    module procedure base_potential_init_copy
  end interface base_potential_init

  interface base_potential_set
    module procedure base_potential_set_info
  end interface base_potential_set

  interface base_potential_gets
    module procedure base_potential_gets_config
    module procedure base_potential_gets_name
  end interface base_potential_gets

  interface base_potential_get
    module procedure base_potential_get_info
    module procedure base_potential_get_config
    module procedure base_potential_get_system
    module procedure base_potential_get_density
    module procedure base_potential_get_simulation
    module procedure base_potential_get_storage
    module procedure base_potential_get_potential_1d
    module procedure base_potential_get_potential_md
  end interface base_potential_get

  interface base_potential_copy
    module procedure base_potential_copy_type
  end interface base_potential_copy

  interface base_potential_end
    module procedure base_potential_end_type
  end interface base_potential_end

#define TEMPLATE_PREFIX base_potential
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_potential
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_potential_new(this, that)
    type(base_potential_t),  target, intent(inout) :: this
    type(base_potential_t), pointer                :: that

    PUSH_SUB(base_potential_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call base_potential_list_push(this%list, that)

    POP_SUB(base_potential_new)
  end subroutine base_potential_new

  ! ---------------------------------------------------------
  subroutine base_potential__idel__(this)
    type(base_potential_t), pointer :: this

    PUSH_SUB(base_potential__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_potential__idel__)
  end subroutine base_potential__idel__

  ! ---------------------------------------------------------
  subroutine base_potential_del(this)
    type(base_potential_t), pointer :: this

    PUSH_SUB(base_potential_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_potential_list_del(this%prnt%list, this)
        call base_potential_end(this)
        call base_potential__idel__(this)
      end if
    end if

    POP_SUB(base_potential_del)
  end subroutine base_potential_del

  ! ---------------------------------------------------------
  subroutine base_potential__init__type(this, sys, config)
    type(base_potential_t),      intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    integer :: nspin, ierr
    logical :: uspin, alloc

    PUSH_SUB(base_potential__init__type)

    this%config => config
    this%sys => sys
    nspin = 1
    call json_get(this%config, "spin", uspin, ierr)
    if(ierr/=JSON_OK) uspin = .false.
    if(uspin) call base_system_get(this%sys, nspin=nspin)
    ASSERT(nspin>0)
    ASSERT(nspin<3)
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK) alloc = .true.
    call storage_init(this%data, nspin, full=.false., allocate=alloc)
    call config_dict_init(this%dict)
    call base_potential_hash_init(this%hash)
    call base_potential_list_init(this%list)

    POP_SUB(base_potential__init__type)
  end subroutine base_potential__init__type

  ! ---------------------------------------------------------
  subroutine base_potential__init__copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(base_potential__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_potential__init__(this, that%sys, that%config)

    POP_SUB(base_potential__init__copy)
  end subroutine base_potential__init__copy

  ! ---------------------------------------------------------
  subroutine base_potential_init_type(this, sys, config)
    type(base_potential_t), intent(out) :: this
    type(base_system_t),    intent(in)  :: sys
    type(json_object_t),    intent(in)  :: config

    PUSH_SUB(base_potential_init_type)

    call base_potential__init__(this, sys, config)

    POP_SUB(base_potential_init_type)
  end subroutine base_potential_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_potential_init_copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: osub, isub
    type(json_object_t),    pointer :: cnfg
    integer                         :: ierr

    PUSH_SUB(base_potential_init_copy)

    nullify(cnfg, osub, isub)
    call base_potential__init__(this, that)
    call base_potential_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_potential_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_new(this, osub)
      call base_potential_init(osub, isub)
      call base_potential_sets(this, osub, cnfg)
    end do
    call base_potential_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_potential_init_copy)
  end subroutine base_potential_init_copy

  ! ---------------------------------------------------------
  subroutine base_potential__istart__(this, sim)
    type(base_potential_t),     intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_potential__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call storage_start(this%data, sim)

    POP_SUB(base_potential__istart__)
  end subroutine base_potential__istart__

  ! ---------------------------------------------------------
  subroutine base_potential__start__(this, sim)
    type(base_potential_t),       intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(base_potential__start__)

    if(present(sim))then
      call base_potential__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sim))
        call base_potential__istart__(this, this%prnt%sim)
      end if
    end if

    POP_SUB(base_potential__start__)
  end subroutine base_potential__start__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_start(this, sim)
    type(base_potential_t),       intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr

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
  end subroutine base_potential_start

  ! ---------------------------------------------------------
  subroutine base_potential__update__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)

    POP_SUB(base_potential__update__)
  end subroutine base_potential__update__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_update(this)
    type(base_potential_t), intent(inout) :: this

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr

    PUSH_SUB(base_potential_update)

    nullify(subs)
    call base_potential_init(iter, this)
    do
      nullify(subs)
      call base_potential_next(iter, subs, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call base_potential_update(subs)
    end do
    call base_potential_end(iter)
    nullify(subs)
    call base_potential__update__(this)

    POP_SUB(base_potential_update)
  end subroutine base_potential_update

  ! ---------------------------------------------------------
  subroutine base_potential__stop__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)

    POP_SUB(base_potential__stop__)
  end subroutine base_potential__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_stop(this)
    type(base_potential_t), intent(inout) :: this

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                         :: ierr

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
  end subroutine base_potential_stop

  ! ---------------------------------------------------------
  subroutine base_potential__reset__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy=0.0_wp
    call storage_reset(this%data)

    POP_SUB(base_potential__reset__)
  end subroutine base_potential__reset__

  ! ---------------------------------------------------------
  subroutine base_potential__acc__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(base_potential__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = this%energy + that%energy
    call storage_add(this%data, that%data)

    POP_SUB(base_potential__acc__)
  end subroutine base_potential__acc__

  ! ---------------------------------------------------------
  subroutine base_potential__sub__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(base_potential__sub__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%energy = this%energy - that%energy
    call storage_sub(this%data, that%data)

    POP_SUB(base_potential__sub__)
  end subroutine base_potential__sub__

  ! ---------------------------------------------------------
  subroutine base_potential_sets(this, that, config)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_potential_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_potential_hash_set(this%hash, config, that)

    POP_SUB(base_potential_sets)
  end subroutine base_potential_sets

  ! ---------------------------------------------------------
  subroutine base_potential_gets_config(this, config, that)
    type(base_potential_t),  intent(in) :: this
    type(json_object_t),     intent(in) :: config
    type(base_potential_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_potential_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_potential_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_POTENTIAL_OK) nullify(that)

    POP_SUB(base_potential_gets_config)
  end subroutine base_potential_gets_config

  ! ---------------------------------------------------------
  subroutine base_potential_gets_name(this, name, that)
    type(base_potential_t),  intent(in) :: this
    character(len=*),        intent(in) :: name
    type(base_potential_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_potential_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_potential_gets(this, config, that)

    POP_SUB(base_potential_gets_name)
  end subroutine base_potential_gets_name

  ! ---------------------------------------------------------
  subroutine base_potential_set_info(this, energy)
    type(base_potential_t),  intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_potential_set_info)

    if(present(energy)) this%energy = energy

    POP_SUB(base_potential_set_info)
  end subroutine base_potential_set_info
 
  ! ---------------------------------------------------------
  subroutine base_potential_get_info(this, size, nspin, energy, use)
    type(base_potential_t),  intent(in)  :: this
    integer,       optional, intent(out) :: size
    integer,       optional, intent(out) :: nspin
    real(kind=wp), optional, intent(out) :: energy
    logical,       optional, intent(out) :: use

    PUSH_SUB(base_potential_get_info)

    if(present(energy)) energy = this%energy
    call storage_get(this%data, size=size, dim=nspin, alloc=use)

    POP_SUB(base_potential_get_info)
  end subroutine base_potential_get_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_config(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(json_object_t),   pointer             :: that

    PUSH_SUB(base_potential_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_potential_get_config)
  end subroutine base_potential_get_config

  ! ---------------------------------------------------------
  subroutine base_potential_get_system(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(base_system_t),   pointer             :: that

    PUSH_SUB(base_potential_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_potential_get_system)
  end subroutine base_potential_get_system

  ! ---------------------------------------------------------
  subroutine base_potential_get_density(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(base_density_t),  pointer             :: that

    type(base_states_t), pointer :: st

    PUSH_SUB(base_potential_get_)

    nullify(that)
    if(associated(this%sys))then
      call base_system_get(this%sys, st)
      if(associated(st)) call base_states_get(st, that)
    end if

    POP_SUB(base_potential_get_)
  end subroutine base_potential_get_density

  ! ---------------------------------------------------------
  subroutine base_potential_get_simulation(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(simulation_t),    pointer             :: that

    PUSH_SUB(base_potential_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_potential_get_simulation)
  end subroutine base_potential_get_simulation

  ! ---------------------------------------------------------
  subroutine base_potential_get_storage(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(storage_t),       pointer             :: that

    PUSH_SUB(base_potential_get_storage)

    that => this%data

    POP_SUB(base_potential_get_storage)
  end subroutine base_potential_get_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_1d(this, that)
    type(base_potential_t),       intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_potential_get_potential_1d)

    call storage_get(this%data, that)

    POP_SUB(base_potential_get_potential_1d)
  end subroutine base_potential_get_potential_1d

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_md(this, that)
    type(base_potential_t),         intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(base_potential_get_potential_md)

    call storage_get(this%data, that)

    POP_SUB(base_potential_get_potential_md)
  end subroutine base_potential_get_potential_md

  ! ---------------------------------------------------------
  subroutine base_potential__copy__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(base_potential__copy__)

    call base_potential__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_potential__init__(this, that%sys, that%config)
      this%energy=that%energy
      if(associated(that%sim))then
        call base_potential__start__(this, that%sim)
        call storage_copy(this%data, that%data)
      end if
    end if

    POP_SUB(base_potential__copy__)
  end subroutine base_potential__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_copy_type(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: osub, isub
    type(json_object_t),    pointer :: cnfg
    integer                         :: ierr

    PUSH_SUB(base_potential_copy_type)

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
      call base_potential_sets(this, osub, cnfg)
    end do
    call base_potential_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_potential_copy_type)
  end subroutine base_potential_copy_type

  ! ---------------------------------------------------------
  subroutine base_potential__end__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__end__)

    nullify(this%config, this%sys, this%sim, this%prnt)
    this%energy = 0.0_wp
    call storage_end(this%data)
    call config_dict_end(this%dict)
    call base_potential_hash_end(this%hash)
    call base_potential_list_end(this%list)

    POP_SUB(base_potential__end__)
  end subroutine base_potential__end__

  ! ---------------------------------------------------------
  recursive subroutine base_potential_end_type(this)
    type(base_potential_t), intent(inout) :: this

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_end_type)

    do
      nullify(subs)
      call base_potential_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_potential_end(subs)
      call base_potential__idel__(subs)
    end do
    nullify(subs)
    call base_potential__end__(this)

    POP_SUB(base_potential_end_type)
  end subroutine base_potential_end_type

#define TEMPLATE_PREFIX base_potential
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_potential_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
