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

#define HASH_TEMPLATE_NAME base_density
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_density

module base_density_m

  use config_dict_m    
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use storage_m

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_density
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                   &
    BASE_DENSITY_OK,          &
    BASE_DENSITY_KEY_ERROR,   &
    BASE_DENSITY_EMPTY_ERROR

  public ::         &
    base_density_t

  public ::                 &
    base_density__init__,   &
    base_density__start__,  &
    base_density__update__, &
    base_density__stop__,   &
    base_density__reset__,  &
    base_density__acc__,    &
    base_density__copy__,   &
    base_density__end__

  public ::              &
    base_density_new,    &
    base_density_del,    &
    base_density_init,   &
    base_density_start,  &
    base_density_update, &
    base_density_stop,   &
    base_density_sets,   &
    base_density_gets,   &
    base_density_set,    &
    base_density_get,    &
    base_density_copy,   &
    base_density_end

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_DENSITY_OK          = BASE_DENSITY_HASH_OK
  integer, parameter :: BASE_DENSITY_KEY_ERROR   = BASE_DENSITY_HASH_KEY_ERROR
  integer, parameter :: BASE_DENSITY_EMPTY_ERROR = BASE_DENSITY_HASH_EMPTY_ERROR

  integer, parameter :: default_nspin = 1

  type :: base_density_t
    private
    type(json_object_t),             pointer :: config =>null()
    type(simulation_t),              pointer :: sim    =>null()
    type(storage_t),                 pointer :: total  =>null()
    type(base_density_t),            pointer :: prnt   =>null()
    integer                                  :: nspin  = 0
    real(kind=wp), dimension(:), allocatable :: charge
    type(storage_t)                          :: data
    type(config_dict_t)                      :: dict
    type(base_density_hash_t)                :: hash
    type(base_density_list_t)                :: list
  end type base_density_t

  interface base_density__init__
    module procedure base_density__init__type
    module procedure base_density__init__copy
  end interface base_density__init__

  interface base_density_init
    module procedure base_density_init_type
    module procedure base_density_init_copy
  end interface base_density_init

  interface base_density_set
    module procedure base_density_set_charge
  end interface base_density_set

  interface base_density_gets
    module procedure base_density_gets_config
    module procedure base_density_gets_name
  end interface base_density_gets

  interface base_density_get
    module procedure base_density_get_info
    module procedure base_density_get_charge
    module procedure base_density_get_config
    module procedure base_density_get_simulation
    module procedure base_density_get_storage
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
  end interface base_density_get

  interface base_density_copy
    module procedure base_density_copy_type
  end interface base_density_copy

  interface base_density_end
    module procedure base_density_end_type
  end interface base_density_end

#define TEMPLATE_PREFIX base_density
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_density_new(this, that)
    type(base_density_t),  target, intent(inout) :: this
    type(base_density_t), pointer                :: that

    PUSH_SUB(base_density_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call base_density_list_push(this%list, that)

    POP_SUB(base_density_new)
  end subroutine base_density_new

  ! ---------------------------------------------------------
  subroutine base_density__idel__(this)
    type(base_density_t), pointer :: this

    PUSH_SUB(base_density__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_density__idel__)
  end subroutine base_density__idel__

  ! ---------------------------------------------------------
  subroutine base_density_del(this)
    type(base_density_t), pointer :: this

    PUSH_SUB(base_density_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_density_list_del(this%prnt%list, this)
        call base_density_end(this)
        call base_density__idel__(this)
      end if
    end if

    POP_SUB(base_density_del)
  end subroutine base_density_del

  ! ---------------------------------------------------------
  subroutine base_density__init__type(this, config)
    type(base_density_t), target, intent(out) :: this
    type(json_object_t),  target, intent(in)  :: config

    integer :: ierr
    logical :: alloc

    PUSH_SUB(base_density__init__type)

    this%config => config
    call json_get(this%config, "nspin", this%nspin, ierr)
    if(ierr/=JSON_OK) this%nspin = default_nspin
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    SAFE_ALLOCATE(this%charge(1:this%nspin))
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK) this%charge = 0.0_wp
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK) alloc = .true.
    if(this%nspin>1)then
      SAFE_ALLOCATE(this%total)
      call storage_init(this%total, allocate=alloc)
    else
      this%total => this%data
    end if
    call storage_init(this%data, this%nspin, allocate=alloc)
    call config_dict_init(this%dict)
    call base_density_hash_init(this%hash)
    call base_density_list_init(this%list)

    POP_SUB(base_density__init__type)
  end subroutine base_density__init__type

  ! ---------------------------------------------------------
  subroutine base_density__init__copy(this, that)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that

    PUSH_SUB(base_density__init__copy)

    ASSERT(associated(that%config))
    call base_density__init__(this, that%config)

    POP_SUB(base_density__init__copy)
  end subroutine base_density__init__copy

  ! ---------------------------------------------------------
  subroutine base_density_init_type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(base_density_init_type)

    call base_density__init__(this, config)

    POP_SUB(base_density_init_type)
  end subroutine base_density_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_density_init_copy(this, that)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: osub, isub
    type(json_object_t),  pointer :: cnfg
    integer                       :: ierr

    PUSH_SUB(base_density_init_copy)

    nullify(cnfg, osub, isub)
    call base_density__init__(this, that)
    call base_density_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_density_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density_new(this, osub)
      call base_density_init(osub, isub)
      call base_density_sets(this, osub, cnfg)
    end do
    call base_density_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_density_init_copy)
  end subroutine base_density_init_copy

  ! ---------------------------------------------------------
  subroutine base_density__istart__(this, sim)
    type(base_density_t),       intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_density__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    if(this%nspin>1) call storage_start(this%total, this%sim, fine=.true.)
    call storage_start(this%data, this%sim, fine=.true.)

    POP_SUB(base_density__istart__)
  end subroutine base_density__istart__

  ! ---------------------------------------------------------
  subroutine base_density__start__(this, sim)
    type(base_density_t),         intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(base_density__start__)

    if(present(sim))then
      call base_density__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sim))
        call base_density__istart__(this, this%prnt%sim)
      end if
    end if

    POP_SUB(base_density__start__)
  end subroutine base_density__start__

  ! ---------------------------------------------------------
  recursive subroutine base_density_start(this, sim)
    type(base_density_t),         intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs
    integer                       :: ierr

    PUSH_SUB(base_density_start)

    nullify(subs)
    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, subs, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density_start(subs, sim)
    end do
    call base_density_end(iter)
    nullify(subs)
    call base_density__start__(this, sim)

    POP_SUB(base_density_start)
  end subroutine base_density_start

  ! ---------------------------------------------------------
  subroutine base_density__norm__(this)
    type(base_density_t), intent(inout) :: this

    real(kind=wp) :: intg
    integer       :: ispn

    PUSH_SUB(base_density__norm__)

    do ispn = 1, this%nspin
      if(this%charge(ispn)>0.0_wp)then
        call storage_integrate(this%data, ispn, intg)
        ASSERT(.not.(intg<0.0_wp))
        if(intg>0.0_wp)call storage_mlt(this%data, ispn, this%charge(ispn)/intg)
      end if
    end do

    POP_SUB(base_density__norm__)
  end subroutine base_density__norm__

  ! ---------------------------------------------------------
  subroutine base_density__update__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_update(this%data)
    call base_density__norm__(this)
    if(this%nspin>1) call storage_reduce(this%total, this%data)

    POP_SUB(base_density__update__)
  end subroutine base_density__update__

  ! ---------------------------------------------------------
  recursive subroutine base_density_update(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs
    integer                       :: ierr

    PUSH_SUB(base_density_update)

    nullify(subs)
    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, subs, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density_update(subs)
    end do
    call base_density_end(iter)
    nullify(subs)
    call base_density__update__(this)

    POP_SUB(base_density_update)
  end subroutine base_density_update

  ! ---------------------------------------------------------
  subroutine base_density__stop__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%nspin>1) call storage_stop(this%total)
    call storage_stop(this%data)

    POP_SUB(base_density__stop__)
  end subroutine base_density__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_density_stop(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs
    integer                       :: ierr

    PUSH_SUB(base_density_stop)

    nullify(subs)
    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, subs, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density_stop(subs)
    end do
    call base_density_end(iter)
    nullify(subs)
    call base_density__stop__(this)

    POP_SUB(base_density_stop)
  end subroutine base_density_stop

  ! ---------------------------------------------------------
  subroutine base_density__reset__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%charge = 0.0_wp
    if(this%nspin>1) call storage_reset(this%total)
    call storage_reset(this%data)

    POP_SUB(base_density__reset__)
  end subroutine base_density__reset__

  ! ---------------------------------------------------------
  subroutine base_density__acc__(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(this%nspin==that%nspin)
    this%charge(:) = this%charge(:) + that%charge(:)
    call storage_add(this%data, that%data)
    call base_density__update__(this)

    POP_SUB(base_density__acc__)
  end subroutine base_density__acc__

  ! ---------------------------------------------------------
  subroutine base_density_sets(this, that, config)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that
    type(json_object_t),  intent(in)    :: config

    character(len=BASE_DENSITY_NAME_LEN) :: name
    integer                              :: ierr

    PUSH_SUB(base_density_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_density_hash_set(this%hash, config, that)

    POP_SUB(base_density_sets)
  end subroutine base_density_sets

  ! ---------------------------------------------------------
  subroutine base_density_gets_config(this, config, that)
    type(base_density_t),  intent(in) :: this
    type(json_object_t),   intent(in) :: config
    type(base_density_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_density_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_density_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_DENSITY_OK) nullify(that)

    POP_SUB(base_density_gets_config)
  end subroutine base_density_gets_config

  ! ---------------------------------------------------------
  subroutine base_density_gets_name(this, name, that)
    type(base_density_t),  intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_density_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_density_gets(this, config, that)

    POP_SUB(base_density_gets_name)
  end subroutine base_density_gets_name

  ! ---------------------------------------------------------
  subroutine base_density_set_charge(this, charge, spin)
    type(base_density_t), intent(inout) :: this
    real(kind=wp),        intent(in)    :: charge
    integer,    optional, intent(in)    :: spin

    PUSH_SUB(base_density_set_charge)

    if(present(spin))then
      this%charge(spin) = charge
    else
      this%charge = charge / real(this%nspin, kind=wp)
    end if

    POP_SUB(base_density_set_charge)
  end subroutine base_density_set_charge

  ! ---------------------------------------------------------
  subroutine base_density_get_info(this, size, nspin, fine, use)
    type(base_density_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin
    logical,    optional, intent(out) :: fine
    logical,    optional, intent(out) :: use

    PUSH_SUB(base_density_get_info)

    call storage_get(this%data, size=size, dim=nspin, fine=fine, alloc=use)

    POP_SUB(base_density_get_info)
  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_charge(this, charge, spin)
    type(base_density_t), intent(in)  :: this
    real(kind=wp),        intent(out) :: charge
    integer,    optional, intent(in)  :: spin

    PUSH_SUB(base_density_get_charge)

    if(present(spin))then
      charge = this%charge(spin)
    else
      charge = sum(this%charge)
    end if

    POP_SUB(base_density_get_charge)
  end subroutine base_density_get_charge

  ! ---------------------------------------------------------
  subroutine base_density_get_config(this, that)
    type(base_density_t), target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_density_get_config)

    nullify(that)
    if(associated(this%config)) that=>this%config

    POP_SUB(base_density_get_config)
  end subroutine base_density_get_config

  ! ---------------------------------------------------------
  subroutine base_density_get_simulation(this, that)
    type(base_density_t), target, intent(in) :: this
    type(simulation_t),  pointer             :: that

    PUSH_SUB(base_density_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_density_get_simulation)
  end subroutine base_density_get_simulation

  ! ---------------------------------------------------------
  subroutine base_density_get_storage(this, that, total)
    type(base_density_t), target, intent(in) :: this
    type(storage_t),     pointer             :: that
    logical,            optional, intent(in) :: total

    logical :: itotal

    PUSH_SUB(base_density_get_storage)

    nullify(that)
    itotal = .false.
    if(present(total)) itotal = total
    if(itotal)then
      that => this%total
    else
      that => this%data
    end if

    POP_SUB(base_density_get_storage)
  end subroutine base_density_get_storage

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, total)
    type(base_density_t),         intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    logical,            optional, intent(in) :: total

    logical :: itotal

    PUSH_SUB(base_density_get_density_1d)

    nullify(that)
    itotal = .false.
    if(present(total)) itotal = total
    if(itotal)then
      call storage_get(this%total, that)
    else
      call storage_get(this%data, that)
    end if

    POP_SUB(base_density_get_density_1d)
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that, total)
    type(base_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    logical :: itotal

    PUSH_SUB(density_get_base_density_2d)

    nullify(that)
    itotal = .false.
    if(present(total)) itotal = total
    if(itotal)then
      call storage_get(this%total, that)
    else
      call storage_get(this%data, that)
    end if

    POP_SUB(base_density_get_density_2d)
  end subroutine base_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_density__copy__(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density__copy__)

    call base_density__end__(this)
    if(associated(that%config))then
      call base_density__init__(this, that%config)
      this%charge(:) = that%charge(:)
      if(associated(that%sim)) then
        call base_density__start__(this, that%sim)
        call storage_copy(this%data, that%data)
        call base_density__update__(this)
      end if
    end if

    POP_SUB(base_density__copy__)
  end subroutine base_density__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_density_copy_type(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: osub, isub
    type(json_object_t),  pointer :: cnfg
    integer                       :: ierr

    PUSH_SUB(base_density_copy_type)

    nullify(cnfg, osub, isub)
    call base_density_end(this)
    call base_density__copy__(this, that)
    call base_density_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_density_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call base_density_new(this, osub)
      call base_density_copy(osub, isub)
      call base_density_sets(this, osub, cnfg)
    end do
    call base_density_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_density_copy_type)
  end subroutine base_density_copy_type

  ! ---------------------------------------------------------
  subroutine base_density__end__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__end__)

    if(this%nspin>1)then
      call storage_end(this%total)
      SAFE_DEALLOCATE_P(this%total)
    end if
    call storage_end(this%data)
    nullify(this%config, this%sim, this%total, this%prnt)
    this%nspin = 0
    SAFE_DEALLOCATE_A(this%charge)
    call config_dict_end(this%dict)
    call base_density_hash_end(this%hash)
    call base_density_list_end(this%list)

    POP_SUB(base_density__end__)
  end subroutine base_density__end__

  ! ---------------------------------------------------------
  recursive subroutine base_density_end_type(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_end_type)

    do
      nullify(subs)
      call base_density_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_density_end(subs)
      call base_density__idel__(subs)
    end do
    nullify(subs)
    call base_density__end__(this)

    POP_SUB(base_density_end_type)
  end subroutine base_density_end_type

#define TEMPLATE_PREFIX base_density
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_density_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

