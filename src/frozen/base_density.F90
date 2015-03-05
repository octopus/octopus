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

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_PREFIX
#include "tlist.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_TEMPLATE_NAME base_density
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_density

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_density_m
  use config_dict_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m
  use storage_m

  implicit none

  private
  public ::                 &
    base_density__init__,   &
    base_density__start__,  &
    base_density__update__, &
    base_density__stop__,   &
    base_density__reset__,  &
    base_density__acc__,    &
    base_density__add__,    &
    base_density__copy__,   &
    base_density__end__

  public ::              &
    base_density_new,    &
    base_density_del,    &
    base_density_init,   &
    base_density_start,  &
    base_density_update, &
    base_density_stop,   &
    base_density_next,   &
    base_density_eval,   &
    base_density_get,    &
    base_density_copy,   &
    base_density_end

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_HEADER
#include "tlist.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_density_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(storage_t),      pointer :: total  =>null()
    type(base_density_t), pointer :: prnt   =>null()
    type(storage_t)               :: data
    type(config_dict_t)           :: dict
    type(base_density_hash_t)     :: hash
    type(base_density_list_t)     :: list
  end type base_density_t

  type, public :: base_density_iterator_t
    private
    type(base_density_t),      pointer :: self =>null()
    type(base_density_hash_iterator_t) :: iter
  end type base_density_iterator_t

  type, public :: base_density_intrpl_t
    private
    type(base_density_t), pointer :: self =>null()
    type(storage_intrpl_t)        :: intrp
  end type base_density_intrpl_t
 
  interface base_density__init__
    module procedure base_density__init__density
    module procedure base_density_intrpl_init
  end interface base_density__init__

  interface base_density__end__
    module procedure base_density__end__density
    module procedure base_density_intrpl_end
  end interface base_density__end__

  interface base_density_init
    module procedure base_density_init_density
    module procedure base_density_iterator_init
    module procedure base_density_intrpl_init
  end interface base_density_init

  interface base_density_get
    module procedure base_density_get_density
    module procedure base_density_get_info
    module procedure base_density_get_config
    module procedure base_density_get_simulation
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
  end interface base_density_get

  interface base_density_next
    module procedure base_density_iterator_next_config_density
    module procedure base_density_iterator_next_config
    module procedure base_density_iterator_next_density
  end interface base_density_next

  interface base_density_eval
    module procedure base_density_intrpl_eval_1d
    module procedure base_density_intrpl_eval_2d
  end interface base_density_eval

  interface base_density_copy
    module procedure base_density_copy_density
    module procedure base_density_iterator_copy
    module procedure base_density_intrpl_copy
  end interface base_density_copy

  interface base_density_end
    module procedure base_density_end_density
    module procedure base_density_iterator_end
    module procedure base_density_intrpl_end
  end interface base_density_end

  integer, parameter :: default_nspin = 1

  integer, public, parameter :: BASE_DENSITY_OK          = BASE_DENSITY_HASH_OK
  integer, public, parameter :: BASE_DENSITY_KEY_ERROR   = BASE_DENSITY_HASH_KEY_ERROR
  integer, public, parameter :: BASE_DENSITY_EMPTY_ERROR = BASE_DENSITY_HASH_EMPTY_ERROR

contains

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_BODY
#include "tlist.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash.F90"
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
  subroutine base_density__inull__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__inull__)

    nullify(this%config, this%sim, this%total, this%prnt)

    POP_SUB(base_density__inull__)
  end subroutine base_density__inull__

  ! ---------------------------------------------------------
  subroutine base_density__init__density(this, config)
    type(base_density_t), target, intent(out) :: this
    type(json_object_t),  target, intent(in)  :: config

    integer :: nspin, ierr
    logical :: alloc

    PUSH_SUB(base_density__init__density)

    call base_density__inull__(this)
    this%config=>config
    call json_get(this%config, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=default_nspin
    ASSERT(nspin>0)
    call json_get(this%config, "allocate", alloc, ierr)
    if(ierr/=JSON_OK)alloc=.true.
    if(nspin>1)then
      SAFE_ALLOCATE(this%total)
      call storage_init(this%total, allocate=alloc)
    else
      this%total=>this%data
    end if
    call storage_init(this%data, nspin, allocate=alloc)
    call config_dict_init(this%dict)
    call base_density_hash_init(this%hash)
    call base_density_list_init(this%list)

    POP_SUB(base_density__init__density)
  end subroutine base_density__init__density

  ! ---------------------------------------------------------
  subroutine base_density_init_density(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(base_density_init_density)

    call base_density__init__(this, config)

    POP_SUB(base_density_init_density)
  end subroutine base_density_init_density

  ! ---------------------------------------------------------
  subroutine base_density__start__(this, sim)
    type(base_density_t),       intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_density__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    if(base_density_get_nspin(this)>1)&
      call storage_start(this%total, sim, fine=.true.)
    call storage_start(this%data, sim, fine=.true.)

    POP_SUB(base_density__start__)
  end subroutine base_density__start__

  ! ---------------------------------------------------------
  recursive subroutine base_density_start(this, sim)
    type(base_density_t), intent(inout) :: this
    type(simulation_t),   intent(in)    :: sim

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
  subroutine base_density__update__(this)
    type(base_density_t), intent(inout) :: this

    real(kind=wp), dimension(:,:), pointer :: prho
    real(kind=wp), dimension(:),   pointer :: trho
    integer                                :: indx

    PUSH_SUB(base_density__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(prho, trho)
    call storage_update(this%data)
    if(base_density_get_nspin(this)>1)then
      call storage_get(this%data, prho)
      call storage_get(this%total, trho)
      do indx = 1, base_density_get_size(this)
        trho(indx)=sum(prho(indx,:))
      end do
      call storage_update(this%total)
      nullify(prho, trho)
    end if

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
    if(base_density_get_nspin(this)>1)&
      call storage_stop(this%total)
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
    call storage_accumulate(this%data, that%data)

    POP_SUB(base_density__acc__)
  end subroutine base_density__acc__

  ! ---------------------------------------------------------
  subroutine base_density__add__(this, that, config)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that
    type(json_object_t),  intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_density__add__)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_density_hash_set(this%hash, config, that)

    POP_SUB(base_density__add__)
  end subroutine base_density__add__

  ! ---------------------------------------------------------
  subroutine base_density_get_density(this, name, that)
    type(base_density_t),  intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_density_get_density)

    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_density_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_DENSITY_OK)nullify(that)
    end if

    POP_SUB(base_density_get_density)
  end subroutine base_density_get_density

  ! ---------------------------------------------------------
  elemental function base_density_get_size(this) result(that)
    type(base_density_t), intent(in) :: this

    integer :: that

    that=storage_get_size(this%data)
  end function base_density_get_size

  ! ---------------------------------------------------------
  elemental function base_density_get_nspin(this) result(that)
    type(base_density_t), intent(in) :: this

    integer :: that

    that=storage_get_dimension(this%data)
  end function base_density_get_nspin

  ! ---------------------------------------------------------
  !elemental 
  subroutine base_density_get_info(this, size, nspin)
    type(base_density_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin

    if(present(size)) size = base_density_get_size(this)
    if(present(nspin)) nspin = base_density_get_nspin(this)

  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_config(this, that)
    type(base_density_t), target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_density_get_config)

    that=>null()
    if(associated(this%config)) that=>this%config

    POP_SUB(base_density_get_config)
  end subroutine base_density_get_config

  ! ---------------------------------------------------------
  subroutine base_density_get_simulation(this, that)
    type(base_density_t), target, intent(in) :: this
    type(simulation_t),  pointer             :: that

    PUSH_SUB(base_density_get_simulation)

    that=>null()
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_density_get_simulation)
  end subroutine base_density_get_simulation

  ! ---------------------------------------------------------
  subroutine base_density_get_storage(this, that)
    type(base_density_t), target, intent(in) :: this
    type(storage_t),     pointer             :: that

    PUSH_SUB(base_density_get_storage)

    that=>this%data

    POP_SUB(base_density_get_storage)
  end subroutine base_density_get_storage

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, total)
    type(base_density_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(base_density_get_density_1d)

    call storage_get(this%data, that)
    if(present(total))call base_density_get_total_density(this, total)

    POP_SUB(base_density_get_density_1d)
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that)
    type(base_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_base_density_2d)

    call storage_get(this%data, that)

    POP_SUB(base_density_get_density_2d)
  end subroutine base_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_density_get_total_density(this, that)
    type(base_density_t),         intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(base_density_get_total_density)

    call storage_get(this%total, that)

    POP_SUB(base_density_get_total_density)
  end subroutine base_density_get_total_density

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_1_n(this, that)
    real(kind=wp),               intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    PUSH_SUB(base_density_adjust_spin_1_n)

    select case(size(that))
    case(1)
      this=that(1)
    case(2)
      this=sum(that)
    case(4)
      ASSERT(.false.)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(base_density_adjust_spin_1_n)
  end subroutine base_density_adjust_spin_1_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_2_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    PUSH_SUB(base_density_adjust_spin_2_n)

    select case(size(that))
    case(1)
      this=0.5_wp*that(1)
    case(2)
      this=that
    case(4)
      ASSERT(.false.)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(base_density_adjust_spin_2_n)
  end subroutine base_density_adjust_spin_2_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_4_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    PUSH_SUB(base_density_adjust_spin_4_n)

    select case(size(that))
    case(1)
      ASSERT(.false.)
    case(2)
      ASSERT(.false.)
    case(4)
      this=that
    case default
      ASSERT(.false.)
    end select

    POP_SUB(base_density_adjust_spin_4_n)
  end subroutine base_density_adjust_spin_4_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    PUSH_SUB(base_density_adjust_spin)

    select case(size(this))
    case(1)
      call base_density_adjust_spin_1_n(this(1), that)
    case(2)
      call base_density_adjust_spin_2_n(this, that)
    case(4)
      call base_density_adjust_spin_4_n(this, that)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(base_density_adjust_spin)
  end subroutine base_density_adjust_spin

  ! ---------------------------------------------------------
  subroutine base_density__copy__(this, that)
    type(base_density_t), target, intent(inout) :: this
    type(base_density_t),         intent(in)    :: that

    PUSH_SUB(base_density__copy__)

    call base_density__end__(this)
    if(associated(that%config))then
      call base_density__init__(this, that%config)
      if(associated(that%sim))then
        call base_density__start__(this, that%sim)
        call storage_copy(this%data, that%data)
        if(base_density_get_nspin(this)>1)then
          call storage_copy(this%total, that%total)
        else
          this%total=>this%data
        end if
      end if
    end if

    POP_SUB(base_density__copy__)
  end subroutine base_density__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_density_copy_density(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: osub, isub
    type(json_object_t),  pointer :: cnfg
    integer                       :: ierr

    PUSH_SUB(base_density_copy_density)

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
      call base_density__add__(this, osub, cnfg)
    end do
    call base_density_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_density_copy_density)
  end subroutine base_density_copy_density

  ! ---------------------------------------------------------
  subroutine base_density__end__density(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__end__density)

    if(base_density_get_nspin(this)>1)then
      call storage_end(this%total)
      SAFE_DEALLOCATE_P(this%total)
    end if
    call storage_end(this%data)
    call base_density__inull__(this)
    call config_dict_end(this%dict)
    call base_density_hash_end(this%hash)
    call base_density_list_end(this%list)

    POP_SUB(base_density__end__density)
  end subroutine base_density__end__density

  ! ---------------------------------------------------------
  recursive subroutine base_density_end_density(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_end_density)

    do
      nullify(subs)
      call base_density_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_density_end(subs)
      call base_density__idel__(subs)
    end do
    nullify(subs)
    call base_density__end__(this)

    POP_SUB(base_density_end_density)
  end subroutine base_density_end_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_init(this, that)
    type(base_density_iterator_t), intent(out) :: this
    type(base_density_t),  target, intent(in)  :: that

    PUSH_SUB(base_density_iterator_init)

    this%self => that
    call base_density_hash_init(this%iter, that%hash)

    POP_SUB(base_density_iterator_init)
  end subroutine base_density_iterator_init

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_config_density(this, config, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(json_object_t),          pointer        :: config
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)   :: ierr

    PUSH_SUB(base_density_iterator_next_config_density)

    call base_density_hash_next(this%iter, config, that, ierr)

    POP_SUB(base_density_iterator_next_config_density)
  end subroutine base_density_iterator_next_config_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_config(this, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(json_object_t),          pointer        :: that
    integer,             optional, intent(out)   :: ierr

    PUSH_SUB(base_density_iterator_next_config)

    call base_density_hash_next(this%iter, that, ierr)

    POP_SUB(base_density_iterator_next_config)
  end subroutine base_density_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_density(this, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)    :: ierr

    PUSH_SUB(base_density_iterator_next_density)

    call base_density_hash_next(this%iter, that, ierr)

    POP_SUB(base_density_iterator_next_density)
  end subroutine base_density_iterator_next_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_copy(this, that)
    type(base_density_iterator_t), intent(inout) :: this
    type(base_density_iterator_t), intent(in)    :: that

    PUSH_SUB(base_density_iterator_copy)

    this%self=>that%self
    call base_density_hash_copy(this%iter, that%iter)

    POP_SUB(base_density_iterator_copy)
  end subroutine base_density_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_density_iterator_end(this)
    type(base_density_iterator_t), intent(inout) :: this

    PUSH_SUB(base_density_iterator_end)

    nullify(this%self)
    call base_density_hash_end(this%iter)

    POP_SUB(base_density_iterator_end)
  end subroutine base_density_iterator_end

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_init(this, that, type)
    type(base_density_intrpl_t),  intent(out) :: this
    type(base_density_t), target, intent(in)  :: that
    integer,            optional, intent(in)  :: type

    PUSH_SUB(base_density_intrpl_init)

    this%self=>that
    call storage_init(this%intrp, that%data, type)

    POP_SUB(base_density_intrpl_init)
  end subroutine base_density_intrpl_init

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_eval_1d(this, x, val)
    type(base_density_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    real(kind=wp), dimension(base_density_get_nspin(this%self)) :: tvl
    real(kind=wp), dimension(1)                                 :: tv1
    integer                                                     :: ierr

    PUSH_SUB(base_density_intrpl_eval_1d)

    if(base_density_get_nspin(this%self)==1)then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call base_density_adjust_spin(tv1, tvl)
      val=tv1(1)
    end if

    POP_SUB(base_density_intrpl_eval_1d)
  end subroutine base_density_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_eval_2d(this, x, val)
    type(base_density_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val

    real(kind=wp), dimension(base_density_get_nspin(this%self)) :: tvl
    integer                                                     :: ierr

    PUSH_SUB(base_density_intrpl_eval_2d)

    if(base_density_get_nspin(this%self)==size(val))then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call base_density_adjust_spin(val, tvl)
    end if

    POP_SUB(base_density_intrpl_eval_2d)
  end subroutine base_density_intrpl_eval_2d

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_copy(this, that)
    type(base_density_intrpl_t), intent(out) :: this
    type(base_density_intrpl_t), intent(in)  :: that

    PUSH_SUB(base_density_intrpl_copy)

    this%self=>that%self
    call storage_copy(this%intrp, that%intrp)

    POP_SUB(base_density_intrpl_copy)
  end subroutine base_density_intrpl_copy

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_end(this)
    type(base_density_intrpl_t), intent(inout) :: this

    PUSH_SUB(base_density_intrpl_end)

    nullify(this%self)
    call storage_end(this%intrp)

    POP_SUB(base_density_intrpl_end)
  end subroutine base_density_intrpl_end

end module base_density_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

