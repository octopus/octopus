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

#define HASH_TEMPLATE_NAME base_states
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_states

module base_states_oct_m

  use base_density_oct_m
  use config_dict_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m

#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_states
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                  &
    BASE_STATES_OK,          &
    BASE_STATES_KEY_ERROR,   &
    BASE_STATES_EMPTY_ERROR

  public ::        &
    base_states_t

  public ::                &
    base_states__init__,   &
    base_states__start__,  &
    base_states__update__, &
    base_states__stop__,   &
    base_states__reset__,  &
    base_states__acc__,    &
    base_states__copy__,   &
    base_states__end__

  public ::             &
    base_states_new,    &
    base_states_del,    &
    base_states_init,   &
    base_states_start,  &
    base_states_update, &
    base_states_stop,   &
    base_states_sets,   &
    base_states_gets,   &
    base_states_set,    &
    base_states_get,    &
    base_states_copy,   &
    base_states_end

#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_STATES_OK          = BASE_STATES_HASH_OK
  integer, parameter :: BASE_STATES_KEY_ERROR   = BASE_STATES_HASH_KEY_ERROR
  integer, parameter :: BASE_STATES_EMPTY_ERROR = BASE_STATES_HASH_EMPTY_ERROR

  type :: base_states_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_states_t), pointer :: prnt   =>null()
    real(kind=wp)                :: charge = 0.0_wp
    type(base_density_t)         :: density
    type(config_dict_t)          :: dict
    type(base_states_hash_t)     :: hash
    type(base_states_list_t)     :: list
  end type base_states_t

  interface base_states__init__
    module procedure base_states__init__type
    module procedure base_states__init__copy
  end interface base_states__init__

  interface base_states_init
    module procedure base_states_init_type
    module procedure base_states_init_copy
  end interface base_states_init

  interface base_states_set
    module procedure base_states_set_info
  end interface base_states_set

  interface base_states_gets
    module procedure base_states_gets_config
    module procedure base_states_gets_type
    module procedure base_states_gets_density
    module procedure base_states_gets_density_1d
    module procedure base_states_gets_density_2d
  end interface base_states_gets

  interface base_states_get
    module procedure base_states_get_info
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_density
    module procedure base_states_get_density_1d
    module procedure base_states_get_density_2d
  end interface base_states_get

  interface base_states_copy
    module procedure base_states_copy_type
  end interface base_states_copy

  interface base_states_end
    module procedure base_states_end_type
  end interface base_states_end

#define TEMPLATE_PREFIX base_states
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains
    
#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_states__new__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_states__new__)
  end subroutine base_states__new__

  ! ---------------------------------------------------------
  subroutine base_states__del__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_states__del__)
  end subroutine base_states__del__

  ! ---------------------------------------------------------
  subroutine base_states_new(this, that)
    type(base_states_t),  target, intent(inout) :: this
    type(base_states_t), pointer                :: that

    PUSH_SUB(base_states_new)

    nullify(that)
    call base_states__new__(that)
    that%prnt => this
    call base_states_list_push(this%list, that)

    POP_SUB(base_states_new)
  end subroutine base_states_new

  ! ---------------------------------------------------------
  subroutine base_states_del(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_states_list_del(this%prnt%list, this)
        call base_states_end(this)
        call base_states__del__(this)
      end if
    end if
    nullify(this)

    POP_SUB(base_states_del)
  end subroutine base_states_del

  ! ---------------------------------------------------------
  subroutine base_states__init__type(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states__init__type)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK) this%charge = 0.0_wp
    call json_get(this%config, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_density__init__(this%density, cnfg)
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_states_hash_init(this%hash)
    call base_states_list_init(this%list)

    POP_SUB(base_states__init__type)
  end subroutine base_states__init__type
    
  ! ---------------------------------------------------------
  subroutine base_states__init__copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    PUSH_SUB(base_states__init__copy)

    ASSERT(associated(that%config))
    call base_states__init__(this, that%config)
    if(associated(that%sim)) call base_states__start__(this, that%sim)

    POP_SUB(base_states__init__copy)
  end subroutine base_states__init__copy
    
  ! ---------------------------------------------------------
  subroutine base_states_init_type(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_states_init_type)

    call base_states__init__(this, config)

    POP_SUB(base_states_init_type)
  end subroutine base_states_init_type
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_init_copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states_init_copy)

    nullify(cnfg, osub, isub)
    call base_states__init__(this, that)
    call base_states_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_states_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_new(this, osub)
      call base_states_init(osub, isub)
      call base_states_sets(this, osub, cnfg)
    end do
    call base_states_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_states_init_copy)
  end subroutine base_states_init_copy

  ! ---------------------------------------------------------
  subroutine base_states__start__(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_states__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call base_density__start__(this%density, sim)

    POP_SUB(base_states__start__)
  end subroutine base_states__start__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_start(this, sim)
    type(base_states_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_states_start)

    nullify(subs)
    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_start(subs, sim)
    end do
    call base_states_end(iter)
    nullify(subs)
    call base_states__start__(this, sim)

    POP_SUB(base_states_start)
  end subroutine base_states_start

  ! ---------------------------------------------------------
  subroutine base_states__update__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__update__(this%density)

    POP_SUB(base_states__update__)
  end subroutine base_states__update__

  ! ---------------------------------------------------------
  recursive subroutine base_states_update(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_states_update)

    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_update(subs)
    end do
    call base_states_end(iter)
    nullify(subs)
    if(associated(this%sim)) call base_states__update__(this)

    POP_SUB(base_states_update)
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states__stop__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_density__stop__(this%density)

    POP_SUB(base_states__stop__)
  end subroutine base_states__stop__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_stop(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr

    PUSH_SUB(base_states_stop)

    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_stop(subs)
    end do
    call base_states_end(iter)
    nullify(subs)
    call base_states__stop__(this)

    POP_SUB(base_states_stop)
  end subroutine base_states_stop

  ! ---------------------------------------------------------
  subroutine base_states__reset__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    this%charge = 0.0_wp
    call base_density__reset__(this%density)

    POP_SUB(base_states__reset__)
  end subroutine base_states__reset__
    
  ! ---------------------------------------------------------
  subroutine base_states__acc__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    this%charge = this%charge + that%charge
    call base_density__acc__(this%density, that%density)

    POP_SUB(base_states__acc__)
  end subroutine base_states__acc__

  ! ---------------------------------------------------------
  subroutine base_states__sets__(this, that, config)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(base_states__sets__)

    call base_density_sets(this%density, that%density, config)

    POP_SUB(base_states__sets__)
  end subroutine base_states__sets__
    
  ! ---------------------------------------------------------
  subroutine base_states_sets(this, that, config)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_states_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_states_hash_set(this%hash, config, that)
    call base_states__sets__(this, that, config)

    POP_SUB(base_states_sets)
  end subroutine base_states_sets
    
  ! ---------------------------------------------------------
  subroutine base_states_gets_config(this, config, that)
    type(base_states_t),  intent(in) :: this
    type(json_object_t),  intent(in) :: config
    type(base_states_t), pointer     :: that

    integer :: ierr

    PUSH_SUB(base_states_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_states_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_STATES_OK) nullify(that)

    POP_SUB(base_states_gets_config)
  end subroutine base_states_gets_config

  ! ---------------------------------------------------------
  subroutine base_states_gets_type(this, name, that)
    type(base_states_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_states_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_states_gets_type)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_states_gets(this, config, that)

    POP_SUB(base_states_gets_type)
  end subroutine base_states_gets_type

  ! ---------------------------------------------------------
  subroutine base_states_gets_density(this, name, that)
    type(base_states_t),   intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_gets_density)

    nullify(that, subs)
    call base_states_gets(this, name, subs)
    if(associated(subs)) call base_states_get(subs, that)

    POP_SUB(base_states_gets_density)
  end subroutine base_states_gets_density
    
  ! ---------------------------------------------------------
  subroutine base_states_gets_density_1d(this, name, that, total)
    type(base_states_t),          intent(in) :: this
    character(len=*),             intent(in) :: name
    real(kind=wp), dimension(:), pointer     :: that
    logical,            optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_gets_density_1d)

    nullify(that, dnst)
    call base_states_gets(this, name, dnst)
    if(associated(dnst)) call base_density_get(dnst, that, total)

    POP_SUB(base_states_gets_density_1d)
  end subroutine base_states_gets_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_gets_density_2d(this, name, that, total)
    type(base_states_t),            intent(in) :: this
    character(len=*),               intent(in) :: name
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_gets_density_2d)

    nullify(that, dnst)
    call base_states_gets(this, name, dnst)
    if(associated(dnst)) call base_density_get(dnst, that, total)

    POP_SUB(base_states_gets_density_2d)
  end subroutine base_states_gets_density_2d

  ! ---------------------------------------------------------
  subroutine base_states_set_info(this, charge)
    type(base_states_t),     intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: charge

    PUSH_SUB(base_states_set_info)

    if(present(charge)) this%charge = charge

    POP_SUB(base_states_set_info)
  end subroutine base_states_set_info
    
  ! ---------------------------------------------------------
  subroutine base_states_get_info(this, charge, nspin)
    type(base_states_t),     intent(in)  :: this
    real(kind=wp), optional, intent(out) :: charge
    integer,       optional, intent(out) :: nspin

    PUSH_SUB(base_states_get_info)

    if(present(charge)) charge = this%charge
    call base_density_get(this%density, nspin=nspin)

    POP_SUB(base_states_get_info)
  end subroutine base_states_get_info
    
  ! ---------------------------------------------------------
  subroutine base_states_get_config(this, that)
    type(base_states_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_states_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_states_get_config)
  end subroutine base_states_get_config
    
  ! ---------------------------------------------------------
  subroutine base_states_get_simulation(this, that)
    type(base_states_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_states_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_states_get_simulation)
  end subroutine base_states_get_simulation
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density(this, that)
    type(base_states_t),   target, intent(in) :: this
    type(base_density_t), pointer             :: that

    PUSH_SUB(base_states_get_density)

    that => this%density

    POP_SUB(base_states_get_density)
  end subroutine base_states_get_density
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density_1d(this, that, total)
    type(base_states_t),          intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    logical,            optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_1d)

    nullify(that)
    call base_density_get(this%density, that, total)

    POP_SUB(base_states_get_density_1d)
  end subroutine base_states_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_get_density_2d(this, that, total)
    type(base_states_t),            intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_2d)

    nullify(that)
    call base_density_get(this%density, that, total)

    POP_SUB(base_states_get_density_2d)
  end subroutine base_states_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_states__copy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__copy__)

    call base_states__end__(this)
    if(associated(that%config))then
      call base_states__init__(this, that)
      this%charge = that%charge
      if(associated(that%sim)) call base_density__copy__(this%density, that%density)
    end if

    POP_SUB(base_states__copy__)
  end subroutine base_states__copy__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_copy_type(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states_copy_type)

    nullify(cnfg, osub, isub)
    call base_states_end(this)
    call base_states__copy__(this, that)
    call base_states_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_states_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_new(this, osub)
      call base_states_copy(osub, isub)
      call base_states_sets(this, osub, cnfg)
    end do
    call base_states_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_states_copy_type)
  end subroutine base_states_copy_type

  ! ---------------------------------------------------------
  subroutine base_states__end__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__end__)

    nullify(this%config, this%sim, this%prnt)
    this%charge = 0.0_wp
    call base_density__end__(this%density)
    call config_dict_end(this%dict)
    call base_states_hash_end(this%hash)
    call base_states_list_end(this%list)

    POP_SUB(base_states__end__)
  end subroutine base_states__end__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_end_type(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_end_type)

    do
      nullify(subs)
      call base_states_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_states_end(subs)
      call base_states__del__(subs)
    end do
    nullify(subs)
    call base_states__end__(this)

    POP_SUB(base_states_end_type)
  end subroutine base_states_end_type

#define TEMPLATE_PREFIX base_states
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_states_oct_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
