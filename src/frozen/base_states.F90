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

module base_states_m

  use base_density_m
  use config_dict_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use simulation_m

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
    module procedure base_states__init__states
    module procedure base_states__init__copy
  end interface base_states__init__

  interface base_states_init
    module procedure base_states_init_states
    module procedure base_states_init_copy
  end interface base_states_init

  interface base_states_set
    module procedure base_states_set_info
    module procedure base_states_set_simulation
  end interface base_states_set

  interface base_states_gets
    module procedure base_states_gets_config
    module procedure base_states_gets_name
  end interface base_states_gets

  interface base_states_get
    module procedure base_states_get_info
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_density
  end interface base_states_get

  interface base_states_copy
    module procedure base_states_copy_states
  end interface base_states_copy

  interface base_states_end
    module procedure base_states_end_states
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
  subroutine base_states_new(this, that)
    type(base_states_t),  target, intent(inout) :: this
    type(base_states_t), pointer                :: that

    PUSH_SUB(base_states_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt => this
    call base_states_list_push(this%list, that)

    POP_SUB(base_states_new)
  end subroutine base_states_new

  ! ---------------------------------------------------------
  subroutine base_states__idel__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_states__idel__)
  end subroutine base_states__idel__

  ! ---------------------------------------------------------
  subroutine base_states_del(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_states_list_del(this%prnt%list, this)
        call base_states_end(this)
        call base_states__idel__(this)
      end if
    end if

    POP_SUB(base_states_del)
  end subroutine base_states_del

  ! ---------------------------------------------------------
  subroutine base_states__inull__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__inull__)

    nullify(this%config, this%sim, this%prnt)
    this%charge = 0.0_wp

    POP_SUB(base_states__inull__)
  end subroutine base_states__inull__
    
  ! ---------------------------------------------------------
  subroutine base_states__iinit__(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    integer :: ierr

    PUSH_SUB(base_states__iinit__)

    call base_states__inull__(this)
    this%config => config
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK) this%charge = 0.0_wp
    call config_dict_init(this%dict)
    call base_states_hash_init(this%hash)
    call base_states_list_init(this%list)

    POP_SUB(base_states__iinit__)
  end subroutine base_states__iinit__
    
  ! ---------------------------------------------------------
  subroutine base_states__init__states(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states__init__states)

    nullify(cnfg)
    call base_states__iinit__(this, config)
    call json_get(this%config, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_density__init__(this%density, cnfg)
    nullify(cnfg)

    POP_SUB(base_states__init__states)
  end subroutine base_states__init__states
    
  ! ---------------------------------------------------------
  subroutine base_states__init__copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    PUSH_SUB(base_states__init__copy)

    ASSERT(associated(that%config))
    call base_states__iinit__(this, that%config)
    call base_density__init__(this%density, that%density)

    POP_SUB(base_states__init__copy)
  end subroutine base_states__init__copy
    
  ! ---------------------------------------------------------
  subroutine base_states_init_states(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_states_init_states)

    call base_states__init__(this, config)

    POP_SUB(base_states_init_states)
  end subroutine base_states_init_states
    
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
  subroutine base_states__istart__(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_states__istart__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim

    POP_SUB(base_states__istart__)
  end subroutine base_states__istart__
    
  ! ---------------------------------------------------------
  subroutine base_states__start__(this, sim)
    type(base_states_t),          intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

    PUSH_SUB(base_states__start__)

    if(present(sim))then
      call base_states__istart__(this, sim)
    else
      if(.not.associated(this%sim))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sim))
        call base_states__istart__(this, this%prnt%sim)
      end if
    end if
    call base_density__start__(this%density, sim)

    POP_SUB(base_states__start__)
  end subroutine base_states__start__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_start(this, sim)
    type(base_states_t),          intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: sim

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

    nullify(subs)
    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call base_states_update(subs)
    end do
    call base_states_end(iter)
    nullify(subs)
    call base_states__update__(this)

    POP_SUB(base_states_update)
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states__stop__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
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

    nullify(subs)
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
    ASSERT(associated(this%sim))
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
  subroutine base_states_gets_name(this, name, that)
    type(base_states_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_states_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_states_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_states_gets(this, config, that)

    POP_SUB(base_states_gets_name)
  end subroutine base_states_gets_name

  ! ---------------------------------------------------------
  subroutine base_states_set_info(this, charge)
    type(base_states_t),     intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: charge

    PUSH_SUB(base_states_set_info)

    if(present(charge)) this%charge = charge

    POP_SUB(base_states_set_info)
  end subroutine base_states_set_info
    
  ! ---------------------------------------------------------
  subroutine base_states_set_simulation(this, that)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that

    PUSH_SUB(base_states_set_simulation)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => that

    POP_SUB(base_states_set_simulation)
  end subroutine base_states_set_simulation
    
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
  subroutine base_states__icopy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__icopy__)

    call base_states__iend__(this)
    if(associated(that%config))then
      call base_states__iinit__(this, that%config)
      this%charge = that%charge
      if(associated(that%sim)) call base_states__istart__(this, that%sim)
    end if

    POP_SUB(base_states__icopy__)
  end subroutine base_states__icopy__
    
  ! ---------------------------------------------------------
  subroutine base_states__copy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__copy__)

    call base_states__icopy__(this, that)
    call base_density__copy__(this%density, that%density)

    POP_SUB(base_states__copy__)
  end subroutine base_states__copy__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_copy_states(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states_copy_states)

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

    POP_SUB(base_states_copy_states)
  end subroutine base_states_copy_states

  ! ---------------------------------------------------------
  subroutine base_states__iend__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__iend__)

    call base_states__inull__(this)
    call config_dict_end(this%dict)
    call base_states_hash_end(this%hash)
    call base_states_list_end(this%list)

    POP_SUB(base_states__iend__)
  end subroutine base_states__iend__
    
  ! ---------------------------------------------------------
  subroutine base_states__end__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__end__)

    call base_states__iend__(this)
    call base_density__end__(this%density)

    POP_SUB(base_states__end__)
  end subroutine base_states__end__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_end_states(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_end_states)

    do
      nullify(subs)
      call base_states_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_states_end(subs)
      call base_states__idel__(subs)
    end do
    nullify(subs)
    call base_states__end__(this)

    POP_SUB(base_states_end_states)
  end subroutine base_states_end_states

#define TEMPLATE_PREFIX base_states
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_states_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
