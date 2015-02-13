#include "global.h"

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

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_states_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

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

  use simulation_m, only: &
    simulation_t

  use base_density_m, only: &
    base_density__start__,  &
    base_density__update__, &
    base_density__stop__,   &
    base_density__add__

  use base_density_m, only: &
    base_density_t,         &
    base_density_init,      &
    base_density_get,       &
    base_density_copy,      &
    base_density_end

  implicit none

  private
  public ::                &
    base_states__start__,  &
    base_states__update__, &
    base_states__stop__,   &
    base_states__add__

  public ::             &
    base_states_init,   &
    base_states_start,  &
    base_states_update, &
    base_states_stop,   &
    base_states_next,   &
    base_states_set,    &
    base_states_get,    &
    base_states_copy,   &
    base_states_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_states_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    real(kind=wp)                :: charge = 0.0_wp
    type(base_density_t)         :: density
    type(config_dict_t)          :: dict
    type(base_states_hash_t)     :: hash
  end type base_states_t

  type, public :: base_states_iterator_t
    private
    type(base_states_t),      pointer :: self =>null()
    type(base_states_hash_iterator_t) :: iter
  end type base_states_iterator_t

  interface base_states_init
    module procedure base_states_init_states
    module procedure base_states_iterator_init
  end interface base_states_init

  interface base_states_next
    module procedure base_states_iterator_next_config_states
    module procedure base_states_iterator_next_config
    module procedure base_states_iterator_next_states
  end interface base_states_next

  interface base_states_set
    module procedure base_states_set_charge
    module procedure base_states_set_simulation
  end interface base_states_set

  interface base_states_get
    module procedure base_states_get_info
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_density
  end interface base_states_get

  interface base_states_copy
    module procedure base_states_copy_states
    module procedure base_states_iterator_copy
  end interface base_states_copy

  interface base_states_end
    module procedure base_states_end_states
    module procedure base_states_iterator_end
  end interface base_states_end

  integer, public, parameter :: BASE_STATES_OK          = BASE_STATES_HASH_OK
  integer, public, parameter :: BASE_STATES_KEY_ERROR   = BASE_STATES_HASH_KEY_ERROR
  integer, public, parameter :: BASE_STATES_EMPTY_ERROR = BASE_STATES_HASH_EMPTY_ERROR

contains
    
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_states_init_states(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_states_init_states)
    this%config=>config
    nullify(this%sim, cnfg)
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK)this%charge=0.0_wp
    nullify(cnfg)
    call json_get(config, "density", cnfg, ierr)
    if(ierr==JSON_OK)call base_density_init(this%density, cnfg)
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_states_hash_init(this%hash)
    POP_SUB(base_states_init_states)
    return
  end subroutine base_states_init_states
    
  ! ---------------------------------------------------------
  subroutine base_states__start__(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_states__start__)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call base_density__start__(this%density, sim)
    POP_SUB(base_states__start__)
    return
  end subroutine base_states__start__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_start(this, sim)
    type(base_states_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim
    !
    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr
    !
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
    return
  end subroutine base_states_start

  ! ---------------------------------------------------------
  subroutine base_states__update__(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states__update__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__update__(this%density)
    POP_SUB(base_states__update__)
    return
  end subroutine base_states__update__

  ! ---------------------------------------------------------
  recursive subroutine base_states_update(this)
    type(base_states_t), intent(inout) :: this
    !
    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr
    !
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
    return
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states__stop__(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states__stop__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__stop__(this%density)
    POP_SUB(base_states__stop__)
    return
  end subroutine base_states__stop__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_stop(this)
    type(base_states_t), intent(inout) :: this
    !
    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    integer                      :: ierr
    !
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
    return
  end subroutine base_states_stop

  ! ---------------------------------------------------------
  subroutine base_states__add__(this, that, config)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(base_states__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_states_hash_set(this%hash, config, that)
    call base_density__add__(this%density, that%density, config)
    POP_SUB(base_states__add__)
    return
  end subroutine base_states__add__
    
  ! ---------------------------------------------------------
  subroutine base_states__get__(this, name, that)
    type(base_states_t),  intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(base_states_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_states__get__)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_states_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_STATES_OK)nullify(that)
    end if
    POP_SUB(base_states__get__)
    return
  end subroutine base_states__get__

  ! ---------------------------------------------------------
  elemental subroutine base_states_set_charge(this, that)
    type(base_states_t), intent(inout) :: this
    real(kind=wp),       intent(in)    :: that
    !
    this%charge=that
    return
  end subroutine base_states_set_charge
    
  ! ---------------------------------------------------------
  subroutine base_states_set_simulation(this, that)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_states_set_simulation)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>that
    POP_SUB(base_states_set_simulation)
    return
  end subroutine base_states_set_simulation
    
  ! ---------------------------------------------------------
  elemental function base_states_get_charge(this) result(that)
    type(base_states_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%charge
    return
  end function base_states_get_charge
    
  ! ---------------------------------------------------------
  elemental subroutine base_states_get_info(this, charge)
    type(base_states_t),     intent(in)  :: this
    real(kind=wp), optional, intent(out) :: charge
    !
    if(present(charge))&
      charge=base_states_get_charge(this)
    return
  end subroutine base_states_get_info
    
  ! ---------------------------------------------------------
  subroutine base_states_get_config(this, that)
    type(base_states_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_states_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_states_get_config)
    return
  end subroutine base_states_get_config
    
  ! ---------------------------------------------------------
  subroutine base_states_get_simulation(this, that)
    type(base_states_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(base_states_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_states_get_simulation)
    return
  end subroutine base_states_get_simulation
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density(this, that)
    type(base_states_t),   target, intent(in) :: this
    type(base_density_t), pointer             :: that
    !
    PUSH_SUB(base_states_get_density)
    that=>this%density
    POP_SUB(base_states_get_density)
    return
  end subroutine base_states_get_density
    
  ! ---------------------------------------------------------
  subroutine base_states_copy_states(this, that)
    type(base_states_t),         intent(out) :: this
    type(base_states_t), target, intent(in)  :: that
    !
    PUSH_SUB(base_states_copy_states)
    this%config=>that%config
    this%sim=>that%sim
    this%charge=that%charge
    call base_density_copy(this%density, that%density)
    call config_dict_copy(this%dict, that%dict)
    call base_states_hash_copy(this%hash, that%hash)
    POP_SUB(base_states_copy_states)
    return
  end subroutine base_states_copy_states
    
  ! ---------------------------------------------------------
  subroutine base_states_end_states(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_end_states)
    nullify(this%config, this%sim)
    this%charge=0.0_wp
    call base_density_end(this%density)
    call config_dict_end(this%dict)
    call base_states_hash_end(this%hash)
    POP_SUB(base_states_end_states)
    return
  end subroutine base_states_end_states
    
  ! ---------------------------------------------------------
  subroutine base_states_iterator_init(this, that)
    type(base_states_iterator_t), intent(out) :: this
    type(base_states_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_states_iterator_init)
    this%self=>that
    call base_states_hash_init(this%iter, that%hash)
    POP_SUB(base_states_iterator_init)
    return
  end subroutine base_states_iterator_init

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_config_states(this, config, system, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: config
    type(base_states_t),         pointer        :: system
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_config_states)
    call base_states_hash_next(this%iter, config, system, ierr)
    POP_SUB(base_states_iterator_next_config_states)
    return
  end subroutine base_states_iterator_next_config_states

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_config(this, that, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_config)
    call base_states_hash_next(this%iter, that, ierr)
    POP_SUB(base_states_iterator_next_config)
    return
  end subroutine base_states_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_states(this, that, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(base_states_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_states)
    call base_states_hash_next(this%iter, that, ierr)
    POP_SUB(base_states_iterator_next_states)
    return
  end subroutine base_states_iterator_next_states

  ! ---------------------------------------------------------
  subroutine base_states_iterator_copy(this, that)
    type(base_states_iterator_t), intent(inout) :: this
    type(base_states_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_states_iterator_copy)
    this%self=>that%self
    call base_states_hash_copy(this%iter, that%iter)
    POP_SUB(base_states_iterator_copy)
    return
  end subroutine base_states_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_states_iterator_end(this)
    type(base_states_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_iterator_end)
    nullify(this%self)
    call base_states_hash_end(this%iter)
    POP_SUB(base_states_iterator_end)
    return
  end subroutine base_states_iterator_end

end module base_states_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
