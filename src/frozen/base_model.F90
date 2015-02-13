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

#define HASH_TEMPLATE_NAME base_model
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_model

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_model_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m, only: JSON_OK, json_object_t, json_get

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

  use space_m, only: &
    space_t

  use geometry_m, only: &
    geometry_t

  use igrid_m, only: &
    grid_t

  use simulation_m, only: &
    simulation__add__

  use simulation_m, only: &
    simulation_t,         &
    simulation_init,      &
    simulation_start,     &
    simulation_set,       &
    simulation_copy,      &
    simulation_end

  use base_geom_m, only: &
    base_geom_t,         &
    base_geom_get

  use base_system_m, only: &
    base_system__init__,   &
    base_system__start__,  &
    base_system__update__, &
    base_system__stop__,   &
    base_system__add__

  use base_system_m, only: &
    base_system_t,         &
    base_system_init,      &
    base_system_set,       &
    base_system_get,       &
    base_system_copy,      &
    base_system_end

  use base_hamiltonian_m, only: &
    base_hamiltonian__start__,  &
    base_hamiltonian__update__, &
    base_hamiltonian__stop__,   &
    base_hamiltonian__add__

  use base_hamiltonian_m, only: &
    base_hamiltonian_t,         &
    base_hamiltonian_init,      &
    base_hamiltonian_get,       &
    base_hamiltonian_copy,      &
    base_hamiltonian_end

  implicit none

  private
  public ::               &
    base_model__init__,   &
    base_model__start__,  &
    base_model__update__, &
    base_model__stop__,   &
    base_model__add__,    &
    base_model__get__

  public ::            &
    base_model_init,   &
    base_model_start,  &
    base_model_update, &
    base_model_stop,   &
    base_model_next,   &
    base_model_get,    &
    base_model_copy,   &
    base_model_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_model_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t)           :: sim
    type(base_system_t)          :: sys
    type(base_hamiltonian_t)     :: hm
    type(config_dict_t)          :: dict
    type(base_model_hash_t)      :: hash
  end type base_model_t

  type, public :: base_model_iterator_t
    private
    type(base_model_t),      pointer :: self =>null()
    type(base_model_hash_iterator_t) :: iter
  end type base_model_iterator_t

  interface base_model__init__
    module procedure base_model__init__begin
    module procedure base_model__init__finish
 end interface base_model__init__

  interface base_model_init
    module procedure base_model_init_model
    module procedure base_model_iterator_init
  end interface base_model_init

  interface base_model_get
    module procedure base_model_get_config
    module procedure base_model_get_simulation
    module procedure base_model_get_system
    module procedure base_model_get_hamiltonian
  end interface base_model_get

  interface base_model_next
    module procedure base_model_iterator_next_config_model
    module procedure base_model_iterator_next_config
    module procedure base_model_iterator_next_model
  end interface base_model_next

  interface base_model_copy
    module procedure base_model_copy_model
    module procedure base_model_iterator_copy
  end interface base_model_copy

  interface base_model_end
    module procedure base_model_end_model
    module procedure base_model_iterator_end
  end interface base_model_end

  integer, public, parameter :: BASE_MODEL_OK          = BASE_MODEL_HASH_OK
  integer, public, parameter :: BASE_MODEL_KEY_ERROR   = BASE_MODEL_HASH_KEY_ERROR
  integer, public, parameter :: BASE_MODEL_EMPTY_ERROR = BASE_MODEL_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_model__init__begin(this, config)
    type(base_model_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    type(space_t),       pointer :: space
    integer                      :: ierr
    !
    PUSH_SUB(base_model__init__begin)
    nullify(cnfg, space)
    this%config=>config
    call json_get(this%config, "system", cnfg, ierr)
    if(ierr==JSON_OK)then
       call base_system__init__(this%sys, cnfg)
       nullify(cnfg)
       call base_system_get(this%sys, space)
       ASSERT(associated(space))
       call json_get(this%config, "simulation", cnfg, ierr)
       if(ierr==JSON_OK)call simulation_init(this%sim, space, cnfg)
       nullify(cnfg, space)
       call json_get(this%config, "hamiltonian", cnfg, ierr)
       if(ierr==JSON_OK)call base_hamiltonian_init(this%hm, this%sys, cnfg)
    end if
    nullify(cnfg)
    call config_dict_init(this%dict)
    call base_model_hash_init(this%hash)
    POP_SUB(base_model__init__begin)
    return
  end subroutine base_model__init__begin

  ! ---------------------------------------------------------
  subroutine base_model__init__finish(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model__init__finish)
    call base_system__init__(this%sys)
    POP_SUB(base_model__init__finish)
    return
  end subroutine base_model__init__finish

  ! ---------------------------------------------------------
  subroutine base_model_init_model(this, config)
    type(base_model_t),  intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(base_model_init_model)
    call base_model__init__begin(this, config)
    call base_model__init__finish(this)
    POP_SUB(base_model_init_model)
    return
  end subroutine base_model_init_model

  ! ---------------------------------------------------------
  subroutine base_model__start__(this, grid)
    type(base_model_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid
    !
    type(base_geom_t), pointer :: geom
    type(geometry_t),  pointer :: geo
    !
    PUSH_SUB(base_model__start__)
    call base_system_get(this%sys, geom)
    ASSERT(associated(geom))
    call base_geom_get(geom, geo)
    ASSERT(associated(geo))
    call simulation_start(this%sim, grid, geo)
    call base_system__start__(this%sys, this%sim)
    call base_hamiltonian__start__(this%hm, this%sim)
    POP_SUB(base_model__start__)
    return
  end subroutine base_model__start__

  ! ---------------------------------------------------------
  recursive subroutine base_model_start(this, grid)
    type(base_model_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid
    !
    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr
    !
    PUSH_SUB(base_model_start)
    nullify(subs)
    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_start(subs, grid)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__start__(this, grid)
    POP_SUB(base_model_start)
    return
  end subroutine base_model_start

  ! ---------------------------------------------------------
  subroutine base_model__update__(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model__update__)
    call base_system__update__(this%sys)
    call base_hamiltonian__update__(this%hm)
    POP_SUB(base_model__update__)
    return
  end subroutine base_model__update__

  ! ---------------------------------------------------------
  recursive subroutine base_model_update(this)
    type(base_model_t), intent(inout) :: this
    !
    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr
    !
    PUSH_SUB(base_model_update)
    nullify(subs)
    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_update(subs)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__update__(this)
    POP_SUB(base_model_update)
    return
  end subroutine base_model_update

  ! ---------------------------------------------------------
  subroutine base_model__stop__(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model__stop__)
    call base_system__stop__(this%sys)
    call base_hamiltonian__stop__(this%hm)
    POP_SUB(base_model__stop__)
    return
  end subroutine base_model__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_model_stop(this)
    type(base_model_t), intent(inout) :: this
    !
    type(base_model_iterator_t) :: iter
    type(base_model_t), pointer :: subs
    integer                     :: ierr
    !
    PUSH_SUB(base_model_stop)
    nullify(subs)
    call base_model_init(iter, this)
    do
      nullify(subs)
      call base_model_next(iter, subs, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_stop(subs)
    end do
    call base_model_end(iter)
    nullify(subs)
    call base_model__stop__(this)
    POP_SUB(base_model_stop)
    return
  end subroutine base_model_stop

  ! ---------------------------------------------------------
  subroutine base_model__add__(this, that, config)
    type(base_model_t),  intent(inout) :: this
    type(base_model_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(base_model__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_model_hash_set(this%hash, config, that)
    call simulation__add__(this%sim, that%sim, config)
    call base_system__add__(this%sys, that%sys, config)
    call base_hamiltonian__add__(this%hm, that%hm, config)
    POP_SUB(base_model__add__)
    return
  end subroutine base_model__add__

  ! ---------------------------------------------------------
  subroutine base_model__get__(this, name, that)
    type(base_model_t),  intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_model_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_model__get__)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_model_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_MODEL_OK)nullify(that)
    end if
    POP_SUB(base_model__get__)
    return
  end subroutine base_model__get__

  ! ---------------------------------------------------------
  subroutine base_model_get_config(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_model_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_model_get_config)
    return
  end subroutine base_model_get_config

  ! ---------------------------------------------------------
  subroutine base_model_get_simulation(this, that)
    type(base_model_t),  target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(base_model_get_simulation)
    that=>this%sim
    POP_SUB(base_model_get_simulation)
    return
  end subroutine base_model_get_simulation

  ! ---------------------------------------------------------
  subroutine base_model_get_system(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(base_system_t), pointer             :: that
    !
    PUSH_SUB(base_model_get_system)
    that=>this%sys
    POP_SUB(base_model_get_system)
    return
  end subroutine base_model_get_system

  ! ---------------------------------------------------------
  subroutine base_model_get_hamiltonian(this, that)
    type(base_model_t),        target, intent(in) :: this
    type(base_hamiltonian_t), pointer             :: that
    !
    PUSH_SUB(base_model_get_hamiltonian)
    that=>this%hm
    POP_SUB(base_model_get_hamiltonian)
    return
  end subroutine base_model_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_model_copy_model(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that
    !
    PUSH_SUB(base_model_copy_model)
    this%config=>that%config
    call simulation_copy(this%sim, that%sim)
    call base_system_copy(this%sys, that%sys)
    call base_hamiltonian_copy(this%hm, that%hm)
    call config_dict_copy(this%dict, that%dict)
    call base_model_hash_copy(this%hash, that%hash)
    POP_SUB(base_model_copy_model)
    return
  end subroutine base_model_copy_model

  ! ---------------------------------------------------------
  subroutine base_model_end_model(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_end_model)
    nullify(this%config)
    call base_hamiltonian_end(this%hm)
    call base_system_end(this%sys)
    call simulation_end(this%sim)
    call config_dict_end(this%dict)
    call base_model_hash_end(this%hash)
    POP_SUB(base_model_end_model)
    return
  end subroutine base_model_end_model

  ! ---------------------------------------------------------
  subroutine base_model_iterator_init(this, that)
    type(base_model_iterator_t), intent(out) :: this
    type(base_model_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_model_iterator_init)
    this%self=>that
    call base_model_hash_init(this%iter, that%hash)
    POP_SUB(base_model_iterator_init)
    return
  end subroutine base_model_iterator_init

  ! ---------------------------------------------------------
  subroutine base_model_iterator_next_config_model(this, config, that, ierr)
    type(base_model_iterator_t), intent(inout) :: this
    type(json_object_t),        pointer        :: config
    type(base_model_t),         pointer        :: that
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_model_iterator_next_config_model)
    call base_model_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_model_iterator_next_config_model)
    return
  end subroutine base_model_iterator_next_config_model

  ! ---------------------------------------------------------
  subroutine base_model_iterator_next_config(this, that, ierr)
    type(base_model_iterator_t), intent(inout) :: this
    type(json_object_t),        pointer        :: that
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_model_iterator_next_config)
    call base_model_hash_next(this%iter, that, ierr)
    POP_SUB(base_model_iterator_next_config)
    return
  end subroutine base_model_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_model_iterator_next_model(this, that, ierr)
    type(base_model_iterator_t), intent(inout) :: this
    type(base_model_t),         pointer        :: that
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_model_iterator_next_model)
    call base_model_hash_next(this%iter, that, ierr)
    POP_SUB(base_model_iterator_next_model)
    return
  end subroutine base_model_iterator_next_model

  ! ---------------------------------------------------------
  subroutine base_model_iterator_copy(this, that)
    type(base_model_iterator_t), intent(inout) :: this
    type(base_model_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_model_iterator_copy)
    this%self=>that%self
    call base_model_hash_copy(this%iter, that%iter)
    POP_SUB(base_model_iterator_copy)
    return
  end subroutine base_model_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_model_iterator_end(this)
    type(base_model_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_iterator_end)
    nullify(this%self)
    call base_model_hash_end(this%iter)
    POP_SUB(base_model_iterator_end)
    return
  end subroutine base_model_iterator_end

end module base_model_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
