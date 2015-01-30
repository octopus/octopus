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

  use space_m, only: &
    space_t

  use geometry_m, only: &
    geometry_t

  use grid_m, only: &
    grid_t

  use simulation_m, only: &
    simulation_t,         &
    simulation_init,      &
    simulation_start,     &
    simulation_set,       &
    simulation_copy,      &
    simulation_end

  use base_geom_m, only:       &
    geom_t   => base_geom_t,   &
    geom_get => base_geom_get

  use base_system_m, only:                  &
    system_t      => base_system_t,         &
    system_init   => base_system_init,      &
    system_start  => base_system_start,     &
    system_update => base_system_update,    &
    system_stop   => base_system_stop,      &
    system_get    => base_system_get,       &
    system_copy   => base_system_copy,      &
    system_end    => base_system_end

  use base_hamiltonian_m, only:                       &
    hamiltonian_t      => base_hamiltonian_t,         &
    hamiltonian_init   => base_hamiltonian_init,      &
    hamiltonian_start  => base_hamiltonian_start,     &
    hamiltonian_update => base_hamiltonian_update,    &
    hamiltonian_stop   => base_hamiltonian_stop,      &
    hamiltonian_get    => base_hamiltonian_get,       &
    hamiltonian_copy   => base_hamiltonian_copy,      &
    hamiltonian_end    => base_hamiltonian_end

  implicit none

  private
  public ::            &
    base_model_init,   &
    base_model_start,  &
    base_model_update, &
    base_model_stop,   &
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
    type(system_t)               :: sys
    type(hamiltonian_t)          :: hm
    type(base_model_hash_t)      :: hash
  end type base_model_t

  type, public :: base_model_iterator_t
    private
    type(base_model_t),      pointer :: self =>null()
    type(base_model_hash_iterator_t) :: iter
  end type base_model_iterator_t

  interface base_model_init
    module procedure base_model_init_begin
    module procedure base_model_init_build
    module procedure base_model_init_finish
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

  integer, parameter :: BASE_MODEL_OK          = BASE_MODEL_HASH_OK
  integer, parameter :: BASE_MODEL_KEY_ERROR   = BASE_MODEL_HASH_KEY_ERROR
  integer, parameter :: BASE_MODEL_EMPTY_ERROR = BASE_MODEL_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_model_init_begin(this, config)
    type(base_model_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_model_init_begin)
    this%config=>config
    nullify(cnfg)
    call json_get(this%config, "simulation", cnfg, ierr)
    if(ierr==JSON_OK)call simulation_init(this%sim, cnfg)
    nullify(cnfg)
    call json_get(this%config, "system", cnfg, ierr)
    if(ierr==JSON_OK)call system_init(this%sys, cnfg)
    nullify(cnfg)
    call json_get(this%config, "hamiltonian", cnfg, ierr)
    if(ierr==JSON_OK)call hamiltonian_init(this%hm, this%sys, cnfg)
    nullify(cnfg)
    call base_model_hash_init(this%hash)
    POP_SUB(base_model_init_begin)
    return
  end subroutine base_model_init_begin

  ! ---------------------------------------------------------
  subroutine base_model_init_build(this, that, config)
    type(base_model_t),  intent(inout) :: this
    type(base_model_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(base_model_init_build)
    call base_model_hash_set(this%hash, config, that)
    call simulation_init(this%sim, that%sim, config)
    call system_init(this%sys, that%sys, config)
    call hamiltonian_init(this%hm, that%hm, config)
    POP_SUB(base_model_init_build)
    return
  end subroutine base_model_init_build

  ! ---------------------------------------------------------
  subroutine base_model_init_finish(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_init_finish)
    call system_init(this%sys)
    POP_SUB(base_model_init_finish)
    return
  end subroutine base_model_init_finish

  ! ---------------------------------------------------------
  subroutine base_model_start(this, grid)
    type(base_model_t),     intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid
    !
    type(geom_t),     pointer :: geom
    type(geometry_t), pointer :: geo
    type(space_t),    pointer :: space
    !
    PUSH_SUB(base_model_start)
    nullify(geom, geo, space)
    call system_get(this%sys, geom)
    ASSERT(associated(geom))
    call system_get(this%sys, space)
    ASSERT(associated(space))
    call geom_get(geom, geo)
    ASSERT(associated(geo))
    call simulation_start(this%sim, grid, geo, space)
    nullify(geom, geo, space)
    call system_start(this%sys, this%sim)
    call hamiltonian_start(this%hm, this%sim)
    POP_SUB(base_model_start)
    return
  end subroutine base_model_start

  ! ---------------------------------------------------------
  subroutine base_model_update(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_update)
    call system_update(this%sys)
    call hamiltonian_update(this%hm)
    POP_SUB(base_model_update)
    return
  end subroutine base_model_update

  ! ---------------------------------------------------------
  subroutine base_model_stop(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_stop)
    call system_stop(this%sys)
    call hamiltonian_stop(this%hm)
    POP_SUB(base_model_stop)
    return
  end subroutine base_model_stop

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
    type(base_model_t), target, intent(in) :: this
    type(system_t),    pointer             :: that
    !
    PUSH_SUB(base_model_get_system)
    that=>this%sys
    POP_SUB(base_model_get_system)
    return
  end subroutine base_model_get_system

  ! ---------------------------------------------------------
  subroutine base_model_get_hamiltonian(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(hamiltonian_t), pointer             :: that
    !
    PUSH_SUB(base_model_get_hamiltonian)
    that=>this%hm
    POP_SUB(base_model_get_hamiltonian)
    return
  end subroutine base_model_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_model_copy_model(this, that)
    type(base_model_t), intent(out) :: this
    type(base_model_t), intent(in)  :: that
    !
    PUSH_SUB(base_model_copy_model)
    this%config=>that%config
    call simulation_copy(this%sim, that%sim)
    call system_copy(this%sys, that%sys)
    call hamiltonian_copy(this%hm, that%hm)
    call base_model_hash_copy(this%hash, that%hash)
    POP_SUB(base_model_copy_model)
    return
  end subroutine base_model_copy_model

  ! ---------------------------------------------------------
  subroutine base_model_end_model(this)
    type(base_model_t), intent(inout) :: this
    !
    PUSH_SUB(base_model_end_model)
    call hamiltonian_end(this%hm)
    call system_end(this%sys)
    call simulation_end(this%sim)
    nullify(this%config)
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
