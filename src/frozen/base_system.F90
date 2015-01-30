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

#define HASH_TEMPLATE_NAME base_system
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_system

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_system_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

  use simulation_m, only: &
    simulation_t

  use space_m, only: &
    operator(==),    &
    space_t,         &
    space_init,      &
    space_copy,      &
    space_end

  use geometry_m, only: &
    geometry_t

  use base_geom_m, only:           &
    geom_t     => base_geom_t,     &
    geom_init  => base_geom_init,  &
    geom_get   => base_geom_get,   &
    geom_copy  => base_geom_copy,  &
    geom_end   => base_geom_end

  use base_states_m, only:               &
    states_t      => base_states_t,      &
    states_init   => base_states_init,   &
    states_start  => base_states_start,  &
    states_update => base_states_update, &
    states_stop   => base_states_stop,   &
    states_get    => base_states_get,    &
    states_copy   => base_states_copy,   &
    states_end    => base_states_end

  implicit none

  private
  public ::             &
    base_system_init,   &
    base_system_start,  &
    base_system_update, &
    base_system_stop,   &
    base_system_set,    &
    base_system_get,    &
    base_system_copy,   &
    base_system_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_system_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(space_t)                :: space
    type(geom_t)                 :: geom
    type(states_t)               :: st
    type(base_system_hash_t)     :: hash
  end type base_system_t

  type, public :: base_system_iterator_t
    private
    type(base_system_t),      pointer :: self =>null()
    type(base_system_hash_iterator_t) :: iter
  end type base_system_iterator_t

  interface base_system_init
    module procedure base_system_init_begin
    module procedure base_system_init_build
    module procedure base_system_init_finish
    module procedure base_system_iterator_init
  end interface base_system_init

  interface base_system_next
    module procedure base_system_iterator_next_config_system
    module procedure base_system_iterator_next_config
    module procedure base_system_iterator_next_system
  end interface base_system_next

  interface base_system_set
    module procedure base_system_set_simulation
  end interface base_system_set

  interface base_system_get
    module procedure base_system_get_config
    module procedure base_system_get_simulation
    module procedure base_system_get_space
    module procedure base_system_get_geom
    module procedure base_system_get_states
  end interface base_system_get

  interface base_system_copy
    module procedure base_system_copy_system
    module procedure base_system_iterator_copy
  end interface base_system_copy

  interface base_system_end
    module procedure base_system_end_system
    module procedure base_system_iterator_end
  end interface base_system_end

  integer, parameter :: BASE_SYSTEM_OK          = BASE_SYSTEM_HASH_OK
  integer, parameter :: BASE_SYSTEM_KEY_ERROR   = BASE_SYSTEM_HASH_KEY_ERROR
  integer, parameter :: BASE_SYSTEM_EMPTY_ERROR = BASE_SYSTEM_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_system_init_begin(this, config)
    type(base_system_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_system_init_begin)
    this%config=>config
    nullify(this%sim, cnfg)
    call json_get(this%config, "space", cnfg, ierr)
    if(ierr==JSON_OK)call space_init(this%space, cnfg)
    nullify(cnfg)
    call json_get(this%config, "geometry", cnfg, ierr)
    if(ierr==JSON_OK)call geom_init(this%geom, this%space, cnfg)
    nullify(cnfg)
    call json_get(this%config, "states", cnfg, ierr)
    if(ierr==JSON_OK)call states_init(this%st, cnfg)
    nullify(cnfg)
    call base_system_hash_init(this%hash)
    POP_SUB(base_system_init_begin)
    return
  end subroutine base_system_init_begin

  ! ---------------------------------------------------------
  subroutine base_system_init_build(this, that, config)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(base_system_init_build)
    ASSERT(this%space==that%space)
    call geom_init(this%geom, that%geom, config)
    call states_init(this%st, that%st, config)
    call base_system_hash_set(this%hash, config, that)
    POP_SUB(base_system_init_build)
    return
  end subroutine base_system_init_build

  ! ---------------------------------------------------------
  subroutine base_system_init_finish(this)
    type(base_system_t), intent(inout) :: this
    !
    PUSH_SUB(base_system_init_finish)
    call geom_init(this%geom)
    POP_SUB(base_system_init_finish)
    return
  end subroutine base_system_init_finish

  ! ---------------------------------------------------------
  subroutine base_system_start(this, sim)
    type(base_system_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_system_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call states_start(this%st, sim)
    POP_SUB(base_system_start)
    return
  end subroutine base_system_start

  ! ---------------------------------------------------------
  subroutine base_system_update(this)
    type(base_system_t), intent(inout) :: this
    !
    PUSH_SUB(base_system_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call states_update(this%st)
    POP_SUB(base_system_update)
    return
  end subroutine base_system_update

  ! ---------------------------------------------------------
  subroutine base_system_stop(this)
    type(base_system_t), intent(inout) :: this
    !
    PUSH_SUB(base_system_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call states_stop(this%st)
    POP_SUB(base_system_stop)
    return
  end subroutine base_system_stop

  ! ---------------------------------------------------------
  subroutine base_system_set_simulation(this, that)
    type(base_system_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_system_set_simulation)
    ASSERT(.not.associated(this%sim))
    this%sim=>that
    POP_SUB(base_system_set_simulation)
    return
  end subroutine base_system_set_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_config(this, that)
    type(base_system_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_system_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_system_get_config)
    return
  end subroutine base_system_get_config

  ! ---------------------------------------------------------
  subroutine base_system_get_simulation(this, that)
    type(base_system_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(base_system_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_system_get_simulation)
    return
  end subroutine base_system_get_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_space(this, that)
    type(base_system_t), target, intent(in) :: this
    type(space_t),      pointer             :: that
    !
    PUSH_SUB(base_system_get_space)
    that=>this%space
    POP_SUB(base_system_get_space)
    return
  end subroutine base_system_get_space

  ! ---------------------------------------------------------
  subroutine base_system_get_geom(this, that)
    type(base_system_t), target, intent(in) :: this
    type(geom_t),       pointer             :: that
    !
    PUSH_SUB(base_system_get_geom)
    that=>this%geom
    POP_SUB(base_system_get_geom)
    return
  end subroutine base_system_get_geom

  ! ---------------------------------------------------------
  subroutine base_system_get_states(this, that)
    type(base_system_t), target, intent(in) :: this
    type(states_t),     pointer             :: that
    !
    PUSH_SUB(base_system_get_states)
    that=>this%st
    POP_SUB(base_system_get_states)
    return
  end subroutine base_system_get_states

  ! ---------------------------------------------------------
  subroutine base_system_copy_system(this, that)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that
    !
    PUSH_SUB(base_system_copy_system)
    this%config=>that%config
    this%sim=>that%sim
    call geom_copy(this%geom, that%geom)
    call space_copy(this%space, that%space)
    call states_copy(this%st, that%st)
    call base_system_hash_copy(this%hash, that%hash)
    POP_SUB(base_system_copy_system)
    return
  end subroutine base_system_copy_system

  ! ---------------------------------------------------------
  subroutine base_system_end_system(this)
    type(base_system_t), intent(inout) :: this
    !
    PUSH_SUB(base_system_end_system)
    nullify(this%config, this%sim)
    call space_end(this%space)
    call geom_end(this%geom)
    call states_end(this%st)
    call base_system_hash_end(this%hash)
    POP_SUB(base_system_end_system)
    return
  end subroutine base_system_end_system

  ! ---------------------------------------------------------
  subroutine base_system_iterator_init(this, that)
    type(base_system_iterator_t), intent(out) :: this
    type(base_system_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_system_iterator_init)
    this%self=>that
    call base_system_hash_init(this%iter, that%hash)
    POP_SUB(base_system_iterator_init)
    return
  end subroutine base_system_iterator_init

  ! ---------------------------------------------------------
  subroutine base_system_iterator_next_config_system(this, config, system, ierr)
    type(base_system_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: config
    type(base_system_t),         pointer        :: system
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_system_iterator_next_config_system)
    call base_system_hash_next(this%iter, config, system, ierr)
    POP_SUB(base_system_iterator_next_config_system)
    return
  end subroutine base_system_iterator_next_config_system

  ! ---------------------------------------------------------
  subroutine base_system_iterator_next_config(this, that, ierr)
    type(base_system_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_system_iterator_next_config)
    call base_system_hash_next(this%iter, that, ierr)
    POP_SUB(base_system_iterator_next_config)
    return
  end subroutine base_system_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_system_iterator_next_system(this, that, ierr)
    type(base_system_iterator_t), intent(inout) :: this
    type(base_system_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_system_iterator_next_system)
    call base_system_hash_next(this%iter, that, ierr)
    POP_SUB(base_system_iterator_next_system)
    return
  end subroutine base_system_iterator_next_system

  ! ---------------------------------------------------------
  subroutine base_system_iterator_copy(this, that)
    type(base_system_iterator_t), intent(inout) :: this
    type(base_system_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_system_iterator_copy)
    this%self=>that%self
    call base_system_hash_copy(this%iter, that%iter)
    POP_SUB(base_system_iterator_copy)
    return
  end subroutine base_system_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_system_iterator_end(this)
    type(base_system_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_system_iterator_end)
    nullify(this%self)
    call base_system_hash_end(this%iter)
    POP_SUB(base_system_iterator_end)
    return
  end subroutine base_system_iterator_end

end module base_system_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
