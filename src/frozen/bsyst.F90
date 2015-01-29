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

#define HASH_TEMPLATE_NAME bsyst
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bsyst

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bsyst_m

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

  use bgeom_m, only:           &
    geom_t     => bgeom_t,     &
    geom_init  => bgeom_init,  &
    geom_get   => bgeom_get,   &
    geom_copy  => bgeom_copy,  &
    geom_end   => bgeom_end

  use bstts_m, only:               &
    states_t      => bstts_t,      &
    states_init   => bstts_init,   &
    states_start  => bstts_start,  &
    states_update => bstts_update, &
    states_stop   => bstts_stop,   &
    states_get    => bstts_get,    &
    states_copy   => bstts_copy,   &
    states_end    => bstts_end

  implicit none

  private
  public ::       &
    bsyst_init,   &
    bsyst_start,  &
    bsyst_update, &
    bsyst_stop,   &
    bsyst_set,    &
    bsyst_get,    &
    bsyst_copy,   &
    bsyst_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bsyst_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(space_t)                :: space
    type(geom_t)                 :: geom
    type(states_t)               :: st
    type(bsyst_hash_t)           :: hash
  end type bsyst_t

  type, public :: bsyst_iterator_t
    private
    type(bsyst_t),      pointer :: self =>null()
    type(bsyst_hash_iterator_t) :: iter
  end type bsyst_iterator_t

  interface bsyst_init
    module procedure bsyst_init_begin
    module procedure bsyst_init_build
    module procedure bsyst_init_finish
    module procedure bsyst_iterator_init
  end interface bsyst_init

  interface bsyst_next
    module procedure bsyst_iterator_next_config_bsyst
    module procedure bsyst_iterator_next_config
    module procedure bsyst_iterator_next_bsyst
  end interface bsyst_next

  interface bsyst_set
    module procedure bsyst_set_simulation
  end interface bsyst_set

  interface bsyst_get
    module procedure bsyst_get_config
    module procedure bsyst_get_simulation
    module procedure bsyst_get_space
    module procedure bsyst_get_geom
    module procedure bsyst_get_states
  end interface bsyst_get

  interface bsyst_copy
    module procedure bsyst_copy
    module procedure bsyst_iterator_copy
  end interface bsyst_copy

  interface bsyst_end
    module procedure bsyst_end
    module procedure bsyst_iterator_end
  end interface bsyst_end

  integer, parameter :: BSYST_OK          = BSYST_HASH_OK
  integer, parameter :: BSYST_KEY_ERROR   = BSYST_HASH_KEY_ERROR
  integer, parameter :: BSYST_EMPTY_ERROR = BSYST_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bsyst_init_begin(this, config)
    type(bsyst_t),               intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(bsyst_init_begin)
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
    call bsyst_hash_init(this%hash)
    POP_SUB(bsyst_init_begin)
    return
  end subroutine bsyst_init_begin

  ! ---------------------------------------------------------
  subroutine bsyst_init_build(this, that, config)
    type(bsyst_t),       intent(inout) :: this
    type(bsyst_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bsyst_init_build)
    ASSERT(this%space==that%space)
    call geom_init(this%geom, that%geom, config)
    call states_init(this%st, that%st, config)
    call bsyst_hash_set(this%hash, config, that)
    POP_SUB(bsyst_init_build)
    return
  end subroutine bsyst_init_build

  ! ---------------------------------------------------------
  subroutine bsyst_init_finish(this)
    type(bsyst_t), intent(inout) :: this
    !
    PUSH_SUB(bsyst_init_finish)
    call geom_init(this%geom)
    POP_SUB(bsyst_init_finish)
    return
  end subroutine bsyst_init_finish

  ! ---------------------------------------------------------
  subroutine bsyst_start(this, sim)
    type(bsyst_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(bsyst_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call states_start(this%st, sim)
    POP_SUB(bsyst_start)
    return
  end subroutine bsyst_start

  ! ---------------------------------------------------------
  subroutine bsyst_update(this)
    type(bsyst_t), intent(inout) :: this
    !
    PUSH_SUB(bsyst_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call states_update(this%st)
    POP_SUB(bsyst_update)
    return
  end subroutine bsyst_update

  ! ---------------------------------------------------------
  subroutine bsyst_stop(this)
    type(bsyst_t), intent(inout) :: this
    !
    PUSH_SUB(bsyst_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call states_stop(this%st)
    POP_SUB(bsyst_stop)
    return
  end subroutine bsyst_stop

  ! ---------------------------------------------------------
  subroutine bsyst_set_simulation(this, that)
    type(bsyst_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that
    !
    PUSH_SUB(bsyst_set_simulation)
    ASSERT(.not.associated(this%sim))
    this%sim=>that
    POP_SUB(bsyst_set_simulation)
    return
  end subroutine bsyst_set_simulation

  ! ---------------------------------------------------------
  subroutine bsyst_get_config(this, that)
    type(bsyst_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bsyst_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bsyst_get_config)
    return
  end subroutine bsyst_get_config

  ! ---------------------------------------------------------
  subroutine bsyst_get_simulation(this, that)
    type(bsyst_t),       target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(bsyst_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(bsyst_get_simulation)
    return
  end subroutine bsyst_get_simulation

  ! ---------------------------------------------------------
  subroutine bsyst_get_space(this, that)
    type(bsyst_t),  target, intent(in) :: this
    type(space_t), pointer             :: that
    !
    PUSH_SUB(bsyst_get_space)
    that=>this%space
    POP_SUB(bsyst_get_space)
    return
  end subroutine bsyst_get_space

  ! ---------------------------------------------------------
  subroutine bsyst_get_geom(this, that)
    type(bsyst_t), target, intent(in) :: this
    type(geom_t), pointer             :: that
    !
    PUSH_SUB(bsyst_get_geom)
    that=>this%geom
    POP_SUB(bsyst_get_geom)
    return
  end subroutine bsyst_get_geom

  ! ---------------------------------------------------------
  subroutine bsyst_get_states(this, that)
    type(bsyst_t),   target, intent(in) :: this
    type(states_t), pointer             :: that
    !
    PUSH_SUB(bsyst_get_states)
    that=>this%st
    POP_SUB(bsyst_get_states)
    return
  end subroutine bsyst_get_states

  ! ---------------------------------------------------------
  subroutine bsyst_copy(this, that)
    type(bsyst_t), intent(out) :: this
    type(bsyst_t), intent(in)  :: that
    !
    PUSH_SUB(bsyst_copy)
    this%config=>that%config
    this%sim=>that%sim
    call geom_copy(this%geom, that%geom)
    call space_copy(this%space, that%space)
    call states_copy(this%st, that%st)
    call bsyst_hash_copy(this%hash, that%hash)
    POP_SUB(bsyst_copy)
    return
  end subroutine bsyst_copy

  ! ---------------------------------------------------------
  subroutine bsyst_end(this)
    type(bsyst_t), intent(inout) :: this
    !
    PUSH_SUB(bsyst_end)
    nullify(this%config, this%sim)
    call space_end(this%space)
    call geom_end(this%geom)
    call states_end(this%st)
    call bsyst_hash_end(this%hash)
    POP_SUB(bsyst_end)
    return
  end subroutine bsyst_end

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_init(this, that)
    type(bsyst_iterator_t), intent(out) :: this
    type(bsyst_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bsyst_iterator_init)
    this%self=>that
    call bsyst_hash_init(this%iter, that%hash)
    POP_SUB(bsyst_iterator_init)
    return
  end subroutine bsyst_iterator_init

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_next_config_bsyst(this, config, system, ierr)
    type(bsyst_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bsyst_t),         pointer        :: system
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bsyst_iterator_next_config_bsyst)
    call bsyst_hash_next(this%iter, config, system, ierr)
    POP_SUB(bsyst_iterator_next_config_bsyst)
    return
  end subroutine bsyst_iterator_next_config_bsyst

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_next_config(this, that, ierr)
    type(bsyst_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bsyst_iterator_next_config)
    call bsyst_hash_next(this%iter, that, ierr)
    POP_SUB(bsyst_iterator_next_config)
    return
  end subroutine bsyst_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_next_bsyst(this, that, ierr)
    type(bsyst_iterator_t), intent(inout) :: this
    type(bsyst_t),         pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bsyst_iterator_next_bsyst)
    call bsyst_hash_next(this%iter, that, ierr)
    POP_SUB(bsyst_iterator_next_bsyst)
    return
  end subroutine bsyst_iterator_next_bsyst

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_copy(this, that)
    type(bsyst_iterator_t), intent(inout) :: this
    type(bsyst_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bsyst_iterator_copy)
    this%self=>that%self
    call bsyst_hash_copy(this%iter, that%iter)
    POP_SUB(bsyst_iterator_copy)
    return
  end subroutine bsyst_iterator_copy

  ! ---------------------------------------------------------
  subroutine bsyst_iterator_end(this)
    type(bsyst_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bsyst_iterator_end)
    nullify(this%self)
    call bsyst_hash_end(this%iter)
    POP_SUB(bsyst_iterator_end)
    return
  end subroutine bsyst_iterator_end

end module bsyst_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
