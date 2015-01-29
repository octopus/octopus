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

#define HASH_TEMPLATE_NAME bmodl
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bmodl

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bmodl_m

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

  use bgeom_m, only:       &
    geom_t   => bgeom_t,   &
    geom_get => bgeom_get

  use bsyst_m, only:                  &
    system_t      => bsyst_t,         &
    system_init   => bsyst_init,      &
    system_start  => bsyst_start,     &
    system_update => bsyst_update,    &
    system_stop   => bsyst_stop,      &
    system_get    => bsyst_get,       &
    system_copy   => bsyst_copy,      &
    system_end    => bsyst_end

  use bhmlt_m, only:                       &
    hamiltonian_t      => bhmlt_t,         &
    hamiltonian_init   => bhmlt_init,      &
    hamiltonian_start  => bhmlt_start,     &
    hamiltonian_update => bhmlt_update,    &
    hamiltonian_stop   => bhmlt_stop,      &
    hamiltonian_copy   => bhmlt_copy,      &
    hamiltonian_end    => bhmlt_end

  implicit none

  private
  public ::       &
    bmodl_init,   &
    bmodl_start,  &
    bmodl_update, &
    bmodl_stop,   &
    bmodl_get,    &
    bmodl_copy,   &
    bmodl_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bmodl_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t)           :: sim
    type(system_t)               :: sys
    type(hamiltonian_t)          :: hm
    type(bmodl_hash_t)           :: hash
  end type bmodl_t

  type, public :: bmodl_iterator_t
    private
    type(bmodl_t),      pointer :: self =>null()
    type(bmodl_hash_iterator_t) :: iter
  end type bmodl_iterator_t

  interface bmodl_init
    module procedure bmodl_init_begin
    module procedure bmodl_init_build
    module procedure bmodl_init_finish
    module procedure bmodl_iterator_init
  end interface bmodl_init

  interface bmodl_get
    module procedure bmodl_get_config
    module procedure bmodl_get_simulation
    module procedure bmodl_get_system
    module procedure bmodl_get_hamiltonian
  end interface bmodl_get

  interface bmodl_next
    module procedure bmodl_iterator_next_config_bmodl
    module procedure bmodl_iterator_next_config
    module procedure bmodl_iterator_next_bmodl
  end interface bmodl_next

  interface bmodl_copy
    module procedure bmodl_copy_bmodl
    module procedure bmodl_iterator_copy
  end interface bmodl_copy

  interface bmodl_end
    module procedure bmodl_end_bmodl
    module procedure bmodl_iterator_end
  end interface bmodl_end

  integer, parameter :: BMODL_OK          = BMODL_HASH_OK
  integer, parameter :: BMODL_KEY_ERROR   = BMODL_HASH_KEY_ERROR
  integer, parameter :: BMODL_EMPTY_ERROR = BMODL_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bmodl_init_begin(this, config)
    type(bmodl_t),               intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(bmodl_init_begin)
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
    call bmodl_hash_init(this%hash)
    POP_SUB(bmodl_init_begin)
    return
  end subroutine bmodl_init_begin

  ! ---------------------------------------------------------
  subroutine bmodl_init_build(this, that, config)
    type(bmodl_t),       intent(inout) :: this
    type(bmodl_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bmodl_init_build)
    call bmodl_hash_set(this%hash, config, that)
    call simulation_init(this%sim, that%sim, config)
    call system_init(this%sys, that%sys, config)
    call hamiltonian_init(this%hm, that%hm, config)
    POP_SUB(bmodl_init_build)
    return
  end subroutine bmodl_init_build

  ! ---------------------------------------------------------
  subroutine bmodl_init_finish(this)
    type(bmodl_t), intent(inout) :: this
    !
    PUSH_SUB(bmodl_init_finish)
    call system_init(this%sys)
    POP_SUB(bmodl_init_finish)
    return
  end subroutine bmodl_init_finish

  ! ---------------------------------------------------------
  subroutine bmodl_start(this, grid)
    type(bmodl_t), intent(inout) :: this
    type(grid_t),  intent(in)    :: grid
    !
    type(geom_t),     pointer :: geom
    type(geometry_t), pointer :: geo
    type(space_t),    pointer :: space
    !
    PUSH_SUB(bmodl_start)
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
    POP_SUB(bmodl_start)
    return
  end subroutine bmodl_start

  ! ---------------------------------------------------------
  subroutine bmodl_update(this)
    type(bmodl_t), intent(inout) :: this
    !
    PUSH_SUB(bmodl_update)
    call system_update(this%sys)
    call hamiltonian_update(this%hm)
    POP_SUB(bmodl_update)
    return
  end subroutine bmodl_update

  ! ---------------------------------------------------------
  subroutine bmodl_stop(this)
    type(bmodl_t), intent(inout) :: this
    !
    PUSH_SUB(bmodl_stop)
    call system_stop(this%sys)
    call hamiltonian_stop(this%hm)
    POP_SUB(bmodl_stop)
    return
  end subroutine bmodl_stop

  ! ---------------------------------------------------------
  subroutine bmodl_get_config(this, that)
    type(bmodl_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bmodl_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bmodl_get_config)
    return
  end subroutine bmodl_get_config

  ! ---------------------------------------------------------
  subroutine bmodl_get_simulation(this, that)
    type(bmodl_t),       target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(bmodl_get_simulation)
    that=>this%sim
    POP_SUB(bmodl_get_simulation)
    return
  end subroutine bmodl_get_simulation

  ! ---------------------------------------------------------
  subroutine bmodl_get_system(this, that)
    type(bmodl_t),   target, intent(in) :: this
    type(system_t), pointer             :: that
    !
    PUSH_SUB(bmodl_get_system)
    that=>this%sys
    POP_SUB(bmodl_get_system)
    return
  end subroutine bmodl_get_system

  ! ---------------------------------------------------------
  subroutine bmodl_get_hamiltonian(this, that)
    type(bmodl_t),        target, intent(in) :: this
    type(hamiltonian_t), pointer             :: that
    !
    PUSH_SUB(bmodl_get_hamiltonian)
    that=>this%hm
    POP_SUB(bmodl_get_hamiltonian)
    return
  end subroutine bmodl_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine bmodl_copy_bmodl(this, that)
    type(bmodl_t), intent(out) :: this
    type(bmodl_t), intent(in)  :: that
    !
    PUSH_SUB(bmodl_copy_bmodl)
    this%config=>that%config
    call simulation_copy(this%sim, that%sim)
    call system_copy(this%sys, that%sys)
    call hamiltonian_copy(this%hm, that%hm)
    call bmodl_hash_copy(this%hash, that%hash)
    POP_SUB(bmodl_copy_bmodl)
    return
  end subroutine bmodl_copy_bmodl

  ! ---------------------------------------------------------
  subroutine bmodl_end_bmodl(this)
    type(bmodl_t), intent(inout) :: this
    !
    PUSH_SUB(bmodl_end_bmodl)
    call hamiltonian_end(this%hm)
    call system_end(this%sys)
    call simulation_end(this%sim)
    nullify(this%config)
    call bmodl_hash_end(this%hash)
    POP_SUB(bmodl_end_bmodl)
    return
  end subroutine bmodl_end_bmodl

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_init(this, that)
    type(bmodl_iterator_t), intent(out) :: this
    type(bmodl_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bmodl_iterator_init)
    this%self=>that
    call bmodl_hash_init(this%iter, that%hash)
    POP_SUB(bmodl_iterator_init)
    return
  end subroutine bmodl_iterator_init

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_next_config_bmodl(this, config, bmodl, ierr)
    type(bmodl_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bmodl_t),         pointer        :: bmodl
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bmodl_iterator_next_config_bmodl)
    call bmodl_hash_next(this%iter, config, bmodl, ierr)
    POP_SUB(bmodl_iterator_next_config_bmodl)
    return
  end subroutine bmodl_iterator_next_config_bmodl

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_next_config(this, that, ierr)
    type(bmodl_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bmodl_iterator_next_config)
    call bmodl_hash_next(this%iter, that, ierr)
    POP_SUB(bmodl_iterator_next_config)
    return
  end subroutine bmodl_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_next_bmodl(this, that, ierr)
    type(bmodl_iterator_t), intent(inout) :: this
    type(bmodl_t),         pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bmodl_iterator_next_bmodl)
    call bmodl_hash_next(this%iter, that, ierr)
    POP_SUB(bmodl_iterator_next_bmodl)
    return
  end subroutine bmodl_iterator_next_bmodl

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_copy(this, that)
    type(bmodl_iterator_t), intent(inout) :: this
    type(bmodl_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bmodl_iterator_copy)
    this%self=>that%self
    call bmodl_hash_copy(this%iter, that%iter)
    POP_SUB(bmodl_iterator_copy)
    return
  end subroutine bmodl_iterator_copy

  ! ---------------------------------------------------------
  subroutine bmodl_iterator_end(this)
    type(bmodl_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bmodl_iterator_end)
    nullify(this%self)
    call bmodl_hash_end(this%iter)
    POP_SUB(bmodl_iterator_end)
    return
  end subroutine bmodl_iterator_end

end module bmodl_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
