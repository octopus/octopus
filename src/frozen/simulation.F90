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

#define HASH_TEMPLATE_NAME simulation
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME simulation

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module simulation_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: operator(==), json_hash
  use kinds_m, only: wp

  use json_m,   only: json_object_t, json_get

  use space_m, only: &
    space_t

  use geometry_m, only: &
    geometry_t

  use simul_box_m, only: &
    simul_box_t

  use mesh_m, only: &
    mesh_t

  use igrid_m, only: &
    grid_t,          &
    grid_get

  use domain_m, only: &
    domain_t,         &
    domain_init,      &
    domain_start,     &
    domain_copy,      &
    domain_end

  implicit none

  private
  public ::           &
    simulation_init,  &
    simulation_start, &
    simulation_get,   &
    simulation_set,   &
    simulation_copy,  &
    simulation_end
  
#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: simulation_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(grid_t),        pointer :: gr     =>null()
    type(domain_t)               :: domain
    type(simulation_hash_t)      :: hash
  end type simulation_t

  type, public :: simulation_iterator_t
    private
    type(simulation_t),      pointer :: self =>null()
    type(simulation_hash_iterator_t) :: iter
  end type simulation_iterator_t

  interface simulation_init
    module procedure simulation_init_start
    module procedure simulation_init_build
    module procedure simulation_iterator_init
  end interface simulation_init

  interface simulation_set
    module procedure simulation_set_space
    module procedure simulation_set_geometry
    module procedure simulation_set_grid
  end interface simulation_set

  interface simulation_get
    module procedure simulation_get_config
    module procedure simulation_get_space
    module procedure simulation_get_geometry
    module procedure simulation_get_simul_box
    module procedure simulation_get_mesh
    module procedure simulation_get_grid
    module procedure simulation_get_domain
  end interface simulation_get

  interface simulation_next
    module procedure simulation_iterator_next_config_simulation
    module procedure simulation_iterator_next_config
    module procedure simulation_iterator_next_simulation
  end interface simulation_next

  interface simulation_copy
    module procedure simulation_copy_simulation
    module procedure simulation_iterator_copy
  end interface simulation_copy

  interface simulation_end
    module procedure simulation_end_simulation
    module procedure simulation_iterator_end
  end interface simulation_end

  integer, parameter :: SIMULATION_OK          = SIMULATION_HASH_OK
  integer, parameter :: SIMULATION_KEY_ERROR   = SIMULATION_HASH_KEY_ERROR
  integer, parameter :: SIMULATION_EMPTY_ERROR = SIMULATION_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine simulation_init_start(this, config)
    type(simulation_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(simulation_init_start)
    nullify(this%config, this%gr)
    this%config=>config
    call domain_init(this%domain)
    call simulation_hash_init(this%hash)
    POP_SUB(simulation_init_start)
    return
  end subroutine simulation_init_start

  ! ---------------------------------------------------------
  subroutine simulation_init_build(this, that, config)
    type(simulation_t),  intent(inout) :: this
    type(simulation_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(simulation_init_build)
    ASSERT(associated(this%config))
    call domain_init(this%domain, that%domain)
    call simulation_hash_set(this%hash, config, that)
    POP_SUB(simulation_init_build)
    return
  end subroutine simulation_init_build

  ! ---------------------------------------------------------
  subroutine simulation_start(this, grid, geo, space)
    type(simulation_t),       intent(inout) :: this
    type(grid_t),     target, intent(in)    :: grid
    type(geometry_t), target, intent(in)    :: geo
    type(space_t),    target, intent(in)    :: space
    !
    PUSH_SUB(simulation_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%gr))
    ASSERT(.not.associated(this%geo))
    ASSERT(.not.associated(this%space))
    this%gr=>grid
    this%geo=>geo
    this%space=>space
    ASSERT(this%gr%sb%dim==this%space%dim)
    call domain_start(this%domain, this%gr%sb, this%geo)
    POP_SUB(simulation_start)
    return
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation_set_space(this, that)
    type(simulation_t),    intent(inout) :: this
    type(space_t), target, intent(in)    :: that
    !
    PUSH_SUB(simulation_set_space)
    ASSERT(.not.associated(this%space))
    this%space=>that
    POP_SUB(simulation_set_space)
    return
  end subroutine simulation_set_space

  ! ---------------------------------------------------------
  subroutine simulation_set_geometry(this, that)
    type(simulation_t),       intent(inout) :: this
    type(geometry_t), target, intent(in)    :: that
    !
    PUSH_SUB(simulation_set_geometry)
    ASSERT(.not.associated(this%geo))
    this%geo=>that
    POP_SUB(simulation_set_geometry)
    return
  end subroutine simulation_set_geometry

  ! ---------------------------------------------------------
  subroutine simulation_set_grid(this, that)
    type(simulation_t),   intent(inout) :: this
    type(grid_t), target, intent(in)    :: that
    !
    PUSH_SUB(simulation_set_grid)
    ASSERT(.not.associated(this%gr))
    this%gr=>that
    POP_SUB(simulation_set_grid)
    return
  end subroutine simulation_set_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_config(this, that)
    type(simulation_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(simulation_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(simulation_get_config)
    return
  end subroutine simulation_get_config

  ! ---------------------------------------------------------
  subroutine simulation_get_space(this, that)
    type(simulation_t), target, intent(in) :: this
    type(space_t),     pointer             :: that
    !
    PUSH_SUB(simulation_get_space)
    nullify(that)
    if(associated(this%space))&
      that=>this%space
    POP_SUB(simulation_get_space)
    return
  end subroutine simulation_get_space

  ! ---------------------------------------------------------
  subroutine simulation_get_geometry(this, that)
    type(simulation_t), target, intent(in) :: this
    type(geometry_t),  pointer             :: that
    !
    PUSH_SUB(simulation_get_geometry)
    nullify(that)
    if(associated(this%geo))&
      that=>this%geo
    POP_SUB(simulation_get_geometry)
    return
  end subroutine simulation_get_geometry

  ! ---------------------------------------------------------
  subroutine simulation_get_simul_box(this, that)
    type(simulation_t), target, intent(in) :: this
    type(simul_box_t), pointer             :: that
    !
    PUSH_SUB(simulation_get_simul_box)
    nullify(that)
    if(associated(this%gr))&
      call grid_get(this%gr, that)
    POP_SUB(simulation_get_simul_box)
    return
  end subroutine simulation_get_simul_box

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that)
    type(simulation_t), target, intent(in) :: this
    type(mesh_t),      pointer             :: that
    !
    PUSH_SUB(simulation_get_mesh)
    nullify(that)
    if(associated(this%gr))&
      call grid_get(this%gr, that)
    POP_SUB(simulation_get_mesh)
    return
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that
    !
    PUSH_SUB(simulation_get_grid)
    nullify(that)
    if(associated(this%gr))&
      that=>this%gr
    POP_SUB(simulation_get_grid)
    return
  end subroutine simulation_get_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_domain(this, that)
    type(simulation_t), target, intent(in) :: this
    type(domain_t),    pointer             :: that
    !
    PUSH_SUB(simulation_get_domain)
    that=>this%domain
    POP_SUB(simulation_get_domain)
    return
  end subroutine simulation_get_domain

  ! ---------------------------------------------------------
  subroutine simulation_copy_simulation(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that
    !
    PUSH_SUB(simulation_copy_simulation)
    this%config=>that%config
    this%space=>that%space
    this%geo=>that%geo
    this%gr=>that%gr
    call domain_copy(this%domain, that%domain)
    call simulation_hash_copy(this%hash, that%hash)
    POP_SUB(simulation_copy_simulation)
    return
  end subroutine simulation_copy_simulation

  ! ---------------------------------------------------------
  subroutine simulation_end_simulation(this)
    type(simulation_t), intent(inout) :: this
    !
    PUSH_SUB(simulation_end_simulation)
    nullify(this%config, this%space, this%geo, this%gr)
    call domain_end(this%domain)
    call simulation_hash_end(this%hash)
    POP_SUB(simulation_end_simulation)
    return
  end subroutine simulation_end_simulation

  ! ---------------------------------------------------------
  subroutine simulation_iterator_init(this, that)
    type(simulation_iterator_t), intent(out) :: this
    type(simulation_t),  target, intent(in)  :: that
    !
    PUSH_SUB(simulation_iterator_init)
    this%self=>that
    call simulation_hash_init(this%iter, that%hash)
    POP_SUB(simulation_iterator_init)
    return
  end subroutine simulation_iterator_init

  ! ---------------------------------------------------------
  subroutine simulation_iterator_next_config_simulation(this, config, sim, ierr)
    type(simulation_iterator_t), intent(inout) :: this
    type(json_object_t),        pointer        :: config
    type(simulation_t),         pointer        :: sim
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(simulation_iterator_next_config_simulation)
    call simulation_hash_next(this%iter, config, sim, ierr)
    POP_SUB(simulation_iterator_next_config_simulation)
    return
  end subroutine simulation_iterator_next_config_simulation

  ! ---------------------------------------------------------
  subroutine simulation_iterator_next_config(this, that, ierr)
    type(simulation_iterator_t), intent(inout) :: this
    type(json_object_t),        pointer        :: that
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(simulation_iterator_next_config)
    call simulation_hash_next(this%iter, that, ierr)
    POP_SUB(simulation_iterator_next_config)
    return
  end subroutine simulation_iterator_next_config

  ! ---------------------------------------------------------
  subroutine simulation_iterator_next_simulation(this, that, ierr)
    type(simulation_iterator_t), intent(inout) :: this
    type(simulation_t),         pointer        :: that
    integer,           optional, intent(out)   :: ierr
    !
    PUSH_SUB(simulation_iterator_next_simulation)
    call simulation_hash_next(this%iter, that, ierr)
    POP_SUB(simulation_iterator_next_simulation)
    return
  end subroutine simulation_iterator_next_simulation

  ! ---------------------------------------------------------
  subroutine simulation_iterator_copy(this, that)
    type(simulation_iterator_t), intent(inout) :: this
    type(simulation_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(simulation_iterator_copy)
    this%self=>that%self
    call simulation_hash_copy(this%iter, that%iter)
    POP_SUB(simulation_iterator_copy)
    return
  end subroutine simulation_iterator_copy

  ! ---------------------------------------------------------
  subroutine simulation_iterator_end(this)
    type(simulation_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(simulation_iterator_end)
    nullify(this%self)
    call simulation_hash_end(this%iter)
    POP_SUB(simulation_iterator_end)
    return
  end subroutine simulation_iterator_end

end module simulation_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
