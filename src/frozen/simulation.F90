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

#define HASH_TEMPLATE_NAME simulation
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME simulation

module simulation_m

  use global_m
  use messages_m
  use profiling_m

#define LIST_TEMPLATE_NAME simulation
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

  use derivatives_m
  use grid_m
  use json_m
  use space_m
  use geometry_m
  use simul_box_m
  use mesh_m
  use config_dict_m
  use domain_m

#define TEMPLATE_PREFIX simulation
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private
  public ::              &
    simulation__init__,  &
    simulation__start__, &
    simulation__add__,   &
    simulation__copy__,  &
    simulation__end__

  public ::           &
    simulation_new,   &
    simulation_del,   &
    simulation_init,  &
    simulation_start, &
    simulation_set,   &
    simulation_get,   &
    simulation_copy,  & 
    simulation_end
  
#define LIST_TEMPLATE_NAME simulation
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: simulation_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(grid_t),        pointer :: grid   =>null()
    type(simulation_t),  pointer :: prnt   =>null()
    type(domain_t)               :: domain
    type(config_dict_t)          :: dict
    type(simulation_hash_t)      :: hash
    type(simulation_list_t)      :: list
  end type simulation_t

  interface simulation__init__
    module procedure simulation__init__simulation
    module procedure simulation__init__copy
  end interface simulation__init__

  interface simulation_init
    module procedure simulation_init_simulation
    module procedure simulation_init_copy
  end interface simulation_init

  interface simulation_set
    module procedure simulation_set_space
    module procedure simulation_set_geometry
    module procedure simulation_set_grid
  end interface simulation_set

  interface simulation_get
    module procedure simulation_get_simulation_by_config
    module procedure simulation_get_simulation_by_name
    module procedure simulation_get_config
    module procedure simulation_get_info
    module procedure simulation_get_space
    module procedure simulation_get_geometry
    module procedure simulation_get_simul_box
    module procedure simulation_get_mesh
    module procedure simulation_get_derivatives
    module procedure simulation_get_grid
    module procedure simulation_get_domain
  end interface simulation_get

  interface simulation_copy
    module procedure simulation_copy_simulation
  end interface simulation_copy

  interface simulation_end
    module procedure simulation_end_simulation
  end interface simulation_end

  integer, parameter :: SIMULATION_OK          = SIMULATION_HASH_OK
  integer, parameter :: SIMULATION_KEY_ERROR   = SIMULATION_HASH_KEY_ERROR
  integer, parameter :: SIMULATION_EMPTY_ERROR = SIMULATION_HASH_EMPTY_ERROR

#define TEMPLATE_PREFIX simulation
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME simulation
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine simulation_new(this, that)
    type(simulation_t),  target, intent(inout) :: this
    type(simulation_t), pointer                :: that
    !
    PUSH_SUB(simulation_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call simulation_list_push(this%list, that)
    POP_SUB(simulation_new)
    return
  end subroutine simulation_new

  ! ---------------------------------------------------------
  subroutine simulation_del(this)
    type(simulation_t), pointer :: this
    !
    PUSH_SUB(simulation_del)
    if(associated(this))then
      if(associated(this%prnt))then
        call simulation_list_del(this%prnt%list, this)
        call simulation_end(this)
        SAFE_DEALLOCATE_P(this)
      end if
    end if
    POP_SUB(simulation_del)
    return
  end subroutine simulation_del

  ! ---------------------------------------------------------
  subroutine simulation__inull__(this)
    type(simulation_t), intent(inout) :: this
    !
    PUSH_SUB(simulation__inull__)
    nullify(this%config, this%space, this%geo, this%grid, this%prnt)
    POP_SUB(simulation__inull__)
    return
  end subroutine simulation__inull__

  ! ---------------------------------------------------------
  subroutine simulation__iinit__(this, space, config)
    type(simulation_t),          intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(simulation__iinit__)
    call simulation__inull__(this)
    this%config=>config
    this%space=>space
    call config_dict_init(this%dict)
    call simulation_hash_init(this%hash)
    call simulation_list_init(this%list)
    POP_SUB(simulation__iinit__)
    return
  end subroutine simulation__iinit__

  ! ---------------------------------------------------------
  subroutine simulation__init__simulation(this, space, config)
    type(simulation_t),  intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(simulation__init__simulation)
    call simulation__iinit__(this, space, config)
    call domain__init__(this%domain)
    POP_SUB(simulation__init__simulation)
    return
  end subroutine simulation__init__simulation

  ! ---------------------------------------------------------
  subroutine simulation__init__copy(this, that)
    type(simulation_t),  intent(out) :: this
    type(simulation_t),  intent(in)  :: that
    !
    PUSH_SUB(simulation__init__copy)
    if(associated(that%config).and.associated(that%space))then
      call simulation__iinit__(this, that%space, that%config)
      call domain__init__(this%domain)
    end if
    POP_SUB(simulation__init__copy)
    return
  end subroutine simulation__init__copy

  ! ---------------------------------------------------------
  subroutine simulation_init_simulation(this, space, config)
    type(simulation_t),  intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(simulation_init_simulation)
    call simulation__init__(this, space, config)
    POP_SUB(simulation_init_simulation)
    return
  end subroutine simulation_init_simulation

  ! ---------------------------------------------------------
  recursive subroutine simulation_init_copy(this, that)
    type(simulation_t), intent(out) :: this
    type(simulation_t), intent(in)  :: that
    !
    type(simulation_iterator_t)  :: iter
    type(simulation_t),  pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(simulation_init_copy)
    nullify(cnfg, osub, isub)
    call simulation__init__(this, that)
    call simulation_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call simulation_next(iter, cnfg, isub, ierr)
      if(ierr/=SIMULATION_OK)exit
      call simulation_new(this, osub)
      call simulation_init(osub, isub)
      call simulation__add__(this, osub, cnfg)
    end do
    call simulation_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(simulation_init_copy)
    return
  end subroutine simulation_init_copy

  ! ---------------------------------------------------------
  subroutine simulation__istart__(this, grid, geo)
    type(simulation_t),       intent(inout) :: this
    type(grid_t),     target, intent(in)    :: grid
    type(geometry_t), target, intent(in)    :: geo
    !
    PUSH_SUB(simulation__istart__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.associated(this%geo))
    ASSERT(.not.associated(this%grid))
    this%grid=>grid
    this%geo=>geo
    ASSERT(this%grid%sb%dim==this%space%dim)
    POP_SUB(simulation__istart__)
    return
  end subroutine simulation__istart__

  ! ---------------------------------------------------------
  subroutine simulation__start__(this, grid, geo)
    type(simulation_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid
    type(geometry_t),   intent(in)    :: geo
    !
    PUSH_SUB(simulation__start__)
    call simulation__istart__(this, grid, geo)
    call domain__start__(this%domain, grid%sb, geo)
    POP_SUB(simulation__start__)
    return
  end subroutine simulation__start__

  ! ---------------------------------------------------------
  subroutine simulation_start(this, grid, geo)
    type(simulation_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid
    type(geometry_t),   intent(in)    :: geo
    !
    PUSH_SUB(simulation_start)
    call simulation__start__(this, grid, geo)
    POP_SUB(simulation_start)
    return
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation__add__(this, that, config)
    type(simulation_t),  intent(inout) :: this
    type(simulation_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(simulation__add__)
    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call simulation_hash_set(this%hash, config, that)
    call domain__add__(this%domain, that%domain, config)
    POP_SUB(simulation__add__)
    return
  end subroutine simulation__add__

  ! ---------------------------------------------------------
  subroutine simulation_get_simulation_by_config(this, config, that)
    type(simulation_t),  intent(in) :: this
    type(json_object_t), intent(in) :: config
    type(simulation_t), pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(simulation_get_simulation_by_config)
    nullify(that)
    ASSERT(associated(this%config))
    call simulation_hash_get(this%hash, config, that, ierr)
    if(ierr/=SIMULATION_OK)nullify(that)
    POP_SUB(simulation_get_simulation_by_config)
    return
  end subroutine simulation_get_simulation_by_config

  ! ---------------------------------------------------------
  subroutine simulation_get_simulation_by_name(this, name, that)
    type(simulation_t),  intent(in) :: this
    character(len=*),    intent(in) :: name
    type(simulation_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(simulation_get_simulation_by_name)
    nullify(config, that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)&
      call simulation_get(this, config, that)
    POP_SUB(simulation_get_simulation_by_name)
    return
  end subroutine simulation_get_simulation_by_name

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
    ASSERT(.not.associated(this%grid))
    this%grid=>that
    POP_SUB(simulation_set_grid)
    return
  end subroutine simulation_set_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_simulation(this, name, that)
    type(simulation_t),  intent(in) :: this
    character(len=*),    intent(in) :: name
    type(simulation_t), pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(simulation_get_simulation)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call simulation_hash_get(this%hash, config, that, ierr)
      if(ierr/=SIMULATION_OK)nullify(that)
    end if
    POP_SUB(simulation_get_simulation)
    return
  end subroutine simulation_get_simulation

  ! ---------------------------------------------------------
  subroutine simulation_get_info(this, ndim)
    type(simulation_t), intent(in)  :: this
    integer,  optional, intent(out) :: ndim
    !
    PUSH_SUB(simulation_get_info)
    if(present(ndim))then
      ndim=0
      if(associated(this%space))&
        ndim=this%space%dim
    end if
    POP_SUB(simulation_get_info)
    return
  end subroutine simulation_get_info

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
    type(simulation_t), intent(in) :: this
    type(simul_box_t), pointer     :: that
    !
    PUSH_SUB(simulation_get_simul_box)
    nullify(that)
    if(associated(this%grid))&
      that=>this%grid%sb
    POP_SUB(simulation_get_simul_box)
    return
  end subroutine simulation_get_simul_box

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that, fine)
    type(simulation_t), intent(in) :: this
    type(mesh_t),      pointer     :: that
    logical,  optional, intent(in) :: fine
    !
    logical :: fn
    !
    PUSH_SUB(simulation_get_mesh)
    fn=.false.
    nullify(that)
    if(present(fine))fn=fine
    if(associated(this%grid))then
      if(fn)then
        that=>this%grid%fine%mesh
      else
        that=>this%grid%mesh
      end if
    end if
    POP_SUB(simulation_get_mesh)
    return
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_derivatives(this, that, fine)
    type(simulation_t),   intent(in) :: this
    type(derivatives_t), pointer     :: that
    logical,    optional, intent(in) :: fine
    !
    logical :: fn
    !
    PUSH_SUB(simulation_get_derivatives)
    fn=.false.
    nullify(that)
    if(present(fine))fn=fine
    if(associated(this%grid))then
      if(fn)then
        that=>this%grid%fine%der
      else
        that=>this%grid%der
      end if
    end if
    POP_SUB(simulation_get_derivatives)
    return
  end subroutine simulation_get_derivatives

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that
    !
    PUSH_SUB(simulation_get_grid)
    nullify(that)
    if(associated(this%grid))&
      that=>this%grid
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
  subroutine simulation__icopy__(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that
    !
    PUSH_SUB(simulation__icopy__)
    call simulation__iend__(this)
    if(associated(that%config).and.associated(that%space))then
      call simulation__iinit__(this, that%space, that%config)
      if(associated(that%grid).and.associated(that%geo))&
        call simulation__istart__(this, that%grid, that%geo)
    end if
    POP_SUB(simulation__icopy__)
    return
  end subroutine simulation__icopy__

  ! ---------------------------------------------------------
  subroutine simulation__copy__(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that
    !
    PUSH_SUB(simulation__copy__)
    call simulation__icopy__(this, that)
    call domain__copy__(this%domain, that%domain)
    POP_SUB(simulation__copy__)
    return
  end subroutine simulation__copy__

  ! ---------------------------------------------------------
  recursive subroutine simulation_copy_simulation(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that
    !
    type(simulation_iterator_t) :: iter
    type(simulation_t), pointer :: osub, isub
    type(json_object_t),  pointer :: cnfg
    integer                       :: ierr
    !
    PUSH_SUB(simulation_copy_simulation)
    nullify(cnfg, osub, isub)
    call simulation_end(this)
    call simulation__copy__(this, that)
    call simulation_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call simulation_next(iter, cnfg, isub, ierr)
      if(ierr/=SIMULATION_OK)exit
      call simulation_new(this, osub)
      call simulation_copy(osub, isub)
      call simulation__add__(this, osub, cnfg)
    end do
    call simulation_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(simulation_copy_simulation)
    return
  end subroutine simulation_copy_simulation

  ! ---------------------------------------------------------
  subroutine simulation__iend__(this)
    type(simulation_t), intent(inout) :: this
    !
    PUSH_SUB(simulation__iend__)
    call simulation__inull__(this)
    call config_dict_end(this%dict)
    call simulation_hash_end(this%hash)
    call simulation_list_end(this%list)
    POP_SUB(simulation__iend__)
    return
  end subroutine simulation__iend__

  ! ---------------------------------------------------------
  subroutine simulation__end__(this)
    type(simulation_t), intent(inout) :: this
    !
    PUSH_SUB(simulation__end__)
    call simulation__iend__(this)
    call domain__end__(this%domain)
    POP_SUB(simulation__end__)
    return
  end subroutine simulation__end__

  ! ---------------------------------------------------------
  recursive subroutine simulation_end_simulation(this)
    type(simulation_t), intent(inout) :: this
    !
    type(simulation_t), pointer :: subs
    !
    PUSH_SUB(simulation_end_simulation)
    do
      nullify(subs)
      call simulation_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call simulation_end(subs)
      SAFE_DEALLOCATE_P(subs)
    end do
    nullify(subs)
    call simulation__end__(this)
    POP_SUB(simulation_end_simulation)
    return
  end subroutine simulation_end_simulation

#define TEMPLATE_PREFIX simulation
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module simulation_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
