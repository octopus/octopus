#include "global.h"

module simulation_m

  use global_m
  use messages_m
  use profiling_m

  use basis_m,     only: basis_t
  use domain_m,    only: domain_t, domain_init, domain_copy, domain_end
  use geometry_m,  only: geometry_t
  use igrid_m,     only: grid_t, grid_init, grid_copy, grid_end
  use json_m,      only: JSON_OK, json_object_t, json_get
  use mesh_m,      only: mesh_t 
  use simul_box_m, only: simul_box_t
  use space_m,     only: space_t

  implicit none

  private
  public ::              &
    simulation_init,     &
    simulation_start,    &
    simulation_extend,   &
    simulation_get,      &
    simulation_get_ndim, &
    simulation_copy,     &
    simulation_end
  
  type, public :: simulation_t
    private
    type(json_object_t), pointer :: config =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(space_t),       pointer :: space  =>null()
    type(grid_t),        pointer :: gr     =>null()
    logical                      :: init   = .false.
    logical                      :: start  = .false.
    logical                      :: alloc  = .false.
    type(domain_t)               :: domain
  end type simulation_t

  interface simulation_init
    module procedure simulation_init_grid
    module procedure simulation_init_config
  end interface simulation_init

  interface simulation_get
    module procedure simulation_get_config
    module procedure simulation_get_simul_box
    module procedure simulation_get_grid
    module procedure simulation_get_mesh
    module procedure simulation_get_geometry
    module procedure simulation_get_domain
  end interface simulation_get

contains

  ! ---------------------------------------------------------
  subroutine simulation_init_common(this, geo, space)
    type(simulation_t),       intent(out) :: this
    type(geometry_t), target, intent(in)  :: geo
    type(space_t),    target, intent(in)  :: space
    !
    PUSH_SUB(simulation_init_common)
    this%config=>null()
    this%geo=>geo
    this%space=>space
    this%gr=>null()
    this%init=.true.
    this%start=.false.
    this%alloc=.false.
    POP_SUB(simulation_init_common)
    return
  end subroutine simulation_init_common

  ! ---------------------------------------------------------
  subroutine simulation_init_grid(this, geo, space, grid)
    type(simulation_t),   intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(space_t),        intent(in)  :: space
    type(grid_t), target, intent(in)  :: grid
    !
    PUSH_SUB(simulation_init_grid)
    call simulation_init_common(this, geo, space)
    this%gr=>grid
    this%alloc=.false.
    POP_SUB(simulation_init_grid)
    return
  end subroutine simulation_init_grid

  ! ---------------------------------------------------------
  subroutine simulation_init_config(this, geo, space, config)
    type(simulation_t),          intent(out) :: this
    type(geometry_t),            intent(in)  :: geo
    type(space_t),               intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(simulation_init_config)
    call simulation_init_common(this, geo, space)
    this%config=>config
    this%alloc=.true.
    POP_SUB(simulation_init_config)
    return
  end subroutine simulation_init_config

  ! ---------------------------------------------------------
  subroutine simulation_start(this)
    type(simulation_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(simulation_start)
    ASSERT(this%init)
    ASSERT(.not.this%start)
    this%start=.true.
    if(this%alloc)then
      ASSERT(.not.associated(this%gr))
      SAFE_ALLOCATE(this%gr)
      call json_get(this%config, "grid", cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call grid_init(this%gr, this%geo, cnfg)
      nullify(cnfg)
    end if
    call domain_init(this%domain, this%gr%sb, this%geo)
    POP_SUB(simulation_start)
    return
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation_extend(this, that, basis)
    type(simulation_t),           intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: that
    type(basis_t),      optional, intent(in)    :: basis
    !
    PUSH_SUB(simulation_extend)
    ASSERT(this%init)
    ASSERT(.not.this%start)
    ASSERT(this%alloc)
    !call simul_box_extend(this%gr%sb, that, basis)
    POP_SUB(simulation_extend)
    return
  end subroutine simulation_extend

  ! ---------------------------------------------------------
  elemental function simulation_get_ndim(this) result(that)
    type(simulation_t),   target, intent(in) :: this
    !
    integer :: that
    !
    that=this%space%dim
    return
  end function simulation_get_ndim

  ! ---------------------------------------------------------
  subroutine simulation_get_config(this, that)
    type(simulation_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(simulation_get_config)
    that=>null()
    if(associated(this%config))&
      that=>this%config
    POP_SUB(simulation_get_config)
    return
  end subroutine simulation_get_config

  ! ---------------------------------------------------------
  subroutine simulation_get_simul_box(this, that)
    type(simulation_t), target, intent(in) :: this
    type(simul_box_t), pointer             :: that
    !
    PUSH_SUB(simulation_get_simul_box)
    that=>null()
    if(associated(this%gr))&
      that=>this%gr%sb
    POP_SUB(simulation_get_simul_box)
    return
  end subroutine simulation_get_simul_box

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that
    !
    PUSH_SUB(simulation_get_grid)
    that=>null()
    if(associated(this%gr))&
      that=>this%gr
    POP_SUB(simulation_get_grid)
    return
  end subroutine simulation_get_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that)
    type(simulation_t), target, intent(in) :: this
    type(mesh_t),      pointer             :: that
    !
    PUSH_SUB(simulation_get_mesh)
    that=>null()
    if(associated(this%gr))&
      that=>this%gr%mesh
    POP_SUB(simulation_get_mesh)
    return
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_geometry(this, that)
    type(simulation_t), target, intent(in) :: this
    type(geometry_t),  pointer             :: that
    !
    PUSH_SUB(simulation_get_geometry)
    that=>null()
    if(associated(this%geo))&
      that=>this%geo
    POP_SUB(simulation_get_geometry)
    return
  end subroutine simulation_get_geometry

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
  subroutine simulation_copy(this, that)
    type(simulation_t), intent(out) :: this
    type(simulation_t), intent(in)  :: that
    !
    PUSH_SUB(simulation_copy)
    this%config=>that%config
    this%geo=>that%geo
    this%init=that%init
    this%alloc=that%alloc
    if(this%alloc)then
      this%gr=>null()
      if(associated(this%gr))then
        SAFE_ALLOCATE(this%gr)
        call grid_copy(this%gr, that%gr)
      end if
    else
      this%gr=>that%gr
    end if
    call domain_copy(this%domain, that%domain)
    POP_SUB(simulation_copy)
    return
  end subroutine simulation_copy

  ! ---------------------------------------------------------
  subroutine simulation_end(this)
    type(simulation_t), intent(inout) :: this
    !
    PUSH_SUB(simulation_end)
    call domain_end(this%domain)
    if(this%init)then
      if(this%alloc)then
        if(associated(this%gr))then
          call grid_end(this%gr)
          SAFE_DEALLOCATE_P(this%gr)
        end if
      end if
      nullify(this%gr)
      this%alloc=.false.
      this%init=.false.
    end if
    nullify(this%config, this%geo)
    POP_SUB(simulation_end)
    return
  end subroutine simulation_end

end module simulation_m

!! Local Variables:
!! mode: f90
!! End:
