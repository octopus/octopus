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
  use map_m,       onlY: map_t, map_init, map_end
  use mesh_m,      only: mesh_t 
  use simul_box_m, only: simul_box_t
  use space_m,     only: space_t

  use wrap_mesh_m, only: wrap_mesh_t, wrap_mesh_init, wrap_mesh_end

  use meshmap_table_m, only:          &
    table_t    => meshmap_table_t,    &
    table_init => meshmap_table_init, &
    table_pop  => meshmap_table_pop,  &
    table_get  => meshmap_table_get,  &
    table_set  => meshmap_table_set,  &
    table_copy => meshmap_table_copy, &
    table_end  => meshmap_table_end

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
  
  integer, parameter :: TABLE_INIT_LEN = 7

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
    type(table_t)                :: table
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
    module procedure simulation_get_map
  end interface simulation_get

contains

  ! ---------------------------------------------------------
  subroutine simulation_init_common(this, geo, space)
    type(simulation_t),       intent(out) :: this
    type(geometry_t), target, intent(in)  :: geo
    type(space_t),    target, intent(in)  :: space
    !
    this%config=>null()
    this%geo=>geo
    this%space=>space
    this%gr=>null()
    this%init=.true.
    this%start=.false.
    this%alloc=.false.
    return
  end subroutine simulation_init_common

  ! ---------------------------------------------------------
  subroutine simulation_init_grid(this, geo, space, grid)
    type(simulation_t),   intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(space_t),        intent(in)  :: space
    type(grid_t), target, intent(in)  :: grid
    !
    call simulation_init_common(this, geo, space)
    this%gr=>grid
    this%alloc=.false.
    return
  end subroutine simulation_init_grid

  ! ---------------------------------------------------------
  subroutine simulation_init_config(this, geo, space, config)
    type(simulation_t),          intent(out) :: this
    type(geometry_t),            intent(in)  :: geo
    type(space_t),               intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config
    !
    print *, "***: simulation_init_config"
    call simulation_init_common(this, geo, space)
    this%config=>config
    this%alloc=.true.
    return
  end subroutine simulation_init_config

  ! ---------------------------------------------------------
  subroutine simulation_start(this)
    type(simulation_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
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
    call table_init(this%table, TABLE_INIT_LEN)
    return
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation_extend(this, that, basis)
    type(simulation_t),           intent(inout) :: this
    type(simulation_t), optional, intent(in)    :: that
    type(basis_t),      optional, intent(in)    :: basis
    !
    print *, "***: simulation_extend"
    ASSERT(this%init)
    ASSERT(.not.this%start)
    ASSERT(this%alloc)
    !call simul_box_extend(this%gr%sb, that, basis)
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
    that=>null()
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine simulation_get_config

  ! ---------------------------------------------------------
  subroutine simulation_get_simul_box(this, that)
    type(simulation_t), target, intent(in) :: this
    type(simul_box_t), pointer             :: that
    !
    that=>null()
    if(associated(this%gr))&
      that=>this%gr%sb
    return
  end subroutine simulation_get_simul_box

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that
    !
    that=>null()
    if(associated(this%gr))&
      that=>this%gr
    return
  end subroutine simulation_get_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that)
    type(simulation_t), target, intent(in) :: this
    type(mesh_t),      pointer             :: that
    !
    that=>null()
    if(associated(this%gr))&
      that=>this%gr%mesh
    return
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_geometry(this, that)
    type(simulation_t), target, intent(in) :: this
    type(geometry_t),  pointer             :: that
    !
    that=>null()
    if(associated(this%geo))&
      that=>this%geo
    return
  end subroutine simulation_get_geometry

  ! ---------------------------------------------------------
  subroutine simulation_get_domain(this, that)
    type(simulation_t), target, intent(in) :: this
    type(domain_t),    pointer             :: that
    !
    that=>this%domain
    return
  end subroutine simulation_get_domain

  ! ---------------------------------------------------------
  subroutine simulation_get_map(this, mesh, that)
    type(simulation_t), target, intent(inout) :: this
    type(mesh_t),               intent(in)    :: mesh
    type(map_t),       pointer                :: that
    !
    type(wrap_mesh_t), pointer :: wmsh
    type(map_t),       pointer :: map
    !
    nullify(that, wmsh, map)
    SAFE_ALLOCATE(wmsh)
    call wrap_mesh_init(wmsh, mesh)
    call table_get(this%table, wmsh, that)
    if(associated(that))then
      call wrap_mesh_end(wmsh)
      SAFE_DEALLOCATE_P(wmsh)
    else
      SAFE_ALLOCATE(map)
      call map_init(map, this%gr%mesh, mesh)
      call table_set(this%table, wmsh, map)
      that=>map
    end if
    nullify(wmsh, map)
    return
  end subroutine simulation_get_map

  ! ---------------------------------------------------------
  subroutine simulation_copy(this, that)
    type(simulation_t), intent(out) :: this
    type(simulation_t), intent(in)  :: that
    !
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
    call table_copy(this%table, that%table)
    return
  end subroutine simulation_copy

  ! ---------------------------------------------------------
  subroutine simulation_end(this)
    type(simulation_t), intent(inout) :: this
    !
    type(wrap_mesh_t), pointer :: wmsh
    type(map_t),       pointer :: map
    !
    nullify(wmsh, map)
    do
      call table_pop(this%table, wmsh, map)
      if((.not.associated(wmsh)).or.(.not.associated(map)))exit
      call wrap_mesh_end(wmsh)
      SAFE_DEALLOCATE_P(wmsh)
      call map_end(map)
      SAFE_DEALLOCATE_P(map)
      nullify(wmsh, map)
    end do
    call table_end(this%table)
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
    return
  end subroutine simulation_end

end module simulation_m

!! Local Variables:
!! mode: f90
!! End:
