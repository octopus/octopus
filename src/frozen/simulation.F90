#include "global.h"

module simulation_m

  use derivatives_m
  use domain_m
  use geometry_m
  use global_m
  use grid_m
  use grid_intrf_m
  use json_m
  use mesh_m
  use messages_m
  use profiling_m
  use space_m

  implicit none

  private

  public ::       &
    simulation_t

  public ::            &
    simulation_init,   &
    simulation_start,  &
    simulation_extend, &
    simulation_get,    &
    simulation_copy,   & 
    simulation_end

  type :: simulation_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(grid_intrf_t)           :: igrid
    type(domain_t)               :: domain
  end type simulation_t

  interface simulation_get
    module procedure simulation_get_config
    module procedure simulation_get_info
    module procedure simulation_get_space
    module procedure simulation_get_mesh
    module procedure simulation_get_derivatives
    module procedure simulation_get_grid
    module procedure simulation_get_grid_intrf
    module procedure simulation_get_domain
  end interface simulation_get

contains

  ! ---------------------------------------------------------
  subroutine simulation__init__(this, space, config)
    type(simulation_t),          intent(out) :: this
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(simulation__init__)

    this%config => config
    this%space => space

    POP_SUB(simulation__init__)
  end subroutine simulation__init__

  ! ---------------------------------------------------------
  subroutine simulation_init(this, geo, space, config)
    type(simulation_t),  intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(simulation_init)

    nullify(cnfg)
    call simulation__init__(this, space, config)
    call json_get(config, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call grid_intrf_init(this%igrid, geo, space, cnfg)
    nullify(cnfg)
    call domain_init(this%domain, space)

    POP_SUB(simulation_init)
  end subroutine simulation_init

  ! ---------------------------------------------------------
  subroutine simulation__start__(this, grid)
    type(simulation_t),     intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(simulation__start__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    if(present(grid)) call grid_intrf_set(this%igrid, grid)

    POP_SUB(simulation__start__)
  end subroutine simulation__start__


  ! ---------------------------------------------------------
  subroutine simulation_start(this, grid)
    type(simulation_t),     intent(inout) :: this
    type(grid_t), optional, intent(in)    :: grid

    PUSH_SUB(simulation_start)

    call simulation__start__(this, grid)
    call grid_intrf_start(this%igrid)
    call domain_start(this%domain, this%igrid)

    POP_SUB(simulation_start)
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation_extend(this, that, config)
    type(simulation_t),  intent(inout) :: this
    type(simulation_t),  intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(simulation_extend)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    nullify(cnfg, list)
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(associated(list))
    if(json_len(list)>0)then
      call json_init(iter, list)
      do
        nullify(cnfg)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call domain_extend(this%domain, that%domain, cnfg)
      end do
      call json_end(iter)
      nullify(cnfg, list)
    else
      call domain_extend(this%domain, that%domain)
    end if

    POP_SUB(simulation_extend)
  end subroutine simulation_extend

  ! ---------------------------------------------------------
  subroutine simulation_get_info(this, ndim)
    type(simulation_t), intent(in)  :: this
    integer,            intent(out) :: ndim

    PUSH_SUB(simulation_get_info)

    ndim = 0
    if(associated(this%space)) ndim = this%space%dim

    POP_SUB(simulation_get_info)
  end subroutine simulation_get_info

  ! ---------------------------------------------------------
  subroutine simulation_get_config(this, that)
    type(simulation_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(simulation_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(simulation_get_config)
  end subroutine simulation_get_config

  ! ---------------------------------------------------------
  subroutine simulation_get_space(this, that)
    type(simulation_t), target, intent(in) :: this
    type(space_t),     pointer             :: that

    PUSH_SUB(simulation_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(simulation_get_space)
  end subroutine simulation_get_space

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that, fine)
    type(simulation_t), intent(in) :: this
    type(mesh_t),      pointer     :: that
    logical,  optional, intent(in) :: fine

    PUSH_SUB(simulation_get_mesh)

    call grid_intrf_get(this%igrid, that, fine)

    POP_SUB(simulation_get_mesh)
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_derivatives(this, that, fine)
    type(simulation_t),   intent(in) :: this
    type(derivatives_t), pointer     :: that
    logical,    optional, intent(in) :: fine

    PUSH_SUB(simulation_get_derivatives)

    call grid_intrf_get(this%igrid, that, fine)

    POP_SUB(simulation_get_derivatives)
  end subroutine simulation_get_derivatives

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), intent(in) :: this
    type(grid_t),      pointer     :: that

    PUSH_SUB(simulation_get_grid)

    call grid_intrf_get(this%igrid, that)

    POP_SUB(simulation_get_grid)
  end subroutine simulation_get_grid

  ! ---------------------------------------------------------
  subroutine simulation_get_grid_intrf(this, that)
    type(simulation_t),  target, intent(in) :: this
    type(grid_intrf_t), pointer             :: that

    PUSH_SUB(simulation_get_grid_intrf)

    that => this%igrid

    POP_SUB(simulation_get_grid_intrf)
  end subroutine simulation_get_grid_intrf

  ! ---------------------------------------------------------
  subroutine simulation_get_domain(this, that)
    type(simulation_t), target, intent(in) :: this
    type(domain_t),    pointer             :: that

    PUSH_SUB(simulation_get_domain)

    that => this%domain

    POP_SUB(simulation_get_domain)
  end subroutine simulation_get_domain

  ! ---------------------------------------------------------
  subroutine simulation_copy(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that

    PUSH_SUB(simulation_copy)

    call simulation_end(this)
    if(associated(that%config).and.associated(that%space))&
      call simulation__init__(this, that%space, that%config)
    call grid_intrf_copy(this%igrid, that%igrid)
    call domain_copy(this%domain, that%domain)

    POP_SUB(simulation_copy)
  end subroutine simulation_copy

  ! ---------------------------------------------------------
  subroutine simulation_end(this)
    type(simulation_t), intent(inout) :: this

    PUSH_SUB(simulation_end)

    nullify(this%config, this%space)
    call grid_intrf_end(this%igrid)
    call domain_end(this%domain)

    POP_SUB(simulation_end)
  end subroutine simulation_end

end module simulation_m

!! Local Variables:
!! mode: f90
!! End:
