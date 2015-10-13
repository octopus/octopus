#include "global.h"

module simulation_m

  use derivatives_m
  use domain_m
  use geometry_m
  use global_m
  use grid_m
  use json_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m
  use space_m

  implicit none

  private

  public ::            &
    simulation_init,   &
    simulation_start,  &
    simulation_extend, &
    simulation_get,    &
    simulation_copy,   & 
    simulation_end

  type, public :: simulation_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(grid_t),        pointer :: grid   =>null()
    type(domain_t)               :: domain
  end type simulation_t

  interface simulation_init
    module procedure simulation_init_simulation
    module procedure simulation_init_copy
  end interface simulation_init

  interface simulation_get
    module procedure simulation_get_config
    module procedure simulation_get_info
    module procedure simulation_get_space
    module procedure simulation_get_simul_box
    module procedure simulation_get_mesh
    module procedure simulation_get_derivatives
    module procedure simulation_get_grid
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
  subroutine simulation_init_simulation(this, space, config)
    type(simulation_t),  intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(simulation_init_simulation)

    call simulation__init__(this, space, config)
    call domain_init(this%domain, space)

    POP_SUB(simulation_init_simulation)
  end subroutine simulation_init_simulation

  ! ---------------------------------------------------------
  subroutine simulation_init_copy(this, that)
    type(simulation_t), intent(out) :: this
    type(simulation_t), intent(in)  :: that

    PUSH_SUB(simulation_init_copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    call simulation_init_simulation(this, that%space, that%config)

    POP_SUB(simulation_init_copy)
  end subroutine simulation_init_copy

  ! ---------------------------------------------------------
  subroutine simulation__start__(this, grid)
    type(simulation_t),   intent(inout) :: this
    type(grid_t), target, intent(in)    :: grid

    PUSH_SUB(simulation__start__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.associated(this%grid))
    this%grid => grid
    ASSERT(this%grid%sb%dim==this%space%dim)

    POP_SUB(simulation__start__)
  end subroutine simulation__start__


  ! ---------------------------------------------------------
  subroutine simulation_start(this, grid, geo)
    type(simulation_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid
    type(geometry_t),   intent(in)    :: geo

    PUSH_SUB(simulation_start)

    call simulation__start__(this, grid)
    call domain_start(this%domain, this%grid%sb, geo)

    POP_SUB(simulation_start)
  end subroutine simulation_start

  ! ---------------------------------------------------------
  subroutine simulation_extend(this, that, config)
    type(simulation_t),  intent(inout) :: this
    type(domain_t),      intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(simulation_extend)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.associated(this%grid))
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
        call domain_extend(this%domain, that, cnfg)
      end do
      call json_end(iter)
      nullify(cnfg, list)
    else
      call domain_extend(this%domain, that)
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
  subroutine simulation_get_simul_box(this, that)
    type(simulation_t), intent(in) :: this
    type(simul_box_t), pointer     :: that

    PUSH_SUB(simulation_get_simul_box)

    nullify(that)
    if(associated(this%grid)) that => this%grid%sb

    POP_SUB(simulation_get_simul_box)
  end subroutine simulation_get_simul_box

  ! ---------------------------------------------------------
  subroutine simulation_get_mesh(this, that, fine)
    type(simulation_t), intent(in) :: this
    type(mesh_t),      pointer     :: that
    logical,  optional, intent(in) :: fine

    logical :: fn

    PUSH_SUB(simulation_get_mesh)

    fn = .false.
    nullify(that)
    if(present(fine)) fn = fine
    if(associated(this%grid))then
      if(fn)then
        that => this%grid%fine%mesh
      else
        that => this%grid%mesh
      end if
    end if

    POP_SUB(simulation_get_mesh)
  end subroutine simulation_get_mesh

  ! ---------------------------------------------------------
  subroutine simulation_get_derivatives(this, that, fine)
    type(simulation_t),   intent(in) :: this
    type(derivatives_t), pointer     :: that
    logical,    optional, intent(in) :: fine

    logical :: fn

    PUSH_SUB(simulation_get_derivatives)

    fn = .false.
    nullify(that)
    if(present(fine)) fn = fine
    if(associated(this%grid))then
      if(fn)then
        that => this%grid%fine%der
      else
        that => this%grid%der
      end if
    end if

    POP_SUB(simulation_get_derivatives)
  end subroutine simulation_get_derivatives

  ! ---------------------------------------------------------
  subroutine simulation_get_grid(this, that)
    type(simulation_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that

    PUSH_SUB(simulation_get_grid)

    nullify(that)
    if(associated(this%grid)) that => this%grid

    POP_SUB(simulation_get_grid)
  end subroutine simulation_get_grid

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
    if(associated(that%config).and.associated(that%space))then
      call simulation__init__(this, that%space, that%config)
      if(associated(that%grid)) call simulation__start__(this, that%grid)
    end if
    call domain_copy(this%domain, that%domain)

    POP_SUB(simulation_copy)
  end subroutine simulation_copy

  ! ---------------------------------------------------------
  subroutine simulation_end(this)
    type(simulation_t), intent(inout) :: this

    PUSH_SUB(simulation_end)

    nullify(this%config, this%space, this%grid)
    call domain_end(this%domain)

    POP_SUB(simulation_end)
  end subroutine simulation_end

end module simulation_m

!! Local Variables:
!! mode: f90
!! End:
