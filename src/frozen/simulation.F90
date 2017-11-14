#include "global.h"

module simulation_oct_m

  use derivatives_oct_m
  use domain_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use grid_intrf_oct_m
  use json_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use space_oct_m

  implicit none

  private

  public ::       &
    simulation_t

  public ::            &
    simulation_assoc,  &
    simulation_init,   &
    simulation_start,  &
    simulation_stop,   &
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

  interface simulation_init
    module procedure simulation_init_cnfg
    module procedure simulation_init_type
  end interface simulation_init

  interface simulation_start
    module procedure simulation_start_type
    module procedure simulation_start_pass
  end interface simulation_start

  interface simulation_stop
    module procedure simulation_stop_type
    module procedure simulation_stop_pass
  end interface simulation_stop

  interface simulation_get
    module procedure simulation_get_config
    module procedure simulation_get_info
    module procedure simulation_get_space
    module procedure simulation_get_mesh
    module procedure simulation_get_group
    module procedure simulation_get_derivatives
    module procedure simulation_get_grid
    module procedure simulation_get_domain
  end interface simulation_get

  interface simulation_copy
    module procedure simulation_copy_type
    module procedure simulation_copy_pass
  end interface simulation_copy

  interface simulation_end
    module procedure simulation_end_type
    module procedure simulation_end_pass
  end interface simulation_end

contains

  ! ---------------------------------------------------------
  function simulation_assoc(this) result(that)
    type(simulation_t), intent(in) :: this

    logical :: that

    PUSH_SUB(simulation_assoc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    that = grid_intrf_assoc(this%igrid)

    POP_SUB(simulation_assoc)
  end function simulation_assoc

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
  subroutine simulation_init_cnfg(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(simulation_init_cnfg)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call grid_intrf_init(cnfg)
    call json_set(this, "grid", cnfg)
    nullify(cnfg)

    POP_SUB(simulation_init_cnfg)
  end subroutine simulation_init_cnfg

  ! ---------------------------------------------------------
  subroutine simulation_init_type(this, geo, space, config)
    type(simulation_t),  intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(simulation_init_type)

    nullify(cnfg)
    call simulation__init__(this, space, config)
    call json_get(config, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call grid_intrf_init(this%igrid, geo, space, cnfg)
    nullify(cnfg)
    call domain_init(this%domain, space)

    POP_SUB(simulation_init_type)
  end subroutine simulation_init_type

  ! ---------------------------------------------------------
  subroutine simulation_start_type(this, grid)
    type(simulation_t), intent(inout) :: this
    type(grid_t),       intent(in)    :: grid

    PUSH_SUB(simulation_start_type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.grid_intrf_assoc(this%igrid))
    call grid_intrf_set(this%igrid, grid)
    ASSERT(grid_intrf_assoc(this%igrid))
    call domain_start(this%domain, this%igrid)

    POP_SUB(simulation_start_type)
  end subroutine simulation_start_type

  ! ---------------------------------------------------------
  subroutine simulation_start_pass(this, init)
    type(simulation_t), intent(inout) :: this

    interface
      subroutine init(this, geo, space, config)
        use geometry_oct_m
        use grid_oct_m
        use json_oct_m
        use space_oct_m
        type(grid_t),        intent(out) :: this
        type(geometry_t),    intent(in)  :: geo
        type(space_t),       intent(in)  :: space
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    PUSH_SUB(simulation_start_pass)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(.not.grid_intrf_assoc(this%igrid))
    call grid_intrf_new(this%igrid, init)
    ASSERT(grid_intrf_assoc(this%igrid))
    call domain_start(this%domain, this%igrid)

    POP_SUB(simulation_start_pass)
  end subroutine simulation_start_pass

  ! ---------------------------------------------------------
  subroutine simulation_stop_type(this)
    type(simulation_t), intent(inout) :: this

    PUSH_SUB(simulation_stop_type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(grid_intrf_assoc(this%igrid))
    call grid_intrf_del(this%igrid)
    ASSERT(.not.grid_intrf_assoc(this%igrid))
    call domain_stop(this%domain)

    POP_SUB(simulation_stop_type)
  end subroutine simulation_stop_type

  ! ---------------------------------------------------------
  subroutine simulation_stop_pass(this, finis)
    type(simulation_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(simulation_stop_pass)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(grid_intrf_assoc(this%igrid))
    call grid_intrf_del(this%igrid, finis)
    ASSERT(.not.grid_intrf_assoc(this%igrid))
    call domain_stop(this%domain)

    POP_SUB(simulation_stop_pass)
  end subroutine simulation_stop_pass

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
  subroutine simulation_get_group(this, that)
    type(simulation_t), intent(in) :: this
    type(mpi_grp_t),   pointer     :: that

    PUSH_SUB(simulation_get_group)

    call grid_intrf_get(this%igrid, that)

    POP_SUB(simulation_get_group)
  end subroutine simulation_get_group

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
  subroutine simulation_get_domain(this, that)
    type(simulation_t), target, intent(in) :: this
    type(domain_t),    pointer             :: that

    PUSH_SUB(simulation_get_domain)

    that => this%domain

    POP_SUB(simulation_get_domain)
  end subroutine simulation_get_domain

  ! ---------------------------------------------------------
  subroutine simulation_copy_type(this, that)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that

    PUSH_SUB(simulation_copy_type)

    call simulation_end(this)
    if(associated(that%config).and.associated(that%space))then
      call simulation__init__(this, that%space, that%config)
      call grid_intrf_copy(this%igrid, that%igrid)
      call domain_copy(this%domain, that%domain)
    end if

    POP_SUB(simulation_copy_type)
  end subroutine simulation_copy_type

  ! ---------------------------------------------------------
  subroutine simulation_copy_pass(this, that, copy)
    type(simulation_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: that

    interface
      subroutine copy(this, that)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
        type(grid_t), intent(in)    :: that
      end subroutine copy
    end interface

    PUSH_SUB(simulation_copy_pass)

    call simulation_end(this)
    if(associated(that%config).and.associated(that%space))then
      call simulation__init__(this, that%space, that%config)
      call grid_intrf_copy(this%igrid, that%igrid, copy)
      call domain_copy(this%domain, that%domain)
    end if

    POP_SUB(simulation_copy_pass)
  end subroutine simulation_copy_pass

  ! ---------------------------------------------------------
  subroutine simulation_end_type(this)
    type(simulation_t), intent(inout) :: this

    PUSH_SUB(simulation_end_type)

    nullify(this%config, this%space)
    call grid_intrf_end(this%igrid)
    call domain_end(this%domain)

    POP_SUB(simulation_end_type)
  end subroutine simulation_end_type

  ! ---------------------------------------------------------
  subroutine simulation_end_pass(this, finis)
    type(simulation_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(simulation_end_pass)

    nullify(this%config, this%space)
    call grid_intrf_end(this%igrid, finis)
    call domain_end(this%domain)

    POP_SUB(simulation_end_pass)
  end subroutine simulation_end_pass

end module simulation_oct_m

!! Local Variables:
!! mode: f90
!! End:
