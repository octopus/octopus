#include "global.h"

module grid_intrf_oct_m

  use derivatives_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m

  implicit none

  private

  public ::       &
    grid_intrf_t

  public ::            &
    grid_intrf_new,    &
    grid_intrf_del,    &
    grid_intrf_assoc,  &
    grid_intrf_init,   &
    grid_intrf_set,    &
    grid_intrf_get,    &
    grid_intrf_copy,   & 
    grid_intrf_end

  integer, parameter :: GRID_NULL = 0
  integer, parameter :: GRID_ASSC = 1
  integer, parameter :: GRID_ALLC = 2

  type :: grid_intrf_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(grid_t),        pointer :: grid   =>null()
    integer                      :: type   = GRID_NULL
  end type grid_intrf_t

  interface grid_intrf_new
    module procedure grid_intrf_new_grid
    module procedure grid_intrf_new_pass
  end interface grid_intrf_new

  interface grid_intrf_del
    module procedure grid_intrf_del_grid
    module procedure grid_intrf_del_pass
  end interface grid_intrf_del

  interface grid_intrf_init
    module procedure grid_intrf_init_config
    module procedure grid_intrf_init_type
    module procedure grid_intrf_init_copy
  end interface grid_intrf_init

  interface grid_intrf_set
    module procedure grid_intrf_set_grid
  end interface grid_intrf_set

  interface grid_intrf_get
    module procedure grid_intrf_get_config
    module procedure grid_intrf_get_space
    module procedure grid_intrf_get_geometry
    module procedure grid_intrf_get_simul_box
    module procedure grid_intrf_get_mesh
    module procedure grid_intrf_get_derivatives
    module procedure grid_intrf_get_grid
  end interface grid_intrf_get

contains

  ! ---------------------------------------------------------
  subroutine grid_intrf__new__(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf__new__)
    
    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(.not.associated(this%grid))
    nullify(this%grid)
    SAFE_ALLOCATE(this%grid)
    !call grid_nullify(this%grid)
    this%type = GRID_ALLC

    POP_SUB(grid_intrf__new__)
  end subroutine grid_intrf__new__

  ! ---------------------------------------------------------
  subroutine grid_intrf_new_grid(this, that)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_t),      pointer        :: that

    PUSH_SUB(grid_intrf_new_grid)
    
    call grid_intrf__new__(this)
    call grid_intrf_get(this, that)

    POP_SUB(grid_intrf_new_grid)
  end subroutine grid_intrf_new_grid

  ! ---------------------------------------------------------
  subroutine grid_intrf_new_pass(this, that, grid_type_init)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_t),      pointer        :: that

    interface
      subroutine grid_type_init(this, geo, space, config)
        use geometry_oct_m
        use grid_oct_m
        use json_oct_m
        use space_oct_m
        type(grid_t),        intent(out) :: this
        type(geometry_t),    intent(in)  :: geo
        type(space_t),       intent(in)  :: space
        type(json_object_t), intent(in)  :: config
      end subroutine grid_type_init
    end interface

    PUSH_SUB(grid_intrf_new_pass)
    
    nullify(that)
    call grid_intrf_new(this, that)
    call grid_type_init(that, this%geo, this%space, this%config)

    POP_SUB(grid_intrf_new_pass)
  end subroutine grid_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine grid_intrf_del_grid(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf_del_grid)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(associated(this%grid))
    ASSERT(this%type==GRID_ALLC)
    !call grid_end(this%grid)
    SAFE_DEALLOCATE_P(this%grid)
    nullify(this%grid)
    this%type = GRID_NULL
    
    POP_SUB(grid_intrf_del_grid)
  end subroutine grid_intrf_del_grid

  ! ---------------------------------------------------------
  subroutine grid_intrf_del_pass(this, grid_type_end)
    type(grid_intrf_t), intent(inout) :: this

    interface
      subroutine grid_type_end(this)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
      end subroutine grid_type_end
    end interface

    PUSH_SUB(grid_intrf_del_pass)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(associated(this%grid))
    ASSERT(this%type==GRID_ALLC)
    call grid_type_end(this%grid)
    call grid_intrf_del(this)
    
    POP_SUB(grid_intrf_del_pass)
  end subroutine grid_intrf_del_pass

  ! ---------------------------------------------------------
  function grid_intrf_assoc(this) result(that)
    type(grid_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(grid_intrf_assoc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    select case(this%type)
    case(GRID_NULL)
      ASSERT(.not.associated(this%grid))
      that = .false.
    case(GRID_ASSC,GRID_ALLC)
      ASSERT(associated(this%grid))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(grid_intrf_assoc)
  end function grid_intrf_assoc

  ! ---------------------------------------------------------
  subroutine grid_intrf_init_config(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(grid_intrf_init_config)

    call json_init(this)

    POP_SUB(grid_intrf_init_config)
  end subroutine grid_intrf_init_config

  ! ---------------------------------------------------------
  subroutine grid_intrf_init_type(this, geo, space, config)
    type(grid_intrf_t),          intent(out) :: this
    type(geometry_t),    target, intent(in)  :: geo
    type(space_t),       target, intent(in)  :: space
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(grid_intrf_init_type)

    this%config => config
    this%space => space
    this%geo => geo
    nullify(this%grid)
    this%type = GRID_NULL

    POP_SUB(grid_intrf_init_type)
  end subroutine grid_intrf_init_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_init_copy(this, that)
    type(grid_intrf_t), intent(out) :: this
    type(grid_intrf_t), intent(in)  :: that

    PUSH_SUB(grid_intrf_init_copy)
    
    ASSERT(associated(that%config))
    ASSERT(associated(that%space))
    ASSERT(associated(that%geo))
    call grid_intrf_init(this, that%geo, that%space, that%config)

    POP_SUB(grid_intrf_init_copy)
  end subroutine grid_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine grid_intrf_set_grid(this, that)
    type(grid_intrf_t),   intent(inout) :: this
    type(grid_t), target, intent(in)    :: that

    PUSH_SUB(grid_intrf_set_grid)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(this%space%dim==that%sb%dim)
    ASSERT(associated(this%geo))
    ASSERT(.not.associated(this%grid))
    ASSERT(this%type==GRID_NULL)
    this%grid => that
    this%type = GRID_ASSC

    POP_SUB(grid_intrf_set_grid)
  end subroutine grid_intrf_set_grid

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_config(this, that)
    type(grid_intrf_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(grid_intrf_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(grid_intrf_get_config)
  end subroutine grid_intrf_get_config

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_space(this, that)
    type(grid_intrf_t), target, intent(in) :: this
    type(space_t),     pointer             :: that

    PUSH_SUB(grid_intrf_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(grid_intrf_get_space)
  end subroutine grid_intrf_get_space

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_geometry(this, that)
    type(grid_intrf_t), target, intent(in) :: this
    type(geometry_t),  pointer             :: that

    PUSH_SUB(grid_intrf_get_geometry)

    nullify(that)
    if(associated(this%geo)) that => this%geo

    POP_SUB(grid_intrf_get_geometry)
  end subroutine grid_intrf_get_geometry

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_simul_box(this, that)
    type(grid_intrf_t), target, intent(in) :: this
    type(simul_box_t), pointer             :: that

    PUSH_SUB(grid_intrf_get_simul_box)

    nullify(that)
    if(associated(this%grid)) that => this%grid%sb

    POP_SUB(grid_intrf_get_simul_box)
  end subroutine grid_intrf_get_simul_box

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_mesh(this, that, fine)
    type(grid_intrf_t), intent(in) :: this
    type(mesh_t),      pointer     :: that
    logical,  optional, intent(in) :: fine

    logical :: fn

    PUSH_SUB(grid_intrf_get_mesh)

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

    POP_SUB(grid_intrf_get_mesh)
  end subroutine grid_intrf_get_mesh

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_derivatives(this, that, fine)
    type(grid_intrf_t),   intent(in) :: this
    type(derivatives_t), pointer     :: that
    logical,    optional, intent(in) :: fine

    logical :: fn

    PUSH_SUB(grid_intrf_get_derivatives)

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

    POP_SUB(grid_intrf_get_derivatives)
  end subroutine grid_intrf_get_derivatives

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_grid(this, that)
    type(grid_intrf_t), target, intent(in) :: this
    type(grid_t),      pointer             :: that

    PUSH_SUB(grid_intrf_get_grid)

    nullify(that)
    if(associated(this%grid)) that => this%grid

    POP_SUB(grid_intrf_get_grid)
  end subroutine grid_intrf_get_grid

  ! ---------------------------------------------------------
  subroutine grid_intrf_copy(this, that)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_intrf_t), intent(in)    :: that

    PUSH_SUB(grid_intrf_copy)

    call grid_intrf_end(this)
    this%config => that%config
    this%space => that%space
    this%geo => that%geo
    select case(that%type)
    case(GRID_NULL)
      nullify(this%grid)
    case(GRID_ASSC)
      this%grid => that%grid
    case(GRID_ALLC)
      call grid_intrf__new__(this)
      !call grid_copy(this%grid, that%grid)
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(grid_intrf_copy)
  end subroutine grid_intrf_copy

  ! ---------------------------------------------------------
  subroutine grid_intrf_end(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf_end)

    if(this%type==GRID_ALLC) call grid_intrf_del(this)
    nullify(this%config, this%space, this%geo, this%grid)
    this%type = GRID_NULL

    POP_SUB(grid_intrf_end)
  end subroutine grid_intrf_end

end module grid_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:
