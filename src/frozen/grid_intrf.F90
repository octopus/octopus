#include "global.h"

module grid_intrf_oct_m

  use derivatives_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
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
    grid_intrf_alloc,  &
    grid_intrf_init,   &
    grid_intrf_set,    &
    grid_intrf_get,    &
    grid_intrf_copy,   & 
    grid_intrf_end

  integer, parameter :: GRID_INTRF_DISA = 0
  integer, parameter :: GRID_INTRF_NULL = 1
  integer, parameter :: GRID_INTRF_ASSC = 2
  integer, parameter :: GRID_INTRF_ALLC = 3

  type :: grid_intrf_t
    private
    type(json_object_t), pointer :: config =>null()
    type(space_t),       pointer :: space  =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(grid_t),        pointer :: self   =>null()
    integer                      :: type   = GRID_INTRF_DISA
  end type grid_intrf_t

  interface grid_intrf_new
    module procedure grid_intrf_new_pass
    module procedure grid_intrf_new_copy
  end interface grid_intrf_new

  interface grid_intrf_del
    module procedure grid_intrf_del_type
    module procedure grid_intrf_del_pass
  end interface grid_intrf_del

  interface grid_intrf_init
    module procedure grid_intrf_init_cnfg
    module procedure grid_intrf_init_type
    module procedure grid_intrf_init_copy
  end interface grid_intrf_init

  interface grid_intrf_set
    module procedure grid_intrf_set_type
  end interface grid_intrf_set

  interface grid_intrf_get
    module procedure grid_intrf_get_config
    module procedure grid_intrf_get_type
    module procedure grid_intrf_get_space
    module procedure grid_intrf_get_geometry
    module procedure grid_intrf_get_simul_box
    module procedure grid_intrf_get_mesh
    module procedure grid_intrf_get_group
    module procedure grid_intrf_get_derivatives
  end interface grid_intrf_get

  interface grid_intrf_copy
    module procedure grid_intrf_copy_type
    module procedure grid_intrf_copy_pass
  end interface grid_intrf_copy

  interface grid_intrf_end
    module procedure grid_intrf_end_type
    module procedure grid_intrf_end_pass
  end interface grid_intrf_end

contains

  ! ---------------------------------------------------------
  subroutine grid_intrf__new__(this, that)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_t),      pointer        :: that

    PUSH_SUB(grid_intrf__new__)
    
    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_isnull(this))
    nullify(that)
    SAFE_ALLOCATE(that)
    call grid_intrf_set(this, that)
    this%type = GRID_INTRF_ALLC

    POP_SUB(grid_intrf__new__)
  end subroutine grid_intrf__new__

  ! ---------------------------------------------------------
  subroutine grid_intrf__del__(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf__del__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_assoc(this))
    if(grid_intrf_alloc(this))then
      SAFE_DEALLOCATE_P(this%self)
    end if
    nullify(this%self)
    this%type = GRID_INTRF_NULL

    POP_SUB(grid_intrf__del__)
  end subroutine grid_intrf__del__

  ! ---------------------------------------------------------
  subroutine grid_intrf_new_pass(this, init)
    type(grid_intrf_t), intent(inout) :: this

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

    type(grid_t), pointer :: self

    PUSH_SUB(grid_intrf_new_pass)
    
    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_isnull(this))
    call grid_intrf__new__(this, self)
    call init(self, this%geo, this%space, this%config)

    POP_SUB(grid_intrf_new_pass)
  end subroutine grid_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine grid_intrf_new_copy(this, that)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_intrf_t), intent(in)    :: that

    type(grid_t), pointer :: self

    PUSH_SUB(grid_intrf_new_copy)

    ASSERT(associated(this%config))
    ASSERT(grid_intrf_assoc(that))
    call grid_intrf__new__(this, self)
    ASSERT(.false.)
    !call grid_copy(self, that%self)
    nullify(self)

    POP_SUB(grid_intrf_new_copy)
  end subroutine grid_intrf_new_copy

  ! ---------------------------------------------------------
  subroutine grid_intrf_del_type(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf_del_type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_assoc(this))
    call grid_intrf_del(this, grid_end)
    
    POP_SUB(grid_intrf_del_type)
  end subroutine grid_intrf_del_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_del_pass(this, finis)
    type(grid_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(grid_intrf_del_pass)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_assoc(this))
    if(grid_intrf_alloc(this)) call finis(this%self)
    call grid_intrf__del__(this)
    
    POP_SUB(grid_intrf_del_pass)
  end subroutine grid_intrf_del_pass

  ! ---------------------------------------------------------
  function grid_intrf_isnull(this) result(that)
    type(grid_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(grid_intrf_isnull)

    select case(this%type)
    case(GRID_INTRF_DISA)
      that = .false.
    case(GRID_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .true.
    case(GRID_INTRF_ASSC, GRID_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .false.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(grid_intrf_isnull)
  end function grid_intrf_isnull

  ! ---------------------------------------------------------
  function grid_intrf_assoc(this) result(that)
    type(grid_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(grid_intrf_assoc)

    select case(this%type)
    case(GRID_INTRF_DISA)
      that = .false.
    case(GRID_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(GRID_INTRF_ASSC, GRID_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(grid_intrf_assoc)
  end function grid_intrf_assoc

  ! ---------------------------------------------------------
  function grid_intrf_alloc(this) result(that)
    type(grid_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(grid_intrf_alloc)

    select case(this%type)
    case(GRID_INTRF_DISA)
      that = .false.
    case(GRID_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(GRID_INTRF_ASSC)
      ASSERT(associated(this%self))
      that = .false.
    case(GRID_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(grid_intrf_alloc)
  end function grid_intrf_alloc

  ! ---------------------------------------------------------
  subroutine grid_intrf_init_cnfg(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(grid_intrf_init_cnfg)

    call json_init(this)

    POP_SUB(grid_intrf_init_cnfg)
  end subroutine grid_intrf_init_cnfg

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
    nullify(this%self)
    this%type = GRID_INTRF_NULL

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
  subroutine grid_intrf_set_type(this, that)
    type(grid_intrf_t),   intent(inout) :: this
    type(grid_t), target, intent(in)    :: that

    PUSH_SUB(grid_intrf_set_type)

    ASSERT(associated(this%config))
    ASSERT(associated(this%space))
    ASSERT(associated(this%geo))
    ASSERT(grid_intrf_isnull(this))
    this%self => that
    this%type = GRID_INTRF_ASSC

    POP_SUB(grid_intrf_set_type)
  end subroutine grid_intrf_set_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_config(this, that)
    type(grid_intrf_t),   intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(grid_intrf_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(grid_intrf_get_config)
  end subroutine grid_intrf_get_config

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_space(this, that)
    type(grid_intrf_t), intent(in) :: this
    type(space_t),     pointer     :: that

    PUSH_SUB(grid_intrf_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(grid_intrf_get_space)
  end subroutine grid_intrf_get_space

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_geometry(this, that)
    type(grid_intrf_t), intent(in) :: this
    type(geometry_t),  pointer     :: that

    PUSH_SUB(grid_intrf_get_geometry)

    nullify(that)
    if(associated(this%geo)) that => this%geo

    POP_SUB(grid_intrf_get_geometry)
  end subroutine grid_intrf_get_geometry

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_simul_box(this, that)
    type(grid_intrf_t), intent(in) :: this
    type(simul_box_t), pointer     :: that

    PUSH_SUB(grid_intrf_get_simul_box)

    nullify(that)
    if(grid_intrf_assoc(this)) that => this%self%sb

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
    if(grid_intrf_assoc(this))then
      if(fn)then
        that => this%self%fine%mesh
      else
        that => this%self%mesh
      end if
    end if

    POP_SUB(grid_intrf_get_mesh)
  end subroutine grid_intrf_get_mesh

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_group(this, that)
    type(grid_intrf_t), intent(in) :: this
    type(mpi_grp_t),   pointer     :: that

    PUSH_SUB(grid_intrf_get_group)

    nullify(that)
    if(grid_intrf_assoc(this)) that => this%self%mesh%mpi_grp

    POP_SUB(grid_intrf_get_group)
  end subroutine grid_intrf_get_group

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
    if(grid_intrf_assoc(this))then
      if(fn)then
        that => this%self%fine%der
      else
        that => this%self%der
      end if
    end if

    POP_SUB(grid_intrf_get_derivatives)
  end subroutine grid_intrf_get_derivatives

  ! ---------------------------------------------------------
  subroutine grid_intrf_get_type(this, that)
    type(grid_intrf_t), intent(in) :: this
    type(grid_t),      pointer     :: that

    PUSH_SUB(grid_intrf_get_type)

    ASSERT(this%type>GRID_INTRF_DISA)
    nullify(that)
    if(grid_intrf_assoc(this)) that => this%self

    POP_SUB(grid_intrf_get_type)
  end subroutine grid_intrf_get_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_copy_type(this, that)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_intrf_t), intent(in)    :: that

    PUSH_SUB(grid_intrf_copy_type)

    call grid_intrf_end(this)
    if(associated(that%config))then
      call grid_intrf_init(this, that)
      select case(that%type)
      case(GRID_INTRF_NULL)
      case(GRID_INTRF_ASSC)
        call grid_intrf_set(this, that%self)
      case(GRID_INTRF_ALLC)
        call grid_intrf_new(this, that)
      case default
        ASSERT(.false.)
      end select
      this%type = that%type
    end if

    POP_SUB(grid_intrf_copy_type)
  end subroutine grid_intrf_copy_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_copy_pass(this, that, copy)
    type(grid_intrf_t), intent(inout) :: this
    type(grid_intrf_t), intent(in) :: that

    interface
      subroutine copy(this, that)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
        type(grid_t), intent(in)    :: that
      end subroutine copy
    end interface

    type(grid_t), pointer :: self
    integer               :: type

    PUSH_SUB(grid_intrf_copy_pass)

    type = this%type
    self => this%self
    call grid_intrf__end__(this)
    ASSERT(.not.(associated(that%config).or.(type==grid_INTRF_ALLC)))
    if(associated(that%config))then
      call grid_intrf_init(this, that)
      select case(that%type)
      case(grid_INTRF_NULL)
        ASSERT(type/=grid_INTRF_ALLC)
      case(grid_INTRF_ASSC)
        ASSERT(type/=grid_INTRF_ALLC)
        call grid_intrf_set(this, that%self)
      case(grid_INTRF_ALLC)
        if(type/=grid_INTRF_ALLC) call grid_intrf__new__(this, self)
        call copy(self, that%self)
      case default
        ASSERT(.false.)
      end select
    end if
    nullify(self)

    POP_SUB(grid_intrf_copy_pass)
  end subroutine grid_intrf_copy_pass

  ! ---------------------------------------------------------
  subroutine grid_intrf__end__(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf__end__)

    nullify(this%config, this%space, this%geo, this%self)
    this%type = GRID_INTRF_DISA

    POP_SUB(grid_intrf__end__)
  end subroutine grid_intrf__end__

  ! ---------------------------------------------------------
  subroutine grid_intrf_end_type(this)
    type(grid_intrf_t), intent(inout) :: this

    PUSH_SUB(grid_intrf_end_type)

    call grid_intrf_end(this, grid_end)
    
    POP_SUB(grid_intrf_end_type)
  end subroutine grid_intrf_end_type

  ! ---------------------------------------------------------
  subroutine grid_intrf_end_pass(this, finis)
    type(grid_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use grid_oct_m
        type(grid_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(grid_intrf_end_pass)

    if(grid_intrf_alloc(this)) call grid_intrf_del(this, finis)
    call grid_intrf__end__(this)

    POP_SUB(grid_intrf_end_pass)
  end subroutine grid_intrf_end_pass

end module grid_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:
