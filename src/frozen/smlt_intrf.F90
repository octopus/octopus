#include "global.h"

module smlt_intrf_oct_m

  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m
  use uuid_oct_m

  implicit none

  private
  public ::       &
    smlt_intrf_t

  public ::           &
    smlt_intrf_new,   &
    smlt_intrf_del,   &
    smlt_intrf_assoc, &
    smlt_intrf_init,  &
    smlt_intrf_apply, &
    smlt_intrf_set,   &
    smlt_intrf_get,   &
    smlt_intrf_copy,  &
    smlt_intrf_end

  integer, parameter :: SMLT_INTRF_DISA = 0
  integer, parameter :: SMLT_INTRF_NULL = 1
  integer, parameter :: SMLT_INTRF_ASSC = 2

  type :: smlt_intrf_t
    private
    type(simulation_t), pointer :: self =>null()
    integer                     :: type = SMLT_INTRF_DISA
    type(uuid_t)                :: uuid
  end type smlt_intrf_t

  interface smlt_intrf_new
    module procedure smlt_intrf_new_type
    module procedure smlt_intrf_new_pass
  end interface smlt_intrf_new

  interface smlt_intrf_del
    module procedure smlt_intrf_del_type
    module procedure smlt_intrf_del_pass
  end interface smlt_intrf_del

  interface smlt_intrf_init
    module procedure smlt_intrf_init_type
    module procedure smlt_intrf_init_copy
  end interface smlt_intrf_init

  interface smlt_intrf_end
    module procedure smlt_intrf_end_type
    module procedure smlt_intrf_end_pass
  end interface smlt_intrf_end

contains

  ! ---------------------------------------------------------
  subroutine smlt_intrf_new_type(this, geo, space, config)
    type(smlt_intrf_t),  intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(space_t),       intent(in)    :: space
    type(json_object_t), intent(in)    :: config

    type(simulation_t), pointer :: sim

    PUSH_SUB(smlt_intrf_new_type)

    ASSERT(smlt_intrf_isnull(this))
    nullify(sim)
    sim => simulation_new(geo, space, config)
    call smlt_intrf_set(this, sim)
    nullify(sim)

    POP_SUB(smlt_intrf_new_type)
  end subroutine smlt_intrf_new_type

  ! ---------------------------------------------------------
  subroutine smlt_intrf_new_pass(this, geo, space, config, init)
    type(smlt_intrf_t),  intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(space_t),       intent(in)    :: space
    type(json_object_t), intent(in)    :: config

    interface
      subroutine init(this, geo, space, config)
        use simulation_oct_m
        use geometry_oct_m
        use space_oct_m
        use json_oct_m
        type(simulation_t),  intent(out) :: this
        type(geometry_t),    intent(in)  :: geo
        type(space_t),       intent(in)  :: space
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    type(simulation_t), pointer :: sim

    PUSH_SUB(smlt_intrf_new_pass)

    ASSERT(smlt_intrf_isnull(this))
    nullify(sim)
    sim => simulation_new(geo, space, config, init)
    call smlt_intrf_set(this, sim)
    nullify(sim)

    POP_SUB(smlt_intrf_new_pass)
  end subroutine smlt_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine smlt_intrf_del_type(this)
    type(smlt_intrf_t), intent(inout) :: this

    PUSH_SUB(smlt_intrf_del_type)

    ASSERT(smlt_intrf_assoc(this))
    call smlt_intrf_del(this, finis)

    POP_SUB(smlt_intrf_del_type)

  contains

    subroutine finis(this)
      type(simulation_t), intent(inout) :: this

      PUSH_SUB(smlt_intrf_del_type.finis)

      call simulation_end(this)

      POP_SUB(smlt_intrf_del_type.finis)
    end subroutine finis

  end subroutine smlt_intrf_del_type

  ! ---------------------------------------------------------
  subroutine smlt_intrf_del_pass(this, finis)
    type(smlt_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use simulation_oct_m
        type(simulation_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(smlt_intrf_del_pass)

    ASSERT(smlt_intrf_assoc(this))
    call simulation_detach(this%self, this%uuid)
    call simulation_del(this%self, finis)
    nullify(this%self)
    this%type = SMLT_INTRF_NULL

    POP_SUB(smlt_intrf_del_pass)
  end subroutine smlt_intrf_del_pass

  ! ---------------------------------------------------------
  function smlt_intrf_isnull(this) result(that)
    type(smlt_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(smlt_intrf_isnull)

    select case(this%type)
    case(SMLT_INTRF_DISA)
      that = .false.
    case(SMLT_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .true.
    case(SMLT_INTRF_ASSC)
      ASSERT(associated(this%self))
      that = .false.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(smlt_intrf_isnull)
  end function smlt_intrf_isnull

  ! ---------------------------------------------------------
  function smlt_intrf_assoc(this) result(that)
    type(smlt_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(smlt_intrf_assoc)

    select case(this%type)
    case(SMLT_INTRF_DISA)
      that = .false.
    case(SMLT_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(SMLT_INTRF_ASSC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(smlt_intrf_assoc)
  end function smlt_intrf_assoc

  ! ---------------------------------------------------------
  subroutine smlt_intrf_init_type(this, that)
    type(smlt_intrf_t),           intent(out) :: this
    type(simulation_t), optional, intent(in)  :: that

    PUSH_SUB(smlt_intrf_init_type)

    nullify(this%self)
    this%type = SMLT_INTRF_NULL
    call uuid_init(this%uuid)
    if(present(that)) call smlt_intrf_set(this, that)

    POP_SUB(smlt_intrf_init_type)
  end subroutine smlt_intrf_init_type

  ! ---------------------------------------------------------
  subroutine smlt_intrf_init_copy(this, that)
    type(smlt_intrf_t), intent(out) :: this
    type(smlt_intrf_t), intent(in)  :: that

    type(simulation_t), pointer :: sim

    PUSH_SUB(smlt_intrf_init_copy)

    ASSERT(smlt_intrf_assoc(that))
    nullify(sim)
    call smlt_intrf_get(that, sim)
    ASSERT(associated(sim))
    call smlt_intrf_init(this, sim)
    nullify(sim)

    POP_SUB(smlt_intrf_init_copy)
  end subroutine smlt_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine smlt_intrf_apply(this, operation)
    type(smlt_intrf_t), intent(inout) :: this

    interface
      subroutine operation(this)
        use simulation_oct_m
        type(simulation_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(simulation_t), pointer :: sim

    PUSH_SUB(smlt_intrf_apply)

    ASSERT(smlt_intrf_assoc(this))
    nullify(sim)
    call smlt_intrf_get(this, sim)
    ASSERT(associated(sim))
    call operation(sim)
    nullify(sim)

    POP_SUB(smlt_intrf_apply)
  end subroutine smlt_intrf_apply

  ! ---------------------------------------------------------
  subroutine smlt_intrf_set(this, that)
    type(smlt_intrf_t),         intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that

    PUSH_SUB(smlt_intrf_set)

    ASSERT(smlt_intrf_isnull(this))
    this%self => that
    call simulation_attach(this%self, this%uuid)
    this%type = SMLT_INTRF_ASSC

    POP_SUB(smlt_intrf_set)
  end subroutine smlt_intrf_set

  ! ---------------------------------------------------------
  subroutine smlt_intrf_get(this, that)
    type(smlt_intrf_t),          intent(in)  :: this
    type(simulation_t), pointer, intent(out) :: that

    PUSH_SUB(smlt_intrf_get)

    ASSERT(this%type>SMLT_INTRF_DISA)
    nullify(that)
    if(smlt_intrf_assoc(this)) that => this%self

    POP_SUB(smlt_intrf_get)
  end subroutine smlt_intrf_get

  ! ---------------------------------------------------------
  subroutine smlt_intrf_copy(this, that)
    type(smlt_intrf_t), intent(inout) :: this
    type(smlt_intrf_t), intent(in) :: that

    PUSH_SUB(smlt_intrf_copy)

    call smlt_intrf_end(this)
    this%type = that%type
    call uuid_init(this%uuid)
    if(smlt_intrf_assoc(that)) call smlt_intrf_set(this, that%self)

    POP_SUB(smlt_intrf_copy)
  end subroutine smlt_intrf_copy

  ! ---------------------------------------------------------
  subroutine smlt_intrf__end__(this)
    type(smlt_intrf_t), intent(inout) :: this

    PUSH_SUB(smlt_intrf__end__)

    nullify(this%self)
    this%type = SMLT_INTRF_DISA
    call uuid_end(this%uuid)

    POP_SUB(smlt_intrf__end__)
  end subroutine smlt_intrf__end__

  ! ---------------------------------------------------------
  subroutine smlt_intrf_end_type(this)
    type(smlt_intrf_t), intent(inout) :: this

    PUSH_SUB(smlt_intrf_end_type)

    call smlt_intrf_end(this, finis)

    POP_SUB(smlt_intrf_end_type)

  contains

    subroutine finis(this)
      type(simulation_t), intent(inout) :: this

      PUSH_SUB(smlt_intrf_end_type.finis)

      call simulation_end(this)

      POP_SUB(smlt_intrf_end_type.finis)
    end subroutine finis

  end subroutine smlt_intrf_end_type

  ! ---------------------------------------------------------
  subroutine smlt_intrf_end_pass(this, finis)
    type(smlt_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use simulation_oct_m
        type(simulation_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(smlt_intrf_end_pass)

    if(smlt_intrf_assoc(this)) call smlt_intrf_del(this, finis)
    call smlt_intrf__end__(this)

    POP_SUB(smlt_intrf_end_pass)
  end subroutine smlt_intrf_end_pass

end module smlt_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:
