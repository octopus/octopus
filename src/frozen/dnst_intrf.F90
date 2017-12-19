#include "global.h"

module dnst_intrf_oct_m

  use dnst_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::       &
    dnst_intrf_t

  public ::           &
    dnst_intrf_new,   &
    dnst_intrf_del,   &
    dnst_intrf_assoc, &
    dnst_intrf_alloc, &
    dnst_intrf_init,  &
    dnst_intrf_set,   &
    dnst_intrf_get,   &
    dnst_intrf_copy,  &
    dnst_intrf_end

  integer, parameter :: DNST_INTRF_DISA = 0
  integer, parameter :: DNST_INTRF_NULL = 1
  integer, parameter :: DNST_INTRF_ASSC = 2
  integer, parameter :: DNST_INTRF_ALLC = 3

  type :: dnst_intrf_t
    private
    type(dnst_t), pointer :: self =>null()
    integer               :: type = DNST_INTRF_DISA
  end type dnst_intrf_t

  interface dnst_intrf_new
    module procedure dnst_intrf_new_type
    module procedure dnst_intrf_new_pass
    module procedure dnst_intrf_new_copy
  end interface dnst_intrf_new

  interface dnst_intrf_del
    module procedure dnst_intrf_del_type
    module procedure dnst_intrf_del_pass
  end interface dnst_intrf_del

  interface dnst_intrf_init
    module procedure dnst_intrf_init_type
    module procedure dnst_intrf_init_copy
  end interface dnst_intrf_init

  interface dnst_intrf_end
    module procedure dnst_intrf_end_type
    module procedure dnst_intrf_end_pass
  end interface dnst_intrf_end

contains

  ! ---------------------------------------------------------
  subroutine dnst_intrf__new__(this, that)
    type(dnst_intrf_t), intent(inout) :: this
    type(dnst_t),      pointer        :: that

    PUSH_SUB(dnst_intrf__new__)
    
    ASSERT(dnst_intrf_isnull(this))
    nullify(that)
    SAFE_ALLOCATE(that)
    call dnst_intrf_set(this, that)
    this%type = DNST_INTRF_ALLC

    POP_SUB(dnst_intrf__new__)
  end subroutine dnst_intrf__new__

  ! ---------------------------------------------------------
  subroutine dnst_intrf__del__(this)
    type(dnst_intrf_t), intent(inout) :: this

    PUSH_SUB(dnst_intrf__del__)

    ASSERT(dnst_intrf_assoc(this))
    if(dnst_intrf_alloc(this))then
      SAFE_DEALLOCATE_P(this%self)
    end if
    nullify(this%self)
    this%type = DNST_INTRF_NULL

    POP_SUB(dnst_intrf__del__)
  end subroutine dnst_intrf__del__

  ! ---------------------------------------------------------
  subroutine dnst_intrf_new_type(this, config)
    type(dnst_intrf_t),  intent(inout) :: this
    type(json_object_t), intent(in)    :: config

    PUSH_SUB(dnst_intrf_new_type)
    
    ASSERT(dnst_intrf_isnull(this))
    call dnst_intrf_new(this, init)

    POP_SUB(dnst_intrf_new_type)

  contains

    subroutine init(this)
      type(dnst_t), intent(out) :: this

      PUSH_SUB(dnst_intrf_new_type.init)

      call dnst_init(this, config)

      POP_SUB(dnst_intrf_new_type.init)
    end subroutine init

  end subroutine dnst_intrf_new_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_new_pass(this, init)
    type(dnst_intrf_t), intent(inout) :: this
    
    interface
      subroutine init(this)
        use dnst_oct_m
        type(dnst_t), intent(out) :: this
      end subroutine init
    end interface

    type(dnst_t), pointer :: self

    PUSH_SUB(dnst_intrf_new_pass)
    
    ASSERT(dnst_intrf_isnull(this))
    nullify(self)
    call dnst_intrf__new__(this, self)
    call init(self)
    nullify(self)

    POP_SUB(dnst_intrf_new_pass)
  end subroutine dnst_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine dnst_intrf_new_copy(this, source, mold)
    type(dnst_intrf_t),           intent(inout) :: this
    type(dnst_intrf_t), optional, intent(in)    :: source
    type(dnst_intrf_t), optional, intent(in)    :: mold

    type(dnst_t), pointer :: self

    PUSH_SUB(dnst_intrf_new_copy)

    ASSERT(present(source).or.present(mold))
    ASSERT(.not.(present(source).and.present(mold)))
    nullify(self)
    call dnst_intrf__new__(this, self)
    if(present(source))then
      ASSERT(dnst_intrf_assoc(source))
      call dnst_copy(self, source%self)
    elseif(present(mold))then
      ASSERT(dnst_intrf_assoc(mold))
      call dnst_init(self, mold%self)
    else
      ASSERT(.FALSE.)
    end if
    nullify(self)

    POP_SUB(dnst_intrf_new_copy)
  end subroutine dnst_intrf_new_copy

  ! ---------------------------------------------------------
  subroutine dnst_intrf_del_type(this)
    type(dnst_intrf_t), intent(inout) :: this

    PUSH_SUB(dnst_intrf_del_type)

    ASSERT(dnst_intrf_assoc(this))
    call dnst_intrf_del(this, dnst_end)

    POP_SUB(dnst_intrf_del_type)
  end subroutine dnst_intrf_del_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_del_pass(this, finis)
    type(dnst_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use dnst_oct_m
        type(dnst_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(dnst_intrf_del_pass)

    ASSERT(dnst_intrf_assoc(this))
    if(dnst_intrf_alloc(this)) call finis(this%self)
    call dnst_intrf__del__(this)
    
    POP_SUB(dnst_intrf_del_pass)
  end subroutine dnst_intrf_del_pass

  ! ---------------------------------------------------------
  function dnst_intrf_isnull(this) result(that)
    type(dnst_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(dnst_intrf_isnull)

    select case(this%type)
    case(DNST_INTRF_DISA)
      that = .false.
    case(DNST_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .true.
    case(DNST_INTRF_ASSC, DNST_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .false.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(dnst_intrf_isnull)
  end function dnst_intrf_isnull

  ! ---------------------------------------------------------
  function dnst_intrf_assoc(this) result(that)
    type(dnst_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(dnst_intrf_assoc)

    select case(this%type)
    case(DNST_INTRF_DISA)
      that = .false.
    case(DNST_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(DNST_INTRF_ASSC, DNST_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(dnst_intrf_assoc)
  end function dnst_intrf_assoc

  ! ---------------------------------------------------------
  function dnst_intrf_alloc(this) result(that)
    type(dnst_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(dnst_intrf_alloc)

    select case(this%type)
    case(DNST_INTRF_DISA)
      that = .false.
    case(DNST_INTRF_NULL)
      ASSERT(.not.associated(this%self))
      that = .false.
    case(DNST_INTRF_ASSC)
      ASSERT(associated(this%self))
      that = .false.
    case(DNST_INTRF_ALLC)
      ASSERT(associated(this%self))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(dnst_intrf_alloc)
  end function dnst_intrf_alloc

  ! ---------------------------------------------------------
  subroutine dnst_intrf_init_type(this)
    type(dnst_intrf_t), intent(out) :: this

    PUSH_SUB(dnst_intrf_init_type)
    
    nullify(this%self)
    this%type = DNST_INTRF_NULL

    POP_SUB(dnst_intrf_init_type)
  end subroutine dnst_intrf_init_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_init_copy(this, that)
    type(dnst_intrf_t), intent(out) :: this
    type(dnst_intrf_t), intent(in)  :: that

    PUSH_SUB(dnst_intrf_init_copy)
    
    call dnst_intrf_init(this)
    select case(that%type)
    case(DNST_INTRF_DISA)
    case(DNST_INTRF_NULL)
    case(DNST_INTRF_ASSC)
      call dnst_intrf_set(this, that%self)
    case(DNST_INTRF_ALLC)
      call dnst_intrf_new(this, mold=that)
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(dnst_intrf_init_copy)
  end subroutine dnst_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine dnst_intrf_set(this, that)
    type(dnst_intrf_t),   intent(inout) :: this
    type(dnst_t), target, intent(in)    :: that

    PUSH_SUB(dnst_intrf_set)

    ASSERT(dnst_intrf_isnull(this))
    this%self => that
    this%type = DNST_INTRF_ASSC

    POP_SUB(dnst_intrf_set)
  end subroutine dnst_intrf_set

  ! ---------------------------------------------------------
  subroutine dnst_intrf_get(this, that)
    type(dnst_intrf_t), intent(in) :: this
    type(dnst_t),      pointer     :: that

    PUSH_SUB(dnst_intrf_get)

    ASSERT(this%type/=DNST_INTRF_DISA)
    nullify(that)
    if(dnst_intrf_assoc(this)) that => this%self

    POP_SUB(dnst_intrf_get)
  end subroutine dnst_intrf_get

  ! ---------------------------------------------------------
  subroutine dnst_intrf_copy(this, that)
    type(dnst_intrf_t), intent(inout) :: this
    type(dnst_intrf_t), intent(in) :: that

    PUSH_SUB(dnst_intrf_copy)

    call dnst_intrf_end(this)
    select case(that%type)
    case(DNST_INTRF_DISA)
    case(DNST_INTRF_NULL)
    case(DNST_INTRF_ASSC)
      call dnst_intrf_set(this, that%self)
    case(DNST_INTRF_ALLC)
      call dnst_intrf_new(this, source=that)
    case default
      ASSERT(.false.)
    end select
    this%type = that%type

    POP_SUB(dnst_intrf_copy)
  end subroutine dnst_intrf_copy

  ! ---------------------------------------------------------
  subroutine dnst_intrf__end__(this)
    type(dnst_intrf_t), intent(inout) :: this

    PUSH_SUB(dnst_intrf__end__)

    nullify(this%self)
    this%type = DNST_INTRF_DISA

    POP_SUB(dnst_intrf__end__)
  end subroutine dnst_intrf__end__

  ! ---------------------------------------------------------
  subroutine dnst_intrf_end_type(this)
    type(dnst_intrf_t), intent(inout) :: this

    PUSH_SUB(dnst_intrf_end_type)

    call dnst_intrf_end(this, dnst_end)

    POP_SUB(dnst_intrf_end_type)
  end subroutine dnst_intrf_end_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_end_pass(this, finis)
    type(dnst_intrf_t), intent(inout) :: this

    interface
      subroutine finis(this)
        use dnst_oct_m
        type(dnst_t), intent(inout) :: this
      end subroutine finis
    end interface

    PUSH_SUB(dnst_intrf_end_pass)

    if(dnst_intrf_assoc(this)) call dnst_intrf_del(this, finis)
    call dnst_intrf__end__(this)

    POP_SUB(dnst_intrf_end_pass)
  end subroutine dnst_intrf_end_pass

end module dnst_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:

