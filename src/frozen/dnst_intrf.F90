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

  integer, parameter :: DNST_NULL = 0
  integer, parameter :: DNST_ASSC = 1
  integer, parameter :: DNST_ALLC = 2

  type :: dnst_intrf_t
    private
    type(json_object_t), pointer :: config =>null()
    type(dnst_t),        pointer :: pdns   =>null()
    integer                      :: type   = DNST_NULL
  end type dnst_intrf_t

  interface dnst_intrf_new
    module procedure dnst_intrf_new_type
    module procedure dnst_intrf_new_pass
  end interface dnst_intrf_new

  interface dnst_intrf_init
    module procedure dnst_intrf_init_type
    module procedure dnst_intrf_init_copy
  end interface dnst_intrf_init

  interface dnst_intrf_get
    module procedure dnst_intrf_get_config
    module procedure dnst_intrf_get_dnst
  end interface dnst_intrf_get

contains

  ! ---------------------------------------------------------
  subroutine dnst_intrf__new__(this, that)
    type(dnst_intrf_t), intent(inout) :: this
    type(dnst_t),      pointer        :: that

    PUSH_SUB(dnst_intrf__new__)
    
    ASSERT(associated(this%config))
    nullify(that)
    SAFE_ALLOCATE(that)
    call dnst_intrf_set(this, that)
    this%type = DNST_ALLC

    POP_SUB(dnst_intrf__new__)
  end subroutine dnst_intrf__new__

  ! ---------------------------------------------------------
  subroutine dnst_intrf_new_type(this, that)
    type(dnst_intrf_t), intent(inout) :: this
    type(dnst_t),      pointer        :: that

    PUSH_SUB(dnst_intrf_new_type)
    
    nullify(that)
    call dnst_intrf__new__(this, that)

    POP_SUB(dnst_intrf_new_type)
  end subroutine dnst_intrf_new_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_new_pass(this, type_init)
    type(dnst_intrf_t), intent(inout) :: this
    
    interface
      subroutine type_init(this, config)
        use dnst_oct_m
        use json_oct_m
        type(dnst_t),        intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine type_init
    end interface

    type(dnst_t), pointer :: dnst

    PUSH_SUB(dnst_intrf_new_pass)
    
    nullify(dnst)
    call dnst_intrf__new__(this, dnst)
    call type_init(dnst, this%config)

    POP_SUB(dnst_intrf_new_pass)
  end subroutine dnst_intrf_new_pass

  ! ---------------------------------------------------------
  subroutine dnst_intrf_del(this, type_end)
    type(dnst_intrf_t), intent(inout) :: this

    interface
      subroutine type_end(this)
        use dnst_oct_m
        type(dnst_t), intent(inout) :: this
      end subroutine type_end
    end interface

    PUSH_SUB(dnst_intrf_del)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this))
    ASSERT(this%type==DNST_ALLC)
    call type_end(this%pdns)
    SAFE_DEALLOCATE_P(this%pdns)
    nullify(this%pdns)
    this%type = DNST_NULL
    
    POP_SUB(dnst_intrf_del)
  end subroutine dnst_intrf_del

  ! ---------------------------------------------------------
  function dnst_intrf_assoc(this) result(that)
    type(dnst_intrf_t), intent(in) :: this

    logical :: that

    PUSH_SUB(dnst_intrf_assoc)

    ASSERT(associated(this%config))
    select case(this%type)
    case(DNST_NULL)
      ASSERT(.not.associated(this%pdns))
      that = .false.
    case(DNST_ASSC,DNST_ALLC)
      ASSERT(associated(this%pdns))
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

    ASSERT(associated(this%config))
    select case(this%type)
    case(DNST_NULL)
      ASSERT(.not.associated(this%pdns))
      that = .false.
    case(DNST_ASSC)
      ASSERT(associated(this%pdns))
      that = .false.
    case(DNST_ALLC)
      ASSERT(associated(this%pdns))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(dnst_intrf_alloc)
  end function dnst_intrf_alloc

  ! ---------------------------------------------------------
  subroutine dnst_intrf_init_type(this, config)
    type(dnst_intrf_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(dnst_intrf_init_type)
    
    this%config => config
    nullify(this%pdns)
    this%type = DNST_NULL

    POP_SUB(dnst_intrf_init_type)
  end subroutine dnst_intrf_init_type

  ! ---------------------------------------------------------
  subroutine dnst_intrf_init_copy(this, that)
    type(dnst_intrf_t), intent(out) :: this
    type(dnst_intrf_t), intent(in)  :: that

    PUSH_SUB(dnst_intrf_init_copy)
    
    ASSERT(associated(that%config))
    call dnst_intrf_init(this, that%config)

    POP_SUB(dnst_intrf_init_copy)
  end subroutine dnst_intrf_init_copy

  ! ---------------------------------------------------------
  subroutine dnst_intrf_set(this, that)
    type(dnst_intrf_t),   intent(inout) :: this
    type(dnst_t), target, intent(in)    :: that

    PUSH_SUB(dnst_intrf_set)

    ASSERT(associated(this%config))
    ASSERT(.not.dnst_intrf_assoc(this))
    ASSERT(this%type==DNST_NULL)
    this%pdns => that
    this%type = DNST_ASSC

    POP_SUB(dnst_intrf_set)
  end subroutine dnst_intrf_set

  ! ---------------------------------------------------------
  subroutine dnst_intrf_get_config(this, that)
    type(dnst_intrf_t),   intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(dnst_intrf_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(dnst_intrf_get_config)
  end subroutine dnst_intrf_get_config

  ! ---------------------------------------------------------
  subroutine dnst_intrf_get_dnst(this, that)
    type(dnst_intrf_t), intent(in) :: this
    type(dnst_t),      pointer     :: that

    PUSH_SUB(dnst_intrf_get_dnst)

    nullify(that)
    if(dnst_intrf_assoc(this)) that => this%pdns

    POP_SUB(dnst_intrf_get_dnst)
  end subroutine dnst_intrf_get_dnst

  ! ---------------------------------------------------------
  subroutine dnst_intrf_copy(this, that, type_copy)
    type(dnst_intrf_t), intent(inout) :: this
    type(dnst_intrf_t), intent(in)    :: that

    interface
      subroutine type_copy(this, that)
        use dnst_oct_m
        type(dnst_t), intent(inout) :: this
        type(dnst_t), intent(in)    :: that
      end subroutine type_copy
    end interface

    type(dnst_t), pointer :: pdns
    integer               :: type

    PUSH_SUB(dnst_intrf_copy)

    nullify(pdns)
    type = this%type
    pdns => this%pdns
    nullify(this%config, this%pdns)
    this%type = DNST_NULL
    if(associated(that%config))then
      call dnst_intrf_init(this, that)
      select case(that%type)
      case(DNST_NULL)
        ASSERT(type/=DNST_ALLC)
      case(DNST_ASSC)
        ASSERT(type/=DNST_ALLC)
        call dnst_intrf_set(this, that%pdns)
      case(DNST_ALLC)
        if(type/=DNST_ALLC) call dnst_intrf__new__(this, pdns)
        call type_copy(pdns, that%pdns)
      case default
        ASSERT(.false.)
      end select
    else
      ASSERT(type/=DNST_ALLC)
    end if
    nullify(pdns)

    POP_SUB(dnst_intrf_copy)
  end subroutine dnst_intrf_copy

  ! ---------------------------------------------------------
  subroutine dnst_intrf_end(this, type_end)
    type(dnst_intrf_t), intent(inout) :: this

    interface
      subroutine type_end(this)
        use dnst_oct_m
        type(dnst_t), intent(inout) :: this
      end subroutine type_end
    end interface

    PUSH_SUB(dnst_intrf_end)

    if(this%type==DNST_ALLC) call dnst_intrf_del(this, type_end)
    nullify(this%config, this%pdns)
    this%type = DNST_NULL

    POP_SUB(dnst_intrf_end)
  end subroutine dnst_intrf_end

end module dnst_intrf_oct_m

!! Local Variables:
!! mode: f90
!! End:

