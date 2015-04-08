#include "template.h"

#if defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  use intrpl_m, only: &
    intrpl_t,         &
    intrpl_init,      &
    intrpl_eval,      &
    intrpl_copy,      &
    intrpl_end

  use intrpl_m, only: &
    INTRPL_OK,        &
    INTRPL_OD,        &
    INTRPL_NI

#endif
#if !defined(INCLUDE_PREFIX) && defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  public ::             &
    TEMPLATE(intrpl_t)

  public ::         &
    TEMPLATE(eval)

  public ::              &
    TEMPLATE(INTRPL_OK), &
    TEMPLATE(INTRPL_OD), &
    TEMPLATE(INTRPL_NI)

  integer, parameter :: TEMPLATE(INTRPL_OK) = INTRPL_OK
  integer, parameter :: TEMPLATE(INTRPL_OD) = INTRPL_OD
  integer, parameter :: TEMPLATE(INTRPL_NI) = INTRPL_NI

  type :: TEMPLATE(intrpl_t)
    private
    type(TEMPLATE(t)), pointer :: self =>null()
    type(intrpl_t)             :: intrp
  end type TEMPLATE(intrpl_t)

  interface TEMPLATE(_init__)
    module procedure TEMPLATE(intrpl_init)
  end interface TEMPLATE(_init__)

#if 0

  interface TEMPLATE(_copy__)
    module procedure TEMPLATE(intrpl_copy)
  end interface TEMPLATE(_copy__)

#endif
  
  interface TEMPLATE(_end__)
    module procedure TEMPLATE(intrpl_end)
  end interface TEMPLATE(_end__)

  interface TEMPLATE(init)
    module procedure TEMPLATE(intrpl_init)
  end interface TEMPLATE(init)

  interface TEMPLATE(get)
    module procedure TEMPLATE(intrpl_get)
  end interface TEMPLATE(get)

  interface TEMPLATE(eval)
    module procedure TEMPLATE(intrpl_eval_1d)
    module procedure TEMPLATE(intrpl_eval_md)
  end interface TEMPLATE(eval)

  interface TEMPLATE(copy)
    module procedure TEMPLATE(intrpl_copy)
  end interface TEMPLATE(copy)

  interface TEMPLATE(end)
    module procedure TEMPLATE(intrpl_end)
  end interface TEMPLATE(end)

#endif
#if !defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && defined(INCLUDE_BODY)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_init)(this, that, type, default)
    type(TEMPLATE(intrpl_t)),  intent(out) :: this
    type(TEMPLATE(t)), target, intent(in)  :: that
    integer,         optional, intent(in)  :: type
    real(kind=wp),   optional, intent(in)  :: default

    real(kind=wp), dimension(:,:), pointer :: data

    PUSH_SUB(TEMPLATE(intrpl_init))

    nullify(data)
    this%self=>that
    call TEMPLATE(get)(that, data)
    ASSERT(associated(data))
    call intrpl_init(this%intrp, that%sim, data, type=type, default=default)

    POP_SUB(TEMPLATE(intrpl_init))
  end subroutine TEMPLATE(intrpl_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_eval_1d)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    integer,                     intent(out) :: ierr

    PUSH_SUB(TEMPLATE(intrpl_eval_1d))

    ierr=TEMPLATE(INTRPL_NI)
    if(associated(this%self))call intrpl_eval(this%intrp, x, v, ierr)

    POP_SUB(TEMPLATE(intrpl_eval_1d))
  end subroutine TEMPLATE(intrpl_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_eval_md)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: v
    integer,                     intent(out) :: ierr

    PUSH_SUB(TEMPLATE(intrpl_eval_md))

    ierr=TEMPLATE(INTRPL_NI)
    if(associated(this%self))call intrpl_eval(this%intrp, x, v, ierr)

    POP_SUB(TEMPLATE(intrpl_eval_md))
  end subroutine TEMPLATE(intrpl_eval_md)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_get)(this, that)
    type(TEMPLATE(intrpl_t)), intent(in) :: this
    type(TEMPLATE(t)),       pointer     :: that

    PUSH_SUB(TEMPLATE(intrpl_get))

    nullify(that)
    if(associated(this%self))&
      that=>this%self

    POP_SUB(TEMPLATE(intrpl_get))
  end subroutine TEMPLATE(intrpl_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_copy)(this, that)
    type(TEMPLATE(intrpl_t)), intent(out) :: this
    type(TEMPLATE(intrpl_t)), intent(in)  :: that

    PUSH_SUB(TEMPLATE(intrpl_copy))

    call TEMPLATE(intrpl_end)(this)
    if(associated(that%self))then
      this%self=>that%self
      call intrpl_copy(this%intrp, that%intrp)
    end if

    POP_SUB(TEMPLATE(intrpl_copy))
  end subroutine TEMPLATE(intrpl_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_end)(this)
    type(TEMPLATE(intrpl_t)), intent(inout) :: this

    PUSH_SUB(TEMPLATE(intrpl_end))

    if(associated(this%self))&
      call intrpl_end(this%intrp)
    nullify(this%self)

    POP_SUB(TEMPLATE(intrpl_end))
  end subroutine TEMPLATE(intrpl_end)

#endif

!! Local Variables:
!! mode: f90
!! End:
