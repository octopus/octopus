#include "template.h"

#if defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  use storage_m, only: &
    storage_init,      &
    storage_eval,      &
    storage_copy,      &
    storage_end

  use storage_m, only: &
    storage_intrpl_t

  use storage_m, only:                        &
    TEMPLATE(INTRPL_OK) => STORAGE_INTRPL_OK, &
    TEMPLATE(INTRPL_OD) => STORAGE_INTRPL_OD, &
    TEMPLATE(INTRPL_NI) => STORAGE_INTRPL_NI

#endif
#if !defined(INCLUDE_PREFIX) && defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  public ::         &
    TEMPLATE(eval)

  public ::              &
    TEMPLATE(INTRPL_OK), &
    TEMPLATE(INTRPL_OD), &
    TEMPLATE(INTRPL_NI)

  type, public :: TEMPLATE(intrpl_t)
    private
    type(TEMPLATE(t)), pointer :: self =>null()
    type(storage_intrpl_t)     :: intrp
  end type TEMPLATE(intrpl_t)

  interface TEMPLATE(_init__)
    module procedure TEMPLATE(intrpl_init)
  end interface TEMPLATE(_init__)

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
  subroutine TEMPLATE(intrpl_init)(this, that, type)
    type(TEMPLATE(intrpl_t)),  intent(out) :: this
    type(TEMPLATE(t)), target, intent(in)  :: that
    integer,         optional, intent(in)  :: type
    !
    PUSH_SUB(TEMPLATE(intrpl_init))
    this%self=>that
    call storage_init(this%intrp, that%data, type)
    POP_SUB(TEMPLATE(intrpl_init))
    return
  end subroutine TEMPLATE(intrpl_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_eval_1d)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    integer,                     intent(out) :: ierr
    !
    PUSH_SUB(TEMPLATE(intrpl_eval_1d))
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(TEMPLATE(intrpl_eval_1d))
    return
  end subroutine TEMPLATE(intrpl_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_eval_md)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: v
    integer,                     intent(out) :: ierr
    !
    PUSH_SUB(TEMPLATE(intrpl_eval_md))
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(TEMPLATE(intrpl_eval_md))
    return
  end subroutine TEMPLATE(intrpl_eval_md)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_get)(this, that)
    type(TEMPLATE(intrpl_t)), intent(in) :: this
    type(TEMPLATE(t)),       pointer     :: that
    !
    PUSH_SUB(TEMPLATE(intrpl_get))
    nullify(that)
    if(associated(this%self))&
      that=>this%self
    POP_SUB(TEMPLATE(intrpl_get))
    return
  end subroutine TEMPLATE(intrpl_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_copy)(this, that)
    type(TEMPLATE(intrpl_t)), intent(out) :: this
    type(TEMPLATE(intrpl_t)), intent(in)  :: that
    !
    PUSH_SUB(TEMPLATE(intrpl_copy))
    call TEMPLATE(intrpl_end)(this)
    if(associated(that%self))then
      this%self=>that%self
      call storage_copy(this%intrp, that%intrp)
    end if
    POP_SUB(TEMPLATE(intrpl_copy))
    return
  end subroutine TEMPLATE(intrpl_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_end)(this)
    type(TEMPLATE(intrpl_t)), intent(inout) :: this
    !
    PUSH_SUB(TEMPLATE(intrpl_end))
    if(associated(this%self))&
      call storage_end(this%intrp)
    nullify(this%self)
    POP_SUB(TEMPLATE(intrpl_end))
    return
  end subroutine TEMPLATE(intrpl_end)

#endif

!! Local Variables:
!! mode: f90
!! End:
