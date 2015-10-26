#include "template.h"

#if defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  use intrpl_m
  use storage_m

#endif
#if !defined(INCLUDE_PREFIX) && defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  public ::             &
    TEMPLATE(intrpl_t)

  integer, parameter :: TEMPLATE(INTRPL_OK) = INTRPL_OK
  integer, parameter :: TEMPLATE(INTRPL_OD) = INTRPL_OD
  integer, parameter :: TEMPLATE(INTRPL_NI) = INTRPL_NI

  type :: TEMPLATE(intrpl_t)
    private
    type(EXTERNAL(t)), pointer :: self =>null()
    type(intrpl_t)             :: intrp
  end type TEMPLATE(intrpl_t)

  interface TEMPLATE(intrpl__eval__)
    module procedure INTERNAL(intrpl__eval__1d)
    module procedure INTERNAL(intrpl__eval__md)
  end interface TEMPLATE(intrpl__eval__)

  interface TEMPLATE(init)
    module procedure TEMPLATE(intrpl_init)
  end interface TEMPLATE(init)

  interface TEMPLATE(get)
    module procedure TEMPLATE(intrpl_get_info)
    module procedure TEMPLATE(intrpl_get_type)
  end interface TEMPLATE(get)

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
    type(EXTERNAL(t)), target, intent(in)  :: that
    integer,         optional, intent(in)  :: type

    type(storage_t), pointer :: data

    PUSH_SUB(TEMPLATE(intrpl_init))

    nullify(data)
    this%self => that
    call EXTERNAL(get)(that, data)
    ASSERT(associated(data))
    call intrpl_init(this%intrp, data, type=type)
    nullify(data)

    POP_SUB(TEMPLATE(intrpl_init))
  end subroutine TEMPLATE(intrpl_init)

  ! ---------------------------------------------------------
  subroutine INTERNAL(intrpl__eval__1d)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    integer,                     intent(out) :: ierr

    PUSH_SUB(INTERNAL(intrpl__eval__1d))

    ierr = TEMPLATE(INTRPL_NI)
    if(associated(this%self)) call intrpl_eval(this%intrp, x, v, ierr)

    POP_SUB(INTERNAL(intrpl__eval__1d))
  end subroutine INTERNAL(intrpl__eval__1d)

  ! ---------------------------------------------------------
  subroutine INTERNAL(intrpl__eval__md)(this, x, v, ierr)
    type(TEMPLATE(intrpl_t)),    intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: v
    integer,                     intent(out) :: ierr

    PUSH_SUB(INTERNAL(intrpl__eval__md))

    ierr = TEMPLATE(INTRPL_NI)
    if(associated(this%self)) call intrpl_eval(this%intrp, x, v, ierr)

    POP_SUB(INTERNAL(intrpl__eval__md))
  end subroutine INTERNAL(intrpl__eval__md)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_get_info)(this, type, dim, default)
    type(TEMPLATE(intrpl_t)), intent(in)  :: this
    integer,        optional, intent(out) :: type
    integer,        optional, intent(out) :: dim
    real(kind=wp),  optional, intent(out) :: default

    PUSH_SUB(TEMPLATE(intrpl_get_info))

    call intrpl_get(this%intrp, type, dim, default)

    POP_SUB(TEMPLATE(intrpl_get_info))
  end subroutine TEMPLATE(intrpl_get_info)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_get_type)(this, that)
    type(TEMPLATE(intrpl_t)), intent(in) :: this
    type(EXTERNAL(t)),       pointer     :: that

    PUSH_SUB(TEMPLATE(intrpl_get_type))

    nullify(that)
    if(associated(this%self)) that => this%self

    POP_SUB(TEMPLATE(intrpl_get_type))
  end subroutine TEMPLATE(intrpl_get_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(intrpl_copy)(this, that)
    type(TEMPLATE(intrpl_t)), intent(out) :: this
    type(TEMPLATE(intrpl_t)), intent(in)  :: that

    PUSH_SUB(TEMPLATE(intrpl_copy))

    call TEMPLATE(intrpl_end)(this)
    if(associated(that%self))then
      this%self => that%self
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

#undef TEMPLATE_PREFIX

!! Local Variables:
!! mode: f90
!! End:
