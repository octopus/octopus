#include "global.h"

module intrpl_m

  use curvilinear_m
  use domain_m
  use global_m
  use grid_m
  use index_m
  use kinds_m
  use mesh_m
  use messages_m
  use profiling_m
  use qshep_m
  use simulation_m
  use storage_m

  implicit none

  private

  public ::   &
    NONE,     &
    NEAREST,  &
    QSHEP

  public ::    &
    INTRPL_OK, &
    INTRPL_OD, &
    INTRPL_NI

  public ::   &
    intrpl_t

  public ::      &
    intrpl_init, &
    intrpl_eval, &
    intrpl_get,  &
    intrpl_copy, &
    intrpl_end

  integer, parameter :: NONE    = 0
  integer, parameter :: NEAREST = 1
  integer, parameter :: QSHEP   = 2

  integer, parameter :: INTRPL_OK = 0
  integer, parameter :: INTRPL_OD = 1
  integer, parameter :: INTRPL_NI = 2

  type :: intrp_t
    private
    integer                              :: type  = NONE
    type(qshep_t),               pointer :: qshep =>null()
    real(kind=wp), dimension(:), pointer :: vals  =>null()
  end type intrp_t

  type :: intrpl_t
    private
    type(simulation_t),              pointer :: sim     =>null()
    type(mesh_t),                    pointer :: mesh    =>null()
    type(domain_t),                  pointer :: domain  =>null()
    integer                                  :: type    = NONE
    integer                                  :: nint    = 0
    real(kind=wp)                            :: default = 0.0_wp
    type(intrp_t), dimension(:), allocatable :: intr
  end type intrpl_t

  interface intrpl_eval
    module procedure intrpl_eval_1d
    module procedure intrpl_eval_md
  end interface intrpl_eval

  interface intrpl_get
    module procedure intrpl_get_info
  end interface intrpl_get

contains

  ! ---------------------------------------------------------
  subroutine intrp_init(this, mesh, vals, type)
    type(intrp_t),                       intent(out) :: this
    type(mesh_t),                        intent(in)  :: mesh
    real(kind=wp), dimension(:), target, intent(in)  :: vals
    integer,                             intent(in)  :: type

    integer :: nsz

    PUSH_SUB(intrp_init)

    this%type = type
    this%vals => vals
    ASSERT(associated(this%vals))
    this%qshep => null()
    select case(this%type)
    case(QSHEP)
      nsz = size(vals)
      ASSERT(nsz>10)
      SAFE_ALLOCATE(this%qshep)
      select case(mesh%sb%dim)
      case(2)
        call qshep_init(this%qshep, nsz, vals, mesh%x(1:nsz,1), mesh%x(1:nsz,2))
      case(3)
        call qshep_init(this%qshep, nsz, vals, mesh%x(1:nsz,1), mesh%x(1:nsz,2), mesh%x(1:nsz,3))
      case default
        message(1)='Quadratic Shepard interpolation only works in two or three dimensions.'
        call messages_fatal(1)
      end select
    end select

    POP_SUB(intrp_init)
  end subroutine intrp_init

  ! ---------------------------------------------------------
  subroutine intrp_qshep(this, x, val)
    type(intrp_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    PUSH_SUB(intrp_qshep)

    val = qshep_interpolate(this%qshep, this%vals, x)

    POP_SUB(intrp_qshep)
  end subroutine intrp_qshep

  ! ---------------------------------------------------------
  subroutine intrp_eval(this, x, val)
    type(intrp_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    PUSH_SUB(intrp_eval)

    val = 0.0_wp
    select case(this%type)
    case(QSHEP)
      call intrp_qshep(this, x, val)
    end select

    POP_SUB(intrp_eval)
  end subroutine intrp_eval

  ! ---------------------------------------------------------
  subroutine intrp_copy(this, that)
    type(intrp_t), intent(out) :: this
    type(intrp_t), intent(in)  :: that

    PUSH_SUB(intrp_copy)

    this%type = that%type
    this%vals => that%vals
    this%qshep => that%qshep

    POP_SUB(intrp_copy)
  end subroutine intrp_copy

  ! ---------------------------------------------------------
  subroutine intrp_end(this)
    type(intrp_t), intent(inout) :: this

    PUSH_SUB(intrp_end)

    select case(this%type)
    case(QSHEP)
      call qshep_end(this%qshep)
      SAFE_DEALLOCATE_P(this%qshep)
    end select
    nullify(this%qshep, this%vals)
    this%type = NONE

    POP_SUB(intrp_end)
  end subroutine intrp_end

  ! ---------------------------------------------------------
  subroutine intrpl_init(this, that, type)
    type(intrpl_t),    intent(out) :: this
    type(storage_t),   intent(in)  :: that
    integer, optional, intent(in)  :: type

    real(kind=wp), dimension(:,:), pointer :: data
    integer                                :: indx
    logical                                :: allc

    PUSH_SUB(intrpl_init)

    nullify(data)
    ASSERT(.not.associated(this%sim))
    ASSERT(.not.associated(this%mesh))
    ASSERT(.not.associated(this%domain))
    call storage_get(that, this%sim)
    ASSERT(associated(this%sim))
    call storage_get(that, this%mesh)
    ASSERT(associated(this%mesh))
    call simulation_get(this%sim, this%domain)
    ASSERT(associated(this%domain))
    call storage_get(that, dim=this%nint, default=this%default, alloc=allc)
    ASSERT(this%nint>0)
    if(allc)then
      this%type = NEAREST
      if(present(type)) this%type = type
    else
      this%type = NONE
    end if
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
      call storage_get(that, data)
      ASSERT(associated(data))
      do indx = 1, this%nint
        call intrp_init(this%intr(indx), this%mesh, data(:,indx), this%type)
      end do
    end if

    POP_SUB(intrpl_init)
  end subroutine intrpl_init

  ! ---------------------------------------------------------
  function intrpl_nearest_index(this, x) result(n)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x

    integer :: n

    real(kind=wp), dimension(MAX_DIM) :: xp, chi
    integer,       dimension(MAX_DIM) :: ix
    integer                           :: dm

    PUSH_SUB(intrpl_nearest_index)

    dm = this%mesh%sb%dim
    xp(1:dm) = x(1:dm)
    xp(dm+1:MAX_DIM) = 0.0_wp
    call curvilinear_x2chi(this%mesh%sb, this%mesh%cv, xp, chi)
    ix(1:dm) = nint(chi(1:dm)/this%mesh%spacing(1:dm))
    ix(dm+1:MAX_DIM) = 0
    n = index_from_coords(this%mesh%idx, ix)

    POP_SUB(intrpl_nearest_index)
  end function intrpl_nearest_index

  ! ---------------------------------------------------------
  function intrpl_in_domain(this, x) result(in)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x

    logical :: in

    integer :: n

    PUSH_SUB(intrpl_in_domain)

    in = .true.
    if(associated(this%domain))&
      in = domain_in(this%domain, x)
    if(in)then
      n = intrpl_nearest_index(this, x)
      in = ( (0<n) .and. (n<=size(this%intr(1)%vals)) )
    end if

    POP_SUB(intrpl_in_domain)
  end function intrpl_in_domain

  ! ---------------------------------------------------------
  subroutine intrpl_nearest(this, x, val)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val

    real(kind=wp) :: tol, dlt
    integer       :: ix, np, dm

    PUSH_SUB(intrpl_nearest)

    dm = this%mesh%sb%dim
    np = intrpl_nearest_index(this, x)
    tol = 0.51_wp * sqrt(sum(this%mesh%spacing(1:dm)**2))
    dlt = sqrt(sum((x(1:dm)-this%mesh%x(np,1:dm))**2))
    ASSERT(dlt<tol)
    forall(ix=1:this%nint) val(ix) = this%intr(ix)%vals(np)

    POP_SUB(intrpl_nearest)
  end subroutine intrpl_nearest

  ! ---------------------------------------------------------
  subroutine intrpl_eval_internal(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr

    integer :: idx

    PUSH_SUB(intrpl_eval_internal)

    if(intrpl_in_domain(this, x))then
      ierr = INTRPL_OK
      select case(this%type)
      case(NEAREST)
        call intrpl_nearest(this, x, val)
      case(QSHEP)
        do idx = 1, this%nint
          call intrp_eval(this%intr(idx), x, val(idx))
        end do
      case default
        val = this%default
        ierr = INTRPL_NI
      end select
    else
      val = this%default
      ierr = INTRPL_OD
    end if

    POP_SUB(intrpl_eval_internal)
  end subroutine intrpl_eval_internal

  ! ---------------------------------------------------------
  subroutine intrpl_eval_1d(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: val
    integer,                      intent(out) :: ierr

    real(kind=wp), dimension(1) :: tvl

    PUSH_SUB(intrpl_eval_1d)

    if(this%type>NONE)then
      call intrpl_eval_internal(this, x, tvl, ierr)
      val = tvl(1)
    else
      val = this%default
      ierr = INTRPL_NI
    end if

    POP_SUB(intrpl_eval_1d)
  end subroutine intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine intrpl_eval_md(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr

    PUSH_SUB(intrpl_eval_md)

    if(this%type>NONE)then
      call intrpl_eval_internal(this, x, val, ierr)
    else
      val = this%default
      ierr = INTRPL_NI
    end if

    POP_SUB(intrpl_eval_md)
  end subroutine intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine intrpl_get_info(this, type, dim, default)
    type(intrpl_t),          intent(in)  :: this
    integer,       optional, intent(out) :: type
    integer,       optional, intent(out) :: dim
    real(kind=wp), optional, intent(out) :: default

    PUSH_SUB(intrpl_get_info)

    if(present(type)) type = this%type
    if(present(dim)) dim = this%nint
    if(present(default)) default = this%default

    POP_SUB(intrpl_get_info)
  end subroutine intrpl_get_info

  ! ---------------------------------------------------------
  subroutine intrpl_copy(this, that)
    type(intrpl_t),         intent(out) :: this
    type(intrpl_t), target, intent(in)  :: that

    integer :: idx

    PUSH_SUB(intrpl_copy)
    this%sim => that%sim
    this%mesh => that%mesh
    this%domain => that%domain
    this%type = that%type
    this%nint = that%nint
    this%default = that%default
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
      do idx = 1, this%nint
        call intrp_copy(this%intr(idx), that%intr(idx))
      end do
    end if

    POP_SUB(intrpl_copy)
  end subroutine intrpl_copy

  ! ---------------------------------------------------------
  subroutine intrpl_end(this)
    type(intrpl_t), intent(inout) :: this

    integer :: idx

    PUSH_SUB(intrpl_end)

    if(this%type>NONE)then
      do idx = 1, this%nint
        call intrp_end(this%intr(idx))
      end do
      SAFE_DEALLOCATE_A(this%intr)
    end if
    this%default = 0.0_wp
    this%nint = 0
    this%type = NONE
    nullify(this%domain, this%mesh, this%sim)

    POP_SUB(intrpl_end)
  end subroutine intrpl_end

end module intrpl_m

!! Local Variables:
!! mode: f90
!! End:
