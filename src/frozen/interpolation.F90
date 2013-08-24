#include "global.h"

module interpolation_m

  use global_m
  use messages_m
  use profiling_m

  use basis_m,       only: basis_t, basis_to_internal
  use curvilinear_m, only: curvilinear_x2chi
  use domain_m,      only: domain_t, domain_in_domain
  use index_m,       only: index_from_coords
  use kinds_m,       only: wp
  use mesh_m,        only: mesh_t
  use qshepmod_m,    only: qshep_t, init_qshep, qshep_interpolate, kill_qshep
  use simulation_m,  only: simulation_t, simulation_get

  implicit none

  private
  public ::                      &
    interpolation_init,          &
    interpolation_get_dimension, &
    interpolation_eval,          &
    interpolation_copy,          &
    interpolation_end

  integer, public, parameter :: NONE    = 0
  integer, public, parameter :: NEAREST = 1
  integer, public, parameter :: QSHEP   = 2

  integer, public, parameter :: INTRP_OK = 0
  integer, public, parameter :: INTRP_OD = 1
  integer, public, parameter :: INTRP_NI = 2

  type, private :: intrp_t
    private
    integer                              :: type  = NONE
    type(qshep_t),               pointer :: qshep =>null()
    real(kind=wp), dimension(:), pointer :: vals  =>null()
  end type intrp_t

  type, public :: interpolation_t
    private
    integer                              :: type    = NONE
    integer                              :: nint    = 0
    type(simulation_t),          pointer :: sim     =>null()
    type(mesh_t),                pointer :: mesh    =>null()
    type(domain_t),              pointer :: domain  =>null()
    type(basis_t),               pointer :: basis   =>null()
    type(intrp_t), dimension(:), pointer :: intr    =>null()
    real(kind=wp)                        :: default = 0.0_wp
  end type interpolation_t

  interface interpolation_init
    module procedure interpolation_init_1d
    module procedure interpolation_init_2d
  end interface interpolation_init

  interface interpolation_eval
    module procedure interpolation_eval_1d
    module procedure interpolation_eval_2d
  end interface interpolation_eval

contains

  ! ---------------------------------------------------------
  subroutine intrp_init(this, mesh, vals, type)
    type(intrp_t),                       intent(out) :: this
    type(mesh_t),                        intent(in)  :: mesh
    real(kind=wp), dimension(:), target, intent(in)  :: vals
    integer,                             intent(in)  :: type
    !
    integer :: n
    !
    PUSH_SUB(intrp_init)
    this%type  = type
    this%vals  =>vals
    ASSERT(associated(this%vals))
    this%qshep =>null()
    select case(this%type)
    case(QSHEP)
      n=size(vals)
      ASSERT(n>10)
      SAFE_ALLOCATE(this%qshep)
      select case(mesh%sb%dim)
      case(2)
        call init_qshep(this%qshep, n, vals, mesh%x(1:n,1), mesh%x(1:n,2))
      case(3)
        call init_qshep(this%qshep, n, vals, mesh%x(1:n,1), mesh%x(1:n,2), mesh%x(1:n,3))
      case default
        message(1)='Quadratic Shepard interpolation only works in two or three dimensions.'
        call messages_fatal(1)
      end select
    end select
    POP_SUB(intrp_init)
    return
  end subroutine intrp_init

  ! ---------------------------------------------------------
  subroutine intrp_qshep(this, x, val)
    type(intrp_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val
    !
    PUSH_SUB(intrp_qshep)
    val=qshep_interpolate(this%qshep, this%vals, x)
    POP_SUB(intrp_qshep)
    return
  end subroutine intrp_qshep

  ! ---------------------------------------------------------
  subroutine intrp_eval(this, x, val)
    type(intrp_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val
    !
    PUSH_SUB(intrp_eval)
    val=0.0_wp
    select case(this%type)
    case(QSHEP)
      call intrp_qshep(this, x, val)
    end select
    POP_SUB(intrp_eval)
    return
  end subroutine intrp_eval

  ! ---------------------------------------------------------
  subroutine intrp_copy(this, that)
    type(intrp_t), intent(out) :: this
    type(intrp_t), intent(in)  :: that
    !
    PUSH_SUB(intrp_copy)
    this%type  = that%type
    this%vals  =>that%vals
    this%qshep =>that%qshep
    POP_SUB(intrp_copy)
    return
  end subroutine intrp_copy

  ! ---------------------------------------------------------
  subroutine intrp_end(this)
    type(intrp_t), intent(inout) :: this
    !
    PUSH_SUB(intrp_end)
    select case(this%type)
    case(QSHEP)
      call kill_qshep(this%qshep)
      SAFE_DEALLOCATE_P(this%qshep)
    end select
    this%qshep =>null()
    this%vals  =>null()
    this%type  = NONE
    POP_SUB(intrp_end)
    return
  end subroutine intrp_end

  ! ---------------------------------------------------------
  subroutine interpolation_init_common(this, sim, nint, basis, type, default)
    type(interpolation_t),           intent(out) :: this
    type(simulation_t),      target, intent(in)  :: sim
    integer,                         intent(in)  :: nint
    type(basis_t), optional, target, intent(in)  :: basis
    integer,       optional,         intent(in)  :: type
    real(kind=wp), optional,         intent(in)  :: default
    !
    PUSH_SUB(interpolation_init_common)
    this%type=NEAREST
    if(present(type))this%type=type
    this%sim=>sim
    this%nint=nint
    call simulation_get(sim, this%mesh)
    ASSERT(associated(this%mesh))
    call simulation_get(sim, this%domain)
    ASSERT(associated(this%domain))
    if(present(basis))then
      this%basis=>basis
    else
      this%basis=>null()
    end if
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
    end if
    this%default=0.0_wp
    if(present(default))this%default=default
    POP_SUB(interpolation_init_common)
    return
  end subroutine interpolation_init_common

  ! ---------------------------------------------------------
  subroutine interpolation_init_1d(this, sim, vals, basis, type, default)
    type(interpolation_t),               intent(out) :: this
    type(simulation_t),          target, intent(in)  :: sim
    real(kind=wp), dimension(:), target, intent(in)  :: vals
    type(basis_t),     optional, target, intent(in)  :: basis
    integer,           optional,         intent(in)  :: type
    real(kind=wp),     optional,         intent(in)  :: default
    !
    PUSH_SUB(interpolation_init_1d)
    call interpolation_init_common(this, sim, 1, basis, type, default)
    if(type>NONE)&
      call intrp_init(this%intr(1), this%mesh, vals, this%type)
    POP_SUB(interpolation_init_1d)
    return
  end subroutine interpolation_init_1d

  ! ---------------------------------------------------------
  subroutine interpolation_init_2d(this, sim, vals, basis, type, default)
    type(interpolation_t),                 intent(out) :: this
    type(simulation_t),            target, intent(in)  :: sim
    real(kind=wp), dimension(:,:), target, intent(in)  :: vals
    type(basis_t),       optional, target, intent(in)  :: basis
    integer,             optional,         intent(in)  :: type
    real(kind=wp),       optional,         intent(in)  :: default
    !
    integer :: i
    !
    PUSH_SUB(interpolation_init_2d)
    call interpolation_init_common(this, sim, size(vals,dim=2), basis, type, default)
    if(type>NONE)then
      do i = 1, this%nint
        call intrp_init(this%intr(i), this%mesh, vals(:,i), this%type)
      end do
    end if
    POP_SUB(interpolation_init_2d)
    return
  end subroutine interpolation_init_2d

  ! ---------------------------------------------------------
  function interpolation_nearest_index(this, x) result(n)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    !
    integer :: n
    !
    real(kind=wp), dimension(MAX_DIM) :: xp, chi
    integer,       dimension(MAX_DIM) :: ix
    integer                           :: i, dm
    !
    PUSH_SUB(interpolation_nearest_index)
    dm=this%mesh%sb%dim
    xp=(/x(1:dm),(0.0_wp, i=dm+1,MAX_DIM)/)
    call curvilinear_x2chi(this%mesh%sb, this%mesh%cv, xp, chi)
    ix(1:dm)=nint(chi(1:dm)/this%mesh%spacing(1:dm))
    ix(dm+1:MAX_DIM)=0
    n=index_from_coords(this%mesh%idx, dm, ix)
    POP_SUB(interpolation_nearest_index)
    return
  end function interpolation_nearest_index

  ! ---------------------------------------------------------
  elemental function interpolation_get_dimension(this) result(dim)
    type(interpolation_t),  intent(in) :: this
    !
    integer :: dim
    !
    dim=this%nint
    return
  end function interpolation_get_dimension

  ! ---------------------------------------------------------
  function interpolation_in_domain(this, x) result(in)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    !
    logical :: in
    !
    integer :: n
    !
    PUSH_SUB(interpolation_in_domain)
    in=domain_in_domain(this%domain, x)
    if(in)then
      n=interpolation_nearest_index(this, x)
      in=((0<n).and.(n<=size(this%intr(1)%vals)))
    end if
    POP_SUB(interpolation_in_domain)
    return
  end function interpolation_in_domain

  ! ---------------------------------------------------------
  subroutine interpolation_nearest(this, x, val)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    !
    real(kind=wp) :: tol, dlt
    integer       :: i, n, dm
    !
    PUSH_SUB(interpolation_nearest)
    dm=this%mesh%sb%dim
    n=interpolation_nearest_index(this, x)
    tol=CNST(0.51)*sqrt(sum(this%mesh%spacing(1:dm)**2))
    dlt=sqrt(sum((x(1:dm)-this%mesh%x(n,1:dm))**2))
    ASSERT(dlt<tol)
    forall(i=1:this%nint)val(i)=this%intr(i)%vals(n)
    POP_SUB(interpolation_nearest)
    return
  end subroutine interpolation_nearest

  ! ---------------------------------------------------------
  subroutine interpolation_eval_internal(this, x, val, ierr)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    integer :: i
    !
    PUSH_SUB(interpolation_eval_internal)
    if(interpolation_in_domain(this, x))then
      ierr=INTRP_OK
      select case(this%type)
      case(NEAREST)
        call interpolation_nearest(this, x, val)
      case(QSHEP)
        do i = 1, this%nint
          call intrp_eval(this%intr(i), x, val(i))
        end do
      case default
        val=this%default
        ierr=INTRP_NI
      end select
    else
      val=this%default
      ierr=INTRP_OD
    end if
    POP_SUB(interpolation_eval_internal)
    return
  end subroutine interpolation_eval_internal

  ! ---------------------------------------------------------
  subroutine interpolation_eval_1d(this, x, val, ierr)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    real(kind=wp), dimension(size(x)) :: y
    real(kind=wp), dimension(1)       :: tvl
    !
    PUSH_SUB(interpolation_eval_1d)
    if(associated(this%basis))then
      call basis_to_internal(this%basis, x, y)
      call interpolation_eval_internal(this, y, tvl, ierr)
    else
      call interpolation_eval_internal(this, x, tvl, ierr)
    end if
    val=tvl(1)
    POP_SUB(interpolation_eval_1d)
    return
  end subroutine interpolation_eval_1d

  ! ---------------------------------------------------------
  subroutine interpolation_eval_2d(this, x, val, ierr)
    type(interpolation_t),        intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    real(kind=wp), dimension(size(x)) :: y
    !
    PUSH_SUB(interpolation_eval_2d)
    if(associated(this%basis))then
      call basis_to_internal(this%basis, x, y)
      call interpolation_eval_internal(this, y, val, ierr)
    else
      call interpolation_eval_internal(this, x, val, ierr)
    end if
    POP_SUB(interpolation_eval_2d)
    return
  end subroutine interpolation_eval_2d

  ! ---------------------------------------------------------
  subroutine interpolation_copy(this, that)
    type(interpolation_t),         intent(out) :: this
    type(interpolation_t), target, intent(in)  :: that
    !
    integer :: i
    !
    PUSH_SUB(interpolation_copy)
    this%type  = that%type
    this%nint  = that%nint
    this%sim   =>that%sim
    this%mesh  =>that%mesh
    this%domain=>that%domain
    this%basis =>that%basis
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
      do i = 1, this%nint
        call intrp_copy(this%intr(i), that%intr(i))
      end do
    end if
    POP_SUB(interpolation_copy)
    return
  end subroutine interpolation_copy

  ! ---------------------------------------------------------
  subroutine interpolation_end(this)
    type(interpolation_t), intent(inout) :: this
    !
    integer :: i
    !
    PUSH_SUB(interpolation_end)
    this%type=0
    do i = 1, this%nint
      call intrp_end(this%intr(i))
    end do
    if(this%type>NONE)then
      SAFE_DEALLOCATE_P(this%intr)
    end if
    nullify(this%sim, this%mesh, this%domain, this%basis, this%intr)
    this%nint=0
    POP_SUB(interpolation_end)
    return
  end subroutine interpolation_end

end module interpolation_m

!! Local Variables:
!! mode: f90
!! End:
