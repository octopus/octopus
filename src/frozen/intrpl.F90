#include "global.h"

module intrpl_m

  use global_m
  use messages_m
  use profiling_m

  use basis_m,       only: basis_t, basis_to_internal
  use curvilinear_m, only: curvilinear_x2chi
  use domain_m,      only: domain_t, domain_in_domain
  use grid_m,        only: grid_t
  use index_m,       only: index_from_coords
  use kinds_m,       only: wp
  use mesh_m,        only: mesh_t
  use qshep_m,       only: qshep_t, qshep_init, qshep_interpolate, qshep_end

  use simulation_m, only: &
    simulation_t,         &
    simulation_get

  implicit none

  private
  public ::                      &
    intrpl_init,          &
    intrpl_inited,        &
    intrpl_get_dimension, &
    intrpl_set,           &
    intrpl_eval,          &
    intrpl_copy,          &
    intrpl_end

  integer, public, parameter :: NONE    = 0
  integer, public, parameter :: NEAREST = 1
  integer, public, parameter :: QSHEP   = 2

  integer, public, parameter :: INTRPL_OK = 0
  integer, public, parameter :: INTRPL_OD = 1
  integer, public, parameter :: INTRPL_NI = 2

  type, private :: intrp_t
    private
    integer                              :: type  = NONE
    type(qshep_t),               pointer :: qshep =>null()
    real(kind=wp), dimension(:), pointer :: vals  =>null()
  end type intrp_t

  type, public :: intrpl_t
    private
    type(simulation_t),              pointer :: sim     =>null()
    type(mesh_t),                    pointer :: mesh    =>null()
    type(domain_t),                  pointer :: domain  =>null()
    type(basis_t),                   pointer :: basis   =>null()
    integer                                  :: type    = NONE
    integer                                  :: nint    = 0
    real(kind=wp)                            :: default = 0.0_wp
    type(intrp_t), dimension(:), allocatable :: intr
  end type intrpl_t

  interface intrpl_init
    module procedure intrpl_init_1d
    module procedure intrpl_init_2d
    module procedure intrpl_init_null
  end interface intrpl_init

  interface intrpl_set
    module procedure intrpl_set_basis
  end interface intrpl_set

  interface intrpl_eval
    module procedure intrpl_eval_1d
    module procedure intrpl_eval_md
  end interface intrpl_eval

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
        call qshep_init(this%qshep, n, vals, mesh%x(1:n,1), mesh%x(1:n,2))
      case(3)
        call qshep_init(this%qshep, n, vals, mesh%x(1:n,1), mesh%x(1:n,2), mesh%x(1:n,3))
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
      call qshep_end(this%qshep)
      SAFE_DEALLOCATE_P(this%qshep)
    end select
    nullify(this%qshep, this%vals)
    this%type=NONE
    POP_SUB(intrp_end)
    return
  end subroutine intrp_end

  ! ---------------------------------------------------------
  subroutine intrpl_init_common(this, sim, type, default)
    type(intrpl_t),             intent(out) :: this
    type(simulation_t), target, intent(in)  :: sim
    integer,          optional, intent(in)  :: type
    real(kind=wp),    optional, intent(in)  :: default
    !
    type(grid_t), pointer :: grid
    !
    PUSH_SUB(intrpl_init_common)
    nullify(grid)
    ASSERT(.not.associated(this%sim))
    ASSERT(.not.associated(this%mesh))
    ASSERT(.not.associated(this%domain))
    this%sim=>sim
    call simulation_get(sim, grid)
    ASSERT(associated(grid))
    this%mesh=>grid%mesh
    nullify(grid)
    call simulation_get(sim, this%domain)
    ASSERT(associated(this%domain))
    nullify(this%basis)
    this%type=NEAREST
    if(present(type))this%type=type
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
    end if
    this%default=0.0_wp
    if(present(default))this%default=default
    POP_SUB(intrpl_init_common)
    return
  end subroutine intrpl_init_common

  ! ---------------------------------------------------------
  subroutine intrpl_init_1d(this, sim, vals, type, default)
    type(intrpl_t),                      intent(out) :: this
    type(simulation_t),          target, intent(in)  :: sim
    real(kind=wp), dimension(:), target, intent(in)  :: vals
    integer,                   optional, intent(in)  :: type
    real(kind=wp),             optional, intent(in)  :: default
    !
    PUSH_SUB(intrpl_init_1d)
    this%nint=1
    call intrpl_init_common(this, sim, type, default)
    if(this%type>NONE)&
      call intrp_init(this%intr(1), this%mesh, vals, this%type)
    POP_SUB(intrpl_init_1d)
    return
  end subroutine intrpl_init_1d

  ! ---------------------------------------------------------
  subroutine intrpl_init_2d(this, sim, vals, type, default)
    type(intrpl_t),                        intent(out) :: this
    type(simulation_t),            target, intent(in)  :: sim
    real(kind=wp), dimension(:,:), target, intent(in)  :: vals
    integer,                     optional, intent(in)  :: type
    real(kind=wp),               optional, intent(in)  :: default
    !
    integer :: i
    !
    PUSH_SUB(intrpl_init_2d)
    this%nint=size(vals,dim=2)
    call intrpl_init_common(this, sim, type, default)
    if(this%type>NONE)then
      do i = 1, this%nint
        call intrp_init(this%intr(i), this%mesh, vals(:,i), this%type)
      end do
    end if
    POP_SUB(intrpl_init_2d)
    return
  end subroutine intrpl_init_2d

  ! ---------------------------------------------------------
  subroutine intrpl_init_null(this, default)
    type(intrpl_t),          intent(out) :: this
    real(kind=wp), optional, intent(in)  :: default
    !
    nullify(this%mesh, this%domain, this%basis)
    this%type=NONE
    this%nint=0
    this%default=0.0_wp
    if(present(default))this%default=default
    return
  end subroutine intrpl_init_null

  ! ---------------------------------------------------------
  elemental function intrpl_inited(this) result(that)
    type(intrpl_t), intent(in) :: this
    !
    logical :: that
    !
    that=(allocated(this%intr))
    return
  end function intrpl_inited

  ! ---------------------------------------------------------
  function intrpl_nearest_index(this, x) result(n)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    !
    integer :: n
    !
    real(kind=wp), dimension(MAX_DIM) :: xp, chi
    integer,       dimension(MAX_DIM) :: ix
    integer                           :: i, dm
    !
    PUSH_SUB(intrpl_nearest_index)
    dm=this%mesh%sb%dim
    xp=(/x(1:dm),(0.0_wp, i=dm+1,MAX_DIM)/)
    call curvilinear_x2chi(this%mesh%sb, this%mesh%cv, xp, chi)
    ix(1:dm)=nint(chi(1:dm)/this%mesh%spacing(1:dm))
    ix(dm+1:MAX_DIM)=0
    n=index_from_coords(this%mesh%idx, ix)
    POP_SUB(intrpl_nearest_index)
    return
  end function intrpl_nearest_index

  ! ---------------------------------------------------------
  elemental function intrpl_get_dimension(this) result(dim)
    type(intrpl_t), intent(in) :: this
    !
    integer :: dim
    !
    dim=this%nint
    return
  end function intrpl_get_dimension

  ! ---------------------------------------------------------
  subroutine intrpl_set_basis(this, that)
    type(intrpl_t),        intent(inout) :: this
    type(basis_t), target, intent(in)    :: that
    !
    PUSH_SUB(intrpl_set_basis)
    this%basis=>that
    POP_SUB(intrpl_set_basis)
    return
  end subroutine intrpl_set_basis

  ! ---------------------------------------------------------
  function intrpl_in_domain(this, x) result(in)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    !
    logical :: in
    !
    integer :: n
    !
    PUSH_SUB(intrpl_in_domain)
    in=.true.
    if(associated(this%domain))&
      in=domain_in_domain(this%domain, x)
    if(in)then
      n=intrpl_nearest_index(this, x)
      in=((0<n).and.(n<=size(this%intr(1)%vals)))
    end if
    POP_SUB(intrpl_in_domain)
    return
  end function intrpl_in_domain

  ! ---------------------------------------------------------
  subroutine intrpl_nearest(this, x, val)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    !
    real(kind=wp) :: tol, dlt
    integer       :: i, n, dm
    !
    PUSH_SUB(intrpl_nearest)
    dm=this%mesh%sb%dim
    n=intrpl_nearest_index(this, x)
    tol=CNST(0.51)*sqrt(sum(this%mesh%spacing(1:dm)**2))
    dlt=sqrt(sum((x(1:dm)-this%mesh%x(n,1:dm))**2))
    ASSERT(dlt<tol)
    forall(i=1:this%nint)val(i)=this%intr(i)%vals(n)
    POP_SUB(intrpl_nearest)
    return
  end subroutine intrpl_nearest

  ! ---------------------------------------------------------
  subroutine intrpl_eval_internal(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    integer :: i
    !
    PUSH_SUB(intrpl_eval_internal)
    if(intrpl_in_domain(this, x))then
      ierr=INTRPL_OK
      select case(this%type)
      case(NEAREST)
        call intrpl_nearest(this, x, val)
      case(QSHEP)
        do i = 1, this%nint
          call intrp_eval(this%intr(i), x, val(i))
        end do
      case default
        val=this%default
        ierr=INTRPL_NI
      end select
    else
      val=this%default
      ierr=INTRPL_OD
    end if
    POP_SUB(intrpl_eval_internal)
    return
  end subroutine intrpl_eval_internal

  ! ---------------------------------------------------------
  subroutine intrpl_eval_1d(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    real(kind=wp), dimension(size(x)) :: y
    real(kind=wp), dimension(1)       :: tvl
    !
    PUSH_SUB(intrpl_eval_1d)
    if(this%type>NONE)then
      if(associated(this%basis))then
        call basis_to_internal(this%basis, x, y)
        call intrpl_eval_internal(this, y, tvl, ierr)
      else
        call intrpl_eval_internal(this, x, tvl, ierr)
      end if
      val=tvl(1)
    else
      val=this%default
      ierr=INTRPL_NI
    end if
    POP_SUB(intrpl_eval_1d)
    return
  end subroutine intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine intrpl_eval_md(this, x, val, ierr)
    type(intrpl_t),               intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    real(kind=wp), dimension(size(x)) :: y
    !
    PUSH_SUB(intrpl_eval_md)
    if(this%type>NONE)then
      if(associated(this%basis))then
        call basis_to_internal(this%basis, x, y)
        call intrpl_eval_internal(this, y, val, ierr)
      else
        call intrpl_eval_internal(this, x, val, ierr)
      end if
    else
      val=this%default
      ierr=INTRPL_NI
    end if
    POP_SUB(intrpl_eval_md)
    return
  end subroutine intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine intrpl_copy(this, that)
    type(intrpl_t),         intent(out) :: this
    type(intrpl_t), target, intent(in)  :: that
    !
    integer :: i
    !
    PUSH_SUB(intrpl_copy)
    this%sim   =>that%sim
    this%mesh  =>that%mesh
    this%domain=>that%domain
    this%basis =>that%basis
    this%type  = that%type
    this%nint  = that%nint
    this%default=that%default
    if(this%type>NONE)then
      SAFE_ALLOCATE(this%intr(this%nint))
      do i = 1, this%nint
        call intrp_copy(this%intr(i), that%intr(i))
      end do
    end if
    POP_SUB(intrpl_copy)
    return
  end subroutine intrpl_copy

  ! ---------------------------------------------------------
  subroutine intrpl_end(this)
    type(intrpl_t), intent(inout) :: this
    !
    integer :: i
    !
    PUSH_SUB(intrpl_end)
    if(this%type>NONE)then
      do i = 1, this%nint
        call intrp_end(this%intr(i))
      end do
      SAFE_DEALLOCATE_A(this%intr)
    end if
    this%default=0.0_wp
    this%nint=0
    this%type=NONE
    nullify(this%basis, this%domain, this%mesh, this%sim)
    POP_SUB(intrpl_end)
    return
  end subroutine intrpl_end

end module intrpl_m

!! Local Variables:
!! mode: f90
!! End:
