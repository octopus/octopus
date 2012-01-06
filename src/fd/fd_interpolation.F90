#include "global.h"

module fd_interpolation_m
  use global_m
  use mesh_m, only: mesh_t
  use messages_m
  use profiling_m

  implicit none

  private
  public ::                &
    fd_interpolation_t,    &
    fd_interpolation_init, &
    fd_interpolation,      &
    fd_interpolation_copy, &
    fd_interpolation_end

  integer, public, parameter :: &
    NEAREST = 1

  type fd_intrp_t
    private
    type(mesh_t), pointer               :: mesh
    FLOAT,        pointer, dimension(:) :: vals
  end type fd_intrp_t

  type fd_interpolation_t
    private
    integer                                 :: type
    integer                                 :: nint
    type(mesh_t),     pointer               :: mesh
    type(fd_intrp_t), pointer, dimension(:) :: intr
  end type fd_interpolation_t

  interface fd_interpolation_init
    module procedure fd_intrp_init_1d
    module procedure fd_intrp_init_md
  end interface

  interface fd_interpolation
    module procedure fd_intrp_1d
    module procedure fd_intrp_md
  end interface

contains

  subroutine fd_intrp_init_1d(this, type, mesh, vals)
    type(fd_interpolation_t),           intent(inout) :: this
    integer,                            intent(in)    :: type
    type(mesh_t), target,               intent(in)    :: mesh
    FLOAT,        target, dimension(:), intent(in)    :: vals
    !
    this%type=type
    this%nint=1
    this%mesh=>mesh
    SAFE_ALLOCATE(this%intr(1))
    this%intr(1)%mesh=>mesh
    this%intr(1)%vals=>vals
    return
  end subroutine fd_intrp_init_1d

  subroutine fd_intrp_init_md(this, type, mesh, vals)
    type(fd_interpolation_t),             intent(inout) :: this
    integer,                              intent(in)    :: type
    type(mesh_t), target,                 intent(in)    :: mesh
    FLOAT,        target, dimension(:,:), intent(in)    :: vals
    !
    integer :: i
    !
    this%type=type
    this%nint=size(vals,dim=2)
    this%mesh=>mesh
    SAFE_ALLOCATE(this%intr(this%nint))
    do i = 1, this%nint
      this%intr(i)%mesh=>mesh
      this%intr(i)%vals=>vals(:,i)
    end do
    return
  end subroutine fd_intrp_init_md

  subroutine fd_intrp_nearest(this, x, val)
    type(fd_interpolation_t),         intent(in)  :: this
    FLOAT,      dimension(:),         intent(in)  :: x
    FLOAT,      dimension(this%nint), intent(out) :: val
    !
    integer, dimension(MAX_DIM) :: p
    integer                     :: i, n, m
    FLOAT                       :: tol, dlt, mnv
    !
    n=size(x)
    if(this%mesh%use_curvilinear)then
      m=-1
      mnv=huge(mnv)
      do i = 1, this%mesh%np_part_global
        dlt=sqrt(sum((x-this%mesh%x(i,1:n))**2))
        if(dlt<mnv)then
          mnv=dlt
          m=i
        end if
      end do
    else
      tol=0.5*sqrt(sum(this%mesh%spacing(1:n)**2))+M_EPSILON
      p=nint((/x,(this%mesh%spacing(i), i=n+1, MAX_DIM)/)/this%mesh%spacing)
      m=this%mesh%idx%lxyz_inv(p(1),p(2),p(3))
      dlt=sqrt(sum((x-this%mesh%x(m,1:n))**2))
      ASSERT(dlt<tol)
    end if
    if(m>this%mesh%np_global)then
      val=0.0
    else
      forall(i=1:this%nint) val(i)=this%intr(i)%vals(m)
    end if
    return
  end subroutine fd_intrp_nearest

  subroutine fd_intrp_1d(this, x, val)
    type(fd_interpolation_t), intent(in)  :: this
    FLOAT,      dimension(:), intent(in)  :: x
    FLOAT,                    intent(out) :: val
    !
    FLOAT, dimension(1) :: tvl
    !
    call fd_intrp_md(this, x, tvl)
    val=tvl(1)
    return
  end subroutine fd_intrp_1d

  subroutine fd_intrp_md(this, x, val)
    type(fd_interpolation_t),         intent(in)  :: this
    FLOAT,      dimension(:),         intent(in)  :: x
    FLOAT,      dimension(this%nint), intent(out) :: val
    !
    select case(this%type)
    case(NEAREST)
      call fd_intrp_nearest(this, x, val)
    case default
      write(message(1), '(a,i2)') "Unknown interpolation type: ", this%type
      call messages_fatal(1)
    end select
    return
  end subroutine fd_intrp_md

  subroutine fd_interpolation_copy(this_out, this_in)
    type(fd_interpolation_t),         intent(inout) :: this_out
    type(fd_interpolation_t), target, intent(in)    :: this_in
    !
    integer :: i
    !
    this_out%type=this_in%type
    this_out%nint=this_in%nint
    this_out%mesh=>this_in%mesh
    SAFE_ALLOCATE(this_out%intr(this_in%nint))
    do i = 1, this_in%nint
      this_out%intr(i)%mesh=>this_in%intr(i)%mesh
      this_out%intr(i)%vals=>this_in%intr(i)%vals
    end do
    return
  end subroutine fd_interpolation_copy

  subroutine fd_interpolation_end(this)
    type(fd_interpolation_t), intent(inout) :: this
    !
    integer :: i
    !
    this%type=0
    nullify(this%mesh)
    do i = 1, this%nint
      nullify(this%intr(i)%mesh)
      nullify(this%intr(i)%vals)
    end do
    SAFE_DEALLOCATE_P(this%intr)
    this%nint=0
    return
  end subroutine fd_interpolation_end

end module fd_interpolation_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
