#include "global.h"

module frozen_basis_m
  use curvilinear_m, only: curvilinear_t, &
    curvilinear_init_from_dump, curvilinear_copy, curvilinear_end
  use global_m
  use messages_m
  use profiling_m
  use simul_box_m, only: simul_box_t, &
    simul_box_init_from_dump, simul_box_copy, simul_box_end
  
  implicit none

  private
  public ::                          &
    frozen_basis_t,                  &
    frozen_basis_init,               &
    frozen_basis_rotate_to_internal, &
    frozen_basis_rotate_to_external, &
    frozen_basis_to_internal,        &
    frozen_basis_to_external,        &
    frozen_basis_copy,               &
    frozen_basis_end

  type frozen_basis_t
    private
    type(simul_box_t)         :: sb
    type(curvilinear_t)       :: cv
    FLOAT, dimension(MAX_DIM) :: r
    FLOAT, dimension(MAX_DIM) :: theta
  end type frozen_basis_t

contains
  
  ! ---------------------------------------------------------
  subroutine frozen_basis_init(this, r, theta, iunit)
    type(frozen_basis_t),      intent(inout) :: this
    FLOAT, dimension(MAX_DIM), intent(in)    :: r
    FLOAT, dimension(MAX_DIM), intent(in)    :: theta
    integer,                   intent(in)    :: iunit
    !
    integer :: gb
    !
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    call simul_box_init_from_dump(this%sb, iunit)
    call curvilinear_init_from_dump(this%cv, iunit)
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    this%r=r
    this%theta=theta
    return
  end subroutine frozen_basis_init

  ! ---------------------------------------------------------
  subroutine frozen_basis_rot(a, b, o)
    FLOAT, intent(inout) :: a, b
    FLOAT, intent(in)    :: o
    !
    FLOAT :: c, s, n, t
    !
    c=cos(o)
    s=sin(o)
    n=sqrt(c**2+s**2)
    c=c/n
    s=s/n
    t=a
    a= c*a + s*b
    b=-s*t + c*b
    return
  end subroutine frozen_basis_rot

  ! ---------------------------------------------------------
  subroutine frozen_basis_rotate_to_internal(this, x)
    type(frozen_basis_t),          intent(in)    :: this
    FLOAT, dimension(this%sb%dim), intent(inout) :: x
    !
    integer :: i, j, k
    !
    k=0
    do i = 1, this%sb%dim-1
      do j = i+1, this%sb%dim
        k=k+1
        call frozen_basis_rot( x(i), x(j), -this%theta(k))
      end do
    end do
    return
  end subroutine frozen_basis_rotate_to_internal

  ! ---------------------------------------------------------
  subroutine frozen_basis_rotate_to_external(this, y)
    type(frozen_basis_t),          intent(in)    :: this
    FLOAT, dimension(this%sb%dim), intent(inout) :: y
    !
    integer :: i, j, k
    !
    k=0
    do i = 1, this%sb%dim-1
      do j = i+1, this%sb%dim
        k=k+1
        call frozen_basis_rot( y(i), y(j), this%theta(k))
      end do
    end do
    return
  end subroutine frozen_basis_rotate_to_external

  ! ---------------------------------------------------------
  subroutine frozen_basis_to_internal(this, x, y)
    type(frozen_basis_t),          intent(in)  :: this
    FLOAT, dimension(this%sb%dim), intent(in)  :: x
    FLOAT, dimension(this%sb%dim), intent(out) :: y
    !
    y=x-this%r
    call frozen_basis_rotate_to_internal(this, y)
    return
  end subroutine frozen_basis_to_internal

  ! ---------------------------------------------------------
  subroutine frozen_basis_to_external(this, y, x)
    type(frozen_basis_t),          intent(in)  :: this
    FLOAT, dimension(this%sb%dim), intent(in)  :: y
    FLOAT, dimension(this%sb%dim), intent(out) :: x
    !
    x=y
    call frozen_basis_rotate_to_external(this, x)
    x=x+this%r
    return
  end subroutine frozen_basis_to_external

  ! ---------------------------------------------------------
  subroutine frozen_basis_copy(this_out, this_in)
    type(frozen_basis_t), intent(inout) :: this_out
    type(frozen_basis_t), intent(in)    :: this_in
    !
    call simul_box_copy(this_out%sb, this_in%sb)
    call curvilinear_copy(this_out%cv, this_in%cv)
    this_out%r=this_in%r
    this_out%theta=this_in%theta
    return
  end subroutine frozen_basis_copy

  ! ---------------------------------------------------------
  subroutine frozen_basis_end(this)
    type(frozen_basis_t), intent(inout) :: this
    !
    call simul_box_end(this%sb)
    call curvilinear_end(this%cv)
    this%r=0.0
    this%theta=0.0
    return
  end subroutine frozen_basis_end

end module frozen_basis_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
