#include "global.h"

module basis_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_len, json_get
  use kinds_m, only: wp, operator(.equal.)

  implicit none

  private
  public ::       &
    operator(==), &
    operator(/=)

  public ::            &
    basis_init,        &
    basis_use,         &
    basis_to_internal, &
    basis_to_external, &
    basis_copy,        &
    basis_end

  integer,       parameter :: ndim = 3
  real(kind=wp), parameter :: eps  = epsilon(1.0_wp)

  type, private :: translation_t
    private
    logical                                  :: do = .false.
    integer                                  :: n  = 0
    real(kind=wp), allocatable, dimension(:) :: r
  end type translation_t

  type, private :: rotation_t
    private
    logical       :: do = .false.
    real(kind=wp) :: c  = 1.0_wp
    real(kind=wp) :: s  = 0.0_wp
  end type rotation_t

  type, private :: rotations_t
    private
    logical                                     :: do = .false.
    integer                                     :: n  = 0
    type(rotation_t), allocatable, dimension(:) :: r
  end type rotations_t

  type, public :: basis_t
    private
    logical             :: do = .false.
    integer             :: dim
    type(translation_t) :: trn
    type(rotations_t)   :: rot
  end type basis_t

  interface operator(==)
    module procedure basis_equal
  end interface operator(==)

  interface operator(/=)
    module procedure basis_not_equal
  end interface operator(/=)

  interface basis_init
    module procedure basis_init_array
    module procedure basis_init_json
  end interface basis_init

contains
  
  ! ---------------------------------------------------------
  subroutine translation_init(this, r)
    type(translation_t),                   intent(out) :: this
    real(kind=wp), dimension(:), optional, intent(in)  :: r
    !
    PUSH_SUB(translation_init)
    this%do=.false.
    this%n=0
    if(present(r))then
      if(any(r>eps))then
        this%do=.true.
        this%n=size(r)
        SAFE_ALLOCATE(this%r(this%n))
        this%r=r
      end if
    end if
    POP_SUB(translation_init)
    return
  end subroutine translation_init

  ! ---------------------------------------------------------
  elemental function translation_equal(this, that) result(eqv)
    type(translation_t), intent(in) :: this
    type(translation_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(this%do.and.that%do)then
      if(this%n==that%n)&
        eqv=all(this%r.equal.that%r)
    end if
    return
  end function translation_equal

  ! ---------------------------------------------------------
  subroutine translation_copy(this_out, this_in)
    type(translation_t), intent(out) :: this_out
    type(translation_t), intent(in)  :: this_in
    !
    PUSH_SUB(translation_copy)
    call translation_init(this_out)
    if(this_in%do)then
      this_out%n=this_in%n
      SAFE_ALLOCATE(this_out%r(this_out%n))
      this_out%r=this_in%r
      this_out%do=.true.
    end if
    POP_SUB(translation_copy)
    return
  end subroutine translation_copy

  ! ---------------------------------------------------------
  subroutine translation_end(this)
    type(translation_t), intent(inout) :: this
    !
    PUSH_SUB(translation_end)
    call translation_init(this)
    if(allocated(this%r))then
      SAFE_DEALLOCATE_A(this%r)
    end if
    POP_SUB(translation_end)
    return
  end subroutine translation_end

  ! ---------------------------------------------------------
  elemental subroutine rotation_init(this, theta)
    type(rotation_t),        intent(out) :: this
    real(kind=wp), optional, intent(in)  :: theta
    !
    real(kind=wp) :: o
    !
    call rotation_end(this)
    if(present(theta))then
      o=modulo(theta,2.0_wp*M_PI)
      if(o<0.0_wp)o=o+2.0_wp*M_PI
      if(o>eps)then
        this%do=.true.
        this%c=cos(theta)
        this%s=sin(theta)
      end if
    end if
    return
  end subroutine rotation_init

  ! ---------------------------------------------------------
  elemental function rotation_equal(this, that) result(eqv)
    type(rotation_t), intent(in) :: this
    type(rotation_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(this%do.and.that%do)then
      eqv=((this%c.equal.that%c).and.(this%s.equal.that%s))
    end if
    return
  end function rotation_equal

  ! ---------------------------------------------------------
  elemental subroutine rotation_copy(this, that)
    type(rotation_t), intent(out) :: this
    type(rotation_t), intent(in)  :: that
    !
    call rotation_end(this)
    if(that%do)then
      this%c=that%c
      this%s=that%s
      this%do=that%do
    end if
    return
  end subroutine rotation_copy

  ! ---------------------------------------------------------
  elemental subroutine rotation_end(this)
    type(rotation_t), intent(inout) :: this
    !
    this%do=.false.
    this%c=1.0_wp
    this%s=0.0_wp
    return
  end subroutine rotation_end

  ! ---------------------------------------------------------
  subroutine rotations_init(this, theta)
    type(rotations_t),                     intent(out) :: this
    real(kind=wp), dimension(:), optional, intent(in)  :: theta
    !
    integer :: i
    !
    PUSH_SUB(rotations_init)
    this%do=.false.
    this%n=0
    if(present(theta))then
      this%n=size(theta)
      SAFE_ALLOCATE(this%r(this%n))
      do i = 1, this%n
        call rotation_init(this%r(i), theta(i))
      end do
      this%do=any(this%r%do)
    end if
    POP_SUB(rotations_init)
    return
  end subroutine rotations_init

  ! ---------------------------------------------------------
  elemental function rotations_equal(this, that) result(eqv)
    type(rotations_t), intent(in) :: this
    type(rotations_t), intent(in) :: that
    !
    logical :: eqv
    !
    integer :: i
    !
    eqv=.false.
    if(this%do.and.that%do)then
      if(this%n==that%n)then
        eqv=.true.
        do i = 1, this%n
          eqv=eqv.and.rotation_equal(this%r(i), that%r(i))
        end do
      end if
    end if
    return
  end function rotations_equal

  ! ---------------------------------------------------------
  subroutine rotations_copy(this_out, this_in)
    type(rotations_t), intent(out) :: this_out
    type(rotations_t), intent(in)  :: this_in
    !
    integer :: i
    !
    PUSH_SUB(rotations_copy)
    call rotations_init(this_out)
    if(this_in%do)then
      this_out%n=this_in%n
      SAFE_ALLOCATE(this_out%r(this_out%n))
      do i = 1, this_out%n
        call rotation_copy(this_out%r(i), this_in%r(i))
      end do
      this_out%do=.true.
    end if
    POP_SUB(rotations_copy)
    return
  end subroutine rotations_copy

  ! ---------------------------------------------------------
  subroutine rotations_end(this)
    type(rotations_t), intent(inout) :: this
    !
    integer :: i
    !
    PUSH_SUB(rotations_end)
    this%do=.false.
    if(allocated(this%r))then
      do i = 1, this%n
        call rotation_end(this%r(i))
      end do
      SAFE_DEALLOCATE_A(this%r)
    end if
    this%n=0
    POP_SUB(rotations_end)
    return
  end subroutine rotations_end

  ! ---------------------------------------------------------
  subroutine basis_init_array(this, r, theta)
    type(basis_t),                         intent(out) :: this
    real(kind=wp), dimension(:), optional, intent(in)  :: r
    real(kind=wp), dimension(:), optional, intent(in)  :: theta
    !
    PUSH_SUB(basis_init_array)
    ASSERT(size(r)==size(theta))
    call translation_init(this%trn, r)
    call rotations_init(this%rot, theta)
    this%do=(this%trn%do.or.this%rot%do)
    POP_SUB(basis_init_array)
    return
  end subroutine basis_init_array

  ! ---------------------------------------------------------
  subroutine basis_init_json(this, config)
    type(basis_t),       intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    real(kind=wp), allocatable, dimension(:) :: r, theta
    !
    integer :: n, ierr
    !
    PUSH_SUB(basis_init_json)
    call json_get(config, "dim", this%dim, ierr)
    if(ierr/=JSON_OK)this%dim=ndim
    n=json_len(config, "r")
    if(n>0)then
      SAFE_ALLOCATE(r(max(n,this%dim)))
      call json_get(config, "r", r, ierr)
      if(ierr/=JSON_OK)r=0.0_wp
      if(n<this%dim)r(n+1:)=0.0_wp
    end if
    n=json_len(config, "theta")
    if(n>0)then
      SAFE_ALLOCATE(theta(max(n,this%dim)))
      call json_get(config, "theta", theta, ierr)
      if(ierr/=JSON_OK)theta=0.0_wp
      if(n<this%dim)theta(n+1:)=0.0_wp
    end if
    if(allocated(r))then
      if(allocated(theta))then
        call basis_init_array(this, r(1:this%dim), theta(1:this%dim))
      else
        call basis_init_array(this, r(1:this%dim))
      end if
    else
      if(allocated(theta))then
        call basis_init_array(this, theta=theta(1:this%dim))
      else
        call basis_init_array(this)
      end if
    end if
    POP_SUB(basis_init_json)
    return
  end subroutine basis_init_json

  ! -----------------------------------------------------
  elemental function basis_equal(this, that) result(eqv)
    type(basis_t), intent(in) :: this
    type(basis_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=.false.
    if(this%do.and.that%do)then
      eqv=translation_equal(this%trn, that%trn)
      if(eqv)eqv=rotations_equal(this%rot, that%rot)
    end if
    return
  end function basis_equal
  
  ! -----------------------------------------------------
  elemental function basis_not_equal(this, that) result(neqv)
    type(basis_t), intent(in) :: this
    type(basis_t), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=.not.basis_equal(this, that)
    return
  end function basis_not_equal

  ! -----------------------------------------------------
  elemental function basis_use(this) result(that)
    type(basis_t), intent(in) :: this
    !
    logical :: that
    !
    that=this%do
    return
  end function basis_use

  ! ---------------------------------------------------------
  elemental subroutine basis_pos_rot(a, b, rot)
    real(kind=wp),    intent(inout) :: a, b
    type(rotation_t), intent(in)    :: rot
    !
    real(kind=wp) :: t
    !
    if(rot%do)then
      t=a
      a= rot%c*a + rot%s*b
      b=-rot%s*t + rot%c*b
    end if
    return
  end subroutine basis_pos_rot

  ! ---------------------------------------------------------
  elemental subroutine basis_neg_rot(a, b, rot)
    real(kind=wp),    intent(inout) :: a, b
    type(rotation_t), intent(in)    :: rot
    !
    real(kind=wp) :: t
    !
    if(rot%do)then
      t=a
      a= rot%c*a - rot%s*b
      b= rot%s*t + rot%c*b
    end if
    return
  end subroutine basis_neg_rot

  ! ---------------------------------------------------------
  pure subroutine basis_to_internal(this, x, y)
    type(basis_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: y
    !
    integer :: i, j, k
    !
    k=0
    if(this%trn%do)then
      y=x-this%trn%r
    else
      y=x
    end if
    do i = 1, this%rot%n-1
      do j = i+1, this%rot%n
        k=k+1
        call basis_neg_rot( y(i), y(j), this%rot%r(k))
      end do
    end do
    return
  end subroutine basis_to_internal

  ! ---------------------------------------------------------
  pure subroutine basis_to_external(this, y, x)
    type(basis_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp), dimension(:), intent(out) :: x
    !
    integer :: i, j, k
    !
    k=0
    x=y
    do i = 1, this%rot%n-1
      do j = i+1, this%rot%n
        k=k+1
        call basis_pos_rot( x(i), x(j), this%rot%r(k))
      end do
    end do
    if(this%trn%do)x=x+this%trn%r
    return
  end subroutine basis_to_external

  ! ---------------------------------------------------------
  subroutine basis_copy(this_out, this_in)
    type(basis_t), intent(out) :: this_out
    type(basis_t), intent(in)  :: this_in
    !
    PUSH_SUB(basis_copy)
    call basis_init_array(this_out)
    if(this_in%do)then
      call translation_copy(this_out%trn, this_in%trn)
      call rotations_copy(this_out%rot, this_in%rot)
      this_out%do=.true.
    end if
    POP_SUB(basis_copy)
    return
  end subroutine basis_copy

  ! ---------------------------------------------------------
  subroutine basis_end(this)
    type(basis_t), intent(inout) :: this
    !
    PUSH_SUB(basis_end)
    if(this%do)then
      this%do=.false.
      call translation_end(this%trn)
      call rotations_end(this%rot)
    end if
    POP_SUB(basis_end)
    return
  end subroutine basis_end

end module basis_m

!! Local Variables:
!! mode: f90
!! End:
