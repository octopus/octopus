#include "global.h"

module basis_m

  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m
  use space_m
  
  implicit none

  private

  public ::  &
    basis_t

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

  real(kind=wp), parameter :: eps = epsilon(1.0_wp)

  type :: translation_t
    private
    logical                                  :: do = .false.
    integer                                  :: n  = 0
    real(kind=wp), allocatable, dimension(:) :: r
  end type translation_t

  type :: rotation_t
    private
    logical       :: do = .false.
    real(kind=wp) :: c  = 1.0_wp
    real(kind=wp) :: s  = 0.0_wp
  end type rotation_t

  type :: rotations_t
    private
    logical                                     :: do = .false.
    integer                                     :: n  = 0
    type(rotation_t), allocatable, dimension(:) :: r
  end type rotations_t

  type :: basis_t
    private
    logical             :: do  = .false.
    integer             :: dim = 0
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

    PUSH_SUB(translation_init)

    if(present(r))then
      if(any(abs(r)>eps))then
        this%do = .true.
        this%n = size(r)
        SAFE_ALLOCATE(this%r(this%n))
        this%r(1:this%n) = r
      end if
    end if

    POP_SUB(translation_init)
  end subroutine translation_init

  ! ---------------------------------------------------------
  elemental function translation_equal(this, that) result(eqv)
    type(translation_t), intent(in) :: this
    type(translation_t), intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(this%do.and.that%do)then
      if(this%n==that%n)&
        eqv = all(this%r.equal.that%r)
    end if

  end function translation_equal

  ! ---------------------------------------------------------
  subroutine translation_copy(this, that)
    type(translation_t), intent(out) :: this
    type(translation_t), intent(in)  :: that

    PUSH_SUB(translation_copy)

    call translation_end(this)
    if(that%do)then
      this%n = that%n
      SAFE_ALLOCATE(this%r(this%n))
      this%r(:) = that%r(:)
      this%do = .true.
    end if

    POP_SUB(translation_copy)
  end subroutine translation_copy

  ! ---------------------------------------------------------
  subroutine translation_end(this)
    type(translation_t), intent(inout) :: this

    PUSH_SUB(translation_end)

    this%do = .false.
    this%n = 0
    if(allocated(this%r))then
      SAFE_DEALLOCATE_A(this%r)
    end if

    POP_SUB(translation_end)
  end subroutine translation_end

  ! ---------------------------------------------------------
  elemental subroutine rotation_init(this, theta)
    type(rotation_t),        intent(out) :: this
    real(kind=wp), optional, intent(in)  :: theta

    real(kind=wp) :: ang

    call rotation_end(this)
    if(present(theta))then
      ang = modulo(theta, 2.0_wp*M_PI)
      if(ang<0.0_wp) ang = ang + 2.0_wp * M_PI
      if(ang>eps)then
        this%do = .true.
        this%c = cos(theta)
        this%s = sin(theta)
      end if
    end if

  end subroutine rotation_init

  ! ---------------------------------------------------------
  elemental function rotation_equal(this, that) result(eqv)
    type(rotation_t), intent(in) :: this
    type(rotation_t), intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(this%do.and.that%do)then
      eqv = ((this%c.equal.that%c).and.(this%s.equal.that%s))
    end if

  end function rotation_equal

  ! ---------------------------------------------------------
  elemental subroutine rotation_copy(this, that)
    type(rotation_t), intent(out) :: this
    type(rotation_t), intent(in)  :: that

    call rotation_end(this)
    if(that%do)then
      this%c = that%c
      this%s = that%s
      this%do = that%do
    end if

  end subroutine rotation_copy

  ! ---------------------------------------------------------
  elemental subroutine rotation_end(this)
    type(rotation_t), intent(inout) :: this

    this%do = .false.
    this%c = 1.0_wp
    this%s = 0.0_wp

  end subroutine rotation_end

  ! ---------------------------------------------------------
  subroutine rotations_init(this, theta)
    type(rotations_t),                     intent(out) :: this
    real(kind=wp), dimension(:), optional, intent(in)  :: theta

    integer :: idx

    PUSH_SUB(rotations_init)

    this%do = .false.
    this%n = 0
    if(present(theta))then
      this%n = size(theta)
      SAFE_ALLOCATE(this%r(this%n))
      do idx = 1, this%n
        call rotation_init(this%r(idx), theta(idx))
      end do
      this%do = any(this%r%do)
    end if

    POP_SUB(rotations_init)
  end subroutine rotations_init

  ! ---------------------------------------------------------
  elemental function rotations_equal(this, that) result(eqv)
    type(rotations_t), intent(in) :: this
    type(rotations_t), intent(in) :: that

    logical :: eqv

    integer :: idx

    eqv = .false.
    if(this%do.and.that%do)then
      if(this%n==that%n)then
        eqv = .true.
        do idx = 1, this%n
          eqv = eqv .and. rotation_equal(this%r(idx), that%r(idx))
        end do
      end if
    end if

  end function rotations_equal

  ! ---------------------------------------------------------
  subroutine rotations_copy(this, that)
    type(rotations_t), intent(out) :: this
    type(rotations_t), intent(in)  :: that

    integer :: idx

    PUSH_SUB(rotations_copy)

    call rotations_end(this)
    if(that%do)then
      this%n = that%n
      SAFE_ALLOCATE(this%r(this%n))
      do idx = 1, this%n
        call rotation_copy(this%r(idx), that%r(idx))
      end do
      this%do = .true.
    end if

    POP_SUB(rotations_copy)
  end subroutine rotations_copy

  ! ---------------------------------------------------------
  subroutine rotations_end(this)
    type(rotations_t), intent(inout) :: this

    integer :: idx

    PUSH_SUB(rotations_end)

    this%do = .false.
    if(allocated(this%r))then
      do idx = 1, this%n
        call rotation_end(this%r(idx))
      end do
      SAFE_DEALLOCATE_A(this%r)
    end if
    this%n = 0

    POP_SUB(rotations_end)
  end subroutine rotations_end

  ! ---------------------------------------------------------
  subroutine basis_init_array(this, space, r, theta)
    type(basis_t),                         intent(out) :: this
    type(space_t),                         intent(in)  :: space
    real(kind=wp), dimension(:), optional, intent(in)  :: r
    real(kind=wp), dimension(:), optional, intent(in)  :: theta

    PUSH_SUB(basis_init_array)

    this%dim = space%dim
    if(present(r))then
      if(this%dim<size(r))then
        call translation_init(this%trn, r(1:this%dim))
      else
        call translation_init(this%trn, r)
      end if
    else
      call translation_init(this%trn)
    end if
    if(present(theta))then
      if(this%dim<size(theta))then
        call rotations_init(this%rot, theta(1:this%dim))
      else
        call rotations_init(this%rot, theta)
      end if
    else
      call rotations_init(this%rot)
    end if
    this%do = ( this%trn%do .or. this%rot%do )

    POP_SUB(basis_init_array)
  end subroutine basis_init_array

  ! ---------------------------------------------------------
  subroutine basis_init_json(this, space, config)
    type(basis_t),       intent(out) :: this
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    real(kind=wp), allocatable, dimension(:) :: r, theta

    integer :: nd, nr, nt, ierr

    PUSH_SUB(basis_init_json)

    call json_get(config, "dim", nd, ierr)
    if(ierr==JSON_OK)then
      ASSERT(nd<=space%dim)
    else
      nd = space%dim
    end if
    nr = json_len(config, "r")
    if(nr>0)then
      SAFE_ALLOCATE(r(nr))
      r = 0.0_wp
      call json_get(config, "r", r, ierr)
      if(ierr/=JSON_OK)then
        nr = 0
        SAFE_DEALLOCATE_A(r)
      end if
    end if
    nt = json_len(config, "theta")
    if(nt>0)then
      SAFE_ALLOCATE(theta(nt))
      theta = 0.0_wp
      call json_get(config, "theta", theta, ierr)
      if(ierr/=JSON_OK)then
        nt = 0
        SAFE_DEALLOCATE_A(theta)
      end if
    end if
    nr = min(nr, nd)
    nt = min(nt, nd)
    if(nr>0)then
      if(nt>0)then
        call basis_init_array(this, space, r=r(1:nr), theta=theta(1:nt))
      else
        call basis_init_array(this, space, r=r(1:nr))
      end if
    else
      if(nt>0)then
        call basis_init_array(this, space, theta=theta(1:nt))
      else
        call basis_init_array(this, space)
      end if
    end if

    POP_SUB(basis_init_json)
  end subroutine basis_init_json

  ! -----------------------------------------------------
  elemental function basis_equal(this, that) result(eqv)
    type(basis_t), intent(in) :: this
    type(basis_t), intent(in) :: that

    logical :: eqv

    eqv = .false.
    if(this%do.and.that%do)then
      eqv = translation_equal(this%trn, that%trn)
      if(eqv) eqv = rotations_equal(this%rot, that%rot)
    end if

  end function basis_equal
  
  ! -----------------------------------------------------
  elemental function basis_not_equal(this, that) result(neqv)
    type(basis_t), intent(in) :: this
    type(basis_t), intent(in) :: that

    logical :: neqv

    neqv = .not. basis_equal(this, that)

  end function basis_not_equal

  ! -----------------------------------------------------
  elemental function basis_use(this) result(that)
    type(basis_t), intent(in) :: this

    logical :: that

    that = this%do

  end function basis_use

  ! ---------------------------------------------------------
  elemental subroutine basis_pos_rot(avl, bvl, rot)
    real(kind=wp),    intent(inout) :: avl
    real(kind=wp),    intent(inout) :: bvl
    type(rotation_t), intent(in)    :: rot

    real(kind=wp) :: tmp

    if(rot%do)then
      tmp =  avl
      avl =  rot%c * avl + rot%s * bvl
      bvl = -rot%s * tmp + rot%c * bvl
    end if

  end subroutine basis_pos_rot

  ! ---------------------------------------------------------
  elemental subroutine basis_neg_rot(avl, bvl, rot)
    real(kind=wp),    intent(inout) :: avl
    real(kind=wp),    intent(inout) :: bvl
    type(rotation_t), intent(in)    :: rot

    real(kind=wp) :: tmp

    if(rot%do)then
      tmp = avl
      avl = rot%c * avl - rot%s * bvl
      bvl = rot%s * tmp + rot%c * bvl
    end if

  end subroutine basis_neg_rot

  ! ---------------------------------------------------------
  pure subroutine basis_to_internal(this, x, y)
    type(basis_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: y

    integer :: idx, jdx, kdx

    kdx = 0
    if(this%trn%do)then
      y = x - this%trn%r
    else
      y = x
    end if
    do idx = 1, this%rot%n-1
      do jdx = idx+1, this%rot%n
        kdx = kdx + 1
        call basis_neg_rot( y(idx), y(jdx), this%rot%r(kdx))
      end do
    end do

  end subroutine basis_to_internal

  ! ---------------------------------------------------------
  pure subroutine basis_to_external(this, y, x)
    type(basis_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp), dimension(:), intent(out) :: x

    integer :: idx, jdx, kdx

    kdx = 0
    x = y
    do idx = 1, this%rot%n-1
      do jdx = idx+1, this%rot%n
        kdx = kdx + 1
        call basis_pos_rot( x(idx), x(jdx), this%rot%r(kdx))
      end do
    end do
    if(this%trn%do) x = x + this%trn%r

  end subroutine basis_to_external

  ! ---------------------------------------------------------
  subroutine basis_copy(this, that)
    type(basis_t), intent(out) :: this
    type(basis_t), intent(in)  :: that

    PUSH_SUB(basis_copy)

    call basis_end(this)
    if(that%do)then
      this%dim = that%dim
      call translation_copy(this%trn, that%trn)
      call rotations_copy(this%rot, that%rot)
      this%do = .true.
    end if

    POP_SUB(basis_copy)
  end subroutine basis_copy

  ! ---------------------------------------------------------
  subroutine basis_end(this)
    type(basis_t), intent(inout) :: this

    PUSH_SUB(basis_end)

    if(this%do)then
      this%do = .false.
      call translation_end(this%trn)
      call rotations_end(this%rot)
    end if
    this%dim = 0

    POP_SUB(basis_end)
  end subroutine basis_end

end module basis_m

!! Local Variables:
!! mode: f90
!! End:
