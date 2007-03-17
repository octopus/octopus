module qshepmod_m

  implicit none

  private
  public :: qshep_t, &
            init_qshep, &
            qshep_interpolate, &
            kill_qshep

  type qshepr_t
   private
   integer :: npoints, nq, nw, nr, dim
   integer, pointer :: lcell(:, :, :), lnext(:)
   real(8), pointer :: rsq(:), a(:, :)
   real(8) :: xmin, ymin, dx, dy, rmax, xyzmin(3), xyzdel(3)

   real(8), pointer :: x(:), y(:), z(:)
  end type qshepr_t

  type qshep_t
   private
   integer :: kind ! 0 for real functions (im is not used);  for complex ones.
   type(qshepr_t) :: re, im
  end type qshep_t

  interface init_qshep
    module procedure dinit_qshep, zinit_qshep
  end interface

  interface qshep_interpolate
    module procedure dqshep_interpolate, zqshep_interpolate
  end interface qshep_interpolate

  contains

  subroutine dinit_qshep(interp, npoints, f, x, y, z)
    type(qshep_t), intent(out) :: interp
    integer, intent(in) :: npoints
    real(8), intent(in) :: f(:)
    real(8), target :: x(:), y(:)
    real(8), target, optional :: z(:)

    interp%kind = 0
    if(present(z)) then
      call init_qshepr(interp%re, npoints, f, x, y, z)
    else
      call init_qshepr(interp%re, npoints, f, x, y)
    end if

  end subroutine dinit_qshep

  subroutine zinit_qshep(interp, npoints, f, x, y, z)
    type(qshep_t), intent(out) :: interp
    integer, intent(in) :: npoints
    complex(8), intent(in) :: f(:)
    real(8), target :: x(:), y(:)
    real(8), target, optional :: z(:)

    interp%kind = 1
    if(present(z)) then
      call init_qshepr(interp%re, npoints, real(f), x, y, z)
      call init_qshepr(interp%im, npoints, aimag(f), x, y, z)
    else
      call init_qshepr(interp%re, npoints, real(f), x, y)
      call init_qshepr(interp%im, npoints, aimag(f), x, y)
    end if

  end subroutine zinit_qshep

  subroutine init_qshepr(interp, npoints, f, x, y, z)
    type(qshepr_t), intent(out) :: interp
    integer, intent(in) :: npoints
    real(8), intent(in) :: f(:)
    real(8), target :: x(:), y(:)
    real(8), target, optional :: z(:)

    integer :: ier

    interp%npoints = npoints

    if(present(z)) then
      interp%dim = 3
    else
      interp%dim = 2
    end if

    select case(interp%dim)
    case(2)
      interp%nq = 13
      interp%nw = 19
      interp%nr = nint(sqrt(interp%npoints/3.0_8))
      allocate(interp%lcell(interp%nr, interp%nr, 1))
      allocate(interp%a(5, npoints))
    case(3) 
      interp%nq = 17 ! This is the recommended value in qshep3d.f90
      interp%nw = 16 ! The recommended value in qshep3d.f90 is 32, but this speeds up things.
      interp%nr = nint((interp%npoints/3.0_8)**(1.0_8/3.0_8))
      allocate(interp%lcell(interp%nr, interp%nr, interp%nr))
      allocate(interp%a(9, interp%npoints))
    end select

    allocate(interp%lnext(npoints))
    allocate(interp%rsq(npoints))

    select case(interp%dim)
    case(2)
      call qshep2 ( npoints, x, y, f, interp%nq, interp%nw, interp%nr, interp%lcell(:, :, 1), &
                    interp%lnext, interp%xmin, interp%ymin, interp%dx, interp%dy, &
                    interp%rmax, interp%rsq, interp%a, ier )
      interp%x => x
      interp%y => y
    case(3)
      call qshep3 ( npoints, x, y, z, f, interp%nq, interp%nw, interp%nr, interp%lcell, &
                    interp%lnext, interp%xyzmin, interp%xyzdel, interp%rmax, interp%rsq, &
                    interp%a, ier )
      interp%x => x
      interp%y => y
      interp%z => z
    end select

  end subroutine init_qshepr


  real(8) function qshep_interpolater(interp, f, px, py, pz) result(v)
    type(qshepr_t), intent(in) :: interp
    real(8), intent(in) :: f(:)
    real(8), intent(in) :: px, py
    real(8), optional, intent(in) :: pz

    real(8), external :: qs2val, qs3val

    select case(interp%dim)
    case(2)
      v = qs2val ( px, py, interp%npoints, interp%x, interp%y, &
                   f, interp%nr, interp%lcell(:, :, 1), interp%lnext, interp%xmin, &
                   interp%ymin, interp%dx, interp%dy, interp%rmax, interp%rsq, interp%a )
    case(3)
      v = qs3val ( px, py, pz, interp%npoints, interp%x, interp%y, interp%z, &
                   f, interp%nr, interp%lcell, interp%lnext, interp%xyzmin, &
                   interp%xyzdel, interp%rmax, interp%rsq, interp%a )
    end select

  end function qshep_interpolater

  real(8) function dqshep_interpolate(interp, f, px, py, pz) result(v)
    type(qshep_t), intent(in) :: interp
    real(8), intent(in) :: f(:)
    real(8), intent(in) :: px, py
    real(8), optional, intent(in) :: pz

    if(present(pz)) then
      v = qshep_interpolater(interp%re, f, px, py, pz)
    else
      v = qshep_interpolater(interp%re, f, px, py)
    end if

  end function dqshep_interpolate

  complex(8) function zqshep_interpolate(interp, f, px, py, pz) result(v)
    type(qshep_t), intent(in) :: interp
    complex(8), intent(in) :: f(:)
    real(8), intent(in) :: px, py
    real(8), optional, intent(in) :: pz

    if(present(pz)) then
      v = cmplx(qshep_interpolater(interp%re, real(f), px, py, pz), &
               qshep_interpolater(interp%im, aimag(f), px, py, pz), 8)
    else
      v = cmplx(qshep_interpolater(interp%re, real(f), px, py), &
               qshep_interpolater(interp%im, aimag(f), px, py), 8)
    end if

  end function zqshep_interpolate

  subroutine kill_qshep(interp)
    type(qshep_t), intent(inout) :: interp

    call kill_qshepr(interp%re)
    if(interp%kind == 1) call kill_qshepr(interp%im)

  end subroutine kill_qshep

  subroutine kill_qshepr(interp)
    type(qshepr_t), intent(inout) :: interp

    if(associated(interp%lcell)) then
      nullify(interp%lcell, interp%lnext, interp%rsq, interp%a, interp%x, interp%y)
      if(interp%dim .eq. 3) nullify(interp%z)
    end if

  end subroutine kill_qshepr

end module qshepmod_m
