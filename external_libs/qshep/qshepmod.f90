module qshepmod_m

  implicit none

  private
  public :: qshep_t, &
            init_qshep, &
            qshep_interpolate, &
            kill_qshep

  type qshep_t
   private
   integer :: npoints, nq, nw, nr, dim
   integer, pointer :: lcell(:, :, :), lnext(:)
   real(8), pointer :: rsq(:), a(:, :)
   real(8) :: xmin, ymin, dx, dy, rmax, xyzmin(3), xyzdel(3)

   real(8), pointer :: x(:), y(:), z(:)
  end type qshep_t

  contains

  subroutine init_qshep(interp, npoints, f, x, y, z)
    type(qshep_t), intent(out) :: interp
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

  end subroutine init_qshep

  real(8) function qshep_interpolate(interp, f, px, py, pz) result(v)
    type(qshep_t), intent(in) :: interp
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

  end function qshep_interpolate


  subroutine kill_qshep(interp)
    type(qshep_t), intent(inout) :: interp

    if(associated(interp%lcell)) then
      nullify(interp%lcell, interp%lnext, interp%rsq, interp%a, interp%x, interp%y)
      if(interp%dim .eq. 3) nullify(interp%z)
    end if

  end subroutine kill_qshep

end module qshepmod_m
