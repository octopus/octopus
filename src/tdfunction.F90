!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: lasers.F90 2898 2007-05-08 19:13:46Z acastro $

#include "global.h"

!--------------------------------------------------------------
! This module defines "time-dependent functions", to be used by
! the lasers module, or in the future in order to define time-dependent
! magnetic fields.
!--------------------------------------------------------------
module tdf_m
  use datasets_m
  use global_m
  use io_m
  use messages_m
  use lib_oct_gsl_spline_m

  implicit none

  private
  public :: tdf_t,                &
            tdf_init_cw,          &
            tdf_init_gaussian,    &
            tdf_init_cosinoidal,  &
            tdf_init_trapezoidal, &
            tdf_init_fromfile,    &
            tdf_init_numerical,   &
            tdf_set_numerical,    &
            tdf,                  &
            tdf_end,              &
            assignment(=)


  integer, parameter ::      &
    TDF_CW            =  0,  &
    TDF_GAUSSIAN      =  1,  &
    TDF_COSINOIDAL    =  2,  &
    TDF_TRAPEZOIDAL   =  3,  &
    TDF_FROM_FILE     = 10,  &
    TDF_NUMERICAL     = 99

  type tdf_t
    private
    integer :: mode

    FLOAT   :: t0     ! the time at the maximum of the pulse
    FLOAT   :: tau0   ! the width of the pulse
    FLOAT   :: tau1   ! for the ramped shape, the length of the "ramping" intervals

    FLOAT   :: a0
    FLOAT   :: omega0

    type(loct_spline_t) :: amplitude

    FLOAT   :: dt     ! the time-discretization value.
    integer :: niter
    CMPLX, pointer :: val(:)
  end type tdf_t

  interface assignment (=)
    module procedure tdf_copy
  end interface

  contains

  !------------------------------------------------------------
  subroutine tdf_init_cw(f, a0, omega0)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: a0, omega0

    call push_sub("tdfunction.tdf_init_cw")

    f%mode = TDF_CW
    f%a0 = a0
    f%omega0 = omega0

    call pop_sub()
  end subroutine tdf_init_cw


  !------------------------------------------------------------
  subroutine tdf_init_gaussian(f, a0, omega0, t0, tau0)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: a0, omega0, t0, tau0

    call push_sub("tdfunction.tdf_init_gaussian")

    f%mode = TDF_GAUSSIAN
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0

    call pop_sub()
  end subroutine tdf_init_gaussian


  !------------------------------------------------------------
  subroutine tdf_init_cosinoidal(f, a0, omega0, t0, tau0)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: a0, omega0, t0, tau0

    call push_sub("tdfunction.tdf_init_cosinoidal")

    f%mode = TDF_COSINOIDAL
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0

    call pop_sub()
  end subroutine tdf_init_cosinoidal


  !------------------------------------------------------------
  subroutine tdf_init_trapezoidal(f, a0, omega0, t0, tau0, tau1)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: a0, omega0, t0, tau0, tau1

    call push_sub("tdfunction.tdf_init_trapezoidal")

    f%mode = TDF_TRAPEZOIDAL
    f%a0 = a0
    f%omega0 = omega0
    f%t0 = t0
    f%tau0 = tau0
    f%tau1 = tau1

    call pop_sub()
  end subroutine tdf_init_trapezoidal


  !------------------------------------------------------------
  subroutine tdf_init_fromfile(f, filename, ierr)
    type(tdf_t), intent(inout) :: f
    character(len=*)           :: filename
    integer, intent(out)       :: ierr

    integer :: iunit, lines, j
    FLOAT :: dummy
    FLOAT, allocatable :: t(:), am(:)

    call push_sub("tdfunction.tdf_init_fromfile")

    f%mode = TDF_FROM_FILE
    ierr = 0

    iunit = io_open(trim(filename), action='read', status='old')

    ! count lines in file
    ! WARNING: This job should be done by loct_number_of_lines
    ! WARNING: We should allow for the possibility of a header.
    lines = 0
    do
      read(iunit, *, err=100, end=100) dummy, dummy, dummy
      lines = lines + 1
    end do
100 continue
    rewind(iunit)

    ! allocate and read info
    ALLOCATE( t(lines), lines)
    ALLOCATE(am(lines), lines)
    do j = 1, lines
      read(iunit, *) t(j), am(j)
    end do
    call io_close(iunit)

    call loct_spline_init(f%amplitude)
    call loct_spline_fit(lines, t, am, f%amplitude)

    deallocate(t, am)
    call pop_sub()
  end subroutine tdf_init_fromfile


  !------------------------------------------------------------
  subroutine tdf_init_numerical(f, niter, dt)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: niter
    FLOAT, intent(in)   :: dt

    call push_sub("tdfunction.tdf_init_numerical")

    f%mode = TDF_NUMERICAL
    f%niter = niter
    ALLOCATE(f%val(niter), niter)
    f%val = M_z0

    call pop_sub()
  end subroutine tdf_init_numerical


  !------------------------------------------------------------
  subroutine tdf_set_numerical(f, values)
    type(tdf_t), intent(inout) :: f
    CMPLX,       intent(in) :: values(:)
    call push_sub("tdfunction.tdf_set_numerical") 

    f%val(1:f%niter) = values(1:f%niter)

    call pop_sub()
  end subroutine tdf_set_numerical



  !------------------------------------------------------------
  CMPLX function tdf(f, t) result(y)
    type(tdf_t), intent(in) :: f
    FLOAT, intent(in)       :: t

    FLOAT :: r
    integer :: il, iu

    call push_sub('tdfunction.tdf')

    select case(f%mode)

    case(TDF_CW)

      y = f%a0 * exp(M_zI * f%omega0*t )

    case(TDF_GAUSSIAN)

      r = exp(-(t - f%t0)**2 / (M_TWO*f%tau0**2))
      y = f%a0 * r * exp(M_zI * (f%omega0*t))

    case(TDF_COSINOIDAL)

      if(abs(t - f%t0) <= f%tau0) then
        r = cos( (M_Pi/2)*((t - 2*f%tau0 - f%t0)/f%tau0) )
      end if
      y = f%a0 * r * exp(M_zI * (f%omega0*t))

    case(TDF_TRAPEZOIDAL)

      if(t > f%t0-f%tau0/M_TWO-f%tau1 .and. t <= f%t0-f%tau0/M_TWO) then
        r = (t - (f%t0 - f%tau0/M_TWO - f%tau1))/f%tau1
      elseif(t>f%t0-f%tau0/M_TWO .and. t <=f%t0+f%tau0/M_TWO) then
        r = M_ONE
      elseif(t>f%t0+f%tau0/M_TWO .and. t <=f%t0+f%tau0/M_TWO+f%tau1) then
        r = (f%t0 + f%tau0/M_TWO + f%tau1 - t)/f%tau1
      end if
      y = f%a0 * r * exp(M_zI * (f%omega0*t))

    case(TDF_FROM_FILE)

      y = loct_splint(f%amplitude, t)

    case(TDF_NUMERICAL)

      il = int(t/f%dt); iu = int(t/f%dt) + 1
      y = f%val(il) + ((f%val(iu)-f%val(il))/f%dt)*(t-(il-1)*f%dt)

    end select

    call pop_sub()
  end function tdf


  !------------------------------------------------------------
  subroutine tdf_end(f)
    type(tdf_t), intent(inout) :: f
    call push_sub('tdfunction.tdf_init')

    select case(f%mode)
    case(TDF_FROM_FILE)
      call loct_spline_end(f%amplitude)
    case(TDF_NUMERICAL)
      deallocate(f%val); nullify(f%val)
    end select

    call pop_sub()
  end subroutine tdf_end


  !------------------------------------------------------------
  subroutine tdf_copy(fout, fin)
    type(tdf_t), intent(inout) :: fout
    type(tdf_t), intent(in)  :: fin

    fout%mode   = fin%mode
    fout%t0     = fin%t0  
    fout%tau0   = fin%tau0
    fout%tau1   = fin%tau1
    fout%dt     = fin%dt 
    fout%a0     = fin%a0
    fout%omega0 = fin%omega0 
    fout%niter  = fin%niter
    if(fin%mode .eq. TDF_FROM_FILE) then
      fout%amplitude = fin%amplitude
    end if

  end subroutine tdf_copy


  !------------------------------------------------------------
  ! Nullifies pointers (without deallocation); sets numerical values to scalar variables.
  subroutine tdf_null(f)
    type(tdf_t), intent(inout) :: f

    f%t0     = CNST(0.0)
    f%tau0   = CNST(0.0)
    f%tau1   = CNST(0.0)
    f%dt     = CNST(0.0)
    f%a0     = CNST(0.0)
    f%omega0 = CNST(0.0)
    f%niter  = 0
    f%mode = -1

  end subroutine tdf_null



end module tdf_m
