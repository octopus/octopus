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
  use units_m

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
            tdf_to_numerical,     &
            tdf,                  &
            tdf_scalar_multiply,  &
            tdf_write,            &
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

  interface tdf_set_numerical
    module procedure tdf_set_numericalr, tdf_set_numericalc, &
                     tdf_set_numericalr1, tdf_set_numericalc1
  end interface

  interface tdf
    module procedure tdfi, tdft
  end interface

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
    ! WARNING: this weird allocataion scheme should be changed in the future:
    ALLOCATE(f%val(niter+1), niter+1)
    f%val = M_z0
    f%dt = dt

    call pop_sub()
  end subroutine tdf_init_numerical


  !------------------------------------------------------------
  subroutine tdf_set_numericalc(f, values)
    type(tdf_t), intent(inout) :: f
    CMPLX,       intent(in) :: values(:)
    call push_sub("tdfunction.tdf_set_numerical") 

    f%val(1:f%niter+1) = values(1:f%niter+1)

    call pop_sub()
  end subroutine tdf_set_numericalc


  !------------------------------------------------------------
  subroutine tdf_set_numericalr(f, values)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in) :: values(:)
    call push_sub("tdfunction.tdf_set_numerical") 

    f%val(1:f%niter+1) = values(1:f%niter+1)

    call pop_sub()
  end subroutine tdf_set_numericalr


  !------------------------------------------------------------
  subroutine tdf_set_numericalc1(f, index, value)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: index
    CMPLX,       intent(in) :: value
    f%val(index) = value
  end subroutine tdf_set_numericalc1


  !------------------------------------------------------------
  subroutine tdf_set_numericalr1(f, index, value)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: index
    FLOAT,       intent(in) :: value
    f%val(index) = value
  end subroutine tdf_set_numericalr1


  !------------------------------------------------------------
  subroutine tdf_to_numerical(f, niter, dt)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: niter
    FLOAT :: dt

    FLOAT :: t
    integer :: j

    if(f%mode .eq. TDF_NUMERICAL) return

    ALLOCATE(f%val(niter+1), niter+1)    

    do j = 1, niter + 1
      t = (j-1)*dt
      f%val(j) = tdf(f, t)
    end do

    f%dt = dt
    f%niter = niter
    f%mode = TDF_NUMERICAL

  end subroutine tdf_to_numerical



  !------------------------------------------------------------
  CMPLX function tdfi(f, i) result(y)
    type(tdf_t), intent(in) :: f
    integer, intent(in)     :: i

    ! Maybe there should be a grid for any kind of function, so
    ! that a meaningul number is produced in any case.
    if(f%mode.ne.TDF_NUMERICAL) then
      y = M_z0
      return
    else
      y = f%val(i)
    end if
    
  end function tdfi


  !------------------------------------------------------------
  CMPLX function tdft(f, t) result(y)
    type(tdf_t), intent(in) :: f
    FLOAT, intent(in)       :: t

    FLOAT :: r
    integer :: il, iu

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

      il = int(t/f%dt)+1; iu = il+1
      if(iu>f%niter+1) then
        y = f%val(il)
      else
        y = f%val(il) + ((f%val(iu)-f%val(il))/f%dt)*(t-(il-1)*f%dt)
      end if

    end select

  end function tdft


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
  ! WARNING: this should be improved to make it more robust.
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
    if(fin%mode .eq. TDF_NUMERICAL) then
      fout%val  = fin%val
    end if

  end subroutine tdf_copy


  !------------------------------------------------------------
  subroutine tdf_scalar_multiply(alpha, f) 
    FLOAT, intent(in) :: alpha
    type(tdf_t), intent(inout) :: f

    select case(f%mode)
    case(TDF_CW, TDF_GAUSSIAN, TDF_COSINOIDAL, TDF_TRAPEZOIDAL)
      f%a0 = alpha*f%a0
    case(TDF_NUMERICAL)
      f%val = alpha*f%val
    case(TDF_FROM_FILE)
      call loct_spline_times(alpha, f%amplitude)
    end select

  end subroutine tdf_scalar_multiply


  !------------------------------------------------------------
  subroutine tdf_write(f, iunit)
    type(tdf_t), intent(in) :: f
    integer, intent(in) :: iunit

    select case(f%mode)
    case(TDF_CW)
      write(iunit,'(3x,a)')          'Mode: continuous wave.'
      write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', f%omega0/units_inp%energy%factor, &
        ' [', trim(units_inp%energy%abbrev), ']'
      write(iunit,'(3x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
    case(TDF_GAUSSIAN)
      write(iunit,'(3x,a)')          'Mode: Gaussian envelope.'
      write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', f%omega0/units_inp%energy%factor, &
        ' [', trim(units_inp%energy%abbrev), ']'
      write(iunit,'(3x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_COSINOIDAL)
      write(iunit,'(3x,a)') 'Mode: cosinoidal envelope.'
      write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', f%omega0/units_inp%energy%factor, &
        ' [', trim(units_inp%energy%abbrev), ']'
      write(iunit,'(3x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_TRAPEZOIDAL)
      write(iunit,'(3x,a)') 'Mode: trapezoidal envelope.'
      write(iunit,'(3x,a,f10.4,3a)') 'Frequency: ', f%omega0/units_inp%energy%factor, &
        ' [', trim(units_inp%energy%abbrev), ']'
      write(iunit,'(3x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(3x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(3x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(3x,a,f10.4,3a)') 'Ramp time: ', f%tau1/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_FROM_FILE)
      write(iunit,'(3x,a)') 'Mode: time-dependent function read from file.'
    case(TDF_NUMERICAL)
      write(iunit,'(3x,a)') 'Mode: time-dependent function stored in a numerical array.'
    end select

  end subroutine tdf_write


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
    nullify(f%val)

  end subroutine tdf_null



end module tdf_m
