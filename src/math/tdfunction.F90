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
  use loct_parser_m
  use splines_m
  use units_m
  use fft_m

  implicit none

  private
  public :: tdf_t,                       &
            tdf_init,                    &
            tdf_init_cw,                 &
            tdf_init_gaussian,           &
            tdf_init_cosinoidal,         &
            tdf_init_trapezoidal,        &
            tdf_init_fromfile,           &
            tdf_init_fromexpr,           &
            tdf_init_numerical,          &
            tdf_set_numerical,           &
            tdf_to_numerical,            &
            tdf_numerical_keep_real,     &
            tdf,                         &
            tdf_dot_product,             &
            tdf_scalar_multiply,         &
            tdf_cosine_multiply,         &
            tdf_cosine_divide,           &
            tdf_fft_forward,             &
            tdf_fft_backward,            &
            tdf_sineseries_to_numerical, &
            tdf_numerical_to_sineseries, &
            tdf_fourier_grid,            &
            tdf_write,                   &
            tdf_niter,                   &
            tdf_nfreqs,                  &
            tdf_dt,                      &
            tdf_copy,                    &
            tdf_read,                    &
            tdf_is_empty,                &
            tdf_end


  integer, parameter ::      &
    TDF_EMPTY         =  10001,  &
    TDF_CW            =  10002,  &
    TDF_GAUSSIAN      =  10003,  &
    TDF_COSINOIDAL    =  10004,  &
    TDF_TRAPEZOIDAL   =  10005,  &
    TDF_FROM_FILE     =  10006,  &
    TDF_NUMERICAL     =  10007,  &
    TDF_FROM_EXPR     =  10008,  &
    TDF_SINE_SERIES   =  10009

  type tdf_t
    private
    integer :: mode        = TDF_EMPTY
    FLOAT   :: t0          = M_ZERO  ! the time at the maximum of the pulse
    FLOAT   :: tau0        = M_ZERO  ! the width of the pulse
    FLOAT   :: tau1        = M_ZERO  ! for the ramped shape, the length of the "ramping" intervals
    FLOAT   :: a0          = M_ZERO
    FLOAT   :: omega0      = M_ZERO
    FLOAT   :: dt          = M_ZERO ! the time-discretization value.
    FLOAT   :: init_time   = M_ZERO
    FLOAT   :: final_time  = M_ZERO
    integer :: niter       = 0
    integer :: nfreqs      = 0

    type(spline_t)         :: amplitude
    character(len=200)     :: expression
    CMPLX, pointer :: val(:)    => NULL()
    FLOAT, pointer :: coeffs(:) => NULL()
  end type tdf_t

  interface tdf_set_numerical
    module procedure tdf_set_numericalr, tdf_set_numericalc, &
                     tdf_set_numericalr1, tdf_set_numericalc1
  end interface

  interface tdf
    module procedure tdfi, tdft
  end interface

  contains



  !------------------------------------------------------------
  subroutine tdf_read(f, function_name, ierr)
    type(tdf_t),      intent(inout) :: f
    character(len=*), intent(in)    :: function_name
    integer,          intent(out)   :: ierr

    type(block_t) :: blk
    integer :: nrows, i, function_type
    character(len=100) :: row_name, filename, function_expression
    FLOAT :: a0, tau0, t0, tau1

    !%Variable TDExternalFields
    !%Type block
    !%Section Time Dependent
    !%Description
    !%
    !%    (1) tdf_cw
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; tdf_cw | amplitude 
    !% <br>%</tt>
    !%
    !% The function is just a constant of value "amplitude".
    !%
    !% <math> f(t) = amplitude
    !%
    !%    (A.2) tdf_gaussian
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; tdf_gaussian | amplitude | tau0 | t0
    !% <br>%</tt>
    !% 
    !% The function is a gaussian:
    !%
    !% <math> f(t) = F_0 exp( - (t-t_0)/(2\tau_0^2) ) </math>
    !%
    !% <math>F_0</math> = amplitude.
    !%
    !%    (A.3) tdf_cosinoidal
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; tdf_cosinoidal | amplitude | tau0 | t0
    !% <br>%</tt>
    !%
    !% <math> f(t) =  F_0 cos( \pi/2 \frac{t-2\tau_0-t_0}{\tau0} )  </math>
    !%
    !% If <math> | t - t_0 | > \tau_0 </math>, then <math> f(t) = 0 </math>.
    !%
    !%    (A.4) tdf_trapezoidal
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; envelope_trapezoidal | amplitude | tau0 | t0 | tau1
    !% <br>%</tt>
    !%
    !% The function ramps linearly during <math>tau_1</math> time units, stays constant for
    !% <math>tau_0</math> timu units, and the decays to zero linearly again for <math>tau_1</math>
    !% time units.
    !%
    !%    (A.5) tdf_from_file
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; envelope_from_file | "filename"
    !% <br>%</tt>
    !%
    !% The temporal shape of the function is contained in a file called "filename". This file
    !% should contain three columns: first column is time, second and third column are the
    !% real part and the imaginary part of the temporal function f(t).
    !%
    !%    (A.6) envelope_from_expr
    !%
    !% <tt>%TDExternalFields
    !% <br>&nbsp;&nbsp; envelope_from_expr | "expression"
    !% <br>%</tt>
    !%
    !% The temporal shape of the field is given as an expression (e.g., "cos(2.0*t)". The 
    !% letter "t" means time, obviously. The expression is used to construct the function f
    !% that defines the field:
    !%
    !%Option tdf_cw 10002
    !% Explained above.
    !%Option tdf_gaussian 10003
    !% Explained above.
    !%Option tdf_cosinoidal 10004
    !% Explained above.
    !%Option tdf_trapezoidal 10005
    !% Explained above.
    !%Option tdf_from_file 10006
    !% Explained above.
    !%Option tdf_from_expr 10008
    !% Explained above.
    !%End
    if(loct_parse_block(check_inp('TDFunctions'), blk) .ne. 0) then
      ierr = -1
      return
    end if

    nrows = loct_parse_block_n(blk)
    row_loop: do i = 1, nrows
      call loct_parse_block_string(blk, i-1, 0, row_name)
      if(trim(row_name).eq.trim(function_name)) then

        call loct_parse_block_int  (blk, i-1, 1, function_type)

        select case(function_type)
          case(TDF_CW)
            call loct_parse_block_float(blk, i-1, 2, a0)
          case(TDF_GAUSSIAN)
            call loct_parse_block_float(blk, i-1, 2, a0)
            call loct_parse_block_float(blk, i-1, 3, tau0)
          case(TDF_COSINOIDAL)
            call loct_parse_block_float(blk, i-1, 2, a0)
            call loct_parse_block_float(blk, i-1, 3, tau0)
            call loct_parse_block_float(blk, i-1, 4, t0)
          case(TDF_TRAPEZOIDAL)
            call loct_parse_block_float(blk, i-1, 2, a0)
            call loct_parse_block_float(blk, i-1, 3, tau0)
            call loct_parse_block_float(blk, i-1, 4, t0)
            call loct_parse_block_float(blk, i-1, 5, tau1)
          case(TDF_FROM_FILE)
            call loct_parse_block_string(blk, i-1, 2, filename)
          case(TDF_FROM_EXPR)
            call loct_parse_block_string(blk, i-1, 2, function_expression)
          case default
            ierr = -2
            call loct_parse_block_end(blk)
            return
        end select

        a0     = a0 * units_inp%energy%factor / units_inp%length%factor
        tau0 = tau0 * units_inp%time%factor
        t0   = t0   * units_inp%time%factor
        tau1 = tau1 * units_inp%time%factor

        select case(function_type)
        case(TDF_CW)
          call tdf_init_cw(f, a0, M_ZERO)
        case(TDF_GAUSSIAN)
          call tdf_init_gaussian(f, a0, M_ZERO, t0, tau0)
        case(TDF_COSINOIDAL)
          call tdf_init_cosinoidal(f, a0, M_ZERO, t0, tau0)
        case(TDF_TRAPEZOIDAL)
          call tdf_init_trapezoidal(f, a0, M_ZERO, t0, tau0, tau1)
        case(TDF_FROM_FILE)
          call tdf_init_fromfile(f, trim(filename), ierr)
        case(TDF_FROM_EXPR)
          call tdf_init_fromexpr(f, trim(function_expression))
        end select
        ierr = 0
        exit row_loop
      end if
    end do row_loop

    call loct_parse_block_end(blk)
  end subroutine tdf_read
  !------------------------------------------------------------


  !------------------------------------------------------------
  integer function tdf_niter(f)
    type(tdf_t), intent(in) :: f
    tdf_niter = f%niter
  end function tdf_niter
  !------------------------------------------------------------


  !------------------------------------------------------------
  integer function tdf_nfreqs(f)
    type(tdf_t), intent(in) :: f
    tdf_nfreqs = f%nfreqs
  end function tdf_nfreqs
  !------------------------------------------------------------


  !------------------------------------------------------------
  FLOAT function tdf_dt(f)
    type(tdf_t), intent(in) :: f
    tdf_dt = f%dt
  end function tdf_dt
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init(f)
    type(tdf_t), intent(inout) :: f
    f%mode = TDF_EMPTY
    f%niter = 0
    f%dt = M_ZERO
    nullify(f%val)
    nullify(f%coeffs)
  end subroutine tdf_init
  !------------------------------------------------------------


  !------------------------------------------------------------
  logical function tdf_is_empty(f)
    type(tdf_t), intent(in) :: f
    tdf_is_empty = (f%mode == TDF_EMPTY)
  end function tdf_is_empty
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_cw(f, a0, omega0)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: a0, omega0

    call push_sub("tdfunction.tdf_init_cw")

    f%mode = TDF_CW
    f%a0 = a0
    f%omega0 = omega0
    nullify(f%val)
    nullify(f%coeffs)

    call pop_sub()
  end subroutine tdf_init_cw
  !------------------------------------------------------------


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
    nullify(f%val)
    nullify(f%coeffs)

    call pop_sub()
  end subroutine tdf_init_gaussian
  !------------------------------------------------------------


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
    nullify(f%val)
    nullify(f%coeffs)

    call pop_sub()
  end subroutine tdf_init_cosinoidal
  !------------------------------------------------------------


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
    nullify(f%val)
    nullify(f%coeffs)

    call pop_sub()
  end subroutine tdf_init_trapezoidal
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_fromexpr(f, expression)
    type(tdf_t), intent(inout)   :: f
    character(len=*), intent(in) :: expression

    call push_sub('tdfunction.tdf_init_fromexpr')

    f%mode = TDF_FROM_EXPR
    f%expression = trim(expression)
    nullify(f%val)
    nullify(f%coeffs)
    
    call pop_sub()
  end subroutine tdf_init_fromexpr
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_fromfile(f, filename, ierr)
    type(tdf_t),      intent(inout) :: f
    character(len=*), intent(in)    :: filename
    integer,          intent(out)   :: ierr

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

    f%init_time  = t(1)
    f%final_time = t(lines)

    call spline_init(f%amplitude)
    call spline_fit(lines, t, am, f%amplitude)

    nullify(f%val)
    nullify(f%coeffs)
    deallocate(t, am)
    call pop_sub()
  end subroutine tdf_init_fromfile
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_init_numerical(f, niter, dt, initval, omegamax)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: niter
    FLOAT, intent(in)   :: dt
    CMPLX, intent(in), optional :: initval
    FLOAT, intent(in), optional :: omegamax

    FLOAT :: bigt

    call push_sub("tdfunction.tdf_init_numerical")

    f%mode = TDF_NUMERICAL
    f%niter = niter
    ! WARNING: this weird allocation scheme should be changed in the future:
    ALLOCATE(f%val(niter+1), niter+1)
    if(present(initval)) then
      f%val = initval
    else
      f%val = M_z0
    end if
    f%dt = dt

    f%init_time = M_ZERO
    f%final_time = f%dt * f%niter

    f%nfreqs = f%niter
    if(present(omegamax)) then
      if(omegamax > M_ZERO) then
        bigt = f%final_time - f%init_time
        f%nfreqs = int(bigt * omegamax / M_PI) + 1
      end if
    end if

    nullify(f%coeffs)
    call pop_sub()
  end subroutine tdf_init_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_fourier_grid(f, wgrid)
    type(tdf_t), intent(in) :: f
    FLOAT, intent(inout)    :: wgrid(:)
    integer :: i, steps
    FLOAT   :: df

    wgrid = M_ZERO
    steps = f%niter + 1

    df = M_ONE/(real(steps, REAL_PRECISION) * f%dt)
    do i = 1, int((steps)/2)
       wgrid(i)   = real(i-1, REAL_PRECISION) * df * M_TWO * M_PI
    enddo
    do i = int((steps)/2), f%niter-1 
       wgrid(i+1) = real(i-steps, REAL_PRECISION) * df * M_TWO * M_PI
    enddo   

  end subroutine tdf_fourier_grid
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_fft_forward(f)
    type(tdf_t), intent(inout) :: f
    integer :: steps, n(3)
    type(fft_t)    :: fft_handler
    CMPLX, allocatable :: tmp(:, :, :), tmp2(:, :, :)

    steps = f%niter
    ALLOCATE(tmp(steps+1, 1, 1), steps+1)
    ALLOCATE(tmp2(steps+1, 1, 1), steps+1)

    n(1:3) = (/ steps+1, 1, 1 /)
    call fft_init(n, fft_complex, fft_handler, optimize = .false.)
    tmp2(:, 1, 1) = f%val(:)
    call zfft_forward(fft_handler, tmp2, tmp)
    f%val(:) = tmp(:, 1, 1)
    deallocate(tmp, tmp2)
    call fft_end(fft_handler)

  end subroutine tdf_fft_forward
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_fft_backward(f)
    type(tdf_t), intent(inout) :: f
    integer :: steps, n(3)
    type(fft_t)    :: fft_handler
    CMPLX, allocatable :: tmp(:, :, :), tmp2(:, :, :)

    steps = f%niter
    ALLOCATE(tmp(steps+1, 1, 1), steps+1)
    ALLOCATE(tmp2(steps+1, 1, 1), steps+1)

    n(1:3) = (/ steps+1, 1, 1 /)
    call fft_init(n, fft_complex, fft_handler, optimize = .false.)
    tmp2(:, 1, 1) = f%val(:)
    call zfft_backward(fft_handler, tmp2, tmp)
    f%val(:) = tmp(:, 1, 1)
    deallocate(tmp, tmp2)
    call fft_end(fft_handler)

  end subroutine tdf_fft_backward
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalc(f, values)
    type(tdf_t), intent(inout) :: f
    CMPLX,       intent(in) :: values(:)
    call push_sub("tdfunction.tdf_set_numerical") 

    select case(f%mode)
    case(TDF_NUMERICAL)
      f%val(1:f%niter+1) = values(1:f%niter+1)
    case(TDF_SINE_SERIES)
      f%coeffs(1:f%nfreqs) = values(1:f%nfreqs)
    end select

    call pop_sub()
  end subroutine tdf_set_numericalc
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalr(f, values)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in) :: values(:)
    call push_sub("tdfunction.tdf_set_numerical") 

    select case(f%mode)
    case(TDF_NUMERICAL)
      f%val(1:f%niter+1) = values(1:f%niter+1)
    case(TDF_SINE_SERIES)
      f%coeffs(1:f%nfreqs) = values(1:f%nfreqs)
    end select

    call pop_sub()
  end subroutine tdf_set_numericalr
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalc1(f, index, value)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: index
    CMPLX,       intent(in) :: value
    select case(f%mode)
    case(TDF_NUMERICAL)
      f%val(index) = value
    case(TDF_SINE_SERIES)
      f%coeffs(index) = value
    end select
  end subroutine tdf_set_numericalc1
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_set_numericalr1(f, index, value)
    type(tdf_t), intent(inout) :: f
    integer, intent(in) :: index
    FLOAT,       intent(in) :: value
    select case(f%mode)
    case(TDF_NUMERICAL)
      f%val(index) = value
    case(TDF_SINE_SERIES)
      f%coeffs(index) = value
    end select
  end subroutine tdf_set_numericalr1
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_numerical_keep_real(f)
    type(tdf_t), intent(inout) :: f

    if(f%mode .ne. TDF_NUMERICAL) return

    f%val(:) = real(f%val(:), REAL_PRECISION)

  end subroutine tdf_numerical_keep_real
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_to_numerical(f, t0, niter, dt)
    type(tdf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: t0
    integer,     intent(in)    :: niter
    FLOAT,       intent(in)    :: dt

    FLOAT :: t
    integer :: j

    if(f%mode .eq. TDF_NUMERICAL) return

    ALLOCATE(f%val(niter+1), niter+1)

    f%init_time = t0
    do j = 1, niter + 1
      t = t0 + (j-1)*dt
      f%val(j) = tdf(f, t)
    end do
    f%final_time = t

    f%dt = dt
    f%niter = niter
    f%mode = TDF_NUMERICAL
    if (f%mode .eq. TDF_FROM_FILE) call spline_end(f%amplitude)

  end subroutine tdf_to_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_numerical_to_sineseries(f, omegamax)
    type(tdf_t), intent(inout) :: f
    FLOAT, intent(in) :: omegamax

    integer :: j, k, nfreqs
    FLOAT :: bigt, t, omega

    ASSERT(f%mode .eq. TDF_NUMERICAL)

    bigt = f%final_time - f%init_time

    if(omegamax > M_ZERO) then
      nfreqs = int(bigt * omegamax / M_PI) + 1
    else
      nfreqs = f%niter
    end if

    f%nfreqs = nfreqs
    ALLOCATE(f%coeffs(nfreqs), nfreqs)

    do j = 1, nfreqs
      omega = (M_PI/bigt) * j
      f%coeffs(j) = M_ZERO
      do k = 2, f%niter
        t = (k-1)*f%dt
        f%coeffs(j) = f%coeffs(j) + sin(omega*t) * tdf(f, k)
      end do
      f%coeffs(j) = sqrt(M_TWO/bigt) * f%coeffs(j) * f%dt
    end do

    deallocate(f%val); nullify(f%val)
    f%mode = TDF_SINE_SERIES

  end subroutine tdf_numerical_to_sineseries
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_sineseries_to_numerical(f)
    type(tdf_t), intent(inout) :: f

    FLOAT :: bigt, t, omega
    integer :: j, k

    ASSERT(f%mode .eq. TDF_SINE_SERIES)

    ALLOCATE(f%val(f%niter+1), f%niter+1)

    bigt = f%final_time - f%init_time

    do k = 1, f%niter + 1
      t = (k-1)*f%dt
      f%val(k) = M_ZERO
      do j = 1, f%nfreqs
        omega = (M_PI/bigt) * j
        f%val(k) = f%val(k) + sin(omega*t) * f%coeffs(j)
      end do
      f%val(k) = f%val(k) * sqrt(M_TWO/bigt)
    end do

    deallocate(f%coeffs); nullify(f%coeffs)
    f%mode = TDF_NUMERICAL

  end subroutine tdf_sineseries_to_numerical
  !------------------------------------------------------------


  !------------------------------------------------------------
  CMPLX function tdfi(f, i) result(y)
    type(tdf_t), intent(in) :: f
    integer, intent(in)     :: i

    ! Maybe there should be a grid for any kind of function, so
    ! that a meaningul number is produced in any case.
    y = M_z0
    select case(f%mode)
    case(TDF_NUMERICAL)
      y = f%val(i)
    case(TDF_SINE_SERIES)
      y = f%coeffs(i)
    end select
    
  end function tdfi
  !------------------------------------------------------------


  !------------------------------------------------------------
  CMPLX function tdft(f, t) result(y)
    type(tdf_t), intent(in) :: f
    FLOAT, intent(in)       :: t

    FLOAT :: r, fre, fim
    integer :: il, iu

    select case(f%mode)

    case(TDF_CW)

      y = f%a0 * exp(M_zI * f%omega0*t )

    case(TDF_GAUSSIAN)

      r = exp(-(t - f%t0)**2 / (M_TWO*f%tau0**2))
      y = f%a0 * r * exp(M_zI * (f%omega0*t))

    case(TDF_COSINOIDAL)

      r = M_ZERO
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
      else
        r = M_ZERO
      end if
      y = f%a0 * r * exp(M_zI * (f%omega0*t))
    case(TDF_FROM_FILE)

      if( t >= f%init_time .and. t <= f%final_time) then
        y = spline_eval(f%amplitude, t)
      else
        y = M_ZERO
      end if

    case(TDF_NUMERICAL)

      il = int(t/f%dt)+1; iu = il+1
      if(iu>f%niter+1) then
        y = f%val(il)
      else
        y = f%val(il) + ((f%val(iu)-f%val(il))/f%dt)*(t-(il-1)*f%dt)
      end if

    case(TDF_FROM_EXPR)
      call loct_parse_expression(fre, fim, 't', t, f%expression)
      y = cmplx(fre, fim)

    case default
      y = M_z0

    end select

  end function tdft
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_end(f)
    type(tdf_t), intent(inout) :: f
    call push_sub('tdfunction.tdf_init')

    select case(f%mode)
    case(TDF_FROM_FILE)
      call spline_end(f%amplitude)
    case(TDF_NUMERICAL)
      deallocate(f%val); nullify(f%val)
    case(TDF_SINE_SERIES)
      deallocate(f%coeffs); nullify(f%coeffs)
    end select
    f%mode = TDF_EMPTY

    call pop_sub()
  end subroutine tdf_end
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_copy(fout, fin)
    type(tdf_t), intent(inout) :: fout
    type(tdf_t), intent(in)  :: fin

    ASSERT( (fin%mode >= TDF_EMPTY)  .and. (fin%mode <= TDF_SINE_SERIES) )
    ASSERT( (fout%mode >= TDF_EMPTY)  .and. (fout%mode <= TDF_SINE_SERIES) )

    call tdf_end(fout)
    call tdf_init(fout)

    fout%t0     = fin%t0  
    fout%tau0   = fin%tau0
    fout%tau1   = fin%tau1
    fout%dt     = fin%dt 
    fout%a0     = fin%a0
    fout%omega0 = fin%omega0 
    fout%niter  = fin%niter
    fout%final_time = fin%final_time
    fout%init_time  = fin%init_time
    fout%expression = fin%expression
    fout%nfreqs     = fin%nfreqs
    if(fin%mode .eq. TDF_FROM_FILE) then
      fout%amplitude = fin%amplitude
    end if
    if(fin%mode .eq. TDF_NUMERICAL) then
      ALLOCATE(fout%val(fout%niter+1), fout%niter+1)
      fout%val  = fin%val
    end if
    if(fin%mode .eq. TDF_SINE_SERIES) then
      ALLOCATE(fout%coeffs(fout%nfreqs), fout%nfreqs)
      fout%coeffs = fin%coeffs
    end if
    fout%mode   = fin%mode

  end subroutine tdf_copy
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_scalar_multiply(alpha, f) 
    FLOAT, intent(in) :: alpha
    type(tdf_t), intent(inout) :: f

    select case(f%mode)
    case(TDF_CW, TDF_GAUSSIAN, TDF_COSINOIDAL, TDF_TRAPEZOIDAL)
      f%a0 = alpha*f%a0
    case(TDF_NUMERICAL)
      f%val = alpha*f%val
    case(TDF_SINE_SERIES)
      f%coeffs = alpha*f%coeffs
    case(TDF_FROM_FILE)
      call spline_times(alpha, f%amplitude)
    end select

  end subroutine tdf_scalar_multiply
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_cosine_multiply(omega, f) 
    FLOAT, intent(in) :: omega
    type(tdf_t), intent(inout) :: f

    integer :: j
    FLOAT :: t

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode .eq. TDF_NUMERICAL)

    do j = 1, f%niter + 1
      t = f%init_time + (j-1)*f%dt
      f%val(j) = f%val(j) * cos(omega*t)
    end do

  end subroutine tdf_cosine_multiply
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_cosine_divide(omega, f) 
    FLOAT, intent(in) :: omega
    type(tdf_t), intent(inout) :: f

    integer :: j
    FLOAT :: t

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode .eq. TDF_NUMERICAL)


    ! WARNING: no check is done for the case cos(omega*t) = 0
    do j = 1, f%niter + 1
      t = f%init_time + (j-1)*f%dt
      f%val(j) = f%val(j) / cos(omega*t)
    end do

  end subroutine tdf_cosine_divide
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine tdf_write(f, iunit)
    type(tdf_t), intent(in) :: f
    integer, intent(in) :: iunit

    select case(f%mode)
    case(TDF_CW)
      write(iunit,'(6x,a)')          'Mode: continuous wave.'
      write(iunit,'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
    case(TDF_GAUSSIAN)
      write(iunit,'(6x,a)')          'Mode: Gaussian envelope.'
      write(iunit,'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(6x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(6x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_COSINOIDAL)
      write(iunit,'(6x,a)') 'Mode: cosinoidal envelope.'
      write(iunit,'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(6x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(6x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_TRAPEZOIDAL)
      write(iunit,'(6x,a)') 'Mode: trapezoidal envelope.'
      write(iunit,'(6x,a,f10.4,a)')  'Amplitude: ', f%a0, ' [a.u]'
      write(iunit,'(6x,a,f10.4,3a)') 'Width:     ', f%tau0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(6x,a,f10.4,3a)') 'Middle t:  ', f%t0/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
      write(iunit,'(6x,a,f10.4,3a)') 'Ramp time: ', f%tau1/units_inp%time%factor, &
        ' [', trim(units_inp%time%abbrev), ']'
    case(TDF_FROM_FILE)
      write(iunit,'(6x,a)') 'Mode: time-dependent function read from file.'
    case(TDF_NUMERICAL)
      write(iunit,'(6x,a)') 'Mode: time-dependent function stored in a numerical array.'
    case(TDF_FROM_EXPR)
      write(iunit,'(6x,a)') 'Mode: time-dependent function parsed from the expression:'
      write(iunit,'(6x,a)') '      f(t) = '//trim(f%expression)
    end select
    if(f%omega0 .ne. M_ZERO) then
      write(iunit,'(6x,a,f10.4,3a)') 'Frequency: ', f%omega0/units_inp%energy%factor, &
        ' [', trim(units_inp%energy%abbrev), ']'
    end if

  end subroutine tdf_write
  !------------------------------------------------------------


  !------------------------------------------------------------
  ! Returns the dot product of f and g, defined as:
  !    < f | g > = \int_0^T dt f^*(t) g(t)
  ! It assumes that both f and m are in the same mode, otherwise
  ! it will fail and stop the code.
  !------------------------------------------------------------
  CMPLX function tdf_dot_product(f, g) result (fg)
    type(tdf_t), intent(in) :: f, g
    integer :: i
    FLOAT :: t

    fg = M_z0

    ! For the moment, we will just assume that f and g are of the same type.
    ASSERT(f%mode .eq. g%mode)

    select case(f%mode)
    case(TDF_NUMERICAL)
      ! We assume that the grid is the same for both functions.
      do i = 1, f%niter + 1
        fg = fg + conjg( f%val(i) ) * g%val(i)
      end do
      fg = fg * f%dt

    case(TDF_SINE_SERIES)
      ! We assume that the frequencies grid is the same for both functions
      do i = 1, f%nfreqs
        fg = fg + f%coeffs(i) * g%coeffs(i)
      end do
    case default
      do i = 1, f%niter + 1
        t = (i-1) * f%dt
        fg = fg + conjg( tdf(f, i) ) * tdf(g, i)
      end do

    end select

  end function tdf_dot_product
  !------------------------------------------------------------


end module tdf_m
