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
!! $Id$

#include "global.h"

module filter_m  
  use datasets_m
  use fft_m
  use global_m
  use io_m 
  use loct_m
  use parser_m
  use messages_m
  use profiling_m
  use string_m
  use tdf_m

  implicit none

  private
  
  public ::             &
       filter_t,        &
       filter_number,   &
       filter_init,     &
       filter_apply,    &
       filter_write,    &
       filter_end

  type filter_t
    private
    integer :: no_filters
    type(tdf_t), pointer :: f(:) => null()
    character(len=1024), pointer :: expression(:) => null()
    integer, pointer :: domain(:) => null()
  end type filter_t

  integer, parameter, public  ::  &
    filter_freq = 1,              &
    filter_time = 2
  
contains


  ! ---------------------------------------------------------
  subroutine filter_init(steps, dt, filter)    
    integer,          intent(in)  :: steps 
    FLOAT,            intent(in)  :: dt
    type(filter_t), intent(inout) :: filter

    type(block_t) :: blk
    integer :: no_c, i, no_f

    PUSH_SUB(filter_init)

    filter%no_filters = 0
    nullify(filter%domain)
    nullify(filter%expression)
        
    !%Variable OCTFilter
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% The block <tt>OCTFilter</tt> describes the type and shape of the filter function
    !% that are applied to the optimized laser field in each iteration.
    !% The filter forces the laser field to obtain the given form in frequency space.
    !% Each line of the block describes a filter; this way you can actually have more
    !% than one filter function (<i>e.g.</i> a filter in time and two in frequency space). 
    !% The filters are applied in the given order, <i>i.e.</i>, first the filter specified 
    !% by the first line is applied, then second line. 
    !% The syntax of each line is, then:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;domain | function
    !% <br>%</tt>
    !%  
    !%
    !% Possible arguments for domain are:
    !%  
    !% (i) <tt>frequency_filter</tt>: Specifies a spectral filter.
    !% 
    !% (ii) <tt>time_filter</tt>: DISABLED IN THIS VERSION.
    !%
    !% Example:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;time | "exp(-80*( w + 0.1567 )^2  ) + exp(-80*( w - 0.1567 )^2  )"
    !% <br>%</tt>
    !%
    !% Be careful that also the negative-frequency component is filtered since the resulting 
    !% field has to be real-valued.
    !%
    !%Option frequency_filter 1
    !% The filter is applied in the frequency domain.
    !%End
    if( parse_block(datasets_check('OCTFilter'),blk) == 0 ) then
      no_f = parse_block_n(blk)

      if(no_f <= 0) then
        POP_SUB(filter_init)
        return
      end if

      filter%no_filters = no_f
      SAFE_ALLOCATE(filter%f(1:no_f))
      SAFE_ALLOCATE(filter%expression(1:no_f))
      SAFE_ALLOCATE(filter%domain(1:no_f))

      do i=1, no_f
        no_c = parse_block_cols(blk, i-1)       
        call parse_block_integer(blk, i-1, 0, filter%domain(i))
        call parse_block_string(blk, i-1, 1, filter%expression(i))
        call conv_to_C_string(filter%expression(i))
        call tdf_init_numerical(filter%f(i), steps, dt, -M_ONE)
        call tdf_numerical_to_fourier(filter%f(i))
      end do
      call build_filter(filter)

      call parse_block_end(blk)
    end if

    POP_SUB(filter_init)
  end subroutine filter_init
  ! ---------------------------------------------------------

  
  ! ---------------------------------------------------------
  subroutine filter_apply(f, filter)  
    type(tdf_t), intent(inout) :: f
    type(filter_t), intent(in) :: filter

    integer        :: no_f, i, j, nfreqs

    PUSH_SUB(filter_apply)

    no_f = filter%no_filters
    if(no_f <= 0) then
      POP_SUB(filter_apply)
      return
    end if
    
    do i = 1, no_f
      write(message(1),'(a,i2)') "Info: Applying filter "         
      call messages_info(1)

      nfreqs = tdf_nfreqs(f)

      select case(filter%domain(i))
      case(filter_freq)
       call tdf_numerical_to_fourier(f)
       call tdf_set_numerical(f, 1, tdf(f, 1)*tdf(filter%f(i), 1))
       do j = 2, nfreqs !+ 1
         call tdf_set_numerical(f, j, tdf(f, j)*tdf(filter%f(i), j) / sqrt(M_TWO))
       end do
       do j = nfreqs + 1, 2*nfreqs - 1
         call tdf_set_numerical(f, j, tdf(f, j)*tdf(filter%f(i), j-nfreqs+1) / sqrt(M_TWO) )
       end do
       call tdf_fourier_to_numerical(f)

      case default
        message(1) = "...I don't know this filter type..."
        call messages_fatal(1)
      end select

    end do

    POP_SUB(filter_apply)
  end subroutine filter_apply
  ! ---------------------------------------------------------


  !------------------------------------------------
  subroutine build_filter(filter)
    type(filter_t), intent(inout) :: filter

    integer :: i, ip, j, nfreqs
    real(8) :: f_re, f_im
    CMPLX, allocatable :: ff(:)
    FLOAT, allocatable :: grid(:)

    PUSH_SUB(build_filter)

    do i = 1, filter%no_filters

      nfreqs = tdf_nfreqs(filter%f(i))

      SAFE_ALLOCATE(ff(1:tdf_nfreqs(filter%f(i))))
      SAFE_ALLOCATE(grid(1:tdf_nfreqs(filter%f(i))))
      ff = M_z0
      grid = M_ZERO

      select case(filter%domain(i))
      case(filter_freq)
        call tdf_fourier_grid(filter%f(i), grid)
        ff = M_z1
        do ip = 1, tdf_nfreqs(filter%f(i))
          call parse_expression(f_re, f_im, "w", real(grid(ip), 8), filter%expression(i))
          ff(ip) = f_re + M_zI*f_im
        end do
       
      case default
        write(message(1),'(a)') "Unknown choice of domain for filter."
        write(message(2),'(a)') "Choose: time or freq ."
        call messages_fatal(2)
      end select

      call tdf_set_numerical(filter%f(i), 1, real(ff(1), REAL_PRECISION))
      do j = 2, nfreqs !+ 1
         call tdf_set_numerical(filter%f(i), j, sqrt(M_TWO)*real(ff(j), REAL_PRECISION))
      end do
      ! WARNING: all the sine coefficients (imaginary part of ff) should be null.
      do j = nfreqs + 1, 2*nfreqs - 1
         call tdf_set_numerical(filter%f(i), j, -sqrt(M_TWO)*aimag(ff(j-nfreqs+1)))
      end do

      SAFE_DEALLOCATE_A(ff)
      SAFE_DEALLOCATE_A(grid)

    end do
    
    POP_SUB(build_filter)
  end subroutine build_filter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine filter_write(filter)
    type(filter_t), intent(in) :: filter

    integer :: kk, iunit, i, max_iter
    FLOAT :: dt
    character(len=80) :: filename
    FLOAT, allocatable :: wgrid(:)

    PUSH_SUB(filter_write)

    if(filter%no_filters <= 0) then
      POP_SUB(filter_write)
      return
    end if

    do kk = 1, filter%no_filters
      write(filename,'(a,i2.2)') OCT_DIR//'filter', kk
      max_iter = tdf_niter(filter%f(kk))
      dt = tdf_dt(filter%f(kk))
      iunit = io_open(filename, action='write')
      SAFE_ALLOCATE(wgrid(1:max_iter/2+1))
      call tdf_fourier_grid(filter%f(kk), wgrid)
      do i = 1, max_iter/2+1
        write(iunit, '(3es30.16e4)') wgrid(i), tdf(filter%f(kk), i)
      end do
      SAFE_DEALLOCATE_A(wgrid)
      call io_close(iunit)
    end do

    POP_SUB(filter_write)
  end subroutine filter_write
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------------
  subroutine filter_end(filter)
    type(filter_t), intent(inout) :: filter
    integer :: i

    PUSH_SUB(filter_end)

    if(filter%no_filters <= 0) then
      POP_SUB(filter_end)
      return
    end if

    do i = 1, filter%no_filters
      call tdf_end(filter%f(i))
    end do
    SAFE_DEALLOCATE_P(filter%f)
    filter%no_filters = 0
    SAFE_DEALLOCATE_P(filter%expression)
    SAFE_DEALLOCATE_P(filter%domain)

    POP_SUB(filter_end)
  end subroutine filter_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function filter_number(filter)
    type(filter_t), intent(in) :: filter

    PUSH_SUB(filter_number)
    filter_number = filter%no_filters

    POP_SUB(filter_number)
  end function filter_number
  ! ---------------------------------------------------------


end module filter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
