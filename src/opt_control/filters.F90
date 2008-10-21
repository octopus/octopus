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
  use global_m
  use messages_m
  use datasets_m
  use io_m 
  use loct_m
  use string_m
  use loct_parser_m
  use fft_m
  use units_m
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
    type(tdf_t), pointer :: f(:)
    character(len=1024), pointer :: expression(:)
    integer, pointer :: domain(:)
  end type filter_t

  integer, parameter, public  :: &
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

    call push_sub('filters.filter_init')

    filter%no_filters = 0
    nullify(filter%domain)
    nullify(filter%expression)
        
    !%Variable OCTFilter
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% The block OCTFilter describes the type and shape of the filter function
    !% that are applied to the optimized laser field in each iteration.
    !% The filter forces the laser field to obtain the given form in frequency space.
    !% Each line of the block describes a filter; this way you can actually have more
    !% than one filter function (e.g. a filter in time and two in frequency space). 
    !% The filters are applied in the given order, i.e., first the filter specified 
    !% by the first line is applied, then second line. This order is important if 
    !% the filters are conjugated, like time and frequency. If they are conjugated 
    !% the second filter can lift the action of the first one. Use with care.
    !% The syntax of each line is, then:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;domain | function
    !% <br>%</tt>
    !%  
    !%
    !% Possible arguments for domain are:
    !%  
    !% (i) frequency_filter : Specifies a spectral filter.
    !% 
    !% (ii) time_filter: Specifies a time-dependent envelope similar to a 
    !% time-dependent penalty. Use it in the case of a fixed fluence where a 
    !% time-dependent penalty is not possible.
    !%  
    !% An example for the block is:
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;time_filter | "exp(-gamma*( t - stime/4 )^2  )" 
    !% <br>&nbsp;&nbsp;time_filter | "exp(-gamma*( t - stime/4 )^2  )" 
    !% <br>%</tt>
    !% 
    !% Or a spectral Gaussian Filter at w=0.1567:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;time | "exp(-80*( w + 0.1567 )^2  ) + exp(-80*( w - 0.1567 )^2  )"
    !% <br>%</tt>
    !%
    !% Be careful that also the negative frequency component is filtered since the resulting 
    !% field has to be real valued.
    !%
    !%Option frequency_filter 1
    !% The filter is applied in the frequency domain
    !%Option time_filter 2
    !% The filter is applied in the time domain.
    !%End
    if( loct_parse_block(check_inp('OCTFilter'),blk) == 0 ) then
      no_f = loct_parse_block_n(blk)

      if(no_f <= 0) then
        call pop_sub(); return
      end if

      filter%no_filters = no_f
      ALLOCATE(filter%f(no_f), no_f)
      ALLOCATE(filter%expression(no_f), no_f)
      ALLOCATE(filter%domain(no_f), no_f)

      do i=1, no_f
        no_c = loct_parse_block_cols(blk, i-1)       
        call loct_parse_block_int(blk, i-1, 0, filter%domain(i))
        call loct_parse_block_string(blk, i-1, 1, filter%expression(i))
        call conv_to_C_string(filter%expression(i))
        call tdf_init_numerical(filter%f(i), steps, dt, -M_ONE)
        call tdf_fft_forward(filter%f(i))
      end do
      call build_filter(filter)

      call loct_parse_block_end(blk)
    end if

    call pop_sub()
  end subroutine filter_init
  ! ---------------------------------------------------------

  
  ! ---------------------------------------------------------
  subroutine filter_apply(f, filter)  
    type(tdf_t), intent(inout) :: f
    type(filter_t), intent(in) :: filter

    integer        :: steps, no_f, i, j

    call push_sub('filters.filter_apply')

    no_f = filter%no_filters
    if(no_f <= 0) then
      call pop_sub(); return
    end if
    
    do i = 1, no_f
      write(message(1),'(a,i2)') "Info: Applying filter "         
      call write_info(1)

      steps = tdf_niter(filter%f(i)) 

      select case(filter%domain(i))
      case(filter_freq)
       call tdf_fft_forward(f)
       do j = 1, steps/2+1
         call tdf_set_fourier(f, j, tdfw(f, j)* tdfw(filter%f(i), j) )
       end do
       call tdf_fft_backward(f)
      case(filter_time)
       do j = 1, steps + 1
         call tdf_set_numerical(f, j, tdf(f, j)* tdf(filter%f(i), j) )
       end do
      case default
        message(1) = "...I don't know this filter type..."
        call write_fatal(1)
      end select

    end do

    call pop_sub()
  end subroutine filter_apply
  ! ---------------------------------------------------------


  !------------------------------------------------
  subroutine build_filter(filter)
    type(filter_t), intent(inout) :: filter

    integer :: i, ip, steps
    real(8) :: f_re, f_im
    FLOAT   :: t, dt
    CMPLX, allocatable :: ff(:)
    FLOAT, allocatable :: grid(:)

    call push_sub('filter.build_filter')

    do i = 1, filter%no_filters

      steps = tdf_niter(filter%f(i))
      dt = tdf_dt(filter%f(i))

      ALLOCATE(ff((steps+1)/2+1), (steps+1)/2+1)
      ALLOCATE(grid((steps+1)/2+1), (steps+1)/2+1)
      ff = M_z0
      grid = M_ZERO

      select case(filter%domain(i))
      case(filter_time)
        do ip = 1, steps + 1
          t = (ip-1)*dt
          call loct_parse_expression(f_re, f_im, "t", real(t, 8), filter%expression(i))
          ff(ip) = f_re + M_zI*f_im
        end do
      
      case(filter_freq)
        call tdf_fourier_grid(filter%f(i), grid)
        ff = M_z1
        do ip = 1, steps/2+1
          call loct_parse_expression(f_re, f_im, "w", real(grid(ip), 8), filter%expression(i))
          ff(ip) = f_re + M_zI*f_im
        end do
       
      case default
        write(message(1),'(a)') "Unknown choice of domain for filter."
        write(message(2),'(a)') "Choose: time or freq ."
        call write_fatal(2)
      end select

      call tdf_set_fourier(filter%f(i), ff)

      deallocate(ff)
      deallocate(grid)

    end do
    
    call pop_sub()
  end subroutine build_filter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine filter_write(filter)
    type(filter_t), intent(in) :: filter

    integer :: kk, iunit, i, max_iter
    FLOAT :: dt
    character(len=80) :: filename
    FLOAT, allocatable :: wgrid(:)

    call push_sub('filters.filter_write')

    if(filter%no_filters <= 0) then
      call pop_sub(); return
    end if

    do kk = 1, filter%no_filters
      write(filename,'(a,i2.2)') 'opt-control/filter', kk

      max_iter = tdf_niter(filter%f(kk))
      dt = tdf_dt(filter%f(kk))

      if(filter%domain(kk) .eq. filter_freq) then
        iunit = io_open(filename, action='write')
        ALLOCATE(wgrid(max_iter/2+1), max_iter/2+1)
        call tdf_fourier_grid(filter%f(kk), wgrid)
        do i = 1, max_iter/2+1
          write(iunit, '(3es30.16e4)') wgrid(i), tdfw(filter%f(kk), i)
        end do
        deallocate(wgrid)
        call io_close(iunit)
      else
        iunit = io_open(filename, action='write')
        do i = 1, max_iter + 1
          write(iunit, '(4ES30.16E4)') (i-1)*dt, tdf(filter%f(kk), i)
        end do
        call io_close(iunit)
      end if

    end do

    call pop_sub()
  end subroutine filter_write
  ! ---------------------------------------------------------


  ! ----------------------------------------------------------------
  subroutine filter_end(filter)
    type(filter_t), intent(inout) :: filter
    integer :: i

    call push_sub('filter.filter_end')

    if(filter%no_filters <= 0) then
      call pop_sub(); return
    end if

    do i = 1, filter%no_filters
      call tdf_end(filter%f(i))
    end do
    deallocate(filter%f)
    filter%no_filters = 0
    deallocate(filter%expression)
    deallocate(filter%domain)

    call pop_sub()
  end subroutine filter_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function filter_number(filter)
    type(filter_t), intent(in) :: filter
    filter_number = filter%no_filters
  end function filter_number
  ! ---------------------------------------------------------


end module filter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
