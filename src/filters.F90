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
  use lib_oct_m
  use string_m
  use lib_oct_parser_m
#if defined(HAVE_FFT)
  use fft_m
#endif
  use units_m
  use tdf_m

  implicit none

  private
  
  public ::             &
       filter_t,        &
       def_tdpenalty,   &
       t_lookup,        &
       w_lookup,        &
       def_filter,      &
       build_filter,    &
       apply_filter,    &
       filter_write,    &
       filter_end

  
  type filter_t
    FLOAT :: weight       ! relative weight of filter
    type(tdf_t) :: f
    character(len=1024) :: expression
    integer :: domain, ftype
  end type filter_t

  integer, parameter, private  :: &
    timefreq = 1,                 &
    phase    = 2 
  
  integer, parameter, public  :: &
    filter_freq = 1,              &
    filter_time = 2,              &
    filter_phase= 3
  
  integer, parameter, private  :: &
    filter_pol_x  = 1,            &
    filter_pol_y  = 2,            &
    filter_pol_xy = 3
  
  integer, parameter, private  :: &
    filter_type_gauss     = 1,    &
    filter_type_gaussFWHM = 2,    &
    filter_type_sin2      = 3
  
  integer, parameter, private  :: &
    gauss      = 1,               &
    gaussFWHM  = 2,               &
    sin2       = 3
  
contains
  
  ! ---------------------------------------------------------
  subroutine def_tdpenalty(tdpenalty, steps, dt, mode_tdpenalty)
    FLOAT,        intent(out):: tdpenalty(:,:)
    integer,      intent(in) :: steps
    FLOAT,        intent(in) :: dt
    logical,      intent(out):: mode_tdpenalty

    C_POINTER                :: blk
    integer                  :: no_lines, no_col, i, j
    type(filter_t),  pointer :: tdp(:)

    call push_sub('filters.def_tdpenalty')
    mode_tdpenalty = .FALSE.

    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
      no_lines = loct_parse_block_n(blk)
      ALLOCATE(tdp(no_lines), no_lines)
      mode_tdpenalty = .TRUE.
      
      ! The structure of the block is:
      ! function | weight 
      do i = 1, no_lines
        no_col = loct_parse_block_cols(blk, i-1)

        ! parse formula string
        call loct_parse_block_string(                            &
          blk, i-1, 0, tdp(i)%expression)  
        ! convert to C string
        call conv_to_C_string(tdp(i)%expression)
        call loct_parse_block_float(blk, i-1, 1, tdp(i)%weight)
        !call loct_parse_block_int(blk, i-1, 1, tdp(i)%ftype)
        
        tdp(i)%domain = filter_time
        
        call build_filter(tdp(i), steps, dt)
      end do

      do i = 1, no_lines      
        tdp(i)%weight = tdp(i)%weight/SUM(tdp(1:no_lines)%weight)
      end do

      tdpenalty = M_ZERO
      do i=1,no_lines
        do j = 1, 2*steps+1
          tdpenalty(1, j-1) =  tdpenalty(1, j-1) + tdp(i)%weight * tdf(tdp(i)%f, j)
        end do
      end do

      tdpenalty = M_ONE / (tdpenalty + real(0.0000001,REAL_PRECISION )) 
      
      ! all we want to know is tdpenalty
      do i = 1, no_lines
        call tdf_end(tdp(i)%f)
      end do
      call loct_parse_block_end(blk)
    end if

    call pop_sub()
  end subroutine def_tdpenalty
  

  ! ---------------------------------------------------------
  subroutine apply_filter(steps, filtermode, filter, f)  
    type(filter_t),   pointer      :: filter(:)
    integer,          intent(in)   :: filtermode, steps
    type(tdf_t), intent(inout) :: f

    FLOAT, allocatable :: cfunction(:)

#if defined(HAVE_FFT)

    integer        :: no_f, i, kk, n(3), last, first, grouplength, ii, j
    CMPLX          :: tmp(0:2*steps, 1, 1), tmp2(0:2*steps, 1, 1) ! Have to be three-dimensional to use the fft_m modul
    CMPLX, allocatable :: filt(:), numerical(:)
    type(fft_t)    :: fft_handler

    call push_sub('filters.apply_filter')

    ALLOCATE(cfunction(0:2*steps), 2*steps+1)
    do i = 1, 2*steps+1
      cfunction(i-1) = tdf(f, i)
    end do

    n(1:3) = (/ 2*steps+1, 1, 1 /)
    
    call fft_init(n, fft_complex, fft_handler, optimize = .false.)
    
    no_f  = size(filter)
    i     = 0
    first = 0
    last  = 0

    ALLOCATE(filt(0:2*steps), 2*steps+1)
    ALLOCATE(numerical(0:2*steps), 2*steps+1)
    filt  = M_z0
    numerical = M_z0

    !do i=1, no_f
    do while(i.lt.no_f )
      i = i + 1
      ! decide time or freq
      write(message(1),'(a,i2)') "Info: Applying filter "         
      call write_info(1)

      do j = 0, 2*steps
        numerical(j) = tdf(filter(i)%f, j+1)
      end do  
      
      select case(filter(i)%domain)
        
      case(filter_freq)
        ! transform to freq space
        
        ! group filters
        ii = i + 1
        grouplength = 0
        filt(:) = numerical(:)
        do while(ii.lt.no_f+1) 
          if(filter(ii)%domain.eq.filter_freq) then
            grouplength = grouplength + 1
            do j = 0, 2*steps
              filt(j) = filt(j) + tdf(filter(ii)%f, j+1)
            end do
          end if
          ii = ii + 1
          first = i
        end do
        !this is out of the loop as a workaround for a bug in sun studio express
        last  = i + grouplength

        write(message(1),'(a,i2,a,i2)') 'Adding filters from: ',first,'to: ',last
        call write_info(1)

        i = ii - 1
          !write(6,*) 'in:', SUM(laser_inout(kk,:)**2)
          tmp2(:, 1, 1) = cmplx(cfunction(:))
          
          call zfft_forward(fft_handler, tmp2, tmp)
          !write(6,*) SUM(abs(tmp)**2)/real((2*steps+1), REAL_PRECISION)
          !
          tmp(:, 1, 1) = tmp(:, 1, 1)*filt(:)
          call zfft_backward(fft_handler, tmp, tmp2)
          
          cfunction(:) = real(tmp2(:, 1, 1), REAL_PRECISION)!/real((2*steps+1), REAL_PRECISION)
          !write(6,*) 'out:', SUM(laser_inout(kk,:)**2)
          !write(6,*) SUM(laser_inout(kk,:)**2)          
        
      case(filter_time)
        ! group filters
        ii = i + 1
        grouplength = 0
        filt(:) = numerical(:)
        do while(ii.lt.no_f+1) 
          if(filter(ii)%domain.eq.filter_time) then
            grouplength = grouplength + 1
            do j = 0, 2*steps
              filt(j) = filt(j) + tdf(filter(ii)%f, j+1)
            end do
          end if
          first = i
          ii = ii + 1
        end do
        !this is out of the loop as a workaround for a bug in sun studio express
        last  = i + grouplength
        i = ii - 1

        write(message(1),*) 'Adding filters from: ',first,'to: ',last
        call write_info(1)

        cfunction(:) = cfunction(:)*filt(:)
        
      case(filter_phase)
        ! phase only optimization
        message(1) = "to be implemented soon"
        call write_fatal(1)
      case default
        message(1) = "...I don't know this filter type..."
        call write_fatal(1)
      end select
    end do

    do i = 1, 2*steps+1
      call tdf_set_numerical(f, i, cfunction(i-1))
    end do

    deallocate(cfunction)    
    deallocate(filt)
    call fft_end(fft_handler)
    call pop_sub()
#else
    write(message(1),'(a)') "Internal error in filters.apply_filter"
    call write_fatal(1)
#endif
  end subroutine apply_filter
  
  
  ! ---------------------------------------------------------
  subroutine def_filter(steps, dt, filtermode, f)    
    integer,          intent(in)  :: filtermode,steps 
    type(filter_t),   pointer     :: f(:)
    FLOAT,            intent(in)  :: dt

    C_POINTER :: blk
    integer :: no_c, i, no_f

    call push_sub('filters.def_filter')
        
    ! filter function
    !%Variable OCTFilter
    !%Type block
    !%Section Optimal Control
    !%Description
    !% The block OCTFilter describes the type and shape of the filter function
    !% that are applied to the optimized laser field in each iteration.
    !% The filter forces the laser field to obtain the given form in frequency space.
    !% Each line of the block describes a filter; this way you can actually have more
    !% than one filter function (e.g. a filter in time and two in frequency space). The filters are 
    !% applied in the given order, i.e., first the filter specified by the first line is applied, then 
    !% second line. This order is important if the filters are conjugated, like time and frequency. 
    !% If they are conjugated the second filter can lift the action of the first one. Use with care.
    !% The syntax of each line is, then:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;domain | pol | pol | pol | function
    !% <br>%</tt>
    !%  
    !%
    !% Possible arguments for domain are:
    !%  
    !% *frequency_filter* - Specifies a spectral filter.
    !% 
    !% *time_filter*  - 
    !%
    !% Specifies a time-dependent envelope similar to a time-dependent penalty.
    !% Use it in the case of a fixed fluence where a time-dependent penalty is not possible.
    !% Characterize the polarization vector with pol_x, pol_y, pol_z.
    !% Finally define the filter function by specifying the "function".
    !% Some examples for "function" are:
    !%
    !%  
    !% An example for the block is:
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;time | 1 | 0 | 0 | "exp(-gamma*( t - stime/4 )^2  )" 
    !% <br>&nbsp;&nbsp;time | 0 | 1 | 0 | "exp(-gamma*( t - stime/4 )^2  )" 
    !% <br>&nbsp;&nbsp;
    !% <br>&nbsp;&nbsp;
    !% <br>%</tt>
    !%   
    !% 
    !% Or a spectral Gaussian Filter at w=0.1567:
    !%
    !% <tt>%OCTFilter
    !% <br>&nbsp;&nbsp;time | 1 | 0 | 0 | "exp(-80*( w + 0.1567 )^2  ) + exp(-80*( w - 0.1567 )^2  )"
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

    if((filtermode==timefreq).AND.(loct_parse_block(check_inp('OCTFilter'),blk)==0)) then
      no_f = loct_parse_block_n(blk)
      ALLOCATE(f(no_f), no_f)
      do i=1, no_f
        ! The structure of the block is:
        ! domain | function_type | center | width | weight 
        no_c = loct_parse_block_cols(blk, i-1)       
        call loct_parse_block_int(blk, i-1, 0, f(i)%domain)

        ! parse formula string
        call loct_parse_block_string(blk, i-1, 1, f(i)%expression)

        ! convert to C string
        call conv_to_C_string(f(i)%expression)
        
        ! WARNING: One should actually read the weight?
        !call loct_parse_block_int(blk, i-1, 2, f(i)%ftype)
        !call loct_parse_block_float(blk, i-1, 5, f(i)%weight)

        call tdf_init_numerical(f(i)%f, 2*steps, dt*M_HALF)
        call build_filter(f(i), steps, dt)
      end do
      
      call loct_parse_block_end(blk)
      
    end if

    call pop_sub()
  end subroutine def_filter


  !------------------------------------------------
  subroutine build_filter(fp, steps, dt)
    type(filter_t), intent(inout) :: fp
    integer,        intent(in)    :: steps
    FLOAT,          intent(in)    :: dt

    integer :: i, iunit
    integer :: ip, pol
    FLOAT   :: grid(0:2*steps), f_re, f_im
    CMPLX, allocatable :: ff(:)


    call push_sub('filter.build_filter_')

    ALLOCATE(ff(0:2*steps), 2*steps+1)

    ff = M_z0
    select case(fp%domain)
    case(filter_time)
      call t_lookup(2*steps+1,dt/real(2.0,REAL_PRECISION),grid)
      do ip=0, 2*steps
        call loct_parse_expression(f_re, f_im, "t", grid(ip), fp%expression)
        ff(ip) = f_re + M_zI*f_im
      end do
      
    case(filter_freq)
      call w_lookup(2*steps+1,dt,grid)
      ff = M_z1
      do ip=0, 2*steps
        call loct_parse_expression(f_re, f_im, "w", grid(ip), fp%expression)      
        ff(ip) = f_re + M_zI*f_im
      end do
      
    case default
      write(message(1),'(a)') "Unknown choice of domain for filter."
      write(message(2),'(a)') "Choose: time or freq ."
      call write_fatal(2)
    end select
    
    
    iunit = io_open('opt-control/filtertest', action='write')
    do i = 0, steps
      write(iunit, '(4es20.12)') grid(i), ff(i)
    end do
    call io_close(iunit)
    
    call tdf_set_numerical(fp%f, ff(0:2*steps))
    
    deallocate(ff)
    call pop_sub()
  end subroutine build_filter


  ! ---------------------------------------------------------
  subroutine t_lookup(steps, dt, tgrid)
    integer, intent(in)  :: steps
    FLOAT,   intent(in)  :: dt
    FLOAT,   intent(out) :: tgrid(:)

    integer  :: i
 
    call push_sub('filter.t_lookup_')

    do i=1, steps
       tgrid(i) = ( real(i, REAL_PRECISION) - M_HALF ) * dt 
    enddo

    call pop_sub()
  end subroutine t_lookup


  ! ---------------------------------------------------------
  subroutine w_lookup(steps, dt, wgrid)
    integer, intent(in)  :: steps
    FLOAT,   intent(in)  :: dt
    FLOAT,   intent(out) :: wgrid(:)
 
    integer :: i
    FLOAT   :: df

    call push_sub('filter.w_lookup')

    wgrid = M_ZERO

    df = M_TWO/(real(steps,REAL_PRECISION) * dt)
    do i=1, int((steps)/2)
       wgrid(i)   = real(i-1, REAL_PRECISION) * df * M_TWO * M_PI
    enddo
    do i=int((steps)/2), steps-1 
       wgrid(i+1) = real(i-steps, REAL_PRECISION) * df * M_TWO * M_PI
    enddo   
    
    call pop_sub()
  end subroutine w_lookup


  ! ----------------------------------------------------------------
  subroutine filter_write(f, ndim, dt, max_iter)
    type(filter_t), pointer :: f(:)
    integer, intent(in)     :: ndim
    FLOAT, intent(in)       :: dt
    integer, intent(in)     :: max_iter

    integer :: kk, no_f, iunit, i
    character(len=80) :: filename
    CMPLX, allocatable :: numerical(:)

    call push_sub('filters.filter_write')

    if(.not.associated(f)) then
      call pop_sub(); return
    endif
    no_f = size(f)
    if(no_f .eq. 0) then
      call pop_sub(); return
    end if

    ALLOCATE(numerical(0:2*max_iter), 2*max_iter+1)
    do kk=1, size(f)
       write(filename,'(a,i2.2)') 'opt-control/filter', kk


       if(f(kk)%domain.eq.1) then
         do i = 1, 2*max_iter+1
           numerical(i-1) = tdf(f(kk)%f, i)
         end do
         call write_fieldw(filename, ndim, max_iter, real(numerical, REAL_PRECISION), dt)
       else
         !write(6,*) f(kk)%numerical(:,0)
         iunit = io_open(filename, action='write')
         do i = 0, 2*max_iter
           write(iunit, '(4ES30.16E4)') i*dt*M_HALF, tdf(f(kk)%f, i+1)
         end do
         call io_close(iunit)
       end if
    end do

    deallocate(numerical)
    call pop_sub()
  end subroutine filter_write


  ! ----------------------------------------------------------------
  subroutine filter_end(f)
    type(filter_t), pointer    :: f(:)

    integer :: i, no_f
    
    if(.not.associated(f)) return

    call push_sub('filter.filter_end_')

    no_f = size(f)
    if(no_f <= 0) return

    do i = 1, no_f
      call tdf_end(f(i)%f)
    end do

    deallocate(f)
    nullify(f)

    call pop_sub()
  end subroutine filter_end


  ! ---------------------------------------------------------
  subroutine write_fieldw(filename, ndim, steps, las, dt)
    ! in w=(2pi f) space
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: ndim
    integer,          intent(in) :: steps
    FLOAT,            intent(in) :: las(0:2*steps)
    FLOAT,            intent(in) :: dt
  
    integer :: i, iunit
    FLOAT   :: wgrid(0:2*steps)

    call push_sub('opt_control.write_fieldw')

    call w_lookup(2*steps+1, dt, wgrid)

    iunit = io_open(filename, action='write')
    do i = 0, 2*steps
       write(iunit, '(4es30.16e4)') wgrid(i), las(i)
    end do
    call io_close(iunit)
    
    call pop_sub()
  end subroutine write_fieldw

end module filter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
