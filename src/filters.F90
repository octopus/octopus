!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: filters.F90,v 1.42 2006/03/18 17:56:50 appel Exp $
#include "global.h"

module filter_m  
  use global_m
  use grid_m
  !use varinfo_m
  use messages_m
  use datasets_m
  use io_m 
  use lib_oct_m
  use string_m
  use lib_oct_parser_m
  use fft_m
  !use lib_oct_gsl_spline_m
  use units_m
  use mesh_m
  !use simul_box_m
  !use output_m

  implicit none

  type filter_t
     FLOAT :: weight       ! relative weight of filter
     FLOAT :: center       ! center in freq or time 
     FLOAT :: width        ! width in freq or time

     CMPLX :: fpol(MAX_DIM)
     CMPLX, pointer :: numerical(:,:) ! (ndim,0:2*steps) store values of filter
     character(len=1024) :: expression

     integer :: domain, ftype

  end type filter_t

  integer, parameter, private  ::  &
       timefreq = 1, &
       phase    = 2 

  integer, parameter, private  ::  &
       filter_freq = 1,            &
       filter_time = 2,            &
       filter_phase= 3

  integer, parameter, private  ::  &
       filter_pol_x  = 1,           &
       filter_pol_y  = 2,           &
       filter_pol_xy = 3
  
  integer, parameter, private  ::  &
       filter_type_gauss      = 1,           &
       filter_type_gaussFWHM  = 2,           &
       filter_type_sin2       = 3
 
 integer, parameter, private  ::  &
       gauss      = 1,           &
       gaussFWHM  = 2,           &
       sin2       = 3

contains
  
  subroutine def_tdpenalty(gr, tdpenalty, steps, dt, mode_tdpenalty)
    type(grid_t), pointer    :: gr
    FLOAT,        intent(out):: tdpenalty(:,:)
    integer,      intent(in) :: steps
    FLOAT,        intent(in) :: dt
    logical,      intent(out):: mode_tdpenalty

    integer(POINTER_SIZE)    :: blk
    integer                  :: no_lines, no_col, i
    type(filter_t),   pointer :: tdp(:)

    call push_sub('filter.def_tdpenalty_')
    mode_tdpenalty = .FALSE.
    !%Variable OCTLaserEnvelope
    !%Type Block
    !%Section Optimal Control
    !%Description
    !% Often a predefined time-dependent envelope on the laser field is required. 
    !% This can be achieved by making the penalty factor time-dependent. 
    !% Here, you may specify the required time dependent envelope.
    !% The code allows for more than one enevelope, all of them will be added together and can weighted against each other. 
    !% It is possible to choose different envelopes for different polarization directions. 
    !% The block is similar to the time filter.
    !%End
        
    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
       no_lines = loct_parse_block_n(blk)
       ALLOCATE(tdp(no_lines), no_lines)
       mode_tdpenalty = .TRUE.

       ! The structure of the block is:
       ! polarization | function | weight 
       do i = 1, no_lines
          no_col = loct_parse_block_cols(blk, i-1)
          call loct_parse_block_cmplx(blk, i-1, 0, tdp(i)%fpol(1)) 
          call loct_parse_block_cmplx(blk, i-1, 1, tdp(i)%fpol(2))
          call loct_parse_block_cmplx(blk, i-1, 2, tdp(i)%fpol(3))
          ! parse formula string
          call loct_parse_block_string(                            &
               blk, i-1, 3, tdp(i)%expression)  
          ! convert to C string
          call conv_to_C_string(tdp(i)%expression)
          call loct_parse_block_float(blk, i-1, 4, tdp(i)%weight)
          !call loct_parse_block_int(blk, i-1, 1, tdp(i)%ftype)
          !call loct_parse_block_float(blk, i-1, 2, tdp(i)%center)
          !call loct_parse_block_float(blk, i-1, 3, tdp(i)%width)

          tdp(i)%domain = filter_time

          ALLOCATE(tdp(i)%numerical(NDIM,0:2*steps), NDIM*(2*steps+1))
          call build_filter(gr, tdp(i), steps, dt)
       end do
       
       tdp(i)%weight = tdp(i)%weight/SUM(tdp(:)%weight)

       tdpenalty = M_ZERO
       do i=1,no_lines
          
          tdpenalty =  tdpenalty &
               + tdp(i)%weight * real(tdp(i)%numerical,PRECISION)
       end do

       tdpenalty = M_ONE / (tdpenalty + real(0.0000001,PRECISION )) 

       ! all we want to know is tdpenalty
       do i = 1, no_lines
          if(associated(tdp(i)%numerical)) then
             deallocate(tdp(i)%numerical) 
             nullify(tdp(i)%numerical)
          end if
       end do
       call loct_parse_block_end(blk)
    end if
    call pop_sub()

  end subroutine def_tdpenalty
  
  ! ---------------------------------------------------------
  subroutine apply_filter(gr, steps, filtermode, filter, laser_inout)  
    type(grid_t),     pointer      :: gr
    type(filter_t),   pointer      :: filter(:)
    integer,          intent(in)   :: filtermode, steps
    FLOAT,            intent(inout):: laser_inout(:,:)

    integer        :: no_f, i, kk, n(3)
    CMPLX          :: tmp(0:2*steps, 1, 1), tmp2(0:2*steps, 1, 1) ! Have to be three-dimensional to use the fft_m module
    type(fft_t)    :: fft_handler

    call push_sub('filter.apply_filter_')

    n(1:3) = (/ 2*steps+1, 1, 1 /)
    call fft_init(n, fft_complex, fft_handler, optimize = .false.)
    
    no_f = size(filter)
    !kk   = M_ZERO
    ! apply filters line by line
    do i=1, no_f
    !do while(kk.lt.no_f)
    !   kk = kk + 1
    !   groupmembers = 0
    !   do i=kk+1, no_f-1
    !      if(filter(i)%group.eq.filter(kk)%group).AND.(filter(i)%domain.eq.filter(kk)%domain) then
    !      groupmembers = groupmembers + 1 
    !   end do
            
       ! decide time or freq
       write(message(1),'(a,i2)') "Info: Applying filter ",i         
       call write_info(1)  
       
       select case(filter(i)%domain)
    
       case(filter_freq)
          ! transform to freq space
          do kk=1, NDIM
             !write(6,*) 'in:', SUM(laser_inout(kk,:)**2)
             tmp2(:, 1, 1) = cmplx(laser_inout(kk,:))
             call zfft_forward(fft_handler, tmp2, tmp)
             !write(6,*) SUM(abs(tmp)**2)/real((2*steps+1), PRECISION)
             !
             tmp(:, 1, 1) = tmp(:, 1, 1)*filter(i)%numerical(kk,:)
             call zfft_backward(fft_handler, tmp, tmp2)
    
             laser_inout(kk,:) = real(tmp2(:, 1, 1), PRECISION)/real((2*steps+1), PRECISION)
             !write(6,*) 'out:', SUM(laser_inout(kk,:)**2)
             !write(6,*) SUM(laser_inout(kk,:)**2)          
          enddo
          
       case(filter_time)
          ! filter in time
          !if(filter(i)%group)
          do kk=1, NDIM
             laser_inout(kk,:) = laser_inout(kk,:)*filter(i)%numerical(kk,:)
          enddo
          
       case(filter_phase)
          ! phase only optimization
          message(1) = "to be implemented soon"
          call write_fatal(1)
       case default
          message(1) = "...I don't know this filter type..."
          call write_fatal(1)
       end select
    end do
    
    call fft_end(fft_handler)
    call pop_sub()
  end subroutine apply_filter

  ! ---------------------------------------------------------
  subroutine def_filter(gr, steps, dt, filtermode, f)    
    type(grid_t),     pointer     :: gr
    integer,          intent(in)  :: filtermode,steps 
    type(filter_t),   pointer     :: f(:)
    FLOAT,            intent(in)  :: dt

    integer(POINTER_SIZE) :: blk
    integer :: no_c, i, no_f, ii

    call push_sub('filter.def_filter')
    
    
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
    !%An example for the block is:
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
    !% Be careful that also the negative frequency component is filtered since the resulting field has to be real valued.
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
          call loct_parse_block_cmplx(blk, i-1, 1, f(i)%fpol(1))
          call loct_parse_block_cmplx(blk, i-1, 2, f(i)%fpol(2))
          call loct_parse_block_cmplx(blk, i-1, 3, f(i)%fpol(3))
          ! parse formula string
          call loct_parse_block_string(                            &
               blk, i-1, 4, f(i)%expression)
          ! convert to C string
          call conv_to_C_string(f(i)%expression)

!          call loct_parse_block_int(blk, i-1, 2, f(i)%ftype)
!          call loct_parse_block_float(blk, i-1, 3, f(i)%center)
!          call loct_parse_block_float(blk, i-1, 4, f(i)%width)
!          call loct_parse_block_float(blk, i-1, 5, f(i)%weight)
          ALLOCATE(f(i)%numerical(NDIM,0:2*steps), NDIM*(2*steps+1))
          call build_filter(gr, f(i), steps, dt)
       end do
           
       call loct_parse_block_end(blk)

       end if
    call pop_sub()
    
  end subroutine def_filter

  !------------------------------------------------
  subroutine build_filter(gr, fp, steps, dt)
    implicit none
    type(grid_t),   pointer       :: gr
    type(filter_t), intent(inout) :: fp
    integer,        intent(in)    :: steps
    FLOAT,          intent(in)    :: dt
    integer :: i, iunit
    integer :: kk, ip
    FLOAT   :: grid(0:2*steps)
    FLOAT   :: width, f_re, f_im
    CMPLX   :: ff(NDIM,0:2*steps)

    call push_sub('filter.build_filter_')

    ff = M_z0
    select case(fp%domain)
    case(filter_time)
       call t_lookup(2*steps+1,dt/real(2.0,PRECISION),grid)
       do ip=0, 2*steps
          call loct_parse_expression(f_re, f_im, "t", grid(ip), fp%expression)
          ! FIXME: polarization   
          ff(:,ip) = f_re + M_zI*f_im 
          !write(6,*) f_re, f_im
       end do
   

    case(filter_freq)
       call w_lookup(2*steps+1,dt,grid)
       do ip=0, 2*steps
          call loct_parse_expression(f_re, f_im, "w", grid(ip), fp%expression)      
          ! FIXME: polarization
          ff(:,ip) = f_re + M_zI*f_im 
       end do
       
    case default
       write(message(1),'(a)') "Unknown choice of domain for filter."
       write(message(2),'(a)') "Choose: time or freq ."
       call write_fatal(2)
    end select

    

    !select case (fp%ftype)
    !case(filter_type_gauss)
    !   write(message(1), '(a)') 'Info: Building Gaussian'
    !   call write_info(1)
    !   do kk=1, NDIM
    !      ff(kk,:) = m_z0
    !      ff(kk,:) =  exp(-(grid - fp%center)**2 / (M_TWO*fp%width**2))
    !   end do
    !   if(fp%domain.eq.filter_freq) then
    !     do kk=1, NDIM
    !        ff(kk,:) =  ff(kk,:) + &
    !             exp(-(grid + fp%center)**2 / (M_TWO*fp%width**2))
    !     end do
    !  end if
 
    iunit = io_open('opt-control/filtertest', action='write')
    do i = 0, steps
       write(iunit, '(4es20.12)') grid(i), ff(:, i)
    end do
    call io_close(iunit)
    
    !case(filter_type_gaussFWHM)
    !   write(message(1), '(a)') 'Info: Building Gaussian FWHM'
    !   call write_info(1)
    !   do kk=1, NDIM
     !     width = sqrt((fp%width-fp%center)**2/(M_TWO*log(M_TWO)))
    !      ff(kk,:) = exp(-(grid/M_TWO - fp%center)**2 / (M_TWO*width**2))
    !   enddo
  

    !case default
    !   write(message(1),'(a)') "Unknown filter function."
    !   call write_fatal(1)
    !end select
    
      
    fp%numerical(:,:) = ff(:,:)
    
    call pop_sub()

  end subroutine build_filter

  ! ---------------------------------------------------------
  subroutine t_lookup(steps,dt,tgrid)
    implicit none
    
    integer, intent(in)  :: steps
    FLOAT,   intent(in)  :: dt
    FLOAT,   intent(out) :: tgrid(:)

    integer  :: i
    call push_sub('filter.t_lookup_')

    do i=1, steps
       tgrid(i) = ( real(i, PRECISION) - M_HALF ) * dt 
    enddo

    call pop_sub()

  end subroutine t_lookup

  ! ---------------------------------------------------------
  subroutine w_lookup(steps,dt,wgrid)
    implicit none
   
    integer, intent(in)  :: steps
    FLOAT,   intent(in)  :: dt
    FLOAT,   intent(out) :: wgrid(:)
 
    integer :: i
    FLOAT   :: df

    call push_sub('filter.w_lookup')

    wgrid = M_ZERO

    df = M_TWO/(real(steps,PRECISION) * dt)
    do i=1, int((steps)/2)
       wgrid(i)   = real(i-1, PRECISION) * df * M_TWO * M_PI
    enddo
    do i=int((steps)/2), steps-1 
       wgrid(i+1) = real(i-steps, PRECISION) * df * M_TWO * M_PI
    enddo   


    call pop_sub()
    
  end subroutine w_lookup


! ----------------------------------------------------------------
  subroutine filter_end(f)
    type(filter_t), pointer    :: f(:)

    integer :: i, no_f
    
    if(.not.associated(f)) return

    call push_sub('filter.filter_end_')

    no_f = size(f)
    if(no_f <= 0) return

    do i = 1, no_f
       
       if(associated(f(i)%numerical)) then
          deallocate(f(i)%numerical) 
          nullify(f(i)%numerical)
       end if
       
    end do

    deallocate(f)
    nullify(f)

    call pop_sub()

  end subroutine filter_end

end module filter_m
