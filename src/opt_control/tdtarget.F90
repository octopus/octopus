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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"

module opt_control_tdtarget_m
  use datasets_m
  use varinfo_m
  use messages_m
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use global_m
  use grid_m
  use mesh_function_m
  use states_m
  use hamiltonian_m
  use timedep_m
  use td_write_m
  use opt_control_constants_m

  implicit none

  private
  public :: td_target_t,      &
            tdtarget_end, &
            build_tdshape, &
            build_tdtarget

  type td_target_t
    integer :: type       ! local or state 
    FLOAT   :: par 
    FLOAT   :: width      ! width parameter
    FLOAT   :: weight     ! relative weight of filter
    character(len=1024) :: expression(MAX_DIM) ! expression of user def func
    integer             :: ftype               ! which function: gaussian etc
    CMPLX, pointer :: tdshape(:,:) ! evaluate user defined functions
  end type td_target_t

  contains

  !----------------------------------------------------------
  subroutine tdtarget_end(tdtg)
    type(td_target_t), intent(inout) :: tdtg
    call push_sub('opt_control_tdtarget.tdtarget_end')
    deallocate(tdtg%tdshape); nullify(tdtg%tdshape)
    call pop_sub()
  end subroutine tdtarget_end


  !----------------------------------------------------------
  ! 
  !----------------------------------------------------------
  subroutine build_tdshape(tdtg, gr, steps, dt)
    type(td_target_t), intent(inout) :: tdtg
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: steps
    FLOAT,             intent(in)    :: dt

    FLOAT   :: f_re, f_im, t
    integer :: kk, pol
    call push_sub('opt_control_tdtarget.build_tdshape')

    do kk=0, steps
      t = (kk-1)*steps
      do pol=1, NDIM
        f_re = M_ZERO
        f_im = M_ZERO
        call loct_parse_expression(f_re, f_im, & 	
          M_ZERO, M_ZERO, M_ZERO, M_ZERO, t, &
          tdtg%expression(pol)) 
        tdtg%tdshape(pol,kk) = f_re + M_zI*f_im
      end do
    end do

    call pop_sub()
  end subroutine build_tdshape


  ! ---------------------------------------------------------
  subroutine build_tdtarget(tdtg, gr, tgt, iter)
    type(td_target_t), intent(in)   :: tdtg
    type(grid_t),      intent(in)   :: gr
    integer,           intent(in)   :: iter
    CMPLX  ,           intent(out)  :: tgt(:)

    integer :: pol
    call push_sub('opt_control_tdtarget.build_tdtarget')
    
    select case(tdtg%ftype)
    case(1)
      ! gauss
      ! Product x = gr%m%x(ip, :) 
      tgt = M_ONE
      do pol=1, NDIM
        tgt = tgt * exp( - (gr%m%x(:, pol) - &
          real(tdtg%tdshape(pol,iter),REAL_PRECISION ))**2/(M_TWO*tdtg%width**2))
      enddo
      
    case default
      message(1)="Unknown TD Target Operator: Choose Gaussian(1) or Step(2)"
      call write_fatal(1)
    end select
    
    call pop_sub()
  end subroutine build_tdtarget


end module opt_control_tdtarget_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

