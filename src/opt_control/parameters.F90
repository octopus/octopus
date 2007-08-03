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

module opt_control_parameters_m
  use string_m
  use datasets_m
  use varinfo_m
  use global_m
  use messages_m
  use io_m
  use lib_oct_parser_m
  use lib_oct_m
  use filter_m
  use external_pot_m
  use tdf_m

  implicit none

  private
  public :: oct_control_parameters_t,     &
            parameters_init,              &
            parameters_set,               &
            parameters_end,               &
            parameters_to_h,              &
            parameters_write,             &
            laser_fluence,                &
            j2_functional


  type oct_control_parameters_t
    integer :: no_parameters
    FLOAT   :: dt
    integer :: ntiter
    type(tdf_t), pointer :: f(:)
    CMPLX, pointer  :: laser_pol(:, :)

    FLOAT, pointer   :: a_penalty(:)
    type(tdf_t), pointer :: td_penalty(:)
  end type oct_control_parameters_t


  contains

  ! ---------------------------------------------------------
  subroutine parameters_init(cp, no_parameters, dt, ntiter, octiter)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer, intent(in) :: no_parameters   
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter, octiter
    integer :: j

    call push_sub('opt_control_parameters.parameters_init')

    cp%no_parameters = no_parameters
    cp%dt = dt
    cp%ntiter = ntiter
    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    do j = 1, cp%no_parameters
      call tdf_init_numerical(cp%f(j), cp%ntiter, cp%dt)
    end do

    ALLOCATE(cp%laser_pol(MAX_DIM, cp%no_parameters), MAX_DIM*cp%no_parameters)
    cp%laser_pol = M_z0

    call parameters_penalty_init(cp, octiter)

    call pop_sub()
  end subroutine parameters_init


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(epot_t), intent(in) :: ep
    integer :: j

    call push_sub('opt_control_parameters.parameters_set')

    do j = 1, cp%no_parameters
      cp%f(j) = ep%lasers(j)%f
      cp%laser_pol(:, j) = ep%lasers(j)%pol(:)
    end do

    call pop_sub()
  end subroutine parameters_set


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep

    integer :: j

    call push_sub('opt_control_paramters.parameters_to_h')

    do j = 1, cp%no_parameters
      ep%lasers(j)%f = cp%f(j)
      ep%lasers(j)%pol(:) = cp%laser_pol(:, j)
    end do

    call pop_sub()
  end subroutine parameters_to_h


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j

    call push_sub('opt_control_parameters.parameters_end')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
      call tdf_end(cp%td_penalty(j))
    end do
    deallocate(cp%f); nullify(cp%f)
    deallocate(cp%laser_pol)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, cp)
    character(len=*), intent(in) :: filename
    type(oct_control_parameters_t), intent(in) :: cp

    integer :: i, j, iunit
    FLOAT :: t
    character(len=2) :: digit

    call push_sub('opt_control_parameters.parameters_write')

    call io_mkdir(trim(filename))
    do j = 1, cp%no_parameters
      if(cp%no_parameters > 1) then
        write(digit,'(i2.2)') j
        iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
      else
        iunit = io_open(trim(filename)//'/cp', action='write')
      end if
      do i = 1, cp%ntiter + 1
        t = (i-1)*cp%dt
        write(iunit, '(3es20.8e3)') t, tdf(cp%f(j), t)
      end do
      call io_close(iunit)
    end do

    call pop_sub()
  end subroutine parameters_write


  ! ---------------------------------------------------------
  subroutine parameters_penalty_init(par, ctr_iter_max)
    type(oct_control_parameters_t), intent(inout) :: par
    integer,      intent(in)                   :: ctr_iter_max

    character(len=1024)      :: expression
    FLOAT                    :: weight, t, f_re, f_im, octpenalty, dt
    C_POINTER                :: blk
    integer                  :: no_lines, i, j, steps
    call push_sub('opt_control_penalty.parameters_penalty_init')

    dt = par%dt
    steps = par%ntiter

    !%Variable OCTPenalty
    !%Type float
    !%Section Optimal Control
    !%Default 1.0
    !%Description
    !% The variable specificies the value of the penalty factor for the 
    !% integrated field strength (fluence). Large value - small fluence. The value is 
    !% always positive. A transient shape can be specified using the block OCTLaserEnvelope.
    !% In this case OCTPenalty is multiplied with time-dependent function. 
    !% The value depends on the coupling between the states. A good start might be a 
    !% value from 0.1 (strong fields) to 10 (weak fields). 
    !%End
    call loct_parse_float(check_inp('OCTPenalty'), M_ONE, octpenalty)
    ! penalty array for fixed fluence run 
    ! the array is only interesting for the development of new algorithms

    ALLOCATE(par%a_penalty(0:ctr_iter_max+1), ctr_iter_max+2)
    par%a_penalty = octpenalty

    ALLOCATE(par%td_penalty(par%no_parameters), par%no_parameters)
    do i = 1, par%no_parameters
      call tdf_init_numerical(par%td_penalty(i), steps, dt, initval = M_z1)
    end do

    !%Variable OCTLaserEnvelope
    !%Type block
    !%Section Optimal Control
    !%Description
    !% Often a predefined time-dependent envelope on the control parameter is desired. 
    !% This can be achieved by making the penalty factor time-dependent. 
    !% Here, you may specify the required time dependent envelope.
    !%
    !% It is possible to choose different envelopes for different control parameters.
    !% There should be one line for each control parameter.
    !%End
        
    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
      no_lines = loct_parse_block_n(blk)
      if(no_lines .ne.par%no_parameters) call input_error('OCTLaserEnvelope')

      ! The structure of the block is:
      ! function | weight 
      do i = 1, no_lines
        call loct_parse_block_string(blk, i-1, 0, expression)  
        call conv_to_C_string(expression)
        call loct_parse_block_float(blk, i-1, 1, weight)
        do j = 1, steps+1
          t = (j-1)*dt
          call loct_parse_expression(f_re, f_im, "t", t, expression)
          call tdf_set_numerical(par%td_penalty(i), j, M_ONE /(f_re + CNST(1.0e-7))  )
        end do
      end do

      call loct_parse_block_end(blk)
    end if

    do i = 1, par%no_parameters
      call tdf_scalar_multiply(octpenalty, par%td_penalty(i))
    end do

    call pop_sub()
  end subroutine parameters_penalty_init


  ! ---------------------------------------------------------
  ! Gets the fluence of the laser field, defined as:
  ! laser_fluence = \sum_{pol} \integrate_0^T
  !                 laserin(t, pol)^2 dt
  ! ---------------------------------------------------------
  FLOAT function laser_fluence(par)
    type(oct_control_parameters_t), intent(in) :: par
    integer :: i, j
    FLOAT :: t
    call push_sub('opt_control_aux.calc_fluence')

    ! WARNING: This is probably very inefficient; there should be functions in
    ! the tdf module taken care of integrating functions.
    laser_fluence = M_ZERO
    do j = 1, par%no_parameters
      do i = 1, par%ntiter+1
        t = (i-1) * par%dt
        laser_fluence = laser_fluence + abs(tdf(par%f(j), i))**2 
      end do
    end do
    laser_fluence = laser_fluence * par%dt

    call pop_sub()
  end function laser_fluence


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function j2_functional(par) result(j2)
    type(oct_control_parameters_t), intent(in) :: par

    integer :: i, j
    FLOAT :: t
    j2 = M_ZERO
    do j = 1, par%no_parameters
      do i = 1, par%ntiter + 1
        t = (i-1) * par%dt
        j2 = j2 + tdf(par%td_penalty(j), i)*abs(tdf(par%f(j), i))**2 
      end do
    end do
    j2 = j2 * par%dt
  end function j2_functional


end module opt_control_parameters_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
