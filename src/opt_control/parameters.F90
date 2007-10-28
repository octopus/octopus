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
  use mesh_m
  use mix_m

  implicit none

  private
  public :: oct_control_parameters_t,     &
            parameters_init,              &
            parameters_set,               &
            parameters_end,               &
            parameters_copy,              &
            parameters_to_h,              &
            parameters_to_h_val,          &
            parameters_write,             &
            parameters_mixing,            &
            laser_fluence,                &
            j2_functional

  type oct_control_parameters_t
    integer :: no_parameters
    FLOAT   :: dt
    integer :: ntiter
    type(tdf_t), pointer :: f(:)
    FLOAT, pointer :: alpha(:)
    type(tdf_t), pointer :: td_penalty(:)
  end type oct_control_parameters_t

  type(mix_t), public :: parameters_mix

contains

  ! ---------------------------------------------------------
  FLOAT function parameters_dotp(x, y) result(res)
    FLOAT, intent(in) :: x(:)
    FLOAT, intent(in) :: y(:)
    
    res = sum(x(:)*y(:))
  end function parameters_dotp


  ! ---------------------------------------------------------
  subroutine parameters_mixing(iter, par_in, par_out, par_new)
    integer, intent(in) :: iter
    type(oct_control_parameters_t), intent(in) :: par_in, par_out
    type(oct_control_parameters_t), intent(inout) :: par_new

    !FLOAT, external :: oct_control_parameters_dotp

    integer :: i, j
    FLOAT, allocatable :: e_in(:, :, :), e_out(:, :, :), e_new(:, :, :)

    ALLOCATE(e_in (par_in%ntiter+1, par_in%no_parameters, 1), (par_in%ntiter+1)*par_in%no_parameters)
    ALLOCATE(e_out(par_in%ntiter+1, par_in%no_parameters, 1), (par_in%ntiter+1)*par_in%no_parameters)
    ALLOCATE(e_new(par_in%ntiter+1, par_in%no_parameters, 1), (par_in%ntiter+1)*par_in%no_parameters)

    do i = 1, par_in%no_parameters
      do j = 1, par_in%ntiter + 1
        e_in (j, i, 1) = tdf(par_in%f(i), j)
        e_out(j, i, 1) = tdf(par_out%f(i), j)
      end do
    end do
    e_new = M_ZERO

    call dmixing(parameters_mix, iter, e_in, e_out, e_new, parameters_dotp)

    do i = 1, par_out%no_parameters
      call tdf_set_numerical(par_new%f(i), e_new(:, i, 1))
    end do

    deallocate(e_in, e_out, e_new)
  end subroutine parameters_mixing


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, no_parameters, dt, ntiter)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer, intent(in) :: no_parameters   
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter
    integer :: j

    call push_sub('opt_control_parameters.parameters_init')

    cp%no_parameters = no_parameters
    cp%dt = dt
    cp%ntiter = ntiter
    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    ALLOCATE(cp%alpha(cp%no_parameters), cp%no_parameters)
    cp%alpha = M_ZERO
    do j = 1, cp%no_parameters
      call tdf_init_numerical(cp%f(j), cp%ntiter, cp%dt)
    end do

    call parameters_penalty_init(cp)

    call pop_sub()
  end subroutine parameters_init


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(epot_t), intent(in) :: ep
    integer :: j

    call push_sub('opt_control_parameters.parameters_set')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
      cp%f(j) = ep%lasers(j)%f
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
      call tdf_end(ep%lasers(j)%f)
      ep%lasers(j)%f = cp%f(j)
    end do
    call pop_sub()
  end subroutine parameters_to_h


  ! ---------------------------------------------------------
  subroutine parameters_to_h_val(cp, ep, val)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep
    integer, intent(in) :: val
    integer :: j
    do j = 1, cp%no_parameters
        call tdf_set_numerical(ep%lasers(j)%f, val, tdf(cp%f(j), val))
    end do
  end subroutine parameters_to_h_val


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j

    call push_sub('opt_control_parameters.parameters_end')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
      call tdf_end(cp%td_penalty(j))
    end do
    deallocate(cp%f)
    nullify(cp%f)
    deallocate(cp%alpha)
    nullify(cp%alpha)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, cp, fourier)
    character(len=*), intent(in) :: filename
    type(oct_control_parameters_t), intent(in) :: cp
    logical, optional, intent(in) :: fourier

    type(tdf_t) :: g
    integer :: i, j, iunit
    FLOAT :: t
    FLOAT, allocatable :: wgrid(:)
    character(len=2) :: digit

    call push_sub('opt_control_parameters.parameters_write')

    call io_mkdir(trim(filename))

    iunit = io_open(trim(filename)//'/Fluence'//digit, action='write')
    write(iunit, '(a,es20.8e3)') 'Flunce = ', laser_fluence(cp)
    call io_close(iunit)

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

    if(present(fourier)) then
    if(fourier) then
      do j = 1, cp%no_parameters
        if(cp%no_parameters > 1) then
          write(digit,'(i2.2)') j
          iunit = io_open(trim(filename)//'/cpw-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cpw', action='write')
        end if
        g = cp%f(j)
        call tdf_fft_forward(g)
        ALLOCATE(wgrid(0:cp%ntiter), cp%ntiter+1)
        call tdf_fourier_grid(g, wgrid)
        do i = 0, cp%ntiter
          write(iunit, '(3es30.16e4)') wgrid(i), tdf(g, i+1)
        end do
        deallocate(wgrid)
        call tdf_end(g)
        call io_close(iunit)
      end do
    end if
    end if

    call pop_sub()
  end subroutine parameters_write


  ! ---------------------------------------------------------
  subroutine parameters_penalty_init(par)
    type(oct_control_parameters_t), intent(inout) :: par

    character(len=1024)      :: expression
    FLOAT                    :: t, octpenalty, dt
    real(8)                  :: f_re, f_im
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
    par%alpha(1:par%no_parameters) = octpenalty

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
    !% There should be one line for each control parameter. Each line should
    !% have only one element: a string with the function that defines the
    !% *inverse* of the time-dependent penalty, which is then defined as
    !% 1 upon this function + 1.0e-7 (to avoid possible singularities).
    !%
    !% The usual choices should be functions between zero and one.
    !%
    !% If, instead of defining a function, the string is "default", then
    !% the program will use the function:
    !%
    !% <math> \frac{1}{\alpha(t)} = \frac{1}{2}( erf(t-T/20)+ erf(-(t-T+T/20)) </math>
    !%End
        
    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then
      no_lines = loct_parse_block_n(blk)
      if(no_lines .ne.par%no_parameters) call input_error('OCTLaserEnvelope')

      do i = 1, no_lines
        call loct_parse_block_string(blk, i-1, 0, expression)
        if(trim(expression)=='default') then
          do j = 1, steps + 1
            t = (j-1)*dt
            f_re = M_HALF * (loct_erf(t-CNST(0.05)*steps*dt) + &
                             loct_erf(-(t-steps*dt+CNST(0.05)*steps*dt)) )
            call tdf_set_numerical(par%td_penalty(i), j, &
              TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        else
          call conv_to_C_string(expression)
          do j = 1, steps+1
            t = (j-1)*dt
            call loct_parse_expression(f_re, f_im, "t", real(t, 8), expression)
            call tdf_set_numerical(par%td_penalty(i), j, TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        end if
      end do

      call loct_parse_block_end(blk)
    end if

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
        laser_fluence = laser_fluence + tdf(par%td_penalty(j), i) * abs(tdf(par%f(j), i))**2 
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
        j2 = j2 + par%alpha(j) * tdf(par%td_penalty(j), i) * abs(tdf(par%f(j), i))**2 
      end do
    end do
    j2 = j2 * par%dt
  end function j2_functional


  ! ---------------------------------------------------------
  subroutine parameters_copy(cp_out, cp_in)
    type(oct_control_parameters_t), intent(inout) :: cp_out
    type(oct_control_parameters_t), intent(in)    :: cp_in
    integer :: j

    cp_out%no_parameters = cp_in%no_parameters
    cp_out%dt = cp_in%dt
    cp_out%ntiter = cp_in%ntiter
    ALLOCATE(cp_out%f(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%alpha(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%td_penalty(cp_out%no_parameters), cp_out%no_parameters)
    do j = 1, cp_in%no_parameters
      cp_out%alpha(j) = cp_in%alpha(j)
      cp_out%f(j) = cp_in%f(j)
      cp_out%td_penalty(j) = cp_in%td_penalty(j)
    end do

  end subroutine parameters_copy

end module opt_control_parameters_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
