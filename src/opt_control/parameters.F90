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
  use loct_parser_m
  use loct_math_m
  use filter_m
  use lasers_m
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
            parameters_diff,              &
            parameters_apply_envelope,    &
            parameters_set_fluence,       &
            parameters_change_rep,        &
            parameters_set_initial,       &
            laser_fluence,                &
            j2_functional

  integer, parameter ::                 &
    ctr_parameter_real_space       = 1, &
    ctr_parameter_frequency_space  = 2


  type oct_control_parameters_t
    integer :: no_parameters
    FLOAT   :: dt
    integer :: ntiter
    FLOAT   :: targetfluence
    type(tdf_t), pointer :: f(:)
    CMPLX, pointer :: pol(:, :) ! the polarization of the field, this is
                                ! necessary to calculate the fluence.
    FLOAT, pointer :: alpha(:)
    type(tdf_t), pointer :: td_penalty(:)

    integer :: representation
    integer :: current_representation
    FLOAT   :: omegamax
  end type oct_control_parameters_t

  type(mix_t), public :: parameters_mix

contains


  ! ---------------------------------------------------------
  subroutine parameters_set_initial(par, ep, m, dt, max_iter)
    type(oct_control_parameters_t), intent(inout) :: par
    type(epot_t), intent(inout)                   :: ep
    type(mesh_t), intent(inout)                   :: m
    FLOAT, intent(in)                             :: dt
    integer, intent(in)                           :: max_iter

    integer :: i
    FLOAT :: targetfluence, omegamax
    logical :: fix_initial_fluence

    !%Variable OCTFixFluenceTo
    !%Type float
    !%Section Optimal Control
    !%Default 0.0
    !%Description
    !% The algorithm tries to obtain the specified fluence for the laser field. 
    !% This works only in conjunction with either the WG05 or the straight iteration scheme.
    !%
    !% If this variable is not present in the input file, by default the code will not
    !% attempt a fixed-fluence QOCT run. The same holds if the value given to this
    !% variable is exactly zero.
    !%
    !% If this variable is given a negative value, then the target fluence will be that of
    !% the initial laser pulse given as guess in the input file. Note, however, that
    !% first the code applies the envelope provided by the "OCTLaserEnvelope" input
    !% option, and afterwards it calculates the fluence.
    !%End
    call loct_parse_float(check_inp('OCTFixFluenceTo'), M_ZERO, targetfluence)


    !%Variable OCTParameterOmegaMax
    !%Type float
    !%Section Optimal Control
    !%Default -1.0
    !%Description
    !%
    !%End
    call loct_parse_float(check_inp('OCTParameterOmegaMax'), -M_ONE, omegamax)

    !%Variable OCTParameterRepresentation
    !%Type integer
    !%Section Optimal Control
    !%Default control_parameter_real_space
    !%Description
    !%
    !%Option control_parameters_real_space 1
    !%
    !%Option control_parameters_fourier_space 2
    !%
    !%End
    call loct_parse_int(check_inp('OCTParameterRepresentation'), ctr_parameter_real_space, par%representation)
    if(.not.varinfo_valid_option('OCTParameterRepresentation', par%representation)) &
      call input_error('OCTParameterRepresentation')

    !%Variable OCTFixInitialFluence
    !%Type logical
    !%Section Optimal Control
    !%Default yes
    !%Description
    !% By default, when asking for a fixed-fluence optimization ("OCTFixFluenceTo = whatever"), 
    !% the initial laser guess provided in the input file is scaled to match this
    !% fluence. However, you can force the program to use that initial laser as the initial
    !% guess, no matter the fluence, by setting "OCTFixInitialFluence = no".
    !%End
    call loct_parse_logical(check_inp('OCTFixInitialFluence'), .true., fix_initial_fluence)


    ! Initial guess for the laser: read from the input file.
    call laser_init(ep%no_lasers, ep%lasers, m)

    ! Check that the laser polarizations are not imaginary.
    do i = 1, ep%no_lasers
      if(any(aimag(ep%lasers(i)%pol(:)) .ne. M_ZERO)) then
        write(message(1),'(a)') 'For optimal control runs, you cannot specify an initial field guess'
        write(message(2),'(a)') 'with complex polarization direction'
        call write_fatal(2)
      end if
    end do

    do i = 1, ep%no_lasers
      call laser_to_numerical(ep%lasers(i), dt, max_iter)
    end do

    call parameters_init(par, ep%no_lasers, dt, max_iter, targetfluence, omegamax)
    call parameters_set(par, ep)
    call parameters_apply_envelope(par)

    ! Moving to the sine-Fourier space.
    if(par%representation .eq. ctr_parameter_frequency_space) then
      message(1) = "Info: The control function will be represented as a sine-Fourier series."
      call write_info(1)
      call parameters_change_rep(par)
    end if

    if(par%targetfluence .ne. M_ZERO) then
      if(fix_initial_fluence) then
        if(par%targetfluence < M_ZERO) then
          par%targetfluence = laser_fluence(par) 
          write(0, *) 'LASER_FLUENCE 2', laser_fluence(par)
        else
          call parameters_set_fluence(par, par%targetfluence)
        end if
      end if
    end if

    if(par%representation .eq. ctr_parameter_frequency_space) then
      call parameters_change_rep(par)
    end if

  end subroutine parameters_set_initial
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_change_rep(par)
    type(oct_control_parameters_t), intent(inout) :: par

    integer :: j

    if(par%current_representation .eq. ctr_parameter_real_space) then
      do j = 1, par%no_parameters
        call tdf_numerical_to_sineseries(par%f(j), par%omegamax)
      end do
      par%current_representation = ctr_parameter_frequency_space
    else
      do j = 1, par%no_parameters
        call tdf_sineseries_to_numerical(par%f(j))
      end do
      par%current_representation = ctr_parameter_real_space
    end if

  end subroutine parameters_change_rep
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_diff(p, q) result(res)
    type(oct_control_parameters_t), intent(in) :: p, q
    integer :: i, j
    FLOAT, allocatable :: delta(:)

    res = M_ZERO
    ALLOCATE(delta(p%ntiter + 1), p%ntiter +1)
    do i = 1, p%no_parameters
      do j = 1, p%ntiter + 1
        delta(j) = abs(tdf(p%f(i), j) - tdf(q%f(i), j))**2
      end do
      res = res + sum(delta)*p%dt
    end do

    deallocate(delta)
  end function parameters_diff


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

    integer :: i, j
    FLOAT, allocatable :: e_in(:, :, :), e_out(:, :, :), e_new(:, :, :)
    call push_sub('parameters.parameters_mixing')


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
    call pop_sub()
  end subroutine parameters_mixing


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, no_parameters, dt, ntiter, targetfluence, omegamax)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer, intent(in) :: no_parameters   
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter
    FLOAT, intent(in) :: targetfluence
    FLOAT, intent(in) :: omegamax

    integer :: j

    call push_sub('parameters.parameters_init')

    cp%current_representation = ctr_parameter_real_space
    cp%omegamax = omegamax

    cp%no_parameters = no_parameters
    cp%dt = dt
    cp%ntiter = ntiter
    cp%targetfluence = targetfluence
    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    ALLOCATE(cp%alpha(cp%no_parameters), cp%no_parameters)
    ALLOCATE(cp%pol(MAX_DIM, cp%no_parameters), MAX_DIM*cp%no_parameters)
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

    call push_sub('parameters.parameters_set')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
      cp%f(j) = ep%lasers(j)%f
      cp%pol(1:MAX_DIM, j) = ep%lasers(j)%pol(1:MAX_DIM)
    end do

    call pop_sub()
  end subroutine parameters_set


  ! ---------------------------------------------------------
  subroutine parameters_apply_envelope(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j, i

    call push_sub('parameters.parameters_apply_envelope')

    ! Do not apply the envelope if the parameters are represented as a sine Fourier series.
    if(cp%representation .eq. ctr_parameter_real_space) then
      do j = 1, cp%no_parameters
        do i = 1, cp%ntiter + 1
          call tdf_set_numerical(cp%f(j), i, tdf(cp%f(j), i) / tdf(cp%td_penalty(j), i) )
        end do
      end do
    end if

    call pop_sub()
  end subroutine parameters_apply_envelope


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep
    integer :: j
    call push_sub('parameters.parameters_to_h')
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

    call push_sub('parameters.parameters_end')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
      call tdf_end(cp%td_penalty(j))
    end do
    deallocate(cp%f)
    nullify(cp%f)
    deallocate(cp%td_penalty)
    nullify(cp%td_penalty)
    deallocate(cp%alpha)
    nullify(cp%alpha)
    deallocate(cp%pol)
    nullify(cp%pol)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, cp, fourier)
    character(len=*), intent(in) :: filename
    type(oct_control_parameters_t), intent(in), target :: cp
    logical, optional, intent(in) :: fourier

    type(tdf_t) :: g
    integer :: i, j, iunit
    logical :: change_rep
    FLOAT :: t
    FLOAT, allocatable :: wgrid(:)
    character(len=2) :: digit
    type(oct_control_parameters_t), pointer :: par

    call push_sub('parameters.parameters_write')

    call io_mkdir(trim(filename))

    ! First, check that we are in the real space representation.
    if(cp%current_representation .eq. ctr_parameter_frequency_space) then
      ALLOCATE(par, 1)
      call parameters_copy(par, cp)
      call parameters_change_rep(par)
      change_rep = .true.
    else
      par => cp
      change_rep = .false.
    end if

    iunit = io_open(trim(filename)//'/Fluence', action='write')
    write(iunit, '(a,es20.8e3)') 'Fluence = ', laser_fluence(par)
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
        write(iunit, '(3es20.8e3)') t, tdf(par%f(j), t)
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
        call tdf_init(g)
        call tdf_copy(g, par%f(j))
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

    if(change_rep) then 
      call parameters_end(par)
    else
      nullify(par)
    end if

    call pop_sub()
  end subroutine parameters_write


  ! ---------------------------------------------------------
  subroutine parameters_penalty_init(par)
    type(oct_control_parameters_t), intent(inout) :: par

    character(len=1024)      :: expression
    FLOAT                    :: t, octpenalty, dt
    real(8)                  :: f_re, f_im
    type(block_t)            :: blk
    integer                  :: no_lines, i, j, steps
    call push_sub('parameters.parameters_penalty_init')

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
  ! laser_fluence = \sum_i^{no_parameters} \integrate_0^T |epsilon(t)|^2
  ! ---------------------------------------------------------
  FLOAT function laser_fluence(par)
    type(oct_control_parameters_t), intent(in) :: par
    integer :: i, j, k, nfreqs
    FLOAT :: t
    call push_sub('parameters.laser_fluence')

    laser_fluence = M_ZERO
    do j = 1, par%no_parameters
      laser_fluence = laser_fluence + tdf_dot_product(par%f(j), par%f(j))
    end do

    call pop_sub()
  end function laser_fluence


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function j2_functional(par) result(j2)
    type(oct_control_parameters_t), intent(in) :: par
    integer :: i, j, k

    FLOAT :: t
    j2 = M_ZERO
    do j = 1, par%no_parameters
      do i = 1, par%ntiter + 1
        t = (i-1) * par%dt
        do k = 1, MAX_DIM
          j2 = j2 + tdf(par%td_penalty(j), i) * real( par%pol(k, j) * tdf(par%f(j), i) )**2 
        end do
      end do
    end do
    j2 = - par%alpha(1) * (j2 * par%dt - par%targetfluence)
  end function j2_functional
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_set_fluence(par, fluence)
    type(oct_control_parameters_t), intent(inout) :: par
    FLOAT, intent(in)                             :: fluence

    FLOAT   :: old_fluence
    integer :: j

    old_fluence = laser_fluence(par) 
    do j = 1, par%no_parameters
      call tdf_scalar_multiply( sqrt(fluence/old_fluence) , par%f(j) )
    end do

  end subroutine parameters_set_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_copy(cp_out, cp_in)
    type(oct_control_parameters_t), intent(inout) :: cp_out
    type(oct_control_parameters_t), intent(in)    :: cp_in
    integer :: j

    cp_out%targetfluence = cp_in%targetfluence
    cp_out%no_parameters = cp_in%no_parameters
    cp_out%dt = cp_in%dt
    cp_out%ntiter = cp_in%ntiter
    cp_out%targetfluence = cp_in%targetfluence
    cp_out%representation = cp_in%representation
    cp_out%current_representation = cp_in%current_representation
    cp_out%omegamax = cp_in%omegamax
    ALLOCATE(cp_out%f(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%alpha(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%td_penalty(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%pol(MAX_DIM, cp_out%no_parameters), MAX_DIM*cp_out%no_parameters)
    do j = 1, cp_in%no_parameters
      cp_out%alpha(j) = cp_in%alpha(j)
      call tdf_init(cp_out%f(j))
      call tdf_copy(cp_out%f(j), cp_in%f(j))
      call tdf_init(cp_out%td_penalty(j))
      call tdf_copy(cp_out%td_penalty(j), cp_in%td_penalty(j))
      cp_out%pol(1:MAX_DIM, j) = cp_in%pol(1:MAX_DIM, j)
    end do

  end subroutine parameters_copy

end module opt_control_parameters_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
