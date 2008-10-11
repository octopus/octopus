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
  use lalg_adv_m
  use string_m
  use datasets_m
  use varinfo_m
  use global_m
  use messages_m
  use io_m
  use units_m
  use loct_parser_m
  use loct_math_m
  use filter_m
  use lasers_m
  use external_pot_m
  use tdf_m
  use mesh_m
  use mix_m
  use filter_m

  implicit none

  private
  public :: oct_control_parameters_t,     &
            parameters_read,              &
            parameters_init,              &
            parameters_set,               &
            parameters_end,               &
            parameters_copy,              &
            parameters_to_h,              &
            parameters_to_h_val,          &
            parameters_write,             &
            parameters_mixing,            &
            parameters_mixing_init,       &
            parameters_mixing_end,        &
            parameters_diff,              &
            parameters_apply_envelope,    &
            parameters_set_fluence,       &
            parameters_set_alpha,         &
            parameters_set_rep,           &
            parameters_to_realtime,       &
            parameters_prepare_initial,   &
            parameters_fluence,           &
            parameters_j2,                &
            parameters_par_to_x,          &
            parameters_x_to_par,          &
            parameters_update,            &
            parameters_number,            &
            parameters_nfreqs,            &
            parameters_w0,                &
            parameters_alpha,             &
            parameters_targetfluence,     &
            parameters_filter


  integer, parameter ::                 &
    ctr_parameter_real_space       = 1, &
    ctr_parameter_frequency_space  = 2


  integer, parameter :: parameter_mode_none      = 0, &
                        parameter_mode_epsilon   = 1, &
                        parameter_mode_f         = 2, &
                        parameter_mode_phi       = 3, &
                        parameter_mode_f_and_phi = 4

  type oct_parameters_common_t
    private
    integer :: representation      = 0
    FLOAT   :: omegamax            = M_ZERO
    FLOAT   :: targetfluence       = M_ZERO
    logical :: fix_initial_fluence = .false.
    FLOAT   :: w0                  = M_ZERO
    integer :: mode                = parameter_mode_none
    integer :: no_parameters       = 0
    FLOAT,       pointer :: alpha(:)      => NULL()
    type(tdf_t), pointer :: td_penalty(:) => NULL()

    CMPLX, pointer :: pol(:, :)           => NULL() ! the polarization of the field, this is
                                                    ! necessary to calculate the fluence.
    type(tdf_t)    :: f ! This is the envelope of the laser field, only used in the phase-only
                        ! optimization (necessary to store it in order to calculate fluences)
  end type oct_parameters_common_t


  type oct_control_parameters_t
    private
    integer :: no_parameters = 0
    FLOAT   :: dt            = M_ZERO
    integer :: ntiter        = 0
    integer :: nfreqs        = 0
    FLOAT   :: targetfluence = M_ZERO
    FLOAT   :: intphi        = M_ZERO
    type(tdf_t), pointer :: f(:) => NULL()
    FLOAT, pointer :: alpha(:)   => NULL()

    integer :: representation         = 0
    integer :: current_representation = 0
    FLOAT   :: omegamax               = M_ZERO

    FLOAT   :: w0       = M_ZERO
    FLOAT, pointer :: utransf(:, :)  => NULL()
    FLOAT, pointer :: utransfi(:, :) => NULL()
  end type oct_control_parameters_t

  type(oct_parameters_common_t), save :: par_common
  type(mix_t) :: parameters_mix

contains


  ! ---------------------------------------------------------
  subroutine parameters_read(ep, dt, max_iter, mode_fixed_fluence, mode_basis_set)
    type(epot_t), intent(inout)                   :: ep
    FLOAT, intent(in)                             :: dt
    integer, intent(in)                           :: max_iter
    logical, intent(out)                          :: mode_fixed_fluence
    logical, intent(out)                          :: mode_basis_set

    character(len=1024)      :: expression
    integer :: i, j, no_lines, steps
    FLOAT   :: octpenalty, t, f_re, f_im
    type(block_t)            :: blk

    call push_sub('parameters.parameters_read')

    !%Variable OCTParameterRepresentation
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default control_parameter_real_space
    !%Description
    !% The control functions can be represented in real space (default), or expanded
    !% in a finite basis set (now, the only basis set defined in the code for this is
    !% a sine Fourier series). 
    !%Option control_parameters_real_space 1
    !%
    !%Option control_parameters_fourier_space 2
    !%
    !%End
    call loct_parse_int(check_inp('OCTParameterRepresentation'), &
      ctr_parameter_real_space, par_common%representation)
    if(.not.varinfo_valid_option('OCTParameterRepresentation', par_common%representation)) &
      call input_error('OCTParameterRepresentation')
    select case(par_common%representation)
    case(ctr_parameter_real_space)
      write(message(1), '(a)') 'Info: The OCT control functions will be represented in real time.'
      call write_info(1)
    case(ctr_parameter_frequency_space)
      write(message(1), '(a)') 'Info: The OCT control functions will be represented as a sine '
      write(message(2), '(a)') '      Fourier series.'
      call write_info(2)
    end select

    !%Variable OCTParameterOmegaMax
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 1.0
    !%Description
    !%
    !%End
    call loct_parse_float(check_inp('OCTParameterOmegaMax'), M_ONE / units_inp%energy%factor, par_common%omegamax)
    if(par_common%representation .eq. ctr_parameter_frequency_space) then
      write(message(1), '(a)')         'Info: The representation of the OCT control parameters as a sine'
      write(message(2), '(a,f10.5,a)') '      Fourier series will be done with a cut-off of ', &
        par_common%omegamax / units_out%energy%factor, ' ['//trim(units_out%energy%abbrev) // ']'
      call write_info(2)
    end if

    !%Variable OCTFixFluenceTo
    !%Type float
    !%Section Calculation Modes::Optimal Control
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
    call loct_parse_float(check_inp('OCTFixFluenceTo'), M_ZERO, par_common%targetfluence)

    !%Variable OCTFixInitialFluence
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default yes
    !%Description
    !% By default, when asking for a fixed-fluence optimization ("OCTFixFluenceTo = whatever"), 
    !% the initial laser guess provided in the input file is scaled to match this
    !% fluence. However, you can force the program to use that initial laser as the initial
    !% guess, no matter the fluence, by setting "OCTFixInitialFluence = no".
    !%End
    call loct_parse_logical(check_inp('OCTFixInitialFluence'), .true., par_common%fix_initial_fluence)

    !%Variable OCTControlFunctionType
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 1
    !%Description
    !% 
    !%Option parameter_mode_epsilon   1
    !%
    !%Option parameter_mode_f         2
    !%
    !%Option parameter_mode_phi       3
    !%
    !%Option parameter_mode_f_and_phi 4
    !%
    !%End
    call loct_parse_int(check_inp('OCTControlFunctionType'), parameter_mode_epsilon, par_common%mode)
    if(.not.varinfo_valid_option('OCTControlFunctionType', par_common%mode)) call input_error('OCTControlFunctionType')
     select case(par_common%mode)
    case(parameter_mode_f, parameter_mode_phi, parameter_mode_f_and_phi)
      if(par_common%representation .ne. ctr_parameter_frequency_space) then
        write(message(1),'(a)') 'If "OCTControlFunctionType = parameter_mode_f", '
        write(message(2),'(a)') '   "OCTControlFunctionType = parameter_mode_phi", '
        write(message(3),'(a)') 'or "OCTControlFunctionType = parameter_mode_f_and_phi", then'
        write(message(4),'(a)') 'the control parameter representation cannot be real time, i.e. you must have'
        write(message(5),'(a)') '"OCTParameterRepresentation = control_parameters_fourier_space".'
        call write_fatal(5)
      end if
     end select

    ! Check that the laser polarizations are not imaginary.
    ALLOCATE(par_common%pol(MAX_DIM, ep%no_lasers), MAX_DIM*ep%no_lasers)
    do i = 1, ep%no_lasers

      par_common%pol(1:MAX_DIM, i) = laser_polarization(ep%lasers(i))

      ! WARNING: Check if this is working....
      if(any(aimag(laser_polarization(ep%lasers(i))) .ne. M_ZERO)) then
        write(message(1),'(a)') 'For optimal control runs, you cannot specify an initial field guess'
        write(message(2),'(a)') 'with complex polarization direction'
        call write_fatal(2)
      end if
    end do

    ! Note that in QOCT runs, it is not acceptable to have complex time-dependent functions.
    do i = 1, ep%no_lasers
      select case(par_common%mode)
      case(parameter_mode_epsilon)
        call laser_to_numerical_all(ep%lasers(i), dt, max_iter, real_part = .true.)
      case default
        call laser_to_numerical(ep%lasers(i), dt, max_iter, real_part = .true.)
      end select
    end do

    ! For phase-only optimization, we need to store the envelope, in order to be able
    ! to calculate the fluence.
    if(par_common%mode .eq. parameter_mode_phi) then
      call laser_get_f(ep%lasers(1), par_common%f)
    end if

    ! Fix the carrier frequency
    call obsolete_variable('OCTCarrierFrequency')
    par_common%w0 = laser_carrier_frequency(ep%lasers(1))

    ! Fix the number of parameters
    select case(par_common%mode)
      case(parameter_mode_epsilon)
        par_common%no_parameters = ep%no_lasers
      case(parameter_mode_f)
        par_common%no_parameters = 1
      case(parameter_mode_phi)
        par_common%no_parameters = 1
      case(parameter_mode_f_and_phi)
        par_common%no_parameters = 2
    end select

    mode_fixed_fluence = .false.
    if (par_common%targetfluence .ne. M_ZERO) mode_fixed_fluence = .true.

    mode_basis_set = .false.
    if(par_common%representation .ne. ctr_parameter_real_space ) mode_basis_set = .true.


    !%Variable OCTPenalty
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 1.0
    !%Description
    !% The variable specificies the value of the penalty factor for the 
    !% integrated field strength (fluence). Large value - small fluence.
    !% A transient shape can be specified using the block OCTLaserEnvelope.
    !% In this case OCTPenalty is multiplied with time-dependent function. 
    !% The value depends on the coupling between the states. A good start might be a 
    !% value from 0.1 (strong fields) to 10 (weak fields). 
    !%
    !% Note that if there are several control functions, one can specify this
    !% variable as a one-line code, each column being the penalty factor for each
    !% of the control functions. Make sure that the number of columns is equal to the
    !% number of control functions. If it is not a block, all control functions will
    !% have the same penalty factor. 
    !%
    !% All penalty factors must be positive. 
    !%End
    ALLOCATE(par_common%alpha(par_common%no_parameters), par_common%no_parameters)
    par_common%alpha = M_ZERO
    if(loct_parse_block('OCTPenalty', blk) == 0) then
      ! We have a block
      i = loct_parse_block_cols(blk, 0)
      if(i.ne.par_common%no_parameters) then
        call input_error('OCTPenalty')
      else
        do j = 1, i
          call loct_parse_block_float(blk, 0, j-1, par_common%alpha(j))
          if(par_common%alpha(j) <= M_ZERO) call input_error('OCTPenalty')
        end do
      end if
    else
      ! We have the same penalty for all the control functions.
      call loct_parse_float(check_inp('OCTPenalty'), M_ONE, octpenalty)
      par_common%alpha(1:par_common%no_parameters) = octpenalty
    end if


    !%Variable OCTLaserEnvelope
    !%Type block
    !%Section Calculation Modes::Optimal Control
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
    steps = max_iter
    ALLOCATE(par_common%td_penalty(par_common%no_parameters), par_common%no_parameters)
    do i = 1, par_common%no_parameters
      call tdf_init_numerical(par_common%td_penalty(i), steps, dt, initval = M_z1)
    end do

    if (loct_parse_block(check_inp('OCTLaserEnvelope'), blk)==0) then

      ! Cannot have this unless we have the "usual" parameter_mode_epsilon.
      if(par_common%mode .ne. parameter_mode_epsilon) then
        write(message(1),'(a)') 'The block "OCTLaserEnvelope" is only compatible with the option'
        write(message(2),'(a)') '"OCTControlFunctionType = parameter_mode_epsilon".'
        call write_fatal(2)
      end if

      no_lines = loct_parse_block_n(blk)
      if(no_lines .ne.par_common%no_parameters) call input_error('OCTLaserEnvelope')

      do i = 1, no_lines
        call loct_parse_block_string(blk, i-1, 0, expression)
        if(trim(expression)=='default') then
          do j = 1, steps + 1
            t = (j-1)*dt
            f_re = M_HALF * (loct_erf(t-CNST(0.05)*steps*dt) + &
                             loct_erf(-(t-steps*dt+CNST(0.05)*steps*dt)) )
            call tdf_set_numerical(par_common%td_penalty(i), j, &
              TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        else
          call conv_to_C_string(expression)
          do j = 1, steps+1
            t = (j-1)*dt
            call loct_parse_expression(f_re, f_im, "t", real(t, 8), expression)
            call tdf_set_numerical(par_common%td_penalty(i), j, TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        end if
      end do

      call loct_parse_block_end(blk)
    end if

    call pop_sub()
  end subroutine parameters_read
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_prepare_initial(par)
    type(oct_control_parameters_t), intent(inout) :: par

    integer :: i, mm, nn
    FLOAT :: t, det
    type(tdf_t) :: fn, fm
    FLOAT, allocatable :: eigenvec(:, :), eigenval(:)

    call push_sub('parameters.parameters_set_initial')

    call parameters_apply_envelope(par)

    ! Move to the sine-Fourier space if required.
    call parameters_set_rep(par)
    if(par%representation .eq. ctr_parameter_frequency_space) then
      if(par%nfreqs <= 1) then
        write(message(1), '(a)')    'Error: The dimension of the basis set used to represent the control'
        write(message(2), '(a)')    '       functions must be larger than one. The input options that you'
        write(message(3), '(a)')    '       supply do not meet this criterion.'
        call write_fatal(3)
      end if
      write(message(1), '(a)')      'Info: The expansion of the control parameters in a sine Fourier series'
      write(message(2), '(a,i6,a)') '      expansion implies the use of ', par%nfreqs, ' basis set functions.'
      call write_info(2)
    end if

    if(par%targetfluence .ne. M_ZERO) then
      if(par%targetfluence < M_ZERO) then
        par%targetfluence = parameters_fluence(par) 
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with the same'
        write(message(2), '(a,f10.5,a)') '      fluence as the input external fields: F = ', par%targetfluence, ' a.u.'
      else
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with a predefined'
        write(message(2), '(a,f10.5,a)') '      fluence: F = ', par%targetfluence, ' a.u.2'
      end if
      call write_info(2)
      if(par_common%fix_initial_fluence) call parameters_set_fluence(par)
    end if

    ! Now we ave to find the "fluence" of the phase, in order to keep it constant.
    if(par_common%mode .eq. parameter_mode_phi) then
      par%intphi = tdf_dot_product(par%f(1), par%f(1))
      if(par%intphi <= M_ZERO) then
        write(message(1), '(a)') 'You must supply a non-null initial-guess phase.'
        call write_fatal(1)
      end if
    end if

    if(par_common%mode .eq. parameter_mode_f) then
      ! If par%envelope = .true., the object to optimize is the envelope of the
      ! the laser pulse. Being e(t) the laser pulse, it is assumed that it
      ! has the form:
      !   e(t) = f(t) cos(w0*t),
      ! where f(t) is the envelope. This is then expanded in a basis set:
      !   f(t) = sum_{n=1}^N f_n g_n(t).
      ! The fluence F[e] is then given by:
      !   F[e] = sum_{m=1}^N sum_{n=1}^N f_m f_n S_{nm},
      !   S_{nm} = \int_{0}^{T} dt g_n(t) g_m(t) cos^2(w0*t).
      ! The following lines of code calculate a matrix U, placed in par%utransf,
      ! that performs the transformation of variables that takes {f_n} to
      ! {h_n}, (\vec{h} = U\vec{f}), such that
      !   F[e] = sum-{m=1}^N h_n^2.
      ! The inverse matrix U^{-1} is placed in par%utransfi.
      !
      ! This scan probably be optimized in some way?
      ALLOCATE(eigenvec(par%nfreqs, par%nfreqs), par%nfreqs*par%nfreqs)
      ALLOCATE(eigenval(par%nfreqs), par%nfreqs)
      ALLOCATE(par%utransf (par%nfreqs, par%nfreqs), par%nfreqs*par%nfreqs)
      ALLOCATE(par%utransfi(par%nfreqs, par%nfreqs), par%nfreqs*par%nfreqs)
      par%utransf  = M_ZERO
      par%utransfi = M_ZERO

      do mm = 1, par%nfreqs
        do nn = mm, par%nfreqs
          call tdf_init_numerical(fn, par%ntiter, par%dt, initval = M_z0, omegamax = par%omegamax)
          call tdf_init_numerical(fm, par%ntiter, par%dt, initval = M_z0, omegamax = par%omegamax)
          call tdf_numerical_to_sineseries(fn, par%omegamax)
          call tdf_numerical_to_sineseries(fm, par%omegamax)
          call tdf_set_numerical(fm, mm, M_ONE)
          call tdf_set_numerical(fn, nn, M_ONE)

          call tdf_sineseries_to_numerical(fm)
          call tdf_sineseries_to_numerical(fn)

          do i = 1, par%ntiter + 1
            t = (i-1)*par%dt
            call tdf_set_numerical(fm, i, tdf(fm, i)*cos(par%w0*t))
            call tdf_set_numerical(fn, i, tdf(fn, i)*cos(par%w0*t))
          end do
          par%utransf(mm, nn) = tdf_dot_product(fm, fn)
          call tdf_numerical_to_sineseries(fn, par%omegamax)
          call tdf_numerical_to_sineseries(fm, par%omegamax)

          call tdf_end(fm)
          call tdf_end(fn)
        end do
      end do
      do mm = 1, par%nfreqs
        do nn = 1, mm - 1
          par%utransf(mm, nn) = par%utransf(nn, mm)
        end do
      end do

      call lalg_eigensolve(par%nfreqs, par%utransf, eigenvec, eigenval)
      do mm = 1, par%nfreqs
        do nn = 1, par%nfreqs
          eigenvec(mm, nn) = eigenvec(mm, nn) * sqrt(eigenval(nn))
        end do
      end do

      par%utransf = transpose(eigenvec)
      par%utransfi = par%utransf
      det =  lalg_inverter(par%nfreqs, par%utransfi)

    end if

    call pop_sub()
  end subroutine parameters_prepare_initial
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Transforms the control function to frequency space, if
  ! this is the space in which the functions are defined (and it
  ! is necessary to perform the transformation). 
  ! And, transforms the control function to real-time space, if
  ! this is the space in which the functions are defined (and it
  ! is necessary to perform the transformation). 
  ! ---------------------------------------------------------
  subroutine parameters_set_rep(par)
    type(oct_control_parameters_t), intent(inout) :: par
    integer :: j

    call push_sub('parameters.parameters_set_rep')

    if(par%current_representation .eq. par%representation) return

    select case(par%current_representation)
    case(ctr_parameter_frequency_space)
      do j = 1, par%no_parameters
        call tdf_sineseries_to_numerical(par%f(j))
      end do
    case(ctr_parameter_real_space)
      do j = 1, par%no_parameters
        call tdf_numerical_to_sineseries(par%f(j), par%omegamax)
      end do
    end select
    par%current_representation = par%representation

    call pop_sub()
  end subroutine parameters_set_rep
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_realtime(par)
    type(oct_control_parameters_t), intent(inout) :: par
    integer :: j
    call push_sub('parameters.parameters_to_realtime')

    if(par%current_representation .eq. ctr_parameter_real_space) return

    do j = 1, par%no_parameters
      call tdf_sineseries_to_numerical(par%f(j))
    end do
    par%current_representation = ctr_parameter_real_space

    call pop_sub()
  end subroutine parameters_to_realtime
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_diff(p, q) result(res)
    type(oct_control_parameters_t), intent(in) :: p, q
    integer :: i, j
    FLOAT, allocatable :: delta(:)
    integer :: nfreqs

    call push_sub('parameters.parameters_diff')

    ASSERT(p%current_representation .eq. q%current_representation)

    select case(p%current_representation)
    case(ctr_parameter_real_space)
      res = M_ZERO
      ALLOCATE(delta(p%ntiter + 1), p%ntiter +1)
      do i = 1, p%no_parameters
        do j = 1, p%ntiter + 1
          delta(j) = abs(tdf(p%f(i), j) - tdf(q%f(i), j))**2
        end do
        res = res + sum(delta)*p%dt
      end do

    case(ctr_parameter_frequency_space)
      nfreqs = tdf_nfreqs(p%f(1))
      res = M_ZERO
      ALLOCATE(delta(nfreqs), nfreqs +1)
      do i = 1, p%no_parameters
        do j = 1, nfreqs
          delta(j) = abs(tdf(p%f(i), j) - tdf(q%f(i), j))**2
        end do
        res = res + sum(delta)
      end do

    end select

    deallocate(delta)
    call pop_sub()
  end function parameters_diff
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_dotp(x, y) result(res)
    FLOAT, intent(in) :: x(:)
    FLOAT, intent(in) :: y(:)
    
    res = sum(x(:)*y(:))
  end function parameters_dotp
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_mixing_init(par)
    type(oct_control_parameters_t), intent(in) :: par

    integer :: nfreqs

    call push_sub('parameters.parameters_mixing_init')

    select case(par%representation)
    case(ctr_parameter_real_space)
      call mix_init(parameters_mix, par%ntiter + 1, par%no_parameters, 1)
    case(ctr_parameter_frequency_space)
      nfreqs = tdf_nfreqs(par%f(1))
      call mix_init(parameters_mix, nfreqs, par%no_parameters, 1)
    end select

    call pop_sub()
  end subroutine parameters_mixing_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_mixing_end
    call mix_end(parameters_mix)
  end subroutine parameters_mixing_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_mixing(iter, par_in, par_out, par_new)
    integer, intent(in) :: iter
    type(oct_control_parameters_t), intent(in) :: par_in, par_out
    type(oct_control_parameters_t), intent(inout) :: par_new

    integer :: i, j, nfreqs
    FLOAT, allocatable :: e_in(:, :, :), e_out(:, :, :), e_new(:, :, :)
    call push_sub('parameters.parameters_mixing')

    ! First, some sanity checks:
    ASSERT(par_in%representation .eq. par_out%representation)
    ASSERT(par_in%representation .eq. par_new%representation)
    ASSERT(par_in%representation .eq. par_in%current_representation)
    ASSERT(par_out%representation .eq. par_out%current_representation)
    ASSERT(par_new%representation .eq. par_new%current_representation)


    select case(par_in%representation)
    case(ctr_parameter_real_space)

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

    case(ctr_parameter_frequency_space)

      nfreqs = tdf_nfreqs(par_in%f(1))
      ALLOCATE(e_in (nfreqs, par_in%no_parameters, 1), nfreqs*par_in%no_parameters)
      ALLOCATE(e_out(nfreqs, par_in%no_parameters, 1), nfreqs*par_in%no_parameters)
      ALLOCATE(e_new(nfreqs, par_in%no_parameters, 1), nfreqs*par_in%no_parameters)
      do i = 1, par_in%no_parameters
        do j = 1, nfreqs
          e_in (j, i, 1) = tdf(par_in%f(i), j)
          e_out(j, i, 1) = tdf(par_out%f(i), j)
        end do
      end do
      e_new = M_ZERO
      call dmixing(parameters_mix, iter, e_in, e_out, e_new, parameters_dotp)
      do i = 1, par_out%no_parameters
        call tdf_set_numerical(par_new%f(i), e_new(:, i, 1))
      end do

    end select

    deallocate(e_in, e_out, e_new)
    call pop_sub()
  end subroutine parameters_mixing
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, dt, ntiter)
    type(oct_control_parameters_t), intent(inout) :: cp
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter

    integer :: j

    call push_sub('parameters.parameters_init')


    cp%representation  = par_common%representation
    cp%w0              = par_common%w0
    cp%omegamax        = par_common%omegamax
    cp%no_parameters   = par_common%no_parameters
    cp%dt              = dt
    cp%ntiter          = ntiter

    cp%current_representation = ctr_parameter_real_space

    cp%targetfluence = par_common%targetfluence
    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    ALLOCATE(cp%alpha(cp%no_parameters), cp%no_parameters)
    cp%alpha = M_ZERO
    do j = 1, cp%no_parameters
      call tdf_init_numerical(cp%f(j), cp%ntiter, cp%dt, omegamax = par_common%omegamax)
    end do

    cp%alpha = par_common%alpha

    cp%nfreqs = 0
    if(cp%representation .eq. ctr_parameter_frequency_space) then
      cp%nfreqs = tdf_nfreqs(cp%f(1))
    end if

    call pop_sub()
  end subroutine parameters_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(epot_t), intent(in) :: ep
    integer :: j

    call push_sub('parameters.parameters_set')

    select case(par_common%mode)
    case(parameter_mode_epsilon, parameter_mode_f)
      do j = 1, cp%no_parameters
        call tdf_end(cp%f(j))
        call laser_get_f(ep%lasers(j), cp%f(j))
      end do
    case(parameter_mode_phi)
      call tdf_end(cp%f(1))
      call laser_get_phi(ep%lasers(1), cp%f(1))
    case(parameter_mode_f_and_phi)
      call tdf_end(cp%f(1))
      call laser_get_f(ep%lasers(1), cp%f(1))
      call tdf_end(cp%f(2))
      call laser_get_phi(ep%lasers(1), cp%f(2))
    end select

    call pop_sub()
  end subroutine parameters_set
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_apply_envelope(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j, i

    call push_sub('parameters.parameters_apply_envelope')

    ! Do not apply the envelope if the parameters are represented as a sine Fourier series.
    if(cp%representation .eq. ctr_parameter_real_space) then
      do j = 1, cp%no_parameters
        do i = 1, cp%ntiter + 1
          call tdf_set_numerical(cp%f(j), i, tdf(cp%f(j), i) / tdf(par_common%td_penalty(j), i) )
        end do
      end do
    end if

    call pop_sub()
  end subroutine parameters_apply_envelope
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, ep)
    type(oct_control_parameters_t), target, intent(in) :: cp
    type(epot_t), intent(inout) :: ep

    integer :: j
    logical :: change_rep
    type(oct_control_parameters_t), pointer :: par
    call push_sub('parameters.parameters_to_h')

    ! First, check that we are in the real space representation.
    if(cp%current_representation .eq. ctr_parameter_frequency_space) then
      ALLOCATE(par, 1)
      call parameters_copy(par, cp)
      call parameters_to_realtime(par)
      change_rep = .true.
    else
      par => cp
      change_rep = .false.
    end if

    select case(par_common%mode)
    case(parameter_mode_epsilon, parameter_mode_f)
      do j = 1, cp%no_parameters
        call laser_set_f(ep%lasers(j), par%f(j))
      end do
    case(parameter_mode_phi)
      call laser_set_phi(ep%lasers(1), par%f(1))
    end select

    if(change_rep) then 
      call parameters_end(par)
    else
      nullify(par)
    end if
    call pop_sub()
  end subroutine parameters_to_h
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_h_val(cp, ep, val)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep
    integer, intent(in) :: val

    integer :: j

    do j = 1, cp%no_parameters
      call laser_set_f_value(ep%lasers(j), val, real(tdf(cp%f(j), val)) )
    end do

  end subroutine parameters_to_h_val
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j

    call push_sub('parameters.parameters_end')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
    end do
    deallocate(cp%f)
    nullify(cp%f)
    deallocate(cp%alpha)
    nullify(cp%alpha)
    if(associated(cp%utransf)) then
      deallocate(cp%utransf)
      nullify(cp%utransf)
    end if
    if(associated(cp%utransfi)) then
      deallocate(cp%utransfi)
      nullify(cp%utransfi)
    end if

    call pop_sub()
  end subroutine parameters_end
  ! ---------------------------------------------------------


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
      call parameters_to_realtime(par)
      change_rep = .true.
    else
      par => cp
      change_rep = .false.
    end if

    iunit = io_open(trim(filename)//'/Fluence', action='write')
    write(iunit, '(a,es20.8e3)') 'Fluence = ', parameters_fluence(par)
    call io_close(iunit)


    select case(par_common%mode)
    case(parameter_mode_epsilon)

      do j = 1, cp%no_parameters
        if(cp%no_parameters > 1) then
          write(digit,'(i2.2)') j
          iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cp', action='write')
        end if
        write(iunit,'(2a20)') '#       t [a.u]      ', '        f(t)         '
        do i = 1, cp%ntiter + 1
          t = (i-1)*cp%dt
          write(iunit, '(2es20.8e3)') t, real(tdf(par%f(j), t))
        end do
        call io_close(iunit)
      end do

    case(parameter_mode_f)

      do j = 1, cp%no_parameters
        if(cp%no_parameters > 1) then
          write(digit,'(i2.2)') j
          iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cp', action='write')
        end if
        write(iunit,'(3a20)') '#       t [a.u]      ', '        f(t)         ', '        E(t)         '
        do i = 1, cp%ntiter + 1
          t = (i-1)*cp%dt
          write(iunit, '(3es20.8e3)') t, real(tdf(par%f(j), t)), real(tdf(par%f(j), t)) * cos(par%w0*t)
        end do
        call io_close(iunit)
      end do

    case(parameter_mode_phi)

      ! In this case, there is only one parameter (for the moment)
      iunit = io_open(trim(filename)//'/cp', action='write')
      write(iunit,'(4a20)') '#       t [a.u]      ', '        f(t)         ', &
                            '        phi(t)       ', '        E(t)         ' 
      do i = 1, cp%ntiter + 1
        t = (i-1)*cp%dt
        write(iunit, '(4es20.8e3)') t, real(tdf(par_common%f, t)), real(tdf(par%f(1), t)), real(tdf(par_common%f, t)) * &
          cos(par%w0*t + real(tdf(par%f(1), t)) )
      end do
      call io_close(iunit)

    case(parameter_mode_f_and_phi)

      ! In this case, there is only one parameter (for the moment)
      iunit = io_open(trim(filename)//'/cp', action='write')
      write(iunit,'(3a20)') '#       t [a.u]      ', '        f(t)         ', &
                            '        phi(t)       ', '        E(t)         '
      do i = 1, cp%ntiter + 1
        t = (i-1)*cp%dt
        write(iunit, '(4es20.8e3)') t, real(tdf(par%f(1), t)), real(tdf(par%f(2), t)), &
          real(tdf(par%f(1), t)) * cos(par%w0*t + real(tdf(par%f(2), t)) )
      end do
      call io_close(iunit)

    end select

    ! This Fourier analysis is only done for the case of optimization of the full electric field.
    ! TODO: extend this to the envelope and/or phaseo optimization.
    if(present(fourier) .and. (par_common%mode .eq. parameter_mode_epsilon) ) then
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


  ! ---------------------------------------------------------
  ! Gets the fluence of the laser field, defined as:
  ! parameters_fluence = \sum_i^{no_parameters} \integrate_0^T |epsilon(t)|^2
  ! ---------------------------------------------------------
  FLOAT function parameters_fluence(par)
    type(oct_control_parameters_t), target, intent(in) :: par
    type(oct_control_parameters_t), pointer :: par_
    integer :: j, i
    FLOAT :: t, fi, phi
    type(tdf_t) :: f
    logical :: change_rep
    call push_sub('parameters.parameters_fluence')

    if(par%current_representation .eq. ctr_parameter_frequency_space) then
      ALLOCATE(par_, 1)
      call parameters_copy(par_, par)
      call parameters_to_realtime(par_)
      change_rep = .true.
    else
      par_ => par
      change_rep = .false.
    end if

    parameters_fluence = M_ZERO

    select case(par_common%mode)
    case(parameter_mode_epsilon)
      do j = 1, par%no_parameters
        parameters_fluence = parameters_fluence + tdf_dot_product(par%f(j), par%f(j))
      end do
    case(parameter_mode_f)
      do j = 1, par%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par%f(j))
        if(par%current_representation .eq. ctr_parameter_frequency_space) then
          call tdf_sineseries_to_numerical(f)
        end if
        call tdf_cosine_multiply(par%w0, f)
        parameters_fluence = parameters_fluence + tdf_dot_product(f, f)
        call tdf_end(f)
      end do
    case(parameter_mode_phi)
      call tdf_init(f)
      call tdf_copy(f, par%f(1))
      if(par%current_representation .eq. ctr_parameter_frequency_space) then
        call tdf_sineseries_to_numerical(f)
      end if
      do i = 1, par%ntiter + 1
        t = (i-1)*par%dt
        fi = tdf(par_common%f, i)
        phi = real(tdf(f, i)) 
        call tdf_set_numerical(f, i, fi *cos(par%w0*t+phi))
      end do
      parameters_fluence = tdf_dot_product(f, f)
      call tdf_end(f)
    case(parameter_mode_f_and_phi)
      call tdf_init(f)
      call tdf_copy(f, par_%f(1))
      do i = 1, par_%ntiter + 1
        t = (i-1)*par_%dt
        fi = tdf(par_%f(1), i)
        phi = real(tdf(par_%f(2), i))
        call tdf_set_numerical(f, i, fi*cos(par%w0*t+phi))
      end do
      parameters_fluence = tdf_dot_product(f, f)
      call tdf_end(f)
    end select

    if(change_rep) then
      call parameters_end(par_)
    else
      nullify(par_)
    end if

    call pop_sub()
  end function parameters_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function parameters_j2(par) result(j2)
    type(oct_control_parameters_t), target, intent(in) :: par
    type(oct_control_parameters_t), pointer :: par_
    integer :: i, j, k
    FLOAT   :: t, integral, fi, phi, tdp
    type(tdf_t) :: f
    logical :: change_rep

    call push_sub('parameters.parameters_j2')

    ASSERT(par%current_representation .eq. par%representation)

    if(par%current_representation .eq. ctr_parameter_frequency_space) then
      ALLOCATE(par_, 1)
      call parameters_copy(par_, par)
      call parameters_to_realtime(par_)
      change_rep = .true.
    else
      par_ => par
      change_rep = .false.
    end if

    integral = M_ZERO
    select case(par_common%mode)
    case(parameter_mode_epsilon)
      do j = 1, par_%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par_%f(j))
        do i = 1, par_%ntiter + 1
          t = (i-1)*par_%dt
          fi = tdf(par_%f(j), i)
          tdp = sqrt(real(tdf(par_common%td_penalty(j), i)))
          call tdf_set_numerical(f, i, fi*tdp)
        end do
        integral = integral + tdf_dot_product(par_%f(j), par_%f(j))
        call tdf_end(f)
      end do
    case(parameter_mode_f)
      do j = 1, par%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par_%f(j))
        if(par_%current_representation .eq. ctr_parameter_frequency_space) then
          call tdf_sineseries_to_numerical(f)
        end if
        do i = 1, par_%ntiter + 1
          t = (i-1)*par_%dt
          fi = tdf(par_%f(j), i)
          tdp = sqrt(real(tdf(par_common%td_penalty(j), i)))
          call tdf_set_numerical(f, i, fi*tdp*cos(par_%w0*t))
        end do
        integral = integral + tdf_dot_product(f, f)
        call tdf_end(f)
      end do
    case(parameter_mode_phi)
      call tdf_init(f)
      call tdf_copy(f, par_%f(1))
      do i = 1, par_%ntiter + 1
        t = (i-1)*par_%dt
        fi = tdf(par_common%f, i)
        phi = real(tdf(par_%f(1), i))
        tdp = sqrt(real(tdf(par_common%td_penalty(1), i)))
        call tdf_set_numerical(f, i, fi*cos(par_%w0*t+phi))
      end do
      integral = tdf_dot_product(f, f)
      call tdf_end(f)
    case(parameter_mode_f_and_phi)
      call tdf_init(f)
      call tdf_copy(f, par_%f(1))
      do i = 1, par_%ntiter + 1
        t = (i-1)*par_%dt
        fi = tdf(par_%f(1), i)
        phi = real(tdf(par_%f(2), i))
        tdp = sqrt(real(tdf(par_common%td_penalty(1), i)))
        call tdf_set_numerical(f, i, tdp*fi*cos(par_%w0*t+phi))
      end do
      integral = tdf_dot_product(f, f)
      call tdf_end(f)
    end select

    j2 = - par_%alpha(1) * (integral - par_%targetfluence)

    if(change_rep) then
      call parameters_end(par_)
    else
      nullify(par_)
    end if

    call pop_sub()
  end function parameters_j2
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_set_fluence(par)
    type(oct_control_parameters_t), intent(inout) :: par
    FLOAT   :: old_fluence
    integer :: j

    call push_sub('parameters.parameters_set_fluence')

    old_fluence = parameters_fluence(par) 
    do j = 1, par%no_parameters
      call tdf_scalar_multiply( sqrt(par%targetfluence/old_fluence) , par%f(j) )
    end do

    call pop_sub()
  end subroutine parameters_set_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_set_alpha(par, alpha)
    type(oct_control_parameters_t), intent(inout) :: par
    FLOAT, intent(in) :: alpha
    par%alpha(:) = alpha
  end subroutine parameters_set_alpha
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
    cp_out%nfreqs = cp_in%nfreqs
    cp_out%targetfluence = cp_in%targetfluence
    cp_out%intphi = cp_in%intphi
    cp_out%representation = cp_in%representation
    cp_out%current_representation = cp_in%current_representation
    cp_out%omegamax = cp_in%omegamax
    cp_out%w0 = cp_in%w0
    ALLOCATE(cp_out%f(cp_out%no_parameters), cp_out%no_parameters)
    ALLOCATE(cp_out%alpha(cp_out%no_parameters), cp_out%no_parameters)
    do j = 1, cp_in%no_parameters
      cp_out%alpha(j) = cp_in%alpha(j)
      call tdf_init(cp_out%f(j))
      call tdf_copy(cp_out%f(j), cp_in%f(j))
    end do
    if(associated(cp_in%utransf)) then
      ALLOCATE(cp_out%utransf(cp_out%nfreqs, cp_out%nfreqs), cp_out%nfreqs**2)
      cp_out%utransf = cp_in%utransf
    end if
    if(associated(cp_in%utransfi)) then
      ALLOCATE(cp_out%utransfi(cp_out%nfreqs, cp_out%nfreqs), cp_out%nfreqs**2)
      cp_out%utransfi = cp_in%utransfi
    end if

  end subroutine parameters_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_par_to_x(par, x)
    type(oct_control_parameters_t), intent(in)    :: par
    REAL_DOUBLE,                    intent(inout) :: x(:)
    integer :: j, k, n
    FLOAT :: sumx2
    FLOAT, allocatable :: ep(:), e(:)

    n = par%nfreqs
    ALLOCATE(e(n), n)
    ALLOCATE(ep(n), n)

    ASSERT(n-1 .eq. size(x))
    ASSERT(par%current_representation .eq. ctr_parameter_frequency_space)

    do j =  1, n
      ep(j) = tdf(par%f(1), j)
    end do
    if(par_common%mode .eq. parameter_mode_f) then
      e = matmul(par%utransf, ep)
    else
      e = ep
    end if
    x(n-1) = atan2(e(n), e(n-1))
    do k = n-2, 1, -1
      sumx2 = M_ZERO
      do j = n, k+1, -1
        sumx2 = sumx2 + e(j)**2
      end do
      x(k) = atan2(sqrt(sumx2), e(k))
    end do

    deallocate(e, ep)
  end subroutine parameters_par_to_x
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_x_to_par(par, x)
    type(oct_control_parameters_t), intent(inout) :: par
    REAL_DOUBLE,                          intent(in)    :: x(:)

    integer :: j, k, n
    FLOAT, allocatable :: e(:), ep(:)
    call push_sub('parameters.parameters_x_to_par')

    n = par%nfreqs
    ASSERT(n-1 .eq. size(x))
    ASSERT(par%current_representation .eq. ctr_parameter_frequency_space)

    ALLOCATE(e(n), n)
    ALLOCATE(ep(n), n)
    e = M_ZERO; ep = M_ZERO

    if(n.eq.2) then
      e(1) = cos(x(1))
      e(2) = sin(x(1))
      call pop_sub()
      return
    elseif(n.eq.3) then
      e(1) = cos(x(1))
      e(2) = sin(x(1))*cos(x(2))
      e(3) = sin(x(1))*sin(x(2))
      call pop_sub()
      return
    end if

    e(1) = cos(x(1))
    e(2) = sin(x(1))*cos(x(2))
    e(3) = sin(x(1))*sin(x(2))*cos(x(3))
    do j = 4, n - 1
      e(j) = M_ONE
      do k = 1, j - 1
        e(j) = e(j) * sin(x(k))
      end do
      e(j) = e(j) * cos(x(j))
    end do
    e(n) = M_ONE
    do k = 1, n - 2
      e(n) = e(n) * sin(x(k))
    end do
    e(n) = e(n) * sin(x(n-1))

    select case(par_common%mode)
    case(parameter_mode_epsilon, parameter_mode_f)
      e = sqrt(par%targetfluence) * e
    case(parameter_mode_phi)
      e = sqrt(par%intphi)*e
    end select

    if(par_common%mode .eq. parameter_mode_f) then
      ep = matmul(par%utransfi, e)
    else
      ep = e
    end if
    call tdf_set_numerical(par%f(1), ep)

    deallocate(e, ep)
    call pop_sub()
  end subroutine parameters_x_to_par
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_update(cp, cpp, dir, iter, delta, eta, d1, dl, dq)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(oct_control_parameters_t), intent(in)    :: cpp
    character(len=1),               intent(in)    :: dir
    integer,                        intent(in)    :: iter
    FLOAT,                          intent(in)    :: delta, eta
    CMPLX,                          intent(in)    :: d1
    CMPLX,                          intent(in)    :: dl(:), dq(:)

    FLOAT :: value
    integer :: j

    call push_sub('parameters.parameters_update')

    select case(dir)
      case('f')
        do j = 1, cp%no_parameters
          value = (M_ONE / parameters_alpha(cp, j)) * aimag(d1*dl(j)) / &
           ( tdf(par_common%td_penalty(j), iter) - M_TWO*aimag(dq(j)) )
          value = (M_ONE - delta)*tdf(cpp%f(j), iter) + delta * value
          call tdf_set_numerical(cp%f(j), iter, value)
          if(iter+1 <= cp%ntiter + 1)  call tdf_set_numerical(cp%f(j), iter+1, value)
          if(iter+2 <= cp%ntiter + 1)  call tdf_set_numerical(cp%f(j), iter+2, value)
        end do

      case('b')
        do j = 1, cp%no_parameters
          value = (M_ONE / parameters_alpha(cp, j)) * aimag(d1*dl(j)) / &
           ( tdf(par_common%td_penalty(j), iter+1) - M_TWO*aimag(dq(j)) ) 
          value = (M_ONE - eta)*tdf(cpp%f(j), iter+1) + eta * value
          call tdf_set_numerical(cp%f(j), iter+1, value)
          if(iter > 0) call tdf_set_numerical(cp%f(j), iter, value)
          if(iter-1 > 0) call tdf_set_numerical(cp%f(j), iter-1, value)
        end do
    end select

    call pop_sub()
  end subroutine parameters_update
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_alpha(par, j)
    type(oct_control_parameters_t), intent(in) :: par
    integer,                        intent(in) :: j
    parameters_alpha = par%alpha(j)
  end function parameters_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_targetfluence(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_targetfluence = par%targetfluence
  end function parameters_targetfluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function parameters_number(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_number = par%no_parameters
  end function parameters_number
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function parameters_nfreqs(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_nfreqs = par%nfreqs
  end function parameters_nfreqs
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_w0(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_w0 = par%w0
  end function parameters_w0
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_filter(par, filter)
    type(oct_control_parameters_t), intent(inout) :: par
    type(filter_t),                 intent(inout) :: filter
    integer :: j

    do j = 1, par%no_parameters
      call filter_apply(par%f(j), filter)
    end do

  end subroutine parameters_filter
  ! ---------------------------------------------------------

end module opt_control_parameters_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
