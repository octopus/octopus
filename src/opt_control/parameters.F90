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
  use loct_m
  use lalg_adv_m
  use string_m
  use datasets_m
  use varinfo_m
  use global_m
  use mpi_m
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
  use math_m
  use profiling_m

  implicit none

  private
  public :: oct_control_parameters_t,     &
            parameters_mod_init,          &
            parameters_mod_close,         &
            parameters_init,              &
            parameters_representation,    &
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
            parameters_to_basis,          &
            parameters_prepare_initial,   &
            parameters_fluence,           &
            parameters_j2,                &
            parameters_basis_to_theta,    &
            parameters_theta_to_basis,    &
            parameters_get_theta,         &
            parameters_set_theta,         &
            parameters_randomize,         &
            parameters_update,            &
            parameters_number,            &
            parameters_bounds,            &
            parameters_dof,               &
            parameters_w0,                &
            parameters_alpha,             &
            parameters_targetfluence,     &
            parameters_filter


  integer, public, parameter ::   &
    ctr_real_time            = 1, &
    ctr_sine_fourier_series  = 2, &
    ctr_fourier_series       = 3, &
    ctr_zero_fourier_series  = 4


  integer, parameter :: parameter_mode_none      = 0, &
                        parameter_mode_epsilon   = 1, &
                        parameter_mode_f         = 2, &
                        parameter_mode_phi       = 3

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


    type(tdf_t)    :: f ! This is the envelope of the laser field, only used in the phase-only
                        ! optimization (necessary to store it in order to calculate fluences)

  end type oct_parameters_common_t


  type oct_control_parameters_t
    private
    integer :: no_parameters = 0
    integer :: dim           = 0
    integer :: dof           = 0
    FLOAT   :: intphi        = M_ZERO
    type(tdf_t), pointer :: f(:) => NULL()
    FLOAT, pointer :: alpha(:)   => NULL()

    integer :: current_representation = 0

    FLOAT   :: w0       = M_ZERO
    FLOAT, pointer :: utransf(:, :)  => NULL()
    FLOAT, pointer :: utransfi(:, :) => NULL()

    REAL_DOUBLE, pointer :: theta(:) => NULL()
  end type oct_control_parameters_t
  
  ! the next variable has to be a pointer to avoid a bug in the IBM compiler
  type(oct_parameters_common_t), pointer :: par_common
  type(mix_t) :: parameters_mix

contains


  ! ---------------------------------------------------------
  ! Initializes the module, should be the first subroutine to be called (the last one
  ! should be parameters_mod_close, when the module is no longer to be used.
  !
  ! It fills the module variable "par_common", whose type is par_common_t, with 
  ! information obtained from the inp file.
  !
  ! Output argument "mode_fixed_fluence" is also given a value, depending on whether
  ! the user requires a fixed-fluence run (.true.) or not (.false.).
  ! ---------------------------------------------------------
  subroutine parameters_mod_init(ep, dt, max_iter, mode_fixed_fluence, parametrized_controls)
    type(epot_t), intent(inout)                   :: ep
    FLOAT, intent(in)                             :: dt
    integer, intent(in)                           :: max_iter
    logical, intent(out)                          :: mode_fixed_fluence
    logical, intent(in)                           :: parametrized_controls

    character(len=1024)      :: expression
    integer :: i, j, no_lines, steps, iunit
    FLOAT   :: octpenalty, t, f_re, f_im, total_time
    type(block_t)            :: blk

    call push_sub('parameters.parameters_read')

    if(.not. associated(par_common)) ALLOCATE(par_common, 1) 

    call messages_print_stress(stdout, "OCT: Info about control functions")

    !%Variable OCTParameterRepresentation
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default control_parameter_fourier_series
    !%Description
    !% If "OCTControlRepresentation = control_function_parametrized", one must 
    !% specify the kind of parameters that determine the control function.
    !% If "OCTControlRepresentation = control_function_real_time", then this variable
    !% is ignored, and the control function is handled directly in real time.
    !%Option control_sine_fourier_series 2
    !% The control function is expanded in a sine Fourier series (which implies that it
    !% it starts and ends at zero).
    !%Option control_fourier_series 3
    !% The control function is expanded as a full Fourier series (although it must, of 
    !% course, be a real function).
    !%Option control_zero_fourier_series 4
    !% The control function is expanded as a Fourier series, but assuming (1) that the zero
    !% frequency component is zero, and (2) the control function, integrated in time, adds
    !% up to zero (this essentially means that the sum of all the cosine coefficients is zero).
    !%End
    if (parametrized_controls) then
      call loct_parse_int(datasets_check('OCTParameterRepresentation'), &
        ctr_fourier_series, par_common%representation)
      if(.not.varinfo_valid_option('OCTParameterRepresentation', par_common%representation)) &
        call input_error('OCTParameterRepresentation')
      select case(par_common%representation)
      case(ctr_sine_fourier_series)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a sine '
        write(message(2), '(a)') '      Fourier series.'
        call write_info(2)
      case(ctr_fourier_series)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series.'
        call write_info(1)
      case(ctr_zero_fourier_series)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series,'
        write(message(2), '(a)') '      in which the zero-frequency component is assumed to be zero.'
        call write_info(1)
      end select
    else
      par_common%representation = ctr_real_time
      write(message(1), '(a)') 'Info: The OCT control functions will be represented in real time.'
      call write_info(1)
    end if

    !%Variable OCTParameterOmegaMax
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default -1.0
    !%Description
    !% The Fourier series that can be used to represent the control functions must be truncated;
    !% the truncation is given by a cut-off frequency which is determined by this variable.
    !%End
    call loct_parse_float(datasets_check('OCTParameterOmegaMax'), -M_ONE, par_common%omegamax)
    if(par_common%representation .ne. ctr_real_time) then
      write(message(1), '(a)')         'Info: The representation of the OCT control parameters will be restricted'
      write(message(2), '(a,f10.5,a)') '      with an energy cut-off of ', &
        par_common%omegamax / units_out%energy%factor, ' ['//trim(units_out%energy%abbrev) // ']'
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
    call loct_parse_float(datasets_check('OCTFixFluenceTo'), M_ZERO, par_common%targetfluence)

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
    call loct_parse_logical(datasets_check('OCTFixInitialFluence'), .true., &
      par_common%fix_initial_fluence)

    !%Variable OCTControlFunctionType
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 1
    !%Description
    !% The control function may fully determine the time-dependent form of the 
    !% external field, or only the envelope function of this external field, or its phase. 
    !% Or, we may have two different control functions, one of them providing the phase 
    !% and the other one, the envelope.
    !%
    !% Note that, if "OCTControlRepresentation = control_function_real_time", then the control
    !% function must *always* determine the full external field.
    !%Option parameter_mode_epsilon   1
    !% In this case, the control function determines the full control function. That is, 
    !% if we are considering the electric field of a laser, the time-dependent electric field.
    !%Option parameter_mode_f         2
    !% The optimization process attempts to find the best possible envelope. The full 
    !% control field is this envelope times a cosine function with a "carrier" frequency. 
    !% This carrier frequencey is given by the carrier frequency of the "TDExternalFields" 
    !% in the inp file.
    !%Option parameter_mode_phi       3
    !% The optimization process attempts to find the best possible time-dependent phase. That is,
    !% the external field would be given by a function in the form e(t) = f(t)*cos(w0*t+phi(t)), 
    !% where f(t) is an "envelope", w0 a carrier frequency, and phi(t) the td phase that we 
    !% wish to optimize.
    !%End
    call loct_parse_int(datasets_check('OCTControlFunctionType'), parameter_mode_epsilon, par_common%mode)
    if(.not.varinfo_valid_option('OCTControlFunctionType', par_common%mode)) &
      call input_error('OCTControlFunctionType')
    if ( (.not.parametrized_controls)  .and.  (par_common%mode .ne. parameter_mode_epsilon) ) &
      call input_error('OCTControlFunctionType')
    call messages_print_var_option(stdout, 'OCTControlFunctionType', par_common%mode)


    ! The laser field is defined by "td functions", as implemented in module "tdf_m". At this point, they
    ! can be in "non-numerical" representation (i.e. described with a set of parameters, e.g. frequency, width, etc).
    ! We need them to be in numerical form (i.e. time grid, values at the time grid). Here we do the transformation.
    ! It cannot be done before calling parameters_mod_init because we need to pass the omegamax value.
    do i = 1, ep%no_lasers
      select case(par_common%mode)
      case(parameter_mode_epsilon)
        call laser_to_numerical_all(ep%lasers(i), dt, max_iter, par_common%omegamax)
      case default
        call laser_to_numerical(ep%lasers(i), dt, max_iter, par_common%omegamax)
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

    ! Fix the number of control functions: if we have "traditional" QOCT (i.e. the control functions
    ! are represented directly in real time, then the number of control functions can be larger than
    ! one; it will be the number of lasers found in the input file. Otherwise, if the control function(s)
    ! are parametrized ("OCTControlRepresentation = control_function_parametrized"), we only have one
    ! control function. If there is more than one laser field in the input file, the program stops.
    if(parametrized_controls) then
      if(ep%no_lasers > 1) then
        write(message(1), '(a)') 'If "OCTControlRepresentation = control_function_parametrized", you can'
        write(message(2), '(a)') 'have only one external field in the input file.'
        call write_fatal(2)
      end if
      par_common%no_parameters = 1
    else
      par_common%no_parameters = ep%no_lasers
    end if


    mode_fixed_fluence = .false.
    if (par_common%targetfluence .ne. M_ZERO) mode_fixed_fluence = .true.


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
      call loct_parse_float(datasets_check('OCTPenalty'), M_ONE, octpenalty)
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
    !% <math> \frac{1}{\alpha(t)} = \frac{1}{2}( erf((100/T)*(t-T/20))+ erf(-(100/T)*(t-T+T/20)) </math>
    !%End
    steps = max_iter
    ALLOCATE(par_common%td_penalty(par_common%no_parameters), par_common%no_parameters)
    do i = 1, par_common%no_parameters
      call tdf_init_numerical(par_common%td_penalty(i), steps, dt, -M_ONE, initval = M_ONE)
    end do

    if (loct_parse_block(datasets_check('OCTLaserEnvelope'), blk)==0) then

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
        total_time = steps*dt
        if(trim(expression)=='default') then
          do j = 1, steps + 1
            t = (j-1)*dt
            f_re = M_HALF * (loct_erf((CNST(100.0)/total_time)*(t-CNST(0.05)*total_time)) + &
              loct_erf(-(CNST(100.0)/total_time)*(t-total_time+CNST(0.05)*total_time)) )
            call tdf_set_numerical(par_common%td_penalty(i), j, &
              TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        else
          call conv_to_C_string(expression)
          do j = 1, steps+1
            t = (j-1)*dt
            call loct_parse_expression(f_re, f_im, "t", t, expression)
            call tdf_set_numerical(par_common%td_penalty(i), j, &
              TOFLOAT(M_ONE /(f_re + CNST(1.0e-7)))  )
          end do
        end if
      end do

      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open('opt-control/td_penalty', action='write' )
        do j = 1, steps+1
          t = (j-1)*dt
          write(iunit, '(f14.8)', advance='no') t
          do i = 1, par_common%no_parameters - 1
            write(iunit, '(es20.8e3)') M_ONE/tdf(par_common%td_penalty(i), j)
          end do
          write(iunit, '(es20.8e3)') M_ONE/tdf(par_common%td_penalty(par_common%no_parameters), j)
        end do
        write(iunit,'()')
        call io_close(iunit)
      end if

      call loct_parse_block_end(blk)
    end if

    call messages_print_stress(stdout)
    call pop_sub()
  end subroutine parameters_mod_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Before using an oct_control_parameters_t variable, it needs
  ! to be initialized, either by call ing parameters_init, or
  ! by copying anther initialized variable through
  ! parameters_copy.
  ! ---------------------------------------------------------
  subroutine parameters_init(cp, dt, ntiter)
    type(oct_control_parameters_t), intent(inout) :: cp
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter

    integer :: j

    call push_sub('parameters.parameters_init')

    cp%w0              = par_common%w0
    cp%no_parameters   = par_common%no_parameters
    cp%current_representation = ctr_real_time
    call loct_pointer_copy(cp%alpha, par_common%alpha)

    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    do j = 1, cp%no_parameters
      call tdf_init_numerical(cp%f(j), ntiter, dt, par_common%omegamax)
    end do

    ! If the control function is represented directly in real time, the "dimension" (cp%dim) is
    ! the number of values that represent the function in the discretized time axis.
    !
    ! If the control function is parametrized, up to now (in the future this might change), all 
    ! parametrizations are based on a previous basis set expansion (sine-Fourier series, or "normal"
    ! Fourier series with or without the zero term). The parameters are not directly the coefficients
    ! of the control function in this basis set expansion, but are constructed from them (e.g. 
    ! by performing a coordinate transformation to hyperspherical coordinates). The "dimension" (cp%dim)
    ! is the dimension of this basis set.
    select case(par_common%representation)
    case(ctr_real_time)
      cp%dim = ntiter + 1
    case(ctr_sine_fourier_series)
      ! cp%dim is directly the number of frequencies in the sine-Fourier expansion
      cp%dim = tdf_sine_nfreqs(cp%f(1))
    case(ctr_fourier_series)
      ! If nf is the number of frequencies, we will have nf-1 non-zero "sines", nf-1 non-zero "cosines",
      ! and the zero frequency component. Total, 2*(nf-1)+1
      cp%dim = 2*(tdf_nfreqs(cp%f(1))-1)+1
    case(ctr_zero_fourier_series)
      ! If nf is the number of frequencies, we will have nf-1 non-zero "sines", nf-1 non-zero "cosines",
      ! but no zero frequency component. Total, 2*(nf-1)
      cp%dim = 2*(tdf_nfreqs(cp%f(1))-1)
    end select


    ! The "degrees of freedom" cp%dof is the number of parameters that define the control function.
    ! (if it is represented directly in real time, this would be meaningless, but we put the number of 
    ! control functions, times the "dimension", which in this case is the number of time discretization points).
    ! This is not equal to the dimension of the basis set employed (cp%dim), because we may add further
    ! constrains, and do a coordinate transformation to account for them.
    select case(par_common%representation)
    case(ctr_real_time)
      cp%dof = cp%no_parameters * cp%dim
    case(ctr_sine_fourier_series, ctr_fourier_series)
      ! The number of degrees of freedom is one less than the number of basis coefficients, since we
      ! add the constrain of fixed fluence.
      cp%dof = cp%dim - 1
    case(ctr_zero_fourier_series)
      ! The number of degrees of freedom is one less than the number of basis coefficients, since we
      ! add (1) the constrain of fixed fluence, and (2) the constrain of zero-integral.
      cp%dof = cp%dim - 2
    end select


    if(cp%dof <= 0) then
      write(message(1),'(a)') 'Error: The number of degrees of freedom used to describe the control function'
      write(message(2),'(a)') '       is zero or less. This should not happen. Please review your input file.'
      call write_fatal(2)
    else
      if(par_common%representation .ne. ctr_real_time) then
        write(message(1), '(a)')      'Info: The expansion of the control parameters in a Fourier series'
        write(message(2), '(a,i6,a)') '      expansion implies the use of ', cp%dim, ' basis set functions.'
        write(message(3), '(a,i6,a)') '      The number of degrees of freedom is ', cp%dof,'.'
        call write_info(3)

        ALLOCATE(cp%theta(cp%dof), cp%dof)
        cp%theta = M_ZERO
      end if
    end if

    ! Construct the transformation matrix, if needed.
    call parameters_trans_matrix(cp)

    call pop_sub()
  end subroutine parameters_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! The external fields defined in epot_t "ep" are transferred to
  ! the control functions described in "cp". This shhould have been
  ! initialized previously.
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
    end select

    call pop_sub()
  end subroutine parameters_set
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function parameters_representation()
    parameters_representation = par_common%representation
  end function parameters_representation
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_prepare_initial(par)
    type(oct_control_parameters_t), intent(inout) :: par

    FLOAT   :: dt
    integer :: ntiter

    call push_sub('parameters.parameters_prepare_initial')

    call parameters_apply_envelope(par)

    if(par_common%targetfluence .ne. M_ZERO) then
      if(par_common%targetfluence < M_ZERO) then
        par_common%targetfluence = parameters_fluence(par) 
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with the same'
        write(message(2), '(a,f10.5,a)') '      fluence as the input external fields: F = ', &
          par_common%targetfluence, ' a.u.'
      else
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with a predefined'
        write(message(2), '(a,f10.5,a)') '      fluence: F = ', par_common%targetfluence, ' a.u.'
      end if
      call write_info(2)
      if(par_common%fix_initial_fluence) call parameters_set_fluence(par)
    end if

    ! Now we have to find the "fluence" of the phase, in order to keep it constant.
    select case(par_common%mode)
    case(parameter_mode_phi)
      par%intphi = tdf_dot_product(par%f(1), par%f(1))
      if(par%intphi <= M_ZERO) then
        dt = tdf_dt(par%f(1))
        ntiter = tdf_niter(par%f(1))

        par%intphi = CNST(0.1)*(M_PI/M_TWO)**2*dt*ntiter
      else
        par%intphi = tdf_dot_product(par%f(1), par%f(1))
      end if
    end select

    ! Move to the "native" representation, if necessary.
    call parameters_set_rep(par)

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
    call push_sub('parameters.parameters_set_rep')

    if(par%current_representation .ne. par_common%representation) then
      if(par_common%representation .eq. ctr_real_time) then
        call parameters_to_realtime(par)
      else
        call parameters_to_basis(par)
      end if
    end if

    call pop_sub()
  end subroutine parameters_set_rep
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_basis(par)
    type(oct_control_parameters_t), intent(inout) :: par
    integer :: j
    call push_sub('parameters.parameters_to_basis')

    if(par%current_representation.eq.ctr_real_time) then
      select case(par_common%representation)
      case(ctr_sine_fourier_series)
        do j = 1, par%no_parameters
          call tdf_numerical_to_sineseries(par%f(j))
        end do
        par%current_representation = ctr_sine_fourier_series
        call parameters_basis_to_theta(par)
      case(ctr_fourier_series)
        do j = 1, par%no_parameters
          call tdf_numerical_to_fourier(par%f(j))
        end do
        par%current_representation = ctr_fourier_series
        call parameters_basis_to_theta(par)
      case(ctr_zero_fourier_series)
        do j = 1, par%no_parameters
          call tdf_numerical_to_zerofourier(par%f(j))
        end do
        par%current_representation = ctr_zero_fourier_series
        call parameters_basis_to_theta(par)
      end select
    end if

    call pop_sub()
  end subroutine parameters_to_basis
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_realtime(par)
    type(oct_control_parameters_t), intent(inout) :: par
    integer :: j
    call push_sub('parameters.parameters_to_realtime')

    select case(par%current_representation)
    case(ctr_real_time)
      call pop_sub(); return
    case(ctr_sine_fourier_series)
      call parameters_theta_to_basis(par)
      do j = 1, par%no_parameters
        call tdf_sineseries_to_numerical(par%f(j))
      end do
    case(ctr_fourier_series)
      call parameters_theta_to_basis(par)
      do j = 1, par%no_parameters
        call tdf_fourier_to_numerical(par%f(j))
      end do
    case(ctr_zero_fourier_series)
      call parameters_theta_to_basis(par)
      do j = 1, par%no_parameters
        call tdf_zerofourier_to_numerical(par%f(j))
      end do
    end select

    par%current_representation = ctr_real_time
    call pop_sub()
  end subroutine parameters_to_realtime
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function parameters_diff(p, q) result(res)
    type(oct_control_parameters_t), intent(in) :: p, q
    integer :: i

    call push_sub('parameters.parameters_diff')

    ASSERT(p%current_representation .eq. q%current_representation)

    res = M_ZERO
    do i = 1, p%no_parameters
      res = res + tdf_diff(p%f(i), q%f(i))
    end do

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
    call push_sub('parameters.parameters_mixing_init')
    call mix_init(parameters_mix, par%dim, par%no_parameters, 1)
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

    integer :: i, j, dim!, ntiter
    FLOAT, allocatable :: e_in(:, :, :), e_out(:, :, :), e_new(:, :, :)
    call push_sub('parameters.parameters_mixing')

    dim = par_in%dim
    ALLOCATE(e_in (dim, par_in%no_parameters, 1), dim*par_in%no_parameters)
    ALLOCATE(e_out(dim, par_in%no_parameters, 1), dim*par_in%no_parameters)
    ALLOCATE(e_new(dim, par_in%no_parameters, 1), dim*par_in%no_parameters)
    do i = 1, par_in%no_parameters
      do j = 1, dim
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


  ! ---------------------------------------------------------
  subroutine parameters_apply_envelope(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j, i

    call push_sub('parameters.parameters_apply_envelope')

    ! Do not apply the envelope if the parameters are represented as a sine Fourier series.
    if(par_common%representation .eq. ctr_real_time) then
      do j = 1, cp%no_parameters
        do i = 1, tdf_niter(cp%f(j)) + 1
          call tdf_set_numerical(cp%f(j), i, tdf(cp%f(j), i) / tdf(par_common%td_penalty(j), i) )
        end do
      end do
    end if

    call pop_sub()
  end subroutine parameters_apply_envelope
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep

    integer :: j
    type(oct_control_parameters_t) :: par
    call push_sub('parameters.parameters_to_h')

    call parameters_copy(par, cp)
    call parameters_to_realtime(par)

    select case(par_common%mode)
    case(parameter_mode_epsilon, parameter_mode_f)
      do j = 1, cp%no_parameters
        call laser_set_f(ep%lasers(j), par%f(j))
      end do
    case(parameter_mode_phi)
      call laser_set_phi(ep%lasers(1), par%f(1))
    end select

    call parameters_end(par)
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
      call laser_set_f_value(ep%lasers(j), val, tdf(cp%f(j), val) )
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
    if(associated(cp%theta)) then
      deallocate(cp%theta)
      nullify(cp%theta)
    end if

    call pop_sub()
  end subroutine parameters_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, cp)
    character(len=*), intent(in) :: filename
    type(oct_control_parameters_t), intent(in) :: cp

    integer :: i, j, k, iunit, niter, nfreqs
    FLOAT :: t, wmax, dw, w, wa, wb, dt
    FLOAT, allocatable :: func(:, :)
    CMPLX :: ft, ez, ezdt
    character(len=2) :: digit
    type(oct_control_parameters_t) :: par

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('parameters.parameters_write')

    call io_mkdir(trim(filename))

    call parameters_copy(par, cp)
    call parameters_to_realtime(par)

    iunit = io_open(trim(filename)//'/Fluence', action='write')
    write(iunit, '(a,es20.8e3)') 'Fluence = ', parameters_fluence(par)
    call io_close(iunit)

    niter = tdf_niter(par%f(1))
    ALLOCATE(func(niter+1, cp%no_parameters), (niter+1)*cp%no_parameters)

    select case(par_common%mode)
    case(parameter_mode_epsilon)

      do j = 1, cp%no_parameters
        if(cp%no_parameters > 1) then
          write(digit,'(i2.2)') j
          iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cp', action='write')
        end if
        write(iunit,'(2a20)') '#       t [a.u]      ', '        e(t)         '
        do i = 1, tdf_niter(par%f(j)) + 1
          t = (i-1)*tdf_dt(par%f(j))
          write(iunit, '(2es20.8e3)') t, tdf(par%f(j), t)
          func(i, j) = tdf(par%f(j), t)
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
        write(iunit,'(3a20)') '#       t [a.u]      ', '        e(t)         ', '        f(t)         '
        do i = 1, tdf_niter(par%f(j)) + 1
          t = (i-1)*tdf_dt(par%f(j))
          write(iunit, '(3es20.8e3)') t, tdf(par%f(j), t) * cos(par%w0*t), tdf(par%f(j), t)
          func(i, j) = tdf(par%f(j), t) * cos(par%w0*t)
        end do
        call io_close(iunit)
      end do

    case(parameter_mode_phi)

      ! In this case, there is only one parameter (for the moment)
      iunit = io_open(trim(filename)//'/cp', action='write')
      write(iunit,'(4a20)') '#       t [a.u]      ', '        e(t)         ', &
                            '         f(t)        ', '       phi(t)        ' 
      do i = 1, tdf_niter(par%f(j)) + 1
        t = (i-1)*tdf_dt(par%f(j))
        write(iunit, '(4es20.8e3)') t, tdf(par_common%f, t) * &
          cos(par%w0*t + tdf(par%f(1), t) ), tdf(par_common%f, t), tdf(par%f(1), t)
        func(i, 1) = tdf(par_common%f, t) * cos(par%w0*t + tdf(par%f(1), t) )
      end do
      call io_close(iunit)

    end select


    !Now, the Fourier transforms.
    select case(par_common%mode)
    case(parameter_mode_epsilon)

      do j = 1, cp%no_parameters
        if(cp%no_parameters > 1) then
          write(digit,'(i2.2)') j
          iunit = io_open(trim(filename)//'/cpw-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cpw', action='write')
        end if
        write(iunit,'(3a20)') '#       w [a.u]      ', '      Re[e(w)]       ', &
                              '      Im[e(w)]       '

        nfreqs = 1000
        wa = M_ZERO
        wb = M_THREE ! hard coded to three atomic units... this should be improved.
        wmax = wb
        dw = wmax/(nfreqs-1)
        dt = tdf_dt(par%f(1))

        do k = 1, nfreqs
          w = wa + (k-1)*dw
          ft = M_z0
          ez = M_z1
          ezdt = exp(M_zI*w*tdf_dt(par%f(j)))
          do i = 1, niter + 1
            t = (i-1)*dt
            ft = ft + func(i, j)*ez
            ez = ez*ezdt
          end do
          ft = ft*dt
          write(iunit,'(3es20.8e3)') w, real(ft), aimag(ft)
        end do
      end do

    case(parameter_mode_f, parameter_mode_phi)
      iunit = io_open(trim(filename)//'/cpw', action='write')
      write(iunit,'(3a20)') '#       w [a.u]      ', '      Re[e(w)]       ', &
                            '      Im[e(w)]       '
      
      nfreqs = 1000
      wa = cp%w0 - M_THREE * par_common%omegamax
      wb = cp%w0 + M_THREE * par_common%omegamax
      wmax = CNST(6.0)*par_common%omegamax
      dw = wmax/(nfreqs-1)
      dt = tdf_dt(par%f(1))

      do j = 1, nfreqs
        w = wa + (j-1)*dw
        ft = M_z0
        ez = M_z1
        ezdt = exp(M_zI*w*tdf_dt(par%f(1)))
        do i = 1, niter + 1
          t = (i-1)*dt
          ft = ft + func(i, 1)*ez
          ez = ez*ezdt
        end do
        ft = ft*dt
        write(iunit,'(3es20.8e3)') w, real(ft), aimag(ft)
      end do

      call io_close(iunit)
    end select

    call parameters_end(par)
    call pop_sub()
  end subroutine parameters_write
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Gets the fluence of the laser field, defined as:
  ! parameters_fluence = \sum_i^{no_parameters} \integrate_0^T |epsilon(t)|^2
  ! ---------------------------------------------------------
  FLOAT function parameters_fluence(par)
    type(oct_control_parameters_t), intent(in) :: par
    type(oct_control_parameters_t)             :: par_
    integer :: j, i
    FLOAT :: t, fi, phi
    type(tdf_t) :: f
    call push_sub('parameters.parameters_fluence')

    call parameters_copy(par_, par)
    call parameters_to_realtime(par_)

    parameters_fluence = M_ZERO

    select case(par_common%mode)
    case(parameter_mode_epsilon)
      do j = 1, par_%no_parameters
        parameters_fluence = parameters_fluence + tdf_dot_product(par_%f(j), par_%f(j))
      end do
    case(parameter_mode_f)
      do j = 1, par%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par_%f(j))
        call tdf_cosine_multiply(par%w0, f)
        parameters_fluence = parameters_fluence + tdf_dot_product(f, f)
        call tdf_end(f)
      end do
    case(parameter_mode_phi)
      call tdf_init(f)
      call tdf_copy(f, par_%f(1))
      do i = 1, tdf_niter(f) + 1
        t = (i-1)*tdf_dt(f)
        fi = tdf(par_common%f, i)
        phi = real(tdf(f, i)) 
        call tdf_set_numerical(f, i, fi *cos(par%w0*t+phi))
      end do
      parameters_fluence = tdf_dot_product(f, f)
      call tdf_end(f)
    end select

    call parameters_end(par_)
    call pop_sub()
  end function parameters_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function parameters_j2(par) result(j2)
    type(oct_control_parameters_t), intent(in) :: par
    type(oct_control_parameters_t)             :: par_
    integer :: i, j
    FLOAT   :: t, integral, fi, phi, tdp
    type(tdf_t) :: f

    call push_sub('parameters.parameters_j2')

    ASSERT(par%current_representation .eq. par_common%representation)

    call parameters_copy(par_, par)
    call parameters_to_realtime(par_)

    integral = M_ZERO
    select case(par_common%mode)
    case(parameter_mode_epsilon)
      do j = 1, par_%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par_%f(j))
        do i = 1, tdf_niter(f) + 1
          t = (i-1)*tdf_dt(f)
          fi = tdf(par_%f(j), i)
          tdp = sqrt(real(tdf(par_common%td_penalty(j), i)))
          call tdf_set_numerical(f, i, fi*tdp)
        end do
        integral = integral + tdf_dot_product(f, f)
        call tdf_end(f)
      end do
    case(parameter_mode_f)
      do j = 1, par_%no_parameters
        call tdf_init(f)
        call tdf_copy(f, par_%f(j))
        if(par_%current_representation .eq. ctr_sine_fourier_series) then
          call tdf_sineseries_to_numerical(f)
        end if
        do i = 1, tdf_niter(f) + 1
          t = (i-1)*tdf_dt(f)
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
      do i = 1, tdf_niter(f) + 1
        t = (i-1)*tdf_dt(f)
        fi = tdf(par_common%f, i)
        phi = real(tdf(par_%f(1), i))
        tdp = sqrt(real(tdf(par_common%td_penalty(1), i)))
        call tdf_set_numerical(f, i, fi*cos(par_%w0*t+phi))
      end do
      integral = tdf_dot_product(f, f)
      call tdf_end(f)
    end select

    j2 = - par_%alpha(1) * (integral - par_common%targetfluence)

    call parameters_end(par_)
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
      call tdf_scalar_multiply( sqrt(par_common%targetfluence/old_fluence) , par%f(j) )
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
    call push_sub('parameters.parameters_copy')

    cp_out%no_parameters = cp_in%no_parameters
    cp_out%dim = cp_in%dim
    cp_out%dof = cp_in%dof
    cp_out%intphi = cp_in%intphi
    cp_out%current_representation = cp_in%current_representation
    cp_out%w0 = cp_in%w0
    call loct_pointer_copy(cp_out%alpha, cp_in%alpha)
    ALLOCATE(cp_out%f(cp_out%no_parameters), cp_out%no_parameters)
    do j = 1, cp_in%no_parameters
      call tdf_init(cp_out%f(j))
      call tdf_copy(cp_out%f(j), cp_in%f(j))
    end do
    call loct_pointer_copy(cp_out%utransf, cp_in%utransf)
    call loct_pointer_copy(cp_out%utransfi, cp_in%utransfi)
    call loct_pointer_copy(cp_out%theta, cp_in%theta)

    call pop_sub()
  end subroutine parameters_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_randomize(par)
    type(oct_control_parameters_t), intent(inout) :: par

     integer :: j
     call push_sub('parameters.parameters_randomize')

     ASSERT(par_common%representation .ne. ctr_real_time)

     call parameters_set_rep(par)

     select case(par_common%mode)
     case(parameter_mode_epsilon, parameter_mode_f)
       do j = 1, par%no_parameters
         call tdf_set_random(par%f(j))
       end do
     case(parameter_mode_phi)
       call tdf_set_random(par%f(1))
     end select

     call pop_sub()
  end subroutine parameters_randomize
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
          if(iter+1 <= tdf_niter(cp%f(j)) + 1)  call tdf_set_numerical(cp%f(j), iter+1, value)
          if(iter+2 <= tdf_niter(cp%f(j)) + 1)  call tdf_set_numerical(cp%f(j), iter+2, value)
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
  FLOAT pure function parameters_alpha(par, j)
    type(oct_control_parameters_t), intent(in) :: par
    integer,                        intent(in) :: j
    parameters_alpha = par%alpha(j)
  end function parameters_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function parameters_targetfluence()
    parameters_targetfluence = par_common%targetfluence
  end function parameters_targetfluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function parameters_number(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_number = par%no_parameters
  end function parameters_number
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine parameters_bounds(par, lower_bounds, upper_bounds)
    type(oct_control_parameters_t), intent(in)  :: par
    FLOAT,                          intent(out) :: lower_bounds(:)
    FLOAT,                          intent(out) :: upper_bounds(:)
    integer :: dog

    call push_sub('parameters.parameters_bounds')

    upper_bounds = M_PI
    dog = parameters_dof(par)

    select case(par_common%mode)
    case(parameter_mode_epsilon, parameter_mode_f, parameter_mode_phi)
      lower_bounds(1:dog-1) = M_ZERO
      lower_bounds(dog)     = -M_PI
    end select

    call pop_sub()
  end subroutine parameters_bounds
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function parameters_dof(par)
    type(oct_control_parameters_t), intent(in) :: par
    parameters_dof = par%dof
  end function parameters_dof
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function parameters_w0(par)
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


  ! ---------------------------------------------------------
  subroutine parameters_mod_close()
    integer :: j

    deallocate(par_common%alpha); nullify(par_common%alpha)
    do j = 1, par_common%no_parameters
      call tdf_end(par_common%td_penalty(j))
    end do
    deallocate(par_common%td_penalty); nullify(par_common%td_penalty)
    if(associated(par_common)) deallocate(par_common)

  end subroutine parameters_mod_close
  ! ---------------------------------------------------------

#include "parameters_trans.F90"

end module opt_control_parameters_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
