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

!> This module contains the definition of the data type that holds a "control function"
!! used for OCT runs. 
!!
!! In addition, the module also contains the necessary procedures to manipulate these objects.
module controlfunction_m
  use datasets_m
  use epot_m
  use filter_m
  use global_m
  use io_m
  use lalg_adv_m
  use lasers_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_m
  use messages_m
  use mix_m
  use mpi_m
  use parser_m
  use profiling_m
  use states_m
  use string_m
  use tdfunction_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public :: controlfunction_t,     &
            controlfunction_mod_init,          &
            controlfunction_mod_close,         &
            controlfunction_init,              &
            controlfunction_representation,    &
            controlfunction_mode,              &
            controlfunction_set,               &
            controlfunction_end,               &
            controlfunction_copy,              &
            controlfunction_to_h,              &
            controlfunction_to_h_val,          &
            controlfunction_write,             &
            controlfunction_mixing,            &
            controlfunction_mixing_init,       &
            controlfunction_mixing_end,        &
            controlfunction_diff,              &
            controlfunction_apply_envelope,    &
            controlfunction_set_fluence,       &
            controlfunction_set_alpha,         &
            controlfunction_set_rep,           &
            controlfunction_to_realtime,       &
            controlfunction_to_basis,          &
            controlfunction_prepare_initial,   &
            controlfunction_fluence,           &
            controlfunction_j2,                &
            controlfunction_basis_to_theta,    &
            controlfunction_theta_to_basis,    &
            controlfunction_get_theta,         &
            controlfunction_set_theta,         &
            controlfunction_randomize,         &
            controlfunction_update,            &
            controlfunction_number,            &
            controlfunction_bounds,            &
            controlfunction_dof,               &
            controlfunction_w0,                &
            controlfunction_alpha,             &
            controlfunction_targetfluence,     &
            controlfunction_filter,            &
            controlfunction_gradient,          &
            controlfunction_cosine_multiply


  integer, public, parameter ::     &
    ctr_real_time              = 1, &
    ctr_fourier_series_h       = 3, &
    ctr_zero_fourier_series_h  = 4, &
    ctr_fourier_series         = 5, &
    ctr_zero_fourier_series    = 6


  integer, parameter, public :: controlfunction_mode_none      = 0, &
                                controlfunction_mode_epsilon   = 1, &
                                controlfunction_mode_f         = 2, &
                                controlfunction_mode_phi       = 3

  !> This data type contains information that is filled when the module
  !! is initialized ("controlfunction_mod_init"), and stored while the module
  !! is in use (until "controlfunction_mod_close" is called). It is information
  !! more or less common to all control functions.
  type controlfunction_common_t
    private
    integer :: representation      = 0
    FLOAT   :: omegamax            = M_ZERO
    FLOAT   :: targetfluence       = M_ZERO
    logical :: fix_initial_fluence = .false.
    FLOAT   :: w0                  = M_ZERO
    integer :: mode                = controlfunction_mode_none
    integer :: no_controlfunctions = 0
    FLOAT,       pointer :: alpha(:) => NULL()
    type(tdf_t), pointer :: td_penalty(:) => NULL()


    type(tdf_t)    :: f ! This is the envelope of the laser field, only used in the phase-only
                        ! optimization (necessary to store it in order to calculate fluences)

  end type controlfunction_common_t

  !> This is the data type used to hold a control function.
  type controlfunction_t
    private
    integer :: no_controlfunctions = 0           ! In fact, not only one control function may be used
                                                 ! to control the propagation, but several. This is the variable that holds
                                                 ! this number.
    integer :: dim           = 0                 ! If the control function is not represented directly in real time, but through
                                                 ! a number of control functions, it will actually be expanded in a basis set. This is the 
                                                 ! the dimension of the basis set. However, this does not mean necessarily that the 
                                                 ! parameters are the coefficients of the basis set.
    integer :: dof           = 0                 ! This is the number of degrees of freedom, or number of parameters, used to represent
                                                 ! a control function (this may be different -- smaller -- than "dim").
    FLOAT   :: intphi        = M_ZERO
    type(tdf_t), pointer :: f(:) => NULL()
    FLOAT, pointer :: alpha(:) => NULL()

    integer :: current_representation = 0

    FLOAT   :: w0       = M_ZERO
    FLOAT, pointer :: u(:, :) => NULL()
    FLOAT, pointer :: utransf(:, :) => NULL()
    FLOAT, pointer :: utransfi(:, :) => NULL()

    FLOAT, pointer :: theta(:) => NULL()
  end type controlfunction_t
  
  ! the next variable has to be a pointer to avoid a bug in the IBM compiler
  ! and it can not be properly initialized thanks to a bug in the PGI compiler
  logical                                 :: cf_common_initialized=.false.
  type(controlfunction_common_t), pointer :: cf_common => NULL()
  type(mix_t) :: controlfunction_mix

contains



  elemental subroutine controlfunction_common_nullify(this)
    type(controlfunction_common_t), intent(out) :: this
    !
    this%representation      = 0
    this%omegamax            = M_ZERO
    this%targetfluence       = M_ZERO
    this%fix_initial_fluence = .false.
    this%w0                  = M_ZERO
    this%mode                = controlfunction_mode_none
    this%no_controlfunctions = 0
    this%alpha               =>NULL()
    this%td_penalty          =>NULL()
    !call tdf_nullify(this%f)
    return
  end subroutine controlfunction_common_nullify

  !> Initializes the module, should be the first subroutine to be called (the last one
  !! should be controlfunction_mod_close, when the module is no longer to be used.
  !!
  !! It fills the module variable "cf_common", whose type is controlfunction_common_t, with 
  !! information obtained from the inp file.
  !!
  !! Output argument "mode_fixed_fluence" is also given a value, depending on whether
  !! the user requires a fixed-fluence run (.true.) or not (.false.).
  subroutine controlfunction_mod_init(ep, dt, max_iter, mode_fixed_fluence, parametrized_controls)
    type(epot_t), intent(inout)                   :: ep
    FLOAT, intent(in)                             :: dt
    integer, intent(in)                           :: max_iter
    logical, intent(out)                          :: mode_fixed_fluence
    logical, intent(in)                           :: parametrized_controls

    character(len=1024) :: expression
    integer :: no_lines, steps, iunit, il, idir, ncols, ipar, irow, istep
    FLOAT   :: octpenalty, f_re, f_im, total_time, time
    CMPLX   :: pol(MAX_DIM)
    type(block_t) :: blk

    PUSH_SUB(controlfunction_mod_init)

    if(.not.cf_common_initialized)then
      cf_common => NULL()
      cf_common_initialized=.true.
    end if

    if(.not. associated(cf_common)) then
      SAFE_ALLOCATE(cf_common)
      call controlfunction_common_nullify(cf_common)
    end if

    call messages_print_stress(stdout, "OCT: Info about control functions")

    !%Variable OCTControlFunctionRepresentation
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default control_fourier_series_h
    !%Description
    !% If <tt>OCTControlRepresentation = control_function_parametrized</tt>, one must 
    !% specify the kind of parameters that determine the control function.
    !% If <tt>OCTControlRepresentation = control_function_real_time</tt>, then this variable
    !% is ignored, and the control function is handled directly in real time.
    !%Option control_fourier_series_h 3
    !% The control function is expanded as a full Fourier series (although it must, of 
    !% course, be a real function). Then, the total fluence is fixed, and a transformation
    !% to hyperspherical coordinates is done; the parameters to optimize are the hyperspherical
    !% angles.
    !%Option control_zero_fourier_series_h 4
    !% The control function is expanded as a Fourier series, but assuming (1) that the zero
    !% frequency component is zero, and (2) the control function, integrated in time, adds
    !% up to zero (this essentially means that the sum of all the cosine coefficients is zero).
    !% Then, the total fluence is fixed, and a transformation to hyperspherical coordinates is 
    !% done; the parameters to optimize are the hyperspherical angles.
    !%Option control_fourier_series 5
    !% The control function is expanded as a full Fourier series (although it must, of 
    !% course, be a real function). The control parameters are the coefficients of this
    !% basis-set expansion.
    !%Option control_zero_fourier_series 6
    !% The control function is expanded as a full Fourier series (although it must, of 
    !% course, be a real function). The control parameters are the coefficients of this
    !% basis-set expansion. The difference with the option <tt>control_fourier_series</tt> is that
    !% (1) that the zero-frequency component is zero, and (2) the control function, integrated 
    !% in time, adds up to zero (this essentially means that the sum of all the cosine 
    !% coefficients is zero).
    !%End


    if (parametrized_controls) then
      call parse_integer(datasets_check('OCTControlFunctionRepresentation'), &
        ctr_fourier_series_h, cf_common%representation)
      if(.not.varinfo_valid_option('OCTControlFunctionRepresentation', cf_common%representation)) &
        call input_error('OCTControlFunctionRepresentation')
      select case(cf_common%representation)
      case(ctr_fourier_series_h)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series,'
        write(message(2), '(a)') '      and then a transformation to hyperspherical coordinates will be made.'
        call messages_info(2)
      case(ctr_zero_fourier_series_h)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series,'
        write(message(2), '(a)') '      in which (i) the zero-frequency component is assumed to be zero,'
        write(message(3), '(a)') '      and  (ii) the sum of all the cosine coefficients are zero, so that'
        write(message(4), '(a)') '      the control function starts and ends at zero.'
        write(message(5), '(a)') '      Then, a transformation to hyperspherical coordinates will be made.'
        call messages_info(6)
      case(ctr_fourier_series)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series.'
        call messages_info(1)
      case(ctr_zero_fourier_series)
        write(message(1), '(a)') 'Info: The OCT control functions will be represented as a Fourier series,'
        write(message(2), '(a)') '      in which the zero-frequency component is assumed to be zero,'
        write(message(3), '(a)') '      and  (ii) the sum of all the cosine coefficients are zero, so that'
        write(message(4), '(a)') '      the control function starts and ends at zero.'
        call messages_info(4)
      end select
    else
      cf_common%representation = ctr_real_time
      write(message(1), '(a)') 'Info: The OCT control functions will be represented in real time.'
      call messages_info(1)
    end if

    !%Variable OCTControlFunctionOmegaMax
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default -1.0
    !%Description
    !% The Fourier series that can be used to represent the control functions must be truncated;
    !% the truncation is given by a cut-off frequency which is determined by this variable.
    !%End
    call parse_float(datasets_check('OCTControlFunctionOmegaMax'), -M_ONE, cf_common%omegamax)
    if(cf_common%representation .ne. ctr_real_time) then
      write(message(1), '(a)')         'Info: The representation of the OCT control parameters will be restricted'
      write(message(2), '(a,f10.5,a)') '      with an energy cut-off of ', &
        units_from_atomic(units_out%energy, cf_common%omegamax), ' ['//trim(units_abbrev(units_out%energy)) // ']'
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
    !% first the code applies the envelope provided by the <tt>OCTLaserEnvelope</tt> input
    !% option, and afterwards it calculates the fluence.
    !%End
    call parse_float(datasets_check('OCTFixFluenceTo'), M_ZERO, cf_common%targetfluence)

    !%Variable OCTFixInitialFluence
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default yes
    !%Description
    !% By default, when asking for a fixed-fluence optimization (<tt>OCTFixFluenceTo = whatever</tt>), 
    !% the initial laser guess provided in the input file is scaled to match this
    !% fluence. However, you can force the program to use that initial laser as the initial
    !% guess, no matter the fluence, by setting <tt>OCTFixInitialFluence = no</tt>.
    !%End
    call parse_logical(datasets_check('OCTFixInitialFluence'), .true., &
      cf_common%fix_initial_fluence)

    !%Variable OCTControlFunctionType
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default controlfunction_mode_epsilon
    !%Description
    !% The control function may fully determine the time-dependent form of the 
    !% external field, or only the envelope function of this external field, or its phase. 
    !% Or, we may have two different control functions, one of them providing the phase 
    !% and the other one, the envelope.
    !%
    !% Note that, if <tt>OCTControlRepresentation = control_function_real_time</tt>, then the control
    !% function must <b>always</b> determine the full external field.
    !%Option controlfunction_mode_epsilon   1
    !% In this case, the control function determines the full control function: namely,
    !% if we are considering the electric field of a laser, the time-dependent electric field.
    !%Option controlfunction_mode_f         2
    !% The optimization process attempts to find the best possible envelope. The full 
    !% control field is this envelope times a cosine function with a "carrier" frequency. 
    !% This carrier frequency is given by the carrier frequency of the <tt>TDExternalFields</tt> 
    !% in the <tt>inp</tt> file.
    !%Option controlfunction_mode_phi       3
    !% The optimization process attempts to find the best possible time-dependent phase. That is,
    !% the external field would be given by a function in the form e(t) = f(t)*cos(w0*t+phi(t)), 
    !% where f(t) is an "envelope", w0 a carrier frequency, and phi(t) the td phase that we 
    !% wish to optimize.
    !%End
    call parse_integer(datasets_check('OCTControlFunctionType'), controlfunction_mode_epsilon, cf_common%mode)
    if(.not.varinfo_valid_option('OCTControlFunctionType', cf_common%mode)) &
      call input_error('OCTControlFunctionType')
    if ( (.not.parametrized_controls)  .and.  (cf_common%mode .ne. controlfunction_mode_epsilon) ) &
      call input_error('OCTControlFunctionType')
    call messages_print_var_option(stdout, 'OCTControlFunctionType', cf_common%mode)


    ! Check that there are no complex polarization vectors.
    do il = 1, ep%no_lasers
      pol(1:MAX_DIM) = laser_polarization(ep%lasers(il))
      do idir = 1, MAX_DIM
        if( aimag(pol(idir))**2 > CNST(1.0e-20) ) then
          write(message(1), '(a)') 'In QOCT runs, the polarization vector cannot be complex. Complex'
          write(message(2), '(a)') 'polarization vectors are only truly necessary if one wants a'  
          write(message(3), '(a)') 'circularly / elliptically polarized laser. This concepts assumes'
          write(message(4), '(a)') 'the existence of a well defined carrier frequency (otherwise it'
          write(message(5), '(a)') 'would not make sense to speak of a fixed phase difference). So in'
          write(message(6), '(a)') 'QOCT runs it would only make sense for envelope-only optimizations.'
          write(message(7), '(a)') 'This possibility should be implemented in the future.'
          call messages_fatal(7)
        end if
      end do
    end do


    ! The laser field is defined by "td functions", as implemented in module "tdfunction_m". At this point, they
    ! can be in "non-numerical" representation (i.e. described with a set of parameters, e.g. frequency, 
    ! width, etc). We need them to be in numerical form (i.e. time grid, values at the time grid). 
    ! Here we do the transformation.
    ! It cannot be done before calling controlfunction_mod_init because we need to pass the omegamax value.
    do il = 1, ep%no_lasers
      select case(cf_common%mode)
      case(controlfunction_mode_epsilon)
        call laser_to_numerical_all(ep%lasers(il), dt, max_iter, cf_common%omegamax)
      case default
        call laser_to_numerical(ep%lasers(il), dt, max_iter, cf_common%omegamax)
      end select
    end do

    ! For phase-only optimization, we need to store the envelope, in order to be able
    ! to calculate the fluence.
    if(cf_common%mode .eq. controlfunction_mode_phi) then
      call laser_get_f(ep%lasers(1), cf_common%f)
    end if 

    ! Fix the carrier frequency
    call messages_obsolete_variable('OCTCarrierFrequency')
    cf_common%w0 = laser_carrier_frequency(ep%lasers(1))

    ! Fix the number of control functions: if we have "traditional" QOCT (i.e. the control functions
    ! are represented directly in real time, then the number of control functions can be larger than
    ! one; it will be the number of lasers found in the input file. Otherwise, if the control function(s)
    ! are parametrized ("OCTControlRepresentation = control_function_parametrized"), we only have one
    ! control function. If there is more than one laser field in the input file, the program stops.
    if(parametrized_controls) then
      if(ep%no_lasers > 1) then
        write(message(1), '(a)') 'If "OCTControlRepresentation = control_function_parametrized", you can'
        write(message(2), '(a)') 'have only one external field in the input file.'
        call messages_fatal(2)
      end if
      cf_common%no_controlfunctions = 1
    else
      cf_common%no_controlfunctions = ep%no_lasers
    end if

    mode_fixed_fluence = .false.
    select case(cf_common%representation)
    case(ctr_fourier_series_h, ctr_zero_fourier_series_h)
      if(cf_common%targetfluence .eq. M_ZERO) then
        write(message(1), '(a)') 'Error: If you set "OCTControlFunctionRepresentation" to either'
        write(message(2), '(a)') '       "control_fourier_series_h", or "control_zero_fourier_series_h", then the run'
        write(message(3), '(a)') '       must be done in fixed fluence mode.'
        call messages_fatal(3)
      end if
      mode_fixed_fluence = .true.
    case(ctr_fourier_series, ctr_zero_fourier_series)
      if(cf_common%targetfluence .ne. M_ZERO) then
        write(message(1), '(a)') 'Error: If you set "OCTControlFunctionRepresentation" to "control_fourier_series",'
        write(message(2), '(a)') '       then you cannot run in fixed fluence mode.'
        call messages_fatal(2)
      end if
      mode_fixed_fluence = .false.
    case default
      if (cf_common%targetfluence .ne. M_ZERO) mode_fixed_fluence = .true.
    end select


    !%Variable OCTPenalty
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 1.0
    !%Description
    !% The variable specifies the value of the penalty factor for the 
    !% integrated field strength (fluence). Large value = small fluence.
    !% A transient shape can be specified using the block <tt>OCTLaserEnvelope</tt>.
    !% In this case <tt>OCTPenalty</tt> is multiplied with time-dependent function. 
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
    SAFE_ALLOCATE(cf_common%alpha(1:cf_common%no_controlfunctions))
    cf_common%alpha = M_ZERO
    if(parse_block('OCTPenalty', blk) == 0) then
      ! We have a block
      ncols = parse_block_cols(blk, 0)
      if(ncols .ne. cf_common%no_controlfunctions) then
        call input_error('OCTPenalty')
      else
        do ipar = 1, ncols
          call parse_block_float(blk, 0, ipar - 1, cf_common%alpha(ipar))
          if(cf_common%alpha(ipar) <= M_ZERO) call input_error('OCTPenalty')
        end do
      end if
    else
      ! We have the same penalty for all the control functions.
      call parse_float(datasets_check('OCTPenalty'), M_ONE, octpenalty)
      cf_common%alpha(1:cf_common%no_controlfunctions) = octpenalty
    end if


    !%Variable OCTLaserEnvelope
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% Often a pre-defined time-dependent envelope on the control function is desired. 
    !% This can be achieved by making the penalty factor time-dependent. 
    !% Here, you may specify the required time-dependent envelope.
    !%
    !% It is possible to choose different envelopes for different control functions.
    !% There should be one line for each control function. Each line should
    !% have only one element: a string with the function that defines the
    !% <b>inverse</b> of the time-dependent penalty, which is then defined as
    !% 1 divided by (this function + 1.0e-7) (to avoid possible singularities).
    !%
    !% The usual choices should be functions between zero and one.
    !%
    !% If, instead of defining a function, the string is <tt>default</tt>, then
    !% the program will use the function:
    !%
    !% <math> \frac{1}{\alpha(t)} = \frac{1}{2}( erf((100/T)*(t-T/20))+ erf(-(100/T)*(t-T+T/20)) </math>
    !%End
    steps = max_iter
    SAFE_ALLOCATE(cf_common%td_penalty(1:cf_common%no_controlfunctions))
    do ipar = 1, cf_common%no_controlfunctions
      call tdf_init_numerical(cf_common%td_penalty(ipar), steps, dt, -M_ONE, initval = M_ONE)
    end do

    if (parse_block(datasets_check('OCTLaserEnvelope'), blk)==0) then

      ! Cannot have this unless we have the "usual" controlfunction_mode_epsilon.
      if(cf_common%mode .ne. controlfunction_mode_epsilon) then
        write(message(1),'(a)') 'The block "OCTLaserEnvelope" is only compatible with the option'
        write(message(2),'(a)') '"OCTControlFunctionType = controlfunction_mode_epsilon".'
        call messages_fatal(2)
      end if

      no_lines = parse_block_n(blk)
      if(no_lines .ne. cf_common%no_controlfunctions) call input_error('OCTLaserEnvelope')

      do irow = 1, no_lines
        call parse_block_string(blk, irow - 1, 0, expression)
        total_time = steps * dt
        if(trim(expression) == 'default') then
          do istep = 1, steps + 1
            time = (istep - 1) * dt
            f_re = M_HALF * (loct_erf((CNST(100.0) / total_time) * (time - CNST(0.05) * total_time)) + &
              loct_erf(-(CNST(100.0) / total_time) * (time - total_time + CNST(0.05) * total_time)) )
            call tdf_set_numerical(cf_common%td_penalty(irow), istep, &
              TOFLOAT(M_ONE / (f_re + CNST(1.0e-7)))  )
          end do
        else
          call conv_to_C_string(expression)
          do istep = 1, steps + 1
            time = (istep - 1) * dt
            call parse_expression(f_re, f_im, "t", time, expression)
            call tdf_set_numerical(cf_common%td_penalty(irow), istep, &
              TOFLOAT(M_ONE / (f_re + CNST(1.0e-7)))  )
          end do
        end if
      end do

      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(OCT_DIR//'td_penalty', action='write' )
        do istep = 1, steps + 1
          time = (istep - 1) * dt
          write(iunit, '(f14.8)', advance='no') time
          do ipar = 1, cf_common%no_controlfunctions - 1
            write(iunit, '(es20.8e3)') M_ONE / tdf(cf_common%td_penalty(ipar), istep)
          end do
          write(iunit, '(es20.8e3)') M_ONE / tdf(cf_common%td_penalty(cf_common%no_controlfunctions), istep)
        end do
        write(iunit,'()')
        call io_close(iunit)
      end if

      call parse_block_end(blk)
    end if

    call messages_print_stress(stdout)
    POP_SUB(controlfunction_mod_init)
  end subroutine controlfunction_mod_init
  ! ---------------------------------------------------------

  elemental subroutine controlfunction_nullify(this)
    type(controlfunction_t), intent(out) :: this
    !
    this%no_controlfunctions    = 0
    this%dim                    = 0
    this%dof                    = 0
    this%intphi                 = M_ZERO
    this%f                      =>NULL()
    this%alpha                  =>NULL()
    this%current_representation = 0
    this%w0                     = M_ZERO
    this%u                      =>NULL()
    this%utransf                =>NULL()
    this%utransfi               =>NULL()
    this%theta                  =>NULL()
    return
  end subroutine controlfunction_nullify
  
  !> Before using an controlfunction_t variable, it needs
  !! to be initialized, either by calling controlfunction_init, or
  !! by copying another initialized variable through
  !! controlfunction_copy.
  subroutine controlfunction_init(cp, dt, ntiter)
    type(controlfunction_t), intent(inout) :: cp
    FLOAT, intent(in) :: dt
    integer, intent(in) :: ntiter

    integer :: ipar

    PUSH_SUB(controlfunction_init)

    call controlfunction_nullify(cp)

    cp%w0                  = cf_common%w0
    cp%no_controlfunctions = cf_common%no_controlfunctions
    cp%current_representation = ctr_real_time
    call loct_pointer_copy(cp%alpha, cf_common%alpha)

    SAFE_ALLOCATE(cp%f(1:cp%no_controlfunctions))
    do ipar = 1, cp%no_controlfunctions
      call tdf_init_numerical(cp%f(ipar), ntiter, dt, cf_common%omegamax)
    end do

    ! If the control function is represented directly in real time, the "dimension" (cp%dim) is
    ! the number of values that represent the function on the discretized time-axis.
    !
    ! If the control function is parametrized, up to now (in the future this might change), all 
    ! parametrizations are based on a previous basis-set expansion (sine-Fourier series, or "normal"
    ! Fourier series with or without the zero term). For the representations whose name ends in "_h", 
    ! the parameters are not directly the coefficients of the control function in this basis-set 
    ! expansion, but are constructed from them (e.g. by performing a coordinate transformation to 
    ! hyperspherical coordinates). The "dimension" (cp%dim) is the dimension of this basis set.
    select case(cf_common%representation)
    case(ctr_real_time)
      cp%dim = ntiter + 1
    case(ctr_fourier_series_h, ctr_fourier_series)
      ! If nf is the number of frequencies, we will have nf-1 non-zero "sines", nf-1 non-zero "cosines",
      ! and the zero-frequency component. Total, 2*(nf-1)+1
      cp%dim = 2 * (tdf_nfreqs(cp%f(1)) - 1) + 1
    case(ctr_zero_fourier_series, ctr_zero_fourier_series_h)
      ! If nf is the number of frequencies, we will have nf-1 non-zero "sines", nf-1 non-zero "cosines",
      ! but no zero-frequency component. Total, 2*(nf-1)
      cp%dim = 2 * (tdf_nfreqs(cp%f(1)) - 1)
    case default
      message(1) = "Internal error: invalid representation."
      call messages_fatal(1)
    end select

    ! The "degrees of freedom" cp%dof is the number of parameters that define the control function.
    ! (if it is represented directly in real time, this would be meaningless, but we put the number of 
    ! control functions, times the "dimension", which in this case is the number of time discretization 
    ! points). This is not equal to the dimension of the basis set employed (cp%dim), because we may 
    ! add further constraints, and do a coordinate transformation to account for them.
    select case(cf_common%representation)
    case(ctr_real_time)
      cp%dof = cp%no_controlfunctions * cp%dim
    case(ctr_fourier_series_h)
      ! The number of degrees of freedom is one fewer than the number of basis coefficients, since we
      ! add the constraint of fixed fluence.
      cp%dof = cp%dim - 1
    case(ctr_zero_fourier_series_h)
      ! The number of degrees of freedom is one fewer than the number of basis coefficients, since we
      ! add (1) the constraint of fixed fluence, and (2) the constraint of the field starting and
      ! ending at zero, which amounts to having all the cosine coefficients summing up to zero.
      cp%dof = cp%dim - 2
    case(ctr_fourier_series)
      ! In this case, we have no constraints: the dof is equal to the dimension of the basis set, since
      ! the parameters are directly the coefficients of the basis-set expansion.
      cp%dof = cp%dim
    case(ctr_zero_fourier_series)
      ! The number of degrees of freedom is reduced by one, since we add the constraint forcing the
      ! the field to start and end at zero, which amounts to having all the cosine coefficients 
      ! summing up to zero.
      cp%dof = cp%dim - 1
    end select


    if(cp%dof <= 0) then
      write(message(1),'(a)') 'Error: The number of degrees of freedom used to describe the control function'
      write(message(2),'(a)') '       is less than or equal to zero. This should not happen. Please review your input file.'
      call messages_fatal(2)
    else
      if(cf_common%representation .ne. ctr_real_time) then
        write(message(1), '(a)')      'Info: The expansion of the control functions in a Fourier series'
        write(message(2), '(a,i6,a)') '      expansion implies the use of ', cp%dim, ' basis-set functions.'
        write(message(3), '(a,i6,a)') '      The number of degrees of freedom is ', cp%dof,'.'
        call messages_info(3)

        SAFE_ALLOCATE(cp%theta(1:cp%dof))
        cp%theta = M_ZERO
      end if
    end if

    ! Construct the transformation matrix, if needed.
    call controlfunction_trans_matrix(cp)

    POP_SUB(controlfunction_init)
  end subroutine controlfunction_init
  ! ---------------------------------------------------------



  !> The external fields defined in epot_t "ep" are transferred to
  !! the control functions described in "cp". This should have been
  !! initialized previously.
  subroutine controlfunction_set(cp, ep)
    type(controlfunction_t), intent(inout) :: cp
    type(epot_t), intent(in) :: ep

    integer :: ipar

    PUSH_SUB(controlfunction_set)

    select case(cf_common%mode)
    case(controlfunction_mode_epsilon, controlfunction_mode_f)
      do ipar = 1, cp%no_controlfunctions
        call tdf_end(cp%f(ipar))
        call laser_get_f(ep%lasers(ipar), cp%f(ipar))
      end do
    case(controlfunction_mode_phi)
      call tdf_end(cp%f(1))
      call laser_get_phi(ep%lasers(1), cp%f(1))
    end select

    POP_SUB(controlfunction_set)
  end subroutine controlfunction_set
  ! ---------------------------------------------------------


  !> Returns the representation type for the control functions used in the OCT run.
  integer pure function controlfunction_representation()
    controlfunction_representation = cf_common%representation
  end function controlfunction_representation
  ! ---------------------------------------------------------


  !> Returns the "mode" of the control function, i.e. if it is the full pulse, the envelope, 
  !! or the phase.
  integer pure function controlfunction_mode()
    controlfunction_mode = cf_common%mode
  end function controlfunction_mode
  ! ---------------------------------------------------------


  !> "Prepares" the initial guess control field: maybe it has to be normalized to
  !! a certain fluence, maybe it should be randomized, etc.
  subroutine controlfunction_prepare_initial(par)
    type(controlfunction_t), intent(inout) :: par

    FLOAT   :: dt
    integer :: ntiter

    PUSH_SUB(controlfunction_prepare_initial)

    call controlfunction_apply_envelope(par)

    if(cf_common%targetfluence .ne. M_ZERO) then
      if(cf_common%targetfluence < M_ZERO) then
        cf_common%targetfluence = controlfunction_fluence(par) 
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with the same'
        write(message(2), '(a,f10.5,a)') '      fluence as the input external fields: F = ', &
          cf_common%targetfluence, ' a.u.'
      else
        write(message(1), '(a)')         'Info: The QOCT run will attempt to find a solution with a predefined'
        write(message(2), '(a,f10.5,a)') '      fluence: F = ', cf_common%targetfluence, ' a.u.'
      end if
      call messages_info(2)
      if(cf_common%fix_initial_fluence) call controlfunction_set_fluence(par)
    end if

    ! Now we have to find the "fluence" of the phase, in order to keep it constant.
    select case(cf_common%mode)
    case(controlfunction_mode_phi)
      par%intphi = tdf_dot_product(par%f(1), par%f(1))
      if(par%intphi <= M_ZERO) then
        dt = tdf_dt(par%f(1))
        ntiter = tdf_niter(par%f(1))

        par%intphi = CNST(0.1) * (M_PI / M_TWO)**2 * dt * ntiter
      else
        par%intphi = tdf_dot_product(par%f(1), par%f(1))
      end if
    end select

    ! Move to the "native" representation, if necessary.
    call controlfunction_set_rep(par)

    POP_SUB(controlfunction_prepare_initial)
  end subroutine controlfunction_prepare_initial
  ! ---------------------------------------------------------


  !> Transforms the control function to frequency space, if
  !! this is the space in which the functions are defined (and it
  !! is necessary to perform the transformation). 
  !! And, transforms the control function to real-time space, if
  !! this is the space in which the functions are defined (and it
  !! is necessary to perform the transformation). 
  subroutine controlfunction_set_rep(par)
    type(controlfunction_t), intent(inout) :: par

    PUSH_SUB(controlfunction_set_rep)

    if(par%current_representation .ne. cf_common%representation) then
      if(cf_common%representation .eq. ctr_real_time) then
        call controlfunction_to_realtime(par)
      else
        call controlfunction_to_basis(par)
      end if
    end if

    POP_SUB(controlfunction_set_rep)
  end subroutine controlfunction_set_rep
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_to_basis(par)
    type(controlfunction_t), intent(inout) :: par

    integer :: ipar

    PUSH_SUB(controlfunction_to_basis)

    if(par%current_representation.eq.ctr_real_time) then
      select case(cf_common%representation)
      case(ctr_fourier_series_h)
        do ipar = 1, par%no_controlfunctions
          call tdf_numerical_to_fourier(par%f(ipar))
        end do
        par%current_representation = ctr_fourier_series_h
        call controlfunction_basis_to_theta(par)
      case(ctr_zero_fourier_series_h)
        do ipar = 1, par%no_controlfunctions
          call tdf_numerical_to_zerofourier(par%f(ipar))
        end do
        par%current_representation = ctr_zero_fourier_series_h
        call controlfunction_basis_to_theta(par)
      case(ctr_fourier_series)
        do ipar = 1, par%no_controlfunctions
          call tdf_numerical_to_fourier(par%f(ipar))
        end do
        par%current_representation = ctr_fourier_series
        call controlfunction_basis_to_theta(par)
      case(ctr_zero_fourier_series)
        do ipar = 1, par%no_controlfunctions
          call tdf_numerical_to_zerofourier(par%f(ipar))
        end do
        par%current_representation = ctr_zero_fourier_series
        call controlfunction_basis_to_theta(par)
      end select
    end if

    POP_SUB(controlfunction_to_basis)
  end subroutine controlfunction_to_basis
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_to_realtime(par)
    type(controlfunction_t), intent(inout) :: par

    integer :: ipar

    PUSH_SUB(controlfunction_to_realtime)

    select case(par%current_representation)
    case(ctr_real_time)
      POP_SUB(controlfunction_to_realtime)
      return
    case(ctr_fourier_series_h)
      call controlfunction_theta_to_basis(par)
      do ipar = 1, par%no_controlfunctions
        call tdf_fourier_to_numerical(par%f(ipar))
      end do
    case(ctr_zero_fourier_series_h)
      call controlfunction_theta_to_basis(par)
      do ipar = 1, par%no_controlfunctions
        call tdf_zerofourier_to_numerical(par%f(ipar))
      end do
    case(ctr_fourier_series)
      call controlfunction_theta_to_basis(par)
      do ipar = 1, par%no_controlfunctions
        call tdf_fourier_to_numerical(par%f(ipar))
      end do
    case(ctr_zero_fourier_series)
      call controlfunction_theta_to_basis(par)
      do ipar = 1, par%no_controlfunctions
        call tdf_zerofourier_to_numerical(par%f(ipar))
      end do
    end select

    par%current_representation = ctr_real_time
    POP_SUB(controlfunction_to_realtime)
  end subroutine controlfunction_to_realtime
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function controlfunction_diff(pp, qq) result(res)
    type(controlfunction_t), intent(in) :: pp, qq

    integer :: ipar

    PUSH_SUB(controlfunction_diff)

    ASSERT(pp%current_representation .eq. qq%current_representation)

    res = M_ZERO
    do ipar = 1, pp%no_controlfunctions
      res = res + tdf_diff(pp%f(ipar), qq%f(ipar))
    end do

    POP_SUB(controlfunction_diff)
  end function controlfunction_diff
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function controlfunction_dotp(xx, yy) result(res)
    FLOAT, intent(in) :: xx(:)
    FLOAT, intent(in) :: yy(:)

    PUSH_SUB(controlfunction_dotp)
    res = sum(xx(:) * yy(:))

    POP_SUB(controlfunction_dotp)
  end function controlfunction_dotp
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_mixing_init(par)
    type(controlfunction_t), intent(in) :: par
    PUSH_SUB(controlfunction_mixing_init)
    call mix_init(controlfunction_mix, par%dim, par%no_controlfunctions, 1)
    POP_SUB(controlfunction_mixing_init)
  end subroutine controlfunction_mixing_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_mixing_end
    PUSH_SUB(controlfunction_mixing_end)
    call mix_end(controlfunction_mix)
    POP_SUB(controlfunction_mixing_end)
  end subroutine controlfunction_mixing_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_mixing(iter, par_in, par_out, par_new)
    integer, intent(in) :: iter
    type(controlfunction_t), intent(in) :: par_in, par_out
    type(controlfunction_t), intent(inout) :: par_new

    integer :: ipar, idir, dim
    FLOAT, allocatable :: e_in(:, :, :), e_out(:, :, :), e_new(:, :, :)
    PUSH_SUB(controlfunction_mixing)

    dim = par_in%dim
    SAFE_ALLOCATE(e_in (1:dim, 1:par_in%no_controlfunctions, 1:1))
    SAFE_ALLOCATE(e_out(1:dim, 1:par_in%no_controlfunctions, 1:1))
    SAFE_ALLOCATE(e_new(1:dim, 1:par_in%no_controlfunctions, 1:1))

    do ipar = 1, par_in%no_controlfunctions
      do idir = 1, dim
        e_in (idir, ipar, 1) = tdf(par_in%f(ipar), idir)
        e_out(idir, ipar, 1) = tdf(par_out%f(ipar), idir)
      end do
    end do

    e_new = M_ZERO
    call dmixing(controlfunction_mix, iter, e_in, e_out, e_new, controlfunction_dotp)
    do ipar = 1, par_out%no_controlfunctions
      call tdf_set_numerical(par_new%f(ipar), e_new(:, ipar, 1))
    end do

    SAFE_DEALLOCATE_A(e_in)
    SAFE_DEALLOCATE_A(e_out)
    SAFE_DEALLOCATE_A(e_new)
    POP_SUB(controlfunction_mixing)
  end subroutine controlfunction_mixing
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_apply_envelope(cp)
    type(controlfunction_t), intent(inout) :: cp
    integer :: ipar, iter

    PUSH_SUB(controlfunction_apply_envelope)

    ! Do not apply the envelope if the control functions are represented as a sine-Fourier series.
    if(cf_common%representation .eq. ctr_real_time) then
      do ipar = 1, cp%no_controlfunctions
        do iter = 1, tdf_niter(cp%f(ipar)) + 1
          call tdf_set_numerical(cp%f(ipar), iter, tdf(cp%f(ipar), iter) / tdf(cf_common%td_penalty(ipar), iter) )
        end do
      end do
    end if

    POP_SUB(controlfunction_apply_envelope)
  end subroutine controlfunction_apply_envelope
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_to_h(cp, ep)
    type(controlfunction_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep

    integer :: ipar
    type(controlfunction_t) :: par
    PUSH_SUB(controlfunction_to_h)

    call controlfunction_copy(par, cp)
    call controlfunction_to_realtime(par)

    select case(cf_common%mode)
    case(controlfunction_mode_epsilon, controlfunction_mode_f)
      do ipar = 1, cp%no_controlfunctions
        call laser_set_f(ep%lasers(ipar), par%f(ipar))
      end do
    case(controlfunction_mode_phi)
      call laser_set_phi(ep%lasers(1), par%f(1))
    end select

    call controlfunction_end(par)
    POP_SUB(controlfunction_to_h)
  end subroutine controlfunction_to_h
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_to_h_val(cp, ep, val)
    type(controlfunction_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep
    integer, intent(in) :: val

    integer :: ipar

    PUSH_SUB(controlfunction_to_h_val)

    do ipar = 1, cp%no_controlfunctions
      call laser_set_f_value(ep%lasers(ipar), val, tdf(cp%f(ipar), val) )
    end do

    POP_SUB(controlfunction_to_h_val)
  end subroutine controlfunction_to_h_val
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_end(cp)
    type(controlfunction_t), intent(inout) :: cp
    integer :: ipar

    PUSH_SUB(controlfunction_end)

    do ipar = 1, cp%no_controlfunctions
      call tdf_end(cp%f(ipar))
    end do
    SAFE_DEALLOCATE_P(cp%f)
    SAFE_DEALLOCATE_P(cp%alpha)
    SAFE_DEALLOCATE_P(cp%u)
    SAFE_DEALLOCATE_P(cp%utransf)
    SAFE_DEALLOCATE_P(cp%utransfi)
    SAFE_DEALLOCATE_P(cp%theta)

    POP_SUB(controlfunction_end)
  end subroutine controlfunction_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_write(filename, cp)
    character(len=*), intent(in) :: filename
    type(controlfunction_t), intent(in) :: cp

    integer :: iter, ipar, ifreq, iunit, niter, nfreqs, idof
    FLOAT :: time, wmax, dw, ww, wa, wb, dt
    FLOAT, allocatable :: func(:, :)
    CMPLX :: ft, ez, ezdt
    character(len=2) :: digit
    type(controlfunction_t) :: par

    if(.not.mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(controlfunction_write)

    call io_mkdir(trim(filename))

    call controlfunction_copy(par, cp)
    call controlfunction_to_realtime(par)

    iunit = io_open(trim(filename)//'/Fluence', action='write')
    write(iunit, '(a,es20.8e3)') 'Fluence = ', controlfunction_fluence(par)
    call io_close(iunit)

    niter = tdf_niter(par%f(1))
    SAFE_ALLOCATE(func(1:niter + 1, 1:cp%no_controlfunctions))

    select case(cf_common%mode)
    case(controlfunction_mode_epsilon)

      do ipar = 1, cp%no_controlfunctions
        if(cp%no_controlfunctions > 1) then
          write(digit,'(i2.2)') ipar
          iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cp', action='write')
        end if
        write(iunit,'(2a20)') '#       t [a.u]      ', '        e(t)         '
        do iter = 1, tdf_niter(par%f(ipar)) + 1
          time = (iter - 1) * tdf_dt(par%f(ipar))
          write(iunit, '(2es20.8e3)') time, tdf(par%f(ipar), iter)
          func(iter, ipar) = tdf(par%f(ipar), time)
        end do
        call io_close(iunit)
      end do

    case(controlfunction_mode_f)

      do ipar = 1, cp%no_controlfunctions
        if(cp%no_controlfunctions > 1) then
          write(digit,'(i2.2)') ipar
          iunit = io_open(trim(filename)//'/cp-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cp', action='write')
        end if
        write(iunit,'(3a20)') '#       t [a.u]      ', '        e(t)         ', '        f(t)         '
        do iter = 1, tdf_niter(par%f(ipar)) + 1
          time = (iter - 1) * tdf_dt(par%f(ipar))
          write(iunit, '(3es20.8e3)') time, tdf(par%f(ipar), time) * cos(par%w0 * time), tdf(par%f(ipar), time)
          func(iter, ipar) = tdf(par%f(ipar), time) * cos(par%w0 * time)
        end do
        call io_close(iunit)
      end do

    case(controlfunction_mode_phi)

      ! In this case, there is only one control function (for the moment)
      iunit = io_open(trim(filename)//'/cp', action='write')
      write(iunit,'(4a20)') '#       t [a.u]      ', '        e(t)         ', &
                            '         f(t)        ', '       phi(t)        ' 
      do iter = 1, tdf_niter(par%f(ipar)) + 1
        time = (iter - 1) * tdf_dt(par%f(ipar))
        write(iunit, '(4es20.8e3)') time, tdf(cf_common%f, time) * &
          cos(par%w0 * time + tdf(par%f(1), time) ), tdf(cf_common%f, time), tdf(par%f(1), time)
        func(iter, 1) = tdf(cf_common%f, time) * cos(par%w0 * time + tdf(par%f(1), time) )
      end do
      call io_close(iunit)

    end select


    !Now, the Fourier transforms.
    select case(cf_common%mode)
    case(controlfunction_mode_epsilon)

      do ipar = 1, cp%no_controlfunctions
        if(cp%no_controlfunctions > 1) then
          write(digit,'(i2.2)') ipar
          iunit = io_open(trim(filename)//'/cpw-'//digit, action='write')
        else
          iunit = io_open(trim(filename)//'/cpw', action='write')
        end if
        write(iunit,'(3a20)') '#       w [a.u]      ', '      Re[e(w)]       ', &
                              '      Im[e(w)]       '

        nfreqs = 1000
        wa = M_ZERO
        wb = M_THREE ! hard-coded to three atomic units... this should be improved.
        wmax = wb
        dw = wmax / (nfreqs - 1)
        dt = tdf_dt(par%f(1))

        do ifreq = 1, nfreqs
          ww = wa + (ifreq - 1) * dw
          ft = M_z0
          ez = M_z1
          ezdt = exp(M_zI * ww * tdf_dt(par%f(ipar)))
          do iter = 1, niter + 1
            time = (iter - 1) * dt
            ft = ft + func(iter, ipar) * ez
            ez = ez * ezdt
          end do
          ft = ft * dt
          write(iunit,'(3es20.8e3)') ww, real(ft), aimag(ft)
        end do
        call io_close(iunit)
      end do
      

    case(controlfunction_mode_f, controlfunction_mode_phi)
      iunit = io_open(trim(filename)//'/cpw', action='write')
      write(iunit,'(3a20)') '#       w [a.u]      ', '      Re[e(w)]       ', &
                            '      Im[e(w)]       '
      
      nfreqs = 1000
      wa = cp%w0 - M_THREE * cf_common%omegamax
      wb = cp%w0 + M_THREE * cf_common%omegamax
      wmax = CNST(6.0)*cf_common%omegamax
      dw = wmax/(nfreqs-1)
      dt = tdf_dt(par%f(1))

      do ifreq = 1, nfreqs
        ww = wa + (ifreq - 1) * dw
        ft = M_z0
        ez = M_z1
        ezdt = exp(M_zI * ww * tdf_dt(par%f(1)))
        do iter = 1, niter + 1
          time = (iter - 1) * dt
          ft = ft + func(iter, 1) * ez
          ez = ez * ezdt
        end do
        ft = ft * dt
        write(iunit,'(3es20.8e3)') ww, real(ft), aimag(ft)
      end do

      call io_close(iunit)
    end select

    ! Now, in case of a parametrized control function, the parameters.
    if(cf_common%representation .ne. ctr_real_time) then
      iunit = io_open(trim(filename)//'/theta', action='write')
      do idof = 1, par%dof
        write(iunit,'(i5,es20.8e3)') idof, par%theta(idof)
      end do
      call io_close(iunit)
    end if

    call controlfunction_end(par)
    POP_SUB(controlfunction_write)
  end subroutine controlfunction_write
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Gets the fluence of the laser field, defined as:
  ! controlfunction_fluence = \sum_i^{no_controlfunctions} \integrate_0^T |epsilon(t)|^2
  ! ---------------------------------------------------------
  FLOAT function controlfunction_fluence(par)
    type(controlfunction_t), intent(in) :: par
    type(controlfunction_t)             :: par_
    integer :: iter, ipar
    FLOAT :: time, fi, phi
    type(tdf_t) :: ff
    PUSH_SUB(controlfunction_fluence)

    call controlfunction_copy(par_, par)
    call controlfunction_to_realtime(par_)

    controlfunction_fluence = M_ZERO

    select case(cf_common%mode)
    case(controlfunction_mode_epsilon)
      do ipar = 1, par_%no_controlfunctions
        controlfunction_fluence = controlfunction_fluence + tdf_dot_product(par_%f(ipar), par_%f(ipar))
      end do
    case(controlfunction_mode_f)
      do ipar = 1, par%no_controlfunctions
        call tdf_init(ff)
        call tdf_copy(ff, par_%f(ipar))
        call tdf_cosine_multiply(par%w0, ff)
        controlfunction_fluence = controlfunction_fluence + tdf_dot_product(ff, ff)
        call tdf_end(ff)
      end do
    case(controlfunction_mode_phi)
      call tdf_init(ff)
      call tdf_copy(ff, par_%f(1))
      do iter = 1, tdf_niter(ff) + 1
        time = (iter - 1) * tdf_dt(ff)
        fi = tdf(cf_common%f, iter)
        phi = real(tdf(ff, iter)) 
        call tdf_set_numerical(ff, iter, fi * cos(par%w0 * time + phi))
      end do
      controlfunction_fluence = tdf_dot_product(ff, ff)
      call tdf_end(ff)
    end select

    call controlfunction_end(par_)
    POP_SUB(controlfunction_fluence)
  end function controlfunction_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function controlfunction_j2(par) result(j2)
    type(controlfunction_t), intent(in) :: par
    type(controlfunction_t)             :: par_
    integer :: iter, ipar
    FLOAT   :: time, integral, fi, phi, tdp
    type(tdf_t) :: ff

    PUSH_SUB(controlfunction_j2)

    ASSERT(par%current_representation .eq. cf_common%representation)

    call controlfunction_copy(par_, par)
    call controlfunction_to_realtime(par_)

    integral = M_ZERO
    select case(cf_common%mode)
    case(controlfunction_mode_epsilon)
      do ipar = 1, par_%no_controlfunctions
        call tdf_init(ff)
        call tdf_copy(ff, par_%f(ipar))
        do iter = 1, tdf_niter(ff) + 1
          time = (iter - 1) * tdf_dt(ff)
          fi = tdf(par_%f(ipar), iter)
          tdp = sqrt(real(tdf(cf_common%td_penalty(ipar), iter), kind=REAL_PRECISION))
          call tdf_set_numerical(ff, iter, fi * tdp)
        end do
        integral = integral + tdf_dot_product(ff, ff)
        call tdf_end(ff)
      end do
    case(controlfunction_mode_f)
      do ipar = 1, par_%no_controlfunctions
        call tdf_init(ff)
        call tdf_copy(ff, par_%f(ipar))
        do iter = 1, tdf_niter(ff) + 1
          time = (iter - 1) * tdf_dt(ff)
          fi = tdf(par_%f(ipar), iter)
          tdp = sqrt(real(tdf(cf_common%td_penalty(ipar), iter)))
          call tdf_set_numerical(ff, iter, fi * tdp * cos(par_%w0 * time))
        end do
        integral = integral + tdf_dot_product(ff, ff)
        call tdf_end(ff)
      end do
    case(controlfunction_mode_phi)
      call tdf_init(ff)
      call tdf_copy(ff, par_%f(1))
      do iter = 1, tdf_niter(ff) + 1
        time = (iter - 1) * tdf_dt(ff)
        fi = tdf(cf_common%f, iter)
        phi = real(tdf(par_%f(1), iter), kind=REAL_PRECISION)
        tdp = sqrt(real(tdf(cf_common%td_penalty(1), iter), kind=REAL_PRECISION))
        call tdf_set_numerical(ff, iter, fi * cos(par_%w0 * time + phi))
      end do
      integral = tdf_dot_product(ff, ff)
      call tdf_end(ff)
    end select

    j2 = - par_%alpha(1) * (integral - cf_common%targetfluence)

    call controlfunction_end(par_)
    POP_SUB(controlfunction_j2)
  end function controlfunction_j2
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_set_fluence(par)
    type(controlfunction_t), intent(inout) :: par
    FLOAT   :: old_fluence
    integer :: ipar

    PUSH_SUB(controlfunction_set_fluence)

    old_fluence = controlfunction_fluence(par) 
    do ipar = 1, par%no_controlfunctions
      call tdf_scalar_multiply( sqrt(cf_common%targetfluence / old_fluence), par%f(ipar) )
    end do

    POP_SUB(controlfunction_set_fluence)
  end subroutine controlfunction_set_fluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_set_alpha(par, alpha)
    type(controlfunction_t), intent(inout) :: par

    FLOAT, intent(in) :: alpha

    PUSH_SUB(controlfunction_set_alpha)

    par%alpha(:) = alpha

    POP_SUB(controlfunction_set_alpha)
  end subroutine controlfunction_set_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_copy(cp_out, cp_in)
    type(controlfunction_t), intent(inout) :: cp_out
    type(controlfunction_t), intent(in)    :: cp_in

    integer :: ipar

    PUSH_SUB(controlfunction_copy)

    call controlfunction_nullify(cp_out)

    cp_out%no_controlfunctions = cp_in%no_controlfunctions
    cp_out%dim = cp_in%dim
    cp_out%dof = cp_in%dof
    cp_out%intphi = cp_in%intphi
    cp_out%current_representation = cp_in%current_representation
    cp_out%w0 = cp_in%w0

    call loct_pointer_copy(cp_out%alpha, cp_in%alpha)
    SAFE_ALLOCATE(cp_out%f(1:cp_out%no_controlfunctions))

    do ipar = 1, cp_in%no_controlfunctions
      call tdf_init(cp_out%f(ipar))
      call tdf_copy(cp_out%f(ipar), cp_in%f(ipar))
    end do

    call loct_pointer_copy(cp_out%u, cp_in%u)
    call loct_pointer_copy(cp_out%utransf, cp_in%utransf)
    call loct_pointer_copy(cp_out%utransfi, cp_in%utransfi)
    call loct_pointer_copy(cp_out%theta, cp_in%theta)

    POP_SUB(controlfunction_copy)
  end subroutine controlfunction_copy
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_randomize(par)
    type(controlfunction_t), intent(inout) :: par

     integer :: ipar

     PUSH_SUB(controlfunction_randomize)

     ASSERT(cf_common%representation .ne. ctr_real_time)

     call controlfunction_set_rep(par)

     select case(cf_common%mode)
     case(controlfunction_mode_epsilon, controlfunction_mode_f)
       do ipar = 1, par%no_controlfunctions
         call tdf_set_random(par%f(ipar))
       end do
     case(controlfunction_mode_phi)
       call tdf_set_random(par%f(1))
     end select

     POP_SUB(controlfunction_randomize)
  end subroutine controlfunction_randomize
  ! ---------------------------------------------------------



  !> Update the control function(s) given in "cp", according to the formula
  !! cp = (1 - mu) * cpp + mu * dd / (td_penalty - 2 * dq)
  subroutine controlfunction_update(cp, cpp, dir, iter, mu, dd, dq)
    type(controlfunction_t), intent(inout) :: cp
    type(controlfunction_t), intent(in)    :: cpp
    character(len=1),        intent(in)    :: dir
    integer,                 intent(in)    :: iter
    FLOAT,                   intent(in)    :: mu
    FLOAT,                   intent(in)    :: dd(:)
    CMPLX,                   intent(in)    :: dq(:)

    FLOAT :: value
    integer :: ipar
    
    PUSH_SUB(controlfunction_update)

    select case(dir)
      case('f')
        do ipar = 1, cp%no_controlfunctions
          value = dd(ipar) / ( tdf(cf_common%td_penalty(ipar), iter) - M_TWO * aimag(dq(ipar)) )
          value = (M_ONE - mu) * tdf(cpp%f(ipar), iter) + mu * value
          call tdf_set_numerical(cp%f(ipar), iter, value)
          if(iter + 1 <= tdf_niter(cp%f(ipar)) + 1)  call tdf_set_numerical(cp%f(ipar), iter+1, value)
          if(iter + 2 <= tdf_niter(cp%f(ipar)) + 1)  call tdf_set_numerical(cp%f(ipar), iter+2, value)
        end do

      case('b')
        do ipar = 1, cp%no_controlfunctions
          value = dd(ipar) / ( tdf(cf_common%td_penalty(ipar), iter + 1) - M_TWO * aimag(dq(ipar)) )
          value = (M_ONE - mu) * tdf(cpp%f(ipar), iter + 1) + mu * value
          call tdf_set_numerical(cp%f(ipar), iter + 1, value)
          if(iter > 0) call tdf_set_numerical(cp%f(ipar), iter, value)
          if(iter - 1 > 0) call tdf_set_numerical(cp%f(ipar), iter-1, value)
        end do
    end select

    POP_SUB(controlfunction_update)
  end subroutine controlfunction_update
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function controlfunction_alpha(par, ipar)
    type(controlfunction_t), intent(in) :: par
    integer,                 intent(in) :: ipar
    controlfunction_alpha = par%alpha(ipar)
  end function controlfunction_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function controlfunction_targetfluence()
    controlfunction_targetfluence = cf_common%targetfluence
  end function controlfunction_targetfluence
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function controlfunction_number(par)
    type(controlfunction_t), intent(in) :: par
    controlfunction_number = par%no_controlfunctions
  end function controlfunction_number
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_bounds(par, lower_bounds, upper_bounds)
    type(controlfunction_t), intent(in)  :: par
    FLOAT,                   intent(out) :: lower_bounds(:)
    FLOAT,                   intent(out) :: upper_bounds(:)
    integer :: dog

    PUSH_SUB(controlfunction_bounds)

    upper_bounds = M_PI
    dog = controlfunction_dof(par)

    select case(cf_common%mode)
    case(controlfunction_mode_epsilon, controlfunction_mode_f, controlfunction_mode_phi)
      lower_bounds(1:dog - 1) = M_ZERO
      lower_bounds(dog)       = -M_PI
    end select

    POP_SUB(controlfunction_bounds)
  end subroutine controlfunction_bounds
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function controlfunction_dof(par)
    type(controlfunction_t), intent(in) :: par
    controlfunction_dof = par%dof
  end function controlfunction_dof
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function controlfunction_w0(par)
    type(controlfunction_t), intent(in) :: par
    controlfunction_w0 = par%w0
  end function controlfunction_w0
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_filter(par, filter)
    type(controlfunction_t), intent(inout) :: par
    type(filter_t),          intent(inout) :: filter

    integer :: ipar

    PUSH_SUB(controlfunction_filter)

    do ipar = 1, par%no_controlfunctions
      call filter_apply(par%f(ipar), filter)
    end do

    POP_SUB(controlfunction_filter)
  end subroutine controlfunction_filter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine controlfunction_mod_close()
    integer :: ipar

    PUSH_SUB(controlfunction_mod_close)

    SAFE_DEALLOCATE_P(cf_common%alpha)

    do ipar = 1, cf_common%no_controlfunctions
      call tdf_end(cf_common%td_penalty(ipar))
    end do

    SAFE_DEALLOCATE_P(cf_common%td_penalty)
    SAFE_DEALLOCATE_P(cf_common)

    POP_SUB(controlfunction_mod_close)
  end subroutine controlfunction_mod_close
  ! ---------------------------------------------------------




  !> controlfunction_gradient computes the gradient with respect
  !! to the theta basis (dim=dof)
  subroutine controlfunction_gradient(xx, par, par_output, grad)
    FLOAT,                   intent(in)    :: xx(:)
    type(controlfunction_t), intent(in)    :: par, par_output
    FLOAT,                   intent(inout) :: grad(:)

    integer :: dim, jj, mm, kk, ss, tt
    FLOAT :: rr
    FLOAT, allocatable :: theta(:), grad_matrix(:,:), gradb(:)

    PUSH_SUB(controlfunction_gradient)


    select case(par%current_representation)
    case(ctr_fourier_series)

       dim = par%dim
       SAFE_ALLOCATE(theta(1:dim)) ! dim = dof for fourier-series
       call controlfunction_get_theta(par_output, theta)
       forall(jj = 1:dim) 
         grad(jj) = M_TWO * controlfunction_alpha(par, 1) * sum(par%u(jj, :)*xx(:)) - M_TWO * theta(jj)
       end forall

    case(ctr_zero_fourier_series)

       dim = par%dof
       forall(jj = 1:dim) 
         grad(jj) = M_TWO * controlfunction_alpha(par, 1) * sum(par%u(jj, :)*xx(:))
       end forall
       SAFE_ALLOCATE(theta(1:dim+1))
       forall(jj = 1:dim+1) theta(jj) = tdf(par_output%f(1), jj)
       do jj = 1, dim/2
         grad(jj) = grad(jj) - M_TWO*(theta(jj+1) - theta(1))
       end do
       do jj = dim/2+1, dim
         grad(jj) = grad(jj) - M_TWO*theta(jj+1)
       end do


    case(ctr_fourier_series_h)

       dim = par%dim
       SAFE_ALLOCATE(theta(1:dim))   ! dim = dof + 1 for fourier-series-h
       SAFE_ALLOCATE(grad_matrix(1:dim - 1, 1:dim))

       forall(jj = 1:dim) theta(jj) = M_TWO * tdf(par_output%f(1), jj) ! get the projection on my basis set of function (theta=b)
       rr = sqrt(cf_common%targetfluence)
       call hypersphere_grad_matrix(grad_matrix, rr, xx)
       grad = matmul(grad_matrix, matmul(transpose(par%utransfi), theta))
       grad = -grad  ! the CG algorithm minimizes, so we need to give the negative gradient for maximization
       
       SAFE_DEALLOCATE_A(grad_matrix)

    case(ctr_zero_fourier_series_h)

       dim = par%dim
       SAFE_ALLOCATE(theta(1:dim))
       SAFE_ALLOCATE(gradb(1:dim-1))

       forall(jj = 1:dim) theta(jj) = tdf(par_output%f(1), jj)

       do jj = 1, dim/2-1
         gradb(jj) = - M_TWO*(theta(jj+1) - theta(1))
       end do
       do jj = dim/2, dim-1
         gradb(jj) = - M_TWO*theta(jj+1)
       end do

       SAFE_ALLOCATE(grad_matrix(1:dim - 2, 1:dim - 1))
       rr = sqrt(cf_common%targetfluence)
       call hypersphere_grad_matrix(grad_matrix, rr, xx)

       grad = matmul(grad_matrix, matmul(transpose(par%utransfi), gradb))

       SAFE_DEALLOCATE_A(grad_matrix)
       SAFE_DEALLOCATE_A(theta)
       SAFE_DEALLOCATE_A(gradb)

    end select

    SAFE_DEALLOCATE_A(theta)
    POP_SUB(controlfunction_gradient)
  end subroutine controlfunction_gradient
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Multiplies all the control function by cos(w0*t), where
  !! w0 is the carrier frequency.
  subroutine controlfunction_cosine_multiply(par)
    type(controlfunction_t), intent(inout) :: par

    integer :: i

    PUSH_SUB(controlfunction_cosine_multiply)

    call controlfunction_to_realtime(par)

    do i = 1, par%no_controlfunctions
      call tdf_cosine_multiply(par%w0, par%f(i))
    end do

    POP_SUB(controlfunction_cosine_multiply)
  end subroutine controlfunction_cosine_multiply

#include "controlfunction_trans_inc.F90"

end module controlfunction_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
