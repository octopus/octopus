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
!! $Id: opt_control_global.F90 3099 2007-07-23 14:21:35Z lorenzen $

#include "global.h"

!> This module contains the definition of the oct_t data type, which
!! contains some of the basic information about how the OCT run will
!! be done.
!! The "oct" variable (whose datatype is oct_t) is declared in the main
!! OCT module, opt_control_m (as a module variable). It is initialized
!! by calling oct_read_inp, defined below in this module.
module opt_control_global_m
  use datasets_m
  use global_m
  use messages_m
  use parser_m
  use varinfo_m

  implicit none

  !> The oct_t datatype stores the basic information about how the OCT run
  !! is done: which algorithm, how the control funtion is stored, should the
  !! intermediate results be stored for debugging, etc.
  type oct_t
    integer :: ctr_function_rep
    logical :: mode_fixed_fluence
    integer :: algorithm
    FLOAT   :: eta, delta
    FLOAT   :: direct_step
    logical :: use_mixing
    logical :: oct_double_check
    FLOAT   :: check_gradient
    logical :: dump_intermediate
    integer :: number_checkpoints
    logical :: random_initial_guess
  end type oct_t

  integer, parameter :: &
    oct_ctr_function_real_time       = 1,       &
    oct_ctr_function_parametrized    = 2

  ! These are the possible QOCT schemes or algorithms; the component "algorithm"
  ! of the oct_t datatype can get any of these values.
  integer, parameter ::  &
    oct_algorithm_zbr98              = 1,       &
    oct_algorithm_zr98               = 2,       &
    oct_algorithm_wg05               = 3,       &
    oct_algorithm_mt03               = 4,       &
    oct_algorithm_krotov             = 5,       &
    oct_algorithm_str_iter           = 6,       &
    oct_algorithm_cg                 = 7,       &
    oct_algorithm_direct             = 8,       &
    oct_algorithm_newuoa             = 9

  contains



  !> Reads, from the inp file, some global information about how the QOCT run
  !! should be. It uses this information to fill the "oct" variable. All the components
  !! of oct are filled, except for mode_fixed_fluence, which is filled when the control
  !! parameters module is initialized.
  subroutine oct_read_inp(oct)
    type(oct_t), intent(inout) :: oct
    PUSH_SUB(oct_read_inp)

    call messages_print_stress(stdout, "OCT run mode")


    !%Variable OCTControlRepresentation
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default control_function_real_time
    !%Description
    !% Optimal Control Theory can be performed with <tt>Octopus</tt> in two different modes:
    !% either considering the control function to be described in full in real time,
    !% or to be represented by a set of parameters (which may, or may not be,
    !% the coefficients of its expansion in a given basis). The particular choice
    !% for these parameters is specified by variable <tt>OCTParameterRepresentation</tt>
    !% (this variable will be ignored if the control function is to be represented
    !% directly in real time).
    !%Option control_function_real_time 1
    !% The control functions are represented directly in real time.
    !%Option control_function_parametrized 2
    !% The control functions are specified by a set of parameters.
    !%End
    call parse_integer(datasets_check('OCTControlRepresentation'), &
      oct_ctr_function_real_time, oct%ctr_function_rep)
    if(.not.varinfo_valid_option('OCTControlRepresentation', oct%ctr_function_rep)) &
      call input_error('OCTControlRepresentation')
    call messages_print_var_option(stdout, &
      'OCTControlRepresentation', oct%ctr_function_rep)


    !%Variable OCTScheme
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default oct_algorithm_zbr98
    !%Description
    !% Optimal Control Theory can be performed with <tt>Octopus</tt> with a variety of different
    !% algorithms. Not all of them can be used with any choice of target or control function
    !% representation. For example, some algorithms cannot be used if 
    !% <tt>OCTControlRepresentation = control_function_real_time</tt>    
    !% (<tt>OCTScheme</tt> .gt. <tt>oct_algorithm_straight_iteration</tt>), and others cannot be used 
    !% if <tt>OCTControlRepresentation = control_function_parametrized</tt>
    !% (<tt>OCTScheme</tt> .lt. <tt>oct_algorithm_straight_iteration</tt>).
    !%Option oct_algorithm_zbr98 1 
    !% Backward-Forward-Backward scheme described in <i>JCP</i> <b>108</b>, 1953 (1998).
    !% Only possible if target operator is a projection operator.
    !% Provides fast, stable and monotonic convergence.
    !%Option oct_algorithm_zr98  2
    !% Forward-Backward-Forward scheme described in <i>JCP</i> <b>109</b>, 385 (1998).
    !% Works for projection and more general target operators also. The convergence is 
    !% stable but slower than ZBR98. 
    !% Note that local operators show an extremely slow convergence. It ensures monotonic 
    !% convergence.
    !%Option oct_algorithm_wg05  3
    !% Forward-Backward scheme described in <i>J. Opt. B.</i> <b>7</b>, 300 (2005).
    !% Works for all kinds of target operators, can be used with all kinds of filters, and 
    !% allows a fixed fluence.
    !% The price is a rather unstable convergence. 
    !% If the restrictions set by the filter and fluence are reasonable, a good overlap can be 
    !% expected within 20 iterations.
    !% No monotonic convergence.
    !%Option oct_algorithm_mt03 4
    !% Basically an improved and generalized scheme. 
    !% Comparable to ZBR98/ZR98. See [Y. Maday and G. Turinici, <i>J. Chem. Phys.</i> <b>118</b>, 8191 (2003)].
    !%Option oct_algorithm_krotov 5
    !% The procedure reported in [D. Tannor, V. Kazakov and V.
    !% Orlov, in <i>Time-Dependent Quantum Molecular Dynamics</i>, edited by J. Broeckhove
    !% and L. Lathouweres (Plenum, New York, 1992), pp. 347-360].
    !%Option oct_algorithm_straight_iteration 6
    !% Straight iteration: one forward and one backward propagation is performed at each
    !% iteration, both with the same control field. An output field is calculated with the
    !% resulting wavefunctions. Note that this scheme typically does not converge, unless
    !% some mixing (<tt>OCTMixing = yes</tt>) is used.
    !%Option oct_algorithm_cg 7
    !% Conjugate-gradients, as implemented in the GNU GSL library.
    !%Option oct_algorithm_direct 8
    !% This is a "direct" optimization scheme. This means that we do not make use of the
    !% "usual" QOCT equations (backward-forward propagations, etc), but we use some gradient-free
    !% maximization algorithm for the function that we want to optimize. In this case, the
    !% maximization algorithm is the Nelder-Mead algorithm as implemeted in the GSL. The function
    !% values are obtained by successive forward propagations.
    !%Option oct_algorithm_newuoa 9
    !% This is exactly the same as <tt>oct_algorithm_direct</tt>, except in this case the maximization
    !% algorithm is the so-called NEWUOA algorithm [M. J. D. Powell, <i>IMA J. Numer. Analysis</i>
    !% <b>28</b>, 649-664 (2008)].
    !%End
    call parse_integer(datasets_check('OCTScheme'), oct_algorithm_zr98, oct%algorithm)
    if(.not.varinfo_valid_option('OCTScheme', oct%algorithm)) call input_error('OCTScheme')
    ! We must check that the algorithm is consistent with OCTControlRepresentation, i.e.
    ! some algorithms only make sense if the control functions are handled directly in real
    ! time, some others make only sense if the control functions are parameterized.
    if(oct%ctr_function_rep .eq. oct_ctr_function_real_time) then
      if(oct%algorithm > oct_algorithm_str_iter) call input_error('OCTScheme')
    else
      if(oct%algorithm < oct_algorithm_str_iter) call input_error('OCTScheme')
    end if
    call messages_print_var_option(stdout, "OCTScheme", oct%algorithm)
    select case(oct%algorithm)
    case(oct_algorithm_mt03)
      oct%delta = M_TWO; oct%eta = M_ZERO
    case(oct_algorithm_zr98)
      oct%delta = M_ONE; oct%eta = M_ONE
    case(oct_algorithm_krotov)
      oct%delta = M_ONE; oct%eta = M_ZERO
    case(oct_algorithm_str_iter)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(oct_algorithm_cg)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(oct_algorithm_direct)
      call parse_float(datasets_check('OCTEta'), M_ONE, oct%eta)
      call parse_float(datasets_check('OCTDelta'), M_ZERO, oct%delta)
    case(oct_algorithm_newuoa)
#if defined(HAVE_NEWUOA)
      call parse_float(datasets_check('OCTEta'), M_ONE, oct%eta)
      call parse_float(datasets_check('OCTDelta'), M_ZERO, oct%delta)
#else
      write(message(1), '(a)') '"OCTScheme = oct_algorithm_newuoa" is only possible if the newuoa'
      write(message(2), '(a)') 'code has been compiled. You must configure octopus passing the'
      write(message(3), '(a)') 'the "--enable-newuoa" switch.'
      call messages_fatal(3)
#endif
    case default
      oct%delta = M_ONE; oct%eta = M_ONE
    end select

    !%Variable OCTDoubleCheck
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default true
    !%Description 
    !% In order to make sure that the optimized field indeed does its job, the code 
    !% may run a normal propagation after the optimization using the optimized field.
    !%End
    call parse_logical(datasets_check('OCTDoubleCheck'), .true., oct%oct_double_check)
    call messages_print_var_value(stdout, "OCTDoubleCheck", oct%oct_double_check)


    !%Variable OCTCheckGradient
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 0.0
    !%Description 
    !% When doing QOCT with the conjugate-gradient optimization scheme, the gradient is
    !% computed thanks to a forward-backwards propagation. For debugging purposes, this
    !% gradient can be compared with the value obtained "numerically" (<i>i.e.</i> by doing
    !% successive forward propagations with control fields separated by small finite
    !% differences).
    !%
    !% In order to activate this feature, set <tt>OCTCheckGradient</tt> to some non-zero value,
    !% which will be the finite difference used to numerically compute the gradient.
    !%End
    call parse_float(datasets_check('OCTCheckGradient'), CNST(0.0), oct%check_gradient)
    call messages_print_var_value(stdout, "OCTCheckGradient", oct%check_gradient)


    !%Variable OCTMixing
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default false
    !%Description 
    !% Use mixing algorithms to create the input fields in the iterative OCT schemes.
    !% Note that this idea is still a little bit experimental, and depending on the
    !% kind of mixing that you use, and the parameters that you set, it may or may
    !% not accelerate the convergence, or even spoil the convergence.
    !%
    !% Using <tt>TypeOfMixing = broyden</tt>, <tt>Mixing = 0.1</tt> and <tt>MixNumberSteps = 3</tt> seems
    !% to work in many cases, but your mileage may vary.
    !%
    !% Note that mixing does not make sense (and is therefore not done, this variable
    !% being ignored), for some OCT algorithms (in particular, if <tt>OCTScheme</tt> is
    !% <tt>oct_algorithm_direct</tt> or <tt>oct_algorithm_newuoa</tt>).
    !%End
    call parse_logical(datasets_check('OCTMixing'), .false., oct%use_mixing)
    call messages_print_var_value(stdout, "OCTMixing", oct%use_mixing)

    !%Variable OCTDirectStep
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 0.25
    !%Description 
    !% If you choose <tt>OCTScheme = oct_algorithm_direct</tt> or <tt>OCTScheme = oct_algorithm_newuoa</tt>,
    !% the algorithms necessitate an initial "step" to perform the direct search for the
    !% optimal value. The precise meaning of this "step" differs.
    !%End
    call parse_float(datasets_check('OCTDirectStep'), CNST(0.25), oct%direct_step)
    call messages_print_var_value(stdout, "OCTDirectStep", oct%direct_step)

    !%Variable OCTDumpIntermediate
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default true
    !%Description 
    !% Writes to disk some data during the OCT algorithm at intermediate steps.
    !% This is rather technical and it should be considered only for debugging
    !% purposes. Nevertheless, since the whole OCT infrastructure is at a very
    !% preliminary stage of development, it is set to true by default.
    !%End
    call parse_logical(datasets_check('OCTDumpIntermediate'), .true., oct%dump_intermediate)
    call messages_print_var_value(stdout, "OCTDumpIntermediate", oct%dump_intermediate)

    !%Variable OCTNumberCheckPoints
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 0
    !%Description 
    !% During an OCT propagation, the code may write the wavefunctions at some time steps (the
    !% "check points"). When the inverse backward or forward propagation
    !% is performed in a following step, the wavefunction should reverse its path
    !% (almost) exactly. This can be checked to make sure that it is the case -- otherwise
    !% one should try reducing the time-step, or altering in some other way the
    !% variables that control the propagation.
    !%
    !% If the backward (or forward) propagation is not retracing the steps of the previous
    !% forward (or backward) propagation, the code will write a warning.
    !%End
    call parse_integer(datasets_check('OCTNumberCheckPoints'), 0, oct%number_checkpoints)
    call messages_print_var_value(stdout, "OCTNumberCheckPoints", oct%number_checkpoints)

    !%Variable OCTRandomInitialGuess
    !%Type logical
    !%Section Calculation Modes::Optimal Control
    !%Default false
    !%Description 
    !% The initial field to start the optimization search is usually given in the <tt>inp</tt> file,
    !% through a <tt>TDExternalFields</tt> block. However, you can start from a random guess if you
    !% set this variable to true.
    !%
    !% Note, however, that this is only valid for the "direct" optimization schemes; moreover
    !% you still need to provide a <tt>TDExternalFields</tt> block.
    !%End
    call parse_logical(datasets_check('OCTRandomInitialGuess'), &
      .false., oct%random_initial_guess)
    call messages_print_var_value(stdout, "OCTRandomInitialGuess", oct%random_initial_guess)

    call messages_print_stress(stdout)
    POP_SUB(oct_read_inp)
  end subroutine oct_read_inp
  ! ---------------------------------------------------------



  !> Returns .true. if the algorithm to be used is one of the "direct" or "gradient-less"
  !! algorithms -- the ones that do not require backwards propagations. Returns .false. otherwise
  logical pure function oct_algorithm_is_direct(oct)
    type(oct_t), intent(in) :: oct
    oct_algorithm_is_direct = (oct%algorithm >= oct_algorithm_direct)
  end function oct_algorithm_is_direct
  ! ---------------------------------------------------------

 
end module opt_control_global_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
