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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module contains the definition of the oct_t data type, which
!! contains some of the basic information about how the OCT run will
!! be done.
!! The "oct" variable (whose datatype is oct_t) is declared in the main
!! OCT module, opt_control_oct_m (as a module variable). It is initialized
!! by calling oct_read_inp, defined below in this module.
module opt_control_global_oct_m
  use global_oct_m
  use messages_oct_m
  use parser_oct_m
  use varinfo_oct_m

  implicit none

  private
  public :: oct_t,                    &
            oct_read_inp,             &
            oct_algorithm_is_direct  

  !> The oct_t datatype stores the basic information about how the OCT run
  !! is done: which algorithm, how the control funtion is stored, should the
  !! intermediate results be stored for debugging, etc.
  type oct_t
    ! Components are public by default
    integer :: algorithm            !< The algorithm to optimize depends on whether the control function is
                                    !! represented in real time, or is parametrized. Filled by the OCTScheme input variable.
    logical :: mode_fixed_fluence   !< Whether or not the optimization is performed in the subspace of external fields
                                    !! of equal fluence. This is filled by the controlfunction_mod_init subroutine.
    FLOAT   :: eta, delta           !< "Technical" variables, that complete the definition of some algorithms.
    FLOAT   :: direct_step          !< The "initial step" of the optimization search, used by some algorithms. Filled
                                    ! by the OCTDirectStep input variable.
    logical :: oct_double_check     !< At the end of the optimization, a final run can be performed in order to make sure
                                    !! that, indeed, the optimized field produces the optimal value.
    FLOAT   :: check_gradient       !< If using the conjugate gradients algorithm, one may make sure that the forward-backward
                                    !! propagation is indeed computing the gradient of the functional, by computing this
                                    !! gradient numerically. This is sent by the OCTCheckGradient input variable.
    integer :: number_checkpoints   !< When propagating backwards, the code may check that the evolution is preserving
                                    !! time-reversal symmetry by checking that the state is equal to a number of previously
                                    !! stored "check-points", saved during the forward propagation.
    logical :: random_initial_guess !< Can be used only with some algorithms; instead of using the field described in the input
                                    !! file as initial guess, the code may generate a random field.
  end type oct_t

contains

  !> Reads, from the inp file, some global information about how the QOCT run
  !! should be. It uses this information to fill the "oct" variable. All the components
  !! of oct are filled, except for mode_fixed_fluence, which is filled when the control
  !! parameters module is initialized.
  subroutine oct_read_inp(oct, parser)
    type(oct_t),    intent(inout) :: oct
    type(parser_t), intent(in)    :: parser

    PUSH_SUB(oct_read_inp)

    call messages_print_stress(stdout, "OCT run mode")
    call messages_obsolete_variable(parser, 'OCTControlRepresentation')

    !%Variable OCTScheme
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default OPTION__OCTSCHEME__ZBR98
    !%Description
    !% Optimal Control Theory can be performed with <tt>Octopus</tt> with a variety of different
    !% algorithms. Not all of them can be used with any choice of target or control function
    !% representation. For example, some algorithms cannot be used if 
    !% <tt>OCTControlRepresentation = control_function_real_time</tt>    
    !% (<tt>OCTScheme</tt> > <tt>oct_straight_iteration</tt>), and others cannot be used 
    !% if <tt>OCTControlRepresentation = control_function_parametrized</tt>
    !% (<tt>OCTScheme</tt>  <  <tt>oct_straight_iteration</tt>).
    !%Option oct_zbr98 1 
    !% Backward-Forward-Backward scheme described in <i>JCP</i> <b>108</b>, 1953 (1998).
    !% Only possible if target operator is a projection operator.
    !% Provides fast, stable and monotonic convergence.
    !%Option oct_zr98  2
    !% Forward-Backward-Forward scheme described in <i>JCP</i> <b>109</b>, 385 (1998).
    !% Works for projection and more general target operators also. The convergence is 
    !% stable but slower than ZBR98. 
    !% Note that local operators show an extremely slow convergence. It ensures monotonic 
    !% convergence.
    !%Option oct_wg05  3
    !% Forward-Backward scheme described in <i>J. Opt. B.</i> <b>7</b>, 300 (2005).
    !% Works for all kinds of target operators, can be used with all kinds of filters, and 
    !% allows a fixed fluence.
    !% The price is a rather unstable convergence. 
    !% If the restrictions set by the filter and fluence are reasonable, a good overlap can be 
    !% expected within 20 iterations.
    !% No monotonic convergence.
    !%Option oct_mt03 4
    !% Basically an improved and generalized scheme. 
    !% Comparable to ZBR98/ZR98. See [Y. Maday and G. Turinici, <i>J. Chem. Phys.</i> <b>118</b>, 8191 (2003)].
    !%Option oct_krotov 5
    !% The procedure reported in [D. Tannor, V. Kazakov and V.
    !% Orlov, in <i>Time-Dependent Quantum Molecular Dynamics</i>, edited by J. Broeckhove
    !% and L. Lathouweres (Plenum, New York, 1992), pp. 347-360].
    !%Option oct_straight_iteration 6
    !% Straight iteration: one forward and one backward propagation is performed at each
    !% iteration, both with the same control field. An output field is calculated with the
    !% resulting wavefunctions. 
    !%Option oct_cg 7
    !% Conjugate-gradients, as implemented in the GNU GSL library. In particular, the
    !% Fletcher-Reeves version.
    !%Option oct_bfgs 8
    !% The methods use the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.  
    !% Also, it calls the GNU GSL library version of the algorithm. It is a quasi-Newton 
    !% method which builds up an approximation to the second derivatives of the function using
    !% the difference between successive gradient vectors.  By combining the first and second 
    !% derivatives the algorithm is able to take Newton-type steps towards the function minimum, 
    !% assuming quadratic behavior in that region. We have chosen to implement the "bfgs2" version,
    !% as GSL calls it, which is supposed to be the most efficient version available, and a faithful 
    !% implementation of the line minimization scheme described in "Practical Methods of Optimization",
    !% (Fletcher), Algorithms 2.6.2 and 2.6.4.  
    !%Option oct_direct 9
    !% This is a "direct" optimization scheme. This means that we do not make use of the
    !% "usual" QOCT equations (backward-forward propagations, etc), but we use some gradient-free
    !% maximization algorithm for the function that we want to optimize. In this case, the
    !% maximization algorithm is the Nelder-Mead algorithm as implemeted in the GSL. The function
    !% values are obtained by successive forward propagations.
    !%Option oct_nlopt_bobyqa 11
    !% The BOBYQA algorithm, as implemented in the NLOPT library -- therefore, octopus has to
    !% be compiled with it in order to be able to use this option.
    !%Option oct_nlopt_lbfgs 12
    !% The local BFGS, as implemented in the NLOPT library -- therefore, octopus has to
    !% be compiled with it in order to be able to use this option.
    !%End
    call parse_variable(dummy_parser, 'OCTScheme', OPTION__OCTSCHEME__OCT_ZR98, oct%algorithm)
    if(.not.varinfo_valid_option('OCTScheme', oct%algorithm)) call messages_input_error('OCTScheme')
    ! We must check that the algorithm is consistent with OCTControlRepresentation, i.e.
    ! some algorithms only make sense if the control functions are handled directly in real
    ! time, some others make only sense if the control functions are parameterized.
    ! This check cannot be here any more, and it should be placed somewhere else.
    call messages_print_var_option(stdout, "OCTScheme", oct%algorithm)
    select case(oct%algorithm)
    case(OPTION__OCTSCHEME__OCT_MT03)
      oct%delta = M_TWO; oct%eta = M_ZERO
      !%Variable OCTEta
      !%Type float
      !%Section Calculation Modes::Optimal Control
      !%Default 1.0
      !%Description 
      !% If <tt>OCTScheme = oct_mt03</tt>, then you can supply the "eta" and "delta" parameters
      !% described in [Y. Maday and G. Turinici, <i>J. Chem. Phys.</i> <b>118</b>, 8191 (2003)], using the
      !% <tt>OCTEta</tt> and <tt>OCTDelta</tt> variables.
      !%End
      call parse_variable(dummy_parser, 'OCTEta', M_ONE, oct%eta)
      !%Variable OCTDelta
      !%Type float
      !%Section Calculation Modes::Optimal Control
      !%Default 0.0
      !%Description 
      !% If <tt>OCTScheme = oct_mt03</tt>, then you can supply the "eta" and "delta" parameters
      !% described in [Y. Maday and G. Turinici, <i>J. Chem. Phys.</i> <b>118</b>, 8191 (2003)], using the
      !% <tt>OCTEta</tt> and <tt>OCTDelta</tt> variables.
      !%End
      call parse_variable(dummy_parser, 'OCTDelta', M_ZERO, oct%delta)

    case(OPTION__OCTSCHEME__OCT_ZR98)
      oct%delta = M_ONE; oct%eta = M_ONE
    case(OPTION__OCTSCHEME__OCT_KROTOV)
      oct%delta = M_ONE; oct%eta = M_ZERO
    case(OPTION__OCTSCHEME__OCT_STRAIGHT_ITERATION)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(OPTION__OCTSCHEME__OCT_CG)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(OPTION__OCTSCHEME__OCT_BFGS)
      oct%delta = M_ZERO; oct%eta = M_ONE
    case(OPTION__OCTSCHEME__OCT_DIRECT)
      ! The use of these variables for the direct and bobyqa schemes remain undocumented for the moment.
      call parse_variable(dummy_parser, 'OCTEta', M_ONE, oct%eta)
      call parse_variable(dummy_parser, 'OCTDelta', M_ZERO, oct%delta)
    case(OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA, OPTION__OCTSCHEME__OCT_NLOPT_LBFGS)
#if defined(HAVE_NLOPT)
      !WARNING: not clear if this is needed, probably not.
      call parse_variable(dummy_parser, 'OCTEta', M_ONE, oct%eta)
      call parse_variable(dummy_parser, 'OCTDelta', M_ZERO, oct%delta)
#else
      write(message(1), '(a)') '"OCTScheme = oct_nlopt_bobyqa" or "OCTScheme = oct_nlopt_lbfgs" are'
      write(message(2), '(a)') ' only possible if the nlopt library has been compiled.'
      call messages_fatal(2)
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
    call parse_variable(dummy_parser, 'OCTDoubleCheck', .true., oct%oct_double_check)
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
    call parse_variable(dummy_parser, 'OCTCheckGradient', CNST(0.0), oct%check_gradient)
    call messages_print_var_value(stdout, "OCTCheckGradient", oct%check_gradient)


    !%Variable OCTDirectStep
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 0.25
    !%Description 
    !% If you choose <tt>OCTScheme = oct_direct</tt> or <tt>OCTScheme = oct_nlopt_bobyqa</tt>,
    !% the algorithms necessitate an initial "step" to perform the direct search for the
    !% optimal value. The precise meaning of this "step" differs.
    !%End
    call parse_variable(dummy_parser, 'OCTDirectStep', CNST(0.25), oct%direct_step)
    call messages_print_var_value(stdout, "OCTDirectStep", oct%direct_step)

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
    call parse_variable(dummy_parser, 'OCTNumberCheckPoints', 0, oct%number_checkpoints)
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
    call parse_variable(dummy_parser, 'OCTRandomInitialGuess', .false., oct%random_initial_guess)
    call messages_print_var_value(stdout, "OCTRandomInitialGuess", oct%random_initial_guess)

    call messages_print_stress(stdout)
    POP_SUB(oct_read_inp)
  end subroutine oct_read_inp
  ! ---------------------------------------------------------



  !> Returns .true. if the algorithm to be used is one of the "direct" or "gradient-less"
  !! algorithms -- the ones that do not require backwards propagations. Returns .false. otherwise
  logical pure function oct_algorithm_is_direct(oct)
    type(oct_t), intent(in) :: oct
    oct_algorithm_is_direct = (oct%algorithm == OPTION__OCTSCHEME__OCT_DIRECT) .or. &
                              (oct%algorithm == OPTION__OCTSCHEME__OCT_NLOPT_BOBYQA)
  end function oct_algorithm_is_direct
  ! ---------------------------------------------------------

 
end module opt_control_global_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
