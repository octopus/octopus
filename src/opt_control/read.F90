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
!! $Id: opt_control.F90 2862 2007-04-25 23:34:20Z acastro $


  ! ---------------------------------------------------------
  subroutine oct_read_inp(oct)
    type(oct_t), intent(inout) :: oct

    !read the parameters for the optimal control run     
    integer :: kk

    call push_sub('opt_control_read.oct_read_inp')  

    !%Variable OCTFixFluenceTo
    !%Type float
    !%Section Optimal Control
    !%Default -1.0
    !%Description
    !% The algorithm tries to obtain the specified fluence for the laser field. 
    !% This works only in conjunction with the WG05 scheme.
    !%End
    oct%mode_fixed_fluence = .false.
    call loct_parse_float(check_inp('OCTFixFluenceTo'), -M_ONE, oct%targetfluence)
    if (oct%targetfluence.ne.-M_ONE) oct%mode_fixed_fluence = .true.
      
    !%Variable OCTTargetMode
    !%Type integer
    !%Section Optimal Control
    !%Description
    !%Option oct_targetmode_static  0
    !% Static or time-independent targets
    !%Option oct_targetmode_td      1
    !% Time-dependent targets, specify block OCTTdTarget
    !%End
    call loct_parse_int(check_inp('OCTTargetMode'), oct_targetmode_static, oct%targetmode)
    if(.not.varinfo_valid_option('OCTTargetMode', oct%targetmode)) call input_error('OCTTargetMode')
       
    !%Variable OCTScheme
    !%Type integer
    !%Section Optimal Control
    !%Default oct_algorithm_zbr98
    !%Description
    !% In order to find the optimal laser field for a given task, e.g., the excitation from an
    !% initial state to a predefined final state at the final time, optimal control theory can 
    !% be applied to quantum mechanics. The mathematical derivation leads a set of equations 
    !% which require the propagation of the wavefunction and a lagrange multiplier (sometimes
    !% comparable to a wavefunction). Several schemes have been sought to solve these control 
    !% equations which boils down to forward and backward propagations. However, the order in
    !% which these equations are solved makes a huge difference. Some schemes can be proven 
    !% to increase the value of the target functional (merit function) in each step. (In 
    !% practice this can be violated if the accuracy of the numerical time propagation is
    !% small. Most likely in 3D.)
    !%Option oct_algorithm_zbr98 1 
    !% Backward-Forward-Backward scheme described in JCP 108, 1953 (1998).
    !% Only possible if targetoperator is a projection operator
    !% Provides the fastest and most stable convergence.
    !% Monotonic convergence.
    !%Option oct_algorithm_zr98  2
    !% Forward-Backward-Forward scheme described in JCP 109,385 (1998).
    !% Works for projection and local target operators
    !% Convergence is stable but slower than ZBR98. 
    !% Note that local operators show an extremely slow convergence.
    !% Monotonic convergence.
    !%Option oct_algorithm_wg05  3
    !% Forward-Backward scheme described in J. Opt. B. 7 300 (2005).
    !% Works for all kind target operators and 
    !% can be used with all kind of filters and allows a fixed fluence.
    !% The price is a rather instable convergence. 
    !% If the restrictions set by the filter and fluence are reasonable, a good overlap can be 
    !% expected with 20 iterations.
    !% No monotonic convergence.
    !%Option oct_algorithm_krotov 4
    !% Yet to be implemented and tested, some people swear on it.
    !% Described in Zhao & Rice: Optical Control of Molecules. (Appendix A)
    !%Option oct_algorithm_mt03 5
    !% Yet to be implemented and tested. Basically an improved and generalized scheme. 
    !% Comparable to ZBR98/ZR98. (see JCP 118,8191 (2003))
    !%End
    call loct_parse_int(check_inp('OCTScheme'), oct_algorithm_zr98, oct%algorithm_type)
    if(.not.varinfo_valid_option('OCTScheme', oct%algorithm_type)) call input_error('OCTScheme')

    !%Variable OCTDoubleCheck
    !%Type logical
    !%Section Optimal Control
    !%Default true
    !%Description 
    !% Run a normal propagation after the optimization using the optimized field.
    !%End
    call loct_parse_logical(check_inp('OCTDoubleCheck'), .TRUE., oct%oct_double_check)


    !%Variable OCTDumpIntermediate
    !%Type logical
    !%Section Optimal Control
    !%Default true
    !%Description 
    !% Writes to disk some data during the OCT algorithm at intermediate steps.
    !% This is rather technical and it should be considered only for debugging
    !% purposes. Nevertheless, since the whole OCT infrastructure is at a very
    !% preliminary developing stage, it is set to true by default
    !%End
    call loct_parse_logical(check_inp('OCTDumpIntermediate'), .true., oct%dump_intermediate)

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Optimal Control
    !%Default 2
    !%Description
    !% The string OCTTargetOperator describes the initial state of the quantum system
    !% Possible arguments are:
    !%Option oct_tg_groundstate 1 
    !% Targetoperator is a projection operator on the ground state
    !%Option oct_tg_excited 2
    !% Currently not in use.
    !%Option oct_tg_gstransformation 3
    !% Targetoperator is a projection operator on a transformation of the ground state 
    !% orbitals defined by the block OCTTargetTransformStates
    !%Option oct_tg_userdefined 4
    !% Targetoperator is a projection operator on a user defined state
    !%Option oct_tg_local 5
    !% Targetoperator is a local operator. Specify the shape and position within the block OCTLocalTarget.
    !%Option oct_tg_td_local 6
    !% Target operator is time-dependent, please specify block OCTTdTarget
    !%End
    call loct_parse_int(check_inp('OCTTargetOperator'),oct_tg_excited, oct%totype)
    if(.not.varinfo_valid_option('OCTTargetOperator', oct%totype)) call input_error('OCTTargetOperator')    

    call pop_sub()
  end subroutine oct_read_inp


  ! ---------------------------------------------------------
  ! Tries to avoid ill defined combinations of run modes
  ! be careful with the order !!
  subroutine check_faulty_runmodes(oct)
    type(oct_t), intent(inout) :: oct

    integer :: jj
    call push_sub('opt_control.check_faulty_runmodes')

    ! FixedFluence and Filter only work with WG05
    if((oct%mode_fixed_fluence).and.(oct%algorithm_type.ne.oct_algorithm_wg05)) then
      write(message(1),'(a)') "Cannot optimize to a given fluence with the chosen algorithm."
      write(message(2),'(a)') "Switching to scheme WG05."         
      call write_info(2)
      oct%algorithm_type = oct_algorithm_wg05
    end if
      
    ! WARNING filters can only be used with WG05, and this is not checked any more.
      
    ! WARNING: tdpenalty and fixed fluence do not work, and this is not checked any more.
      
    ! local targets only in ZR98 and WG05
    if((oct%totype.eq.oct_tg_local) & 
      .AND.(oct%algorithm_type.eq.oct_algorithm_zbr98)) then
      write(message(1),'(a)') "Warning: Local targets work" &
        // " only with ZR98 and WG05."
      write(message(2),'(a)') "Warning: Switching to ZR98."
      call write_info(2)
      oct%algorithm_type = oct_algorithm_zr98
    end if
      
    ! tdtargets only in ZR98 and WG05
    if((oct%targetmode.eq.oct_targetmode_td) & 
      .AND.(oct%algorithm_type.eq.oct_algorithm_zbr98)) then
      write(message(1),'(a)') "Warning: Time-dependent targets work" &
        // " only with ZR98 and WG05."
      write(message(2),'(a)') "Warning: Please change algorithm type."
      call write_fatal(2)
    end if

    call pop_sub()      
  end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
