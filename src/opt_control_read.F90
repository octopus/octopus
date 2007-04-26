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

#include "global.h"


  ! ---------------------------------------------------------
  subroutine oct_read_inp(oct)
    type(oct_t), intent(inout) :: oct

    !read the parameters for the optimal control run     
    integer :: kk

    call push_sub('opt_control_read.oct_read_inp')  

    !%Variable OCTEps
    !%Type float
    !%Section Optimal Control
    !%Default 0.001
    !%Description
    !% Define the convergence threshold.
    !% For the monotonically convergent scheme: If the increase of the 
    !% target functional is less then OCTEps the iteration is stopped.
    !% Example
    !% OCTEps = 0.00001
    !%End
    call loct_parse_float(check_inp('OCTEps'), CNST(1.0e-3), oct%eps)

    !%Variable OCTMaxIter
    !%Type integer
    !%Section Optimal Control
    !%Default 10
    !%Description
    !% OCTMaxIter defines the maximum number of iterations.
    !% Typical values range from 10-100.
    !%End
    call loct_parse_int(check_inp('OCTMaxIter'), 10, oct%ctr_iter_max)
      
    if( oct%ctr_iter_max < 0 .and. oct%eps < M_ZERO ) then
      message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
      call write_fatal(1)
    end if
    if(oct%ctr_iter_max < 0) oct%ctr_iter_max = huge(oct%ctr_iter_max)
   
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
    ALLOCATE(oct%a_penalty(0:oct%ctr_iter_max+1), oct%ctr_iter_max+2)
    call loct_parse_float(check_inp('OCTPenalty'), M_ONE, oct%penalty)
    ! penalty array for fixed fluence run 
    ! the array is only interesting for the development of new algorithms
    oct%a_penalty = oct%penalty

    ! read in laser polarization and degress of freedom
    call oct_read_laserpol(oct%laser_pol, oct%dof)

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
       
    !%Variable OCTFilterMode
    !%Type integer
    !%Section Optimal Control
    !%Default 0
    !%Description     
    !%End
    call loct_parse_int(check_inp('OCTFilterMode'), 0, oct%filtermode)

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
      
    call pop_sub()
  end subroutine oct_read_inp


  ! ---------------------------------------------------------
  subroutine oct_read_laserpol(laserpol, dof)
    ! read the parameters for the laser polarization
    CMPLX,   intent(out) :: laserpol(MAX_DIM)
    integer, intent(out) :: dof
  
    integer   :: i, no_blk, no_c
    C_POINTER :: blk

    call push_sub('opt_control_read.oct_read_laserpol') 

    !%Variable OCTPolarization
    !%Type block
    !%Section Optimal Control
    !%Description
    !% Define how many degress of freedom the laser has and how it is polarized.
    !% The different examples below will explain this.
    !% First of all, the syntax of each line is:
    !%
    !% <tt>%OCTPolarization
    !% <br>&nbsp;&nbsp;dof | pol1 | pol2 | pol3 
    !% <br>%</tt>
    !% The variable defines the degress of freedom which is either 1 or 2.
    !% pol1, pol2, pol3 define the polarization.
    !%
    !% Some examples:
    !%
    !% <tt>%OCTPolarization
    !% <br>&nbsp;&nbsp;1 | 1 | 0 | 0 
    !% <br>%</tt> 
    !% Here we try to optimize the x-polarized laser.
    !%
    !% <tt>%OCTPolarization
    !% <br>&nbsp;&nbsp;1 | 1 | 1 | 0 
    !% <br>%</tt> 
    !% The polarization is linear and lies in the x-y plane, only one laser field is optimized.
    !%
    !% <tt>%OCTPolarization
    !% <br>&nbsp;&nbsp;2 | 1 | 1 | 0 
    !% <br>%</tt> 
    !% The polarization lies in the x-y plane, but this time two components of the laser
    !% field are optimized. This may lead to linear, circular, or ellipitically polarized 
    !% fields, dependening on the problem.
    !%
    !% <tt>%OCTPolarization
    !% <br>&nbsp;&nbsp;1 | 1 | i | 0 
    !% <br>%</tt> 
    !% If we know that the answer must be a circular polarized field we can also fix it 
    !% and optimize only one component.  
    !%End
    if(loct_parse_block(check_inp('OCTPolarization'),blk)==0) then
      no_blk = loct_parse_block_n(blk)
      ! TODO: for each orbital            
      do i=1, no_blk
        no_c = loct_parse_block_cols(blk, i-1)
        call loct_parse_block_int(blk,i-1, 0, dof)
        if((dof.gt.3).OR.(dof.lt.1) ) then
          message(1) = "OCTPolarization: Choose degrees of freedom between 1 and 3."
          call write_fatal(1)
          
        end if
        call loct_parse_block_cmplx(blk,i-1, 1, laserpol(1))
        call loct_parse_block_cmplx(blk,i-1, 2, laserpol(2))
        call loct_parse_block_cmplx(blk,i-1, 3, laserpol(3))
        
      end do
      call loct_parse_block_end(blk)
    else
      message(1) = 'Input: No Polarization defined and degrees of freedom defined'
      message(2) = 'Input: Using default: x-polarized.'
      laserpol(1) = m_z1
      laserpol(2) = m_z0
      laserpol(3) = m_z0
      dof  = 1
      call write_info(2)
    end if

    call pop_sub()
  end subroutine oct_read_laserpol


  ! ---------------------------------------------------------
  subroutine read_state(st, m, filename)
    type(states_t),   intent(inout) :: st
    type(mesh_t),     intent(in)  :: m
    character(len=*), intent(in)  :: filename

    integer :: ierr

    call push_sub('opt_control.read_state')

    call zinput_function('tmp/restart_gs/'//trim(filename), m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)
    ! if we do not succeed try obf
    if(ierr>0) call zinput_function ('tmp/restart_gs/'//trim(filename)//'.obf', &
      m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)
    ! if we do not succeed try NetCDF
    if(ierr>0) call zinput_function('tmp/restart_gs/'//trim(filename)//'.ncdf', &
      m, st%zpsi(:, 1, 1, 1), ierr, is_tmp=.TRUE.)

    if(ierr > 0) then
       message(1) = "Unsuccesfull read of states in 'tmp/restart_gs/" // trim(filename) // "'"
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine read_state

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
      
    ! Filters work only with WG05
    if((oct%filtermode.gt.0).and.(oct%algorithm_type.ne.oct_algorithm_wg05)) then
      write(message(1),'(a)') "Warning: Cannot use filters with the chosen algorithm."
      write(message(2),'(a)') "Warning: Switching to scheme WG05."
      call write_info(2)
      oct%algorithm_type = oct_algorithm_wg05
    end if
      
    ! tdpenalty and fixed fluence do not work !
    if((oct%mode_tdpenalty).AND.(oct%mode_fixed_fluence)) then
      write(message(1),'(a)') "Warning: Cannot use fixed fluence and" &
        //" td penalty."
      write(message(2),'(a)') "Warning: Disabling td penalty."
      call write_info(2)
      oct%tdpenalty = oct%penalty
    end if
      
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
!!$      
!!$      ! td target with states works only for 1 electron
!!$      ! TODO: HOW TO RETRIEVE NUMBER OF ELECTRONS
!!$      !       UNCOMMENT LAST LINES
!!$      if(associated(td_tg)) then
!!$        do jj=1, size(td_tg)
!!$          if(td_tg(jj)%type.eq.oct_tgtype_state) &
!!$            td_tg_state = .TRUE.
!!$        end do
!!$      end if
!!$      ! occ in states_t
!!$      !write(6,*) 'DEBUG: occ: ', shape(psi_i%occ), size(psi_i%occ)
!!$      if( (td_tg_state) .AND. (SUM(psi_i%occ(1,:)).gt.1) ) then
!!$      !if(td_tg_state) then
!!$      !     if(psi_i%occ.gt.1) then
!!$        write(message(1),'(a)') "Warning: Time-dependent state targets" &
!!$          // " work only with one electron."
!!$        write(message(2),'(a)') "Warning: Please change input file."
!!$        call write_fatal(2)
!!$      end if
!!$
      call pop_sub()      
    end subroutine check_faulty_runmodes


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
