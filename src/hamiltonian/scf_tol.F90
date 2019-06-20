!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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

module scf_tol_oct_m
  use global_oct_m
  use messages_oct_m
  use parser_oct_m
  use varinfo_oct_m

  implicit none
  
  private

  integer, public, parameter :: &
       SCF_TOL_FIXED    = 0,    &
       SCF_TOL_ADAPTIVE = 1,    &
       SCF_TOL_LINEAR   = 2,    &
       SCF_TOL_EXP      = 3

  public ::                       &
       scf_tol_t,                 &
       scf_tol_init,              & 
       scf_tol_end,               & 
       scf_tol_stop,              & 
       scf_tol_step,              &
       scf_tol_final,             &
       scf_tol_obsolete_variables

  type scf_tol_t
    private
     integer, public   :: max_iter
     integer           :: scheme
     FLOAT,   public   :: conv_rel_dens
     FLOAT,   public   :: conv_abs_dens
     FLOAT             :: conv_threshold_use
     FLOAT             :: dynamic_tol_factor
     FLOAT             :: current_tol
     FLOAT             :: initial_tol
     FLOAT,   public   :: final_tol
     integer           :: iter_window
  end type scf_tol_t

contains

  !-----------------------------------------------------------------
  subroutine scf_tol_init(this, parser, qtot, def_maximumiter, tol_scheme)
    type(scf_tol_t),    intent(out) :: this
    type(parser_t),     intent(in)  :: parser
    FLOAT,              intent(in)  :: qtot
    integer, optional,  intent(in)  :: def_maximumiter
    integer, optional,  intent(in)  :: tol_scheme

    integer :: def_maximumiter_

    PUSH_SUB(scf_tol_init)
    
    !%Variable LRMaximumIter
    !%Type integer
    !%Default 200
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The maximum number of SCF iterations to calculate response.
    !%End
    def_maximumiter_ = 200
    if(present(def_maximumiter)) def_maximumiter_ = def_maximumiter

    call parse_variable(dummy_parser, 'LRMaximumIter', def_maximumiter_, this%max_iter)
    
    !%Variable LRConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The tolerance in the absolute variation of the density response, to determine if
    !% the SCF for linear response is converged.
    !% <math>\varepsilon = \int {\rm d}^3r \left| \rho^{out}(\bf r) -\rho^{inp}(\bf r) \right|</math>.
    !% A zero value means do not use this criterion.
    !%End

    call parse_variable(dummy_parser, 'LRConvAbsDens', CNST(1e-5), this%conv_abs_dens)

    !%Variable LRConvRelDens
    !%Type float
    !%Default 0.0
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The tolerance in the relative variation of the density response, to determine if
    !% the SCF for linear response is converged.
    !% <math>\varepsilon = \frac{1}{N_{\rm elec}}</math> <tt>LRConvAbsDens</tt>.
    !% A zero value means do not use this criterion.
    !%End
    
    call parse_variable(dummy_parser, 'LRConvRelDens', M_ZERO, this%conv_rel_dens)

    ! value to use for adaptive tol scheme
    if(this%conv_abs_dens <= M_ZERO) then
      this%conv_threshold_use = this%conv_rel_dens * qtot
    else if(this%conv_rel_dens <= M_ZERO) then
      this%conv_threshold_use = this%conv_abs_dens
    else
      this%conv_threshold_use = min(this%conv_abs_dens, this%conv_rel_dens * qtot)
    end if

    if(this%conv_abs_dens <= M_ZERO .and. this%conv_rel_dens <= M_ZERO) then
      message(1) = "Input: Not all convergence criteria can be <= 0"
      message(2) = "Please set one of the following to a nonzero value:"
      message(3) = "LRConvAbsDens | LRConvRelDens"
      call messages_fatal(3)
    end if

    !%Variable LRTolScheme
    !%Type integer
    !%Default tol_adaptive
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The scheme used to adjust the tolerance of the solver during
    !% the SCF iteration. For <tt>kdotp</tt> and magnetic <tt>em_resp</tt> modes, or
    !% whenever <tt>HamiltonianVariation = V_ext_only</tt>, the
    !% scheme is set to <tt>tol_fixed</tt>, and this variable is ignored.
    !%Option tol_fixed 0
    !% The solver tolerance is fixed for all the iterations; this
    !% improves convergence but increases the computational cost
    !%Option tol_adaptive 1 
    !% The tolerance is increased according to the level of
    !% convergence of the SCF.
    !%Option tol_linear 2
    !% The tolerance decreases linearly for the first <tt>LRTolIterWindow</tt> iterations.
    !%Option tol_exp 3
    !% The tolerance decreases exponentially for the first <tt>LRTolIterWindow</tt> iterations.
    !%End
    if(present(tol_scheme)) then
      this%scheme = tol_scheme
    else
      call parse_variable(dummy_parser, 'LRTolScheme', SCF_TOL_ADAPTIVE, this%scheme)
    end if
    if(.not.varinfo_valid_option('LRTolScheme', this%scheme)) &
      call messages_input_error('LRTolScheme')

    !%Variable LRTolInitTol
    !%Type float
    !%Default 1e-2
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged,
    !% for the first SCF iteration. Ignored if <tt>LRTolScheme = fixed</tt>.
    !%End

    call parse_variable(dummy_parser, 'LRTolInitTol', CNST(1e-2), this%initial_tol)
    this%current_tol = this%initial_tol

    !%Variable LRTolFinalTol
    !%Type float
    !%Default 1e-6
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged.
    !%End

    call parse_variable(dummy_parser, 'LRTolFinalTol', CNST(1e-6), this%final_tol)

    if(this%scheme == SCF_TOL_ADAPTIVE) then 
      !%Variable LRTolAdaptiveFactor
      !%Type float
      !%Default 0.1
      !%Section Linear Response::SCF in LR calculations
      !%Description
      !% This factor controls how much the tolerance is decreased
      !% during the self-consistency process. Smaller values mean that
      !% tolerance is decreased faster.
      !%End

      call parse_variable(dummy_parser, 'LRTolAdaptiveFactor', CNST(0.1), this%dynamic_tol_factor)
    end if

    if(this%scheme==SCF_TOL_LINEAR.or.this%scheme==SCF_TOL_EXP) then
      !%Variable LRTolIterWindow
      !%Type float
      !%Default 10
      !%Section Linear Response::SCF in LR calculations
      !%Description
      !% Number of iterations necessary to reach the final tolerance.
      !%End

      call parse_variable(dummy_parser, 'LRTolIterWindow', 10, this%iter_window)
    end if

    POP_SUB(scf_tol_init)

  end subroutine scf_tol_init
    

  !-----------------------------------------------------------------
  FLOAT function scf_tol_step(this, iter, scf_res) result(r)
    type(scf_tol_t), intent(inout) :: this
    integer,         intent(in)    :: iter
    FLOAT,           intent(in)    :: scf_res

    FLOAT :: logi, logf

    PUSH_SUB(scf_tol_step)

    if(iter == 0) this%current_tol = M_HUGE

    select case(this%scheme)
    case(SCF_TOL_FIXED)
      r = this%final_tol

    case(SCF_TOL_ADAPTIVE)
      if(iter == 0) then
        r = this%initial_tol
      else
        r = this%dynamic_tol_factor * (this%final_tol/this%conv_threshold_use)*scf_res
      end if

    case(SCF_TOL_LINEAR)
      r = this%initial_tol + (this%final_tol - this%initial_tol) * &
        real(iter, REAL_PRECISION) / real(this%iter_window, REAL_PRECISION)

    case(SCF_TOL_EXP)
      logi = log(this%initial_tol)
      logf = log(this%final_tol)
      r = logi + (logf - logi) * &
        real(iter, REAL_PRECISION) / real(this%iter_window, REAL_PRECISION)
      r = exp(r)
    end select
      
    ! tolerance can never be larger than final tolerance
    r = max(r, this%final_tol)
    ! tolerance always has to decrease
    r = min(r, this%current_tol)

    this%current_tol = r

    POP_SUB(scf_tol_step)
  end function scf_tol_step


  !-----------------------------------------------------------------
  subroutine scf_tol_stop(this)
    type(scf_tol_t),     intent(inout) :: this

    PUSH_SUB(scf_tol_stop)
    this%current_tol = M_ZERO

    POP_SUB(scf_tol_stop)
  end subroutine scf_tol_stop


  !-----------------------------------------------------------------
  subroutine scf_tol_end(this)
    type(scf_tol_t),     intent(inout) :: this

    PUSH_SUB(scf_tol_end)
    this%current_tol = M_ZERO

    POP_SUB(scf_tol_end)
  end subroutine scf_tol_end


  !-----------------------------------------------------------------
  FLOAT function scf_tol_final(this)
    type(scf_tol_t), intent(in) :: this

    scf_tol_final = this%final_tol

  end function scf_tol_final

  !-----------------------------------------------------------------

  subroutine scf_tol_obsolete_variables(parser, old_prefix, new_prefix)
    type(parser_t),      intent(in)    :: parser
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix

    PUSH_SUB(scf_tol_obsolete_variables)

    call messages_obsolete_variable(parser, trim(old_prefix)//'LRMaximumIter', trim(new_prefix)//'LRMaximumIter')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRConvAbsDens', trim(new_prefix)//'LRConvAbsDens')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRTolScheme', trim(new_prefix)//'LRTolScheme')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRTolInitTol', trim(new_prefix)//'LRTolInitTol')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRTolFinalTol', trim(new_prefix)//'LRTolFinalTol')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRTolAdaptiveFactor', trim(new_prefix)//'LRTolAdaptiveFactor')
    call messages_obsolete_variable(parser, trim(old_prefix)//'LRTolIterWindow', trim(new_prefix)//'LRTolIterWindow')

    POP_SUB(scf_tol_obsolete_variables)
  end subroutine scf_tol_obsolete_variables

end module scf_tol_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
