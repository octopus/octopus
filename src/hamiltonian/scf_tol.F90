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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: linear_response.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module scf_tol_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use varinfo_m

  implicit none
  
  private

  integer, public, parameter :: &
       SCF_TOL_FIXED    = 0, &
       SCF_TOL_ADAPTIVE = 1, &
       SCF_TOL_LINEAR   = 2, &
       SCF_TOL_EXP      = 3

  public :: &
       scf_tol_t, &
       scf_tol_init, & 
       scf_tol_end, & 
       scf_tol_stop, & 
       scf_tol_step
       
  type scf_tol_t
     integer           :: max_iter
     integer           :: scheme
     FLOAT             :: conv_abs_dens
     FLOAT             :: dynamic_tol_factor
     FLOAT             :: current_tol
     FLOAT             :: initial_tol
     FLOAT             :: final_tol
     integer           :: iter_window
  end type scf_tol_t

contains

  !-----------------------------------------------------------------
  subroutine scf_tol_init(this, prefix, def_maximumiter, tol_scheme)
    type(scf_tol_t),    intent(out) :: this
    character(len=*),   intent(in)  :: prefix
    integer, optional,  intent(in)  :: def_maximumiter
    integer, optional,  intent(in)  :: tol_scheme

    integer :: def_maximumiter_
    character(len=256) :: str

    call push_sub('scf_tol.scf_tol_init')
    
    !%Variable LRMaximumIter
    !%Type integer
    !%Default 200
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The maximum number of SCF iterations to calculate response.
    !%End
    def_maximumiter_ = 200
    if(present(def_maximumiter)) def_maximumiter_ = def_maximumiter

    str = 'LRMaximumIter'
    if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
         str = trim(prefix)//trim(str)
    call parse_integer(datasets_check(str), def_maximumiter_, this%max_iter)
    
    !%Variable LRConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The tolerance in the variation of the density, to determine if
    !% the SCF for linear response is converged.
    !%End
    str = 'LRConvAnsDens'
    if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
         str = trim(prefix)//trim(str)
    call parse_float(datasets_check(str), CNST(1e-5), this%conv_abs_dens)

    !%Variable LRTolScheme
    !%Type integer
    !%Default adaptive
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The scheme used to adjust the tolerance of the solver during
    !% the SCF iteration. For kdotp and magnetic em_resp modes, or
    !% whenever HamiltonianVariation = V_ext_only, the
    !% scheme is set to fixed, and this variable is ignored.
    !%Option tol_fixed 0
    !% The solver tolerance is fixed for all the iterations; this
    !% improves convergence but increases the computational cost
    !%Option tol_adaptive 1 
    !% The tolerance is increased according to the level of
    !% convergence of the SCF.
    !%Option tol_linear 2
    !% The tolerance decreases linearly for the first LRTolIterWindow iterations
    !%Option tol_exp 3
    !% The tolerance decreases exponentially for the first LRTolIterWindow iterations
    !%End
    if(present(tol_scheme)) then
      this%scheme = tol_scheme
    else
      str = 'LRTolScheme'
      if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
           str = trim(prefix)//trim(str)
      call parse_integer(datasets_check(str), SCF_TOL_ADAPTIVE, this%scheme)
    end if
    if(.not.varinfo_valid_option('LRTolScheme', this%scheme)) &
         call input_error('LRTolScheme')

    !%Variable LRTolInitTol
    !%Type float
    !%Default 1e-2
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged,
    !% for the first SCF iteration. Ignored if LRTolScheme = fixed.
    !%End
    str = 'LRTolInitTol'
    if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
         str = trim(prefix)//trim(str)
    call parse_float(datasets_check(str), CNST(1e-2), this%initial_tol)
    this%current_tol = this%initial_tol

    !%Variable LinearSolverTol
    !%Type float
    !%Default 1e-6
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged.
    !%End
    str = 'LRTolFinalTol'
    if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
         str = trim(prefix)//trim(str)
    call parse_float(datasets_check(str), CNST(1e-6), this%final_tol)

    if(this%scheme == SCF_TOL_ADAPTIVE) then 
      !%Variable LRTolAdaptiveFactor
      !%Type float
      !%Default 0.9
      !%Section Linear Response::SCF in LR calculations
      !%Description
      !% This factor controls how much the tolerance is decreased
      !% during the self-consistency process. Smaller values mean that
      !% tolerance is decreased faster. The default is 0.9.
      !%End
      str = 'LRAdaptiveTolFactor'
      if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
           str = trim(prefix)//trim(str)
      call parse_float(datasets_check(str), CNST(0.1), this%dynamic_tol_factor)
    end if

    if(this%scheme==SCF_TOL_LINEAR.or.this%scheme==SCF_TOL_EXP) then
      !%Variable LRTolIterWindow
      !%Type float
      !%Default 10
      !%Section Linear Response::SCF in LR calculations
      !%Description
      !% Number of iterations necessary to reach the final tolerance
      !%End
      str = 'LRTolIterWindow'
      if(parse_isdef(datasets_check(trim(prefix)//trim(str))) /= 0) &
           str = trim(prefix)//trim(str)
      call parse_integer(datasets_check(str), 10, this%iter_window)
    end if

    call pop_sub('scf_tol.scf_tol_init')

  end subroutine scf_tol_init
    

  !-----------------------------------------------------------------
  FLOAT function scf_tol_step(this, iter, scf_res) result(r)
    type(scf_tol_t), intent(inout) :: this
    integer,         intent(in)    :: iter
    FLOAT,           intent(in)    :: scf_res

    FLOAT :: logi, logf

    if(iter == 0) this%current_tol = M_HUGE

    select case(this%scheme)
    case(SCF_TOL_FIXED)
      r = this%final_tol

    case(SCF_TOL_ADAPTIVE)
      if(iter == 0) then
        r = this%initial_tol
      else
        r = this%dynamic_tol_factor * (this%final_tol/this%conv_abs_dens)*scf_res
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

  end function scf_tol_step


  !-----------------------------------------------------------------
  subroutine scf_tol_stop(this)
    type(scf_tol_t),     intent(inout) :: this
    this%current_tol = M_ZERO
  end subroutine scf_tol_stop


  !-----------------------------------------------------------------
  subroutine scf_tol_end(this)
    type(scf_tol_t),     intent(inout) :: this
    this%current_tol = M_ZERO
  end subroutine scf_tol_end

end module scf_tol_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
