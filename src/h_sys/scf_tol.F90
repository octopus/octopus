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
  use loct_parser_m
  use messages_m

  implicit none
  
  private

  integer, public, parameter :: SCF_FIXED = 0
  integer, public, parameter :: SCF_ADAPTIVE = 1

  public :: &
       scf_tol_t, &
       scf_tol_init, & 
       scf_tol_end, & 
       scf_tol_start, & 
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
  end type scf_tol_t

contains

  subroutine scf_tol_init(this, prefix, def_maximumiter)
    type(scf_tol_t),    intent(out) :: this
    character(len=*),   intent(in)  :: prefix
    integer, optional,  intent(in)  :: def_maximumiter

    integer :: def_maximumiter_

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

    if (loct_parse_isdef(check_inp(trim(prefix)//'LRMaximumIter')) /= 0) then 
      call loct_parse_int(check_inp(trim(prefix)//'LRMaximumIter'), def_maximumiter_, this%max_iter)
    else
      call loct_parse_int(check_inp('LRMaximumIter'), def_maximumiter_, this%max_iter)
    end if
    
    !%Variable LRConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The tolerance in the variation of the density, to determine if
    !% the SCF for linear response is converged.
    !%End
    if (loct_parse_isdef(check_inp(trim(prefix)//'LRConvAnsDens')) /= 0 ) then 
      call loct_parse_float(check_inp(trim(prefix)//"LRConvAbsDens"), &
        CNST(1e-5), this%conv_abs_dens)
    else
      call loct_parse_float(check_inp("LRConvAbsDens"), &
        CNST(1e-5), this%conv_abs_dens)
    end if

    !%Variable LRTolScheme
    !%Type integer
    !%Default true
    !%Section Linear Response::SCF in LR calculations
    !%Description
    !% The scheme used to adjust the tolerance of the Solver during
    !% the SCF iteration.
    !%Option fixed 0
    !% The solver tolerance is fixed during all the iteration, this
    !% improves convergency but increases the computational cost
    !%Option adaptive 1 
    !% The tolerance is increased according to the level of
    !% convergency of the SCF.
    !%End
    if (loct_parse_isdef(check_inp(trim(prefix)//'LRTolScheme')) /= 0 ) then
      call loct_parse_int(check_inp(trim(prefix)//'LRTolScheme'), SCF_ADAPTIVE, this%scheme)
    else
      call loct_parse_int(check_inp('LRTolScheme'), SCF_ADAPTIVE, this%scheme)
    end if

    if( this%scheme /= SCF_ADAPTIVE .and. this%scheme /= SCF_FIXED ) then 
      call input_error('LRTolScheme')
    end if

    if( this%scheme == SCF_ADAPTIVE) then 
      
      !%Variable LRAdaptiveTolFactor
      !%Type float
      !%Default 0.5
      !%Section Linear Response::SCF in LR calculations
      !%Description
      !% This factor controls how much the tolerance is increased at
      !% first iterations. Larger values means larger tolerance.
      !%End
      if (loct_parse_isdef(check_inp(trim(prefix)//'LRAdaptiveTolFactor')) /=0 ) then 
        call loct_parse_float(check_inp(trim(prefix)//'LRAdaptiveTolFactor'), &
          CNST(0.5), this%dynamic_tol_factor)
      else
        call loct_parse_float(check_inp('LRAdaptiveTolFactor'), &
          CNST(0.5), this%dynamic_tol_factor)
      end if

    end if

    call pop_sub()

  end subroutine scf_tol_init
    
  FLOAT function scf_tol_start(this, initial_tol, final_tol) result (r)
    type(scf_tol_t), intent(out) :: this
    FLOAT, intent(in) :: initial_tol
    FLOAT, intent(in) :: final_tol

    this%initial_tol = initial_tol
    this%final_tol = final_tol
    
    select case(this%scheme)
      
    case(SCF_FIXED)
      r = this%final_tol
    case(SCF_ADAPTIVE)
      this%current_tol = this%initial_tol
      r =  this%initial_tol
    end select
    
  end function scf_tol_start

  FLOAT function scf_tol_step(this, iter, scf_res) result(r)
    type(scf_tol_t),     intent(inout) :: this
    integer,        intent(in)    :: iter
    FLOAT,          intent(in)    :: scf_res

    select case(this%scheme)

    case(SCF_FIXED)
      r = this%final_tol

    case(SCF_ADAPTIVE)

      r = this%dynamic_tol_factor * (this%final_tol/this%conv_abs_dens)*scf_res 
      r = max(r, this%final_tol)
      r = min(r, this%current_tol)
      
      this%current_tol = r

    end select
      
  end function scf_tol_step


  subroutine scf_tol_stop(this)
    type(scf_tol_t),     intent(inout) :: this
    this%current_tol = M_ZERO
  end subroutine scf_tol_stop

  subroutine scf_tol_end(this)
    type(scf_tol_t),     intent(inout) :: this
    this%current_tol = M_ZERO
  end subroutine scf_tol_end

end module scf_tol_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
