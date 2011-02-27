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
!! $Id$

#include "global.h"

! This module is intended to contain "only mathematical" functions
! and procedures.

module solvers_m
  use blas_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::                     &
    dconjugate_gradients,       &
    zconjugate_gradients,       &
    zqmr_sym,                   &
    dqmr_sym,                   &
    zqmr,                       &
    dqmr

  ! ---------------------------------------------------------
  ! QMR (quasi-minimal residual) algorithm for complex symmetric matrices
  ! algorithm taken from:
  ! Parallel implementation of efficient preconditioned linear solver for
  ! grid-based applications in chemical physics. II: QMR linear solver
  ! Appendix A. Simplified QMR algorithm
  ! W Chen and B Poirier, J Comput Phys 219, 198-209 (2006)
  interface zqmr_sym
    module procedure zqmr_sym_spec_dotu, zqmr_sym_gen_dotu
  end interface

  interface dqmr_sym
    module procedure dqmr_sym_spec_dotu, dqmr_sym_gen_dotu
  end interface
  integer, pointer :: np_p

  ! ---------------------------------------------------------
  ! QMR (quasi-minimal residual) algorithm for complex matrices
  ! algorithm taken from: An Implementation of the QMR Method based on
  ! Coupled Two-Term Recurrences by R. W. Freund and N. M. Nachtigal (page 25)
  ! http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950017192_1995117192.pdf
  interface zqmr
    module procedure zqmr_spec_dotu, zqmr_gen_dotu
  end interface

  interface dqmr
    module procedure dqmr_spec_dotu, dqmr_gen_dotu
  end interface

  interface dconjugate_gradients
    module procedure dsym_conjugate_gradients, dbi_conjugate_gradients
  end interface

  interface zconjugate_gradients
    module procedure zsym_conjugate_gradients, zbi_conjugate_gradients
  end interface

contains

  ! ---------------------------------------------------------
  ! the complex dot product without conjugated vector
  CMPLX function zdotu_qmr(x, y)
    CMPLX, intent(in) :: x(:)
    CMPLX, intent(in) :: y(:)

    PUSH_SUB(zdotu_qmr)

    ASSERT(ubound(x, dim = 1) >= np_p)
    ASSERT(ubound(y, dim = 1) >= np_p)

    zdotu_qmr = blas_dotu(np_p, x(1), 1, y(1), 1)

    POP_SUB(zdotu_qmr)
  end function zdotu_qmr

  FLOAT function ddotu_qmr(x, y)
    FLOAT, intent(in) :: x(:)
    FLOAT, intent(in) :: y(:)

    PUSH_SUB(ddotu_qmr)

    ASSERT(ubound(x, dim = 1) >= np_p)
    ASSERT(ubound(y, dim = 1) >= np_p)

    ddotu_qmr = blas_dot(np_p, x(1), 1, y(1), 1)

    POP_SUB(ddotu_qmr)
  end function ddotu_qmr


#include "undef.F90"
#include "complex.F90"
#include "solvers_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "solvers_inc.F90"

end module solvers_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
