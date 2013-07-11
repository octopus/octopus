!! Copyright (C) 2005-2006 Heiko Appel
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
!! $Id$

#include "global.h"

module sparskit_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use profiling_m

  implicit none

  private

  integer, public, parameter ::  &
    SK_CG      =  1,             &  !< Conjugate Gradient Method
    SK_CGNR    =  2,             &  !< Conjugate Gradient Method (Normal Residual equation)
    SK_BCG     =  3,             &  !< Bi-Conjugate Gradient Method
    SK_DBCG    =  4,             &  !< BCG with partial pivoting
    SK_BCGSTAB =  5,             &  !< BCG stabilized
    SK_TFQMR   =  6,             &  !< Transpose-Free Quasi-Minimum Residual method
    SK_FOM     =  7,             &  !< Full Orthogonalization Method
    SK_GMRES   =  8,             &  !< Generalized Minimum Residual method
    SK_FGMRES  =  9,             &  !< Flexible version of Generalized Minimum Residual method
    SK_DQGMRES = 10,             &  !< Direct versions of Quasi Generalized Minimum Residual method
    SK_MINVAL  = SK_CG,          &
    SK_MAXVAL  = SK_DQGMRES

#ifdef HAVE_SPARSKIT

  public ::                      &
    sparskit_solver_t,           &
    dsparskit_solver_init,       &
    dsparskit_solver_run,        &
    dsparskit_solver_end,        &
    zsparskit_solver_init,       &
    zsparskit_solver_run,        &
    zsparskit_solver_end

  type sparskit_solver_t
    integer :: size                 !< size of the linear system
    integer :: solver_type          !< which solver to use
    integer :: krylov_size          !< size of the Krylov subspace (used for some solvers)
    integer :: preconditioning      !< what kind of preconditioning to use
    integer :: maxiter              !< maximum number of iterations
    integer :: used_iter            !< number of performed iterations
    integer :: iter_out             !< determines how often status info of the solver is printed
    FLOAT   :: residual_norm        !< used store current error norm
    FLOAT   :: rel_tolerance        !< relative tolerance
    FLOAT   :: abs_tolerance        !< absolute tolerance

    FLOAT, allocatable :: sk_work(:), sk_b(:), sk_y(:)

    integer :: ipar(16)             !< integer parameter array for the reverse communication protocol
    FLOAT   :: fpar(16)             !< floating-point parameter array for the reverse communication protocol
    logical :: verbose              !< if .true. then the solver will write more details
  end type sparskit_solver_t


contains


#include "undef.F90"
#include "real.F90"
#include "sparskit_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "sparskit_inc.F90"

#endif /* HAVE_SPARSKIT */

! distdot function for dot products is defined in mesh_function_m

end module sparskit_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
