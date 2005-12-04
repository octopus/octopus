!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module sparskit
  use global
  use messages
  use datasets_mod
  use lib_oct_parser

  implicit none

  private
  public ::                      &
    sparskit_solver_type,        &
    dsparskit_solver_init,       &
    dsparskit_solver_run,        &
    dsparskit_solver_end,        &
    zsparskit_solver_init,       &
    zsparskit_solver_run,        &
    zsparskit_solver_end

  ! the drivers (we make them public for direct calls)
  public ::                      &
    dsk_driver_cg,               &
    dsk_driver_cgnr,             &
    dsk_driver_bcg,              &
    dsk_driver_dbcg,             &
    dsk_driver_bcgstab,          &
    dsk_driver_tfqmr,            &
    dsk_driver_fom,              &
    dsk_driver_gmres,            &
    dsk_driver_fgmres,           &
    dsk_driver_dqgmres,          &
    zsk_driver_cg,               &
    zsk_driver_cgnr,             &
    zsk_driver_bcg,              &
    zsk_driver_dbcg,             &
    zsk_driver_bcgstab,          &
    zsk_driver_tfqmr,            &
    zsk_driver_fom,              &
    zsk_driver_gmres,            &
    zsk_driver_fgmres,           &
    zsk_driver_dqgmres

  integer, private, parameter :: &
    SK_CG      =  1,             &  ! Conjugate Gradient Method
    SK_CGNR    =  2,             &  ! Conjugate Gradient Method (Normal Residual equation)
    SK_BCG     =  3,             &  ! Bi-Conjugate Gradient Method
    SK_DBCG    =  4,             &  ! BCG with partial pivoting
    SK_BCGSTAB =  5,             &  ! BCG stabilized
    SK_TFQMR   =  6,             &  ! Transpose-Free Quasi-Minimum Residual method
    SK_FOM     =  7,             &  ! Full Orthogonalization Method
    SK_GMRES   =  8,             &  ! Generalized Minimum Residual method
    SK_FGMRES  =  9,             &  ! Flexible version of Generalized Minimum Residual method
    SK_DQGMRES = 10,             &  ! Direct versions of Quasi Generalized Minimum Residual method
    SK_MINVAL  = SK_CG,          &
    SK_MAXVAL  = SK_DQGMRES

  FLOAT, allocatable   :: sk_work(:), sk_b(:), sk_y(:)

  type sparskit_solver_type
    integer :: size                 ! size of the linear system
    integer :: solver_type          ! which solver to use
    integer :: krylov_size          ! size of the Krylov subspace (used for some solvers)
    integer :: preconditioning      ! what kind of preconditioning to use
    integer :: maxiter              ! maximum number of iterations
    integer :: used_iter            ! number of performed iterations
    FLOAT   :: rel_tolerance        ! relative tolerance
    FLOAT   :: abs_tolerance        ! absolute tolerance

    integer :: ipar(16)             ! integer parameter array for the reverse communication protocol
    FLOAT   :: fpar(16)             ! floating-point parameter array for the reverse communication protocol
  end type sparskit_solver_type


contains


#include "undef.F90"
#include "real.F90"
#include "sparskit_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "sparskit_inc.F90"


end module sparskit


! ---------------------------------------------------------
FLOAT function distdot(n,x,ix,y,iy)
  use blas
  !  use lib_basic_alg

  integer, intent(in) :: n, ix, iy
  FLOAT, intent(in)   :: x, y

  distdot = ddot(n,x,ix,y,iy)
  !  distdot = lalg_dot(n, x, y)

end function distdot
