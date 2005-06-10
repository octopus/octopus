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
  use syslabels

  implicit none

  private
  public :: sparskit_solver_type, &
       sparskit_solver_init, &
       sparskit_solver_run,  & 
       sparskit_solver_end

  ! the drivers
  public :: sk_driver_cg, sk_driver_cgnr, sk_driver_bcg, sk_driver_dbcg, &
       sk_driver_bcgstab, sk_driver_tfqmr, sk_driver_fom, sk_driver_gmres, &
       sk_driver_fgmres, sk_driver_dqgmres


  integer, public, parameter :: &
       SK_CG      =  1,   &  ! Conjugate Gradient Method
       SK_CGNR    =  2,   &  ! Conjugate Gradient Method (Normal Residual equation)
       SK_BCG     =  3,   &  ! Bi-Conjugate Gradient Method
       SK_DBCG    =  4,   &  ! BCG with partial pivoting 
       SK_BCGSTAB =  5,   &  ! BCG stabilized
       SK_TFQMR   =  6,   &  ! Transpose-Free Quasi-Minimum Residual method
       SK_FOM     =  7,   &  ! Full Orthogonalization Method
       SK_GMRES   =  8,   &  ! Generalized Minimum Residual method
       SK_FGMRES  =  9,   &  ! Flexible version of Generalized Minimum Residual method
       SK_DQGMRES = 10       ! Direct versions of Quasi Generalize Minimum Residual method

  FLOAT, allocatable   :: work(:)

  type sparskit_solver_type
     integer :: size            ! size of the linear system
     integer :: solver_type     ! which solver to use
     integer :: krylov_size     ! size of the Krylov subspace (used for some solvers)
     integer :: preconditioning ! what kind of preconditioning to use
     integer :: maxiter         ! maximum number of iterations
     integer :: used_iter       ! number of performed iterations 
     FLOAT   :: rel_tolerance   ! relative tolerance
     FLOAT   :: abs_tolerance   ! absolute tolerance

     integer :: ipar(16)        ! integer parameter array for the reverse communication protocol
     FLOAT   :: fpar(16)        ! floating-point parameter array for the reverse communication protocol
  end type sparskit_solver_type


contains

  ! ---------------------------------------------------------
  subroutine sparskit_solver_init(n, sk)
    type(sparskit_solver_type), intent(out) :: sk
    integer, intent(in)  :: n

    integer :: workspace_size, m 
    call push_sub('sparskit_solver_init')


    !%Variable SparskitSolver
    !%Type integer
    !%Section 1 Generalities
    !%Description
    !% Specifies what kind of linear solver will be used
    !%Option sk_cg 1
    !% Conjugate Gradient Method
    !%Option sk_cgnr 2
    !% Conjugate Gradient Method (Normal Residual equation)
    !%Option sk_bcg 3
    !% Bi-Conjugate Gradient Method
    !%Option sk_dbcg 4
    !% BCG with partial pivoting 
    !%Option sk_bcgstab 5
    !% BCG stabilized
    !%Option sk_tfqmr 6
    !% Transpose-Free Quasi-Minimum Residual method
    !%Option sk_fom 7
    !% Full Orthogonalization Method
    !%Option sk_gmres 8
    !% Generalized Minimum Residual method
    !%Option sk_fgmres 9
    !% Flexible version of Generalized Minimum Residual method
    !%Option sk_dqgmres 10
    !% Direct versions of Quasi Generalize Minimum Residual method
    !%End
    call loct_parse_int(check_inp('SparskitSolver'),          SK_CG, sk%solver_type)
    call loct_parse_int(check_inp('SparskitKrylovSubspaceSize'), 15, sk%krylov_size)
    call loct_parse_int(check_inp('SparskitPreconditioning'),     0, sk%preconditioning)
    call loct_parse_int(check_inp('SparskitMaxIter'),          1000, sk%maxiter)
    if (sk%preconditioning.ne.0) then
       message(1) = 'Error: Preconditioning not implemented yet ...'
       call write_fatal(1)
    endif
    call loct_parse_float(check_inp('SparskitRelTolerance'),    CNST(1e-4), sk%rel_tolerance)
    call loct_parse_float(check_inp('SparskitAbsTolerance'),    CNST(1e-8), sk%abs_tolerance)

    ! size of the problem
    sk%size = n
    ! Krylov subspace size
    m = sk%krylov_size
    if (mod(m, 2).ne.0) m = m + 1

    select case(sk%solver_type)
    case(SK_CG)
       message(1) = 'Info: Sparskit solver type: Conjugate Gradient Method'
       workspace_size = 5*n
    case(SK_CGNR)
       message(1) = 'Info: Sparskit solver type: Conjugate Gradient Method (Normal Residual equation)'
       workspace_size = 5*n
    case(SK_BCG)
       message(1) = 'Info: Sparskit solver type: Bi-Conjugate Gradient Method'
       workspace_size = 7*n
    case(SK_DBCG)
       message(1) = 'Info: Sparskit solver type: BCG with partial pivoting'
       workspace_size = 11*n
    case(SK_BCGSTAB)
       message(1) = 'Info: Sparskit solver type: BCG stabilized'
       workspace_size = 8*n
    case(SK_TFQMR)
       message(1) = 'Info: Sparskit solver type: Transpose-Free Quasi-Minimum Residual method'
       workspace_size = 11*n
    case(SK_FOM)
       message(1) = 'Info: Sparskit solver type: Full Orthogonalization Method'
       workspace_size = (n+3)*(m+2) + (m+1)*m/2
    case(SK_GMRES)
       message(1) = 'Info: Sparskit solver type: Generalized Minimum Residual method'
       workspace_size = (n+3)*(m+2) + (m+1)*m/2       
    case(SK_FGMRES)
       message(1) = 'Info: Sparskit solver type: Flexible version of Generalized Minimum Residual method'
       workspace_size =  2*n*(m+1) + (m+1)*m/2 + 3*m + 2
    case(SK_DQGMRES)
       message(1) = 'Info: Sparskit solver type: Direct versions of Quasi Generalize Minimum Residual method'
       workspace_size = n + (m+1) * (2*n+4)
    case default
       write(message(1), '(a,i4,a)') "Input: '", sk%solver_type, &
            "' is not a valid Sparskit Solver"
       message(2) = '( Sparskit Solver =  cg | cgnr | bcg | dbcg | bcgstab | tfqmr | fom | gmres | fgmres | dqgmres )'
       call write_fatal(2)
    end select
    call write_info(1)

    ! Now we initialize the arrays for the reverse communication protocol
    sk%ipar = 0
    sk%fpar = 0

    ! A call to the solver with ipar(1) == 0 will initialize the iterative solver.
    sk%ipar(1) = 0

    ! Stopping criteria; use convergence test scheme 1
    sk%ipar(3) = 1 

    sk%ipar(4) = workspace_size
    sk%ipar(5) = sk%krylov_size

    ! Maximum number of matrix-vector multiplies
    sk%ipar(6) = sk%maxiter

    ! Relative tolerance
    sk%fpar(1) = sk%rel_tolerance

    ! Absolute tolerance
    sk%fpar(2) = sk%abs_tolerance


    allocate(work(1:workspace_size))
    work = 0

    call pop_sub
  end subroutine sparskit_solver_init


  ! ---------------------------------------------------------
  subroutine sparskit_solver_run(sk, rhs, sol, op)
    type(sparskit_solver_type), intent(inout) :: sk
    FLOAT, intent(in)  :: rhs(:)
    FLOAT, intent(out) :: sol(:)
    interface
       subroutine op(x, y)
         FLOAT, intent(in)  :: x(:)
         FLOAT, intent(out) :: y(:)
       end subroutine op
    end interface

    call push_sub('sparskit_solver_run')

    select case(sk%solver_type)
    case(SK_CG)
       call sk_driver_cg(sk, rhs, sol, work, op)
    case(SK_CGNR)
       call sk_driver_cgnr(sk, rhs, sol, work, op)
    case(SK_BCG)
       call sk_driver_bcg(sk, rhs, sol, work, op)
    case(SK_DBCG)
       call sk_driver_dbcg(sk, rhs, sol, work, op)
    case(SK_BCGSTAB)
       call sk_driver_bcgstab(sk, rhs, sol, work, op)
    case(SK_TFQMR)
       call sk_driver_tfqmr(sk, rhs, sol, work, op)
    case(SK_FOM)
       call sk_driver_fom(sk, rhs, sol, work, op)
    case(SK_GMRES)
       call sk_driver_gmres(sk, rhs, sol, work, op)
    case(SK_FGMRES)
       call sk_driver_fgmres(sk, rhs, sol, work, op)
    case(SK_DQGMRES)
       call sk_driver_dqgmres(sk, rhs, sol, work, op)
    case default
       write(message(1), '(a,i4,a)') "Input: '", sk%solver_type, &
            "' is not a valid Sparsekit Solver"
       message(2) = '( Sparsekit Solver =  cg | cgnr | bcg | dbcg | bcgstab | tfqmr | fom | gmres | fgmres | dqgmres )'
       call write_fatal(2)
    end select

    call pop_sub
  end subroutine sparskit_solver_run


  ! ---------------------------------------------------------
  subroutine sparskit_solver_end()
    call push_sub('sparskit_solver_end')

    deallocate(work)

    call pop_sub
  end subroutine sparskit_solver_end



#define SOLVER cg
#define DRIVER sk_driver_cg
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER cgnr
#define DRIVER sk_driver_cgnr
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER bcg
#define DRIVER sk_driver_bcg
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER dbcg
#define DRIVER sk_driver_dbcg
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER bcgstab
#define DRIVER sk_driver_bcgstab
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER tfqmr
#define DRIVER sk_driver_tfqmr
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER fom
#define DRIVER sk_driver_fom
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER gmres
#define DRIVER sk_driver_gmres
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER fgmres
#define DRIVER sk_driver_fgmres
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER
#define SOLVER dqgmres
#define DRIVER sk_driver_dqgmres
#include "sparskit_inc.F90"
#undef  SOLVER
#undef  DRIVER

end module sparskit


! ---------------------------------------------------------
FLOAT function distdot(n,x,ix,y,iy)
  use blas
!  use lib_basic_alg

  integer, intent(in) :: n, ix, iy
  FLOAT, intent(inout)   :: x, y
  
  distdot = ddot(n,x,ix,y,iy)
!  distdot = lalg_dot(n, x, y)
  
end function distdot

