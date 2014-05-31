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
    sparskit_solver_init,        &
    dsparskit_solver_run,        &
    zsparskit_solver_run,        &
    sparskit_solver_end

  type sparskit_solver_t
    logical :: is_complex           !< whether set up for complex (otherwise real)
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

  ! ---------------------------------------------------------
  subroutine sparskit_solver_init(n, sk, is_complex)
    integer,                 intent(in)  :: n
    type(sparskit_solver_t), intent(out) :: sk
    logical,                 intent(in)  :: is_complex

    integer :: workspace_size, m

    PUSH_SUB(sparskit_solver_init)

    sk%is_complex = is_complex
    ! there might be some incompatibilities to check of real/complex and available methods?

    !%Variable SPARSKITSolver
    !%Type integer
    !%Default sk_bcg
    !%Section Math::SPARSKIT
    !%Description
    !% Specifies what kind of linear solver will be used.
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
    !% Direct versions of the Quasi-Generalized Minimum Residual method
    !%End
    call parse_integer(datasets_check('SPARSKITSolver'), SK_BCG, sk%solver_type)
    if ( sk%solver_type < SK_MINVAL.or.sk%solver_type > SK_MAXVAL ) then
      call input_error('SPARSKITSolver')
    end if

    !%Variable SPARSKITKrylovSubspaceSize
    !%Type integer
    !%Default 15
    !%Section Math::SPARSKIT
    !%Description
    !% Some of the SPARSKIT solvers are Krylov subspace methods.
    !% This variable determines what size the solver will use 
    !% for the subspace.
    !%End
    call parse_integer(datasets_check('SPARSKITKrylovSubspaceSize'), 15, sk%krylov_size)

    ! preconditioner not implemented
    sk%preconditioning = 0

    !%Variable SPARSKITMaxIter
    !%Type integer
    !%Default 50000
    !%Section Math::SPARSKIT
    !%Description
    !% This variable controls the maximum number of iteration steps that
    !% will be performed by the (iterative) linear solver.
    !%End
    call parse_integer(datasets_check('SPARSKITMaxIter'),          5000, sk%maxiter)
    
    !%Variable SPARSKITIterOut
    !%Type integer
    !%Default -1
    !%Section Math::SPARSKIT
    !%Description
    !% Determines how often status info of the solver is printed.
    !% If <= 0, will never be printed.
    !%End
    call parse_integer(datasets_check('SPARSKITIterOut'),            -1, sk%iter_out)

    !%Variable SPARSKITRelTolerance
    !%Type float
    !%Default 1e-8
    !%Section Math::SPARSKIT
    !%Description
    !% Some SPARSKIT solvers use a relative tolerance as a stopping criterion 
    !% for the iterative solution process. This variable can be used to 
    !% specify the tolerance.
    !%End
    call parse_float(datasets_check('SPARSKITRelTolerance'), CNST(1e-8), sk%rel_tolerance)
    
    !%Variable SPARSKITAbsTolerance
    !%Type float
    !%Default 1e-10
    !%Section Math::SPARSKIT
    !%Description
    !% Some SPARSKIT solvers use an absolute tolerance as a stopping criterion 
    !% for the iterative solution process. This variable can be used to 
    !% specify the tolerance.
    !%End
    call parse_float(datasets_check('SPARSKITAbsTolerance'), CNST(1e-10), sk%abs_tolerance)

    !%Variable SPARSKITVerboseSolver
    !%Type logical
    !%Default no
    !%Section Math::SPARSKIT
    !%Description
    !% When set to yes, the SPARSKIT solver will write more detailed output.
    !%End
    call parse_logical(datasets_check('SPARSKITVerboseSolver'), .false., sk%verbose)

    ! size of the problem
    if(is_complex) then
      sk%size = 2*n
    else
      sk%size = n
    endif

    ! initialize workspace size
    workspace_size = 0 

    ! Krylov subspace size
    m = sk%krylov_size
    if (mod(m, 2) /= 0) m = m + 1
    
    select case(sk%solver_type)
    case(SK_CG)
      message(1) = 'Info: SPARSKIT solver type: Conjugate Gradient Method'
      workspace_size = 5*sk%size
    case(SK_CGNR)
      message(1) = 'Info: SPARSKIT solver type: Conjugate Gradient Method (Normal Residual equation)'
      workspace_size = 5*sk%size
    case(SK_BCG)
      message(1) = 'Info: SPARSKIT solver type: Bi-Conjugate Gradient Method'
      workspace_size = 7*sk%size
    case(SK_DBCG)
      message(1) = 'Info: SPARSKIT solver type: BCG with partial pivoting'
      workspace_size = 11*sk%size
    case(SK_BCGSTAB)
      message(1) = 'Info: SPARSKIT solver type: BCG stabilized'
      workspace_size = 8*sk%size
    case(SK_TFQMR)
      message(1) = 'Info: SPARSKIT solver type: Transpose-Free Quasi-Minimum Residual method'
      workspace_size = 11*sk%size
    case(SK_FOM)
      message(1) = 'Info: SPARSKIT solver type: Full Orthogonalization Method'
      workspace_size = (sk%size+3)*(m+2) + (m+1)*m/2
    case(SK_GMRES)
      message(1) = 'Info: SPARSKIT solver type: Generalized Minimum Residual method'
      workspace_size = (sk%size+3)*(m+2) + (m+1)*m/2
    case(SK_FGMRES)
      message(1) = 'Info: SPARSKIT solver type: Flexible version of Generalized Minimum Residual method'
      workspace_size =  2*sk%size*(m+1) + (m+1)*m/2 + 3*m + 2
    case(SK_DQGMRES)
      message(1) = 'Info: SPARSKIT solver type: Direct versions of Quasi-Generalized Minimum Residual method'
      workspace_size = sk%size + (m+1) * (2*sk%size+4)
    case default
      write(message(1), '(a,i4,a)') "Input: '", sk%solver_type, &
        "' is not a valid SPARSKIT Solver"
      message(2) = '( SPARSKIT Solver =  cg | cgnr | bcg | dbcg | bcgstab | tfqmr | fom | gmres | fgmres | dqgmres )'
      call messages_fatal(2)
    end select
    call messages_info(1)

    ! Now we initialize the arrays for the reverse communication protocol
    sk%ipar = 0
    sk%fpar = 0

    ! A call to the solver with ipar(1) == 0 will initialize the iterative solver.
    sk%ipar(1) = 0

    ! Stopping criteria; use convergence test scheme 2
    sk%ipar(3) = 2

    sk%ipar(4) = workspace_size
    sk%ipar(5) = sk%krylov_size

    ! Maximum number of matrix-vector multiplies
    sk%ipar(6) = sk%maxiter

    ! Relative tolerance
    sk%fpar(1) = sk%rel_tolerance

    ! Absolute tolerance
    sk%fpar(2) = sk%abs_tolerance

    ! allocate and initialize work arrays
    SAFE_ALLOCATE(sk%sk_b(1:sk%size))
    SAFE_ALLOCATE(sk%sk_y(1:sk%size))
    SAFE_ALLOCATE(sk%sk_work(1:workspace_size))
    sk%sk_work = M_ZERO
    sk%sk_y    = M_ZERO
    sk%sk_b    = M_ZERO

    POP_SUB(sparskit_solver_init)
  end subroutine sparskit_solver_init

  ! ---------------------------------------------------------
  subroutine sparskit_solver_end(sk)
    type(sparskit_solver_t), intent(inout) :: sk
    
    PUSH_SUB(sparskit_solver_end)
    
    SAFE_DEALLOCATE_A(sk%sk_b)
    SAFE_DEALLOCATE_A(sk%sk_y)
    SAFE_DEALLOCATE_A(sk%sk_work)
    
    POP_SUB(sparskit_solver_end)
  end subroutine sparskit_solver_end

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
