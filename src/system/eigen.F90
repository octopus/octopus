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
!! $Id: eigen.F90 6342 2010-03-04 01:25:36Z dstrubbe $

#include "global.h"

module eigensolver_m
  use batch_m
  use datasets_m
  use derivatives_m
  use eigen_cg_m
  use eigen_lobpcg_m
  use eigen_rmmdiis_m
  use exponential_m
  use global_m
  use grid_m
  use hamiltonian_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_calc_m
  use states_dim_m
  use subspace_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::            &
    eigensolver_t,     &
    eigensolver_init,  &
    eigensolver_end,   &
    eigensolver_run

  type eigensolver_t
    integer :: es_type    !< which eigensolver to use

    logical :: subspace_diag
    
    FLOAT   :: tolerance
    integer :: es_maxiter

    integer :: arnoldi_vectors
    FLOAT   :: imag_time

    !> Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer        :: matvec
    integer, pointer :: converged(:)

    !> Stores information about the preconditioning.
    type(preconditioner_t) :: pre

    type(subspace_t) :: sdiag

    integer :: rmmdiis_minimization_iter
  end type eigensolver_t


  integer, public, parameter :: &
       RS_PLAN    = 11,         &
       RS_CG      =  5,         &
       RS_MG      =  7,         &
       RS_CG_NEW  =  6,         &
       RS_EVO     =  9,         &
       RS_LOBPCG  =  8,         &
       RS_RMMDIIS = 10,         &
       RS_ARPACK  = 12,         &  
       RS_DIRECT  = 14

contains

  ! ---------------------------------------------------------
  subroutine eigensolver_init(eigens, gr, st)
    type(eigensolver_t), intent(out)   :: eigens
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st

    integer :: default_iter, default_es

    PUSH_SUB(eigensolver_init)

    !%Variable Eigensolver
    !%Type integer
    !%Section SCF::Eigensolver
    !%Description
    !% Which eigensolver to use to obtain the lowest eigenvalues and
    !% eigenfunctions of the Kohn-Sham Hamiltonian. The default is
    !% conjugate gradients (<tt>cg</tt>); when parallelization in states is
    !% enabled, the default is <tt>lobpcg</tt>.
    !%Option cg 5
    !% Conjugate-gradients algorithm.
    !%Option plan 11
    !% Preconditioned Lanczos scheme.
    !%Option cg_new 6
    !% An alternative conjugate-gradients eigensolver, faster for
    !% larger systems but less mature.
    !%Option evolution 9
    !% Propagation in imaginary time. WARNING: Sometimes it misbehaves. Use with 
    !% caution.
    !%Option lobpcg 8
    !% (Experimental) Locally optimal block-preconditioned
    !% conjugate-gradient algorithm. Ref: A. Knyazev, Toward the
    !% Optimal Preconditioned Eigensolver: Locally Optimal Block
    !% Preconditioned Conjugate Gradient Method, <i>SIAM Journal on
    !% Scientific Computing</i>, 23(2):517-541, 2001.  
    !%Option rmmdiis 10 
    !% Residual minimization scheme, direct inversion in the
    !% iterative subspace eigensolver, based on the implementation of
    !% Kresse and Furthmüller [<i>Phys. Rev. B</i> <b>54</b>, 11169
    !% (1996)]. This eigensolver requires almost no orthogonalization
    !% so it can be considerably faster than the other options for
    !% large systems; however it might suffer stability problems. To
    !% improve its performance a large number of <tt>ExtraStates</tt>
    !% are required (around 10-20% of the number of occupied states).
    !%Option multigrid 7
    !% (Experimental) Multigrid eigensolver.
    !%Option arpack 12
    !% Implicitly Restarted Arnoldi Method. Requires the ARPACK package. 
    !% method.
    !%End

    if(st%parallel_in_states) then
      default_es = RS_LOBPCG
    else
      default_es = RS_CG
    endif

    call parse_integer(datasets_check('Eigensolver'), default_es, eigens%es_type)

    if(st%parallel_in_states .and. .not. eigensolver_parallel_in_states(eigens)) then
      message(1) = "The selected eigensolver is not parallel in states."
      message(2) = "Please use the lobpcg or rmmdiis eigensolvers."
      call messages_fatal(2)
    end if

    call messages_obsolete_variable('EigensolverVerbose')

    !%Variable EigensolverSubspaceDiag
    !%Type logical
    !%Default yes
    !%Section SCF::Eigensolver
    !%Description
    !% Allows you to turn off subspace diagonalization during the diagonalization of
    !% the Hamiltonian. Subspace diagonalization sometimes creates problems when restarting
    !% unoccupied-states calculations with a larger number of unoccupied states.
    !%End
    call parse_logical(datasets_check('EigensolverSubspaceDiag'), .true., eigens%subspace_diag)

    default_iter = 25

    select case(eigens%es_type)
    case(RS_CG_NEW)
    case(RS_MG)
    case(RS_CG)
    case(RS_PLAN)
    case(RS_DIRECT)
    case(RS_EVO)
      !%Variable EigensolverImaginaryTime
      !%Type float
      !%Default 10.0
      !%Section SCF::Eigensolver
      !%Description
      !% The imaginary-time step that is used in the imaginary-time evolution
      !% method (<tt>Eigensolver = evolution</tt>) to obtain the lowest eigenvalues/eigenvectors.
      !% It must satisfy <tt>EigensolverImaginaryTime > 0</tt>.
      !%End
      call parse_float(datasets_check('EigensolverImaginaryTime'), CNST(10.0), eigens%imag_time)
      if(eigens%imag_time <= M_ZERO) call input_error('EigensolverImaginaryTime')
    case(RS_LOBPCG)
    case(RS_RMMDIIS)
      default_iter = 3

      !%Variable EigensolverMinimizationIter
      !%Type integer
      !%Default 5
      !%Section SCF::Eigensolver
      !%Description
      !% During the first iterations, the RMMDIIS eigensolver requires
      !% some steepest-descent minimizations to improve
      !% convergence. This variable determines the number of those
      !% minimizations. The default is 5.
      !%End

      call parse_integer(datasets_check('EigensolverMinimizationIter'), 5, eigens%rmmdiis_minimization_iter)

      if(gr%mesh%use_curvilinear) call messages_experimental("RMMDIIS eigensolver for curvilinear coordinates")

#if defined(HAVE_ARPACK) 
    case(RS_ARPACK) 
      	 	 
      !%Variable EigenSolverArnoldiVectors 
      !%Type integer 
      !%Default 20 
      !%Section SCF::EigenSolver 
      !%Description 
      !% This indicates how many Arnoldi vectors are generated 
      !% It must satisfy EigenSolverArnoldiVectors - Number Of Eigenvectors >= 2. 
      !% See the ARPACK documentation for more details. It will default to  
      !% twice the number of eigenvectors (which is the number of states) 
      !%End 
      call parse_integer(datasets_check('EigenSolverArnoldiVectors'), 2*st%nst, eigens%arnoldi_vectors) 
      if(eigens%arnoldi_vectors-st%nst < (M_TWO - st%nst)) call input_error('EigenSolverArnoldiVectors') 
      	 	 
      ! Arpack is not working in some cases, so let us check. 
      if(st%d%ispin .eq. SPINORS) then 
        write(message(1), '(a)') 'The ARPACK diagonalizer does not handle spinors (yet).' 
        write(message(2), '(a)') 'Please provide a different EigenSolver.' 
        call messages_fatal(2) 
      end if 
#endif 

    case default
      call input_error('Eigensolver')
    end select
    call messages_print_var_option(stdout, "Eigensolver", eigens%es_type)

    call messages_obsolete_variable('EigensolverInitTolerance', 'EigensolverTolerance')
    call messages_obsolete_variable('EigensolverFinalTolerance', 'EigensolverTolerance')
    call messages_obsolete_variable('EigensolverFinalToleranceIteration')

    !%Variable EigensolverTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::Eigensolver
    !%Description
    !% This is the tolerance for the eigenvectors. The default is 1e-6.
    !%End
    call parse_float(datasets_check('EigensolverTolerance'), CNST(1.0e-6), eigens%tolerance)

    !%Variable EigensolverMaxIter
    !%Type integer
    !%Section SCF::Eigensolver
    !%Description
    !% Determines the maximum number of iterations that the
    !% eigensolver will perform if the desired tolerance is not
    !% achieved. The default is 25 iterations for all eigensolvers
    !% except for <tt>rmdiis</tt>, which performs only 3 iterations (only
    !% increase it if you know what you are doing).
    !%End
    call parse_integer(datasets_check('EigensolverMaxIter'), default_iter, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) call input_error('EigensolverMaxIter')

    if(eigens%es_maxiter > default_iter) then
      call messages_write('You have specified a large number of eigensolver iterations (')
      call messages_write(eigens%es_maxiter)
      call messages_write(').', new_line = .true.)
      call messages_write('This is not a good idea as it might slow down convergence, even for', new_line = .true.)
      call messages_write('independent particles, as subspace diagonalization will not be used', new_line = .true.)
      call messages_write('often enough.')
      call messages_warning()
    end if
    
    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_DIRECT, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_init(eigens%pre, gr)
    case default
      call preconditioner_null(eigens%pre)
    end select

    nullify(eigens%diff)
    SAFE_ALLOCATE(eigens%diff(1:st%nst, 1:st%d%nik))

    SAFE_ALLOCATE(eigens%converged(1:st%d%nik))
    eigens%converged(1:st%d%nik) = 0
    eigens%matvec    = 0
    
    call subspace_init(eigens%sdiag, st)
        
    POP_SUB(eigensolver_init)

  end subroutine eigensolver_init


  ! ---------------------------------------------------------
  subroutine eigensolver_end(eigens)
    type(eigensolver_t), intent(inout) :: eigens

    PUSH_SUB(eigensolver_end)

    call subspace_end(eigens%sdiag)

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_end(eigens%pre)
    end select

    SAFE_DEALLOCATE_P(eigens%converged)
    SAFE_DEALLOCATE_P(eigens%diff)

    POP_SUB(eigensolver_end)
  end subroutine eigensolver_end


  ! ---------------------------------------------------------
  subroutine eigensolver_run(eigens, gr, st, hm, iter, conv)
    type(eigensolver_t),  intent(inout) :: eigens
    type(grid_t),         intent(in)    :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    integer,              intent(in)    :: iter
    logical,    optional, intent(inout) :: conv

    integer :: maxiter, ik, ns, idir
    FLOAT :: kpoint(1:MAX_DIM)
#ifdef HAVE_MPI
    logical :: conv_reduced
    integer :: outcount, ist
    FLOAT, allocatable :: ldiff(:), leigenval(:)
#endif
    type(profile_t), save :: prof
    character(len=100) :: str_tmp

    call profiling_in(prof, "EIGEN_SOLVER")
    PUSH_SUB(eigensolver_run)

    if(present(conv)) conv = .false.

    eigens%matvec = 0

    ns = 1

    if(st%d%nspin == 2) ns = 2

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens)) then
      call loct_progress_bar(-1, st%nst*st%d%nik)
    end if

    ik_loop: do ik = st%d%kpt%start, st%d%kpt%end
      maxiter = eigens%es_maxiter
      
      if(eigens%subspace_diag .and. eigens%converged(ik) == 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES &
       .and. eigens%es_type /= RS_DIRECT) then
        if (states_are_real(st)) then
          call dsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        else
          call zsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        end if
      end if

      if (states_are_real(st)) then
        
        select case(eigens%es_type)
        case(RS_CG_NEW)
          call deigensolver_cg2_new(gr, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_CG)
          call deigensolver_cg2(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_PLAN)
          call deigensolver_plan(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_EVO)
          call deigensolver_evolution(gr, st, hm, eigens%tolerance, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
        case(RS_LOBPCG)
          call deigensolver_lobpcg(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
        case(RS_MG)
          call deigensolver_mg(gr%der, st, hm, eigens%sdiag, eigens%tolerance, maxiter, &
            eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter <= eigens%rmmdiis_minimization_iter) then
            call deigensolver_rmmdiis_min(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik)
          else
            call deigensolver_rmmdiis(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik, eigens%diff(:, ik))
          end if
#if defined(HAVE_ARPACK) 
 	      case(RS_ARPACK) 
 	        call deigen_solver_arpack(gr, st, hm, eigens%tolerance, maxiter, eigens%arnoldi_vectors, & 
 	             eigens%converged(ik), ik, eigens%diff(:,ik)) 
#endif 
        end select

        if(eigens%subspace_diag.and.eigens%es_type /= RS_RMMDIIS .and. eigens%es_type /= RS_ARPACK) then
          call dsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        end if

      else

        select case(eigens%es_type)
        case(RS_CG_NEW)
          call zeigensolver_cg2_new(gr, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_CG)
          call zeigensolver_cg2(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_DIRECT)
           call eigensolver_direct(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_PLAN)
          call zeigensolver_plan(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_EVO)
          call zeigensolver_evolution(gr, st, hm, eigens%tolerance, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
        case(RS_LOBPCG)
          call zeigensolver_lobpcg(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
            eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
        case(RS_MG)
          call zeigensolver_mg(gr%der, st, hm, eigens%sdiag, eigens%tolerance, maxiter, &
            eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter <= eigens%rmmdiis_minimization_iter) then
            call zeigensolver_rmmdiis_min(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik)
          else
            call zeigensolver_rmmdiis(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik,  eigens%diff(:, ik))
          end if         
#if defined(HAVE_ARPACK) 
       	case(RS_ARPACK) 
 	        call zeigen_solver_arpack(gr, st, hm, eigens%tolerance, maxiter, eigens%arnoldi_vectors, & 
 	             eigens%converged(ik), ik, eigens%diff(:,ik)) 
#endif 
        end select

        if(eigens%subspace_diag.and.eigens%es_type /= RS_RMMDIIS .and.eigens%es_type /= RS_ARPACK &
          .and.eigens%es_type /= RS_DIRECT ) then
          call zsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
         
        end if

      end if

      eigens%matvec = eigens%matvec + maxiter
    end do ik_loop
    
    ! If we complex scale H the eigenstates need to be otrhonormalized with respect to the c-product.
    ! Moreover the eigenvalues ordering need to be imposed as there is no eigensolver 
    ! supporting this ordering (yet).
    if(hm%cmplxscl) then !cmplxscl
      call states_sort_complex(gr%mesh, st, eigens%diff )
      call states_orthogonalize_cproduct(st, gr%mesh)
    end if
    
    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens)) then
      write(stdout, '(1x)')
    end if

    if(present(conv)) conv = all(eigens%converged(st%d%kpt%start:st%d%kpt%end) == st%nst)

#ifdef HAVE_MPI
    if(st%d%kpt%parallel) then
      if(present(conv)) then
        call MPI_Allreduce(conv, conv_reduced, 1, MPI_LOGICAL, MPI_LAND, st%d%kpt%mpi_grp%comm, mpi_err)
        conv = conv_reduced
      end if

      ! every node needs to know all eigenvalues (and diff)
      SAFE_ALLOCATE(ldiff(1:st%d%kpt%nlocal))
      SAFE_ALLOCATE(leigenval(1:st%d%kpt%nlocal))
      do ist = st%st_start, st%st_end
        ldiff(1:st%d%kpt%nlocal) = eigens%diff(ist, st%d%kpt%start:st%d%kpt%end)
        leigenval(1:st%d%kpt%nlocal) = st%eigenval(ist, st%d%kpt%start:st%d%kpt%end)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, ldiff, outcount, &
                                 eigens%diff(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount.eq.st%d%nik)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, leigenval, outcount, &
                                 st%eigenval(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount.eq.st%d%nik)
      end do
      SAFE_DEALLOCATE_A(ldiff)
      SAFE_DEALLOCATE_A(leigenval)
    end if
#endif

    POP_SUB(eigensolver_run)
    call profiling_out(prof)
  end subroutine eigensolver_run


  ! ---------------------------------------------------------
  logical function eigensolver_parallel_in_states(this) result(par_stat)
    type(eigensolver_t), intent(in) :: this
    
    PUSH_SUB(eigensolver_parallel_in_states)

    par_stat = .false.

    select case(this%es_type)
    case(RS_RMMDIIS, RS_LOBPCG)
      par_stat = .true.
    end select
    
    POP_SUB(eigensolver_parallel_in_states)
  end function eigensolver_parallel_in_states
    

  ! ---------------------------------------------------------
  logical function eigensolver_has_progress_bar(this) result(has)
    type(eigensolver_t), intent(in) :: this

    PUSH_SUB(eigensolver_has_progress_bar)

    has = .false.

    select case(this%es_type)
    case(RS_RMMDIIS, RS_CG, RS_CG_NEW, RS_LOBPCG)
      has = .true.
    end select

    POP_SUB(eigensolver_has_progress_bar)
  end function eigensolver_has_progress_bar
  
#include "undef.F90"
#include "real.F90"
#include "eigen_mg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "eigen_mg_inc.F90"
#include "eigen_plan_inc.F90"
#include "eigen_evolution_inc.F90"

#if defined(HAVE_ARPACK) 
#include "undef.F90" 
#include "real.F90" 
#include "eigen_arpack_inc.F90" 
#include "undef.F90" 
#include "complex.F90" 
#include "eigen_arpack_inc.F90" 
#endif 

end module eigensolver_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
