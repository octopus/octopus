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

module eigensolver_m
  use datasets_m
  use eigen_cg_m
  use eigen_lobpcg_m
  use eigen_rmmdiis_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use loct_parser_m
  use varinfo_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use mpi_lib_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_dim_m
  use states_calc_m
  use subspace_m
  use units_m
  use exponential_m

  implicit none

  private
  public ::            &
    eigensolver_t,     &
    eigensolver_init,  &
    eigensolver_end,   &
    eigensolver_run

  type eigensolver_t
    integer :: es_type    ! which eigensolver to use
    logical :: verbose    ! If true, the solver prints additional information.

    FLOAT   :: init_tol
    FLOAT   :: final_tol
    integer :: final_tol_iter
    integer :: es_maxiter

    integer :: arnoldi_vectors
    FLOAT   :: imag_time

    ! Stores information about how well it performed.
    FLOAT, pointer :: diff(:, :)
    integer        :: matvec
    integer, pointer :: converged(:)

    ! Stores information about the preconditioning.
    type(preconditioner_t) :: pre
  end type eigensolver_t


  integer, public, parameter :: &
       RS_PLAN    = 11,         &
       RS_CG      =  5,         &
       RS_MG      =  7,         &
       RS_CG_NEW  =  6,         &
       RS_EVO     =  9,         &
       RS_LOBPCG  =  8,         &
       RS_RMMDIIS = 10

contains

  ! ---------------------------------------------------------
  subroutine eigensolver_init(gr, eigens, st)
    type(grid_t),         intent(inout) :: gr
    type(eigensolver_t), intent(out)   :: eigens
    type(states_t),       intent(in)    :: st

    integer :: default_iter

    call push_sub('eigen.eigensolver_init')

    !%Variable Eigensolver
    !%Type integer
    !%Default cg
    !%Section SCF::Eigensolver
    !%Description
    !% Decides the eigensolver that obtains the lowest eigenvalues and
    !% eigenfunctions of the Kohn-Sham Hamiltonian. The default is
    !% conjugate gradients (cg), when parallelization in states is
    !% enabled the default is lobpcg.
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
    !% Locally optimal block-preconditioned conjugate-gradient algorithm
    !% (only available if DevelVersion = yes),
    !% see: A. Knyazev. Toward the Optimal Preconditioned Eigensolver: Locally
    !% Optimal Block Preconditioned Conjugate Gradient Method. SIAM
    !% Journal on Scientific Computing, 23(2):517??541, 2001.
    !%Option rmmdiis 10
    !% Residual minimization scheme, direct inversion in the iterative
    !% subspace eigensolver, based on the implementation of Kresse and
    !% Furthmüller [Phys. Rev. B 54, 11169 (1996)]. This eigensolver
    !% requires almost no orthogonalization so it can be considerably
    !% faster than the other options for large systems, however it
    !% might suffer stability problems. To improve its performance a
    !% large number of ExtraStates are required (around 10-20% of the
    !% number of occupied states).
    !%Option multigrid 7
    !% Multigrid eigensolver (experimental).
    !%End
    call loct_parse_int(datasets_check('Eigensolver'), RS_CG, eigens%es_type)

    if(st%parallel_in_states .and. .not. eigensolver_parallel_in_states(eigens)) then
      message(1) = "The selected eigensolver is not parallel in states."
      message(2) = "Please use the lobpcg or rmmdiis eigensolvers."
      call write_fatal(2)
    end if

    !%Variable EigensolverVerbose
    !%Type logical
    !%Default no
    !%Section SCF::Eigensolver
    !%Description
    !% If enabled the eigensolver prints additional information.
    !%End
    call loct_parse_logical(datasets_check('EigensolverVerbose'), .false., eigens%verbose)

    default_iter = 25

    select case(eigens%es_type)
    case(RS_CG_NEW)
    case(RS_MG)
    case(RS_CG)
    case(RS_PLAN)
    case(RS_EVO)
      !%Variable EigensolverImaginaryTime
      !%Type float
      !%Default 10.0
      !%Section SCF::Eigensolver
      !%Description
      !% The imaginary-time step that is used in the imaginary-time evolution
      !% method to obtain the lowest eigenvalues/eigenvectors.
      !% It must satisfy EigensolverImaginaryTime > 0.
      !%End
      call loct_parse_float(datasets_check('EigensolverImaginaryTime'), CNST(10.0), eigens%imag_time)
      if(eigens%imag_time <= M_ZERO) call input_error('EigensolverImaginaryTime')
    case(RS_LOBPCG)
    case(RS_RMMDIIS)
      default_iter = 3
      call messages_devel_version("RMMDIIS eigensolver")
    case default
      call input_error('Eigensolver')
    end select
    call messages_print_var_option(stdout, "Eigensolver", eigens%es_type)

    !%Variable EigensolverInitTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::Eigensolver
    !%Description
    !% This is the initial tolerance for the eigenvectors.
    !%End
    call loct_parse_float(datasets_check('EigensolverInitTolerance'), CNST(1.0e-6), eigens%init_tol)
    if(eigens%init_tol < 0) call input_error('EigensolverInitTolerance')

    !%Variable EigensolverFinalTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::Eigensolver
    !%Description
    !% This is the final tolerance for the eigenvectors. Must be smaller than <tt>EigensolverInitTolerance</tt>.
    !%End
    call loct_parse_float(datasets_check('EigensolverFinalTolerance'), CNST(1.0e-6), eigens%final_tol)
    if(eigens%final_tol < 0 .or. eigens%final_tol > eigens%init_tol) call input_error('EigensolverFinalTolerance')

    !%Variable EigensolverFinalToleranceIteration
    !%Type integer
    !%Default 7
    !%Section SCF::Eigensolver
    !%Description
    !% Determines how many iterations are needed 
    !% to go from <tt>EigensolverInitTolerance</tt> to <tt>EigensolverFinalTolerance</tt>.
    !% Must be larger than 1.
    !%End
    call loct_parse_int(datasets_check('EigensolverFinalToleranceIteration'), 7, eigens%final_tol_iter)
    if(eigens%final_tol_iter <= 1) call input_error('EigensolverFinalToleranceIteration')

    !%Variable EigensolverMaxIter
    !%Type integer
    !%Default 25
    !%Section SCF::Eigensolver
    !%Description
    !% Determines the maximum number of iterations that the
    !% eigensolver will perform if the desired tolerance is not
    !% achieved. The default is 25 iterations for all eigensolvers
    !% except for the rmmdiis that only performs one iteration (only
    !% increase it if you know what you are doing).
    !%End
    call loct_parse_int(datasets_check('EigensolverMaxIter'), default_iter, eigens%es_maxiter)
    if(eigens%es_maxiter < 1) call input_error('EigensolverMaxIter')

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_init(eigens%pre, gr)
    case default
      call preconditioner_null(eigens%pre)
    end select

    nullify(eigens%diff)
    SAFE_ALLOCATE(eigens%diff(1:st%nst, 1:st%d%nik))

    SAFE_ALLOCATE(eigens%converged(1:st%d%nik))
    eigens%converged(1:st%d%nik) = 0
    eigens%matvec    = 0

    call pop_sub()

  end subroutine eigensolver_init


  ! ---------------------------------------------------------
  subroutine eigensolver_end(eigens)
    type(eigensolver_t), intent(inout) :: eigens

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_end(eigens%pre)
    end select

    SAFE_DEALLOCATE_P(eigens%converged)
    SAFE_DEALLOCATE_P(eigens%diff)
    nullify(eigens%diff)
  end subroutine eigensolver_end


  ! ---------------------------------------------------------
  subroutine eigensolver_run(eigens, gr, st, hm, iter, conv, verbose)
    type(eigensolver_t), intent(inout) :: eigens
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    integer,              intent(in)    :: iter
    logical,    optional, intent(inout) :: conv
    logical,    optional, intent(in)    :: verbose

    logical :: verbose_
    integer :: maxiter, ik, ns
    FLOAT :: tol
#ifdef HAVE_MPI
    logical :: conv_reduced
    integer :: outcount, ist
    FLOAT, allocatable :: ldiff(:), leigenval(:)
#endif
    type(profile_t), save :: prof

    call profiling_in(prof, "EIGEN_SOLVER")
    call push_sub('eigen.eigensolver_run')

    verbose_ = eigens%verbose; if(present(verbose)) verbose_ = verbose

    if(iter < eigens%final_tol_iter) then
      tol = log(eigens%final_tol/eigens%init_tol)/(eigens%final_tol_iter - 1)*(iter - 1) + &
           log(eigens%init_tol)
      tol = exp(tol)
    else
      tol = eigens%final_tol
    end if

    if(present(conv)) conv = .false.

    eigens%matvec = 0

    ns = 1

    if(st%d%nspin == 2) ns = 2

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens) .and. .not.verbose_) then
      call loct_progress_bar(-1, st%nst*st%d%nik)
    end if

    ik_loop: do ik = st%d%kpt%start, st%d%kpt%end
      maxiter = eigens%es_maxiter
      if(verbose_) then
        if(st%d%nik > ns) then
          write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
               st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
               st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
               st%d%kpoints(3, ik)*units_out%length%factor, ')'
          call write_info(1)
        end if
      end if
      
      if(eigens%converged(ik) == 0) then
        if (states_are_real(st)) then
          call dsubspace_diag(gr, st, hm, ik, st%eigenval(:, ik), st%dpsi(:, :, :, ik), eigens%diff(:, ik))
        else
          call zsubspace_diag(gr, st, hm, ik, st%eigenval(:, ik), st%zpsi(:, :, :, ik), eigens%diff(:, ik))
        end if
      end if

      if (states_are_real(st)) then
        
        select case(eigens%es_type)
        case(RS_CG_NEW)
          call deigensolver_cg2_new(gr, st, hm, tol, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik), verbose = verbose_)
        case(RS_CG)
          call deigensolver_cg2(gr, st, hm, eigens%pre, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), verbose = verbose_)
        case(RS_PLAN)
          call deigensolver_plan(gr, st, hm, eigens%pre, tol, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_EVO)
          call deigensolver_evolution(gr, st, hm, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
        case(RS_LOBPCG)
          call deigensolver_lobpcg(gr, st, hm, eigens%pre, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size, verbose = verbose_)
        case(RS_MG)
          call deigensolver_mg(gr, st, hm, tol, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter == 1) then
            call deigensolver_rmmdiis_start(gr, st, hm, eigens%pre, tol, maxiter, &
                 eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
          else
            call deigensolver_rmmdiis(gr, st, hm, eigens%pre, tol, maxiter, &
                 eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
          end if
        end select

        if(eigens%es_type /= RS_RMMDIIS) then
          call dsubspace_diag(gr, st, hm, ik, st%eigenval(:, ik), st%dpsi(:, :, :, ik), eigens%diff(:, ik))
        end if

      else

        select case(eigens%es_type)
        case(RS_CG_NEW)
          call zeigensolver_cg2_new(gr, st, hm, tol, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik), verbose = verbose_)
        case(RS_CG)
          call zeigensolver_cg2(gr, st, hm, eigens%pre, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), verbose = verbose_)
        case(RS_PLAN)
          call zeigensolver_plan(gr, st, hm, eigens%pre, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_EVO)
          call zeigensolver_evolution(gr, st, hm, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
        case(RS_LOBPCG)
          call zeigensolver_lobpcg(gr, st, hm, eigens%pre, tol, maxiter, &
               eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size, verbose = verbose_)
        case(RS_MG)
          call zeigensolver_mg(gr, st, hm, tol, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter == 1) then
            call zeigensolver_rmmdiis_start(gr, st, hm, eigens%pre, tol, maxiter, &
                 eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
          else
            call zeigensolver_rmmdiis(gr, st, hm, eigens%pre, tol, maxiter, &
                 eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
          end if
        end select

        if(eigens%es_type /= RS_RMMDIIS) then
          call zsubspace_diag(gr, st, hm, ik, st%eigenval(:, ik), st%zpsi(:, :, :, ik), eigens%diff(:, ik))
        end if

      end if

      eigens%matvec = eigens%matvec + maxiter
    end do ik_loop

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

    call pop_sub()
    call profiling_out(prof)
  end subroutine eigensolver_run

  logical function eigensolver_parallel_in_states(this) result(par_stat)
    type(eigensolver_t), intent(in) :: this
    
    par_stat = .false.

    select case(this%es_type)
    case(RS_RMMDIIS, RS_LOBPCG)
      par_stat = .true.
    end select
    
  end function eigensolver_parallel_in_states
    
  logical function eigensolver_has_progress_bar(this) result(has)
    type(eigensolver_t), intent(in) :: this

    has = .false.

    select case(this%es_type)
    case(RS_RMMDIIS, RS_CG, RS_CG_NEW)
      has = .true.
    end select

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

end module eigensolver_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
