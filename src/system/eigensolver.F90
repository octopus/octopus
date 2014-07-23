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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module eigensolver_m
  use batch_m
  use datasets_m
  use derivatives_m
  use eigen_arpack_m
  use eigen_cg_m
  use eigen_feast_m
  use eigen_lobpcg_m
  use eigen_rmmdiis_m
  use energy_calc_m
  use exponential_m
  use global_m
  use grid_m
  use hamiltonian_m
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

    FLOAT   :: tolerance
    integer :: es_maxiter

    type(eigen_arpack_t) :: arpack !< arpack solver !!!
    integer :: arnoldi_vectors
    FLOAT   :: current_rel_dens_error !< for the arpack solver it is important to base precision on how well density is converged
    FLOAT   :: imag_time

    !> Stores information about how well it performed.
    FLOAT, pointer   :: diff(:, :)
    integer          :: matvec
    integer, pointer :: converged(:)

    !> Stores information about the preconditioning.
    type(preconditioner_t) :: pre

    type(subspace_t) :: sdiag

    type(feast_t) :: feast !< Contains data used by FEAST solver (at least a linear solver)

    integer :: rmmdiis_minimization_iter

    logical :: save_mem
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
       RS_FEAST   = 13

contains

  subroutine cmplxscl_choose_state_order(eigens, st, gr)
    type(eigensolver_t), intent(inout) :: eigens
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    
    !CMPLX, allocatable :: kinetic_elements(:, :) ! XXXXXX
    integer :: ist
    FLOAT   :: limitvalue
    CMPLX, allocatable :: boundary_norms(:)
    FLOAT, allocatable :: buf(:), score(:)

    PUSH_SUB(cmplxscl_choose_state_order)

    call states_orthogonalize_cproduct(st, gr%mesh)

    !SAF!E_ALL!OCATE(kinetic_elements(1:st%nst, 1))
    SAFE_ALLOCATE(score(1:st%nst))
    SAFE_ALLOCATE(buf(1:st%nst))
    SAFE_ALLOCATE(boundary_norms(1:st%nst))

    call identify_continuum_states(gr%mesh, st, boundary_norms)
        
    buf(:) = abs(boundary_norms(:))
    call sort(buf)
    if(st%cmplxscl%nlocalizedstates == 0) then
      limitvalue = -M_ONE ! all states will be above limitvalue
    else
      limitvalue = min(buf(st%cmplxscl%nlocalizedstates) * CNST(1.01), st%cmplxscl%localizationthreshold)
    end if
    
    do ist=1, st%nst
      if(abs(boundary_norms(ist)) <= limitvalue) then ! XXX spin
        score(ist) = -M_PI + atan(st%zeigenval%Re(ist, 1)) ! all negative values
      else
        score(ist) = M_PI + atan(abs(boundary_norms(ist))) ! all positive values
      end if
    end do
    
    write(message(1), *) 'Partial norm'
    call messages_info(1)
    do ist=1, st%nst
      write(message(1), *) ist, boundary_norms(ist)
      call messages_info(1)
    end do
    write(message(1), *) 'state Re(eps) Im(eps) bnorm^2 score'
    call messages_info(1)
    do ist=1, st%nst
      write(message(1), '(i4,2x,f11.6,f11.6,4x,f9.4,2x,f9.4)') ist, st%zeigenval%Re(ist, 1), st%zeigenval%Im(ist, 1), &
        abs(boundary_norms(ist)), score(ist)
      !'(i4,1x,2f7.3)') ist, st%zeigenval%Re(ist), st%zeigenval%Im(ist)
      !write(message(1), *) ist, score(ist), abs(boundary_norms(ist)), (abs(boundary_norms(ist)) <= limitvalue)
      call messages_info(1)
    end do

    call states_sort_complex(gr%mesh, st, eigens%diff, score)

    call identify_continuum_states(gr%mesh, st, boundary_norms) ! XXX This is just to get correct ordering again
    write(message(1), *) 'and again: Score / abs.part.norm'
    call messages_info(1)
    do ist=1, st%nst
      write(message(1), '(i4,2x,f11.6,f11.6,4x,f9.4,2x)') ist, st%zeigenval%Re(ist, 1), st%zeigenval%Im(ist, 1), &
        abs(boundary_norms(ist))
      !write(message(1), ist, 
      write(message(1), *) ist, score(ist), abs(boundary_norms(ist)), (abs(boundary_norms(ist)) <= limitvalue)
      call messages_info(1)
    end do

    SAFE_DEALLOCATE_A(boundary_norms)
    SAFE_DEALLOCATE_A(buf)
    SAFE_DEALLOCATE_A(score)
    !call cmplxscl_get_kinetic_elements(st, hm, gr%der, kinetic_elements)

    !print*, 'kinetic elements'
    !do ist = 1, st%nst
    !  print*, kinetic_elements(ist, 1)
    !end do

    !print*, 'kinetic elements / eps'
    !do ist = 1, st%nst
    !  print*, kinetic_elements(ist, 1) / (st%zeigenval%Re(ist, 1) + M_zI * st%zeigenval%Im(ist, 1)) !&
    !* exp(-M_TWO * M_zI * st%cmplxscl%theta)
    !end do
    !SAF!E_DEA!LLOC!ATE_A(kinetic_elements)
    POP_SUB(cmplxscl_choose_state_order)
  end subroutine cmplxscl_choose_state_order



  ! XXXX complex scaling.  Move to some better place
  ! This is a hack to identify continuum states.
  ! Continuum states are the only states that do not approach 
  ! zero asymptotically.  Thus they will be very nonzero even far
  ! from the center of the system.  To identify them we therefore
  ! look at a state and check whether a very large fraction of its norm
  ! is contributed within some radius.
  ! If it isn't, then it's considered a continuum state.
  ! Really this method should consider only boundary points, because some continuum
  ! states could potentially also have a large fraction of the norm in the center.
  ! But this is difficult to implement when the box shape is user defined.
  ! Hence we use a fixed radius related to total system size.
  subroutine identify_continuum_states(mesh, st, boundary_norms)
    type(mesh_t),   intent(in)   :: mesh
    type(states_t), intent(in)   :: st
    CMPLX,          intent(out)  :: boundary_norms(:)

    FLOAT              :: r2
    integer            :: ist, ik, ii
    CMPLX, allocatable :: psi(:, :)
    FLOAT, allocatable :: mask(:)
    CMPLX              :: cnorm2

    PUSH_SUB(identify_continuum_states)
    !print*, size(mesh%x, 1)
    !print*, size(mesh%x, 2)
    !print*, mesh%x(1:5, 1)
    !print*, mesh%x(1:5, 2)
    !print*, mesh%x(1:5, 3)

    ! Well, we hardcode that everything is centered on 0.
    ! Then we calculate distances to 0 and use those to weight
    ! the norm of each state.  Ugly ugly...
    SAFE_ALLOCATE(mask(1:mesh%np_part))
    SAFE_ALLOCATE(psi(1:mesh%np_part, 1))

    r2 = st%cmplxscl%localizationradius**2
    do ii=1, mesh%np_part
      if(sum(mesh%x(ii, :)**2) < r2) then
        mask(ii) = M_ZERO
      else
        mask(ii) = M_ONE
      end if
    end do
    !print*, 'mask', sum(mask)

    ASSERT(st%d%dim == 1)
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist=1, st%nst
        call states_get_state(st, mesh, ist, ik, psi)
        psi(:, 1) = psi(:, 1) * mask(:)
        cnorm2 = zmf_dotp(mesh, 1, psi, psi, dotu = .true.)
        !print*, 'cnorm2', ist, cnorm2
        boundary_norms(ist) = cnorm2
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(mask)

    POP_SUB(identify_continuum_states)
  end subroutine identify_continuum_states



  ! ---------------------------------------------------------
  subroutine eigensolver_init(eigens, gr, st)
    type(eigensolver_t), intent(out)   :: eigens
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st

    integer :: default_iter, default_es
    FLOAT   :: default_tol
    real(8) :: mem

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
    !% Preconditioned Lanczos scheme. Ref: Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu and S. Ogut,
    !% "Solution of Large Eigenvalue Problems in Electronic Structure Calculations", <i>BIT</i> <b>36</b>, 1 (1996).
    !%Option cg_new 6
    !% An alternative conjugate-gradients eigensolver, faster for
    !% larger systems but less mature.
    !% Ref: Jiang et al., <i>Phys. Rev. B</i> <b>68</b>, 165337 (2003)
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
    !%Option feast 13
    !% (Experimental) Non-Hermitian FEAST eigensolver. Requires the FEAST
    !% package.
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

    if(eigens%es_type == RS_LOBPCG .and. st%group%block_start /= st%group%block_end) &
      call messages_experimental("lobpcg eigensolver with more than one block per node")

    call messages_obsolete_variable('EigensolverVerbose')
    call messages_obsolete_variable('EigensolverSubspaceDiag', 'SubspaceDiagonalization')

    default_iter = 25
    default_tol = CNST(1e-6)

    select case(eigens%es_type)
    case(RS_CG_NEW)
    case(RS_MG)
      call messages_experimental("multigrid eigensolver")
    case(RS_CG)
    case(RS_PLAN)
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

      !%Variable EigensolverSaveMemory
      !%Type logical
      !%Default no
      !%Section SCF::Eigensolver
      !%Description
      !% The RMMDIIS eigensolver may require a considerable amount of
      !% extra memory. When this variable is set to yes, the
      !% eigensolver will use less memory at the expense of some
      !% performance. This is especially useful for GPUs. The default is no.
      !%End

      call parse_logical(datasets_check('EigensolverSaveMemory'), .false., eigens%save_mem)

      if(gr%mesh%use_curvilinear) call messages_experimental("RMMDIIS eigensolver for curvilinear coordinates")

#if defined(HAVE_ARPACK) 
    case(RS_ARPACK) 

      !%Variable EigensolverArnoldiVectors 
      !%Type integer 
      !%Default 20 
      !%Section SCF::Eigensolver
      !%Description 
      !% For <tt>Eigensolver = arpack</tt>, this indicates how many Arnoldi vectors are generated.
      !% It must satisfy <tt>EigensolverArnoldiVectors</tt> - Number Of Eigenvectors >= 2. 
      !% See the ARPACK documentation for more details. It will default to  
      !% twice the number of eigenvectors (which is the number of states) 
      !%End 
      call parse_integer(datasets_check('EigensolverArnoldiVectors'), 2*st%nst, eigens%arnoldi_vectors) 
      if(eigens%arnoldi_vectors-st%nst < (M_TWO - st%nst)) call input_error('EigensolverArnoldiVectors') 

      eigens%current_rel_dens_error = -M_ONE ! Negative initial value to signify that no value has been assigned yet

      ! Arpack is not working in some cases, so let us check. 
      if(st%d%ispin  ==  SPINORS) then 
        write(message(1), '(a)') 'The ARPACK diagonalizer does not handle spinors (yet).' 
        write(message(2), '(a)') 'Please provide a different Eigensolver.' 
        call messages_fatal(2) 
      end if

      call arpack_init(eigens%arpack, gr, st%nst)

      !Some default values 
      default_iter = 500  ! empirical value based upon experience
      default_tol = M_ZERO ! default is machine precision   
#endif 
    case(RS_FEAST)
      call messages_experimental("FEAST eigensolver")
      call feast_init(eigens%feast, gr, st%nst)

    case default
      call input_error('Eigensolver')
    end select

    call messages_print_stress(stdout, 'Eigensolver')

    call messages_print_var_option(stdout, "Eigensolver", eigens%es_type)

    call messages_obsolete_variable('EigensolverInitTolerance', 'EigensolverTolerance')
    call messages_obsolete_variable('EigensolverFinalTolerance', 'EigensolverTolerance')
    call messages_obsolete_variable('EigensolverFinalToleranceIteration')

    !%Variable EigensolverTolerance
    !%Type float
    !%Default 1.0e-6
    !%Section SCF::Eigensolver
    !%Description
    !% This is the tolerance for the eigenvectors. The default is 1e-6
    !% except for the ARPACK solver for which it is 0.
    !%End
    call parse_float(datasets_check('EigensolverTolerance'), default_tol, eigens%tolerance)

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
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_init(eigens%pre, gr)
    case default
      call preconditioner_null(eigens%pre)
    end select

    nullify(eigens%diff)
    SAFE_ALLOCATE(eigens%diff(1:st%nst, 1:st%d%nik))
    eigens%diff(1:st%nst, 1:st%d%nik) = 0

    SAFE_ALLOCATE(eigens%converged(1:st%d%nik))
    eigens%converged(1:st%d%nik) = 0
    eigens%matvec = 0

    ! FEAST: subspace diagonalization or not?  I guess not.
    ! But perhaps something could be gained by changing this.
    call subspace_init(eigens%sdiag, st, no_sd = (eigens%es_type == RS_ARPACK).or.(eigens%es_type == RS_FEAST))

    ! print memory requirements
    select case(eigens%es_type)
    case(RS_RMMDIIS)
      call messages_write('Info: The rmmdiis eigensolver requires ')
      mem = (2.0_8*eigens%es_maxiter - 1.0_8)*st%d%block_size*dble(gr%mesh%np_part)
      if(states_are_real(st)) then
        mem = mem*CNST(8.0)
      else
        mem = mem*CNST(16.0)
      end if
      call messages_write(mem, units = unit_megabytes, fmt = '(f9.1)')
      call messages_write(' of additional')
      call messages_new_line()
      call messages_write('      memory.  This amount can be reduced by decreasing the value')
      call messages_new_line()
      call messages_write('      of the variable StatesBlockSize (currently set to ')
      call messages_write(st%d%block_size)
      call messages_write(').')
      call messages_info()
    end select

    call messages_print_stress(stdout)

    POP_SUB(eigensolver_init)
  end subroutine eigensolver_init


  ! ---------------------------------------------------------
  subroutine eigensolver_end(eigens)
    type(eigensolver_t), intent(inout) :: eigens

    PUSH_SUB(eigensolver_end)

    select case(eigens%es_type)
    case(RS_PLAN, RS_CG, RS_LOBPCG, RS_RMMDIIS)
      call preconditioner_end(eigens%pre)

    case(RS_FEAST)
      call feast_end(eigens%feast)
    
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
    logical,    optional, intent(out)   :: conv

    integer :: maxiter, ik, ns, ist
#ifdef HAVE_MPI
    logical :: conv_reduced
    integer :: outcount, lmatvec
    FLOAT, allocatable :: ldiff(:), leigenval(:)
    integer, allocatable :: lconv(:)
#endif
    type(profile_t), save :: prof

    call profiling_in(prof, "EIGEN_SOLVER")
    PUSH_SUB(eigensolver_run)

    if(present(conv)) conv = .false.

    eigens%matvec = 0

    ns = 1

    if(st%d%nspin == 2) ns = 2

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens) .and. .not. in_debug_mode) then
      call loct_progress_bar(-1, st%nst*st%d%nik)
    end if

    ik_loop: do ik = st%d%kpt%start, st%d%kpt%end
      maxiter = eigens%es_maxiter

      if(eigens%converged(ik) == 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
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
          call deigensolver_mg(gr%der, st, hm, eigens%sdiag, maxiter, ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter <= eigens%rmmdiis_minimization_iter) then
            call deigensolver_rmmdiis_min(gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
          else
            call deigensolver_rmmdiis(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik, eigens%diff(:, ik), eigens%save_mem)
          end if
        case(RS_ARPACK)
          ! We don`t have any tests of this presently and would not like
          ! to guarantee that it works right now
          call messages_not_implemented('ARPACK solver for Hermitian problems')
        case(RS_FEAST)
          ! same as for ARPACK
          call messages_not_implemented('FEAST solver for Hermitian problems')
        end select

        ! FEAST: subspace diag or not?
        if(eigens%es_type /= RS_RMMDIIS .and. eigens%es_type /= RS_ARPACK) then
          call dsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        end if

      else

        select case(eigens%es_type)
        case(RS_CG_NEW)
          call zeigensolver_cg2_new(gr, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
        case(RS_CG)
          call zeigensolver_cg2(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
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
          call zeigensolver_mg(gr%der, st, hm, eigens%sdiag, maxiter, ik, eigens%diff(:, ik))
        case(RS_RMMDIIS)
          if(iter <= eigens%rmmdiis_minimization_iter) then
            call zeigensolver_rmmdiis_min(gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
          else
            call zeigensolver_rmmdiis(gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
              eigens%converged(ik), ik,  eigens%diff(:, ik), eigens%save_mem)
          end if
#if defined(HAVE_ARPACK) 
       	case(RS_ARPACK)
          call zeigensolver_arpack(eigens%arpack, gr, st, hm, eigens%tolerance, eigens%current_rel_dens_error, maxiter, & 
            eigens%converged(ik), ik, eigens%diff(:,ik))
#endif 
        case(RS_FEAST)
          call zeigensolver_feast(eigens%feast, gr, st, hm, eigens%converged(ik), ik, eigens%diff(:, ik))
        end select

        if(eigens%es_type /= RS_RMMDIIS .and.eigens%es_type /= RS_ARPACK) then
          call zsubspace_diag(eigens%sdiag, gr%der, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        end if

      end if

      ! recheck convergence after subspace diagonalization, since states may have reordered
      eigens%converged(ik) = 0
      do ist = 1, st%nst
        if(eigens%diff(ist, ik) < eigens%tolerance) then
          eigens%converged(ik) = ist
        else
          exit
        endif
      enddo

      eigens%matvec = eigens%matvec + maxiter
    end do ik_loop

    ! If we complex scale H the eigenstates need to be orthonormalized with respect to the c-product.
    ! Moreover the eigenvalues ordering need to be imposed as there is no eigensolver 
    ! supporting this ordering (yet).
    if(st%cmplxscl%space) then
      call cmplxscl_choose_state_order(eigens, st, gr)
    end if

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens) .and. .not. in_debug_mode) then
      write(stdout, '(1x)')
    end if

    if(present(conv)) conv = all(eigens%converged(st%d%kpt%start:st%d%kpt%end) == st%nst)

#ifdef HAVE_MPI
    if(st%d%kpt%parallel) then
      if(present(conv)) then
        call MPI_Allreduce(conv, conv_reduced, 1, MPI_LOGICAL, MPI_LAND, st%d%kpt%mpi_grp%comm, mpi_err)
        conv = conv_reduced
      end if

      lmatvec = eigens%matvec
      call MPI_Allreduce(lmatvec, eigens%matvec, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)

      SAFE_ALLOCATE(lconv(1:st%d%kpt%nlocal))
      lconv(1:st%d%kpt%nlocal) = eigens%converged(st%d%kpt%start:st%d%kpt%end)
      call lmpi_gen_allgatherv(st%d%kpt%nlocal, lconv, outcount, eigens%converged, st%d%kpt%mpi_grp)
      ASSERT(outcount == st%d%nik)
      SAFE_DEALLOCATE_A(lconv)

      ! every node needs to know all eigenvalues (and diff)
      SAFE_ALLOCATE(ldiff(1:st%d%kpt%nlocal))
      SAFE_ALLOCATE(leigenval(1:st%d%kpt%nlocal))
      do ist = st%st_start, st%st_end
        ldiff(1:st%d%kpt%nlocal) = eigens%diff(ist, st%d%kpt%start:st%d%kpt%end)
        leigenval(1:st%d%kpt%nlocal) = st%eigenval(ist, st%d%kpt%start:st%d%kpt%end)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, ldiff, outcount, &
          eigens%diff(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, leigenval, outcount, &
          st%eigenval(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
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

  !#if defined(HAVE_ARPACK) 
  !#include "undef.F90" 
  !#include "real.F90" 
  !#include "eigen_arpack_inc.F90" 
  !#include "undef.F90" 
  !#include "complex.F90" 
  !#include "eigen_arpack_inc.F90" 
  !#endif 

end module eigensolver_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
