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

#include "global.h"

module test_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use calc_mode_par_oct_m
  use clock_oct_m
  use density_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use exponential_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ion_interaction_oct_m
  use io_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use mesh_interpolation_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use simul_box_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use subspace_oct_m
  use electrons_oct_m
  use types_oct_m
  use v_ks_oct_m
  use wfs_elec_oct_m
  use XC_F90(lib_m)
  use xc_oct_m

  implicit none

  type test_parameters_t
    private
    integer :: type
    integer :: repetitions
    integer :: min_blocksize
    integer :: max_blocksize
  end type test_parameters_t

  public :: test_run

contains

  ! ---------------------------------------------------------
  subroutine test_run(namespace)
    type(namespace_t),       intent(in)    :: namespace

    type(test_parameters_t) :: param
    integer :: test_mode

    PUSH_SUB(test_run)

    call messages_obsolete_variable(namespace, 'WhichTest', 'TestMode')

    !%Variable TestMode
    !%Type integer
    !%Default hartree
    !%Section Utilities::oct-test
    !%Description
    !% Decides what kind of test should be performed.
    !%Option hartree 1
    !% Tests the Poisson solvers used to calculate the Hartree potential.
    !%Option derivatives 2
    !% Tests and benchmarks the implementation of the finite-difference operators, used to calculate derivatives.
    !%Option orthogonalization 3
    !% Tests the implementation of the orthogonalization routines.
    !%Option interpolation 4
    !% Test the interpolation routines.
    !%Option ion_interaction 5
    !% Tests the ion-ion interaction routines.
    !%Option projector 6
    !% Tests the code that applies the nonlocal part of the pseudopotentials
    !% in case of spin-orbit coupling
    !%Option dft_u 7
    !% Tests the DFT+U part of the code for projections on the basis.
    !%Option hamiltonian_apply 8
    !% Tests the application of the Hamiltonian, or a part of it
    !%Option density_calc 9
    !% Calculation of the density.
    !%Option exp_apply 10
    !% Tests the exponential of the Hamiltonian
    !%Option boundaries 11
    !% Tests the boundaries conditions
    !%Option subspace_diag 12
    !% Tests the subspace diagonalization
    !%Option batch_ops 13
    !% Tests the batch operations
    !%Option clock 18
    !% Tests for clock
    !%Option density_grad 19
    !% Calculation of the density gradient.
    !%End
    call parse_variable(namespace, 'TestMode', OPTION__TESTMODE__HARTREE, test_mode)

    call messages_obsolete_variable(namespace, 'TestDerivatives', 'TestType')
    call messages_obsolete_variable(namespace, 'TestOrthogonalization', 'TestType')

    !%Variable TestType
    !%Type integer
    !%Default all
    !%Section Utilities::oct-test
    !%Description
    !% Decides on what type of values the test should be performed.
    !%Option real 1
    !% Test for double-precision real functions.
    !%Option complex 2
    !%Option all 3
    !% Tests for double-precision real and complex functions.
    !%End
    call parse_variable(namespace, 'TestType', OPTION__TESTTYPE__ALL, param%type)
    if(param%type < 1 .or. param%type > 5) then
      message(1) = "Invalid option for TestType."
      call messages_fatal(1, only_root_writes = .true.)
    endif

    !%Variable TestRepetitions
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-test
    !%Description
    !% This variable controls the behavior of oct-test for performance
    !% benchmarking purposes. It sets the number of times the
    !% computational kernel of a test will be executed, in order to
    !% provide more accurate timings.
    !%
    !% Currently this variable is used by the <tt>hartree_test</tt>,
    !% <tt>derivatives</tt>, and <tt>projector</tt> tests.
    !%End
    call parse_variable(namespace, 'TestRepetitions', 1, param%repetitions)

    !%Variable TestMinBlockSize
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-test
    !%Description
    !% Some tests can work with multiple blocksizes, in this case of
    !% range of blocksizes will be tested. This variable sets the lower
    !% bound of that range.
    !%
    !% Currently this variable is only used by the derivatives test.
    !%End
    call parse_variable(namespace, 'TestMinBlockSize', 1, param%min_blocksize)

    !%Variable TestMaxBlockSize
    !%Type integer
    !%Default 128
    !%Section Utilities::oct-test
    !%Description
    !% Some tests can work with multiple blocksizes, in this case of
    !% range of blocksizes will be tested. This variable sets the lower
    !% bound of that range.
    !%
    !% Currently this variable is only used by the derivatives test.
    !%End
    call parse_variable(namespace, 'TestMaxBlockSize', 128, param%max_blocksize)

    call messages_print_stress(stdout, "Test mode")
    call messages_print_var_option(stdout, "TestMode", test_mode)
    call messages_print_var_option(stdout, "TestType", param%type)
    call messages_print_var_value(stdout, "TestRepetitions", param%repetitions)
    call messages_print_var_value(stdout, "TestMinBlockSize", param%min_blocksize)
    call messages_print_var_value(stdout, "TestMaxBlockSize", param%max_blocksize)
    call messages_print_stress(stdout)

    select case(test_mode)
    case(OPTION__TESTMODE__HARTREE)
      call test_hartree(param, namespace)
    case(OPTION__TESTMODE__DERIVATIVES)
      call test_derivatives(param, namespace)
    case(OPTION__TESTMODE__ORTHOGONALIZATION)
      call test_orthogonalization(param, namespace)
    case(OPTION__TESTMODE__INTERPOLATION)
      call test_interpolation(param, namespace)
    case(OPTION__TESTMODE__ION_INTERACTION)
      call test_ion_interaction(namespace)
    case(OPTION__TESTMODE__PROJECTOR)
      call test_projector(param, namespace)
    case(OPTION__TESTMODE__DFT_U)
      call test_dft_u(param, namespace)
    case(OPTION__TESTMODE__HAMILTONIAN_APPLY)
      call test_hamiltonian(param, namespace)
    case(OPTION__TESTMODE__DENSITY_CALC)
      call test_density_calc(param, namespace)
    case(OPTION__TESTMODE__DENSITY_GRAD)
      call test_density_grad(param, namespace)
    case(OPTION__TESTMODE__EXP_APPLY)
      call test_exponential(param, namespace)
    case(OPTION__TESTMODE__BOUNDARIES)
      call test_boundaries(param, namespace)
    case(OPTION__TESTMODE__SUBSPACE_DIAG)
      call test_subspace_diagonalization(param, namespace)
    case(OPTION__TESTMODE__BATCH_OPS)
      call test_batch_ops(param, namespace)
    case(OPTION__TESTMODE__CLOCK)
      call test_clock()
    end select

    POP_SUB(test_run)
  end subroutine test_run

  ! ---------------------------------------------------------
  subroutine test_hartree(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    sys => electrons_t(namespace)
    call poisson_test(sys%hm%psolver, sys%gr%mesh, namespace, param%repetitions)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree

 ! ---------------------------------------------------------
  subroutine test_projector(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t), pointer :: epsib
    integer :: itime
    CMPLX, allocatable :: psi(:, :)

    PUSH_SUB(test_projector)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the nonlocal part of the pseudopotential with SOC')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    !Initialize external potential
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)


    !Initialize external potential
    SAFE_ALLOCATE(epsib)
    call sys%st%group%psib(1, 1)%copy_to(epsib)

    call batch_set_zero(epsib)

    do itime = 1, param%repetitions
      call zproject_psi_batch(sys%gr%mesh, sys%gr%der%boundaries, sys%hm%ep%proj,  &
                              sys%hm%ep%natoms, 2, sys%st%group%psib(1, 1), epsib)
    end do

    SAFE_ALLOCATE(psi(sys%gr%mesh%np, sys%st%d%dim))
    do itime = 1, epsib%nst
      call batch_get_state(epsib, itime, sys%gr%mesh%np, psi)
      write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, zmf_nrm2(sys%gr%mesh, 2, psi)
      call messages_info(1)
    end do
    SAFE_DEALLOCATE_A(psi)

    call epsib%end()
    SAFE_DEALLOCATE_P(epsib)
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_projector)
  end subroutine test_projector

  ! ---------------------------------------------------------
  subroutine test_dft_u(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t), pointer :: epsib
    integer :: itime
    type(orbitalbasis_t) :: basis
    FLOAT, allocatable :: ddot(:,:,:), dweight(:,:)
    CMPLX, allocatable :: zdot(:,:,:), zweight(:,:)

    PUSH_SUB(test_dft_u)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing some DFT+U routines')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    SAFE_ALLOCATE(epsib)
    call sys%st%group%psib(1, 1)%copy_to(epsib, copy_data = .true.)

    !Initialize the orbital basis
    call orbitalbasis_init(basis, sys%namespace)
    if (states_are_real(sys%st)) then
      call dorbitalbasis_build(basis, sys%geo, sys%gr%mesh, sys%st%d%kpt, sys%st%d%dim, .false., .false.)
      SAFE_ALLOCATE(dweight(1:basis%orbsets(1)%sphere%np,1:epsib%nst_linear))
      SAFE_ALLOCATE(ddot(1:sys%st%d%dim,1:basis%orbsets(1)%norbs, 1:epsib%nst))
    else
      call zorbitalbasis_build(basis, sys%geo, sys%gr%mesh, sys%st%d%kpt, sys%st%d%dim, .false., .false.)
      call orbitalset_update_phase(basis%orbsets(1), sys%gr%sb, sys%st%d%kpt, (sys%st%d%ispin==SPIN_POLARIZED))
      SAFE_ALLOCATE(zweight(1:basis%orbsets(1)%sphere%np,1:epsib%nst_linear))
      SAFE_ALLOCATE(zdot(1:sys%st%d%dim,1:basis%orbsets(1)%norbs, 1:epsib%nst))
    end if

    do itime = 1, param%repetitions
      call batch_set_zero(epsib)
      if(states_are_real(sys%st)) then
        dweight = M_ONE
        ddot = M_ZERO
        call dorbitalset_get_coeff_batch(basis%orbsets(1), 1, sys%st%group%psib(1, 1), .false., ddot)
        call dorbitalset_add_to_batch(basis%orbsets(1), 1, epsib, .false., dweight)
      else
        zweight = M_ONE
        zdot = M_ZERO
        call zorbitalset_get_coeff_batch(basis%orbsets(1), sys%st%d%dim, sys%st%group%psib(1, 1), .false., zdot)
        call zorbitalset_add_to_batch(basis%orbsets(1), sys%st%d%dim, epsib, .false., zweight)
      end if
    end do

    if(epsib%is_packed()) then
      call epsib%do_unpack(force = .true.)
    end if

    call test_prints_info_batch(sys%st, sys%gr, epsib)

    SAFE_DEALLOCATE_A(dweight)
    SAFE_DEALLOCATE_A(zweight)
    SAFE_DEALLOCATE_A(ddot)
    SAFE_DEALLOCATE_A(zdot)

    call epsib%end()
    SAFE_DEALLOCATE_P(epsib)
    call orbitalbasis_end(basis)
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_dft_u)
  end subroutine test_dft_u

  ! ---------------------------------------------------------
  subroutine test_hamiltonian(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(wfs_elec_t), pointer :: hpsib
    integer :: itime, terms

    PUSH_SUB(test_hamiltonian)

    !%Variable TestHamiltonianApply
    !%Type integer
    !%Default term_all
    !%Section Utilities::oct-test
    !%Description
    !% Decides which part of the Hamiltonian is applied.
    !%Option term_all 0
    !% Apply the full Hamiltonian.
    !%Option term_kinetic 1
    !% Apply only the kinetic operator
    !%Option term_local_potential 2
    !% Apply only the local potential.
    !%Option term_non_local_potential 4
    !% Apply only the non_local potential.
    !%End
    call parse_variable(namespace, 'TestHamiltonianApply', OPTION__TESTHAMILTONIANAPPLY__TERM_ALL, terms)
    if(terms==0) terms = huge(1)


    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the application of the Hamiltonian')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    !Initialize external potential
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%hm, sys%st, sys%geo)

    call boundaries_set(sys%gr%der%boundaries, sys%st%group%psib(1, 1))

    SAFE_ALLOCATE(hpsib)
    call sys%st%group%psib(1, 1)%copy_to(hpsib)

    if(hamiltonian_elec_apply_packed(sys%hm)) then
      call sys%st%group%psib(1, 1)%do_pack()
      call hpsib%do_pack(copy = .false.)
    end if

    do itime = 1, param%repetitions
      if(states_are_real(sys%st)) then
        call dhamiltonian_elec_apply_batch(sys%hm, sys%namespace, sys%gr%mesh, sys%st%group%psib(1, 1), hpsib, terms = terms, &
          set_bc = .false.)
      else
        call zhamiltonian_elec_apply_batch(sys%hm, sys%namespace, sys%gr%mesh, sys%st%group%psib(1, 1), hpsib, terms = terms, &
          set_bc = .false.)
      end if
    end do

    if(hpsib%is_packed()) then
      call hpsib%do_unpack(force = .true.)
    end if

    call test_prints_info_batch(sys%st, sys%gr, hpsib)

    call hpsib%end(copy = .false.)
    SAFE_DEALLOCATE_P(hpsib)
    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_hamiltonian)
  end subroutine test_hamiltonian


  ! ---------------------------------------------------------
  subroutine test_density_calc(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing density calculation')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call density_calc(sys%st, sys%gr, sys%st%rho)
    end do

    write(message(1),'(a,3x, f12.6)') "Norm density  ", dmf_nrm2(sys%gr%mesh, sys%st%rho(:,1))
    call messages_info(1)

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_density_calc)
  end subroutine test_density_calc

  ! ---------------------------------------------------------
  subroutine test_density_grad(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    integer :: ib, iq, idir, np, np_part
    type(wfs_elec_t) :: psib, grad_psib(1:MAX_DIM)
  
    FLOAT,  allocatable :: grad_psi(:, :, :)
    FLOAT,  allocatable :: grad_rho(:, :)
  

    PUSH_SUB(test_density_grad)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing density gradient calculation')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    np = sys%gr%mesh%np
    np_part = sys%gr%mesh%np_part
  
    SAFE_ALLOCATE(grad_psi(1:np, 1:sys%gr%mesh%sb%dim, 1:sys%st%d%dim))
    SAFE_ALLOCATE(grad_rho(1:np, 1:sys%gr%mesh%sb%dim))
    grad_rho = M_ZERO
  

    do iq = sys%st%d%kpt%start, sys%st%d%kpt%end

      do ib = sys%st%group%block_start, sys%st%group%block_end

        call sys%st%group%psib(ib, iq)%copy_to(psib, copy_data = .true.)

        do idir = 1, sys%gr%mesh%sb%dim
          call psib%copy_to(grad_psib(idir))
          if (states_are_real(sys%st) ) then
            call dderivatives_batch_perform(sys%gr%der%grad(idir), sys%gr%der, psib, grad_psib(idir), set_bc = .false.)
          else
            call zderivatives_batch_perform(sys%gr%der%grad(idir), sys%gr%der, psib, grad_psib(idir), set_bc = .false.)
          end if
        end do

        if (states_are_real(sys%st) ) then
          call ddensity_accumulate_grad(sys%gr, sys%st, psib, grad_psib, grad_rho)
        else
          call zdensity_accumulate_grad(sys%gr, sys%st, psib, grad_psib, grad_rho)
        end if

      end do

    end do


    write(message(1),'(a,3x, f12.6)') "Norm grad(1)  ", dmf_nrm2(sys%gr%mesh, grad_rho(:,1))
    write(message(2),'(a,3x, f12.6)') "Norm grad(2)  ", dmf_nrm2(sys%gr%mesh, grad_rho(:,2))
    write(message(3),'(a,3x, f12.6)') "Norm grad(3)  ", dmf_nrm2(sys%gr%mesh, grad_rho(:,3))
    call messages_info(3)

    SAFE_DEALLOCATE_A(grad_psi)
    SAFE_DEALLOCATE_A(grad_rho)


    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_density_grad)
  end subroutine test_density_grad


  ! ---------------------------------------------------------
  subroutine test_boundaries(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call boundaries_set(sys%gr%der%boundaries, sys%st%group%psib(1, 1))
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_density_calc)
  end subroutine test_boundaries


   ! ---------------------------------------------------------
  subroutine test_exponential(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    type(exponential_t) :: te
    integer :: itime

    PUSH_SUB(test_exponential)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing exponential')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh, wfs_type=TYPE_CMPLX)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    !Initialize external potential
    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%hm, sys%st, sys%geo)

    call exponential_init(te, namespace)

    if(hamiltonian_elec_apply_packed(sys%hm)) then
      call sys%st%group%psib(1, 1)%do_pack()
    end if

    do itime = 1, param%repetitions
      call exponential_apply_batch(te, sys%namespace, sys%gr%mesh, sys%hm, sys%st%group%psib(1, 1), CNST(1.0))
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call exponential_end(te)

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_exponential)
  end subroutine test_exponential


  ! ---------------------------------------------------------
  subroutine test_subspace_diagonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime
    type(subspace_t) :: sdiag

    PUSH_SUB(test_subspace_diagonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%pack()
    call hamiltonian_elec_epot_generate(sys%hm, sys%namespace, sys%gr, sys%geo, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%namespace, sys%hm, sys%st, sys%geo)

    call subspace_init(sdiag, sys%namespace, sys%st, no_sd = .false.)

    do itime = 1, param%repetitions
      if(states_are_real(sys%st)) then
        call dsubspace_diag(sdiag, sys%namespace, sys%gr%mesh, sys%st, sys%hm, 1, sys%st%eigenval(:, 1))
      else
        call zsubspace_diag(sdiag, sys%namespace, sys%gr%mesh, sys%st, sys%hm, 1, sys%st%eigenval(:, 1))
      end if
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_subspace_diagonalization)
  end subroutine test_subspace_diagonalization


  ! ---------------------------------------------------------
  subroutine test_batch_ops(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime, ops, ops_default, ist, jst, nst
    type(wfs_elec_t) :: xx, yy
    FLOAT, allocatable :: tmp(:)
    FLOAT, allocatable :: ddotv(:)
    CMPLX, allocatable :: zdotv(:)
    FLOAT, allocatable :: ddot(:,:)
    CMPLX, allocatable :: zdot(:,:)


    PUSH_SUB(test_batch_ops)

    !%Variable TestBatchOps
    !%Type flag
    !%Default ops_axpy + ops_scal + ops_nrm2
    !%Section Utilities::oct-test
    !%Description
    !% Decides which part of the Hamiltonian is applied.
    !%Option ops_axpy bit(1)
    !% Tests batch_axpy operation
    !%Option ops_scal bit(2)
    !% Tests batch_scal operation
    !%Option ops_nrm2 bit(3)
    !% Tests batch_nrm2 operation
    !%Option ops_dotp_matrix bit(4)
    !% Tests X(mesh_batch_dotp_matrix)
    !%Option ops_dotp_self bit(5)
    !% Tests X(mesh_batch_dotp_self)
    !%Option ops_dotp_vector bit(6)
    !% Tests X(mesh_batch_dotp_vector)
    !%End
    ops_default = OPTION__TESTBATCHOPS__OPS_AXPY &
                + OPTION__TESTBATCHOPS__OPS_SCAL &
                + OPTION__TESTBATCHOPS__OPS_NRM2 &
                + OPTION__TESTBATCHOPS__OPS_DOTP_MATRIX &
                + OPTION__TESTBATCHOPS__OPS_DOTP_SELF &
                + OPTION__TESTBATCHOPS__OPS_DOTP_VECTOR

    call parse_variable(namespace, 'TestBatchOps', ops_default, ops)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing batch operations')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    sys => electrons_t(namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_AXPY) /= 0) then
      message(1) = 'Info: Testing axpy'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      do itime = 1, param%repetitions
        call batch_axpy(sys%gr%mesh%np, CNST(0.1), xx, yy)
      end do
      call test_prints_info_batch(sys%st, sys%gr, yy, string = "axpy")

      call xx%end()
      call yy%end()
    end if

    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_SCAL) /= 0) then
      message(1) = 'Info: Testing scal'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      do itime = 1, param%repetitions
        call batch_scal(sys%gr%mesh%np, CNST(0.1), yy)
      end do
      call test_prints_info_batch(sys%st, sys%gr, yy, string="scal")

      call xx%end()
      call yy%end()
    end if

    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_NRM2) /= 0) then
      message(1) = 'Info: Testing nrm2'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      SAFE_ALLOCATE(tmp(1:xx%nst))

      do itime = 1, param%repetitions
        call mesh_batch_nrm2(sys%gr%mesh, yy, tmp)
      end do
      do itime = 1, xx%nst
        write(message(1),'(a,i1,3x,e13.6)') "Nrm2 norm state  ", itime, tmp(itime)
        call messages_info(1)
      end do

      SAFE_DEALLOCATE_A(tmp)

      call xx%end()
      call yy%end()
    end if

    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_MATRIX) /= 0) then
    
      message(1) = 'Info: Testing dotp_matrix'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if(states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddot(nst, nst))
        call dmesh_batch_dotp_matrix(sys%gr%mesh, xx, yy, ddot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,e13.6)') 'Dotp_matrix states', ist, jst, ddot(ist,jst)
          end do
          call messages_info(nst)
        end do
        SAFE_DEALLOCATE_A(ddot)
      else
        SAFE_ALLOCATE(zdot(nst, nst))
        call zmesh_batch_dotp_matrix(sys%gr%mesh, xx, yy, zdot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,2e14.6)') 'Dotp_matrix states', ist, jst, zdot(ist,jst)
          end do
          call messages_info(nst)
        end do
        SAFE_DEALLOCATE_A(zdot)
      end if
  
      call xx%end()
      call yy%end()    
    end if
  
    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_VECTOR) /= 0) then
    
      message(1) = 'Info: Testing dotp_vector'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)
      call sys%st%group%psib(1, 1)%copy_to(yy, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if(states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddotv(nst))
        call dmesh_batch_dotp_vector(sys%gr%mesh, xx, yy, ddotv)

        do ist = 1, nst
          write(message(ist), '(a,i3,3x,e13.6)') 'Dotp_vector state', ist, ddotv(ist)
        end do
        call messages_info(nst)
        SAFE_DEALLOCATE_A(ddotv)
      else
        SAFE_ALLOCATE(zdotv(nst))
        call zmesh_batch_dotp_vector(sys%gr%mesh, xx, yy, zdotv)
        do ist = 1, nst
          write(message(ist), '(a,i3,3x,2e14.6)') 'Dotp_vector state', ist, zdotv(ist)
        end do
        call messages_info(nst)
        SAFE_DEALLOCATE_A(zdotv)
      end if
  
      call xx%end()
      call yy%end()    
    end if

    if(bitand(ops, OPTION__TESTBATCHOPS__OPS_DOTP_SELF) /= 0) then
    
      message(1) = 'Info: Testing dotp_self'
      call messages_info(1)

      call sys%st%group%psib(1, 1)%copy_to(xx, copy_data = .true.)

      nst = sys%st%group%psib(1, 1)%nst

      if(states_are_real(sys%st)) then
        SAFE_ALLOCATE(ddot(nst, nst))
        call dmesh_batch_dotp_self(sys%gr%mesh, xx, ddot)

        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,e13.6)') 'Dotp_self states', ist, jst, ddot(ist,jst)
          end do
          call messages_info(nst*nst)
        end do
        SAFE_DEALLOCATE_A(ddot)
      else
        SAFE_ALLOCATE(zdot(nst, nst))
        call zmesh_batch_dotp_self(sys%gr%mesh, xx, zdot)
        do ist = 1, nst
          do jst = 1, nst
            write(message(jst), '(a,2i3,3x,2e14.6)') 'Dotp_self states', ist, jst, zdot(ist,jst)
          end do
          call messages_info(nst*nst)
        end do
        SAFE_DEALLOCATE_A(zdot)
      end if
  
      call xx%end()
    end if

    call states_elec_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_batch_ops)
  end subroutine test_batch_ops


! ---------------------------------------------------------
  subroutine test_derivatives(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_derivatives)

    sys => electrons_t(namespace)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call dderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call zderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys
    integer :: itime

    PUSH_SUB(test_orthogonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_par_set_scalapack_compat()

    sys => electrons_t(namespace)

    message(1) = 'Info: Testing orthogonalization.'
    message(2) = ''
    call messages_info(2)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      message(1) = 'Info: Real wave-functions.'
      call messages_info(1)
      do itime = 1, param%repetitions
        call dstates_elec_calc_orth_test(sys%st, sys%namespace, sys%gr%mesh, sys%gr%sb)
      end do
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      message(1) = 'Info: Complex wave-functions.'
      call messages_info(1)
      do itime = 1, param%repetitions
        call zstates_elec_calc_orth_test(sys%st, sys%namespace, sys%gr%mesh, sys%gr%sb)
      end do
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

  ! ---------------------------------------------------------

  subroutine test_interpolation(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_interpolation)

    sys => electrons_t(namespace)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call messages_write('Info: Testing real interpolation routines')
      call messages_new_line()
      call messages_new_line()
      call messages_info()

      call dmesh_interpolation_test(sys%gr%mesh)
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call messages_new_line()
      call messages_write('Info: Testing complex interpolation routines')
      call messages_new_line()
      call messages_new_line()
      call messages_info()

      call zmesh_interpolation_test(sys%gr%mesh)
    end if

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_interpolation)
  end subroutine test_interpolation


  ! ---------------------------------------------------------

  subroutine test_ion_interaction(namespace)
    type(namespace_t),        intent(in) :: namespace

    type(electrons_t), pointer :: sys

    PUSH_SUB(test_ion_interaction)

    sys => electrons_t(namespace)

    call ion_interaction_test(sys%geo, sys%namespace, sys%gr%sb, sys%mc)

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(test_ion_interaction)
  end subroutine test_ion_interaction

  ! ---------------------------------------------------------

  subroutine test_prints_info_batch(st, gr, psib, string)
    type(states_elec_t), intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    class(batch_t),      intent(inout) :: psib
    character(*), optional,  intent(in)    :: string      

    integer :: itime
    CMPLX, allocatable :: zpsi(:, :)
    FLOAT, allocatable :: dpsi(:, :)

    character(80)      :: string_

    string_ = optional_default(string, "")

    PUSH_SUB(test_prints_info_batch)

    if(states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(gr%mesh%np, st%d%dim))
    else
      SAFE_ALLOCATE(zpsi(gr%mesh%np, st%d%dim))
    end if

    do itime = 1, psib%nst
      if(states_are_real(st)) then
        call batch_get_state(psib, itime, gr%mesh%np, dpsi)
        write(message(1),'(a,i1,3x,e13.6)') "Norm state "//trim(string_)//" ", itime, dmf_nrm2(gr%mesh, st%d%dim, dpsi)
      else
        call batch_get_state(psib, itime, gr%mesh%np, zpsi)
        write(message(1),'(a,i1,3x,e13.6)') "Norm state "//trim(string_)//" ", itime, zmf_nrm2(gr%mesh, st%d%dim, zpsi)
      end if
      call messages_info(1)
    end do

    if(states_are_real(st)) then
      SAFE_DEALLOCATE_A(dpsi)
    else
      SAFE_DEALLOCATE_A(zpsi)
    end if

    POP_SUB(test_prints_info_batch)

  end subroutine test_prints_info_batch


  ! ---------------------------------------------------------
  subroutine test_clock()

    type(clock_t) :: test_clock_a, test_clock_b

    PUSH_SUB(test_clock)

    test_clock_a = clock_t('test_clock_a', CNST(2.0), CNST(1.0), 100)
    test_clock_b = clock_t('test_clock_b', CNST(1.0), CNST(1.0))
    call test_clock_a%print()
    call test_clock_b%print()

    call test_clock_a%set_time(test_clock_b)
    call test_clock_a%print()
    call test_clock_a%increment()
    call test_clock_a%print()
    call test_clock_a%decrement()
    call test_clock_a%print()
    call test_clock_a%increment()
    call test_clock_a%print()
    call test_clock_a%reset()
    call test_clock_a%print()
    call test_clock_a%increment(3)
    call test_clock_a%print()
    call test_clock_a%decrement(2)
    call test_clock_a%print()
    message(1) = test_clock_a%print_str()
    call messages_info(1)

    write(message(1),'(A,x,I10.10)') &
	'clock_get_tick', test_clock_a%get_tick()
    write(message(2),'(A,x,F15.10)') &
	'clock_get_sim_time', test_clock_a%get_sim_time()
    write(message(3),'(A,x,I1)')     &
	'clock_is_earlier', abs(transfer(test_clock_a .lt. test_clock_b, 0))
    write(message(4),'(A,x,I1)')     &
	'clock_is_equal_or_earlier', abs(transfer(test_clock_a .le. test_clock_b, 0))
    write(message(5),'(A,x,I1)')     &
	'clock_is_later', abs(transfer(test_clock_a .gt. test_clock_b, 0))
    write(message(6),'(A,x,I1)')     &
	'clock_is_equal_or_later', abs(transfer(test_clock_a .ge. test_clock_b, 0))
    write(message(7),'(A,x,I1)')     &
	'clock_is_equal', abs(transfer(test_clock_a .eq. test_clock_b, 0))
    write(message(8),'(A,x,I1)')     &
	'clock_is_later_with_step', abs(transfer(test_clock_a%is_later_with_step(test_clock_b), 0))
    call messages_info(8)


    POP_SUB(test_clock)
  end subroutine test_clock

end module test_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
