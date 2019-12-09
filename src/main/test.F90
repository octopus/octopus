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
  use celestial_body_oct_m
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
  use propagator_verlet_oct_m
  use simul_box_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use subspace_oct_m
  use system_oct_m
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
    !%Calculation of the density.
    !%Option celestial_dynamics 17
    !% Test of celestial dynamics using multisystems
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
    case(OPTION__TESTMODE__EXP_APPLY)
      call test_exponential(param, namespace)
    case(OPTION__TESTMODE__BOUNDARIES)
      call test_boundaries(param, namespace)
    case(OPTION__TESTMODE__SUBSPACE_DIAG)
      call test_subspace_diagonalization(param, namespace)
    case(OPTION__TESTMODE__BATCH_OPS)
      call test_batch_ops(param, namespace)
    case(OPTION__TESTMODE__CELESTIAL_DYNAMICS)
      call test_celestial_dynamics(param)
    end select

    POP_SUB(test_run)
  end subroutine test_run

  ! ---------------------------------------------------------
  subroutine test_hartree(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call system_init(sys, namespace)
    call poisson_test(sys%hm%psolver, sys%gr%mesh, namespace, param%repetitions)
    call system_end(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree

 ! ---------------------------------------------------------
  subroutine test_projector(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    type(wfs_elec_t), pointer :: epsib
    integer :: itime
    CMPLX, allocatable :: psi(:, :)

    PUSH_SUB(test_projector)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the nonlocal part of the pseudopotential with SOC')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

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
    call system_end(sys)

    POP_SUB(test_projector)
  end subroutine test_projector

  ! ---------------------------------------------------------
  subroutine test_dft_u(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
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

    call system_init(sys, namespace)

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
    call system_end(sys)

    POP_SUB(test_dft_u)
  end subroutine test_dft_u

  ! ---------------------------------------------------------
  subroutine test_hamiltonian(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    type(wfs_elec_t), pointer :: hpsib
    integer :: itime, terms
    type(simul_box_t) :: sb

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

    call system_init(sys, namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    !Initialize external potential
    call simul_box_init(sb, sys%namespace, sys%geo, sys%space)
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
    call simul_box_end(sb)
    call states_elec_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_hamiltonian)
  end subroutine test_hamiltonian


  ! ---------------------------------------------------------
  subroutine test_density_calc(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing density calculation')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call density_calc(sys%st, sys%gr, sys%st%rho)
    end do

    write(message(1),'(a,3x, f12.6)') "Norm density  ", dmf_nrm2(sys%gr%mesh, sys%st%rho(:,1))
    call messages_info(1)

    call states_elec_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_density_calc)
  end subroutine test_density_calc


  ! ---------------------------------------------------------
  subroutine test_boundaries(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

    call states_elec_allocate_wfns(sys%st, sys%gr%mesh)
    call states_elec_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call sys%st%pack()

    do itime = 1, param%repetitions
      call boundaries_set(sys%gr%der%boundaries, sys%st%group%psib(1, 1))
    end do

    call test_prints_info_batch(sys%st, sys%gr, sys%st%group%psib(1, 1))

    call states_elec_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_density_calc)
  end subroutine test_boundaries


   ! ---------------------------------------------------------
  subroutine test_exponential(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    type(exponential_t) :: te
    integer :: itime

    PUSH_SUB(test_exponential)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing exponential')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

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
    call system_end(sys)

    POP_SUB(test_exponential)
  end subroutine test_exponential


  ! ---------------------------------------------------------
  subroutine test_subspace_diagonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    integer :: itime
    type(subspace_t) :: sdiag

    PUSH_SUB(test_subspace_diagonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing boundary conditions')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

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
    call system_end(sys)

    POP_SUB(test_subspace_diagonalization)
  end subroutine test_subspace_diagonalization


  ! ---------------------------------------------------------
  subroutine test_batch_ops(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    integer :: itime, ops
    type(wfs_elec_t) :: xx, yy
    FLOAT, allocatable :: tmp(:)

    PUSH_SUB(test_density_calc)

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
    !%End
    ops = OPTION__TESTBATCHOPS__OPS_AXPY &
        + OPTION__TESTBATCHOPS__OPS_SCAL &
        + OPTION__TESTBATCHOPS__OPS_NRM2
    call parse_variable(namespace, 'TestBatchOps', ops, ops)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing batch operations')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, namespace)

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
      call test_prints_info_batch(sys%st, sys%gr, yy)

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
      call test_prints_info_batch(sys%st, sys%gr, yy)

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

    call states_elec_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_density_calc)
  end subroutine test_batch_ops


! ---------------------------------------------------------
  subroutine test_derivatives(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys

    PUSH_SUB(test_derivatives)

    call system_init(sys, namespace)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call dderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call zderivatives_test(sys%gr%der, sys%namespace, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    call system_end(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys
    integer :: itime

    PUSH_SUB(test_orthogonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_par_set_scalapack_compat()

    call system_init(sys, namespace)

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

    call system_end(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

  ! ---------------------------------------------------------

  subroutine test_interpolation(param, namespace)
    type(test_parameters_t), intent(in) :: param
    type(namespace_t),       intent(in) :: namespace

    type(system_t) :: sys

    PUSH_SUB(test_interpolation)

    call system_init(sys, namespace)

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

    call system_end(sys)

    POP_SUB(test_interpolation)
  end subroutine test_interpolation


  ! ---------------------------------------------------------

  subroutine test_ion_interaction(namespace)
    type(namespace_t),        intent(in) :: namespace

    type(system_t) :: sys

    PUSH_SUB(test_ion_interaction)

    call system_init(sys, namespace)

    call ion_interaction_test(sys%geo, sys%namespace, sys%gr%sb)

    call system_end(sys)

    POP_SUB(test_ion_interaction)
  end subroutine test_ion_interaction

  ! ---------------------------------------------------------

  subroutine test_prints_info_batch(st, gr, psib)
    type(states_elec_t), intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    class(batch_t),      intent(inout) :: psib

    integer :: itime
    CMPLX, allocatable :: zpsi(:, :)
    FLOAT, allocatable :: dpsi(:, :)

    PUSH_SUB(test_prints_info_batch)

    if(states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(gr%mesh%np, st%d%dim))
    else
      SAFE_ALLOCATE(zpsi(gr%mesh%np, st%d%dim))
    end if

    do itime = 1, psib%nst
      if(states_are_real(st)) then
        call batch_get_state(psib, itime, gr%mesh%np, dpsi)
        write(message(1),'(a,i1,3x,e13.6)') "Norm state  ", itime, dmf_nrm2(gr%mesh, st%d%dim, dpsi)
      else
        call batch_get_state(psib, itime, gr%mesh%np, zpsi)
        write(message(1),'(a,i1,3x,e13.6)') "Norm state  ", itime, zmf_nrm2(gr%mesh, st%d%dim, zpsi)
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

  subroutine test_celestial_dynamics(param)
    type(test_parameters_t), intent(in) :: param

    type(namespace_t) :: global_namespace, earth_namespace, moon_namespace, sun_namespace
    type(celestial_body_t) :: sun, earth, moon
    type(propagator_verlet_t) :: prop_sun, prop_earth, prop_moon
    integer :: it, Nstep, internal_loop, iunit
    logical :: all_done
    FLOAT :: dt

    PUSH_SUB(test_celestial_dynamics)

    call messages_write('Info: Testing celestial dynamics using multisystems')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    global_namespace = namespace_t("")
    earth_namespace = namespace_t("Earth")
    moon_namespace = namespace_t("Moon")
    sun_namespace = namespace_t("Sun")

    !Initialize subsystems
    call sun%init(sun_namespace)
    call earth%init(earth_namespace)
    call moon%init(moon_namespace)

    !Define interactions manually
    call sun%add_interaction_partner(earth)
    call sun%add_interaction_partner(moon)
    call earth%add_interaction_partner(moon)
    call earth%add_interaction_partner(sun)
    call moon%add_interaction_partner(earth)
    call moon%add_interaction_partner(sun)

    !Allocate receivers
    call sun%allocate_receiv_structure()
    call earth%allocate_receiv_structure()
    call moon%allocate_receiv_structure()

    !Creates Verlet propagators
    call parse_variable(global_namespace, 'TDTimeStep', CNST(10.0), dt)
    call prop_sun%init(M_ZERO, dt)
    call prop_earth%init(M_ZERO, dt)
    call prop_moon%init(M_ZERO, dt)

    !Associate them to subsystems
    call sun%set_propagator(prop_sun)
    call earth%set_propagator(prop_earth)
    call moon%set_propagator(prop_moon)

    iunit = io_open('celestial_dynamics.dat', global_namespace, action='write')
    write(iunit, '(i5,6(e13.6,1x))') 0, sun%pos(1:2), earth%pos(1:2), moon%pos(1:2)

    call parse_variable(global_namespace, 'TDMaxSteps', 1000, Nstep)
    do it = 1, Nstep

      all_done = .false.
      internal_loop = 1

      call prop_sun%rewind()
      call prop_earth%rewind()
      call prop_moon%rewind()

      do while(.not. all_done .and. internal_loop < 1000)

        call sun%system_dt(prop_sun)
        call earth%system_dt(prop_earth)
        call moon%system_dt(prop_moon)

        !We check the exit condition
        all_done = prop_sun%step_is_done() .and. prop_earth%step_is_done() .and. prop_moon%step_is_done()
        INCR(internal_loop, 1)
      end do

      write(stdout,'(a,i5)') 'Iteraction : ', it
      write(iunit, '(i5,6(1x,e13.6))') it, sun%pos(1:2), earth%pos(1:2), moon%pos(1:2)
    end do

    call io_close(iunit)


    POP_SUB(test_celestial_dynamics)
  end subroutine test_celestial_dynamics

end module test_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
