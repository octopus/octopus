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
  use density_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use ion_interaction_oct_m
  use mesh_function_oct_m
  use mesh_interpolation_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use system_oct_m
  use types_oct_m
  use v_ks_oct_m
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
  subroutine test_run(parser)
    type(parser_t),       intent(in)    :: parser
        
    type(test_parameters_t) :: param
    integer :: test_mode
    
    PUSH_SUB(test_run)

    call messages_obsolete_variable(parser, 'WhichTest', 'TestMode')

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
    !%Calculation of the density.
    !%End
    call parse_variable(dummy_parser, 'TestMode', OPTION__TESTMODE__HARTREE, test_mode)

    call messages_obsolete_variable(parser, 'TestDerivatives', 'TestType')
    call messages_obsolete_variable(parser, 'TestOrthogonalization', 'TestType')
  
    !%Variable TestType
    !%Type integer
    !%Default all
    !%Section Utilities::oct-test
    !%Description
    !% Decides on what type of values the test should be performed.
    !%Option real 1
    !% Test for double-precision real functions.
    !%Option complex 2
    !% Test for double-precision complex functions.
    !%Option real_single 4
    !% Test for single-precision real functions. (Only implemented for derivatives.)
    !%Option complex_single 5
    !% Test for single-precision complex functions. (Only implemented for derivatives.)
    !%Option all 3
    !% Tests for double-precision real and complex functions.
    !%End
    call parse_variable(dummy_parser, 'TestType', OPTION__TESTTYPE__ALL, param%type)
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
    call parse_variable(dummy_parser, 'TestRepetitions', 1, param%repetitions)
  
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
    call parse_variable(dummy_parser, 'TestMinBlockSize', 1, param%min_blocksize)

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
    call parse_variable(dummy_parser, 'TestMaxBlockSize', 128, param%max_blocksize)

    call messages_print_stress(stdout, "Test mode")
    call messages_print_var_option(stdout, "TestMode", test_mode)
    call messages_print_var_option(stdout, "TestType", param%type)
    call messages_print_var_value(stdout, "TestRepetitions", param%repetitions)
    call messages_print_var_value(stdout, "TestMinBlockSize", param%min_blocksize)
    call messages_print_var_value(stdout, "TestMaxBlockSize", param%max_blocksize)
    call messages_print_stress(stdout)
  
    select case(test_mode)
    case(OPTION__TESTMODE__HARTREE)
      call test_hartree(param, parser)
    case(OPTION__TESTMODE__DERIVATIVES)
      call test_derivatives(param, parser)
    case(OPTION__TESTMODE__ORTHOGONALIZATION)
      call test_orthogonalization(param, parser)
    case(OPTION__TESTMODE__INTERPOLATION)
      call test_interpolation(param, parser)
    case(OPTION__TESTMODE__ION_INTERACTION)
      call test_ion_interaction(parser)
    case(OPTION__TESTMODE__PROJECTOR)
      call test_projector(param, parser)
    case(OPTION__TESTMODE__DFT_U)
      call test_dft_u(param, parser)
    case(OPTION__TESTMODE__HAMILTONIAN_APPLY)
      call test_hamiltonian(param, parser)
    case(OPTION__TESTMODE__DENSITY_CALC)
      call test_density_calc(param, parser)
    end select
  
    POP_SUB(test_run)
  end subroutine test_run

  ! ---------------------------------------------------------
  subroutine test_hartree(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser

    type(system_t) :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call system_init(sys, parser)
    call poisson_test(sys%gr%mesh, param%repetitions)
    call system_end(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree

 ! ---------------------------------------------------------
  subroutine test_projector(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys
    type(epot_t) :: ep
    type(batch_t), pointer :: epsib
    integer :: itime

    PUSH_SUB(test_projector)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the nonlocal part of the pseudopotential with SOC')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, parser)

    call states_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX)
    call states_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
  
    !Initialize external potential
    call epot_init(ep, sys%parser, sys%gr, sys%geo, SPINORS, 1, XC_FAMILY_NONE)
    call epot_generate(ep, sys%parser, sys%gr, sys%geo, sys%st)
   
    !Initialize external potential
    SAFE_ALLOCATE(epsib)
    call batch_copy(sys%st%group%psib(1, 1), epsib)

    call batch_set_zero(epsib)
    
    do itime = 1, param%repetitions
      call zproject_psi_batch(sys%gr%mesh, ep%proj, ep%natoms, 2, sys%st%group%psib(1, 1), epsib, 1)
    end do
    
    do itime = 1, epsib%nst
      write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, zmf_nrm2(sys%gr%mesh, 2, epsib%states(itime)%zpsi)
      call messages_info(1)
    end do

    call batch_end(epsib)
    SAFE_DEALLOCATE_P(epsib)
    call epot_end(ep)
    call states_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_projector)
  end subroutine test_projector

  ! ---------------------------------------------------------
  subroutine test_dft_u(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys
    type(batch_t), pointer :: epsib
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

    call system_init(sys, parser)

    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call states_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call states_pack(sys%st)

    SAFE_ALLOCATE(epsib)
    call batch_copy(sys%st%group%psib(1, 1), epsib, copy_data = .true.)

    !Initialize the orbital basis
    call orbitalbasis_init(basis)
    if (states_are_real(sys%st)) then
      call dorbitalbasis_build(basis, sys%geo, sys%gr%mesh, sys%st%d%kpt, sys%st%d%dim, &
                                .false., .false.)
      SAFE_ALLOCATE(dweight(1:basis%orbsets(1)%sphere%np,1:epsib%nst_linear))
      SAFE_ALLOCATE(ddot(1:sys%st%d%dim,1:basis%orbsets(1)%norbs, 1:epsib%nst))
    else
      call zorbitalbasis_build(basis, sys%geo, sys%gr%mesh, sys%st%d%kpt, sys%st%d%dim, &
                                .false., .false.)
      call orbitalset_update_phase(basis%orbsets(1), sys%gr%sb, sys%st%d%kpt, (sys%st%d%ispin==SPIN_POLARIZED))
      SAFE_ALLOCATE(zweight(1:basis%orbsets(1)%sphere%np,1:epsib%nst_linear))
      SAFE_ALLOCATE(zdot(1:sys%st%d%dim,1:basis%orbsets(1)%norbs, 1:epsib%nst))
    end if

    do itime = 1, param%repetitions
      call batch_set_zero(epsib)
      if(states_are_real(sys%st)) then
        dweight = M_ONE
        ddot = M_ZERO
        call dorbitalset_get_coeff_batch(basis%orbsets(1), 1, sys%st%group%psib(1, 1), 1, .false., &
                                           .false., ddot)
        call dorbitalset_add_to_batch(basis%orbsets(1), 1, epsib, 1, .false., .false., dweight)
      else
        zweight = M_ONE
        zdot = M_ZERO
        call zorbitalset_get_coeff_batch(basis%orbsets(1), sys%st%d%dim, sys%st%group%psib(1, 1), 1, &
                                           .true., .false., zdot)
        call zorbitalset_add_to_batch(basis%orbsets(1), sys%st%d%dim, epsib, 1, .true., .false., zweight)
      end if
    end do

    if(batch_is_packed(epsib)) then
      call batch_unpack(epsib, force = .true.)
    end if

    do itime = 1, epsib%nst
      if(states_are_real(sys%st)) then 
        write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, dmf_nrm2(sys%gr%mesh, sys%st%d%dim, epsib%states(itime)%dpsi)
      else
        write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, zmf_nrm2(sys%gr%mesh, sys%st%d%dim, epsib%states(itime)%zpsi)
      end if
      call messages_info(1)
    end do

    SAFE_DEALLOCATE_A(dweight)
    SAFE_DEALLOCATE_A(zweight)
    SAFE_DEALLOCATE_A(ddot)
    SAFE_DEALLOCATE_A(zdot)

    call batch_end(epsib)
    SAFE_DEALLOCATE_P(epsib)
    call orbitalbasis_end(basis)
    call states_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_dft_u)
  end subroutine test_dft_u

  ! ---------------------------------------------------------
  subroutine test_hamiltonian(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys
    type(batch_t), pointer :: hpsib
    integer :: itime, terms
    type(hamiltonian_t) :: hm
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
    call parse_variable(dummy_parser, 'TestHamiltonianApply', OPTION__TESTMODE__HARTREE, terms)
    if(terms==0) terms = huge(1)


    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the application of the Hamiltonian')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, parser)

    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call states_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)

    !Initialize external potential
    call simul_box_init(sb, sys%parser, sys%geo, sys%space)
    call hamiltonian_init(hm, sys%parser, sys%gr, sys%geo, sys%st, sys%ks%theory_level, sys%ks%xc_family, &
             family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
    if(sys%st%d%pack_states .and. hamiltonian_apply_packed(hm, sys%gr%mesh)) call states_pack(sys%st)
    call hamiltonian_epot_generate(hm, sys%parser, sys%gr, sys%geo, sys%st)
    call density_calc(sys%st, sys%gr, sys%st%rho)
    call v_ks_calc(sys%ks, sys%parser, hm, sys%st, sys%geo)

    call boundaries_set(sys%gr%der%boundaries, sys%st%group%psib(1, 1)) 

    SAFE_ALLOCATE(hpsib)
    call batch_copy(sys%st%group%psib(1, 1), hpsib)

    if(hamiltonian_apply_packed(hm, sys%gr%der%mesh)) then
      call batch_pack(sys%st%group%psib(1, 1))
      call batch_pack(hpsib, copy = .false.)
    end if

    do itime = 1, param%repetitions
      if(states_are_real(sys%st)) then
        call dhamiltonian_apply_batch(hm, sys%gr%der, sys%st%group%psib(1, 1), hpsib, 1, terms = terms, set_bc = .false.)
      else
        call zhamiltonian_apply_batch(hm, sys%gr%der, sys%st%group%psib(1, 1), hpsib, 1, terms = terms, set_bc = .false.)
      end if
    end do

    if(batch_is_packed(hpsib)) then
      call batch_unpack(hpsib, force = .true.)
    end if
    
    do itime = 1, hpsib%nst
      if(states_are_real(sys%st)) then 
        write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, dmf_nrm2(sys%gr%mesh, sys%st%d%dim, hpsib%states(itime)%dpsi)
      else
        write(message(1),'(a,i1,3x, f12.6)') "Norm state  ", itime, zmf_nrm2(sys%gr%mesh, sys%st%d%dim, hpsib%states(itime)%zpsi)
      end if
      call messages_info(1)
    end do

    call batch_end(hpsib, copy = .false.)
    SAFE_DEALLOCATE_P(hpsib)
    call hamiltonian_end(hm)
    call simul_box_end(sb)
    call states_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_hamiltonian)
  end subroutine test_hamiltonian


  ! ---------------------------------------------------------
  subroutine test_density_calc(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys
    integer :: itime

    PUSH_SUB(test_density_calc)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing density calculation')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys, parser)

    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call states_generate_random(sys%st, sys%gr%mesh, sys%gr%sb)
    if(sys%st%d%pack_states) call states_pack(sys%st)

    do itime = 1, param%repetitions
      call density_calc(sys%st, sys%gr, sys%st%rho)
    end do

    write(message(1),'(a,3x, f12.6)') "Norm density  ", dmf_nrm2(sys%gr%mesh, sys%st%rho(:,1))
    call messages_info(1)

    call states_deallocate_wfns(sys%st)
    call system_end(sys)

    POP_SUB(test_density_calc)
  end subroutine test_density_calc



! ---------------------------------------------------------
  subroutine test_derivatives(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys

    PUSH_SUB(test_derivatives)

    call system_init(sys, parser)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      call dderivatives_test(sys%gr%der, sys%parser, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      call zderivatives_test(sys%gr%der, sys%parser, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    if(param%type == OPTION__TESTTYPE__REAL_SINGLE) then
      call sderivatives_test(sys%gr%der, sys%parser, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if
   
    if(param%type == OPTION__TESTTYPE__COMPLEX_SINGLE) then
      call cderivatives_test(sys%gr%der, sys%parser, param%repetitions, param%min_blocksize, param%max_blocksize)
    end if

    call system_end(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys

    PUSH_SUB(test_orthogonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_par_set_scalapack_compat()

    call system_init(sys, parser)

    message(1) = 'Info: Testing orthogonalization.'
    message(2) = ''
    call messages_info(2)

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__REAL) then
      message(1) = 'Info: Real wave-functions.'
      call messages_info(1)
      call dstates_calc_orth_test(sys%st, sys%gr%mesh, sys%gr%sb)
    end if

    if(param%type == OPTION__TESTTYPE__ALL .or. param%type == OPTION__TESTTYPE__COMPLEX) then
      message(1) = 'Info: Complex wave-functions.'
      call messages_info(1)
      call zstates_calc_orth_test(sys%st, sys%gr%mesh, sys%gr%sb)
    end if

    call system_end(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

  ! ---------------------------------------------------------

  subroutine test_interpolation(param, parser)
    type(test_parameters_t), intent(in) :: param
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys

    PUSH_SUB(test_interpolation)

    call system_init(sys, parser)

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

  subroutine test_ion_interaction(parser)
    type(parser_t),          intent(in) :: parser
    
    type(system_t) :: sys

    PUSH_SUB(test_ion_interaction)

    call system_init(sys, parser)

    call ion_interaction_test(sys%geo, sys%gr%sb)

    call system_end(sys)

    POP_SUB(test_ion_interaction)
  end subroutine test_ion_interaction
  

end module test_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
