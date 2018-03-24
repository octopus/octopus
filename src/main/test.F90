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

program oct_test
  use global_oct_m
  use batch_oct_m
  use base_hamiltonian_oct_m
  use calc_mode_par_oct_m
  use command_line_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use fft_oct_m
  use io_oct_m
  use ion_interaction_oct_m
  use mesh_interpolation_oct_m
  use mesh_function_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use restart_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use system_oct_m
  use test_parameters_oct_m
  use types_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use XC_F90(lib_m)

  implicit none

  character(len=256) :: config_str
  integer :: test_type
  integer :: test_mode
  type(test_parameters_t) :: test_param
  integer :: ierr
 
  integer, parameter ::              &
    XC_FLAGS_NONE = 0

  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call global_init()
  call calc_mode_par_init()
  call messages_init()

  call messages_obsolete_variable('WhichTest', 'TestMode')

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
  !%End
  call parse_variable('TestMode', OPTION__TESTMODE__HARTREE, test_mode)

  call messages_obsolete_variable('TestDerivatives', 'TestType')
  call messages_obsolete_variable('TestOrthogonalization', 'TestType')
  
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
  call parse_variable('TestType', OPTION__TESTTYPE__ALL, test_type)
  if(test_type < 1 .or. test_type > 5) then
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
  call parse_variable('TestRepetitions', 1, test_param%repetitions)

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
  call parse_variable('TestMinBlockSize', 1, test_param%min_blocksize)

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
  call parse_variable('TestMaxBlockSize', 128, test_param%max_blocksize)
  
  call io_init()
  call profiling_init()

  call print_header()

  call messages_print_stress(stdout, "Test mode")
  call messages_print_var_option(stdout, "TestMode", test_mode)
  call messages_print_var_option(stdout, "TestType", test_type)
  call messages_print_var_value(stdout, "TestRepetitions", test_param%repetitions)
  call messages_print_var_value(stdout, "TestMinBlockSize", test_param%min_blocksize)
  call messages_print_var_value(stdout, "TestMaxBlockSize", test_param%max_blocksize)
  call messages_print_stress(stdout)

  call restart_module_init()
  call fft_all_init()
  call unit_system_init()

  select case(test_mode)
  case(OPTION__TESTMODE__HARTREE)
    call test_hartree()
  case(OPTION__TESTMODE__DERIVATIVES)
    call test_derivatives()
  case(OPTION__TESTMODE__ORTHOGONALIZATION)
    call test_orthogonalization()
  case(OPTION__TESTMODE__INTERPOLATION)
    call test_interpolation()
  case(OPTION__TESTMODE__ION_INTERACTION)
    call test_ion_interaction() 
  case(OPTION__TESTMODE__PROJECTOR)
    call test_projector()
  end select

  call fft_all_end()
  call profiling_output()
  call profiling_end()
  call io_end()
  call print_date("Calculation ended on ")
  call messages_end()
  call calc_mode_par_end()
  call global_end()

  contains

! ---------------------------------------------------------
  subroutine test_hartree
    type(system_t) :: sys

    PUSH_SUB(test_hartree)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call system_init(sys)
    call poisson_test(sys%gr%mesh, test_param)
    call system_end(sys)

    POP_SUB(test_hartree)
  end subroutine test_hartree

 ! ---------------------------------------------------------
  subroutine test_projector
    type(system_t)      :: sys
    type(epot_t) :: ep
    type(batch_t), pointer :: epsib
    integer :: itime
    type(base_hamiltonian_t), pointer :: subsys_hm

    PUSH_SUB(test_projector)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)

    call messages_write('Info: Testing the nonlocal part of the pseudopotential with SOC')
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    call system_init(sys)

    call states_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX)
    call states_generate_random(sys%st, sys%gr%mesh)
  
    !Initialize external potential
    nullify(subsys_hm)
    call epot_init(ep, sys%gr, sys%geo, SPINORS, 1, .false., subsys_hm, XC_FAMILY_NONE)
    call epot_generate(ep, sys%gr, sys%geo, sys%st, .false.)
   
    !Initialize external potential
    SAFE_ALLOCATE(epsib)
    call batch_copy(sys%st%group%psib(1, 1), epsib)

    do itime = 1, test_param%repetitions
      call zproject_psi_batch(sys%gr%mesh, ep%proj, ep%natoms, 2, &
                               sys%st%group%psib(1, 1), epsib, 1)
    end do
    do itime = 1, epsib%nst
      call messages_write(zmf_nrm2(sys%gr%mesh, 2, epsib%states(itime)%zpsi))
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
  subroutine test_derivatives()
    type(system_t) :: sys

    PUSH_SUB(test_derivatives)

    call system_init(sys)

    message(1) = 'Info: Testing the finite-differences derivatives.'
    message(2) = ''
    call messages_info(2)

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__REAL) then
      call dderivatives_test(sys%gr%der, test_param)
    end if

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__COMPLEX) then
      call zderivatives_test(sys%gr%der, test_param)
    end if

    if(test_type == OPTION__TESTTYPE__REAL_SINGLE) then
      call sderivatives_test(sys%gr%der, test_param)
    end if
   
    if(test_type == OPTION__TESTTYPE__COMPLEX_SINGLE) then
      call cderivatives_test(sys%gr%der, test_param)
    end if

    call system_end(sys)

    POP_SUB(test_derivatives)
  end subroutine test_derivatives

  ! ---------------------------------------------------------

  subroutine test_orthogonalization()
    type(system_t) :: sys

    PUSH_SUB(test_orthogonalization)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call calc_mode_par_set_scalapack_compat()

    call system_init(sys)

    message(1) = 'Info: Testing orthogonalization.'
    message(2) = ''
    call messages_info(2)

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__REAL) then
      message(1) = 'Info: Real wave-functions.'
      call messages_info(1)
      call dstates_calc_orth_test(sys%st, sys%gr%mesh)
    end if

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__COMPLEX) then
      message(1) = 'Info: Complex wave-functions.'
      call messages_info(1)
      call zstates_calc_orth_test(sys%st, sys%gr%mesh)
    end if

    call system_end(sys)

    POP_SUB(test_orthogonalization)
  end subroutine test_orthogonalization

  ! ---------------------------------------------------------

  subroutine test_interpolation()
    type(system_t) :: sys

    PUSH_SUB(test_interpolation)

    call system_init(sys)

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__REAL) then
      call messages_write('Info: Testing real interpolation routines')
      call messages_new_line()
      call messages_new_line()
      call messages_info()

      call dmesh_interpolation_test(sys%gr%mesh)
    end if

    if(test_type == OPTION__TESTTYPE__ALL .or. test_type == OPTION__TESTTYPE__COMPLEX) then
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

  subroutine test_ion_interaction()
    type(system_t) :: sys

    PUSH_SUB(test_ion_interaction)

    call system_init(sys)

    call ion_interaction_test(sys%geo, sys%gr%sb)

    call system_end(sys)

    POP_SUB(test_ion_interaction)
  end subroutine test_ion_interaction
  

end program oct_test

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
