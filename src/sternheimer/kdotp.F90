!!! Copyright (C) 2008-2010 David Strubbe
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

module kdotp_oct_m
  use global_oct_m
  use grid_oct_m
  use output_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use ions_oct_m
  use kdotp_calc_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use linear_response_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use sternheimer_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m

  implicit none

  private

  public :: &
       kdotp_lr_run,       &
       int2str      

  type kdotp_t
    private
    type(pert_t) :: perturbation
    type(pert_t) :: perturbation2

    FLOAT, allocatable :: eff_mass_inv(:,:,:,:)  !< inverse effective-mass tensor
                                                 !! (idir1, idir2, ist, ik)
    FLOAT, allocatable :: velocity(:,:,:) !< group velocity vector (idir, ist, ik)

    type(lr_t), allocatable :: lr(:,:) !< linear response for (periodic_dim,1)
                                       !! second index is dummy; should only be 1
                                       !! for compatibility with em_resp routines

    type(lr_t), allocatable :: lr2(:,:,:) !< second-order response for 
                                          !! (periodic_dim,periodic_dim,1)

    logical :: ok                   !< is converged?
    integer :: occ_solution_method  !< how to get occupied components of response
    FLOAT   :: degen_thres          !< maximum energy difference to be considered degenerate
    FLOAT   :: eta                  !< imaginary freq. added to Sternheimer eqn.
  end type kdotp_t

contains

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(kdotp_lr_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = kdotp not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call kdotp_lr_run_legacy(system, from_scratch)
    end select

    POP_SUB(kdotp_lr_run)
  end subroutine kdotp_lr_run

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run_legacy(sys, fromScratch)
    type(electrons_t),   intent(inout) :: sys
    logical,             intent(in)    :: fromScratch

    type(kdotp_t)           :: kdotp_vars
    type(sternheimer_t)     :: sh, sh2
    logical                 :: calc_eff_mass, calc_2nd_order, complex_response

    integer              :: idir, idir2, ierr, pdim, ispin
    character(len=100)   :: str_tmp
    FLOAT                :: errornorm
    type(restart_t)      :: restart_load, restart_dump
	
    type(pert_t)            :: pert2  ! for the second direction in second-order kdotp

    PUSH_SUB(kdotp_lr_run_legacy)

    call messages_experimental("k.p perturbation and calculation of effective masses")

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    !TODO: This test belongs to the pert.F90 file, in the case of the velocity operator not been defined 
    ! from the Hamiltonian
    ! In this case, there are other terms missing (MGGA, DFT+U for instance).
    if(sys%hm%theory_level == HARTREE_FOCK) then
      call messages_not_implemented('Commutator of Fock operator')
    end if
    if(sys%hm%theory_level == GENERALIZED_KOHN_SHAM_DFT) then
      call messages_not_implemented('k.p with generalized Kohn-Sham DFT')
    end if


    if (sys%kpoints%use_symmetries) then
      call messages_experimental("KPoints symmetries with CalculationMode = kdotp")
    end if

    pdim = sys%space%periodic_dim

    if(.not. sys%space%is_periodic()) then
       message(1) = "k.p perturbation cannot be used for a finite system."
       call messages_fatal(1)
    end if

    call pert_init(kdotp_vars%perturbation, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%ions)
    SAFE_ALLOCATE(kdotp_vars%lr(1:1, 1:pdim))

    call parse_input()

    if(calc_2nd_order) then
      call pert_init(kdotp_vars%perturbation2, sys%namespace, PERTURBATION_NONE, sys%gr, sys%ions)
      call pert_init(pert2, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%ions)
      call pert_setup_dir(kdotp_vars%perturbation2, 1) ! direction is irrelevant
      SAFE_ALLOCATE(kdotp_vars%lr2(1:1, 1:pdim, 1:pdim))
    end if

    !Read ground-state wavefunctions
    complex_response = (kdotp_vars%eta /= M_ZERO ) .or. states_are_complex(sys%st)
    call restart_init(restart_load, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_look_and_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                                       is_complex = complex_response)
      call restart_end(restart_load)
    else
      message(1) = "A previous gs calculation is required."
      call messages_fatal(1)
    end if

    ! Use of ForceComplex will make this true after states_elec_look_and_load even if it was not before.
    ! Otherwise, this line is a tautology.
    complex_response = states_are_complex(sys%st)

    ! Start restart. Note: we are going to use the same directory to read and write.
    ! Therefore, restart_dump must be initialized first to make sure the directory
    ! exists when we initialize restart_load.
    call restart_init(restart_dump, sys%namespace, RESTART_KDOTP, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)
    ! no problem if this fails
    call restart_init(restart_load, sys%namespace, RESTART_KDOTP, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call messages_info(1)
    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
    
    if(states_are_real(sys%st)) then
      message(1) = 'Info: Using real wavefunctions.'
    else
      message(1) = 'Info: Using complex wavefunctions.'
    end if
    call messages_info(1)

    message(1) = 'Calculating band velocities.'
    call messages_info(1)
 
    SAFE_ALLOCATE(kdotp_vars%velocity(1:pdim, 1:sys%st%nst, 1:sys%st%d%nik))
    kdotp_vars%velocity(:,:,:) = M_ZERO
    if(states_are_complex(sys%st)) then
      call zcalc_band_velocity(sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ions, kdotp_vars%perturbation, &
        kdotp_vars%velocity(:,:,:))
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(KDOTP_DIR, sys%namespace) ! data output
      call kdotp_write_band_velocity(sys%st, pdim, kdotp_vars%velocity(:,:,:), sys%namespace)
    end if
    SAFE_DEALLOCATE_A(kdotp_vars%velocity)

    call sternheimer_obsolete_variables(sys%namespace, 'KdotP_', 'KdotP')
    call sternheimer_init(sh, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, complex_response, &
      set_ham_var = 0, set_occ_response = (kdotp_vars%occ_solution_method == 0), &
      set_last_occ_response = (kdotp_vars%occ_solution_method == 0), occ_response_by_sternheimer = .true.)
    ! ham_var_set = 0 results in HamiltonianVariation = V_ext_only
    if(calc_2nd_order) then
      call sternheimer_init(sh2, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, complex_response, &
        set_ham_var = 0, set_occ_response = .false., set_last_occ_response = .false.)
    end if

    do idir = 1, pdim
      call lr_init(kdotp_vars%lr(1, idir))
      call lr_allocate(kdotp_vars%lr(1, idir), sys%st, sys%gr%mesh)

      if(calc_2nd_order) then
        do idir2 = idir, pdim
          call lr_init(kdotp_vars%lr2(1, idir, idir2))
          call lr_allocate(kdotp_vars%lr2(1, idir, idir2), sys%st, sys%gr%mesh)
        end do
      end if

      ! load wavefunctions
      if(.not. fromScratch) then
        str_tmp = kdotp_wfs_tag(idir)
        call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, 1), ierr)
        if (ierr == 0) call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                                 ierr, lr=kdotp_vars%lr(1, idir))
        call restart_close_dir(restart_load)
          
        if(ierr /= 0) then
          message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'."
          call messages_warning(1)
        end if

        if(calc_2nd_order) then
          do idir2 = idir, pdim
            str_tmp = kdotp_wfs_tag(idir, idir2)
            call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, 1), ierr)
            if (ierr == 0) then
              call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                           ierr, lr=kdotp_vars%lr2(1, idir, idir2))
            end if
            call restart_close_dir(restart_load)
          
            if(ierr /= 0) then
              message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'."
              call messages_warning(1)
            end if
          end do
        end if
      end if
    end do

    call info()
    message(1) = "Info: Calculating k.p linear response of ground-state wavefunctions."
    call messages_info(1)
    kdotp_vars%ok = .true.

    ! solve the Sternheimer equation
    do idir = 1, pdim
      write(message(1), '(3a)') 'Info: Calculating response for the ', index2axis(idir), &
                                '-direction.' 
      call messages_info(1)
      call pert_setup_dir(kdotp_vars%perturbation, idir)

      if(states_are_real(sys%st)) then
        call dsternheimer_solve(sh, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, sys%ions, &
          kdotp_vars%lr(1:1, idir), 1, M_ZERO, kdotp_vars%perturbation, restart_dump, "", kdotp_wfs_tag(idir), &
          have_restart_rho = .false.)
        if (kdotp_vars%occ_solution_method == 1) then
          call dkdotp_add_occ(sys%namespace, sys%gr, sys%st, sys%hm, sys%ions, kdotp_vars%perturbation, &
            kdotp_vars%lr(1, idir), kdotp_vars%degen_thres)
        end if
      else
        call zsternheimer_solve(sh, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, sys%ions, &
          kdotp_vars%lr(1:1, idir), 1, M_zI * kdotp_vars%eta, kdotp_vars%perturbation, restart_dump, "", &
          kdotp_wfs_tag(idir), have_restart_rho = .false.)
        if (kdotp_vars%occ_solution_method == 1) then
          call zkdotp_add_occ(sys%namespace, sys%gr, sys%st, sys%hm, sys%ions, kdotp_vars%perturbation, kdotp_vars%lr(1, idir), &
            kdotp_vars%degen_thres)
        end if
      end if

      kdotp_vars%ok = kdotp_vars%ok .and. sternheimer_has_converged(sh)         

      errornorm = M_ZERO
      if(states_are_real(sys%st)) then 
        call doutput_lr(sys%outp, sys%namespace, sys%space, KDOTP_DIR, sys%st, sys%gr%mesh, kdotp_vars%lr(1, idir), idir, 1, &
          sys%ions, units_out%force)

        do ispin = 1, sys%st%d%nspin
          errornorm = hypot(errornorm, TOFLOAT(dmf_nrm2(sys%gr%mesh, kdotp_vars%lr(1, idir)%ddl_rho(:, ispin))))
        end do
      else
        call zoutput_lr(sys%outp, sys%namespace, sys%space, KDOTP_DIR, sys%st, sys%gr%mesh, kdotp_vars%lr(1, idir), idir, 1, &
          sys%ions, units_out%force)

        do ispin = 1, sys%st%d%nspin
          errornorm = hypot(errornorm, TOFLOAT(zmf_nrm2(sys%gr%mesh, kdotp_vars%lr(1, idir)%zdl_rho(:, ispin))))
        end do
      end if

      write(message(1),'(a,g12.6)') "Norm of relative density variation = ", errornorm / sys%st%qtot
      call messages_info(1)

      if(calc_2nd_order) then
        ! by equality of mixed partial derivatives, kdotp_vars%lr2(idir, idir2) = kdotp_vars%lr2(idir2, idir)
        do idir2 = idir, pdim
          write(message(1), '(3a)') 'Info: Calculating second-order response in the ', index2axis(idir2), &
            '-direction.' 
          call messages_info(1)
		  
          call pert_setup_dir(pert2, idir2)

          if(states_are_real(sys%st)) then
            call dsternheimer_solve_order2(sh, sh, sh2, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, &
              sys%ks%xc, sys%mc, sys%ions, kdotp_vars%lr(1:1, idir), kdotp_vars%lr(1:1, idir2), &
              1, M_ZERO, M_ZERO, kdotp_vars%perturbation, pert2, &
              kdotp_vars%lr2(1:1, idir, idir2), kdotp_vars%perturbation2, restart_dump, "", kdotp_wfs_tag(idir, idir2), &
              have_restart_rho = .false., have_exact_freq = .true.)
          else
            call zsternheimer_solve_order2(sh, sh, sh2, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, &
              sys%ks%xc, sys%mc, sys%ions, kdotp_vars%lr(1:1, idir), kdotp_vars%lr(1:1, idir2), &
              1, M_zI * kdotp_vars%eta, M_zI * kdotp_vars%eta, kdotp_vars%perturbation, pert2, &
              kdotp_vars%lr2(1:1, idir, idir2), kdotp_vars%perturbation2, restart_dump, "", kdotp_wfs_tag(idir, idir2), &
              have_restart_rho = .false., have_exact_freq = .true.)
          end if

        end do
        message(1) = ""
        call messages_info(1)
      end if
    end do ! idir

    ! calculate effective masses
    if (calc_eff_mass) then
      message(1) = "Info: Calculating effective masses."
      call messages_info(1)

      SAFE_ALLOCATE(kdotp_vars%eff_mass_inv(1:pdim, 1:pdim, 1:sys%st%nst, 1:sys%st%d%nik))
      kdotp_vars%eff_mass_inv(:,:,:,:) = M_ZERO


      if(states_are_real(sys%st)) then
        call dcalc_eff_mass_inv(sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ions, kdotp_vars%lr, &
          kdotp_vars%perturbation, kdotp_vars%eff_mass_inv, kdotp_vars%degen_thres)
      else
        call zcalc_eff_mass_inv(sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ions, kdotp_vars%lr, &
          kdotp_vars%perturbation, kdotp_vars%eff_mass_inv, kdotp_vars%degen_thres)
      end if

      call kdotp_write_degeneracies(sys%st, kdotp_vars%degen_thres)
      call kdotp_write_eff_mass(sys%st, sys%kpoints, kdotp_vars, sys%namespace, sys%space%periodic_dim)

      SAFE_DEALLOCATE_A(kdotp_vars%eff_mass_inv)
    end if

    ! clean up some things
    do idir = 1, pdim
      call lr_dealloc(kdotp_vars%lr(1, idir))

      if(calc_2nd_order) then
        do idir2 = idir, pdim
          call lr_dealloc(kdotp_vars%lr2(1, idir, idir2))
        end do
      end if
    end do

    call restart_end(restart_load)
    call restart_end(restart_dump)
    call sternheimer_end(sh)
    call pert_end(kdotp_vars%perturbation)
    SAFE_DEALLOCATE_A(kdotp_vars%lr)

    if(calc_2nd_order) then
      call sternheimer_end(sh2)
      call pert_end(pert2)
      call pert_end(kdotp_vars%perturbation2)
      SAFE_DEALLOCATE_A(kdotp_vars%lr2)
    end if

    call states_elec_deallocate_wfns(sys%st)

    POP_SUB(kdotp_lr_run_legacy)

  contains

    ! ---------------------------------------------------------

    subroutine parse_input()

      PUSH_SUB(kdotp_lr_run_legacy.parse_input)

      !%Variable KdotPOccupiedSolutionMethod
      !%Type integer
      !%Default sternheimer_eqn
      !%Section Linear Response::KdotP
      !%Description
      !% Method of calculating the contribution of the projection of the
      !% linear-response wavefunctions in the occupied subspace.
      !%Option sternheimer_eqn 0
      !% The Sternheimer equation is solved including the occupied subspace,
      !% to get the full linear-response wavefunctions.
      !%Option sum_over_states 1
      !% The Sternheimer equation is solved only in the unoccupied subspace,
      !% and a sum-over-states perturbation-theory expression is used to
      !% evaluate the contributions in the occupied subspace.
      !%End

      call messages_obsolete_variable(sys%namespace, 'KdotP_OccupiedSolutionMethod', 'KdotPOccupiedSolutionMethod')

      call parse_variable(sys%namespace, 'KdotPOccupiedSolutionMethod', 0, kdotp_vars%occ_solution_method)
      if(kdotp_vars%occ_solution_method == 1 .and. .not. smear_is_semiconducting(sys%st%smear)) then
        call messages_not_implemented('KdotPOccupiedSolutionMethod = sum_over_states for non-semiconducting smearing')
      end if

      !%Variable DegeneracyThreshold
      !%Type float
      !%Default 1e-5
      !%Section States
      !%Description
      !% States with energy <math>E_i</math> and <math>E_j</math> will be considered degenerate
      !% if <math> \left| E_i - E_j \right| < </math><tt>DegeneracyThreshold</tt>.
      !%End
      call parse_variable(sys%namespace, 'DegeneracyThreshold', units_from_atomic(units_inp%energy, CNST(1e-5)), kdotp_vars%degen_thres)
      kdotp_vars%degen_thres = units_to_atomic(units_inp%energy, kdotp_vars%degen_thres)

      !%Variable KdotPEta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::KdotP
      !%Description
      !% Imaginary frequency added to Sternheimer equation which may improve convergence.
      !% Not recommended.
      !%End
      call messages_obsolete_variable(sys%namespace, 'KdotP_Eta', 'KdotPEta')
      call parse_variable(sys%namespace, 'KdotPEta', M_ZERO, kdotp_vars%eta)
      kdotp_vars%eta = units_to_atomic(units_inp%energy, kdotp_vars%eta)

      !%Variable KdotPCalculateEffectiveMasses
      !%Type logical
      !%Default true
      !%Section Linear Response::KdotP
      !%Description
      !% If true, uses <tt>kdotp</tt> perturbations of ground-state wavefunctions
      !% to calculate effective masses. It is not correct for degenerate states.
      !%End      
      call messages_obsolete_variable(sys%namespace, 'KdotP_CalculateEffectiveMasses', 'KdotPCalculateEffectiveMasses')
      call parse_variable(sys%namespace, 'KdotPCalculateEffectiveMasses', .true., calc_eff_mass)

      !%Variable KdotPCalcSecondOrder
      !%Type logical
      !%Default false
      !%Section Linear Response::KdotP
      !%Description
      !% If true, calculates second-order response of wavefunctions as well as first-order response.
      !% Note that the second derivative of the Hamiltonian is NOT included in this calculation.
      !% This is needed for a subsequent run in <tt>CalculationMode = em_resp</tt> with <tt>EMHyperpol</tt>.
      !%End      
      call parse_variable(sys%namespace, 'KdotPCalcSecondOrder', .false., calc_2nd_order)

      POP_SUB(kdotp_lr_run_legacy.parse_input)

   end subroutine parse_input

    ! ---------------------------------------------------------
    subroutine info()

      PUSH_SUB(kdotp_lr_run_legacy.info)

      call pert_info(kdotp_vars%perturbation)

      write(message(1),'(a)') 'k.p perturbation theory'
      call messages_print_stress(stdout, trim(message(1)))

      if (kdotp_vars%occ_solution_method == 0) then
        message(1) = 'Occupied solution method: Sternheimer equation.'
      else
        message(1) = 'Occupied solution method: sum over states.'
      end if

      call messages_info(1)

      call messages_print_stress(stdout)
      
      POP_SUB(kdotp_lr_run_legacy.info)

    end subroutine info

  end subroutine kdotp_lr_run_legacy

  ! ---------------------------------------------------------
  subroutine kdotp_write_band_velocity(st, periodic_dim, velocity, namespace)
    type(states_elec_t), intent(inout) :: st
    integer,             intent(in)    :: periodic_dim
    FLOAT,               intent(in)    :: velocity(:,:,:)
    type(namespace_t),   intent(inout) :: namespace

    character(len=80) :: filename, tmp
    integer :: iunit, ik, ist, ik2, ispin, idir

    PUSH_SUB(kdotp_write_band_velocity)

    write(filename, '(a)') KDOTP_DIR//'velocity'
    iunit = io_open(trim(filename), namespace, action='write')
    write(iunit,'(a)') '# Band velocities'

    do ik = 1, st%d%nik
      ispin = st%d%get_spin_index(ik)
      ik2 = st%d%get_kpoint_index(ik)
      tmp = int2str(ik2)

      write(iunit,'(a,i1,a,a)') '# spin = ', ispin, ', k-point = ', trim(tmp)

      write(iunit,'(a)',advance='no') '# state    energy       '
      do idir = 1, periodic_dim
        write(iunit,'(3a)',advance='no') 'vg(', trim(index2axis(idir)), ')       '
      end do
      write(iunit,'(a)')

      write(iunit,'(3a)',advance='no')       '#           [', trim(units_abbrev(units_out%energy)), ']     '
      do idir = 1, periodic_dim
        write(iunit,'(3a)',advance='no') '[', trim(units_abbrev(units_out%velocity)), '] '
      end do
      write(iunit,'(a)')

      do ist = 1, st%nst
        write(iunit,'(i5,f12.5,3f12.5)') ist, units_from_atomic(units_out%energy, st%eigenval(ist, ik)), &
          velocity(1:periodic_dim, ist, ik)
      end do
    end do

    call io_close(iunit)
    POP_SUB(kdotp_write_band_velocity)
  end subroutine kdotp_write_band_velocity

  ! ---------------------------------------------------------
  subroutine kdotp_write_eff_mass(st, kpoints, kdotp_vars, namespace, periodic_dim)
    type(states_elec_t),  intent(inout) :: st
    type(kpoints_t),      intent(in)    :: kpoints
    type(kdotp_t),        intent(inout) :: kdotp_vars
    type(namespace_t),    intent(in)    :: namespace
    integer,              intent(in)    :: periodic_dim

    character(len=80) :: filename, tmp
    integer :: iunit, ik, ist, ik2, ispin

    PUSH_SUB(kdotp_write_eff_mass)

    do ik = 1, st%d%nik
      ispin = st%d%get_spin_index(ik)
      ik2 = st%d%get_kpoint_index(ik)

      tmp = int2str(ik2)
      write(filename, '(3a, i1)') KDOTP_DIR//'kpoint_', trim(tmp), '_', ispin
      iunit = io_open(trim(filename), namespace, action='write')

      write(iunit,'(a, i10)')    '# spin    index = ', ispin
      write(iunit,'(a, i10)')    '# k-point index = ', ik2
      write(iunit,'(a, 99f12.8)') '# k-point coordinates = ', kpoints%get_point(ik2)
      if (.not. kdotp_vars%ok) write(iunit, '(a)') "# WARNING: not converged"      
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Inverse effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        call output_tensor(iunit, kdotp_vars%eff_mass_inv(:, :, ist, ik), periodic_dim, unit_one)
      end do
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        call lalg_inverter(periodic_dim, kdotp_vars%eff_mass_inv(:, :, ist, ik))
        call output_tensor(iunit, kdotp_vars%eff_mass_inv(:, :, ist, ik), periodic_dim, unit_one)
      end do

      call io_close(iunit)
    end do

    POP_SUB(kdotp_write_eff_mass)
  end subroutine kdotp_write_eff_mass

  ! ---------------------------------------------------------
  subroutine kdotp_write_degeneracies(st, threshold)
    type(states_elec_t), intent(inout) :: st
    FLOAT,               intent(in)    :: threshold

    character(len=80) :: tmp
    integer :: ik, ist, ist2, ik2, ispin

    PUSH_SUB(kdotp_write_degeneracies)

    call messages_print_stress(stdout, 'Degenerate subspaces')

    do ik = 1, st%d%nik
      ispin = st%d%get_spin_index(ik)
      ik2 = st%d%get_kpoint_index(ik)

      tmp = int2str(ik2)
      write(message(1), '(3a, i1)') 'k-point ', trim(tmp), ', spin ', ispin 
      call messages_info(1)

      ist = 1
      do while (ist <= st%nst)
      ! test for degeneracies
         write(message(1),'(a)') '===='
         tmp = int2str(ist)
         write(message(2),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
           units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
         call messages_info(2)

         ist2 = ist + 1
         do while (ist2 <= st%nst .and. &
           abs(st%eigenval(min(ist2, st%nst), ik) - st%eigenval(ist, ik)) < threshold)
           tmp = int2str(ist2)
           write(message(1),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
             units_from_atomic(units_out%energy, st%eigenval(ist2, ik)), ' ', units_abbrev(units_out%energy)
           call messages_info(1)
           ist2 = ist2 + 1
         end do

         ist = ist2
      end do

      write(message(1),'()')
      call messages_info(1)
    end do

    message(1) = "Velocities and effective masses are not correct within degenerate subspaces."
    call messages_warning(1)

    POP_SUB(kdotp_write_degeneracies)
  end subroutine kdotp_write_degeneracies

  ! ---------------------------------------------------------
  character(len=12) pure function int2str(ii) result(str)
    integer, intent(in) :: ii
    
    write(str, '(i11)') ii
    str = trim(adjustl(str))

  end function int2str
            
end module kdotp_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
