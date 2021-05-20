!! Copyright (C) 2004-2012 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca), David Strubbe
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

module em_resp_oct_m
  use born_charges_oct_m
  use electrons_oct_m
  use em_resp_calc_oct_m
  use forces_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use output_oct_m
  use io_oct_m
  use ions_oct_m
  use kdotp_oct_m
  use kdotp_calc_oct_m
  use kpoints_oct_m
  use linear_response_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pcm_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use sort_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use sternheimer_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  private

  public :: &
       em_resp_run,             &
       out_hyperpolarizability

  type em_resp_t
    private
    type(pert_t) :: perturbation

    integer :: nsigma !< 1: consider only positive values of the frequency
                      !! 2: consider both positive and negative
    integer :: nfactor!< 1: only one frequency needed
                      !! 2: two frequencies (complex conjugate for magneto-optics)
                      !! 3: three frequencies (for the hyperpolarizabilities)
    integer :: nomega !< number of frequencies to consider

    FLOAT :: eta                     !< small imaginary part to add to the frequency
    FLOAT :: freq_factor(3)
    FLOAT,      allocatable :: omega(:)  !< the frequencies to consider
    type(lr_t), allocatable :: lr(:,:,:) !< linear response for (space%dim, nsigma, nfactor)
    CMPLX,      allocatable :: alpha_k(:, :, :, :)    !< contributions of k-points to 
                                                      !! the linear polarizability
    CMPLX,      allocatable :: alpha_be_k(:, :, :, :) !< contributions of k-points to 
                                                      !! the magneto-optical response
    logical :: calc_hyperpol
    CMPLX   :: alpha(MAX_DIM, MAX_DIM, 3)        !< the linear polarizability
    CMPLX   :: alpha_be(MAX_DIM, MAX_DIM, MAX_DIM) !< the magneto-optical response
    CMPLX   :: alpha0(MAX_DIM, MAX_DIM, 3)        !< the linear polarizability without G = G` = 0 term
                                                  !! of the LRC kernel
    CMPLX   :: alpha_be0(MAX_DIM, MAX_DIM, MAX_DIM) !< the magneto-optical response without G = G` = 0
                                                  !! term of the LRC kernel
    CMPLX   :: beta (MAX_DIM, MAX_DIM, MAX_DIM)  !< first hyperpolarizability

    CMPLX   :: chi_para(MAX_DIM, MAX_DIM)     !< The paramagnetic part of the susceptibility
    CMPLX   :: chi_dia (MAX_DIM, MAX_DIM)     !< The diamagnetic  part of the susceptibility
    CMPLX   :: magn(MAX_DIM)                     !< The orbital magnetization

    logical :: ok(1:3)                           !< whether calculation is converged
    logical :: force_no_kdotp                    !< whether to use kdotp run for periodic system

    logical :: calc_rotatory                     !< whether to calculate rotatory response
    logical :: calc_Born                         !< whether to calculate Born effective charges
    type(Born_charges_t) :: Born_charges(3)      !< one set for each frequency factor
    logical :: occ_response                      !< whether to calculate full response in Sternheimer eqn.
    logical :: wfns_from_scratch                 !< whether to ignore restart LR wfns and initialize to zero
    logical :: calc_magnetooptics                !< whether to calculate magneto-optical response
    logical :: magnetooptics_nohvar              !< whether to consider corrections to exchange-correlation 
                                                 !! and Hartree terms for magnetic perturbations in magneto-optics
    logical :: kpt_output                        !< whether to include in the output contributions of different 
                                                 !! k-points to dielectric constant  
    logical :: lrc_kernel                        !< whether the LRC kernel is used   

  end type em_resp_t

contains

  ! ---------------------------------------------------------
  subroutine em_resp_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(em_resp_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = em_resp not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call em_resp_run_legacy(system, from_scratch)
    end select

    POP_SUB(em_resp_run)
  end subroutine em_resp_run

  ! ---------------------------------------------------------
  subroutine em_resp_run_legacy(sys, fromScratch)
    type(electrons_t), intent(inout) :: sys
    logical,           intent(in)    :: fromScratch

    type(em_resp_t)         :: em_vars
    type(sternheimer_t)     :: sh, sh_kdotp, sh2, sh_kmo, sh_mo
    type(lr_t)              :: kdotp_lr(MAX_DIM, 1)
    type(lr_t), allocatable :: kdotp_em_lr2(:, :, :, :)
    type(lr_t), allocatable :: b_lr(:, :)
    type(lr_t), allocatable :: kb_lr(:, :, :), k2_lr(:, :, :)
    type(lr_t), allocatable :: ke_lr(:, :, :, :)
    type(pert_t)            :: pert_kdotp, pert2_none, pert_b

    integer :: sigma, idir, idir2, ierr, iomega, ifactor
    integer :: ierr_e(3), ierr_e2(3), nfactor_ke
    character(len=100) :: str_tmp
    logical :: complex_response, have_to_calculate, use_kdotp, opp_freq, &
      exact_freq(3), complex_wfs, allocate_rho_em, allocate_rho_mo

    FLOAT :: last_omega, frequency
    FLOAT, allocatable :: dl_eig(:,:,:)
    CMPLX :: frequency_eta, frequency_zero, lrc_coef(MAX_DIM, MAX_DIM)
    type(restart_t) :: gs_restart, kdotp_restart

    PUSH_SUB(em_resp_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%kpoints%use_symmetries) then
      call messages_experimental("em_resp with k-points symmetries")
    end if

    if(sys%kpoints%reduced%npoints /= sys%kpoints%full%npoints) then
      call messages_experimental('em_resp with reduced k-grid')
    end if

    call parse_input()

    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC .and. &
      any(abs(em_vars%omega(1:em_vars%nomega)) > M_EPSILON)) then
      call messages_not_implemented('Dynamical magnetic response')
    end if

    em_vars%lrc_kernel = .false.
    if(abs(sys%ks%xc%kernel_lrc_alpha) > M_EPSILON) em_vars%lrc_kernel = .true.

    if(em_vars%lrc_kernel .and. sys%space%periodic_dim < sys%space%dim) then
      message(1) = 'The use of the LRC kernel for non-periodic dimensions makes no sense.'
      call messages_warning(1)
    end if

    complex_wfs = states_are_complex(sys%st)
    complex_response = (em_vars%eta > M_EPSILON) .or. states_are_complex(sys%st)
    call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_look_and_load(gs_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                    is_complex = complex_response)
      call restart_end(gs_restart)
    else
      message(1) = "Previous gs calculation is required."
      call messages_fatal(1)
    end if

    ! Use of ForceComplex will make this true after states_elec_look_and_load even if it was not before.
    ! Otherwise, this line is a tautology.
    complex_response = states_are_complex(sys%st)

    if (states_are_real(sys%st)) then
      message(1) = 'Info: Using real wavefunctions.'
    else
      message(1) = 'Info: Using complex wavefunctions.'
    end if
    call messages_info(1)

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call messages_info(1)
    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)

    use_kdotp = sys%space%is_periodic() .and. .not. em_vars%force_no_kdotp

    if(use_kdotp .and. .not. smear_is_semiconducting(sys%st%smear)) then
      ! there needs to be a gap.
      message(1) = "em_resp with kdotp can only be used with semiconducting smearing"
      call messages_fatal(1)
    end if

    ! read kdotp wavefunctions if necessary
    if (use_kdotp) then
      message(1) = "Reading kdotp wavefunctions for periodic directions."
      call messages_info(1)

      call restart_init(kdotp_restart, sys%namespace, RESTART_KDOTP, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)
      if(ierr /= 0) then
        message(1) = "Unable to read kdotp wavefunctions."
        message(2) = "Previous kdotp calculation required."
        call messages_fatal(2)
      end if
      
      do idir = 1, sys%space%periodic_dim
        call lr_init(kdotp_lr(idir, 1))
        call lr_allocate(kdotp_lr(idir, 1), sys%st, sys%gr%mesh, allocate_rho = .false.)

        ! load wavefunctions
        str_tmp = kdotp_wfs_tag(idir)
        ! 1 is the sigma index which is used in em_resp
        call restart_open_dir(kdotp_restart, wfs_tag_sigma(str_tmp, 1), ierr)
        if (ierr == 0) then
          call states_elec_load(kdotp_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr, &
            lr=kdotp_lr(idir, 1))
        end if
        call restart_close_dir(kdotp_restart)

        if(ierr /= 0) then
          message(1) = "Could not load kdotp wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'"
          message(2) = "Previous kdotp calculation required."
          call messages_fatal(2)
        end if
      end do

      call restart_end(kdotp_restart)
    end if

    em_vars%nfactor = 1
    if(em_vars%calc_hyperpol) em_vars%nfactor = 3

    ! in effect, nsigma = 1 only if hyperpol not being calculated, and the only frequency is zero
    if(em_vars%calc_hyperpol .or. any(abs(em_vars%omega(1:em_vars%nomega)) > M_EPSILON)) then
      em_vars%nsigma = 2
      ! positive and negative values of the frequency must be considered
    else
      em_vars%nsigma = 1
      ! only considering positive values
    end if

    if(em_vars%calc_hyperpol .and. use_kdotp) then
      call pert_init(pert_kdotp, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%ions)
      call pert_init(pert2_none, sys%namespace, PERTURBATION_NONE,  sys%gr, sys%ions)
      call messages_experimental("Second-order Sternheimer equation")
      call pert_setup_dir(pert2_none, 1)  ! direction is irrelevant
      SAFE_ALLOCATE(kdotp_em_lr2(1:sys%space%periodic_dim, 1:sys%space%dim, 1:em_vars%nsigma, 1:em_vars%nfactor))
      do ifactor = 1, em_vars%nfactor
        do sigma = 1, em_vars%nsigma
          do idir = 1, sys%space%periodic_dim
            do idir2 = 1, sys%space%dim
              call lr_init(kdotp_em_lr2(idir, idir2, sigma, ifactor))
              call lr_allocate(kdotp_em_lr2(idir, idir2, sigma, ifactor), sys%st, sys%gr%mesh, allocate_rho = .false.)
            end do
          end do
        end do
      end do
      call sternheimer_init(sh2, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
        complex_response, set_ham_var = 0, set_last_occ_response = .false.)
      call sternheimer_init(sh_kdotp, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
        complex_response, set_ham_var = 0, set_last_occ_response = .true.)
      em_vars%occ_response = .true.
      SAFE_ALLOCATE(dl_eig(1:sys%st%nst, 1:sys%st%d%nik, 1:sys%space%periodic_dim))
    end if

    ! Hyperpolarizability requires full corrections to wavefunctions (with projections on occupied states).
    ! Magnetooptics is implemented only for projections of corrections to wavefunctions on unoccupied states.
    if(em_vars%calc_magnetooptics) then 
      if(em_vars%calc_hyperpol .and. use_kdotp) then
        message(1) = "Hyperpolarizability and magnetooptics with kdotp are not compatible."
        message(2) = "Only calculation of hyperpolarizability will be performed."
        call messages_warning(2)
        em_vars%calc_magnetooptics = .false.
      else
        em_vars%nfactor = 2    
        em_vars%freq_factor(1) = M_ONE
        em_vars%freq_factor(2) = -M_ONE
      end if
    end if

    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
      em_vars%nsigma = 1
      if(use_kdotp) call messages_experimental("Magnetic perturbation for periodic systems")
    end if  

    if(em_vars%calc_magnetooptics .or. &
      (pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC)) then
      frequency_zero = M_ZERO
      em_vars%occ_response = .false.
   
      if(use_kdotp) then
        call pert_init(pert2_none, sys%namespace, PERTURBATION_NONE,  sys%gr, sys%ions)
        call pert_setup_dir(pert2_none, 1) 

        SAFE_ALLOCATE(k2_lr(1:sys%space%dim, 1:sys%space%dim, 1:1))
        SAFE_ALLOCATE(kb_lr(1:sys%space%dim, 1:sys%space%dim, 1:1))
        do idir = 1, sys%space%dim
          do idir2 = 1, sys%space%dim
            call lr_init(kb_lr(idir, idir2, 1))
            call lr_allocate(kb_lr(idir, idir2, 1), sys%st, sys%gr%mesh, allocate_rho = .false.)
            if(idir2 <= idir) then
              call lr_init(k2_lr(idir, idir2, 1))
              call lr_allocate(k2_lr(idir, idir2, 1), sys%st, sys%gr%mesh, allocate_rho = .false.)
            end if 
          end do
        end do

        if(sys%space%periodic_dim < sys%space%dim) then
          if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
            message(1) = "All directions should be periodic for magnetic perturbations with kdotp."
          else
            message(1) = "All directions should be periodic for magnetooptics with kdotp."
          end if
          call messages_fatal(1)
        end if
        if(.not. complex_response) then
          do idir = 1, sys%space%dim
            call dlr_orth_response(sys%gr%mesh, sys%st, kdotp_lr(idir, 1), M_ZERO)
          end do
        else
          do idir = 1, sys%space%dim
            call zlr_orth_response(sys%gr%mesh, sys%st, kdotp_lr(idir, 1), frequency_zero)
          end do
        end if
        call sternheimer_init(sh_kmo, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
          complex_response, set_ham_var = 0, set_last_occ_response = em_vars%occ_response)  
      end if
    end if

    SAFE_ALLOCATE(em_vars%lr(1:sys%space%dim, 1:em_vars%nsigma, 1:em_vars%nfactor))
    do ifactor = 1, em_vars%nfactor
      call born_charges_init(em_vars%Born_charges(ifactor), sys%namespace, sys%ions, sys%st, sys%space%dim)
    end do

    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC &
      .and. sys%st%d%nspin == 1 .and. states_are_real(sys%st)) then
      ! first-order response is zero if there is time-reversal symmetry. F Mauri and SG Louie, PRL 76, 4246 (1996)
      call sternheimer_init(sh, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
        complex_response, set_ham_var = 0, set_last_occ_response = em_vars%occ_response)
      ! set HamiltonianVariation to V_ext_only, in magnetic case
    else
      call sternheimer_init(sh, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
        complex_response, set_last_occ_response = em_vars%occ_response)
      ! otherwise, use default, which is hartree + fxc
    end if

    if(em_vars%lrc_kernel .and. (.not. sternheimer_add_hartree(sh)) &
      .and. (.not. sternheimer_add_fxc(sh))) then
      message(1) = "Only the G = G'= 0 term of the LRC kernel is taken into account."
      call messages_warning(1)
    end if   
        

    if(mpi_grp_is_root(mpi_world)) then
      call info()
      call io_mkdir(EM_RESP_DIR, sys%namespace) ! output
    end if

    allocate_rho_em = sternheimer_add_fxc(sh) .or. sternheimer_add_hartree(sh)
    do ifactor = 1, em_vars%nfactor
      do idir = 1, sys%space%dim
        do sigma = 1, em_vars%nsigma
          call lr_init(em_vars%lr(idir, sigma, ifactor))
          call lr_allocate(em_vars%lr(idir, sigma, ifactor), sys%st, sys%gr%mesh, allocate_rho = allocate_rho_em)
        end do
      end do
    end do
    
    if((pert_type(em_vars%perturbation) .ne. PERTURBATION_ELECTRIC) .or. (.not. use_kdotp) &
      .or. sys%st%d%nik == 1) em_vars%kpt_output = .false.

    if(em_vars%kpt_output) then
      SAFE_ALLOCATE(em_vars%alpha_k(1:MAX_DIM, 1:MAX_DIM, 1:em_vars%nfactor, 1:sys%st%d%nik))
    end if

    if(em_vars%calc_magnetooptics) then
      if(em_vars%magnetooptics_nohvar) then
        call sternheimer_init(sh_mo, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
          complex_response, set_ham_var = 0, set_last_occ_response = em_vars%occ_response) 
      else
        call sternheimer_init(sh_mo, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
          complex_response, set_last_occ_response = em_vars%occ_response) 
        call sternheimer_build_kxc(sh_mo, sys%namespace, sys%gr%mesh, sys%st, sys%ks%xc)
      end if
      call messages_experimental("Magneto-optical response")
      allocate_rho_mo = sternheimer_add_fxc(sh_mo) .or. sternheimer_add_hartree(sh_mo)
      SAFE_ALLOCATE(b_lr(1:sys%space%dim, 1:1))
      do idir = 1, sys%space%dim
        call lr_init(b_lr(idir, 1))
        call lr_allocate(b_lr(idir, 1), sys%st, sys%gr%mesh, allocate_rho = allocate_rho_mo)
      end do
      
      if(use_kdotp) then
        if(em_vars%kpt_output) then 
          SAFE_ALLOCATE(em_vars%alpha_be_k(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM, 1:sys%st%d%nik))
        end if
        nfactor_ke = 1
        if(sys%kpoints%use_time_reversal .and. sys%kpoints%full%npoints > 1) nfactor_ke = em_vars%nfactor
        SAFE_ALLOCATE(ke_lr(1:sys%space%dim, 1:sys%space%dim, 1:em_vars%nsigma, 1:nfactor_ke))
        do idir = 1, sys%space%dim
          do idir2 = 1, sys%space%dim
            do sigma = 1, em_vars%nsigma
              do ifactor = 1, nfactor_ke              
                call lr_init(ke_lr(idir, idir2, sigma, ifactor))
                call lr_allocate(ke_lr(idir, idir2, sigma, ifactor), sys%st, sys%gr%mesh, allocate_rho = .false.)
              end do
            end do
          end do
        end do
      else
        call pert_init(pert_b, sys%namespace, PERTURBATION_MAGNETIC,  sys%gr, sys%ions)
      end if
    end if


    last_omega = M_HUGE
    do iomega = 1, em_vars%nomega

      em_vars%ok(1:3) = .true.

      do ifactor = 1, em_vars%nfactor
        frequency = em_vars%freq_factor(ifactor)*em_vars%omega(iomega)
        frequency_eta = frequency + M_zI * em_vars%eta
        if(em_vars%calc_magnetooptics .and. ifactor == 2) frequency_eta = frequency - M_zI * em_vars%eta

        if(abs(frequency) < M_EPSILON .and. em_vars%calc_magnetooptics .and. use_kdotp) then
          message(1) = "Magnetooptical response with kdotp requires non-zero frequency."
          call messages_warning(1)       
        end if

        ierr = 0
        ierr_e(:) = 0
        ierr_e2(:) = 0

        have_to_calculate = .true.
        opp_freq = .false.

        ! if this frequency is zero and this is not the first
        ! iteration we do not have to do anything
        if(iomega > 1 .and. abs(em_vars%freq_factor(ifactor)) <= M_EPSILON) have_to_calculate = .false. 

        if(ifactor > 1 .and. (.not. em_vars%calc_magnetooptics)) then 

          ! if this frequency is the same as the previous one, just copy it
          if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor - 1) * em_vars%omega(iomega) &
                                            - frequency) < M_EPSILON ) then

            do idir = 1, sys%space%dim
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor - 1), em_vars%lr(idir, 1, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor - 1), em_vars%lr(idir, 2, ifactor))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, sys%space%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 1, ifactor))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 2, ifactor))
                end do
              end if
            end do

            have_to_calculate = .false.

          end if

          ! if this frequency is minus the previous one, copy it inverted
          if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor - 1) * em_vars%omega(iomega) &
                                            + frequency) < M_EPSILON ) then 

            do idir = 1, sys%space%dim
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor - 1), em_vars%lr(idir, 2, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor - 1), em_vars%lr(idir, 1, ifactor))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, sys%space%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 2, ifactor))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 1, ifactor))
                end do
              end if
            end do

            have_to_calculate = .false.

          end if

        end if

        if(iomega > 1 .and. ifactor == 1 .and. (.not. em_vars%calc_magnetooptics)) then 

          ! if this frequency is the same as the previous one, just copy it
          if( have_to_calculate .and. abs(frequency - last_omega) < M_EPSILON ) then

            do idir = 1, sys%space%dim
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, em_vars%nfactor), em_vars%lr(idir, 1, 1))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, em_vars%nfactor), em_vars%lr(idir, 2, 1))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, sys%space%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 1, 1))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 2, 1))
                end do
              end if
            end do

            have_to_calculate = .false.

          end if

          ! if this frequency is minus the previous one, copy it inverted
          if( have_to_calculate .and. abs(frequency + last_omega) < M_EPSILON ) then

            do idir = 1, sys%space%dim
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, em_vars%nfactor), em_vars%lr(idir, 2, 1))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, em_vars%nfactor), em_vars%lr(idir, 1, 1))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, sys%space%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 2, 1))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 1, 1))
                end do
              end if
            end do

            have_to_calculate = .false.

          end if

        end if

        if(have_to_calculate) then 

          exact_freq(:) = .false.

          if(states_are_real(sys%st)) then
            call drun_sternheimer(em_vars, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, &
              sys%ions)
          else
            call zrun_sternheimer(em_vars, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, &
              sys%ions)
          end if

        end if ! have_to_calculate

        if(.not. have_to_calculate) cycle

        if(states_are_real(sys%st)) then
          call dcalc_properties_linear(em_vars, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, &
            sys%ions, sys%outp)
        else
          call zcalc_properties_linear(em_vars, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, &
            sys%ions, sys%outp)
        end if

      end do ! ifactor

      if(states_are_real(sys%st)) then
        call dcalc_properties_nonlinear(em_vars, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%ions)
      else
        call zcalc_properties_nonlinear(em_vars, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%ions)
      end if

      last_omega = em_vars%freq_factor(em_vars%nfactor) * em_vars%omega(iomega)

    end do ! iomega

    do idir = 1, sys%space%dim
      do sigma = 1, em_vars%nsigma
        do ifactor = 1, em_vars%nfactor
          call lr_dealloc(em_vars%lr(idir, sigma, ifactor))
        end do
      end do
    end do

    call sternheimer_end(sh)
    call pert_end(em_vars%perturbation)

    if(use_kdotp) then
      do idir = 1, sys%space%periodic_dim
        call lr_dealloc(kdotp_lr(idir, 1))
      end do
    end if

    if(em_vars%calc_hyperpol .and. use_kdotp) then
      call sternheimer_end(sh_kdotp)
      call sternheimer_end(sh2)
      call pert_end(pert_kdotp)
      call pert_end(pert2_none)
      do idir = 1, sys%space%periodic_dim
        do idir2 = 1, sys%space%periodic_dim
          do sigma = 1, em_vars%nsigma
            do ifactor = 1, em_vars%nfactor
              call lr_dealloc(kdotp_em_lr2(idir, idir2, sigma, ifactor))
            end do
          end do
        end do
      end do
      SAFE_DEALLOCATE_A(kdotp_em_lr2)
      SAFE_DEALLOCATE_A(dl_eig)
    end if

    if(em_vars%kpt_output) then
      SAFE_DEALLOCATE_A(em_vars%alpha_k)
    end	if

    if(em_vars%calc_magnetooptics .or. &
      (pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC)) then
      if(use_kdotp) then  
        call pert_end(pert2_none) 
        call sternheimer_end(sh_kmo)
        do idir = 1, sys%space%dim
          do idir2 = 1, sys%space%dim 
            call lr_dealloc(kb_lr(idir, idir2, 1))
            if(idir2 <= idir) call lr_dealloc(k2_lr(idir, idir2, 1))
          end do
        end do
        SAFE_DEALLOCATE_A(k2_lr)
        SAFE_DEALLOCATE_A(kb_lr)
      end if
    end if

    if(em_vars%calc_magnetooptics) then
      if(.not. em_vars%magnetooptics_nohvar) call sternheimer_unset_kxc(sh_mo)
      call sternheimer_end(sh_mo)
      do idir = 1, sys%space%dim
        call lr_dealloc(b_lr(idir, 1))
      end do
      SAFE_DEALLOCATE_A(b_lr)
      
      if(use_kdotp) then
        do idir = 1, sys%space%dim 
          do idir2 = 1, sys%space%dim
            do sigma = 1, em_vars%nsigma
              do ifactor = 1, nfactor_ke
                call lr_dealloc(ke_lr(idir, idir2, sigma, ifactor))
              end do
            end do
          end do
        end do
        SAFE_DEALLOCATE_A(ke_lr)
        if(em_vars%kpt_output) then
          SAFE_DEALLOCATE_A(em_vars%alpha_be_k)
        end if        
      else
        call pert_end(pert_b)
      end if
    end if

    SAFE_DEALLOCATE_A(em_vars%omega)
    SAFE_DEALLOCATE_A(em_vars%lr)
    do ifactor = 1, em_vars%nfactor
      call Born_charges_end(em_vars%Born_charges(ifactor))
    end do
    call states_elec_deallocate_wfns(sys%st)

    POP_SUB(em_resp_run_legacy)

  contains

    ! ---------------------------------------------------------
    subroutine parse_input()

      type(block_t) :: blk
      integer :: nrow, irow, nfreqs_in_row, ifreq, istep, perturb_type
      FLOAT :: omega_ini, omega_fin, domega
      logical :: freq_sort  

      PUSH_SUB(em_resp_run_legacy.parse_input)

      call messages_obsolete_variable(sys%namespace, 'PolFreqs               ', 'EMFreqs             ')
      call messages_obsolete_variable(sys%namespace, 'PolHyper               ', 'EMHyperpol          ')
      call messages_obsolete_variable(sys%namespace, 'PolEta                 ', 'EMEta               ')
      call messages_obsolete_variable(sys%namespace, 'PolHamiltonianVariation', 'HamiltonianVariation')

      !%Variable EMFreqs
      !%Type block
      !%Section Linear Response::Polarizabilities
      !%Description
      !% This block defines for which frequencies the polarizabilities
      !% will be calculated. If it is not present, the static (<math>\omega = 0</math>) response
      !% is calculated.
      !%
      !% Each row of the block indicates a sequence of frequency values, the
      !% first column is an integer that indicates the number of steps, the
      !% second number is the initial frequency, and the third number the final
      !% frequency. If the first number is one, then only the initial value is
      !% considered. The block can have any number of rows. Consider the next example:
      !%
      !% <tt>%EMFreqs
      !% <br>31 | 0.0 | 1.0
      !% <br> 1 | 0.32
      !% <br>%</tt>
      !%
      !%End

      if(parse_block(sys%namespace, 'EMFreqs', blk) == 0) then 

        nrow = parse_block_n(blk)
        em_vars%nomega = 0

        !count the number of frequencies
        do irow = 0, nrow-1
          call parse_block_integer(blk, irow, 0, nfreqs_in_row)
          if(nfreqs_in_row < 1) then
            message(1) = "EMFreqs: invalid number of frequencies."
            call messages_fatal(1)
          end if
          em_vars%nomega = em_vars%nomega + nfreqs_in_row
        end do

        SAFE_ALLOCATE(em_vars%omega(1:em_vars%nomega))

        !read frequencies
        ifreq = 1
        do irow = 0, nrow-1
          call parse_block_integer(blk, irow, 0, nfreqs_in_row)
          call parse_block_float(blk, irow, 1, omega_ini)
          if(nfreqs_in_row > 1) then 
            call parse_block_float(blk, irow, 2, omega_fin)
            domega = (omega_fin - omega_ini)/(nfreqs_in_row - M_ONE)
            do istep = 0, nfreqs_in_row-1
              em_vars%omega(ifreq + istep) = units_to_atomic(units_inp%energy, omega_ini + domega*istep)
            end do
            ifreq = ifreq + nfreqs_in_row
          else
            em_vars%omega(ifreq) = units_to_atomic(units_inp%energy, omega_ini)
            ifreq = ifreq + 1
          end if
        end do

        call parse_block_end(blk)

        !%Variable EMFreqsSort
        !%Type logical
        !%Default true
        !%Section Linear Response::Polarizabilities
        !%Description
        !% If true, the frequencies specified by the <tt>EMFreqs</tt> block are sorted, so that
        !% they are calculated in increasing order. Can be set to false to use the order as stated,
        !% in case this makes better use of available restart information.
        !%End
        call parse_variable(sys%namespace, 'EMFreqsSort', .true., freq_sort)

        if(freq_sort) call sort(em_vars%omega)

      else
        !there is no frequency block, we calculate response for w = 0.0
        em_vars%nomega = 1
        SAFE_ALLOCATE(em_vars%omega(1:em_vars%nomega))
        em_vars%omega(1) = M_ZERO
      end if

      !%Variable EMEta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::Polarizabilities
      !%Description
      !% The imaginary part of the frequency, effectively a Lorentzian broadening
      !% for peaks in the spectrum. It can help convergence of the SCF cycle for the
      !% Sternheimer equation when on a resonance, and it can be used as a positive
      !% infinitesimal to get the imaginary parts of response functions at poles.
      !% In units of energy. Cannot be negative.
      !%End

      call parse_variable(sys%namespace, 'EMEta', M_ZERO, em_vars%eta, units_inp%energy)
      if(em_vars%eta < -M_EPSILON) then
        message(1) = "EMEta cannot be negative."
        call messages_fatal(1)
      end if

      ! reset the values of these variables
      em_vars%calc_hyperpol = .false.
      em_vars%freq_factor(1:3) = M_ONE
      em_vars%calc_magnetooptics = .false.
      em_vars%magnetooptics_nohvar = .true.
      em_vars%kpt_output = .false.

      !%Variable EMPerturbationType
      !%Type integer
      !%Default electric
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Which perturbation to consider for electromagnetic linear response.
      !%Option electric 1
      !% Electric perturbation used to calculate electric polarizabilities
      !% and hyperpolarizabilities.
      !%Option magnetic 2
      !% Magnetic perturbation used to calculate magnetic susceptibilities.
      !%Option none 0
      !% Zero perturbation, for use in testing.
      !%End 
      call parse_variable(sys%namespace, 'EMPerturbationType', PERTURBATION_ELECTRIC, perturb_type)
      call messages_print_var_option(stdout, 'EMPerturbationType', perturb_type)
      
      call pert_init(em_vars%perturbation, sys%namespace, perturb_type, sys%gr, sys%ions)

      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        !%Variable EMCalcRotatoryResponse
        !%Type logical
        !%Default false
        !%Section Linear Response::Polarizabilities
        !%Description
        !% Calculate circular-dichroism spectrum from electric perturbation,
        !% and write to file <tt>rotatory_strength</tt>.
        !%End

        call parse_variable(sys%namespace, 'EMCalcRotatoryResponse', .false., em_vars%calc_rotatory)

        !%Variable EMHyperpol
        !%Type block
        !%Section Linear Response::Polarizabilities
        !%Description
        !% This block describes the multiples of the frequency used for
        !% the dynamic hyperpolarizability. The results are written to the
        !% file <tt>beta</tt> in the directory for the first multiple.
        !% There must be three factors, summing to zero: <math>\omega_1 + \omega_2 + \omega_3 = 0</math>.
        !% For example, for second-harmonic generation, you could use
        !% <tt>1 | 1 | -2</tt>.
        !%End

        if (parse_block(sys%namespace, 'EMHyperpol', blk) == 0) then 
          call parse_block_float(blk, 0, 0, em_vars%freq_factor(1))
          call parse_block_float(blk, 0, 1, em_vars%freq_factor(2))
          call parse_block_float(blk, 0, 2, em_vars%freq_factor(3))

          call parse_block_end(blk)

          if(abs(sum(em_vars%freq_factor(1:3))) > M_EPSILON) then
            message(1) = "Frequency factors specified by EMHyperpol must sum to zero."
            call messages_fatal(1)
          end if

          em_vars%calc_hyperpol = .true.
        end if

        !%Variable EMCalcMagnetooptics
        !%Type logical
        !%Default false
        !%Section Linear Response::Polarizabilities
        !%Description
        !% Calculate magneto-optical response.
        !%End
        call parse_variable(sys%namespace, 'EMCalcMagnetooptics', .false., em_vars%calc_magnetooptics)

        !%Variable EMMagnetoopticsNoHVar
        !%Type logical
        !%Default true
        !%Section Linear Response::Polarizabilities
        !%Description
        !% Exclude corrections to the exchange-correlation and Hartree terms 
        !% from consideration of perturbations induced by a magnetic field
        !%End
        call parse_variable(sys%namespace, 'EMMagnetoopticsNoHVar', .true., em_vars%magnetooptics_nohvar)

        !%Variable EMKPointOutput
        !%Type logical
        !%Default false
        !%Section Linear Response::Polarizabilities
        !%Description
        !% Give in the output contributions of different k-points to the dielectric constant.
        !% Can be also used for magneto-optical effects.
        !%End

        call parse_variable(sys%namespace, 'EMKPointOutput', .false., em_vars%kpt_output)

      end if

      !%Variable EMForceNoKdotP
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% If the system is periodic, by default wavefunctions from a previous <tt>kdotp</tt> run will
      !% be read, to be used in the formulas for the polarizability and
      !% hyperpolarizability in the quantum theory of polarization. For testing purposes,
      !% you can set this variable to true to disregard the <tt>kdotp</tt> run, and use the formulas
      !% for the finite system. This variable has no effect for a finite system.
      !%End

      call parse_variable(sys%namespace, 'EMForceNoKdotP', .false., em_vars%force_no_kdotp)

      !%Variable EMCalcBornCharges
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Calculate linear-response Born effective charges from electric perturbation (experimental).
      !%End

      call parse_variable(sys%namespace, 'EMCalcBornCharges', .false., em_vars%calc_Born)
      if (em_vars%calc_Born) call messages_experimental("Calculation of Born effective charges")

      !%Variable EMOccupiedResponse
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Solve for full response without projector into unoccupied subspace.
      !% Not possible if there are partial occupations.
      !% When <tt>EMHyperpol</tt> is set for a periodic system, this variable is ignored and
      !% the full response is always calculated.
      !%End

      call parse_variable(sys%namespace, 'EMOccupiedResponse', .false., em_vars%occ_response)
      if(em_vars%occ_response .and. .not. (smear_is_semiconducting(sys%st%smear) .or. sys%st%smear%method == SMEAR_FIXED_OCC)) then
        message(1) = "EMOccupiedResponse cannot be used if there are partial occupations."
        call messages_fatal(1)
      end if

      !%Variable EMWavefunctionsFromScratch
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Do not use saved linear-response wavefunctions from a previous run as starting guess.
      !% Instead initialize to zero as in <tt>FromScratch</tt>, but restart densities will still
      !% be used. Restart wavefunctions from a very different frequency can hinder convergence.
      !%End

      call parse_variable(sys%namespace, 'EMWavefunctionsFromScratch', .false., em_vars%wfns_from_scratch)

      POP_SUB(em_resp_run_legacy.parse_input)

    end subroutine parse_input


    ! ---------------------------------------------------------
    subroutine info()

      PUSH_SUB(em_resp_run_legacy.info)

      call pert_info(em_vars%perturbation)
      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        if(em_vars%calc_hyperpol) then 
          write(message(1),'(a)') 'Linear-Response First-Order Hyperpolarizabilities'
          call messages_print_stress(stdout, trim(message(1)))
        else 
          write(message(1),'(a)') 'Linear-Response Polarizabilities'
          call messages_print_stress(stdout, trim(message(1)))
        end if
      else
        write(message(1),'(a)') 'Magnetic Susceptibilities'
        call messages_print_stress(stdout, trim(message(1)))
      end if

      if (states_are_real(sys%st)) then 
        message(1) = 'Wavefunctions type: Real'
      else
        message(1) = 'Wavefunctions type: Complex'
      end if
      call messages_info(1)

      write(message(1),'(a,i3,a)') 'Calculating response for ', em_vars%nomega, ' frequencies.'
      call messages_info(1)

      call messages_print_stress(stdout)

      POP_SUB(em_resp_run_legacy.info)

    end subroutine info

! Note: unlike the typical usage, here the templates make 'internal procedures'
#include "undef.F90"
#include "real.F90"
#include "em_resp_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "em_resp_inc.F90"

#include "undef.F90"
  end subroutine em_resp_run_legacy


  ! ---------------------------------------------------------
  subroutine em_resp_output(st, namespace, space, gr, hm, ions, outp, sh, em_vars, iomega, ifactor)
    type(states_elec_t),      intent(inout) :: st
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(ions_t),             intent(in)    :: ions
    type(output_t),           intent(in)    :: outp
    type(sternheimer_t),      intent(in)    :: sh
    type(em_resp_t),          intent(inout) :: em_vars
    integer,                  intent(in)    :: iomega
    integer,                  intent(in)    :: ifactor
    
    integer :: iunit
    character(len=80) :: dirname, str_tmp
    logical :: use_kdotp
    CMPLX :: epsilon(MAX_DIM, MAX_DIM) 

    PUSH_SUB(em_resp_output)

    use_kdotp = space%is_periodic() .and. .not. em_vars%force_no_kdotp

    str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(ifactor)*em_vars%omega(iomega)))
    if(em_vars%calc_magnetooptics) str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%omega(iomega)))
    write(dirname, '(a, a)') EM_RESP_DIR//'freq_', trim(str_tmp)
    call io_mkdir(trim(dirname), namespace)

    if (sh%enable_el_pt_coupling) then
      iunit = io_open(trim(dirname)//'/photon_coord_q', namespace, action='write')
      write(iunit, '(a)') 'Photon coordinate Q [', trim(units_abbrev(units_out%energy)), ']'
      write(iunit, '(a)') '                 Re                Im'
      write(iunit, '(f20.6,f20.6)') units_from_atomic(units_out%energy, real(sh%zphoton_coord_q)), &
                                    units_from_atomic(units_out%energy, aimag(sh%zphoton_coord_q))
      call io_close(iunit)
    end if

    call write_eta()

    if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
      if((.not. em_vars%calc_magnetooptics) .or. ifactor == 1) then
        call out_polarizability()
        if(em_vars%calc_Born) then
          call out_Born_charges(em_vars%Born_charges(ifactor), ions, namespace, space%dim, dirname, &
            write_real = em_vars%eta < M_EPSILON)
        end if

        if (space%periodic_dim  ==  space%dim) then
          call out_dielectric_constant()
        end if

        if ((.not. space%is_periodic() .or. em_vars%force_no_kdotp) .and. em_vars%calc_rotatory) then
          call out_circular_dichroism()
        end if

      else 
        call out_magnetooptics() 
        if(iomega == 1) call out_susceptibility()
      end if

    else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
      call out_susceptibility()
    end if

    call out_wfn_and_densities()

    POP_SUB(em_resp_output)

  contains


    subroutine write_eta()

      PUSH_SUB(em_resp_output.write_eta)

      iunit = io_open(trim(dirname)//'/eta', namespace, action='write')

      write(iunit, '(3a)') 'Imaginary part of frequency [', trim(units_abbrev(units_out%energy)), ']'
      write(iunit, '(f20.6)') units_from_atomic(units_out%energy, em_vars%eta)

      call io_close(iunit)

      POP_SUB(em_resp_output.write_eta)
    end subroutine write_eta


    ! ---------------------------------------------------------
    !> Note: this should be in spectrum.F90
    subroutine cross_section_header(out_file)
      integer, intent(in) :: out_file

      character(len=80) :: header_string
      integer :: ii, idir, kdir

      PUSH_SUB(em_resp_output.cross_section_header)

      !this header is the same as spectrum.F90
      write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
      write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
      write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)

      do idir = 1, space%dim
        do kdir = 1, space%dim
          write(header_string,'(a6,i1,a1,i1,a1)') 'sigma(', idir, ',', kdir, ')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do

      write(out_file, *)
      write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
      do ii = 1, 2 + space%dim**2
        write(out_file, '(a20)', advance = 'no')  str_center('['//trim(units_abbrev(units_out%length**2)) // ']', 20)
      end do
      write(out_file,*)

      POP_SUB(em_resp_output.cross_section_header)
    end subroutine cross_section_header


    ! ---------------------------------------------------------
    subroutine out_polarizability()
      FLOAT :: cross(MAX_DIM, MAX_DIM), crossp(MAX_DIM, MAX_DIM)
      FLOAT :: cross_sum, crossp_sum, anisotropy
      integer :: idir, idir2
      
      PUSH_SUB(em_resp_output.out_polarizability)
  
      iunit = io_open(trim(dirname)//'/alpha', namespace, action='write')
  
      if (.not. em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
      write(iunit, '(3a)') '# Polarizability tensor [', trim(units_abbrev(units_out%polarizability)), ']'
      call output_tensor(iunit, TOFLOAT(em_vars%alpha(:, :, ifactor)), space%dim, units_out%polarizability)
  
      call io_close(iunit)
  
      ! CROSS SECTION (THE IMAGINARY PART OF POLARIZABILITY)
      if(em_vars%eta > M_EPSILON) then 
        cross(1:space%dim, 1:space%dim) = aimag(em_vars%alpha(1:space%dim, 1:space%dim, ifactor)) * &
          em_vars%freq_factor(ifactor) * em_vars%omega(iomega) * (M_FOUR * M_PI / P_c)

        do idir = 1, space%dim
          do idir2 = 1, space%dim
            cross(idir, idir2) = units_from_atomic(units_out%length**2, cross(idir, idir2))
          end do
        end do

        iunit = io_open(trim(dirname)//'/cross_section', namespace, action='write')
        if (.not. em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
        crossp(1:space%dim, 1:space%dim) = matmul(cross(1:space%dim, 1:space%dim), cross(1:space%dim, 1:space%dim))

        cross_sum = M_ZERO
        crossp_sum = M_ZERO
        do idir = 1, space%dim
          cross_sum = cross_sum + cross(idir, idir)
          crossp_sum = crossp_sum + crossp(idir, idir)
        end do

        anisotropy = crossp_sum - M_THIRD * cross_sum**2
            
        call cross_section_header(iunit)
        write(iunit,'(3e20.8)', advance = 'no') &
          units_from_atomic(units_out%energy, em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), &
          cross_sum * M_THIRD, sqrt(max(anisotropy, M_ZERO))
        do idir = 1, space%dim
          do idir2 = 1, space%dim
            write(iunit,'(e20.8)', advance = 'no') cross(idir, idir2)
          end do
        end do
        write(iunit,'(a)', advance = 'yes')
  
        call io_close(iunit)
      end if
      
      POP_SUB(em_resp_output.out_polarizability)
    end subroutine out_polarizability


    ! ---------------------------------------------------------
    !> epsilon = 1 + 4 * pi * alpha/volume
    subroutine out_dielectric_constant()
      integer :: idir, idir1, ik
      character(len=80) :: header_string
      CMPLX, allocatable :: epsilon_k(:, :, :)

      PUSH_SUB(em_resp_output.out_dielectric_constant)
  
      iunit = io_open(trim(dirname)//'/epsilon', namespace, action='write')
      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
      epsilon(1:space%dim, 1:space%dim) = &
        4 * M_PI * em_vars%alpha(1:space%dim, 1:space%dim, ifactor) / ions%latt%rcell_volume
      do idir = 1, space%dim
        epsilon(idir, idir) = epsilon(idir, idir) + M_ONE
      end do

      write(iunit, '(a)') '# Real part of dielectric constant'
      call output_tensor(iunit, TOFLOAT(epsilon(1:space%dim, 1:space%dim)), space%dim, unit_one)
      write(iunit, '(a)')
      write(iunit, '(a)') '# Imaginary part of dielectric constant'
      call output_tensor(iunit, aimag(epsilon(1:space%dim, 1:space%dim)), space%dim, unit_one)

      if(em_vars%lrc_kernel) then
        write(iunit, '(a)')
        write(iunit, '(a)') '# Without G = G'' = 0 term of the LRC kernel'

        epsilon(1:space%dim, 1:space%dim) = &
          4 * M_PI * em_vars%alpha0(1:space%dim, 1:space%dim, ifactor) / ions%latt%rcell_volume
        do idir = 1, space%dim
          epsilon(idir, idir) = epsilon(idir, idir) + M_ONE
        end do

        write(iunit, '(a)') '# Real part of dielectric constant'
        call output_tensor(iunit, TOFLOAT(epsilon(1:space%dim, 1:space%dim)), space%dim, unit_one)
        write(iunit, '(a)')
        write(iunit, '(a)') '# Imaginary part of dielectric constant'
        call output_tensor(iunit, aimag(epsilon(1:space%dim, 1:space%dim)), space%dim, unit_one)
      end if

      call io_close(iunit)

      if(em_vars%kpt_output) then
        SAFE_ALLOCATE(epsilon_k(1:space%dim, 1:space%dim, 1:hm%kpoints%reduced%npoints))
        do ik = 1, hm%kpoints%reduced%npoints
          do idir = 1, space%dim
            do idir1 = 1, space%dim
              epsilon_k(idir, idir1, ik) = M_FOUR * M_PI * em_vars%alpha_k(idir, idir1, ifactor, ik) / ions%latt%rcell_volume
            end do
          end do
        end do
        iunit = io_open(trim(dirname)//'/epsilon_k_re', namespace, action='write')

        write(iunit, '(a)') '# Real part of dielectric constant'
        write(iunit, '(a10)', advance = 'no') '#  index  '
        write(iunit, '(a20)', advance = 'no') str_center("weight", 20)
        write(iunit, '(a20)', advance = 'no') str_center("kx", 20)        
        write(iunit, '(a20)', advance = 'no') str_center("ky", 20)     
        write(iunit, '(a20)', advance = 'no') str_center("kz", 20)

        do idir = 1, space%dim
          do idir1 = 1, space%dim
            write(header_string,'(a7,i1,a1,i1,a1)') 'Re eps(', idir, ',', idir1, ')'
            write(iunit, '(a20)', advance = 'no') str_center(trim(header_string), 20)
          end do
        end do
        write(iunit, *)

        do ik = 1, hm%kpoints%reduced%npoints
          write(iunit, '(i8)', advance = 'no') ik
          write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%weight(ik)
          do idir = 1, space%dim
            write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%red_point(idir, ik)
          end do
          do idir = 1, space%dim
            do idir1 = 1, space%dim
              write(iunit, '(e20.8)', advance = 'no') TOFLOAT(epsilon_k(idir, idir1, ik))
            end do
          end do
          write(iunit, *)
        end do
        call io_close(iunit)

        iunit = io_open(trim(dirname)//'/epsilon_k_im', namespace, action='write')

        write(iunit, '(a)') '# Imaginary part of dielectric constant'
        write(iunit, '(a10)', advance = 'no') '#  index  '
        write(iunit, '(a20)', advance = 'no') str_center("weight", 20)
        write(iunit, '(a20)', advance = 'no') str_center("kx", 20)        
        write(iunit, '(a20)', advance = 'no') str_center("ky", 20)     
        write(iunit, '(a20)', advance = 'no') str_center("kz", 20)

        do idir = 1, space%dim
          do idir1 = 1, space%dim
            write(header_string,'(a7,i1,a1,i1,a1)') 'Im eps(', idir, ',', idir1,')'
            write(iunit, '(a20)', advance = 'no') str_center(trim(header_string), 20)
          end do
        end do
        write(iunit, *)

        do ik = 1, hm%kpoints%reduced%npoints                              
          write(iunit, '(i8)', advance = 'no') ik
          write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%weight(ik)
          do idir = 1, space%dim
            write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%red_point(idir, ik)
          end do
          do idir = 1, space%dim
            do idir1 = 1, space%dim
              write(iunit, '(e20.8)', advance = 'no') aimag(epsilon_k(idir, idir1, ik))
            end do
          end do
          write(iunit, *)          
        end do
        call io_close(iunit)
        SAFE_DEALLOCATE_A(epsilon_k)
      end if

      POP_SUB(em_resp_output.out_dielectric_constant)
    end subroutine out_dielectric_constant


    ! ---------------------------------------------------------
    subroutine out_susceptibility()

      character(len=80) :: dirname1

      PUSH_SUB(em_resp_output.out_susceptibility)

      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        write(dirname1, '(a)') EM_RESP_DIR//'freq_0.0000'
        call io_mkdir(trim(dirname1), namespace)
        iunit = io_open(trim(dirname1)//'/susceptibility', namespace, action='write')
      else  
        iunit = io_open(trim(dirname)//'/susceptibility', namespace, action='write')
      end if

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      ! There is no separation into the diamagnetic and paramagnetic terms in the expression 
      ! for periodic systems 
      if(.not. use_kdotp) then
        write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm a.u.]'
        call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :)), space%dim, unit_ppm)
        write(iunit, '(1x)')

        write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm a.u.]'
        call output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :)), space%dim, unit_ppm)
        write(iunit, '(1x)')
      end if

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm a.u.]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :) + em_vars%chi_dia(:,:)), &
        space%dim, unit_ppm)
      write(iunit, '(1x)')

      write(iunit, '(a)') hyphens

      if(.not. use_kdotp) then
        write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
        call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :)), space%dim, unit_susc_ppm_cgs)
        write(iunit, '(1x)')

        write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
        call output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :)), space%dim, unit_susc_ppm_cgs)
        write(iunit, '(1x)')
      end if

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm cgs / mol]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :) + em_vars%chi_dia(:,:)), &
           space%dim, unit_susc_ppm_cgs)
      write(iunit, '(1x)')

      if(use_kdotp) then
        write(iunit, '(a)') hyphens
        write(iunit, '(1a)') '# Magnetization [ppm a.u.]'
        write(iunit, '(3f20.8)') units_from_atomic(unit_ppm, TOFLOAT(em_vars%magn(1))), &
          units_from_atomic(unit_ppm, TOFLOAT(em_vars%magn(2))), &
          units_from_atomic(unit_ppm, TOFLOAT(em_vars%magn(3)))
      end if

      call io_close(iunit)      
      POP_SUB(em_resp_output.out_susceptibility)
    end subroutine out_susceptibility


    ! ---------------------------------------------------------
    subroutine out_wfn_and_densities()
      integer :: idir, isigma

      PUSH_SUB(em_resp_output.out_wfn_and_densities)

      do idir = 1, space%dim
        if (states_are_complex(st)) then 

          if (space%dim == 3 .and. outp%what(OPTION__OUTPUT__ELF)) then
            if (em_vars%nsigma == 1) then
              call zlr_calc_elf(st, gr, hm%kpoints, em_vars%lr(idir, 1, ifactor))
            else
              call zlr_calc_elf(st, gr, hm%kpoints, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
            end if
          end if
          do isigma = 1, em_vars%nsigma
            call zoutput_lr(outp, namespace, space, dirname, st, gr%mesh, em_vars%lr(idir, isigma, ifactor), idir, isigma, ions, &
              units_out%force)
          end do
        else

          if (space%dim == 3 .and. outp%what(OPTION__OUTPUT__ELF)) then
            if (em_vars%nsigma == 1) then
              call dlr_calc_elf(st, gr, hm%kpoints, em_vars%lr(idir, 1, ifactor))
            else
              call dlr_calc_elf(st, gr, hm%kpoints, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
            end if
          end if

          do isigma = 1, em_vars%nsigma
            call doutput_lr(outp, namespace, space, dirname, st, gr%mesh, em_vars%lr(idir, isigma, ifactor), idir, isigma, ions, &
              units_out%force)
          end do

        end if
      end do

      POP_SUB(em_resp_output.out_wfn_and_densities)

    end subroutine out_wfn_and_densities
    

  ! ---------------------------------------------------------
  !> See D Varsano, LA Espinosa Leal, Xavier Andrade, MAL Marques, Rosa di Felice, Angel Rubio,
  !! Phys. Chem. Chem. Phys. 11, 4481 (2009)
    subroutine out_circular_dichroism()
      
      type(pert_t) :: angular_momentum
      integer :: idir
      FLOAT :: ff
      CMPLX :: dic
      CMPLX, allocatable :: psi(:, :, :, :)
  
      PUSH_SUB(em_resp_output.out_circular_dichroism)

      if(states_are_complex(st) .and. em_vars%nsigma == 2) then       

        message(1) = "Info: Calculating rotatory response."
        call messages_info(1)

        call pert_init(angular_momentum, namespace, PERTURBATION_MAGNETIC, gr, ions)
        
        SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

        call states_elec_get_state(st, gr%mesh, psi)

        dic = M_ZERO
        do idir = 1, space%dim
          call pert_setup_dir(angular_momentum, idir)
          dic = dic &
            + zpert_expectation_value(angular_momentum, namespace, gr, ions, hm, st, &
            psi, em_vars%lr(idir, 1, ifactor)%zdl_psi) &
            + zpert_expectation_value(angular_momentum, namespace, gr, ions, hm, st, &
            em_vars%lr(idir, 2, ifactor)%zdl_psi, psi)
        end do

        SAFE_DEALLOCATE_A(psi)
        
        call pert_end(angular_momentum)
        
        dic = dic*M_zI*M_HALF

        iunit = io_open(trim(dirname)//'/rotatory_strength', namespace, action='write')

        ! print header
        write(iunit, '(a1,a20,a20,a20)') '#', str_center("Energy", 20), str_center("R", 20), str_center("Re[beta]", 20)
        write(iunit, '(a1,a20,a20,a20)') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
             str_center('['//trim(units_abbrev(units_out%length**3)) //']', 20), &
             str_center('['//trim(units_abbrev(units_out%length**4)) //']', 20)

        ff = M_ZERO
        if(em_vars%omega(iomega) /= 0) ff = TOFLOAT(dic)/(M_THREE*em_vars%omega(iomega))

        write(iunit, '(3e20.8)') units_from_atomic(units_out%energy, em_vars%omega(iomega)), &
             units_from_atomic(units_out%length**3, aimag(dic)/(P_C*M_PI)), units_from_atomic(units_out%length**4, ff)

        call io_close(iunit)
      end if
      
      POP_SUB(em_resp_output.out_circular_dichroism)

    end subroutine out_circular_dichroism
    
   ! ---------------------------------------------------------
    subroutine out_magnetooptics   
      integer :: idir, ik
      CMPLX :: epsilon_m(4), diff(4), eps_mk(MAX_DIM)
      
      PUSH_SUB(em_resp_output.out_magnetooptics)
      
      ! This code assumes 3D
      ASSERT(space%dim == 3)

      diff(:) = M_ZERO
      do idir = 1, space%dim
        diff(idir) = M_HALF * (em_vars%alpha_be(magn_dir(idir, 1), magn_dir(idir, 2), idir) - &
          em_vars%alpha_be(magn_dir(idir, 2), magn_dir(idir, 1), idir))
      end do
      diff(4) = (diff(1) + diff(2) + diff(3)) / M_THREE

      iunit = io_open(trim(dirname)//'/alpha_mo', namespace, action='write')

      if (.not. em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      write(iunit, '(a1, a25)', advance = 'no') '#', str_center(" ", 25)
      write(iunit, '(a20)', advance = 'no') str_center("   yz,x = -zy,x", 20)
      write(iunit, '(a20)', advance = 'no') str_center("   zx,y = -xz,y", 20)
      write(iunit, '(a20)', advance = 'no') str_center("   xy,z = -yx,z", 20)
      write(iunit, '(a20)', advance = 'no') str_center(" Average", 20)
      write(iunit, *)
 
      write(iunit, '(a25)', advance = 'no') str_center("Re alpha [a.u.]", 25)
      do idir = 1, space%dim + 1 
        write(iunit, '(e20.8)', advance = 'no') TOFLOAT(diff(idir))
      end do
      write(iunit, *)

      write(iunit, '(a25)', advance = 'no') str_center("Im alpha [a.u.]", 25)
      do idir = 1, space%dim + 1
        write(iunit, '(e20.8)', advance = 'no') aimag(diff(idir))
      end do
      write(iunit, *)

      if (space%is_periodic()) then
        ! This code assumes 3D periodic
        ASSERT(space%periodic_dim == 3)

        do idir = 1, space%dim
          epsilon_m(idir) = 4 * M_PI * diff(idir) / ions%latt%rcell_volume
        end do
        epsilon_m(4) = 4 * M_PI * diff(4) / ions%latt%rcell_volume

        write(iunit, '(a25)', advance = 'no') str_center("Re epsilon (B = 1 a.u.)", 25)
        do idir = 1, space%dim + 1
          write(iunit, '(e20.8)', advance = 'no') TOFLOAT(epsilon_m(idir))
        end do
        write(iunit, *)

        write(iunit, '(a25)', advance = 'no') str_center("Im epsilon (B = 1 a.u.)", 25)
        do idir = 1, space%dim + 1
          write(iunit, '(e20.8)', advance = 'no') aimag(epsilon_m(idir))
        end do
        write(iunit, *)

        if (em_vars%lrc_kernel) then
          write(iunit, '(a)')
          write(iunit, '(a)') '# Without the G = G'' = 0 term of the LRC kernel'

          diff(:) = M_ZERO
          epsilon_m(:) = M_ZERO
          do idir = 1, space%dim
            diff(idir) = M_HALF * (em_vars%alpha_be0(magn_dir(idir, 1), magn_dir(idir, 2), idir) - &
              em_vars%alpha_be0(magn_dir(idir, 2), magn_dir(idir, 1), idir))

            epsilon_m(idir) = 4 * M_PI * diff(idir) / ions%latt%rcell_volume
          end do
          diff(4) =(diff(1) + diff(2) + diff(3)) / M_THREE
          epsilon_m(4) = 4 * M_PI * diff(4) / ions%latt%rcell_volume

          write(iunit, '(a1, a25)', advance = 'no') '#', str_center(" ", 25)
          write(iunit, '(a20)', advance = 'no') str_center("   yz,x = -zy,x", 20)
          write(iunit, '(a20)', advance = 'no') str_center("   zx,y = -xz,y", 20)
          write(iunit, '(a20)', advance = 'no') str_center("   xy,z = -yx,z", 20)
          write(iunit, '(a20)', advance = 'no') str_center(" Average", 20)
          write(iunit, *)

          write(iunit, '(a25)', advance = 'no') str_center("Re alpha [a.u.]", 25)
          do idir = 1, space%dim + 1
            write(iunit, '(e20.8)', advance = 'no') TOFLOAT(diff(idir))
          end do
          write(iunit, *)

          write(iunit, '(a25)', advance = 'no') str_center("Im alpha [a.u.]", 25)
          do idir = 1, space%dim + 1
            write(iunit, '(e20.8)', advance = 'no') aimag(diff(idir))
          end do
          write(iunit, *)

          write(iunit, '(a25)', advance = 'no') str_center("Re epsilon (B = 1 a.u.)", 25)
          do idir = 1, space%dim + 1
            write(iunit, '(e20.8)', advance = 'no') TOFLOAT(epsilon_m(idir))
          end do
          write(iunit, *)

          write(iunit, '(a25)', advance = 'no') str_center("Im epsilon (B = 1 a.u.)", 25)
          do idir = 1, space%dim + 1
            write(iunit, '(e20.8)', advance = 'no') aimag(epsilon_m(idir))
          end do
          write(iunit, *)
        end if
      end if
      call io_close(iunit)

      if (space%is_periodic() .and. em_vars%kpt_output) then
	iunit = io_open(trim(dirname)//'/epsilon_mo_k', namespace, action='write')

        write(iunit, '(a)') '# Contribution to dielectric tensor for B = 1 a.u.'
        write(iunit, '(a10)', advance = 'no') '#  index  '
        write(iunit, '(a20)', advance = 'no') str_center("weight", 20)
        write(iunit, '(a20)', advance = 'no') str_center("kx", 20)
        write(iunit, '(a20)', advance = 'no') str_center("ky", 20)
        write(iunit, '(a20)', advance = 'no') str_center("kz", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Re eps_yz,x", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Re eps_zx,y", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Re eps_xy,z", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Im eps_yz,x", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Im eps_zx,y", 20)
        write(iunit, '(a20)', advance = 'no') str_center("Im eps_xy,z", 20)
        write(iunit, *)

        do ik = 1, hm%kpoints%reduced%npoints
          write(iunit, '(i8)', advance = 'no') ik
          write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%weight(ik)
          do idir = 1, space%dim
            eps_mk(idir) = M_TWO * M_PI * (em_vars%alpha_be_k(magn_dir(idir, 1), magn_dir(idir, 2), idir, ik) - &
              em_vars%alpha_be_k(magn_dir(idir, 2), magn_dir(idir, 1), idir, ik)) / ions%latt%rcell_volume
          end do

          do idir = 1, space%dim
            write(iunit, '(e20.8)', advance = 'no') hm%kpoints%reduced%red_point(idir, ik)
          end do
          do idir = 1, space%dim
            write(iunit, '(e20.8)', advance = 'no') TOFLOAT(eps_mk(idir))
          end do
          do idir = 1, space%dim
            write(iunit, '(e20.8)', advance = 'no') aimag(eps_mk(idir))
          end do
          write(iunit, *)
        end do
        call io_close(iunit)
      end if

      POP_SUB(em_resp_output.out_magnetooptics)
    end subroutine out_magnetooptics

  end subroutine em_resp_output

  ! ---------------------------------------------------------
  !> Ref: David M Bishop, Rev Mod Phys 62, 343 (1990)
  !! beta // and _L are eqn (154), beta  k is eqn (155)
  !! generalized to lack of Kleinman symmetry
  subroutine out_hyperpolarizability(sb, beta, freq_factor, converged, dirname, namespace)
    type(simul_box_t),  intent(in) :: sb
    CMPLX,              intent(in) :: beta(:, :, :)
    FLOAT,              intent(in) :: freq_factor(:)
    logical,            intent(in) :: converged
    character(len=*),   intent(in) :: dirname
    type(namespace_t),  intent(in) :: namespace

    CMPLX :: bpar(1:MAX_DIM), bper(1:MAX_DIM), bk(1:MAX_DIM)
    CMPLX :: HRS_VV, HRS_HV
    integer :: ii, jj, kk, iunit

    PUSH_SUB(out_hyperpolarizability)

    ! Output first hyperpolarizability (beta)
    iunit = io_open(trim(dirname)//'/beta', namespace, action='write')

    if (.not. converged) write(iunit, '(a)') "# WARNING: not converged"

    write(iunit, '(a,3(f4.1,a),2a)', advance='no') 'First hyperpolarizability tensor: beta(', &
         freq_factor(1), ', ', freq_factor(2), ', ', freq_factor(3), ') [', &
         trim(units_abbrev(units_out%hyperpolarizability)), ']'

    write(iunit, '()')

    do ii = 1, sb%dim
      do jj = 1, sb%dim
        do kk = 1, sb%dim
          write(iunit,'(a,e20.8,e20.8)') 'beta '// &
               index2axis(ii)//index2axis(jj)//index2axis(kk)//' ', &
               units_from_atomic(units_out%hyperpolarizability, TOFLOAT( beta(ii, jj, kk))), &
               units_from_atomic(units_out%hyperpolarizability, aimag(beta(ii, jj, kk)))
        end do
      end do
    end do

    if (sb%dim == 3) then 
      bpar = M_ZERO
      bper = M_ZERO

      do ii = 1, sb%dim
        do jj = 1, sb%dim
          bpar(ii) = bpar(ii) + beta(ii, jj, jj) + beta(jj, ii, jj) + beta(jj, jj, ii)
          bper(ii) = bper(ii) + M_TWO*beta(ii, jj, jj) - M_THREE*beta(jj, ii, jj) + M_TWO*beta(jj, jj, ii)
        end do
      end do

      write(iunit, '()')

      bpar = bpar / M_FIVE
      bper = bper / M_FIVE
      bk(1:sb%dim) = M_THREE*M_HALF*(bpar(1:sb%dim) - bper(1:sb%dim))

      do ii = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta // '//index2axis(ii), &
          units_from_atomic(units_out%hyperpolarizability, TOFLOAT(bpar(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bpar(ii)))
      end do

      write(iunit, '()')

      do ii = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta _L '//index2axis(ii), &
          units_from_atomic(units_out%hyperpolarizability, TOFLOAT(bper(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bper(ii)))
      end do

      write(iunit, '()')

      do ii = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta  k '//index2axis(ii), &
          units_from_atomic(units_out%hyperpolarizability, TOFLOAT(bk(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bk(ii)))
      end do

      call calc_beta_HRS(sb, beta, HRS_VV, HRS_HV)

      write(iunit, '()')
      write(iunit, '(a)') 'beta for liquid- or gas-phase hyper-Rayleigh scattering:'
      write(iunit, '(a, 2e20.8)') 'VV polarization ', &
         units_from_atomic(units_out%hyperpolarizability, TOFLOAT(sqrt(HRS_VV))), &
         units_from_atomic(units_out%hyperpolarizability, aimag(sqrt(HRS_VV)))
      write(iunit, '(a, 2e20.8)') 'HV polarization ', &
         units_from_atomic(units_out%hyperpolarizability, TOFLOAT(sqrt(HRS_HV))), &
         units_from_atomic(units_out%hyperpolarizability, aimag(sqrt(HRS_HV)))
    end if

    call io_close(iunit)
    POP_SUB(out_hyperpolarizability)

    contains

      ! ---------------------------------------------------------
      !> calculate hyper-Rayleigh scattering hyperpolarizabilities
      !! SJ Cyvin, JE Rauch, and JC Decius, J Chem Phys 43, 4083 (1965)
      !! generalized to avoid assumption of Kleinman symmetry (permutation of indices)
      !! as in R Bersohn, Y-H Pao, and HL Frisch, J Chem Phys 45, 3184 (1966)
      subroutine calc_beta_HRS(sb, beta, HRS_VV, HRS_HV)
        type(simul_box_t), intent(in)  :: sb
        CMPLX,             intent(in)  :: beta(:, :, :)
        CMPLX,             intent(out) :: HRS_VV, HRS_HV

        CMPLX :: HRS_A, HRS_B, HRS_C, HRS_D, HRS_E
        CMPLX :: HRS_B1, HRS_B2, HRS_C1, HRS_C2, HRS_C3, HRS_D1, HRS_D2, HRS_D3, HRS_E1, HRS_E2
        integer :: ii, jj
        
        PUSH_SUB(out_hyperpolarizability.calc_beta_HRS)
        
        ! first calculate VV (vertical-vertical) polarization, FFF in Decius et al.
        HRS_A = M_ZERO
        do ii = 1, sb%dim
          HRS_A = HRS_A + beta(ii,ii,ii)**2
        end do

        HRS_B = M_ZERO
        HRS_C = M_ZERO
        do ii = 1, sb%dim
          do jj = 1, sb%dim
            if (ii /= jj) then
              HRS_B = HRS_B + beta(ii,ii,ii) * (beta(ii,jj,jj) + beta(jj,ii,jj) + beta(jj,jj,ii))
              HRS_C = HRS_C + (beta(ii,ii,jj) + beta(ii,jj,ii) + beta(jj,ii,ii))**2
            end if
          end do
        end do

        HRS_D = (beta(1,1,2) + beta(1,2,1) + beta(2,1,1)) * (beta(2,3,3) + beta(3,2,3) + beta(3,3,2)) &
              + (beta(2,2,3) + beta(2,3,2) + beta(3,2,2)) * (beta(3,1,1) + beta(1,3,1) + beta(1,1,3)) &
              + (beta(3,3,1) + beta(3,1,3) + beta(1,3,3)) * (beta(1,2,2) + beta(2,1,2) + beta(2,2,1))
    
        HRS_E = (beta(1,2,3) + beta(1,3,2) + beta(2,1,3) + beta(2,3,1) + beta(3,1,2) + beta(3,2,1))**2
    
        HRS_VV = (M_ONE / CNST(7.0))     * HRS_A &
               + (M_TWO / CNST(35.0))  * HRS_B &
               + (M_ONE / CNST(35.0))  * HRS_C &
               + (M_TWO / CNST(105.0)) * HRS_D &
               + (M_ONE / CNST(105.0)) * HRS_E

        ! now calculate HV (horizontal-vertical) polarization, FGG in Decius et al.
        HRS_B1 = M_ZERO
        HRS_B2 = M_ZERO
        HRS_C1 = M_ZERO
        HRS_C2 = M_ZERO
        HRS_C3 = M_ZERO
        do ii = 1, sb%dim
          do jj = 1, sb%dim
            if (ii /= jj) then
              HRS_B1 = HRS_B1 + beta(ii,ii,ii) * beta(ii,jj,jj)
              HRS_B2 = HRS_B2 + beta(ii,ii,ii) * (beta(jj,ii,jj) + beta(jj,jj,ii))
              HRS_C1 = HRS_C1 + (beta(ii,ii,jj) + beta(ii,jj,ii))**2
              HRS_C2 = HRS_C2 + beta(jj,ii,ii) * (beta(ii,ii,jj) + beta(ii,jj,ii))
              HRS_C3 = HRS_C3 + beta(jj,ii,ii)**2
            end if
          end do
        end do
  
        HRS_D1 = (beta(1,1,2) + beta(1,2,1) + beta(2,1,1)) * (beta(3,2,3) + beta(3,3,2)) &
               + (beta(2,2,3) + beta(2,3,2) + beta(3,2,2)) * (beta(1,3,1) + beta(1,1,3)) &
               + (beta(3,3,1) + beta(3,1,3) + beta(1,3,3)) * (beta(2,1,2) + beta(2,2,1))
        HRS_D2 = (beta(1,1,2) + beta(1,2,1)) * beta(2,3,3) &
               + (beta(2,2,3) + beta(2,3,2)) * beta(3,1,1) &
               + (beta(3,3,1) + beta(3,1,3)) * beta(1,2,2)
        HRS_D3 = beta(2,1,1) * beta(2,3,3) &
               + beta(3,2,2) * beta(3,1,1) &
               + beta(1,3,3) * beta(1,2,2)
  
        HRS_E1 = (beta(1,2,3) + beta(1,3,2))**2 &
               + (beta(2,1,3) + beta(2,3,1))**2 &
               + (beta(3,1,2) + beta(3,2,1))**2
  
        HRS_E2 = (beta(1,2,3) + beta(1,3,2)) * (beta(2,1,3) + beta(2,3,1)) &
               + (beta(2,1,3) + beta(2,3,1)) * (beta(3,1,2) + beta(3,2,1)) &
               + (beta(3,1,2) + beta(3,2,1)) * (beta(1,2,3) + beta(1,3,2))
  
        HRS_HV = (M_ONE   / CNST(35.0))  * HRS_A &
               + (M_FOUR  / CNST(105.0)) * HRS_B1 &
               - (M_ONE   / CNST(35.0))  * HRS_B2 &
               + (M_TWO   / CNST(105.0)) * HRS_C1 &
               - (M_ONE   / CNST(35.0))  * HRS_C2 &
               + (M_THREE / CNST(35.0))  * HRS_C3 &
               - (M_ONE   / CNST(105.0)) * HRS_D1 &
               - (M_ONE   / CNST(105.0)) * HRS_D2 &
               + (M_TWO   / CNST(35.0))  * HRS_D3 &
               + (M_ONE   / CNST(35.0))  * HRS_E1 &
               - (M_ONE   / CNST(105.0)) * HRS_E2
  
        POP_SUB(out_hyperpolarizability.calc_beta_HRS)
    end subroutine calc_beta_HRS

  end subroutine out_hyperpolarizability

end module em_resp_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
