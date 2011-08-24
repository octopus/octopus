!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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

module em_resp_m
  use born_charges_m
  use datasets_m
  use em_resp_calc_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use output_m
  use io_m
  use kdotp_m
  use kdotp_calc_m
  use lalg_basic_m
  use linear_response_m
  use linear_solver_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use mpi_m
  use parser_m
  use pert_m
  use profiling_m
  use restart_m
  use simul_box_m
  use smear_m
  use species_m
  use states_m
  use states_dim_m
  use sternheimer_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       em_resp_run,             &
       out_hyperpolarizability

  type em_resp_t
    type(pert_t) :: perturbation

    integer :: nsigma ! 1: consider only positive values of the frequency
                      ! 2: consider both positive and negative
    integer :: nfactor! 1: only one frequency needed
                      ! 3: three frequencies (for the hyperpolarizabilities)
    integer :: nomega ! number of frequencies to consider

    FLOAT :: eta                     ! small imaginary part to add to the frequency
    FLOAT :: freq_factor(3)
    FLOAT,      pointer :: omega(:)  ! the frequencies to consider
    type(lr_t), pointer :: lr(:,:,:) ! linear response for (gr%sb%dim, nsigma, nfactor)

    logical :: calc_hyperpol
    CMPLX   :: alpha(MAX_DIM, MAX_DIM, 3)        ! the linear polarizability
    CMPLX   :: beta (MAX_DIM, MAX_DIM, MAX_DIM)  ! first hyperpolarizability

    CMPLX   :: chi_para(MAX_DIM, MAX_DIM, 3)     ! The paramagnetic part of the susceptibility
    CMPLX   :: chi_dia (MAX_DIM, MAX_DIM, 3)     ! The diamagnetic  part of the susceptibility

    logical :: ok(1:3)                           ! whether calculation is converged
    logical :: force_no_kdotp                    ! whether to use kdotp run for periodic system

    logical :: calc_rotatory                     ! whether to calculate rotatory response
    logical :: calc_Born                         ! whether to calculate Born effective charges
    type(Born_charges_t) :: Born_charges(3)      ! one set for each frequency factor
    logical :: occ_response                      ! whether to calculate full response in Sternheimer eqn.
    logical :: wfns_from_scratch                 ! whether to ignore restart LR wfns and initialize to zero
    
  end type em_resp_t

contains

  ! ---------------------------------------------------------
  subroutine em_resp_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(em_resp_t)         :: em_vars
    type(sternheimer_t)     :: sh, sh_kdotp, sh2
    type(lr_t)              :: kdotp_lr(MAX_DIM, 1)
    type(lr_t), allocatable :: kdotp_em_lr2(:, :, :, :)
    type(pert_t)            :: pert_kdotp, pert2_none

    integer :: sigma, sigma_alt, ndim, idir, idir2, ierr, iomega, ifactor
    character(len=100) :: dirname_restart, dirname_output, str_tmp
    logical :: complex_response, have_to_calculate, use_kdotp, opp_freq, exact_freq

    FLOAT :: closest_omega, last_omega

    PUSH_SUB(em_resp_run)

    gr => sys%gr
    ndim = sys%gr%sb%dim

    call parse_input()

    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC .and. &
      any(abs(em_vars%omega(1:em_vars%nomega)) > M_EPSILON)) then
      call messages_not_implemented('Dynamical magnetic response')
    endif

    complex_response = (em_vars%eta > M_EPSILON) .or. states_are_complex(sys%st)
    call restart_look_and_read(sys%st, sys%gr, sys%geo, is_complex = complex_response)

    if (states_are_real(sys%st)) then
      message(1) = 'Info: SCF using real wavefunctions.'
    else
      message(1) = 'Info: SCF using complex wavefunctions.'
    end if
    call messages_info(1)

    use_kdotp = simul_box_is_periodic(gr%sb) .and. .not. em_vars%force_no_kdotp

    if(use_kdotp .and. pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
      call messages_not_implemented('Magnetic perturbation in periodic system')
    endif

    ! read kdotp wavefunctions if necessary
    if (use_kdotp) then
      message(1) = "Reading kdotp wavefunctions for periodic directions."
      call messages_info(1)

      do idir = 1, gr%sb%periodic_dim
        call lr_init(kdotp_lr(idir, 1))
        call lr_allocate(kdotp_lr(idir, 1), sys%st, sys%gr%mesh)

        ! load wavefunctions
        str_tmp = kdotp_wfs_tag(idir)
        write(dirname_restart,'(2a)') KDOTP_DIR, trim(wfs_tag_sigma(str_tmp, 1))
        ! 1 is the sigma index which is used in em_resp
        call restart_read(trim(tmpdir)//dirname_restart, sys%st, sys%gr, sys%geo, &
          ierr, lr=kdotp_lr(idir, 1))

        if(ierr .ne. 0) then
          message(1) = "Could not load kdotp wavefunctions from '"//trim(tmpdir)//trim(dirname_restart)//"'"
          message(2) = "Previous kdotp calculation required."
          call messages_fatal(2)
        end if
      end do
    endif

    em_vars%nfactor = 1
    if(em_vars%calc_hyperpol) em_vars%nfactor = 3

    ! in effect, nsigma = 1 only if hyperpol not being calculated, and the only frequency is zero
    if(em_vars%calc_hyperpol .or. any(abs(em_vars%omega(1:em_vars%nomega)) > M_EPSILON)) then
      em_vars%nsigma = 2
      ! positive and negative values of the frequency must be considered
    else
      em_vars%nsigma = 1
      ! only considering positive values
    endif

    if(em_vars%calc_hyperpol .and. use_kdotp) then
      call pert_init(pert_kdotp, PERTURBATION_KDOTP, sys%gr, sys%geo)
      call pert_init(pert2_none, PERTURBATION_NONE,  sys%gr, sys%geo)
      call messages_experimental("Second-order Sternheimer equation")
      call pert_setup_dir(pert2_none, 1)  ! direction is irrelevant
      SAFE_ALLOCATE(kdotp_em_lr2(1:gr%sb%periodic_dim, 1:gr%sb%dim, 1:em_vars%nsigma, 1:em_vars%nfactor))
      do ifactor = 1, em_vars%nfactor
        do sigma = 1, em_vars%nsigma
          do idir = 1, gr%sb%periodic_dim
            do idir2 = 1, gr%sb%dim
              call lr_init(kdotp_em_lr2(idir, idir2, sigma, ifactor))
              call lr_allocate(kdotp_em_lr2(idir, idir2, sigma, ifactor), sys%st, sys%gr%mesh)
            enddo
          enddo
        enddo
      enddo
      call sternheimer_init(sh2, sys, hm, "EM", complex_response, set_ham_var = 0, set_last_occ_response = .true.)
      call sternheimer_init(sh_kdotp, sys, hm, "EM", complex_response, set_ham_var = 0, &
        set_last_occ_response = .true.)
      em_vars%occ_response = .true.
    endif

    SAFE_ALLOCATE(em_vars%lr(1:gr%sb%dim, 1:em_vars%nsigma, 1:em_vars%nfactor))
    do ifactor = 1, em_vars%nfactor
      call Born_charges_init(em_vars%Born_charges(ifactor), sys%geo, sys%st, gr%sb%dim)
    enddo

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call messages_info(1)
    call system_h_setup(sys, hm)

    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC &
      .and. sys%st%d%nspin == 1 .and. states_are_real(sys%st)) then
      ! first-order response is zero if there is time-reversal symmetry. F Mauri and SG Louie, PRL 76, 4246 (1996)
      call sternheimer_init(sh, sys, hm, "EM", complex_response, set_ham_var = 0, set_last_occ_response = em_vars%occ_response)
      ! set HamiltonianVariation to V_ext_only, in magnetic case
    else
      call sternheimer_init(sh, sys, hm, "EM", complex_response, set_last_occ_response = em_vars%occ_response)
      ! otherwise, use default, which is hartree + fxc
    endif

    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(trim(tmpdir)//EM_RESP_DIR, is_tmp=.true.) ! restart
      call info()
      call io_mkdir(EM_RESP_DIR) ! output
    endif

    do ifactor = 1, em_vars%nfactor
      do idir = 1, sys%gr%sb%dim
        do sigma = 1, em_vars%nsigma
          call lr_init(em_vars%lr(idir, sigma, ifactor))
          call lr_allocate(em_vars%lr(idir, sigma, ifactor), sys%st, sys%gr%mesh)
        end do
      end do
    end do

    last_omega = M_HUGE
    do iomega = 1, em_vars%nomega

      em_vars%ok(1:3) = .true.

      do ifactor = 1, em_vars%nfactor
        do idir = 1, sys%gr%sb%dim

          ierr = 0

          have_to_calculate = .true.
          opp_freq = .false.

          ! if this frequency is zero and this is not the first
          ! iteration we do not have to do anything
          if(iomega > 1 .and. em_vars%freq_factor(ifactor) == M_ZERO) have_to_calculate = .false. 

          if(ifactor > 1) then 

            ! if this frequency is the same as the previous one, just copy it
            if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor - 1) * em_vars%omega(iomega) &
                                            - em_vars%freq_factor(ifactor) * em_vars%omega(iomega)) < M_EPSILON ) then

              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor - 1), em_vars%lr(idir, 1, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor - 1), em_vars%lr(idir, 2, ifactor))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, gr%sb%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 1, ifactor))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 2, ifactor))
                enddo
              endif

              have_to_calculate = .false.

            end if

            ! if this frequency is minus the previous one, copy it inverted
            if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor - 1) * em_vars%omega(iomega) &
                                            + em_vars%freq_factor(ifactor) * em_vars%omega(iomega)) < M_EPSILON ) then 

              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor - 1), em_vars%lr(idir, 2, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor - 1), em_vars%lr(idir, 1, ifactor))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, gr%sb%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 2, ifactor))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, ifactor - 1), &
                    kdotp_em_lr2(idir, idir2, 1, ifactor))
                enddo
              endif

              have_to_calculate = .false.

            end if

          end if

          if(iomega > 1 .and. ifactor == 1) then 

            ! if this frequency is the same as the previous one, just copy it
            if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor) * em_vars%omega(iomega) - last_omega) < M_EPSILON ) then

              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, em_vars%nfactor), em_vars%lr(idir, 1, 1))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, em_vars%nfactor), em_vars%lr(idir, 2, 1))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, gr%sb%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 1, 1))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 2, 1))
                enddo
              endif

              have_to_calculate = .false.

            end if

            ! if this frequency is minus the previous one, copy it inverted
            if( have_to_calculate .and. abs(em_vars%freq_factor(ifactor) * em_vars%omega(iomega) + last_omega) < M_EPSILON ) then

              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, em_vars%nfactor), em_vars%lr(idir, 2, 1))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, em_vars%nfactor), em_vars%lr(idir, 1, 1))

              if(em_vars%calc_hyperpol .and. use_kdotp) then
                do idir2 = 1, gr%sb%periodic_dim
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 1, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 2, 1))
                  call lr_copy(sys%st, sys%gr%mesh, kdotp_em_lr2(idir, idir2, 2, em_vars%nfactor), &
                    kdotp_em_lr2(idir, idir2, 1, 1))
                enddo
              endif

              have_to_calculate = .false.

            end if

          end if

          if(have_to_calculate) then 

            str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(ifactor) * em_vars%omega(iomega)))
            write(message(1), '(5a)') 'Info: Calculating response for the ', index2axis(idir), &
              '-direction and frequency ', trim(str_tmp), '.'
            call messages_info(1)

            exact_freq = .false.
            if(.not. fromscratch) then 

              ! try to load wavefunctions, if first frequency; otherwise will already be initialized
              if(iomega == 1 .and. .not. em_vars%wfns_from_scratch) then
                do sigma = 1, em_vars%nsigma
                  if(sigma == 2 .and. abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)) < M_EPSILON) then
                    if(states_are_real(sys%st)) then
                      em_vars%lr(idir, 2, ifactor)%ddl_psi = em_vars%lr(idir, 1, ifactor)%ddl_psi
                    else
                      em_vars%lr(idir, 2, ifactor)%zdl_psi = em_vars%lr(idir, 1, ifactor)%zdl_psi
                    endif

                    if(em_vars%calc_hyperpol .and. use_kdotp) then
                      do idir2 = 1, gr%sb%periodic_dim
                        if(states_are_real(sys%st)) then
                          kdotp_em_lr2(idir2, idir, 2, ifactor)%ddl_psi = kdotp_em_lr2(idir2, idir, 1, ifactor)%ddl_psi
                        else
                          kdotp_em_lr2(idir2, idir, 2, ifactor)%zdl_psi = kdotp_em_lr2(idir2, idir, 1, ifactor)%zdl_psi
                        endif
                      enddo
                    endif
                  else
                    sigma_alt = sigma
                    if(em_vars%freq_factor(ifactor) * em_vars%omega(iomega) < -M_EPSILON .and. em_vars%nsigma == 2) &
                      sigma_alt = swap_sigma(sigma)
                    
                    str_tmp = em_wfs_tag(idir, ifactor)
                    write(dirname_restart,'(2a)') EM_RESP_DIR, trim(wfs_tag_sigma(str_tmp, sigma))
                    call restart_read(trim(tmpdir)//dirname_restart, sys%st, sys%gr, sys%geo, &
                      ierr, lr=em_vars%lr(idir, sigma_alt, ifactor))
                    
                    if(ierr .ne. 0) then
                      message(1) = "Initializing to zero, could not load response wavefunctions from '" &
                        //trim(tmpdir)//trim(dirname_restart)//"'"
                      call messages_warning(1)
                    end if

                    if(em_vars%calc_hyperpol .and. use_kdotp) then
                      do idir2 = 1, gr%sb%periodic_dim
                        str_tmp = em_wfs_tag(idir, ifactor, idir2)
                        write(dirname_restart,'(2a)') EM_RESP_DIR, trim(wfs_tag_sigma(str_tmp, sigma))
                        call restart_read(trim(tmpdir)//dirname_restart, sys%st, sys%gr, sys%geo, &
                          ierr, lr=kdotp_em_lr2(idir2, idir, sigma_alt, ifactor))
                        
                        if(ierr .ne. 0) then
                          message(1) = "Initializing to zero, could not load second-order response wavefunctions from '" &
                            //trim(tmpdir)//trim(dirname_restart)//"'"
                          call messages_warning(1)
                        end if
                      enddo
                    endif
                  endif
                end do
              endif

              ! if opposite sign from last frequency, swap signs to get a closer frequency
              if(iomega > 1 .and. em_vars%nsigma == 2) then
                if(em_vars%omega(iomega - 1) * em_vars%omega(iomega) < M_ZERO) then
                  if(states_are_complex(sys%st)) then
                    call zlr_swap_sigma(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
                  else
                    call dlr_swap_sigma(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
                  endif
                endif
              end if
                
              !search for the density of the closest frequency, including negative
              closest_omega = em_vars%freq_factor(ifactor) * em_vars%omega(iomega)
              call oct_search_file_lr(closest_omega, idir, ierr, trim(tmpdir)//EM_RESP_DIR)
              sigma_alt = 1
              if(closest_omega * em_vars%freq_factor(ifactor) * em_vars%omega(iomega) < M_ZERO) opp_freq = .true.
              if(opp_freq .and. em_vars%nsigma == 2) sigma_alt = 2

              !attempt to read 
              if(ierr == 0) then 
                if(states_are_complex(sys%st)) then 
                  call zrestart_read_lr_rho(em_vars%lr(idir, sigma_alt, ifactor), sys%gr, sys%st%d%nspin, &
                    EM_RESP_DIR, em_rho_tag(closest_omega, idir), ierr)
                else 
                  call drestart_read_lr_rho(em_vars%lr(idir, sigma_alt, ifactor), sys%gr, sys%st%d%nspin, &
                    EM_RESP_DIR, em_rho_tag(closest_omega, idir), ierr)
                end if

                if(ierr == 0 .and. &
                  abs(abs(closest_omega) - abs(em_vars%freq_factor(ifactor) * em_vars%omega(iomega))) <= CNST(1e-4)) then
                  ! the frequencies are written to four decimals in the restart directory, so we cannot expect higher precision
                  exact_freq = .true.
                endif
              end if

              if(ierr == 0 .and. em_vars%nsigma == 2) then 
                sigma_alt = 1
                if(opp_freq) sigma_alt = 2

                if(states_are_complex(sys%st)) then 
                  em_vars%lr(idir, swap_sigma(sigma_alt), ifactor)%zdl_rho = conjg(em_vars%lr(idir, sigma_alt, ifactor)%zdl_rho)
                else 
                  em_vars%lr(idir, swap_sigma(sigma_alt), ifactor)%ddl_rho =       em_vars%lr(idir, sigma_alt, ifactor)%ddl_rho
                end if
              end if

            end if ! .not. fromscratch
            
            call pert_setup_dir(em_vars%perturbation, idir)

            if(use_kdotp .and. idir <= gr%sb%periodic_dim) then
              if (states_are_complex(sys%st)) then
                call zsternheimer_set_rhs(sh, kdotp_lr(idir, 1)%zdl_psi)
              else
                call dsternheimer_set_rhs(sh, kdotp_lr(idir, 1)%ddl_psi)
              endif
            end if

            ! if the frequency is zero, we do not need to calculate both responses
            if(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)) < M_EPSILON .and. em_vars%nsigma == 2) then

              if (states_are_complex(sys%st)) then
                call zsternheimer_solve(sh, sys, hm, em_vars%lr(idir, 1:1, ifactor), 1, &
                  em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, &
                  em_vars%perturbation, EM_RESP_DIR, &
                  em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
                  em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0), have_exact_freq = exact_freq)

                em_vars%lr(idir, 2, ifactor)%zdl_psi = em_vars%lr(idir, 1, ifactor)%zdl_psi
                em_vars%lr(idir, 2, ifactor)%zdl_rho = conjg(em_vars%lr(idir, 1, ifactor)%zdl_rho)
              else
                call dsternheimer_solve(sh, sys, hm, em_vars%lr(idir, 1:1, ifactor), 1, &
                  em_vars%freq_factor(ifactor)*em_vars%omega(iomega), &
                  em_vars%perturbation, EM_RESP_DIR, &
                  em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
                  em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0), have_exact_freq = exact_freq)

                em_vars%lr(idir, 2, ifactor)%ddl_psi = em_vars%lr(idir, 1, ifactor)%ddl_psi
                em_vars%lr(idir, 2, ifactor)%ddl_rho = em_vars%lr(idir, 1, ifactor)%ddl_rho
              end if

            else

              if (states_are_complex(sys%st)) then
                call zsternheimer_solve(sh, sys, hm, em_vars%lr(idir, 1:em_vars%nsigma, ifactor), em_vars%nsigma, &
                  em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, &
                  em_vars%perturbation, EM_RESP_DIR, &
                  em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
                  em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0), have_exact_freq = exact_freq)
              else
                call dsternheimer_solve(sh, sys, hm, em_vars%lr(idir, 1:em_vars%nsigma, ifactor), em_vars%nsigma, &
                  em_vars%freq_factor(ifactor)*em_vars%omega(iomega), &
                  em_vars%perturbation, EM_RESP_DIR, &
                  em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
                  em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0), have_exact_freq = exact_freq)
              end if

            end if

            if(use_kdotp) then
              call sternheimer_unset_rhs(sh)
            end if

            em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)

            if(em_vars%calc_hyperpol .and. use_kdotp) then
              do idir2 = 1, gr%sb%periodic_dim
                write(message(1), '(a,a,a)') 'Info: Calculating kdotp response in ', index2axis(idir2), '-direction.'
                call messages_info(1)
                call pert_setup_dir(pert_kdotp, idir2)
                
                ! need nsigma
                ! need to give a proper name to the restart files

                ! if the frequency is zero, we do not need to calculate both responses
                if(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)) < M_EPSILON .and. em_vars%nsigma == 2) then
                  if (states_are_complex(sys%st)) then
                    call zsternheimer_solve_order2(sh, sh_kdotp, sh2, sys, hm, em_vars%lr(idir, 1:1, ifactor), &
                      kdotp_lr(idir2, 1:1), 1, &
                      em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, M_z0, &
                      em_vars%perturbation, pert_kdotp, kdotp_em_lr2(idir2, idir, 1:1, ifactor), &
                      pert2_none, EM_RESP_DIR, &
                      "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true.)
                    kdotp_em_lr2(idir2, idir, 2, ifactor)%zdl_psi = kdotp_em_lr2(idir2, idir, 1, ifactor)%zdl_psi
                  else
                    call dsternheimer_solve_order2(sh, sh_kdotp, sh2, sys, hm, em_vars%lr(idir, 1:1, ifactor), &
                      kdotp_lr(idir2, 1:1), 1, &
                      em_vars%freq_factor(ifactor)*em_vars%omega(iomega), M_ZERO, &
                      em_vars%perturbation, pert_kdotp, kdotp_em_lr2(idir2, idir, 1:1, ifactor), &
                      pert2_none, EM_RESP_DIR, &
                      "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true.)
                    kdotp_em_lr2(idir2, idir, 2, ifactor)%ddl_psi = kdotp_em_lr2(idir2, idir, 1, ifactor)%ddl_psi
                  end if
                else
                  if (states_are_complex(sys%st)) then
                    call zsternheimer_solve_order2(sh, sh_kdotp, sh2, sys, hm, em_vars%lr(idir, 1:em_vars%nsigma, ifactor), &
                      kdotp_lr(idir2, 1:1), em_vars%nsigma, &
                      em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, M_z0, &
                      em_vars%perturbation, pert_kdotp, kdotp_em_lr2(idir2, idir, 1:em_vars%nsigma, ifactor), &
                      pert2_none, EM_RESP_DIR, &
                      "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true.)
                  else
                    call dsternheimer_solve_order2(sh, sh_kdotp, sh2, sys, hm, em_vars%lr(idir, 1:em_vars%nsigma, ifactor), &
                      kdotp_lr(idir2, 1:1), em_vars%nsigma, &
                      em_vars%freq_factor(ifactor)*em_vars%omega(iomega), M_ZERO, &
                      em_vars%perturbation, pert_kdotp, kdotp_em_lr2(idir2, idir, 1:em_vars%nsigma, ifactor), &
                      pert2_none, EM_RESP_DIR, &
                      "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true.)
                  end if
                endif

                em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
              enddo
              write(message(1), '(a)') ''
              call messages_info(1)
            endif
          end if ! have_to_calculate

        end do ! idir

        if(.not. have_to_calculate) cycle

        if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
    
          ! calculate polarizability
          message(1) = "Info: Calculating polarizabilities."
          call messages_info(1)
    
          if(use_kdotp) then
            if(states_are_complex(sys%st)) then
              call zcalc_polarizability_periodic(sys, em_vars%lr(:, :, ifactor), kdotp_lr(:, 1), &
                em_vars%nsigma, em_vars%alpha(:, :, ifactor))
            else
              call dcalc_polarizability_periodic(sys, em_vars%lr(:, :, ifactor), kdotp_lr(:, 1), &
                em_vars%nsigma, em_vars%alpha(:, :, ifactor))
            endif
          endif

          if(states_are_complex(sys%st)) then
            call zcalc_polarizability_finite(sys, hm, em_vars%lr(:, :, ifactor), em_vars%nsigma, &
              em_vars%perturbation, em_vars%alpha(:, :, ifactor), doalldirs = .not. use_kdotp)
          else
            call dcalc_polarizability_finite(sys, hm, em_vars%lr(:, :, ifactor), em_vars%nsigma, &
              em_vars%perturbation, em_vars%alpha(:, :, ifactor), doalldirs = .not. use_kdotp)
          end if
    
          if(em_vars%calc_Born) then
            ! calculate Born effective charges
            message(1) = "Info: Calculating (frequency-dependent) Born effective charges."
            call messages_info(1)
    
            do idir = 1, sys%gr%sb%dim
              ! time = M_ZERO
              if(states_are_complex(sys%st)) then
                if(em_vars%nsigma == 2) then
                  call zforces_born_charges(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                    lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 2, ifactor), &
                    lr_dir = idir, Born_charges = em_vars%Born_charges(ifactor))
                else
                  call zforces_born_charges(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                    lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 1, ifactor), &
                    lr_dir = idir, Born_charges = em_vars%Born_charges(ifactor))
                endif
              else
                if(em_vars%nsigma == 2) then
                  call dforces_born_charges(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                    lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 2, ifactor), &
                    lr_dir = idir, Born_charges = em_vars%Born_charges(ifactor))
                else
                  call dforces_born_charges(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                    lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 1, ifactor), &
                    lr_dir = idir, Born_charges = em_vars%Born_charges(ifactor))
                endif
              endif
            enddo
          endif

        else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
          message(1) = "Info: Calculating magnetic susceptibilities."
          call messages_info(1)
    
          if(states_are_complex(sys%st)) then 
            call zlr_calc_susceptibility(sys, hm, em_vars%lr(:,:, ifactor), em_vars%nsigma, em_vars%perturbation, &
              em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          else
            call dlr_calc_susceptibility(sys, hm, em_vars%lr(:,:, ifactor), em_vars%nsigma, em_vars%perturbation, &
              em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          end if
        end if

        call em_resp_output(sys%st, sys%gr, hm, sys%geo, sys%outp, em_vars, iomega, ifactor)

      end do ! ifactor

      ! calculate hyperpolarizability
      if(em_vars%calc_hyperpol) then
        write(message(1), '(a)') 'Info: Calculating hyperpolarizabilities.'
        call messages_info(1)

        if(use_kdotp) then
          if(states_are_complex(sys%st)) then
            call zpost_orthogonalize(sys, em_vars%nfactor, em_vars%nsigma, em_vars%freq_factor(:), &
              em_vars%omega(iomega) + M_zI * em_vars%eta, kdotp_lr(:, 1), em_vars%lr, kdotp_em_lr2)
            call zlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, &
              kdotp_lr = kdotp_lr(:, 1), kdotp_em_lr = kdotp_em_lr2, occ_response = .false.)
          else
            call dpost_orthogonalize(sys, em_vars%nfactor, em_vars%nsigma, em_vars%freq_factor(:), &
              em_vars%omega(iomega), kdotp_lr(:, 1), em_vars%lr, kdotp_em_lr2)
            call dlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, &
              kdotp_lr = kdotp_lr(:, 1), kdotp_em_lr = kdotp_em_lr2, occ_response = .false.)
          end if
        else
          if(states_are_complex(sys%st)) then
            call zlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, occ_response = em_vars%occ_response)
          else
            call dlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, occ_response = em_vars%occ_response)
          end if
        endif

        str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(1)*em_vars%omega(iomega)))
        write(dirname_output, '(a, a)') EM_RESP_DIR//'freq_', trim(str_tmp)
        call io_mkdir(trim(dirname_output))
        call out_hyperpolarizability(gr%sb, em_vars%beta, em_vars%freq_factor(1:3), em_vars%ok(1), dirname_output)
      end if

      last_omega = em_vars%freq_factor(em_vars%nfactor) * em_vars%omega(iomega)

    end do ! iomega

    do idir = 1, ndim
      do sigma = 1, em_vars%nsigma
        do ifactor = 1, em_vars%nfactor
          call lr_dealloc(em_vars%lr(idir, sigma, ifactor))
        end do
      end do
    end do

    call sternheimer_end(sh)
    call pert_end(em_vars%perturbation)

    SAFE_DEALLOCATE_P(em_vars%omega)
    SAFE_DEALLOCATE_P(em_vars%lr)
    do ifactor = 1, em_vars%nfactor
      call Born_charges_end(em_vars%Born_charges(ifactor))
    enddo
    call states_deallocate_wfns(sys%st)

    POP_SUB(em_resp_run)

  contains

    ! ---------------------------------------------------------
    subroutine parse_input()

      type(block_t) :: blk
      integer :: nrow, irow, nfreqs_in_row, ifreq, istep, perturb_type
      FLOAT :: omega_ini, omega_fin, domega
      logical :: freq_sort  

      PUSH_SUB(em_resp_run.parse_input)

      call messages_obsolete_variable('PolFreqs               ', 'EMFreqs             ')
      call messages_obsolete_variable('PolHyper               ', 'EMHyperpol          ')
      call messages_obsolete_variable('PolEta                 ', 'EMEta               ')
      call messages_obsolete_variable('PolConvAbsDens         ', 'LRConvAbsDens       ')
      call messages_obsolete_variable('PolHamiltonianVariation', 'HamiltonianVariation')

      !%Variable EMFreqs
      !%Type block
      !%Section Linear Response::Polarizabilities
      !%Description
      !% This block defines for which frequencies the polarizabilities
      !% will be calculated. If it is not present, the static (omega = 0) response
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

      if (parse_block(datasets_check('EMFreqs'), blk) == 0) then 

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
        call parse_logical(datasets_check('EMFreqsSort'), .true., freq_sort)

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

      call parse_float(datasets_check('EMEta'), M_ZERO, em_vars%eta, units_inp%energy)
      if(em_vars%eta < -M_EPSILON) then
        message(1) = "EMEta cannot be negative."
        call messages_fatal(1)
      endif

      ! reset the values of these variables
      em_vars%calc_hyperpol = .false.
      em_vars%freq_factor(1:3) = M_ONE

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
      call parse_integer(datasets_check('EMPerturbationType'), PERTURBATION_ELECTRIC, perturb_type)

      call pert_init(em_vars%perturbation, perturb_type, sys%gr, sys%geo)

      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        !%Variable EMHyperpol
        !%Type block
        !%Section Linear Response::Polarizabilities
        !%Description
        !% This block describes the multiples of the frequency used for
        !% the dynamic hyperpolarizability. The results are written to the
        !% file <tt>beta</tt> in the directory for the first multiple.
        !% There must be three factors, summing to zero. For example,
        !% for second-harmonic generation, you could use
        !% <tt>1 | 1 | -2</tt>.
        !%End

        if (parse_block(datasets_check('EMHyperpol'), blk) == 0) then 
          call parse_block_float(blk, 0, 0, em_vars%freq_factor(1))
          call parse_block_float(blk, 0, 1, em_vars%freq_factor(2))
          call parse_block_float(blk, 0, 2, em_vars%freq_factor(3))

          call parse_block_end(blk)

          if(abs(sum(em_vars%freq_factor(1:3))) > M_EPSILON) then
            message(1) = "Frequency factors specified by EMHyperpol must sum to zero."
            call messages_fatal(1)
          endif

          em_vars%calc_hyperpol = .true.
        end if
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

      call parse_logical(datasets_check('EMForceNoKdotP'), .false., em_vars%force_no_kdotp)

      !%Variable EMCalcBornCharges
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Calculate linear-response Born effective charges from electric perturbation (experimental).
      !%End

      call parse_logical(datasets_check('EMCalcBornCharges'), .false., em_vars%calc_Born)
      if (em_vars%calc_Born) call messages_experimental("Calculation of Born effective charges")

      !%Variable EMCalcRotatoryResponse
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Calculate circular-dichroism spectrum from electric perturbation,
      !% and write to file <tt>rotatory_strength</tt>.
      !%End

      call parse_logical(datasets_check('EMCalcRotatoryResponse'), .false., em_vars%calc_rotatory)

      !%Variable EMOccupiedResponse
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Solve for full response without projector into unoccupied subspace.
      !% Not possible if there are partial occupations.
      !%End

      call parse_logical(datasets_check('EMOccupiedResponse'), .false., em_vars%occ_response)
      if(em_vars%occ_response .and. .not. (smear_is_semiconducting(sys%st%smear) .or. sys%st%smear%method == SMEAR_FIXED_OCC)) then
        message(1) = "EMOccupiedResponse cannot be used if there are partial occupations."
        call messages_fatal(1)
      endif

      !%Variable EMWavefunctionsFromScratch
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Do not use saved linear-response wavefunctions from a previous run as starting guess.
      !% Instead initialize to zero as in <tt>FromScratch</tt>, but restart densities will still
      !% be used. Restart wavefunctions from a very different frequency can hinder convergence.
      !%End

      call parse_logical(datasets_check('EMWavefunctionsFromScratch'), .false., em_vars%wfns_from_scratch)

      POP_SUB(em_resp_run.parse_input)

    end subroutine parse_input


    ! ---------------------------------------------------------
    subroutine info()

      PUSH_SUB(em_resp_run.info)

      call pert_info(em_vars%perturbation, stdout)
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

      POP_SUB(em_resp_run.info)

    end subroutine info

  end subroutine em_resp_run


  ! ---------------------------------------------------------
  subroutine em_resp_output(st, gr, hm, geo, outp, em_vars, iomega, ifactor)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(geometry_t),     intent(inout) :: geo
    type(output_t),       intent(in)    :: outp
    type(em_resp_t),      intent(inout) :: em_vars
    integer,              intent(in)    :: iomega
    integer,              intent(in)    :: ifactor
    
    integer :: iunit
    character(len=80) :: dirname, str_tmp

    PUSH_SUB(em_resp_output)

    str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(ifactor)*em_vars%omega(iomega)))
    write(dirname, '(a, a)') EM_RESP_DIR//'freq_', trim(str_tmp)
    call io_mkdir(trim(dirname))

    call write_eta()

    if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
      call out_polarizability()
      if(em_vars%calc_Born) then
        call out_Born_charges(em_vars%Born_charges(ifactor), geo, gr%sb%dim, dirname, &
          write_real = em_vars%eta < M_EPSILON)
      endif
    else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
      call out_susceptibility()
    end if

!      call out_projections()
!      This routine does not give any useful information.

    call out_wfn_and_densities()

    if(gr%sb%periodic_dim .eq. gr%sb%dim) then
      call out_dielectric_constant()
    endif

    if((.not. simul_box_is_periodic(gr%sb) .or. em_vars%force_no_kdotp) .and. em_vars%calc_rotatory) then
      call out_circular_dichroism()
    endif

    POP_SUB(em_resp_output)

  contains


    subroutine write_eta()

      PUSH_SUB(em_resp_output.write_eta)

      iunit = io_open(trim(dirname)//'/eta', action='write')

      write(iunit, '(3a)') 'Imaginary part of frequency [', trim(units_abbrev(units_out%energy)), ']'
      write(iunit, '(f20.6)') units_from_atomic(units_out%energy, em_vars%eta)

      call io_close(iunit)

      POP_SUB(em_resp_output.write_eta)
    end subroutine write_eta


    ! ---------------------------------------------------------
    ! Note: this should be in spectrum.F90
    subroutine cross_section_header(out_file)
      integer, intent(in) :: out_file

      character(len=80) :: header_string
      integer :: ii, idir, kdir

      PUSH_SUB(em_resp_output.cross_section_header)

      !this header is the same as spectrum.F90
      write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
      write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
      write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)

      do idir = 1, gr%sb%dim
        do kdir = 1, gr%sb%dim
          write(header_string,'(a6,i1,a1,i1,a1)') 'sigma(', idir, ',', kdir, ')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do

      write(out_file, *)
      write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20)
      do ii = 1, 2 + gr%sb%dim**2
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
  
      iunit = io_open(trim(dirname)//'/alpha', action='write')
  
      if (.not. em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
      write(iunit, '(3a)') '# Polarizability tensor [', trim(units_abbrev(units_out%polarizability)), ']'
      call output_tensor(iunit, TOFLOAT(em_vars%alpha(:, :, ifactor)), &
        gr%sb%dim, units_out%polarizability)
  
      call io_close(iunit)
  
      ! CROSS SECTION (THE IMAGINARY PART OF POLARIZABILITY)
      if(em_vars%eta > M_EPSILON) then 
        cross(1:gr%sb%dim, 1:gr%sb%dim) = aimag(em_vars%alpha(1:gr%sb%dim, 1:gr%sb%dim, ifactor)) * &
          em_vars%freq_factor(ifactor) * em_vars%omega(iomega) * (M_FOUR * M_PI / P_c)

        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            cross(idir, idir2) = units_from_atomic(units_out%length**2, cross(idir, idir2))
          enddo
        enddo

        iunit = io_open(trim(dirname)//'/cross_section', action='write')
        if (.not. em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
        crossp(1:gr%sb%dim, 1:gr%sb%dim) = matmul(cross(1:gr%sb%dim, 1:gr%sb%dim), cross(1:gr%sb%dim, 1:gr%sb%dim))

        cross_sum = M_ZERO
        crossp_sum = M_ZERO
        do idir = 1, gr%sb%dim
          cross_sum = cross_sum + cross(idir, idir)
          crossp_sum = crossp_sum + crossp(idir, idir)
        enddo

        anisotropy = crossp_sum - M_THIRD * cross_sum**2
            
        call cross_section_header(iunit)
        write(iunit,'(3e20.8)', advance = 'no') &
          units_from_atomic(units_out%energy, em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), &
          cross_sum * M_THIRD, sqrt(max(anisotropy, M_ZERO))
        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            write(iunit,'(e20.8)', advance = 'no') cross(idir, idir2)
          enddo
        enddo
        write(iunit,'(a)', advance = 'yes')
  
        call io_close(iunit)
      end if
      
      POP_SUB(em_resp_output.out_polarizability)
    end subroutine out_polarizability


    ! ---------------------------------------------------------
    ! epsilon = 1 + 4 * pi * alpha/volume
    subroutine out_dielectric_constant()
      CMPLX epsilon(MAX_DIM, MAX_DIM) 
      integer idir

      PUSH_SUB(em_resp_output.out_dielectric_constant)
  
      iunit = io_open(trim(dirname)//'/epsilon', action='write')
      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
  
      epsilon(1:gr%sb%dim, 1:gr%sb%dim) = &
        4 * M_PI * em_vars%alpha(1:gr%sb%dim, 1:gr%sb%dim, ifactor) / gr%sb%rcell_volume
      do idir = 1, gr%sb%dim
        epsilon(idir, idir) = epsilon(idir, idir) + M_ONE
      enddo

      write(iunit, '(a)') '# Real part of dielectric constant'
      call output_tensor(iunit, real(epsilon(1:gr%sb%dim, 1:gr%mesh%sb%dim)), gr%sb%dim, unit_one)
      write(iunit, '(a)')
      write(iunit, '(a)') '# Imaginary part of dielectric constant'
      call output_tensor(iunit, aimag(epsilon(1:gr%sb%dim, 1:gr%mesh%sb%dim)), gr%sb%dim, unit_one)
  
      call io_close(iunit)
      POP_SUB(em_resp_output.out_dielectric_constant)
    end subroutine out_dielectric_constant


    ! ---------------------------------------------------------
    subroutine out_susceptibility()

      PUSH_SUB(em_resp_output.out_susceptibility)

      iunit = io_open(trim(dirname)//'/susceptibility', action='write')

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor)), gr%sb%dim, unit_ppm)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :, ifactor)), gr%sb%dim, unit_ppm)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm a.u.]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor) + em_vars%chi_dia(:,:, ifactor)), &
        gr%sb%dim, unit_ppm)
      write(iunit, '(1x)')

      write(iunit, '(a)') hyphens

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor)), gr%sb%dim, unit_susc_ppm_cgs)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :, ifactor)), gr%sb%dim, unit_susc_ppm_cgs)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm cgs / mol]'
      call output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor) + em_vars%chi_dia(:,:, ifactor)), &
           gr%sb%dim, unit_susc_ppm_cgs)
      write(iunit, '(1x)')

      call io_close(iunit)      
      POP_SUB(em_resp_output.out_susceptibility)
    end subroutine out_susceptibility

    ! ---------------------------------------------------------
    subroutine out_projections()
      CMPLX   :: proj
      integer :: ist, ivar, ik, idir, sigma
      character(len=80) :: fname

      PUSH_SUB(em_resp_output.out_projections)

      do ik = st%d%kpt%start, st%d%kpt%end
        do idir = 1, gr%sb%dim

          write(fname, '(2a,i1,2a)') trim(dirname), '/projection-k', ik, '-', index2axis(idir)
          iunit = io_open(trim(fname), action='write')

          if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

          write(iunit, '(a)', advance='no') '# state '
          do ivar = 1, st%nst
            do sigma = 1, em_vars%nsigma

              if( sigma == em_vars%nsigma .and. ivar == st%nst) then 
                write(iunit, '(i3)', advance='yes') (3 - 2*sigma)*ivar
              else 
                write(iunit, '(i3)', advance='no') (3 - 2*sigma)*ivar
              end if

            end do
          end do

          do ist = 1, st%nst
            write(iunit, '(i3)', advance='no') ist

            do ivar = 1, st%nst
              do sigma = 1, em_vars%nsigma

                if(states_are_complex(st)) then
                  proj = &
                       zmf_dotp(gr%mesh, st%d%dim, st%zpsi(:, :, ist, ik), em_vars%lr(idir, sigma, ifactor)%zdl_psi(:, :, ivar, ik))
                else
                  proj = &
                       dmf_dotp(gr%mesh, st%d%dim, st%dpsi(:, :, ist, ik), em_vars%lr(idir, sigma, ifactor)%ddl_psi(:, :, ivar, ik))
                end if
                  
                if( sigma == em_vars%nsigma .and. ivar == st%nst) then 
                  write(iunit, '(f12.6)', advance='yes') abs(proj)
                else 
                  write(iunit, '(f12.6,a)', advance='no') abs(proj), ' '
                end if

              end do
            end do

          end do
          call io_close(iunit)

        end do ! dir
      end do !ik

      POP_SUB(em_resp_output.out_projections)

    end subroutine out_projections


    ! ---------------------------------------------------------
    subroutine out_wfn_and_densities()
      integer :: idir, isigma

      PUSH_SUB(em_resp_output.out_wfn_and_densities)

      do idir = 1, gr%sb%dim
        if(states_are_complex(st)) then 

          if(gr%sb%dim == 3) then
            if(iand(outp%what, C_OUTPUT_ELF).ne.0) &
              call zlr_calc_elf(st, gr, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
          end if
          do isigma = 1, em_vars%nsigma
            call zoutput_lr(st, gr, em_vars%lr(idir, isigma, ifactor), dirname, idir, isigma, outp, geo, units_out%force)
          end do
        else

          if(gr%sb%dim == 3) then
            if(iand(outp%what, C_OUTPUT_ELF) .ne. 0) &
              call dlr_calc_elf(st, gr, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
          end if

          do isigma = 1, em_vars%nsigma
            call doutput_lr(st, gr, em_vars%lr(idir, isigma, ifactor), dirname, idir, isigma, outp, geo, units_out%force)
          end do

        end if
      end do

      POP_SUB(em_resp_output.out_wfn_and_densities)

    end subroutine out_wfn_and_densities
    

  ! ---------------------------------------------------------
  ! See D Varsano, LA Espinosa Leal, Xavier Andrade, MAL Marques, Rosa di Felice, Angel Rubio,
  ! Phys. Chem. Chem. Phys. 11, 4481 (2009)
    subroutine out_circular_dichroism
      type(pert_t) :: angular_momentum
      integer :: idir
      FLOAT :: ff
      CMPLX :: dic

      PUSH_SUB(em_resp_output.out_circular_dichroism)

      if(states_are_complex(st) .and. em_vars%nsigma == 2) then       

        message(1) = "Info: Calculating rotatory response."
        call messages_info(1)

        call pert_init(angular_momentum, PERTURBATION_MAGNETIC, gr, geo)
        
        dic = M_ZERO
        do idir = 1, gr%sb%dim
          call pert_setup_dir(angular_momentum, idir)
          dic = dic &
               + zpert_expectation_value(angular_momentum, gr, geo, hm, st, st%zpsi, em_vars%lr(idir, 1, ifactor)%zdl_psi) &
               + zpert_expectation_value(angular_momentum, gr, geo, hm, st, em_vars%lr(idir, 2, ifactor)%zdl_psi, st%zpsi)
        end do
        
        call pert_end(angular_momentum)
        
        dic = dic*M_zI*M_HALF

        iunit = io_open(trim(dirname)//'/rotatory_strength', action='write')

        ! print header
        write(iunit, '(a1,a20,a20,a20)') '#', str_center("Energy", 20), str_center("R", 20), str_center("Re[beta]", 20)
        write(iunit, '(a1,a20,a20,a20)') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
             str_center('['//trim(units_abbrev(units_out%length**3)) //']', 20), &
             str_center('['//trim(units_abbrev(units_out%length**4)) //']', 20)

        ff = M_ZERO
        if(em_vars%omega(iomega) .ne. 0) ff = real(dic)/(M_THREE*em_vars%omega(iomega))

        write(iunit, '(3e20.8)') units_from_atomic(units_out%energy, em_vars%omega(iomega)), &
             units_from_atomic(units_out%length**3, aimag(dic)/(P_C*M_PI)), units_from_atomic(units_out%length**4, ff)

        call io_close(iunit)
      end if
      
      POP_SUB(em_resp_output.out_circular_dichroism)

    end subroutine out_circular_dichroism
    
  end subroutine em_resp_output

  ! ---------------------------------------------------------
  ! Ref: David M Bishop, Rev Mod Phys 62, 343 (1990)
  ! beta // and _L are eqn (154), beta  k is eqn (155)
  ! generalized to lack of Kleinman symmetry
  subroutine out_hyperpolarizability(sb, beta, freq_factor, converged, dirname)
    type(simul_box_t),  intent(in) :: sb
    CMPLX,              intent(in) :: beta(:, :, :)
    FLOAT,              intent(in) :: freq_factor(:)
    logical,            intent(in) :: converged
    character(len=*),   intent(in) :: dirname

    CMPLX :: bpar(1:MAX_DIM), bper(1:MAX_DIM), bk(1:MAX_DIM)
    CMPLX :: HRS_VV, HRS_HV
    integer :: ii, jj, kk, iunit

    PUSH_SUB(out_hyperpolarizability)

    ! Output first hyperpolarizability (beta)
    iunit = io_open(trim(dirname)//'/beta', action='write')

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
               units_from_atomic(units_out%hyperpolarizability, real( beta(ii, jj, kk))), &
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
          units_from_atomic(units_out%hyperpolarizability, real(bpar(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bpar(ii)))
      end do

      write(iunit, '()')

      do ii = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta _L '//index2axis(ii), &
          units_from_atomic(units_out%hyperpolarizability, real(bper(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bper(ii)))
      end do

      write(iunit, '()')

      do ii = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta  k '//index2axis(ii), &
          units_from_atomic(units_out%hyperpolarizability, real(bk(ii))), &
          units_from_atomic(units_out%hyperpolarizability, aimag(bk(ii)))
      end do

      call calc_beta_HRS(sb, beta, HRS_VV, HRS_HV)

      write(iunit, '()')
      write(iunit, '(a)') 'beta for liquid- or gas-phase hyper-Rayleigh scattering:'
      write(iunit, '(a, 2e20.8)') 'VV polarization ', &
         units_from_atomic(units_out%hyperpolarizability, real(sqrt(HRS_VV))), &
         units_from_atomic(units_out%hyperpolarizability, aimag(sqrt(HRS_VV)))
      write(iunit, '(a, 2e20.8)') 'HV polarization ', &
         units_from_atomic(units_out%hyperpolarizability, real(sqrt(HRS_HV))), &
         units_from_atomic(units_out%hyperpolarizability, aimag(sqrt(HRS_HV)))
    endif

    call io_close(iunit)
    POP_SUB(out_hyperpolarizability)

    contains

      ! ---------------------------------------------------------
      ! calculate hyper-Rayleigh scattering hyperpolarizabilities
      ! SJ Cyvin, JE Rauch, and JC Decius, J Chem Phys 43, 4083 (1965)
      ! generalized to avoid assumption of Kleinman symmetry (permutation of indices)
      ! as in R Bersohn, Y-H Pao, and HL Frisch, J Chem Phys 45, 3184 (1966)
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
        enddo

        HRS_B = M_ZERO
        HRS_C = M_ZERO
        do ii = 1, sb%dim
          do jj = 1, sb%dim
            if (ii .ne. jj) then
              HRS_B = HRS_B + beta(ii,ii,ii) * (beta(ii,jj,jj) + beta(jj,ii,jj) + beta(jj,jj,ii))
              HRS_C = HRS_C + (beta(ii,ii,jj) + beta(ii,jj,ii) + beta(jj,ii,ii))**2
            endif
          enddo
        enddo

        HRS_D = (beta(1,1,2) + beta(1,2,1) + beta(2,1,1)) * (beta(2,3,3) + beta(3,2,3) + beta(3,3,2)) &
              + (beta(2,2,3) + beta(2,3,2) + beta(3,2,2)) * (beta(3,1,1) + beta(1,3,1) + beta(1,1,3)) &
              + (beta(3,3,1) + beta(3,1,3) + beta(1,3,3)) * (beta(1,2,2) + beta(2,1,2) + beta(2,2,1))
    
        HRS_E = (beta(1,2,3) + beta(1,3,2) + beta(2,1,3) + beta(2,3,1) + beta(3,1,2) + beta(3,2,1))**2
    
        HRS_VV = (M_ONE / M_SEVEN)     * HRS_A &
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
            if (ii .ne. jj) then
              HRS_B1 = HRS_B1 + beta(ii,ii,ii) * beta(ii,jj,jj)
              HRS_B2 = HRS_B2 + beta(ii,ii,ii) * (beta(jj,ii,jj) + beta(jj,jj,ii))
              HRS_C1 = HRS_C1 + (beta(ii,ii,jj) + beta(ii,jj,ii))**2
              HRS_C2 = HRS_C2 + beta(jj,ii,ii) * (beta(ii,ii,jj) + beta(ii,jj,ii))
              HRS_C3 = HRS_C3 + beta(jj,ii,ii)**2
            endif
          enddo
        enddo
  
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

end module em_resp_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
