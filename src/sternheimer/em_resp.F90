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
! defines EM_RESP_RESTART_DIR
#define OUTPUT_DIR "linear/"

module em_resp_m
  use datasets_m
  use em_resp_calc_m
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kdotp_m
  use kdotp_calc_m
  use lalg_basic_m
  use loct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use h_sys_output_m
  use pert_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       em_resp_run,            &
       out_polarizability,     &
       out_hyperpolarizability

  type em_resp_t
    type(pert_t) :: perturbation

    integer :: nsigma ! 1: consider only positive values of the frequency
                      ! 2: consider both positive and negative
    integer :: nfactor! 1: only one frequency needed
                      ! 3: three frequencies (for the hyperpolarizabilities)
    integer :: nomega ! number of frequencies to consider

    FLOAT :: eta                     ! small imaginary part to add to the frequency
    FLOAT :: freq_factor(MAX_DIM)    !
    FLOAT,      pointer :: omega(:)  ! the frequencies to consider
    type(lr_t), pointer :: lr(:,:,:) ! linear response for (gr%mesh%sb%dim, nsigma, nfactor)

    logical :: calc_hyperpol
    CMPLX   :: alpha(MAX_DIM, MAX_DIM, 3)        ! the linear polarizability
    CMPLX   :: beta (MAX_DIM, MAX_DIM, MAX_DIM)  ! first hyperpolarizability

    CMPLX   :: chi_para(MAX_DIM, MAX_DIM, 3)     ! The paramagnetic part of the susceptibility
    CMPLX   :: chi_dia (MAX_DIM, MAX_DIM, 3)     ! The diamagnetic  part of the susceptibility

    logical :: ok(1:3)                           ! whether calculation is converged
    logical :: force_no_kdotp                    ! whether to use kdotp run for periodic system
    CMPLX :: Born_sum(MAX_DIM, MAX_DIM)          ! sum over atoms of Born charge tensors
  end type em_resp_t

contains

  ! ---------------------------------------------------------
  subroutine em_resp_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(em_resp_t)         :: em_vars
    type(sternheimer_t)     :: sh
    type(lr_t)              :: kdotp_lr(MAX_DIM, 1)

    integer :: sigma, ndim, i, idir, ierr, iomega, ifactor
    character(len=100) :: dirname, str_tmp
    logical :: complex_response, have_to_calculate, use_kdotp

    FLOAT :: closest_omega

    call push_sub('em_resp.em_resp_run')

    gr => sys%gr
    ndim = sys%gr%sb%dim

    call parse_input()

    complex_response = (em_vars%eta /= M_ZERO ) .or. states_are_complex(sys%st)
    call restart_look_and_read(sys%st, sys%gr, sys%geo, is_complex = complex_response)

    use_kdotp = simul_box_is_periodic(gr%sb) .and. .not. em_vars%force_no_kdotp
    ! read kdotp wavefunctions if necessary
    if (use_kdotp) then
      message(1) = "Reading kdotp wavefunctions since system is periodic."
      call write_info(1)

      do idir = 1, gr%mesh%sb%dim
        call lr_init(kdotp_lr(idir, 1))
        call lr_allocate(kdotp_lr(idir, 1), sys%st, sys%gr%mesh)

        ! load wave-functions
        str_tmp = kdotp_wfs_tag(idir)
        write(dirname,'(3a)') "kdotp/", trim(str_tmp), '_1'
        ! 1 is the sigma index which is used in em_resp
        call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
          ierr, lr=kdotp_lr(idir, 1))

        if(ierr.ne.0) then
          message(1) = "Could not load kdotp wavefunctions from '"//trim(tmpdir)//trim(dirname)//"'"
          message(2) = "Previous kdotp calculation required."
          call write_fatal(2)
        end if
      end do
    endif

    em_vars%nfactor = 1
    if(em_vars%calc_hyperpol) em_vars%nfactor = 3

    ! in effect, nsigma = 1 only if hyperpol not being calculated, and the only frequency is zero
    if(em_vars%calc_hyperpol .or. (em_vars%nomega > 1) .or. (abs(em_vars%omega(1)) >= M_EPSILON) ) then
       em_vars%nsigma = 2
       ! positive and negative values of the frequency must be considered
    else
       em_vars%nsigma = 1
       ! only considering positive values
    endif
    
    ALLOCATE(em_vars%lr(1:gr%mesh%sb%dim, 1:em_vars%nsigma, 1:em_vars%nfactor), gr%mesh%sb%dim*em_vars%nsigma*em_vars%nfactor)
    em_vars%lr(1:gr%mesh%sb%dim, 1:em_vars%nsigma, 1:em_vars%nfactor)%nst = sys%st%nst

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    call system_h_setup(sys, hm)
    
    if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
       call sternheimer_init(sh, sys, hm, "EM", hermitian = states_are_real(sys%st), set_ham_var = 0)
       ! set HamiltonVariation to V_ext_only, in magnetic case
    else
       call sternheimer_init(sh, sys, hm, "EM", hermitian = states_are_real(sys%st))
       ! otherwise, use default, which is hartree + fxc
    endif

    call io_mkdir(trim(tmpdir)//EM_RESP_RESTART_DIR)
    call info()

    call io_mkdir(OUTPUT_DIR)

    do ifactor = 1, em_vars%nfactor
      do idir = 1, sys%gr%sb%dim
        do sigma = 1, em_vars%nsigma
          call lr_init(em_vars%lr(idir, sigma, ifactor))
          call lr_allocate(em_vars%lr(idir, sigma, ifactor), sys%st, sys%gr%mesh)
        end do
      end do
    end do

    do iomega = 1, em_vars%nomega

      em_vars%ok(1:3) = .true.

      do ifactor = 1, em_vars%nfactor
        do idir = 1, sys%gr%sb%dim

          ierr = 0

          have_to_calculate = .true.
          ! if this frequency is zero and this is not the first
          ! iteration we do not have to do anything
          if( iomega > 1 .and. em_vars%freq_factor(ifactor) == M_ZERO) have_to_calculate = .false. 
            
          if(ifactor > 1) then 

            ! if this frequency is the same as the previous one, just copy it
            if( have_to_calculate .and. &
                 em_vars%freq_factor(ifactor)*em_vars%omega(iomega) == & 
                 em_vars%freq_factor(ifactor-1)*em_vars%omega(iomega) ) then
              
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor-1), em_vars%lr(idir, 1, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor-1), em_vars%lr(idir, 2, ifactor))
              
              have_to_calculate = .false.
              
            end if

            ! if this frequency is minus the previous one, copy it inverted
            if( have_to_calculate .and. & 
                 em_vars%freq_factor(ifactor) == -em_vars%freq_factor(ifactor-1) ) then 
              
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor-1), em_vars%lr(idir, 2, ifactor))
              call lr_copy(sys%st, sys%gr%mesh, em_vars%lr(idir, 2, ifactor-1), em_vars%lr(idir, 1, ifactor))
              
              have_to_calculate = .false.
              
            end if

          end if

          if(have_to_calculate) then 

            str_tmp = freq2str(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor)
            write(message(1), '(a,i1,2a)') 'Info: Calculating response for direction ', idir, &
              ' and frequency ' , trim(str_tmp)
            call write_info(1)

            if(.not. fromscratch) then 

               ! try to load wavefunctions, if first frequency; otherwise will already be initialized
               if(iomega == 1) then
                  do sigma = 1, em_vars%nsigma
                     str_tmp =  em_wfs_tag(idir, ifactor)
                     write(dirname,'(3a, i1)') EM_RESP_RESTART_DIR, trim(str_tmp), '_', sigma
                     call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
                          ierr, lr=em_vars%lr(idir, sigma, ifactor))

                     if(ierr.ne.0) then
                        message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//trim(dirname)//"'"
                        call write_warning(1)
                     end if
                  end do
               end if

              !try to load restart density
              if (states_are_complex(sys%st)) then 
                call zrestart_read_lr_rho(em_vars%lr(idir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                     EM_RESP_RESTART_DIR, &
                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), idir), ierr)
              else 
                call drestart_read_lr_rho(em_vars%lr(idir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                     EM_RESP_RESTART_DIR, &
                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), idir), ierr)
              end if

              !search for the density of the closest frequency
              if(ierr /= 0) then 
                
                closest_omega = em_vars%freq_factor(ifactor)*em_vars%omega(iomega)
                call oct_search_file_lr(closest_omega, idir, ierr, trim(tmpdir)//EM_RESP_RESTART_DIR)
                
                !attempt to read 
                if(ierr == 0 ) then 
                  if (states_are_complex(sys%st)) then 
                    call zrestart_read_lr_rho(em_vars%lr(idir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                         EM_RESP_RESTART_DIR, em_rho_tag(closest_omega, idir), ierr)
                  else 
                    call drestart_read_lr_rho(em_vars%lr(idir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                         EM_RESP_RESTART_DIR, em_rho_tag(closest_omega, idir), ierr)
                  end if
                end if
                
              end if

              if(ierr == 0 .and. em_vars%nsigma == 2 ) then 
                if (states_are_complex(sys%st)) then 
                  em_vars%lr(idir, 2, ifactor)%zdl_rho = conjg(em_vars%lr(idir, 1, ifactor)%zdl_rho)
                else 
                  em_vars%lr(idir, 2, ifactor)%ddl_rho = em_vars%lr(idir, 1, ifactor)%ddl_rho
                end if
              end if

            end if ! .not. fromscratch
            
            call pert_setup_dir(em_vars%perturbation, idir)

            if(use_kdotp) then
              kdotp_lr(idir, 1)%zdl_psi = M_ZI*kdotp_lr(idir, 1)%zdl_psi
              call zsternheimer_set_rhs(sh, kdotp_lr(idir, 1)%zdl_psi)
            end if

            if (states_are_complex(sys%st)) then 
              call zsternheimer_solve(sh, sys, hm, em_vars%lr(idir, :, ifactor), em_vars%nsigma, &
                em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, &
                em_vars%perturbation, EM_RESP_RESTART_DIR, &
                em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), idir), &
                em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0))
            else
              call dsternheimer_solve(sh, sys, hm, em_vars%lr(idir, :, ifactor), em_vars%nsigma, &
                em_vars%freq_factor(ifactor)*em_vars%omega(iomega), &
                em_vars%perturbation, EM_RESP_RESTART_DIR, &
                em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), idir), &
                em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0))
            end if
            
            if(use_kdotp) then
              kdotp_lr(idir, 1)%zdl_psi = -M_ZI*kdotp_lr(idir, 1)%zdl_psi
              call sternheimer_unset_rhs(sh)
            end if

            em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
            
          end if

        end do ! idir
      end do ! ifactor
      
      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then

        ! calculate polarizability
         message(1) = "Info: Calculating polarizabilities."
         call write_info(1)

        do ifactor = 1, em_vars%nfactor
          if(use_kdotp) then
            call zcalc_polarizability_periodic(sys, em_vars%lr(:, :, ifactor), kdotp_lr(:, 1), &
              em_vars%nsigma, em_vars%alpha(:, :, ifactor))
          else if(states_are_complex(sys%st)) then
            call zcalc_polarizability_finite(sys, hm, em_vars%lr(:, :, ifactor), em_vars%nsigma, &
              em_vars%perturbation, em_vars%alpha(:, :, ifactor))
          else
            call dcalc_polarizability_finite(sys, hm, em_vars%lr(:, :, ifactor), em_vars%nsigma, &
              em_vars%perturbation, em_vars%alpha(:, :, ifactor))
          end if
        end do

        ! calculate Born effective charges
         message(1) = "Info: Calculating (frequency-dependent) Born effective charges."
         call write_info(1)

        do ifactor = 1, em_vars%nfactor
          do idir = 1, sys%gr%sb%dim
            ! time = M_ZERO
            if(states_are_complex(sys%st)) then
              if(em_vars%nsigma == 2) then
                call zcalc_forces_from_potential(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                  lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 2, ifactor), &
                  lr_dir = idir, Born_sum = em_vars%Born_sum(idir, :))
              else
                call zcalc_forces_from_potential(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                  lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 1, ifactor), &
                  lr_dir = idir, Born_sum = em_vars%Born_sum(idir, :))
              endif
            else
              call dcalc_forces_from_potential(sys%gr, sys%geo, hm%ep, sys%st, M_ZERO, &
                lr = em_vars%lr(idir, 1, ifactor), lr2 = em_vars%lr(idir, 1, ifactor), &
                lr_dir = idir, Born_sum = em_vars%Born_sum(idir, :))
            endif
          enddo
        enddo

        ! calculate hyperpolarizability
        if(em_vars%calc_hyperpol) then
          if(states_are_complex(sys%st)) then
            call zlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta)
          else
            call dlr_calc_beta(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta)
          end if
        end if
      
      else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
        message(1) = "Info: Calculating magnetic susceptibilities."
        call write_info(1)

        do ifactor = 1, em_vars%nfactor
          if(states_are_complex(sys%st)) then 
            call zlr_calc_susceptibility(sys, hm, em_vars%lr(:,:, ifactor), em_vars%nsigma, em_vars%perturbation, &
               em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          else
            call dlr_calc_susceptibility(sys, hm, em_vars%lr(:,:, ifactor), em_vars%nsigma, em_vars%perturbation, &
               em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          end if
        end do
      end if

      call em_resp_output(sys%st, sys%gr, hm, sys%geo, sys%outp, em_vars, iomega)

    end do

    do idir = 1, ndim
      do sigma = 1, em_vars%nsigma
        do ifactor = 1, em_vars%nfactor
          call lr_dealloc(em_vars%lr(idir, sigma, ifactor))
        end do
      end do
    end do

    call sternheimer_end(sh)
    call pert_end(em_vars%perturbation)

    deallocate(em_vars%omega, em_vars%lr)
    call states_deallocate_wfns(sys%st)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine parse_input()
      type(block_t) :: blk
      integer   :: nrow
      integer   :: number, j, k
      FLOAT     :: omega_ini, omega_fin, domega

      call push_sub('em_resp.em_resp_run.parse_input')

      call obsolete_variable('PolFreqs               ', 'EMFreqs             ')
      call obsolete_variable('PolHyper               ', 'EMHyperpol          ')
      call obsolete_variable('PolEta                 ', 'EMEta               ')
      call obsolete_variable('PolConvAbsDens         ', 'LRConvAbsDens       ')
      call obsolete_variable('PolHamiltonianVariation', 'HamiltonianVariation')

      !%Variable EMFreqs
      !%Type block
      !%Section Linear Response::Polarizabilities
      !%Description
      !% This block defines for which frequencies the polarizabilities
      !% will be calculated. If is not present the static (omega = 0) response
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

      if (loct_parse_block(datasets_check('EMFreqs'), blk) == 0) then 

        nrow = loct_parse_block_n(blk)
        em_vars%nomega = 0

        !count the number of frequencies
        do i = 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          if(number < 1) then
            message(1) = "EMFreqs: invalid number of frequencies"
            call write_fatal(1)
          end if
          em_vars%nomega = em_vars%nomega + number
        end do

        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega)

        !read frequencies
        j = 1
        do i = 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          call loct_parse_block_float(blk, i, 1, omega_ini)
          if(number > 1) then 
            call loct_parse_block_float(blk, i, 2, omega_fin)
            domega = (omega_fin - omega_ini)/(number - M_ONE)
            do k = 0, number-1
              em_vars%omega(j + k) = (omega_ini + domega*k)*units_inp%energy%factor
            end do
            j = j + number
          else
            em_vars%omega(j) = omega_ini*units_inp%energy%factor
            j = j + 1
          end if
        end do

        call loct_parse_block_end(blk)

        call sort(em_vars%omega)

      else
        !there is no frequency block, we calculate response for w = 0.0
        em_vars%nomega = 1
        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega)
        em_vars%omega(1) = M_ZERO
      end if

      !%Variable EMEta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Imaginary part of the frequency.
      !%End

      call loct_parse_float(datasets_check('EMEta'), M_ZERO, em_vars%eta)
      em_vars%eta = em_vars%eta*units_inp%energy%factor

      ! reset the values of these variables
      em_vars%calc_hyperpol = .false.
      em_vars%freq_factor(1:MAX_DIM) = M_ONE

      call pert_init(em_vars%perturbation, sys%gr, sys%geo)

      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        !%Variable EMHyperpol
        !%Type block
        !%Section Linear Response::Polarizabilities
        !%Description
        !% This blocks describes the multiples of the frequency used for
        !% the dynamic hyperpolarizability.
        !%End

        if (loct_parse_block(datasets_check('EMHyperpol'), blk) == 0) then 
          call loct_parse_block_float(blk, 0, 0, em_vars%freq_factor(1))
          call loct_parse_block_float(blk, 0, 1, em_vars%freq_factor(2))
          call loct_parse_block_float(blk, 0, 2, em_vars%freq_factor(3))
          
          call loct_parse_block_end(blk)
          
          em_vars%calc_hyperpol = .true.
        end if
      end if

      !%Variable EMForceNoKdotP
      !%Type logical
      !%Default false
      !%Section Linear Response::Polarizabilities
      !%Description
      !% If the system is periodic, by default wavefunctions from a previous kdotp run will
      !% be read, to be used in the formulas for the polarizability and
      !% hyperpolarizability in the quantum theory of polarization. For testing purposes,
      !% you can set this variable to true to disregard the kdotp run, and use the formulas
      !% for the finite system. This variable has no effect for a finite system.
      !%End

      call loct_parse_logical(datasets_check('EMForceNoKdotP'), .false., em_vars%force_no_kdotp)

      call pop_sub()

    end subroutine parse_input


    ! ---------------------------------------------------------
    subroutine info()

      call push_sub('em_resp.em_resp_run.info')

      call pert_info(em_vars%perturbation, stdout)
      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        if(em_vars%calc_hyperpol) then 
          write(message(1),'(a)') 'Linear-Reponse First-Order Hyperpolarizabilities'
          call messages_print_stress(stdout, trim(message(1)))
        else 
          write(message(1),'(a)') 'Linear-Reponse Polarizabilities'
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
      call write_info(1)

      write(message(1),'(a,i3,a)') 'Calculating response for ', em_vars%nomega, ' frequencies.'
      call write_info(1)

      call messages_print_stress(stdout)

      call pop_sub()

    end subroutine info

  end subroutine em_resp_run


  ! ---------------------------------------------------------
  subroutine em_resp_output(st, gr, hm, geo, outp, em_vars, iomega)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(geometry_t),     intent(inout) :: geo
    type(h_sys_output_t), intent(in)    :: outp
    type(em_resp_t),      intent(inout) :: em_vars
    integer,              intent(in)    :: iomega
    
    integer :: iunit, ifactor
    character(len=80) :: dirname, str_tmp

    call push_sub('em_resp.em_resp_output')

    do ifactor = 1, em_vars%nfactor
      str_tmp = freq2str(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor)
      write(dirname, '(a, a)') OUTPUT_DIR//'freq_', trim(str_tmp)
      call io_mkdir(trim(dirname))

! TODO: make out_polarizability and out_susceptibility subroutines of this subroutine again
      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
        call out_polarizability(st, gr, em_vars%alpha(1:MAX_DIM, 1:MAX_DIM, ifactor), &
          em_vars%freq_factor(ifactor) * em_vars%omega(iomega), em_vars%ok(ifactor), dirname)
        call out_Born_charges()
        if(em_vars%calc_hyperpol) call out_hyperpolarizability(gr%sb, em_vars%beta, em_vars%ok(ifactor), dirname)
      else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
        call out_susceptibility()
      end if

      call out_projections()

      str_tmp = freq2str(em_vars%omega(iomega)/units_out%energy%factor)
      write(dirname, '(a, a)') OUTPUT_DIR//'freq_', trim(str_tmp)
      call io_mkdir(dirname)
      call out_wavefunctions()

      call out_circular_dichroism
    end do

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine out_Born_charges()
      integer iatom
      FLOAT :: phase(1:MAX_DIM, 1:MAX_DIM)

      call push_sub('em_resp.em_resp_output.out_Born_charges')
      iunit = io_open(trim(dirname)//'/Born_charges', action='write')

      write(iunit,'(a)') '# (Frequency-dependent) Born effective charges'

      do iatom = 1, geo%natoms
         write(iunit,'(a)')
         write(iunit,'(i5,a10)') iatom, trim(geo%atom(iatom)%spec%label)

         write(iunit,'(a)') 'Magnitude:'
         call io_output_tensor(iunit, TOFLOAT(abs(geo%atom(iatom)%Born_charge(:,:))), gr%mesh%sb%dim, M_ONE)

         write(iunit,'(a)') 'Phase:'
         phase = atan2(aimag(geo%atom(iatom)%Born_charge(:,:)),real(geo%atom(iatom)%Born_charge(:,:)))
         call io_output_tensor(iunit, phase, gr%mesh%sb%dim, M_ONE)
      enddo

      write(iunit,'(a)')
      write(iunit,'(a)') '# Sum of Born effective charges before correction to satisfy acoustic sum rule.' 
      write(iunit,'(a)') 'Real:'
      call io_output_tensor(iunit, real(em_vars%Born_sum(:, :)), gr%mesh%sb%dim, M_ONE)
      write(iunit,'(a)') 'Imaginary:'
      call io_output_tensor(iunit, aimag(em_vars%Born_sum(:, :)), gr%mesh%sb%dim, M_ONE)

      call io_close(iunit)
      call pop_sub()
    end subroutine out_Born_charges

    ! ---------------------------------------------------------
    subroutine out_susceptibility()
      FLOAT :: to_ppmcgs

      call push_sub('em_resp.em_resp_output.out_susceptibility')

      iunit = io_open(trim(dirname)//'/susceptibility', action='write')

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor)), gr%mesh%sb%dim, CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :, ifactor)), gr%mesh%sb%dim, CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm a.u.]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor) + em_vars%chi_dia(:,:, ifactor)), &
        gr%mesh%sb%dim, CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(a)') hyphens
      to_ppmcgs = M_ONE/CNST(8.9238878e-2)*CNST(1e-6)

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor)), gr%mesh%sb%dim, to_ppmcgs)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_dia(:, :, ifactor)), gr%mesh%sb%dim, to_ppmcgs)
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm cgs / mol]'
      call io_output_tensor(iunit, TOFLOAT(em_vars%chi_para(:, :, ifactor) + em_vars%chi_dia(:,:, ifactor)), &
           gr%mesh%sb%dim, to_ppmcgs)
      write(iunit, '(1x)')

      call io_close(iunit)      
      call pop_sub()
    end subroutine out_susceptibility

    ! ---------------------------------------------------------
    subroutine out_projections()
      CMPLX   :: proj
      integer :: ist, ivar, ik, dir, sigma
      character(len=80) :: fname

      call push_sub('em_resp.em_resp_output.out_projections')

      do ik = st%d%kpt%start, st%d%kpt%end
        do dir = 1, gr%mesh%sb%dim

          write(fname, '(2a,i1,a,i1)') trim(dirname), '/projection-', ik, '-', dir
          iunit = io_open(trim(fname), action='write')

          if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

          write(iunit, '(a)', advance='no') '# state '
          do ivar = 1, em_vars%lr(dir, 1, 1)%nst
            do sigma = 1, em_vars%nsigma

              if( sigma == em_vars%nsigma .and. ivar == em_vars%lr(dir, 1, 1)%nst) then 
                write(iunit, '(i3)', advance='yes') (3 - 2*sigma)*ivar
              else 
                write(iunit, '(i3)', advance='no') (3 - 2*sigma)*ivar
              end if

            end do
          end do

          do ist = 1, st%nst
            write(iunit, '(i3)', advance='no') ist

            do ivar = 1, em_vars%lr(dir, 1, 1)%nst
              do sigma = 1, em_vars%nsigma

                if(states_are_complex(st)) then
                  proj = &
                       zmf_dotp(gr%mesh, st%d%dim, st%zpsi(:, :, ist, ik), em_vars%lr(dir, sigma, ifactor)%zdl_psi(:, :, ivar, ik))
                else
                  proj = &
                       dmf_dotp(gr%mesh, st%d%dim, st%dpsi(:, :, ist, ik), em_vars%lr(dir, sigma, ifactor)%ddl_psi(:, :, ivar, ik))
                end if
                  
                if( sigma == em_vars%nsigma .and. ivar == em_vars%lr(dir, 1, 1)%nst) then 
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

      call pop_sub()

    end subroutine out_projections


    ! ---------------------------------------------------------
    subroutine out_wavefunctions()
      integer :: dir, isigma

      call push_sub('em_resp.em_resp_output.out_wavefunctions')

      do dir = 1, gr%mesh%sb%dim
        if(states_are_complex(st)) then 

          if(gr%mesh%sb%dim==3) then
            if(iand(outp%what, output_elf).ne.0) &
              call zlr_calc_elf(st, gr, em_vars%lr(dir, 1, ifactor), em_vars%lr(dir, 2, ifactor))
          end if
          do isigma = 1, em_vars%nsigma
            call zh_sys_output_lr(st, gr, em_vars%lr(dir, isigma, ifactor), dirname, dir, isigma, outp, geo)
          end do
        else

          if(gr%mesh%sb%dim==3) then
            if(iand(outp%what, output_elf) .ne. 0) &
              call dlr_calc_elf(st, gr, em_vars%lr(dir, 1, ifactor), em_vars%lr(dir, 2, ifactor))
          end if

          do isigma = 1, em_vars%nsigma
            call dh_sys_output_lr(st, gr, em_vars%lr(dir, isigma, ifactor), dirname, dir, isigma, outp, geo)
          end do

        end if
      end do

      call pop_sub()

    end subroutine out_wavefunctions
    

  ! ---------------------------------------------------------
    subroutine out_circular_dichroism
      type(pert_t) :: angular_momentum
      integer :: idir
      FLOAT :: ff
      CMPLX :: dic

      call push_sub('em_resp.em_resp_output.out_circular_dichroism')

      if(states_are_complex(st) .and. em_vars%nsigma == 2) then       

        message(1) = "Info: Calculating rotatory response."
        call write_info(1)

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
        write(iunit, '(a1,a20,a20,a20)') '#', str_center('['//trim(units_out%energy%abbrev) // ']', 20), &
             str_center('['//trim(units_out%length%abbrev) //'^3]', 20), &
             str_center('['//trim(units_out%length%abbrev) //'^4]', 20)

        ff = M_ZERO
        if(em_vars%omega(iomega) .ne. 0) ff = real(dic)/(M_THREE*em_vars%omega(iomega))

        write(iunit, '(3e20.8)') em_vars%omega(iomega)/units_inp%energy%factor, &
             aimag(dic)/(P_C*M_PI*units_out%length%factor**3), ff/units_out%length%factor**4

        call io_close(iunit)
      end if
      
      call pop_sub()

    end subroutine out_circular_dichroism
    
  end subroutine em_resp_output


  ! ---------------------------------------------------------
  subroutine out_polarizability(st, gr, alpha, omega, converged, dirname)
    type(states_t),    intent(in) :: st
    type(grid_t),      intent(in) :: gr
    CMPLX,             intent(in) :: alpha(:, :)
    FLOAT,             intent(in) :: omega
    logical,           intent(in) :: converged
    character(len=*),  intent(in) :: dirname

    FLOAT :: cross(MAX_DIM, MAX_DIM), crossp(MAX_DIM, MAX_DIM)
    FLOAT :: average, anisotropy
    integer iunit
    
    call push_sub('em_resp.out_polarizability')

    iunit = io_open(trim(dirname)//'/alpha', action='write')

    if (.not.converged) write(iunit, '(a)') "# WARNING: not converged"

    write(iunit, '(2a)', advance='no') '# Polarizability tensor [', &
      trim(units_out%length%abbrev)
    if(gr%mesh%sb%dim.ne.1) write(iunit, '(a,i1)', advance='no') '^', gr%mesh%sb%dim
    write(iunit, '(a)') ']'

    call io_output_tensor(iunit, TOFLOAT(alpha(1:MAX_DIM, 1:MAX_DIM)), gr%mesh%sb%dim, units_out%length%factor**gr%mesh%sb%dim)

    call io_close(iunit)

    ! CROSS SECTION (THE IMAGINARY PART OF POLARIZABILITY)
    if(states_are_complex(st)) then 
      cross(1:MAX_DIM, 1:MAX_DIM) = aimag(alpha(1:MAX_DIM, 1:MAX_DIM)) * &
        omega / units_out%energy%factor * M_FOUR * M_PI / P_c 
        
      iunit = io_open(trim(dirname)//'/cross_section', action='write')
      if (.not. converged) write(iunit, '(a)') "# WARNING: not converged"

      average = M_THIRD * (cross(1, 1) + cross(2, 2) + cross(3, 3))
      crossp(:, :) = matmul(cross(:, :), cross(:, :))
      anisotropy = &
        M_THIRD*(M_THREE*(crossp(1, 1) + crossp(2, 2) + crossp(3, 3)) - (cross(1, 1) + cross(2, 2) + cross(3, 3))**2)
          
      call cross_section_header(iunit)
      write(iunit,'(3e20.8)', advance = 'no') &
        omega / units_out%energy%factor, average, sqrt(max(anisotropy, M_ZERO))
      write(iunit,'(9e20.8)', advance = 'no') cross(1:3, 1:3)
      write(iunit,'(a)', advance = 'yes')

      call io_close(iunit)
    end if
    
    call pop_sub()

  contains

    ! ---------------------------------------------------------
    ! Note: this should be in spectrum.F90
    subroutine cross_section_header(out_file)
      integer, intent(in) :: out_file

      character(len=80) :: header_string
      integer :: i, k

      call push_sub('em_resp.out_polarizability.cross_section_header')

      !this header is the same as spectrum.F90
      write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
      write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
      write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)

      do i = 1, 3
        do k = 1, 3
          write(header_string,'(a6,i1,a1,i1,a1)') 'sigma(',i,',',k,')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do

      write(out_file, *)
      write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_out%energy%abbrev) // ']', 20)
      do i = 1, 11
        write(out_file, '(a20)', advance = 'no')  str_center('['//trim(units_out%length%abbrev) //'^2]', 20)
      end do
      write(out_file,*)

      call pop_sub()
    end subroutine cross_section_header

  end subroutine out_polarizability


  ! ---------------------------------------------------------
  subroutine out_hyperpolarizability(sb, beta, converged, dirname)
    type(simul_box_t),  intent(in) :: sb
    CMPLX,              intent(in) :: beta(:, :, :)
    logical,            intent(in) :: converged
    character(len=*),   intent(in) :: dirname

    character, parameter :: axis(1:3) = (/ 'x', 'y', 'z' /)
    CMPLX :: bpar(1:MAX_DIM), bper(1:MAX_DIM), bk(1:MAX_DIM)
    FLOAT :: HRS_VV, HRS_HV, HRS_A, HRS_B, HRS_C, HRS_D, HRS_E
    integer :: i, j, k, iunit

    call push_sub('em_resp_out_hyperpolarizability')

    ! Output first hyperpolarizability (beta)
    iunit = io_open(trim(dirname)//'/beta', action='write')

    if (.not. converged) write(iunit, '(a)') "# WARNING: not converged"

    write(iunit, '(2a)', advance='no') 'First hyperpolarizability tensor: beta [', trim(units_out%length%abbrev)
    write(iunit, '(a,i1)', advance='no') '^', 5
    write(iunit, '(a)') ']'

    write(iunit, '()')

    do i = 1, sb%dim
      do j = 1, sb%dim
        do k = 1, sb%dim
          write(iunit,'(a,e20.8,e20.8)') 'beta '//axis(i)//axis(j)//axis(k)//' ', &
               real( beta(i, j, k))/units_out%length%factor**(5), &
               aimag(beta(i, j, k))/units_out%length%factor**(5)
        end do
      end do
    end do

    if (sb%dim == 3) then 
      bpar = M_ZERO
      bper = M_ZERO

      do i = 1, sb%dim
        do j = 1, sb%dim
          bpar(i) = bpar(i) + beta(i, j, j) + beta(j, i, j) + beta(j, j, i)
          bper(i) = bper(i) + M_TWO*beta(i, j, j) - M_THREE*beta(j, i, j) + M_TWO*beta(j, j, i)
        end do
      end do

      write(iunit, '()')

      bpar = bpar/(M_FIVE * units_out%length%factor**(5))
      bper = bper/(M_FIVE * units_out%length%factor**(5))
      bk(1:sb%dim) = M_THREE*M_HALF*(bpar(1:sb%dim) - bper(1:sb%dim))

      do i = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta // '//axis(i), real(bpar(i)), aimag(bpar(i))
      end do

      write(iunit, '()')

      do i = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta _L '//axis(i), real(bper(i)), aimag(bper(i))
      end do

      write(iunit, '()')

      do i = 1, sb%dim
        write(iunit, '(a, 2e20.8)') 'beta  k '//axis(i), real(bk(i)), aimag(bk(i))
      end do

      ! calculate hyper-Rayleigh scattering hyperpolarizabilities
      ! SJ Cyvin, JE Rauch, and JC Decius, J Chem Phys 43, 4083 (1965)
      HRS_VV = M_ZERO
      HRS_HV = M_ZERO

      HRS_A = M_ZERO
      do i = 1, sb%dim
        HRS_A = HRS_A + beta(i, i, i)**2
      enddo

      HRS_B = M_ZERO
      HRS_C = M_ZERO
      do i = 1, sb%dim
        do j = 1, sb%dim
          if (i .ne. j) then
            HRS_B = HRS_B + beta(i, i, i) * beta(i, j, j)
            HRS_C = HRS_C + beta(i, i, j)**2
          endif
        enddo
      enddo

      HRS_D = beta(1, 1, 2) * beta(2, 3, 3) &
            + beta(2, 2, 3) * beta(3, 1, 1) &
            + beta(3, 3, 1) * beta(1, 2, 2)

      HRS_E = beta(1, 2, 3)**2

      HRS_VV = (M_ONE      / M_SEVEN)    * HRS_A &
             + (M_SIX      / CNST(35.0)) * HRS_B &
             + (M_NINE     / CNST(35.0)) * HRS_C &
             + (M_SIX      / CNST(35.0)) * HRS_D &
             + (CNST(12.0) / CNST(35.0)) * HRS_E

      HRS_HV = (M_ONE      / CNST(35.0))  * HRS_A &
             - (M_TWO      / CNST(105.0)) * HRS_B &
             + (CNST(11.0) / CNST(105.0)) * HRS_C &
             - (M_TWO      / CNST(105.0)) * HRS_D &
             + (M_EIGHT    / CNST(35.0))  * HRS_E

      write(iunit, '()')
      write(iunit, '(a)') 'beta hyper-Rayleigh scattering:'
      write(iunit, '(a, e20.8)') 'VV polarization ', sqrt(HRS_VV)
      write(iunit, '(a, e20.8)') 'HV polarization ', sqrt(HRS_HV)
    endif

    call io_close(iunit)
    call pop_sub()

  end subroutine out_hyperpolarizability

end module em_resp_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
