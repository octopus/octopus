!! Copyright (C) 2004-2014 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca), David Strubbe
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

subroutine X(run_sternheimer)()

  PUSH_SUB(em_resp_run.X(run_sternheimer))

  call restart_init(restart_dump, RESTART_EM_RESP, RESTART_TYPE_DUMP, sys%st%dom_st_kpt_mpi_grp, &
                    mesh=sys%gr%mesh, sb=sys%gr%sb)

  if(.not. fromscratch) then

    call restart_init(restart_load, RESTART_EM_RESP, RESTART_TYPE_LOAD, sys%st%dom_st_kpt_mpi_grp, &
                      mesh=sys%gr%mesh, sb=sys%gr%sb)

    ! try to load wavefunctions, if first frequency; otherwise will already be initialized
    if(iomega == 1 .and. .not. em_vars%wfns_from_scratch) then
      do sigma = 1, em_vars%nsigma
        if(sigma == 2 .and. abs(frequency) < M_EPSILON) then
          em_vars%lr(idir, 2, ifactor)%X(dl_psi) = em_vars%lr(idir, 1, ifactor)%X(dl_psi)
          
          if(em_vars%calc_hyperpol .and. use_kdotp) then
            do idir2 = 1, gr%sb%periodic_dim
              kdotp_em_lr2(idir2, idir, 2, ifactor)%X(dl_psi) = kdotp_em_lr2(idir2, idir, 1, ifactor)%X(dl_psi)
            enddo
          endif
        else
          sigma_alt = sigma
          if(frequency < -M_EPSILON .and. em_vars%nsigma == 2) &
            sigma_alt = swap_sigma(sigma)
          
          str_tmp = em_wfs_tag(idir, ifactor)
          call restart_cd(restart_load, dirname=wfs_tag_sigma(str_tmp, sigma))
          call states_load(restart_load, sys%st, sys%gr, ierr, lr=em_vars%lr(idir, sigma_alt, ifactor))
          call restart_cd(restart_load)

          if(ierr /= 0) then
            message(1) = "Initializing to zero, could not load response wavefunctions from '" &
              //trim(wfs_tag_sigma(str_tmp, sigma))//"'"
            call messages_warning(1)
          end if
          
          if(em_vars%calc_hyperpol .and. use_kdotp) then
            do idir2 = 1, gr%sb%periodic_dim
              str_tmp = em_wfs_tag(idir, ifactor, idir2)              
              call restart_cd(restart_load, dirname=wfs_tag_sigma(str_tmp, sigma))
              call states_load(restart_load, sys%st, sys%gr, ierr, &
                lr=kdotp_em_lr2(idir2, idir, sigma_alt, ifactor))
              call restart_cd(restart_load)
              
              if(ierr /= 0) then
                message(1) = "Initializing to zero, could not load second-order response wavefunctions from '" &
                  //trim(wfs_tag_sigma(str_tmp, sigma))//"'"
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
        call X(lr_swap_sigma)(sys%st, sys%gr%mesh, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
      endif
    end if
    
    !search for the density of the closest frequency, including negative
    closest_omega = em_vars%freq_factor(ifactor) * em_vars%omega(iomega)
    call loct_search_file_lr(closest_omega, idir, ierr, trim(restart_dir(restart_load)))
    sigma_alt = 1
    if(closest_omega * frequency < M_ZERO) opp_freq = .true.
    if(opp_freq .and. em_vars%nsigma == 2) sigma_alt = 2
    
    !attempt to read 
    if(ierr == 0) then 
      call X(lr_load_rho)(em_vars%lr(idir, sigma_alt, ifactor)%X(dl_rho), sys%gr%mesh, sys%st%d%nspin, &
        restart_load, em_rho_tag(closest_omega, idir), ierr)
      
      if(ierr == 0 .and. &
        abs(abs(closest_omega) - abs(frequency)) <= CNST(1e-4)) then
        ! the frequencies are written to four decimals in the restart directory, so we cannot expect higher precision
        exact_freq = .true.
      endif
    end if
    
    if(ierr == 0 .and. em_vars%nsigma == 2) then 
      sigma_alt = 1
      if(opp_freq) sigma_alt = 2
      
      em_vars%lr(idir, swap_sigma(sigma_alt), ifactor)%X(dl_rho) = R_CONJ(em_vars%lr(idir, sigma_alt, ifactor)%X(dl_rho))
    end if
    
    call restart_end(restart_load)
  end if ! .not. fromscratch
  
  call pert_setup_dir(em_vars%perturbation, idir)
  
  if(use_kdotp .and. idir <= gr%sb%periodic_dim) then
    call X(sternheimer_set_rhs)(sh, kdotp_lr(idir, 1)%X(dl_psi))
  end if
  
  ! if the frequency is zero, we do not need to calculate both responses
  if(abs(frequency) < M_EPSILON .and. em_vars%nsigma == 2) then
    nsigma_eff = 1
  else
    nsigma_eff = em_vars%nsigma
  endif
  
  call X(sternheimer_solve)(sh, sys, hm, em_vars%lr(idir, 1:nsigma_eff, ifactor), nsigma_eff, &
    R_TOPREC(frequency_eta), em_vars%perturbation, restart_dump, &
    em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
    em_wfs_tag(idir, ifactor), have_restart_rho=(ierr==0), have_exact_freq = exact_freq)
  
  if(nsigma_eff == 1 .and. em_vars%nsigma == 2) then
    em_vars%lr(idir, 2, ifactor)%X(dl_psi) = em_vars%lr(idir, 1, ifactor)%X(dl_psi)
    em_vars%lr(idir, 2, ifactor)%X(dl_rho) = R_CONJ(em_vars%lr(idir, 1, ifactor)%X(dl_rho))
  endif
  
  if(use_kdotp) then
    call sternheimer_unset_rhs(sh)
  end if
  
  em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
  
  if(em_vars%calc_hyperpol .and. use_kdotp) then
    call X(em_resp_calc_eigenvalues)(sys, dl_eig)

    call restart_init(kdotp_restart, RESTART_KDOTP, RESTART_TYPE_LOAD, sys%st%dom_st_kpt_mpi_grp, &
                       mesh=sys%gr%mesh, sb=sys%gr%sb)

    do idir2 = 1, gr%sb%periodic_dim
      write(message(1), '(a,a,a)') 'Info: Calculating kdotp response in ', index2axis(idir2), '-direction.'
      call messages_info(1)
      call pert_setup_dir(pert_kdotp, idir2)
      
      message(1) = "Reading 2nd-order kdotp wavefunction."
      call messages_info(1)

      ! load wavefunctions
      str_tmp = kdotp_wfs_tag(min(idir, idir2), max(idir, idir2))
      ! 1 is the sigma index which is used in em_resp
      call restart_cd(kdotp_restart, dirname=wfs_tag_sigma(str_tmp, 1))
      call states_load(kdotp_restart, sys%st, sys%gr, ierr, lr=kdotp_lr2)
      call restart_cd(kdotp_restart)
      if(ierr /= 0) then
        message(1) = "Could not load 2nd-order kdotp wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'"
        message(2) = "Previous kdotp calculation (with KdotPCalcSecondOrder) required."
        call messages_fatal(2)
      end if
          
      call X(sternheimer_solve_order2)(sh, sh_kdotp, sh2, sys, hm, em_vars%lr(idir, 1:nsigma_eff, ifactor), &
        kdotp_lr(idir2, 1:1), nsigma_eff, R_TOPREC(frequency_eta), R_TOTYPE(M_ZERO), &
        em_vars%perturbation, pert_kdotp, kdotp_em_lr2(idir2, idir, 1:nsigma_eff, ifactor), &
        pert2_none, restart_dump, &
        "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true., &
        give_pert1psi2 = kdotp_lr2%X(dl_psi), give_dl_eig1 = dl_eig(:, :, idir2))
      
      ! if the frequency is zero, we do not need to calculate both responses
      if(nsigma_eff == 1 .and. em_vars%nsigma == 2) then
        kdotp_em_lr2(idir2, idir, 2, ifactor)%X(dl_psi) = kdotp_em_lr2(idir2, idir, 1, ifactor)%X(dl_psi)
      endif
      
      em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
    end do
    write(message(1), '(a)') ''
    call messages_info(1)

    call restart_end(kdotp_restart)
  end if

  call restart_end(restart_dump)

  POP_SUB(em_resp_run.X(run_sternheimer))
end subroutine X(run_sternheimer)

! ---------------------------------------------------------
subroutine X(calc_properties_linear)()

  PUSH_SUB(em_resp_run.X(calc_properties_linear))
  
  if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
    
    ! calculate polarizability
    message(1) = "Info: Calculating polarizabilities."
    call messages_info(1)
    
    if(use_kdotp) then
      call X(calc_polarizability_periodic)(sys, em_vars%lr(:, :, ifactor), kdotp_lr(:, 1), &
        em_vars%nsigma, em_vars%alpha(:, :, ifactor))
    endif

    call X(calc_polarizability_finite)(sys, hm, em_vars%lr(:, :, ifactor), em_vars%nsigma, &
      em_vars%perturbation, em_vars%alpha(:, :, ifactor), doalldirs = .not. use_kdotp)
    
    if(em_vars%calc_Born) then
      ! calculate Born effective charges
      message(1) = "Info: Calculating (frequency-dependent) Born effective charges."
      call messages_info(1)
      
      call X(forces_born_charges)(sys%gr, sys%geo, hm%ep, sys%st, &
        lr = em_vars%lr(:, 1, ifactor), lr2 = em_vars%lr(:, em_vars%nsigma, ifactor), &
        Born_charges = em_vars%Born_charges(ifactor))
    endif
    
  else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
    message(1) = "Info: Calculating magnetic susceptibilities."
    call messages_info(1)
    
    call X(lr_calc_susceptibility)(sys, hm, em_vars%lr(:,:, ifactor), em_vars%nsigma, em_vars%perturbation, &
      em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
  end if
  
  call em_resp_output(sys%st, sys%gr, hm, sys%geo, sys%outp, em_vars, iomega, ifactor)
  
  POP_SUB(em_resp_run.X(calc_properties_linear))
end subroutine X(calc_properties_linear)

! ---------------------------------------------------------
subroutine X(calc_properties_nonlinear)()

  PUSH_SUB(em_resp_run.X(calc_properties_nonlinear))

  ! calculate hyperpolarizability
  if(em_vars%calc_hyperpol) then
    write(message(1), '(a)') 'Info: Calculating hyperpolarizabilities.'
    call messages_info(1)
    
    if(use_kdotp) then
      call X(post_orthogonalize)(sys, em_vars%nfactor, em_vars%nsigma, em_vars%freq_factor(:), &
        em_vars%omega(iomega), em_vars%eta, em_vars%lr, kdotp_em_lr2)
      call X(lr_calc_beta)(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, &
        kdotp_lr = kdotp_lr(:, 1), kdotp_em_lr = kdotp_em_lr2, dl_eig = dl_eig, occ_response = .false.)
    else
      call X(lr_calc_beta)(sh, sys, hm, em_vars%lr, em_vars%perturbation, em_vars%beta, occ_response = em_vars%occ_response)
    endif
    
    str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(1)*em_vars%omega(iomega)))
    write(dirname_output, '(a, a)') EM_RESP_DIR//'freq_', trim(str_tmp)
    call io_mkdir(trim(dirname_output))
    call out_hyperpolarizability(gr%sb, em_vars%beta, em_vars%freq_factor(1:3), em_vars%ok(1), dirname_output)
  end if

  POP_SUB(em_resp_run.X(calc_properties_nonlinear))
end subroutine X(calc_properties_nonlinear)
