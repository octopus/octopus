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

subroutine X(run_sternheimer)(em_vars, namespace, space, gr, kpoints, st, hm, xc, mc, ions)
  type(em_resp_t),          intent(inout) :: em_vars
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(kpoints_t),          intent(in)    :: kpoints
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(inout) :: hm
  type(xc_t),               intent(in)    :: xc
  type(multicomm_t),        intent(in)    :: mc
  type(ions_t),             intent(in)    :: ions

  integer, parameter :: PB = 1, PK2 = 2, PKB = 3, PKE = 4, PE = 5
  integer :: ik, ik0, ist
  integer :: sigma, sigma_alt, idir, idir2, ierr, nsigma_eff, ipert
  FLOAT :: closest_omega
  character(len=100) :: str_tmp
  type(restart_t) :: restart_load, restart_dump
  type(lr_t) :: kdotp_lr2
  R_TYPE, allocatable :: inhomog(:,:,:,:,:)

  PUSH_SUB(em_resp_run_legacy.X(run_sternheimer))

  call restart_init(restart_dump, namespace, RESTART_EM_RESP, RESTART_TYPE_DUMP, mc, ierr, mesh=gr%mesh)

  if(.not. fromscratch) then

    call restart_init(restart_load, namespace, RESTART_EM_RESP, RESTART_TYPE_LOAD, mc, ierr, mesh=gr%mesh)

    do idir = 1, space%dim

      ! try to load wavefunctions, if first frequency; otherwise will already be initialized
      if(iomega == 1 .and. .not. em_vars%wfns_from_scratch) then
        do sigma = 1, em_vars%nsigma
          if(sigma == 2 .and. abs(frequency) < M_EPSILON) then
            em_vars%lr(idir, 2, ifactor)%X(dl_psi) = em_vars%lr(idir, 1, ifactor)%X(dl_psi)
          
            if(em_vars%calc_hyperpol .and. use_kdotp) then
              do idir2 = 1, space%periodic_dim
                kdotp_em_lr2(idir2, idir, 2, ifactor)%X(dl_psi) = kdotp_em_lr2(idir2, idir, 1, ifactor)%X(dl_psi)
              end do
            end if
          else
            sigma_alt = sigma
            if(frequency < -M_EPSILON .and. em_vars%nsigma == 2) &
              sigma_alt = swap_sigma(sigma)
          
            if((.not. em_vars%calc_magnetooptics) .or. complex_wfs .or. ifactor == 1) then
              str_tmp = em_wfs_tag(idir, ifactor)
              call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, sigma), ierr)
              if (ierr == 0) then
                call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, &
                  lr=em_vars%lr(idir, sigma_alt, ifactor))
              end if
              call restart_close_dir(restart_load)

              if(ierr /= 0) then
                message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, sigma))//&
                   "': Initializing to zero."
                call messages_warning(1)
              end if
            end if
          
            if(em_vars%calc_hyperpol .and. use_kdotp) then
              do idir2 = 1, space%periodic_dim
                str_tmp = em_wfs_tag(idir, ifactor, idir2)              
                call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, sigma), ierr)
                if (ierr == 0) then
                  call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, &
                    lr=kdotp_em_lr2(idir2, idir, sigma_alt, ifactor))
                end if
                call restart_close_dir(restart_load)
              
                if(ierr /= 0) then
                  message(1) = "Unable to read second-order response wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, sigma))//&
                     "': Initializing to zero."
                  call messages_warning(1)
                end if
              end do
            end if
            if((pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC .or. em_vars%calc_magnetooptics) .and. use_kdotp) then
              if (sigma == 1 .and. ifactor == 1) then 
                do ipert = PK2, PKB
                  do idir2 = 1, space%dim 
                    if((ipert == PKB) .or. (idir2 <= idir)) then
                      str_tmp = em_wfs_tag(idir, ifactor, idir2, ipert)              
                      call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, sigma), ierr)
                      if (ierr == 0) then
                        select case(ipert)
                        case(PK2)
                          call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, &
                            lr = k2_lr(idir, idir2, sigma))
                        case(PKB)
                          call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, &
                            lr = kb_lr(idir, idir2, sigma))
                        end select
                      end if
                      call restart_close_dir(restart_load)
                    end if
                    if(ierr /= 0) then
                      select case(ipert)
                      case(PK2)
                        message(1) = "Unable to read magneto-optics response wavefunctions (K2) from '"&
                          //trim(wfs_tag_sigma(str_tmp, sigma))// "': Initializing to zero."
                      case(PKB)
                        message(1) = "Unable to read magneto-optics response wavefunctions (KB) from '"&
                          //trim(wfs_tag_sigma(str_tmp, sigma))// "': Initializing to zero."
                      end select
                      call messages_warning(1)
                    end if
                  end do
                end do
              end if
            end if
            if(em_vars%calc_magnetooptics) then
              if(sigma == 1 .and. ifactor == 1) then
                str_tmp = em_wfs_tag(idir, ifactor, ipert = PB)  
                call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, sigma), ierr)
                if(ierr == 0) then
                  call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, lr = b_lr(idir, sigma))
                end if
                call restart_close_dir(restart_load)
                if(ierr /= 0) then
                  message(1) = "Unable to read magneto-optics response wavefunctions (B) from '"&
                    //trim(wfs_tag_sigma(str_tmp, sigma))//"': Initializing to zero."
                  call messages_warning(1)
                end if
              end if
              if(use_kdotp) then
                if(nfactor_ke > 1 .or. ifactor == 1) then
                do idir2 = 1, space%dim
                  str_tmp = em_wfs_tag(idir, ifactor, idir2, ipert = PKE)              
                  call restart_open_dir(restart_load, wfs_tag_sigma(str_tmp, sigma), ierr)
                  if (ierr == 0) then
                    call states_elec_load(restart_load, namespace, space, st, gr%mesh, kpoints, ierr, &
                      lr=ke_lr(idir, idir2, sigma_alt, ifactor))
                  end if
                  call restart_close_dir(restart_load)
                  if(ierr /= 0) then
                    message(1) = "Unable to read magneto-optics response wavefunctions (KE) from '"&
                      //trim(wfs_tag_sigma(str_tmp, sigma))// "': Initializing to zero."
                    call messages_warning(1)
                  end if
                end do
                end if
              end if
            end if
          end if
        end do
      end if
    
      ! if opposite sign from last frequency, swap signs to get a closer frequency
      if(iomega > 1 .and. em_vars%nsigma == 2) then
        if(em_vars%omega(iomega - 1) * em_vars%omega(iomega) < M_ZERO) then
          call X(lr_swap_sigma)(st, gr%mesh, em_vars%lr(idir, 1, ifactor), em_vars%lr(idir, 2, ifactor))
        end if
      end if
    
      if(allocate_rho_em .and. (em_vars%calc_hyperpol .or. (iomega == 1 .and. ifactor == 1))) then
        !search for the density of the closest frequency, including negative
        closest_omega = em_vars%freq_factor(ifactor) * em_vars%omega(iomega)
        call loct_search_file_lr(closest_omega, idir, ierr_e(idir), trim(restart_dir(restart_load)))
    
        !attempt to read 
        if(ierr_e(idir) == 0) then 
          sigma_alt = 1
          if(closest_omega * frequency < M_ZERO) opp_freq = .true.
          if(opp_freq .and. em_vars%nsigma == 2) sigma_alt = 2
      
          call X(lr_load_rho)(em_vars%lr(idir, sigma_alt, ifactor)%X(dl_rho), space, gr%mesh, st%d%nspin, restart_load, &
            em_rho_tag(closest_omega, idir), ierr_e(idir))
      
          if (ierr_e(idir) == 0) then 
            message(1) = "Read response density '"//trim(em_rho_tag(closest_omega, idir))//"'."
          else
            message(1) = "Unable to read response density '"//trim(em_rho_tag(closest_omega, idir))//"'."
          end if
          call messages_info(1)

          if (ierr_e(idir) == 0 .and. &
            abs(abs(closest_omega) - abs(frequency)) <= CNST(1e-4)) then
            ! the frequencies are written to four decimals in the restart directory, so we cannot expect higher precision
            exact_freq(idir) = .true.
          end if
        end if
    
        if (ierr_e(idir) == 0 .and. em_vars%nsigma == 2) then 
          sigma_alt = 1
          if(opp_freq) sigma_alt = 2
      
          em_vars%lr(idir, swap_sigma(sigma_alt), ifactor)%X(dl_rho) = R_CONJ(em_vars%lr(idir, sigma_alt, ifactor)%X(dl_rho))
        end if
       
        if(ierr_e(idir) == 0 .and. em_vars%calc_magnetooptics .and. use_kdotp) then
          if(complex_wfs) then
            call X(lr_load_rho)(em_vars%lr(idir, swap_sigma(sigma_alt), 2)%X(dl_rho), space, gr%mesh, st%d%nspin, &
              restart_load, em_rho_tag(closest_omega, idir, ipert = PE), ierr_e2(idir))
      
            if (ierr_e2(idir) == 0) then 
              message(1) = "Read response density '"//trim(em_rho_tag(closest_omega, idir, ipert = PE))//"'."
            else
              message(1) = "Unable to read response density '"//trim(em_rho_tag(closest_omega, idir, ipert = PE))//"'."
            end if
            call messages_info(1)

            if (ierr_e2(idir) == 0 .and. em_vars%nsigma == 2) then       
              em_vars%lr(idir, sigma_alt, 2)%X(dl_rho) = R_CONJ(em_vars%lr(idir, swap_sigma(sigma_alt), 2)%X(dl_rho))
            end if
          end if
        end if
      end if
    end do

    call restart_end(restart_load)
  end if ! .not. fromscratch
  
  if(em_vars%calc_magnetooptics .and. ifactor == 2) then
    do idir = 1, space%dim
      do sigma = 1, em_vars%nsigma
        sigma_alt = swap_sigma(sigma) 
        if((ierr /= 0) .or. (.not. complex_wfs)) then
          em_vars%lr(idir, sigma_alt, 2)%X(dl_psi) = R_CONJ(em_vars%lr(idir, sigma, 1)%X(dl_psi))
        end if
        if(allocate_rho_em .and. ((ierr_e2(idir) /= 0) .or. (.not. complex_wfs))) then
          em_vars%lr(idir, sigma_alt, 2)%X(dl_rho) = R_CONJ(em_vars%lr(idir, sigma, 1)%X(dl_rho))
        end if
      end do
    end do
  end if
  ! if the frequency is zero, we do not need to calculate both responses
  if(abs(frequency) < M_EPSILON .and. em_vars%nsigma == 2) then
    nsigma_eff = 1
  else
    nsigma_eff = em_vars%nsigma
  end if
  SAFE_ALLOCATE(inhomog(1:gr%mesh%np, 1:hm%d%dim, 1:st%nst, 1:st%d%kpt%nlocal, 1:nsigma_eff))

  if((pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) .and. use_kdotp) then
    do idir = 1, space%dim 
      message(1)="Info: Calculating response for B-perturbation"
      call messages_info(1)
      call X(inhomog_B)(sh, namespace, gr, st, hm, xc, ions, magn_dir(idir,1), magn_dir(idir,2), kdotp_lr(magn_dir(idir, 1), 1:1), &
        kdotp_lr(magn_dir(idir, 2), 1:1), inhomog)
      call X(sternheimer_set_inhomog)(sh, inhomog)   
      call X(sternheimer_solve)(sh, namespace, space, gr, kpoints, st, hm, xc, mc, ions, em_vars%lr(idir, 1:1, ifactor), 1, &
        R_TOPREC(frequency_zero), pert2_none, restart_dump, &
        em_rho_tag(abs(em_vars%freq_factor(ifactor) * em_vars%omega(iomega)), idir), &
        em_wfs_tag(idir, ifactor), have_restart_rho = (ierr_e(idir) == 0), have_exact_freq = exact_freq(idir))
      call sternheimer_unset_inhomog(sh)
      em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
    end do

    do idir = 1, space%dim
      do idir2 = 1, idir
        message(1)="Info: Calculating response for K2-perturbation"
        call messages_info(1)
        call X(inhomog_k2_tot)(namespace, gr, st, hm, ions, idir, idir2, kdotp_lr(idir, 1:1), kdotp_lr(idir2, 1:1), inhomog)
        call X(sternheimer_set_inhomog)(sh_kmo, inhomog)
        call X(sternheimer_solve)(sh_kmo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, k2_lr(idir, idir2, 1:1), 1, &
          R_TOPREC(frequency_zero), pert2_none, restart_dump, "null", &
          em_wfs_tag(idir, ifactor, idir2, PK2), have_restart_rho = .false., have_exact_freq = .false.)
        call sternheimer_unset_inhomog(sh_kmo)
        em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_kmo)
      end do
    end do

    do idir = 1, space%dim
      do idir2 = 1, space%dim    
        message(1)="Info: Calculating response for KB-perturbation"
        call messages_info(1)
        call X(inhomog_kB_tot)(sh, namespace, gr, st, hm, xc, ions, idir, magn_dir(idir2, 1), &
          magn_dir(idir2, 2), kdotp_lr(idir, 1:1), em_vars%lr(idir2, 1:1, ifactor), kdotp_lr(magn_dir(idir2, 1), 1:1), &
          kdotp_lr(magn_dir(idir2, 2), 1:1), & 
          k2_lr(max(magn_dir(idir2, 1), idir), min(magn_dir(idir2,1), idir), 1:1),&
          k2_lr(max(magn_dir(idir2, 2), idir), min(magn_dir(idir2,2), idir), 1:1), inhomog)
        call X(sternheimer_set_inhomog)(sh_kmo, inhomog)   
        call X(sternheimer_solve)(sh_kmo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, &
          kb_lr(idir, idir2, 1:1), 1, R_TOPREC(frequency_zero), pert2_none, restart_dump, "null", &
          em_wfs_tag(idir, ifactor, idir2, PKB), have_restart_rho = .false., have_exact_freq = .false.)
        call sternheimer_unset_inhomog(sh_kmo)
        em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_kmo)
      end do
    end do
  else 
    do idir = 1, space%dim
      if((.not. em_vars%calc_magnetooptics) .or. ifactor == 1 .or. complex_wfs) then
        call pert_setup_dir(em_vars%perturbation, idir)
  
        if(use_kdotp .and. idir <= space%periodic_dim) then
          do ik = st%d%kpt%start, st%d%kpt%end
            ik0 = ik - st%d%kpt%start + 1
            do ist = st%st_start, st%st_end
              inhomog(1:gr%mesh%np, 1:st%d%dim, ist, ik0, 1) = kdotp_lr(idir, 1)%X(dl_psi)(1:gr%mesh%np, 1:st%d%dim, ist, ik)
            end do
          end do
          call X(sternheimer_set_rhs)(sh, inhomog(:,:,:,:,1))
        end if

        str_tmp = freq2str(units_from_atomic(units_out%energy, frequency))
        write(message(1), '(5a)') 'Info: Calculating response for the ', index2axis(idir), &
          '-direction and frequency ', trim(str_tmp), '.'
        call messages_info(1)
  
        if((.not. em_vars%calc_magnetooptics) .or. ifactor==1) then
          call X(sternheimer_solve)(sh, namespace, space, gr, kpoints, st, hm, xc, mc, ions, &
            em_vars%lr(idir, 1:nsigma_eff, ifactor), nsigma_eff, R_TOPREC(frequency_eta), &
            em_vars%perturbation, restart_dump, em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir), &
            em_wfs_tag(idir, ifactor), have_restart_rho=(ierr_e(idir)==0), have_exact_freq = exact_freq(idir))
        else
          call X(sternheimer_solve)(sh, namespace, space, gr, kpoints, st, hm, xc, mc, ions, &
            em_vars%lr(idir, 1:nsigma_eff, ifactor), nsigma_eff, R_TOPREC(frequency_eta), &
            em_vars%perturbation, restart_dump, &
            em_rho_tag(abs(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)), idir, ipert = PE), &
            em_wfs_tag(idir, ifactor), have_restart_rho=.true., have_exact_freq = exact_freq(idir))
        end if
  
        if(nsigma_eff == 1 .and. em_vars%nsigma == 2) then
          em_vars%lr(idir, 2, ifactor)%X(dl_psi) = em_vars%lr(idir, 1, ifactor)%X(dl_psi)
          if(allocate_rho_em) em_vars%lr(idir, 2, ifactor)%X(dl_rho) = R_CONJ(em_vars%lr(idir, 1, ifactor)%X(dl_rho))
        end if
  
        if(use_kdotp) then
          call sternheimer_unset_rhs(sh)
        end if
  
        em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
      end if
  
      if(em_vars%calc_hyperpol .and. use_kdotp) then
        call X(em_resp_calc_eigenvalues)(space, gr%mesh, gr%sb, st, dl_eig)

        call restart_init(kdotp_restart, namespace, RESTART_KDOTP, RESTART_TYPE_LOAD, mc, ierr, mesh=gr%mesh)

        call lr_init(kdotp_lr2)
        call lr_allocate(kdotp_lr2, st, gr%mesh, allocate_rho = .false.)

        do idir2 = 1, space%periodic_dim
          write(message(1), '(a,a,a)') 'Info: Calculating kdotp response in ', index2axis(idir2), '-direction.'
          call messages_info(1)
          call pert_setup_dir(pert_kdotp, idir2)
      
          message(1) = "Reading 2nd-order kdotp wavefunction."
          call messages_info(1)

          ! load wavefunctions
          str_tmp = kdotp_wfs_tag(min(idir, idir2), max(idir, idir2))
          ! 1 is the sigma index which is used in em_resp
          call restart_open_dir(kdotp_restart, wfs_tag_sigma(str_tmp, 1), ierr)
          if (ierr == 0) then
            call states_elec_load(kdotp_restart, namespace, space, st, gr%mesh, kpoints, ierr, lr=kdotp_lr2)
          end if

          call restart_close_dir(kdotp_restart)
          if(ierr /= 0) then
            message(1) = "Unable to read second-order kdotp wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'."
            message(2) = "Previous kdotp calculation (with KdotPCalcSecondOrder) required."
            call messages_fatal(2)
          end if
          
          call X(sternheimer_solve_order2)(sh, sh_kdotp, sh2, namespace, space, gr, kpoints, st, hm, xc, mc, ions, &
            em_vars%lr(idir, 1:nsigma_eff, ifactor), kdotp_lr(idir2, 1:1), nsigma_eff, &
            R_TOPREC(frequency_eta), R_TOTYPE(M_ZERO), em_vars%perturbation, pert_kdotp, &
            kdotp_em_lr2(idir2, idir, 1:nsigma_eff, ifactor), pert2_none, &
            restart_dump, "null", em_wfs_tag(idir, ifactor, idir2), have_restart_rho=.true., have_exact_freq = .true., &
            give_pert1psi2 = kdotp_lr2%X(dl_psi), give_dl_eig1 = dl_eig(:, :, idir2))
      
          ! if the frequency is zero, we do not need to calculate both responses
          if(nsigma_eff == 1 .and. em_vars%nsigma == 2) then
            kdotp_em_lr2(idir2, idir, 2, ifactor)%X(dl_psi) = kdotp_em_lr2(idir2, idir, 1, ifactor)%X(dl_psi)
          end if
      
          em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
        end do
        write(message(1), '(a)') ''
        call messages_info(1)

        call restart_end(kdotp_restart)
        call lr_dealloc(kdotp_lr2)
      end if
    end do
  end if

  if(em_vars%calc_magnetooptics) then 
    if(use_kdotp) then
      if(iomega == 1 .and. ifactor == 1) then
        do idir = 1, space%dim
          message(1) = "Info: Calculating response for B-perturbation"
          call messages_info(1)  
          call X(inhomog_B)(sh_mo, namespace, gr, st, hm, xc, ions, magn_dir(idir, 1), &
            magn_dir(idir, 2), kdotp_lr(magn_dir(idir, 1), 1:1), kdotp_lr(magn_dir(idir, 2), 1:1), inhomog)
          call X(sternheimer_set_inhomog)(sh_mo, inhomog)   
          call X(sternheimer_solve)(sh_mo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, b_lr(idir, 1:1), 1, &
            R_TOPREC(frequency_zero), pert2_none, restart_dump, "null", &
            em_wfs_tag(idir, ifactor, ipert = PB), have_restart_rho = .false., have_exact_freq = .false.)
          call sternheimer_unset_inhomog(sh_mo)
          em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_mo)
        end do
      end if   

      do ipert = PK2, PKE
        do idir2 = 1, space%dim
          do idir = 1, space%dim 
            select case(ipert)
            case(PK2) 
              if(iomega == 1 .and. ifactor == 1) then 
                if(idir2 <= idir) then
                  message(1) = "Info: Calculating response for K2-perturbation"
                  call messages_info(1)
                  call X(inhomog_k2_tot)(namespace, gr, st, hm, ions, idir, idir2, kdotp_lr(idir, 1:1), kdotp_lr(idir2, 1:1), inhomog)
                  call X(sternheimer_set_inhomog)(sh_kmo, inhomog)
                  call X(sternheimer_solve)(sh_kmo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, k2_lr(idir, idir2, 1:1), &
                    1, R_TOPREC(frequency_zero), pert2_none, restart_dump, "null", &
                    em_wfs_tag(idir, ifactor, idir2, ipert), have_restart_rho = .false., have_exact_freq = .false.)
                  call sternheimer_unset_inhomog(sh_kmo)
                  em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_kmo)
                end if
              end if
                
            case(PKB)
              if(iomega == 1 .and. ifactor == 1) then
                message(1) = "Info: Calculating response for KB-perturbation"
                call messages_info(1) 
                call X(inhomog_kB_tot)(sh_mo, namespace, gr, st, hm, xc, ions, idir, magn_dir(idir2, 1), magn_dir(idir2, 2), &
                  kdotp_lr(idir, 1:1), b_lr(idir2, 1:1), kdotp_lr(magn_dir(idir2, 1), 1:1), kdotp_lr(magn_dir(idir2, 2), 1:1), &
                  k2_lr(max(magn_dir(idir2, 1), idir), min(magn_dir(idir2, 1), idir), 1:1), & 
                  k2_lr(max(magn_dir(idir2, 2), idir), min(magn_dir(idir2, 2), idir), 1:1), inhomog)
                call X(sternheimer_set_inhomog)(sh_kmo, inhomog)   
                call X(sternheimer_solve)(sh_kmo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, kb_lr(idir, idir2, 1:1), &
                  1, R_TOPREC(frequency_zero), pert2_none, restart_dump, "null", &
                  em_wfs_tag(idir, ifactor, idir2, ipert), have_restart_rho = .false., have_exact_freq = .false.)
                call sternheimer_unset_inhomog(sh_kmo)
                em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_kmo)
              end if
                  
            case(PKE)
              if (ifactor .le. nfactor_ke) then
                message(1) = "Info: Calculating response for KE-perturbation"
                call messages_info(1) 
                call X(inhomog_kE_tot)(sh, namespace, gr, st, hm, xc, ions, idir, nsigma_eff, kdotp_lr(idir, 1:1), &
                  em_vars%lr(idir2, :, ifactor), k2_lr(max(idir2, idir), min(idir2,idir), 1:1), inhomog)
                call X(sternheimer_set_inhomog)(sh_kmo, inhomog)   
                call X(sternheimer_solve)(sh_kmo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, &
                  ke_lr(idir, idir2, 1:nsigma_eff, ifactor), nsigma_eff, &
                  R_TOPREC(frequency_eta), pert2_none, restart_dump, "null", &
                  em_wfs_tag(idir, ifactor, idir2, ipert), have_restart_rho = .false., have_exact_freq = .false.)
                if(nsigma_eff == 1 .and. em_vars%nsigma == 2) then
                  ke_lr(idir, idir2, 2, ifactor)%X(dl_psi) = ke_lr(idir, idir2, 1, ifactor)%X(dl_psi)
                end if
                call sternheimer_unset_inhomog(sh_kmo)
                em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_kmo)
              end if
            end select
          end do
        end do
      end do
    else
      if(iomega == 1 .and. ifactor == 1) then
        do idir = 1, space%dim
          message(1)="Info: Calculating response for B-perturbation"
          call messages_info(1)  
          call pert_setup_dir(pert_b, idir)
          call X(sternheimer_solve)(sh_mo, namespace, space, gr, kpoints, st, hm, xc, mc, ions, b_lr(idir, 1:1), 1, &
            R_TOPREC(frequency_zero), pert_b, restart_dump, "null", &
            em_wfs_tag(idir, ifactor, ipert = PB), have_restart_rho = .false., have_exact_freq = .false.)
          em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh_mo)
        end do
      end if
    end if
  end if

  call restart_end(restart_dump)
  SAFE_DEALLOCATE_A(inhomog) 

  POP_SUB(em_resp_run_legacy.X(run_sternheimer))
end subroutine X(run_sternheimer)

! ---------------------------------------------------------
subroutine X(calc_properties_linear)(em_vars, namespace, space, gr, kpoints, st, hm, xc, ions, outp)
  type(em_resp_t),          intent(inout) :: em_vars
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(kpoints_t),          intent(in)    :: kpoints
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(inout) :: hm
  type(xc_t),               intent(in)    :: xc
  type(ions_t),             intent(in)    :: ions
  type(output_t),           intent(in)    :: outp
  
  integer :: idir, idir1, idir2

  PUSH_SUB(em_resp_run_legacy.X(calc_properties_linear))
  
  if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
    if((.not. em_vars%calc_magnetooptics) .or. ifactor == 1) then 
      ! calculate polarizability
      message(1) = "Info: Calculating polarizabilities."
      call messages_info(1)
    
      if(use_kdotp) then
        if(.not. em_vars%kpt_output) then
          call X(calc_polarizability_periodic)(space, gr%mesh, gr%symm, st, em_vars%lr(:, :, ifactor), &
            kdotp_lr(:, 1), em_vars%nsigma, em_vars%alpha(:, :, ifactor))
        else
          call X(calc_polarizability_periodic)(space, gr%mesh, gr%symm, st, em_vars%lr(:, :, ifactor), &
            kdotp_lr(:, 1), em_vars%nsigma, em_vars%alpha(:, :, ifactor), zpol_k = em_vars%alpha_k(:, :, ifactor, :))
        end if
        if(em_vars%lrc_kernel) then
          do idir = 1, space%dim
            do idir2 = 1, space%dim
              lrc_coef(idir, idir2) = M_ONE/(M_ONE - M_HALF * xc%kernel_lrc_alpha * &
                (em_vars%alpha(idir, idir, ifactor) + em_vars%alpha(idir2, idir2, ifactor)) / ions%latt%rcell_volume)
            end do
          end do
          do idir = 1, space%dim
            do idir2 = 1, space%dim
              em_vars%alpha0(idir, idir2, ifactor) = em_vars%alpha(idir, idir2, ifactor)
              em_vars%alpha(idir, idir2, ifactor) = em_vars%alpha(idir, idir2, ifactor) * lrc_coef(idir, idir2)
            end do
          end do 
        end if
      end if

      call X(calc_polarizability_finite)(namespace, space, gr, st, hm, ions, em_vars%lr(:, :, ifactor), &
        em_vars%nsigma, em_vars%perturbation, em_vars%alpha(:, :, ifactor), doalldirs = .not. use_kdotp)
    
      if(em_vars%calc_Born) then
        ! calculate Born effective charges
        message(1) = "Info: Calculating (frequency-dependent) Born effective charges."
        call messages_info(1)
      
        call X(forces_born_charges)(gr, namespace, space, ions, hm%ep, st, kpoints, &
          lr = em_vars%lr(:, 1, ifactor), lr2 = em_vars%lr(:, em_vars%nsigma, ifactor), &
          Born_charges = em_vars%Born_charges(ifactor), lda_u_level= hm%lda_u_level)
      end if

    else
    
      write(message(1), '(a)') 'Info: Calculating magneto-optical response.'
      call messages_info(1)
      em_vars%alpha_be(:, :, :) = M_ZERO
      em_vars%chi_para(:, :) = M_ZERO
      em_vars%chi_dia(:, :) = M_ZERO
    
      if(use_kdotp) then
        if(abs(frequency) > M_EPSILON) then
          if(.not. em_vars%kpt_output) then
            call X(lr_calc_magneto_optics_periodic)(sh, sh_mo, namespace, space, gr, st, hm, xc, ions, &
              em_vars%nsigma, em_vars%nfactor, nfactor_ke, em_vars%freq_factor, em_vars%lr(:, :, :), b_lr(:, :), &
              kdotp_lr(:, :), ke_lr(:, :, :, :), kb_lr(:, :, :), -frequency_eta, em_vars%alpha_be(:, :, :))
          else
            call X(lr_calc_magneto_optics_periodic)(sh, sh_mo, namespace, space, gr, st, hm, xc, ions, &
              em_vars%nsigma, em_vars%nfactor, nfactor_ke, em_vars%freq_factor, em_vars%lr(:, :, :), b_lr(:, :), &
              kdotp_lr(:, :), ke_lr(:, :, :, :), kb_lr(:, :, :), -frequency_eta, em_vars%alpha_be(:, :, :), &
              zpol_kout = em_vars%alpha_be_k(:, :, :, :))
          end if
          if(em_vars%lrc_kernel) then
            do idir = 1, space%dim
              do idir1 = 1, space%dim
                do idir2 = 1, space%dim
                  em_vars%alpha_be0(idir, idir1, idir2) = em_vars%alpha_be(idir, idir1, idir2)
                  em_vars%alpha_be(idir, idir1, idir2) = em_vars%alpha_be(idir, idir1, idir2) * lrc_coef(idir, idir1)
                end do
              end do
            end do
          end if
        end if
        if(iomega == 1) then
          message(1) = "Info: Calculating magnetic susceptibilities."
          call messages_info(1)
          call X(lr_calc_susceptibility_periodic)(namespace, space, gr%symm, gr%mesh, st, hm, &
            kdotp_lr(:, 1), b_lr(:, 1), k2_lr(:, :, 1), kb_lr(:, :, 1), em_vars%chi_dia(:, :))
          em_vars%chi_para(:, :) = M_ZERO  
          call X(lr_calc_magnetization_periodic)(namespace, space, gr%mesh, st, hm, kdotp_lr(:, 1), em_vars%magn(:))
        end if
      else
        call X(lr_calc_magneto_optics_finite)(sh, sh_mo, namespace, space, gr, st, hm, xc, ions, &
          em_vars%nsigma, em_vars%nfactor, em_vars%lr(:, :, :), b_lr(:, :), em_vars%alpha_be(:, :, :))
        call X(lr_calc_susceptibility)(namespace, space, gr, st, hm, ions, b_lr(:, :), 1, pert_b, &
          em_vars%chi_para(: ,:), em_vars%chi_dia(:, :))
      end if
    end if
    
  else if(pert_type(em_vars%perturbation) == PERTURBATION_MAGNETIC) then
    message(1) = "Info: Calculating magnetic susceptibilities."
    call messages_info(1)
  
    if(use_kdotp) then
      call X(lr_calc_susceptibility_periodic)(namespace, space, gr%symm, gr%mesh, st, hm, &
        kdotp_lr(:, 1), em_vars%lr(:, 1, ifactor), k2_lr(:, :, 1), kb_lr(:, :, 1), em_vars%chi_dia(:, :))
        em_vars%chi_para(:, :) = M_ZERO
      call X(lr_calc_magnetization_periodic)(namespace, space, gr%mesh, st, hm, kdotp_lr(:, 1), em_vars%magn(:))
    else
      call X(lr_calc_susceptibility)(namespace, space, gr, st, hm, ions, em_vars%lr(:, :, ifactor), &
        em_vars%nsigma, em_vars%perturbation, em_vars%chi_para(:, :), em_vars%chi_dia(:, :))
    end if
  end if
  
  call em_resp_output(st, namespace, space, gr, hm, ions, outp, sh, em_vars, iomega, ifactor)
  
  POP_SUB(em_resp_run_legacy.X(calc_properties_linear))
end subroutine X(calc_properties_linear)

! ---------------------------------------------------------
subroutine X(calc_properties_nonlinear)(em_vars, namespace, space, gr, st, hm, xc, ions)
  type(em_resp_t),          intent(inout) :: em_vars
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(states_elec_t),      intent(in)    :: st
  type(hamiltonian_elec_t), intent(inout) :: hm
  type(xc_t),               intent(in)    :: xc
  type(ions_t),             intent(in)    :: ions

  character(len=100) :: dirname_output, str_tmp

  PUSH_SUB(em_resp_run_legacy.X(calc_properties_nonlinear))

  ! calculate hyperpolarizability
  if(em_vars%calc_hyperpol) then
    write(message(1), '(a)') 'Info: Calculating hyperpolarizabilities.'
    call messages_info(1)
    
    if(use_kdotp) then
      call X(post_orthogonalize)(space, gr%mesh, st, em_vars%nfactor, em_vars%nsigma, em_vars%freq_factor(:), &
        em_vars%omega(iomega), em_vars%eta, em_vars%lr, kdotp_em_lr2)
      call X(lr_calc_beta)(sh, namespace, space, gr, st, hm, xc, ions, em_vars%lr, em_vars%perturbation, em_vars%beta, &
        kdotp_lr = kdotp_lr(:, 1), kdotp_em_lr = kdotp_em_lr2, dl_eig = dl_eig, occ_response = .false.)
    else
      call X(lr_calc_beta)(sh, namespace, space, gr, st, hm, xc, ions, em_vars%lr, em_vars%perturbation, em_vars%beta, &
        occ_response = em_vars%occ_response)
    end if
    
    str_tmp = freq2str(units_from_atomic(units_out%energy, em_vars%freq_factor(1)*em_vars%omega(iomega)))
    write(dirname_output, '(a, a)') EM_RESP_DIR//'freq_', trim(str_tmp)
    call io_mkdir(trim(dirname_output), namespace)
    call out_hyperpolarizability(gr%sb, em_vars%beta, em_vars%freq_factor(1:3), em_vars%ok(1), dirname_output, namespace)
  end if

  POP_SUB(em_resp_run_legacy.X(calc_properties_nonlinear))
end subroutine X(calc_properties_nonlinear)
