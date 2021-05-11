!! Copyright (C) 2005-2012 M. Marques, X. Andrade, David Strubbe
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

!--------------------------------------------------------------
!> This routine calculates the first-order variations of the wavefunctions 
!! for an applied perturbation.
subroutine X(sternheimer_solve)(this, namespace, space, gr, kpoints, st, hm, xc, mc, ions, lr, nsigma, omega, perturbation, restart, &
  rho_tag, wfs_tag, have_restart_rho, have_exact_freq)
  type(sternheimer_t),         intent(inout) :: this
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),       target,  intent(in)    :: gr
  type(kpoints_t),             intent(in)    :: kpoints
  type(states_elec_t),         intent(in)    :: st
  type(hamiltonian_elec_t),    intent(inout) :: hm
  type(xc_t),                  intent(in)    :: xc
  type(multicomm_t),           intent(in)    :: mc
  type(ions_t),                intent(in)    :: ions
  type(lr_t),                  intent(inout) :: lr(:) 
  integer,                     intent(in)    :: nsigma 
  R_TYPE,                      intent(in)    :: omega
  type(pert_t),                intent(in)    :: perturbation
  type(restart_t),             intent(inout) :: restart
  character(len=*),            intent(in)    :: rho_tag
  character(len=*),            intent(in)    :: wfs_tag
  logical,           optional, intent(in)    :: have_restart_rho
  logical,           optional, intent(in)    :: have_exact_freq

  FLOAT :: tol
  FLOAT, allocatable :: dpsimod(:, :), residue(:, :)
  integer, allocatable :: conv_iters(:, :)
  integer :: iter, sigma, sigma_alt, ik, ist, err, sst, est, ii, ierr
  R_TYPE, allocatable :: dl_rhoin(:, :, :), dl_rhonew(:, :, :), dl_rhotmp(:, :, :)
  R_TYPE, allocatable :: rhs(:, :, :), hvar(:, :, :), psi(:, :), rhs_full(:, :, :)
  R_TYPE, allocatable :: tmp(:), rhs_tmp(:, :, :)
  FLOAT  :: abs_dens, rel_dens
  R_TYPE :: omega_sigma, proj
  logical, allocatable :: orth_mask(:)
  type(wfs_elec_t) :: rhsb, dlpsib
  logical :: conv_last, conv, states_conv, have_restart_rho_
  integer :: total_iter, idim, ip, ispin, ib, total_iter_reduced
  logical :: calculate_rho
#ifdef HAVE_MPI
  logical :: states_conv_reduced
#endif
  
  PUSH_SUB(X(sternheimer_solve))
  call profiling_in(prof, TOSTRING(X(STERNHEIMER)))

  ASSERT(nsigma == 1 .or. nsigma == 2)

  call mix_clear(this%mixer)
  
  calculate_rho = this%add_fxc .or. this%add_hartree

  if (st%d%ispin == SPINORS .and. calculate_rho) then
    call messages_not_implemented('linear response density for spinors', namespace=namespace)
  end if

  SAFE_ALLOCATE(dpsimod(1:nsigma, st%st_start:st%st_end))
  SAFE_ALLOCATE(residue(1:nsigma, st%st_start:st%st_end))
  SAFE_ALLOCATE(conv_iters(1:nsigma, st%st_start:st%st_end))
  SAFE_ALLOCATE(tmp(1:gr%mesh%np))
  SAFE_ALLOCATE(rhs(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))
  SAFE_ALLOCATE(rhs_tmp(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))
  if(this%last_occ_response .and. .not. this%occ_response_by_sternheimer) then
    SAFE_ALLOCATE(rhs_full(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))
  end if
  if(calculate_rho) then
    SAFE_ALLOCATE(hvar(1:gr%mesh%np, 1:st%d%nspin, 1:nsigma))
    SAFE_ALLOCATE(dl_rhoin(1:gr%mesh%np, 1:st%d%nspin, 1:1))
    SAFE_ALLOCATE(dl_rhonew(1:gr%mesh%np, 1:st%d%nspin, 1:1))
    SAFE_ALLOCATE(dl_rhotmp(1:gr%mesh%np, 1:st%d%nspin, 1:1))
  end if
  SAFE_ALLOCATE(orth_mask(1:st%nst))

  conv = .false.
  conv_last = .not. calculate_rho .or. .not. this%last_occ_response
  ! otherwise it is not actually SCF, and there can only be one pass through

  have_restart_rho_ = .false.
  if(present(have_restart_rho)) have_restart_rho_ = have_restart_rho
  if((.not. have_restart_rho_) .and. calculate_rho) call X(lr_build_dl_rho)(gr%mesh, st, lr, nsigma)

  message(1)="--------------------------------------------"
  call messages_info(1)

  total_iter = 0

  ! preorthogonalization
  if (this%preorthogonalization) then 
    do sigma = 1, nsigma
      if(sigma == 1) then 
        omega_sigma = omega
      else
        omega_sigma = -R_CONJ(omega)
      end if
      call X(lr_orth_response)(gr%mesh, st, lr(sigma), omega_sigma)
    end do
  end if

  !this call is required to reset the scf_tol object, whether we want its result or not
  tol = scf_tol_step(this%scf_tol, 0, M_ONE)
  if(have_restart_rho_ .and. present(have_exact_freq)) then
    if(have_exact_freq) then
      tol = scf_tol_final(this%scf_tol) * CNST(10.0)
      ! if rho is converged already, then we should try to solve fully for the wavefunctions
    end if
  end if

  !self-consistency iteration for response
  iter_loop: do iter = 1, this%scf_tol%max_iter
    if (calculate_rho) then
      write(message(1), '(a, i3)') "LR SCF Iteration: ", iter
      call messages_info(1)
    end if

    write(message(1), '(a, f20.6, a, f20.6, a, i1)') &
      "Frequency: ", units_from_atomic(units_out%energy,  R_REAL(omega)), &
      " Eta : ",     units_from_atomic(units_out%energy, R_AIMAG(omega))
    write(message(2), '(a)') &
      '   ik  ist                norm   iters            residual'
    call messages_info(2)

    if (calculate_rho) then
      call lalg_copy(gr%mesh%np, st%d%nspin, lr(1)%X(dl_rho)(:, :), dl_rhoin(:, :, 1))

      this%X(omega) = omega
      call X(sternheimer_calc_hvar)(this, gr%mesh, st, hm, xc, lr, nsigma, hvar)
    end if

    SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

    states_conv = .true.

    do ik = st%d%kpt%start, st%d%kpt%end
      !now calculate response for each state
      ispin = st%d%get_spin_index(ik)

      do ib = st%group%block_start, st%group%block_end
        
        sst = states_elec_block_min(st, ib)
        est = states_elec_block_max(st, ib)

        !calculate the RHS of the Sternheimer eq

        call wfs_elec_init(rhsb, st%d%dim, sst, est, rhs_tmp, ik)

        if(sternheimer_have_rhs(this)) then
          do ist = sst, est
            call batch_set_state(rhsb, ist-sst+1, gr%mesh%np, this%X(rhs)(:, :, ist, ik - st%d%kpt%start + 1))
          end do
        else
          call X(pert_apply_batch)(perturbation, namespace, gr, ions, hm, st%group%psib(ib, ik), rhsb)
        end if

        call rhsb%end()

        ii = 0
        do ist = sst, est
          ii = ii + 1

          do idim = 1, st%d%dim
            rhs_tmp(1:gr%mesh%np, idim, ii) = -rhs_tmp(1:gr%mesh%np, idim, ii)
          end do
        end do
        
        do sigma = 1, nsigma

          if(sigma == 1) then 
            omega_sigma = omega
          else 
            omega_sigma = -R_CONJ(omega)
          end if

          ii = 0
          do ist = sst, est
            ii = ii + 1

            call lalg_copy(gr%mesh%np, st%d%dim, rhs_tmp(:,:,ii), rhs(:,:,ii))
            
            if(calculate_rho) then
              call states_elec_get_state(st, gr%mesh, ist, ik, psi)
              do idim = 1, st%d%dim
                rhs(1:gr%mesh%np, idim, ii) = rhs(1:gr%mesh%np, idim, ii) - &
                  hvar(1:gr%mesh%np, ispin, sigma)*psi(1:gr%mesh%np, idim)
              end do
            end if

            if(sternheimer_have_inhomog(this)) then
              call lalg_axpy(gr%mesh%np, st%d%dim, M_ONE, this%X(inhomog)(:, :, ist, ik-st%d%kpt%start + 1, sigma), rhs(:, :, ii))
            end if

            if(conv_last .and. this%last_occ_response .and. .not. this%occ_response_by_sternheimer) then
              call lalg_copy(gr%mesh%np, st%d%dim, rhs(:, :, ii), rhs_full(:, :, ii))
            end if

            call X(lr_orth_vector)(gr%mesh, st, rhs(:, :, ii), ist, ik, omega_sigma, &
              min_proj = conv_last .and. this%last_occ_response .and. this%occ_response_by_sternheimer)

          end do

          !solve the Sternheimer equation
          call X(wfs_elec_init)(dlpsib, st%d%dim, sst, est, gr%mesh%np_part, ik)
          do ist = sst, est
            call batch_set_state(dlpsib, ist-sst+1, gr%mesh%np, lr(sigma)%X(dl_psi)(:, :, ist, ik))
          end do
          call wfs_elec_init(rhsb, st%d%dim, sst, est, rhs, ik)

          call X(linear_solver_solve_HXeY_batch)(this%solver, namespace, hm, gr%mesh, st, dlpsib, rhsb, &
            -st%eigenval(sst:est, ik) + omega_sigma, tol, residue(sigma, sst:est), conv_iters(sigma, sst:est), &
            occ_response = this%occ_response)

          do ist = sst, est
            call batch_get_state(dlpsib, ist-sst+1, gr%mesh%np, lr(sigma)%X(dl_psi)(:, :, ist, ik))
          end do
          call dlpsib%end()
          call rhsb%end()

          !re-orthogonalize the resulting vector
          ii = 0
          do ist = sst, est
            ii = ii + 1

            if (this%preorthogonalization) then 
              ! should remove degenerate states here too
              if (this%occ_response) then
                call states_elec_get_state(st, gr%mesh, ist, ik, psi)
                proj = X(mf_dotp)(gr%mesh, st%d%dim, psi, lr(sigma)%X(dl_psi)(:, :, ist, ik))
                call lalg_axpy(gr%mesh%np, st%d%dim,-proj, psi, lr(sigma)%X(dl_psi)(:, :, ist, ik))
              else
                call X(lr_orth_vector)(gr%mesh, st, lr(sigma)%X(dl_psi)(1:gr%mesh%np_part, 1:st%d%dim, ist, ik), ist, ik, &
                  omega_sigma)
              end if
            end if

            dpsimod(sigma, ist) = X(mf_nrm2)(gr%mesh, st%d%dim, lr(sigma)%X(dl_psi)(:, :, ist, ik))

          end do !ist

          if(conv_last .and. this%last_occ_response .and. .not. this%occ_response_by_sternheimer) then
            call X(sternheimer_add_occ)(gr%mesh, st, lr(sigma)%X(dl_psi)(:, :, :, ik), rhs_full, sst, est, ik, omega_sigma, &
              CNST(1e-5))
          end if

        end do !sigma

        do ist = sst, est
          do sigma = 1, nsigma
            states_conv = states_conv .and. (residue(sigma, ist) < tol)
            total_iter = total_iter + conv_iters(sigma, ist)

            write(message(1), '(i5, i5, f20.6, i8, e20.6)') &
              ik, (3 - 2*sigma)*ist, dpsimod(sigma, ist), conv_iters(sigma, ist), residue(sigma, ist)
            call messages_info(1)
          end do !sigma
        end do !ist

      end do !sst
    end do !ik

    SAFE_DEALLOCATE_A(psi)

    total_iter_reduced = total_iter
#ifdef HAVE_MPI
    if(st%d%kpt%parallel) then
      call MPI_Allreduce(states_conv, states_conv_reduced, 1, MPI_LOGICAL, MPI_LAND, st%d%kpt%mpi_grp%comm, mpi_err)
      states_conv = states_conv_reduced

      call MPI_Allreduce(total_iter, total_iter_reduced, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
    end if
#endif

    if(calculate_rho) then
      call X(lr_build_dl_rho)(gr%mesh, st, lr, nsigma)

      dl_rhonew(1:gr%mesh%np, 1:st%d%nspin, 1) = M_ZERO

      !write restart info
      !save all frequencies as positive
      if(R_REAL(omega) >= M_ZERO) then
        sigma_alt = 1
      else
        sigma_alt = 2
      end if
      if(st%d%kpt%start == 1) then
        call X(lr_dump_rho)(lr(sigma_alt), space, gr%mesh, st%d%nspin, restart, rho_tag, ierr)
      end if
      if (ierr /= 0) then
        message(1) = "Unable to write response density '"//trim(rho_tag)//"'."
        call messages_warning(1)
      end if
    end if

    do sigma = 1, nsigma 
      !save all frequencies as positive
      sigma_alt = sigma
      if(R_REAL(omega) < M_ZERO) sigma_alt = swap_sigma(sigma)

      call restart_open_dir(restart, wfs_tag_sigma(wfs_tag, sigma_alt), err)
      if (err == 0) then
        call states_elec_dump(restart, space, st, gr%mesh, kpoints, err, iter = iter, lr = lr(sigma))
      end if
      if (err /= 0) then
        message(1) = "Unable to write response wavefunctions."
        call messages_warning(1)
      end if
      call restart_close_dir(restart)
    end do

    if (.not. states_conv) then
      message(1) = "Linear solver failed to converge all states."
      call messages_warning(1)
    end if

    if (.not. calculate_rho) then
      ! no need to deal with mixing, SCF iterations, etc.
      ! dealing with restart density above not necessary, but easier to just leave it
      ! convergence criterion is now about individual states, rather than SCF residual
      this%ok = states_conv

      message(1)="--------------------------------------------"
      write(message(2), '(a, i8)') &
        'Info: Total Hamiltonian applications:', total_iter_reduced * linear_solver_ops_per_iter(this%solver)
      call messages_info(2)
      exit
    end if

    ! all the rest is the mixing and checking for convergence

    if(this%scf_tol%max_iter == iter) then 
      message(1) = "Self-consistent iteration for response did not converge."
      this%ok = .false.
      call messages_warning(1)
    end if

    call lalg_copy(gr%mesh%np, st%d%nspin, lr(1)%X(dl_rho)(:, :), dl_rhotmp(:, :, 1))

    call X(mixing)(this%mixer, dl_rhoin, dl_rhotmp, dl_rhonew)

    abs_dens = M_ZERO

    !NTD: This is quite different from the scf criterium.
    ! In the scf routine, we use a norm 1 to evaluate the change in density
    ! whereas here we evaluate a norm 2.
    ! For the spinor case, the present version is most likely not correct, as the off-diagonal term
    ! will not cancel each other at convergence and the normalization by the charge is thus not correct
    do ispin = 1, st%d%nspin
      do ip = 1, gr%mesh%np
        tmp(ip) = dl_rhoin(ip, ispin, 1) - dl_rhotmp(ip, ispin, 1)
      end do
      abs_dens = hypot(abs_dens, TOFLOAT(X(mf_nrm2)(gr%mesh, tmp)))
    end do
    rel_dens = abs_dens / st%qtot

    write(message(1), '(a,e16.6,a,e16.6,a)') "SCF Residual: ", abs_dens, " (abs), ", rel_dens, " (rel)"

    message(2)="--------------------------------------------"
    call messages_info(2)

    if( (abs_dens <= this%scf_tol%conv_abs_dens .or. this%scf_tol%conv_abs_dens <= M_ZERO) .and. &
      (rel_dens <= this%scf_tol%conv_rel_dens .or. this%scf_tol%conv_rel_dens <= M_ZERO)) then 
      if(conv_last .and. states_conv) then 
        ! if not all states are converged, keep working
        conv = .true.
      else
        conv_last = .true.
      end if
    end if

    if(conv) then
      this%ok = .true.

      write(message(1), '(a, i4, a)') &
        'Info: SCF for response converged in ', iter, ' iterations.'
      write(message(2), '(a, i8)') &
        '      Total Hamiltonian applications:', total_iter_reduced * linear_solver_ops_per_iter(this%solver)
      call messages_info(2)
      exit
    else
      ! not quitting if converged allows results to be calculated if possible
      ! before dying on the next direction or frequency
      if(clean_stop(mc%master_comm)) then
        message(1) = "Exiting cleanly."
        call messages_fatal(1, only_root_writes = .true.)
      end if

      call lalg_copy(gr%mesh%np, st%d%nspin, dl_rhonew(:, :, 1), lr(1)%X(dl_rho)(:, :))

      if(nsigma == 2) then
        ! we do it in this way to avoid a bug in ifort 11.1
        dl_rhonew(:, :, 1) = R_CONJ(dl_rhonew(:, :, 1))
        call lalg_copy(gr%mesh%np, st%d%nspin, dl_rhonew(:, :, 1), lr(2)%X(dl_rho)(:, :))
      end if

      tol = scf_tol_step(this%scf_tol, iter, TOFLOAT(abs_dens))
    end if

  end do iter_loop

  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(dpsimod)
  SAFE_DEALLOCATE_A(residue)
  SAFE_DEALLOCATE_A(conv_iters)
  SAFE_DEALLOCATE_A(rhs)
  SAFE_DEALLOCATE_A(rhs_tmp)
  if(this%last_occ_response .and. .not. this%occ_response_by_sternheimer) then
    SAFE_DEALLOCATE_A(rhs_full)
  end if
  if(calculate_rho) then
    SAFE_DEALLOCATE_A(hvar)
    SAFE_DEALLOCATE_A(dl_rhoin)
    SAFE_DEALLOCATE_A(dl_rhonew)
    SAFE_DEALLOCATE_A(dl_rhotmp)
  end if
  SAFE_DEALLOCATE_A(orth_mask)

  call profiling_out(prof)
  POP_SUB(X(sternheimer_solve))
end subroutine X(sternheimer_solve)


! ---------------------------------------------------------
!> add projection onto occupied states, by sum over states
subroutine X(sternheimer_add_occ)(mesh, st, lr_psi, rhs, sst, est, ik, omega_sigma, degen_thres)
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  R_TYPE,              intent(inout) :: lr_psi(:,:,:)
  R_TYPE,              intent(in)    :: rhs(:,:,:) !< (np, ndim, nst)
  integer,             intent(in)    :: sst !< start state
  integer,             intent(in)    :: est !< start state
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: omega_sigma
  FLOAT,               intent(in)    :: degen_thres

  integer :: ist, ist2, ii
  R_TYPE :: mtxel
  R_TYPE, allocatable :: psi(:,:)
  FLOAT :: ediff

  PUSH_SUB(X(sternheimer_add_occ))

  if(st%parallel_in_states) then
    call messages_not_implemented("sternheimer_add_occ parallel in states")
  end if

  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

  ii = 0
  ! iteration within states block
  do ist = sst, est
    ii = ii + 1

    do ist2 = 1, st%nst
      ediff = st%eigenval(ist, ik) - st%eigenval(ist2, ik)

      ! avoid dividing by zero below; these contributions are arbitrary anyway
      if (abs(ediff) < degen_thres) cycle

      ! the unoccupied subspace was handled by the Sternheimer equation
      if (st%occ(ist2, ik) < M_HALF) cycle

      call states_elec_get_state(st, mesh, ist2, ik, psi)
      mtxel = X(mf_dotp)(mesh, st%d%dim, psi, rhs(:, :, ii))

      call lalg_axpy(mesh%np, st%d%dim, mtxel/(ediff - omega_sigma), psi(1:mesh%np, 1:st%d%dim), &
        lr_psi(1:mesh%np, 1:st%d%dim, ist))

     ! need to get psi(ist) to do this. correct for a Hermitian operator, not for kdotp (which would need -mtxel)
!     lr%X(dl_psi)(:, :, ist2, ik) = lr%X(dl_psi)(:, :, ist2, ik) + &
!       st%X(psi)(:, :, ist, ik) * R_CONJ(mtxel) / (st%eigenval(ist2, ik) - st%eigenval(ist, ik) - omega_sigma)

    end do
  end do

  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(sternheimer_add_occ))
end subroutine X(sternheimer_add_occ)


!--------------------------------------------------------------
subroutine X(sternheimer_calc_hvar)(this, mesh, st, hm, xc, lr, nsigma, hvar)
  type(sternheimer_t),      intent(inout) :: this
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(in)    :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(xc_t),               intent(in)    :: xc
  type(lr_t),               intent(inout) :: lr(:) 
  integer,                  intent(in)    :: nsigma 
  R_TYPE,                   intent(out)   :: hvar(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:nsigma)

  PUSH_SUB(X(sternheimer_calc_hvar))

  if(this%add_fxc) then
    call X(calc_hvar)(this%add_hartree, mesh, st, hm, xc, lr(1)%X(dl_rho), nsigma, hvar, fxc = this%fxc)
  else
    call X(calc_hvar)(this%add_hartree, mesh, st, hm, xc, lr(1)%X(dl_rho), nsigma, hvar)
  end if
  if (this%enable_el_pt_coupling) then
#ifdef R_TCOMPLEX    
    call calc_hvar_photons(this, mesh, st, lr(1)%zdl_rho, nsigma, hvar)
#else
    ! Photons only available with complex states
    ASSERT(.false.)
#endif
  end if

  POP_SUB(X(sternheimer_calc_hvar))
end subroutine X(sternheimer_calc_hvar)


!--------------------------------------------------------------
subroutine X(calc_hvar)(add_hartree, mesh, st, hm, xc, lr_rho, nsigma, hvar, fxc)
  logical,                  intent(in)    :: add_hartree
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(in)    :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(xc_t),               intent(in)    :: xc 
  R_TYPE,                   intent(in)    :: lr_rho(:,:) !< (1:mesh%np, 1:st%d%nspin)
  integer,                  intent(in)    :: nsigma 
  R_TYPE,                   intent(out)   :: hvar(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:nsigma)
  FLOAT, optional,          intent(in)    :: fxc(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:st%d%nspin)

  R_TYPE, allocatable :: tmp(:), hartree(:)
  integer :: ip, ispin, ispin2
  FLOAT :: coeff_hartree

  PUSH_SUB(X(calc_hvar))
  call profiling_in(prof_hvar, TOSTRING(X(CALC_HVAR)))

  coeff_hartree = -xc%kernel_lrc_alpha / (M_FOUR * M_PI)
  if(add_hartree) coeff_hartree = coeff_hartree + 1

  if (abs(coeff_hartree) > M_EPSILON) then
    SAFE_ALLOCATE(    tmp(1:mesh%np))
    SAFE_ALLOCATE(hartree(1:mesh%np))
    do ip = 1, mesh%np
      tmp(ip) = sum(lr_rho(ip, 1:st%d%nspin))
    end do
    hartree(1:mesh%np) = R_TOTYPE(M_ZERO)
    call X(poisson_solve)(hm%psolver, hartree, tmp, all_nodes = .false.)

    SAFE_DEALLOCATE_A(tmp)
  end if

  do ispin = 1, st%d%nspin
    !* initialize
    hvar(1:mesh%np, ispin, 1) = M_ZERO

    !* hartree
   if (abs(coeff_hartree) > M_EPSILON) &
     hvar(1:mesh%np, ispin, 1) = hvar(1:mesh%np, ispin, 1) + coeff_hartree * hartree(1:mesh%np)
    
    !* fxc
    if(present(fxc)) then
      do ispin2 = 1, st%d%nspin
        hvar(1:mesh%np, ispin, 1) = hvar(1:mesh%np, ispin, 1) + fxc(1:mesh%np, ispin, ispin2)*lr_rho(1:mesh%np, ispin2)
      end do
    end if
  end do

  if (nsigma == 2) hvar(1:mesh%np, 1:st%d%nspin, 2) = R_CONJ(hvar(1:mesh%np, 1:st%d%nspin, 1))

  if (add_hartree) then
    SAFE_DEALLOCATE_A(hartree)
  end if

  call profiling_out(prof_hvar)
  POP_SUB(X(calc_hvar))
end subroutine X(calc_hvar)

!--------------------------------------------------------------
subroutine X(calc_kvar)(this, mesh, st, lr_rho1, lr_rho2, nsigma, kvar)
  type(sternheimer_t),    intent(inout) :: this
  type(mesh_t),           intent(in)    :: mesh
  type(states_elec_t),    intent(in)    :: st
  R_TYPE,                 intent(in)    :: lr_rho1(:,:) !< (1:mesh%np, 1:st%d%nspin)
  R_TYPE,                 intent(in)    :: lr_rho2(:,:) !< (1:mesh%np, 1:st%d%nspin)
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(out)   :: kvar(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:nsigma)

  integer :: ispin, ispin2, ispin3

  PUSH_SUB(X(calc_kvar))

  do ispin = 1, st%d%nspin
    kvar(1:mesh%np, ispin, 1) = M_ZERO

    if(this%add_fxc) then
      do ispin2 = 1, st%d%nspin
        do ispin3 = 1, st%d%nspin
          kvar(1:mesh%np, ispin, 1) = kvar(1:mesh%np, ispin, 1) + &
            this%kxc(1:mesh%np, ispin, ispin2, ispin3) * lr_rho1(1:mesh%np, ispin2) * lr_rho2(1:mesh%np, ispin3)
        end do
      end do
    end if
  end do
  
  if (nsigma==2) kvar(1:mesh%np,1:st%d%nspin,2) = R_CONJ(kvar(1:mesh%np,1:st%d%nspin,1))

  POP_SUB(X(calc_kvar))
end subroutine X(calc_kvar)

!-----------------------------------------------------------
subroutine X(sternheimer_set_rhs)(this, rhs)
  type(sternheimer_t), intent(inout) :: this
  R_TYPE, target,      intent(in)    :: rhs(:, :, :, :)

  PUSH_SUB(X(sternheimer_set_rhs))
  this%X(rhs) => rhs

  POP_SUB(X(sternheimer_set_rhs))
end subroutine X(sternheimer_set_rhs)


!-----------------------------------------------------------
subroutine X(sternheimer_set_inhomog)(this, inhomog)
  type(sternheimer_t), intent(inout) :: this
  R_TYPE, target,      intent(in)    :: inhomog(:, :, :, :, :)

  PUSH_SUB(X(sternheimer_set_inhomog))
  this%X(inhomog) => inhomog

  POP_SUB(X(sternheimer_set_inhomog))
end subroutine X(sternheimer_set_inhomog)


!--------------------------------------------------------------
subroutine X(sternheimer_solve_order2)(sh1, sh2, sh_2ndorder, namespace, space, gr, kpoints, st, hm, xc, mc, ions, lr1, lr2, &
  nsigma, omega1, omega2, pert1, pert2, lr_2ndorder, pert_2ndorder, restart, rho_tag, wfs_tag, have_restart_rho, have_exact_freq, &
  give_pert1psi2, give_dl_eig1)
  type(sternheimer_t),         intent(inout) :: sh1
  type(sternheimer_t),         intent(inout) :: sh2
  type(sternheimer_t),         intent(inout) :: sh_2ndorder
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(grid_t),        target, intent(in)    :: gr
  type(kpoints_t),             intent(in)    :: kpoints
  type(states_elec_t),         intent(in)    :: st
  type(hamiltonian_elec_t),    intent(inout) :: hm
  type(xc_t),                  intent(in)    :: xc
  type(multicomm_t),           intent(in)    :: mc
  type(ions_t),                intent(in)    :: ions
  type(lr_t),                  intent(inout) :: lr1(:) 
  type(lr_t),                  intent(inout) :: lr2(:) 
  integer,                     intent(in)    :: nsigma 
  R_TYPE,                      intent(in)    :: omega1
  R_TYPE,                      intent(in)    :: omega2
  type(pert_t),                intent(in)    :: pert1
  type(pert_t),                intent(in)    :: pert2
  type(lr_t),                  intent(inout) :: lr_2ndorder(:)
  type(pert_t),                intent(in)    :: pert_2ndorder
  type(restart_t),             intent(inout) :: restart
  character(len=*),            intent(in)    :: rho_tag
  character(len=*),            intent(in)    :: wfs_tag
  logical,           optional, intent(in)    :: have_restart_rho
  logical,           optional, intent(in)    :: have_exact_freq
  R_TYPE,            optional, intent(in)    :: give_pert1psi2(:,:,:,:) !< (np, ndim, ist, ik)
  FLOAT,             optional, intent(in)    :: give_dl_eig1(:,:) !< (nst, nk) expectation values of bare perturbation

  integer :: isigma, ik, ist, idim, ispin, ip
  R_TYPE :: dl_eig1, dl_eig2
  R_TYPE, allocatable :: inhomog(:,:,:,:,:), hvar1(:,:,:), hvar2(:,:,:), &
    pert1psi2(:,:), pert2psi1(:,:), pert1psi(:,:), pert2psi(:,:)
  R_TYPE, allocatable :: psi(:, :), psi2(:, :)
  type(mesh_t), pointer :: mesh

  PUSH_SUB(X(sternheimer_solve_order2))

  ASSERT(nsigma == 1 .or. nsigma == 2)

  mesh => gr%mesh

  ! FIXME: do not allocate ones we can point instead

  SAFE_ALLOCATE(inhomog(1:mesh%np, 1:st%d%dim, 1:st%nst, 1:st%d%kpt%nlocal, 1:nsigma))
  SAFE_ALLOCATE(hvar1(1:mesh%np, 1:st%d%nspin, 1:nsigma))
  SAFE_ALLOCATE(hvar2(1:mesh%np, 1:st%d%nspin, 1:nsigma))
  SAFE_ALLOCATE(pert1psi2(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert2psi1(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert1psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert2psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))

  hvar1 = M_ZERO
  if (sh1%add_fxc .or. sh1%add_hartree) then
    call X(sternheimer_calc_hvar)(sh1, mesh, st, hm, xc, lr1, nsigma, hvar1)
  end if
!  call X(sternheimer_calc_hvar)(sh2, mesh, st, hm, xc, lr2, nsigma, hvar2)
! for kdotp, hvar = 0
  hvar2 = M_ZERO

  inhomog(:,:,:,:,:) = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end

    ispin = st%d%get_spin_index(ik)
    do ist = st%st_start, st%st_end

      call states_elec_get_state(st, mesh, ist, ik, psi)
    
      do idim = 1, st%d%dim
        do ip = 1, mesh%np
          psi2(ip, idim) = R_CONJ(psi(ip, idim))*psi(ip, idim) 
        end do
      end do

      do isigma = 1, nsigma

        call X(pert_apply)(pert1, namespace, gr, ions, hm, ik, psi, pert1psi)
        call X(pert_apply)(pert2, namespace, gr, ions, hm, ik, psi, pert2psi)
        if(present(give_pert1psi2)) then
          pert1psi2(1:mesh%np, 1:st%d%dim) = give_pert1psi2(1:mesh%np, 1:st%d%dim, ist, ik)
        else
          call X(pert_apply)(pert1, namespace, gr, ions, hm, ik, lr2(isigma)%X(dl_psi)(:, :, ist, ik), pert1psi2)
        end if
        call X(pert_apply)(pert2, namespace, gr, ions, hm, ik, lr1(isigma)%X(dl_psi)(:, :, ist, ik), pert2psi1)

        ! derivative of the eigenvalues:
        ! bare perturbation
        if(present(give_dl_eig1)) then
          dl_eig1 = R_TOTYPE(give_dl_eig1(ist, ik))
        else
          dl_eig1 = X(mf_dotp)(mesh, st%d%dim, psi, pert1psi)
        end if
        dl_eig2 = X(mf_dotp)(mesh, st%d%dim, psi, pert2psi)

        ! Hxc perturbation
        do idim = 1, st%d%dim
          dl_eig1 = dl_eig1 + X(mf_dotp)(mesh, psi2(:, idim), hvar1(:, ispin, isigma))
          dl_eig2 = dl_eig2 + X(mf_dotp)(mesh, psi2(:, idim), hvar2(:, ispin, isigma))
        end do

!        write(message(1),*) 'dl_eig1 ist ', ist, 'ik ', ik, dl_eig1
!        write(message(2),*) 'dl_eig2 ist ', ist, 'ik ', ik, dl_eig2
!        call messages_info(2)

        ! FIXME: sort out proper isigma`s when not just freq = 0 for second perturbation?

        do idim = 1, st%d%dim
          inhomog(1:mesh%np, idim, ist, ik - st%d%kpt%start + 1, isigma) = &
            - pert1psi2(1:mesh%np, idim) &
            - (hvar1(1:mesh%np, ispin, isigma) - dl_eig1) * lr2(1)%X(dl_psi)(1:mesh%np, idim, ist, ik) &
            - pert2psi1(1:mesh%np, idim) &
            - (hvar2(1:mesh%np, ispin, isigma) - dl_eig2) * lr1(isigma)%X(dl_psi)(1:mesh%np, idim, ist, ik)
        end do

      end do
    end do
  end do

  inhomog(:,:,:,:,:) = inhomog(:,:,:,:,:) * M_HALF

  ! sum frequency
  call X(sternheimer_set_inhomog)(sh_2ndorder, inhomog)
  call X(sternheimer_solve)(sh_2ndorder, namespace, space, gr, kpoints, st, hm, xc, mc, ions, lr_2ndorder, nsigma, &
    omega1 + omega2, pert_2ndorder, restart, rho_tag, wfs_tag, have_restart_rho = have_restart_rho, &
    have_exact_freq = have_exact_freq)
  call sternheimer_unset_inhomog(sh_2ndorder)
  
  ! call X(sternheimer_solve) for difference frequency
  ! unless of course they are identical (one or both freqs the same)

  SAFE_DEALLOCATE_A(inhomog)
  SAFE_DEALLOCATE_A(hvar1)
  SAFE_DEALLOCATE_A(hvar2)
  SAFE_DEALLOCATE_A(pert1psi2)
  SAFE_DEALLOCATE_A(pert2psi1)
  SAFE_DEALLOCATE_A(pert1psi)
  SAFE_DEALLOCATE_A(pert2psi)
  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(sternheimer_solve_order2))
end subroutine X(sternheimer_solve_order2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
