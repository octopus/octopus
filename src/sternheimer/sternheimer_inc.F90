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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: sternheimer_inc.F90 2663 2007-01-25 09:04:29Z lorenzen $

!--------------------------------------------------------------
!> This routine calculates the first-order variations of the wavefunctions 
!! for an applied perturbation.
subroutine X(sternheimer_solve)(                           &
     this, sys, hm, lr, nsigma, omega, perturbation,       &
     restart_dir, rho_tag, wfs_tag, have_restart_rho, have_exact_freq)
  type(sternheimer_t),    intent(inout) :: this
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr(:) 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(in)    :: omega
  type(pert_t),           intent(in)    :: perturbation
  character(len=*),       intent(in)    :: restart_dir
  character(len=*),       intent(in)    :: rho_tag
  character(len=*),       intent(in)    :: wfs_tag
  logical,      optional, intent(in)    :: have_restart_rho
  logical,      optional, intent(in)    :: have_exact_freq

  FLOAT :: dpsimod, tol
  integer :: iter, sigma, sigma_alt, ik, ist, err
  R_TYPE, allocatable :: dl_rhoin(:, :, :), dl_rhonew(:, :, :), dl_rhotmp(:, :, :)
  R_TYPE, allocatable :: rhs(:, :, :), hvar(:, :, :), psi(:, :)
  R_TYPE, allocatable :: tmp(:)
  real(8):: abs_dens
  R_TYPE :: omega_sigma, proj
  logical, allocatable :: orth_mask(:)

  logical :: conv_last, conv, states_conv, have_restart_rho_
  type(mesh_t), pointer :: mesh
  type(states_t), pointer :: st
  integer :: total_iter, idim, ip, ispin
  character(len=100) :: dirname

  PUSH_SUB(X(sternheimer_solve))
  call profiling_in(prof, "STERNHEIMER")

  ASSERT(nsigma == 1 .or. nsigma == 2)

  mesh => sys%gr%mesh
  st => sys%st
  
  call mix_clear(this%mixer, func_type = R_TYPE_VAL)

  call mesh_init_mesh_aux(sys%gr%mesh)
  
  SAFE_ALLOCATE(tmp(1:mesh%np))
  SAFE_ALLOCATE(rhs(1:mesh%np, 1:st%d%dim, 1:nsigma))
  SAFE_ALLOCATE(hvar(1:mesh%np, 1:st%d%nspin, 1:nsigma))
  SAFE_ALLOCATE(dl_rhoin(1:mesh%np, 1:st%d%nspin, 1:1))
  SAFE_ALLOCATE(dl_rhonew(1:mesh%np, 1:st%d%nspin, 1:1))
  SAFE_ALLOCATE(dl_rhotmp(1:mesh%np, 1:st%d%nspin, 1:1))
  SAFE_ALLOCATE(orth_mask(1:st%nst))

  conv = .false.
  conv_last = .not. (this%add_fxc .or. this%add_hartree) .or. .not. this%last_occ_response
  ! otherwise it is not actually SCF, and there can only be one pass through

  have_restart_rho_ = .false.
  if(present(have_restart_rho)) have_restart_rho_ = have_restart_rho
  if(.not. have_restart_rho_) call X(lr_build_dl_rho)(mesh, st, lr, nsigma)
  
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
      call X(lr_orth_response)(mesh, st, lr(sigma), omega_sigma)
    enddo
  endif

  !this call is required to reset the scf_tol object, whether we want its result or not
  tol = scf_tol_step(this%scf_tol, 0, M_ONE)
  if(have_restart_rho_ .and. present(have_exact_freq)) then
    if(have_exact_freq) then
      tol = scf_tol_final(this%scf_tol) * M_TEN
      ! if rho is converged already, then we should try to solve fully for the wavefunctions
    endif
  endif

  !self-consistency iteration for response
  iter_loop: do iter = 1, this%scf_tol%max_iter
    if (this%add_fxc .or. this%add_hartree) then
      write(message(1), '(a, i3)') "LR SCF Iteration: ", iter
      call messages_info(1)
    endif

    write(message(1), '(a, f20.6, a, f20.6, a, i1)') &
         "Frequency: ", units_from_atomic(units_out%energy,  R_REAL(omega)), &
         " Eta : ",     units_from_atomic(units_out%energy, R_AIMAG(omega))
    write(message(2), '(a)') &
         '   ik  ist                norm   iters            residual'
    call messages_info(2)

    do ispin = 1, st%d%nspin
       call lalg_copy(mesh%np, lr(1)%X(dl_rho)(:, ispin), dl_rhoin(:, ispin, 1))
    end do

    call X(sternheimer_calc_hvar)(this, sys, lr, nsigma, hvar)

    SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

    do ik = st%d%kpt%start, st%d%kpt%end
      !now calculate response for each state
      ispin = states_dim_get_spin_index(sys%st%d, ik)

      states_conv = .true.

      do ist = st%st_start, st%st_end

        call states_get_state(sys%st, sys%gr%mesh, ist, ik, psi)

        do sigma = 1, nsigma
          !calculate the RHS of the Sternheimer eq
          if(sternheimer_have_rhs(this)) then
            ASSERT(associated(this%X(rhs)))
            forall(idim = 1:st%d%dim, ip = 1:mesh%np) rhs(ip, idim, sigma) = this%X(rhs)(ip, idim, ist, ik)
          else
            rhs(1:mesh%np, 1:st%d%dim, sigma) = R_TOTYPE(M_ZERO)
            call X(pert_apply)(perturbation, sys%gr, sys%geo, hm, ik, psi, rhs(:, :, sigma))
          end if

          do idim = 1, st%d%dim
            rhs(1:mesh%np, idim, sigma) = -rhs(1:mesh%np, idim, sigma) &
              - hvar(1:mesh%np, ispin, sigma)*psi(1:mesh%np, idim)
          end do

          if(sternheimer_have_inhomog(this)) then
            forall(idim = 1:st%d%dim, ip = 1:mesh%np)
              rhs(ip, idim, sigma) = rhs(ip, idim, sigma) + this%X(inhomog)(ip, idim, ist, ik, sigma)
            end forall
          endif

          ! Let Pc = projector onto unoccupied states, Pn` = projector that removes state n
          ! For an SCF run, we will apply Pn` for the last step always, since the whole wavefunction is useful for some
          ! things and the extra cost here is small. If occ_response, previous steps will also use Pn`. If !occ_response,
          ! previous steps will use Pc, which generally reduces the number of linear-solver iterations needed. Only the
          ! wavefunctions in the unoccupied subspace are needed to construct the first-order density.
          ! I am not sure what the generalization of this scheme is for metals, so we will just use Pc if there is smearing.

          if(sigma == 1) then 
            omega_sigma = omega
          else 
            omega_sigma = -R_CONJ(omega)
          end if

          if (conv_last .and. this%last_occ_response) then
            ! project out only the component of the unperturbed wavefunction
            proj = X(mf_dotp)(mesh, st%d%dim, psi, rhs(:, :, sigma))
            do idim = 1, st%d%dim
              call lalg_axpy(mesh%np, -proj, psi(:, idim), rhs(:, idim, sigma))
            end do
          else
            ! project RHS onto the unoccupied states
            call X(lr_orth_vector)(mesh, st, rhs(:, :, sigma), ist, ik, omega_sigma)
          endif

          !solve the Sternheimer equation
          call X(solve_HXeY)(this%solver, hm, sys%gr, sys%st, ist, ik, &
               lr(sigma)%X(dl_psi)(1:mesh%np_part, 1:st%d%dim, ist, ik), &
               rhs(:, :, sigma), -sys%st%eigenval(ist, ik) + omega_sigma, tol, this%occ_response)

          if (this%preorthogonalization) then 
            !re-orthogonalize the resulting vector
            if (this%occ_response) then
              proj = X(mf_dotp)(mesh, st%d%dim, psi, lr(sigma)%X(dl_psi)(:, :, ist, ik))
              do idim = 1, st%d%dim
                call lalg_axpy(mesh%np, -proj, psi(:, idim), lr(sigma)%X(dl_psi)(:, idim, ist, ik))
              end do
            else
              call X(lr_orth_vector)(mesh, st, lr(sigma)%X(dl_psi)(1:mesh%np_part, 1:st%d%dim, ist, ik), ist, ik, omega_sigma)
            endif
          end if

          ! print the norm of the variations, and the number of
          ! iterations and residual of the linear solver
          dpsimod = X(mf_nrm2)(mesh, st%d%dim, lr(sigma)%X(dl_psi)(:, :, ist, ik))

          write(message(1), '(i5, i5, f20.6, i8, e20.6)') &
            ik, (3 - 2 * sigma) * ist, dpsimod, this%solver%iter, this%solver%abs_psi 
          call messages_info(1)

          states_conv = states_conv .and. (this%solver%abs_psi < tol)
          total_iter = total_iter + this%solver%iter
          
        end do !sigma
      end do !ist
    end do !ik
    
    SAFE_DEALLOCATE_A(psi)
    
    call X(lr_build_dl_rho)(mesh, st, lr, nsigma)
    
    dl_rhonew(1:mesh%np, 1:st%d%nspin, 1) = M_ZERO
    
    if(this%add_fxc .or. this%add_hartree) then
      !write restart info
      !save all frequencies as positive
      if(R_REAL(omega) >= M_ZERO) then
        sigma_alt = 1
      else
        sigma_alt = 2
      endif
      call X(restart_write_lr_rho)(lr(sigma_alt), sys%gr, st%d%nspin, restart_dir, rho_tag)
    endif

    do sigma = 1, nsigma 
      !save all frequencies as positive
      sigma_alt = sigma
      if(R_REAL(omega) < M_ZERO) sigma_alt = swap_sigma(sigma)

      write(dirname,'(2a)') trim(restart_dir), trim(wfs_tag_sigma(wfs_tag, sigma_alt))
      call restart_write(trim(tmpdir)//dirname, st, sys%gr, sys%geo, err, iter = iter, lr = lr(sigma))
    end do
    
    if (.not. states_conv) then
      message(1) = "Linear solver failed to converge all states."
      call messages_warning(1)
    endif

    if (.not.(this%add_fxc .or. this%add_hartree)) then
    ! no need to deal with mixing, SCF iterations, etc.
    ! dealing with restart density above not necessary, but easier to just leave it
    ! convergence criterion is now about individual states, rather than SCF residual
      this%ok = states_conv

      message(1)="--------------------------------------------"
      write(message(2), '(a, i8)') &
           'Info: Total Hamiltonian applications:', total_iter * linear_solver_ops_per_iter(this%solver)
      call messages_info(2)
      exit
    endif

    ! all the rest is the mixing and checking for convergence

    if(this%scf_tol%max_iter == iter) then 
      message(1) = "Self-consistent iteration for response did not converge."
      this%ok = .false.
      call messages_warning(1)
    end if

    do ispin = 1, st%d%nspin
      call lalg_copy(mesh%np, lr(1)%X(dl_rho)(:, ispin), dl_rhotmp(:, ispin, 1))
    end do

    call X(mixing)(this%mixer, iter, dl_rhoin, dl_rhotmp, dl_rhonew, X(mf_dotp_aux))

    abs_dens = M_ZERO

    do ispin = 1, st%d%nspin
      forall(ip = 1:mesh%np) tmp(ip) = dl_rhoin(ip, ispin, 1) - dl_rhotmp(ip, ispin, 1)
      abs_dens = hypot(abs_dens, real(X(mf_nrm2)(mesh, tmp), 8))
    end do

    write(message(1), '(a, e20.6)') "SCF Residual ", abs_dens

    message(2)="--------------------------------------------"
    call messages_info(2)
      
    if( abs_dens <= this%scf_tol%conv_abs_dens ) then 
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
           '      Total Hamiltonian applications:', total_iter * linear_solver_ops_per_iter(this%solver)
      call messages_info(2)
      exit
    else
      ! not quitting if converged allows results to be calculated if possible
      ! before dying on the next direction or frequency
      if(clean_stop()) then
        message(1) = "Exiting cleanly."
        call messages_fatal(1, only_root_writes = .true.)
      endif

      do ispin = 1, st%d%nspin
        call lalg_copy(mesh%np, dl_rhonew(:, ispin, 1), lr(1)%X(dl_rho)(:, ispin))
      end do
      
      if(nsigma == 2) then
        ! we do it in this way to avoid a bug in ifort 11.1
        dl_rhonew(:, :, 1) = R_CONJ(dl_rhonew(:, :, 1))
        do ispin = 1, st%d%nspin
          call lalg_copy(mesh%np, dl_rhonew(:, ispin, 1), lr(2)%X(dl_rho)(:, ispin))
        end do
      end if
      
      tol = scf_tol_step(this%scf_tol, iter, TOFLOAT(abs_dens))
    end if
    
  end do iter_loop

  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(rhs)
  SAFE_DEALLOCATE_A(hvar)
  SAFE_DEALLOCATE_A(dl_rhoin)
  SAFE_DEALLOCATE_A(dl_rhonew)
  SAFE_DEALLOCATE_A(dl_rhotmp)
  SAFE_DEALLOCATE_A(orth_mask)

  call profiling_out(prof)
  POP_SUB(X(sternheimer_solve))

end subroutine X(sternheimer_solve)


!--------------------------------------------------------------
subroutine X(sternheimer_calc_hvar)(this, sys, lr, nsigma, hvar)
  type(sternheimer_t),    intent(inout) :: this
  type(system_t),         intent(inout) :: sys
  type(lr_t),             intent(inout) :: lr(:) 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(out)   :: hvar(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:nsigma)

  R_TYPE, allocatable :: tmp(:), hartree(:)
  integer :: np, ip, ispin, ispin2

  PUSH_SUB(X(sternheimer_calc_hvar))
  call profiling_in(prof_hvar, "STERNHEIMER_HVAR")

  np = sys%gr%mesh%np

  if (this%add_hartree) then 
    SAFE_ALLOCATE(    tmp(1:np))
    SAFE_ALLOCATE(hartree(1:np))
    do ip = 1, np
      tmp(ip) = sum(lr(1)%X(dl_rho)(ip, 1:sys%st%d%nspin))
    end do
    hartree(1:np) = R_TOTYPE(M_ZERO)
    call X(poisson_solve)(psolver, hartree, tmp, all_nodes = .false.)

    SAFE_DEALLOCATE_A(tmp)
  end if

  do ispin = 1, sys%st%d%nspin
    !* initialize
    hvar(1:np, ispin, 1) = M_ZERO

    !* hartree
    if (this%add_hartree) hvar(1:np, ispin, 1) = hvar(1:np, ispin, 1) + hartree(1:np)
    
    !* fxc
    if(this%add_fxc) then
      do ispin2 = 1, sys%st%d%nspin
        hvar(1:np, ispin, 1) = hvar(1:np, ispin, 1) + this%fxc(1:np, ispin, ispin2)*lr(1)%X(dl_rho)(1:np, ispin2)
      end do
    end if
  end do
  
  if (nsigma == 2) hvar(1:np, 1:sys%st%d%nspin, 2) = R_CONJ(hvar(1:np, 1:sys%st%d%nspin, 1))

  if (this%add_hartree) then
    SAFE_DEALLOCATE_A(hartree)
  end if

  call profiling_out(prof_hvar)
  POP_SUB(X(sternheimer_calc_hvar))
end subroutine X(sternheimer_calc_hvar)


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
subroutine X(sternheimer_solve_order2)( &
     sh1, sh2, sh_2ndorder, sys, hm, lr1, lr2, nsigma, omega1, omega2, pert1, pert2,       &
     lr_2ndorder, pert_2ndorder, restart_dir, rho_tag, wfs_tag, have_restart_rho, have_exact_freq, &
     give_pert1psi2)
  type(sternheimer_t),    intent(inout) :: sh1
  type(sternheimer_t),    intent(inout) :: sh2
  type(sternheimer_t),    intent(inout) :: sh_2ndorder
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr1(:) 
  type(lr_t),             intent(inout) :: lr2(:) 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(in)    :: omega1
  R_TYPE,                 intent(in)    :: omega2
  type(pert_t),           intent(in)    :: pert1
  type(pert_t),           intent(in)    :: pert2
  type(lr_t),             intent(inout) :: lr_2ndorder(:)
  type(pert_t),           intent(in)    :: pert_2ndorder
  character(len=*),       intent(in)    :: restart_dir
  character(len=*),       intent(in)    :: rho_tag
  character(len=*),       intent(in)    :: wfs_tag
  logical,      optional, intent(in)    :: have_restart_rho
  logical,      optional, intent(in)    :: have_exact_freq
  R_TYPE,       optional, intent(in)    :: give_pert1psi2(:,:,:,:)

  integer :: isigma, ik, ist, idim, ispin
  R_TYPE :: dl_eig1, dl_eig2
  R_TYPE, allocatable :: inhomog(:,:,:,:,:), hvar1(:,:,:), hvar2(:,:,:), &
    pert1psi2(:,:), pert2psi1(:,:), pert1psi(:,:), pert2psi(:,:)
  R_TYPE, allocatable :: psi(:, :)
  type(mesh_t), pointer :: mesh
  type(states_t), pointer :: st

  PUSH_SUB(X(sternheimer_solve_order2))

  ASSERT(nsigma == 1 .or. nsigma == 2)

  mesh => sys%gr%mesh
  st => sys%st

  SAFE_ALLOCATE(inhomog(1:mesh%np, 1:st%d%dim, 1:st%nst, 1:st%d%nik, 1:nsigma))
  SAFE_ALLOCATE(hvar1(1:mesh%np, 1:st%d%nspin, 1:nsigma))
  SAFE_ALLOCATE(hvar2(1:mesh%np, 1:st%d%nspin, 1:nsigma))
  SAFE_ALLOCATE(pert1psi2(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert2psi1(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert1psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pert2psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

  call X(sternheimer_calc_hvar)(sh1, sys, lr1, nsigma, hvar1)
!  call X(sternheimer_calc_hvar)(sh2, sys, lr2, nsigma, hvar2)
! for kdotp, hvar = 0
  hvar2 = M_ZERO

  inhomog(:,:,:,:,:) = M_ZERO

  do isigma = 1, nsigma

    do ik = st%d%kpt%start, st%d%kpt%end

      ispin = states_dim_get_spin_index(sys%st%d, ik)
      do ist = st%st_start, st%st_end

        call states_get_state(st, sys%gr%mesh, ist, ik, psi)

        call X(pert_apply)(pert1, sys%gr, sys%geo, hm, ik, psi, pert1psi)
        call X(pert_apply)(pert2, sys%gr, sys%geo, hm, ik, psi, pert2psi)
        if(present(give_pert1psi2)) then
          pert1psi2(:,:) = give_pert1psi2(:, :, ist, ik)
        else
          call X(pert_apply)(pert1, sys%gr, sys%geo, hm, ik, lr2(isigma)%X(dl_psi)(:, :, ist, ik), pert1psi2)
        endif
        call X(pert_apply)(pert2, sys%gr, sys%geo, hm, ik, lr1(isigma)%X(dl_psi)(:, :, ist, ik), pert2psi1)

        ! derivative of the eigenvalues
        dl_eig1 = X(mf_dotp)(mesh, st%d%dim, psi, pert1psi)
        dl_eig2 = X(mf_dotp)(mesh, st%d%dim, psi, pert2psi)

        do idim = 1, st%d%dim
          dl_eig1 = dl_eig1 + X(mf_dotp)(mesh, R_TOTYPE(abs(psi(:, idim))**2), hvar1(:, ispin, isigma))
          dl_eig2 = dl_eig2 + X(mf_dotp)(mesh, R_TOTYPE(abs(psi(:, idim))**2), hvar2(:, ispin, isigma))
        enddo

!        write(message(1),*) 'dl_eig1 ist ', ist, 'ik ', ik, dl_eig1
!        write(message(2),*) 'dl_eig2 ist ', ist, 'ik ', ik, dl_eig2
!        call messages_info(2)
        
        do idim = 1, st%d%dim
          inhomog(1:mesh%np, idim, ist, ik, isigma) = &
            - pert1psi2(1:mesh%np, idim) &
            - (hvar1(1:mesh%np, ispin, isigma) - dl_eig1) * lr2(isigma)%X(dl_psi)(1:mesh%np, idim, ist, ik) &
            - pert2psi1(1:mesh%np, idim) &
            - (hvar2(1:mesh%np, ispin, isigma) - dl_eig2) * lr1(isigma)%X(dl_psi)(1:mesh%np, idim, ist, ik)
        end do

      enddo
    enddo
  enddo

  inhomog(:,:,:,:,:) = inhomog(:,:,:,:,:) * M_HALF

  ! sum frequency
  call X(sternheimer_set_inhomog)(sh_2ndorder, inhomog)
  call X(sternheimer_solve)(sh_2ndorder, sys, hm, lr_2ndorder, nsigma, &
    omega1 + omega2, pert_2ndorder, restart_dir, rho_tag, wfs_tag, &
    have_restart_rho = have_restart_rho, have_exact_freq = have_exact_freq)
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
