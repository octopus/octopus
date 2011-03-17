!! Copyright (C) 2005-2006 M. Marques, X. Andrade
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
! This routine calculates the first-order variations of the wavefunctions 
! for an applied perturbation.
!--------------------------------------------------------------
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
  R_TYPE, allocatable :: rhs(:, :, :), hvar(:, :, :)
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
  conv_last = .false.

  have_restart_rho_ = .false.
  if(present(have_restart_rho)) have_restart_rho_ = have_restart_rho
  if(.not. have_restart_rho_) call X(lr_build_dl_rho)(mesh, st, lr, nsigma)
  
  message(1)="--------------------------------------------"
  call write_info(1)

  total_iter = 0

  ! preorthogonalization
  if (this%preorthogonalization) then 
    do sigma = 1, nsigma
      call X(lr_orth_response)(mesh, st, lr(sigma), omega)
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
      call write_info(1)
    else
      conv_last = .true.
    endif
    ! otherwise it is not actually SCF, and there can only be one pass through

    write(message(1), '(a, f20.6, a, f20.6, a, i1)') &
         "Frequency: ", units_from_atomic(units_out%energy,  R_REAL(omega)), &
         " Eta : ",     units_from_atomic(units_out%energy, R_AIMAG(omega))
    write(message(2), '(a)') &
         '   ik  ist                norm   iters            residual'
    call write_info(2)

    do ispin = 1, st%d%nspin
       call lalg_copy(mesh%np, lr(1)%X(dl_rho)(:, ispin), dl_rhoin(:, ispin, 1))
    end do

    call X(sternheimer_calc_hvar)(this, sys, hm, lr, nsigma, hvar)

    do ik = st%d%kpt%start, st%d%kpt%end
      !now calculate response for each state
      ispin = states_dim_get_spin_index(sys%st%d, ik)

      states_conv = .true.

      do ist = st%st_start, st%st_end
        do sigma = 1, nsigma
          !calculate the RHS of the Sternheimer eq
          if(sternheimer_have_rhs(this)) then
            ASSERT(associated(this%X(rhs)))
            forall(idim = 1:st%d%dim, ip = 1:mesh%np) rhs(ip, idim, sigma) = this%X(rhs)(ip, idim, ist, ik)
          else
            rhs(1:mesh%np, 1, sigma) = R_TOTYPE(M_ZERO)
            call X(pert_apply)(perturbation, sys%gr, sys%geo, hm, ik, st%X(psi)(:, :, ist, ik), rhs(:, :, sigma))
          end if

          do idim = 1, st%d%dim
            rhs(1:mesh%np, idim, sigma) = -rhs(1:mesh%np, idim, sigma) &
              - hvar(1:mesh%np, ispin, sigma)*st%X(psi)(1:mesh%np, idim, ist, ik)
          end do

          ! Let Pc = projector onto unoccupied states, Pn` = projector that removes state n
          ! For an SCF run, we will apply Pn` for the last step always, since the whole wavefunction is useful for some
          ! things and the extra cost here is small. If occ_response, previous steps will also use Pn`. If !occ_response,
          ! previous steps will use Pc, which generally reduces the number of linear-solver iterations needed. Only the
          ! wavefunctions in the unoccupied subspace are needed to construct the first-order density.
          ! I am not sure what the generalization of this scheme is for metals, so we will just use Pc if there is smearing.

          if (conv_last .and. this%last_occ_response) then
            ! project out only the component of the unperturbed wavefunction
            proj = X(mf_dotp)(mesh, st%d%dim, st%X(psi)(:, :, ist, ik), rhs(:, :, sigma))
            do idim = 1, st%d%dim
              call lalg_axpy(mesh%np, -proj, st%X(psi)(:, idim, ist, ik), rhs(:, idim, sigma))
            end do
          else
            ! project RHS onto the unoccupied states
            call X(lr_orth_vector)(mesh, st, rhs(:, :, sigma), ist, ik, omega)
          endif
        
          if(sigma == 1) then 
            omega_sigma = omega
          else 
            omega_sigma = -R_CONJ(omega)
          end if

          !solve the Sternheimer equation
          call X(solve_HXeY)(this%solver, hm, sys%gr, sys%st, ist, ik, &
               lr(sigma)%X(dl_psi)(1:mesh%np_part, 1:st%d%dim, ist, ik), &
               rhs(:, :, sigma), -sys%st%eigenval(ist, ik) + omega_sigma, tol, this%occ_response)

          if (this%preorthogonalization) then 
            !re-orthogonalize the resulting vector
            if (this%occ_response) then
              proj = X(mf_dotp)(mesh, st%d%dim, st%X(psi)(:, :, ist, ik), lr(sigma)%X(dl_psi)(:, :, ist, ik))
              do idim = 1, st%d%dim
                call lalg_axpy(mesh%np, -proj, st%X(psi)(:, idim, ist, ik), lr(sigma)%X(dl_psi)(:, idim, ist, ik))
              end do
            else
              call X(lr_orth_vector)(mesh, st, lr(sigma)%X(dl_psi)(1:mesh%np_part, 1:st%d%dim, ist, ik), ist, ik, omega)
            endif
          end if

          ! print the norm of the variations, and the number of
          ! iterations and residual of the linear solver
          dpsimod = X(mf_nrm2)(mesh, st%d%dim, lr(sigma)%X(dl_psi)(:, :, ist, ik))

          write(message(1), '(i5, i5, f20.6, i8, e20.6)') &
            ik, (3 - 2 * sigma) * ist, dpsimod, this%solver%iter, this%solver%abs_psi 
          call write_info(1)

          states_conv = states_conv .and. (this%solver%abs_psi < tol)
          total_iter = total_iter + this%solver%iter
          
        end do !sigma
      end do !ist
    end do !ik
    
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
      call restart_write(trim(tmpdir)//dirname, st, sys%gr, err, iter = iter, lr = lr(sigma))
    end do
    
    if (.not. states_conv) then
      message(1) = "Linear solver failed to converge all states."
      call write_warning(1)
    endif

    if (.not.(this%add_fxc .or. this%add_hartree)) then
    ! no need to deal with mixing, SCF iterations, etc.
    ! dealing with restart density above not necessary, but easier to just leave it
    ! convergence criterion is now about individual states, rather than SCF residual
      this%ok = states_conv

      message(1)="--------------------------------------------"
      write(message(2), '(a, i8)') &
           'Info: Total Hamiltonian applications:', total_iter * linear_solver_ops_per_iter(this%solver)
      call write_info(2)
      exit
    endif

    ! all the rest is the mixing and checking for convergence

    if(this%scf_tol%max_iter == iter) then 
      message(1) = "Self-consistent iteration for response did not converge."
      this%ok = .false.
      call write_warning(1)
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
    call write_info(2)
      
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
      call write_info(2)
      exit
    else
      ! not quitting if converged allows results to be calculated if possible
      ! before dying on the next direction or frequency
      if(clean_stop()) then
        message(1) = "Exiting cleanly."
        call write_fatal(1, only_root_writes = .true.)
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
subroutine X(sternheimer_calc_hvar)(this, sys, hm, lr, nsigma, hvar)
  type(sternheimer_t),    intent(inout) :: this
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr(:) 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(out)   :: hvar(:,:,:)

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
    do ispin2 = 1, sys%st%d%nspin
      hvar(1:np, ispin, 1) = hvar(1:np, ispin, 1) + this%fxc(1:np, ispin, ispin)*lr(1)%X(dl_rho)(1:np, ispin)
    end do
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


!--------------------------------------------------------------
subroutine X(sternheimer_solve_order2)( &
     this, sys, hm, lr, nsigma, omega, perturbation,       &
     restart_dir, rho_tag, wfs_tag, have_restart_rho)
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

  PUSH_SUB(X(sternheimer_solve_order2))

  ! construct inhomogeneous RHS term from first-order solution
  ! call X(sternheimer_solve) for sum frequency
  ! call X(sternheimer_solve) for difference frequency
  ! unless of course they are identical (one or both freqs the same)

  POP_SUB(X(sternheimer_solve_order2))
end subroutine X(sternheimer_solve_order2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
