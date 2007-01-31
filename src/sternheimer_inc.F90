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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: em_resp_inc.F90 2663 2007-01-25 09:04:29Z lorenzen $


! -------------------------------------------------------------
! this is the central routine of the electromagnetic response
! it calculates the first order variations of the wavefunctions 
! for an electric field
!--------------------------------------------------------------

subroutine X(sternheimer_solve)(this, sys, h, lr, dir, tag, nsigma, omega, restart_dir)
  type(sternheimer_t) :: this
  type(system_t), target, intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: h
  type(lr_t),             intent(inout) :: lr(:,:) !ndim,nsigma
  integer,                intent(in)    :: dir
  integer,                intent(in)    :: tag
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(in)    :: omega
  character(len=*),       intent(in)    :: restart_dir

  FLOAT :: dpsimod
  integer :: iter, sigma, ik, ik2, ist, i, err
  R_TYPE, allocatable :: dl_rhoin(:, :, :), dl_rhonew(:, :, :), dl_rhotmp(:, :, :)
  R_TYPE, allocatable :: Y(:, :, :),dV(:, :, :), tmp(:)
  FLOAT  :: abs_dens
  R_TYPE :: omega_sigma

  logical :: conv_last, conv

  type(mesh_t), pointer :: m
  type(states_t), pointer :: st

  integer total_iter

  character(len=30) :: dirname

  call push_sub('static_pol_lr.get_response_e')

  ASSERT( nsigma==1 .or. nsigma ==2 )
  
  m => sys%gr%m
  st => sys%st
  
  call mix_init(this%mixer, sys%gr%m, sys%st%d%nspin, 1, func_type=sys%st%d%wfs_type)
  
  ALLOCATE(tmp(m%np),m%np)
  ALLOCATE(Y(m%np, 1, nsigma), m%np*1*nsigma)
  ALLOCATE(dV(m%np, st%d%nspin, nsigma), m%np*st%d%nspin*nsigma)
  ALLOCATE(dl_rhoin(m%np, st%d%nspin, 1), 1*m%np*st%d%nspin)
  ALLOCATE(dl_rhonew(m%np, st%d%nspin, 1), 1*m%np*st%d%nspin)
  ALLOCATE(dl_rhotmp(m%np, st%d%nspin, 1),1*m%np*st%d%nspin)

  conv = .false.
  conv_last = .false.

  call init_response_e()

  this%solver%tol = scf_tol_start(this%scftol, this%solver%initial_tol, this%solver%final_tol)

  message(1)="--------------------------------------------"
  call write_info(1)

  !self consistency iteration for response
  iter_loop: do iter=1, this%scftol%max_iter
    write(message(1), '(a, i3)') "LR SCF Iteration: ", iter
    write(message(2), '(a, f20.6, a, f20.6, a, i1)') &
         "Frequency: ", R_REAL(omega), " Eta : ", R_AIMAG(omega), " Dir: ", dir
    call write_info(2)

    total_iter = 0
    
    dl_rhoin(1:m%np, 1:st%d%nspin, 1) = lr(dir, 1)%X(dl_rho)(1:m%np, 1:st%d%nspin)

    !calculate the variation of hartree term
    if (this%add_hartree) then
      do sigma=1,nsigma
        do i = 1, m%np
          tmp(i) = sum(lr(dir, sigma)%X(dl_rho)(i, 1:st%d%nspin))
        end do
        call X(poisson_solve)(sys%gr, lr(dir,sigma)%X(dl_Vhar), tmp, all_nodes=.false.)
      end do
    end if
    

    do ik = 1, st%d%nspin
      do sigma = 1, nsigma

        !Calculate H^(1):

        !* Vext
        dV(1:m%np, ik, sigma) = m%x(1:m%np, dir)

        !* hartree
        if (this%add_hartree) dV(1:m%np, ik, sigma) = dV(1:m%np, ik, sigma) &
             + lr(dir, sigma)%X(dl_Vhar)(1:m%np)

        !* fxc
        if(this%add_fxc) then 
          do ik2 = 1, st%d%nspin
            dV(1:m%np, ik, sigma) = dV(1:m%np, ik, sigma) + &
                 lr(dir, sigma)%dl_Vxc(1:m%np, ik, ik2)*lr(dir, sigma)%X(dl_rho)(1:m%np, ik2)
          end do
        end if
      end do

      !now calculate response for each state
      do ist = 1, st%nst
        do sigma = 1, nsigma
          if(st%occ(ist, ik) <= lr_min_occ) cycle
          
          !calculate the RHS of the Sternheimer eq
          Y(1:m%np, 1, sigma) = -dV(1:m%np, ik, sigma)*st%X(psi)(1:m%np, 1, ist, ik)

          !and project it into the unoccupied states
          if(this%orth_response) then 
            call X(lr_orth_vector)(m, st, Y(:,:, sigma), ik)
          end if
        
          if(sigma==1) then 
            omega_sigma = omega
          else 
            omega_sigma = -R_CONJ(omega)
          end if

          !solve the Sternheimer equation
          call X(solve_HXeY) (this%solver, h, sys%gr, sys%st, ik, lr(dir, sigma)%X(dl_psi)(:,:, ist, ik),&
               Y(:,:, sigma), -sys%st%eigenval(ist, ik) + omega_sigma)
          
          !altough the dl_psi we get should be orthogonal to psi
          !a re-orthogonalization is sometimes necessary 
          if(this%orth_response) then 
            call X(lr_orth_vector)(m, st, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), ik)
          end if

          ! print the norm of the variations, and the number of
          ! iterations and residual of the linear solver
          tmp(1:m%np) = R_ABS(lr(dir, sigma)%X(dl_psi)(1:m%np, 1, ist, ik))**2
          dpsimod = X(mf_integrate)(m, tmp)
          write(message(1), '(i4, f20.6, i5, e20.6)') &
               (3-2*sigma)*ist, dpsimod, iter, this%solver%abs_psi 
          call write_info(1)

          total_iter=total_iter + iter
          
        end do !sigma
      end do !ist
    end do !ik

    call X(lr_build_dl_rho)(m, st, lr(dir,:), omega, nsigma)

    dl_rhonew(1:m%np, 1:st%d%nspin, 1) = M_ZERO

    !write restart info
!    call X(restart_write_lr_density)(sys, lr(dir, 1), R_REAL(omega), dir)
    do sigma=1,nsigma 
      write(dirname,'(a,i1,a,i1,a,i1)') restart_dir//"wfs", dir, "_", tag, "_", sigma
      call restart_write(trim(tmpdir)//dirname, st, sys%gr, err, iter=iter, lr=lr(dir, sigma))
    end do
    
    !all the rest is the mixing and checking for convergency

    if( this%scftol%max_iter == iter  ) then 
      message(1) = "Self-consistent iteration for response did not converge"
      this%ok = .false.
      call write_warning(1);
    end if


    dl_rhotmp(1:m%np, 1:st%d%nspin, 1) = lr(dir,1)%X(dl_rho)(1:m%np, 1:st%d%nspin)

    call X(mixing)(this%mixer, m, iter, st%d%nspin, 1, dl_rhoin, dl_rhotmp, dl_rhonew)

    abs_dens =  M_ZERO

    do ik = 1, st%d%nspin
      tmp(:) = R_REAL(dl_rhoin(:, ik, 1) - dl_rhotmp(:, ik, 1))**2 &
           + M_zI*R_AIMAG(dl_rhoin(:,ik, 1) - dl_rhotmp(:, ik, 1))**2
      abs_dens = abs_dens + abs(X(mf_integrate)(m, tmp))
    end do

    abs_dens = sqrt(abs_dens)

    write(message(1), '(a, e20.6, a, i5)') "SCF Residual ", abs_dens, " Total Iterations ", total_iter

    message(2)="--------------------------------------------"
    call write_info(2)
      
    if( abs_dens <= this%scftol%conv_abs_dens ) then 
      if(conv_last) then 
        conv = .true.
      else
        conv_last = .true.
      end if
    end if

    if(conv) then
      this%ok = .true.

      write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, ' iterations'
      call write_info(1)
      exit

    else
      
      lr(dir,1)%X(dl_rho)(1:m%np, 1:st%d%nspin) = dl_rhonew(1:m%np, 1:st%d%nspin, 1)
      if(nsigma == 2) lr(dir,2)%X(dl_rho)(1:m%np, 1:st%d%nspin) = R_CONJ(dl_rhonew(1:m%np, 1:st%d%nspin, 1))

      this%solver%tol = scf_tol_step(this%scftol, iter, abs_dens)

    end if
    
  end do iter_loop

  call scf_tol_end(this%scftol)

  call mix_end(this%mixer)

  deallocate(tmp)
  deallocate(Y)
  deallocate(dV)
  deallocate(dl_rhoin)
  deallocate(dl_rhonew)
  deallocate(dl_rhotmp)

  call pop_sub()

contains

  !------------------------------------------------------------
  subroutine init_response_e()
    integer :: ierr

    call push_sub('static_pol_lr.init_response_e')

    do sigma=1,nsigma
      lr(dir,sigma)%X(dl_Vhar)(:) = M_ZERO
      if(this%orth_response) then 
        call X(lr_orth_response)(m, st, lr(dir,sigma))
      end if
    end do

    if(.not. this%from_scratch ) then 
      !try to read the density from restart information
!      call X(restart_read_lr_density)(sys, lr(dir, 1), R_REAL(omega), dir, ierr)
      if(nsigma == 2) lr(dir, 2)%X(dl_rho) = R_CONJ(lr(dir, 1)%X(dl_rho))
    else 
      ierr = 1 
    end if

    !if this fails, build density from wavefunctions
    if ( ierr .ne. 0 ) then 
      call X(lr_build_dl_rho)(m, st, lr(dir,:), omega, nsigma)
    end if

    call pop_sub()
  end subroutine init_response_e

end subroutine X(sternheimer_solve)


