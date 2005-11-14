!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! conjugate-gradients method.
subroutine eigen_solver_cg2(gr, st, h, tol, niter, converged, diff, reorder, verbose)
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(inout) :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,      optional, intent(in)    :: reorder
  logical,      optional, intent(in)    :: verbose

  R_TYPE, allocatable :: h_psi(:,:), g(:,:), g0(:,:), &
    cg(:,:), ppsi(:,:)

  R_TYPE :: es(2), a0, b0, gg, gg0, gg1, gamma, theta, norma
  FLOAT :: cg0, e0, res
  integer  :: ik, moved, p, j, iter, maxter, conv, conv_
  logical  :: reord = .true., verbose_

  call push_sub('eigen_cg.eigen_solver_cg2')

  verbose_ = .false.; if(present(verbose)) verbose_ = verbose
  if(verbose_) then
    call messages_print_stress(stdout, "CG Info")
    message(1) = "Diagonalization with the conjugate gradients algorithm."
    write(message(2),'(a,e8.2)') '  Tolerance: ',tol
    write(message(3),'(a,i6)')   '  Maximum number of iterations per eigenstate:', niter
    message(4) = ''
    call write_info(4)
  end if

  if(present(reorder)) reord = reorder
  conv_ = 0
  maxter = niter
  niter = 0
  moved = 0

  allocate(h_psi(1:NP_PART, st%d%dim), ppsi(1:NP_PART, st%d%dim), &
    g(1:NP_PART, st%d%dim), g0(1:NP_PART, st%d%dim), cg(1:NP_PART, st%d%dim))

  ! Start of main loop, which runs over all the eigenvectors searched
  ik_loop: do ik = 1, st%d%nik
    conv = converged
    eigenfunction_loop : do p = conv + 1, st%nst

      if(verbose_) then
        write(message(2),'(a,i4,a)') 'Eigenstate # ',p,':'
      end if

      ! Orthogonalize starting eigenfunctions to those already calculated...
      call X(states_gram_schmidt)(p, gr%m, st%d%dim, &
        st%X(psi)(:, 1:st%d%dim, 1:p, ik), start=p)

      ! Calculate starting gradient: |hpsi> = H|psi>
      call X(Hpsi)(h, gr, st%X(psi)(:,:, p, ik) , h_psi, ik)

      ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
      st%eigenval(p, ik) = R_REAL(X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, p, ik), h_psi))

      ! Starts iteration for this band
      iter_loop: do iter = 1, maxter

        ! if approximate inverse preconditioner....
        !call pre(hpsi%val   , g%val   )
        !call pre(psi(m)%val , ppsi%val)
        g(1:NP,:) = h_psi(1:NP,:)
        ppsi(1:NP,:) = st%X(psi)(1:NP,:, p, ik)

        es(1) = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, p, ik), g)
        es(2) = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, p, ik), ppsi)
        es(1) = es(1)/es(2)
        g(1:NP,:) = g(1:NP,:) - es(1)*ppsi(1:NP,:)

        ! Orthogonalize to lowest eigenvalues (already calculated)
        do j = 1, p - 1
          a0 = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, j, ik), g)
          g(1:NP,:) = g(1:NP,:) - a0 * st%X(psi)(1:NP,:, j, ik)
        end do

        if(iter .ne. 1) gg1 = X(states_dotp) (gr%m, st%d%dim, g, g0)

        ! Approximate inverse preconditioner...
        !call ipre(g, g0)
        g0(1:NP,:) = g(1:NP,:)

        gg = X(states_dotp) (gr%m, st%d%dim, g, g0)
        if( abs(gg) < CNST(1.0e-15) ) then
          conv = conv + 1
          st%eigenval(p, ik) = es(1)
          res = sqrt(abs(gg))
          exit
        end if

        ! Starting or following iterations...
        if(iter .eq. 1) then
          gg0 = gg
          cg (1:NP, :) = g(1:NP, :)
        else
          !gamma = gg/gg0        ! (Fletcher-Reeves)
          gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
          gg0 = gg
          cg(1:NP,:) = gamma*cg(1:NP,:)
          cg(1:NP,:) = cg(1:NP,:) + g(1:NP,:)

          norma = gamma*cg0*sin(theta)
          cg(1:NP, :) = cg(1:NP, :) - norma * st%X(psi)(1:NP, :, p, ik)
        end if

        ! cg contains now the conjugate gradient
        cg0 = X(states_nrm2) (gr%m, st%d%dim, cg(:,:))
        call X(Hpsi) (h, gr, cg, ppsi, ik)

        ! Line minimization.
        a0 = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, p, ik), ppsi)
        a0 = M_TWO * a0 / cg0
        b0 = X(states_dotp) (gr%m, st%d%dim, cg(:,:), ppsi)
        b0 = b0/cg0**2
        e0 = st%eigenval(p, ik)
        theta = atan(R_REAL(a0/(e0 - b0)))/M_TWO
        es(1) = ((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) + e0 + b0) / M_TWO
        es(2) =(-(e0-b0)*cos(M_TWO*theta) - a0*sin(M_TWO*theta) + e0 + b0) / M_TWO

        ! Choose the minimum solutions.
        if (R_REAL(es(2)) < R_REAL(es(1))) then
          theta = theta + M_PI/M_TWO
        end if
        st%eigenval(p, ik) = min(R_REAL(es(1)), R_REAL(es(2)))

        ! Upgrade psi...
        a0 = cos(theta)
        b0 = sin(theta)/cg0
        ! This does the sum: st%X(psi)(:, :, p, ik) = a0*st%X(psi)(:, :, p, ik) + b0*cg(:, :)
        ! It can crash in Intel compiler version 8 otherwise.
        do j = 1, st%d%dim
          call lalg_scal(NP, a0, st%X(psi)(:, j, p, ik))
          call lalg_axpy(NP, b0, cg(:, j), st%X(psi)(:, j, p, ik))
        end do

        ! Calculate H|psi>
        h_psi(1:NP,:) = a0*h_psi(1:NP,:) + b0*ppsi(1:NP,:)

        res = X(states_residue)(gr%m, st%d%dim, h_psi, st%eigenval(p, ik), &
          st%X(psi)(:, :, p, ik))
        ! Test convergence.
        if(res < tol) then
          conv = conv + 1
          exit iter_loop
        end if

      end do iter_loop

      if(verbose_) then
        if(res<tol) then
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2)),"     converged. Iterations:", iter, '   [Res = ',res,']'
        else
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2))," not converged. Iterations:", iter, '   [Res = ',res,']'
        end if
        call write_info(1)
      end if

      niter = niter + iter + 1

      if(present(diff)) then
        diff(p, ik) = res
      end if

      if(p>1 .and. reord) call sort(st%eigenval(1:p, ik), st%X(psi)(1:NP, :, 1:p, ik))

    end do eigenfunction_loop

    conv_ = conv_ + conv

  end do ik_loop

  converged = conv_
  ! Deallocation of variables
  deallocate(h_psi, g, g0, cg, ppsi)

  if(verbose_) call messages_print_stress(stdout)

  call pop_sub()
end subroutine eigen_solver_cg2


! The algorithm is essentially taken from Jiang et al. Phys. Rev. B 68, 165337 (2003).
subroutine eigen_solver_cg2_new(gr, st, h, tol, niter, converged, diff, reorder, verbose)
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(inout) :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,      optional, intent(in)    :: reorder
  logical,      optional, intent(in)    :: verbose

  integer :: nik, nst, dim, ik, ist, maxter, i, k, conv, conv_
  logical :: verbose_, reorder_
  R_TYPE, allocatable :: psi(:,:), phi(:, :), hpsi(:, :), hcgp(:, :), cg(:, :), sd(:, :), cgp(:, :)
  FLOAT :: ctheta, stheta, ctheta2, stheta2, mu, lambda, dump, &
           gamma, sol(2), alpha, beta, theta, theta2, res ! Could be complex?

  call push_sub('eigen_cg.eigen_solver_cg2_new')

  verbose_ = .false.; if(present(verbose)) verbose_ = verbose
  reorder_ = .true. ; if(present(reorder)) reorder_ = reorder
  if(verbose_) then
    call messages_print_stress(stdout, "CG Info")
    message(1) = "Diagonalization with the conjugate gradients algorithm [new]."
    write(message(2),'(a,e8.2)') '  Tolerance: ',tol
    write(message(3),'(a,i6)')   '  Maximum number of iterations per eigenstate:', niter
    message(4) = ""
    call write_info(4)
  end if

  dim = st%d%dim; nik = st%d%nik; nst = st%nst

  conv_ = 0
  maxter = niter
  niter = 0

  allocate(phi(NP_PART, dim), psi(NP_PART, dim), hpsi(NP_PART, dim), &
     cg(NP_PART, dim), hcgp(NP_PART, dim), sd(NP_PART, dim), cgp(NP_PART, dim))

  kpoints: do ik = 1, nik
    conv = converged
    states: do ist = conv + 1, nst

      ! Orthogonalize starting eigenfunctions to those already calculated...
      call X(states_gram_schmidt)(ist, gr%m, dim, st%X(psi)(:, 1:dim, 1:ist, ik), start=ist)
      psi(1:NP, :) = st%X(psi)(1:NP, :, ist, ik)

      ! Calculate starting gradient: |hpsi> = H|psi>
      call X(Hpsi)(h, gr, psi, phi, ik); niter = niter + 1

      ! Initial settings for scalar variables.
      ctheta = M_ONE
      stheta = M_ZERO
      mu     = M_ONE

      ! Initialize to zero the vector variables.
      hcgp = R_TOTYPE(M_ZERO)
      cg   = R_TOTYPE(M_ZERO)

      band: do i = 1, maxter - 1 ! One operation has already been made.

         ! Get H|psi> (through the linear formula)
         phi(1:NP, :) = ctheta*phi(1:NP, :) + stheta*hcgp(1:NP, :)

         ! lambda = <psi|H|psi> = <psi|phi>
         lambda = X(states_dotp)(gr%m, dim, psi, phi)

         ! Check convergence
         res = X(states_residue)(gr%m, dim, phi, lambda, psi)
         if(present(diff)) diff(ist, ik) = res
         if(res < tol) then
           conv = conv + 1
           exit band
         end if

         ! Get steepest descent vector
         sd(1:NP, :) = lambda*psi(1:NP, :) - phi(1:NP, :)
         do k = 1, ist - 1
            dump = X(states_dotp)(gr%m, dim, st%X(psi)(:, :, k, ik), sd(:, :))
            sd(1:NP, :) = sd(1:NP, :) - dump*st%X(psi)(1:NP, :, k, ik)
         end do

         ! Get conjugate-gradient vector
         gamma = X(states_dotp)(gr%m, dim, sd, sd)/mu
         mu    = X(states_dotp)(gr%m, dim, sd, sd)
         cg(1:NP,:) = sd(1:NP,:) + gamma*cg(1:NP,:)

         !
         dump = X(states_dotp)(gr%m, dim, psi, cg)
         cgp(1:NP,:) = cg(1:NP,:) - dump*psi(1:NP,:)
         dump = sqrt(X(states_dotp)(gr%m, dim, cgp, cgp))
         cgp(1:NP,:) = cgp(1:NP,:)/dump

         call X(Hpsi)(h, gr, cgp, hcgp, ik); niter = niter + 1

         alpha = - lambda + X(states_dotp)(gr%m, dim, cgp, hcgp)
         beta  = M_TWO*X(states_dotp)(gr%m, dim, cgp, phi)
         theta = M_HALF*atan(-beta/alpha)
         ctheta = cos(theta)
         stheta = sin(theta)

         ! I am not sure wether this is necessary or not.
         theta2 = theta + M_PI/M_TWO
         ctheta2 = cos(theta2)
         stheta2 = sin(theta2)
         sol(1) = ctheta**2*lambda + stheta**2*X(states_dotp)(gr%m, dim, cgp, hcgp) + &
                                     M_TWO*stheta*ctheta*X(states_dotp)(gr%m, dim, cgp, phi)
         sol(2) = ctheta**2*lambda + stheta2**2*X(states_dotp)(gr%m, dim, cgp, hcgp) + &
                                     M_TWO*stheta2*ctheta2*X(states_dotp)(gr%m, dim, cgp, phi)

         if(sol(2) < sol(1)) then
            theta = theta2
            stheta = stheta2
            ctheta = ctheta2
         end if

         psi(1:NP,:) = ctheta*psi(1:NP,:) + stheta*cgp(1:NP,:)

      end do band

      st%X(psi)(1:NP, :, ist, ik) = psi(1:NP, :)
      st%eigenval(ist, ik) = lambda

      if(verbose_) then
        if(res<tol) then
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2)),"     converged. Iterations:", i, '   [Res = ',res,']'
        else
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2))," not converged. Iterations:", i, '   [Res = ',res,']'
        end if
        call write_info(1)
      end if

      ! Reordering.
      if(ist>1 .and. reorder_) call sort(st%eigenval(1:ist, ik), st%X(psi)(:, :, 1:ist, ik))

    end do states

    conv_ = conv_ + conv

  end do kpoints

  converged = conv_

  deallocate(phi, psi, hpsi, cg, hcgp, sd, cgp)
  if(verbose_) call messages_print_stress(stdout)

  call pop_sub()
end subroutine eigen_solver_cg2_new
