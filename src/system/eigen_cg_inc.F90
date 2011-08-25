!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: eigen_cg_inc.F90 6376 2010-03-23 14:18:30Z joseba $

! ---------------------------------------------------------
!> conjugate-gradients method.
subroutine X(eigensolver_cg2) (gr, st, hm, pre, tol, niter, converged, ik, diff, verbose)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  FLOAT,        optional, intent(out)   :: diff(1:st%nst)
  logical,      optional, intent(in)    :: verbose

  R_TYPE, allocatable :: h_psi(:,:), g(:,:), g0(:,:),  cg(:,:), ppsi(:,:), psi(:, :)
  R_TYPE   :: es(2), a0, b0, gg, gg0, gg1, gamma, theta, norma
  real(8)  :: cg0, e0, res
  integer  :: p, iter, maxter, idim, ip
  logical  :: verbose_
  R_TYPE   :: sb(3)

  PUSH_SUB(X(eigensolver_cg2))

  verbose_ = .false.
  if(present(verbose)) verbose_ = verbose

  if(verbose_) then
    call messages_print_stress(stdout, "CG Info")
    message(1) = "Diagonalization with the conjugate gradients algorithm."
    write(message(2),'(a,e8.2)') '  Tolerance: ',tol
    write(message(3),'(a,i6)')   '  Maximum number of iterations per eigenstate:', niter
    message(4) = ''
    call messages_info(4)
  end if

  maxter = niter
  niter = 0

  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(h_psi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   cg(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(    g(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   g0(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE( ppsi(1:gr%mesh%np, 1:st%d%dim))
  h_psi = R_TOTYPE(M_ZERO)
  cg    = R_TOTYPE(M_ZERO)
  g     = R_TOTYPE(M_ZERO)
  g0    = R_TOTYPE(M_ZERO)
  ppsi  = R_TOTYPE(M_ZERO)

  do idim = 1, st%d%dim
    cg(1:gr%mesh%np_part, idim) = R_TOTYPE(M_ZERO)
  end do

  ! Set the diff to zero, since it is intent(out)
  if(present(diff)) diff(1:st%nst) = M_ZERO

  ! Start of main loop, which runs over all the eigenvectors searched
  ASSERT(converged >= 0)

  eigenfunction_loop : do p = converged + 1, st%nst

    if(verbose_) then
      write(message(2),'(a,i4,a)') ' Eigenstate # ', p, ':'
    end if

    call states_get_state(st, gr%mesh, p, ik, psi)

    ! Orthogonalize starting eigenfunctions to those already calculated...
    if(p > 1) then
      call X(states_orthogonalization)(gr%mesh, p - 1, st%d%dim, st%X(psi)(:, :, :, ik), psi, normalize = .true.)
    end if

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_apply)(hm, gr%der, psi, h_psi, p, ik)

    ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
    st%eigenval(p, ik) = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_psi))

    ! Starts iteration for this band
    iter_loop: do iter = 1, maxter

      ! inverse preconditioner....
      call  X(preconditioner_apply)(pre, gr, hm, ik, h_psi, g)
      call  X(preconditioner_apply)(pre, gr, hm, ik, psi, ppsi)

      es(1) = X(mf_dotp) (gr%mesh, st%d%dim, psi, g, reduce = .false.)
      es(2) = X(mf_dotp) (gr%mesh, st%d%dim, psi, ppsi, reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%vp%comm, es, dim = 2)

      es(1) = es(1)/es(2)

      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, R_TOPREC(-es(1)), ppsi(:, idim), g(:, idim))
      end do

      ! Orthogonalize to lowest eigenvalues (already calculated)
      if(p > 1) call X(states_orthogonalization)(gr%mesh, p - 1, st%d%dim, st%X(psi)(:, :, :, ik), g, normalize = .false.)

      if(iter .ne. 1) then
        gg1 = X(mf_dotp) (gr%mesh, st%d%dim, g, g0, reduce = .false.)
      else
        gg1 = M_ZERO
      end if

      ! Approximate inverse preconditioner...
      call  X(preconditioner_apply)(pre, gr, hm, ik, g(:,:), g0(:,:))

      gg = X(mf_dotp) (gr%mesh, st%d%dim, g, g0, reduce = .false.)

      if(gr%mesh%parallel_in_domains) then
        sb(1) = gg1
        sb(2) = gg
        call comm_allreduce(gr%mesh%vp%comm, sb, dim = 2)
        gg1 = sb(1)
        gg  = sb(2)
      end if

      if( abs(gg) < M_EPSILON ) then
        if(converged == p - 1) converged = p ! only consider the first converged eigenvectors
        st%eigenval(p, ik) = es(1)
        res = sqrt(abs(gg))
        exit
      end if

      ! Starting or following iterations...
      if(iter .eq. 1) then
        gg0 = gg

        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, g(:,idim), cg(:, idim))
        end do
      else
        !gamma = gg/gg0        ! (Fletcher-Reeves)
        gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
        gg0 = gg
        
        norma = gamma*cg0*sin(theta)
        
        forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
          cg(ip, idim) = gamma*cg(ip, idim) + g(ip, idim) - norma*psi(ip, idim)
        end forall
        
        call profiling_count_operations(st%d%dim*gr%mesh%np*(2*R_ADD + 2*R_MUL))

      end if

      ! cg contains now the conjugate gradient
      call X(hamiltonian_apply)(hm, gr%der, cg, ppsi, p, ik)

      ! Line minimization.
      a0 = X(mf_dotp) (gr%mesh, st%d%dim, psi, ppsi, reduce = .false.)
      b0 = X(mf_dotp) (gr%mesh, st%d%dim, cg, ppsi, reduce = .false.)
      cg0 = X(mf_nrm2) (gr%mesh, st%d%dim, cg, reduce = .false.)

      if(gr%mesh%parallel_in_domains) then
        sb(1) = a0
        sb(2) = b0
        sb(3) = cg0**2
        call comm_allreduce(gr%mesh%vp%comm, sb, dim = 3)
        a0 = sb(1)
        b0 = sb(2)
        cg0 = sqrt(sb(3))
      end if

      a0 = M_TWO * a0 / cg0
      b0 = b0/cg0**2
      e0 = st%eigenval(p, ik)
      theta = atan(R_REAL(a0/(e0 - b0)))/M_TWO
      es(1) = M_HALF*((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) + e0 + b0)
      es(2) = -M_HALF*((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) - (e0 + b0))

      ! Choose the minimum solutions.
      if (R_REAL(es(2)) < R_REAL(es(1))) theta = theta + M_PI/M_TWO
      st%eigenval(p, ik) = min(R_REAL(es(1)), R_REAL(es(2)))

      ! Upgrade psi...
      a0 = cos(theta)
      b0 = sin(theta)/cg0

      forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
        psi(ip, idim) = a0*psi(ip, idim) + b0*cg(ip, idim)
        h_psi(ip, idim) = a0*h_psi(ip, idim) + b0*ppsi(ip, idim)
      end forall

      call profiling_count_operations(st%d%dim*gr%mesh%np*(2*R_ADD + 4*R_MUL))

      res = X(states_residue)(gr%mesh, st%d%dim, h_psi, st%eigenval(p, ik), psi)

      if(in_debug_mode) then
        write(message(1), '(a,i4,a,i4,a,i4,a,f12.6)') 'Debug: CG Eigensolver - ik', ik, ' ist ', p, ' iter ', iter, ' res ', res
        call messages_info(1)
      end if

      ! Test convergence.
      if(res < tol) then
        if(converged == p - 1) converged = p ! only consider the first converged eigenvectors
        exit iter_loop
      end if

    end do iter_loop

    call states_set_state(st, gr%mesh, p, ik, psi)

    if(verbose_) then
      if(res<tol) then
        write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2)),"     converged. Iterations:", iter, '   [Res = ',res,']'
      else
        write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2))," not converged. Iterations:", maxter, '   [Res = ',res,']'
        ! if it didn't converge, then iter = maxter + 1, which is one more than the number of iterations actually done
      end if
      call messages_info(1)
    end if

    niter = niter + iter + 1

    if(present(diff)) then
      diff(p) = res
    end if

    if(mpi_grp_is_root(mpi_world).and..not.verbose_) then
      call loct_progress_bar(st%nst*(ik - 1) +  p, st%nst*st%d%nik)
    end if

  end do eigenfunction_loop

  ! Deallocation of variables
  SAFE_DEALLOCATE_A(h_psi)
  SAFE_DEALLOCATE_A(g)
  SAFE_DEALLOCATE_A(g0)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(ppsi)

  if(verbose_) call messages_print_stress(stdout)

  POP_SUB(X(eigensolver_cg2))
end subroutine X(eigensolver_cg2)


! ---------------------------------------------------------
!> The algorithm is essentially taken from Jiang et al. Phys. Rev. B 68, 165337 (2003).
subroutine X(eigensolver_cg2_new) (gr, st, hm, tol, niter, converged, ik, diff, verbose)
  type(grid_t),        intent(in)    :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  FLOAT,               intent(in)    :: tol
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged
  integer,             intent(in)    :: ik
  FLOAT,     optional, intent(out)   :: diff(1:st%nst)
  logical,   optional, intent(in)    :: verbose

  integer :: nst, dim, ist, maxter, i, conv, ip, idim
  logical :: verbose_
  R_TYPE, allocatable :: psi(:,:), phi(:, :), hcgp(:, :), cg(:, :), sd(:, :), cgp(:, :)
  FLOAT :: ctheta, stheta, ctheta2, stheta2, mu, lambda, dump, &
       gamma, sol(2), alpha, beta, theta, theta2, res ! Could be complex?
  R_TYPE :: dot
  logical, allocatable :: orthogonal(:)

  PUSH_SUB(X(eigensolver_cg2_new))

  verbose_ = .false.
  if(present(verbose)) verbose_ = verbose

  if(verbose_) then
    call messages_print_stress(stdout, "CG Info")
    message(1) = "Diagonalization with the conjugate gradients algorithm [new]."
    write(message(2),'(a,e8.2)') '  Tolerance: ', tol
    write(message(3),'(a,i6)')   '  Maximum number of iterations per eigenstate:', niter
    message(4) = ""
    call messages_info(4)
  end if

  dim = st%d%dim
  nst = st%nst

  maxter = niter
  niter = 0

  SAFE_ALLOCATE( phi(1:gr%mesh%np     , 1:dim))
  SAFE_ALLOCATE( psi(1:gr%mesh%np_part, 1:dim))
  SAFE_ALLOCATE(  cg(1:gr%mesh%np     , 1:dim))
  SAFE_ALLOCATE(hcgp(1:gr%mesh%np     , 1:dim))
  SAFE_ALLOCATE(  sd(1:gr%mesh%np     , 1:dim))
  SAFE_ALLOCATE( cgp(1:gr%mesh%np_part, 1:dim))
  SAFE_ALLOCATE(orthogonal(1:nst))

  phi(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  psi(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  cgp(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)

  ! Set the diff to zero, since it is intent(out)
  if(present(diff)) diff(1:st%nst) = M_ZERO

  conv = converged
  states: do ist = conv + 1, nst

    call states_get_state(st, gr%mesh, ist, ik, psi)

    ! Orthogonalize starting eigenfunctions to those already calculated...
    if(ist.gt.1) then
      call X(states_orthogonalization)(gr%mesh, ist - 1, st%d%dim, &
        st%X(psi)(:, :, 1:ist - 1, ik), psi, normalize = .true.)
    end if

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_apply)(hm, gr%der, psi, phi, ist, ik)
    niter = niter + 1

    ! Initial settings for scalar variables.
    ctheta = M_ONE
    stheta = M_ZERO
    mu     = M_ONE

    ! Initialize to zero the vector variables.
    hcgp = R_TOTYPE(M_ZERO)
    cg   = R_TOTYPE(M_ZERO)

    orthogonal = .false.

    band: do i = 1, maxter - 1 ! One operation has already been made.

      if(mod(i, 5).eq.0) orthogonal = .false.

      ! Get H|psi> (through the linear formula)
      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          phi(ip, idim) = ctheta*phi(ip, idim) + stheta*hcgp(ip, idim)
        end do
      end do

      ! lambda = <psi|H|psi> = <psi|phi>
      lambda = X(mf_dotp)(gr%mesh, dim, psi, phi)

      ! Check convergence
      res = X(states_residue)(gr%mesh, dim, phi, lambda, psi)
      if(present(diff)) diff(ist) = res
      if(res < tol) then
        conv = conv + 1
        exit band
      end if

      ! Get steepest descent vector
      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          sd(ip, idim) = lambda*psi(ip, idim) - phi(ip, idim)
        end do
      end do

      if(ist > 1) call X(states_orthogonalization)(gr%mesh, ist - 1, dim, st%X(psi)(:, :, :, ik), sd, &
        normalize = .false., mask = orthogonal)

      ! Get conjugate-gradient vector
      dot = X(mf_nrm2)(gr%mesh, dim, sd)**2
      gamma = dot/mu
      mu    = dot

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          cg(ip, idim) = sd(ip, idim) + gamma*cg(ip, idim)
        end do
      end do

      dump = X(mf_dotp)(gr%mesh, dim, psi, cg)

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          cgp(ip, idim) = cg(ip, idim) - dump*psi(ip, idim)
        end do
      end do

      dump = X(mf_nrm2)(gr%mesh, dim, cgp)

      do idim = 1, st%d%dim
        call lalg_scal(gr%mesh%np, M_ONE/dump, cgp(:, idim))
      end do

      call X(hamiltonian_apply)(hm, gr%der, cgp, hcgp, ist, ik)

      niter = niter + 1

      alpha = -lambda + X(mf_dotp)(gr%mesh, dim, cgp, hcgp)
      beta  = M_TWO*X(mf_dotp)(gr%mesh, dim, cgp, phi)
      theta = M_HALF*atan(-beta/alpha)
      ctheta = cos(theta)
      stheta = sin(theta)

      ! This checks whether we are picking the maximum or the minimum.
      theta2 = theta + M_PI/M_TWO
      ctheta2 = cos(theta2)
      stheta2 = sin(theta2)
      sol(1) = lambda + stheta**2*alpha + beta*stheta*ctheta
      sol(2) = lambda + stheta2**2*alpha + beta*stheta2*ctheta2

      if(sol(2) < sol(1)) then
        theta = theta2
        stheta = stheta2
        ctheta = ctheta2
      end if

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          psi(ip, idim) = ctheta*psi(ip, idim) + stheta*cgp(ip, idim)
        end do
      end do

    end do band

    call states_set_state(st, gr%mesh, ist, ik, psi)

    st%eigenval(ist, ik) = lambda

    if(verbose_) then
      if(res<tol) then
        write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2)),"     converged. Iterations:", i, '   [Res = ',res,']'
      else
        write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2))," not converged. Iterations:", i, '   [Res = ',res,']'
      end if
      call messages_info(1)
    end if

    if(mpi_grp_is_root(mpi_world).and..not.verbose_) then
      call loct_progress_bar(st%nst*(ik - 1) + ist, st%nst*st%d%nik)
    end if

  end do states

  converged = conv

  SAFE_DEALLOCATE_A(phi)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(hcgp)
  SAFE_DEALLOCATE_A(sd)
  SAFE_DEALLOCATE_A(cgp)
  SAFE_DEALLOCATE_A(orthogonal)
  if(verbose_) call messages_print_stress(stdout)

  POP_SUB(X(eigensolver_cg2_new))
end subroutine X(eigensolver_cg2_new)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
