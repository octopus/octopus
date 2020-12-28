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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
!> conjugate-gradients method.
subroutine X(eigensolver_cg2) (namespace, gr, st, hm, xc, pre, tol, niter, converged, ik, diff, orthogonalize_to_all, &
  conjugate_direction, additional_terms, energy_change_threshold, shift)
  type(namespace_t),        intent(in)    :: namespace
  type(grid_t),             intent(in)    :: gr
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(preconditioner_t),   intent(in)    :: pre
  type(xc_t),               intent(in)    :: xc
  FLOAT,                    intent(in)    :: tol
  integer,                  intent(inout) :: niter
  integer,                  intent(inout) :: converged
  integer,                  intent(in)    :: ik
  FLOAT,                    intent(out)   :: diff(:) !< (1:st%nst)
  logical,                  intent(in)    :: orthogonalize_to_all
  integer,                  intent(in)    :: conjugate_direction
  logical,                  intent(in)    :: additional_terms
  FLOAT,                    intent(in)    :: energy_change_threshold
  FLOAT, pointer, optional, intent(in)   :: shift(:,:)

  R_TYPE, allocatable :: h_psi(:,:), g(:,:), g0(:,:),  cg(:,:), h_cg(:,:), psi(:, :), psi2(:, :), g_prev(:,:), psi_j(:,:)
  R_TYPE   :: gg, gg0, gg1, gamma, norma, dot
  FLOAT    :: es(2), cg0, e0, res, alpha, beta, theta, old_energy, lam, lam_conj, cg_phi, g0_norm, g_norm
  FLOAT    :: stheta, stheta2, ctheta, ctheta2
  FLOAT, allocatable :: chi(:, :), omega(:, :), fxc(:, :, :), lam_sym(:)
  FLOAT    :: integral_hartree, integral_xc, tmp
  integer  :: ist, jst, iter, maxter, idim, ip, isp, ixc, ib
  R_TYPE   :: sb(2)
  FLOAT    :: a0, b0, dsb(3)
  logical  :: fold_ ! use folded spectrum operator (H-shift)^2
  logical  :: add_xc_term
  type(states_elec_group_t) :: hpsi_j

  PUSH_SUB(X(eigensolver_cg2))

  ! if the optional shift argument is present, assume we are computing a folded spectrum
  fold_ =  present(shift)

  ! make sure the passed optional pointer is allocated
  if(fold_) then
    ASSERT(associated(shift))
  end if

  ! do we add the XC term? needs derivatives of the XC functional
  add_xc_term = additional_terms
  if(st%d%ispin == UNPOLARIZED) then
    isp = 1
  else
    isp = 2
  end if
  do ixc = 1, 2
    if(bitand(xc%kernel(ixc, isp)%flags, XC_FLAGS_HAVE_FXC) == 0) then
      add_xc_term = .false.
    end if
  end do
  if(bitand(xc%kernel_family, XC_FAMILY_LDA) == 0) then
    add_xc_term = .false.
  end if
  ! TODO: extend to spinors
  if(st%d%ispin == SPINORS) then
    add_xc_term = .false.
  end if

  maxter = niter
  niter = 0

  SAFE_ALLOCATE(  psi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   cg(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(    g(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(   g0(1:gr%mesh%np,      1:st%d%dim))
  SAFE_ALLOCATE(g_prev(1:gr%mesh%np,     1:st%d%dim))
  if(additional_terms) then
    SAFE_ALLOCATE(  chi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(omega(1:gr%mesh%np_part, 1:st%d%dim))
    if(st%d%ispin == UNPOLARIZED) then
      SAFE_ALLOCATE(fxc(1:gr%mesh%np, 1:1, 1:1))
    else if(st%d%ispin == SPIN_POLARIZED) then
      SAFE_ALLOCATE(fxc(1:gr%mesh%np, 1:2, 1:2))
    end if
  end if
  if(fold_) then
    SAFE_ALLOCATE( psi2(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_psi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE( h_cg(1:gr%mesh%np_part, 1:st%d%dim))
  else
    SAFE_ALLOCATE(h_psi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE( h_cg(1:gr%mesh%np, 1:st%d%dim))
  end if

  if(hm%theory_level == RDMFT) then
    SAFE_ALLOCATE(psi_j(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(lam_sym(1:st%nst))
    call states_elec_group_copy(st%d, st%group, hpsi_j, copy_data=.false.)
    do ib = hpsi_j%block_start, hpsi_j%block_end
      call X(hamiltonian_elec_apply_batch) (hm, namespace, gr%mesh, st%group%psib(ib, ik), hpsi_j%psib(ib, ik))
    end do
  end if

  h_psi = R_TOTYPE(M_ZERO)
  cg    = R_TOTYPE(M_ZERO)
  g     = R_TOTYPE(M_ZERO)
  g0    = R_TOTYPE(M_ZERO)
  h_cg  = R_TOTYPE(M_ZERO)
  g_prev = R_TOTYPE(M_ZERO)

  ! get derivative once here -> the density does not change in the loop
  if(add_xc_term) then
    fxc = M_ZERO
    call xc_get_fxc(xc, gr%mesh, namespace, st%rho, st%d%ispin, fxc)
  end if

  ! Set the diff to zero, since it is intent(out)
  diff(1:st%nst) = M_ZERO

  ! Start of main loop, which runs over all the eigenvectors searched
  ASSERT(converged >= 0)

  ! The steps in this loop follow closely the algorithm from
  ! Payne et al. (1992), Rev. Mod. Phys. 64, 4, section V.B
  eigenfunction_loop : do ist = converged + 1, st%nst
    gg1   = R_TOTYPE(M_ZERO)

    call states_elec_get_state(st, gr%mesh, ist, ik, psi)

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, psi, h_psi, ist, ik)

    if(fold_) then
      call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, h_psi, psi2, ist, ik)
      ! h_psi = (H-shift)^2 psi
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, gr%mesh%np
          h_psi(ip, idim) = psi2(ip, idim) - M_TWO*shift(ist,ik)*h_psi(ip, idim) + shift(ist,ik)**2*psi(ip, idim)
        end do
      end do
    end if

    ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
    st%eigenval(ist, ik) = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_psi))
    old_energy = st%eigenval(ist, ik)

    ! Starts iteration for this band
    iter_loop: do iter = 1, maxter
      ! need to save g from previous iteration for Polak-Ribiere method
      if(conjugate_direction == OPTION__CGDIRECTION__POLAK) then
        if(iter /= 1) then
          g_prev = g
        else
          g_prev = M_ZERO
        end if
      end if

      ! PTA92, eq. 5.10
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, gr%mesh%np
          g(ip, idim) = -h_psi(ip, idim) + st%eigenval(ist, ik)*psi(ip, idim)
        end do
      end do

      if (hm%theory_level == RDMFT) then
        ! For RDMFT, the gradient of the total energy functional differs from
        ! the DFT and HF cases, as the lagrange multiplier matrix lambda cannot
        ! be diagonalized together with the Hamiltonian. This is because the
        ! orbitals of the minimization are not the eigenstates of the
        ! single-body Hamiltonian, but of the systems 1RDM.
        ! The functional to be minimized here is:
        !  F = E[psi_i]-sum_ij lam_ij (<psi_i|psi_j> - delta_ij) + const.
        ! The respective gradient reads:
        !   dF/dpsi_i= dE/dphi_i - sum_j lam_ij |psi_j>= H|psi_i> - sum_j lam_ij |psi_j>
        ! We get the expression for lam_ij from the gradient with respect to
        ! psi*: lam_ij = <psi_i|dE/dpsi_j^*> = <psi_i|H|psi_j>
        ! NB: lam_ij != lam_ji until SCF convergence!
        do jst = 1, st%nst
          if (jst == ist) then
            lam_sym(jst) = M_TWO*st%eigenval(jst, ik)
          else
            call states_elec_get_state(st, gr%mesh, jst, ik, psi_j)
            do idim = 1, st%d%dim
              call batch_get_state(hpsi_j%psib(hpsi_j%iblock(jst, ik), ik), (/jst, idim/), gr%mesh%np, h_cg(:, idim))
            end do

            ! calculate <phi_j|H|phi_i> = lam_ji and <phi_i|H|phi_j> = lam_ij
            lam = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi_j, h_psi))
            lam_conj = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_cg))
            lam_sym(jst) = lam + lam_conj

            do idim = 1, st%d%dim
              call lalg_axpy(gr%mesh%np, -lam_conj, psi_j(:, idim), g(:, idim))
            end do
          end if
        end do
      end if

      ! PTA92, eq. 5.12
      ! Orthogonalize to all states
      call X(states_elec_orthogonalize_single_batch)(st, gr%mesh, ist - 1, ik, g, normalize = .false., &
          against_all=orthogonalize_to_all)
      g_norm = X(mf_nrm2) (gr%mesh, st%d%dim, g)
      if(g_norm /= M_ZERO) then
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, M_ONE/g_norm, g(1:gr%mesh%np, idim))
        end do
      end if

      ! PTA92, eq. 5.17
      ! Approximate inverse preconditioner
      call  X(preconditioner_apply)(pre, namespace, gr, hm, g(:,:), g0(:,:), ik)

      ! PTA92, eq. 5.18
      dot = X(mf_dotp) (gr%mesh, st%d%dim, psi, g0)
      ! This needs to be done before the orthogonalization_single call, as psi is not guaranted
      ! to be orthogonal to the other bands here
      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -dot, psi(:, idim), g0(:, idim))
      end do

      ! orthogonalize against previous or all states, depending on the optional argument orthogonalize_to_all
      call X(states_elec_orthogonalize_single_batch)(st, gr%mesh, ist - 1, ik, g0, normalize = .false., &
          against_all=orthogonalize_to_all)

      ! normalize g0 here, needed for consistency
      g0_norm = X(mf_nrm2) (gr%mesh, st%d%dim, g0)
      if(g0_norm /= M_ZERO) then
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, M_ONE/g0_norm, g0(1:gr%mesh%np, idim))
        end do
      end if

      ! dot products needed for conjugate gradient
      gg = X(mf_dotp) (gr%mesh, st%d%dim, g0, g, reduce = .false.)
      if(iter /= 1 .and. conjugate_direction == OPTION__CGDIRECTION__POLAK) then
        ! only needed for Polak-Ribiere
        gg1 = X(mf_dotp) (gr%mesh, st%d%dim, g0, g_prev, reduce = .false.)
      end if

      if(gr%mesh%parallel_in_domains) then
        sb(1) = gg1
        sb(2) = gg
        call comm_allreduce(gr%mesh%vp%comm, sb, dim = 2)
        gg1 = sb(1)
        gg  = sb(2)
      end if

      ! Starting or following iterations...
      if(iter  ==  1) then
        gg0 = gg
        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, g0(:, idim), cg(:, idim))
        end do
      else
        select case (conjugate_direction)
        case (OPTION__CGDIRECTION__FLETCHER)
          ! PTA eq. 5.20
          gamma = gg/gg0        ! (Fletcher-Reeves)
        case (OPTION__CGDIRECTION__POLAK)
          gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
        case default
          call messages_input_error(namespace, 'Conjugate Direction')
        end select
        ! save for next iteration
        gg0 = gg

        ! PTA92, eq. 5.19
        do idim = 1, st%d%dim
          !$omp parallel do simd schedule(static)
          do ip = 1, gr%mesh%np
            cg(ip, idim) = gamma*cg(ip, idim) + g0(ip, idim)
          end do
        end do
        call profiling_count_operations(st%d%dim*gr%mesh%np*(2*R_ADD + 2*R_MUL))
      end if

      ! PTA92, eq. 5.21
      norma =  X(mf_dotp) (gr%mesh, st%d%dim, psi, cg)
      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -norma, psi(1:gr%mesh%np, idim), cg(:, idim))
      end do

      ! normalize cg here (PTA92, eq. 5.22)
      cg0 = X(mf_nrm2) (gr%mesh, st%d%dim, cg)
      if(cg0 /= M_ZERO) then
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, M_ONE/cg0, cg(1:gr%mesh%np, idim))
        end do
      end if

      ! cg contains now the conjugate gradient
      call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, cg, h_cg, ist, ik)

      if(fold_) then
        call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, h_cg, psi2, ist, ik)
        ! h_psi = (H-shift)^2 psi
        do idim = 1, st%d%dim
          h_cg(1:gr%mesh%np, idim) = psi2(1:gr%mesh%np, idim) - M_TWO*shift(ist,ik)*h_cg(1:gr%mesh%np, idim) &
                                 + shift(ist,ik)**2*cg(1:gr%mesh%np, idim)
        end do
      end if

      ! Line minimization (eq. 5.23 to 5.38)
      a0 = M_TWO*R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_cg, reduce = .false.)) !Eq. 5.26
      b0  = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, cg, h_cg, reduce = .false.))

      if(gr%mesh%parallel_in_domains) then
        dsb(1) = a0
        dsb(2) = b0
        call comm_allreduce(gr%mesh%vp%comm, dsb, dim = 2)
        a0 = dsb(1)
        b0 = dsb(2)
      end if

      ! compare eq. 5.31
      e0 = st%eigenval(ist, ik)
      alpha = M_TWO * (e0 - b0)

      if (additional_terms) then
        ! more terms here, see PTA92 eqs 5.31, 5.32, 5.33, 5.36
        ! Hartree term
        tmp = M_TWO
        do idim = 1, st%d%dim
          !$omp parallel do simd schedule(static)
          do ip = 1, gr%mesh%np
            chi(ip, idim) = tmp * R_REAL(R_CONJ(cg(ip, idim)) * psi(ip, idim))
          end do
        end do
        call dpoisson_solve(hm%psolver, omega(:, 1), chi(:, 1), all_nodes = .false.)
        integral_hartree = dmf_dotp(gr%mesh, st%d%dim, chi, omega)

        ! exchange term
        ! TODO: adapt to different spin cases
        if(add_xc_term) then
          integral_xc = dmf_dotp(gr%mesh, st%d%dim, fxc(:, :, 1), chi(:, :)**2)
        else
          integral_xc = M_ZERO
        end if

        ! add additional terms to alpha (alpha is -d2e/dtheta2 from eq. 5.31)
        alpha = alpha - st%d%kweights(ik)*st%occ(ist, 1)/st%smear%el_per_state * &
          (integral_hartree + integral_xc) / gr%sb%rcell_volume**2
      end if

      beta = a0 * M_TWO

      ! For RDMFT, we get a different formula for the line minimization, which turns out to
      ! only change the beta of the original expression.
      ! beta -> beta + beta_rdmft, with beta_rdmft= - sum_j (lam_ji <cg_i|phi_k> + c.c.)
      if(hm%theory_level == RDMFT) then
        do jst = 1, st%nst
          call states_elec_get_state(st, gr%mesh, jst, ik, psi_j)
          cg_phi = M_TWO*R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi_j, cg))
          beta = beta - cg_phi * lam_sym(jst)
        end do
      end if

      ! Eq. 5.37
      theta = atan(beta/alpha)*M_HALF
      stheta = sin(theta)
      ctheta = cos(theta)
      es(1) = alpha * (M_HALF - stheta**2) + beta*M_TWO*stheta*ctheta
      stheta2 = sin(theta + M_PI*M_HALF)
      ctheta2 = cos(theta + M_PI*M_HALF)
      es(2) = alpha * (M_HALF - stheta2**2) + beta*M_TWO*stheta2*ctheta2

      ! Choose the minimum solutions.
      if (es(2) < es(1)) then
        theta = theta + M_PI*M_HALF
        a0 = ctheta2
        b0 = stheta2
      else
        a0 = ctheta
        b0 = stheta
      end if

      ! PTA92, eq. 5.38
      do idim = 1, st%d%dim
        !$omp parallel do simd schedule(static)
        do ip = 1, gr%mesh%np
          psi(ip, idim) = a0*psi(ip, idim) + b0*cg(ip, idim)
          h_psi(ip, idim) = a0*h_psi(ip, idim) + b0*h_cg(ip, idim)
        end do
      end do

      call profiling_count_operations(st%d%dim*gr%mesh%np*(2*R_ADD + 4*R_MUL))

      st%eigenval(ist, ik) = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_psi))
      res = X(states_elec_residue)(gr%mesh, st%d%dim, h_psi, st%eigenval(ist, ik), psi)

      if (hm%theory_level == RDMFT) then
        do idim = 1, st%d%dim
          call batch_set_state(hpsi_j%psib(hpsi_j%iblock(ist, ik), ik), (/ist, idim/), gr%mesh%np, h_psi(:, idim))
        end do
      end if

      if(debug%info) then
        write(message(1), '(a,i4,a,i4,a,i4,a,es12.5,a,i4,a,es12.5,a,es12.5,a,es12.5,es12.5)') 'Debug: CG Eigensolver - ik', ik, ' ist ', ist, &
             ' iter ', iter, " max ", maxter, &
             ' deltae ', abs(st%eigenval(ist, ik) - old_energy), &
             ' residue ', res
        call messages_info(1)
      end if

      old_energy = st%eigenval(ist, ik)

      ! Test convergence of residual
      if(res < tol .or. abs(gg0) < M_EPSILON) then
        if(converged == ist - 1) converged = ist ! only consider the first converged eigenvectors
        exit iter_loop
      end if
    end do iter_loop

    ! if the folded operator was used, compute the actual eigenvalue
    if(fold_) then
      call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, psi, h_psi, ist, ik)
      st%eigenval(ist, ik) = R_REAL(X(mf_dotp) (gr%mesh, st%d%dim, psi, h_psi, reduce = .true.))
      res = X(states_elec_residue)(gr%mesh, st%d%dim, h_psi, st%eigenval(ist, ik), psi)
    end if

    call states_elec_set_state(st, gr%mesh, ist, ik, psi)

    niter = niter + iter + 1

    diff(ist) = res

    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) +  ist, st%lnst*st%d%kpt%nlocal)
    end if

  end do eigenfunction_loop

  ! Deallocation of variables
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(h_psi)
  SAFE_DEALLOCATE_A(g)
  SAFE_DEALLOCATE_A(g0)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(h_cg)
  SAFE_DEALLOCATE_A(chi)
  SAFE_DEALLOCATE_A(omega)
  SAFE_DEALLOCATE_A(fxc)
  SAFE_DEALLOCATE_A(psi2)
  if(hm%theory_level == RDMFT) then
    SAFE_DEALLOCATE_A(psi_j)
    SAFE_DEALLOCATE_A(lam_sym)
    call states_elec_group_end(hpsi_j, st%d)
  end if

  POP_SUB(X(eigensolver_cg2))
end subroutine X(eigensolver_cg2)


! ---------------------------------------------------------
!> The algorithm is essentially taken from Jiang et al. Phys. Rev. B 68, 165337 (2003).
subroutine X(eigensolver_cg_jiang) (namespace, gr, st, hm, tol, niter, converged, ik, diff)
  type(namespace_t),        intent(in)    :: namespace
  type(grid_t),             intent(in)    :: gr
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  FLOAT,                    intent(in)    :: tol
  integer,                  intent(inout) :: niter
  integer,                  intent(inout) :: converged
  integer,                  intent(in)    :: ik
  FLOAT,          optional, intent(out)   :: diff(:) !< (1:st%nst)

  integer :: nst, dim, ist, maxter, i, conv, ip, idim
  R_TYPE, allocatable :: psi(:,:), phi(:, :), hcgp(:, :), cg(:, :), sd(:, :), cgp(:, :)
  FLOAT :: ctheta, stheta, ctheta2, stheta2, mu, lambda, dump, &
    gamma, sol(2), alpha, beta, theta, theta2, res, norm
  R_TYPE :: dot

  PUSH_SUB(X(eigensolver_cg_jiang))

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

  phi(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  psi(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)
  cgp(1:gr%mesh%np, 1:dim) = R_TOTYPE(M_ZERO)

  ! Set the diff to zero, since it is intent(out)
  if(present(diff)) diff(1:st%nst) = M_ZERO

  conv = converged
  states: do ist = conv + 1, nst

    call states_elec_get_state(st, gr%mesh, ist, ik, psi)

    ! Orthogonalize starting eigenfunctions to those already calculated...
    if(ist > 1) call X(states_elec_orthogonalize_single)(st, gr%mesh, ist - 1, ik, psi, normalize = .true.)

    ! Calculate starting gradient: |hpsi> = H|psi>
    call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, psi, phi, ist, ik)
    niter = niter + 1

    ! Initial settings for scalar variables.
    ctheta = M_ONE
    stheta = M_ZERO
    mu     = M_ONE

    ! Initialize to zero the vector variables.
    hcgp = R_TOTYPE(M_ZERO)
    cg   = R_TOTYPE(M_ZERO)


    band: do i = 1, maxter - 1 ! One operation has already been made.

      if( i >1 ) then ! Get H|psi> (through the linear formula)
        do idim = 1, st%d%dim
          do ip = 1, gr%mesh%np
            phi(ip, idim) = ctheta*phi(ip, idim) + stheta*hcgp(ip, idim)
          end do
        end do
      end if

      ! lambda = <psi|H|psi> = <psi|phi>
      lambda = R_REAL(X(mf_dotp)(gr%mesh, dim, psi, phi))

      ! Check convergence
      res = X(states_elec_residue)(gr%mesh, dim, phi, lambda, psi)

      if(debug%info) then
        norm = X(mf_nrm2)(gr%mesh, dim, phi)
        write(message(1), '(a,i4,a,i4,a,i4,a,es13.6,a,es13.6)') 'Debug: CG New Eigensolver - ik', ik, &
          ' ist ', ist, ' iter ', i + 1, ' res ', res, ' ', res/norm
        call messages_info(1)
      end if

      if(present(diff)) diff(ist) = res
      if(res < tol) then
        if(conv == ist - 1) conv = ist
        exit band
      end if

      ! Get steepest descent vector
      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          sd(ip, idim) = lambda*psi(ip, idim) - phi(ip, idim)
        end do
      end do

      if(ist > 1) call X(states_elec_orthogonalize_single)(st, gr%mesh, ist - 1, ik, sd, normalize = .false.)

      ! Get conjugate-gradient vector
      dump = X(mf_nrm2)(gr%mesh, dim, sd)**2
      gamma = dump/mu
      mu    = dump

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          cg(ip, idim) = sd(ip, idim) + gamma*cg(ip, idim)
        end do
      end do

      dot = X(mf_dotp)(gr%mesh, dim, psi, cg)

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          cgp(ip, idim) = cg(ip, idim) - dot*psi(ip, idim)
        end do
      end do

      norm = X(mf_nrm2)(gr%mesh, dim, cgp)

      call X(hamiltonian_elec_apply_single)(hm, namespace, gr%mesh, cgp, hcgp, ist, ik)

      niter = niter + 1

      alpha = -lambda + R_REAL(X(mf_dotp)(gr%mesh, dim, cgp, hcgp))/norm**2
      beta  = M_TWO*R_REAL(X(mf_dotp)(gr%mesh, dim, cgp, phi))/norm
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
      stheta = stheta/norm

      do idim = 1, st%d%dim
        do ip = 1, gr%mesh%np
          psi(ip, idim) = ctheta*psi(ip, idim) + stheta*cgp(ip, idim)
        end do
      end do

    end do band

    call states_elec_set_state(st, gr%mesh, ist, ik, psi)

    st%eigenval(ist, ik) = lambda

    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) +  ist, st%lnst*st%d%kpt%nlocal)
    end if

  end do states

  converged = conv

  SAFE_DEALLOCATE_A(phi)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(cg)
  SAFE_DEALLOCATE_A(hcgp)
  SAFE_DEALLOCATE_A(sd)
  SAFE_DEALLOCATE_A(cgp)

  POP_SUB(X(eigensolver_cg_jiang))
end subroutine X(eigensolver_cg_jiang)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
