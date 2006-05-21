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

!! This subroutine implements the preconditioned Lanczos eigensolver as
!! described in the paper:
!!
!! Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu and S. Ogut,
!! "Solution of Large Eigenvalue Problems in Electronic Structure Calculations",
!! BIT 36, 1 (1996).
!!
!! We also implement the "smoothing" preconditioning described in that paper.


subroutine X(eigen_solver_plan) (gr, st, hamilt, tol, niter, converged, diff)
  type(grid_t),        target, intent(inout) :: gr
  type(states_t),      target, intent(inout) :: st
  type(hamiltonian_t), target, intent(inout) :: hamilt
  FLOAT,                       intent(in)    :: tol
  integer,                     intent(inout) :: niter
  integer,                     intent(out)   :: converged
  FLOAT,             optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)


  ! Local stuff
  !  integer :: n          ! Dimension of the problem.
  integer :: ned        ! Number of smallest eigenpairs desired
  integer :: nec        ! number of eigen-pairs converged, if initially
                        ! nec > 0, the first nec elements of eigenval, res and
                        ! first nec columns of eigenvec are assumed to have converged
                        ! eigen-pairs and corresponding residual norms.
  integer :: maxmatvecs ! Maximum number of matrix-vectors applications allowed.
                        ! On exit reset to actual number of MATVECs used
  integer :: me         ! array size of eigenval, res and number of columns in eigenvec.
  FLOAT,  allocatable :: eigenval(:)     ! The eigenvalues
  R_TYPE, allocatable :: eigenvec(:,:,:) ! The eigenvectors
  FLOAT,  allocatable :: res(:)          ! The residuals
  R_TYPE, allocatable :: v(:,:,:)        ! The Krylov subspace basis vectors
  R_TYPE, allocatable :: av(:,:,:)       ! Workspace: W = A V
  FLOAT,  allocatable :: tmp(:)          ! Workspace.
  R_TYPE, allocatable :: h(:,:)          ! Projection of the hamiltonian onto Krylov subspace.
  R_TYPE, allocatable :: hevec(:,:)
  R_TYPE, allocatable :: aux(:,:)

  integer  :: blk, i, ii, idim, dim, j, d1, d2, matvec, nconv, ik, knec
  FLOAT :: x

  ! Some hard coded parameters.
  integer, parameter  :: winsiz = 5  ! window size, number of eigenvalues computed simultaneously
  integer, parameter  :: krylov = 15 ! The Krylov subspace size.
  FLOAT,   parameter  :: eps    = CNST(1e-15)


  call push_sub('eigen_plan.eigen_solver_plan')

  !  n          = m%np*st%d%dim
  dim        = st%d%dim
  ned        = st%nst
  nec        = 0
  maxmatvecs = niter*st%d%nik*st%nst
  me         = ned + winsiz - 1

  ! Allocate memory
  ! Careful: aux has to range from 1 to NP_PART because it is input to
  ! hpsi. In parallel the space NP+1:NP_PART is needed for ghost points
  ! in the non local operator.
  ALLOCATE(eigenvec(NP_PART, dim, me), NP_PART*dim*me)
  ALLOCATE(aux(NP_PART, dim),          NP_PART*dim)
  ALLOCATE(tmp(krylov),                krylov)
  ALLOCATE(v(NP_PART, dim, krylov),    NP_PART*dim*krylov)
  ALLOCATE(h(krylov, krylov),          krylov*krylov)
  ALLOCATE(eigenval(me),               me)
  ALLOCATE(av(NP_PART, dim, krylov),   NP_PART*dim*krylov)
  ALLOCATE(hevec(krylov, krylov),      krylov*krylov)
  ALLOCATE(res(me),                    me)

  eigenval = M_ZERO
  eigenvec = R_TOTYPE(M_ZERO)
  res      = M_ZERO
  v        = R_TOTYPE(M_ZERO)
  av       = R_TOTYPE(M_ZERO)
  tmp      = M_ZERO
  h        = R_TOTYPE(M_ZERO)
  hevec    = R_TOTYPE(M_ZERO)
  aux      = R_TOTYPE(M_ZERO)

  knec  = 0 ! Initialize the total (including all irreps.) converged eigenvector counter.
  niter = 0 ! Initialize the total matrix-vector multiplication counter.

  ! Main loop: runs over the irreducible subspaces.
  k_points : do ik = 1, st%d%nik

    ! First of all, copy the initial estimates.
    do i = 1, st%nst
      call lalg_copy(NP_PART, dim, st%X(psi)(:, :, i, ik), eigenvec(:, :, i))
      eigenval(i) = st%eigenval(i, ik)
    end do

    ! Initialization of counters...
    matvec = 0 ! Set the matrix-multiplication counter to zero.
    nec    = 0 ! Sets the converged vectors counter to zero.
    d1     = 0 ! index for inner loop
    nconv  = 0 ! number of eigen-pairs converged.

    ! Sets the projected hamiltonian matrix to zero.
    h = R_TOTYPE(M_ZERO)

    ! Beginning of the outer loop; start/restart
    outer_loop : do

      if(nec >= ned)           exit outer_loop ! :)   Already converged!
      if(matvec >= maxmatvecs) exit outer_loop ! :(   Maximum number of mat-vec operation surpassed...

      if (d1.le.winsiz) then !start from beginning
        blk = winsiz
      else                    !restart to work on another set of eigen-pairs
        blk = min(krylov/2, d1)
      end if

      !copy next set of Ritz vector/initial guesses to V
      do i = 1, winsiz
        call lalg_copy(NP_PART, dim, eigenvec(:, :, nec+i), v(:, :, i))
      end do

      ! Beginning of the inner loop.
      d1 = 0
      inner_loop: do
        d2 = d1 + blk

        ! Orthonormalization. The vectors in the v are orthonormalized against the converged
        ! eigenvectors, and among themselves. A check is done for the case of linear dependence,
        ! and random vectors are created in that case.
        i = d1 + 1
        ortho: do
          if(i>d2) exit ortho
          do ii = 1, nec
            av(ii, 1, d1 + 1) = X(states_dotp)(gr%m, dim, eigenvec(:,:,ii), v(:,:,i))
            call lalg_axpy(NP_PART, dim, -av(ii, 1, d1 + 1), eigenvec(:, :, ii), v(:, :, i))
          end do
          do ii = 1, i - 1
            av(ii, 1, d1 + 1) = X(states_dotp)(gr%m, dim, v(:, :, ii), v(:, :, i))
            call lalg_axpy(NP_PART, dim, -av(ii, 1, d1 + 1), v(:, :, ii), v(:, :, i))
          end do
          x = X(states_nrm2)(gr%m, dim, v(:, :, i))
          if(x .le. eps) then
            call X(mf_random)(gr%m, v(:, 1, i))
          else
            call lalg_scal(NP_PART, dim, R_TOTYPE(M_ONE/x), v(:, :, i))
            i = i + 1
          end if
        end do ortho

        ! matrix-vector multiplication
        do i = 1, blk
          call lalg_copy(NP_PART, dim, v(:, :, d1 + i), aux(:, :))
          av(:, :, d1 + i) = R_TOTYPE(M_ZERO)
          call X(Hpsi)(hamilt, gr, aux, av(:, :, d1 + i), ik)
        end do
        matvec = matvec + blk

        ! Here we calculate the last blk columns of H = V^T A V. We do not need the lower
        ! part of  the matrix since it is symmetric (LAPACK routine only need the upper triangle)
        do i = d1 + 1, d2
          do ii = 1, i
            h(ii, i) = X(states_dotp)(gr%m, dim, v(:, :, ii), av(:, :, i))
          end do
        end do

        ! Diagonalization in the subspace, by using LAPACK.
        call lalg_eigensolve(d2, h(1:d2, 1:d2), hevec(1:d2, 1:d2), tmp(1:d2))

        ! Store the Ritz values as approximate eigenvalues.
        call lalg_copy(winsiz, tmp, eigenval(nec+1:nec+winsiz))

        if ( d2+1.le.krylov .and. matvec.lt.maxmatvecs) then
          ! In this case, compute only the lowest Ritz eigenpair.
          call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, 1), &
            R_TOTYPE(M_ZERO), eigenvec(:, :, nec + 1))
          call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, 1), &
            R_TOTYPE(M_ZERO), av(:, :, d2 + 1))
          call residual(dim, av(:, :, d2+1), eigenvec(:, :, nec+1), tmp(1), av(:, :, d2+1), res(nec+1))

          ! If the first Ritz eigen-pair converged, compute all
          ! Ritz vectors and the residual norms.
          if(res(nec+1)<tol) then
            do i = 2, winsiz
              call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, i), &
                R_TOTYPE(M_ZERO), eigenvec(:, :, nec+i))
            end do
            do i = 2, winsiz
              call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, i), &
                R_TOTYPE(M_ZERO), v  (:, :, i))
            end do
            do i = 2, winsiz
              call residual(dim, v(:, :, i), eigenvec(:, :, nec+i), tmp(i), av(:, :, i), res(nec+i))
            end do
          end if
          d1 = d2
        else
          do i = 1, winsiz
            call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, i), &
              R_TOTYPE(M_ZERO), eigenvec(:, :, nec+i))
          end do
          do i = 1, winsiz
            call lalg_gemv(NP_PART, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, i), &
              R_TOTYPE(M_ZERO), v(:, :, i))
          end do
          do i = 1, winsiz
            call lalg_copy(NP_PART, dim, v(:, :, i), av(:, :, i))
            call lalg_copy(NP_PART, dim, eigenvec(:, :, nec + i), v(:, :, i))
            call residual(dim, av(:, :, i), v(:, :, i), tmp(i), av(:, :, winsiz+i), res(nec+i))
          end do

          ! Forms the first winsiz rows of H = V^T A V
          do i = 1, winsiz
            do ii = 1, i
              h(ii, i) = X(states_dotp)(gr%m, dim, v(:, :, ii), av(:, :, i))
            end do
          end do
          d1 = winsiz
        end if
        blk = 1

        ! Convergence test, and reordering of the eigenpairs. Starts checking
        ! the convergece of the eigenpairs of the window, and stops checking
        ! whenever finds one not converged. Then, for each converged eigenpair,
        ! compares its eigenvalue to the previous one, swapping them if
        ! necessary.
        nconv = 0
        ordering: do i = nec + 1, nec + winsiz - 1
          if(res(i) >= tol) exit ordering
          nconv = nconv + 1
          do j = i, 2, -1
            if (eigenval(j-1) <= eigenval(j)) exit
            x = eigenval(j-1); eigenval(j-1) = eigenval(j); eigenval(j) = x
            x = res(j-1); res(j-1) = res(j); res(j) = x
            call lalg_swap(NP_PART, dim, eigenvec(:, :, j), eigenvec(:, :, j-1))
          end do
        end do ordering

        ! If the maximum mat-vecs is surpassed, get out of here.
        if(matvec > maxmatvecs) exit outer_loop

        ! Restart if: any eigenpair is converged.
        if (nconv > 0) then
          nec = nec + nconv
          if (d2+1 > krylov) d1 = d2
          cycle outer_loop
        end if

        ! Preconditioning
        do idim = 1, dim
          call lalg_copy(NP_PART, av(:, idim, d1 + 1), aux(:, idim))
          call X(nl_operator_operate) (filter, aux(:, idim), v(:, idim, d1+1))
        end do

      end do inner_loop
    end do outer_loop

    do i = 1, st%nst
      call lalg_copy(NP_PART, dim, eigenvec(:, :, i), st%X(psi)(:, :, i, ik))
      st%eigenval(i, ik) = eigenval(i)
      diff(i, ik) = res(i)
    end do

    knec = knec + nec
    niter = niter + matvec

  end do k_points

  converged = knec
  deallocate(eigenval, eigenvec, res,  v, av, tmp, h, hevec, aux)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine residual(dim, hv, v, e, res, r)
    integer, intent(in)    :: dim
    R_TYPE,  intent(inout) :: hv(:,:)
    R_TYPE,  intent(inout) :: v(:,:)
    FLOAT,   intent(in)    :: e
    R_TYPE,  intent(inout) :: res(:,:)
    FLOAT,   intent(out)   :: r

    call push_sub('eigen_plan.residual')

    res = hv - e*v
    r = X(states_nrm2)(gr%m, dim, res)

    call pop_sub()
  end subroutine residual

end subroutine X(eigen_solver_plan)
