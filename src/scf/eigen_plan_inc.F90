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
!! $Id$

!! This subroutine implements the preconditioned Lanczos eigensolver as
!! described in the paper:
!!
!! Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu and S. Ogut,
!! "Solution of Large Eigenvalue Problems in Electronic Structure Calculations",
!! BIT 36, 1 (1996).
!!
!! We also implement the "smoothing" preconditioning described in that paper.


subroutine X(eigensolver_plan) (gr, st, hm, pre, tol, niter, converged, ik, diff)
  type(grid_t),                intent(inout) :: gr
  type(states_t),              intent(inout) :: st
  type(hamiltonian_t),         intent(inout) :: hm
  type(preconditioner_t),      intent(in)    :: pre
  FLOAT,                       intent(in)    :: tol
  integer,                     intent(inout) :: niter
  integer,                     intent(out)   :: converged
  integer,                     intent(in)    :: ik
  FLOAT,             optional, intent(out)   :: diff(1:st%nst)


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

  integer  :: blk, i, ii, idim, dim, j, d1, d2, matvec, nconv
  FLOAT :: x

  ! Some hard coded parameters.
  integer, parameter  :: winsiz = 5  ! window size, number of eigenvalues computed simultaneously
  integer, parameter  :: krylov = 15 ! The Krylov subspace size.

  call push_sub('eigen_plan.eigensolver_plan')

  !  n          = m%np*st%d%dim
  dim        = st%d%dim
  ned        = st%nst
  nec        = 0
  maxmatvecs = niter*st%d%nik*st%nst
  me         = ned + winsiz - 1

  ! Allocate memory
  ! Careful: aux has to range from 1 to gr%mesh%np_part because it is input to
  ! hpsi. In parallel the space NP+1:gr%mesh%np_part is needed for ghost points
  ! in the non local operator.
  ALLOCATE(eigenvec(gr%mesh%np, dim, me),      gr%mesh%np*dim*me)
  ALLOCATE(aux(gr%mesh%np_part, dim),          gr%mesh%np_part*dim)
  ALLOCATE(tmp(krylov),                krylov)
  ALLOCATE(v(gr%mesh%np, dim, krylov),    gr%mesh%np*dim*krylov)
  ALLOCATE(h(krylov, krylov),          krylov*krylov)
  ALLOCATE(eigenval(me),               me)
  ALLOCATE(av(gr%mesh%np, dim, krylov),   gr%mesh%np*dim*krylov)
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

  niter = 0 ! Initialize the total matrix-vector multiplication counter.

  ! First of all, copy the initial estimates.
  do i = 1, st%nst
    do idim = 1, dim
      call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, i, ik), eigenvec(:, idim, i))
    end do
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
      do idim = 1, dim
        call lalg_copy(gr%mesh%np, eigenvec(:, idim, nec+i), v(:, idim, i))
      end do
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
          av(ii, 1, d1 + 1) = X(mf_dotp)(gr%mesh, dim, eigenvec(:,:,ii), v(:,:,i))
          do idim = 1, dim
            call lalg_axpy(gr%mesh%np, -av(ii, 1, d1 + 1), eigenvec(:, idim, ii), v(:, idim, i))
          end do
        end do
        do ii = 1, i - 1
          av(ii, 1, d1 + 1) = X(mf_dotp)(gr%mesh, dim, v(:, :, ii), v(:, :, i))
          do idim = 1, dim
            call lalg_axpy(gr%mesh%np, -av(ii, 1, d1 + 1), v(:, idim, ii), v(:, idim, i))
          end do
        end do
        x = X(mf_nrm2)(gr%mesh, dim, v(:, :, i))
        if(x .le. M_EPSILON) then
          call X(mf_random)(gr%mesh, v(:, 1, i))
        else
          do idim = 1, dim
            call lalg_scal(gr%mesh%np, R_TOTYPE(M_ONE/x), v(:, idim, i))
          end do
          i = i + 1
        end if
      end do ortho

      ! matrix-vector multiplication
      do i = 1, blk
        do idim = 1, dim
          call lalg_copy(gr%mesh%np, v(:, idim, d1 + i), aux(:, idim))
        end do
        av(:, :, d1 + i) = R_TOTYPE(M_ZERO)
        call X(hamiltonian_apply)(hm, gr, aux, av(:, :, d1 + i), d1+i, ik)
      end do
      matvec = matvec + blk

      ! Here we calculate the last blk columns of H = V^T A V. We do not need the lower
      ! part of  the matrix since it is symmetric (LAPACK routine only need the upper triangle)
      do i = d1 + 1, d2
        do ii = 1, i
          h(ii, i) = X(mf_dotp)(gr%mesh, dim, v(:, :, ii), av(:, :, i))
        end do
      end do

      ! Diagonalization in the subspace, by using LAPACK.
      call lalg_eigensolve(d2, h(1:d2, 1:d2), hevec(1:d2, 1:d2), tmp(1:d2))

      ! Store the Ritz values as approximate eigenvalues.
      call lalg_copy(winsiz, tmp, eigenval(nec+1:nec+winsiz))

      if ( d2+1.le.krylov .and. matvec.lt.maxmatvecs) then
        ! In this case, compute only the lowest Ritz eigenpair.
        call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, 1), &
             R_TOTYPE(M_ZERO), eigenvec(:, :, nec + 1))
        call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, 1), &
             R_TOTYPE(M_ZERO), av(:, :, d2 + 1))
        call residual(dim, av(:, :, d2+1), eigenvec(:, :, nec+1), tmp(1), av(:, :, d2+1), res(nec+1))

        ! If the first Ritz eigen-pair converged, compute all
        ! Ritz vectors and the residual norms.
        if(res(nec+1)<tol) then
          do i = 2, winsiz
            call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, i), &
                 R_TOTYPE(M_ZERO), eigenvec(:, :, nec+i))
          end do
          do i = 2, winsiz
            call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, i), &
                 R_TOTYPE(M_ZERO), v  (:, :, i))
          end do
          do i = 2, winsiz
            call residual(dim, v(:, :, i), eigenvec(:, :, nec+i), tmp(i), av(:, :, i), res(nec+i))
          end do
        end if
        d1 = d2
      else
        do i = 1, winsiz
          call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), v(:, :, 1:d2), hevec(1:d2, i), &
               R_TOTYPE(M_ZERO), eigenvec(:, :, nec+i))
        end do
        do i = 1, winsiz
          call lalg_gemv(gr%mesh%np, dim, d2, R_TOTYPE(M_ONE), av(:, :, 1:d2), hevec(1:d2, i), &
               R_TOTYPE(M_ZERO), v(:, :, i))
        end do
        do i = 1, winsiz
          do idim = 1, dim
            call lalg_copy(gr%mesh%np, v(:, idim, i), av(:, idim, i))
            call lalg_copy(gr%mesh%np, eigenvec(:, idim, nec + i), v(:, idim, i))
          end do
          call residual(dim, av(:, :, i), v(:, :, i), tmp(i), av(:, :, winsiz+i), res(nec+i))
        end do

        ! Forms the first winsiz rows of H = V^T A V
        do i = 1, winsiz
          do ii = 1, i
            h(ii, i) = X(mf_dotp)(gr%mesh, dim, v(:, :, ii), av(:, :, i))
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
          do idim = 1, dim
            call lalg_swap(gr%mesh%np, eigenvec(:, idim, j), eigenvec(:, idim, j-1))
          end do
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
        call lalg_copy(gr%mesh%np, av(:, idim, d1 + 1), aux(:, idim))
      end do
      call X(preconditioner_apply)(pre, gr, hm, aux(:,:), v(:,:, d1+1))

    end do inner_loop
  end do outer_loop

  do i = 1, st%nst
    do idim = 1, dim
      call lalg_copy(gr%mesh%np, eigenvec(:, idim, i), st%X(psi)(:, idim, i, ik))
    end do
    st%eigenval(i, ik) = eigenval(i)
    diff(i) = res(i)
  end do

  converged = nec
  niter = niter + matvec

  SAFE_DEALLOCATE_A(eigenval)
  SAFE_DEALLOCATE_A(eigenvec)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(av)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(h)
  SAFE_DEALLOCATE_A(hevec)
  SAFE_DEALLOCATE_A(aux)
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
    r = X(mf_nrm2)(gr%mesh, dim, res)

    call pop_sub()
  end subroutine residual

end subroutine X(eigensolver_plan)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
