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

!! This subroutine implements the preconditioned Lanczos eigensolver aas
!! described in the paper:
!!
!! Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu and S. Ogut,
!! "Solution of Large Eigenvalue Problems in Electronic Structure Calculations",
!! BIT 36, 1 (1996).
!!
!! We also implement the "smoothing" preconditioning described in that paper.

#include "config_F90.h"

subroutine eigen_solver_plan(st, sys, hamilt, tol, niter, converged, diff)
  implicit none
  type(states_type), target, intent(inout)   :: st
  type(system_type), target, intent(IN)      :: sys
  type(hamiltonian_type), target, intent(IN) :: hamilt
  real(r8), intent(in)               :: tol
  integer, intent(inout)             :: niter
  integer, intent(out)               :: converged
  real(r8), intent(out), optional    :: diff(1:st%nst,1:st%nik)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Local stuff
  integer :: n          ! Dimension of the problem.
  integer :: ned        ! Number of smallest eigenpairs desired
  integer :: nec        ! number of eigen-pairs converged, if initially
                        ! nec > 0, the first nec elements of eigenval, res and
                        ! first nec columns of eigenvec are assumed to have converged
                        ! eigen-pairs and corresponding residual norms.
  integer :: maxmatvecs ! Maximum number of matrix-vectors applications allowed.
                        ! On exit reset to actual number of MATVECs used
  integer :: me         ! array size of eigenval, res and number of columns in eigenvec.
  real(r8), allocatable :: eigenval(:)  ! The eigenvalues
  R_TYPE, allocatable :: eigenvec(:, :) ! The eigenvectors
  real(r8), allocatable :: res(:)       ! The residuals
  R_TYPE, allocatable :: v(:, :)        ! The Krylov subspace basis vectors
  R_TYPE, allocatable :: av(:, :)       ! Workspace: W = A V
  real(r8), allocatable :: tmp(:)       ! Workspace.
  R_TYPE, allocatable :: h(:, :)        ! Projection of the hamiltonian onto Krylov subspace.
  R_TYPE, allocatable :: hevec(:, :)
  R_TYPE, allocatable :: aux(:, :)
  
  integer  :: blk, i, ii, idim, j, jj, d1, d2, matvec, nconv, ik, knec, np
  real(r8) :: x

  ! Some hard coded parameters.
  integer, parameter  :: winsiz = 5 ! window size, number of eigenvalues computed simultaneously
  integer, parameter  :: krylov = 15 ! The Krylov subspace size.
  real(r8), parameter :: eps    = 1e-15_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sub_name = 'eigen_solver_blan'; call push_sub()
  
  n          = sys%m%np*st%dim
  np         = sys%m%np
  ned        = st%nst
  nec        = 0
  maxmatvecs = niter*st%nik*st%nst
  me         = ned + winsiz - 1
  
  allocate(eigenval(me),          &
           eigenvec(n, me),       &
           res(me),               &
           v(n, krylov),          &
           av(n, krylov),         &
           tmp(krylov),           &
           h(krylov, krylov),     &
           hevec(krylov, krylov), &
           aux(0:sys%m%np, st%dim))

  knec  = 0 ! Initialize the total (including all irreps.) converged eigenvector counter.
  niter = 0 ! Initialize the total matrix-vector multiplication counter.

  ! Main loop: runs over the irreducible subspaces.
  k_points : do ik = 1, st%nik
    
    ! First of all, copy the initial estimates.
    do i = 1, sys%st%nst
      do idim = 1, st%dim
        call R_FUNC(copy)(np, st%R_FUNC(psi)(1, idim, i, ik), 1, eigenvec((idim-1)*np+1, i), 1)
      enddo
      eigenval(i) = st%eigenval(i, ik)
    enddo
    
    ! Initialization of counters...
    matvec = 0 ! Set the matrix-multiplication counter to zero.
    nec    = 0 ! Sets the converged vectors counter to zero.
    d1     = 0 ! index for inner loop
    nconv  = 0 ! number of eigen-pairs converged.
    
    ! Sets the projected hamiltonian matrix to zero.
    h = R_TOTYPE(M_ZERO)
    
    ! Beginning of the outer loop; start/restart
    outer_loop : do
      
      if(nec>=ned)           exit outer_loop ! :)   Already converged!
      if(matvec>=maxmatvecs) exit outer_loop ! :(   Maximum number of mat-vec operation surpassed...
      
      if (d1.le.winsiz) then !start from beginning
        blk = winsiz
      else                    !restart to work on another set of eigen-pairs
        blk = min(krylov/2, d1)
      endif
      
      !copy next set of Ritz vector/initial guesses to V
      call R_FUNC(copy)(n*winsiz, eigenvec(1, nec + 1), 1, v(1, 1), 1)
      
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
            av(ii, d1 + 1) = R_FUNC(states_dotp)(sys%m, st%dim, eigenvec(1, ii), v(1, i))
            call R_FUNC(axpy)(n, -av(ii, d1 + 1), eigenvec(1, ii), 1, v(1, i), 1)
          enddo
          do ii = 1, i - 1
            av(ii, d1 + 1) = R_FUNC(states_dotp)(sys%m, st%dim, v(1, ii), v(1, i))!/sys%m%vol_pp
            call R_FUNC(axpy)(n, -av(ii, d1 + 1), v(1, ii), 1, v(1, i), 1)
          enddo
          x = R_FUNC(states_nrm2)(sys%m, st%dim, v(1, i))
          if(x .le. eps) then
            call R_FUNC(mesh_random)(sys%m, v(1, i))
          else
            call R_FUNC(scal)(n, R_TOTYPE(M_ONE/x), v(1, i), 1)
            i = i + 1
          endif
        enddo ortho
        
        !matrix-vector multiplication
        do i = 1, blk
          do idim = 1, st%dim
            aux(1:sys%m%np, idim) = v((idim-1)*sys%m%np+1:idim*sys%m%np, d1 + i)
          enddo
          aux(0, :) = R_TOTYPE(M_ZERO)
          call R_FUNC(Hpsi)(hamilt, sys%m, st, sys, ik, aux, av(1, d1 + i))
        enddo
        matvec = matvec + blk
        
        ! Here we calculate the last blk columns of H = V^T A V. We do not need the lower
        ! part of  the matrix since it is symmetric (LAPACK routine only need the upper triangle)
        do i = d1 + 1, d2
          do ii = 1, i
            h(ii, i) = R_FUNC(states_dotp)(sys%m, st%dim, v(1, ii), av(1, i))
          enddo
        enddo
        
        ! Diagonalization in the subspace, by using LAPACK.
        call R_FUNC(iagonalise)(d2, h(1:d2, 1:d2), hevec(1:d2, 1:d2), tmp(1:d2))
        
        ! Store the Ritz values as approximate eigenvalues.
        call dcopy(winsiz, tmp(1), 1, eigenval(nec+1), 1)
       
        if ( d2+1.le.krylov .and. matvec.lt.maxmatvecs) then
          ! In this case, compute only the lowes Ritz eigenpair.
          call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), v, n, hevec, 1, &
               R_TOTYPE(M_ZERO), eigenvec(1, nec + 1), 1)
          call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), av, n, hevec, 1, &
               R_TOTYPE(M_ZERO), av(1, d2 + 1), 1)
          call residual(av(1:n, d2+1), eigenvec(1:n, nec+1), tmp(1), av(1:n, d2+1), res(nec+1))
          
          ! If the first Ritz eigen-pair converged, compute all 
          ! Ritz vectors and the residual norms.
          if(res(nec+1)<tol) then
            do i = 2, winsiz
              call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), v, n, hevec(1, i), 1, &
                   R_TOTYPE(M_ZERO), eigenvec(1, nec+i), 1)
            enddo
            do i = 2, winsiz
              call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), av, n, hevec(1, i), 1, &
                   R_TOTYPE(M_ZERO), v  (1, i),     1) 
            enddo
            do i = 2, winsiz
              call residual(v(1, i), eigenvec(1, nec+i), tmp(i), av(1, i), res(nec+i))
            enddo
          endif
          d1 = d2
        else
          do i = 1, winsiz
            call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), v, n, hevec(1, i), 1, &
                 R_TOTYPE(M_ZERO), eigenvec(1, nec+i), 1)
          enddo
          do i = 1, winsiz
            call R_FUNC(gemv)('N', n, d2, R_TOTYPE(M_ONE), av, n, hevec(1, i), 1, &
                 R_TOTYPE(M_ZERO), v(1, i), 1)
          enddo
          do i = 1, winsiz
            av(:, i) = v(:, i)
            v(:, i) = eigenvec(:, nec + i)
            call residual(av(1, i), v(1, i), tmp(i), av(1, winsiz+i), res(nec+i))
          enddo
          
          ! Forms the first winsiz rows of H = V^T A V
          do i = 1, winsiz
            do ii = 1, i
              h(ii, i) = R_FUNC(states_dotp)(sys%m, st%dim, v(1, ii), av(1, i))
            enddo
          enddo
          
          d1 = winsiz
        endif
        blk = 1
        
        ! Convergence test, and reordering of the eigenpairs. Starts checking
        ! the convergece of the eigenpairs of the window, and stops checking
        ! whenever finds one not converged. Then, for each converged eigenpair,
        ! compares its eigenvalue to the previous one, swapping them if
        ! necessary.
        nconv = 0
        ordering: do i = nec + 1, nec+winsiz-1
          if(res(i) >= tol) exit ordering
          nconv = nconv + 1
          do j = i, 2, -1
            if (eigenval(j-1) <= eigenval(j)) exit
            call dswap(1, eigenval(j-1), 1, eigenval(j), 1)
            call dswap(1, res(j-1),  1, res(j),  1)
            call R_FUNC(swap)(n, eigenvec(1,j), 1, eigenvec(1,j-1), 1)
          enddo
        enddo ordering
        
        ! If the maximum mat-vecs is surpassed, get out of here.
        if(matvec > maxmatvecs) exit outer_loop
        
        ! Restart if: any eigenpair is converged.
        if (nconv > 0) then
          nec = nec + nconv
          if (d2+1 > krylov) d1 = d2
          cycle outer_loop
        endif
        
        ! Preconditioning
        do idim = 1, st%dim
          call R_FUNC(copy) (np, av((idim-1)*np+1, d1 + 1), 1, aux(1, idim), 1)
          aux(0, idim) = R_TOTYPE(M_ZERO)
          call R_FUNC(low_frequency)(sys%m, aux(:, idim), v((idim-1)*sys%m%np+1, d1 + 1))
        enddo
        
      enddo inner_loop
      
    enddo outer_loop
    
    do i = 1, sys%st%nst
      do idim = 1, st%dim
        call R_FUNC(copy)(sys%m%np, eigenvec((idim-1)*sys%m%np+1, i), 1, st%R_FUNC(psi)(1, idim, i, ik), 1)
        st%R_FUNC(psi)(0, 1, i, ik) = R_TOTYPE(M_ZERO)
      enddo
      st%eigenval(i, ik) = eigenval(i)
      diff(i, ik) = res(i)
    enddo
    
    knec = knec + nec
    niter = niter + matvec
    
  enddo k_points
  
  converged = knec
  deallocate(eigenval, eigenvec, res,  v, av, tmp, h, hevec, aux)
  call pop_sub(); return
  
contains
  
  subroutine residual(hv, v, e, res, r)
    implicit none
    R_TYPE, intent(inout), dimension(n) :: hv
    R_TYPE, intent(inout), dimension(n) :: v
    real(r8), intent(in)  :: e
    R_TYPE, intent(inout), dimension(n) :: res
    real(r8), intent(out) :: r
    
    res(1:n) = hv(1:n) - e*v(1:n)
    r = R_FUNC(states_nrm2)(sys%m, st%dim, res)
    
  end subroutine residual
  
  subroutine diagonalise(dim, a, b, e)
    integer  :: dim
    real(r8) :: a(dim, dim), b(dim, dim)
    real(r8) :: e(dim)
    
    integer :: info, lwork
    real(r8), allocatable :: work(:)

    lwork = 6*dim
    allocate(work(lwork))
    b = a
    call dsyev('V', 'U', dim, b, dim, e, work, lwork, info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'dsyev returned error message ', info
      call write_fatal(1)
    endif
    deallocate(work)
  end subroutine diagonalise

  subroutine ziagonalise(dim, a, b, e)
    integer     :: dim
    complex(r8) :: a(dim, dim), b(dim, dim)
    real(r8)    :: e(dim)
    
    integer :: info, lwork
    complex(r8), allocatable :: work(:)
    real(r8), allocatable :: rwork(:)

    lwork = 6*dim
    allocate(work(lwork), rwork(lwork))
    b = a
    call zheev('V','U', dim, b, dim, e, work, lwork, rwork, info)
    if(info.ne.0) then
      write(message(1),'(a,i5)') 'dsyev returned error message ', info
      call write_fatal(1)
    endif
    deallocate(work, rwork)
  end subroutine ziagonalise

end subroutine eigen_solver_plan
