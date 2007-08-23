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
!! $Id: eigen.F90 3030 2007-06-25 16:45:05Z marques $

  ! Implementation of the locally optimal block preconditioned conjugate
  ! gradients algorithm.

#include "global.h"

! Index set of unconverged eigenvectors.
#define UC uc(1:nuc)

  ! ---------------------------------------------------------
  ! Locally optimal block preconditioned conjugate gradient algorithm.
  ! For details, see:
  !
  ! A. Knyazev. Toward the Optimal Preconditioned Eigensolver: Locally
  ! Optimal Block Preconditioned Conjugate Gradient Method. SIAM
  ! Journal on Scientific Computing, 23(2):517­541, 2001.
  !
  ! A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchin-
  ! nikov. Block Locally Optimal Preconditioned Eigenvalue Xolvers
  ! (BLOPEX) in hypre and PETSc. SIAM Journal of Scientific Computing,
  ! 2007.
  subroutine X(eigen_solver_lobpcg)(gr, st, h, pre, tol, niter, converged, diff, verbose)
    type(grid_t),           intent(inout) :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(inout) :: h
    type(preconditioner_t), intent(in)    :: pre
    FLOAT,                  intent(in)    :: tol
    integer,                intent(inout) :: niter
    integer,                intent(inout) :: converged
    FLOAT, optional,        intent(out)   :: diff(1:st%nst, 1:st%d%nik)
    logical, optional,      intent(in)    :: verbose

    integer :: np            ! Number of points per state.
    integer :: nep           ! Number of eigenpairs.
    integer :: nst           ! Number of eigenstates (i. e. the blocksize).
    integer :: ik, ist
    integer :: i, j, k, maxiter
    integer :: conv
    integer :: nuc, uc(st%nst)
    logical :: verbose_

    FLOAT,  allocatable :: diffs(:, :)
    R_TYPE, allocatable :: tmp(:, :, :)     ! Temporary storage of wavefunction size.
    R_TYPE, allocatable :: nuc_tmp(:, :)    ! Temporary storage of Gram block size.
    R_TYPE, allocatable :: res(:, :, :)     ! Residuals.
    R_TYPE, allocatable :: h_res(:, :, :)   ! H res.
    R_TYPE, allocatable :: dir(:, :, :)     ! Conjugate directions.
    R_TYPE, allocatable :: h_dir(:, :, :)   ! H dir.
    R_TYPE, allocatable :: h_psi(:, :, :)   ! H |psi>.
    R_TYPE, pointer     :: psi(:, :, :)     ! Pointer to the eigenstates.
    FLOAT, pointer      :: eval(:, :)       ! Pointer to the eigenvalues.
    R_TYPE, allocatable :: gram_h(:, :)     ! Gram matrix for Hamiltonian.
    R_TYPE, allocatable :: gram_i(:, :)     ! Gram matrix for unit matrix.
    R_TYPE, allocatable :: gram_block(:, :) ! Space to construct the Gram matrix blocks.
    R_TYPE, allocatable :: ritz_vec(:, :)   ! Ritz-vectors.
    FLOAT,  allocatable :: ritz_val(:)      ! Ritz-values.
    R_TYPE, allocatable :: ritz_psi(:, :)   ! Block of ritz_vec to calculate new psi.
    R_TYPE, allocatable :: ritz_dir(:, :)   ! Block of ritz_vec to calculate new dir.
    R_TYPE, allocatable :: ritz_res(:, :)   ! Block of ritz_vec to calculate new res.

    call push_sub('eigen_lobpcg.Xeigen_solver_lobpcg')

    ! Check if LOBPCG can be run.
    ! LOBPCG does not work with domain parallelization for the moment
    ! (because we do things like <psi|psi>).
    if(gr%m%parallel_in_domains) then
      message(1) = 'The LOBPCG eigenvsolver cannot be used with domain parallelization'
      message(2) = 'at the moment. Choose a different eigensolver.'
      call write_fatal(3)
    end if

    ! Some abbreviations.
    eval => st%eigenval
    nep  =  st%nst * st%d%nik
    nst  =  st%nst
    np   =  NP_PART*st%d%dim

    if(present(diff)) then
      diff = R_TOTYPE(M_ZERO)
    end if

    verbose_ = .false.
    if(present(verbose)) verbose_ = verbose
    if(verbose_) then
      call messages_print_stress(stdout, "LOBPCG Info")
      message(1) = 'Diagonalization with the locally optimal block preconditioned'
      message(2) = 'conjugate gradient algorithm.'
      write(message(3),'(a,e8.2)') '  Tolerance: ', tol
      write(message(4),'(a,i6)')   '  Maximum number of iterations per block of eigenstates:', &
        niter
      message(5) = ''
      call write_info(5)
    end if

    ALLOCATE(diffs(nst, st%d%nik), nep)
    ALLOCATE(tmp(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(res(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(h_res(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(dir(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(h_dir(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(h_psi(NP_PART, st%d%dim, nst), NP_PART*st%d%dim*nst)
    ALLOCATE(gram_block(nst, nst), nst**2)
    ALLOCATE(ritz_val(nst), nst)

    maxiter = niter
    niter   = 0

    k_loop: do ik = 1, st%d%nik
      ! At the beginning, all eigenvectors are considered unconverged.
      nuc = nst
      do ist = 1, nst
        uc(ist) = ist
      end do

      ! Set pointer to current block, so we can drop the k index in the following.
      psi => st%X(psi)(:, :, :, ik)

      ! Orthonormalize initial vectors.
      call X(lobpcg_orth)(gr%m, st%d%dim, st%d%wfs_type, nst, np, psi, nuc, uc)

      ! Get initial Ritz-values and -vectors.
      do ist = 1, nst
        call X(hpsi)(h, gr, psi(:, :, ist), h_psi(:, :, ist), ik)
      end do

write(*, *) '**1'
      call states_blockt_mul(gr%m, st%d%dim, psi, h_psi, gram_block)
write(*, *) '**1_'
      call X(lobpcg_conj)(st%d%wfs_type, nst, nst, gram_block)

      ALLOCATE(ritz_vec(nst, nst), nst**2)
      call lalg_eigensolve(nst, gram_block, ritz_vec, eval(:, ik))

write(*, *) '**2'
      call states_block_matr_mul(gr%m, st%d%dim, psi, ritz_vec, tmp)
write(*, *) '**2_'
      call lalg_copy(np*nst, tmp(:, 1, 1), psi(:, 1, 1))

write(*, *) '**3'
      call states_block_matr_mul(gr%m, st%d%dim, h_psi, ritz_vec, tmp)
write(*, *) '**3_'
      call lalg_copy(np*nst, tmp(:, 1, 1), h_psi(:, 1, 1))
      deallocate(ritz_vec)

      ! First iteration is special because Gram matrices are smaller than in
      ! following iterations.

      ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
      call X(lobpcg_res)(nst, ik, np, h_psi, psi, eval, res)

      ! Check for convergence.
      call X(lobpcg_unconv_ev)(ik, gr%m, st%d%dim, tol, res, diffs, nuc, uc)
      conv = nst-nuc

      if(conv.eq.nst) then
        if(present(diff)) then
          diff = diffs
        end if
        converged = converged + conv
        niter     = niter+1
        cycle k_loop
      end if

      ! Apply preconditioner.
      do i = 1, nuc
        ist = uc(i)
        call X(preconditioner_apply)(pre, gr, h, res(:, :, ist), tmp(:, :, ist))
        call lalg_copy(np, tmp(:, 1, ist), res(:, 1, ist))
      end do

      ! Make residuals orthogonal to eigenstates.
      call X(lobpcg_orth_res)(nst, gr%m, st%d%dim, psi, res, nuc, uc)

      ! Orthonormalize residuals.
write(*, *) '**4'
      call X(lobpcg_orth)(gr%m, st%d%dim, st%d%wfs_type, nst, np, res, nuc, uc)
write(*, *) '**4_'

      ! Apply Hamiltonian to residuals.
      do i = 1, nuc
        ist = uc(i)
        call X(hpsi)(h, gr, res(:, :, ist), h_res(:, :, ist), ik)
      end do
      niter = niter+nuc

      ! Rayleigh-Ritz procedure.
      ALLOCATE(ritz_psi(nst, nst), nst**2)
      ALLOCATE(ritz_dir(nuc, nst), nst*nuc)
      ALLOCATE(ritz_res(nuc, nst), nst*nuc)
      ALLOCATE(gram_h(nst+nuc, nst+nuc), (nst+nuc)**2)
      ALLOCATE(gram_i(nst+nuc, nst+nuc), (nst+nuc)**2)
      ALLOCATE(ritz_vec(nst+nuc, nst), (nst+nuc)*nst)

      ! gram_h matrix.
      gram_h(1:nst, 1:nst) = R_TOTYPE(M_ZERO)
      ! (1, 1)-block: eigenvalues in gram_h.
      do ist = 1, nst
        gram_h(ist, ist) = eval(ist, ik)
      end do

      ! (1, 2)-block: (H |psi>)^T res.
write(*, *) '**5'
      call states_blockt_mul(gr%m, st%d%dim, h_psi, res, gram_h(1:nst, nst+1:nst+nuc), idx2=UC)
write(*, *) '**5_'

      ! (2, 2)-block: res^T (H res).
write(*, *) '**7'
      call states_blockt_mul(gr%m, st%d%dim, res, h_res, &
        gram_h(nst+1:nst+nuc, nst+1:nst+nuc), idx1=UC, idx2=UC)
write(*, *) '**7_'
      call X(lobpcg_conj)(st%d%wfs_type, nuc, nuc, gram_h(nst+1:nst+nuc, nst+1:nst+nuc))

      ! gram_i matrix.
      gram_i(1:nst, 1:nst)                 = R_TOTYPE(M_ZERO)
      gram_i(nst+1:nst+nuc, nst+1:nst+nuc) = R_TOTYPE(M_ZERO)

      ! Unit matrices on diagonal blocks.
      do j = 1, nst+nuc
        gram_i(j, j) = 1
      end do

      ! (1, 2)-block: <psi| res.
write(*, *) '**8'
      call states_blockt_mul(gr%m, st%d%dim, psi, res, gram_i(1:nst, nst+1:nst+nuc), idx2=UC)
write(*, *) '**8_'

      call lalg_lowest_geneigensolve(nst, nst+nuc, gram_h, gram_i, ritz_val, ritz_vec)

      ritz_psi = ritz_vec(1:nst, 1:nst)
      ritz_res = ritz_vec(nst+1:nst+nuc, 1:nst)

      ! Calculate new conjugate directions.
write(*, *) '**9'
      call states_block_matr_mul(gr%m, st%d%dim, res, ritz_res, dir, idxp=UC)
write(*, *) '**9_'
write(*, *) '**10'
      call states_block_matr_mul(gr%m, st%d%dim, h_res, ritz_res, h_dir, idxp=UC)
write(*, *) '**10_'

      ! Calculate new eigenstates and update H |psi>
write(*, *) '**11'
      call states_block_matr_mul(gr%m, st%d%dim, psi, ritz_psi, tmp)
write(*, *) '**11_'
      call lalg_copy(np*nst, tmp(:, 1, 1), psi(:, 1, 1))

      call lalg_axpy(np*nst, R_TOTYPE(M_ONE), dir(:, 1, 1), psi(:, 1, 1))

write(*, *) '**12'
      call states_block_matr_mul(gr%m, st%d%dim, h_psi, ritz_psi, tmp)
write(*, *) '**12_'
      call lalg_copy(np*nst, tmp(:, 1, 1), h_psi(:, 1, 1))
      call lalg_axpy(np*nst, R_TOTYPE(M_ONE), h_dir(:, 1, 1), h_psi(:, 1, 1))

      ! Gram matrices have to be reallocated later (because nuc changes).
      deallocate(ritz_vec, gram_h, gram_i, ritz_psi, ritz_res, ritz_dir)

      ! Copy new eigenvalues.
      call lalg_copy(nst, ritz_val, eval(:, ik))

      iter: do k = 2, maxiter
        ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
        call X(lobpcg_res)(nst, ik, np, h_psi, psi, eval, res)

        ! Check for convergence.
        call X(lobpcg_unconv_ev)(ik, gr%m, st%d%dim, tol, res, diffs, nuc, uc)
        if(nuc.eq.0) then
          exit iter
        end if

        ALLOCATE(nuc_tmp(nuc, nuc), nuc**2)
        ALLOCATE(ritz_psi(nst, nst), nst**2)
        ALLOCATE(ritz_dir(nuc, nst), nst*nuc)
        ALLOCATE(ritz_res(nuc, nst), nst*nuc)
        ! Allocate space for Gram matrices in this iterations.
        ALLOCATE(ritz_vec(nst+2*nuc, nst), nst**2+2*nst*nuc)
        ALLOCATE(gram_h(nst+2*nuc, nst+2*nuc), (nst+2*nuc)**2)
        ALLOCATE(gram_i(nst+2*nuc, nst+2*nuc), (nst+2*nuc)**2)

        ! Apply the preconditioner.
        do i = 1, nuc
          ist = uc(i)
          call X(preconditioner_apply)(pre, gr, h, res(:, :, ist), tmp(:, :, ist))
          call lalg_copy(np, tmp(:, 1, ist), res(:, 1, ist))
        end do

        ! Make residuals orthogonal to eigenstates.
        call X(lobpcg_orth_res)(nst, gr%m, st%d%dim, psi, res, nuc, uc)

        ! Orthonormalize residuals.
        call X(lobpcg_orth)(gr%m, st%d%dim, st%d%wfs_type, nst, np, res, nuc, uc)

        ! Apply Hamiltonian to residual.
        do i = 1, nuc
          ist = uc(i)
          call X(hpsi)(h, gr, res(:, :, ist), h_res(:, :, ist), ik)
        end do
        niter = niter+nuc

        ! Orthonormalize conjugate directions.
        ! Since h_dir also has to be modified (to avoid a full calculation of
        ! H dir with the new dir), we cannot use lobpcg_orth at this point.
        !call lalg_herk(nst, np, 'C', R_TOTYPE(M_ONE), dir(:, :, 1), R_TOTYPE(M_ZERO), nst_tmp)
        call states_blockt_mul(gr%m, st%d%dim, dir, dir, nuc_tmp, idx1=UC, idx2=UC)
        call X(lobpcg_conj)(st%d%wfs_type, nuc, nuc, nuc_tmp)
        call lalg_cholesky(nuc, nuc_tmp)
        call lalg_invert_upper_triangular(nuc, nuc_tmp)
        ! Fill lower triangle of nuc_tmp with zeros.
        do i = 2, nuc
          nuc_tmp(i, 1:i-1) = R_TOTYPE(M_ZERO)
        end do
        call states_block_matr_mul(gr%m, st%d%dim, dir, nuc_tmp, tmp, idxp=UC, idxr=UC)
        do i = 1, nuc
          call lalg_copy(np, tmp(:, 1, uc(i)), dir(:, 1, uc(i)))
        end do
        call states_block_matr_mul(gr%m, st%d%dim, h_dir, nuc_tmp, tmp, idxp=UC, idxr=UC)
        do i = 1, nuc
          call lalg_copy(np, tmp(:, 1, uc(i)), h_dir(:, 1, uc(i)))
        end do

        ! Rayleigh-Ritz procedure.
        ! gram_h matrix.
        gram_h(1:nst, 1:nst) = R_TOTYPE(M_ZERO)
        ! (1, 1)-block: eigenvalues on diagonal.
        do ist = 1, nst
          gram_h(ist, ist) = eval(ist, ik)
        end do

        ! (1, 2)-block: (H |psi>)^T res.
        call states_blockt_mul(gr%m, st%d%dim, h_psi, res, gram_h(1:nst, nst+1:nst+nuc), idx2=UC)

        ! (1, 3)-block: (H |psi>)^T dir.
        call states_blockt_mul(gr%m, st%d%dim, h_psi, dir, &
          gram_h(1:nst, nst+nuc+1:nst+2*nuc), idx2=UC)

        ! (2, 2)-block: res^T (H res).
        call states_blockt_mul(gr%m, st%d%dim, res, h_res, &
          gram_h(nst+1:nst+nuc, nst+1:nst+nuc), idx1=UC, idx2=UC)
        call X(lobpcg_conj)(st%d%wfs_type, nuc, nuc, gram_h(nst+1:nst+nuc, nst+1:nst+nuc))

        ! (2, 3)-block: (H res)^T dir.
        call states_blockt_mul(gr%m, st%d%dim, h_res, dir, &
          gram_h(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), idx1=UC, idx2=UC)

        ! (3, 3)-block: dir^T (H dir)
        call states_blockt_mul(gr%m, st%d%dim, dir, h_dir, &
          gram_h(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc), idx1=UC, idx2=UC)
        call X(lobpcg_conj)(st%d%wfs_type, nuc, nuc, &
          gram_h(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc))

        ! gram_i matrix.
        gram_i = R_TOTYPE(M_ZERO)
        ! Unit matrices on diagonal blocks.
        do j = 1, nst+2*nuc
          gram_i(j, j) = 1
        end do

        ! (1, 2)-block: <psi| res.
        call states_blockt_mul(gr%m, st%d%dim, psi, res, gram_i(1:nst, nst+1:nst+nuc), idx2=UC)

        ! (1, 3)-block: <psi| dir.
        call states_blockt_mul(gr%m, st%d%dim, psi, dir, &
          gram_i(1:nst, nst+nuc+1:nst+2*nuc), idx2=UC)

        ! (2, 3)-block: res^T dir
        call states_blockt_mul(gr%m, st%d%dim, res, dir, &
          gram_i(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), idx1=UC, idx2=UC)

        call lalg_lowest_geneigensolve(nst, nst+2*nuc, gram_h, gram_i, ritz_val, ritz_vec)

        ritz_psi = ritz_vec(1:nst, 1:nst)
        ritz_res = ritz_vec(nst+1:nst+nuc, 1:nst)
        ritz_dir = ritz_vec(nst+nuc+1:nst+2*nuc, 1:nst)

        ! Calculate new conjugate directions:
        ! dir <- dir ritz_dir + res ritz_res
        ! h_dir <- (H res) ritz_res + (H dir) ritz_dir
        call states_block_matr_mul(gr%m, st%d%dim, dir, ritz_dir, tmp, idxp=UC)
        call lalg_copy(np*nst, tmp(:, 1, 1), dir(:, 1, 1))
        call states_block_matr_mul_add(gr%m, st%d%dim, R_TOTYPE(M_ONE), &
          res, ritz_res, R_TOTYPE(M_ONE), dir, idxp=UC)
        call states_block_matr_mul(gr%m, st%d%dim, h_dir, ritz_dir, tmp, idxp=UC)
        call lalg_copy(np*nst, tmp(:, 1, 1), h_dir(:, 1, 1))
        call states_block_matr_mul_add(gr%m, st%d%dim, R_TOTYPE(M_ONE), &
          h_res, ritz_res, R_TOTYPE(M_ONE), h_dir, idxp=UC)

        ! Calculate new eigenstates:
        ! |psi> <- |psi> ritz_psi + dir
        ! h_psi <- (H |psi>) ritz_psi + H dir
        call states_block_matr_mul(gr%m, st%d%dim, psi, ritz_psi, tmp)
        call lalg_copy(np*nst, tmp(:, 1, 1), psi(:, 1, 1))
        call lalg_axpy(np*nst, R_TOTYPE(M_ONE), dir(:, 1, 1), psi(:, 1, 1))
        call states_block_matr_mul(gr%m, st%d%dim, h_psi, ritz_psi, tmp)
        call lalg_copy(np*nst, tmp(:, 1, 1), h_psi(:, 1, 1))
        call lalg_axpy(np*nst, R_TOTYPE(M_ONE), h_dir(:, 1, 1), h_psi(:, 1, 1))

        ! Gram matrices have to be reallocated later (because nuc changes).
        deallocate(nuc_tmp, ritz_vec, gram_h, gram_i, ritz_psi, ritz_res, ritz_dir)

        ! Copy new eigenvalues.
        call lalg_copy(nst, ritz_val, eval(:, ik))
      end do iter

      ! Check, which eigenvectors converged.
      ! Calculate latest residuals first if necessary.
      if(k.ge.maxiter) then
        call X(lobpcg_res)(nst, ik, np, h_psi, psi, eval, res)
      end if
      call X(lobpcg_unconv_ev)(ik, gr%m, st%d%dim, tol, res, diffs, nuc, uc)
      conv = nst-nuc
      if(present(diff)) then
        diff = diffs
      end if
      converged = converged + conv
    end do k_loop

    deallocate(tmp, res, h_res, dir, h_dir, h_psi, gram_block, ritz_val, diffs)

    if(verbose_) call messages_print_stress(stdout)

    call pop_sub()
  end subroutine X(eigen_solver_lobpcg)


  ! ---------------------------------------------------------
  ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
  subroutine X(lobpcg_res)(nst, ik, np, h_psi, psi, eval, res)
    integer, intent(in)    :: nst
    integer, intent(in)    :: ik
    integer, intent(in)    :: np
    R_TYPE,  intent(in)    :: h_psi(:, :, :)
    R_TYPE,  intent(in)    :: psi(:, :, :)
    FLOAT,   intent(in)    :: eval(:, :)
    R_TYPE,  intent(inout) :: res(:, :, :)

    integer :: i

    call push_sub('eigen_lobpcg_inc.Xlobpcg_res')

    do i = 1, nst
      call lalg_copy(np, h_psi(:, 1, i), res(:, 1, i))
      call lalg_axpy(np, -eval(i, ik), psi(:, 1, i), res(:, 1, i))
    end do

    call pop_sub()
  end subroutine X(lobpcg_res)


  ! ---------------------------------------------------------
  ! Recalculate set of unconverged eigenvectors.
  subroutine X(lobpcg_unconv_ev)(ik, m, dim, tol, res, diff, nuc, uc)
    integer,      intent(in)    :: ik
    type(mesh_t), intent(in)    :: m
    integer,      intent(in)    :: dim
    FLOAT,        intent(in)    :: tol
    R_TYPE,       intent(in)    :: res(:, :, :)
    FLOAT,        intent(inout) :: diff(:, :)
    integer,      intent(inout) :: nuc
    integer,      intent(inout) :: uc(:)

    integer :: i, ist, j, new_nuc
    integer :: new_uc(nuc)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_unconv_ev')

    j       = 1
    new_nuc = 0
    do i = 1, nuc
      ist = uc(i)
      diff(ist, ik) = X(states_nrm2)(m, dim, res(:, :, ist))
      if(diff(ist, ik).ge.tol) then
        new_uc(j) = ist
        new_nuc   = new_nuc+1
        j         = j+1
        write(*, *) 'UNCONVERGED EV', ist, diff(ist, ik), tol
      else
        write(*, *) '  CONVERGED EV', ist, diff(ist, ik), tol
      end if
    end do
    nuc       = new_nuc
    uc(1:nuc) = new_uc(1:nuc)

    write(*, *) 'UNCONVERGED EV:', uc(1:nuc)

    call pop_sub()
  end subroutine X(lobpcg_unconv_ev)


  ! ---------------------------------------------------------
  ! Check for convergence.
  integer function X(lobpcg_converged)(nst, ik, m, dim, tol, res, diff)
    integer,      intent(in)  :: nst
    integer,      intent(in)  :: ik
    type(mesh_t), intent(in)  :: m
    integer,      intent(in)  :: dim
    FLOAT,        intent(in)  :: tol
    R_TYPE,       intent(in)  :: res(:, :, :)
    FLOAT,        intent(out) :: diff(:, :)

    integer :: ist, conv

    call push_sub('eigen_lobpcg_inc.Xlobpcg_converged')

    conv = 0

    do ist = 1, nst
      diff(ist, ik) = X(states_nrm2)(m, dim, res(:, :, ist))
      if(diff(ist, ik).lt.tol) then
        conv = conv + 1
      end if
    end do

    X(lobpcg_converged) = conv

    call pop_sub()
  end function X(lobpcg_converged)


  ! ---------------------------------------------------------
  ! Orthonormalize the column vectors of vs.
  subroutine X(lobpcg_orth)(m, dim, wfs_type, nst, np, vs, nuc, uc)
    type(mesh_t), intent(in)    :: m
    integer,      intent(in)    :: dim
    integer,      intent(in)    :: wfs_type
    integer,      intent(in)    :: nst
    integer,      intent(in)    :: np
    R_TYPE,       intent(inout) :: vs(:, :, :)
    integer,      intent(in)    :: nuc
    integer,      intent(in)    :: uc(:)

    integer :: i

    R_TYPE  :: vv(nuc, nuc), tmp(m%np_part, dim, nst)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_orth')

    call states_blockt_mul(m, dim, vs, vs, vv, idx1=UC, idx2=UC)
    call X(lobpcg_conj)(wfs_type, nuc, nuc, vv)
    call X(write_matrix)(nuc, nuc, vv)
    call lalg_cholesky(nuc, vv)
    call lalg_invert_upper_triangular(nuc, vv)
    ! Fill lower triangle of vv with zeros.
    do i = 2, nuc
      vv(i, 1:i-1) = R_TOTYPE(M_ZERO)
    end do
    call states_block_matr_mul(m, dim, vs, vv, tmp, idxp=UC, idxr=UC)
    do i = 1, nuc
      call lalg_copy(np, tmp(:, 1, uc(i)), vs(:, 1, uc(i)))
    end do

    call pop_sub()
  end subroutine X(lobpcg_orth)


  ! ---------------------------------------------------------
  ! Make residuals orthogonal to eigenstates.
  subroutine X(lobpcg_orth_res)(nst, m, dim, psi, res, nuc, uc)
    integer,      intent(in)    :: nst
    type(mesh_t), intent(in)    :: m
    integer,      intent(in)    :: dim
    R_TYPE,       intent(in)    :: psi(:, :, :)
    R_TYPE,       intent(inout) :: res(:, :, :)
    integer,      intent(in)    :: nuc
    integer,      intent(in)    :: uc(:)

    R_TYPE :: tmp(nst, nuc)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_orth_res')

    call states_blockt_mul(m, dim, psi, res, tmp, idx2=UC)
    call states_block_matr_mul_add(m, dim, -R_TOTYPE(M_ONE), psi, &
      tmp, R_TOTYPE(M_ONE), res, idxr=UC)

    call pop_sub()
  end subroutine X(lobpcg_orth_res)


  ! ---------------------------------------------------------
  ! Calculate a <- (a + a^+)/2 (for complex wavefunctions).
  subroutine X(lobpcg_conj)(wfs_type, m, n, a)
    integer, intent(in)    :: wfs_type
    integer, intent(in)    :: m
    integer, intent(in)    :: n
    R_TYPE,  intent(inout) :: a(:, :)

    integer :: i, j
    R_TYPE  :: at(n, m)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_conj')

    if(wfs_type.eq.M_CMPLX) then
      do i = 1, m
        do j = 1, n
          at(j, i) = R_CONJ(a(i, j))
        end do
      end do
      a = (a + at)/R_TOTYPE(M_TWO)
    end if

    call pop_sub()
  end subroutine X(lobpcg_conj)


#ifdef R_TREAL
  ! Debug vv.
  subroutine zwrite_matrix(m, n, a)
    integer, intent(in) :: m, n
    CMPLX,   intent(in) :: a(:, :)

    integer   :: i, j
    character(len=20) :: fmt
    character(len=500) :: outp

    fmt = '(a,f10.3,a,f10.3,a)'

    do i = 1, m
      outp = ''
      do j = 1, n
        write(outp, fmt) trim(outp), real(a(i, j)), '+i', aimag(a(i, j)), '  '
      end do
      write(23, '(a)') trim(outp)
    end do
  end subroutine zwrite_matrix

  subroutine dwrite_matrix(m, n, a)
    integer, intent(in) :: m, n
    FLOAT,   intent(in) :: a(:, :)

    integer   :: i
    character(len=20) :: fmt

    write(fmt, *) n
    fmt = '('//trim(fmt)//'f20.8)'

    do i = 1, m
      write(23, fmt) real(a(i, :))
    end do
  end subroutine dwrite_matrix
  ! Debug ^^.
#endif


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
