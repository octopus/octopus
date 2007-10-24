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

  ! FIXME: this file needs cleanup. It contains definitely too much
  ! copy and paste.

  #include "global.h"

  ! If Fortran had subarray pointers, I would use those...
  ! Block of ritz_vec to calculate new psi.
  #define RITZ_PSI ritz_vec(1:nst, 1:nst)
  ! Block of ritz_vec to calculate new dir.
  #define RITZ_RES ritz_vec(nst+1:nst+nuc, 1:nst)
  ! Block of ritz_vec to calculate new res.
  #define RITZ_DIR ritz_vec(nst+nuc+1:nst+2*nuc, 1:nst)

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

    integer :: np               ! Number of points per state.
    integer :: nep              ! Number of eigenpairs.
    integer :: nst              ! Number of eigenstates (i. e. the blocksize).
    integer :: lnst             ! Number of local eigenstates.
    integer :: ik, ist
    integer :: st_start, st_end
    integer :: i, j, k
    integer :: conv, maxiter

    integer, target  :: nuc, uc(st%nst) ! Index set of unconverged eigenpairs.
    integer, pointer :: lnuc, luc(:)    ! Index set of local unconverged eigenpairs.

    integer           :: iunit
    logical           :: verbose_, no_bof
    logical           :: explicit_gram
    character(len=3)  :: rank_number
    character(len=10) :: file_name

    FLOAT,  allocatable :: diffs(:, :)
    R_TYPE, allocatable :: tmp(:, :, :)     ! Temporary storage of wavefunction size.
    R_TYPE, allocatable :: nuc_tmp(:, :)    ! Temporary storage of Gram block size.
    R_TYPE, allocatable :: res(:, :, :)     ! Residuals.
    R_TYPE, allocatable :: h_res(:, :, :)   ! H res.
    R_TYPE, allocatable :: dir(:, :, :)     ! Conjugate directions.
    R_TYPE, allocatable :: h_dir(:, :, :)   ! H dir.
    R_TYPE, allocatable :: h_psi(:, :, :)   ! H |psi>.
    FLOAT, pointer      :: eval(:, :)       ! Pointer to the eigenvalues.
    R_TYPE, allocatable :: gram_h(:, :)     ! Gram matrix for Hamiltonian.
    R_TYPE, allocatable :: gram_i(:, :)     ! Gram matrix for unit matrix.
    R_TYPE, allocatable :: gram_block(:, :) ! Space to construct the Gram matrix blocks.
    R_TYPE, allocatable :: ritz_vec(:, :)   ! Ritz-vectors.
    FLOAT,  allocatable :: ritz_val(:)      ! Ritz-values.
#if defined(HAVE_MPI)
    FLOAT, allocatable  :: ldiffs(:)
#endif

    call push_sub('eigen_lobpcg.Xeigen_solver_lobpcg')

    if(in_profiling_mode.and.st%parallel_in_states) then
      write(rank_number, '(i3.3)') st%mpi_grp%rank
      file_name = 'num_es.'//trim(rank_number)
      iunit = io_open(file_name, action='write', position='append')
      if(iunit.ge.0) then
        write(iunit, '(a,i3,a)') '=== Eigensolver started on node ', st%mpi_grp%rank, ' ==='
        call io_close(iunit)
      end if
    end if
      
    ! The results with explicit Gram diagonal blocks were not better, so it is switched of.
    explicit_gram = .false.

    ! Some abbreviations.
    eval     => st%eigenval
    nst      =  st%nst
    nep      =  nst*st%d%nik
    np       =  NP_PART*st%d%dim
    st_start =  st%st_start
    st_end   =  st%st_end

    ! If not running parallel, the luc and lnuc just point to the global set.
    if(st%parallel_in_states) then
      lnst = st%lnst
      ALLOCATE(luc(lnst), lnst)
      ALLOCATE(lnuc, 1)
    else
      luc  => uc
      lnuc => nuc
      lnst =  nst
    end if

    if(present(diff)) then
      diff = R_TOTYPE(M_ZERO)
    end if

    ! Some verbose output.
    verbose_ = .false.
    if(present(verbose)) verbose_ = verbose
    if(verbose_) then
      call messages_print_stress(stdout, "LOBPCG Info")
      message(1) = 'Diagonalization with the locally optimal block preconditioned'
      message(2) = 'conjugate gradient algorithm.'
      write(message(3),'(a,e8.2)') '  Tolerance: ', tol
      write(message(4),'(a,i6)')   '  Maximum number of iterations per block of eigenstates: ', &
        niter
      message(5) = ''
      call write_info(5)
    end if

    ! FIXME: it should be possible to allocate only vectors of length NP.
    ALLOCATE(diffs(nst, st%d%nik), nep)
    ALLOCATE(tmp(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(res(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(h_res(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(dir(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(h_dir(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(h_psi(NP_PART, st%d%dim, st_start:st_end), NP_PART*st%d%dim*lnst)
    ALLOCATE(gram_block(nst, nst), nst**2)
    ALLOCATE(ritz_val(nst), nst)

    maxiter = niter
    niter   = 0

    ! Iterate over all k-points.
    k_loop: do ik = 1, st%d%nik
      ! At the beginning, all eigenvectors are considered unconverged.
      nuc = nst
      do ist = 1, nst
        uc(ist) = ist
      end do
      if(st%parallel_in_states) then
        lnuc = lnst
        do ist = 1, lnuc
          luc(ist) = ist+st_start-1
        end do
      end if

      ! Orthonormalize initial vectors.
      no_bof = .false.
      call X(lobpcg_orth)(gr%m, st, st%X(psi)(:, :, :, ik), nuc, uc, lnuc, luc, tmp, no_bof)
      if(no_bof) then
        message(1) = 'Bad problem: orthonormalization of initial vectors failed.'
        call write_warning(1)
      end if

      ! Get initial Ritz-values and -vectors.
      do ist = st_start, st_end
        call X(hpsi)(h, gr, st%X(psi)(:, :, ist, ik), h_psi(:, :, ist), ik)
      end do
      niter = niter+lnst

      call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), h_psi, gram_block, symm=.true.)

      ALLOCATE(ritz_vec(nst, nst), nst**2)
      no_bof = .false.
      call lalg_eigensolve(nst, gram_block, ritz_vec, eval(:, ik), bof=no_bof)
      if(no_bof) then
        message(1) = 'Bad problem: Rayleigh-Ritz procedure for initial vectors failed.'
        call write_warning(1)
      end if
      call states_block_matr_mul(gr%m, st, st%X(psi)(:, :, :, ik), ritz_vec, tmp)
      call lalg_copy(np*lnst, tmp(:, 1, st_start), st%X(psi)(:, 1, st_start, ik))
      call states_block_matr_mul(gr%m, st, h_psi, ritz_vec, tmp)
      call lalg_copy(np*lnst, tmp(:, 1, st_start), h_psi(:, 1, st_start))
      deallocate(ritz_vec)

      ! First iteration is special because Gram matrices are smaller than in
      ! following iterations.

      ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
      call X(lobpcg_res)(st, ik, gr%m, h_psi, st%X(psi)(:, :, :, ik), eval, res)

      ! Check for convergence.
      call X(lobpcg_unconv_ev)(ik, gr%m, st, tol, res, diffs, nuc, uc, lnuc, luc)
      conv = nst-nuc

      ! If converged, go to the next k-point.
      if(conv.eq.nst) then
        if(present(diff)) then
#if defined(HAVE_MPI)
          ! Exchange differences.
          if(st%parallel_in_states) then
            ALLOCATE(ldiffs(lnst), lnst)
            ldiffs = diffs(st_start:st_end, ik)
            call lmpi_gen_alltoallv(lnst, ldiffs, i, diffs(:, ik), st%mpi_grp)
            deallocate(ldiffs)
          end if
#endif
          diff = diffs
        end if
        converged = converged + conv
        call X(lobpcg_info)(st, verbose_, ik, 1, nuc, uc, diffs(:, ik))
        cycle k_loop
      end if

      ! Apply preconditioner.
      do i = 1, lnuc
        ist = luc(i)
        call X(preconditioner_apply)(pre, gr, h, res(:, :, ist), tmp(:, :, ist))
        call lalg_copy(np, tmp(:, 1, ist), res(:, 1, ist))
      end do

      ! Orthonormalize residuals.
      no_bof = .false.
      call X(lobpcg_orth)(gr%m, st, res, nuc, uc, lnuc, luc, tmp, no_bof)
      if(no_bof) then
        message(1) = 'Bad problem: orthonormalization of residuals failed.'
        write(message(2), '(a)') 'in iteration #     1'
        call write_warning(2)
        cycle k_loop
      end if

      ! Apply Hamiltonian to residuals.
      do i = 1, lnuc
        ist = luc(i)
        call X(hpsi)(h, gr, res(:, :, ist), h_res(:, :, ist), ik)
      end do
      niter = niter+lnuc

      ! Rayleigh-Ritz procedure.
      ALLOCATE(gram_h(nst+nuc, nst+nuc), (nst+nuc)**2)
      ALLOCATE(gram_i(nst+nuc, nst+nuc), (nst+nuc)**2)
      ALLOCATE(ritz_vec(nst+nuc, nst), (nst+nuc)*nst)

      ! gram_h matrix.
      ! (1, 1)-block:
      if(explicit_gram) then
        call states_blockt_mul(gr%m, st, h_psi, st%X(psi)(:, :, :, ik), gram_h(1:nst, 1:nst))
      else
        ! Eigenvalues in gram_h.
        gram_h(1:nst, 1:nst) = R_TOTYPE(M_ZERO)
        do ist = 1, nst
          gram_h(ist, ist) = eval(ist, ik)
        end do
      end if

      ! (1, 2)-block: (H |psi>)^+ res.
      call states_blockt_mul(gr%m, st, h_psi, res, gram_h(1:nst, nst+1:nst+nuc), xpsi2=UC)

      ! (2, 2)-block: res^+ (H res).
      call states_blockt_mul(gr%m, st, res, h_res, &
        gram_h(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC, symm=.true.)

      ! gram_i matrix.
      ! Diagonal blocks:
      if(explicit_gram) then
        call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), st%X(psi)(:, :, :, ik), gram_i(1:nst, 1:nst))
        call states_blockt_mul(gr%m, st, res, res, &
          gram_i(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC)
      else
        ! Unit matrices on diagonal blocks.
        gram_i(1:nst, 1:nst)                 = R_TOTYPE(M_ZERO)
        gram_i(nst+1:nst+nuc, nst+1:nst+nuc) = R_TOTYPE(M_ZERO)
        do j = 1, nst+nuc
          gram_i(j, j) = 1
        end do
      end if

      ! (1, 2)-block: <psi| res.
      call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), res, gram_i(1:nst, nst+1:nst+nuc), xpsi2=UC)


      no_bof = .false.
      call profiling_in(C_PROFILING_LOBPCG_ESOLVE)
      call lalg_lowest_geneigensolve(nst, nst+nuc, gram_h, gram_i, ritz_val, ritz_vec, bof=no_bof)
      call profiling_out(C_PROFILING_LOBPCG_ESOLVE)
      if(no_bof) then
        message(1) = 'Bad problem: Rayleigh-Ritz procedure failed'
        write(message(2), '(a)') 'in iteration #     1'
        call write_warning(2)
        cycle k_loop
      end if

      ! Calculate new conjugate directions.
      call states_block_matr_mul(gr%m, st, res, RITZ_RES, dir, xpsi=UC)
      call states_block_matr_mul(gr%m, st, h_res, RITZ_RES, h_dir, xpsi=UC)

      ! Calculate new eigenstates and update H |psi>
      call states_block_matr_mul(gr%m, st, st%X(psi)(:, :, :, ik), RITZ_PSI, tmp)
      call lalg_copy(np*lnst, tmp(:, 1, st_start), st%X(psi)(:, 1, st_start, ik))

      do ist = st_start, st_end ! Leave this loop, otherwise xlf90 crashes.
        call lalg_axpy(np, R_TOTYPE(M_ONE), dir(:, 1, ist), st%X(psi)(:, 1, ist, ik))
      end do

      call states_block_matr_mul(gr%m, st, h_psi, RITZ_PSI, tmp)
      call lalg_copy(np*lnst, tmp(:, 1, st_start), h_psi(:, 1, st_start))
      call lalg_axpy(np*lnst, R_TOTYPE(M_ONE), h_dir(:, 1, st_start), h_psi(:, 1, st_start))

      ! Gram matrices have to be reallocated later (because nuc changes).
      deallocate(ritz_vec, gram_h, gram_i) !, ritz_psi, ritz_res, ritz_dir)

      ! Copy new eigenvalues.
      call lalg_copy(nst, ritz_val, eval(:, ik))

      ! This is the big iteration loop.
      iter: do k = 2, maxiter-1 ! One iteration was performed to get initial Ritz-vectors.
        ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
        call X(lobpcg_res)(st, ik, gr%m, h_psi, st%X(psi)(:, :, :, ik), eval, res)

        ! Check for convergence. If converged, quit the eigenpair iteration loop.
        call X(lobpcg_unconv_ev)(ik, gr%m, st, tol, res, diffs, nuc, uc, lnuc, luc)
        if(nuc.eq.0) then
          exit iter
        end if

        ALLOCATE(nuc_tmp(nuc, nuc), nuc**2)
        ! Allocate space for Gram matrices in this iterations.
        ALLOCATE(ritz_vec(nst+2*nuc, nst), nst**2+2*nst*nuc)
        ALLOCATE(gram_h(nst+2*nuc, nst+2*nuc), (nst+2*nuc)**2)
        ALLOCATE(gram_i(nst+2*nuc, nst+2*nuc), (nst+2*nuc)**2)

        ! Apply the preconditioner.
        do i = 1, lnuc
          ist = luc(i)
          call X(preconditioner_apply)(pre, gr, h, res(:, :, ist), tmp(:, :, ist))
          call lalg_copy(np, tmp(:, 1, ist), res(:, 1, ist))
        end do

        ! Orthonormalize residuals.
        no_bof = .false.
        call X(lobpcg_orth)(gr%m, st, res, nuc, uc, lnuc, luc, tmp, no_bof)
        if(no_bof) then
          message(1) = 'Bad problem: orthonormalization of residuals failed.'
        write(message(2), '(a)') 'in iteration #     1'
          call write_warning(2)
          exit iter
        end if

        ! Apply Hamiltonian to residual.
        do i = 1, lnuc
          ist = luc(i)
          call X(hpsi)(h, gr, res(:, :, ist), h_res(:, :, ist), ik)
        end do
        niter = niter+lnuc

        ! Orthonormalize conjugate directions.
        ! Since h_dir also has to be modified (to avoid a full calculation of
        ! H dir with the new dir), we cannot use lobpcg_orth at this point.
        call states_blockt_mul(gr%m, st, dir, dir, nuc_tmp, xpsi1=UC, xpsi2=UC, symm=.true.)
        call profiling_in(C_PROFILING_LOBPCG_CHOL)
        no_bof = .false.
        call lalg_cholesky(nuc, nuc_tmp, bof=no_bof)
        call profiling_out(C_PROFILING_LOBPCG_CHOL)
        if(no_bof) then
          message(1) = 'Bad problem: orthonormalization of conjugate directions failed'
          write(message(2), '(a,i6)') 'in iteration #', k
          call write_warning(2)
          ! Set directions to zero.
          ! FIXME: they should not be included in the subspace at all in this case.
          ! (the code has to be cleaned up anyway, so this can be done then).
          dir   = R_TOTYPE(M_ZERO)
          h_dir = R_TOTYPE(M_ZERO)
        else
          call profiling_in(C_PROFILING_LOBPCG_INV)
          call lalg_invert_upper_triangular(nuc, nuc_tmp)
          call profiling_out(C_PROFILING_LOBPCG_INV)
          ! Fill lower triangle of nuc_tmp with zeros.
          do i = 2, nuc
            nuc_tmp(i, 1:i-1) = R_TOTYPE(M_ZERO)
          end do
          call states_block_matr_mul(gr%m, st, dir, nuc_tmp, tmp, xpsi=UC, xres=UC)
          do i = 1, lnuc
            call lalg_copy(np, tmp(:, 1, luc(i)), dir(:, 1, luc(i)))
          end do
          call states_block_matr_mul(gr%m, st, h_dir, nuc_tmp, tmp, xpsi=UC, xres=UC)
          do i = 1, lnuc
            call lalg_copy(np, tmp(:, 1, luc(i)), h_dir(:, 1, luc(i)))
          end do
        end if

        ! Rayleigh-Ritz procedure.
        ! gram_h matrix.
        if(explicit_gram) then
          call states_blockt_mul(gr%m, st, h_psi, st%X(psi)(:, :, :, ik), gram_h(1:nst, 1:nst))
        else
          ! (1, 1)-block: eigenvalues on diagonal.
          gram_h(1:nst, 1:nst) = R_TOTYPE(M_ZERO)
          do ist = 1, nst
            gram_h(ist, ist) = eval(ist, ik)
          end do
        end if

        ! (1, 2)-block: (H |psi>)^+ res.
        call states_blockt_mul(gr%m, st, h_psi, res, gram_h(1:nst, nst+1:nst+nuc), xpsi2=UC)

        ! (1, 3)-block: (H |psi>)^+ dir.
        call states_blockt_mul(gr%m, st, h_psi, dir, &
          gram_h(1:nst, nst+nuc+1:nst+2*nuc), xpsi2=UC)

        ! (2, 2)-block: res^+ (H res).
        call states_blockt_mul(gr%m, st, res, h_res, &
          gram_h(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC, symm=.true.)

        ! (2, 3)-block: (H res)^+ dir.
        call states_blockt_mul(gr%m, st, h_res, dir, &
          gram_h(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)

        ! (3, 3)-block: dir^+ (H dir)
        call states_blockt_mul(gr%m, st, dir, h_dir, &
          gram_h(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC, symm=.true.)

        ! gram_i matrix.
        ! Diagonal blocks.
        if(explicit_gram) then
          call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), st%X(psi)(:, :, :, ik), gram_i(1:nst, 1:nst))
          call states_blockt_mul(gr%m, st, res, res, &
            gram_i(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC)
          call states_blockt_mul(gr%m, st, dir, dir, &
            gram_i(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)
        else
          ! Unit matrices on diagonal blocks.
          gram_i = R_TOTYPE(M_ZERO)
          do j = 1, nst+2*nuc
            gram_i(j, j) = 1
          end do
        end if

        ! (1, 2)-block: <psi| res.
        call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), res, gram_i(1:nst, nst+1:nst+nuc), xpsi2=UC)

        ! (1, 3)-block: <psi| dir.
        call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), dir, &
          gram_i(1:nst, nst+nuc+1:nst+2*nuc), xpsi2=UC)

        ! (2, 3)-block: res^+ dir.
        call states_blockt_mul(gr%m, st, res, dir, &
          gram_i(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)
        call profiling_in(C_PROFILING_LOBPCG_ESOLVE)
        call lalg_lowest_geneigensolve(nst, nst+2*nuc, gram_h, gram_i, ritz_val, ritz_vec, bof=no_bof)
        call profiling_out(C_PROFILING_LOBPCG_ESOLVE)
        if(no_bof) then
          message(1) = 'Bad problem: Rayleigh-Ritz procedure failed'
          write(message(2), '(a,i6)') 'in iteration #', k
          call write_warning(2)
          exit iter
        end if

        ! Calculate new conjugate directions:
        ! dir <- dir ritz_dir + res ritz_res
        ! h_dir <- (H res) ritz_res + (H dir) ritz_dir
        call states_block_matr_mul(gr%m, st, dir, RITZ_DIR, tmp, xpsi=UC)
        call lalg_copy(np*lnst, tmp(:, 1, st_start), dir(:, 1, st_start))
        call states_block_matr_mul_add(gr%m, st, R_TOTYPE(M_ONE), &
          res, RITZ_RES, R_TOTYPE(M_ONE), dir, xpsi=UC)
        call states_block_matr_mul(gr%m, st, h_dir, RITZ_DIR, tmp, xpsi=UC)
        call lalg_copy(np*lnst, tmp(:, 1, st_start), h_dir(:, 1, st_start))
        call states_block_matr_mul_add(gr%m, st, R_TOTYPE(M_ONE), &
          h_res, RITZ_RES, R_TOTYPE(M_ONE), h_dir, xpsi=UC)

        ! Calculate new eigenstates:
        ! |psi> <- |psi> ritz_psi + dir
        ! h_psi <- (H |psi>) ritz_psi + H dir
        call states_block_matr_mul(gr%m, st, st%X(psi)(:, :, :, ik), RITZ_PSI, tmp)
        call lalg_copy(np*lnst, tmp(:, 1, st_start), st%X(psi)(:, 1, st_start, ik))
        do ist = st_start, st_end ! Leave this loop, otherwise xlf90 crashes.
          call lalg_axpy(np, R_TOTYPE(M_ONE), dir(:, 1, ist), st%X(psi)(:, 1, ist, ik))
        end do
        call states_block_matr_mul(gr%m, st, h_psi, RITZ_PSI, tmp)
        call lalg_copy(np*lnst, tmp(:, 1, st_start), h_psi(:, 1, st_start))
        call lalg_axpy(np*lnst, R_TOTYPE(M_ONE), h_dir(:, 1, st_start), h_psi(:, 1, st_start))

        ! Gram matrices have to be reallocated later (because nuc changes).
        deallocate(nuc_tmp, ritz_vec, gram_h, gram_i)

        !, ritz_psi, ritz_res, ritz_dir)

        ! Copy new eigenvalues.
        call lalg_copy(nst, ritz_val, eval(:, ik))
      end do iter

      ! Check, which eigenvectors converged.
      ! Calculate latest residuals first if necessary.
      if(k.ge.maxiter) then
        call X(lobpcg_res)(st, ik, gr%m, h_psi, st%X(psi)(:, :, :, ik), eval, res)
      end if
      call X(lobpcg_unconv_ev)(ik, gr%m, st, tol, res, diffs, nuc, uc, lnuc, luc)

      conv = nst-nuc
      if(present(diff)) then
        ! Exchange differences.
#if defined(HAVE_MPI)
        if(st%parallel_in_states) then
          ALLOCATE(ldiffs(lnst), lnst)
          ldiffs = diffs(st_start:st_end, ik)
          call lmpi_gen_alltoallv(lnst, ldiffs, i, diffs(:, ik), st%mpi_grp)
          deallocate(ldiffs)
        end if
#endif
        diff = diffs
      end if
      converged = converged + conv
      call X(lobpcg_info)(st, verbose_, ik, k, nuc, uc, diffs(:, ik))
    end do k_loop

#if defined(HAVE_MPI)
      ! Exchange number of matrix-vector operations.
      if(st%parallel_in_states) then
        i = niter
        call MPI_Debug_In(st%mpi_grp%comm, C_MPI_ALLREDUCE)
        call MPI_Allreduce(i, niter, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
        call MPI_Debug_Out(st%mpi_grp%comm, C_MPI_ALLREDUCE)
      end if
#endif

    if(verbose_) call messages_print_stress(stdout)

    deallocate(tmp, res, h_res, dir, h_dir, h_psi, gram_block, ritz_val, diffs)

    if(st%parallel_in_states) then
      deallocate(luc, lnuc)
    end if

    call pop_sub()
  end subroutine X(eigen_solver_lobpcg)


  ! ---------------------------------------------------------
  ! In verbose mode: write convergence information for k-point ik.
  subroutine X(lobpcg_info)(st, verbose, ik, niter, nuc, uc, diff)
    type(states_t), intent(in) :: st
    logical,        intent(in) :: verbose
    integer,        intent(in) :: ik
    integer,        intent(in) :: niter
    integer,        intent(in) :: nuc
    integer,        intent(in) :: uc(:)
    FLOAT,          intent(in) :: diff(:)

    integer :: i, niter_total
    logical :: mask(st%nst)

    call push_sub('eigen_lobpcg_inc.X(lobpcg_info)')

    if(verbose) then
#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        call MPI_Debug_In(st%mpi_grp%comm, C_MPI_ALLREDUCE)
        call MPI_Allreduce(niter, niter_total, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
        call MPI_Debug_Out(st%mpi_grp%comm, C_MPI_ALLREDUCE)
      end if
#else
      niter_total = niter
#endif
      call X(lobpcg_conv_mask)(nuc, uc, mask)
      write(message(1), '(a,i5,a,i5,a)') 'Result for k = ', ik, ' after ', &
        niter_total, ' block iterations:'
      call write_info(1)
      do i = 1, st%nst
        if(mask(i)) then
          write(message(1), '(a,i5,a,e8.2,a)') '  Eigenstate #', i, &
            ':     converged [Res = ', diff(i), ']'
          call write_info(1)
        else
          write(message(1), '(a,i5,a,e8.2,a)') '  Eigenstate #', i, &
            ': not converged [Res = ', diff(i), ']'
          call write_info(1)
        end if
      end do
    end if

    call pop_sub()
  end subroutine X(lobpcg_info)


  ! ---------------------------------------------------------
  ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
  subroutine X(lobpcg_res)(st, ik, m, h_psi, psi, eval, res)
    type(states_t), intent(in)    :: st
    integer,        intent(in)    :: ik
    type(mesh_t),   intent(in)    :: m
    R_TYPE,         intent(in)    :: h_psi(m%np_part, st%d%dim, st%st_start:st%st_end)
    R_TYPE,         intent(in)    :: psi(m%np_part, st%d%dim, st%st_start:st%st_end)
    FLOAT,          intent(in)    :: eval(:, :)
    R_TYPE,         intent(inout) :: res(m%np_part, st%d%dim, st%st_start:st%st_end)

    integer :: ist, np

    call push_sub('eigen_lobpcg_inc.Xlobpcg_res')

    np = m%np_part*st%d%dim
    do ist = st%st_start, st%st_end
      call lalg_copy(np, h_psi(:, 1, ist), res(:, 1, ist))
      call lalg_axpy(np, -eval(ist, ik), psi(:, 1, ist), res(:, 1, ist))
    end do

    call pop_sub()
  end subroutine X(lobpcg_res)


  ! ---------------------------------------------------------
  ! Recalculate set of unconverged eigenvectors.
  subroutine X(lobpcg_unconv_ev)(ik, m, st, tol, res, diff, nuc, uc, lnuc, luc)
    integer,        intent(in)    :: ik
    type(mesh_t),   intent(in)    :: m
    type(states_t), intent(in)    :: st
    FLOAT,          intent(in)    :: tol
    R_TYPE,         intent(in)    :: res(m%np_part, st%d%dim, st%st_start:st%st_end)
    FLOAT,          intent(inout) :: diff(:, :)
    integer,        intent(inout) :: nuc
    integer,        intent(inout) :: uc(:)
    integer,        intent(inout) :: lnuc
    integer,        intent(inout) :: luc(:)

    integer           :: i, ist, j, new_nuc, iunit
    integer           :: new_uc(nuc)
    character(len=3)  :: rank_number
    character(len=10) :: file_name

    call push_sub('eigen_lobpcg_inc.Xlobpcg_unconv_ev')

    j       = 1
    new_nuc = 0
    do i = 1, lnuc
      ist = luc(i)
      diff(ist, ik) = X(states_nrm2)(m, st%d%dim, res(:, :, ist))
      if(diff(ist, ik).ge.tol) then
        new_uc(j) = ist
        new_nuc   = new_nuc+1
        j         = j+1
      end if
    end do
    lnuc        = new_nuc
    luc(1:lnuc) = new_uc(1:lnuc)

#if defined(HAVE_MPI)
    ! Update set of unconverged vectors an all nodes.
    if(st%parallel_in_states) then
      call lmpi_gen_alltoallv(lnuc, luc, nuc, uc, st%mpi_grp)
      if(in_profiling_mode) then
        write(rank_number, '(i3.3)') st%mpi_grp%rank
        file_name = 'num_es.'//trim(rank_number)
        iunit = io_open(file_name, action='write', position='append')
        if(iunit.ge.0) then
          write(iunit, '(a,i6,a,i6)') 'k = ', ik, '       nuc = ', lnuc
          call io_close(iunit)
        end if
      end if
    end if
#else
    uc = uc ! Avoid unused variable warning.
#endif

    call pop_sub()
  end subroutine X(lobpcg_unconv_ev)


  ! ---------------------------------------------------------
  ! Returns a mask with mask(i) = .false. for eigenvector i unconverged.
  subroutine X(lobpcg_conv_mask)(nuc, uc, mask)
    integer, intent(in)  :: nuc
    integer, intent(in)  :: uc(:)
    logical, intent(out) :: mask(:)

    call push_sub('eigen_lobpcg_inc.X(lobpcg_conv_mask)')

    mask     = .true.
    mask(UC) = .false.

    call pop_sub()
  end subroutine X(lobpcg_conv_mask)


  ! ---------------------------------------------------------
  ! Orthonormalize the column vectors of vs.
  subroutine X(lobpcg_orth)(m, st, vs, nuc, uc, lnuc, luc, tmp, chol_failure)
    type(mesh_t),   intent(in)    :: m
    type(states_t), intent(in)    :: st
    R_TYPE,         intent(inout) :: vs(m%np_part, st%d%dim, st%st_start:st%st_end)
    integer,        intent(in)    :: nuc, lnuc
    integer,        intent(in)    :: uc(:), luc(:)
    R_TYPE,         intent(out)   :: tmp(m%np_part, st%d%dim, st%st_start:st%st_end)
    logical,        intent(out)   :: chol_failure

    integer             :: i
    R_TYPE, allocatable :: vv(:, :)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_orth')

    chol_failure = .false.
    ALLOCATE(vv(nuc, nuc), nuc**2)

    call states_blockt_mul(m, st, vs, vs, vv, xpsi1=UC, xpsi2=UC, symm=.true.)
    call profiling_in(C_PROFILING_LOBPCG_CHOL)
    call lalg_cholesky(nuc, vv, bof=chol_failure)
    call profiling_out(C_PROFILING_LOBPCG_CHOL)
    if(chol_failure) then ! Failure in Choleksy decomposition.
      return
    end if
    call profiling_in(C_PROFILING_LOBPCG_INV)
    call lalg_invert_upper_triangular(nuc, vv)
    call profiling_out(C_PROFILING_LOBPCG_INV)
    ! Fill lower triangle of vv with zeros.
    do i = 2, nuc
      vv(i, 1:i-1) = R_TOTYPE(M_ZERO)
    end do
    call states_block_matr_mul(m, st, vs, vv, tmp, xpsi=UC, xres=UC)
    do i = 1, lnuc
      call lalg_copy(m%np_part*st%d%dim, tmp(:, 1, luc(i)), vs(:, 1, luc(i)))
    end do

    deallocate(vv)

    call pop_sub()
  end subroutine X(lobpcg_orth)


  ! ---------------------------------------------------------
  ! Make residuals orthogonal to eigenstates (this routine seems to
  ! be unnessecary).
  subroutine X(lobpcg_orth_res)(nst, m, st, psi, res, nuc, uc)
    integer,        intent(in)    :: nst
    type(mesh_t),   intent(in)    :: m
    type(states_t), intent(in)    :: st
    R_TYPE,         intent(in)    :: psi(:, :, :)
    R_TYPE,         intent(inout) :: res(:, :, :)
    integer,        intent(in)    :: nuc
    integer,        intent(in)    :: uc(:)

    R_TYPE :: tmp(nst, nuc)

    call push_sub('eigen_lobpcg_inc.Xlobpcg_orth_res')

    call states_blockt_mul(m, st, psi, res, tmp, xpsi2=UC)
    call states_block_matr_mul_add(m, st, -R_TOTYPE(M_ONE), psi, &
      tmp, R_TOTYPE(M_ONE), res, xres=UC)

    call pop_sub()
  end subroutine X(lobpcg_orth_res)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
