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
  
! Implementation of the locally optimal block preconditioned conjugate-
! gradients algorithm.

! Index set of unconverged eigenvectors.
#define UC uc(1:nuc)

! ---------------------------------------------------------
!> Driver for the LOBPCG eigensolver that performs a per-block,
!! per-k-point iteration.
subroutine X(eigensolver_lobpcg)(gr, st, hm, pre, tol, niter, converged, ik, diff, block_size)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(in)    :: ik
  integer,                intent(inout) :: converged
  FLOAT,                  intent(out)   :: diff(:) !< (1:st%nst)
  integer,                intent(in)    :: block_size
  
  integer            :: ib, psi_start, psi_end, constr_start, constr_end, bs
  integer            :: n_matvec, conv, maxiter, iblock
#ifdef HAVE_MPI
  integer            :: outcount
  FLOAT, allocatable :: ldiff(:)
#endif
  
  PUSH_SUB(X(eigensolver_lobpcg))
  
  bs = block_size
  
  maxiter = niter
  niter   = 0
  
  diff(1:st%nst) = M_ZERO
  
  iblock = 0
  
  if(mpi_grp_is_root(mpi_world)) then
    call loct_progress_bar(st%nst*(ik - 1), st%nst*st%d%nik)
  end if
  
  ! Iterate over all blocks.
  do ib = st%st_start, st%st_end, bs
    iblock    = iblock+1
    psi_start = ib
    psi_end   = ib+bs-1
    
    if(psi_end > st%st_end) then
      psi_end = st%st_end
    end if
    constr_start = st%st_start
    constr_end   = ib-1
    
    n_matvec = maxiter
    
    if(constr_end >= constr_start) then
      call X(lobpcg)(gr, st, hm, psi_start, psi_end, st%X(psi)(:, :, psi_start:psi_end, ik), &
        constr_start, constr_end, &
        ik, pre, tol, n_matvec, conv, diff, &
        constr = st%X(psi)(:, :, constr_start:constr_end, ik))
    else
      call X(lobpcg)(gr, st, hm, psi_start, psi_end, st%X(psi)(:, :, psi_start:psi_end, ik), &
        constr_start, constr_end, ik, pre, tol, n_matvec, conv, diff)
    end if
    
    niter     = niter + n_matvec
    converged = converged + conv  
    
    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) + psi_end, st%nst*st%d%nik)
    end if
  end do
  
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    SAFE_ALLOCATE(ldiff(1:st%lnst))
    ldiff(1:st%lnst) = diff(st%st_start:st%st_end)
    call lmpi_gen_allgatherv(st%lnst, ldiff, outcount, diff, st%mpi_grp)
    SAFE_DEALLOCATE_A(ldiff)
  end if
#endif
  
  POP_SUB(X(eigensolver_lobpcg))
end subroutine X(eigensolver_lobpcg)


! ---------------------------------------------------------
!> Locally optimal block preconditioned conjugate gradient algorithm.
!! For details, see:
!!
!! A. Knyazev. Toward the Optimal Preconditioned Eigensolver: Locally
!! Optimal Block Preconditioned Conjugate Gradient Method. SIAM
!! Journal on Scientific Computing, 23(2):517-541, 2001.
!!
!! A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchin-
!! nikov. Block Locally Optimal Preconditioned Eigenvalue Xolvers
!! (BLOPEX) in hypre and PETSc. SIAM Journal of Scientific Computing,
!! 2007.
!!
!! There is also a wiki page at
!! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:LOBPCG
subroutine X(lobpcg)(gr, st, hm, st_start, st_end, psi, constr_start, constr_end,  &
  ik, pre, tol, niter, converged, diff, constr)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: st_start
  integer,                intent(in)    :: st_end
  R_TYPE, target,         intent(inout) :: psi(gr%mesh%np_part, st%d%dim, st_start:st_end)
  integer,                intent(in)    :: constr_start
  integer,                intent(in)    :: constr_end
  integer,                intent(in)    :: ik
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(out)   :: converged
  FLOAT,                  intent(inout) :: diff(:) !< (1:st%nst)
  R_TYPE, optional,       intent(in)    :: constr(gr%mesh%np_part, st%d%dim, constr_start:constr_end)

  integer :: nps   !< Number of points per state.
  integer :: nst   !< Number of eigenstates (i.e. the blocksize).
  integer :: lnst  !< Number of local eigenstates.

  integer :: ist, i, j, iter, blks, maxiter, nconstr

  integer, target      :: nuc                 !< Index set of unconverged eigenpairs.
  integer, pointer     :: uc(:), lnuc, luc(:) !< Index set of local unconverged eigenpairs.
  integer, allocatable :: all_ev(:), all_constr(:)

#ifdef HAVE_MPI
  integer :: lnconstr
  integer, allocatable :: lall_constr(:)
#endif

  integer           :: hash_table_size
  logical           :: no_bof, found
  logical           :: explicit_gram
  R_TYPE            :: beta
  type(iihash_t)    :: all_ev_inv

  R_TYPE, pointer             :: ritz_psi(:, :), ritz_res(:, :), ritz_dir(:, :)
  R_TYPE, allocatable         :: tmp(:, :, :)     !< Temporary storage of wavefunction size.
  R_TYPE, allocatable         :: nuc_tmp(:, :)    !< Temporary storage of Gram block size.
  R_TYPE, allocatable         :: res(:, :, :)     !< Residuals.
  R_TYPE, allocatable         :: h_res(:, :, :)   !< H res.
  R_TYPE, allocatable         :: dir(:, :, :)     !< Conjugate directions.
  R_TYPE, allocatable         :: h_dir(:, :, :)   !< H dir.
  R_TYPE, allocatable         :: h_psi(:, :, :)   !< H |psi>.
  FLOAT,  allocatable          :: eval(:)         !< The eigenvalues of the current block.
  R_TYPE, allocatable         :: gram_h(:, :)     !< Gram matrix for Hamiltonian.
  R_TYPE, allocatable         :: gram_i(:, :)     !< Gram matrix for unit matrix.
  R_TYPE, allocatable         :: gram_block(:, :) !< Space to construct the Gram matrix blocks.
  R_TYPE, allocatable, target :: ritz_vec(:, :)   !< Ritz-vectors.
  type(batch_t) :: psib, hpsib
  logical :: there_are_constraints
  
  PUSH_SUB(X(lobpcg))

  if(constr_end >= constr_start) then
    there_are_constraints = .true.
    ASSERT(present(constr) .eqv. there_are_constraints)
  else
    there_are_constraints = .false.
  end if

  ! The results with explicit Gram diagonal blocks were not better, so it is switched off.
  explicit_gram = .false.

  ! Abbreviations.
  nps = gr%mesh%np_part*st%d%dim

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    lnst     = st_end-st_start+1
    lnconstr = constr_end-constr_start+1
    SAFE_ALLOCATE(luc(1:lnst))
    SAFE_ALLOCATE(lnuc)
    call MPI_Allreduce(lnst, nst, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
    call MPI_Allreduce(lnconstr, nconstr, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
#endif
  else
    nconstr =  constr_end-constr_start+1
    nst     =  st_end-st_start+1
    lnuc    => nuc
    lnst    =  nst
  end if
  ASSERT( (nconstr > 0) .eqv. there_are_constraints)
  SAFE_ALLOCATE(uc(1:nst))
  if(.not.st%parallel_in_states) then
    nconstr = 0
    luc => uc
  end if

  if(there_are_constraints) then
    SAFE_ALLOCATE(all_constr(1:nconstr))
    if(st%parallel_in_states) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(lall_constr(1:lnconstr))
      do i = 1, lnconstr
        lall_constr(i) = i + constr_start-1
      end do

      call lmpi_gen_allgatherv(lnconstr, lall_constr, nconstr, all_constr, st%mpi_grp)
      SAFE_DEALLOCATE_A(lall_constr)
#endif
    else
      do i = 1, nconstr
        all_constr(i) = i+constr_start-1
      end do
    end if
  end if

  SAFE_ALLOCATE(  tmp(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(  res(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(h_res(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(  dir(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(h_dir(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(h_psi(1:gr%mesh%np_part, 1:st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(gram_block(1:nst, 1:nst))
  SAFE_ALLOCATE(eval(nst))

  ! Set them to zero, otherwise behaviour may be slightly nondeterministic.
  tmp   = R_TOTYPE(M_ZERO)
  res   = R_TOTYPE(M_ZERO)
  h_res = R_TOTYPE(M_ZERO)
  dir   = R_TOTYPE(M_ZERO)
  h_dir = R_TOTYPE(M_ZERO)
  h_psi = R_TOTYPE(M_ZERO)

  maxiter = niter
  niter   = 0

  ! At the beginning, all eigenvectors are considered unconverged.
  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    lnuc = lnst
    do ist = 1, lnuc
      luc(ist) = ist+st_start-1
    end do
    call lmpi_gen_allgatherv(lnuc, luc, nuc, uc, st%mpi_grp)
#endif
  else
    nuc = nst
    do ist = 1, nuc
      uc(ist) = ist + st_start-1
    end do
  end if

  ! Auxiliary index maps, required because we compact the distributed blocks of each node
  ! in the Rayleigh-Ritz into one big block.
  ! all_ev:     {1, ..., nst} -> {1, ..., st%nst} (blocksize to number of eigenvectors).
  ! all_ev_inv: {1, ..., st%nst} -> {1, ..., nst} (the reverese of all_ev).
  SAFE_ALLOCATE(all_ev(1:nst))
  all_ev = uc
  hash_table_size = max(3, st%nst) ! Minimum size of hash table is 3.
  call iihash_init(all_ev_inv, hash_table_size)
  do ist = 1, nst
    call iihash_insert(all_ev_inv, all_ev(ist), ist)
  end do

  ! Apply the constraints to the initial vectors.
  if(nconstr > 0) then
    call X(lobpcg_apply_constraints)(st_start, st_end, psi, nuc, uc)
  end if

  ! Orthonormalize initial vectors.
  no_bof = .false.
  call X(lobpcg_orth)(st_start, st_end, psi(:, :, :), no_bof)

  if(no_bof) then
    message(1) = 'Problem: orthonormalization of initial vectors failed.'
    call messages_warning(1)
  end if

  ! Get initial Ritz-values and -vectors.
  call batch_init(psib, st%d%dim, st_start, st_end, st%X(psi)(:, :, st_start:, ik))
  call batch_init(hpsib, st%d%dim, st_start, st_end, h_psi(:, :, st_start:))

  call X(hamiltonian_apply_batch)(hm, gr%der, psib, hpsib, ik)
  
  call batch_end(psib)
  call batch_end(hpsib)

  niter = niter+lnst
  call X(blockt_mul)(psi(:, :, :), h_psi, gram_block, xpsi1=all_ev, xpsi2=all_ev, symm=.true.)

  SAFE_ALLOCATE(ritz_vec(1:nst, 1:nst))
  no_bof = .false.
  ritz_vec = gram_block
  call lalg_eigensolve(nst, ritz_vec, eval, bof=no_bof)

  if(no_bof) then
    message(1) = 'Problem: Rayleigh-Ritz procedure for initial vectors failed.'
    call messages_warning(1)
  end if
  call X(block_matr_mul)(psi, ritz_vec, tmp, xpsi=all_ev, xres=all_ev)
  call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), psi(:, :, st_start:))
  call X(block_matr_mul)(h_psi, ritz_vec, tmp, xpsi=all_ev, xres=all_ev)
  call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), h_psi(:, :, st_start:))
  SAFE_DEALLOCATE_A(ritz_vec)

  ! This is the big iteration loop.
  iteration: do iter = 1, maxiter-1 ! One iteration was used up to get initial Ritz-vectors.
    ! Calculate residuals: res(ist, ik) <- H psi(ist, ik) - eval(ist, ik) psi(ist, ik).
    call X(lobpcg_res)()

    ! Check for convergence. If converged, quit the eigenpair iteration loop.
    call X(lobpcg_unconv_ev)

    if(nuc.eq.0) then
      exit iteration
    end if

    SAFE_ALLOCATE(nuc_tmp(1:nuc, 1:nuc))
    ! Allocate space for Gram matrices in this iterations.
    ! blks says if we have one or two additional blocks in the subspace
    ! (i. e. only residuals or residuals (1st iteration) and conjugate
    ! directions (subsequent iterations).
    if(iter > 1) then
      blks = 2
    else
      blks = 1
    end if
    SAFE_ALLOCATE(ritz_vec(1:nst+blks*nuc, 1:nst))
    SAFE_ALLOCATE(  gram_h(1:nst+blks*nuc, 1:nst+blks*nuc))
    SAFE_ALLOCATE(  gram_i(1:nst+blks*nuc, 1:nst+blks*nuc))
    ritz_psi => ritz_vec(1:nst, 1:nst)
    ritz_res => ritz_vec(nst+1:nst+nuc, 1:nst)
    if(iter > 1) then
      ritz_dir => ritz_vec(nst+nuc+1:nst+2*nuc, 1:nst)
    end if

    ! Apply the preconditioner.
    do i = 1, lnuc
      ist = luc(i)
      call X(preconditioner_apply)(pre, gr, hm, ik, res(:, :, ist), tmp(:, :, ist))
      call lalg_copy(gr%mesh%np_part, st%d%dim, tmp(:, :, ist), res(:, :, ist))
    end do

    ! Apply the constraints to the residuals.
    if(nconstr > 0) then
      call X(lobpcg_apply_constraints)(st_start, st_end, res, nuc, uc)
    end if

    ! Orthonormalize residuals.
    no_bof = .false.
    call X(lobpcg_orth)(st_start, st_end, res, no_bof)
    ! FIXME: a proper restart should be initiated here.

    if(no_bof) then
      message(1) = 'Big problem: orthonormalization of residuals failed.'
      message(2) = 'Quitting eigensolver iteration.'
      write(message(3), '(a,i6)') 'in iteration #', iter
      call messages_warning(3)
      exit iteration
    end if

    ! Apply Hamiltonian to residuals.

    if(lnuc > 0) then
      call batch_init(psib, st%d%dim, lnuc)
      call batch_init(hpsib, st%d%dim, lnuc)
    end if
    
    do i = 1, lnuc
      ist = luc(i)
      call batch_add_state(psib, ist, res(:, :, ist))
      call batch_add_state(hpsib, ist, h_res(:, :, ist))
    end do

    if(lnuc > 0) then
      call X(hamiltonian_apply_batch)(hm, gr%der, psib, hpsib, ik)
    end if

    niter = niter + lnuc

    call batch_end(psib)
    call batch_end(hpsib)
      
    ! Orthonormalize conjugate directions in all but the first iteration.
    ! Since h_dir also has to be modified (to avoid a full calculation of
    ! H dir with the new dir), we cannot use lobpcg_orth at this point.
    if(iter > 1) then
      call X(blockt_mul)(dir, dir, nuc_tmp, xpsi1=UC, xpsi2=UC, symm=.true.)
      call profiling_in(C_PROFILING_LOBPCG_CHOL)
      no_bof = .false.
      call lalg_cholesky(nuc, nuc_tmp, bof=no_bof)
      call profiling_out(C_PROFILING_LOBPCG_CHOL)

      if(no_bof) then
        message(1) = 'Problem: orthonormalization of conjugate directions failed'
        write(message(2), '(a,i6)') 'in iteration #', iter
        call messages_warning(2)
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
          nuc_tmp(i, 1:i - 1) = R_TOTYPE(M_ZERO)
        end do
        call X(block_matr_mul)(dir, nuc_tmp, tmp, xpsi=UC, xres=UC)
        do i = 1, lnuc
          call lalg_copy(gr%mesh%np_part, st%d%dim, tmp(:, :, luc(i)), dir(:, :, luc(i)))
        end do
        call X(block_matr_mul)(h_dir, nuc_tmp, tmp, xpsi=UC, xres=UC)
        do i = 1, lnuc
          call lalg_copy(gr%mesh%np_part, st%d%dim, tmp(:, :, luc(i)), h_dir(:, :, luc(i)))
        end do
      end if
    end if

    ! Rayleigh-Ritz procedure.
    ! gram_h matrix.
    if(explicit_gram) then
      call X(blockt_mul)(h_psi, psi(:, :, :), gram_h(1:nst, 1:nst), xpsi1=all_ev, xpsi2=all_ev)
    else
      ! (1, 1)-block: eigenvalues on diagonal.
      gram_h(1:nst, 1:nst) = R_TOTYPE(M_ZERO)
      do ist = 1, nst
        gram_h(ist, ist) = eval(ist)
      end do
    end if

    ! (1, 2)-block: (H |psi>)^+ res.
    call X(blockt_mul)(h_psi, res, gram_h(1:nst, nst+1:nst+nuc), xpsi1=all_ev, xpsi2=UC)

    ! (2, 2)-block: res^+ (H res).
    call X(blockt_mul)(res, h_res, gram_h(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC, symm=.true.)

    if(iter > 1) then
      ! (1, 3)-block: (H |psi>)^+ dir.
      call X(blockt_mul)(h_psi, dir, gram_h(1:nst, nst+nuc+1:nst+2*nuc), xpsi1=all_ev, xpsi2=UC)

      ! (2, 3)-block: (H res)^+ dir.
      call X(blockt_mul)(h_res, dir, gram_h(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)

      ! (3, 3)-block: dir^+ (H dir)
      call X(blockt_mul)(dir, h_dir, gram_h(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC, symm=.true.)
    end if

    ! gram_i matrix.
    ! Diagonal blocks.
    if(explicit_gram) then
      call X(blockt_mul)(psi(:, :, :), psi(:, :, :), gram_i(1:nst, 1:nst), xpsi1=all_ev, xpsi2=all_ev)
      call X(blockt_mul)(res, res, gram_i(nst+1:nst+nuc, nst+1:nst+nuc), xpsi1=UC, xpsi2=UC)
      call X(blockt_mul)(dir, dir, gram_i(nst+nuc+1:nst+2*nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)
    else
      ! Unit matrices on diagonal blocks.
      gram_i = R_TOTYPE(M_ZERO)
      do j = 1, nst+blks*nuc
        gram_i(j, j) = 1
      end do
    end if

    ! (1, 2)-block: <psi| res.
    call X(blockt_mul)(psi(:, :, :), res, gram_i(1:nst, nst+1:nst+nuc), xpsi1=all_ev, xpsi2=UC)

    if(iter > 1) then
      ! (1, 3)-block: <psi| dir.
      call X(blockt_mul)(psi(:, :, :), dir, gram_i(1:nst, nst+nuc+1:nst+2*nuc), xpsi1=all_ev, xpsi2=UC)

      ! (2, 3)-block: res^+ dir.
      call X(blockt_mul)(res, dir, gram_i(nst+1:nst+nuc, nst+nuc+1:nst+2*nuc), xpsi1=UC, xpsi2=UC)
    end if

    call profiling_in(C_PROFILING_LOBPCG_ESOLVE)
    no_bof = .false.
    call lalg_lowest_geneigensolve(nst, nst+blks*nuc, gram_h, gram_i, eval, ritz_vec, bof=no_bof)
    call profiling_out(C_PROFILING_LOBPCG_ESOLVE)

    if(no_bof) then
      message(1) = 'Problem: Rayleigh-Ritz procedure failed'
      write(message(2), '(a,i6)') 'in iteration #', iter
      call messages_warning(2)
      exit iteration
    end if

    ! Calculate new conjugate directions:
    ! dir <- dir ritz_dir + res ritz_res
    ! h_dir <- (H res) ritz_res + (H dir) ritz_dir
    if(iter > 1) then
      call X(block_matr_mul)(dir, ritz_dir, tmp, xpsi=UC, xres=all_ev)
      call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), dir(:, :, st_start:))
      call X(block_matr_mul)(h_dir, ritz_dir, tmp, xpsi=UC, xres=all_ev)
      call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), h_dir(:, :, st_start:))
      beta = R_TOTYPE(M_ONE)
    else
      beta = R_TOTYPE(M_ZERO)
    end if
    call X(block_matr_mul_add)(R_TOTYPE(M_ONE), res, ritz_res, beta, dir, xpsi=UC, xres=all_ev)
    call X(block_matr_mul_add)(R_TOTYPE(M_ONE), h_res, ritz_res, beta, h_dir, xpsi=UC, xres=all_ev)

    ! Calculate new eigenstates:
    ! |psi> <- |psi> ritz_psi + dir
    ! h_psi <- (H |psi>) ritz_psi + H dir
    call X(block_matr_mul)(psi(:, :, :), ritz_psi, tmp, xpsi=all_ev, xres=all_ev)
    call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), psi(:, :, st_start:))
    do ist = st_start, st_end ! Leave this loop, otherwise xlf90 crashes.
      call lalg_axpy(gr%mesh%np_part, st%d%dim, M_ONE, dir(:, :, ist), psi(:, :, ist))
    end do
    call X(block_matr_mul)(h_psi, ritz_psi, tmp, xpsi=all_ev, xres=all_ev)
    call lalg_copy(gr%mesh%np_part, st%d%dim, lnst, tmp(:, :, st_start:), h_psi(:, :, st_start:))
    call lalg_axpy(gr%mesh%np_part, st%d%dim, lnst, R_TOTYPE(M_ONE), h_dir(:, :, st_start:), h_psi(:, :, st_start:))

    ! Gram matrices have to be reallocated later (because nuc changes).
    SAFE_DEALLOCATE_A(nuc_tmp)
    SAFE_DEALLOCATE_A(ritz_vec)
    SAFE_DEALLOCATE_A(gram_h)
    SAFE_DEALLOCATE_A(gram_i)

    ! Copy new eigenvalues.
    do i = 1, nst
      st%eigenval(all_ev(i), ik) = eval(i)
    end do
  end do iteration

  ! Check, which eigenvectors converged.
  ! Calculate latest residuals first if necessary.
  if(iter >= maxiter) then
    call X(lobpcg_res)()
  end if
  call X(lobpcg_unconv_ev)()

  converged = nst - nuc

#if defined(HAVE_MPI)
  ! Exchange number of matrix-vector operations.
  if(st%parallel_in_states) then
    i = niter
    call MPI_Allreduce(i, niter, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
  end if
#endif

  if(there_are_constraints) then
    SAFE_DEALLOCATE_A(all_constr)
  end if
  SAFE_DEALLOCATE_A(all_ev)
  SAFE_DEALLOCATE_P(uc)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(h_res)
  SAFE_DEALLOCATE_A(dir)
  SAFE_DEALLOCATE_A(h_dir)
  SAFE_DEALLOCATE_A(h_psi)
  SAFE_DEALLOCATE_A(gram_block)
  SAFE_DEALLOCATE_A(eval)
  call iihash_end(all_ev_inv)

  if(st%parallel_in_states) then
    SAFE_DEALLOCATE_P(luc)
    SAFE_DEALLOCATE_P(lnuc)
  end if

  POP_SUB(X(lobpcg))

contains

  ! ---------------------------------------------------------
  !> Calculate residuals
  subroutine X(lobpcg_res)()
    integer :: ist, iev
    integer :: idim, ip

    PUSH_SUB(X(lobpcg).X(lobpcg_res))

    do ist = st_start, st_end
      iev = iihash_lookup(all_ev_inv, ist, found)
      ASSERT(found)
     
      forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np) 
        res(ip, idim, ist) = h_psi(ip, idim, ist) - eval(iev)*psi(ip, idim, ist)
      end forall
    end do

    POP_SUB(X(lobpcg).X(lobpcg_res))
  end subroutine X(lobpcg_res)


  ! ---------------------------------------------------------
  !> Recalculate set of unconverged eigenvectors.
  subroutine X(lobpcg_unconv_ev)()
    integer           :: i, ist, j, new_nuc
    integer           :: new_uc(nuc)

    PUSH_SUB(X(lobpcg).X(lobpcg_unconv_ev))

    j       = 1
    new_nuc = 0
    do i = 1, lnuc
      ist = luc(i)
      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, ist))
      if(diff(ist) >= tol) then
        new_uc(j) = ist
        new_nuc   = new_nuc+1
        j         = j+1
      end if
    end do
    lnuc        = new_nuc
    luc(1:lnuc) = new_uc(1:lnuc)

#if defined(HAVE_MPI)
    ! Update set of unconverged vectors on all nodes.
    if(st%parallel_in_states) then
      call lmpi_gen_allgatherv(lnuc, luc, nuc, uc, st%mpi_grp)
    end if
#endif
    POP_SUB(X(lobpcg).X(lobpcg_unconv_ev))
  end subroutine X(lobpcg_unconv_ev)


  ! ---------------------------------------------------------
  !> Returns a mask with mask(i) = .false. for eigenvector i unconverged.
  subroutine X(lobpcg_conv_mask)(mask)
    logical, intent(out) :: mask(:)

    PUSH_SUB(X(lobpcg).X(lobpcg_conv_mask))

    mask     = .true.
    mask(UC) = .false.

    POP_SUB(X(lobpcg).X(lobpcg_conv_mask))
  end subroutine X(lobpcg_conv_mask)


  ! ---------------------------------------------------------
  !> Orthonormalize the column vectors of vs.
  subroutine X(lobpcg_orth)(v_start, v_end, vs, chol_failure)
    integer,        intent(in)    :: v_start
    integer,        intent(in)    :: v_end
    R_TYPE,         intent(inout) :: vs(:, :, v_start:) !< (gr%mesh%np_part, st%d%dim, v_start:v_end)
    logical,        intent(out)   :: chol_failure

    integer             :: i
    R_TYPE, allocatable :: vv(:, :)

    PUSH_SUB(X(lobpcg).X(lobpcg_orth))

    chol_failure = .false.
    SAFE_ALLOCATE(vv(1:nuc, 1:nuc))
    call states_blockt_mul(gr%mesh, st, v_start, v_start, vs, vs, vv, xpsi1=UC, xpsi2=UC, symm=.true.)
    call profiling_in(C_PROFILING_LOBPCG_CHOL)
    call lalg_cholesky(nuc, vv, bof=chol_failure)
    call profiling_out(C_PROFILING_LOBPCG_CHOL)
    if(chol_failure) then ! Failure in Cholesky decomposition.
      POP_SUB(X(lobpcg).X(lobpcg_orth))
      return
    end if
    call profiling_in(C_PROFILING_LOBPCG_INV)
    call lalg_invert_upper_triangular(nuc, vv)
    call profiling_out(C_PROFILING_LOBPCG_INV)
    ! Fill lower triangle of vv with zeros.
    do i = 2, nuc
      vv(i, 1:i-1) = R_TOTYPE(M_ZERO)
    end do

    call X(block_matr_mul)(vs, vv, tmp, xpsi=UC, xres=UC)
    do i = 1, lnuc
      call lalg_copy(gr%mesh%np_part, st%d%dim, tmp(:, :, luc(i)), vs(:, :, luc(i)))
    end do
    SAFE_DEALLOCATE_A(vv)

    POP_SUB(X(lobpcg).X(lobpcg_orth))
  end subroutine X(lobpcg_orth)


  ! ---------------------------------------------------------
  subroutine X(lobpcg_apply_constraints)(vs_start, vs_end, vs, nidx, idx)
    integer, intent(in)    :: vs_start
    integer, intent(in)    :: vs_end
    R_TYPE,  intent(inout) :: vs(:, :, vs_start:) !< (gr%mesh%np_part, st%d%dim, vs_start:vs_end)
    integer, intent(in)    :: nidx
    integer, intent(in)    :: idx(:)
    
    R_TYPE              :: det
    R_TYPE, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "LOBPCG_CONSTRAINTS")
    PUSH_SUB(X(lobpcg).X(lobpcg_apply_constraints))

    SAFE_ALLOCATE(tmp1(1:nconstr, 1:nconstr))
    SAFE_ALLOCATE(tmp2(1:nconstr, 1:nidx))
    SAFE_ALLOCATE(tmp3(1:nconstr, 1:nidx))

    call states_blockt_mul(gr%mesh, st, constr_start, constr_start, &
      constr, constr, tmp1, xpsi1=all_constr, xpsi2=all_constr)
    det = lalg_inverter(nconstr, tmp1, invert=.true.)
    call states_blockt_mul(gr%mesh, st, constr_start, vs_start, &
      constr, vs, tmp2, xpsi1=all_constr, xpsi2=idx(1:nidx))
    call lalg_gemm(nconstr, nidx, nconstr, R_TOTYPE(M_ONE), tmp1, tmp2, R_TOTYPE(M_ZERO), tmp3)
    call states_block_matr_mul_add(gr%mesh, st, -R_TOTYPE(M_ONE), constr_start, vs_start, &
      constr, tmp3, R_TOTYPE(M_ONE), vs, xpsi=all_constr, xres=idx(1:nidx))

    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(tmp3)
    POP_SUB(X(lobpcg).X(lobpcg_apply_constraints))
    call profiling_out(prof)

  end subroutine X(lobpcg_apply_constraints)


  ! ---------------------------------------------------------
  subroutine X(blockt_mul)(psi1, psi2, res, xpsi1, xpsi2, symm)
    R_TYPE,            intent(in)  :: psi1(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    R_TYPE,            intent(in)  :: psi2(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    R_TYPE,            intent(out) :: res(:, :)
    integer,           intent(in)  :: xpsi1(:)
    integer,           intent(in)  :: xpsi2(:)
    logical, optional, intent(in)  :: symm

    PUSH_SUB(X(lobpcg).X(blockt_mul))

    call states_blockt_mul(gr%mesh, st, st_start, st_start, &
      psi1, psi2, res, xpsi1=xpsi1, xpsi2=xpsi2, symm=symm)

    POP_SUB(X(lobpcg).X(blockt_mul))
  end subroutine X(blockt_mul)


  ! ---------------------------------------------------------
  subroutine X(block_matr_mul_add)(alpha, psi, matr, beta, res, xpsi, xres)
    R_TYPE,  intent(in)    :: alpha
    R_TYPE,  intent(in)    :: psi(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    R_TYPE,  intent(in)    :: matr(:, :)
    R_TYPE,  intent(in)    :: beta
    R_TYPE,  intent(inout) :: res(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    integer, intent(in)    :: xpsi(:)
    integer, intent(in)    :: xres(:)

    PUSH_SUB(X(lobpcg).X(block_matr_mul_add))

    call states_block_matr_mul_add(gr%mesh, st, alpha, st_start, st_start, &
      psi, matr, beta, res, xpsi=xpsi, xres=xres)

    POP_SUB(X(lobpcg).X(block_matr_mul_add))
  end subroutine X(block_matr_mul_add)


  ! ---------------------------------------------------------
  subroutine X(block_matr_mul)(psi, matr, res, xpsi, xres)
    R_TYPE,  intent(in)  :: psi(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    R_TYPE,  intent(in)  :: matr(:, :)
    R_TYPE,  intent(out) :: res(:, :, st_start:) !< (gr%mesh%np_part, st%d%dim, st_start:st_end)
    integer, intent(in)  :: xpsi(:)
    integer, intent(in)  :: xres(:)

    PUSH_SUB(X(lobpcg).X(block_matr_mul))

    call states_block_matr_mul_add(gr%mesh, st, R_TOTYPE(M_ONE), st_start, st_start, &
      psi, matr, R_TOTYPE(M_ZERO), res, xpsi=xpsi, xres=xres)

    POP_SUB(X(lobpcg).X(block_matr_mul))
  end subroutine X(block_matr_mul)
end subroutine X(lobpcg)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
