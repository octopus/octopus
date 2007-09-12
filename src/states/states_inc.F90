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


! ---------------------------------------------------------
! Multiplication of two blocks of states:
! res <- psi1(idx1)^+ * psi2(idx2) with the index sets idx1 and idx2.
subroutine X(states_blockt_mul)(mesh, st, psi1, psi2, res, idx1, idx2, symm)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE, target,    intent(in)  :: psi1(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE, target,    intent(in)  :: psi2(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(out) :: res(:, :)
  integer, optional, intent(in)  :: idx1(:)
  integer, optional, intent(in)  :: idx2(:)
  logical, optional, intent(in)  :: symm    ! Indicates if res(j, i) can be calculated as
                                            ! res(i, j)*.

  logical              :: symm_
  integer              :: i, j
  integer              :: m, n
  integer, allocatable :: idx1_(:), idx2_(:)
  R_TYPE, pointer      :: blk1(:, :, :), blk2(:, :, :)

  call profiling_in(C_PROFILING_LOBPCG_BLOCKT)
  call push_sub('states_inc.Xstates_blockt_mul')

  symm_ = .false.
  if(present(symm)) then
    symm_ = symm
  end if

  if(present(idx1)) then
    m = ubound(idx1, 1)
    ALLOCATE(idx1_(m), m)
    idx1_ = idx1
  else
    m = st%nst
    ALLOCATE(idx1_(m), m)
    do i = 1, m
      idx1_(i) = i
    end do
  end if
  if(present(idx2)) then
    n = ubound(idx2, 1)
    ALLOCATE(idx2_(n), n)
    idx2_ = idx2
  else
    n = st%nst
    ALLOCATE(idx2_(n), n)
    do i = 1, n
      idx2_(i) = i
    end do
  end if

  ! FIXME: write something sensible here!
  ! Quick and dirty parallelization. In the parallel case, we simply copy
  ! all states to all nodes. Otherwise, we put pointers to the input blocks.
  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Allocate space for the blocks.
    ALLOCATE(blk1(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)
    ALLOCATE(blk2(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)

    call X(states_gather)(mesh, st, psi1, blk1)
    call X(states_gather)(mesh, st, psi2, blk2)
#endif
  else
    blk1 => psi1
    blk2 => psi2
  end if

  ! FIXME: This works with domain parallelization but the dotps should be blocked.
  if(symm_) then
    do i = 1, m
      res(i, i) = X(states_dotp)(mesh, st%d%dim, blk1(:, :, idx1_(i)), blk2(:, :, idx2_(i)))
      do j = i+1, n
        res(i, j) = X(states_dotp)(mesh, st%d%dim, blk1(:, :, idx1_(i)), blk2(:, :, idx2_(j)))
        res(j, i) = R_CONJ(res(i, j))
      end do
    end do
  else
    do i = 1, m
      do j = 1, n
        res(i, j) = X(states_dotp)(mesh, st%d%dim, blk1(:, :, idx1_(i)), blk2(:, :, idx2_(j)))
      end do
    end do
  end if
  deallocate(idx1_, idx2_)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    deallocate(blk1, blk2)
  end if
#endif

  call pop_sub()
  call profiling_out(C_PROFILING_LOBPCG_BLOCKT)
end subroutine X(states_blockt_mul)


! ---------------------------------------------------------
! Gather all states on all nodes. out has to be of sufficient size.
subroutine X(states_gather)(mesh, st, in, out)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: mesh
  R_TYPE,         intent(in)  :: in(:, :, :)
  R_TYPE,         intent(out) :: out(:, :, :)

  integer              :: i, mpi_err
  integer, allocatable :: sendcnts(:), sdispls(:), recvcnts(:), rdispls(:)

  call push_sub('states_inc.Xstates_gather')

#if defined(HAVE_MPI)
  ALLOCATE(sendcnts(st%mpi_grp%size), st%mpi_grp%size)
  ALLOCATE(sdispls(st%mpi_grp%size), st%mpi_grp%size)
  ALLOCATE(recvcnts(st%mpi_grp%size), st%mpi_grp%size)
  ALLOCATE(rdispls(st%mpi_grp%size), st%mpi_grp%size)

  sendcnts   = mesh%np_part*st%d%dim*st%st_num(st%mpi_grp%rank)
  sdispls    = 0
  recvcnts   = st%st_num*mesh%np_part*st%d%dim
  rdispls(1) = 0
  do i = 2, st%mpi_grp%size
    rdispls(i) = rdispls(i-1) + recvcnts(i-1)
  end do

  call MPI_Debug_In(st%mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(in(:, 1, 1), sendcnts, sdispls, R_MPITYPE, &
    out(:, 1, 1), recvcnts, rdispls, R_MPITYPE, st%mpi_grp%comm, mpi_err)
  call MPI_Debug_Out(st%mpi_grp%comm, C_MPI_ALLTOALLV)

  deallocate(sendcnts, sdispls, recvcnts, rdispls)
#endif

  call pop_sub()
end subroutine X(states_gather)


! ---------------------------------------------------------
! Multiplication of block of states with indices idxp by matrix and
! update columns with idxr in the result.
! res(idxr) <- psi(idxp) * matr.
subroutine X(states_block_matr_mul)(mesh, st, psi, matr, res, idxp, idxr)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE,            intent(in)  :: psi(:, :, :)
  R_TYPE,            intent(in)  :: matr(:, :)
  R_TYPE,            intent(out) :: res(:, :, :)
  integer, optional, intent(in)  :: idxp(:), idxr(:)

  call push_sub('states_inc.Xstates_block_matr_mul')

  call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi, matr, &
    R_TOTYPE(M_ZERO), res, idxp, idxr)

  call pop_sub()
end subroutine X(states_block_matr_mul)


! ---------------------------------------------------------
! Multiplication of block of states by matrix plus block of states
! (with the corresponding column indices):
! res(idxr) <- alpha * psi(idx) * matr + beta * res(idxr).
subroutine X(states_block_matr_mul_add)(mesh, st, alpha, psi, matr, beta, res, idxp, idxr)
  type(mesh_t),      intent(in)    :: mesh
  type(states_t),    intent(in)    :: st
  R_TYPE,            intent(in)    :: alpha
  R_TYPE, target,    intent(in)    :: psi(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(in)    :: matr(:, :)
  R_TYPE,            intent(in)    :: beta
  R_TYPE, target,    intent(inout) :: res(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  integer, optional, intent(in)    :: idxp(:), idxr(:)

  R_TYPE               :: tmp
  integer              :: res_col, psi_col, matr_col, i, j, k, idim
  integer, allocatable :: idxp_(:), idxr_(:)
  R_TYPE, pointer      :: blkp(:, :, :), blkr(:, :, :)

  call profiling_in(C_PROFILING_LOBPCG_BLOCK_MATR)
  call push_sub('states_inc.Xstates_block_matr_add')

  if(present(idxp)) then
    psi_col = ubound(idxp, 1)
    ALLOCATE(idxp_(psi_col), psi_col)
    idxp_ = idxp
  else
    psi_col = st%nst
    ALLOCATE(idxp_(psi_col), psi_col)
    do i = 1, psi_col
      idxp_(i) = i
    end do
  end if
  if(present(idxr)) then
    res_col = ubound(idxr, 1)
    ALLOCATE(idxr_(res_col), res_col)
    idxr_ = idxr
  else
    res_col = st%nst
    ALLOCATE(idxr_(res_col), res_col)
    do i = 1, res_col
      idxr_(i) = i
    end do
  end if

  matr_col = ubound(matr, 2)

  ! FIXME: remove this shit!
  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Allocate space for the block.
    ALLOCATE(blkp(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)
    ALLOCATE(blkr(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)

    call X(states_gather)(mesh, st, res, blkr)
    call X(states_gather)(mesh, st, psi, blkp)
#endif
  else
    blkp => psi
    blkr => res
  end if

  ! FIXME: does not work with domain parallelization.  
  ! We can shortcut for beta being zero.
  if(beta.eq.R_TOTYPE(M_ZERO)) then
    call profiling_in(C_PROFILING_LOBPCG_LOOP)
    do j = 1, matr_col
      do idim = 1, st%d%dim
        do i = 1, mesh%np
          blkr(i, idim, idxr_(j)) = R_TOTYPE(M_ZERO)
          do k = 1, psi_col
            blkr(i, idim, idxr_(j)) = blkr(i, idim, idxr_(j)) + blkp(i, idim, idxp_(k))*matr(k, j)
          end do
          blkr(i, idim, idxr_(j)) = alpha*blkr(i, idim, idxr_(j))
        end do
      end do
    end do
    call profiling_out(C_PROFILING_LOBPCG_LOOP)
 ! And also for alpha=beta=1.
 else if(alpha.eq.R_TOTYPE(M_ONE).and.beta.eq.R_TOTYPE(M_ONE)) then
    call profiling_in(C_PROFILING_LOBPCG_LOOP)
    do j = 1, matr_col
      do idim = 1, st%d%dim
        do i = 1, mesh%np
          do k = 1, psi_col
            blkr(i, idim, idxr_(j)) = blkr(i, idim, idxr_(j)) + blkp(i, idim, idxp_(k))*matr(k, j)
          end do
        end do
      end do
    end do
    call profiling_out(C_PROFILING_LOBPCG_LOOP)
 ! The general case.
  else
    call profiling_in(C_PROFILING_LOBPCG_LOOP)
    do j = 1, matr_col
      do idim = 1, st%d%dim
        do i = 1, mesh%np
          tmp = R_TOTYPE(M_ZERO)
          do k = 1, psi_col
            tmp = tmp + blkp(i, idim, idxp_(k))*matr(k, j)
          end do
          blkr(i, idim, idxr_(j)) = alpha*tmp + beta*blkr(i, idim, idxr_(j))
        end do
      end do
    end do
    call profiling_out(C_PROFILING_LOBPCG_LOOP)
  end if

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Copy result back.
    res(:, :, :) = blkr(:, :, st%st_start:st%st_end)
    deallocate(blkp, blkr)
  end if
#endif

  call pop_sub()
  call profiling_out(C_PROFILING_LOBPCG_BLOCK_MATR)
end subroutine X(states_block_matr_mul_add)


! ---------------------------------------------------------
! Orthonormalizes nst orbital in mesh m
subroutine X(states_gram_schmidt1)(nst, m, dim, psi, start)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  integer, optional, intent(in)    :: start

  integer :: p, q, stst, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT1)
  call push_sub('states_inc.Xstates_gram_schmidt1')

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = X(states_dotp)(m, dim, psi(:,:, q), psi(:,:, p))
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), psi(:, idim, p))
      end do
    end do

    nrm2 = X(states_nrm2)(m, dim, psi(:,:, p))
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, psi(:, idim, p))
    end do
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT1)
end subroutine X(states_gram_schmidt1)

! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits to do only the orthogonalization (no normalization).
! And one can pass an extra optional argument, mask, which:
!  - on input, if mask(p) = .true., the p-orbital is not used.
!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
subroutine X(states_gram_schmidt2)(nst, m, dim, psi, phi, normalize, mask)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     ! phi(m%np_part, dim)
  logical, optional, intent(in)    :: normalize
  logical, optional, intent(inout) :: mask(:)      ! nst

  logical :: normalize_
  integer :: q, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT2)
  call push_sub('states_inc.Xstates_gram_schmidt2')

  do q = 1, nst
    if(present(mask)) then
      if(mask(q)) cycle
    end if
    ss = X(states_dotp)(m, dim, psi(:,:, q), phi)
    if(abs(ss) > M_EPSILON) then
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), phi(:, idim))
      end do
    else
      if(present(mask)) mask(q) = .true.
    end if
  end do

  normalize_ = .false.
  if(present(normalize)) normalize_ = normalize
  if(normalize) then
    nrm2 = X(states_nrm2)(m, dim, phi)
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, phi(:, idim))
    end do
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT2)
end subroutine X(states_gram_schmidt2)


! ---------------------------------------------------------
R_TYPE function X(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f1(:,:), f2(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim))
  end do

  call pop_sub()

end function X(states_dotp)


! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(m, dim, psi)
  type(mesh_t),    intent(in)  :: m
  integer,         intent(in)  :: dim
  R_TYPE,          intent(out) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim

  call push_sub('states_inc.Xstates_normalize_orbital')

  norm = X(states_nrm2) (m, dim, psi)
  norm = sqrt(norm)

  do idim = 1, dim
    !$omp parallel workshare
    psi(1:m%np, idim) = psi(1:m%np, idim)/norm
    !$omp end parallel workshare
  end do

  call pop_sub()
end subroutine X(states_normalize_orbital)


! ---------------------------------------------------------
FLOAT function X(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_nrm2')

  nrm2 = M_ZERO
  do idim = 1, dim
    nrm2 = nrm2 + X(mf_nrm2)(m, f(:, idim))**2
  end do
  nrm2 = sqrt(nrm2)

  call pop_sub()

end function X(states_nrm2)


! ---------------------------------------------------------
FLOAT function X(states_residue)(m, dim, hf, e, f) result(r)
  type(mesh_t),      intent(in)  :: m
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:), f(:,:)
  FLOAT,             intent(in)  :: e

  R_TYPE, allocatable :: res(:,:)

  call push_sub('states_inc.Xstates_residue')

  ALLOCATE(res(m%np_part, dim), m%np_part*dim)

  !$omp parallel workshare
  res(1:m%np, 1:dim) = hf(1:m%np, 1:dim) - e*f(1:m%np, 1:dim)
  !$omp end parallel workshare

  r = X(states_nrm2)(m, dim, res)
  deallocate(res)

  call pop_sub()

end function X(states_residue)


! ---------------------------------------------------------
! The routine calculates the expectation value of the momentum 
! operator
! <p> = < phi*(ist, k) | -i \nabla | phi(ist, ik) >
!
! Note, the blas routines cdotc, zdotc take care of complex 
! conjugation *. Therefore we pass phi directly.
! ---------------------------------------------------------
subroutine X(states_calc_momentum)(gr, st)
  type(grid_t),   intent(inout) :: gr
  type(states_t), intent(inout) :: st

  integer :: idim, ist, ik, i
  CMPLX               :: expect_val_p
  R_TYPE, allocatable :: grad(:,:,:)  

  call push_sub('states_inc.Xstates_calc_momentum')

  ALLOCATE(grad(NP, st%d%dim, NDIM), NP*st%d%dim*NDIM)

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end

      do idim = 1, st%d%dim
        ! compute gradient of st%X(psi)
        call X(f_gradient)(gr%sb, gr%f_der, &
          st%X(psi)(1:NP_PART, idim, ist, ik), grad(1:NP, idim, 1:NDIM))
      end do

      do i = 1, NDIM
        ! since the expectation value of the momentum operator is real
        ! for square integrable wfns this integral should be purely imaginary 
        ! for complex wfns but real for real wfns (see case distinction below)
        expect_val_p = X(states_dotp)(gr%m, st%d%dim, &
          st%X(psi)(1:NP, 1:st%d%dim, ist, ik), grad(1:NP, 1:st%d%dim, i))

        ! In the case of real wave functions we do not include the 
        ! -i prefactor of p = -i \nabla
        if (st%d%wfs_type == M_REAL) then
          st%momentum(i, ist, ik) = real( expect_val_p )
        else
          st%momentum(i, ist, ik) = real( -M_zI*expect_val_p )
        end if
      end do

      ! have to add the momentum vector in the case of periodic systems, 
      ! since st%X(psi) contains only u_k
      do i = 1, gr%sb%periodic_dim
        st%momentum(i, ist, ik) = st%momentum(i, ist, ik) + st%d%kpoints(i, ik)
      end do

    end do
  end do

  deallocate(grad)

  call pop_sub()
end subroutine X(states_calc_momentum)


! ---------------------------------------------------------
! It calculates the expectation value of the angular
! momentum of the state phi. If l2 is passed, it also
! calculates the expectation value of the square of the
! angular momentum of the state phi.
! ---------------------------------------------------------
subroutine X(states_angular_momentum)(gr, phi, l, l2)
  type(grid_t), intent(inout)  :: gr
  R_TYPE,       intent(inout)  :: phi(:, :)
  FLOAT,        intent(out)    :: l(MAX_DIM)
  FLOAT, optional, intent(out) :: l2

  integer :: idim, dim
  R_TYPE, allocatable :: lpsi(:, :)

  call push_sub('states_inc.Xstates_angular_momemtum')

  ASSERT(gr%m%sb%dim .ne.1)

  select case(gr%m%sb%dim)
  case(3)
    ALLOCATE(lpsi(NP_PART, 3), NP_PART*3)
  case(2)
    ALLOCATE(lpsi(NP_PART, 1), NP_PART*1)
  end select

  dim = size(phi, 2)

  l = M_ZERO
  if(present(l2)) l2 = M_ZERO

  do idim = 1, dim
#if defined(R_TREAL)
    l = M_ZERO
#else
    call X(f_angular_momentum)(gr%sb, gr%f_der, phi(:, idim), lpsi)
    select case(gr%m%sb%dim)
    case(3)
      l(1) = l(1) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
      l(2) = l(2) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 2))
      l(3) = l(3) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 3))
    case(2)
      l(3) = l(3) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
    end select
#endif
    if(present(l2)) then
      call X(f_l2)(gr%sb, gr%f_der, phi(:, idim), lpsi(:, 1))
      l2 = l2 + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
    end if
  end do

  deallocate(lpsi)
  call pop_sub()
end subroutine X(states_angular_momentum)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
