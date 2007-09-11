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
! res <- psi1(idx1)^+ * psi2(idx2).
subroutine X(states_blockt_mul)(mesh, st, psi1, psi2, res, idx1, idx2)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE, target,    intent(in)  :: psi1(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE, target,    intent(in)  :: psi2(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(out) :: res(:, :)
  integer, optional, intent(in)  :: idx1(:)
  integer, optional, intent(in)  :: idx2(:)

  integer              :: i, j, idim, mpi_err
  integer              :: m, n, max1, max2
  integer, pointer     :: ix1(:), ix2(:)
  integer, pointer     :: idx1_(:), idx2_(:)
  R_TYPE, pointer      :: bl1(:, :, :), bl2(:, :, :)
#if defined(HAVE_MPI)
  integer, allocatable :: idx1_l(:), idx2_l(:)
  integer, allocatable :: sendcnts(:), sdispls(:), recvcnts(:), rdispls(:)
#endif

  call profiling_in(C_PROFILING_LOBPCG_BLOCKT)
  call push_sub('states_inc.Xstates_blockt_mul')

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


  ! Quick and dirty parallelization. In the parallel case, we simply copy
  ! all states to all nodes.
  if(st%parallel_in_states) then
#if defined(HAVE_MPI)

! write(*, *) 'RANK', st%mpi_grp%rank, 'idx1_', idx1_
! write(*, *) 'RANK', st%mpi_grp%rank, 'idx2_', idx2_

    ! Exchange idx1_, idx2_.
!     ALLOCATE(sendcnts(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(sdispls(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(recvcnts(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(rdispls(st%mpi_grp%size), st%mpi_grp%size)

!     ALLOCATE(idx1_l(st%mpi_grp%size), st%mpi_grp%size)

!     ALLOCATE(idx2_l(st%mpi_grp%size), st%mpi_grp%size)
!     call MPI_Alltoall(m, 1, MPI_INTEGER, idx1_l, 1, MPI_INTEGER, st%mpi_grp%comm, mpi_err)
!     call MPI_Alltoall(n, 1, MPI_INTEGER, idx2_l, 1, MPI_INTEGER, st%mpi_grp%comm, mpi_err)
!     max1 = sum(idx1_l)
!     max2 = sum(idx2_l)
!     ALLOCATE(ix1(max1), max1)
!     ALLOCATE(ix2(max2), max2)

!     sendcnts   = m
!     sdispls    = 0
!     recvcnts   = idx1_l
!     rdispls(1) = 0
!     do i = 2, st%mpi_grp%size
!       rdispls(i) = rdispls(i-1)+recvcnts(i-1)
!     end do
!     call MPI_Alltoallv(idx1_, sendcnts, sdispls, MPI_INTEGER, &
!       ix1, recvcnts, rdispls, MPI_INTEGER, st%mpi_grp%comm, mpi_err)
! write(*, *) 'BAR'

!     sendcnts   = n
!     sdispls    = 0
!     recvcnts   = idx2_l
!     rdispls(1) = 0
!     do i = 2, st%mpi_grp%size
!       rdispls(i) = rdispls(i-1)+recvcnts(i-1)
!     end do
!     call MPI_Alltoallv(idx2_, sendcnts, sdispls, MPI_INTEGER, &
!       ix2, recvcnts, rdispls, MPI_INTEGER, st%mpi_grp%comm, mpi_err)

    ! Allocate space for the blocks.
    ALLOCATE(bl1(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)
    ALLOCATE(bl2(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)

    call X(states_gather)(mesh, st, psi1, bl1)
    call X(states_gather)(mesh, st, psi2, bl2)
#endif
  else
    bl1  => psi1
    bl2  => psi2
  end if

  ix1  => idx1_
  ix2  => idx2_
  max1 =  m
  max2 =  n

! if(st%mpi_grp%rank.eq.0) then
! do i =1, m
! write(*, *) 'BL1', ix1(i), bl1(1, 1, ix1(i))
! write(*, *) 'BL2', ix2(i), bl2(1, 1, ix2(i))
! end do
! end if
  
  ! FIXME: does not work with domain parallelization.
! if(max1.eq.max2.and.all(ix1.eq.ix2)) then
! write(*, *) 'IX1', ix1
! write(*, *) 'IX2'n, ix2
!   do i = 1, max1
!     res(i, i) = X(states_dotp)(mesh, st%d%dim, bl1(:, :, ix1(i)), bl2(:, :, ix2(i)))
!     do j = i+1, max2
!       !res(i, j) = R_TOTYPE(M_ZERO)
!       res(i, j) = X(states_dotp)(mesh, st%d%dim, bl1(:, :, ix1(i)), bl2(:, :, ix2(j)))
!       res(j, i) = R_CONJ(res(i, j))
! !       do idim = 1, st%d%dim
! !         res(i, j) = res(i, j) + &
! !           lalg_dot(mesh%np, bl1(:, idim, ix1(i)), bl2(:, idim, ix2(j)))
! !      end do
!     end do
!   end do
! else
  do i = 1, max1
    do j = 1, max2
      !res(i, j) = R_TOTYPE(M_ZERO)
      res(i, j) = X(states_dotp)(mesh, st%d%dim, bl1(:, :, ix1(i)), bl2(:, :, ix2(j)))
!       do idim = 1, st%d%dim
!         res(i, j) = res(i, j) + &
!           lalg_dot(mesh%np, bl1(:, idim, ix1(i)), bl2(:, idim, ix2(j)))
!      end do
    end do
  end do
!end if
  deallocate(idx1_, idx2_)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    deallocate(bl1, bl2)
!    deallocate(sendcnts, sdispls, recvcnts, rdispls)
  end if
#endif

  call pop_sub()
  call profiling_out(C_PROFILING_LOBPCG_BLOCKT)
end subroutine X(states_blockt_mul)


! ---------------------------------------------------------
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

  call MPI_Alltoallv(in(:, 1, 1), sendcnts, sdispls, R_MPITYPE, &
    out(:, 1, 1), recvcnts, rdispls, R_MPITYPE, st%mpi_grp%comm, mpi_err)

  deallocate(sendcnts, sdispls, recvcnts, rdispls)
#endif

  call pop_sub()
end subroutine X(states_gather)


! ---------------------------------------------------------
! Multiplication of block of states by matrix:
! res <- psi(idxp) * matr.
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
! Multiplication of block of states by matrix plus block of states:
! res <- alpha * psi(idx) * matr + beta * res.
subroutine X(states_block_matr_mul_add)(mesh, st, alpha, psi, matr, beta, res, idxp, idxr)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE,            intent(in)  :: alpha
  R_TYPE, target,    intent(in)  :: psi(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(in)  :: matr(:, :)
  R_TYPE,            intent(in)  :: beta
  R_TYPE, target,    intent(out) :: res(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  integer, optional, intent(in)  :: idxp(:), idxr(:)

  R_TYPE               :: tmp
  integer              :: l, m, n, i, j, k, idim, stat, maxp, maxr, mpi_err
  integer, pointer     :: idxp_(:), idxr_(:)
  integer, pointer     :: ixp(:), ixr(:)
  R_TYPE, pointer      :: blp(:, :, :), blr(:, :, :)
#if defined(HAVE_MPI)
  integer, allocatable :: idxp_l(:), idxr_l(:)
  integer, allocatable :: sendcnts(:), sdispls(:), recvcnts(:), rdispls(:)
#endif

  call profiling_in(C_PROFILING_LOBPCG_BLOCK_MATR)
  call push_sub('states_inc.Xstates_block_matr_add')

  if(present(idxp)) then
    m = ubound(idxp, 1)
    ALLOCATE(idxp_(m), m)
    idxp_ = idxp
  else
    m = st%nst
    ALLOCATE(idxp_(m), m)
    do i = 1, m
      idxp_(i) = i
    end do
  end if
  if(present(idxr)) then
    n = ubound(idxr, 1)
    ALLOCATE(idxr_(n), n)
    idxr_ = idxr
  else
    n = st%nst
    ALLOCATE(idxr_(n), n)
    do i = 1, n
      idxr_(i) = i
    end do
  end if

  l = ubound(matr, 2)

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Exchange idx1_, idx2_.
 !    ALLOCATE(sendcnts(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(sdispls(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(recvcnts(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(rdispls(st%mpi_grp%size), st%mpi_grp%size)

!     ALLOCATE(idxp_l(st%mpi_grp%size), st%mpi_grp%size)
!     ALLOCATE(idxr_l(st%mpi_grp%size), st%mpi_grp%size)
!     call MPI_Alltoall(m, 1, MPI_INTEGER, idxp_l, 1, MPI_INTEGER, st%mpi_grp%comm, mpi_err)
!     call MPI_Alltoall(n, 1, MPI_INTEGER, idxr_l, 1, MPI_INTEGER, st%mpi_grp%comm, mpi_err)
!     maxp = sum(idxp_l)
!     maxr = sum(idxr_l)
!     ALLOCATE(ixp(maxp), maxr)
!     ALLOCATE(ixr(maxp), maxr)

!     sendcnts   = m
!     sdispls    = 0
!     recvcnts   = idxp_l
!     rdispls(1) = 0
!     do i = 2, st%mpi_grp%size
!       rdispls(i) = rdispls(i-1)+recvcnts(i-1)
!     end do
!     call MPI_Alltoallv(idxp_, sendcnts, sdispls, MPI_INTEGER, &
!       ixp, recvcnts, rdispls, MPI_INTEGER, st%mpi_grp%comm, mpi_err)

!     sendcnts   = n
!     sdispls    = 0
!     recvcnts   = idxr_l
!     rdispls(1) = 0
!     do i = 1, st%mpi_grp%size
!       rdispls(i) = rdispls(i-1)+recvcnts(i-1)
!     end do
!     call MPI_Alltoallv(idxr_, sendcnts, sdispls, MPI_INTEGER, &
!       ixr, recvcnts, rdispls, MPI_INTEGER, st%mpi_grp%comm, mpi_err)

    
    ! Allocate space for the block.
    ALLOCATE(blp(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)
    ALLOCATE(blr(mesh%np_part, st%d%dim, st%nst), mesh%np_part*st%d%dim*st%nst)

    call X(states_gather)(mesh, st, psi, blp)
#endif
  else
    blp  => psi
    blr  => res
  end if
! if(st%mpi_grp%rank.eq.0) then
!   do i = 1, st%nst
!     write(*, *) i, blp(1, 1, i)
!   end do
!end if
  ixp  => idxp_
  ixr  => idxr_
  maxp =  m
  maxr =  ubound(matr, 2)
  
  ! FIXME: does not work with domain parallelization.  
  if(beta.eq.R_TOTYPE(M_ZERO)) then
    call profiling_in(C_PROFILING_LOBPCG_LOOP)
    do j = 1, maxr
      do idim = 1, st%d%dim
        do i = 1, mesh%np
          blr(i, idim, ixr(j)) = R_TOTYPE(M_ZERO)
          do k = 1, maxp
            blr(i, idim, ixr(j)) = blr(i, idim, ixr(j)) + blp(i, idim, ixp(k))*matr(k, j)
          end do
          blr(i, idim, ixr(j)) = alpha*blr(i, idim, ixr(j))
        end do
      end do
    end do
    call profiling_out(C_PROFILING_LOBPCG_LOOP)

  else
    call profiling_in(C_PROFILING_LOBPCG_LOOP)
    do j = 1, maxr
      do idim = 1, st%d%dim
        do i = 1, mesh%np
          tmp = R_TOTYPE(M_ZERO)
          do k = 1, maxp
!            blr(i, idim, ixr(j)) = blr(i, idim, ixr(j)) + blp(i, idim, ixp(k))*matr(k, j)
            tmp = tmp + blp(i, idim, ixp(k))*matr(k, j)
          end do
          blr(i, idim, ixr(j)) = alpha*tmp + beta*blr(i, idim, ixr(j))
        end do
      end do
    end do
    call profiling_out(C_PROFILING_LOBPCG_LOOP)
  end if
! if(st%mpi_grp%rank.eq.0) then
!   do i = 1, st%nst
!     write(*, *) 'RES', i, blr(1, 1, i)
!   end do
! end if

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    res(:, :, :) = blr(:, :, st%st_start:st%st_end)
    deallocate(blp, blr)
!    deallocate(sendcnts, sdispls, recvcnts, rdispls)
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
