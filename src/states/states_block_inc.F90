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
!! $Id: states_inc.F90 3262 2007-09-20 21:51:21Z lorenzen $


! ---------------------------------------------------------
!> Multiplication of two blocks of states:
!! res <- psi1(xpsi1)^+ * psi2(xpsi2) with the index sets xpsi1 and xpsi2.
subroutine X(states_blockt_mul)(mesh, st, psi1_start, psi2_start, &
  psi1, psi2, res, xpsi1, xpsi2, symm)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  integer,           intent(in)  :: psi1_start
  integer,           intent(in)  :: psi2_start
  R_TYPE, target,    intent(in)  :: psi1(:, :, psi1_start:)
  R_TYPE, target,    intent(in)  :: psi2(:, :, psi2_start:)
  R_TYPE,            intent(out) :: res(:, :)
  integer, optional, intent(in)  :: xpsi1(:)
  integer, optional, intent(in)  :: xpsi2(:)
  logical, optional, intent(in)  :: symm    !< Indicates if res(j, i) can be calculated as res(i, j)*.

  logical              :: symm_
  integer              :: ii
  integer              :: psi1_col, psi2_col
  integer, pointer     :: xpsi1_(:), xpsi2_(:)
  type(batch_t)        :: psi1b, psi2b
#if defined(HAVE_MPI)
  integer              :: jj
  integer              :: size, rank, round, res_col_offset, res_row_offset
  integer              :: dst, src, kk, ll, recvcnt, sendcnt, left, right, max_count
  integer              :: stats(MPI_STATUS_SIZE, 2), reqs(2)
  integer, pointer     :: xpsi1_count(:), xpsi2_count(:), xpsi1_node(:, :), xpsi2_node(:, :)
  R_TYPE, pointer      :: sendbuf(:, :, :), recvbuf(:, :, :), tmp_ptr(:, :, :)
  R_TYPE, allocatable  :: res_tmp(:, :)
  R_TYPE, allocatable  :: psi1_block(:, :, :), res_local(:, :)
#endif

  call profiling_in(C_PROFILING_BLOCKT)
  PUSH_SUB(X(states_blockt_mul))

  symm_ = .false.
  if(present(symm)) then
    symm_ = symm
  end if

  ! Calculate index sets of state block psi1 and psi2.
  if(present(xpsi1)) then
    call make_idx_set(st%nst, xpsi1_, psi1_col, xpsi1)
  else
    call make_idx_set(st%nst, xpsi1_, psi1_col)
  end if
  if(present(xpsi2)) then
    call make_idx_set(st%nst, xpsi2_, psi2_col, xpsi2)
  else
    call make_idx_set(st%nst, xpsi2_, psi2_col)
  end if

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Shortcuts.
    size = st%mpi_grp%size
    rank = st%mpi_grp%rank

    ! Calculate the index sets per node,
    ! xpsi1_node(1:xpsi1_count(node), node) and
    ! xpsi2_node(1:xpsi2_count(node), node) are the index sets.
    call states_block_local_idx(st, xpsi1_, psi1_col, xpsi1_count, xpsi1_node)
    call states_block_local_idx(st, xpsi2_, psi2_col, xpsi2_count, xpsi2_node)

    ! Has to be zero because we use an allreduce on this at the end.
    res = R_TOTYPE(M_ZERO)

    ! Compact psi1 in order to use BLAS gemm on it.
    if(.not.mesh%use_curvilinear) then
      SAFE_ALLOCATE(psi1_block(1:mesh%np, 1:st%d%dim, 1:xpsi1_count(rank)))
      call profiling_in(C_PROFILING_BLOCKT_CP)
      call X(states_compactify)(st%d%dim, mesh, psi1_start, &
        xpsi1_node(1:xpsi1_count(rank), rank), psi1, psi1_block)
      call profiling_out(C_PROFILING_BLOCKT_CP)
    end if

    ! Allocate send and receive buffers. For some blocks, they are oversized but
    ! buffer reusing saves some overhead in the MPI library.
    max_count = maxval(xpsi2_count)
    SAFE_ALLOCATE(sendbuf(1:mesh%np, 1:st%d%dim, 1:max_count))
    SAFE_ALLOCATE(recvbuf(1:mesh%np, 1:st%d%dim, 1:max_count))

    ! Compact the local block to send away.
    sendcnt = xpsi2_count(rank)
    call X(states_compactify)(st%d%dim, mesh, psi2_start, xpsi2_node(1:xpsi2_count(rank), rank), psi2, sendbuf)

    ! Get neighbours.
    call MPI_Cart_shift(st%dom_st_mpi_grp%comm, P_STRATEGY_STATES-1, -1, src, dst, mpi_err)
    right = lmpi_translate_rank(st%dom_st_mpi_grp%comm, st%mpi_grp%comm, src)
    left  = lmpi_translate_rank(st%dom_st_mpi_grp%comm, st%mpi_grp%comm, dst)

    do round = 0, size-1
      kk = mod(rank+round, size)   ! The column of the block currently being calculated.
      ll = mod(rank+round+1, size) ! The column of the block currently being communicated.
      ! In all but the first rounds we have to wait for the data to arrive and
      ! then swap buffers.
      if(round.gt.0) then
        call MPI_Waitall(2, reqs, stats, mpi_err)
        tmp_ptr => sendbuf
        sendbuf => recvbuf
        recvbuf => tmp_ptr
        sendcnt =  recvcnt
      end if
      ! In all but the last rounds the data has to be passed on to the left neighbour, and
      ! accordingly received from the right neighbour.
      if(round.lt.size-1) then
        recvcnt = xpsi2_count(ll)
        call MPI_Irecv(recvbuf(1, 1, 1), mesh%np*st%d%dim*recvcnt, R_MPITYPE, right, 0, &
          st%mpi_grp%comm, reqs(1), mpi_err)
        call MPI_Isend(sendbuf(1, 1, 1), mesh%np*st%d%dim*sendcnt, R_MPITYPE, left, 0,  & 
          st%mpi_grp%comm, reqs(2), mpi_err)
      end if
      ! Do the matrix multiplication.
      res_row_offset = sum(xpsi1_count(0:rank-1))
      res_col_offset = sum(xpsi2_count(0:kk-1))
      if(.not.mesh%use_curvilinear) then
        if(xpsi1_count(rank).gt.0.and.sendcnt.gt.0) then
          SAFE_ALLOCATE(res_local(1:xpsi1_count(rank), 1:sendcnt))

          call profiling_in(C_PROFILING_BLOCKT_MM)
          call lalg_gemmt(xpsi1_count(rank), sendcnt, mesh%np*st%d%dim, R_TOTYPE(mesh%vol_pp(1)), &
            psi1_block, sendbuf, R_TOTYPE(M_ZERO), res_local)
          call profiling_out(C_PROFILING_BLOCKT_MM)

          call profiling_in(C_PROFILING_BLOCKT_CP)
          res(res_row_offset+1:res_row_offset+xpsi1_count(rank), res_col_offset+1:res_col_offset+sendcnt) = res_local
          call profiling_out(C_PROFILING_BLOCKT_CP)
          SAFE_DEALLOCATE_A(res_local)
        end if
      else ! Curvilinear coordinates.
        do ii = 1, xpsi1_count(rank)
          do jj = 1, sendcnt
            res(ii + res_row_offset, jj + res_col_offset) = &
              X(mf_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(ii, rank)), sendbuf(:, :, jj), reduce=.false.)
          end do
        end do
      end if
    end do

    SAFE_DEALLOCATE_P(sendbuf)
    SAFE_DEALLOCATE_P(recvbuf)
    ! Add up all the individual blocks.
    call profiling_in(C_PROFILING_BLOCKT_AR)
#ifndef HAVE_MPI2
    SAFE_ALLOCATE(res_tmp(1:psi1_col, 1:psi2_col))
    res_tmp = res
#endif
    call MPI_Allreduce(MPI_IN_PLACE_OR(res_tmp), res, psi1_col*psi2_col, R_MPITYPE, MPI_SUM, st%dom_st_mpi_grp%comm, mpi_err)
    call profiling_out(C_PROFILING_BLOCKT_AR)
    SAFE_DEALLOCATE_A(res_tmp)
    SAFE_DEALLOCATE_P(xpsi1_count)
    SAFE_DEALLOCATE_P(xpsi2_count)
    SAFE_DEALLOCATE_P(xpsi1_node)
    SAFE_DEALLOCATE_P(xpsi2_node)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call messages_fatal(1)
#endif
  else ! No states parallelization.

    if(present(xpsi1)) then
      call batch_init(psi1b, st%d%dim, psi1_col)
      do ii = 1, psi1_col
        call batch_add_state(psi1b, ii, psi1(:, :, xpsi1(ii)))
      end do
    else
      call batch_init(psi1b, st%d%dim, 1, psi1_col, psi1(:, :, :))
    end if

    if(present(xpsi2)) then
      call batch_init(psi2b, st%d%dim, psi2_col)
      do ii = 1, psi2_col
        call batch_add_state(psi2b, ii, psi2(:, :, xpsi2(ii)))
      end do
    else
      call batch_init(psi2b, st%d%dim, 1, psi2_col, psi2(:, :, :))
    end if

    ASSERT(batch_is_ok(psi1b))
    ASSERT(batch_is_ok(psi2b))

    call X(mesh_batch_dotp_matrix)(mesh, psi1b, psi2b, res, symm = symm_)
    
    call batch_end(psi1b)
    call batch_end(psi2b)

  end if
  SAFE_DEALLOCATE_P(xpsi1_)
  SAFE_DEALLOCATE_P(xpsi2_)

  POP_SUB(X(states_blockt_mul))
  call profiling_out(C_PROFILING_BLOCKT)
end subroutine X(states_blockt_mul)


! ---------------------------------------------------------
!> Gather all states on all nodes. out has to be of sufficient size.
#if defined(HAVE_MPI)
subroutine X(states_gather)(mesh, st, in, out)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: mesh
  R_TYPE,         intent(in)  :: in(:, :, :)
  R_TYPE,         intent(out) :: out(:, :, :)

  integer              :: ii
  integer, allocatable :: sendcnts(:), sdispls(:), recvcnts(:), rdispls(:)

  PUSH_SUB(X(states_gather))

  SAFE_ALLOCATE(sendcnts(1:st%mpi_grp%size))
  SAFE_ALLOCATE( sdispls(1:st%mpi_grp%size))
  SAFE_ALLOCATE(recvcnts(1:st%mpi_grp%size))
  SAFE_ALLOCATE( rdispls(1:st%mpi_grp%size))

  sendcnts   = mesh%np_part*st%d%dim*st%st_num(st%mpi_grp%rank)
  sdispls    = 0
  recvcnts   = st%st_num*mesh%np_part*st%d%dim
  rdispls(1) = 0
  do ii = 2, st%mpi_grp%size
    rdispls(ii) = rdispls(ii-1) + recvcnts(ii-1)
  end do

  call mpi_debug_in(st%mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(in(:, 1, 1), sendcnts, sdispls, R_MPITYPE, &
    out(:, 1, 1), recvcnts, rdispls, R_MPITYPE, st%mpi_grp%comm, mpi_err)
  call mpi_debug_out(st%mpi_grp%comm, C_MPI_ALLTOALLV)

  SAFE_DEALLOCATE_A(sendcnts)
  SAFE_DEALLOCATE_A(sdispls)
  SAFE_DEALLOCATE_A(recvcnts)
  SAFE_DEALLOCATE_A(rdispls)

  POP_SUB(X(states_gather))
end subroutine X(states_gather)
#endif


! ---------------------------------------------------------
!> Multiplication of block of states with indices idxp by matrix and
!! update columns with idxr in the result.
!! res(xres) <- psi(xpsi) * matr.
subroutine X(states_block_matr_mul)(mesh, st, psi_start, psi_end, res_start, res_end, &
  psi, matr, res, xpsi, xres)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  integer,           intent(in)  :: psi_start
  integer,           intent(in)  :: psi_end
  integer,           intent(in)  :: res_start
  integer,           intent(in)  :: res_end
  R_TYPE,            intent(in)  :: psi(:, :, psi_start:) !< (mesh%np_part, st%d%dim, psi_start:psi_end)
  R_TYPE,            intent(in)  :: matr(:, :)
  R_TYPE,            intent(out) :: res(:, :, :)
  integer, optional, intent(in)  :: xpsi(:), xres(:)

  PUSH_SUB(X(states_block_matr_mul))

  if(present(xpsi).and.present(xres)) then
    call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi_start, psi_end, &
      res_start, res_end, psi, matr, R_TOTYPE(M_ZERO), res, xpsi, xres)
  else if(present(xpsi)) then
    call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi_start, psi_end, &
      res_start, res_end, psi, matr, R_TOTYPE(M_ZERO), res, xpsi=xpsi)
  else if(present(xres)) then
    call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi_start, psi_end, &
      res_start, res_end, psi, matr, R_TOTYPE(M_ZERO), res, xres=xres)
  else
    call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi_start, psi_end, &
      res_start, res_end, psi, matr, R_TOTYPE(M_ZERO), res)
  end if

  POP_SUB(X(states_block_matr_mul))
end subroutine X(states_block_matr_mul)


! ---------------------------------------------------------
!> Multiplication of block of states by matrix plus block of states
!! (with the corresponding column indices):
!! res(xres) <- alpha * psi(xpsi) * matr + beta * res(xres).
subroutine X(states_block_matr_mul_add)(mesh, st, alpha, psi_start, psi_end, res_start, res_end, &
  psi, matr, beta, res, xpsi, xres)
  type(mesh_t),      intent(in)    :: mesh
  type(states_t),    intent(in)    :: st
  R_TYPE,            intent(in)    :: alpha
  integer,           intent(in)    :: psi_start
  integer,           intent(in)    :: psi_end
  integer,           intent(in)    :: res_start
  integer,           intent(in)    :: res_end
  R_TYPE,            intent(in)    :: psi(:, :, psi_start:) !< (mesh%np_part, st%d%dim, psi_start:psi_end)
  R_TYPE,            intent(in)    :: matr(:, :)
  R_TYPE,            intent(in)    :: beta
  R_TYPE,            intent(inout) :: res(:, :, res_start:) !< (mesh%np_part, st%d%dim, res_start:res_end)
  integer, optional, intent(in)    :: xpsi(:), xres(:)

  integer              :: res_col, psi_col, matr_col
  integer, pointer     :: xpsi_(:), xres_(:)
  R_TYPE, allocatable  :: res_block(:, :, :), matr_block(:, :), psi_block(:, :, :)
#if defined(HAVE_MPI)
  integer              :: rank, size, round, matr_row_offset, matr_col_offset, ii
  integer              :: src, dst, left, right, kk, ll, idim, sendcnt, recvcnt, max_count
  integer              :: stats(MPI_STATUS_SIZE, 2), reqs(2)
  integer, pointer     :: xpsi_count(:), xres_count(:), xpsi_node(:, :), xres_node(:, :)
  R_TYPE, pointer      :: sendbuf(:, :, :), recvbuf(:, :, :), tmp_ptr(:, :, :)
#endif

  call profiling_in(C_PROFILING_BLOCK_MATR)
  PUSH_SUB(X(states_block_matr_mul_add))

  ! Calculate global index sets of state block psi and res.
  if(present(xpsi)) then
    call make_idx_set(st%nst, xpsi_, psi_col, xpsi)
  else
    call make_idx_set(st%nst, xpsi_, psi_col)
  end if
  if(present(xres)) then
    call make_idx_set(st%nst, xres_, res_col, xres)
  else
    call make_idx_set(st%nst, xres_, res_col)
  end if

  matr_col = ubound(matr, 2)

  ! There is a little code duplication between the serial and parallel case
  ! but the code is easier to understand having it separated (instead a lot of
  ! conditionals and pointers).
  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Shortcuts.
    size = st%mpi_grp%size
    rank = st%mpi_grp%rank

    ! Calculate the index sets per node.
    ! xpsi_node(1:xpsi_count(node), node) and
    ! xres_node(1:xres_count(node), node) are the index sets.
    call states_block_local_idx(st, xpsi_, psi_col, xpsi_count, xpsi_node)
    call states_block_local_idx(st, xres_, res_col, xres_count, xres_node)

    ! Take care of beta first, if necessary, and compact res to res_block.
    if(beta.ne.R_TOTYPE(M_ZERO)) then
      call profiling_in(C_PROFILING_BLOCK_MATR_CP)
      do ii = 1, xres_count(rank)
        do idim = 1, st%d%dim
          call lalg_scal(mesh%np, beta, res(:, idim, xres_node(ii, rank)))
        end do
      end do
      call profiling_out(C_PROFILING_BLOCK_MATR_CP)
    else
      res = R_TOTYPE(M_ZERO)
    end if
    SAFE_ALLOCATE(res_block(1:mesh%np, 1:st%d%dim, 1:xres_count(rank)))
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st%d%dim, mesh, res_start, xres_node(1:xres_count(rank), rank), res, res_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! Allocate send and receive buffers. For some blocks, they are oversized but
    ! buffer reusing saves some overhead in the MPI library.
    max_count = maxval(xpsi_count)
    SAFE_ALLOCATE(sendbuf(1:mesh%np, 1:st%d%dim, 1:max_count))
    SAFE_ALLOCATE(recvbuf(1:mesh%np, 1:st%d%dim, 1:max_count))

    ! Compact the local block to send away.
    sendcnt = xpsi_count(rank)
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st%d%dim, mesh, psi_start, xpsi_node(1:xpsi_count(rank), rank), psi(:, :, :), sendbuf)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! Get neighbours.
    call MPI_Cart_shift(st%dom_st_mpi_grp%comm, P_STRATEGY_STATES-1, -1, src, dst, mpi_err)
    right = lmpi_translate_rank(st%dom_st_mpi_grp%comm, st%mpi_grp%comm, src)
    left  = lmpi_translate_rank(st%dom_st_mpi_grp%comm, st%mpi_grp%comm, dst)

    ! Asynchronously left-rotate blocks of psi.
    do round = 0, size-1
      kk = mod(rank+round, size)   ! The column of the block currently being calculated.
      ll = mod(rank+round+1, size) ! The column of the block currently being communicated.
      ! In all but the first rounds we have to wait for the data to arrive and
      ! then swap buffers, i.e., what we received in the send buffer in the remainder
      ! of this loop, a bit confusing.
      if(round.gt.0) then
        call MPI_Waitall(2, reqs, stats, mpi_err)
        tmp_ptr => sendbuf
        sendbuf => recvbuf
        recvbuf => tmp_ptr
        sendcnt =  recvcnt
      end if
      ! In all but the last rounds the data has to be passed on to the left neighbour, and
      ! accordingly received from the right neighbour.
      if(round.lt.size-1) then
        recvcnt = xpsi_count(ll)
        call MPI_Irecv(recvbuf, mesh%np*st%d%dim*recvcnt, R_MPITYPE, right, 0, st%mpi_grp%comm, reqs(1), mpi_err)
        call MPI_Isend(sendbuf, mesh%np*st%d%dim*sendcnt, R_MPITYPE, left, 0, st%mpi_grp%comm, reqs(2), mpi_err)
      end if

      ! Do the matrix multiplication.
      matr_row_offset = sum(xpsi_count(0:kk-1))
      matr_col_offset = sum(xres_count(0:rank-1))
      SAFE_ALLOCATE(matr_block(1:xpsi_count(kk), 1:xres_count(rank)))
      if(sendcnt.gt.0.and.xres_count(rank).gt.0) then
        call profiling_in(C_PROFILING_BLOCK_MATR_CP)
        matr_block = matr(matr_row_offset+1:matr_row_offset+xpsi_count(kk), &
          matr_col_offset+1:matr_col_offset+xres_count(rank))
        call profiling_out(C_PROFILING_BLOCK_MATR_CP)
        call profiling_in(C_PROFILING_BLOCK_MATR_MM)
        call lalg_gemm(mesh%np * st%d%dim, xres_count(rank), sendcnt, alpha, &
          sendbuf, matr_block, R_TOTYPE(M_ONE), res_block)
        call profiling_out(C_PROFILING_BLOCK_MATR_MM)
      end if
      SAFE_DEALLOCATE_A(matr_block)
    end do

    ! Copy result.
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_uncompactify)(st%d%dim, mesh, res_start, xres_node(1:xres_count(rank), rank), res_block, res)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    SAFE_DEALLOCATE_A(res_block)
    SAFE_DEALLOCATE_P(sendbuf)
    SAFE_DEALLOCATE_P(recvbuf)
    SAFE_DEALLOCATE_P(xpsi_count)
    SAFE_DEALLOCATE_P(xres_count)
    SAFE_DEALLOCATE_P(xpsi_node)
    SAFE_DEALLOCATE_P(xres_node)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call messages_fatal(1)
#endif
  else ! No states parallelization.
    ! Compact everything to pass it to BLAS.
    SAFE_ALLOCATE(res_block(1:mesh%np, 1:st%d%dim, 1:res_col))
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st%d%dim, mesh, res_start, xres_(1:res_col), res, res_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    SAFE_ALLOCATE(psi_block(1:mesh%np, 1:st%d%dim, 1:psi_col))
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st%d%dim, mesh, psi_start, xpsi_(1:psi_col), psi, psi_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! matr_block is needed because matr may be an assumed-shape array.
    SAFE_ALLOCATE(matr_block(1:psi_col, 1:matr_col))
    matr_block = matr
    call lalg_gemm(mesh%np * st%d%dim, matr_col, psi_col, alpha, &
      psi_block, matr_block, beta, res_block)
    SAFE_DEALLOCATE_A(matr_block)

    ! Copy result.
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_uncompactify)(st%d%dim, mesh, res_start, xres_(1:res_col), res_block, res)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)
    SAFE_DEALLOCATE_A(res_block)
    SAFE_DEALLOCATE_A(psi_block)
  end if

  SAFE_DEALLOCATE_P(xpsi_)
  SAFE_DEALLOCATE_P(xres_)

  POP_SUB(X(states_block_matr_mul_add))
  call profiling_out(C_PROFILING_BLOCK_MATR)
end subroutine X(states_block_matr_mul_add)


! ---------------------------------------------------------
!> Copy in(:, :, idx) to out out(:, :, 1:ubound(idx, 1)).
subroutine X(states_compactify)(dim, mesh, in_start, idx, in, out)
  integer,      intent(in)  :: dim
  type(mesh_t), intent(in)  :: mesh
  integer,      intent(in)  :: in_start
  integer,      intent(in)  :: idx(:)
  R_TYPE,       intent(in)  :: in(:, :, :)
  R_TYPE,       intent(out) :: out(:, :, :)

  integer :: ist, idim, nn
  type(profile_t), save :: prof

  call profiling_in(prof, "STATES_COMPACTIFY")
  PUSH_SUB(X(states_compactify))

  nn = ubound(idx, 1)

  do ist = 1, nn
    do idim = 1, dim
      call lalg_copy(mesh%np, in(:, idim, idx(ist) - in_start + 1), out(:, idim, ist))
    end do
  end do
  
  POP_SUB(X(states_compactify))
  call profiling_out(prof)
end subroutine X(states_compactify)


! ---------------------------------------------------------
!> Undo the effect of X(states_compactify), i.e.
!! X(states_compactify)(st, mesh, idx, in, out) followed by
!! X(states_uncompactify)(st, mesh, idx, out, in) is identity.
subroutine X(states_uncompactify)(dim, mesh, out_start, idx, in, out)
  integer,      intent(in)  :: dim
  type(mesh_t), intent(in)  :: mesh
  integer,      intent(in)  :: out_start
  integer,      intent(in)  :: idx(:)
  R_TYPE,       intent(in)  :: in(:, :, :)
  R_TYPE,       intent(out) :: out(:, :, :)

  integer :: ist, idim, nn
  type(profile_t), save :: prof

  call profiling_in(prof, "STATES_UNCOMPACTIFY")
  PUSH_SUB(X(states_uncompactify))

  nn = ubound(idx, 1)

  do ist = 1, nn
    do idim = 1, dim
      call lalg_copy(mesh%np, in(:, idim, ist), out(:, idim, idx(ist) - out_start + 1))
    end do
  end do
  
  POP_SUB(X(states_uncompactify))
  call profiling_out(prof)
end subroutine X(states_uncompactify)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
