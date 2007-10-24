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
! Multiplication of two blocks of states:
! res <- psi1(xpsi1)^+ * psi2(xpsi2) with the index sets xpsi1 and xpsi2.
subroutine X(states_blockt_mul)(mesh, st, psi1, psi2, res, xpsi1, xpsi2, symm)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE, target,    intent(in)  :: psi1(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE, target,    intent(in)  :: psi2(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(out) :: res(:, :)
  integer, optional, intent(in)  :: xpsi1(:)
  integer, optional, intent(in)  :: xpsi2(:)
  logical, optional, intent(in)  :: symm    ! Indicates if res(j, i) can be calculated as
                                            ! res(i, j)*.

  logical              :: symm_
  integer              :: i, j
  integer              :: psi1_col, psi2_col
  integer, allocatable :: xpsi1_(:), xpsi2_(:)
  R_TYPE, allocatable  :: psi1_block(:, :, :), psi2_block(:, :, :), res_local(:, :)
#if defined(HAVE_MPI)
  integer              :: size, rank, node, round, ist, res_col_offset, res_row_offset
  integer, allocatable :: xpsi1_count(:), xpsi2_count(:), xpsi1_node(:, :), xpsi2_node(:, :)
  integer              :: dst, src, k, l, recvcnt, sendcnt, left, right, max_count
  integer              :: stats(MPI_STATUS_SIZE, 2), reqs(2)
  R_TYPE, pointer      :: sendbuf(:, :, :), recvbuf(:, :, :), tmp_ptr(:, :, :)
  R_TYPE, allocatable  :: res_tmp(:, :)
#endif

  call profiling_in(C_PROFILING_BLOCKT)
  call push_sub('states_inc.Xstates_blockt_mul')

  symm_ = .false.
  if(present(symm)) then
    symm_ = symm
  end if

  ! Calculate index sets.
  ! FIXME: this is the same as in Xstates_block_matr_mul_add and should be
  ! put into an extra routine.
  if(present(xpsi1)) then
    psi1_col = ubound(xpsi1, 1)
    ALLOCATE(xpsi1_(psi1_col), psi1_col)
    xpsi1_ = xpsi1
  else
    psi1_col = st%nst
    ALLOCATE(xpsi1_(psi1_col), psi1_col)
    do i = 1, psi1_col
      xpsi1_(i) = i
    end do
  end if
  if(present(xpsi2)) then
    psi2_col = ubound(xpsi2, 1)
    ALLOCATE(xpsi2_(psi2_col), psi2_col)
    xpsi2_ = xpsi2
  else
    psi2_col = st%nst
    ALLOCATE(xpsi2_(psi2_col), psi2_col)
    do i = 1, psi2_col
      xpsi2_(i) = i
    end do
  end if

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Shortcuts.
    size = st%mpi_grp%size
    rank = st%mpi_grp%rank

    ! Calculate the index sets per node,
    ! xpsi1_node(1:xpsi1_count(node), node) and
    ! xpsi2_node(1:xpsi2_count(node), node) are the index sets.
    ! FIXME: this is the same as in Xstates_block_matr_mul_add and should be
    ! put into an extra routine.
    ALLOCATE(xpsi1_count(0:size-1), size)
    ALLOCATE(xpsi2_count(0:size-1), size)
    ! Count the how many vectors each node has.
    xpsi1_count = 0
    xpsi2_count = 0
    do i = 1, psi1_col
      xpsi1_count(st%node(xpsi1_(i))) = xpsi1_count(st%node(xpsi1_(i))) + 1
    end do
    do i = 1, psi2_col
      xpsi2_count(st%node(xpsi2_(i))) = xpsi2_count(st%node(xpsi2_(i))) + 1
    end do
    ! Allocate space, it is a bit more than really required.
    ALLOCATE(xpsi1_node(maxval(xpsi1_count), 0:size-1), maxval(xpsi1_count)*size)
    ALLOCATE(xpsi2_node(maxval(xpsi2_count), 0:size-1), maxval(xpsi2_count)*size)
    ! Now set up the index sets.
    xpsi1_count = 0
    xpsi2_count = 0
    xpsi1_node  = 0
    xpsi2_node  = 0
    do ist = 1, st%nst
      node = st%node(ist)
      ! A state ist is only included if its in the global index set.
      if(member(ist, xpsi1_)) then
        xpsi1_count(node)                   = xpsi1_count(node) + 1
        xpsi1_node(xpsi1_count(node), node) = ist
      end if
      if(member(ist, xpsi2_)) then
        xpsi2_count(node)                   = xpsi2_count(node) + 1
        xpsi2_node(xpsi2_count(node), node) = ist
      end if
    end do

    ! Has to be zero because we use an allreduce on this at the end.
    res = R_TOTYPE(M_ZERO)

    ! Compact psi1 in order to use BLAS gemm on it.
    if(.not.mesh%use_curvlinear) then
      ALLOCATE(psi1_block(mesh%np, st%d%dim, xpsi1_count(rank)), mesh%np*st%d%dim*xpsi1_count(rank))
      call profiling_in(C_PROFILING_BLOCKT_CP)
      call X(states_compactify)(st, mesh, xpsi1_node(1:xpsi1_count(rank), rank), psi1, psi1_block)
      call profiling_out(C_PROFILING_BLOCKT_CP)
    end if

    ! Allocate send and receive buffers. For some blocks, they are oversized but
    ! buffer reusing saves some overhead in the MPI library.
    max_count = maxval(xpsi2_count)
    ALLOCATE(sendbuf(mesh%np, st%d%dim, max_count), mesh%np*st%d%dim*max_count)
    ALLOCATE(recvbuf(mesh%np, st%d%dim, max_count), mesh%np*st%d%dim*max_count)

    ! Compact the local block to send away.
    sendcnt = xpsi2_count(rank)
    call X(states_compactify)(st, mesh, xpsi2_node(1:xpsi2_count(rank), rank), psi2, sendbuf)

    ! Get neighbours.
    call MPI_Cart_shift(st%dom_st%comm, P_STRATEGY_STATES-1, -1, src, dst, mpi_err)
    right = lmpi_translate_rank(st%dom_st%comm, st%mpi_grp%comm, src)
    left  = lmpi_translate_rank(st%dom_st%comm, st%mpi_grp%comm, dst)

    do round = 0, size-1
      k = mod(rank+round, size)   ! The column of the block currently being calculated.
      l = mod(rank+round+1, size) ! The column of the block currently being communicated.
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
        recvcnt = xpsi2_count(l)
        call MPI_Irecv(recvbuf, mesh%np*st%d%dim*recvcnt, R_MPITYPE, right, 0, st%mpi_grp%comm, reqs(1), mpi_err)
        call MPI_Isend(sendbuf, mesh%np*st%d%dim*sendcnt, R_MPITYPE, left, 0,  st%mpi_grp%comm, reqs(2), mpi_err)
      end if
      
      ! Do the matrix multiplication.
      res_row_offset = sum(xpsi1_count(0:rank-1))
      res_col_offset = sum(xpsi2_count(0:k-1))
      if(.not.mesh%use_curvlinear) then
        if(xpsi1_count(rank).gt.0.and.sendcnt.gt.0) then
          ALLOCATE(res_local(xpsi1_count(rank), sendcnt), xpsi1_count(rank)*sendcnt)

          call profiling_in(C_PROFILING_BLOCKT_MM)
          call lalg_gemmt(xpsi1_count(rank), sendcnt, mesh%np*st%d%dim, R_TOTYPE(mesh%vol_pp(1)), &
            psi1_block(:, :, 1), sendbuf(:, :, 1), R_TOTYPE(M_ZERO), res_local)
          call profiling_out(C_PROFILING_BLOCKT_MM)

          call profiling_in(C_PROFILING_BLOCKT_CP)
          res(res_row_offset+1:res_row_offset+xpsi1_count(rank), res_col_offset+1:res_col_offset+sendcnt) = res_local
          call profiling_out(C_PROFILING_BLOCKT_CP)
          deallocate(res_local)
        end if
      else ! Curvilinear coordinates.
        do i = 1, xpsi1_count(rank)
          do j = 1, sendcnt
            res(i+res_row_offset, j+res_col_offset) = &
              X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(i, rank)), sendbuf(:, :, j), reduce=.false.)
          end do
        end do
      end if
    end do

    deallocate(sendbuf, recvbuf)

    ! Add up all the individual blocks.
    ALLOCATE(res_tmp(psi1_col, psi2_col), psi1_col*psi2_col)
    call profiling_in(C_PROFILING_BLOCKT_AR)
    call MPI_Allreduce(res, res_tmp, psi1_col*psi2_col, R_MPITYPE, MPI_SUM, st%dom_st%comm, mpi_err)
    call profiling_out(C_PROFILING_BLOCKT_AR)
    call profiling_in(C_PROFILING_BLOCKT_CP)
    res = res_tmp
    call profiling_out(C_PROFILING_BLOCKT_CP)
    deallocate(res_tmp)
    deallocate(xpsi1_count, xpsi2_count, xpsi1_node, xpsi2_node, xpsi1_, xpsi2_)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call write_fatal(1)
#endif
  else ! No states parallelization.
    if(.not.mesh%use_curvlinear) then
      ASSERT(psi1_col.ne.0.and.psi2_col.ne.0)
      ALLOCATE(res_local(psi1_col, psi2_col), psi1_col*psi2_col)
      ALLOCATE(psi1_block(mesh%np, st%d%dim, psi1_col), mesh%np*st%d%dim*psi1_col)
      ALLOCATE(psi2_block(mesh%np, st%d%dim, psi2_col), mesh%np*st%d%dim*psi2_col)
      
      call X(states_compactify)(st, mesh, xpsi1_(1:psi1_col), psi1, psi1_block)
      call X(states_compactify)(st, mesh, xpsi2_(1:psi2_col), psi2, psi2_block)
      call lalg_gemmt(psi1_col, psi2_col, mesh%np*st%d%dim, R_TOTYPE(mesh%vol_pp(1)), &
        psi1_block(:, :, 1), psi2_block(:, :, 1), R_TOTYPE(M_ZERO), res_local)
      ! We have to use the intermediate res_local and cannot put res directly into
      ! lalg_gemmt because otherwise our BLAS interface screws up everything because
      ! it puts 1s in all free indices! Cost me several headaches to figure this out.
      ! Perhaps, one day we should write a more comfortable interface to BLAS, or at
      ! least one that is clearer!
      res = res_local
      deallocate(psi1_block, psi2_block, res_local)
    else ! Curvilinear coordinates.
      if(symm_) then
        do i = 1, psi1_col
          res(i, i) = X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_(i)), psi2(:, :, xpsi2_(i)), reduce=.false.)
          do j = i+1, psi2_col
            res(i, j) = X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_(i)), psi2(:, :, xpsi2_(j)), reduce=.false.)
            res(j, i) = R_CONJ(res(i, j))
          end do
        end do
      else
        do i = 1, psi1_col
          do j = 1, psi2_col
            res(i, j) = X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_(i)), psi2(:, :, xpsi2_(j)), reduce=.false.)
          end do
        end do
      end if
    end if
    ! This has to be done because no reduction is done in Xstates_dotp above.
    ! Rationale: this way, all the reduce operations are done at once (hopefully)
    ! saving latency overhead.
    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ALLOCATE(res_tmp(psi1_col, psi2_col), psi1_col*psi2_col)
      call MPI_Allreduce(res, res_tmp, psi1_col*psi2_col, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
      res = res_tmp
      deallocate(res_tmp)
#else
      message(1) = "Running parallel in domain without MPI. This is a bug!"
      call write_fatal(1)
#endif
    end if
  end if
  
  call pop_sub()
  call profiling_out(C_PROFILING_BLOCKT)
end subroutine X(states_blockt_mul)


! ---------------------------------------------------------
! Gather all states on all nodes. out has to be of sufficient size.
#if defined(HAVE_MPI)
subroutine X(states_gather)(mesh, st, in, out)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: mesh
  R_TYPE,         intent(in)  :: in(:, :, :)
  R_TYPE,         intent(out) :: out(:, :, :)

  integer              :: i
  integer, allocatable :: sendcnts(:), sdispls(:), recvcnts(:), rdispls(:)

  call push_sub('states_inc.Xstates_gather')

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

  call mpi_debug_in(st%mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(in(:, 1, 1), sendcnts, sdispls, R_MPITYPE, &
    out(:, 1, 1), recvcnts, rdispls, R_MPITYPE, st%mpi_grp%comm, mpi_err)
  call mpi_debug_out(st%mpi_grp%comm, C_MPI_ALLTOALLV)

  deallocate(sendcnts, sdispls, recvcnts, rdispls)

  call pop_sub()
end subroutine X(states_gather)
#endif


! ---------------------------------------------------------
! Multiplication of block of states with indices idxp by matrix and
! update columns with idxr in the result.
! res(xres) <- psi(xpsi) * matr.
subroutine X(states_block_matr_mul)(mesh, st, psi, matr, res, xpsi, xres)
  type(mesh_t),      intent(in)  :: mesh
  type(states_t),    intent(in)  :: st
  R_TYPE,            intent(in)  :: psi(:, :, :)
  R_TYPE,            intent(in)  :: matr(:, :)
  R_TYPE,            intent(out) :: res(:, :, :)
  integer, optional, intent(in)  :: xpsi(:), xres(:)

  call push_sub('states_inc.Xstates_block_matr_mul')

  call X(states_block_matr_mul_add)(mesh, st, R_TOTYPE(M_ONE), psi, matr, &
    R_TOTYPE(M_ZERO), res, xpsi, xres)

  call pop_sub()
end subroutine X(states_block_matr_mul)


! ---------------------------------------------------------
! Multiplication of block of states by matrix plus block of states
! (with the corresponding column indices):
! res(xres) <- alpha * psi(xpsi) * matr + beta * res(xres).
subroutine X(states_block_matr_mul_add)(mesh, st, alpha, psi, matr, beta, res, xpsi, xres)
  type(mesh_t),      intent(in)    :: mesh
  type(states_t),    intent(in)    :: st
  R_TYPE,            intent(in)    :: alpha
  R_TYPE,            intent(in)    :: psi(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  R_TYPE,            intent(in)    :: matr(:, :)
  R_TYPE,            intent(in)    :: beta
  R_TYPE,            intent(inout) :: res(mesh%np_part, st%d%dim, st%st_start:st%st_end)
  integer, optional, intent(in)    :: xpsi(:), xres(:)

  integer              :: res_col, psi_col, matr_col, i
  integer, allocatable :: xpsi_(:), xres_(:)
  R_TYPE, allocatable  :: res_block(:, :, :), matr_block(:, :), psi_block(:, :, :)
#if defined(HAVE_MPI)
  integer              :: rank, size, node, ist, round, matr_row_offset, matr_col_offset
  integer              :: src, dst, left, right, k, l, idim, sendcnt, recvcnt, max_count
  integer              :: stats(MPI_STATUS_SIZE, 2), reqs(2)
  integer, allocatable :: xpsi_count(:), xres_count(:), xpsi_node(:, :), xres_node(:, :)
  R_TYPE, pointer      :: sendbuf(:, :, :), recvbuf(:, :, :), tmp_ptr(:, :, :)
#endif

  call profiling_in(C_PROFILING_BLOCK_MATR)
  call push_sub('states_inc.Xstates_block_matr_add')

  ! Calculate gobal index sets.
  if(present(xpsi)) then
    psi_col = ubound(xpsi, 1)
    ALLOCATE(xpsi_(psi_col), psi_col)
    xpsi_ = xpsi
  else
    psi_col = st%nst
    ALLOCATE(xpsi_(psi_col), psi_col)
    do i = 1, psi_col
      xpsi_(i) = i
    end do
  end if
  if(present(xres)) then
    res_col = ubound(xres, 1)
    ALLOCATE(xres_(res_col), res_col)
    xres_ = xres
  else
    res_col = st%nst
    ALLOCATE(xres_(res_col), res_col)
    do i = 1, res_col
      xres_(i) = i
    end do
  end if

  matr_col = ubound(matr, 2)

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    ! Shortcuts.
    size = st%mpi_grp%size
    rank = st%mpi_grp%rank

    ! Calculate the index sets per node.
    ! xpsi_node(1:xpsi_count(node), node) and
    ! xres_node(1:xres_count(node), node) are the index sets.
    ALLOCATE(xpsi_count(0:size-1), size)
    ALLOCATE(xres_count(0:size-1), size)
    ! Count the how many vectors each node has.
    xpsi_count = 0
    xres_count = 0
    do i = 1, psi_col
      xpsi_count(st%node(xpsi_(i))) = xpsi_count(st%node(xpsi_(i))) + 1
    end do
    do i = 1, res_col
      xres_count(st%node(xres_(i))) = xres_count(st%node(xres_(i))) + 1
    end do
    ! Allocate space, it is a bit more than really required but makes the code simpler.
    ALLOCATE(xpsi_node(maxval(xpsi_count), 0:size-1), maxval(xpsi_count)*size)
    ALLOCATE(xres_node(maxval(xres_count), 0:size-1), maxval(xres_count)*size)
    ! Now set up the index sets.
    xpsi_count = 0
    xres_count = 0
    xpsi_node  = 0
    xres_node  = 0
    do ist = 1, st%nst
      node = st%node(ist)
      ! A state ist is only included if its in the global index set.
      if(member(ist, xpsi_)) then
        xpsi_count(node)                  = xpsi_count(node) + 1
        xpsi_node(xpsi_count(node), node) = ist
      end if
      if(member(ist, xres_)) then
        xres_count(node)                  = xres_count(node) + 1
        xres_node(xres_count(node), node) = ist
      end if
    end do

    ! Take care of beta first, if necessary, and compact res to res_block.
    if(beta.ne.R_TOTYPE(M_ZERO)) then
      call profiling_in(C_PROFILING_BLOCK_MATR_CP)
      do i = 1, xres_count(rank)
        do idim = 1, st%d%dim
          call lalg_scal(mesh%np, beta, res(:, idim, xres_node(i, rank)))
        end do
      end do
      call profiling_out(C_PROFILING_BLOCK_MATR_CP)
    else
      res = R_TOTYPE(M_ZERO)
    end if
    ALLOCATE(res_block(mesh%np, st%d%dim, xres_count(rank)), mesh%np*st%d%dim*xres_count(rank))
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st, mesh, xres_node(1:xres_count(rank), rank), res, res_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! Allocate send and receive buffers. For some blocks, they are oversized but
    ! buffer reusing saves some overhead in the MPI library.
    max_count = maxval(xpsi_count)
    ALLOCATE(sendbuf(mesh%np, st%d%dim, max_count), mesh%np*st%d%dim*max_count)
    ALLOCATE(recvbuf(mesh%np, st%d%dim, max_count), mesh%np*st%d%dim*max_count)

    ! Compact the local block to send away.
    sendcnt = xpsi_count(rank)
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st, mesh, xpsi_node(1:xpsi_count(rank), rank), psi(:, :, :), sendbuf)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! Get neighbours.
    call MPI_Cart_shift(st%dom_st%comm, P_STRATEGY_STATES-1, -1, src, dst, mpi_err)
    right = lmpi_translate_rank(st%dom_st%comm, st%mpi_grp%comm, src)
    left  = lmpi_translate_rank(st%dom_st%comm, st%mpi_grp%comm, dst)

    do round = 0, size-1
      k = mod(rank+round, size)   ! The column of the block currently being calculated.
      l = mod(rank+round+1, size) ! The column of the block currently being communicated.
      ! In all but the first rounds we have to wait for the data to arrive and
      ! then swap buffers, i. e., what we received in the send buffer in the remainder
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
        recvcnt = xpsi_count(l)
        call MPI_Irecv(recvbuf, mesh%np*st%d%dim*recvcnt, R_MPITYPE, right, 0, st%mpi_grp%comm, reqs(1), mpi_err)
        call MPI_Isend(sendbuf, mesh%np*st%d%dim*sendcnt, R_MPITYPE, left, 0, st%mpi_grp%comm, reqs(2), mpi_err)
      end if

      ! Do the matrix multiplication.
      matr_row_offset = sum(xpsi_count(0:k-1))
      matr_col_offset = sum(xres_count(0:rank-1))
      ALLOCATE(matr_block(xpsi_count(k), xres_count(rank)), xpsi_count(k)*xres_count(rank))

      if(sendcnt.gt.0.and.xres_count(rank).gt.0) then
        call profiling_in(C_PROFILING_BLOCK_MATR_CP)
        matr_block = matr(matr_row_offset+1:matr_row_offset+xpsi_count(k), matr_col_offset+1:matr_col_offset+xres_count(rank))
        call profiling_out(C_PROFILING_BLOCK_MATR_CP)
        call profiling_in(C_PROFILING_BLOCK_MATR_MM)
        call lalg_gemm(mesh%np*st%d%dim, xres_count(rank), sendcnt, alpha, &
          sendbuf(:, :, 1), matr_block, R_TOTYPE(M_ONE), res_block(:, :, 1))
        call profiling_out(C_PROFILING_BLOCK_MATR_MM)
      end if

      deallocate(matr_block)
    end do

    ! Copy result.
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_uncompactify)(st, mesh, xres_node(1:xres_count(rank), rank), res_block, res)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    deallocate(res_block)
    deallocate(sendbuf, recvbuf)
    deallocate(xpsi_count, xres_count, xpsi_node, xres_node)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call write_fatal(1)
#endif
  else
    ALLOCATE(res_block(mesh%np, st%d%dim, res_col), mesh%np*st%d%dim*res_col)
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st, mesh, xres_(1:res_col), res, res_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ALLOCATE(psi_block(mesh%np, st%d%dim, psi_col), mesh%np*st%d%dim*psi_col)
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_compactify)(st, mesh, xpsi_(1:psi_col), psi, psi_block)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)

    ! matr_block is needed because matr may be an assumed shape array.
    ALLOCATE(matr_block(psi_col, matr_col), psi_col*matr_col)
    matr_block = matr
    call lalg_gemm(mesh%np*st%d%dim, matr_col, psi_col, alpha, &
      psi_block(:, :, 1), matr_block, beta, res_block(:, :, 1))
    deallocate(matr_block)

    ! Copy result.
    call profiling_in(C_PROFILING_BLOCK_MATR_CP)
    call X(states_uncompactify)(st, mesh, xres_(1:res_col), res_block, res)
    call profiling_out(C_PROFILING_BLOCK_MATR_CP)
    deallocate(res_block, psi_block)
  end if

  deallocate(xpsi_, xres_)

  call pop_sub()
  call profiling_out(C_PROFILING_BLOCK_MATR)
end subroutine X(states_block_matr_mul_add)


! ---------------------------------------------------------
! Copy in(:, :, idx) to out out(:, :, 1:ubound(idx, 1)).
subroutine X(states_compactify)(st, m, idx, in, out)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: m
  integer,        intent(in)  :: idx(:)
  R_TYPE,         intent(in)  :: in(:, :, :)
  R_TYPE,         intent(out) :: out(:, :, :)

  integer :: ist, idim, n

  call push_sub('states_block_inc.X(states_compactify)')

  n = ubound(idx, 1)

  do ist = 1, n
    do idim = 1, st%d%dim
      call lalg_copy(m%np, in(:, idim, idx(ist)-st%st_start+1), out(:, idim, ist))
    end do
  end do
  
  call pop_sub()
end subroutine X(states_compactify)


! ---------------------------------------------------------
! Undo the effect of X(states_compactify), i. e.
! X(states_compactify)(st, m, idx, in, out) followed by
! X(states_uncompactify)(st, m, idx, out, in) is identity.
subroutine X(states_uncompactify)(st, m, idx, in, out)
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: m
  integer,        intent(in)  :: idx(:)
  R_TYPE,         intent(in)  :: in(:, :, :)
  R_TYPE,         intent(out) :: out(:, :, :)

  integer :: ist, idim, n

  call push_sub('states_block_inc.X(states_uncompactify)')

  n = ubound(idx, 1)

  do ist = 1, n
    do idim = 1, st%d%dim
      call lalg_copy(m%np, in(:, idim, ist), out(:, idim, idx(ist)-st%st_start+1))
    end do
  end do
  
  call pop_sub()
end subroutine X(states_uncompactify)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
