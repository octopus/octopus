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
#if defined(HAVE_MPI)
  integer              :: size, rank, node, round, partner, ist, res_col_offset, res_row_offset
  integer, allocatable :: xpsi1_count(:), xpsi2_count(:), xpsi1_node(:, :), xpsi2_node(:, :)
  integer              :: stat(MPI_STATUS_SIZE)
  R_TYPE, allocatable  :: sendbuf(:, :, :), recvbuf(:, :, :), res_tmp(:, :)
#endif

  call profiling_in(C_PROFILING_LOBPCG_BLOCKT)
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

    ! The local part.
    res_row_offset = sum(xpsi1_count(0:rank-1))
    res_col_offset = sum(xpsi2_count(0:rank-1))
    if(symm_) then
      do i = 1, xpsi1_count(rank)
        res(i+res_row_offset, i+res_col_offset) =                         &
          X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(i, rank)), &
          psi2(:, :, xpsi2_node(i, rank)), reduce=.false.)
        do j = i+1, xpsi2_count(rank)
          res(i+res_row_offset, j+res_col_offset) =                         &
            X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(i, rank)), &
            psi2(:, :, xpsi2_node(j, rank)), reduce=.false.)
          res(j+res_col_offset, i+res_row_offset) = R_CONJ(res(i+res_row_offset, j+res_col_offset))
        end do
      end do
    else
      do i = 1, xpsi1_count(rank)
        do j = 1, xpsi2_count(rank)
          res(i+res_row_offset, j+res_col_offset) =                         &
            X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(i, rank)), &
            psi2(:, :, xpsi2_node(j, rank)), reduce=.false.)
        end do
      end do
    end if
  
    ! The remote parts.
    do round = 1, st%ap%rounds
      partner = st%ap%schedule(rank, round)

      if(partner.eq.rank) then ! We are idle this round.
        ! Just do nothing.
      else ! We are not idle.
        ! Both nodes send psi2 to each other and then calculate one block of the
        ! result.
        ALLOCATE(sendbuf(mesh%np, st%d%dim, xpsi2_count(rank)), mesh%np*st%d%dim*xpsi2_count(rank))
        ALLOCATE(recvbuf(mesh%np, st%d%dim, xpsi2_count(partner)), mesh%np*st%d%dim*xpsi2_count(partner))

        ! Copy block to sendbuffer. This is necessaery because the vectors in the block may not
        ! be contiguous.
        do i = 1, xpsi2_count(rank)
          call lalg_copy(mesh%np*st%d%dim, psi2(:, 1, xpsi2_node(i, rank)), sendbuf(:, 1, i))
        end do

        call MPI_Sendrecv(sendbuf, mesh%np*st%d%dim*xpsi2_count(rank), R_MPITYPE, partner, 12, &
          recvbuf, mesh%np*st%d%dim*xpsi2_count(partner), R_MPITYPE, partner, 12, st%mpi_grp%comm, stat, mpi_err)

        res_row_offset = sum(xpsi1_count(0:rank-1))
        res_col_offset = sum(xpsi2_count(0:partner-1))
        ! Do the matrix multiplication.
        do i = 1, xpsi1_count(rank)
          do j = 1, xpsi2_count(partner)
            res(i+res_row_offset, j+res_col_offset) = &
              X(states_dotp)(mesh, st%d%dim, psi1(:, :, xpsi1_node(i, rank)), recvbuf(:, :, j), reduce=.false.)
          end do
        end do
        deallocate(sendbuf, recvbuf)
      end if
    end do

    ! Add up all the individual blocks.
    ALLOCATE(res_tmp(psi1_col, psi2_col), psi1_col*psi2_col)
    call MPI_Allreduce(res, res_tmp, psi1_col*psi2_col, R_MPITYPE, MPI_SUM, st%mpi_grp%comm, mpi_err)
    ! What is the problem with this BLAS call?
    ! call lalg_copy(psi1_col*psi2_col, res_tmp(:, 1), res(:, 1))
    res = res_tmp
    deallocate(res_tmp)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call write_fatal(1)
#endif
  else ! No states parallelization.
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
  ! FIXME: test this!
  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    ALLOCATE(res_tmp(psi1_col, psi2_col), psi1_col*psi2_col)
    call MPI_Allreduce(res, res_tmp, psi1_col*psi2_col, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    call lalg_copy(psi1_col*psi2_col, res_tmp(:, 1), res(:, 1))
    deallocate(res_tmp)
#endif
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_LOBPCG_BLOCKT)
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

  call MPI_Debug_In(st%mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(in(:, 1, 1), sendcnts, sdispls, R_MPITYPE, &
    out(:, 1, 1), recvcnts, rdispls, R_MPITYPE, st%mpi_grp%comm, mpi_err)
  call MPI_Debug_Out(st%mpi_grp%comm, C_MPI_ALLTOALLV)
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

  integer              :: res_col, psi_col, matr_col, i, j, k, idim
  R_TYPE               :: tmp
  integer, allocatable :: xpsi_(:), xres_(:)
#if defined(HAVE_MPI)
  integer              :: rank, size, node, ist, round, to, from, matr_row_offset, matr_col_offset
  integer              :: req, stat(MPI_STATUS_SIZE)
  integer, allocatable :: xpsi_count(:), xres_count(:), xpsi_node(:, :), xres_node(:, :)
  R_TYPE, allocatable  :: buf(:, :, :), newblock(:, :, :)
#endif

  call profiling_in(C_PROFILING_LOBPCG_BLOCK_MATR)
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
    ! Allocate space, it is a bit more than really required.
    ALLOCATE(xpsi_node(maxval(xpsi_count), 0:size-1), maxval(xpsi_count)*size)
    ALLOCATE(xres_node(maxval(xres_count), 0:size-1), maxval(xres_count)*size)
    ! Now set up the index sets.
    xpsi_count = 0
    xres_count = 0
    xpsi_node  = 0
    xres_node  = 0
    do ist = 1, st%nst
      node = st%node(ist)
      ! A state ist is ony included if its in the global index set.
      if(member(ist, xpsi_)) then
        xpsi_count(node)                  = xpsi_count(node) + 1
        xpsi_node(xpsi_count(node), node) = ist
      end if
      if(member(ist, xres_)) then
        xres_count(node)                  = xres_count(node) + 1
        xres_node(xres_count(node), node) = ist
      end if
    end do

    ! The local block. This calculates the diagonal blocks of the result.
    matr_row_offset = sum(xpsi_count(0:rank-1))
    matr_col_offset = sum(xres_count(0:rank-1))
    if(beta.eq.R_TOTYPE(M_ZERO)) then
      do j = 1, xres_count(rank)
        do idim = 1, st%d%dim
          do i = 1, mesh%np
            tmp = R_TOTYPE(M_ZERO)
            do k = 1, xpsi_count(rank)
              tmp = tmp + psi(i, idim, xpsi_node(k, rank))*matr(k+matr_row_offset, j+matr_col_offset)
            end do
            res(i, idim, xres_node(j, rank)) = alpha*tmp
          end do
        end do
      end do
    else
      do j = 1, xres_count(rank)
        do idim = 1, st%d%dim
          do i = 1, mesh%np
            tmp = R_TOTYPE(M_ZERO)
            do k = 1, xpsi_count(rank)
              tmp = tmp + psi(i, idim, xpsi_node(k, rank))*matr(k+matr_row_offset, j+matr_col_offset)
            end do
            res(i, idim, xres_node(j, rank)) = beta*res(i, idim, xres_node(j, rank)) + alpha*tmp
          end do
        end do
      end do
    end if

    ! The remote blocks (offdiagonal blocks).
    ! Buffer for the blocks received from other nodes.
    ALLOCATE(newblock(mesh%np, st%d%dim, xres_count(rank)), mesh%np*st%d%dim*xres_count(rank))

    ! Communication is done in a ring fashion: 0 sends to 1, 1 to 2, etc. in the first round,
    ! 0 to 2, 1 to 3, etc. in the second etc. In the whole, size-1 blocks have to be
    ! calculated, because the diagonal block has already been done.
    do round = 1, size-1
      to   = mod(rank + round, size)
      from = mod(size + rank - round, size)

      ! MPI_Sendrecv might be interesting here (because we do a shift).
      ! But interleaving computation with non-blocking communication is
      ! probaby better.
      if(xres_count(to).gt.0) then
        matr_row_offset = sum(xpsi_count(0:rank-1))
        matr_col_offset = sum(xres_count(0:to-1))
        ! Allocate sendbuffer.
        ALLOCATE(buf(mesh%np, st%d%dim, xres_count(to)), mesh%np*st%d%dim*xres_count(to))
        do j = 1, xres_count(to)
          do idim = 1, st%d%dim
            do i = 1, mesh%np
              tmp = R_TOTYPE(M_ZERO)
              do k = 1, xpsi_count(rank)
                tmp = tmp + psi(i, idim, xpsi_node(k, rank))*matr(k+matr_row_offset, j+matr_col_offset)
              end do
              buf(i, idim, j) = alpha*tmp
            end do
          end do
        end do
        ! Asynchronous send.
        call MPI_Isend(buf, mesh%np*st%d%dim*xres_count(to), R_MPITYPE, to, 0, &
          st%mpi_grp%comm, req, mpi_err)
      end if
      ! And blocking receive. Otherwise, the addition of the new block to the
      ! local result cannot be done.
      if(xres_count(rank).gt.0) then
        call MPI_Recv(newblock, mesh%np*st%d%dim*xres_count(rank), R_MPITYPE, from, 0, &
          st%mpi_grp%comm, stat, mpi_err)
        do i = 1, xres_count(rank)
          call lalg_axpy(mesh%np, R_TOTYPE(M_ONE), newblock(:, 1, i), res(:, 1, xres_node(i, rank)))
        end do
      end if
      ! Wait for the sending to be complete.
      if(xres_count(to).gt.0) then
        call MPI_Wait(req, stat, mpi_err)
        deallocate(buf)
      end if
    end do
    deallocate(newblock)
#else
    message(1) = 'Running gs parallel in states without MPI. This is a bug!'
    call write_fatal(1)
#endif
  else
    ! We can shortcut for beta being zero.
    if(beta.eq.R_TOTYPE(M_ZERO)) then
      call profiling_in(C_PROFILING_LOBPCG_LOOP)
      do j = 1, matr_col
        do idim = 1, st%d%dim
          do i = 1, mesh%np
            res(i, idim, xres_(j)) = R_TOTYPE(M_ZERO)
            do k = 1, psi_col
              res(i, idim, xres_(j)) = res(i, idim, xres_(j)) + psi(i, idim, xpsi_(k))*matr(k, j)
            end do
            res(i, idim, xres_(j)) = alpha*res(i, idim, xres_(j))
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
              res(i, idim, xres_(j)) = res(i, idim, xres_(j)) + psi(i, idim, xpsi_(k))*matr(k, j)
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
              tmp = tmp + psi(i, idim, xpsi_(k))*matr(k, j)
            end do
            if(beta.eq.R_TOTYPE(M_ZERO)) then
              res(i, idim, xres_(j)) = alpha*tmp
            else
              res(i, idim, xres_(j)) = alpha*tmp + beta*res(i, idim, xres_(j))
            end if
          end do
        end do
      end do
      call profiling_out(C_PROFILING_LOBPCG_LOOP)
    end if
  end if

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    deallocate(xpsi_count, xres_count, xpsi_node, xres_node)
#endif
  end if

  deallocate(xpsi_, xres_)

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
R_TYPE function X(states_dotp)(m, dim, f1, f2, reduce) result(dotp)
  type(mesh_t),      intent(in) :: m
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  call push_sub('states_inc.Xstates_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim), reduce)
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

  integer             :: idim, ist, ik, i
  CMPLX               :: expect_val_p
  R_TYPE, allocatable :: grad(:,:,:)  
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT               :: lmomentum(NDIM, st%lnst)
#endif

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

    ! Exchange momenta in the state parallel case.
#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      lmomentum = st%momentum(:, st%st_start:st%st_end, ik)
      call lmpi_gen_alltoallv(NDIM*st%lnst, lmomentum(:, 1), tmp, &
        st%momentum(:, 1, ik), st%mpi_grp)
    end if
#endif
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
