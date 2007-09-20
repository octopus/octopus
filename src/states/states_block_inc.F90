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
      ! probably better.
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


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
