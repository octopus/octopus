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


! ---------------------------------------------------------
!> Every node has incount (may vary from node to node) items (in
!! array in) to send to everybody else in the group. The total
!! number of items in the out array is given by outcount. out has
!! to be big enough to contain all possible incoming items.
subroutine X(lmpi_gen_allgatherv)(incount, in, outcount, out, mpi_grp)
  integer,         intent(in)  :: incount
  R_TYPE,          intent(in)  :: in(:)
  integer,         intent(out) :: outcount
  R_TYPE,          intent(out) :: out(:)
  type(mpi_grp_t), intent(in)  :: mpi_grp

  integer              :: mpi_err, inode
  integer, allocatable :: rdispls(:), recvbuf(:), recvcnts(:)

  PUSH_SUB(X(lmpi_gen_allgatherv))

  SAFE_ALLOCATE( rdispls(1:mpi_grp%size))
  SAFE_ALLOCATE( recvbuf(1:mpi_grp%size))
  SAFE_ALLOCATE(recvcnts(1:mpi_grp%size))

  ! Query how many elements each node has to contribute.
  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHER)
  call MPI_Allgather(incount, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, mpi_grp%comm, mpi_err)
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHER)

  outcount = sum(recvbuf)

  ! Exchange the data.
  recvcnts = recvbuf
  rdispls(1) = 0
  ! Never try to implement the following loop as a vector assignment like
  ! (caution: wrong code!):
  ! rdispls(2:mpi_grp%size) = rdispls(1:mpi_grp%size-1)+recvcnts(1:mpi_grp%size-1)
  ! (Took me an hour to find the mistake...)
  do inode = 2, mpi_grp%size
    rdispls(inode) = rdispls(inode-1)+recvcnts(inode-1)
  end do

  call mpi_debug_in(mpi_grp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv(in, incount, R_MPITYPE, out, recvcnts, rdispls, R_MPITYPE, mpi_grp%comm, mpi_err)
  call mpi_debug_out(mpi_grp%comm, C_MPI_ALLGATHERV)

  SAFE_DEALLOCATE_A(rdispls)
  SAFE_DEALLOCATE_A(recvbuf)
  SAFE_DEALLOCATE_A(recvcnts)

  POP_SUB(X(lmpi_gen_allgatherv))
end subroutine X(lmpi_gen_allgatherv)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
