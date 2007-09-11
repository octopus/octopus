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
! Every node has incount (may vary from node to node) items (in
! array in) to send to everybody alse in the group. The total
! number of items in the out array is given by outcount. out has
! to be big enough to contain all possible incoming items.
subroutine X(lmpi_gen_alltoallv)(incount, in, outcount, out, mpi_grp)
  integer,         intent(in)  :: incount
  R_TYPE,          intent(in)  :: in(:)
  integer,         intent(out) :: outcount
  R_TYPE,          intent(out) :: out(:)
  type(mpi_grp_t), intent(in)  :: mpi_grp

  integer                          :: mpi_err
  integer, dimension(mpi_grp%size) :: sendcnts, sdispls, rdispls, recvbuf, recvcnts

  call push_sub('mpi_lib_inc.Xlmpi_gen_alltoallv')

  ! Query how many elements each node has to contribute.
  recvcnts                = 0
  sdispls                 = 0
  sendcnts                = 1
  recvcnts                = 1
  rdispls(1)              = 0
  rdispls(2:mpi_grp%size) = rdispls(1:mpi_grp%size)+recvcnts(1:mpi_grp%size)
  call MPI_Debug_In(mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(incount, sendcnts, sdispls, R_MPITYPE, recvbuf, recvcnts, rdispls, &
    R_MPITYPE, mpi_grp%comm, mpi_err)
  call MPI_Debug_Out(mpi_grp%comm, C_MPI_ALLTOALLV)

  outcount = sum(recvbuf)

  ! Exchange the data.
  recvcnts                = recvbuf
  sendcnts                = incount
  sdispls                 = 0
  rdispls(1)              = 0
  rdispls(2:mpi_grp%size) = rdispls(1:mpi_grp%size)+recvcnts(1:mpi_grp%size)
  call MPI_Debug_In(mpi_grp%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(in, sendcnts, sdispls, R_MPITYPE, out, recvcnts, rdispls, &
    R_MPITYPE, mpi_grp%comm, mpi_err)
  call MPI_Debug_Out(mpi_grp%comm, C_MPI_ALLTOALLV)

  call pop_sub()
end subroutine X(lmpi_gen_alltoallv)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
