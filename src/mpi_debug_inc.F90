!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine TS(MPI_Barrier)(comm, ierr)
  integer :: comm, ierr

  call MPI_Debug_IN (comm, C_MPI_BARRIER)

  call MPI_Barrier(comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_BARRIER)

end subroutine TS(MPI_Barrier)


! ---------------------------------------------------------
subroutine TS(MPI_Scatterv)(sendbuf, sendcnts, displs, sendtype, recvbuf, &
  recvcount, recvtype, root, comm, ierr)

  R_TYPE  :: sendbuf(:), recvbuf(:)
  integer :: sendcnts(:), displs(:)
  integer :: recvcount, sendtype, recvtype, root, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_SCATTERV)

  call MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, &
    recvcount, recvtype, root, comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_SCATTERV)

end subroutine TS(MPI_Scatterv)


! ---------------------------------------------------------
subroutine TS(MPI_Gatherv)(sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, root, comm, ierr)

  R_TYPE  :: sendbuf(:), recvbuf(:)
  integer :: recvcount(:), displs(:)
  integer :: sendcnts, sendtype, recvtype, root, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_GATHERV)

  call MPI_Gatherv(sendbuf, sendcnts, sendtype, recvbuf, &
    recvcount, displs, recvtype, root, comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_GATHERV)

end subroutine TS(MPI_Gatherv)


! ---------------------------------------------------------
subroutine TS(MPI_Alltoallv)(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype, comm, ierr)

  R_TYPE  :: sendbuf(:), recvbuf(:)
  integer :: sendcnts(:), sdispls(:), recvcount(:), rdispls(:)
  integer :: sendtype, recvtype, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_ALLTOALLV)

  call MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
    recvcount, rdispls, recvtype, comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_ALLTOALLV)

end subroutine TS(MPI_Alltoallv)


! ---------------------------------------------------------
subroutine TS(MPI_Allgatherv)(sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, comm, ierr)

  R_TYPE  :: sendbuf(:), recvbuf(:)
  integer :: recvcount(:), displs(:)
  integer :: sendcnts, sendtype, recvtype, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_ALLGATHERV)

  call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
    recvcount, displs, recvtype, comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_ALLGATHERV)

end subroutine TS(MPI_Allgatherv)


! ---------------------------------------------------------
subroutine TS(MPI_Bcast)(buf, cnt, sendtype, root, comm, ierr)
  R_TYPE  :: buf(:)
  integer :: cnt, sendtype, root, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_BCAST)

  call MPI_Bcast(buf, cnt, sendtype, root, comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_BCAST)

end subroutine TS(MPI_Bcast)


! ---------------------------------------------------------
subroutine TS(MPI_Allreduce)(sendbuf, recvbuf, count, datatype, op, &
  comm, ierr)

  R_TYPE  :: sendbuf, recvbuf
  integer :: count, datatype, op, comm, ierr

  call MPI_Debug_IN (comm, C_MPI_ALLREDUCE)

  call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
    comm, ierr)

  call MPI_Debug_OUT(comm, C_MPI_ALLREDUCE)

end subroutine TS(MPI_Allreduce)
