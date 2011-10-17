!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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
  integer, intent(in)  :: comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_BARRIER)

  call MPI_Barrier(comm, ierr)

  call mpi_debug_out(comm, C_MPI_BARRIER)
end subroutine TS(MPI_Barrier)


! ---------------------------------------------------------
subroutine TS(MPI_Scatterv)(sendbuf, sendcnts, displs, sendtype, recvbuf, &
  recvcount, recvtype, root, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcnts(:), displs(:), sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount
  integer, intent(in)  :: recvtype, root, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_SCATTERV)

  call MPI_Scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, &
    recvcount, recvtype, root, comm, ierr)

  call mpi_debug_out(comm, C_MPI_SCATTERV)
end subroutine TS(MPI_Scatterv)


! ---------------------------------------------------------
subroutine TS(MPI_Gatherv)(sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, root, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcnts, sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount(:), displs(:)
  integer, intent(in)  :: recvtype, root, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_GATHERV)

  call MPI_Gatherv(sendbuf, sendcnts, sendtype, recvbuf, &
    recvcount, displs, recvtype, root, comm, ierr)

  call mpi_debug_out(comm, C_MPI_GATHERV)
end subroutine TS(MPI_Gatherv)


! ---------------------------------------------------------
subroutine TS(MPI_Alltoallv)(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
  recvcount, rdispls, recvtype, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcnts(:), sdispls(:), sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount(:), rdispls(:), recvtype, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_ALLTOALLV)

  call MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype, recvbuf, &
    recvcount, rdispls, recvtype, comm, ierr)

  call mpi_debug_out(comm, C_MPI_ALLTOALLV)
end subroutine TS(MPI_Alltoallv)


! ---------------------------------------------------------
subroutine TS(MPI_Alltoall)(sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, recvtype, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcnts(:), sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount(:), recvtype, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_ALLTOALLV)

  call MPI_Alltoall(sendbuf, sendcnts, sendtype, recvbuf, &
    recvcount, recvtype, comm, ierr)

  call mpi_debug_out(comm, C_MPI_ALLTOALLV)
end subroutine TS(MPI_Alltoall)


! ---------------------------------------------------------
subroutine TS(MPI_Allgatherv)(sendbuf, sendcnts, sendtype, recvbuf, &
  recvcount, displs, recvtype, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcnts, sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount(:), displs(:)
  integer, intent(in)  :: recvtype, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_ALLGATHERV)

  call MPI_Allgatherv(sendbuf, sendcnts, sendtype, recvbuf, &
    recvcount, displs, recvtype, comm, ierr)

  call mpi_debug_out(comm, C_MPI_ALLGATHERV)
end subroutine TS(MPI_Allgatherv)


! ---------------------------------------------------------
subroutine TS(MPI_Bcast)(buf, cnt, sendtype, root, comm, ierr)
  R_TYPE,  intent(inout) :: buf(:)
  integer, intent(in)    :: cnt, sendtype, root, comm
  integer, intent(out)   :: ierr

  call mpi_debug_in(comm, C_MPI_BCAST)

  call MPI_Bcast(buf, cnt, sendtype, root, comm, ierr)

  call mpi_debug_out(comm, C_MPI_BCAST)
end subroutine TS(MPI_Bcast)


! ---------------------------------------------------------
subroutine TS(MPI_Allreduce)(sendbuf, recvbuf, count, datatype, op, &
  comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf
  R_TYPE,  intent(out) :: recvbuf
  integer, intent(in)  :: count, datatype, op, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_ALLREDUCE)

  call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
    comm, ierr)

  call mpi_debug_out(comm, C_MPI_ALLREDUCE)
end subroutine TS(MPI_Allreduce)


! ---------------------------------------------------------
subroutine TS(MPI_Allgather)(sendbuf, sendcount, sendtype, recvbuf, &
  recvcount, recvtype, comm, ierr)
  R_TYPE,  intent(in)  :: sendbuf(:)
  integer, intent(in)  :: sendcount, sendtype
  R_TYPE,  intent(out) :: recvbuf(:)
  integer, intent(in)  :: recvcount, recvtype, comm
  integer, intent(out) :: ierr

  call mpi_debug_in(comm, C_MPI_ALLGATHER)

  call MPI_Allgather(sendbuf, sendcount, sendtype,&
    recvbuf, recvcount, recvtype, comm, ierr)

  call mpi_debug_out(comm, C_MPI_ALLGATHER)
end subroutine TS(MPI_Allgather)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
