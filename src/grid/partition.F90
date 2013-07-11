!! Copyright (C) 2013 M. Oliveira
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: par_vec.F90 10505 2013-05-03 18:56:16Z dstrubbe $

#include "global.h"
 
module partition_m
  use global_m
  use io_binary_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use profiling_m


  implicit none

  private

  public ::                        &
    partition_t,                   &
    partition_init,                &
    partition_end,                 &
    partition_set,                 &
    partition_write,               &
    partition_read,                &
    partition_get_local_size,      &
    partition_get_global,          &
    partition_get_partition_number


  !> The partition is an array that contains the mapping between some global index 
  !! and a process, such that point ip will be stored in process partition(ip).
  !!
  !! In this module this array is distributed among the processes, such that each process 
  !! only stores a portion of the full array. Because each process needs to know in a 
  !! straighforward way (i.e. without having to perform any kind of communication or 
  !! lengthy operations) which process stores the partition corresponding to any giving
  !! point, all the processes except one (the one with the highest rank) store exactly the 
  !! same number of points.
  !!
  !! Note 1: this module can be a bit confusing as they are in fact two partitions. One is the
  !! partition of some array (in the case of Octopus, this is typically the mesh functions),
  !! which is the main information stored in the partition_t object, and then there is the 
  !! partition of the partition itself, as this is also distributed.
  !! 
  !! Note 2: in principle, the mpi group used by the processes for the partition distribution
  !! does not need to be the same as the mpi group of the processes the partition refers to.
  type partition_t
    private

    ! The following components are the same for all processes:
    type(mpi_grp_t) :: mpi_grp   !< The mpi group use for distributing the partition data.
    integer ::         np_global !< The total number of points in the partition.
    integer ::         nppp      !< Number of points per process. It is different from
                                 !! np_local only for the process with the highest rank.
    integer ::         npart     !< The number of partitions.

    ! The following components are process dependent:
    integer :: np_local          !< The number of points of the partition stored in this process.
    integer :: istart            !< The position of the first point stored in this process.
    integer, pointer :: part(:)  !< The local portion of the partition.

  end type partition_t


contains

  !---------------------------------------------------------
  subroutine partition_init(partition, np_global, mpi_grp)
    type(partition_t), intent(out) :: partition
    integer,           intent(in)  :: np_global

    type(mpi_grp_t),   intent(in)  :: mpi_grp

    PUSH_SUB(partition_init)

    !Global variables
    partition%mpi_grp = mpi_grp
    partition%np_global = np_global
    partition%nppp = ceiling(real(np_global)/real(mpi_grp%size))
    partition%npart = mpi_grp%size

    !Processor dependent
    if (mpi_grp%rank + 1 == partition%npart) then
      partition%np_local = mod(np_global, partition%nppp)
      if (partition%np_local == 0) partition%np_local = partition%nppp
    else
      partition%np_local = partition%nppp
    end if
    partition%istart = partition%nppp*mpi_grp%rank + 1

    !Allocate memory for the partition
    nullify(partition%part)
    SAFE_ALLOCATE(partition%part(partition%np_local))

    POP_SUB(partition_init)
  end subroutine partition_init

  ! ---------------------------------------------------------
  subroutine partition_end(partition)
    type(partition_t), intent(inout) :: partition

    PUSH_SUB(partition_end)

    SAFE_DEALLOCATE_P(partition%part)

    POP_SUB(partition_end)
  end subroutine partition_end

  ! ---------------------------------------------------------
  subroutine partition_set(partition, part)
    type(partition_t), intent(inout) :: partition
    integer,           intent(in)    :: part(:) !< The local portion of the partition.

    PUSH_SUB(partition_set)

    partition%part(1:partition%np_local) = part(1:partition%np_local)

    POP_SUB(partition_set)
  end subroutine partition_set

  ! ---------------------------------------------------------
  subroutine partition_write(partition, filename)
    type(partition_t), intent(in) :: partition
    character(len=*),  intent(in) :: filename

    integer :: ierr
    integer, allocatable :: part_global(:)

    PUSH_SUB(partition_write)

    !Get the global partition in the root node
    if (partition%mpi_grp%rank == 0) then
      SAFE_ALLOCATE(part_global(partition%np_global))
    else
      SAFE_ALLOCATE(part_global(1))
    end if
    call partition_get_global(partition, part_global, 0)
    
    !Only the root node writes
    if (partition%mpi_grp%rank == 0) then
      call io_binary_write(filename, partition%np_global, part_global, ierr)
    end if

    SAFE_DEALLOCATE_A(part_global)

    POP_SUB(partition_write)
  end subroutine partition_write

  ! ---------------------------------------------------------
  subroutine partition_read(partition, filename, ierr)
    type(partition_t), intent(inout) :: partition
    character(len=*),  intent(in)    :: filename
    integer,           intent(out)   :: ierr

    integer :: ipart
    integer, allocatable :: part_global(:)
    integer, allocatable :: scounts(:), sdispls(:)

    PUSH_SUB(partition_read)

    ! The global partition is only read by the root node
    if (partition%mpi_grp%rank == 0) then
      SAFE_ALLOCATE(part_global(partition%np_global))
      call io_binary_read(filename, partition%np_global, part_global, ierr)
    else
      SAFE_ALLOCATE(part_global(1))
    end if
#ifdef HAVE_MPI
    ! All nodes need to know the result
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, partition%mpi_grp%comm, mpi_err)
#endif

    ! If reading was successfull, then scatter the result
    if (ierr == 0) then
      if (partition%mpi_grp%rank == 0) then
        SAFE_ALLOCATE(scounts(partition%npart))
        SAFE_ALLOCATE(sdispls(partition%npart))

        scounts(1:partition%npart-1) = partition%nppp
        scounts(partition%npart) = partition%np_global - partition%nppp*(partition%npart - 1)
        sdispls(1) = 0
        do ipart = 2, partition%npart
          sdispls(ipart) = sdispls(ipart-1) + scounts(ipart-1)
        end do
      else
        SAFE_ALLOCATE(scounts(1))
        SAFE_ALLOCATE(sdispls(1))
      end if

#ifdef HAVE_MPI
      call mpi_debug_in(partition%mpi_grp%comm, C_MPI_SCATTERV)
      call MPI_Scatterv(part_global(1), scounts(1), sdispls(1), MPI_INTEGER, &
                        partition%part(1), partition%np_local, MPI_INTEGER,  &
                        0, partition%mpi_grp%comm, mpi_err)
      call mpi_debug_out(partition%mpi_grp%comm, C_MPI_SCATTERV)
#endif

      SAFE_DEALLOCATE_A(scounts)
      SAFE_DEALLOCATE_A(sdispls)
    end if

    SAFE_DEALLOCATE_A(part_global)

    POP_SUB(partition_read)
  end subroutine partition_read

  ! ---------------------------------------------------------
  subroutine partition_get_local_size(partition, istart, np_local)
    type(partition_t), intent(in)  :: partition    
    integer,           intent(out) :: istart   !< The number of points of the partition stored in this process.
    integer,           intent(out) :: np_local !< The position of the first point stored in this process.

    PUSH_SUB(partition_get_local_size)

    istart = partition%istart
    np_local = partition%np_local

    POP_SUB(partition_get_local_size)
  end subroutine partition_get_local_size

  ! ---------------------------------------------------------
  !> Returns que global partition. If root is present, the partition is
  !! gathered only in that node. Otherwise it is gathered in all nodes.
  subroutine partition_get_global(partition, part_global, root)
    type(partition_t), intent(in)  :: partition
    integer,           intent(out) :: part_global(:)
    integer, optional, intent(in)  :: root

    integer :: ipart
    integer, allocatable :: rdispls(:), rcounts(:)

    PUSH_SUB(partition_get_global)

    SAFE_ALLOCATE(rdispls(partition%npart))
    SAFE_ALLOCATE(rcounts(partition%npart))

    rcounts(1:partition%npart-1) = partition%nppp
    rcounts(partition%npart) = partition%np_global - partition%nppp*(partition%npart - 1)
    rdispls(1) = 0
    do ipart = 2, partition%npart
      rdispls(ipart) = rdispls(ipart-1) + rcounts(ipart-1)
    end do

    if (present(root)) then
#ifdef HAVE_MPI
      call mpi_debug_in(partition%mpi_grp%comm, C_MPI_GATHERV)
      call MPI_Gatherv(partition%part(1), partition%np_local, MPI_INTEGER, &
           part_global(1), rcounts(1), rdispls(1), MPI_INTEGER, &
           root, partition%mpi_grp%comm, mpi_err)
#endif
    else
#ifdef HAVE_MPI
      call mpi_debug_in(partition%mpi_grp%comm, C_MPI_ALLGATHERV)
      call MPI_Allgatherv(partition%part(1), partition%np_local, MPI_INTEGER, &
           part_global(1), rcounts(1), rdispls(1), MPI_INTEGER, &
           partition%mpi_grp%comm, mpi_err)
      call mpi_debug_out(partition%mpi_grp%comm, C_MPI_GATHERV)
#endif
    end if

    SAFE_DEALLOCATE_A(rdispls)
    SAFE_DEALLOCATE_A(rcounts)

    POP_SUB(partition_get_global)
  end subroutine partition_get_global

  ! ---------------------------------------------------------
  !> Given a list of _global_ indexes, it returns the partition number
  !! were those points are stored.
  !! Note that this routine will accept global indexes equal to 0. In that
  !! case it will return 0 as a partition number.
  subroutine partition_get_partition_number(partition, np, points, partno)
    type(partition_t), intent(in)  :: partition
    integer,           intent(in)  :: np
    integer,           intent(in)  :: points(:)
    integer,           intent(out) :: partno(:)

    integer :: ip, nproc, rnp, zero_part
    integer, allocatable :: sbuffer(:), rbuffer(:)
    integer, allocatable :: scounts(:), rcounts(:)
    integer, allocatable :: sdispls(:), rdispls(:)
    integer, allocatable :: ipos(:), order(:)

    PUSH_SUB(partition_get_partition_number)

    SAFE_ALLOCATE(scounts(1:partition%npart))
    SAFE_ALLOCATE(rcounts(1:partition%npart))
    SAFE_ALLOCATE(sdispls(1:partition%npart))
    SAFE_ALLOCATE(rdispls(1:partition%npart))

    ! How many points we will have to send/receive from each process?
    scounts = 0
    zero_part = 1
    do ip = 1, np
      !Who knows where points(ip) is stored?
      if (points(ip) /= 0) then
        nproc = ceiling(real(points(ip))/real(partition%nppp))
      else
        nproc = zero_part
        zero_part = mod(zero_part, partition%npart) + 1
      end if

      !We increase the respective counter
      scounts(nproc) = scounts(nproc) + 1
    end do


    !Tell each process how many points we will need from it
#ifdef HAVE_MPI
    call mpi_debug_in(partition%mpi_grp%comm, C_MPI_ALLTOALL)
    call MPI_Alltoall(scounts(1), 1, MPI_INTEGER, &
                      rcounts(1), 1, MPI_INTEGER, &
                      partition%mpi_grp%comm, mpi_err)
    call mpi_debug_out(partition%mpi_grp%comm, C_MPI_ALLTOALL)
#endif

    !Build displacement arrays
    sdispls(1) = 0
    rdispls(1) = 0
    do ip = 2, partition%npart
      sdispls(ip) = sdispls(ip-1) + scounts(ip-1)
      rdispls(ip) = rdispls(ip-1) + rcounts(ip-1)
    end do


    rnp = sum(rcounts)
    SAFE_ALLOCATE(sbuffer(np))
    SAFE_ALLOCATE(rbuffer(rnp))

    !Put points in correct order for sending
    SAFE_ALLOCATE(ipos(1:partition%npart))
    SAFE_ALLOCATE(order(np))
    ipos = 0
    zero_part = 1
    do ip = 1, np
      !Who knows where points(ip) is stored?
      if (points(ip) /= 0) then
        nproc = ceiling(real(points(ip))/real(partition%nppp))
      else
        nproc = zero_part
        zero_part = mod(zero_part, partition%npart) + 1
      end if

      !We increase the respective counter
      ipos(nproc) = ipos(nproc) + 1

      !Put the point in the correct place
      order(ip) = sdispls(nproc) + ipos(nproc)
      sbuffer(sdispls(nproc) + ipos(nproc)) = points(ip) ! global index of the point
    end do
    SAFE_DEALLOCATE_A(ipos)

    !Send the global indexes of the points to the process that knows what is the corresponding partition
#ifdef HAVE_MPI
    call mpi_debug_in(partition%mpi_grp%comm, C_MPI_ALLTOALLV)
    call MPI_Alltoallv(sbuffer, scounts(1), sdispls(1), MPI_INTEGER, &
                       rbuffer, rcounts(1), rdispls(1), MPI_INTEGER, &
                       partition%mpi_grp%comm, mpi_err)
    call mpi_debug_out(partition%mpi_grp%comm, C_MPI_ALLTOALLV)
#endif

    !We get the partition number from the global index. This will be send back.
    do ip = 1, rnp
      if (rbuffer(ip) == 0) cycle
      rbuffer(ip) = partition%part(rbuffer(ip) - partition%istart + 1)
    end do

    !Now we send the information backwards
#ifdef HAVE_MPI
    call mpi_debug_in(partition%mpi_grp%comm, C_MPI_ALLTOALLV)
    call MPI_Alltoallv(rbuffer, rcounts(1), rdispls(1), MPI_INTEGER, &
                       sbuffer, scounts(1), sdispls(1), MPI_INTEGER,    &
                       partition%mpi_grp%comm, mpi_err)
    call mpi_debug_out(partition%mpi_grp%comm, C_MPI_ALLTOALLV)
#endif

    !Reorder the points
    do ip = 1, np
      partno(ip) = sbuffer(order(ip))
    end do

    !Deallocate memory
    SAFE_DEALLOCATE_A(order)
    SAFE_DEALLOCATE_A(sbuffer)
    SAFE_DEALLOCATE_A(scounts)
    SAFE_DEALLOCATE_A(sdispls)
    SAFE_DEALLOCATE_A(rbuffer)
    SAFE_DEALLOCATE_A(rcounts)
    SAFE_DEALLOCATE_A(rdispls)

    POP_SUB(partition_get_partition_number)
  end subroutine partition_get_partition_number

end module partition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
