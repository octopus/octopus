!! Copyright (C) 2011 M. Oliveira
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

#include "global.h"

module partition_transfer_m
  use global_m
  use iihash_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use par_vec_m
  use profiling_m
  use subarray_m

  implicit none

  private
  public ::                      &
    partition_transfer_t,        &
    partition_transfer_init,     &
    partition_transfer_end,      &
    dpartition_transfer,         &
    zpartition_transfer

  type partition_transfer_t
    private
    integer :: comm

    integer, pointer :: rdispls(:)
    integer, pointer :: sdispls(:)
    integer, pointer :: rcounts(:)
    integer, pointer :: scounts(:)
  end type partition_transfer_t

  type(profile_t), save :: prof_transfer

contains

  ! -----------------------------------------------------------------
  !> \warning  Input and output groups may have a different number of
  !! processes. In that case the transfer will only work if some
  !! further constraints are met (see sanity checks below). One of the
  !! cases where it should work is when one of the groups is a subgroup
  !! of a Cartesian topology created from the other group. This is the
  !! case when one of the groups is mpi_world, the other is the
  !! parallelization in domains group, and there are no slaves.
  subroutine partition_transfer_init(this, np, mpi_grp_in, part_in, mpi_grp_out, part_out, nsend, nrec, order_in, order_out)
    type(partition_transfer_t), intent(out) :: this
    integer,                    intent(in)  :: np 
    type(mpi_grp_t), target,    intent(in)  :: mpi_grp_in
    integer,                    intent(in)  :: part_in(:)  !< point -> partition
    type(mpi_grp_t), target,    intent(in)  :: mpi_grp_out
    integer,                    intent(in)  :: part_out(:) !< point -> partition
    integer,                    intent(out) :: nsend
    integer,                    intent(out) :: nrec
    integer, pointer,           intent(out) :: order_in(:)
    integer, pointer,           intent(out) :: order_out(:)

    logical :: found
    integer :: n12, tmp_partno(2), ipart, opart, ip, pcount, mycolumn, irec, isend, ipos
    type(iihash_t) :: map_out, map_in
    type(mpi_grp_t), pointer :: grp1, grp2
    integer, allocatable :: partno_list(:,:), part_map(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(partition_transfer_init)
   call profiling_in(prof,"P_TRANS_INIT")
   
    ! In order to avoid unnecessary communications, all the data
    ! transfer is going to be made from the point of view of the group
    ! that has more processes.
    if (mpi_grp_in%size >= mpi_grp_out%size) then
      grp1 => mpi_grp_in
      grp2 => mpi_grp_out
    else
      grp1 => mpi_grp_out
      grp2 => mpi_grp_in
    end if
    this%comm = grp1%comm

    ! The number of partitions in group 1 should be a multiple of the
    ! number of partitions in group 2.
    if (mod(grp1%size, grp2%size) /= 0) then
      message(1) = "Incompatible size of mpi groups in partition_transfer_init"
      call messages_fatal(1)
    end if
    n12 = grp1%size/grp2%size

    ! We need to know the partition number of all the processes in
    ! both groups
    SAFE_ALLOCATE(partno_list(1:2, 1:grp1%size))
    tmp_partno(1) = grp1%rank + 1
    tmp_partno(2) = grp2%rank + 1
#ifdef HAVE_MPI
    call MPI_Allgather(tmp_partno(1), 2, MPI_INTEGER, partno_list(1, 1), 2, MPI_INTEGER, &
         grp1%comm, mpi_err)
#endif

    ! Build partition map. This is a matrix with n12 columns and
    ! grp2%size lines. Each line contains the partition numbers of the
    ! processes of group 1 that also store the partition of group 2
    ! with the same number as the line number. The number of columns
    ! for each line should be exactly equal to grp1%size/grp2%size in
    ! all cases and there should be no repeated values.
    SAFE_ALLOCATE(part_map(1:grp2%size, 1:n12))
    part_map = 0
    do ipart = 1, grp2%size
      pcount = 0
      do ip = 1, grp1%size        
        if (partno_list(2, ip) == ipart) then
          pcount = pcount + 1
          
          if (pcount > n12 .or. any(partno_list(1, ip) == part_map(1:ipart,:))) then
            message(1) = "Incompatible mpi groups in partition_transfer_init"
            call messages_fatal(1)
          end if
          part_map(ipart, pcount) = partno_list(1, ip)
          if (ip == grp1%rank + 1) mycolumn = pcount
        end if
      end do
      if (pcount /= n12) then
        message(1) = "Incompatible mpi groups in partition_transfer_init"
        call messages_fatal(1)
      end if
    end do

    ! Build mapping between all the possible receivers and the ouput
    ! group. This map is a hash table, where the keys are the possible
    ! receivers and the values are the output partition these
    ! receivers are responsible for.
    ! If group 1 is the input group, then all members of group 1 are
    ! possible receivers.
    ! If group 1 is the output group, then, in order to avoid
    ! unnecessary communications, each process will only send data to
    ! a subset of all the possible receivers.  We will choose the
    ! processes that are on the same column of the partition map than
    ! the local process. Note that this implies that there are always
    ! mpi_grp_in%size possible receivers.
    call iihash_init(map_out, mpi_grp_in%size)
    do ipart = 1, grp2%size
      if (mpi_grp_in%size >= mpi_grp_out%size) then
        do ip = 1, n12
          call iihash_insert(map_out, part_map(ipart, ip), ipart)
        end do
      else
        call iihash_insert(map_out, part_map(ipart, mycolumn), part_map(ipart, mycolumn))
      end if
    end do

    ! Build mapping between all the possible senders and the input
    ! group. This map is a hash table similar to the map_out one.
    call iihash_init(map_in, mpi_grp_in%size)
    do ipart = 1, grp2%size
      if (mpi_grp_in%size >= mpi_grp_out%size) then
        do ip = 1, n12
          call iihash_insert(map_in, part_map(ipart, ip), part_map(ipart, ip))
        end do
      else
        call iihash_insert(map_in, part_map(ipart, mycolumn), ipart)
      end if
    end do

    ! Total number of points to be sent
    nsend = 0
    nsend = nsend + count(part_in(1:np) == mpi_grp_in%rank + 1)
  
    ! List of points to be send
    SAFE_ALLOCATE(this%sdispls(1:grp1%size))
    SAFE_ALLOCATE(this%scounts(1:grp1%size))
    SAFE_ALLOCATE(order_in(1:nsend))

    ipos = 0
    ! Loop over all possible receivers
    do irec = 1, grp1%size
      this%scounts(irec) = 0
      this%sdispls(irec) = ipos

      opart = iihash_lookup(map_out, irec, found)
      if (.not. found) cycle

      do ip = 1, np
        ! Does point ip belong to partition ipart and should it be sent to partition opart?
        if (part_in(ip) == mpi_grp_in%rank + 1 .and. part_out(ip) == opart) then
          ipos = ipos + 1
          order_in(ipos) = ip
          this%scounts(irec) = this%scounts(irec) + 1
        end if
      end do

    end do

    ! Total number of points to be received
    nrec = 0
    nrec = nrec + count(part_out(1:np) == mpi_grp_out%rank + 1)
    ! Displacements and number of points to be received 
    SAFE_ALLOCATE(this%rdispls(1:grp1%size))
    SAFE_ALLOCATE(this%rcounts(1:grp1%size))
    SAFE_ALLOCATE(order_out(1:nrec))

    ipos = 0
    ! Loop over all possible senders
    do isend = 1, grp1%size
      this%rcounts(isend) = 0
      this%rdispls(isend) = ipos

      ipart = iihash_lookup(map_in, isend, found)
      if (.not. found) cycle

      do ip = 1, np
        ! Is point ip to be received by this output partition and is it stored in partition ipart?
        if (part_in(ip) == ipart .and. part_out(ip) == mpi_grp_out%rank + 1) then
          ipos = ipos + 1
          order_out(ipos) = ip
          this%rcounts(isend) = this%rcounts(isend) + 1
        end if
      end do

    end do

    SAFE_DEALLOCATE_A(part_map)
    SAFE_DEALLOCATE_A(partno_list)
    call iihash_end(map_out)
    call iihash_end(map_in)

    call profiling_out(prof)
    POP_SUB(partition_transfer_init)
  end subroutine partition_transfer_init

  ! -----------------------------------------------------------------
  subroutine partition_transfer_end(this)
    type(partition_transfer_t), intent(inout) :: this

    PUSH_SUB(partition_transfer_end)

    SAFE_DEALLOCATE_P(this%rdispls)
    SAFE_DEALLOCATE_P(this%sdispls)
    SAFE_DEALLOCATE_P(this%rcounts)
    SAFE_DEALLOCATE_P(this%scounts)

    POP_SUB(partition_transfer_end)
  end subroutine partition_transfer_end


#include "undef.F90"
#include "real.F90"
#include "partition_transfer_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "partition_transfer_inc.F90"
  
end module partition_transfer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
