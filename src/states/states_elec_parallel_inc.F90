!! Copyright (C) 2015 X. Andrade
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

subroutine X(states_elec_parallel_gather_3)(st, dims, psi)
  type(states_elec_t), intent(in)  :: st
  integer,           intent(in)    :: dims(2)
  R_TYPE,            intent(inout) :: psi(:, :, :)

  integer :: maxst, ist, i1, i2, irank, ist_local
  R_TYPE, allocatable :: sendpsi(:, :, :), recvpsi(:, :, :)
  
  !no PUSH_SUB, called too often
  
  call profiling_in(prof_gather, 'STATES_GATHER')

  if(st%parallel_in_states) then

    maxst = maxval(st%dist%num(0:st%mpi_grp%size - 1))
    
    SAFE_ALLOCATE(sendpsi(1:dims(1), 1:dims(2), 1:maxst))
    SAFE_ALLOCATE(recvpsi(1:dims(1), 1:dims(2), 1:maxst*st%mpi_grp%size))

    ! We have to use a temporary array to make the data contiguous
    
    do ist = 1, st%lnst
      do i1 = 1, dims(1)
        do i2 = 1, dims(2)
          sendpsi(i1, i2, ist) = psi(st%st_start + ist - 1, i1, i2)
        end do
      end do
    end do

#ifdef HAVE_MPI
    call MPI_Allgather(sendpsi(1, 1, 1), product(dims(1:2))*maxst, R_MPITYPE, &
      recvpsi(1, 1, 1), product(dims(1:2))*maxst, R_MPITYPE, st%mpi_grp%comm, mpi_err)
#endif

    ! now get the correct states from the data of each rank
    ist = 0
    do irank = 0, st%mpi_grp%size - 1
      do ist_local = 1, st%dist%num(irank)
        ist = ist + 1
        do i1 = 1, dims(1)
          do i2 = 1, dims(2)
            psi(ist, i1, i2) = recvpsi(i1, i2, irank*maxst + ist_local)
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(sendpsi)
    SAFE_DEALLOCATE_A(recvpsi)
    
  end if

  call profiling_out(prof_gather)
end subroutine X(states_elec_parallel_gather_3)

!---------------------------------------------------

subroutine X(states_elec_parallel_gather_1)(st, aa)
  type(states_elec_t), intent(in)    :: st
  R_TYPE,              intent(inout) :: aa(:)

  !no PUSH_SUB, called too often

  R_TYPE, allocatable :: sendaa(:)
  integer, allocatable :: displs(:)
  
  call profiling_in(prof_gather, 'STATES_GATHER')

  if(st%parallel_in_states) then

    SAFE_ALLOCATE(sendaa(st%st_start:st%st_end))
    SAFE_ALLOCATE(displs(0:st%mpi_grp%size - 1))
    
    sendaa(st%st_start:st%st_end) = aa(st%st_start:st%st_end)
    displs(0:st%mpi_grp%size - 1) = st%dist%range(1, 0:st%mpi_grp%size - 1) - 1
    
#ifdef HAVE_MPI
    call MPI_Allgatherv(sendaa(st%st_start), st%lnst, R_MPITYPE, &
      aa(1), st%dist%num(0), displs(0), R_MPITYPE, st%mpi_grp%comm, mpi_err)
#endif

    SAFE_DEALLOCATE_A(sendaa)
    SAFE_DEALLOCATE_A(displs)
    
  end if

  call profiling_out(prof_gather)
end subroutine X(states_elec_parallel_gather_1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
