!! Copyright (C) 2008-2016 X. Andrade
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

subroutine X(distributed_allgather)(this, aa)
  type(distributed_t), intent(in)    :: this
  R_TYPE,              intent(inout) :: aa(:)
  
  integer, allocatable :: displs(:)
  
  if(.not. this%parallel) return
  
  PUSH_SUB(X(distributed_allgather))
  
  SAFE_ALLOCATE(displs(0:this%mpi_grp%size - 1))
  
  displs(0:this%mpi_grp%size - 1) = this%range(1, 0:this%mpi_grp%size - 1) - 1
  
#ifdef HAVE_MPI    
  call MPI_Allgatherv(MPI_IN_PLACE, this%nlocal, R_MPITYPE, aa(1), &
    this%num(0), displs(0), R_MPITYPE, this%mpi_grp%comm, mpi_err)
#endif
  
  SAFE_DEALLOCATE_A(displs)
  
  POP_SUB(X(distributed_allgather))
end subroutine X(distributed_allgather)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
