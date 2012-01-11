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

subroutine X(partition_transfer)(this, f_in, f_out)
  type(partition_transfer_t), intent(in)  :: this
  R_TYPE,                     intent(in)  :: f_in(:)
  R_TYPE,                     intent(out) :: f_out(:)
    
  PUSH_SUB(X(partition_transfer))  

  call profiling_in(prof_transfer, "PARTITION_TRANSFER")

#ifdef HAVE_MPI
  call mpi_debug_in(this%comm, C_MPI_ALLTOALLV)
  call MPI_Alltoallv(f_in, this%scounts(1), this%sdispls(1), R_MPITYPE, &
       f_out, this%rcounts(1), this%rdispls(1), R_MPITYPE,    &
       this%comm, mpi_err)
  call mpi_debug_out(this%comm, C_MPI_ALLTOALLV)
#endif

  call profiling_out(prof_transfer)

  POP_SUB(X(partition_transfer))
end subroutine X(partition_transfer)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
