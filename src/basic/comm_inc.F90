!! Copyright (C) 2010 X. Andrade
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
!! $Id: comm_inc.F90 3587 2007-11-22 16:43:00Z xavier $

! -----------------------------------------------------------------------------

subroutine X(comm_allreduce_1)(comm, aa, dim)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa(:)
  integer, optional,                intent(in)    :: dim

  integer :: dim1
#if defined(HAVE_MPI) && !defined(HAVE_MPI2)
  R_TYPE, allocatable :: aac(:)
#endif

  dim1 = ubound(aa, dim = 1)
  if(present(dim)) dim1 = dim  

  ASSERT(ubound(aa, dim = 1) >= dim1)
  
#if defined(HAVE_MPI2)

  call MPI_Allreduce(MPI_IN_PLACE, aa(1), dim1, R_MPITYPE, MPI_SUM, comm, mpi_err)

#elif defined(HAVE_MPI)

  SAFE_ALLOCATE(aac(1:dim1))
  aac(1:dim1) = aa(1:dim1)
  call MPI_Allreduce(aac(1), aa(1), dim1, R_MPITYPE, MPI_SUM, comm, mpi_err)
  SAFE_DEALLOCATE_A(aac)

#endif

end subroutine X(comm_allreduce_1)

! -----------------------------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
