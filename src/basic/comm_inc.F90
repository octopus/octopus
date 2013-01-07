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

subroutine X(comm_allreduce_0)(comm, aa)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa

  R_TYPE :: aac

  !no PUSH SUB, called too often
    
#if defined(HAVE_MPI)
  aac = aa
  call MPI_Allreduce(aac, aa, 1, R_MPITYPE, MPI_SUM, comm, mpi_err)
#endif

end subroutine X(comm_allreduce_0)

! -----------------------------------------------------------------------------

subroutine X(comm_allreduce_1)(comm, aa, dim)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa(:)
  integer, optional,                intent(in)    :: dim

  integer :: dim1
#if defined(HAVE_MPI) && !defined(HAVE_MPI2)
  R_TYPE, allocatable :: aac(:)
#endif

  PUSH_SUB(X(comm_allreduce_1))

  dim1 = ubound(aa, dim = 1)
  if(present(dim)) dim1 = dim  

  ASSERT(ubound(aa, dim = 1) >= dim1)

  if(dim1 > 0) then

#if defined(HAVE_MPI2)
     call MPI_Allreduce(MPI_IN_PLACE, aa(1), dim1, R_MPITYPE, MPI_SUM, comm, mpi_err)
#elif defined(HAVE_MPI)
     
     SAFE_ALLOCATE(aac(1:dim1))
     aac(1:dim1) = aa(1:dim1)
     call MPI_Allreduce(aac(1), aa(1), dim1, R_MPITYPE, MPI_SUM, comm, mpi_err)
     SAFE_DEALLOCATE_A(aac)
     
#endif
  end if

  POP_SUB(X(comm_allreduce_1))
end subroutine X(comm_allreduce_1)

! -----------------------------------------------------------------------------

subroutine X(comm_allreduce_2)(comm, aa, dim)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa(:, :)
  integer, optional,                intent(in)    :: dim(1:2)

  integer :: dim_(1:2), ii
  R_TYPE, allocatable :: aac(:, :)
  
  PUSH_SUB(X(comm_allreduce_2))

  dim_ = ubound(aa)
  if(present(dim)) dim_ = dim  

  ASSERT(all(ubound(aa) >= dim_))

  if(any(dim_(1:2) < 1)) then
    POP_SUB(X(comm_allreduce_2))
    return
  end if

  if(ubound(aa, dim = 1) == dim_(1)) then
    ! the array is contiguous in memory

#if defined(HAVE_MPI2)

    call MPI_Allreduce(MPI_IN_PLACE, aa(1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)

#elif defined(HAVE_MPI)
    SAFE_ALLOCATE(aac(1:dim_(1), 1:dim_(2)))
    aac(1:dim_(1), 1:dim_(2)) = aa(1:dim_(1), 1:dim_(2))
    call MPI_Allreduce(aac(1, 1), aa(1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)
#endif

  else

#if defined(HAVE_MPI2)

    SAFE_ALLOCATE(aac(1:dim_(1), 1:dim_(2)))
    aac(1:dim_(1), 1:dim_(2)) = aa(1:dim_(1), 1:dim_(2))
    call MPI_Allreduce(MPI_IN_PLACE, aac(1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)
    aa(1:dim_(1), 1:dim_(2)) = aac(1:dim_(1), 1:dim_(2))
#elif defined(HAVE_MPI)
   
    SAFE_ALLOCATE(aac(1:dim_(1), 1))
    do ii = 1, dim_(2)
      aac(1:dim_(1), 1) = aa(1:dim_(1), ii)
      call MPI_Allreduce(aac(1, 1), aa(1, ii), dim_(1), R_MPITYPE, MPI_SUM, comm, mpi_err)
    end do
#endif
  end if

  SAFE_DEALLOCATE_A(aac)

  POP_SUB(X(comm_allreduce_2))
end subroutine X(comm_allreduce_2)

! -----------------------------------------------------------------------------

subroutine X(comm_allreduce_3)(comm, aa)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa(:, :, :)

  integer :: dim_(1:3)
  R_TYPE, allocatable :: aac(:, :, :)
  
  PUSH_SUB(X(comm_allreduce_3))

  dim_ = ubound(aa)

#if defined(HAVE_MPI2)

    call MPI_Allreduce(MPI_IN_PLACE, aa(1, 1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)

#elif defined(HAVE_MPI)
    SAFE_ALLOCATE(aac(1:dim_(1), 1:dim_(2), 1:dim_(3)))
    aac(1:dim_(1), 1:dim_(2), 1:dim_(3)) = aa(1:dim_(1), 1:dim_(2), 1:dim_(3))
    call MPI_Allreduce(aac(1, 1, 1), aa(1, 1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)
#endif

  SAFE_DEALLOCATE_A(aac)

  POP_SUB(X(comm_allreduce_3))
end subroutine X(comm_allreduce_3)


! -----------------------------------------------------------------------------

subroutine X(comm_allreduce_4)(comm, aa)
  integer,                          intent(in)    :: comm
  R_TYPE,                           intent(inout) :: aa(:, :, :, :)

  integer :: dim_(1:4)
  R_TYPE, allocatable :: aac(:, :, :, :)
  
  PUSH_SUB(X(comm_allreduce_4))

  dim_ = ubound(aa)

#if defined(HAVE_MPI2)

    call MPI_Allreduce(MPI_IN_PLACE, aa(1, 1, 1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)

#elif defined(HAVE_MPI)
    SAFE_ALLOCATE(aac(1:dim_(1), 1:dim_(2), 1:dim_(3), 1:dim_(4)))
    aac(1:dim_(1), 1:dim_(2), 1:dim_(3), 1:dim_(4)) = aa(1:dim_(1), 1:dim_(2), 1:dim_(3), 1:dim_(4))
    call MPI_Allreduce(aac(1, 1, 1, 1), aa(1, 1, 1, 1), product(dim_), R_MPITYPE, MPI_SUM, comm, mpi_err)
#endif

  SAFE_DEALLOCATE_A(aac)

  POP_SUB(X(comm_allreduce_4))
end subroutine X(comm_allreduce_4)


! -----------------------------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
