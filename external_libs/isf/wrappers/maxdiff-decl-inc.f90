!> @file
!! Include fortran file for maxdiff declarations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  integer, intent(in), optional :: root !<rank of the process retrieving the diff
  integer, intent(in), optional :: comm
  logical, intent(in), optional :: bcast !< if .true. all the proc will have the same maxdiff (default .false.)
  !local variables
  logical :: bcst
  integer :: ndims,nproc,mpi_comm,iroot,i,jproc

  bcst=.false.
  if (present(bcast)) bcst=bcast
  if (present(comm)) then
     mpi_comm=comm
  else
     mpi_comm=MPI_COMM_WORLD 
  end if
  if (present(root)) then
     iroot=root
  else
     iroot=0
  end if
  nproc=mpisize(mpi_comm)
