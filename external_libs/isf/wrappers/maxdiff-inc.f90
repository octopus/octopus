!> @file
!! Include fortran file for maxdiff operations with scalars
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (nproc == 1 .or. ndims == 0) return

  !check that the positions are identical for all the processes
  array_glob=f_malloc((/ndims,nproc/),id='array_glob')

  call mpigather(sendbuf=array,sendcount=ndims,recvbuf=array_glob,&
       root=iroot,comm=mpi_comm)

  include 'maxdiff-end-inc.f90'
