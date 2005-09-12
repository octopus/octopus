!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module mpi_mod
#if defined(HAVE_MPI)
  use varinfo
  use global
  use lib_oct
  use messages

  implicit none

  private

  public :: MPI_Debug_Statistics,               &
       MPI_Debug_IN, MPI_Debug_OUT  

  public :: TSD_MPI_Barrier, TSZ_MPI_Barrier,   &
       TSD_MPI_Scatterv, TSZ_MPI_Scatterv,      &
       TSD_MPI_Gatherv, TSZ_MPI_Gatherv,        &
       TSD_MPI_Alltoallv, TSZ_MPI_Alltoallv,    &
       TSD_MPI_Allreduce, TSZ_MPI_Allreduce

  integer, public, parameter ::  &
       C_MPI_BARRIER   = 1,      &
       C_MPI_SCATTERV  = 2,      &
       C_MPI_GATHERV   = 3,      &
       C_MPI_ALLTOALLV = 4,      &
       C_MPI_ALLREDUCE = 5

  character(len=15), dimension(C_MPI_ALLREDUCE), public :: mpi_rlabel = &
       (/                &
       'MPI_BARRIER  ',  &
       'MPI_SCATTERV ',  &
       'MPI_GATHERV  ',  & 
       'MPI_ALLTOALLV',  &
       'MPI_ALLREDUCE'   &
       /)       

  integer, public :: call_counter(C_MPI_BARRIER:C_MPI_ALLREDUCE) = 0
  integer, public :: sec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)    = 0
  integer, public :: usec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)   = 0

  integer, private :: sec_in, usec_in

contains

  ! ---------------------------------------------------------
  subroutine MPI_Debug_Statistics()

    integer :: j
    integer :: usec_call(C_MPI_BARRIER:C_MPI_ALLREDUCE)

    message(1) = ''
    message(2) = hyphens
    message(3) = ''
    write(message(4), '(23x,a,4x,a,8x,a)') 'total time', 'calls', 'usec/call'
    do j = 1, C_MPI_ALLREDUCE
       if (sec_accum(j).eq.0.and.usec_accum(j).eq.0) then
          usec_call(j) = 0
       else
          usec_call(j) = (sec_accum(j)*1000000+usec_accum(j))/call_counter(j)
       endif

       write(message(j+4),'(a,i6,a,i6.6,6x,i4,6x,i10)')          &
            mpi_rlabel(j)//' : ',                                &
            sec_accum(j), '.', usec_accum(j), call_counter(j),   &
            usec_call(j)
    enddo
    message(C_MPI_ALLREDUCE+5) = ''    
    message(C_MPI_ALLREDUCE+6) = hyphens    
    call write_debug(C_MPI_ALLREDUCE+6)

  end subroutine MPI_Debug_Statistics


  ! ---------------------------------------------------------
  subroutine MPI_Debug_In(comm, index)
    integer, intent(in) :: comm, index

    call_counter(index) = call_counter(index) + 1
    call loct_gettimeofday(sec_in, usec_in)
    call epoch_time_diff(sec_in, usec_in)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i4.4,a,i4.4,a,i6.6)') '* I ',       &
         sec_in, '.', usec_in, ' '//trim(mpi_rlabel(index))//' - ', comm,':',   &
         call_counter(index), ' - ', sec_accum(index), '.', usec_accum(index)
    call write_debug(1)

  end subroutine MPI_Debug_IN


  ! ---------------------------------------------------------
  subroutine MPI_Debug_Out(comm, index)
    integer, intent(in) :: comm, index

    integer :: sec, usec, sec_diff, usec_diff

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)
    call mpi_time_accum(index, sec, usec, sec_diff, usec_diff)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i4.4,a,i4.4,a,i6.6,a,i4.4,a,i6.6)') &
         '* O ',                                                                &
         sec, '.', usec, ' '//trim(mpi_rlabel(index))//' - ', comm, ':',        &
         call_counter(index), ' - ', sec_accum(index), '.', usec_accum(index),  &
         ' - ', sec_diff, '.', usec_diff
    call write_debug(1)

  end subroutine MPI_Debug_Out


  ! ---------------------------------------------------------
  subroutine mpi_time_accum(index, sec, usec, sec_diff, usec_diff)
    integer, intent(in)  :: index, sec, usec  
    integer, intent(out) :: sec_diff, usec_diff

    integer :: sec_tmp, usec_tmp

    sec_tmp  = sec
    usec_tmp = usec

    if (usec_tmp-usec_in .lt. 0) then
       usec_tmp = usec_tmp + 1000000
       sec_tmp  = sec_tmp  - 1
    endif
    usec_tmp = usec_tmp - usec_in
    sec_tmp  = sec_tmp  - sec_in    

    usec_diff = usec_tmp
    sec_diff  = sec_tmp

    ! accumulate values
    if (usec_tmp+usec_accum(index) .gt. 1000000) then
       usec_tmp = usec_tmp - 1000000 
       sec_tmp  = sec_tmp  + 1
    endif
    sec_accum(index)  = sec_accum(index)  + sec_tmp
    usec_accum(index) = usec_accum(index) + usec_tmp

  end subroutine mpi_time_accum


#include "undef.F90"
#include "real.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_inc.F90"

#endif
end module mpi_mod
