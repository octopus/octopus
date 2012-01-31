!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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

! Routines to support MPI debugging.

#include "global.h"

module mpi_debug_m
  use global_m
  use loct_m
  use messages_m
  use mpi_m

  implicit none

  private

#if !defined(HAVE_MPI)
  integer, public :: mpi_debug_dummy !< this avoids compilers complaining about empty module
#else
  public ::               &
    mpi_debug_statistics, &
    mpi_debug_in,         &
    mpi_debug_out

  public ::                      &
    TSD_MPI_Barrier,             &
    TSZ_MPI_Barrier,             &
    TSI_MPI_Barrier,             &
    TSD_MPI_Scatterv,            &
    TSZ_MPI_Scatterv,            &
    TSI_MPI_Scatterv,            &
    TSD_MPI_Gatherv,             &
    TSZ_MPI_Gatherv,             &
    TSI_MPI_Gatherv,             &
    TSD_MPI_Alltoallv,           &
    TSZ_MPI_Alltoallv,           &
    TSI_MPI_Alltoallv,           &
    TSD_MPI_Allgatherv,          &
    TSZ_MPI_Allgatherv,          &
    TSI_MPI_Allgatherv,          &
    TSD_MPI_Bcast,               &
    TSZ_MPI_Bcast,               &
    TSI_MPI_Bcast,               &
    TSD_MPI_Allreduce,           &
    TSZ_MPI_Allreduce,           &
    TSI_MPI_Allreduce,           &
    TSD_MPI_Alltoall,            &
    TSZ_MPI_Alltoall,            &
    TSI_MPI_Alltoall,            &
    TSD_MPI_Allgather,           &
    TSZ_MPI_Allgather,           &
    TSI_MPI_Allgather

  integer, parameter :: C_NUM_MPI_ROUTINES = 9

  integer, public, parameter ::  &
    C_MPI_BARRIER    = 1,        &
    C_MPI_SCATTERV   = 2,        &
    C_MPI_GATHERV    = 3,        &
    C_MPI_ALLTOALLV  = 4,        &
    C_MPI_ALLGATHERV = 5,        &
    C_MPI_BCAST      = 6,        &
    C_MPI_ALLREDUCE  = 7,        &
    C_MPI_ALLTOALL   = 8,        &
    C_MPI_ALLGATHER  = 9

  character(len=14), dimension(C_NUM_MPI_ROUTINES), public :: mpi_rlabel = &
    (/                           &
    'MPI_BARRIER   ',            &
    'MPI_SCATTERV  ',            &
    'MPI_GATHERV   ',            &
    'MPI_ALLTOALLV ',            &
    'MPI_ALLGATHERV',            &
    'MPI_BCAST     ',            &
    'MPI_ALLREDUCE ',            &
    'MPI_ALLTOALL  ',            &
    'MPI_ALLGATHER '             &
    /)

  integer, public :: call_counter(C_NUM_MPI_ROUTINES) = 0
  real(8), public :: sec_accum(C_NUM_MPI_ROUTINES)    = 0

  real(8) :: sec_in
#endif

contains

#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  subroutine mpi_debug_statistics()
    integer :: j
    real(8) :: usec_call(C_NUM_MPI_ROUTINES)

    if(.not.in_debug_mode) return

    message(1) = ''
    message(2) = hyphens
    message(3) = ''
    write(message(4), '(23x,a,4x,a,8x,a)') 'total time', 'calls', 'usec/call'
    do j = 1, C_NUM_MPI_ROUTINES
      if (sec_accum(j).eq.0) then
        usec_call(j) = 0
      else
        usec_call(j) = (sec_accum(j)*1000000)/call_counter(j)
      end if

      write(message(j+4),'(a,f13.6,6x,i4,6x,f13.0)') &
        mpi_rlabel(j)//' : ', sec_accum(j),          &
        call_counter(j), usec_call(j)
    end do
    message(C_NUM_MPI_ROUTINES+5) = ''
    message(C_NUM_MPI_ROUTINES+6) = hyphens
    call messages_debug(C_NUM_MPI_ROUTINES+6)
  end subroutine mpi_debug_statistics


  ! ---------------------------------------------------------
  subroutine mpi_debug_in(comm, index)
    integer, intent(in) :: comm, index

    if(.not.in_debug_mode) return

    call_counter(index) = call_counter(index) + 1
    sec_in              = MPI_Wtime()
    write(message(1),'(a,f18.6,a,z8.8,a,i6.6,a,f13.6)') '* MPI_I ', &
      sec_in, ' '//mpi_rlabel(index)//' : 0x', comm, ' | ',  &
      call_counter(index), ' - ', sec_accum(index)
    call messages_debug(1)
  end subroutine mpi_debug_in


  ! ---------------------------------------------------------
  subroutine mpi_debug_out(comm, index)
    integer, intent(in) :: comm, index

    real(8) :: sec_out, sec_diff

    if(.not.in_debug_mode) return

    sec_out = MPI_Wtime()
    call mpi_time_accum(index, sec_out, sec_diff)
    write(message(1),'(a,f18.6,a,z8.8,a,i6.6,a,f13.6,a,f13.6)')         &
      '* MPI_O ', sec_out, ' '//mpi_rlabel(index)//' : 0x', comm, ' | ', &
      call_counter(index), ' - ', sec_accum(index), ' - ', sec_diff
    call messages_debug(1)
  end subroutine mpi_debug_out


  ! ---------------------------------------------------------
  subroutine mpi_time_accum(index, sec, sec_diff)
    integer, intent(in)  :: index
    real(8), intent(in)  :: sec
    real(8), intent(out) :: sec_diff

    sec_diff         = sec - sec_in
    sec_accum(index) = sec_accum(index) + sec_diff
  end subroutine mpi_time_accum


#include "undef.F90"
#include "real.F90"
#include "mpi_debug_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_debug_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mpi_debug_inc.F90"

#else
  subroutine this_module_is_not_empty()
    integer :: neither_is_this_subroutine
    neither_is_this_subroutine = 0
  end subroutine this_module_is_not_empty
#endif
end module mpi_debug_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
