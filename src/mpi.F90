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

! Notes ragarding the multi communicator part.
!
! This part is for combined parallelization in indices and
! domains.
! Given a number n_index of indices and their ranges
! index_range(1:n_index) communicators for domains and
! indices are created.
!
! Example
! Given to indices j, k with ranges
! index_range(j) = 3, index_range(k) = 4
! and n_node = 48 nodes, we would put
! n_node_domain = n_node/n_domain with
! n_domain = prod(a=j, k)[index_range(a)] = 4
! nodes into each domain parallelization.
! To do collective operations (like sums) over j, k respectively
! n_index_comm = N(n_index) = 7 communicators are introduced
! (for the definition of N see below).
! Each of these communicators contains only the root nodes of
! all domain parallelizations participating in this particular
! index. The reason for this is that only root nodes know
! the complete functions (after vec_gather).
! 
! The example above can visualized as follows:
! 
!   j
!  --->
!  |   (1)  (2)  (3)
! k|   (4)  (5)  (6)
!  |   (7)  (8)  (9)
!  V  (10) (11) (12)
!
! (n) are domain parallelizations of four nodes.
! The communicators for j are as follows:
! index_comm(k=1, j) = {1, 2, 3}
! index_comm(k=2, j) = {4, 5, 6}
! index_comm(k=3, j) = {7, 8, 9}
! index_comm(k=4, j) = {10, 11, 12}
! 
! For k they look like this:
! index_comm(j=1, k) = {1, 4, 7, 10}
! index_comm(j=2, k) = {2, 5, 8, 11}
! index_comm(j=3, k) = {3, 6, 9, 12}
!
! {p} means that the root node of domain parallelization
! p is member of the denoted communicator.
!
! In general the communicators for index a(i) out of
! a(1), ..., a(n_index) is addressed by specifying all
! other indices x, y, ... except i:
! index_comm(a=(x, y, ...), i)
!
! For generality, all index communicators are stored in a
! vector and the adressing is done by a function get_comm(a, i).
! a(1:n_index) specifies the values for all indices except i
! (the value of a(i) is irrelevant) and i is the number of
! the requested index:
! get_comm(a, i) == offset + position
! WHERE
! offset   == sum(x=1, ..., i-1)[prod(y=1, ..., n_index, y|=x)[index_range(y)]]
! position == sum(x=1, ..., n_index, x!=i)[(a(x)-1)*
!             prod(y=x+1, ..., n_index, y|=i)[index_range(y)]] + 1
!
! For each index i the rest if the indices form a
! n_index-1 dimensional array. The offset of this array in
! index_comm is offset in the above function. It is the number of
! communicators all indices 1, ..., i-1 have.
! The position in the array is computed as for any other
! n-dimensional array (for generalities sake it has to be done by hand
! and not by the Fortran compiler).
! With this function, the j communicator for k=2 can be accessed with
! index_comm(get_comm((/0, 2/), 1))
! (with j being the first and k the second index).
!
! Some more stripped down cases:
! (*) Only index parallelization (e. g. k-points and states j):
!     There are no domain communicators. In the sketch above
!     (1), ..., (12) would directly denote nodes and not domain
!     parallelizations
! (*) Only domain parallelization: There would not be any index
!     communicators and just one domain parallelization:
!     n_domain = 1, domain_comm(1) = MPI_COMM_WORLD.
!
!
! ---------- 
! (*) N(1)   = 1
!     N(i+1) = N(i)*index_range(i+1) + prod(a=1, ..., i)[index_range(a)]
!
! (*) prod(x=1, ..., n)[f(x)] = f(1)* ... *f(x)
!

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

  public ::                                                        &
       TSD_MPI_Barrier, TSZ_MPI_Barrier, TSI_MPI_Barrier,          &
       TSD_MPI_Scatterv, TSZ_MPI_Scatterv, TSI_MPI_Scatterv,       &
       TSD_MPI_Gatherv, TSZ_MPI_Gatherv, TSI_MPI_Gatherv,          &
       TSD_MPI_Alltoallv, TSZ_MPI_Alltoallv, TSI_MPI_Alltoallv,    &
       TSD_MPI_Allgatherv, TSZ_MPI_Allgatherv, TSI_MPI_Allgatherv, &
       TSD_MPI_Bcast, TSZ_MPI_Bcast, TSI_MPI_Bcast,                &
       TSD_MPI_Allreduce, TSZ_MPI_Allreduce, TSI_MPI_Allreduce

  integer, public, parameter ::  &
       C_MPI_BARRIER    = 1,     &
       C_MPI_SCATTERV   = 2,     &
       C_MPI_GATHERV    = 3,     &
       C_MPI_ALLTOALLV  = 4,     &
       C_MPI_ALLGATHERV = 5,     &
       C_MPI_BCAST      = 6,     &
       C_MPI_ALLREDUCE  = 7

  character(len=15), dimension(C_MPI_ALLREDUCE), public :: mpi_rlabel = &
       (/                &
       'MPI_BARRIER   ', &
       'MPI_SCATTERV  ', &
       'MPI_GATHERV   ', & 
       'MPI_ALLTOALLV ', &
       'MPI_ALLGATHERV', &
       'MPI_BCAST     ', &
       'MPI_ALLREDUCE '  &
       /)       

  integer, public :: call_counter(C_MPI_BARRIER:C_MPI_ALLREDUCE) = 0
  integer, public :: sec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)    = 0
  integer, public :: usec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)   = 0

  integer, private :: sec_in, usec_in

  ! Stores all communicators.
  type multicomm_type
    integer          :: n_node          ! Total number of nodes.
    integer          :: n_index         ! Number of parallel indices.
    integer          :: n_domain_comm   ! Number of domain communicators.
    integer          :: n_index_comm    ! Number of index communicators.

    integer, pointer :: index_range(:)  ! Range of index i is
                                        ! 1, ..., index_range(i).

    integer, pointer :: domain_comm(:)  ! par_vec communicators.
    integer, pointer :: index_comm(:)   ! Index communicators.

  end type multicomm_type


contains

  ! ---------------------------------------------------------
  subroutine MPI_Debug_Statistics()

    integer :: j
    integer :: usec_call(C_MPI_BARRIER:C_MPI_ALLREDUCE)

    if(.not.in_debug_mode) return

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

    if(.not.in_debug_mode) return

    call_counter(index) = call_counter(index) + 1
    call loct_gettimeofday(sec_in, usec_in)
    call epoch_time_diff(sec_in, usec_in)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i6.6,a,i4.4,a,i6.6)') '* I ',       &
         sec_in, '.', usec_in, ' '//trim(mpi_rlabel(index))//' - ', comm,':',   &
         call_counter(index), ' - ', sec_accum(index), '.', usec_accum(index)
    call write_debug(1)

  end subroutine MPI_Debug_IN


  ! ---------------------------------------------------------
  subroutine MPI_Debug_Out(comm, index)
    integer, intent(in) :: comm, index

    integer :: sec, usec, sec_diff, usec_diff

    if(.not.in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)
    call mpi_time_accum(index, sec, usec, sec_diff, usec_diff)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i6.6,a,i4.4,a,i6.6,a,i4.4,a,i6.6)') &
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

#include "undef.F90"
#include "integer.F90"
#include "mpi_inc.F90"

#endif
end module mpi_mod
