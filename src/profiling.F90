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

! Simple profiling module to count the number of
! times code between
! call profiling_in(tag)
! ...
! call profiling_out(tag)
! is executed and how much time is consumed. The tag
! has to be registered as constant C_PROFILING_[SOMETHING] and
! in tag_label (C_TAG_LENGTH and C_NUM_TAG should be adjusted
! after adding a new tag).
! The routines maintain entry and exit counters (in profiling_in,
! profiling_out respectively). They should be equal, otherwise
! measurements are likely to be incorrect (e. g. because
! of a return inside subroutine without an profiling_out).
! loct_gettimeofday is used to get the times.
! The results are written to .profiling.NNN/profiling.nnn
! for every node with nnn being the node number and NNN the
! total number of nodes.

module profiling_mod

  use io
  use messages
  use global
  use mpi_mod
  use lib_oct

  implicit none
  private

  public ::                            &
    profiling_init,                    &
    profiling_in,                      &
    profiling_out,                     &
    profiling_output

  integer, parameter ::                &
    C_TAG_LENGTH = 17,                 & ! Max. number of characters of tag label.
    C_NUM_TAGS   = 9                     ! Number of tags.

  integer ::                           &
    pass_in(C_NUM_TAGS)           = 0, &
    pass_out(C_NUM_TAGS)          = 0, &
    time_in_sec(C_NUM_TAGS)       = 0, &
    time_in_usec(C_NUM_TAGS)      = 0, &
    time_sec(C_NUM_TAGS)          = 0, &
    time_usec(C_NUM_TAGS)         = 0

  integer, parameter, public ::        &
    C_PROFILING_COMPLETE_SUBSYS   = 1, &
    C_PROFILING_MF_INTEGRATE      = 2, &
    C_PROFILING_MF_DOTP           = 3, &
    C_PROFILING_MF_NRM2           = 4, &
    C_PROFILING_NL_OPERATOR       = 5, &
    C_PROFILING_GHOST_UPDATE      = 6, &
    C_PROFILING_VEC_INTEGRATE     = 7, &
    C_PROFILING_SCF_CYCLE         = 8, &
    C_PROFILING_MF_DOTP_ALLREDUCE = 9

  character(len=C_TAG_LENGTH), dimension(C_NUM_TAGS) :: tag_label = &
    (/                                 &
    'COMPLETE_SUBSYS  ',               &
    'MF_INTEGRATE     ',               &
    'MF_DOTP          ',               &
    'MF_NRM2          ',               &
    'NL_OPERATOR      ',               &
    'GHOST_UPDATE     ',               &
    'VEC_INTEGRATE    ',               &
    'SCF_CYCLE        ',               &
    'MF_DOTP_ALLREDUCE'                &
    /)

contains

  ! Create profiling subdirectory.
  subroutine profiling_init
#if defined(HAVE_MPI)
    character(len=3) :: dirnum
#endif

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_init')

#if defined(HAVE_MPI)
    if(mpi_grp_is_root(mpi_world)) then
      write(dirnum, '(i3.3)') mpi_world%size
      call io_mkdir(trim('profiling.'//dirnum))
    end if
#else
    call io_mkdir('profiling.001')
#endif

    ! initialize counter
    pass_in      = 0
    pass_out     = 0
    time_in_sec  = 0
    time_in_usec = 0
    time_sec     = 0
    time_usec    = 0

    call pop_sub()

  end subroutine profiling_init


  ! Increment in counter and save entry time.
  subroutine profiling_in(tag)
    integer :: tag

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_in')

    pass_in(tag) = pass_in(tag)+1
    call loct_gettimeofday(time_in_sec(tag), time_in_usec(tag))

    call pop_sub()

  end subroutine profiling_in


  ! Increment out counter and sum up difference between entry
  ! and exit time.
  subroutine profiling_out(tag)
    integer :: tag

    integer :: sec, usec

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_out')

    pass_out(tag) = pass_out(tag)+1
    call loct_gettimeofday(sec, usec)
    call time_diff(time_in_sec(tag), time_in_usec(tag), sec, usec)
    call time_sum(sec, usec, time_sec(tag), time_usec(tag))

    call pop_sub()

  end subroutine profiling_out


  ! Write profiling results of each node to profiling.NNN/profifling.nnn
  ! The format of each line is
  ! tag-label    pass_in    pass_out    time   time/pass_in
  !
  ! The last column gives the average time consumed between in and out
  ! (only, if pass_in and pass_out are equal).
  subroutine profiling_output
    integer          :: i
    integer          :: iunit
    character(len=3) :: filenum
    character(len=3) :: dirnum
    real             :: time_per_pass

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_output')

#if defined(HAVE_MPI)
    write(filenum, '(i3.3)') mpi_world%rank
    write(dirnum, '(i3.3)') mpi_world%size
#else
    filenum = '000'
    dirnum  = '001'
#endif

    iunit = io_open('profiling.'//dirnum//'/profiling.'//filenum, action='write')
    if(iunit.lt.0) then
      message(1) = 'Could not write profiling results.'
      call write_warning(1)
      call pop_sub()
      return
    end if
    write(iunit,'(a71)') &
      'TAG                   NUMBER OF CALLS       TOTAL TIME    TIME PER CALL'
    write(iunit,'(a71)') &
      '======================================================================='
    do i = 1, C_NUM_TAGS
      if(pass_in(i).eq.pass_out(i).and.pass_in(i).ne.0) then
        time_per_pass = (real(time_sec(i))+1e-6*real(time_usec(i)))/ &
          real(pass_in(i))
      else
        write(iunit,'(3a,i10,a,i10)') '*** WARNING: Entries and exits for ',trim(tag_label(i)), ' do not coincide:', &
                          pass_in(i), ' .ne.', pass_out(i)
        cycle
      end if
      write(iunit, '(a,i20,i10,a1,i6.6,f17.7)') tag_label(i), pass_in(i), &
        time_sec(i), '.', time_usec(i), time_per_pass
    end do

    call io_close(iunit)
    call pop_sub()
  end subroutine profiling_output

end module profiling_mod
