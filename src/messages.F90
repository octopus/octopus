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

module messages
  use global
  use lib_oct

  implicit none

  private

  public :: write_fatal, write_warning, write_info, print_date
  public :: input_error
  public :: push_sub, pop_sub

  character(len=256), dimension(20), public :: message    ! to be output by fatal, warning
  character(len=70),      parameter, public :: stars =  &
       '**********************************************************************'
  character(len=68),      parameter, public :: hyphens = &
       '--------------------------------------------------------------------'
  character(len=69),      parameter, public :: shyphens = '*'//hyphens

  ! min_lun in io.F90 is equal to 10. We hardwire this here since we 
  ! cannot write "use io" above. Unit 8 and 9 should always be available.
  integer, parameter, private :: iunit_out = 8
  integer, parameter, private :: iunit_err = 9
  character(len=512), private :: msg

contains


  ! ---------------------------------------------------------
  subroutine write_fatal(no_lines)
    integer, intent(in) :: no_lines

    integer :: i
#ifdef HAVE_MPI
    integer :: ierr
#endif

    if(conf%flush_messages.and.mpiv%node.eq.0) then
       open(unit=iunit_err, file='messages.stderr', &
            action='write', position='append')
    endif

    call flush_msg(stderr, '')
    call flush_msg(stderr, stars)
    call flush_msg(stderr, '')
    write(msg, '(a)') '*** Fatal Error (description follows)'
    call flush_msg(stderr, msg)

#ifdef HAVE_MPI
    call flush_msg(stderr, shyphens)
    write(msg, '(a,i4)') "* From node = ", mpiv%node
    call flush_msg(stderr, msg)
#endif
    call flush_msg(stderr, shyphens)
    do i=1,no_lines
       write(msg, '(a,1x,a)') '*', trim(message(i))
       call flush_msg(stderr, msg)
    end do
    call flush_msg(stderr, shyphens)

#ifdef DEBUG
    write(msg, '(a)') '* Stack: '
    call flush_msg(stderr, msg, 'no')
    do i = 1, no_sub_stack
       write(msg, '(a,a)') ' > ', trim(sub_stack(i))
       call flush_msg(stderr, msg, 'no')
    end do
    call flush_msg(stderr, '')
    call flush_msg(stderr, stars)
    call flush_msg(stderr, '')
    ! cannot call this anymore since the move from this routine from global.F90 
    ! to messages.F90
    !    call io_status(stderr)
#endif

    if(conf%flush_messages.and.mpiv%node.eq.0) then
       close(iunit_err)
    endif

#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif

    stop
  end subroutine write_fatal


  ! ---------------------------------------------------------
  subroutine write_warning(no_lines)
    integer, intent(in) :: no_lines
    integer :: i

    if(conf%flush_messages.and.mpiv%node.eq.0) then
       open(unit=iunit_err, file='messages.stderr', &
            action='write', position='append')
    endif

    ! this always writes from ALL nodes

    if(conf%verbose >= VERBOSE_WARNING) then
       call flush_msg(stderr, '')
       write(msg, '(a)') '** Warning:'
       call flush_msg(stderr, msg)
#ifdef HAVE_MPI
       write(msg , '(a,i4)') '** From node = ', mpiv%node
       call flush_msg(stderr, msg)
#endif
       do i=1,no_lines
          write(msg , '(a,3x,a)') '**', trim(message(i))
          call flush_msg(stderr, msg)
       end do
    end if
#ifdef HAVE_FLUSH
    call flush(stderr)
#endif 

    if(conf%flush_messages.and.mpiv%node.eq.0) then
       close(iunit_err)
    endif

    return
  end subroutine write_warning


  ! ---------------------------------------------------------
  subroutine write_info(no_lines, iunit, verbose_limit, stress)
    integer, intent(in) :: no_lines
    integer, intent(in), optional :: iunit
    integer, intent(in), optional :: verbose_limit
    logical, optional, intent(in) :: stress

    integer :: i, iu

#ifdef HAVE_MPI
    if(mpiv%node .ne. 0) return
#endif

    if(conf%flush_messages) then
       open(unit=iunit_out, file='messages.stdout', &
            action='write', position='append')
    endif

    if(present(iunit)) then
       iu = iunit
    else
       iu = stdout
    end if

    if(conf%verbose >= VERBOSE_NORMAL) then
       if(present(stress)) then
          call flush_msg(iu, stars)
       endif
       do i = 1, no_lines
          if(.not.present(verbose_limit)) then
             write(msg, '(a)') trim(message(i))
             call flush_msg(iu, msg)
          else if(conf%verbose>verbose_limit) then
             write(msg, '(a)') trim(message(i))
             call flush_msg(iu, msg)
          endif
       enddo
       if(present(stress)) then
          call flush_msg(iu, stars)
          call flush_msg(iu, '')
       endif
    end if

    if(conf%flush_messages) close(iunit_out)

#ifdef HAVE_FLUSH
    call flush(iu)
#endif
  end subroutine write_info


  ! ---------------------------------------------------------
  subroutine input_error(var)
    character(len=*), intent(in) :: var

#ifdef HAVE_MPI
    integer :: ierr
#endif

    if(mpiv%node == 0) then
       call flush_msg(stderr, '')
       call flush_msg(stderr, stars)
       call flush_msg(stderr, '')
       write(msg, '(a)') '*** Fatal Error in input '
       call flush_msg(stderr, msg)
       call flush_msg(stderr, shyphens)

       call varinfo_print(stderr, var)
       call flush_msg(stderr, shyphens)
    end if

#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif
    stop
  end subroutine input_error


  ! ---------------------------------------------------------
  subroutine flush_msg(iunit, str, adv)
    integer,            intent(in)           :: iunit 
    character(len = *), intent(in)           :: str
    character(len = *), intent(in), optional :: adv

    character(len = 20) :: adv_

    if(present(adv)) then
       adv_ = adv
    else
       adv_ = 'yes'
    endif

    write(iunit, '(a)', advance=trim(adv_)) trim(str)
    if(conf%flush_messages.and.mpiv%node.eq.0) then
       if(iunit.eq.stderr) write(iunit_err, '(a)', advance=trim(adv_)) trim(str)
       if(iunit.eq.stdout) write(iunit_out, '(a)', advance=trim(adv_)) trim(str)
    endif

  end subroutine flush_msg


  ! ---------------------------------------------------------
  subroutine print_date(str)
    character(len = *), intent(in) :: str
    integer :: val(8)
    call date_and_time(values=val)
    message(1) = ""
    write(message(3),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
         str , val(1), "/", val(2), "/", val(3), &
         " at ", val(5), ":", val(6), ":", val(7)
    message(2) = str_center(trim(message(3)), 70)
    message(3) = ""
    call write_info(3)
  end subroutine print_date


  ! ---------------------------------------------------------
#ifdef DEBUG
  subroutine push_sub(sub_name)
    character(len=*), intent(in) :: sub_name
    integer i

    no_sub_stack = no_sub_stack + 1
    if(no_sub_stack > 49) then
       sub_stack(50) = 'push_sub'
       message(1) = 'Too many recursion levels (max=50)'
       call write_fatal(1)
    else
       sub_stack(no_sub_stack) = trim(sub_name)
       time_stack(no_sub_stack) = loct_clock()

       if(conf%verbose >= VERBOSE_DEBUG .and. no_sub_stack <= conf%debug_level .and. mpiv%node == 0) then
          write(stderr,'(a,f10.3,i10, a)', advance='no') "* I ", loct_clock()/CNST(1e6), &
               loct_getmem(), " | "
          do i = no_sub_stack-1, 1, -1
             write(stderr,'(a)', advance='no') "..|"
          end do
          write(stderr,'(a)') trim(sub_name)
       end if
    end if

    return
  end subroutine push_sub


  ! ---------------------------------------------------------
  subroutine pop_sub()
    integer i

    if(no_sub_stack > 0) then
       if(conf%verbose > VERBOSE_DEBUG .and. no_sub_stack <= conf%debug_level .and. mpiv%node == 0) then

          ! It seems in std C libraries the number of clock ticks per second is 1e6...
          write(stderr,'(a,f10.3,i10, a)', advance='no') "* O ", &
               (loct_clock()-time_stack(no_sub_stack))/CNST(1e6), &
               loct_getmem(), " | "
          do i = no_sub_stack-1, 1, -1
             write(stderr,'(a)', advance='no') "..|"
          end do
          write(stderr,'(a)') trim(sub_stack(no_sub_stack))
       end if
       no_sub_stack = no_sub_stack - 1
    else
       no_sub_stack = 1
       sub_stack(1) = 'pop_sub'
       message(1) = 'Too few recursion levels'
       call write_fatal(1)    
    end if

  end subroutine pop_sub
#endif

end module messages
