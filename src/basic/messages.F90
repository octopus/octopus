!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module messages_m
  use datasets_m
  use global_m
  use loct_m
  use loct_parser_m
  use mpi_m
  use string_m
  use varinfo_m

  implicit none

  private

  public ::                   &
    write_fatal,              &
    write_warning,            &
    write_info,               &
    write_debug,              &
    write_debug_marker,       &
    write_debug_newlines,     &
    print_date,               &
    time_diff,                &
    time_sum,                 &
    epoch_time_diff,          &
    alloc_error,              &
    input_error,              &
    push_sub,                 &
    pop_sub,                  &
    messages_print_stress,    &
    messages_print_var_info,  &
    messages_print_var_option,&
    messages_print_var_value, &
    obsolete_variable,        &
    messages_devel_version

  character(len=256), dimension(20), public :: message    ! to be output by fatal, warning
  character(len=68),      parameter, public :: hyphens = &
    '--------------------------------------------------------------------'
  character(len=69),      parameter, public :: shyphens = '*'//hyphens

  logical,                           public :: flush_messages


  ! min_lun in io.F90 is equal to 10. We hardwire this here since we
  ! cannot write "use io" above. Unit 8 and 9 should always be available.
  integer, parameter, public :: iunit_out = 8
  integer, parameter, public :: iunit_err = 9
  ! max_lun is currently 99, i.e. we can hardwire unit_offset above 1000
  integer, parameter, private :: unit_offset = 1000
  character(len=512), private :: msg


  ! ---------------------------------------------------------
  ! Prints out to iunit a message in the form:
  ! ["InputVariable" = value]
  ! where "InputVariable" is given by var.
  ! Since the variable can be integer, real, or logical, we
  ! need a generic interface.
  ! ---------------------------------------------------------
  interface messages_print_var_value
    module procedure messages_print_var_valuei
    module procedure messages_print_var_valuer
    module procedure messages_print_var_valuel
  end interface messages_print_var_value


  integer,    public :: global_alloc_err
  integer(8), public :: global_sizeof
contains

  ! ---------------------------------------------------------
  subroutine write_fatal(no_lines)
    integer, intent(in) :: no_lines

    integer :: i

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_err, file='messages.stderr', &
        action='write', position='append')
    end if

    call messages_print_stress(stderr, "FATAL ERROR")
    write(msg, '(a)') '*** Fatal Error (description follows)'
    call flush_msg(stderr, msg)

#ifdef HAVE_MPI
    call flush_msg(stderr, shyphens)
    write(msg, '(a,i4)') "* From node = ", mpi_world%rank
    call flush_msg(stderr, msg)
#endif
    call flush_msg(stderr, shyphens)
    do i=1,no_lines
      write(msg, '(a,1x,a)') '*', trim(message(i))
      call flush_msg(stderr, msg)
    end do

    ! We only dump the stack in debug mode because subroutine invocations
    ! are only recorded in debug mode (via push_sub/pop_sub). Otherwise,
    ! it is a bit confusing that the stack seems to be empty.
    if(in_debug_mode) then
      call flush_msg(stderr, shyphens)

      write(msg, '(a)') '* Stack: '
      call flush_msg(stderr, msg, 'no')
      do i = 1, no_sub_stack
        write(msg, '(a,a)') ' > ', trim(sub_stack(i))
        call flush_msg(stderr, msg, 'no')
      end do
      call flush_msg(stderr, " ")
    end if

    call messages_print_stress(stderr)

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      close(iunit_err)
    end if

    ! switch file indicator to state aborted
    call switch_status('aborted')

#ifdef HAVE_MPI
    call MPI_Finalize(mpi_err)
#endif

    stop
  end subroutine write_fatal


  ! ---------------------------------------------------------
  subroutine write_warning(no_lines, all_nodes)
    integer,           intent(in) :: no_lines
    logical, optional, intent(in) :: all_nodes

    integer :: i
    logical :: have_to_write
#ifdef HAVE_MPI
    logical :: all_nodes_
#endif

    have_to_write = mpi_grp_is_root(mpi_world)

#ifdef HAVE_MPI
    all_nodes_ = .false.
    if(present(all_nodes)) then
      have_to_write = have_to_write .or. all_nodes
      all_nodes_ = all_nodes
    end if
#endif

    if(have_to_write) then

      if(flush_messages) open(unit=iunit_err, file='messages.stderr', action='write', position='append')
      
      call flush_msg(stderr, '')
      write(msg, '(a)') '** Warning:'
      call flush_msg(stderr, msg)

#ifdef HAVE_MPI
      if(all_nodes_) then
        write(msg , '(a,i4)') '** From node = ', mpi_world%rank
        call flush_msg(stderr, msg)
      end if
#endif

      do i = 1, no_lines
        write(msg , '(a,3x,a)') '**', trim(message(i))
        call flush_msg(stderr, msg)
      end do
      call flush_msg(stderr, '')

#ifdef HAVE_FLUSH
      call flush(stderr)
#endif
      
      if(flush_messages) close(iunit_err)
      
    end if

    return
  end subroutine write_warning


  ! ---------------------------------------------------------
  subroutine write_info(no_lines, iunit, verbose_limit, stress)
    integer, intent(in) :: no_lines
    integer, intent(in), optional :: iunit
    integer, intent(in), optional :: verbose_limit
    logical, optional, intent(in) :: stress

    integer :: i, iu

    if(.not.mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(iunit)) then
      iu = iunit
    else
      iu = stdout
    end if

    if(present(stress)) then
      call messages_print_stress(iu)
    end if

    do i = 1, no_lines
      if(.not.present(verbose_limit)) then
        write(msg, '(a)') trim(message(i))
        call flush_msg(iu, msg)
      else if(in_debug_mode) then
        write(msg, '(a)') trim(message(i))
        call flush_msg(iu, msg)
      end if
    end do
    if(present(stress)) then
      call messages_print_stress(iu)
    end if

    if(flush_messages) close(iunit_out)

#ifdef HAVE_FLUSH
    call flush(iu)
#endif
  end subroutine write_info


  ! ---------------------------------------------------------
  subroutine write_debug(no_lines)
    integer, intent(in) :: no_lines

    integer             :: i, iunit

    if(.not.in_debug_mode) return

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    call open_debug_trace(iunit)
    do i = 1, no_lines
      write(msg, '(a)') trim(message(i))
      call flush_msg(iunit, msg)
    end do
    close(iunit)

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      close(iunit_out)
    end if

  end subroutine write_debug


  ! ---------------------------------------------------------
  subroutine write_debug_newlines(no_lines)
    integer, intent(in) :: no_lines

    integer             :: i, iunit

    if(.not.in_debug_mode) return
    if(mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    call open_debug_trace(iunit)
    do i = 1, no_lines
      write(msg, '(a)') '* -'
      call flush_msg(iunit, msg)
    end do
    close(iunit)

    if(flush_messages) close(iunit_out)

  end subroutine write_debug_newlines


  ! ---------------------------------------------------------
  subroutine write_debug_marker(no)
    integer, intent(in) :: no

    if(.not.in_debug_mode) return

    write(message(1), '(a,i3)') 'debug marker #',no
    call write_debug(1)

  end subroutine write_debug_marker


  ! ---------------------------------------------------------
  subroutine switch_status(status)
    character(len=*), intent(in) :: status

    ! only root node is taking care of file I/O
    if(.not.mpi_grp_is_root(mpi_world)) return     

    ! remove old status files first, before we switch to state aborted
    call loct_rm_status_files(current_label)
    
    ! create empty status file to indicate 'aborted state'
    open(unit=iunit_err, file='exec/'//trim(current_label)//'oct-status-'//trim(status), &
      action='write', status='unknown')
    close(iunit_err)
    
  end subroutine switch_status


  ! ---------------------------------------------------------
  subroutine open_debug_trace(iunit)
    integer, intent(out) :: iunit

    character(len=4) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i3.3)') iunit - unit_offset
    call loct_mkdir(trim(current_label)//'debug')
    open(iunit, file=trim(current_label)//'debug/debug_trace.node.'//filenum, &
      action='write', status='unknown', position='append')

  end subroutine open_debug_trace


  ! ---------------------------------------------------------
  subroutine alloc_error(size, file, line)
    integer(8),       intent(in) :: size
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    write(message(1), '(a,i10,3a,i5)') "Failed to allocate ", size, " words in file '", trim(file), "' line ", line
    call write_fatal(1)

  end subroutine alloc_error


  ! ---------------------------------------------------------
  subroutine input_error(var)
    character(len=*), intent(in) :: var

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call messages_print_stress(stderr, "INPUT ERROR")
      write(msg, '(a)') '*** Fatal Error in input '
      call flush_msg(stderr, msg)
      call flush_msg(stderr, shyphens)

      call varinfo_print(stderr, var)
      call messages_print_stress(stderr)
    end if

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      close(iunit_out)
    end if

    ! switch file indicator to state aborted
    call switch_status('aborted')

#ifdef HAVE_MPI
    call MPI_Finalize(mpi_err)
#endif
    stop
  end subroutine input_error
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuei(iunit, var, value)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: value
    character(len=10) :: intstring
    if(.not.mpi_grp_is_root(mpi_world)) return
    write(intstring,'(i10)') value
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(intstring))//']'
  end subroutine messages_print_var_valuei
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuer(iunit, var, value)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    FLOAT,            intent(in) :: value
    character(len=10) :: floatstring
    if(.not.mpi_grp_is_root(mpi_world)) return
    write(floatstring,'(g10.4)') value
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(floatstring))//']'
  end subroutine messages_print_var_valuer
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuel(iunit, var, value)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    logical,          intent(in) :: value
    character(len=1) :: lstring
    if(.not.mpi_grp_is_root(mpi_world)) return
    write(lstring,'(l1)') value
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//lstring//']'
  end subroutine messages_print_var_valuel
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_info(iunit, var)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var

    if(.not.mpi_grp_is_root(mpi_world)) return

    call varinfo_print(iunit, var)
  end subroutine messages_print_var_info


  ! ---------------------------------------------------------
  subroutine messages_print_var_option(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option
    character(len=*), intent(in), optional :: pre

    if(.not.mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(pre)) then
      call varinfo_print_option(iunit, var, option, pre)
      if(flush_messages) then
        call varinfo_print_option(iunit_out, var, option, pre)
      end if
    else
      call varinfo_print_option(iunit, var, option)
      if(flush_messages) then
        call varinfo_print_option(iunit_out, var, option, pre)
      end if
    end if

    if(flush_messages) then
      close(iunit_out)
    end if
  end subroutine messages_print_var_option


  ! ---------------------------------------------------------
  subroutine messages_print_stress(iunit, msg)
    integer, intent(in) :: iunit
    character(len=*), intent(in), optional :: msg

    integer, parameter :: max_len = 70

    integer :: i, j, l
    character(len=70) :: str

    if(.not.mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(msg)) then
      l   = len(msg)

      str = ''; j = 1

      do i = 1, (max_len - (l + 2))/2
        str(j:j) = '*'; j = j + 1
      end do
 
      str(j:j) = ' '; j = j + 1
     
      do i = 1, l
        str(j:j) = msg(i:i); j = j + 1
      end do

      str(j:j) = ' '; j = j + 1

      do i = j, max_len
        str(j:j) = '*'; j = j + 1
      end do

      call flush_msg(iunit, '')   ! empty line
      call flush_msg(iunit, str)  ! out nice line with the header
    else
      do i = 1, max_len
        str(i:i) = '*'
      end do

      call flush_msg(iunit, str)  ! out nice line with the header
      call flush_msg(iunit, '')   ! empty line
    end if

    if(flush_messages) close(iunit_out)

#ifdef HAVE_FLUSH
    call flush(iunit)
#endif
  end subroutine messages_print_stress


  ! ---------------------------------------------------------
  subroutine flush_msg(iunit, str, adv)
    integer,            intent(in)           :: iunit
    character(len = *), intent(in)           :: str
    character(len = *), intent(in), optional :: adv

    character(len = 20) :: adv_

    adv_ = 'yes'
    if(present(adv)) adv_ = adv

    write(iunit, '(a)', advance=trim(adv_)) trim(str)
    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      if(iunit.eq.stderr) write(iunit_err, '(a)', advance=trim(adv_)) trim(str)
      if(iunit.eq.stdout) write(iunit_out, '(a)', advance=trim(adv_)) trim(str)
    end if

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
  subroutine epoch_time_diff(sec, usec)
    integer, intent(inout) :: sec, usec

    call time_diff(s_epoch_sec, s_epoch_usec, sec, usec)

  end subroutine epoch_time_diff


  ! ---------------------------------------------------------
  ! Computes t2 <- t2-t1. sec1,2 and usec1,2 are
  ! seconds,microseconds of t1,2
  subroutine time_diff(sec1, usec1, sec2, usec2)
    integer, intent(in)    :: sec1, usec1
    integer, intent(inout) :: sec2, usec2

    ! Correct overflow.
    if(usec2-usec1.lt.0) then
      usec2 = 1000000 + usec2
      if(sec2.ge.sec1) then
        sec2 = sec2-1
      end if
    end if

    ! Replace values.
    if(sec2.ge.sec1) then
      sec2 = sec2-sec1
    end if
    usec2 = usec2-usec1

  end subroutine time_diff


  ! ---------------------------------------------------------
  ! Computes t2 <- t1+t2. Parameters as in time_diff
  ! Assert: t1,2 <= 0.
  subroutine time_sum(sec1, usec1, sec2, usec2)
    integer, intent(in)    :: sec1, usec1
    integer, intent(inout) :: sec2, usec2

    sec2  = sec1+sec2
    usec2 = usec1+usec2

    ! Carry?
    if(usec2.ge.1000000) then
      sec2  = sec2+1
      usec2 = usec2-1000000
    end if

  end subroutine time_sum


#ifndef NDEBUG
  ! ---------------------------------------------------------
  subroutine push_sub(sub_name)
    character(len=*), intent(in) :: sub_name

    integer i, iunit, sec, usec

    if(.not.in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)

    no_sub_stack = no_sub_stack + 1
    if(no_sub_stack > 49) then
      sub_stack(50) = 'push_sub'
      message(1) = 'Too many recursion levels (max=50)'
      call write_fatal(1)
    end if

    sub_stack(no_sub_stack)  = trim(sub_name)
    time_stack(no_sub_stack) = loct_clock()

    if(conf%debug_level.ge.100) then
      call open_debug_trace(iunit)
      call push_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    end if

    ! also write to stderr if we are node 0
    if(conf%debug_level.gt.1) then
      if (mpi_grp_is_root(mpi_world)) call push_sub_write(stderr)
    end if

  contains

    subroutine push_sub_write(iunit_out)
      integer,  intent(in) :: iunit_out

      write(iunit_out,'(a,i6,a,i6.6,f20.6,i8,a)', advance='no') "* I ", &
        sec,'.',usec, &
        loct_clock(), &
        get_memory_usage()/1024, " | "
      do i = no_sub_stack-1, 1, -1
        write(iunit_out,'(a)', advance='no') "..|"
      end do
      write(iunit_out,'(a)') trim(sub_name)

    end subroutine push_sub_write

  end subroutine push_sub


  ! ---------------------------------------------------------
  subroutine pop_sub()
    integer i, iunit, sec, usec

    if(.not.in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)

    if(no_sub_stack <= 0) then
      no_sub_stack = 1
      sub_stack(1) = 'pop_sub'
      message(1) = 'Too few recursion levels'
      call write_fatal(1)
    end if

    if(conf%debug_level.ge.100) then
      call open_debug_trace(iunit)
      call pop_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    end if
      
    ! also write to stderr if we are node 0
    if(conf%debug_level.gt.1) then
      if (mpi_grp_is_root(mpi_world)) call pop_sub_write(stderr)
    end if
    
    no_sub_stack = no_sub_stack - 1

  contains

    subroutine pop_sub_write(iunit_out)
      integer, intent(in) :: iunit_out

      write(iunit_out,'(a,i6,a,i6.6,f20.6,i8, a)', advance='no') "* O ", &
        sec,'.',usec, &
        loct_clock()-time_stack(no_sub_stack), &
        get_memory_usage()/1024, " | "
      do i = no_sub_stack-1, 1, -1
        write(iunit_out,'(a)', advance='no') "..|"
      end do
      write(iunit_out,'(a)') trim(sub_stack(no_sub_stack))

    end subroutine pop_sub_write

  end subroutine pop_sub
#endif
  
  subroutine obsolete_variable(name, rep)
    character(len=*),           intent(in) :: name
    character(len=*), optional, intent(in) :: rep
    
    if ( loct_parse_isdef(trim(name)) /= 0 ) then 

      write(message(1), '(a)') 'Input variable '//trim(name)//' is obsolete.'

      if(present(rep)) then
        write(message(2), '(a)') 'Please use variable '//trim(rep)//' instead.'
        call write_fatal(2)
      else
        call write_fatal(1)
      end if

    end if
    
  end subroutine obsolete_variable

  subroutine messages_devel_version(name)
    character(len=*),           intent(in) :: name
    
    if(.not. conf%devel_version) then
      write(message(1), '(a)') 'Error: '//trim(name)//' is under development.'
      write(message(2), '(a)') 'To use it (at your own risk) set variable DevelVersion to yes.'
      call write_fatal(2)
    else
      write(message(1), '(a)') trim(name)//' is under development.'
      write(message(2), '(a)') 'It might not work or produce wrong results.'
      call write_warning(2)
    end if

  end subroutine messages_devel_version
  
end module messages_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
