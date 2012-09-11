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
  use mpi_m
  use parser_m
  use string_m
  use unit_m
  use varinfo_m

  implicit none

  private

  public ::                     &
    messages_init,              &
    messages_end,               &
    messages_fatal,             &
    messages_warning,           &
    messages_info,              &
    messages_debug,             &
    messages_debug_marker,      &
    messages_debug_newlines,    &
    delete_debug_trace,         &
    print_date,                 &
    time_diff,                  &
    time_sum,                   &
    epoch_time_diff,            &
    alloc_error,                &
    dealloc_error,              &
    input_error,                &
    push_sub,                   &
    pop_sub,                    &
    messages_print_stress,      &
    messages_print_var_info,    &
    messages_print_var_option,  &
    messages_print_var_value,   &
    messages_obsolete_variable, &
    messages_experimental,      &
    messages_check_def,         &
    messages_not_implemented,   &
    messages_new_line,          &
    messages_write,             &
    messages_clean_path

  integer, parameter :: max_lines = 20
  character(len=256), dimension(max_lines), public :: message    !< to be output by fatal, warning
  character(len=68),      parameter, public :: hyphens = &
    '--------------------------------------------------------------------'
  character(len=69),      parameter, public :: shyphens = '*'//hyphens

  logical,                           public :: flush_messages


  !> min_lun in io.F90 is equal to 10. We hardwire this here since we
  !! cannot write "use io" above. Unit 8 and 9 should always be available.
  integer, parameter, public :: iunit_out = 8
  integer, parameter, public :: iunit_err = 9
  !> max_lun is currently 99, i.e. we can hardwire unit_offset above 1000
  integer, parameter, private :: unit_offset = 1000
  character(len=512), private :: msg


  ! ---------------------------------------------------------
  !> Prints out to iunit a message in the form:
  !! ["InputVariable" = value]
  !! where "InputVariable" is given by var.
  !! Since the variable can be integer, real, or logical, we
  !! need a generic interface.
  ! ---------------------------------------------------------
  interface messages_print_var_value
    module procedure messages_print_var_valuei
    module procedure messages_print_var_valuer
    module procedure messages_print_var_valuel
  end interface messages_print_var_value

  interface messages_write
    module procedure messages_write_float
    module procedure messages_write_integer
    module procedure messages_write_integer8
    module procedure messages_write_str
    module procedure messages_write_logical
  end interface messages_write

  integer,    public :: global_alloc_err
  integer(8), public :: global_sizeof

  integer :: warnings
  integer :: experimentals
  integer :: current_line

contains

  ! ---------------------------------------------------------

  subroutine messages_init()

    call messages_obsolete_variable('DevelVersion', 'ExperimentalFeatures')

    !%Variable ExperimentalFeatures
    !%Type logical
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% If true, allows the use of certain parts of the code that are
    !% still under development and are not suitable for production
    !% runs. This should not be used unless you know what you are doing.
    !%End
    call parse_logical('ExperimentalFeatures', .false., conf%devel_version)

    !%Variable DebugLevel
    !%Type integer
    !%Default 0
    !%Section Execution::Debug
    !%Description
    !% This variable decides whether or not to enter debug mode.
    !% If it is greater than 0, different amounts of additional information
    !% are written to standard output and additional assertion checks are performed.
    !%Option 0
    !% (default) <tt>Octopus</tt> does not enter debug mode.
    !%Option 1
    !% Moderate amount of debug output; assertion checks enabled.
    !%Option 2
    !% The code prints a stack trace as it enters end exits subroutines.
    !% This is useful for developers and you should include this output when
    !% submitting a bug report.
    !%Option 99
    !% The debug output is additionally written to files in the <tt>debug</tt>
    !% directory. For each node (when running in parallel) there is a file called
    !% <tt>debug_trace.&lt;rank&gt;</tt>. Writing these files slows down the code by a huge factor and
    !% it is usually only necessary for parallel runs. In the serial case all
    !% the information can be obtained from standard out.
    !%End
    call parse_integer('DebugLevel', 0, conf%debug_level)
    if(conf%debug_level>0) then
      in_debug_mode = .true.
    else
      in_debug_mode = .false.
    end if
    
    warnings = 0
    experimentals = 0    

    call messages_reset_lines()

  end subroutine messages_init

  ! ---------------------------------------------------------

  subroutine messages_end()

    if(mpi_grp_is_root(mpi_world)) then
  
      if(experimentals > 0 .or. warnings > 0) then
        message(1) = ''
        call messages_info(1)      
      end if
      
      if(experimentals > 0) then
        call messages_write('Octopus used ')
        call messages_write(experimentals)
        call messages_write(' experimental feature(s).')
        call messages_info()
      end if
      
      if(warnings > 0) then
        call messages_write('Octopus emitted ')
        call messages_write(warnings)
        call messages_write(' warning(s).')
        call messages_info()
      end if
      
      open(unit = iunit_out, file = 'exec/messages', action = 'write')
      write(iunit_out, '(a, i9)') "warnings          = ", warnings
      write(iunit_out, '(a, i9)') "experimental      = ", experimentals
      close(iunit_out)
    
    end if
  
end subroutine messages_end

  ! ---------------------------------------------------------
  subroutine messages_fatal(no_lines, only_root_writes)
    integer, optional, intent(in) :: no_lines
    logical, optional, intent(in) :: only_root_writes

    integer :: ii, no_lines_
    logical :: only_root_writes_, should_write

    no_lines_ = current_line
    if(present(no_lines)) no_lines_ = no_lines

    if(present(only_root_writes)) then
      should_write = mpi_grp_is_root(mpi_world) .or. (.not. only_root_writes)
      only_root_writes_ = only_root_writes
    else
      should_write = .true.
      only_root_writes_ = .false.
    endif

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_err, file='messages.stderr', &
        action='write', position='append')
    end if

    if(should_write) then
      call messages_print_stress(stderr, "FATAL ERROR")
      write(msg, '(a)') '*** Fatal Error (description follows)'
      call flush_msg(stderr, msg)

#ifdef HAVE_MPI
      if(.not. only_root_writes_) then
        call flush_msg(stderr, shyphens)
        write(msg, '(a,i4)') "* From node = ", mpi_world%rank
        call flush_msg(stderr, msg)
      endif
#endif
      call flush_msg(stderr, shyphens)
      do ii = 1, no_lines_
        write(msg, '(a,1x,a)') '*', trim(message(ii))
        call flush_msg(stderr, msg)
      end do
    endif

    ! We only dump the stack in debug mode because subroutine invocations
    ! are only recorded in debug mode (via push_sub/pop_sub). Otherwise,
    ! it is a bit confusing that the stack seems to be empty.
    if(in_debug_mode) then
      call flush_msg(stderr, shyphens)

      write(msg, '(a)') '* Stack: '
      call flush_msg(stderr, msg, 'no')
      do ii = 1, no_sub_stack
        write(msg, '(a,a)') ' > ', trim(sub_stack(ii))
        call flush_msg(stderr, msg, 'no')
      end do
      call flush_msg(stderr, " ")
    end if

    if(should_write) then
      call messages_print_stress(stderr)
    endif

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      close(iunit_err)
    end if

    ! switch file indicator to state aborted
    call switch_status('aborted')

#ifdef HAVE_MPI
    call MPI_Abort(mpi_world%comm, mpi_err)
#endif

    call loct_exit_failure()
  end subroutine messages_fatal


  ! ---------------------------------------------------------
  subroutine messages_warning(no_lines, all_nodes)
    integer, optional, intent(in) :: no_lines
    logical, optional, intent(in) :: all_nodes

    integer :: il, no_lines_
    logical :: have_to_write
#ifdef HAVE_MPI
    logical :: all_nodes_
#endif
    
    no_lines_ = current_line
    if(present(no_lines)) no_lines_ = no_lines

    INCR(warnings, 1)

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

      do il = 1, no_lines_
        write(msg , '(a,3x,a)') '**', trim(message(il))
        call flush_msg(stderr, msg)
      end do
      call flush_msg(stderr, '')

#ifdef HAVE_FLUSH
      call flush(stderr)
#endif
      
      if(flush_messages) close(iunit_err)
      
    end if

    call messages_reset_lines()

  end subroutine messages_warning

  ! ---------------------------------------------------------

  subroutine messages_info(no_lines, iunit, verbose_limit, stress, all_nodes)
    integer, optional, intent(in) :: no_lines
    integer, optional, intent(in) :: iunit
    logical, optional, intent(in) :: verbose_limit
    logical, optional, intent(in) :: stress
    logical, optional, intent(in) :: all_nodes

    integer :: il, iu, no_lines_

    if(.not. mpi_grp_is_root(mpi_world) .and. .not. optional_default(all_nodes, .false.)) then 
      call messages_reset_lines()
      return
    end if

    no_lines_ = current_line
    if(present(no_lines)) no_lines_ = no_lines

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

    do il = 1, no_lines_
      if(.not. present(verbose_limit) .or. in_debug_mode) then
        write(msg, '(a)') trim(message(il))
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

    call messages_reset_lines()

  end subroutine messages_info


  ! ---------------------------------------------------------
  subroutine messages_debug(no_lines)
    integer, intent(in) :: no_lines

    integer             :: il, iunit

    if(.not.in_debug_mode) return

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    call open_debug_trace(iunit)
    do il = 1, no_lines
      write(msg, '(a)') trim(message(il))
      call flush_msg(iunit, msg)
    end do
    close(iunit)

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      close(iunit_out)
    end if

  end subroutine messages_debug


  ! ---------------------------------------------------------
  subroutine messages_debug_newlines(no_lines)
    integer, intent(in) :: no_lines

    integer             :: il, iunit

    if(.not. in_debug_mode) return
    if(mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    call open_debug_trace(iunit)
    do il = 1, no_lines
      write(msg, '(a)') '* -'
      call flush_msg(iunit, msg)
    end do
    close(iunit)

    if(flush_messages) close(iunit_out)

  end subroutine messages_debug_newlines


  ! ---------------------------------------------------------
  subroutine messages_debug_marker(no)
    integer, intent(in) :: no

    if(.not. in_debug_mode) return

    write(message(1), '(a,i3)') 'debug marker #', no
    call messages_debug(1)

  end subroutine messages_debug_marker


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

    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir(trim(current_label)//'debug')
    open(iunit, file=trim(current_label)//'debug/debug_trace.node.'//filenum, &
      action='write', status='unknown', position='append')

  end subroutine open_debug_trace

  ! ---------------------------------------------------------
  subroutine delete_debug_trace()

    integer :: iunit
    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir(trim(current_label)//'debug')
    call loct_rm(trim(current_label)//'debug/debug_trace.node.'//filenum)

  end subroutine delete_debug_trace


  ! ---------------------------------------------------------
  subroutine alloc_error(size, file, line)
    integer(8),       intent(in) :: size
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    write(message(1), '(a,i14,3a,i5)') "Failed to allocate ", size, " words in file '", trim(file), "' line ", line
    call messages_fatal(1)

  end subroutine alloc_error


  ! ---------------------------------------------------------
  subroutine dealloc_error(size, file, line)
    integer(8),       intent(in) :: size
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    write(message(1), '(a,i14,3a,i5)') "Failed to deallocate array of ", size, " words in file '", trim(file), "' line ", line
    call messages_fatal(1)

  end subroutine dealloc_error


  ! ---------------------------------------------------------
  subroutine input_error(var)
    character(len=*), intent(in) :: var

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
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

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      close(iunit_out)
    end if

    ! switch file indicator to state aborted
    call switch_status('aborted')

#ifdef HAVE_MPI
    call MPI_Abort(mpi_world%comm, mpi_err)
#endif

    call loct_exit_failure()
  end subroutine input_error
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuei(iunit, var, value)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: value

    character(len=10) :: intstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    write(intstring,'(i10)') value
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(intstring))//']'

  end subroutine messages_print_var_valuei
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuer(iunit, var, value, unit)
    integer,                intent(in) :: iunit
    character(len=*),       intent(in) :: var
    FLOAT,                  intent(in) :: value
    type(unit_t), optional, intent(in) :: unit

    character(len=10) :: floatstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    if(.not. present(unit)) then
      write(floatstring,'(g10.4)') value
      write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(floatstring))//']'
    else
      write(floatstring,'(g10.4)') units_from_atomic(unit, value)
      write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(floatstring))//' '//trim(units_abbrev(unit))//']'
    end if

  end subroutine messages_print_var_valuer
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuel(iunit, var, value)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    logical,          intent(in) :: value

    character(len=3) :: lstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    if(value) then
      lstring = 'yes'
    else
      lstring = 'no'
    end if
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(lstring)//']'

  end subroutine messages_print_var_valuel
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_info(iunit, var)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var

    if(.not. mpi_grp_is_root(mpi_world)) return

    call varinfo_print(iunit, var)
  end subroutine messages_print_var_info


  ! ---------------------------------------------------------
  subroutine messages_print_var_option(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: option
    character(len=*), intent(in), optional :: pre

    if(.not. mpi_grp_is_root(mpi_world)) return

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
    integer,                    intent(in) :: iunit
    character(len=*), optional, intent(in) :: msg

    integer, parameter :: max_len = 70

    integer :: ii, jj, length
    character(len=70) :: str

    if(.not.mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(msg)) then
      length = len(msg)

      str = ''
      jj = 1

      do ii = 1, (max_len - (length + 2))/2
        str(jj:jj) = '*'
        jj = jj + 1
      end do
 
      str(jj:jj) = ' '
      jj = jj + 1
     
      do ii = 1, length
        str(jj:jj) = msg(ii:ii)
        jj = jj + 1
      end do

      str(jj:jj) = ' '
      jj = jj + 1

      do ii = jj, max_len
        str(jj:jj) = '*'
        jj = jj + 1
      end do

      call flush_msg(iunit, '')   ! empty line
      call flush_msg(iunit, str)  ! out nice line with the header
    else
      do ii = 1, max_len
        str(ii:ii) = '*'
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
    integer,                      intent(in) :: iunit
    character(len = *),           intent(in) :: str
    character(len = *), optional, intent(in) :: adv

    character(len = 20) :: adv_

    adv_ = 'yes'
    if(present(adv)) adv_ = adv

    write(iunit, '(a)', advance=trim(adv_)) trim(str)
    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      if(iunit .eq. stderr) write(iunit_err, '(a)', advance=trim(adv_)) trim(str)
      if(iunit .eq. stdout) write(iunit_out, '(a)', advance=trim(adv_)) trim(str)
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
    call messages_info(3)

  end subroutine print_date

  ! ---------------------------------------------------------
  subroutine epoch_time_diff(sec, usec)
    integer, intent(inout) :: sec
    integer, intent(inout) :: usec

    ! this is called by push/pop so there cannot be a push/pop in this routine

    call time_diff(s_epoch_sec, s_epoch_usec, sec, usec)
  end subroutine epoch_time_diff


  ! ---------------------------------------------------------
  !> Computes t2 <- t2-t1. sec1,2 and usec1,2 are
  !! seconds,microseconds of t1,2
  subroutine time_diff(sec1, usec1, sec2, usec2)
    integer, intent(in)    :: sec1
    integer, intent(in)    :: usec1
    integer, intent(inout) :: sec2
    integer, intent(inout) :: usec2

    ! this is called by push/pop so there cannot be a push/pop in this routine

    ! Correct overflow.
    if(usec2 - usec1 .lt. 0) then
      usec2 = 1000000 + usec2
      if(sec2 .ge. sec1) then
        sec2 = sec2 - 1
      end if
    end if

    ! Replace values.
    if(sec2 .ge. sec1) then
      sec2 = sec2 - sec1
    end if
    usec2 = usec2 - usec1

  end subroutine time_diff

  ! ---------------------------------------------------------
  !> Computes t2 <- t1+t2. Parameters as in time_diff
  !! Assert: t1,2 <= 0.
  subroutine time_sum(sec1, usec1, sec2, usec2)
    integer, intent(in)    :: sec1
    integer, intent(in)    :: usec1
    integer, intent(inout) :: sec2
    integer, intent(inout) :: usec2

    PUSH_SUB(time_sum)

    sec2  = sec1 + sec2
    usec2 = usec1 + usec2

    ! Carry?
    if(usec2 .ge. 1000000) then
      sec2  = sec2 + 1
      usec2 = usec2 - 1000000
    end if

    POP_SUB(time_sum)
  end subroutine time_sum


#ifndef NDEBUG
  ! ---------------------------------------------------------
  subroutine push_sub(sub_name)
    character(len=*), intent(in) :: sub_name

    integer iunit, sec, usec

    if(.not. in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)

    no_sub_stack = no_sub_stack + 1
    if(no_sub_stack > 49) then
      sub_stack(50) = 'push_sub'
      message(1) = 'Too many recursion levels (max=50)'
      call messages_fatal(1)
    end if

    sub_stack(no_sub_stack)  = trim(messages_clean_path(sub_name))
    time_stack(no_sub_stack) = loct_clock()

    if(conf%debug_level .ge. 99) then
      call open_debug_trace(iunit)
      call push_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    else if(conf%debug_level > 1 .and. mpi_grp_is_root(mpi_world)) then
      ! write to stderr if we are node 0
      call push_sub_write(stderr)
    end if

  contains

    subroutine push_sub_write(iunit_out)
      integer,  intent(in) :: iunit_out

      integer :: ii
      character(len=200) :: tmpstr

      write(tmpstr,'(a,i6,a,i6.6,f20.6,i8,a)') "* I ", &
        sec, '.', usec, &
        loct_clock(), &
        get_memory_usage() / 1024, " | "
      do ii = no_sub_stack - 1, 1, -1
        write(tmpstr, '(2a)') trim(tmpstr), "..|"
      end do
      write(tmpstr, '(2a)') trim(tmpstr), trim(messages_clean_path(sub_name))
      call flush_msg(iunit_out, tmpstr)

    end subroutine push_sub_write

  end subroutine push_sub


  ! ---------------------------------------------------------
  subroutine pop_sub(sub_name)
    character(len=*), intent(in) :: sub_name
    
    character(len=80) :: sub_name_short

    integer iunit, sec, usec

    if(.not. in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)

    if(no_sub_stack <= 0) then
      no_sub_stack = 1
      sub_stack(1) = 'pop_sub'
      message(1) = 'Too few recursion levels.'
      call messages_fatal(1)
    end if

    ! the name might be truncated in sub_stack, so we copy to a string
    ! of the same size
    sub_name_short = messages_clean_path(sub_name)

    if(sub_name_short .ne. sub_stack(no_sub_stack)) then
      write (message(1),'(a)') 'Wrong sub name on pop_sub :'
      write (message(2),'(2a)') ' got      : ', sub_name_short
      write (message(3),'(2a)') ' expected : ', sub_stack(no_sub_stack)
      call messages_fatal(3)
    end if

    if(conf%debug_level .ge. 99) then
      call open_debug_trace(iunit)
      call pop_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    else if (conf%debug_level .gt. 1 .and. mpi_grp_is_root(mpi_world)) then
      ! write to stderr if we are node 0
      call pop_sub_write(stderr)
    end if
    
    no_sub_stack = no_sub_stack - 1

  contains

    subroutine pop_sub_write(iunit_out)
      integer, intent(in) :: iunit_out

      integer :: ii
      character(len=200) :: tmpstr

      write(tmpstr,'(a,i6,a,i6.6,f20.6,i8, a)') "* O ", &
        sec, '.', usec, &
        loct_clock() - time_stack(no_sub_stack), &
        get_memory_usage() / 1024, " | "
      do ii = no_sub_stack - 1, 1, -1
        write(tmpstr,'(2a)') trim(tmpstr), "..|"
      end do
      write(tmpstr,'(2a)') trim(tmpstr), trim(sub_stack(no_sub_stack))
      call flush_msg(iunit_out, tmpstr)

    end subroutine pop_sub_write

  end subroutine pop_sub
#endif
  
  ! ---------------------------------------------------------
  subroutine messages_obsolete_variable(name, rep)
    character(len=*),           intent(in) :: name
    character(len=*), optional, intent(in) :: rep
    
    if ( parse_isdef(trim(name)) /= 0 ) then 

      write(message(1), '(a)') 'Input variable '//trim(name)//' is obsolete.'

      if(present(rep)) then
        write(message(2), '(a)') ' '
        write(message(3), '(a)') 'Equivalent functionality can be obtained with the '//trim(rep)
        write(message(4), '(a)') 'variable. Check the documentation for details.'
        write(message(5), '(a)') '(You can use the `oct-help -s '//trim(rep)//'` command).'
        call messages_fatal(5, only_root_writes = .true.)
      else
        call messages_fatal(1, only_root_writes = .true.)
      end if

    end if
    
  end subroutine messages_obsolete_variable

  ! ---------------------------------------------------------
  subroutine messages_experimental(name)
    character(len=*), intent(in) :: name
    
    INCR(experimentals, 1)

    if(.not. conf%devel_version) then
      write(message(1), '(a)') trim(name)//' is under development.'
      write(message(2), '(a)') 'To use it (at your own risk) set the variable ExperimentalFeatures to yes.'
      call messages_fatal(2, only_root_writes = .true.)
    else
      write(message(1), '(a)') trim(name)//' is under development.'
      write(message(2), '(a)') 'It might not work or produce wrong results.'
      call messages_warning(2)

      ! remove this warning from the count
      INCR(warnings, -1)
    end if

  end subroutine messages_experimental
  

  !--------------------------------------------------------------
  subroutine messages_check_def(var, def, text)
    FLOAT,            intent(in) :: var
    FLOAT,            intent(in) :: def
    character(len=*), intent(in) :: text

    PUSH_SUB(messages_check_def)

    if(var > def) then
      write(message(1), '(3a)') "The value for '", text, "' does not match the recommended value."
      write(message(2), '(f8.3,a,f8.3)') var, ' > ', def
      call messages_warning(2)
    end if

    POP_SUB(messages_check_def)
  end subroutine messages_check_def


  ! ------------------------------------------------------------
  subroutine messages_not_implemented(feature)
    character(len=*), intent(in) :: feature

    PUSH_SUB(messages_not_implemented)

    message(1) = trim(feature)//" not implemented."
    call messages_fatal(1, only_root_writes = .true.)

    POP_SUB(messages_not_implemented)
  end subroutine messages_not_implemented

  ! ------------------------------------------------------------

  subroutine messages_reset_lines()

    current_line = 1
    message(1) = ''
    
  end subroutine messages_reset_lines

  ! ------------------------------------------------------------

  subroutine messages_new_line()
    
    current_line = current_line + 1
    message(current_line) = ''

    if(current_line > max_lines) stop 'Too many message lines.'
   
  end subroutine messages_new_line

  ! ------------------------------------------------------------

  subroutine messages_write_float(val, fmt, new_line, units, align_left, print_units)
    FLOAT,                      intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line
    type(unit_t),     optional, intent(in) :: units
    logical,          optional, intent(in) :: align_left
    logical,          optional, intent(in) :: print_units

    character(len=10) :: number
    FLOAT            :: tval

    tval = val
    if(present(units)) tval = units_from_atomic(units, val)
    
    if(present(fmt)) then
      write(number, '('//trim(fmt)//')') tval
    else
      write(number, '(f12.6)') tval
    end if

    if(optional_default(align_left, .false.)) number = ' '//adjustl(number)

    write(message(current_line), '(a, a)') trim(message(current_line)), trim(number)

    if(present(units) .and. optional_default(print_units, .false.)) then
      write(message(current_line), '(a, a, a)') trim(message(current_line)), ' ', trim(units_abbrev(units))
    end if

    if(optional_default(new_line, .false.)) call messages_new_line()

  end subroutine messages_write_float

  ! ------------------------------------------------------------

  subroutine messages_write_integer8(val, fmt, new_line)
    integer(8),                 intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line

    character(len=10) :: number

    if(present(fmt)) then
      write(message(current_line), '(a, '//trim(fmt)//')') trim(message(current_line)), val
    else
      write(number, '(i10)') val
      write(message(current_line), '(3a)') trim(message(current_line)), ' ', trim(adjustl(number))
    end if

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_integer8

  ! ------------------------------------------------------------

  subroutine messages_write_integer(val, fmt, new_line)
    integer,                    intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line

    character(len=10) :: number

    if(present(fmt)) then
      write(message(current_line), '(a, '//trim(fmt)//')') trim(message(current_line)), val
    else
      write(number, '(i10)') val
      write(message(current_line), '(3a)') trim(message(current_line)), ' ', trim(adjustl(number))
    end if

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_integer

  ! ------------------------------------------------------------

  subroutine messages_write_str(val, new_line)
    character(len=*),           intent(in) :: val
    logical,          optional, intent(in) :: new_line

    write(message(current_line), '(2a)') trim(message(current_line)), trim(val)

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_str

  ! ------------------------------------------------------------

  subroutine messages_write_logical(val, new_line)
    logical,           intent(in) :: val
    logical, optional, intent(in) :: new_line

    if(val) then
      write(message(current_line), '(2a)') trim(message(current_line)), ' yes'
    else
      write(message(current_line), '(2a)') trim(message(current_line)), ' no'
    end if

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_logical

  ! -----------------------------------------------------------

  character(len=256) function messages_clean_path(filename) result(clean_path)
    character(len=*), intent(in) :: filename

    integer :: pos, start

    pos = index(filename, 'src/', back = .true.)
    if(pos == 0) then
       ! 'src/' does not occur
       start = pos + 1
    else
       ! remove 'src/'
       start = pos + 4
    endif
    clean_path = filename(start:)
  end function messages_clean_path

end module messages_m

! ---------------------------------------------------------
!> This subroutine is called by the assert macro, it is not in a
!> module so it can be called from any file. The interface is declared
!> in global_m.
subroutine assert_die(s, f, l)
  use global_m
  use messages_m
  use mpi_m

  implicit none

  character(len=*), intent(in) :: s, f
  integer, intent(in) :: l
    
  call messages_write('Node ')
  call messages_write(mpi_world%rank)
  call messages_write(':')
  call messages_new_line()

  call messages_write(' Assertion "'//trim(s)//'"')
  call messages_new_line()

  call messages_write(' failed in line ')
  call messages_write(l)
  call messages_write(' of file "'//trim(messages_clean_path(f))//'".')

  call messages_fatal()

end subroutine assert_die



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
