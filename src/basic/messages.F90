!! Copyright (C) 2002-20016 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module messages_oct_m
  use global_oct_m
  use debug_oct_m
  use loct_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use string_oct_m
  use unit_oct_m
  use varinfo_oct_m

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
    messages_switch_status,     &
    delete_debug_trace,         &
    print_date,                 &
    time_diff,                  &
    time_sum,                   &
    epoch_time_diff,            &
    alloc_error,                &
    dealloc_error,              &
    messages_input_error,       &
    push_sub,                   &
    pop_sub,                    &
    messages_print_stress,      &
    messages_print_var_info,    &
    messages_print_var_option,  &
    messages_print_var_value,   &
    messages_obsolete_variable, &
    messages_variable_is_block, &
    messages_experimental,      &
    messages_check_def,         &
    messages_not_implemented,   &
    messages_new_line,          &
    messages_write,             &
    messages_clean_path,        &
    messages_dump_stack,        &
    debug

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
  integer, parameter, private :: SLEEPYTIME_ALL = 1, SLEEPYTIME_NONWRITERS = 60 !< seconds
  character(len=64),  private :: oct_status = 'undefined' !< start with an undefined status


  ! ---------------------------------------------------------
  !> Prints out to iunit a message in the form:
  !! ["InputVariable" = value]
  !! where "InputVariable" is given by var.
  !! Since the variable can be integer, real, logical, or string we
  !! need a generic interface.
  ! ---------------------------------------------------------
  interface messages_print_var_value
    module procedure messages_print_var_valuei
    module procedure messages_print_var_values
    module procedure messages_print_var_valuer
    module procedure messages_print_var_valuel
    module procedure messages_print_var_valuear
  end interface messages_print_var_value

  interface messages_write
    module procedure messages_write_float
    module procedure messages_write_integer
    module procedure messages_write_integer8
    module procedure messages_write_str
    module procedure messages_write_logical
  end interface messages_write


  interface  messages_print_var_option
    module procedure messages_print_var_option_4
    module procedure messages_print_var_option_8
  end interface messages_print_var_option

  integer :: warnings
  integer :: experimentals
  integer :: current_line

  type(debug_t), save :: debug

  !> from signals.c
  interface
    subroutine get_signal_description(signum, signame)
      implicit none
      integer, intent(in) :: signum
      character(len=*), intent(out) :: signame
    end subroutine get_signal_description

    subroutine trap_segfault()
      implicit none
    end subroutine trap_segfault
  end interface



contains

  ! ---------------------------------------------------------
  subroutine messages_init()

    logical :: trap_signals

    call messages_obsolete_variable(global_namespace, 'DevelVersion', 'ExperimentalFeatures')

    !%Variable ExperimentalFeatures
    !%Type logical
    !%Default no
    !%Section Execution::Debug
    !%Description
    !% If true, allows the use of certain parts of the code that are
    !% still under development and are not suitable for production
    !% runs. This should not be used unless you know what you are doing.
    !% See details on
    !% <a href=http://octopus-code.org/experimental_features>wiki page</a>.
    !%End
    call parse_variable(global_namespace, 'ExperimentalFeatures', .false., conf%devel_version)

    call messages_obsolete_variable(global_namespace, 'DebugLevel', 'Debug')

    call debug_init(debug, global_namespace)

    warnings = 0
    experimentals = 0

    !%Variable DebugTrapSignals
    !%Type logical
    !%Default yes
    !%Section Execution::Debug
    !%Description
    !% If true, trap signals to handle them in octopus itself and
    !% print a custom backtrace. If false, do not trap signals; then,
    !% core dumps can be produced or gdb can be used to stop at the
    !% point a signal was produced (e.g. a segmentation fault).
    !%End
    call parse_variable(global_namespace, 'DebugTrapSignals', .true., trap_signals)

    if (trap_signals) call trap_segfault()

    call messages_reset_lines()

  end subroutine messages_init

  ! ---------------------------------------------------------

  subroutine messages_end()

    if(mpi_grp_is_root(mpi_world)) then

      if(experimentals > 0 .or. warnings > 0) then
        message(1) = ''
        call messages_info(1)
      end if


      if(warnings > 0) then
        call messages_write('Octopus emitted ')
        call messages_write(warnings)
        if(warnings > 1) then
          call messages_write(' warnings.')
        else
          call messages_write(' warning.')
        end if
        call messages_info()
      end if

      if(experimentals > 0) then
        call messages_new_line()
        call messages_write('Octopus used ')
        call messages_write(experimentals)
        if(experimentals > 1) then
          call messages_write(' experimental features:')
        else
          call messages_write(' experimental feature:')
        end if
        call messages_new_line()
        call messages_new_line()
        call messages_write('  Since you used one or more experimental features, results are likely')
        call messages_new_line()
        call messages_write('  wrong and should not  be considered as valid scientific data.  Check')
        call messages_new_line()
        call messages_new_line()
        call messages_write('  http://octopus-code.org/experimental_features')
        call messages_new_line()
        call messages_new_line()
        call messages_write('  or contact the octopus developers for details.')
        call messages_new_line()
        call messages_info()
      end if

      open(unit = iunit_out, file = 'exec/messages', action = 'write')
      write(iunit_out, '(a, i9)') "warnings          = ", warnings
      write(iunit_out, '(a, i9)') "experimental      = ", experimentals
      close(iunit_out)

    end if

  end subroutine messages_end

  ! ---------------------------------------------------------
  subroutine messages_fatal(no_lines, only_root_writes, namespace)
    integer,           optional, intent(in) :: no_lines
    logical,           optional, intent(in) :: only_root_writes
    type(namespace_t), optional, intent(in) :: namespace

    integer :: ii, no_lines_
    logical :: only_root_writes_, should_write
    integer, allocatable :: recv_buf(:), recv_req(:)
#ifdef HAVE_MPI
    integer, parameter :: FATAL_TAG = 32767
    logical :: received
    integer :: send_req
#endif

    no_lines_ = current_line
    if(present(no_lines)) no_lines_ = no_lines

    if(present(only_root_writes)) then
      should_write = mpi_grp_is_root(mpi_world) .or. (.not. only_root_writes)
      only_root_writes_ = only_root_writes
    else
      should_write = .true.
      only_root_writes_ = .false.
    end if

    ! This is to avoid all nodes reporting an error. The root node
    ! post a message reception to all nodes, the rest of the nodes
    ! send a message. If the message is received, the non-root nodes
    ! know that the root node will report the error, so they do not do
    ! anything.

    if(.not. only_root_writes_) then
      if(mpi_world%rank == 0) then

        allocate(recv_buf(1:mpi_world%size - 1))
        allocate(recv_req(1:mpi_world%size - 1))
        do ii = 1, mpi_world%size - 1
#ifdef HAVE_MPI
          call MPI_Recv_init(recv_buf(ii), 1, MPI_INTEGER, ii, FATAL_TAG, mpi_world%comm, recv_req(ii), mpi_err)
#endif
        end do
        deallocate(recv_buf)
        deallocate(recv_req)

      else

#ifdef HAVE_MPI
        call MPI_Send_init(1, 1, MPI_INTEGER, 0, FATAL_TAG, mpi_world%comm, send_req, mpi_err)
#endif
        !sleep for a second and check
        call loct_nanosleep(SLEEPYTIME_ALL, 0)
#ifdef HAVE_MPI
        call MPI_Test(send_req, received, MPI_STATUS_IGNORE, mpi_err)
#endif
        should_write = .false.

      end if
    end if

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_err, file='messages.stderr', &
        action='write', position='append')
    end if

    ! Give a moment for all standard output hopefully to be printed
    call loct_nanosleep(SLEEPYTIME_ALL, 0)

    ! If we are not writing wait for the root node to get here and
    ! write the error message. If the root doesn`t get here, we all print the
    ! error messsage anyways and die. Otherwise, no message might be written.
    if(.not. should_write) call loct_nanosleep(SLEEPYTIME_NONWRITERS, 0)

    call messages_print_stress(stderr, "FATAL ERROR")
    write(msg, '(a)') '*** Fatal Error (description follows)'
    call flush_msg(stderr, msg)

    if(present(namespace)) then
      if(len_trim(namespace%get()) > 0) then
        write(msg, '(3a)') '* In namespace ', trim(namespace%get()), ':'
        call flush_msg(stderr, msg)
      end if
    end if

#ifdef HAVE_MPI
    if(.not. only_root_writes_ .or. .not. mpi_grp_is_root(mpi_world)) then
      call flush_msg(stderr, shyphens)
      write(msg, '(a,i4)') "* From node = ", mpi_world%rank
      call flush_msg(stderr, msg)
    end if
#endif
    call flush_msg(stderr, shyphens)
    do ii = 1, no_lines_
      write(msg, '(a,1x,a)') '*', trim(message(ii))
      call flush_msg(stderr, msg)
    end do

    ! We only dump the stack in debug mode because subroutine invocations
    ! are only recorded in debug mode (via push_sub/pop_sub). Otherwise,
    ! it is a bit confusing that the stack seems to be empty.
    if(debug%trace) then
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
    end if

    if(flush_messages .and. mpi_grp_is_root(mpi_world)) then
      close(iunit_err)
    end if

    ! switch file indicator to state aborted
    call messages_switch_status('aborted')

#ifdef HAVE_MPI
    call MPI_Abort(mpi_world%comm, 999, mpi_err)
#endif

    call loct_exit_failure()
  end subroutine messages_fatal


  ! ---------------------------------------------------------
  subroutine messages_warning(no_lines, all_nodes, namespace)
    integer,           optional, intent(in) :: no_lines
    logical,           optional, intent(in) :: all_nodes
    type(namespace_t), optional, intent(in) :: namespace

    integer :: il, no_lines_
    logical :: have_to_write, all_nodes_

    no_lines_ = current_line
    if(present(no_lines)) no_lines_ = no_lines

    warnings = warnings + 1

    have_to_write = mpi_grp_is_root(mpi_world)

    all_nodes_ = .false.
    if(present(all_nodes)) then
      have_to_write = have_to_write .or. all_nodes
      all_nodes_ = all_nodes
    end if

    if(have_to_write) then

      if(flush_messages) open(unit=iunit_err, file='messages.stderr', action='write', position='append')

      call flush_msg(stderr, '')
      write(msg, '(a)') '** Warning:'
      call flush_msg(stderr, msg)

      if(present(namespace)) then
        if(len_trim(namespace%get()) > 0) then
          write(msg, '(3a)') '** In namespace ', trim(namespace%get()), ':'
          call flush_msg(stderr, msg)
        end if
      end if

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
      if(.not. present(verbose_limit) .or. debug%info) then
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

    if(.not. debug%info) return

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

    if(.not. debug%info) return
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

    if(.not. debug%info) return

    write(message(1), '(a,i3)') 'debug marker #', no
    call messages_debug(1)

  end subroutine messages_debug_marker


  ! ---------------------------------------------------------
  !> create status file for asynchronous communication
  subroutine messages_switch_status(status)
    character(len=*), intent(in) :: status

    ! only root node is taking care of file I/O
    if (.not.mpi_grp_is_root(mpi_world)) return

    ! remove old status files first, before we switch to a new state
    call loct_rm('exec/oct-status-running')
    call loct_rm('exec/oct-status-finished')
    call loct_rm('exec/oct-status-aborted')
    if (oct_status /= 'walltimer-aborted') then
      call loct_rm('exec/oct-status-walltimer-aborted')
    end if

    oct_status = status

    ! create empty status file to indicate new state
    open(unit=iunit_err, file='exec/oct-status-'//trim(status), &
      action='write', status='unknown')
    close(iunit_err)

  end subroutine messages_switch_status


  ! ---------------------------------------------------------
  subroutine open_debug_trace(iunit)
    integer, intent(out) :: iunit

    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir('debug')
    open(iunit, file = 'debug/debug_trace.node.'//filenum, &
      action='write', status='unknown', position='append')

  end subroutine open_debug_trace

  ! ---------------------------------------------------------
  subroutine delete_debug_trace()

    integer :: iunit
    character(len=6) :: filenum

    iunit = mpi_world%rank + unit_offset
    write(filenum, '(i6.6)') iunit - unit_offset
    call loct_mkdir('debug')
    call loct_rm('debug/debug_trace.node.'//filenum)

  end subroutine delete_debug_trace


  ! ---------------------------------------------------------
  subroutine alloc_error(size, file, line)
    integer(8),       intent(in) :: size
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    write(message(1), '(a,i18,3a,i5)') "Failed to allocate ", size, " words in file '", trim(file), "' line ", line
    call messages_fatal(1)

  end subroutine alloc_error


  ! ---------------------------------------------------------
  subroutine dealloc_error(size, file, line)
    integer(8),       intent(in) :: size
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    write(message(1), '(a,i18,3a,i5)') "Failed to deallocate array of ", size, " words in file '", trim(file), "' line ", line
    call messages_fatal(1)

  end subroutine dealloc_error


  ! ---------------------------------------------------------
  subroutine messages_input_error(namespace, var, details, row, column)
    type(namespace_t),          intent(in) :: namespace
    character(len=*),           intent(in) :: var
    character(len=*), optional, intent(in) :: details
    integer,          optional, intent(in) :: row
    integer,          optional, intent(in) :: column

    character(len=10) :: row_str, column_str

    call messages_write('Input error in the input variable '// trim(var))

    if (present(row)) then
      ! Print row and, if available, the column. We add one to both values
      ! in order to translate from the C numbering used by the parser to a
      ! more human-friendly numbering.
      write(row_str, '(I10)') row + 1
      call messages_write(' at row '//adjustl(row_str))
      if (present(column)) then
        write(column_str, '(I10)') column + 1
        call messages_write(', column '//adjustl(column_str))
      end if
    end if
    if(present(details)) then
      call messages_write(':', new_line = .true.)
      call messages_new_line()
      call messages_write('  '//trim(details))
    end if
    call messages_write('.', new_line = .true.)

    call messages_new_line()

    call messages_write('You can get the documentation of the variable with the command:', new_line = .true.)
    call messages_write('  oct-help -p '//trim(var))
    call messages_fatal(namespace=namespace)

  end subroutine messages_input_error
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuei(iunit, var, val)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer,          intent(in) :: val

    character(len=10) :: intstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    write(intstring,'(i10)') val
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(intstring))//']'

  end subroutine messages_print_var_valuei
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine messages_print_var_values(iunit, var, val)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: val

    if(.not. mpi_grp_is_root(mpi_world)) return

    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(val)//']'

  end subroutine messages_print_var_values
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuer(iunit, var, val, unit)
    integer,                intent(in) :: iunit
    character(len=*),       intent(in) :: var
    FLOAT,                  intent(in) :: val
    type(unit_t), optional, intent(in) :: unit

    character(len=10) :: floatstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    if(.not. present(unit)) then
      write(floatstring,'(g10.4)') val
      write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(floatstring))//']'
    else
      write(floatstring,'(g10.4)') units_from_atomic(unit, val)
      write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(adjustl(floatstring))//' '//trim(units_abbrev(unit))//']'
    end if

  end subroutine messages_print_var_valuer
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuel(iunit, var, val)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    logical,          intent(in) :: val

    character(len=3) :: lstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    if(val) then
      lstring = 'yes'
    else
      lstring = 'no'
    end if
    write(iunit,'(a)') 'Input: ['//trim(var)//' = '//trim(lstring)//']'

  end subroutine messages_print_var_valuel
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_valuear(iunit, var, val, unit)
    integer,                intent(in) :: iunit
    character(len=*),       intent(in) :: var
    FLOAT,                  intent(in) :: val(:)
    type(unit_t), optional, intent(in) :: unit

    integer :: ii
    character(len=10) :: floatstring

    if(.not. mpi_grp_is_root(mpi_world)) return

    call messages_write('Input: ['//trim(var)//' = (')
    do ii = 1, size(val)
      write(floatstring,'(g10.4)') val(ii)
      call messages_write(trim(adjustl(floatstring)))
      if(ii < size(val)) call messages_write(', ')
    end do
    call messages_write(')')
    if(present(unit)) then
      call messages_write(' '//trim(units_abbrev(unit))//']')
    else
      call messages_write(']')
    end if
    call messages_info(iunit = iunit)

  end subroutine messages_print_var_valuear
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine messages_print_var_info(iunit, var)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var

    if(.not. mpi_grp_is_root(mpi_world)) return

    call varinfo_print(iunit, var)
  end subroutine messages_print_var_info


  ! ---------------------------------------------------------
  subroutine messages_print_var_option_8(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer(8),       intent(in) :: option
    character(len=*), intent(in), optional :: pre

    integer :: option4

    option4 = int(option)

    if(.not. mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(pre)) then
      call varinfo_print_option(iunit, var, option4, pre)
      if(flush_messages) then
        call varinfo_print_option(iunit_out, var, option4, pre)
      end if
    else
      call varinfo_print_option(iunit, var, option4)
      if(flush_messages) then
        call varinfo_print_option(iunit_out, var, option4, pre)
      end if
    end if

    if(flush_messages) then
      close(iunit_out)
    end if
  end subroutine messages_print_var_option_8

  ! ---------------------------------------------------------
  subroutine messages_print_var_option_4(iunit, var, option, pre)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: var
    integer(4),       intent(in) :: option
    character(len=*), intent(in), optional :: pre

    call messages_print_var_option_8(iunit, var, int(option, 8), pre)

  end subroutine messages_print_var_option_4

  ! ---------------------------------------------------------
  subroutine messages_print_stress(iunit, msg, namespace)
    integer,                     intent(in) :: iunit
    character(len=*),  optional, intent(in) :: msg
    type(namespace_t), optional, intent(in) :: namespace

    integer, parameter :: max_len = 70

    integer :: ii, jj, length
    character(len=70) :: str
    character(len=max_len) :: msg_combined

    if(.not.mpi_grp_is_root(mpi_world)) return

    if(flush_messages) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    if(present(msg)) then
      ! make sure we do not get a segfault for too long messages
      if(len_trim(msg) > max_len) then
        msg_combined = trim(msg(1:max_len))
      else
        msg_combined = trim(msg)
      end if
      if(present(namespace)) then
        ! check if we are below the maximum length
        if(len_trim(msg) + len_trim(namespace%get()) + 1 < max_len) then
          ! only change message if namespace non-empty
          if(len_trim(namespace%get()) > 0) then
            msg_combined = trim(msg) // " " // trim(namespace%get())
          end if
        end if
      end if
      length = len_trim(msg_combined)

      str = ''
      jj = 1

      do ii = 1, (max_len - (length + 2))/2
        str(jj:jj) = '*'
        jj = jj + 1
      end do

      str(jj:jj) = ' '
      jj = jj + 1

      do ii = 1, length
        str(jj:jj) = msg_combined(ii:ii)
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
      if(iunit  ==  stderr) write(iunit_err, '(a)', advance=trim(adv_)) trim(str)
      if(iunit  ==  stdout) write(iunit_out, '(a)', advance=trim(adv_)) trim(str)
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
    if(usec2 - usec1  <  0) then
      usec2 = 1000000 + usec2
      if(sec2 >= sec1) then
        sec2 = sec2 - 1
      end if
    end if

    ! Replace values.
    if(sec2 >= sec1) then
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
    if(usec2 >= 1000000) then
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

    if(.not. debug%trace) return

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

    if(debug%trace_file) then
      call open_debug_trace(iunit)
      call push_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    end if

    if(debug%trace_term .and. mpi_grp_is_root(mpi_world)) then
      ! write to stderr if we are node 0
      call push_sub_write(stderr)
    end if

  contains

    subroutine push_sub_write(iunit_out)
      integer,  intent(in) :: iunit_out

      integer :: ii
      character(len=1000) :: tmpstr

      write(tmpstr,'(a,i6,a,i6.6,f20.6,i8,a)') "* I ", &
        sec, '.', usec, &
        loct_clock(), &
        loct_get_memory_usage() / 1024, " | "
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

    if(.not. debug%trace) return

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
    sub_name_short = trim(messages_clean_path(sub_name))

    if(sub_name_short /= sub_stack(no_sub_stack)) then
      write (message(1),'(a)') 'Wrong sub name on pop_sub :'
      write (message(2),'(2a)') ' got      : ', sub_name_short
      write (message(3),'(2a)') ' expected : ', sub_stack(no_sub_stack)
      call messages_fatal(3)
    end if

    if(debug%trace_file) then
      call open_debug_trace(iunit)
      call pop_sub_write(iunit)
      ! close file to ensure flushing
      close(iunit)
    end if

    if (debug%trace_term .and. mpi_grp_is_root(mpi_world)) then
      ! write to stderr if we are node 0
      call pop_sub_write(stderr)
    end if

    no_sub_stack = no_sub_stack - 1

  contains

    subroutine pop_sub_write(iunit_out)
      integer, intent(in) :: iunit_out

      integer :: ii
      character(len=1000) :: tmpstr

      write(tmpstr,'(a,i6,a,i6.6,f20.6,i8, a)') "* O ", &
        sec, '.', usec, &
        loct_clock() - time_stack(no_sub_stack), &
        loct_get_memory_usage() / 1024, " | "
      do ii = no_sub_stack - 1, 1, -1
        write(tmpstr,'(2a)') trim(tmpstr), "..|"
      end do
      write(tmpstr,'(2a)') trim(tmpstr), trim(sub_stack(no_sub_stack))
      call flush_msg(iunit_out, tmpstr)

    end subroutine pop_sub_write

  end subroutine pop_sub
#endif

  ! ---------------------------------------------------------
  subroutine messages_obsolete_variable(namespace, name, rep)
    type(namespace_t),          intent(in) :: namespace
    character(len=*),           intent(in) :: name
    character(len=*), optional, intent(in) :: rep

    if(parse_is_defined(namespace, trim(name))) then

      write(message(1), '(a)') 'Input variable '//trim(name)//' is obsolete.'

      if(present(rep)) then
        write(message(2), '(a)') ' '
        write(message(3), '(a)') 'Equivalent functionality can be obtained with the '//trim(rep)
        write(message(4), '(a)') 'variable. Check the documentation for details.'
        write(message(5), '(a)') '(You can use the `oct-help -p '//trim(rep)//'` command).'
        call messages_fatal(5, only_root_writes = .true.)
      else
        call messages_fatal(1, only_root_writes = .true.)
      end if

    end if

  end subroutine messages_obsolete_variable

  ! ---------------------------------------------------------
  subroutine messages_variable_is_block(namespace, name)
    type(namespace_t),          intent(in) :: namespace
    character(len=*),           intent(in) :: name

    if(parse_is_defined(namespace, trim(name))) then

      write(message(1), '(a)') 'Input variable `'//trim(name)//'` must be defined as a block.'
      write(message(2), '(a)') 'Please check the documentation for details.'
      write(message(3), '(a)') '(You can use the `oct-help -p '//trim(name)//'` command).'
      call messages_fatal(3, only_root_writes = .true.)

    end if

  end subroutine messages_variable_is_block

  ! ---------------------------------------------------------
  subroutine messages_experimental(name)
    character(len=*), intent(in) :: name

    experimentals = experimentals + 1

    if(.not. conf%devel_version) then
      call messages_write(trim(name)//' is an experimental feature.')
      call messages_new_line()
      call messages_new_line()
      call messages_write('If you still want to use this feature (at your own risk), check:')
      call messages_new_line()
      call messages_new_line()
      call messages_write('http://octopus-code.org/experimental_features')
      call messages_new_line()
      call messages_fatal(only_root_writes = .true.)
    else
      write(message(1), '(a)') trim(name)//' is under development.'
      write(message(2), '(a)') 'It might not work or produce wrong results.'
      call messages_warning(2)

      ! remove this warning from the count
      warnings = warnings - 1
    end if

  end subroutine messages_experimental


  !--------------------------------------------------------------
  subroutine messages_check_def(var, should_be_less, def, name, unit)
    FLOAT,                  intent(in) :: var
    logical,                intent(in) :: should_be_less
    FLOAT,                  intent(in) :: def
    character(len=*),       intent(in) :: name
    type(unit_t), optional, intent(in) :: unit

    logical :: is_bad
    character(len=3) :: op_str

    PUSH_SUB(messages_check_def)

    if(should_be_less) then
      is_bad = var > def
      op_str = ' > '
    else
      is_bad = var < def
      op_str = ' < '
    end if

    if(is_bad) then
      write(message(1), '(3a)') "The value for '", name, "' is inconsistent with the recommended value."
      if(present(unit)) then
        write(message(2), '(a,f8.3,4a,f8.3,a,a)') 'given ', units_from_atomic(unit, var), ' ', trim(units_abbrev(unit)), &
          op_str, 'recommended ', units_from_atomic(unit, def), ' ', trim(units_abbrev(unit))
      else
        write(message(2), '(a,f8.3,2a,f8.3)') 'given ', var, op_str, 'recommended ', def
      end if
      call messages_warning(2)
    end if

    POP_SUB(messages_check_def)
  end subroutine messages_check_def


  ! ------------------------------------------------------------
  subroutine messages_not_implemented(feature, namespace)
    character(len=*),            intent(in) :: feature
    type(namespace_t), optional, intent(in) :: namespace

    PUSH_SUB(messages_not_implemented)

    message(1) = trim(feature)//" not implemented."
    call messages_fatal(1, only_root_writes = .true., namespace=namespace)

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

    character(len=30) :: number
    FLOAT            :: tval

    tval = val
    if(present(units)) tval = units_from_atomic(units, val)

    if(present(fmt)) then
      write(number, '('//trim(fmt)//')') tval
    else
      write(number, '(f12.6)') tval
    end if

    if(optional_default(align_left, .false.)) then
      number = adjustl(number)
      number(1:len(number)) = ' '//number(1:len(number)-1)
    end if

    write(message(current_line), '(a, a)') trim(message(current_line)), trim(number)

    if(present(units) .and. optional_default(print_units, .true.)) then
      write(message(current_line), '(a, a, a)') trim(message(current_line)), ' ', trim(units_abbrev(units))
    end if

    if(optional_default(new_line, .false.)) call messages_new_line()

  end subroutine messages_write_float

  ! ------------------------------------------------------------

  subroutine messages_write_integer8(val, fmt, new_line, units, print_units)
    integer(8),                 intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line
    type(unit_t),     optional, intent(in) :: units
    logical,          optional, intent(in) :: print_units

    character(len=20) :: number
    FLOAT      :: val_conv_float

    if(present(units)) then
      val_conv_float = units_from_atomic(units, dble(val))

      if(present(fmt)) then
        write(message(current_line), '(a, '//trim(fmt)//')') trim(message(current_line)), val_conv_float
      else
        write(number, '(f15.3)') val_conv_float
        write(message(current_line), '(3a)') trim(message(current_line)), ' ', trim(adjustl(number))
      end if

    else

      if(present(fmt)) then
        write(message(current_line), '(a, '//trim(fmt)//')') trim(message(current_line)), val
      else
        write(number, '(i10)') val
        write(message(current_line), '(3a)') trim(message(current_line)), ' ', trim(adjustl(number))
      end if

    end if


    if(present(units) .and. optional_default(print_units, .true.)) then
      write(message(current_line), '(a, a, a)') trim(message(current_line)), ' ', trim(adjustl(units_abbrev(units)))
    end if

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_integer8

  ! ------------------------------------------------------------

  subroutine messages_write_integer(val, fmt, new_line, units, print_units)
    integer(4),                 intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line
    type(unit_t),     optional, intent(in) :: units
    logical,          optional, intent(in) :: print_units

    call messages_write_integer8(int(val, 8), fmt, new_line, units, print_units)

  end subroutine messages_write_integer

  ! ------------------------------------------------------------

  subroutine messages_write_str(val, fmt, new_line)
    character(len=*),           intent(in) :: val
    character(len=*), optional, intent(in) :: fmt
    logical,          optional, intent(in) :: new_line

    character(len=100) :: fmt_

    if(len(trim(message(current_line))) + len(trim(val)) > len(message(current_line))) then
      ! cannot use normal message approach without interfering with message we are trying to write
      ! write directly in case trim(val) is itself too long
      write(0, *) "Exceeded message line length limit, to write string:", trim(val)
    else
      fmt_ = optional_default(fmt, '(a)')
      write(message(current_line), '(a, '//trim(fmt_)//')') trim(message(current_line)), trim(val)
    end if

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_str

  ! ------------------------------------------------------------

  subroutine messages_write_logical(val, new_line)
    logical,           intent(in) :: val
    logical, optional, intent(in) :: new_line

    character(len=3) :: text

    if(val) then
      text = 'yes'
    else
      text = 'no'
    end if

    if(len(trim(message(current_line))) + len(trim(text)) > len(message(current_line))) then
      write(message(current_line + 1), '(3a)') "Exceeded message line length limit, to write logical value '", trim(text), "'"
      call messages_fatal(current_line + 1)
    end if

    write(message(current_line), '(a,1x,a)') trim(message(current_line)), trim(text)

    if(present(new_line)) then
      if(new_line) call messages_new_line()
    end if

  end subroutine messages_write_logical

  ! -----------------------------------------------------------

  character(len=MAX_PATH_LEN) function messages_clean_path(filename) result(clean_path)
    character(len=*), intent(in) :: filename

    integer :: pos, start

    pos = index(filename, 'src/', back = .true.)
    if(pos == 0) then
       ! 'src/' does not occur
       start = pos + 1
    else
       ! remove 'src/'
       start = pos + 4
    end if
    clean_path = filename(start:)
  end function messages_clean_path

  ! -----------------------------------------------------------

  subroutine messages_dump_stack(isignal)
    integer, intent(in) :: isignal

    integer :: ii
    character(len=300) :: description

    call get_signal_description(isignal, description)

    write(msg, '(a,i2)') ''
    call flush_msg(stderr, msg)
    write(msg, '(a,i2)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call flush_msg(stderr, msg)
    write(msg, '(a,i2)') ''
    call flush_msg(stderr, msg)
    write(msg, '(a,i2,2a)') '  Octopus was killed by signal ', isignal, ': ', trim(description)
    call flush_msg(stderr, msg)
    write(msg, '(a,i2)') ''
    call flush_msg(stderr, msg)
    write(msg, '(a)')    '  Note: Octopus is currently trapping signals. This might prevent the'
    call flush_msg(stderr, msg)
    write(msg, '(a)')    '  use of debuggers or the generation of core dumps. To change this'
    call flush_msg(stderr, msg)
    write(msg, '(a)')    '  behavior, use the DebugTrapSignals input option.'
    call flush_msg(stderr, msg)
    write(msg, '(a,i2)') ''
    call flush_msg(stderr, msg)
    write(msg, '(a,i2)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call flush_msg(stderr, msg)

    if(debug%trace) then
      call flush_msg(stderr, shyphens)

      write(msg, '(a)') 'Octopus debug trace: '
      call flush_msg(stderr, msg)
      do ii = 1, no_sub_stack
        write(msg, '(a,a)') ' > ', trim(sub_stack(ii))
        call flush_msg(stderr, msg, 'no')
      end do
      call flush_msg(stderr, " ")
    else
      write(msg, '(a)') " Octopus debug trace not available. You can enable it with 'Debug = trace'."
      call flush_msg(stderr, msg)
    end if

  end subroutine messages_dump_stack

end module messages_oct_m

! ---------------------------------------------------------
!> This subroutine is called by the assert macro, it is not in a
!> module so it can be called from any file. The interface is declared
!> in global_m.
subroutine assert_die(s, f, l)
  use messages_oct_m
  use mpi_oct_m

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

!-------------------------------------------------------

subroutine handle_segv(isignal) bind(c)
  use messages_oct_m
  use iso_c_binding

  implicit none

  integer(c_int), intent(in) :: isignal

  ! Switch status to aborted
  call messages_switch_status('aborted')

  ! Dump stack
  call messages_dump_stack(isignal)

end subroutine handle_segv


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
