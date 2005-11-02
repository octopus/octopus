! --------------------------------------------------------------------------- !
!                                                                             !
! This module implements an interface to the FORTRAN logical unit             !
! system. Based on code by Richard Maine.                                     !
!                                                                             !
! Alberto Garcia, December 30, 1996                                           !
! Rewritten as a single subroutine                                            !
! with multiple entry points, March 7, 1998                                   !
!                                                                             !
! Converted into f90-module, A. Rubio 1999.                                   !
!                                                                             !
! Logical unit management. Units 0 to min_lun-1 are "reserved",               !
! since most of the "typical" files (output, etc) use them.                   !
!                                                                             !
! Logical units min_lun to min_max are managed by this module.                !
!                                                                             !
! --------------------------------------------------------------------------- !
!
! $Id$

#include "global.h"

module io
  use syslabels
  use global
  use messages
  use lib_oct
  use lib_oct_parser


  implicit none

  private
  public ::       &
    io_workpath,  &
    io_open,      &
    io_mkdir,     &
    io_init,      &
    io_end,       &
    io_status,    &
    io_dump_file, &
    io_free,      &
    io_close,     &
    io_assign,    &
    get_extension

  integer, parameter :: min_lun=10, max_lun=99
  logical            :: lun_is_free(min_lun:max_lun)
  character(len=512) :: work_dir    ! name of the output directory

contains

  ! ---------------------------------------------------------
  subroutine io_init()
    character(len=128) :: filename
    character(len=256) :: node_hook
    logical :: file_exists, mpi_debug_hook
    integer :: sec, usec


    lun_is_free(min_lun:max_lun)=.true.

    stdin = 5

    !%Variable stdout
    !%Type string
    !%Default "-"
    !%Section Generalities::IO
    !%Description
    !% The standard output by default goes to, well, to standard output. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call loct_parse_string('stdout', '-', filename)
    stdout = 6
    if(trim(filename).ne.'-') then
      close(stdout)
      open(stdout, file=filename, status='unknown')
    end if

    !%Variable stderr
    !%Type string
    !%Default "-"
    !%Section Generalities::IO
    !%Description
    !% The standard error by default goes to, well, to standard error. This can
    !% be changed by setting this variable: if you give it a name (other than "-")
    !% the output stream is printed in that file instead.
    !%End
    call loct_parse_string('stderr', '-', filename)
    stderr = 0
    if(trim(filename).ne.'-') then
      close(stderr)
      open(stderr, file=filename, status='unknown')
    end if

    !%Variable WorkDir
    !%Type string
    !%Default "."
    !%Section Generalities::IO
    !%Description
    !% By default, all files are written and read from the working directory,
    !% i.e. the directory from which the executable was launched. This behavior can
    !% be changed by setting this variable: if you give it a name (other than ".")
    !% the files are written and read in that directory.
    !%End
    call loct_parse_string('WorkDir', '.', work_dir)
    ! ... and if necessary create workdir (will not harm if work_dir is already there)
    if (work_dir.ne.'.') call loct_mkdir(trim(work_dir))

    !%Variable FlushMessages
    !%Type logical
    !%default no
    !%Section Generalities::IO
    !%Description
    !% In addition to writing to stdout and stderr, the code messages may also be
    !% flushed to "messages.stdout" and "messages.stderr", if this variable is
    !% set to yes.
    !%End
    call loct_parse_logical('FlushMessages', .false., flush_messages)

    ! delete files so that we start writing to empty ones
    if(flush_messages) then
      call loct_rm('messages.stdout')
      call loct_rm('messages.stderr')
    end if

    !%Variable DebugLevel
    !%Type integer
    !%Default 1
    !%Section Generalities::Debug
    !%Description
    !% This variable decides wether or not to enter debug-mode. In debugging mode,
    !% the program prints to standard error when it enters and exits the subroutines,
    !% what is the memory it is using (only, for the moment being, in Linux systems),
    !% and some other information. Useful for developers, and mandatory if you want
    !% to send a bug report to the developers and being considered.
    !% You have two options: (i) setting it to zero -- or less than zero, in which
    !% case you do not run in debugging mode (this is the default), or (ii) setting
    !% it to a positive number. In this case the entries and exits to nested subroutines
    !% are only printed down to the level that is given in this variable.
    !%End
    call loct_parse_int('DebugLevel',0,conf%debug_level)
    if(conf%debug_level>0) then
      in_debug_mode = .true.
    else
      in_debug_mode = .false.
    end if


    if(in_debug_mode) then
      call loct_parse_logical('MPIDebugHook', .false., mpi_debug_hook)
      if (mpi_debug_hook) then
        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* I ',sec,'.',usec,' | MPI debug hook'
        call write_debug(1)

        write(stdout,'(a,i3,a)') 'node:', mpiv%node, ' In debug hook'
        write(node_hook,'(i3.3)') mpiv%node
        file_exists = .false.

        do while (.not.file_exists)
          inquire(file='node_hook.'//node_hook, exist=file_exists)
          call loct_nanosleep(1,0)
          write(stdout,'(a,i3,a)') 'node:', mpiv%node, &
            ' - still sleeping. To release me touch: node_hook.'//trim(node_hook)
        end do

        write(stdout,'(a,i3,a)') 'node:', mpiv%node, ' Leaving debug hook'
        ! remove possible debug hooks
        call loct_rm( 'node_hook.'//trim(node_hook) )

        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* O ',sec,'.',usec,' | MPI debug hook'
        call write_debug(1)
      end if
    end if

    ! create temporary dir (we will need it)
    tmpdir = 'tmp/'
    call io_mkdir(tmpdir)

    ! create debug directory if in debugging mode
    if(in_debug_mode) then
      call io_mkdir('debug')
    end if

    !%Variable ProfilingMode
    !%Default no
    !%Type logical
    !%Section Generalities::Debug
    !%Description
    !% Use this variable to run octopus in profiling mode. In this mode
    !% octopus records time spent in certain areas of the code and
    !% the number of times this code is executed. These numbers
    !% are written in './profiling.NNN/profiling.nnn' with nnn being the
    !% node number (000 in serial) and NNN the number of processors.
    !% This is mainly for development purposes. Note, however, that
    !% octopus should be compiled with --disable-debug to do proper
    !% profiling.
    !%End
    call loct_parse_logical('ProfilingMode', .false., in_profiling_mode)

  end subroutine io_init


  ! ---------------------------------------------------------
  subroutine io_end()
    if(stderr.ne.0) call io_close(stderr)
    if(stdin .ne.5) call io_close(stdin)
    if(stdout.ne.6) call io_close(stdout)
  end subroutine io_end


  ! ---------------------------------------------------------
  subroutine io_assign(lun)
    integer, intent(out) :: lun

    integer :: iostat
    logical :: used

    lun = -1

    ! Looks for a free unit and assigns it to lun
    do lun = min_lun, max_lun
      if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)
        if (iostat .ne. 0) used = .true.
        lun_is_free(lun) = .false.
        if (.not. used) return
      end if
    end do

  end subroutine io_assign


  ! ---------------------------------------------------------
  subroutine io_free(lun)
    integer, intent(in) :: lun

    if (lun .ge. min_lun .and. lun .le. max_lun) &
      lun_is_free(lun) = .true.
  end subroutine io_free


  ! ---------------------------------------------------------
  character(len=512) function io_workpath(path) result(wpath)
    character(len=*), intent(in) :: path

    write(wpath, '(4a)') trim(work_dir), "/", trim(current_label), trim(path)

  end function io_workpath


  ! ---------------------------------------------------------
  subroutine io_mkdir(fname)
    character(len=*), intent(in) :: fname

    call loct_mkdir(trim(io_workpath(fname)))
  end subroutine io_mkdir


  ! ---------------------------------------------------------
  integer function io_open(file, action, status, form, position, die) result(iunit)
    character(len=*), intent(in) :: file, action
    character(len=*), intent(in), optional :: status, form, position
    logical,          intent(in), optional :: die

    character(len=20)  :: status_, form_, position_
    character(len=512) :: file_
    logical            :: die_
    integer            :: iostat

    status_ = 'unknown'
    if(present(status  )) status_   = status
    form_   = 'formatted'
    if(present(form    )) form_     = form
    position_ = 'asis'
    if(present(position)) position_ = position
    die_    = .true.
    if(present(die     )) die_      = die

    call io_assign(iunit)
    if(iunit<0) then
      if(die_) then
        write(stderr, '(a)') '*** IO Error: Too many files open'
        message(1) = 'Error: io_open.'
        call write_fatal(1)
      end if
      return
    end if

    if(file(1:1) .ne. '/') then
      file_ = io_workpath(file)
    else
      file_ = file
    end if

    open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
      action=trim(action), position=trim(position_), iostat=iostat)

    if(iostat.ne.0) then
      call io_free(iunit)
      iunit = -1
      if(die_) then
        write(*, '(5a)') '*** IO Error: Could not open file "', trim(file_), &
          '" for action="', trim(action), '"'
        message(1) = 'Error: io_open.'
        call write_fatal(1)
      end if
    end if

  end function io_open


  ! ---------------------------------------------------------
  subroutine io_close(iunit)
    integer, intent(inout) :: iunit

    close(iunit)
    call io_free(iunit)
    iunit = -1

  end subroutine io_close


  ! ---------------------------------------------------------
  ! Prints a list of the connected logical units and the names of
  ! the associated files
  ! ---------------------------------------------------------
  subroutine io_status(iunit)
    integer, intent(in) :: iunit

    integer :: i, iostat
    logical :: opened, named
    character(len=50) :: filename
    character(len=11) :: form

    write(iunit, '(a)') '******** io_status ********'
    do i = 0, max_lun
      inquire(i, opened=opened, named=named, name=filename, form=form, iostat=iostat)
      if (iostat .eq. 0) then
        if (opened) then
          if(.not.named) filename = 'No name available'
          write(iunit, '(i4,5x,a,5x,a)') i, form, filename
        end if
      else
        write(iunit, '(i4,5x,a)') i, 'Iostat error'
      end if
    end do
    write(iunit,'(a)') '********           ********'

  end subroutine io_status


  ! ---------------------------------------------------------
  subroutine io_dump_file(ounit, filename)
    integer,          intent(in) :: ounit
    character(len=*), intent(in) :: filename

    integer :: iunit, err
    character(len=80) :: s

    call io_assign(iunit)
    open(unit=iunit, file=filename, iostat=err, action='read', status='old')
    do while(err == 0)
      read(iunit, fmt='(a80)', iostat=err) s
      if(err==0) write(ounit, '(a)') s
    end do
    call io_close(iunit)

  end subroutine io_dump_file


  ! ---------------------------------------------------------
  ! Given a path, it returns the extension (if it exists) of the file
  ! (that is, the part of the name that comes after its last point)
  ! If the filename does not have an extension, it returns the empty string.
  ! ---------------------------------------------------------
  character(len=8) function get_extension(path) result(ext)
    character(len = * ), intent(in)  :: path
    integer :: i, j

    i = index(path, ".", back = .true.)
    j = index(path(i+1:), "/")
    if(i.eq.0 .or. j.ne.0) then
      ext = ""
    else
      ext = path(i+1:)
    end if
  end function get_extension

end module io
