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

module io_m
  use datasets_m
  use global_m
  use mpi_m
  use messages_m
  use lib_oct_m
  use lib_oct_parser_m


  implicit none

  private
  public ::              &
    io_workpath,         &
    io_open,             &
    io_mkdir,            &
    io_init,             &
    io_end,              &
    io_status,           &
    io_dump_file,        &
    io_free,             &
    io_close,            &
    io_assign,           &
    io_get_extension,    &
    io_debug_on_the_fly


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
    !%Default no
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

    if(in_debug_mode) then
      !%Variable MPIDebugHook
      !%Type logical
      !%Default no
      !%Section Generalities::Debug
      !%Description
      !% When debugging the code in parallel it is usually difficult to find the origin
      !% of race conditions that appear in MPI communications. This variable introduces 
      !% a facility to control separate MPI processes. If set to yes, all nodes will 
      !% start up, but will get trapped in an endless loop. In every cycle of the loop 
      !% each node is sleeping for one second and is then checking if a file with the 
      !% name node_hook.xxx (where xxx denotes the node number) exists. A given node can 
      !% only be released from the loop if the corresponding file is created. This allows 
      !% to selectively run eg. a compute node first followed by the master node. Or, by
      !% reversing the file creation of the node hooks, to run the master first followed
      !% by a compute node.
      !%End
      call loct_parse_logical('MPIDebugHook', .false., mpi_debug_hook)
      if (mpi_debug_hook) then
        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* I ',sec,'.',usec,' | MPI debug hook'
        call write_debug(1)

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' In debug hook'
        write(node_hook,'(i3.3)') mpi_world%rank
        file_exists = .false.

        do while (.not.file_exists)
          inquire(file='node_hook.'//node_hook, exist=file_exists)
          call loct_nanosleep(1,0)
          write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, &
            ' - still sleeping. To release me touch: node_hook.'//trim(node_hook)
        end do

        write(stdout,'(a,i3,a)') 'node:', mpi_world%rank, ' Leaving debug hook'
        ! remove possible debug hooks
        call loct_rm( 'node_hook.'//trim(node_hook) )

        call loct_gettimeofday(sec, usec)
        call epoch_time_diff(sec,usec)
        write(message(1),'(a,i6,a,i6.6,20x,a)') '* O ', sec, '.', usec,' | MPI debug hook'
        call write_debug(1)
      end if
    end if

    ! create temporary dir (we will need it)
    !%Variable TmpDir
    !%Default "tmp/"
    !%Type string
    !%Section Generalities
    !%Description
    !% The name of the directory where octopus stores binary information
    !% like the wave-functions.
    !%End
    call loct_parse_string('TmpDir', 'tmp/', tmpdir)
    call io_mkdir(tmpdir, is_tmp=.true.)

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
  character(len=512) function io_workpath(path, is_tmp) result(wpath)
    character(len=*),  intent(in) :: path
    logical, optional, intent(in) :: is_tmp

    logical :: is_tmp_

    is_tmp_ = .false.
    if(present(is_tmp)) is_tmp_ = is_tmp

    if(path(1:1) .eq. '/') then
      ! we do not change absolute path names
      wpath = trim(path)
    else
      if(is_tmp_) then
        ! the current label is not added to the path for tmp files
        write(wpath, '(3a)') trim(work_dir), "/", trim(path)
      else
        write(wpath, '(4a)') trim(work_dir), "/", trim(current_label), trim(path)
      end if
    end if

  end function io_workpath


  ! ---------------------------------------------------------
  subroutine io_mkdir(fname, is_tmp)
    character(len=*),  intent(in) :: fname
    logical, optional, intent(in) :: is_tmp

    logical :: is_tmp_

    is_tmp_ = .false.
    if(present(is_tmp)) is_tmp_ = is_tmp

    call loct_mkdir(trim(io_workpath(fname, is_tmp_)))
  end subroutine io_mkdir


  ! ---------------------------------------------------------
  integer function io_open(file, action, status, form, position, die, is_tmp, recl) result(iunit)
    character(len=*), intent(in) :: file, action
    character(len=*), intent(in), optional :: status, form, position
    logical,          intent(in), optional :: die, is_tmp
    integer,          intent(in), optional :: recl

    character(len=20)  :: status_, form_, position_
    character(len=512) :: file_
    logical            :: die_, is_tmp_
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

    is_tmp_ = .false.
    if(present(is_tmp)) is_tmp_ = is_tmp
    file_ = io_workpath(file, is_tmp_)

    if(present(recl)) then
      open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
        recl=recl, action=trim(action), position=trim(position_), iostat=iostat)
    else
      open(unit=iunit, file=trim(file_), status=trim(status_), form=trim(form_), &
        action=trim(action), position=trim(position_), iostat=iostat)
    endif
      

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

    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      open(unit=iunit_out, file='messages.stdout', &
        action='write', position='append')
    end if

    do while(err == 0)
      read(iunit, fmt='(a80)', iostat=err) s
      if(err==0) then
        write(ounit, '(a)') s
        if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
          write(iunit_out, '(a)') s        
        endif
      end if
    end do
    
    if(flush_messages.and.mpi_grp_is_root(mpi_world)) then
      close(iunit_out)
    end if

    call io_close(iunit)

  end subroutine io_dump_file


  ! ---------------------------------------------------------
  ! Given a path, it returns the extension (if it exists) of the file
  ! (that is, the part of the name that comes after its last point)
  ! If the filename does not have an extension, it returns the empty string.
  character(len=8) function io_get_extension(path) result(ext)
    character(len = * ), intent(in)  :: path
    integer :: i, j

    i = index(path, ".", back = .true.)
    j = index(path(i+1:), "/")
    if(i.eq.0 .or. j.ne.0) then
      ext = ""
    else
      ext = path(i+1:)
    end if
  end function io_get_extension


  ! ---------------------------------------------------------
  ! check if debug mode or message flushing should be enabled or 
  ! disabled on the fly
  subroutine io_debug_on_the_fly()

    ! only root node performs the check
    if(.not.mpi_grp_is_root(mpi_world)) return

    if(io_file_exists('enable_debug_mode', 'Enabling DebugMode')) then
      conf%debug_level = 100
      in_debug_mode    = .true.
      ! this call does not hurt if the directory is already there
      ! but is otherwise required
      call io_mkdir('debug')
      ! we have been notified by the user, so we can cleanup the file
      call loct_rm('enable_debug_mode')
      ! artificially increase sub stack to avoid underflow
      no_sub_stack = no_sub_stack + 8
    endif
    if(io_file_exists('enable_flush_messages', 'Enabling flushing of messages')) then
      flush_messages   = .true.
      ! we have been notified by the user, so we can cleanup the file
      call loct_rm('enable_flush_messages')
    endif

    if(io_file_exists('disable_debug_mode', 'Disabling DebugMode')) then
      conf%debug_level = 0
      in_debug_mode    = .false.
      ! we have been notified by the user, so we can cleanup the file
      call loct_rm('disable_debug_mode')
    endif
    if(io_file_exists('disable_flush_messages', 'Disabling flushing of messages')) then
      flush_messages   = .false.
      ! we have been notified by the user, so we can cleanup the file
      call loct_rm('disable_flush_messages')
    endif

  end subroutine io_debug_on_the_fly

  
  ! Returns true if a file with name 'filename' exists
  ! and issues a reminder
  logical function io_file_exists(filename, msg) result(file_exists)
    character(len=*), intent(in)  :: filename, msg

    file_exists = .false.
    inquire(file=trim(filename), exist=file_exists)
    if(file_exists) then
      message(1) = trim(msg)
      call write_warning(1)
    end if

    return
  end function io_file_exists

end module io_m
