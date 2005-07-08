!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  public :: io_workpath, io_open, io_mkdir
  public :: io_init, io_end, io_status, io_dump_file, io_free, io_close, io_assign
  public :: get_extension

  integer, parameter :: min_lun=10, max_lun=99
  logical            :: lun_is_free(min_lun:max_lun)
  character(len=512) :: work_dir    ! name of the ourput directory

contains


  ! ---------------------------------------------------------
  subroutine io_init()
    character(len=128) :: filename

    lun_is_free(min_lun:max_lun)=.true.

    stdin = 5

    ! setup standard output
    call loct_parse_string('stdout', '-', filename)
    stdout = 6
    if(trim(filename).ne.'-') then
       close(stdout)
       open(stdout, file=filename, status='unknown')
    end if

    ! setup standard error
    call loct_parse_string('stderr', '-', filename)
    stderr = 0
    if(trim(filename).ne.'-') then
       close(stderr)
       open(stderr, file=filename, status='unknown')
    end if

    ! check where to output files ...
    call loct_parse_string('WorkDir', '.', work_dir)
    ! ... and if necessary create workdir (will not harm if work_dir is already there)
    if (work_dir.ne.'.') call loct_mkdir(trim(work_dir))

    ! does the user want to flush stdout and stderr to files "messages.{stdout,stderr}" ?
    call loct_parse_logical('FlushMessages', .false., flush_messages)

    ! delete files so that we start writing to empty ones
    if(flush_messages) then
       call loct_rm('messages.stdout')
       call loct_rm('messages.stderr')
    endif

    ! verbosity level
    call loct_parse_int('Verbose', VERBOSE_NORMAL, conf%verbose)
    if(conf%verbose > VERBOSE_DEBUG .and. mpiv%node == 0) then
       call loct_parse_int('DebugLevel', 3, conf%debug_level)
       message(1) = 'Entering DEBUG mode'
       call write_warning(1)
    end if

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
          stop 'io_open'
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
          stop 'io_open'
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
    enddo
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
    endif
  end function get_extension


end module io
