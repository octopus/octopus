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

module io
  use global

  implicit none

  integer, private, parameter :: min_lun=10, max_lun=99
  integer, private, parameter :: nunits=max_lun-min_lun+1
  integer, private   :: i, unit, lun, iostat
  
  logical, private   :: lun_is_free(min_lun:max_lun)=.true.
  logical, private   :: used, named, opened
  
  character(len=50), private :: filename
  character(len=11), private :: form

contains

  ! Next subroutines find or change the standard units...
  subroutine io_seterr(unit)
    integer, intent(inout) :: unit

    stderr = unit
  end subroutine io_seterr

  subroutine io_setout(unit)
    integer, intent(inout) :: unit

    stdout = unit
  end subroutine io_setout

  subroutine io_geterr(unit)
    integer, intent(inout) :: unit
    
    unit = stderr
  end subroutine io_geterr

  subroutine io_getout(unit)
    integer, intent(inout) :: unit
    
    unit = stdout
  end subroutine io_getout

  ! Logical unit management
  subroutine io_assign(lun)
    integer, intent(INOUT) :: lun

    ! Looks for a free unit and assigns it to lun
    do lun = min_lun, max_lun
      if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)
        if (iostat .ne. 0) used = .true.
        lun_is_free(lun) = .false.
        if (.not. used) return
      end if
    end do
    message(1) = 'No luns available in io_assign'
    call write_fatal(1)

  end subroutine io_assign

  ! Useful to specify that one needs to use a particular unit number
  ! For example, assume some legacy code expects to work with unit 15:
  !
  ! call io_reserve(15)   ! this call at the beginning of the program
  ! ...
  ! open(15,....)
  subroutine io_reserve(lun)
    integer, intent(INOUT) :: lun

    inquire(unit=lun, opened=used, iostat=iostat)
    if (iostat .ne. 0) used = .true.
    if (used) then
      write(message(1), '(a,i3,a)') &
           'Cannot reserve unit', lun, '. Already connected!'
      call write_fatal(1)
    end if
    if (lun .ge. min_lun .and. lun .le. max_lun) &
         lun_is_free(lun) = .false.
    
  end subroutine io_reserve

  ! Use this routine instead of a simple close!!
  subroutine io_close(lun)
    integer, intent(INOUT) :: lun

    close(lun)
    if (lun .ge. min_lun .and. lun .le. max_lun) &
         lun_is_free(lun) = .true.
  end subroutine io_close

  ! Prints a list of the connected logical units and the names of
  ! the associated files
  subroutine io_status
    message(1) = '******** io_status ********'
    call write_info(1)
    do i = 0, max_lun
      inquire(i, opened=opened, named=named, name=filename,          &
           form=form, iostat=iostat)
      if (iostat .eq. 0) then
        if (opened) then
          if (named) then
            write(message(1), 9000) i, form, filename
          else
            write(message(1), 9000) i, form, 'No name available'
          end if
        end if
      else
        write(message(1), 9000) i, 'Iostat error'
      end if
      call write_info(1)
    enddo
    message(1) = '********           ********'
    call write_info(1)

9000 format(i4,5x,a,5x,a)
    
  end subroutine io_status

end module io
