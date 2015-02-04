!> @file
!! Manage different low-level operations
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_utils
  use dictionaries, only: f_err_throw,f_err_define
  use yaml_strings, only: yaml_toa
  implicit none

  public 

  integer, save, private :: INPUT_OUTPUT_ERROR

  !preprocessed include file with processor-specific values
  include 'f_utils.inc' !defines recl_kind


contains

  subroutine f_utils_errors()

    call f_err_define('INPUT_OUTPUT_ERROR',&
         'Some of intrinsic I/O fortan routines returned an error code',&
         INPUT_OUTPUT_ERROR,&
         err_action='Check if you have correct file system permission in i/o library or check the fortan runtime library')

  end subroutine f_utils_errors

  !> gives the maximum record length allowed for a given unit
  subroutine f_utils_recl(unt,recl_max,recl)
    implicit none
    integer, intent(in) :: unt !< unit to be checked for record length
    integer, intent(in) :: recl_max !< maximum value for record length
    !> Value for the record length. This corresponds to the minimum between recl_max and the processor-dependent value
    !! provided by inquire statement
    integer, intent(out) :: recl 
    !local variables
    logical :: unit_is_open
    integer :: ierr,ierr_recl
    integer :: recl_file

    !in case of any error, the value is set to recl_max
    recl=recl_max
    ierr_recl=-1
    !initialize the value of recl_file
    recl_file=int(-1234567891,kind=recl_kind)
    inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    if (ierr == 0 .and. .not. unit_is_open) then
       !inquire the record length for the unit
       inquire(unit=unt,recl=recl_file,iostat=ierr_recl)
    end if
    if (ierr_recl == 0) then
       recl=int(min(int(recl_max,kind=recl_kind),recl_file))
    end if
    if (recl <=0) recl=recl_max
  end subroutine f_utils_recl

  !> inquire for the existence of a file
  subroutine f_file_exists(file,exists)
    implicit none
    character(len=*), intent(in) :: file
    logical, intent(out) :: exists
    !local variables
    integer :: ierr 

    exists=.false.
    inquire(file=trim(file),exist=exists,iostat=ierr)
    if (ierr /=0) then
       call f_err_throw('Error in inquiring file='//&
         trim(file)//', iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
    end if
    exists = exists .and. ierr==0

  end subroutine f_file_exists

  !> call the close statement and retrieve the error
  !! do not close the file if the unit is not connected
  subroutine f_close(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    if (unit > 0) then
       close(unit,iostat=ierr)
       if (ierr /= 0) call f_err_throw('Error in closing unit='//&
               trim(yaml_toa(unit))//', iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
    end if
  end subroutine f_close

  !>search the unit associated to a filename.
  !! the unit is -1 if the file does not exists or if the file is
  !! not connected
  subroutine f_file_unit(file,unit)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: unit
    !local variables
    logical ::  exists
    integer :: unt,ierr

    unit=-1
    call f_file_exists(file,exists)
    if (exists) then
       inquire(file=trim(file),number=unt,iostat=ierr)
       if (ierr /= 0) then
          call f_err_throw('Error in inquiring file='//&
               trim(file)//' for number, iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          unit=unt
       end if
    end if
  end subroutine f_file_unit

  !>get a unit which is not opened at present
  !! start the search from the unit
  function f_get_free_unit(unit) result(unt2)
    implicit none
    !> putative free unit. Starts to search from this value
    integer, intent(in), optional :: unit
    integer :: unt2
    !local variables
    logical :: unit_is_open
    integer :: unt,ierr

    unit_is_open=.true.
    unt=7
    if (present(unit)) unt=unit
    do while(unit_is_open)
       unt=unt+1
       inquire(unit=unt,opened=unit_is_open,iostat=ierr)
       if (ierr /=0) then
          call f_err_throw('Error in inquiring unit='//&
               trim(yaml_toa(unt))//', iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
          exit
       end if
    end do
    unt2=unt
  end function f_get_free_unit

  !> delete an existing file. If the file does not exists, it does nothing
  subroutine f_delete_file(file)
    implicit none
    character(len=*), intent(in) :: file
    !local variables
    logical :: exists
    integer :: ierr
    external :: delete

    call f_file_exists(trim(file),exists)
    if (exists) then
       !c-function in utils.c
       call delete(trim(file),len_trim(file),ierr)
       if (ierr /=0) call f_err_throw('Error in deleting file='//&
            trim(file)//'iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    
  end subroutine f_delete_file

  !> get process id
  function f_getpid()
    implicit none
    integer :: f_getpid
    !local variables
    integer :: pid
    external :: getprocid

    call getprocid(pid)
    f_getpid=pid

  end function f_getpid

  !> rewind a unit
  subroutine f_rewind(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    rewind(unit,iostat=ierr)
    if (ierr /=0) call f_err_throw('Error in rewind unit='//&
         trim(yaml_toa(unit))//'iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
    
  end subroutine f_rewind

end module f_utils
