!> @file
!! Include file used in yaml_strings.f90
!! Body of the yaml_toa template.
!! yaml: Yet Another Markup Language (ML for Human)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  character(len=max_value_length) :: str
  character(len=*), optional, intent(in) :: fmt

  !print *,'here',data,fmt
  str=repeat(' ',max_value_length)
  if (present(fmt)) then
     write(str,fmt) data
     !if the format has failed the result is full of stars
     if (trim(str) == repeat('*',len_trim(str))) write(str,cnv_fmt(data)) data
  else
     write(str,cnv_fmt(data)) data
  end if
  !otherwise write it in free format
  if (trim(str) == repeat('*',len_trim(str))) write(str,*) data
  !print *,'hereagain',str,data,fmt
  str=yaml_adjust(str,clean=.not. present(fmt))
