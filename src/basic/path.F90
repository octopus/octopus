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

module path_m

  use global_m
  use messages_m
  use profiling_m

  use loct_m, only: &
    loct_getcwd,    &
    loct_dirname,   &
    loct_basename,  &
    loct_realpath,  &
    loct_getenv

  implicit none

  private
  public ::          &
    path_isabs,      &
    path_join,       &
    path_expanduser, &
    path_basename,   &
    path_dirname,    &
    path_realpath,   &
    path_getcwd
  
  character, private, parameter :: sep='/'

contains

  ! ---------------------------------------------------------
  elemental function startswith(s, c) result(v)
    !intrinsic :: index, trim, adjustl
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: c
    !
    logical :: v
    !
    v=(index(trim(adjustl(s)), trim(adjustl(c)))==1)
    return
  end function startswith

  ! ---------------------------------------------------------
  elemental function endswith(s, c) result(v)
    !intrinsic :: index, len_trim, trim, adjustl
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: c
    !
    logical :: v
    !
    integer :: p
    !
    p=len_trim(adjustl(s))-len_trim(adjustl(c))+1
    v=(index(trim(adjustl(s)), trim(adjustl(c)), back=.true.)==p)
    return
  end function endswith

  ! ---------------------------------------------------------
  elemental function path_isabs(path) result(is)
    character(len=*), intent(in) :: path
    !
    logical :: is
    !
    is=startswith(path, sep)
    return
  end function path_isabs

  ! ---------------------------------------------------------
  elemental subroutine path_join(apath, bpath, rpath)
    !intrinsic :: trim, adjustl
    character(len=*), intent(in)  :: apath
    character(len=*), intent(in)  :: bpath
    character(len=*), intent(out) :: rpath
    !
    rpath=trim(adjustl(apath))
    if(path_isabs(bpath))then
       rpath=trim(adjustl(bpath))
    else if((trim(rpath)=="").or.endswith(rpath, sep))then
       rpath=trim(rpath)//trim(adjustl(bpath))
    else
       rpath=trim(rpath)//sep//trim(adjustl(bpath))
    end if
    return
  end subroutine path_join

  ! ---------------------------------------------------------
  subroutine path_expanduser(path, xpath)
    !intrinsic :: trim, adjustl
    character(len=*), intent(in)  :: path
    character(len=*), intent(out) :: xpath
    !
    character(len=MAX_PATH_LEN) :: home
    !
    PUSH_SUB(path_expanduser)
    xpath=trim(adjustl(path))
    if(startswith(path, "~/"))then
       call loct_getenv("HOME", home)
       call path_join(trim(adjustl(home)), path(3:), xpath)
    end if
    POP_SUB(path_expanduser)
    return
  end subroutine path_expanduser
 
  ! ---------------------------------------------------------
  subroutine path_basename(path, name)
    !intrinsic :: len, trim, adjustl
    character(len=*), intent(in)  :: path
    character(len=*), intent(out) :: name
    !
    PUSH_SUB(path_basename)
    call loct_basename(trim(adjustl(path)), name)
    name=trim(adjustl(name))
    POP_SUB(path_basename)
    return
  end subroutine path_basename
  
  ! ---------------------------------------------------------
  subroutine path_dirname(path, name)
    !intrinsic :: len, trim, adjustl
    character(len=*), intent(in)  :: path
    character(len=*), intent(out) :: name
    !
    PUSH_SUB(path_dirname)
    call loct_dirname(trim(adjustl(path)), name)
    name=trim(adjustl(name))
    POP_SUB(path_dirname)
    return
  end subroutine path_dirname

  ! ---------------------------------------------------------
  subroutine path_realpath(path, resolved_path)
    !intrinsic :: trim, adjustl
    character(len=*), intent(in)  :: path
    character(len=*), intent(out) :: resolved_path
    !
    character(len=MAX_PATH_LEN) :: xpth
    !
    PUSH_SUB(path_realpath)
    call path_expanduser(trim(adjustl(path)), xpth)
    call loct_realpath(trim(adjustl(xpth)), resolved_path)
    resolved_path=trim(adjustl(resolved_path))
    POP_SUB(path_realpath)
    return
  end subroutine path_realpath

  ! ---------------------------------------------------------
  subroutine path_getcwd(path)
    !intrinsic :: trim, adjustl
    character(len=*), intent(out) :: path
    !
    PUSH_SUB(path_getcwd)
    call loct_getcwd(path)
    path=trim(adjustl(path))
    POP_SUB(path_getcwd)
    return
  end subroutine path_getcwd

end module path_m

!! Local Variables:
!! mode: f90
!! End:
