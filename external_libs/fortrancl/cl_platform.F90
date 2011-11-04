  !! Copyright (C) 2010 X. Andrade
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
  !! $Id: cl.F90 3587 2007-11-22 16:43:00Z xavier $

#include "config_F90.h"
 
module cl_platform_m
  use cl_types_m

  implicit none 

  private

  ! the functions
  public ::                          &
    flGetPlatformIDs,                &
    flGetPlatformInfo

  interface flGetPlatformIDs

    subroutine flgetplatformids_num(num_platforms, status)
      use cl_types_m

      implicit none
      integer,              intent(out)  :: num_platforms
      integer,              intent(out)  :: status
    end subroutine flgetplatformids_num

    module procedure flgetplatformids_list

  end interface flGetPlatformIDs

  ! ---------------------------------------------------

  interface

    subroutine flGetPlatformInfo(platform, param_name, param_value, status)
      use cl_types_m

      implicit none
      type(cl_platform_id), intent(in)   :: platform
      integer,              intent(in)   :: param_name
      character(len=*),     intent(out)  :: param_value
      integer,              intent(out)  :: status
    end subroutine flGetPlatformInfo

  end interface

  ! ---------------------------------------------------

contains

  subroutine flgetplatformids_list(num_entries, platforms, num_platforms, status)
    integer,              intent(out)  :: num_entries
    type(cl_platform_id), intent(out)  :: platforms(:)
    integer,              intent(out)  :: num_platforms
    integer,              intent(out)  :: status

#ifdef HAVE_OPENCL
    integer                         :: iplatform
    type(cl_platform_id), allocatable :: plat(:)

    interface
      subroutine flgetplatformids_listall(num_entries, platforms, num_platforms, status)
        use cl_types_m

        implicit none

        integer,              intent(out)  :: num_entries
        type(cl_platform_id), intent(out)  :: platforms
        integer,              intent(out)  :: num_platforms
        integer,              intent(out)  :: status
      end subroutine flgetplatformids_listall

      subroutine flgetplatformids_getplat(allplatforms, iplatform, platform)
        use cl_types_m

        implicit none

        type(cl_platform_id), intent(in)   :: allplatforms
        integer,              intent(in)   :: iplatform
        type(cl_platform_id), intent(out)  :: platform
      end subroutine flgetplatformids_getplat
    end interface

    ! since our cl_platform_id type might be longer than the C
    ! cl_platform_id type we need to get all the values in an array
    ! and the copy them explicitly to the return array

    allocate(plat(1:num_entries))

    call flgetplatformids_listall(num_entries, plat(1), num_platforms, status)

    do iplatform = 1, num_platforms
      call flgetplatformids_getplat(plat(1), iplatform - 1, platforms(iplatform))
    end do

    deallocate(plat)
#endif
  end subroutine flgetplatformids_list

  ! ----------------------------------------------------------

end module cl_platform_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
