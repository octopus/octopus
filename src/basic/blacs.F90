!! Copyright (C) 2011 D. Strubbe
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
!! $Id$

#include "global.h"

! -----------------------------------------------------------------------
!> This module contains interfaces for BLACS routines
!! Interfaces are from http://www.netlib.org/blacs/BLACS/QRef.html (entered manually...)
! -----------------------------------------------------------------------

module blacs_m

  implicit none

  integer, parameter :: BLACS_DLEN = 9

  interface
    subroutine blacs_get(icontxt, what, val)
      implicit none
      integer, intent(in)  :: icontxt
      integer, intent(in)  :: what
      integer, intent(out) :: val
    end subroutine blacs_get
  end interface

  interface
    subroutine blacs_gridinit(icontxt, order, nprow, npcol)
      implicit none
      integer,   intent(inout) :: icontxt
      character, intent(in)    :: order
      integer,   intent(in)    :: nprow
      integer,   intent(in)    :: npcol
    end subroutine blacs_gridinit
  end interface

  interface
    subroutine blacs_gridmap(icontxt, usermap, ldumap, nprow, npcol)
      implicit none
      integer, intent(inout) :: icontxt
      integer, intent(in)    :: usermap
      integer, intent(in)    :: ldumap
      integer, intent(in)    :: nprow
      integer, intent(in)    :: npcol
    end subroutine blacs_gridmap
  end interface

  interface
    subroutine blacs_gridexit(icontxt)
      implicit none
      integer, intent(in)  :: icontxt
    end subroutine blacs_gridexit
  end interface

  interface
    subroutine blacs_exit(icontxt)
      implicit none
      integer, intent(in)  :: icontxt
    end subroutine blacs_exit
  end interface

  interface
    subroutine blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)
      implicit none
      integer, intent(in)  :: icontxt
      integer, intent(out) :: nprow
      integer, intent(out) :: npcol
      integer, intent(out) :: myprow
      integer, intent(out) :: mypcol
    end subroutine blacs_gridinfo
  end interface


  interface
    integer function numroc(n, nb, iproc, isrcproc, nprocs)
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: nb
      integer, intent(in) :: iproc
      integer, intent(in) :: isrcproc
      integer, intent(in) :: nprocs
    end function numroc
  end interface

end module blacs_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
