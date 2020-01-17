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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module is intended to contain "only mathematical" functions
!! and procedures.
module sort_oct_m
  use global_oct_m
  use messages_oct_m

  implicit none

  private
  public ::                     &
    sort

  ! ------------------------------------------------------------------------------
  !> This is the common interface to a sorting routine.
  !! It performs the shell algorithm, not as fast as the quicksort for large numbers,
  !! but it seems that better for moderate numbers (around 100).
  !! Their possible interfaces are:
  !!   subroutine sort(a [, ind] )
  !!     FLOAT_OR_INTEGER, intent(inout) :: a(:)
  !!     integer, intent(inout), optional :: ind(:)
  !!     ! This routine sorts, from smallest to largest, the array a.
  !!     ! If the integer array ind is present, it puts in it the indexing
  !!     ! of the sorting, so that other arrays can be sorted according to
  !!     ! the sorting of a.
  !!   end subroutine sort
  !!
  !!   subroutine sort(a, x)
  !!     FLOAT, intent(inout) :: a(:)
  !!     FLOAT_OR_COMPLEX, intent(inout) :: x(:, : [, :])
  !!     ! This routine sorts, from smallest to largest, the array a.
  !!     ! The real or complex array x, which may be two or three dimensional,
  !!     ! is sorted according to the ordering of a. The last dimension of x
  !!     ! must have the same size as a.
  !!   end subroutine sort
  interface sort
    module procedure dsort, isort
    module procedure dshellsort1, zshellsort1, ishellsort1
    module procedure dshellsort2, zshellsort2, ishellsort2
  end interface sort

contains

  ! ---------------------------------------------------------
  subroutine dsort(a, ind)
    FLOAT,             intent(inout) :: a(:)
    integer, optional, intent(out)   :: ind(:)

    PUSH_SUB(dsort)

    if(size(a) > 0) then

      if(.not. present(ind)) then
        call dsort1(size(a), a(1))
      else
        call dsort2(size(a), a(1), ind(1))
      end if

    end if

    POP_SUB(dsort)
  end subroutine dsort


  ! ---------------------------------------------------------
  !> Shell sort for integer arrays.
  subroutine isort(a, ind)
    integer,           intent(inout) :: a(:)
    integer, optional, intent(out)   :: ind(:)

    PUSH_SUB(isort)

    if(size(a) > 0) then

      if(.not. present(ind)) then
        call isort1(size(a), a(1))
      else
        call isort2(size(a), a(1), ind(1))
      end if

    end if

    POP_SUB(isort)
  end subroutine isort


#include "undef.F90"
#include "complex.F90"
#include "sort_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "sort_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "sort_inc.F90"

end module sort_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
