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
    sort,                       &
    matrix_sort

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
    module procedure sort_complex
  end interface sort

  interface matrix_sort
    module procedure dmatrix_sort, zmatrix_sort, imatrix_sort
  end interface matrix_sort

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


  ! ---------------------------------------------------------
  !> Sort a complex vector vec(:)+i*Imvec(:) and put the ordering in reorder(:)
  !! according to the following order:
  !!
  !! 1. values with zero imaginary part sorted by increasing real part
  !! 2. values with negative imaginary part sorted by decreasing imaginary part
  !! 3. values with positive imaginary part unsorted
  subroutine sort_complex(vec, Imvec, reorder, imthr)
    FLOAT,           intent(inout)  :: vec(:)
    FLOAT,           intent(inout)  :: Imvec(:)
    integer,         intent(out)    :: reorder(:)
    FLOAT, optional, intent(in)     :: imthr !< the threshold for zero imaginary part

    integer              :: dim, n0, n1, n2, i
    integer, allocatable :: table(:),idx0(:)
    FLOAT,   allocatable :: temp(:),tempI(:)
    FLOAT                :: imthr_

    PUSH_SUB(sort_complex)

    dim = size(vec, 1)
    ASSERT(dim == size(Imvec,1) .and. dim == size(reorder,1))

    imthr_ = CNST(1E-6)
    if(present(imthr)) imthr_ = imthr

    allocate(table(1:dim))
    allocate(temp(1:dim))
    allocate(tempI(1:dim))

    tempI = Imvec
    call sort(tempI,table)

    n0 = 0
    n1 = 0
    temp = vec
    do i = 1, dim
      if (abs(Imvec(i)) < imthr_) then
        n0 = n0 + 1 
      else if (Imvec(i) < -imthr_) then 
        n1 = n1 + 1
      end if
      vec(i)     = temp(table(dim - i + 1))
      Imvec(i)   = tempI(dim - i + 1)
      reorder(i) = table(dim - i + 1)
    end do
    n2 = dim - n0 - n1 

    temp = vec
    tempI = Imvec
    table = reorder

    !first zero img parts
    if (n0 > 0) then
      allocate(idx0(1:n0))
      call sort(temp(n2+1:n2+n0),idx0(:))
    end if

    do i = 1, n0
      vec  (i) = temp (n2 + i)
      Imvec  (i) = tempI(n2 + idx0(i))
      reorder(i) = table(n2 + idx0(i))
    end do
    deallocate(idx0)

    !negative Img parts
    do i =  1, n1
      vec    (n0 + i) = temp (n2 + n0 + i)
      Imvec  (n0 + i) = tempI(n2 + n0 + i)
      reorder(n0 + i) = table(n2 + n0 + i)
    end do

    ! positive img parts
    do i = 1, n2
      vec    (n0 + n1 + i) = temp (n2 + 1 -i)
      Imvec  (n0 + n1 + i) = tempI(n2 + 1 -i)
      reorder(n0 + n1 + i) = table(n2 + 1 -i)
    end do

    deallocate(tempI)  
    deallocate(temp)
    deallocate(table)

    POP_SUB(sort_complex)
  end subroutine sort_complex


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
