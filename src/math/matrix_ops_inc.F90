!! Copyright (C) 2015 X. Andrade
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
!! $Id: matrix_inc.F90 14808 2015-11-21 05:43:46Z xavier $


subroutine X(to_submatrix)(this, bsize, proc, nproc, submatrix)
  type(matrix_t),           intent(in)    :: this
  integer,                  intent(in)    :: bsize(:)
  integer,                  intent(in)    :: proc(:)
  integer,                  intent(in)    :: nproc(:)
  R_TYPE,                   intent(inout) :: submatrix(:, :)

  integer :: ic1, ic2, ib1, ib2, i1, i2
  
  ic1 = 1
  do ib1 = 1 + bsize(1)*proc(1), this%dim(1), bsize(1)*nproc(1)
    do i1 = ib1, min(ib1 - 1 + bsize(1), this%dim(1))
      
      ic2 = 1
      do ib2 = 1 + bsize(2)*proc(2), this%dim(2), bsize(2)*nproc(2)
        do i2 = ib2, min(ib2 - 1 + bsize(2), this%dim(2))
          
          submatrix(ic1, ic2) = this%dmat(i1, i2)
          
          ic2 = ic2 + 1
        end do
      end do
      
      ic1 = ic1 + 1
    end do
  end do
  
end subroutine X(to_submatrix)

! ---------------------------------------------------------------

subroutine X(from_submatrix)(this, bsize, proc, nproc, submatrix)
  type(matrix_t),           intent(inout) :: this
  integer,                  intent(in)    :: bsize(:)
  integer,                  intent(in)    :: proc(:)
  integer,                  intent(in)    :: nproc(:)
  R_TYPE,                   intent(in)    :: submatrix(:, :)

  integer :: ic1, ic2, ib1, ib2, i1, i2
  
  ic1 = 1
  do ib1 = 1 + bsize(1)*proc(1), this%dim(1), bsize(1)*nproc(1)
    do i1 = ib1, min(ib1 - 1 + bsize(1), this%dim(1))
      
      ic2 = 1
      do ib2 = 1 + bsize(2)*proc(2), this%dim(2), bsize(2)*nproc(2)
        do i2 = ib2, min(ib2 - 1 + bsize(2), this%dim(2))
          
          this%dmat(i1, i2) = submatrix(ic1, ic2)
          
          ic2 = ic2 + 1
        end do
      end do
      
      ic1 = ic1 + 1
    end do
  end do
  
end subroutine X(from_submatrix)
  
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
