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
!! $Id: matrix.F90 7610 2011-03-25 03:41:01Z xavier $

#include "global.h"

module matrix_ops_m
  use global_m
  use lapack_m
  use matrix_m
  use messages_m
  use mpi_m
  use profiling_m
  use types_m

  implicit none

  private
  public ::                          &
    matrix_diagonalize_hermitian
  
contains

  ! ---------------------------------------------------------
  subroutine matrix_diagonalize_hermitian(this, eigenvalues, metric)
    type(matrix_t),           intent(inout) :: this            !< The matrix to diagonalize
    FLOAT,                    intent(out)   :: eigenvalues(:)  !< The eigenvalues of the matrix
    type(matrix_t), optional, intent(inout) :: metric          !< If present, a generalized eigenvalue problem is solved

    FLOAT :: worksize
    CMPLX :: zworksize
    FLOAT, allocatable :: work(:)
    CMPLX, allocatable :: zwork(:)
    integer :: info
    
    PUSH_SUB(matrix_diagonalize_hermitian)
    
    ASSERT(ubound(eigenvalues, dim = 1) == this%dim(1))
    ASSERT(this%dim(1) == this%dim(2))

    if(present(metric)) then
      ASSERT(all(this%dim(1:2) == metric%dim(1:2)))
    end if


    if(present(metric)) then
      
      if(matrix_type(this) == TYPE_FLOAT) then
        
        call dsygv(1, 'V', 'U', this%dim(1), this%dmat(1, 1), this%dim(1), metric%dmat(1, 1), metric%dim(1), &
          eigenvalues(1), worksize, -1, info)
        
        SAFE_ALLOCATE(work(1:int(worksize)))
        
        call dsygv(1, 'V', 'U', this%dim(1), this%dmat(1, 1), this%dim(1), metric%dmat(1, 1), metric%dim(1), &
          eigenvalues(1), work(1), int(worksize), info)
        
      else

        SAFE_ALLOCATE(work(1:3*this%dim(1) -2))
        
        call zhegv(1, 'V', 'U', this%dim(1), this%zmat(1, 1), this%dim(1), metric%zmat(1, 1), metric%dim(1), &
          eigenvalues(1), zworksize, -1, work(1), info)
        
        SAFE_ALLOCATE(zwork(1:int(zworksize)))
        
        call zhegv(1, 'V', 'U', this%dim(1), this%zmat(1, 1), this%dim(1), metric%zmat(1, 1), metric%dim(1), &
          eigenvalues(1), zwork(1), int(zworksize), work(1), info)
        
      end if

      SAFE_DEALLOCATE_A(work)
      
    end if
    
    POP_SUB(matrix_diagonalize_hermitian)
  end subroutine matrix_diagonalize_hermitian

#include "undef.F90"
#include "real.F90"

#include "matrix_ops_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "matrix_ops_inc.F90"

end module matrix_ops_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
