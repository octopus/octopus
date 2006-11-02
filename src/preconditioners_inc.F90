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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! --------------------------------------------------------- 
subroutine X(preconditioner_apply)(pre, gr, h, a, b, omega)
  type(preconditioner_t), intent(in)    :: pre
  type(grid_t),           intent(inout) :: gr
  type(hamiltonian_t),    intent(inout) :: h
  R_TYPE,                 intent(inout) :: a(:,:)
  R_TYPE,                 intent(out)   :: b(:,:)
  R_TYPE,       optional, intent(in)    :: omega
  
  integer :: i, idim
  FLOAT, allocatable :: diag(:)

  select case(pre%which)
  case(PRECONDITIONER_NONE)
    do idim = 1, h%d%dim
      call lalg_copy(NP, a(:,idim), b(:,idim))
    end do

  case(PRECONDITIONER_SMOOTHING)
    do idim = 1, h%d%dim
      call X(nl_operator_operate) (pre%op, a(:, idim), b(:, idim))
    end do

  case(PRECONDITIONER_JACOBI)
    ALLOCATE(diag(NP), NP)

    do idim = 1, h%d%dim
      diag(:) = pre%diag_lapl(1:NP) + h%ep%vpsl(1:NP) + h%vhxc(1:NP, idim)
      if(present(omega)) then
        diag(:) = diag(:) + omega
      end if

      b(1:NP,idim) = a(1:NP,idim)/diag(1:NP)
    end do

    deallocate(diag)
  end select
  
end subroutine X(preconditioner_apply)
