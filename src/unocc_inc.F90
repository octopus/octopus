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
subroutine X(one_body) (m, st, h)
  type(mesh_t),        intent(in) :: m
  type(states_t),      intent(in) :: st
  type(hamiltonian_t), intent(in) :: h

  integer i, j, iunit
  R_TYPE :: me
  
  call io_assign(iunit)
  iunit = io_open('ME/1-body', action='write')

  do i = 1, st%nst
    do j = 1, st%nst
      if(j > i) cycle
      
      me = st%eigenval(i,1) - X(mf_integrate) (m, R_CONJ(st%X(psi) (:, 1, i, 1)) * &
           h%Vhxc(:, 1) * st%X(psi) (:, 1, j, 1))

      write(iunit, *) i, j, me
    end do
  end do
  
  call io_close(iunit)
end subroutine X(one_body)


! ---------------------------------------------------------
subroutine X(two_body) (gr, st)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st

  integer i, j, k, l, iunit
  R_TYPE :: me
  R_TYPE, allocatable :: n(:), v(:)

  call io_assign(iunit)
  iunit = io_open('ME/2-body', action='write')

  ALLOCATE(n(1:gr%m%np), gr%m%np)
  ALLOCATE(v(1:gr%m%np), gr%m%np)

  do i = 1, st%nst
    do j = 1, st%nst
      if(j > i) cycle

      n(:) = R_CONJ(st%X(psi) (:, 1, i, 1)) * st%X(psi) (:, 1, j, 1)
      call X(poisson_solve) (gr, v, n)

      do k = 1, st%nst
        if(k > i) cycle
        do l = 1, st%nst
          if(l > k) cycle
          if(l > j) cycle

          me = X(mf_integrate) (gr%m, v(:) * &
                st%X(psi) (:, 1, k, 1) *R_CONJ(st%X(psi) (:, 1, l, 1)))

          write(iunit, *) i, j, k, l, me
        end do
      end do
    end do
  end do

  call io_close(iunit)
end subroutine X(two_body)
