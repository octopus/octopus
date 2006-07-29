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
! This routine diagonalises the hamiltonian in the subspace defined by the states.
subroutine X(eigen_diagon_subspace) (gr, st, h)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  
  R_TYPE, allocatable :: h_subspace(:,:), vec(:,:), f(:,:,:)
  integer :: ik, i, j

  ALLOCATE(h_subspace(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(vec(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(f(NP, st%d%dim, st%nst), NP*st%d%dim*st%nst)

  ik_loop: do ik = 1, st%d%nik
    eigenfunction_loop : do i = 1, st%nst
      call X(Hpsi)(h, gr, st%X(psi)(:,:, i, ik) , f(:,:, 1), ik)
      h_subspace(i, i) = st%eigenval(i, ik)
      do j = i+1, st%nst
        h_subspace(i, j) = X(states_dotp) (gr%m, st%d%dim, st%X(psi)(:,:, j, ik), f(:,:, 1))
        h_subspace(j, i) = R_CONJ(h_subspace(i, j))
      end do
    end do eigenfunction_loop

    call lalg_eigensolve(st%nst, h_subspace, vec, st%eigenval(:, ik))

    f(1:NP,1:st%d%dim,1:st%nst) = st%X(psi)(1:NP,1:st%d%dim,1:st%nst, ik)
    do i = 1, st%nst
      ! build new state
      st%X(psi)(1:NP,1:st%d%dim, i, ik) = vec(i, i)*st%X(psi)(1:NP,1:st%d%dim, i, ik)
      do j = 1, st%nst
        if(i.ne.j) st%X(psi)(1:NP,1:st%d%dim,i, ik) = st%X(psi)(1:NP,1:st%d%dim,i, ik) & 
             + vec(j, i)*f(1:NP,1:st%d%dim,j)
      end do

      ! renormalize
      st%X(psi)(1:NP,1:st%d%dim, i, ik) = &
           st%X(psi)(1:NP,1:st%d%dim, i, ik)/X(states_nrm2)(gr%m, st%d%dim, st%X(psi)(:,:, i, ik))
    end do
  end do ik_loop

  deallocate(f, h_subspace, vec)

end subroutine X(eigen_diagon_subspace) 
