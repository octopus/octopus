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


! ---------------------------------------------------------
! This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(eigen_diagon_subspace)(gr, st, h, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  FLOAT, optional,     intent(out)   :: diff(1:st%nst,1:st%d%nik)

  R_TYPE, allocatable :: h_subspace(:,:), vec(:,:), f(:,:,:)
  integer             :: i, ik, tmp
  FLOAT               :: nrm2, ldiff(st%st_end-st%st_start+1)

  call push_sub('eigen_inc.Xeigen_diagon_subspace')

  ALLOCATE(h_subspace(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(vec(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(f(NP_PART, st%d%dim, st%st_start:st%st_end), NP_PART*st%d%dim*(st%st_end-st%st_start+1))

  ik_loop: do ik = 1, st%d%nik
    ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
    do i = st%st_start, st%st_end
      call X(hpsi)(h, gr, st%X(psi)(:, :, i, ik), f(:, :, i), ik)
    end do
    call states_blockt_mul(gr%m, st, st%X(psi)(:, :, :, ik), f, h_subspace)
    
    ! Diagonalize the hamiltonian in the subspace.
    call lalg_eigensolve(st%nst, h_subspace, vec, st%eigenval(:, ik))

    ! The new states are the given by the eigenvectors of the matrix.
    !$omp parallel workshare
    f(1:NP, 1:st%d%dim, st%st_start:st%st_end) = st%X(psi)(1:NP, 1:st%d%dim, st%st_start:st%st_end, ik)
    !$omp end parallel workshare

    call states_block_matr_mul(gr%m, st, f, vec, st%X(psi)(:, :, :, ik))

    ! Renormalize.
    do i = st%st_start, st%st_end
      nrm2 = X(states_nrm2)(gr%m, st%d%dim, st%X(psi)(:, :, i, ik))
      !$omp parallel workshare
      st%X(psi)(1:NP, 1:st%d%dim, i, ik) = st%X(psi)(1:NP, 1:st%d%dim, i, ik)/nrm2
      !$omp end parallel workshare
    end do

    ! Recalculate the residues if requested by the diff argument.
    if(present(diff)) then 
      do i = st%st_start, st%st_end
        call X(Hpsi)(h, gr, st%X(psi)(:, :, i, ik) , f(:, :, 1), ik)
        diff(i, ik) = X(states_residue)(gr%m, st%d%dim, f(:, :, 1), st%eigenval(i, ik), &
          st%X(psi)(:, :, i, ik))
      end do

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ldiff = diff(st%st_start:st%st_end, ik)
        call lmpi_gen_alltoallv(st%st_end-st%st_start+1, ldiff, tmp, diff(:, ik), st%mpi_grp)
      end if
#endif
    end if
  end do ik_loop

  deallocate(f, h_subspace, vec)

  call pop_sub()
end subroutine X(eigen_diagon_subspace) 


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
