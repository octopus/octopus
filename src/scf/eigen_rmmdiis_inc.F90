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
!
! This is an implementation of the RMM-DIIS method, it was taken from
! the GPAW code (revision: 2086 file: gpaw/eigensolvers/rmm_diis2.py).
!
! ---------------------------------------------------------
subroutine X(eigen_solver_rmmdiis) (gr, st, h, pre, tol, niter, converged, diff, verbose)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  type(preconditioner_t), intent(in) :: pre
  FLOAT,               intent(in)    :: tol
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged(:)
  FLOAT,     optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,   optional, intent(in)    :: verbose

  integer :: ik, ist, idim, ip
  R_TYPE, allocatable :: residuals(:, :, :), preres(:, :), resres(:, :)
  R_TYPE :: lambda
  FLOAT :: error

  call push_sub('eigen_rmmdiis_inc.eigen_solver_rmmdiss')

  ALLOCATE(residuals(1:NP, 1:st%d%dim, st%st_start:st%st_end), NP*st%d%dim*st%lnst)
  ALLOCATE(preres(1:NP_PART, 1:st%d%dim), NP_PART*st%d%dim*st%lnst)
  ALLOCATE(resres(1:NP, 1:st%d%dim), NP*st%d%dim*st%lnst)

  call X(subspace_diag)(gr, st, h, diff)

  do ik = 1, st%d%nik

    do ist = st%st_start, st%st_end
      
      call X(Hpsi)(h, gr, st%X(psi)(:,:, ist, ik) , residuals(:, :, ist), ist, ik)
      
      st%eigenval(ist, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:,:, ist, ik) , residuals(:, :, ist))

      do idim = 1, st%d%dim
        call lalg_axpy(NP, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), residuals(:, idim, ist))
      end do
      
    end do

    do ist = st%st_start, st%st_end

      ! dummy preconditioning
      do idim = 1, st%d%dim
        call lalg_copy(NP, residuals(:, idim, ist), preres(:, idim))
      end do
      
      ! calculate the residual of the residual
      call X(Hpsi)(h, gr, preres, resres, ist, ik)
      
      do idim = 1, st%d%dim
        call lalg_axpy(NP, -st%eigenval(ist, ik), preres(:, idim), resres(:, idim))
      end do

      ! the size of the correction
      
      lambda = -X(states_dotp)(gr%m, st%d%dim, residuals(:, :, ist), resres)/X(states_dotp)(gr%m, st%d%dim, resres, resres)

      do idim = 1, st%d%dim
        call lalg_axpy(NP, M_HALF*lambda, resres(:, idim), residuals(:, idim, ist))
      end do
      
      ! dummy preconditioning
      do idim = 1, st%d%dim
        call lalg_copy(NP, residuals(:, idim, ist), preres(:, idim))
      end do

      !now correct psi
      do idim = 1, st%d%dim
        call lalg_axpy(NP, M_TWO*lambda, preres(:, idim), st%X(psi)(:, idim, ist, ik))
      end do

    end do

    call X(states_gram_schmidt_full)(st, st%nst, gr%m, st%d%dim, st%X(psi)(:, :, :, ik))

    do ist = st%st_start, st%st_end
      
      call X(Hpsi)(h, gr, st%X(psi)(:,:, ist, ik) , residuals(:, :, ist), ist, ik)

      st%eigenval(ist, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:,:, ist, ik) , residuals(:, :, ist))
      
      diff(ist, ik) = X(states_residue)(gr%m, st%d%dim, residuals(:, :, ist), st%eigenval(ist, ik), st%X(psi)(:, :, ist, ik))
      
    end do

  end do

  call pop_sub()

end subroutine X(eigen_solver_rmmdiis)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
