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
  integer :: times, ntimes
  R_TYPE, allocatable :: residuals(:, :), preres(:, :), resres(:, :)
  R_TYPE :: lambda
  FLOAT :: error

  call push_sub('eigen_rmmdiis_inc.eigen_solver_rmmdiss')

  ALLOCATE(residuals(1:NP, 1:st%d%dim), NP*st%d%dim)
  ALLOCATE(preres(1:NP_PART, 1:st%d%dim), NP_PART*st%d%dim*st%lnst)
  ALLOCATE(resres(1:NP, 1:st%d%dim), NP*st%d%dim*st%lnst)

  call X(subspace_diag)(gr, st, h, diff)

  ntimes = niter
  niter = 0

  do ik = 1, st%d%nik

    do ist = st%st_start, st%st_end
      
      do times = 1, ntimes

        call X(Hpsi)(h, gr, st%X(psi)(:,:, ist, ik) , residuals, ist, ik)

        st%eigenval(ist, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:,:, ist, ik) , residuals(:, :))

        do idim = 1, st%d%dim
          call lalg_axpy(NP, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), residuals(:, idim))
        end do

        error = X(states_nrm2)(gr%m, st%d%dim, residuals)

        if(error < tol .or. times == ntimes) then
          exit
        end if

        ! dummy preconditioning
        do idim = 1, st%d%dim
          call lalg_copy(NP, residuals(:, idim), preres(:, idim))
        end do

        ! calculate the residual of the residual
        call X(Hpsi)(h, gr, preres, resres, ist, ik)

        do idim = 1, st%d%dim
          call lalg_axpy(NP, -st%eigenval(ist, ik), preres(:, idim), resres(:, idim))
        end do

        ! the size of the correction
        lambda = -X(states_dotp)(gr%m, st%d%dim, residuals, resres)/X(states_dotp)(gr%m, st%d%dim, resres, resres)

        do idim = 1, st%d%dim
          call lalg_axpy(NP, M_HALF*lambda, resres(:, idim), residuals(:, idim))
        end do

        ! dummy preconditioning
        do idim = 1, st%d%dim
          call lalg_copy(NP, residuals(:, idim), preres(:, idim))
        end do

        !now correct psi
        do idim = 1, st%d%dim
          call lalg_axpy(NP, M_TWO*lambda, preres(:, idim), st%X(psi)(:, idim, ist, ik))
        end do

        niter = niter + 2

      end do

    end do

    call X(states_gram_schmidt_full)(st, st%nst, gr%m, st%d%dim, st%X(psi)(:, :, :, ik))

  end do

  call pop_sub()

end subroutine X(eigen_solver_rmmdiis)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
