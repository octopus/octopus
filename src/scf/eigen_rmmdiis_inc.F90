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
! the GPAW code (revision: 2086 file:
! gpaw/eigensolvers/rmm_diis2.py). I believe this is actually a two
! step DIIS. In our case we restart it several times to achieve
! convergency, perhaps it would be interesting to implement the full
! process.
!
! ---------------------------------------------------------
subroutine X(eigensolver_rmmdiis) (gr, st, hm, pre, tol, niter, converged, ik, diff, blocksize)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  type(preconditioner_t), intent(in) :: pre
  FLOAT,               intent(in)    :: tol
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: diff(1:st%nst)
  integer,             intent(in)    :: blocksize

  integer :: ist, idim, ib
  integer :: times, ntimes, iblock, psi_start, psi_end, num_in_block
  R_TYPE, allocatable :: residuals(:, :, :), preres(:, :, :), resres(:, :, :)
  R_TYPE, allocatable :: lambda(:, :)
  logical, allocatable :: conv(:)
#ifdef HAVE_MPI
  R_TYPE, allocatable :: lambda_tmp(:, :)
  FLOAT, allocatable :: ldiff(:)
  integer :: outcount
#endif  
  FLOAT :: error
  type(batch_t) :: psib, hpsib

  call push_sub('eigen_rmmdiis_inc.eigensolver_rmmdiss')

  ALLOCATE(residuals(1:gr%mesh%np_part, 1:st%d%dim, blocksize), gr%mesh%np*st%d%dim*blocksize)
  ALLOCATE(preres(1:gr%mesh%np_part, 1:st%d%dim, blocksize), gr%mesh%np_part*st%d%dim*blocksize)
  ALLOCATE(resres(1:gr%mesh%np, 1:st%d%dim, blocksize), gr%mesh%np*st%d%dim*blocksize)
  ALLOCATE(lambda(1:2, 1:blocksize), 2*blocksize) 

  ALLOCATE(conv(st%st_start:st%st_end), st%lnst)
  conv = .false.

  ntimes = niter
  niter = 0

  converged = 0

  iblock = 0
  do psi_start = st%st_start, st%st_end, blocksize
    iblock  = iblock + 1
    psi_end = min(psi_start + blocksize - 1, st%st_end)

    num_in_block = psi_end - psi_start + 1

    do times = 1, ntimes
      
      ! apply the hamiltonian over the initial vector

      call batch_init(psib, hm%d%dim, num_in_block)
      call batch_init(hpsib, hm%d%dim, num_in_block)

      ib = 0
      do ist = psi_start, psi_end
        ib = ib + 1
        if(conv(ist)) cycle
        call batch_add_state(psib, ist, st%X(psi)(:, :, ist, ik))
        call batch_add_state(hpsib, ist, residuals(:, :, ib))
      end do

      call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

      niter = niter + num_in_block
      
      call batch_end(psib)
      call batch_end(hpsib)
      
      ! calculate the residual

      ib = 0
      do ist = psi_start, psi_end
        ib = ib + 1
        if(conv(ist)) cycle

        st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:,:, ist, ik) , residuals(:, :, ib))
        
        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), residuals(:, idim, ib))
        end do

        error = X(mf_nrm2)(gr%mesh, st%d%dim, residuals(:, :, ib))

        if(error < tol) then
          num_in_block = num_in_block - 1
          conv(ist) = .true.
          converged = converged + 1
          diff(ist) = error
          cycle
        end if

        call  X(preconditioner_apply)(pre, gr, hm, residuals(:, :, ib), preres(:, :, ib))
      end do

      if(num_in_block == 0) exit

      ! apply the hamiltonian to the residuals
      call batch_init(psib, hm%d%dim, num_in_block)
      call batch_init(hpsib, hm%d%dim, num_in_block)

      ib = 0
      do ist = psi_start, psi_end
        ib = ib + 1
        if(conv(ist)) cycle
        call batch_add_state(psib, ist, preres(:, :, ib))
        call batch_add_state(hpsib, ist, resres(:, :, ib))
      end do

      call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)

      niter = niter + num_in_block
      
      call batch_end(psib)
      call batch_end(hpsib)

      ! calculate the correction
      ib = 0
      do ist = psi_start, psi_end
        ib = ib + 1
        if(conv(ist)) cycle

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), preres(:, idim, ib), resres(:, idim, ib))
        end do

        ! the size of the correction
        lambda(1, ib) = X(mf_dotp)(gr%mesh, st%d%dim, residuals(:, :, ib), resres(:, :, ib), reduce = .false.)
        lambda(2, ib) = X(mf_dotp)(gr%mesh, st%d%dim, resres(:, :, ib), resres(:, :, ib), reduce = .false.)
      end do

#ifdef HAVE_MPI
      if(gr%mesh%parallel_in_domains) then
        !reduce the two values together
        ALLOCATE(lambda_tmp(1:2, 1:blocksize), 2*blocksize) 
        call MPI_Allreduce(lambda, lambda_tmp, 2*blocksize, R_MPITYPE, MPI_SUM, gr%mesh%vp%comm, mpi_err)
        lambda(1:2, 1:blocksize) = lambda_tmp(1:2, 1:blocksize)
        deallocate(lambda_tmp)
      end if
#endif 

      ib = 0
      do ist = psi_start, psi_end
        ib = ib + 1
        if(conv(ist)) cycle
        
        lambda(1, ib) = -lambda(1, ib)/lambda(2, ib)

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, M_HALF*lambda(1, ib), resres(:, idim, ib), residuals(:, idim, ib))
        end do

        call X(preconditioner_apply)(pre, gr, hm, residuals(:, :, ib), preres(:, :, ib))

        !now correct psi
        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, M_TWO*lambda(1, ib), preres(:, idim, ib), st%X(psi)(:, idim, ist, ik))
        end do

        niter = niter + 1

      end do

    end do

    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik-1) +  psi_end - 1, st%nst*st%d%nik)
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ALLOCATE(ldiff(st%lnst), st%lnst)
    ldiff(1:st%lnst) = diff(st%st_start:st%st_end)
    call lmpi_gen_allgatherv(st%lnst, ldiff, outcount, diff, st%mpi_grp)
    deallocate(ldiff)
  end if
#endif

  call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))

  call pop_sub()

end subroutine X(eigensolver_rmmdiis)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
