!! Copyright (C) 2009 X. Andrade
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

! See http://prola.aps.org/abstract/PRB/v54/i16/p11169_1

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

  R_TYPE, allocatable :: res(:, :, :, :)
  R_TYPE, allocatable :: psi(:, :, :, :)

  R_TYPE, allocatable :: mm(:, :, :, :), evec(:, :)
  FLOAT,  allocatable :: eval(:)

  FLOAT, allocatable :: lambda(:)
  integer :: ist, jst, idim, ip, ii, iter, nops, maxst, ib, bsize
  R_TYPE :: ca, cb, cc
  R_TYPE, allocatable :: fr(:, :)
  type(profile_t), save :: prof
  type(batch_t), allocatable :: psib(:), resb(:)
  type(batch_t) :: psibit, resbit
  integer, allocatable :: done(:), last(:)
  logical, allocatable :: failed(:)
#ifdef HAVE_MPI
  R_TYPE, allocatable :: fbuff(:)
  R_TYPE, allocatable :: buff(:, :)
  R_TYPE, allocatable :: mmc(:, :, :, :)
#endif

  call push_sub('eigen_rmmdiis_inc.eigensolver_rmmdiis')

  ALLOCATE(psi(gr%mesh%np_part, st%d%dim, niter, blocksize), gr%mesh%np_part*st%d%dim*blocksize*niter)
  ALLOCATE(res(gr%mesh%np_part, st%d%dim, niter, blocksize), gr%mesh%np_part*st%d%dim*blocksize*niter)
  ALLOCATE(lambda(blocksize), blocksize)
  ALLOCATE(psib(niter), niter)
  ALLOCATE(resb(niter), niter)
  ALLOCATE(done(blocksize), blocksize)
  ALLOCATE(last(blocksize), blocksize)
  ALLOCATE(failed(blocksize), blocksize)
  ALLOCATE(fr(4, blocksize), 4*blocksize)

  nops = 0

  call profiling_in(prof, "RMMDIIS")

  do jst = st%st_start, st%st_end, blocksize
    maxst = min(jst + blocksize - 1, st%st_end)
    bsize = maxst - jst + 1

    do ist = jst, maxst
      ib = ist - jst + 1
      do idim = 1, st%d%dim
        call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, ist, ik), psi(:, idim, 1, ib))
      end do
    end do

    call batch_init(psib(1), st%d%dim, jst, maxst, psi(:, :, 1, :))
    call batch_init(resb(1), st%d%dim, jst, maxst, res(:, :, 1, :))

    call X(hamiltonian_apply_batch)(hm, gr, psib(1), resb(1), ik)
    nops = nops + bsize

    call batch_end(psib(1))
    call batch_end(resb(1))

    done = 0

    do ist = jst, maxst
      ib = ist - jst + 1

      st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1, ib), res(:, :, 1, ib), reduce = .false.)
    end do

#ifdef HAVE_MPI
    ! perform all the reductions at once
    if(gr%mesh%parallel_in_domains) then
      ALLOCATE(fbuff(bsize), bsize)
      fbuff(1:bsize) = st%eigenval(jst:maxst, ik)
      call MPI_Allreduce(fbuff, st%eigenval(jst:maxst, ik), bsize, MPI_FLOAT, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
      SAFE_DEALLOCATE_A(fbuff)
    end if
#endif

    do ist = jst, maxst
      ib = ist - jst + 1

      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, 1, ib), res(:, idim, 1, ib))
      end do

      if(X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, 1, ib)) < tol) done(ib) = 1
    end do

    if(any(done(1:bsize) == 0)) then

      ! initialize the batch objects
      do iter = 1, niter

        call batch_init(psib(iter), st%d%dim, maxst - jst + 1 - sum(done))
        call batch_init(resb(iter), st%d%dim, maxst - jst + 1 - sum(done))

        do ist = jst, maxst
          ib = ist - jst + 1
          if(done(ib) /= 0) cycle

          call batch_add_state(psib(iter), ist, psi(:, :, iter, ib))
          call batch_add_state(resb(iter), ist, res(:, :, iter, ib))
        end do
      end do

      ! get lambda 
      call X(preconditioner_apply_batch)(pre, gr, hm, resb(1), psib(2))
      call X(hamiltonian_apply_batch)(hm, gr, psib(2), resb(2), ik)
      nops = nops + bsize - sum(done)

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        fr(1, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 2, ib), psi(:, :, 2, ib), reduce = .false.)
        fr(2, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1, ib), psi(:, :, 2, ib), reduce = .false.)
        fr(3, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 2, ib), res(:, :, 2, ib), reduce = .false.)
        fr(4, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1, ib), res(:, :, 2, ib), reduce = .false.)
      end do

#ifdef HAVE_MPI
      ! perform all the reductions at once
      if(gr%mesh%parallel_in_domains) then
        ALLOCATE(buff(4, blocksize), 4*blocksize)
        buff(1:4, 1:blocksize) = fr(1:4, 1:blocksize)
        call MPI_Allreduce(buff, fr, 4*blocksize, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(buff)
      end if
#endif

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        ca = fr(1, ib)*fr(4, ib) - fr(3, ib)*fr(2, ib)
        cb = fr(3, ib) - st%eigenval(ist, ik)*fr(1, ib)
        cc = st%eigenval(ist, ik)*fr(2, ib) - fr(4, ib)

        lambda(ib) = 2*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

        ! restrict the value of lambda to be between 0.1 and 1.0
        if(abs(lambda(ib)) > CNST(1.0)) lambda(ib) = lambda(ib)/abs(lambda(ib))
        if(abs(lambda(ib)) < CNST(0.1)) lambda(ib) = CNST(0.1)*lambda(ib)/abs(lambda(ib))
      end do

      do iter = 2, niter

        ! for iter == 2 the preconditioning was done already
        if(iter > 2) call X(preconditioner_apply_batch)(pre, gr, hm, resb(iter - 1), psib(iter))

        do ist = jst, maxst
          ib = ist - jst + 1

          if(done(ib) /= 0) cycle

          ! predict by jacobi
          forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
            psi(ip, idim, iter, ib) = lambda(ib)*psi(ip, idim, iter, ib) + psi(ip, idim, iter - 1, ib)
          end forall
        end do

        ! calculate the residual
        call X(hamiltonian_apply_batch)(hm, gr, psib(iter), resb(iter), ik)
        nops = nops + bsize - sum(done)

        ALLOCATE(mm(iter, iter, 2, bsize), iter**2*2*bsize)
        ALLOCATE(evec(iter, 1), iter)
        ALLOCATE(eval(iter), iter)

        do ist = jst, maxst
          ib = ist - jst + 1

          if(done(ib) /= 0 .and. .not. failed(ib)) cycle

          do idim = 1, st%d%dim
            call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter, ib), res(:, idim, iter, ib))
          end do

          ! perform the diis correction

          !  mm_1(i, j) = <res_i|res_j>
          call batch_init(resbit, st%d%dim, 1, iter, res(:, :, :, ib))
          call X(mf_dotp_batch)(gr%mesh, resbit, resbit, mm(:, :, 1, ib), symm = .true., reduce = .false.)
          call batch_end(resbit)

          !  mm_2(i, j) = <psi_i|psi_j>
          call batch_init(psibit, st%d%dim, 1, iter, psi(:, :, :, ib))
          call X(mf_dotp_batch)(gr%mesh, psibit, psibit, mm(:, :, 2, ib), symm = .true., reduce = .false.)
          call batch_end(psibit)

        end do

#ifdef HAVE_MPI
        ! perform all the reductions at once
        if(gr%mesh%parallel_in_domains) then
          ALLOCATE(mmc(iter, iter, 2, bsize), iter**2*2*bsize)
          mmc(1:iter, 1:iter, 1:2, 1:bsize) = mm(1:iter, 1:iter, 1:2, 1:bsize)
          call MPI_Allreduce(mmc, mm, iter**2*2*bsize, R_MPITYPE, MPI_SUM, gr%mesh%mpi_grp%comm, mpi_err)
          SAFE_DEALLOCATE_A(mmc)
        end if
#endif

        do ist = jst, maxst
          ib = ist - jst + 1

          if(done(ib) /= 0 .and. .not. failed(ib)) cycle

          failed(ib) = .false.
          call lalg_lowest_geneigensolve(1, iter, mm(:, :, 1, ib), mm(:, :, 2, ib), eval, evec, bof = failed(ib))
          if(failed(ib)) then
            last(ib) = iter - 1
            cycle
          end if

          ! generate the new vector and the new residual (the residual
          ! might be recalculated instead but that seems to be a bit
          ! slower).
          do idim = 1, st%d%dim
            call lalg_scal(gr%mesh%np, evec(iter, 1), psi(:, idim, iter, ib))
            call lalg_scal(gr%mesh%np, evec(iter, 1), res(:, idim, iter, ib))
          end do

          do ii = 1, iter - 1
            do idim = 1, st%d%dim
              call lalg_axpy(gr%mesh%np, evec(ii, 1), psi(:, idim, ii, ib), psi(:, idim, iter, ib))
              call lalg_axpy(gr%mesh%np, evec(ii, 1), res(:, idim, ii, ib), res(:, idim, iter, ib))
            end do
          end do

        end do

        SAFE_DEALLOCATE_A(mm)
        SAFE_DEALLOCATE_A(eval)      
        SAFE_DEALLOCATE_A(evec)

      end do

      do iter = 1, niter
        call batch_end(psib(iter))
        call batch_end(psib(iter))
      end do

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        if(.not. failed(ib)) then
          ! end with a trial move
          call X(preconditioner_apply)(pre, gr, hm, res(:, :, iter - 1, ib), st%X(psi)(: , :, ist, ik))

          forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
            st%X(psi)(ip, idim, ist, ik) = psi(ip, idim, iter - 1, ib) + lambda(ib)*st%X(psi)(ip, idim, ist, ik)
          end forall
        else
          do idim = 1, st%d%dim
            call lalg_copy(gr%mesh%np, psi(:, idim, last(ib), ib), st%X(psi)(: , idim, ist, ik))
          end do
        end if

        if(mpi_grp_is_root(mpi_world)) then
          call loct_progress_bar(st%nst*(ik - 1) +  ist, st%nst*st%d%nik)
        end if

      end do

    end if
  end do

  call profiling_out(prof)

  call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))

  ! recalculate the eigenvalues and residuals
  do jst = st%st_start, st%st_end, blocksize
    maxst = min(jst + blocksize - 1, st%st_end)

    call batch_init(psib(1), st%d%dim, jst, maxst, st%X(psi)(:, :, jst:, ik))
    call batch_init(resb(1), st%d%dim, jst, maxst, res(:, :, 1, :))

    call X(hamiltonian_apply_batch)(hm, gr, psib(1), resb(1), ik)

    call batch_end(psib(1))
    call batch_end(resb(1))

    do ist = jst, maxst
      ib = ist - jst + 1
      nops = nops + 1

      st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), res(:, :, 1, ib))

      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), res(:, idim, 1, ib))
      end do

      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, 1, ib))
    end do
  end do

  niter = nops

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(psib)
  SAFE_DEALLOCATE_A(resb)
  SAFE_DEALLOCATE_A(done)
  SAFE_DEALLOCATE_A(fr)

  call pop_sub()

end subroutine X(eigensolver_rmmdiis)

subroutine X(eigensolver_rmmdiis_start) (gr, st, hm, pre, tol, niter, converged, ik, diff, blocksize)
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

  integer, parameter :: sweeps = 5
  integer, parameter :: sd_steps = 2

  integer :: isweep, isd, ist, idim
  R_TYPE  :: lambda
  R_TYPE  :: ca, cb, cc, fr, rr, fhr, rhr

  R_TYPE, pointer     :: psi(:, :)
  R_TYPE, allocatable :: res(:, :)
  R_TYPE, allocatable :: kres(:, :)

  FLOAT :: fhf, ff

  call push_sub('eigen_rmmdiis.Xeigensolver_rmmdiis_start')

  ALLOCATE(res(gr%mesh%np_part, st%d%dim),  gr%mesh%np_part*st%d%dim)
  ALLOCATE(kres(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)
  niter = 0

  do isweep = 1, sweeps

    if(isweep > 1) call X(subspace_diag)(gr, st, hm, ik)

    do ist = st%st_start, st%st_end
      do isd = 1, sd_steps
        
        psi => st%X(psi)(:, :, ist, ik)

        call X(hamiltonian_apply)(hm, gr, psi, res, ist, ik)

        fhf = X(mf_dotp)(gr%mesh, st%d%dim, psi, res)
        ff = X(mf_dotp)(gr%mesh, st%d%dim, psi, psi)

        st%eigenval(ist, ik) = fhf/ff

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim), res(:, idim))
        end do
        
        call X(preconditioner_apply)(pre, gr, hm, res, kres)

        call X(hamiltonian_apply)(hm, gr, kres, res, ist, ik)
        niter = niter + 2

        rr  = X(mf_dotp)(gr%mesh, st%d%dim, kres, kres)
        fr  = X(mf_dotp)(gr%mesh, st%d%dim, psi,  kres)
        rhr = X(mf_dotp)(gr%mesh, st%d%dim, kres, res)
        fhr = X(mf_dotp)(gr%mesh, st%d%dim, psi,  res)

        ca = rr*fhr - rhr*fr
        cb = ff*rhr - fhf*rr
        cc = fhf*fr - ff*fhr

        lambda = 2*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, lambda, kres(:, idim), psi(:, idim))
        end do

      end do

    end do
    
    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) +  (ist*(isweep -1))/sweeps, st%nst*st%d%nik)
    end if

    call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))
  end do
  
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(kres)


  call pop_sub()

end subroutine X(eigensolver_rmmdiis_start)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
