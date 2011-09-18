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
!! $Id: eigen_rmmdiis_inc.F90 6342 2010-03-04 01:25:36Z dstrubbe $

! ---------------------------------------------------------
! See http://prola.aps.org/abstract/PRB/v54/i16/p11169_1
subroutine X(eigensolver_rmmdiis) (gr, st, hm, pre, tol, niter, converged, ik, diff, blocksize)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  FLOAT,                  intent(out)   :: diff(1:st%nst)
  integer,                intent(in)    :: blocksize

  R_TYPE, allocatable :: res(:, :, :, :)
  R_TYPE, allocatable :: psi(:, :, :, :)
  R_TYPE, allocatable :: mm(:, :, :, :), evec(:, :)
  R_TYPE, allocatable :: eigen(:)
  FLOAT,  allocatable :: eval(:)
  FLOAT, allocatable :: lambda(:)
  integer :: ist, jst, idim, ip, ii, iter, nops, maxst, ib, bsize
  R_TYPE :: ca, cb, cc
  R_TYPE, allocatable :: fr(:, :)
  type(profile_t), save :: prof
  type(batch_t), allocatable :: psib(:), resb(:)
  type(batch_t) :: psibit, resbit, stpsib
  integer, allocatable :: done(:), last(:)
  logical, allocatable :: failed(:)

  PUSH_SUB(X(eigensolver_rmmdiis))

  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim, 1:niter, 1:blocksize))
  SAFE_ALLOCATE(res(1:gr%mesh%np_part, 1:st%d%dim, 1:niter, 1:blocksize))
  SAFE_ALLOCATE(lambda(1:blocksize))
  SAFE_ALLOCATE(psib(1:niter))
  SAFE_ALLOCATE(resb(1:niter))
  SAFE_ALLOCATE(  done(1:blocksize))
  SAFE_ALLOCATE(  last(1:blocksize))
  SAFE_ALLOCATE(failed(1:blocksize))
  SAFE_ALLOCATE(fr(1:4, 1:blocksize))

  nops = 0

  call profiling_in(prof, "RMMDIIS")

  failed = .false.

  do jst = st%st_start, st%st_end, blocksize
    maxst = min(jst + blocksize - 1, st%st_end)
    bsize = maxst - jst + 1

    call batch_init(stpsib, st%d%dim, jst, maxst, st%X(psi)(:, :, jst:, ik))
    call batch_init(psib(1), st%d%dim, jst, maxst, psi(:, :, 1, :))

    call batch_copy_data(gr%mesh%np, stpsib, psib(1))

    call batch_end(stpsib)

    call batch_init(resb(1), st%d%dim, jst, maxst, res(:, :, 1, :))

    call X(hamiltonian_apply_batch)(hm, gr%der, psib(1), resb(1), ik)
    nops = nops + bsize

    SAFE_ALLOCATE(eigen(1:bsize))

    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1), resb(1), eigen)

    st%eigenval(jst:maxst, ik) = eigen(1:bsize)

    SAFE_DEALLOCATE_A(eigen)

    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(1), resb(1))

    call batch_end(psib(1))
    call batch_end(resb(1))

    done = 0

    do ist = jst, maxst
      ib = ist - jst + 1
      
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
      call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(1), psib(2))
      call X(hamiltonian_apply_batch)(hm, gr%der, psib(2), resb(2), ik)
      nops = nops + bsize - sum(done)

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        fr(1, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 2, ib), psi(:, :, 2, ib), reduce = .false.)
        fr(2, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1, ib), psi(:, :, 2, ib), reduce = .false.)
        fr(3, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 2, ib), res(:, :, 2, ib), reduce = .false.)
        fr(4, ib) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1, ib), res(:, :, 2, ib), reduce = .false.)
      end do

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, fr)

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        ca = fr(1, ib)*fr(4, ib) - fr(3, ib)*fr(2, ib)
        cb = fr(3, ib) - st%eigenval(ist, ik)*fr(1, ib)
        cc = st%eigenval(ist, ik)*fr(2, ib) - fr(4, ib)

        lambda(ib) = CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

        ! restrict the value of lambda to be between 0.1 and 1.0
        if(abs(lambda(ib)) > CNST(1.0)) lambda(ib) = lambda(ib)/abs(lambda(ib))
        if(abs(lambda(ib)) < CNST(0.1)) lambda(ib) = CNST(0.1)*lambda(ib)/abs(lambda(ib))
      end do

      do iter = 2, niter

        ! for iter == 2 the preconditioning was done already
        if(iter > 2) call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(iter - 1), psib(iter))

        do ist = jst, maxst
          ib = ist - jst + 1

          if(done(ib) /= 0) cycle

          ! predict by jacobi
          forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
            psi(ip, idim, iter, ib) = lambda(ib)*psi(ip, idim, iter, ib) + psi(ip, idim, iter - 1, ib)
          end forall
        end do

        ! calculate the residual
        call X(hamiltonian_apply_batch)(hm, gr%der, psib(iter), resb(iter), ik)
        nops = nops + bsize - sum(done)

        SAFE_ALLOCATE(  mm(1:iter, 1:iter, 1:2, 1:bsize))
        SAFE_ALLOCATE(evec(1:iter, 1:1))
        SAFE_ALLOCATE(eval(1:iter))

        do ist = jst, maxst
          ib = ist - jst + 1

          if(done(ib) /= 0 .and. .not. failed(ib)) cycle

          do idim = 1, st%d%dim
            call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter, ib), res(:, idim, iter, ib))
          end do

          ! perform the diis correction

          !  mm_1(i, j) = <res_i|res_j>
          call batch_init(resbit, st%d%dim, 1, iter, res(:, :, :, ib))
          call X(mesh_batch_dotp_matrix)(gr%mesh, resbit, resbit, mm(:, :, 1, ib), symm = .true., reduce = .false.)
          call batch_end(resbit)

          !  mm_2(i, j) = <psi_i|psi_j>
          call batch_init(psibit, st%d%dim, 1, iter, psi(:, :, :, ib))
          call X(mesh_batch_dotp_self)(gr%mesh, psibit, mm(:, :, 2, ib), reduce = .false.)
          call batch_end(psibit)

        end do

        if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, mm)

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
        call batch_end(resb(iter))
      end do

      do ist = jst, maxst
        ib = ist - jst + 1

        if(done(ib) /= 0) cycle

        if(.not. failed(ib)) then
          ! end with a trial move
          call X(preconditioner_apply)(pre, gr, hm, ik, res(:, :, iter - 1, ib), st%X(psi)(: , :, ist, ik))

          forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
            st%X(psi)(ip, idim, ist, ik) = psi(ip, idim, iter - 1, ib) + lambda(ib) * st%X(psi)(ip, idim, ist, ik)
          end forall
        else
          do idim = 1, st%d%dim
            call lalg_copy(gr%mesh%np, psi(:, idim, last(ib), ib), st%X(psi)(: , idim, ist, ik))
          end do
        end if

        if(mpi_grp_is_root(mpi_world)) then
          call loct_progress_bar(st%nst * (ik - 1) +  ist, st%nst*st%d%nik)
        end if

      end do

    end if
  end do

  call profiling_out(prof)

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  ! recalculate the eigenvalues and residuals

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end))

  do jst = st%st_start, st%st_end, blocksize
    maxst = min(jst + blocksize - 1, st%st_end)

    call batch_init(psib(1), st%d%dim, jst, maxst, st%X(psi)(:, :, jst:, ik))
    call batch_init(resb(1), st%d%dim, jst, maxst, res(:, :, 1, :))
    
    call X(hamiltonian_apply_batch)(hm, gr%der, psib(1), resb(1), ik)
    call X(mesh_batch_dotp_vector)(gr%der%mesh, psib(1), resb(1), eigen(jst:maxst))

    st%eigenval(jst:maxst, ik) = eigen(jst:maxst)

    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(1), resb(1))

    call X(mesh_batch_dotp_vector)(gr%der%mesh, resb(1), resb(1), eigen(jst:maxst))

    diff(jst:maxst) = sqrt(abs(eigen(jst:maxst)))

    call batch_end(psib(1))
    call batch_end(resb(1))

    nops = nops + maxst - jst + 1
  end do
  
  SAFE_DEALLOCATE_A(eigen)

  niter = nops

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(psib)
  SAFE_DEALLOCATE_A(resb)
  SAFE_DEALLOCATE_A(done)
  SAFE_DEALLOCATE_A(last)
  SAFE_DEALLOCATE_A(failed)
  SAFE_DEALLOCATE_A(fr)

  POP_SUB(X(eigensolver_rmmdiis))

end subroutine X(eigensolver_rmmdiis)

! ---------------------------------------------------------
subroutine X(eigensolver_rmmdiis_min) (gr, st, hm, pre, tol, niter, converged, ik, blocksize)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  integer,                intent(in)    :: blocksize

  integer, parameter :: sweeps = 5
  integer, parameter :: sd_steps = 2

  integer :: isd, ist, sst, est, bsize, ib, ii
  R_TYPE  :: ca, cb, cc
  R_TYPE, allocatable :: res(:, :, :), lambda(:)
  R_TYPE, allocatable :: kres(:, :, :)
  R_TYPE, allocatable :: me1(:, :), me2(:, :)

  type(batch_t) :: resb, kresb

  PUSH_SUB(X(eigensolver_rmmdiis_min))

  SAFE_ALLOCATE(me1(1:2, 1:st%d%block_size))
  SAFE_ALLOCATE(me2(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(res(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
  SAFE_ALLOCATE(kres(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
  SAFE_ALLOCATE(lambda(1:st%nst))

  niter = 0

  do ib = st%block_start, st%block_end
    sst = st%block_range(ib, 1)
    est = st%block_range(ib, 2)

    call batch_copy(st%psib(ib, ik), resb, reference= .false.)
    call batch_copy(st%psib(ib, ik), kresb, reference= .false.)

    do isd = 1, sd_steps

      call X(hamiltonian_apply_batch)(hm, gr%der, st%psib(ib, ik), resb, ik)

      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik), resb, me1(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik), st%psib(ib, ik), me1(2, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me1, (/2, st%d%block_size/))

      forall(ist = sst:est) st%eigenval(ist, ik) = me1(1, ist - sst + 1)/me1(2, ist - sst + 1)
 
      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), st%psib(ib, ik), resb)

      call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb, kresb)

      call X(hamiltonian_apply_batch)(hm, gr%der, kresb, resb, ik)

      niter = niter + 2*bsize

      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, kresb, me2(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik),  kresb, me2(2, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, resb,  me2(3, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik),  resb,  me2(4, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me2, (/4, st%d%block_size/))

      do ist = sst, est
        ii = ist - sst + 1

        ca = me2(1, ii) * me2(4, ii) - me2(3, ii) * me2(2, ii)
        cb = me1(2, ii) * me2(3, ii) - me1(1, ii) * me2(1, ii)
        cc = me1(1, ii) * me2(2, ii) - me1(2, ii) * me2(4, ii)

        lambda(ist) = CNST(2.0) * cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

      end do

      call batch_axpy(gr%mesh%np, lambda, kresb, st%psib(ib, ik))

    end do

    call batch_end(resb)
    call batch_end(kresb)

    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) +  est, st%nst*st%d%nik)
    end if

  end do

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(me1)
  SAFE_DEALLOCATE_A(me2)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(kres)

  POP_SUB(X(eigensolver_rmmdiis_min))

end subroutine X(eigensolver_rmmdiis_min)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
