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

  R_TYPE, allocatable :: res(:, :, :)
  R_TYPE, allocatable :: psi(:, :, :)

  R_TYPE, allocatable :: aa(:, :), mm(:, :), evec(:, :)
  FLOAT,  allocatable :: eval(:)

  FLOAT :: lambda
  integer :: ist, idim, ip, ii, jj, iter, nops
  R_TYPE :: ca, cb, cc, fr, rr, fhr, rhr
  logical :: fail, did_something

  call push_sub('eigen_rmmdiis_inc.eigensolver_rmmdiis')

  ALLOCATE(psi(gr%mesh%np_part, st%d%dim, niter), gr%mesh%np_part*st%d%dim*niter)
  ALLOCATE(res(gr%mesh%np_part, st%d%dim, niter), gr%mesh%np_part*st%d%dim*niter)

  nops = 0

  do ist = st%st_start, st%st_end

    do idim = 1, st%d%dim
      call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, ist, ik), psi(:, idim, 1))
    end do
    
    call X(hamiltonian_apply)(hm, gr, psi(:, :, 1), res(:, :, 1), ist, ik)
    nops = nops + 1

    st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 1))
    
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, 1), res(:, idim, 1))
    end do
    
    if(X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, 1)) < tol) cycle
    
    ! get lambda 
    call X(hamiltonian_apply)(hm, gr, res(:, :, 1), res(:, :, 2), ist, ik)
    nops = nops + 1

    rr = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, 1), res(:, :, 1))
    fr = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 1))
    rhr = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, 1), res(:, :, 2))
    fhr = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 2))
    
    ca = rr*fhr - rhr*fr
    cb = rhr - st%eigenval(ist, ik)*rr
    cc = fr - fhr
    
    lambda = 2*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))
    ! the article recommends to restrict the value of lambda, but for
    ! us it makes things worst.
    ! lambda = sign(max(min(1.0, abs(lambda)), 0.1), lambda)
    
    do iter = 2, niter
      call X(preconditioner_apply)(pre, gr, hm, res(:, :, iter - 1), psi(:, :, iter))

      ! predict by jacobi
      forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
        psi(ip, idim, iter) = lambda*psi(ip, idim, iter) + psi(ip, idim, iter - 1)
      end forall

      ! calculate the residual
      call X(hamiltonian_apply)(hm, gr, psi(:, :, iter), res(:, :, iter), ist, ik)
      nops = nops + 1

      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter), res(:, idim, iter))
      end do

      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, iter))

      ! perform the diis correction
      ALLOCATE(aa(iter, iter), iter**2)
      ALLOCATE(mm(iter, iter), iter**2)
      ALLOCATE(evec(iter, 1), iter)
      ALLOCATE(eval(iter), iter)

      do ii = 1, iter
        do jj = 1, iter
          aa(ii, jj) = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, ii), res(:, :, jj))
          mm(ii, jj) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, ii), psi(:, :, jj))
        end do
      end do

      fail = .false.
      call lalg_lowest_geneigensolve(1, iter, aa, mm, eval, evec, bof = fail)

      SAFE_DEALLOCATE_A(aa)
      SAFE_DEALLOCATE_A(mm)
      SAFE_DEALLOCATE_A(eval)      
        
      if(fail) then
        SAFE_DEALLOCATE_A(evec)
        exit
      end if
      
      ! generate the new vector and the new residual (the residual
      ! might be recalculated instead but that seems to be a bit
      ! slower).
      do idim = 1, st%d%dim
        call lalg_scal(gr%mesh%np, evec(iter, 1), psi(:, idim, iter))
        call lalg_scal(gr%mesh%np, evec(iter, 1), res(:, idim, iter))
      end do
      
      do ii = 1, iter - 1
        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, evec(ii, 1), psi(:, idim, ii), psi(:, idim, iter))
          call lalg_axpy(gr%mesh%np, evec(ii, 1), res(:, idim, ii), res(:, idim, iter))
        end do
      end do

      SAFE_DEALLOCATE_A(evec)
    end do

    ! end with a trial move
    forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
      st%X(psi)(ip, idim, ist, ik) = psi(ip, idim, iter - 1) + lambda*res(ip, idim, iter - 1)
    end forall

    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) +  ist, st%nst*st%d%nik)
    end if

  end do

  call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))

  ! recalculate the eigenvalues and residuals
  do ist = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, ist, ik), res(:, :, 1), ist, ik)
    nops = nops + 1
    
    st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), res(:, :, 1))
     
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), res(:, idim, 1))
    end do

    diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, 1))
  end do

  niter = nops

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
  R_TYPE, allocatable :: hres(:, :)

  FLOAT :: nrm2

  call push_sub('eigen_rmmdiis.Xeigensolver_rmmdiis_start')

  ALLOCATE(res(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)
  ALLOCATE(hres(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)
  ALLOCATE(kres(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)

  do isweep = 1, sweeps

    if(isweep > 1) call X(subspace_diag)(gr, st, hm, ik)

    do ist = st%st_start, st%st_end
      do isd = 1, sd_steps
        
        psi => st%X(psi)(:, :, ist, ik)

        call X(hamiltonian_apply)(hm, gr, psi, res, ist, ik)
        
        st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, psi, res)

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim), res(:, idim))
        end do
        
        call X(preconditioner_apply)(pre, gr, hm, res, kres)

        call X(hamiltonian_apply)(hm, gr, kres, hres, ist, ik)

        rr  = X(mf_dotp)(gr%mesh, st%d%dim, kres, kres)
        fr  = X(mf_dotp)(gr%mesh, st%d%dim, psi,  kres)
        rhr = X(mf_dotp)(gr%mesh, st%d%dim, kres, hres)
        fhr = X(mf_dotp)(gr%mesh, st%d%dim, psi,  hres)

        ca = rr*fhr - rhr*fr
        cb = rhr - st%eigenval(ist, ik)*rr
        cc = fr - fhr

        lambda = 2*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, lambda, kres(:, idim), psi(:, idim))
        end do

        ! normalize (todo: we can avoid this step by simply generalizing the formula for lambda)
        nrm2 = X(mf_nrm2)(gr%mesh, st%d%dim, psi)
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, M_ONE/nrm2, psi(:, idim))
        end do

      end do

    end do

    call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))
  end do
  
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(hres)
  SAFE_DEALLOCATE_A(kres)


  call pop_sub()

end subroutine X(eigensolver_rmmdiis_start)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
