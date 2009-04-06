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
  integer :: ist, idim, ip, size, ii, jj, iter
  R_TYPE :: ca, cb, cc, fr, rr, fhr, rhr
  logical :: fail

  call push_sub('eigen_rmmdiis_inc.eigensolver_rmmdiis')

  ALLOCATE(psi(gr%mesh%np_part, st%d%dim, niter), gr%mesh%np_part*st%d%dim*niter)
  ALLOCATE(res(gr%mesh%np_part, st%d%dim, niter), gr%mesh%np_part*st%d%dim*niter)

!  ALLOCATE(hpsi(gr%mesh%np, st%d%dim), gr%mesh%np*st%d%dim)
!  ALLOCATE(hres(gr%mesh%np, st%d%dim), gr%mesh%np*st%d%dim)
!  ALLOCATE(resres(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)

  do ist = st%st_start, st%st_end
!    print*, "EV ", ist, st%eigenval(ist, ik), diff(ist)
  end do

  do ist = st%st_start, st%st_end

    do idim = 1, st%d%dim
      call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, ist, ik), psi(:, idim, 1))
    end do

    do iter = 1, niter - 1

!      print*, ist, iter

      call X(hamiltonian_apply)(hm, gr, psi(:, :, iter), res(:, :, iter), ist, ik)
      
      if(iter == 1) st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 1))
      
      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter), res(:, idim, iter))
      end do


      if(X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, iter)) < tol) exit

      if(iter == 1) then
        ! get lambda
        call X(hamiltonian_apply)(hm, gr, res(:, :, 1), res(:, :, 2), ist, ik)

        rr = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, 1), res(:, :, 1))
        fr = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 1))
        rhr = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, 1), res(:, :, 2))
        fhr = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, 1), res(:, :, 2))

        ca = rr*fhr - rhr*fr
        cb = rhr - st%eigenval(ist, ik)*rr
        cc = fr - fhr
        
        lambda = 2*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))
!        print*, "COEFF", ca, cb, cc, cb**2 - CNST(4.0)*ca*cc
!        print*, "LAMBDA", lambda, 2*cc/(cb - sqrt(cb**2 - CNST(4.0)*ca*cc))
!        lambda = sign(max(min(1.0, abs(lambda)), 0.1), lambda)
!        print*, "LAMBDA2", lambda
      end if

      ! predict by jacobi
      forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
        psi(ip, idim, iter + 1) = psi(ip, idim, iter) + lambda*res(ip, idim, iter)
      end forall

      ! calculate the residual
      call X(hamiltonian_apply)(hm, gr, psi(:, :, iter + 1), res(:, :, iter + 1), ist, ik)
      
      do idim = 1, st%d%dim
        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter + 1), res(:, idim, iter + 1))
      end do

      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, iter + 1))

!      print*, "RES 0", diff(ist)

      ! perform the diis correction
      size = iter + 1
      ALLOCATE(aa(size, size), size**2)
      ALLOCATE(mm(size, size), size**2)
      ALLOCATE(evec(size, 1), size)
      ALLOCATE(eval(size), size)

      do ii = 1, size
        do jj = 1, size
          aa(ii, jj) = X(mf_dotp)(gr%mesh, st%d%dim, res(:, :, ii), res(:, :, jj))
          mm(ii, jj) = X(mf_dotp)(gr%mesh, st%d%dim, psi(:, :, ii), psi(:, :, jj))
        end do
      end do

      fail = .false.
      call lalg_lowest_geneigensolve(1, size, aa, mm, eval, evec, bof = fail)
      
      if(fail) then
        deallocate(aa, mm, eval, evec)
        exit
      end if

      !correct the new vector

      do idim = 1, st%d%dim
        call lalg_scal(gr%mesh%np, evec(size, 1), psi(:, idim, size))
        call lalg_scal(gr%mesh%np, evec(size, 1), res(:, idim, size))
      end do

      do ii = 1, size - 1
        do idim = 1, st%d%dim
          call lalg_axpy(gr%mesh%np, evec(ii, 1), psi(:, idim, ii), psi(:, idim, size))
          call lalg_axpy(gr%mesh%np, evec(ii, 1), res(:, idim, ii), res(:, idim, size))
        end do
      end do

      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, iter + 1))

!      print*, "ALPHA", evec(:, 1)
!      print*, "RES 1", diff(ist)

!      call X(hamiltonian_apply)(hm, gr, psi(:, :, iter + 1), res(:, :, iter + 1), ist, ik)
!      
!      do idim = 1, st%d%dim
!        call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), psi(:, idim, iter + 1), res(:, idim, iter + 1))
!      end do

!      diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, iter + 1))
!
!      print*, "RES 2", diff(ist)

      deallocate(aa, mm, eval, evec)

    end do

!    print*, "ITER",  niter, iter
    do idim = 1, st%d%dim
      call lalg_copy(gr%mesh%np, psi(:, idim, iter), st%X(psi)(:, idim, ist, ik))
    end do

  end do

  call X(states_gram_schmidt_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, ik))

  do ist = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, ist, ik), res(:, :, 1), ist, ik)
    
    st%eigenval(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), res(:, :, 1))
     
    do idim = 1, st%d%dim
      call lalg_axpy(gr%mesh%np, -st%eigenval(ist, ik), st%X(psi)(:, idim, ist, ik), res(:, idim, 1))
    end do

    diff(ist) = X(mf_nrm2)(gr%mesh, st%d%dim, res(:, :, 1))

!    print*, "EV ", st%eigenval(ist, ik), diff(ist)
  end do

!  stop
  call pop_sub()

end subroutine X(eigensolver_rmmdiis)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
