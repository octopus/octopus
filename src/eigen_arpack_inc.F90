!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

subroutine X(eigen_solver_arpack)(gr, st, h, tol_, niter, ncv, converged, diff)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  FLOAT,               intent(in)    :: tol_
  integer,             intent(inout) :: niter
  integer,             intent(in)    :: ncv
  integer,             intent(inout) :: converged
  FLOAT,     optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)

  logical, allocatable :: select(:)
  FLOAT, allocatable :: ax(:), d(:, :), resid(:), v(:, :),   &
    workd(:), workev(:), workl(:)
  integer :: ldv, nev, iparam(11), ipntr(14), ido, n, lworkl, info, ierr, &
    i, j, ishfts, maxitr, mode1, ik
  FLOAT :: tol, sigmar, sigmai

!!!!WARNING: No support for spinors, yet. No support for complex wavefunctions.
  call push_sub('eigen_arpack.eigen_solver_arpack')

  kpoints: do ik = 1, st%d%nik

#if defined(HAVE_MPI)
    if(gr%m%parallel_in_domains) then
      message(1) = 'Error: Arpack-Solver not parallelized for domain decomposition.'
      call write_fatal(1)
      !  FIXME: Need to adjust m%x and m%vol_pp occurences in the code below
      !         appropriately for domain decomposition. Also parallelization
      !         of the vectors has to be taken care of.
    end if
#endif

    ldv = NP
    n = NP
    nev = st%nst
    lworkl  = 3*ncv**2+6*ncv
    ALLOCATE(ax(ldv),       ldv)
    ALLOCATE(d(ncv, 3),     ncv*3)
    ALLOCATE(resid(ldv),    ldv)
    ALLOCATE(v(ldv, ncv),   ldv*ncv)
    ALLOCATE(workd(3*ldv),  3*ldv)
    ALLOCATE(workev(3*ncv), 3*ncv)
    ALLOCATE(workl(lworkl), lworkl)
    ALLOCATE(select(ncv),   ncv)

    select = .true.
    tol    = tol_
    ido    = 0
    info = 1

    do i = 1, NP
      resid(i) = sum(st%X(psi)(i, 1, 1:st%nst, ik))*sqrt(gr%m%vol_pp(i))
    end do

    ishfts = 1
    maxitr = niter
    mode1 = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1

    do
      call dnaupd ( ido, 'I', n, 'SR', nev, tol, resid, ncv, &
        v, ldv, iparam, ipntr, workd, workl, lworkl, info )
      if( abs(ido).ne.1) exit
      call av (n, workd(ipntr(1)), workd(ipntr(2)))
    end do
    ! If info is larger than zero, it may not be an error (i.e., not all eigenvectors
    ! were converged)
    if(info .lt. 0) then
      write(message(1),'(a,i5)') 'Error with ARPACK _naupd, info = ', info
      write(message(2),'(a)')    'Check the documentation of _naupd.'
      call write_fatal(2)
    end if

    call dneupd ( .true., 'A', select, d, d(1,2), v, ldv, &
      sigmar, sigmai, workev, 'I', n, 'SR', nev, tol, &
      resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, ierr )

    if(ierr .ne. 0) then
      write(message(1),'(a,i5)') 'Error with ARPACK _neupd, info = ', info
      write(message(2),'(a)')    'Check the documentation of _neupd.'
      call write_fatal(2)
    end if

    ! This sets the number of converged eigenvectors.
    converged =  iparam(5)
    !call dmout(6, converged, 3, d, ncv, -6, 'Ritz values (Real, Imag) and residual residuals')
    ! This sets niter to the number of matrix-vector operations.
    niter = iparam(9)
    do j = 1, min(st%nst, converged)
      write(*, *) 'Modifying the wavefunctions...'
      do i = 1, NP
        st%X(psi)(i, 1, j, ik) = v(i, j)/sqrt(gr%m%vol_pp(i))
      end do
      st%eigenval(j, ik) = d(j, 1)
      if(workl(ipntr(11)+j-1)<CNST(1.0e-99)) then
        diff(j, ik) = M_ZERO
      else
        diff(j, ik) = workl(ipntr(11)+j-1)
      end if
    end do

    deallocate(ax, d, resid, v, workd, workev, workl, select)

  end do kpoints

  call pop_sub()
contains

  ! ---------------------------------------------------------
  subroutine av (n, v, w)
    integer :: n
    FLOAT :: v(n), w(n)
    integer :: i
    R_TYPE, allocatable :: psi(:, :), hpsi(:, :)

    ALLOCATE(psi(NP, 1),  NP*1)
    ALLOCATE(hpsi(NP, 1), NP*1)

    do i = 1, NP
      psi(i, 1) = v(i)/sqrt(gr%m%vol_pp(i))
    end do
    call X(hpsi)(h, gr, psi, hpsi, ik)
    do i = 1, NP
      w(i) = hpsi(i, 1)*sqrt(gr%m%vol_pp(i))
    end do

    deallocate(psi, hpsi)
  end subroutine av

end subroutine X(eigen_solver_arpack)
