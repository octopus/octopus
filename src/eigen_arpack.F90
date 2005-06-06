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

subroutine eigen_solver_arpack(m, f_der, st, h, tol, niter, converged, errorflag, diff, reorder, verbose)
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(IN)    :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(out)   :: errorflag
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,      optional, intent(in)    :: reorder
  logical,      optional, intent(in)    :: verbose

  integer :: ido, n, nev, ncv, ldv, iparam(11), lworkl, info, ist, i, ik, idim
  integer, allocatable :: ipntr(:)
  logical :: rvec, symmetric
  logical, allocatable :: select(:)
  character(len=1) :: bmat, howmny
  character(len=2) :: which
  FLOAT :: tol_, sigma, sigmar, sigmai
  R_TYPE, allocatable :: resid(:), v(:, :), workd(:), workl(:), d(:), &
                         di(:), dr(:), workev(:), psi(:, :), hpsi(:, :)
  FLOAT, allocatable :: z(:, :)
  FLOAT, allocatable :: rwork(:)
  integer, parameter :: ncv_factor = 5

  call push_sub('eigen_solver_arpack')

  symmetric = .not.(m%use_curvlinear)


  kpoints: do ik = 1, st%d%nik


  bmat = 'I'
  n = m%np*st%d%dim
#if defined(R_TREAL)
  if(symmetric) then
    which = 'SA'
  else
    which = 'SR'
  endif
#else
  which = 'SR'
#endif

  nev = st%nst
  ncv = ncv_factor*nev
  ldv = n
#if defined(R_TREAL)
  if(symmetric) then
    lworkl = ncv*(ncv+8)
  else
    lworkl =  ncv*(3*ncv + 6)
  endif
#else
  lworkl = ncv*(3*ncv + 5)
#endif
  tol_ = tol
  ido = 0
  info = 1
  iparam(1) = 1 ! ishfts
  iparam(3) = niter*nev ! maxitr

  iparam(4) = 1
  iparam(5) = 0
  iparam(7) = 1 ! mode1
  rvec = .true.
  howmny = 'A'

  allocate(resid(n), v(n, ncv), workl(lworkl), workd(3*n), psi(m%np, st%d%dim), hpsi(m%np, st%d%dim))
#if defined(R_TREAL)
  if(symmetric) then
    allocate(select(nev), d(nev), ipntr(11))
  else
    allocate(select(ncv), d(nev), dr(nev+1), z(n, nev + 1), di(nev+1), workev(3*ncv), ipntr(14))
  endif
#else
  allocate(rwork(ncv), workev(2*ncv), select(ncv), d(nev+1), ipntr(14))
#endif
  select = .true.

  do idim = 1, st%d%dim
     do i = 1, m%np
        resid((idim-1)*m%np+i) = sum(st%X(psi)(i, idim, 1:st%nst, ik))
     enddo
  enddo

  do
    #if defined(R_TREAL)
    if(symmetric) then
      call dsaupd(ido, bmat, n, which, nev, tol_, resid, &
                  ncv, v, ldv, iparam, ipntr, workd, workl, &
                  lworkl, info)
    else
      call dnaupd(ido, bmat, n, which, nev, tol_, resid, &
                  ncv, v, ldv, iparam(1:11), ipntr(1:14), workd, workl, &
                  lworkl, info )
    endif
    #else
    call znaupd(ido, bmat, n, which, nev, tol_, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, &
                lworkl, rwork, info )
    #endif
    if(info < 0) then
      write(message(1),'(a,i5)') 'Error in ARPACK package routine X([s/n]aupd):', info
      call write_fatal(1)
    endif
    if(   .not.(ido .eq. -1 .or. ido .eq. 1)    ) exit
      call blas_copy(m%np*st%d%dim, workd(ipntr(1)), 1, psi(1, 1), 1)
      call X(hpsi)(h, m, f_der, psi, hpsi, ik)
      call blas_copy(m%np*st%d%dim, hpsi(1, 1), 1, workd(ipntr(2)), 1)
  enddo

  niter = iparam(9)
  converged = iparam(5)

  if(symmetric) then
    iparam(5) = st%nst ! This is a really smart trick that I figured out: it permits
                       ! to retrieve the approximation to the non-converged eigenvectors.
                       ! However, it does not work for the non-symmetric cases.
  endif

#if defined(R_TREAL)
  if(symmetric) then
    call dseupd ( rvec, howmny, select, d, v, ldv, sigma, bmat, n, which, nev, tol_, &
                  resid, ncv, v, ldv, iparam(1:7), ipntr, workd, workl, lworkl, info )
  else
    call dneupd ( rvec, howmny, select, dr, di, z, ldv, sigmar, &
                  sigmai, workev, bmat, n, which, nev, tol_, &
                  resid, ncv, v, ldv, iparam, ipntr, workd, &
                  workl, lworkl, info)
  endif
#else
  call zneupd ( rvec, howmny, select, d, v(1:ldv, 1:nev), ldv, sigma, workev, bmat, &
                n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, &
                workl, lworkl, rwork, info )
#endif
  if(info .ne. 0) then
     write(message(1),'(a,i5)') 'Error in ARPACK package routine X([s/n]eupd):', info
     call write_fatal(1)
  endif

#if defined(R_TREAL)
  if(.not.symmetric) then
    do ist = 1, st%nst
       v(1:n, ist) = z(1:n, ist)
       d(ist) = dr(ist)
    enddo
  endif
#endif

  do ist = 1, st%nst
     call blas_copy(n, v(1, ist), 1, st%X(psi)(1, 1, ist, ik), 1)
     st%X(psi)(:, :, ist, ik) = st%X(psi)(:, :, ist, ik)/sqrt(m%vol_pp(1))
     st%eigenval(ist, ik) = real(d(ist), PRECISION)
     #if defined(R_TREAL)
     if(symmetric) then
       diff(ist, ik) = workl(ipntr(9) + ist - 1)
     else
       diff(ist, ik) = workl(ipntr(11) + ist - 1)
     endif
     #else
     diff(ist, ik) = abs(workl(ipntr(11) + ist - 1))
     #endif
  enddo

#if defined(R_TREAL)
  deallocate(resid, v, workl, workd, psi, hpsi, select, ipntr)
  if(symmetric) deallocate(d)
  if(.not.symmetric) deallocate(z, d, dr, di)
#else
  deallocate(resid, v, workl, workd, rwork, psi, hpsi, select, d, workev, ipntr)
#endif

  enddo kpoints

  call pop_sub()
end subroutine eigen_solver_arpack
