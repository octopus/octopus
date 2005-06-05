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

subroutine eigen_solver_arpack(m, f_der, st, h, tol, niter, converged, diff)
  type(mesh_type),        intent(in)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(in)    :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)

#if defined(R_TREAL)
  integer :: ido, n, nev, ncv, ldv, iparam(11), ipntr(11), &
             lworkl, info, ist, i
  logical :: rvec
  logical, allocatable :: select(:)
  character(len=1) :: bmat, howmny
  character(len=2) :: which
  FLOAT :: tol_, sigma
  FLOAT, allocatable :: resid(:), v(:, :), workd(:), workl(:), &
                        d(:), psi(:, :), hpsi(:, :)
  integer, parameter :: ncv_factor = 5

  call push_sub('eigen_solver_arpack')

  bmat = 'I'
  n = m%np*st%d%dim
  which = 'SA'
  nev = st%nst
  ncv = ncv_factor*nev
  ldv = n
  lworkl = ncv*(ncv+8)
  tol_ = tol
  ido = 0
  info = 1
  iparam(1) = 1 ! ishfts
  iparam(3) = niter*st%nst ! maxitr
  iparam(4) = 1
  iparam(5) = 0
  iparam(7) = 1 ! mode1
  rvec = .true.
  howmny = 'A'

  allocate(resid(n), v(n, ncv), workl(lworkl), workd(3*n))
  allocate(psi(m%np, st%d%dim), hpsi(m%np, st%d%dim))
  allocate(select(nev), d(nev))

  do i = 1, m%np
     resid(i) = sum(st%dpsi(i, 1, 1:st%nst, 1))
  enddo

  do
    call dsaupd(ido, bmat, n, which, nev, tol_, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, &
                lworkl, info)
    if(info < 0) then
      write(message(1),'(a,i5)') 'Error in ARPACK package routine dsaupd:', info
      call write_fatal(1)
    endif
    if(   .not.(ido .eq. -1 .or. ido .eq. 1)    ) exit
      psi(1:m%np, 1) = workd(ipntr(1):ipntr(1)+n-1)
      call dhpsi(h, m, f_der, psi, hpsi, 1)
      workd(ipntr(2):ipntr(2)+n-1) = hpsi(1:m%np, 1)
  enddo

  niter = iparam(9)
  converged = iparam(5)

  iparam(5) = st%nst ! This is a really smart trick that I figured out: it permits
                     ! to retrieve the approximation to the non-converged eigenvectors.
  call dseupd ( rvec, howmny, select, d, v, ldv, sigma, bmat, n, which, nev, tol, &
                resid, ncv, v, ldv, iparam(1:7), ipntr, workd, workl, lworkl, info )
  if(info .ne. 0) then
     write(message(1),'(a,i5)') 'Error in ARPACK package routine dseupd:', info
     call write_fatal(1)
  endif

  do ist = 1, st%nst
     st%dpsi(1:m%np, 1, ist, 1) = v(1:n, ist)/sqrt(m%vol_pp(1))
     st%eigenval(ist, 1) = d(ist)
     diff(ist, 1) = workl(ipntr(9) + ist - 1)
  enddo

  deallocate(resid, v, workl, workd, psi, hpsi, select, d)
  call pop_sub()
#else
  message(1) = 'Not written yet... coming soon.'
  call write_fatal(1)
#endif
end subroutine eigen_solver_arpack
