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
!! $Id: eigen_evolution_inc.F90 6257 2009-12-26 01:13:24Z xavier $

! ---------------------------------------------------------
subroutine X(eigensolver_evolution) (gr, st, hm, tol, niter, converged, ik, diff, tau)
  type(grid_t),        target, intent(in)    :: gr
  type(states_t),              intent(inout) :: st
  type(hamiltonian_t), target, intent(in)    :: hm
  FLOAT,                       intent(in)    :: tol
  integer,                     intent(inout) :: niter
  integer,                     intent(inout) :: converged
  integer,                     intent(in)    :: ik
  FLOAT,                       intent(out)   :: diff(:) !< (1:st%nst)
  FLOAT,                       intent(in)    :: tau

  integer :: ist, iter, maxiter, conv, matvec, i, j
  R_TYPE, allocatable :: hpsi(:, :), m(:, :), c(:, :), phi(:, :, :)
  FLOAT, allocatable :: eig(:)
  type(exponential_t) :: te

  PUSH_SUB(X(eigensolver_evolution))

  maxiter = niter
  matvec = 0

  call exponential_init(te)

  SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(m(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(c(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(eig(1:st%nst))
  SAFE_ALLOCATE(phi(1:gr%mesh%np_part, 1:st%d%dim, 1:st%nst))

  ! Warning: it seems that the algorithm is improved if some extra states are added -- states
  ! whose convergence should not be checked.
  conv = converged

  do iter = 1, maxiter

    do ist = conv + 1, st%nst
      call exponentiate(st%X(psi)(:, :, ist, ik), j)
      matvec = matvec + j
    end do

    ! This is the orthonormalization suggested by Aichinger and Krotschek
    ! [Comp. Mat. Science 34, 188 (2005)]
    do i = 1, st%nst
      do j = i, st%nst
        m(i, j) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, i, ik), st%X(psi)(:, :, j, ik) )
      end do
    end do
    c = m
    call lalg_eigensolve(st%nst, c, eig)
    do i = 1, st%nst
      c(:, i) = c(:, i) / sqrt(eig(i))
    end do
    
    call lalg_gemm(gr%mesh%np_part * st%d%dim, st%nst, st%nst, R_TOTYPE(M_ONE), &
         st%X(psi)(:, :, :, ik), c, R_TOTYPE(M_ZERO), phi)
    do i = 1, st%nst
      st%X(psi)(1:gr%mesh%np, :, i, ik) = phi(1:gr%mesh%np, :, st%nst -i + 1)
    end do

    ! Get the eigenvalues and the residues.
    do ist = conv + 1, st%nst
      call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist, ik), hpsi, ist, ik)
      st%eigenval(ist, ik) = real(X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), hpsi), REAL_PRECISION)
      diff(ist) = X(states_residue)(gr%mesh, st%d%dim, hpsi, st%eigenval(ist, ik), st%X(psi)(:, :, ist, ik))

      if(in_debug_mode) then
        write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') 'Debug: Evolution Eigensolver - ik', ik, &
          ' ist ', ist, ' iter ', iter, ' res ', diff(ist)
        call messages_info(1)
      end if
    end do

    ! Reordering.... (maybe this is unnecessary since the orthonormalization already orders them...)
    if(st%nst > 1) call sort(st%eigenval(1:st%nst, ik), st%X(psi)(:, :, 1:st%nst, ik))

    ! And check for convergence. Note that they must be converged *in order*, so that they can be frozen.
    do ist = conv + 1, st%nst
      if( (diff(ist) < tol) .and. (ist == conv + 1) ) conv = conv + 1
    end do
    if(conv == st%nst) exit
  end do

  converged = conv

  niter = matvec
  call exponential_end(te)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(m)
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(eig)
  SAFE_DEALLOCATE_A(phi)

  POP_SUB(X(eigensolver_evolution))
contains

  ! ---------------------------------------------------------
  subroutine exponentiate(psi, order)
    R_TYPE, target, intent(inout) :: psi(:, :)
    integer,        intent(out)   :: order
    CMPLX,          pointer       :: zpsi(:, :)

    PUSH_SUB(X(eigensolver_evolution).exponentiate)

#if defined(R_TREAL)
    SAFE_ALLOCATE(zpsi(1:gr%mesh%np_part, 1:st%d%dim))
    zpsi(1:gr%mesh%np, 1:st%d%dim) = psi(1:gr%mesh%np, 1:st%d%dim)
#else
    zpsi => psi
#endif
    call exponential_apply(te, gr%der, hm, zpsi, ist, ik, -tau, M_ZERO, order = order, imag_time = .true.)
#if defined(R_TREAL)
    psi(1:gr%mesh%np, 1:st%d%dim) = R_TOTYPE(zpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_DEALLOCATE_P(zpsi)
#else
    nullify(zpsi)
#endif

    POP_SUB(X(eigensolver_evolution).exponentiate)
  end subroutine exponentiate

end subroutine X(eigensolver_evolution)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
