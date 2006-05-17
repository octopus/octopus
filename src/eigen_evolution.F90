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

! ---------------------------------------------------------
subroutine X(eigen_solver_evolution) (gr, st, h, tol, niter, converged, diff, tau)
  type(grid_t), target,intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), target, intent(inout)    :: h
  FLOAT,               intent(in)    :: tol
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged
  FLOAT,               intent(out)   :: diff(1:st%nst,1:st%d%nik)
  FLOAT,               intent(in)    :: tau

  integer :: ik, ist, iter, maxiter, conv, conv_, matvec, i, j
  R_TYPE, allocatable :: hpsi(:, :), m(:, :), c(:, :), phi(:, :, :)
  FLOAT, allocatable :: eig(:)

  call push_sub('eigen_evolution.eigen_solver_evolution')

  maxiter = niter
  conv_ = 0
  matvec = 0

  ALLOCATE(hpsi(gr%m%np, st%d%dim), gr%m%np*st%d%dim)
  ALLOCATE(m(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(c(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(eig(st%nst), st%nst)
  ALLOCATE(phi(gr%m%np, st%d%dim, st%nst), gr%m%np*st%d%dim*st%nst)

  ! Warning: it seems that the algorithm is improved if some extra states are added -- states
  ! whose convergence should not be checked.
  kpoints: do ik = 1, st%d%nik
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
          m(i, j) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, i, ik), st%X(psi)(:, :, j, ik) )
        end do
      end do
      call lalg_eigensolve(st%nst, m, c, eig)
      do i = 1, st%nst
        c(:, i) = c(:, i) / sqrt(eig(i))
      end do
      ! Internally the BLAS does the Winograd-Strassen algorithm?
      call lalg_gemm(gr%m%np * st%d%dim, st%nst, st%nst, R_TOTYPE(M_ONE), &
        st%X(psi)(:, :, :, ik), c, R_TOTYPE(M_ZERO), phi)
      do i = 1, st%nst
        st%X(psi)(:, :, i, ik) = phi(:, :, st%nst -i + 1)
      end do

      ! Get the eigenvalues and the residues.
      do ist = conv + 1, st%nst
        call X(hpsi)(h, gr, st%X(psi)(:, :, ist, ik), hpsi, ik)
        st%eigenval(ist, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), hpsi)
        diff(ist, ik) = X(states_residue)(gr%m, st%d%dim, hpsi, st%eigenval(ist, ik), st%X(psi)(:, :, ist, ik))
      end do

      ! Reordering.... (maybe this is unnecessary since the orthonormalization already orders them...)
      if(st%nst > 1) call sort(st%eigenval(1:st%nst, ik), st%X(psi)(:, :, 1:st%nst, ik))

      ! And check for convergence. Note that they must be converged *in order*, so that they can be frozen.
      do ist = conv + 1, st%nst
        if( (diff(ist, ik) < tol) .and. (ist == conv + 1) ) conv = conv + 1
      end do
      if(conv == st%nst) exit
    end do

    conv_ = conv_ + conv
  end do kpoints

  converged = conv_
  niter = matvec
  deallocate(hpsi, m, c, eig, phi)
  call pop_sub()
contains

  ! ---------------------------------------------------------
  subroutine exponentiate(psi, j)
    R_TYPE, intent(inout) :: psi(:, :)
    integer, intent(out) :: j

    integer :: n, m, itrace, iflag, lwsp, liwsp
    integer, allocatable :: iwsp(:)
    FLOAT :: t, tolerance, anorm
    R_TYPE, allocatable :: w(:, :), wsp(:)


    n = st%d%dim * gr%m%np
    m = 20 ! maximum size for the Krylov basis.
    t = -tau
    tolerance = CNST(1.0e-8)
    anorm = M_ONE
    itrace = 0
    lwsp = n*(m+1)+n+(m+2)**2+4*(m+2)**2+6+1
    liwsp = m + 2
    iflag = 0

    ALLOCATE(wsp(lwsp), lwsp)
    ALLOCATE(iwsp(lwsp), lwsp)
    ALLOCATE(w(gr%m%np, st%d%dim), gr%m%np*st%d%dim)

    h_  => h
    gr_ => gr
    ik_ =  ik
    call X(gexpv)(n, m, t, psi(:, 1), w(:, 1), tolerance, anorm, &
      wsp, lwsp, iwsp, liwsp, X(mv), itrace, iflag)
    nullify(h_)
    nullify(gr_)
    psi = w
    deallocate(w)

    j = iwsp(1)
  end subroutine exponentiate

end subroutine X(eigen_solver_evolution)


! ---------------------------------------------------------
subroutine X(mv) (x, y)
  R_TYPE, intent(in) :: x(:)
  R_TYPE, intent(out) :: y(:)
  integer :: idim
  R_TYPE, allocatable :: psi(:, :), hpsi(:, :)

  ALLOCATE(psi(gr_%m%np_part, h_%d%dim), gr_%m%np_part*h_%d%dim)
  ALLOCATE(hpsi(gr_%m%np, h_%d%dim), gr_%m%np*h_%d%dim)

  do idim = 1, h_%d%dim
    psi(1:gr_%m%np, idim) = x((idim-1)*gr_%m%np+1:idim*gr_%m%np)
  end do
  call X(hpsi)(h_, gr_, psi, hpsi, ik_)
  do idim = 1, h_%d%dim
    y((idim-1)*gr_%m%np+1:idim*gr_%m%np) = hpsi(1:gr_%m%np, idim)
  end do

  deallocate(psi, hpsi)

end subroutine X(mv)

