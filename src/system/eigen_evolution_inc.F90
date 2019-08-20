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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
subroutine X(eigensolver_evolution)(mesh, st, hm, te, tol, niter, converged, ik, diff, tau)
  type(mesh_t),             target, intent(in)    :: mesh
  type(states_elec_t),              intent(inout) :: st
  type(hamiltonian_elec_t), target, intent(inout) :: hm
  type(exponential_t),              intent(inout) :: te
  FLOAT,                            intent(in)    :: tol
  integer,                          intent(inout) :: niter
  integer,                          intent(inout) :: converged
  integer,                          intent(in)    :: ik
  FLOAT,                            intent(out)   :: diff(:) !< (1:st%nst)
  FLOAT,                            intent(in)    :: tau

  integer :: ib, minst, maxst, ist, iter, maxiter, conv, convb, matvec, i
  R_TYPE, allocatable :: c(:, :), zeig(:), res(:)
  FLOAT, allocatable :: eig(:)
  type(batch_t) :: hpsib
#if defined(R_TREAL)
  type(batch_t) :: zpsib
  CMPLX, allocatable :: zpsi(:,:)
  FLOAT, allocatable :: psi(:,:)
#endif

  PUSH_SUB(X(eigensolver_evolution))

  maxiter = niter
  matvec = 0

  SAFE_ALLOCATE(c(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(eig(1:st%nst))
  SAFE_ALLOCATE(zeig(1:st%nst))
  SAFE_ALLOCATE(res(1:st%nst))
#if defined(R_TREAL)
  SAFE_ALLOCATE(zpsi(1:mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:st%d%dim))
#endif

  ! Warning: it seems that the algorithm is improved if some extra states are added -- states
  ! whose convergence should not be checked.
  convb = st%group%block_start - 1
  conv = 0
  do ib = st%group%block_start, st%group%block_end
    if (conv + st%group%psib(ib, ik)%nst <= converged) then
      conv = conv + st%group%psib(ib, ik)%nst
      convb = convb + 1
    else
      exit
    end if
  end do

  do iter = 1, maxiter
#if defined(R_TREAL)
    ! The application of the exponential for the real case is still done one state at a time, as we need
    ! to do a couple of type conversions. To avoid this, we need either a function to convert the
    ! type of a batch, or to modify the batch_copy_data routine to allow the copy between batches of
    ! different types.
    do ist = conv + 1, st%nst
      call batch_init(zpsib, hm%d%dim, 1)
      call states_elec_get_state(st, mesh, ist, ik, zpsi)
      call batch_add_state(zpsib, ist, zpsi)

      call exponential_apply_batch(te, mesh, hm, zpsib, ik, -tau, imag_time = .true.)

      call batch_get_state(zpsib, 1, mesh%np, zpsi)
      psi(1:mesh%np, 1:st%d%dim) = R_TOTYPE(zpsi(1:mesh%np, 1:st%d%dim))
      call states_elec_set_state(st, mesh, ist, ik, psi)
      call batch_end(zpsib)
      matvec = matvec + te%exp_order
    end do
#else
    do ib = convb + 1, st%group%block_end
      call exponential_apply_batch(te, mesh, hm, st%group%psib(ib, ik), ik, -tau, imag_time = .true.)
      matvec = matvec + te%exp_order*(states_elec_block_max(st, ib) - states_elec_block_min(st, ib) + 1)
    end do
#endif

    ! This is the orthonormalization suggested by Aichinger and Krotschek
    ! [Comp. Mat. Science 34, 188 (2005)]
    call X(states_elec_calc_overlap)(st, mesh, ik, c)

    call lalg_eigensolve(st%nst, c, eig)

    do i = 1, st%nst
      c(1:st%nst, i) = c(1:st%nst, i)/sqrt(eig(i))
    end do

    call states_elec_rotate(mesh, st, c, ik)
    
    ! Get the eigenvalues and the residues.
    do ib = convb + 1, st%group%block_end
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      if (hamiltonian_elec_apply_packed(hm, mesh)) call batch_pack(st%group%psib(ib, ik))

      call batch_copy(st%group%psib(ib, ik), hpsib)

      call X(hamiltonian_elec_apply_batch)(hm, mesh, st%group%psib(ib, ik), hpsib, ik)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), hpsib, zeig(minst:maxst))
      call batch_axpy(mesh%np, -zeig(minst:maxst), st%group%psib(ib, ik), hpsib)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), hpsib, res(minst:maxst))
      call batch_end(hpsib)

      do ist = minst, maxst
        diff(ist) = sqrt(abs(res(ist)))
        st%eigenval(ist, ik) = R_REAL(zeig(ist))
      end do

      if (hamiltonian_elec_apply_packed(hm, mesh)) call batch_unpack(st%group%psib(ib, ik), copy=.false.)

      if (debug%info) then
        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') 'Debug: Evolution Eigensolver - ik', ik, &
            ' ist ', ist, ' iter ', iter, ' res ', diff(ist)
          call messages_info(1)
        end do
      end if

    end do

    ! And check for convergence. Note that they must be converged *in order*, so that they can be frozen.
    do ib = convb + 1, st%group%block_end
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)
      if (all(diff(minst:maxst) < tol) .and. ib == convb + 1) then
        convb = convb + 1
        conv = conv + maxst - minst + 1
      end if
    end do
    if (convb == st%group%block_end) exit
  end do

  converged = conv

  niter = matvec

#if defined(R_TREAL)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(zpsi)
#endif
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(eig)
  SAFE_DEALLOCATE_A(zeig)
  SAFE_DEALLOCATE_A(res)

  POP_SUB(X(eigensolver_evolution))
end subroutine X(eigensolver_evolution)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
