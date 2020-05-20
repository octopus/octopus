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
subroutine X(eigensolver_evolution)(namespace, mesh, st, hm, te, tol, niter, converged, ik, diff, tau)
  type(namespace_t),                intent(in)    :: namespace
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
  FLOAT, allocatable :: eig(:), norm(:)
  type(wfs_elec_t) :: hpsib
#if defined(R_TREAL)
  type(wfs_elec_t) :: zpsib
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
  SAFE_ALLOCATE(norm(1:st%nst))
#if defined(R_TREAL)
  SAFE_ALLOCATE(zpsi(1:mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(psi(1:mesh%np_part, 1:st%d%dim))
#endif

  ! Warning: it seems that the algorithm is improved if some extra states are added -- states
  ! whose convergence should not be checked.
  convb = st%group%block_start - 1
  conv = 0
  !In case of independent particles, the Hamiltonian does not depends upon the density and we can freeze
  !the lowest converged states
  if(hm%theory_level == INDEPENDENT_PARTICLES) then
    do ib = st%group%block_start, st%group%block_end
      if (conv + st%group%psib(ib, ik)%nst <= converged) then
        conv = conv + st%group%psib(ib, ik)%nst
        convb = convb + 1
      else
        exit
      end if
    end do
  end if

  !I am leaving the loop as it is in case we want to change the logic of the code.
  ! However, we would need to compute the density and v_ks inside the loop if we do so.
  do iter = 1, maxiter
#if defined(R_TREAL)
    ! The application of the exponential for the real case is still done one state at a time, as we need
    ! to do a couple of type conversions. To avoid this, we need either a function to convert the
    ! type of a batch, or to modify the batch_copy_data_to routine to allow the copy between batches of
    ! different types.
    do ist = conv + 1, st%nst
      call states_elec_get_state(st, mesh, ist, ik, zpsi)
      call wfs_elec_init(zpsib, hm%d%dim, ist, ist, zpsi, ik)

      call hamiltonian_elec_base_set_phase_corr(hm%hm_base, mesh, zpsib)
      call exponential_apply_batch(te, namespace, mesh, hm, zpsib, -tau, imag_time = .true.)
      call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, mesh, zpsib)

      call batch_get_state(zpsib, 1, mesh%np, zpsi)
      psi(1:mesh%np, 1:st%d%dim) = R_TOTYPE(zpsi(1:mesh%np, 1:st%d%dim))
      norm(ist) = X(mf_nrm2)(mesh, st%d%dim, psi)
      if(norm(ist) > CNST(1.1)) then
        message(1) = "Evolution eigensolver: the time evolution seems to be unstable."
        message(2) = "Please reduce the value of EigensolverImaginaryTime."
        call messages_fatal(2)
      end if
      call states_elec_set_state(st, mesh, ist, ik, psi)
      call zpsib%end()
      matvec = matvec + te%exp_order
    end do
#else
    do ib = convb + 1, st%group%block_end
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      if (hamiltonian_elec_apply_packed(hm)) then
        call st%group%psib(ib, ik)%do_pack()
      end if

      call hamiltonian_elec_base_set_phase_corr(hm%hm_base, mesh, st%group%psib(ib, ik))
      call exponential_apply_batch(te, namespace, mesh, hm, st%group%psib(ib, ik), -tau, imag_time = .true.)
      call hamiltonian_elec_base_unset_phase_corr(hm%hm_base, mesh, st%group%psib(ib, ik))
      matvec = matvec + te%exp_order*(states_elec_block_max(st, ib) - states_elec_block_min(st, ib) + 1)

      call mesh_batch_nrm2(mesh, st%group%psib(ib, ik), norm(minst:maxst))
      if(any(norm(minst:maxst) > CNST(1.1))) then
        message(1) = "Evolution eigensolver: the time evolution seems to be unstable."
        message(2) = "Please reduce the value of EigensolverImaginaryTime."
        call messages_fatal(2)
      end if

      if (hamiltonian_elec_apply_packed(hm)) then
        call st%group%psib(ib, ik)%do_unpack()
      end if
    end do
#endif

    ! This is not the orthonormalization suggested by Aichinger and Krotschek
    ! [Comp. Mat. Science 34, 188 (2005)]
    ! However, the current scheme seem to be more robbust (Cholesky).
    ! Otherwise the norms after orthonormalization are not exactly 1, which create a problem
    call X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)

    ! Get the eigenvalues and the residues.
    do ib = convb + 1, st%group%block_end
      if (hamiltonian_elec_apply_packed(hm)) then
        call st%group%psib(ib, ik)%do_pack()
      end if

      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      call st%group%psib(ib, ik)%copy_to(hpsib)

      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), hpsib)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), hpsib, zeig(minst:maxst))
      st%eigenval(minst:maxst, ik) = R_REAL(zeig(minst:maxst))
      call batch_axpy(mesh%np, -st%eigenval(:,ik), st%group%psib(ib, ik), hpsib)
      call mesh_batch_nrm2(mesh, hpsib, diff(minst:maxst))
      call hpsib%end()

      if (debug%info) then
        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') 'Debug: Evolution Eigensolver - ik', ik, &
            ' ist ', ist, ' iter ', iter, ' res ', diff(ist)
          call messages_info(1)
        end do
      end if

      if (hamiltonian_elec_apply_packed(hm)) then
        call st%group%psib(ib, ik)%do_unpack(copy=.false.)
      end if
    end do

    ! And check for convergence. Note that they must be converged *in order*, so that they can be frozen.
    if(hm%theory_level == INDEPENDENT_PARTICLES) then
      do ib = convb + 1, st%group%block_end
        minst = states_elec_block_min(st, ib)
        maxst = states_elec_block_max(st, ib)
        if (all(diff(minst:maxst) < tol) .and. ib == convb + 1) then
          convb = convb + 1
          conv = conv + maxst - minst + 1
        end if
      end do
      if (convb == st%group%block_end) exit
    end if
  end do

  converged = conv

  niter = matvec

#if defined(R_TREAL)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(zpsi)
#endif
  SAFE_DEALLOCATE_A(c)
  SAFE_DEALLOCATE_A(eig)
  SAFE_DEALLOCATE_A(norm)
  SAFE_DEALLOCATE_A(zeig)
  SAFE_DEALLOCATE_A(res)

  POP_SUB(X(eigensolver_evolution))
end subroutine X(eigensolver_evolution)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
