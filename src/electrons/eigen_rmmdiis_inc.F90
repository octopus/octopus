!! Copyright (C) 2009 X. Andrade
!! Copyright (C) 2020 N. Tancogne-Dejean
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
!> See http://prola.aps.org/abstract/PRB/v54/i16/p11169_1
subroutine X(eigensolver_rmmdiis) (namespace, mesh, st, hm, pre, tol, niter, converged, ik, diff)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(states_elec_t), target, intent(inout) :: st
  type(hamiltonian_elec_t),    intent(in)    :: hm
  type(preconditioner_t),      intent(in)    :: pre
  FLOAT,                       intent(in)    :: tol
  integer,                     intent(inout) :: niter
  integer,                     intent(inout) :: converged
  integer,                     intent(in)    :: ik
  FLOAT,                       intent(out)   :: diff(:) !< (1:st%nst)

  R_TYPE, allocatable :: mm(:, :, :, :), evec(:, :, :), finalpsi(:)
  R_TYPE, allocatable :: eigen(:), nrmsq(:), eigen_full(:)
  FLOAT,  allocatable :: eval(:, :)
  FLOAT,  allocatable :: lambda(:)
  integer :: ist, minst, idim, ii, iter, nops, maxst, jj, bsize, ib, jter, kter, prog
  FLOAT :: ca, cb, cc
  R_TYPE, allocatable :: fr(:, :), me(:, :)
  type(batch_pointer_t), allocatable :: psib(:), resb(:)
  integer, allocatable :: done(:), last(:)
  logical, allocatable :: failed(:)
  logical :: pack
  integer :: err

  PUSH_SUB(X(eigensolver_rmmdiis))

  pack = hamiltonian_elec_apply_packed(hm)

  SAFE_ALLOCATE(lambda(1:st%nst))
  SAFE_ALLOCATE(psib(1:niter))
  SAFE_ALLOCATE(resb(1:niter))
  SAFE_ALLOCATE(done(1:st%d%block_size))
  SAFE_ALLOCATE(last(1:st%d%block_size))
  SAFE_ALLOCATE(failed(1:st%d%block_size))
  SAFE_ALLOCATE(me(1:2, 1:st%d%block_size))
  SAFE_ALLOCATE(fr(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(nrmsq(1:st%d%block_size))
  SAFE_ALLOCATE(eigen(1:st%d%block_size))

  do iter = 1, niter
    if(iter /= 1) then
      SAFE_ALLOCATE(psib(iter)%batch)
    end if
    SAFE_ALLOCATE(resb(iter)%batch)
  end do

  nops = 0

  call profiling_in(prof, TOSTRING(X(RMMDIIS)))

  failed = .false.
  prog = 0

  do ib = st%group%block_start, st%group%block_end
    minst = states_elec_block_min(st, ib)
    maxst = states_elec_block_max(st, ib)
    bsize = maxst - minst + 1

    psib(1)%batch => st%group%psib(ib, ik)

    if(pack) call psib(1)%batch%do_pack()

    call psib(1)%batch%copy_to(resb(1)%batch)

    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib(1)%batch, resb(1)%batch)
    nops = nops + bsize

    !me(1) = <psi|H|psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(1)%batch, resb(1)%batch, me(1, :), reduce = .false.)
    !me(2) = <psi|psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(1)%batch, psib(1)%batch, me(2, :), reduce = .false.)
    if(mesh%parallel_in_domains) call mesh%allreduce(me)

    !This is the Rayleigh quotient defined as <psi||H|psi>/<psi|psi>
    do ist = minst, maxst
      st%eigenval(ist, ik) = R_REAL(me(1, ist - minst + 1))/R_REAL(me(2, ist - minst + 1))
    end do

    !resb is now the residue (H-\epsilon)|psi>
    call batch_axpy(mesh%np, -st%eigenval(:, ik), psib(1)%batch, resb(1)%batch)

    done = 0

    !We compute the norm of the residue
    call X(mesh_batch_dotp_vector)(mesh, resb(1)%batch, resb(1)%batch, nrmsq)

    !If the norm is smaller than the tolerance, we are done
    do ii = 1, bsize
      if(sqrt(abs(R_REAL(nrmsq(ii)))) < tol) done(ii) = 1
    end do

    if(all(done(1:bsize) /= 0)) then
      if(pack) call st%group%psib(ib, ik)%do_unpack()
      call resb(1)%batch%end()
      cycle
    end if

    call psib(1)%batch%copy_to(psib(2)%batch)

    ! We now find the new direction to perform the minimization.
    ! In order to do so, we first approximate the new direction |\delta\psi> to do the minimization,
    ! based on the value of the residue.
    ! Then we construct the new approximate state |\phi> = |\psi> + \lambda |\delta\psi> = |\psi> + \lambda K|R>.
    ! For this, we need to find \lambda, which we define as the value minimizing the residue associated 
    ! with \phi.

    !We find the new direction from the residual vector, see Kresse and Furthmueller
    ! Computational Materials Science 6 (1996) 15-50, Eq. 50 and discussion above
    !psib(2) is now \delta\psi> = K|R>
    call X(preconditioner_apply_batch)(pre, namespace, mesh, hm, resb(1)%batch, psib(2)%batch, ik)

    call psib(1)%batch%copy_to(resb(2)%batch)

    !resb(2) is now H\delta\psi> 
    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib(2)%batch, resb(2)%batch)
    nops = nops + bsize

    !Here we have the following quantities
    ! fr(1) is <\delta\psi|\delta\psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(2)%batch, psib(2)%batch, fr(1, :), reduce = .false.)
    ! fr(2) is <\psi|\delta\psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(1)%batch, psib(2)%batch, fr(2, :), reduce = .false.)
    ! fr(3) is <\delta\psi|H|\delta\psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(2)%batch, resb(2)%batch, fr(3, :), reduce = .false.)
    ! fr(4) is <\psi|H|\delta\psi>
    call X(mesh_batch_dotp_vector)(mesh, psib(1)%batch, resb(2)%batch, fr(4, :), reduce = .false.)

    if(mesh%parallel_in_domains) call mesh%allreduce(fr)

    !The residue of the new state \phi reads as
    !  R(\phi) = ( <\psi + \lambda\delta\psi | H | \psi + \lambda\delta\psi >)/
    !            (<\psi + \lambda\delta\psi |\psi + \lambda\delta\psi >)
    ! By taking the derivative of it w.r.t. \lambda, we arrive to second-order polynomial that 
    ! gives the value of \lambda that minimize R(\phi)
    ! The equation as the form of A\lambda^2 + B\lambda + C
    ! where A = fr(3) * fr(2) - fr(4) * fr(1) 
    !       B = fr(3) * me(2) - me(1) * fr(1)
    !       C = fr(4) * me(2) - me(1) * fr(2)

    do ist = minst, maxst
      ii = ist - minst + 1

      ca = R_REAL(fr(3, ii))*R_REAL(fr(2, ii)) - R_REAL(fr(4, ii))*R_REAL(fr(1, ii))
      cb = R_REAL(me(2, ii))*R_REAL(fr(3, ii)) - R_REAL(me(1, ii))*R_REAL(fr(1, ii))
      cc = R_REAL(me(2, ii))*R_REAL(fr(4, ii)) - R_REAL(fr(2, ii))*R_REAL(me(1, ii))

      !This is the solution of ca*x^2+cb*x+cc using Mueller formula
      lambda(ist) = - M_TWO * cc / (cb + sqrt(cb**2 - CNST(4.0) * ca * cc))

      ! restrict the value of lambda to be between 0.1 and 1.0
      if(abs(lambda(ist)) > CNST(1.0)) lambda(ist) = lambda(ist)/abs(lambda(ist))
      if(abs(lambda(ist)) < CNST(0.1)) lambda(ist) = CNST(0.1)*lambda(ist)/abs(lambda(ist))
    end do

    SAFE_ALLOCATE(mm(1:niter, 1:niter, 1:2, 1:bsize))

    do iter = 2, niter

      ! for iter == 2 the preconditioning was done already
      if(iter > 2) then
        call psib(iter - 1)%batch%copy_to(psib(iter)%batch)
        call X(preconditioner_apply_batch)(pre, namespace, mesh, hm, resb(iter - 1)%batch, psib(iter)%batch, ik)
      end if

      ! predict by jacobi
      ! Here we compute the new trial wavefunction using a fixed \lambda found earlier
      call batch_xpay(mesh%np, psib(iter - 1)%batch, lambda, psib(iter)%batch)

      if(iter > 2) then
         call psib(iter)%batch%copy_to(resb(iter)%batch)
      end if

      ! calculate the residual
      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib(iter)%batch, resb(iter)%batch)
      nops = nops + bsize

      ! According to Kresse and Furthmueller, Computational Materials Science 6 (1996) 15-50 
      ! we update st%eigenval here.
      if(iter == 2) then
        call X(mesh_batch_dotp_vector)(mesh, psib(iter)%batch, resb(iter)%batch, me(1, 1:bsize), reduce = .false.)
        call X(mesh_batch_dotp_vector)(mesh, psib(iter)%batch, psib(iter)%batch, me(2, 1:bsize), reduce = .false.)
        if(mesh%parallel_in_domains) call mesh%allreduce(me)
        do ist = minst, maxst
          st%eigenval(ist, ik) = R_REAL(me(1, ist - minst + 1))/R_REAL(me(2, ist - minst + 1))
        end do
      end if

      !We compute the new residue
      call batch_axpy(mesh%np, -st%eigenval(:, ik), psib(iter)%batch, resb(iter)%batch)

      !We now perform the direct inversion in the iterative subspace (DIIS)
      !which is solved as an Hermitian eigenvalue problem
      call profiling_in(prof_iter, TOSTRING(X(RMMDIIS_MATRIX)))
      ! calculate the matrix elements between iterations
      do jter = 1, iter
        do kter = 1, jter
          
          if(jter < iter - 1 .and. kter < iter - 1) then
            ! it was calculated on the previous iteration
            ! in parallel this was already reduced, so we set it to zero in non-root ranks
            if(mesh%parallel_in_domains .and. mesh%mpi_grp%rank /= 0) mm(jter, kter, 1:2, 1:bsize) = CNST(0.0)
            cycle
          end if

          call X(mesh_batch_dotp_vector)(mesh, resb(jter)%batch, resb(kter)%batch, mm(jter, kter, 1, :), reduce = .false.)
          call X(mesh_batch_dotp_vector)(mesh, psib(jter)%batch, psib(kter)%batch, mm(jter, kter, 2, :), reduce = .false.)
        end do
      end do
      call profiling_out(prof_iter)
            
      ! symmetrize
      do jter = 1, iter
        do kter = jter + 1, iter
          mm(jter, kter, 1:2, 1:bsize) = R_CONJ(mm(kter, jter, 1:2, 1:bsize))
        end do
      end do

      if(mesh%parallel_in_domains) call mesh%allreduce(mm)

      SAFE_ALLOCATE(evec(1:iter, 1:1, 1:bsize))
      SAFE_ALLOCATE(eval(1:iter, 1:bsize))

      do ist = minst, maxst
        ii = ist - minst + 1

        failed(ii) = .false.
        call lalg_lowest_geneigensolve(1, iter, mm(:, :, 1, ii), mm(:, :, 2, ii), eval(:, ii),  &
                   evec(:, :, ii), preserve_mat=.true., bof = failed(ii), err_code = err)
        if( err < 0 .or. err > iter ) then

          !Due to numerical errors, lapack sometimes finds the overlap matrix to not be
          !definite positive. Here we impose it to be only ones, as this represent 
          !what the numbers are, up to numerical precision.
          mm(1:iter, 1:iter, 2, 1:bsize) = M_ONE
          do jter = 1, iter
            mm(jter,jter,2,1:bsize) = M_ONE + M_EPSILON
          end do

          failed(ii) = .false.
          call lalg_lowest_geneigensolve(1, iter, mm(:, :, 1, ii), mm(:, :, 2, ii), eval(:, ii),  &
                   evec(:, :, ii), preserve_mat=.true., bof = failed(ii), err_code = err)

          !We should never fall into this case, as we impose B to be definite positive here
          if(err < 0 .or. err > iter) then
            failed(ii) = .true.
            last(ii) = iter - 1
     
            evec(1:iter - 1, 1, ii) = CNST(0.0)
            evec(iter, 1, ii) = CNST(1.0)
            cycle
          else
            failed(ii) = .false.
          end if
        else !In this case we did not reach the tolerance
          failed(ii) = .false.
        end if
      end do

      call resb(iter)%batch%end()

      call profiling_in(prof_lc, TOSTRING(X(RMMDIIS_LC)))

      !Here we construct the linear combination of the states, see Eq. 68 
      call batch_scal(mesh%np, evec(iter, 1, :), psib(iter)%batch, a_start = minst)

      do jj = 1, iter - 1
        if(pack) call psib(jj)%batch%do_pack()
        call batch_axpy(mesh%np, evec(jj, 1, :), psib(jj)%batch, psib(iter)%batch, a_start = minst)
        if(pack) call psib(jj)%batch%do_unpack(copy = .false.)
      end do

      call profiling_out(prof_lc)

      call psib(iter)%batch%copy_to(resb(iter)%batch)

      ! re-calculate the residual
      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib(iter)%batch, resb(iter)%batch)
      nops = nops + bsize
      call batch_axpy(mesh%np, -st%eigenval(:, ik), psib(iter)%batch, resb(iter)%batch)

      ! why not allocate these outside the loop?
      SAFE_DEALLOCATE_A(eval)      
      SAFE_DEALLOCATE_A(evec)

      if(debug%info) then
        call X(mesh_batch_dotp_vector)(mesh, resb(iter)%batch, resb(iter)%batch, eigen)

        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') &
            'Debug: RMMDIIS Eigensolver - ik', ik, ' ist ', ist, ' iter ', iter, ' res ', sqrt(abs(eigen(ist - minst + 1)))
          call messages_info(1)
        end do
      end if

    end do ! iter

    SAFE_DEALLOCATE_A(mm)

    ! end with a trial move
    call X(preconditioner_apply_batch)(pre, namespace, mesh, hm, resb(niter)%batch, resb(niter - 1)%batch, ik)

    call batch_xpay(mesh%np, psib(niter)%batch, lambda, resb(niter - 1)%batch)

    if(any(failed(1:bsize))) then 
      SAFE_ALLOCATE(finalpsi(1:mesh%np))

      do ist = minst, maxst
        ii = ist - minst + 1
        
        if(failed(ii)) then
          do idim = 1, st%d%dim
            call batch_get_state(psib(last(ii))%batch, (/ist, idim/), mesh%np, finalpsi)
            call batch_set_state(resb(niter - 1)%batch, (/ist, idim/), mesh%np, finalpsi)
          end do
        end if
      end do

      SAFE_DEALLOCATE_A(finalpsi)
    end if

    ! we can remove most of the batches
    do iter = 1, niter
      if(iter /= 1) call psib(iter)%batch%end()
      if(iter /= niter - 1) call resb(iter)%batch%end()
    end do

    call resb(niter - 1)%batch%copy_data_to(mesh%np, st%group%psib(ib, ik))

    call resb(niter - 1)%batch%end()

    if(pack) call st%group%psib(ib, ik)%do_unpack()

    prog = prog + bsize
    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) + prog, st%lnst*st%d%kpt%nlocal)
    end if
    
  end do ! ib

  call profiling_out(prof)

  call X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)

  ! recalculate the eigenvalues and residuals
  SAFE_ALLOCATE(eigen_full(1:st%nst))
  eigen_full(1:st%nst) = R_TOTYPE(M_ZERO)
  st%eigenval(1:st%nst, ik) = R_TOTYPE(M_ZERO)

  do ib = st%group%block_start, st%group%block_end
    minst = states_elec_block_min(st, ib)
    maxst = states_elec_block_max(st, ib)

    if(pack) call st%group%psib(ib, ik)%do_pack()
  
    call st%group%psib(ib, ik)%copy_to(resb(1)%batch)
    
    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), resb(1)%batch)
    call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), resb(1)%batch, me(1, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), st%group%psib(ib, ik), me(2, :), reduce = .false.)
    if(mesh%parallel_in_domains) call mesh%allreduce(me)

    !This is the Rayleigh quotient
    do ist = minst, maxst
      st%eigenval(ist, ik) = R_REAL(me(1, ist - minst + 1))/R_REAL(me(2, ist - minst + 1))
    end do

    call batch_axpy(mesh%np, -st%eigenval(:, ik), st%group%psib(ib, ik), resb(1)%batch)

    call X(mesh_batch_dotp_vector)(mesh, resb(1)%batch, resb(1)%batch, eigen_full(minst:maxst), reduce = .false.)

    call resb(1)%batch%end()

    if(pack) call st%group%psib(ib, ik)%do_unpack()
    
    nops = nops + maxst - minst + 1
  end do

  if(mesh%parallel_in_domains) call mesh%allreduce(eigen_full)
  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp, eigen_full)
    call comm_allreduce(st%mpi_grp, st%eigenval(1:st%nst, ik))
  end if

  diff(:) = sqrt(abs(eigen_full(:)))
  SAFE_DEALLOCATE_A(eigen_full)
  
  converged = converged + count(diff(:) <= tol)

  do iter = 1, niter
    if(iter /= 1) then
      SAFE_DEALLOCATE_P(psib(iter)%batch)
    end if
    SAFE_DEALLOCATE_P(resb(iter)%batch)
  end do

  niter = nops

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(psib)
  SAFE_DEALLOCATE_A(resb)
  SAFE_DEALLOCATE_A(done)
  SAFE_DEALLOCATE_A(last)
  SAFE_DEALLOCATE_A(failed)
  SAFE_DEALLOCATE_A(fr)
  SAFE_DEALLOCATE_A(nrmsq)
  SAFE_DEALLOCATE_A(eigen)
  SAFE_DEALLOCATE_A(me)

  POP_SUB(X(eigensolver_rmmdiis))

end subroutine X(eigensolver_rmmdiis)

! ---------------------------------------------------------

subroutine X(eigensolver_rmmdiis_min) (namespace, mesh, st, hm, pre, niter, converged, ik)
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(preconditioner_t),   intent(in)    :: pre
  integer,                  intent(inout) :: niter
  integer,                  intent(inout) :: converged
  integer,                  intent(in)    :: ik

  integer :: sd_steps
  integer :: isd, ist, minst, maxst, ib, ii
  FLOAT  :: ca, cb, cc
  FLOAT, allocatable :: lambda(:)
  R_TYPE, allocatable :: diff(:)
  R_TYPE, allocatable :: me1(:, :), me2(:, :)
  logical :: pack
  type(wfs_elec_t) :: resb, kresb

  PUSH_SUB(X(eigensolver_rmmdiis_min))

  sd_steps = niter
  
  pack = hamiltonian_elec_apply_packed(hm)

  SAFE_ALLOCATE(me1(1:2, 1:st%d%block_size))
  SAFE_ALLOCATE(me2(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(lambda(1:st%nst))

  niter = 0

  if(debug%info) then
    SAFE_ALLOCATE(diff(1:st%d%block_size))
  end if

  do ib = st%group%block_start, st%group%block_end
    minst = states_elec_block_min(st, ib)
    maxst = states_elec_block_max(st, ib)

    if(pack) call st%group%psib(ib, ik)%do_pack()

    call st%group%psib(ib, ik)%copy_to(resb)
    call st%group%psib(ib, ik)%copy_to(kresb)

    do isd = 1, sd_steps

      !We start by computing the Rayleigh quotient
      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), resb)

      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), resb, me1(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik), st%group%psib(ib, ik), me1(2, :), reduce = .false.)

      if(mesh%parallel_in_domains) call mesh%allreduce(me1)

      !This is the Rayleigh quotient
      do ist = minst, maxst
        st%eigenval(ist, ik) = R_REAL(me1(1, ist - minst + 1)/me1(2, ist - minst + 1))
      end do
 
      !We get the residual vector
      call batch_axpy(mesh%np, -st%eigenval(:, ik), st%group%psib(ib, ik), resb)

      if(debug%info) then
        call X(mesh_batch_dotp_vector)(mesh, resb, resb, diff)

        do ist = minst, maxst
          write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') &
            'Debug: RMMDIIS MIN Eigensolver - ik', ik, ' ist ', ist, ' iter ', isd, ' res ', sqrt(abs(diff(ist - minst + 1)))
          call messages_info(1)
        end do
      end if

      call X(preconditioner_apply_batch)(pre, namespace, mesh, hm, resb, kresb, ik)

      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, kresb, resb)

      niter = niter + 2*(maxst - minst + 1)

      call X(mesh_batch_dotp_vector)(mesh, kresb, kresb, me2(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik),  kresb, me2(2, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(mesh, kresb, resb,  me2(3, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(mesh, st%group%psib(ib, ik),  resb,  me2(4, :), reduce = .false.)

      if(mesh%parallel_in_domains) call mesh%allreduce(me2)

      do ist = minst, maxst
        ii = ist - minst + 1

        ca = R_REAL(me2(1, ii))*R_REAL(me2(4, ii)) - R_REAL(me2(3, ii))*R_REAL(me2(2, ii))
        cb = R_REAL(me1(2, ii))*R_REAL(me2(3, ii)) - R_REAL(me1(1, ii))*R_REAL(me2(1, ii))
        cc = R_REAL(me1(1, ii))*R_REAL(me2(2, ii)) - R_REAL(me2(4, ii))*R_REAL(me1(2, ii))

        !This is - the solution of ca*x^2+cb*x+cc
        lambda(ist) = CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))
      end do

      call batch_axpy(mesh%np, lambda, kresb, st%group%psib(ib, ik))
      
    end do

    if(pack) call st%group%psib(ib, ik)%do_unpack()

    call resb%end()
    call kresb%end()

    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(st%lnst*(ik - 1) +  maxst, st%lnst*st%d%kpt%nlocal)
    end if

  end do

  if(debug%info) then
    SAFE_DEALLOCATE_A(diff)
  end if

  call X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)

  converged = 0

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(me1)
  SAFE_DEALLOCATE_A(me2)

  POP_SUB(X(eigensolver_rmmdiis_min))

end subroutine X(eigensolver_rmmdiis_min)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
