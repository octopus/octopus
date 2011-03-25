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
!! $Id: subspace_inc.F90 6258 2009-12-26 20:07:24Z xavier $

! ---------------------------------------------------------
! This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag)(this, der, st, hm, ik, eigenval, psi, diff)
  type(subspace_t),    intent(in)    :: this
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(:, :, :)
  FLOAT, optional,     intent(out)   :: diff(:)

  R_TYPE, allocatable :: h_subspace(:, :), f(:, :, :)
  integer             :: ist, ist2, size, idim
  FLOAT               :: nrm2
  type(profile_t),     save    :: diagon_prof
  type(batch_t) :: psib, hpsib, whole_psib

  PUSH_SUB(X(subspace_diag))
  call profiling_in(diagon_prof, "SUBSPACE_DIAG")

  select case(this%method)
  case(SD_OLD)
    call X(subspace_diag_old)(der, st, hm, ik, eigenval, psi, diff)
  case(SD_SCALAPACK)
    call X(subspace_diag_scalapack)(der, st, hm, ik, eigenval, psi, diff)
  case(SD_STANDARD)
    ASSERT(.not. st%parallel_in_states)

    SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(f(1:der%mesh%np, 1:st%d%dim, 1:st%d%block_size))

    ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
    do ist = st%st_start, st%st_end, st%d%block_size
      size = min(st%d%block_size, st%st_end - ist + 1)

      call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))

      call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)

      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)

      call batch_end(psib)

      call batch_init(whole_psib, hm%d%dim, ist, st%nst, psi(:, :, ist:st%nst))
      call X(mesh_batch_dotp_matrix)(der%mesh, hpsib, whole_psib, h_subspace)
      call batch_end(whole_psib)
      call batch_end(hpsib)
      
    end do

    ! only half of h_subspace has the matrix, but this is what Lapack needs

    ! Diagonalize the Hamiltonian in the subspace.
    call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

    ! Calculate the new eigenfunctions as a linear combination of the
    ! old ones.
    call batch_init(whole_psib, hm%d%dim, 1, st%nst, psi(:, :, :))
    call X(mesh_batch_rotate)(der%mesh, whole_psib, h_subspace)
    call batch_end(whole_psib)

    ! Renormalize.
    do ist = st%st_start, st%st_end
      nrm2 = X(mf_nrm2)(der%mesh, st%d%dim, psi(:, :, ist))
      do idim = 1, st%d%dim
        call lalg_scal(der%mesh%np, M_ONE/nrm2, psi(:, idim, ist))
      end do
    end do

    ! Recalculate the residues if requested by the diff argument.
    if(present(diff)) then 
      
      do ist = st%st_start, st%st_end, st%d%block_size
        size = min(st%d%block_size, st%st_end - ist + 1)

        call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
        call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)
        
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)

        call batch_end(psib)
        call batch_end(hpsib)
        
        do ist2 = ist, ist + size - 1
          diff(ist2) = X(states_residue)(der%mesh, st%d%dim, f(:, :, ist2 - ist + 1), eigenval(ist2), psi(:, :, ist2))
        end do
      end do

    end if

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(h_subspace)

  case default
    ASSERT(.false.)
  end select

  call profiling_out(diagon_prof)
  POP_SUB(X(subspace_diag))

end subroutine X(subspace_diag)

! --------------------------------------------------------- 
! This routine diagonalises the Hamiltonian in the subspace defined by
! the states; this version is aware of parallelization in states but
! consumes more memory.
!
! I leave this function here for the moment. Eventually it will be
! removed. XA
!
subroutine X(subspace_diag_old)(der, st, hm, ik, eigenval, psi, diff)
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(:, :, st%st_start:)
  FLOAT, optional,     intent(out)   :: diff(:)

  R_TYPE, allocatable :: h_subspace(:,:), ff(:,:,:)
  integer             :: ist
  FLOAT               :: nrm2
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT               :: ldiff(st%lnst)
#endif

  PUSH_SUB(X(subspace_diag_old))

  SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(ff(1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))
  
  ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
  do ist = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, der, psi(:, :, ist), ff(:, :, ist), ist, ik)
  end do
  call states_blockt_mul(der%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       psi(:, :, :), ff, h_subspace, symm=.true.)

  ! Diagonalize the Hamiltonian in the subspace.
  call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

  ! The new states are given by the eigenvectors of the matrix.
  ff(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end) = psi(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end)

  call states_block_matr_mul(der%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       ff, h_subspace, psi(:, :, :))

  ! Renormalize.
  do ist = st%st_start, st%st_end
    nrm2 = X(mf_nrm2)(der%mesh, st%d%dim, psi(:, :, ist))
    psi(1:der%mesh%np, 1:st%d%dim, ist) = psi(1:der%mesh%np, 1:st%d%dim, ist)/nrm2
  end do

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do ist = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, der, psi(:, :, ist) , ff(:, :, st%st_start), ist, ik)
      diff(ist) = X(states_residue)(der%mesh, st%d%dim, ff(:, :, st%st_start), eigenval(ist), &
           psi(:, :, ist))
    end do

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      ldiff = diff(st%st_start:st%st_end)
      call lmpi_gen_allgatherv(st%lnst, ldiff, tmp, diff(:), st%mpi_grp)
    end if
#endif
  end if

  SAFE_DEALLOCATE_A(ff)
  SAFE_DEALLOCATE_A(h_subspace)
  
  POP_SUB(X(subspace_diag_old))
  
end subroutine X(subspace_diag_old)

! --------------------------------------------------------- 
! This routine diagonalises the Hamiltonian in the subspace defined by
! the states; this version is aware of parallelization in states but
! consumes more memory.
!
subroutine X(subspace_diag_scalapack)(der, st, hm, ik, eigenval, psi, diff)
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(:, :, st%st_start:)
  FLOAT, optional,     intent(out)   :: diff(:)
 
#ifdef HAVE_SCALAPACK
  R_TYPE, allocatable  :: hs(:, :), hpsi(:, :, :), evectors(:, :), work(:)
  R_TYPE               :: rttmp
  integer              :: tmp, ist, lwork, size
  FLOAT                :: ldiff(st%lnst)
  integer :: psi_block(1:2), total_np, psi_desc(BLACS_DLEN), hs_desc(BLACS_DLEN), info
  integer :: nbl, nrow, ncol
  type(batch_t) :: psib, hpsib
#ifdef R_TCOMPLEX
  integer :: lrwork
  CMPLX, allocatable :: rwork(:)
  CMPLX :: ftmp
#endif

  PUSH_SUB(X(subspace_diag_scalapack))

  SAFE_ALLOCATE(hpsi(1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))
  
  call states_blacs_blocksize(st, der%mesh, psi_block, total_np)

  call descinit(psi_desc(1), total_np, st%nst, psi_block(1), psi_block(2), 0, 0,  st%dom_st_proc_grid%context, &
    st%d%dim*der%mesh%np_part, info)

  ! select the blocksize, we use the division used for state
  ! parallelization but with a maximum of 64
  nbl = min(64, psi_block(2))

  ! calculate the size of the matrix in each node
  nrow = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
  ncol = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

  SAFE_ALLOCATE(hs(1:nrow, 1:ncol))

  call descinit(hs_desc(1), st%nst, st%nst, nbl, nbl, 0, 0, st%dom_st_proc_grid%context, nrow, info)

  ! calculate |hpsi> = H |psi>
  do ist = st%st_start, st%st_end, st%d%block_size
    size = min(st%d%block_size, st%st_end - ist + 1)
    
    call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
    call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, hpsi(: , :, ist:))
    
    call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)
    
    call batch_end(psib)
    call batch_end(hpsib)
  end do

  ! We need to set to zero some extra parts of the array
  if(st%d%dim == 1) then
    psi(der%mesh%np + 1:psi_block(1), 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    hpsi(der%mesh%np + 1:psi_block(1), 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
  else
    psi(der%mesh%np + 1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    hpsi(der%mesh%np + 1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
  end if
  
  ! get the matrix <psi|H|psi> = <psi|hpsi>
  call pblas_gemm('c', 'n', st%nst, st%nst, total_np, &
    R_TOTYPE(der%mesh%vol_pp(1)), psi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    hpsi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    R_TOTYPE(M_ZERO), hs(1, 1), 1, 1, hs_desc(1))

  SAFE_ALLOCATE(evectors(1:nrow, 1:ncol))

  ! now diagonalize
#ifdef R_TCOMPLEX

  call pzheev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), &
    work = rttmp, lwork = -1, rwork = ftmp, lrwork = -1, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pzheev workspace query failure, error code = ", info
    call messages_fatal(1)
  endif

  lwork = nint(abs(rttmp))
  lrwork = nint(real(ftmp, 8))

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:lrwork))

  ! parameter 4 is bad, they say
  call pzheev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), &
    work = work(1), lwork = lwork, rwork = rwork(1), lrwork = lrwork, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pzheev call failure, error code = ", info
    call messages_fatal(1)
  endif

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)

#else

  call pdsyev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), work = rttmp, lwork = -1, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pdsyev workspace query failure, error code = ", info
    call messages_fatal(1)
  endif

  lwork = nint(abs(rttmp))
  SAFE_ALLOCATE(work(1:lwork))

  call pdsyev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), work = work(1), lwork = lwork, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pdsyev call failure, error code = ", info
    call messages_fatal(1)
  endif
  
  SAFE_DEALLOCATE_A(work)

#endif

  SAFE_DEALLOCATE_A(hs)

  hpsi(1:der%mesh%np, 1:st%d%dim,  st%st_start:st%st_end) = psi(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end)

  call pblas_gemm('n', 'n', total_np, st%nst, st%nst, &
    R_TOTYPE(M_ONE), hpsi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    evectors(1, 1), 1, 1, hs_desc(1), &
    R_TOTYPE(M_ZERO), psi(1, 1, st%st_start), 1, 1, psi_desc(1))

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do ist = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, der, psi(:, :, ist) , hpsi(:, :, st%st_start), ist, ik)
      diff(ist) = X(states_residue)(der%mesh, st%d%dim, hpsi(:, :, st%st_start), eigenval(ist), psi(:, :, ist))
    end do

    if(st%parallel_in_states) then
      ldiff = diff(st%st_start:st%st_end)
      call lmpi_gen_allgatherv(st%lnst, ldiff, tmp, diff(:), st%mpi_grp)
    end if

  end if
  
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(hs)
  
  POP_SUB(X(subspace_diag_scalapack))

#endif /* SCALAPACK */  
end subroutine X(subspace_diag_scalapack)

! ------------------------------------------------------

subroutine X(subspace_test)(st, hm, gr)
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  type(grid_t),        intent(in)    :: gr
  
  type(subspace_t) :: sdiag

  PUSH_SUB(X(subspace_test))

  call states_allocate_wfns(st, gr%mesh, wfs_type = R_TYPE_VAL)
  call states_generate_random(st, gr%mesh)

  call subspace_init(sdiag, st)

#ifdef R_TCOMPLEX
  st%zpsi = (M_ONE + M_ZI)*st%zpsi
#endif
  
  hm%vxc = M_ZERO
  
  call hamiltonian_update(hm, gr%mesh)

  message(1) = 'Info: Orthogonalizing random wavefunctions.'
  message(2) = ''
  call messages_info(2)

  call X(states_orthogonalization_full)(st, st%nst, gr%mesh, st%d%dim, st%X(psi)(:, :, :, 1))
  call X(subspace_diag)(sdiag, gr%der, st, hm, 1, st%eigenval(:, 1), st%X(psi)(:, :, :, 1))

  call print_results()

  call subspace_end(sdiag)
  call states_deallocate_wfns(st)
  
  POP_SUB(X(subspace_test))

contains
  subroutine print_results()
    integer :: ist, jst
    R_TYPE :: dd, hme
    R_TYPE, allocatable :: psi1(:, :), psi2(:, :), hpsi(:, :)
#ifdef HAVE_MPI
    integer :: req(1:4), nreq
#endif

    PUSH_SUB(X(subspace_test).print_results)

    SAFE_ALLOCATE(psi1(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(psi2(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))

    message(1) = 'Residuals:'
    call messages_info(1)
    
    do ist = 1, st%nst
      do jst = ist, st%nst
        if(.not. st%parallel_in_states) then
          psi1(1:gr%mesh%np, 1:st%d%dim) = st%X(psi)(1:gr%mesh%np, 1:st%d%dim, ist, 1)
          psi2(1:gr%mesh%np, 1:st%d%dim) = st%X(psi)(1:gr%mesh%np, 1:st%d%dim, jst, 1)
        end if

#ifdef HAVE_MPI
        if(st%parallel_in_states) then

          ! we bring the two vectors to node 0 and calculate the dot
          ! product, this is very simple and very slow, we only do it for testing
          
          nreq = 0

          ! post the receptions
          if(st%mpi_grp%rank == 0) then
            call MPI_Irecv(psi1(1, 1), gr%mesh%np_part*st%d%dim, R_MPITYPE, st%node(ist), ist, &
              st%mpi_grp%comm, req(nreq + 1), mpi_err)
            call MPI_Irecv(psi2(1, 1), gr%mesh%np_part*st%d%dim, R_MPITYPE, st%node(jst), jst, &
              st%mpi_grp%comm, req(nreq + 2), mpi_err)
            INCR(nreq, 2)
          end if

          ! if I have the wave function, I send it (note: a node could be sending to itself, this is by design)
          if(st%node(ist)  == st%mpi_grp%rank) then
            INCR(nreq, 1)
            call MPI_Isend(st%X(psi)(1, 1, ist, 1), gr%mesh%np_part*st%d%dim, R_MPITYPE, 0, ist, &
              st%mpi_grp%comm, req(nreq), mpi_err)
          end if
          
          if(st%node(jst) == st%mpi_grp%rank) then
            INCR(nreq, 1)
            call MPI_Isend(st%X(psi)(1, 1, jst, 1), gr%mesh%np_part*st%d%dim, R_MPITYPE, 0, jst, &
              st%mpi_grp%comm, req(nreq), mpi_err)
          end if

          if(nreq > 0) call MPI_Waitall(nreq, req(1), MPI_STATUSES_IGNORE, mpi_err)

          if(st%mpi_grp%rank /= 0) cycle

        end if
#endif
        dd = X(mf_dotp)(gr%mesh, st%d%dim, psi1, psi2)
        call X(hamiltonian_apply)(hm, gr%der, psi2, hpsi, jst, 1)
        hme = X(mf_dotp)(gr%mesh, st%d%dim, psi1, hpsi)        

        write (message(1), '(2i7, 2e16.6)') ist, jst, abs(dd), abs(hme)
        call messages_info(1)

      end do
    end do
    
    message(1) = ''
    call messages_info(1)

    SAFE_DEALLOCATE_A(psi1)
    SAFE_DEALLOCATE_A(psi2)
    SAFE_DEALLOCATE_A(hpsi)

    POP_SUB(X(subspace_test).print_results)
  end subroutine print_results
end subroutine X(subspace_test)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
