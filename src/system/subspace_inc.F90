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
!> This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag)(this, namespace, mesh, st, hm, ik, eigenval, diff)
  type(subspace_t),            intent(in)    :: this
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(states_elec_t), target, intent(inout) :: st
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  FLOAT,                       intent(out)   :: eigenval(:)
  FLOAT, optional,             intent(out)   :: diff(:)

  integer :: ist
  R_TYPE, allocatable :: psi(:, :, :)
    
  PUSH_SUB(X(subspace_diag))
  call profiling_in(diagon_prof, "SUBSPACE_DIAG")

  select case(this%method)
    
  case(OPTION__SUBSPACEDIAGONALIZATION__SCALAPACK)

    SAFE_ALLOCATE(psi(1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

    do ist = st%st_start, st%st_end
      call states_elec_get_state(st, mesh, ist, ik, psi(:, :, ist))
    end do

    call X(subspace_diag_scalapack)(namespace, mesh, st, hm, ik, eigenval, psi, diff)

    do ist = st%st_start, st%st_end
      call states_elec_set_state(st, mesh, ist, ik, psi(:, :, ist))
    end do
    
    SAFE_DEALLOCATE_A(psi)
    
  case(OPTION__SUBSPACEDIAGONALIZATION__STANDARD)
    call X(subspace_diag_standard)(namespace, mesh, st, hm, ik, eigenval, diff)
    
  case(OPTION__SUBSPACEDIAGONALIZATION__NONE)
    ! do nothing

  case default
    ASSERT(.false.)
    
  end select

  if(present(diff) .and. st%parallel_in_states) then
    call states_elec_parallel_gather(st, diff)
  end if

  call profiling_out(diagon_prof)
  POP_SUB(X(subspace_diag))
end subroutine X(subspace_diag)

! ---------------------------------------------------------
!> This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag_standard)(namespace, mesh, st, hm, ik, eigenval, diff)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(states_elec_t), target, intent(inout) :: st
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  FLOAT,                       intent(out)   :: eigenval(:)
  FLOAT, optional,             intent(out)   :: diff(:)

  R_TYPE, allocatable :: hmss(:, :), rdiff(:)
  integer             :: ib, minst, maxst
  type(batch_t)       :: hpsib
  type(profile_t), save :: prof_diff
  
  PUSH_SUB(X(subspace_diag_standard))

  SAFE_ALLOCATE(hmss(1:st%nst, 1:st%nst))
  
  call X(subspace_diag_hamiltonian)(namespace, mesh, st, hm, ik, hmss)
  
  ! Diagonalize the Hamiltonian in the subspace.
  ! only half of hmss has the matrix, but this is what Lapack needs
  call lalg_eigensolve(st%nst, hmss, eigenval)
  
#ifdef HAVE_MPI
  ! the eigenvectors are not unique due to phases and degenerate subspaces, but
  ! they must be consistent among processors in domain parallelization
  if (mesh%parallel_in_domains) &
      call MPI_Bcast(hmss, st%nst**2, R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif

  ! Calculate the new eigenfunctions as a linear combination of the
  ! old ones.
  call states_elec_rotate(st, namespace, mesh, hmss, ik)
  
  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 

    call profiling_in(prof_diff, 'SUBSPACE_DIFF')
    
    SAFE_ALLOCATE(rdiff(1:st%nst))
    rdiff(1:st%nst) = R_TOTYPE(M_ZERO)
    
    do ib = st%group%block_start, st%group%block_end
      
      minst = states_elec_block_min(st, ib)
      maxst = states_elec_block_max(st, ib)

      if(hamiltonian_elec_apply_packed(hm)) call st%group%psib(ib, ik)%do_pack
      
      call st%group%psib(ib, ik)%copy_to(hpsib)

      call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), hpsib, ik)
      call batch_axpy(mesh%np, -eigenval, st%group%psib(ib, ik), hpsib)
      call X(mesh_batch_dotp_vector)(mesh, hpsib, hpsib, rdiff(minst:maxst), reduce = .false.)

      call hpsib%end()

      if(hamiltonian_elec_apply_packed(hm)) call st%group%psib(ib, ik)%do_unpack(copy = .false.)
      
    end do

    if (mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, rdiff)
    diff(1:st%nst) = sqrt(abs(rdiff(1:st%nst)))

    SAFE_DEALLOCATE_A(rdiff)

    call profiling_out(prof_diff)
    
  end if

  SAFE_DEALLOCATE_A(hmss)

  POP_SUB(X(subspace_diag_standard))

end subroutine X(subspace_diag_standard)

! --------------------------------------------------------- 
!> This routine diagonalises the Hamiltonian in the subspace defined by
!! the states; this version is aware of parallelization in states but
!! consumes more memory.
subroutine X(subspace_diag_scalapack)(namespace, mesh, st, hm, ik, eigenval, psi, diff)
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  type(hamiltonian_elec_t), intent(in)    :: hm
  integer,                  intent(in)    :: ik
  FLOAT,                    intent(out)   :: eigenval(:)
  R_TYPE,                   intent(inout) :: psi(:, :, st%st_start:)
  FLOAT, optional,          intent(out)   :: diff(:)
 
#ifdef HAVE_SCALAPACK
  R_TYPE, allocatable :: hs(:, :), hpsi(:, :, :), evectors(:, :)
  integer             :: ist, size
  integer :: psi_block(1:2), total_np, psi_desc(BLACS_DLEN), hs_desc(BLACS_DLEN), info
  integer :: nbl, nrow, ncol, ip, idim
  type(batch_t) :: psib, hpsib
  type(profile_t), save :: prof_diag, prof_gemm1, prof_gemm2
#ifdef HAVE_ELPA
  class(elpa_t), pointer :: elpa
#else
  integer :: lwork
  R_TYPE :: rttmp
  R_TYPE, allocatable :: work(:)
#ifdef R_TCOMPLEX
  integer :: lrwork
  CMPLX, allocatable :: rwork(:)
  CMPLX :: ftmp
#endif
#endif
  
  PUSH_SUB(X(subspace_diag_scalapack))

  SAFE_ALLOCATE(hpsi(1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))
  
  call states_elec_parallel_blacs_blocksize(st, namespace, mesh, psi_block, total_np)

  call descinit(psi_desc(1), total_np, st%nst, psi_block(1), psi_block(2), 0, 0,  st%dom_st_proc_grid%context, &
    st%d%dim*mesh%np_part, info)

  if(info /= 0) then
    write(message(1), '(a,i6)') "subspace diagonalization descinit for psi failed with error code ", info
    call messages_fatal(1, namespace=namespace)
  end if

  ! select the blocksize, we use the division used for state
  ! parallelization but with a maximum of 64
  nbl = min(64, psi_block(2))

  ! calculate the size of the matrix in each node
  nrow = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
  ncol = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

  SAFE_ALLOCATE(hs(1:nrow, 1:ncol))

  call descinit(hs_desc(1), st%nst, st%nst, nbl, nbl, 0, 0, st%dom_st_proc_grid%context, nrow, info)

  if(info /= 0) then
    write(message(1), '(a,i6)') "subspace diagonalization descinit for Hamiltonian failed with error code ", info
    call messages_fatal(1, namespace=namespace)
  end if

  ! calculate |hpsi> = H |psi>
  do ist = st%st_start, st%st_end, st%d%block_size
    size = min(st%d%block_size, st%st_end - ist + 1)
    
    call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
    call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, hpsi(: , :, ist:))
    
    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, psib, hpsib, ik)
    
    call psib%end()
    call hpsib%end()
  end do

  ! We need to set to zero some extra parts of the array
  if(st%d%dim == 1) then
    psi(mesh%np + 1:psi_block(1), 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    hpsi(mesh%np + 1:psi_block(1), 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
  else
    psi(mesh%np + 1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    hpsi(mesh%np + 1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
  end if
  
  call profiling_in(prof_gemm1, "SCALAPACK_GEMM1")

  ! get the matrix <psi|H|psi> = <psi|hpsi>
  call pblas_gemm('c', 'n', st%nst, st%nst, total_np, &
    R_TOTYPE(mesh%vol_pp(1)), psi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    hpsi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    R_TOTYPE(M_ZERO), hs(1, 1), 1, 1, hs_desc(1))

  SAFE_ALLOCATE(evectors(1:nrow, 1:ncol))
  call profiling_out(prof_gemm1)

  call profiling_in(prof_diag, "SCALAPACK_DIAG")

  ! now diagonalize
#ifdef HAVE_ELPA
  if (elpa_init(20170403) /= elpa_ok) then
    write(message(1),'(a)') "ELPA API version not supported"
    call messages_fatal(1, namespace=namespace)
  endif
  elpa => elpa_allocate()

  ! set parameters describing the matrix
  call elpa%set("na", st%nst, info)
  call elpa%set("nev", st%nst, info)
  call elpa%set("local_nrows", nrow, info)
  call elpa%set("local_ncols", ncol, info)
  call elpa%set("nblk", nbl, info)
  call elpa%set("mpi_comm_parent", st%dom_st_mpi_grp%comm, info)
  call elpa%set("process_row", st%dom_st_proc_grid%myrow, info)
  call elpa%set("process_col", st%dom_st_proc_grid%mycol, info)

  info = elpa%setup()

  ! one stage solver usually shows worse performance than two stage solver
  call elpa%set("solver", elpa_solver_2stage, info)

  ! call eigensolver
  call elpa%eigenvectors(hs, eigenval, evectors, info)

  ! error handling
  if (info /= elpa_ok) then
    write(message(1),'(a,i6,a,a)') "Error in ELPA, code: ", info, ", message: ", &
      elpa_strerr(info)
    call messages_fatal(1, namespace=namespace)
  end if

  call elpa_deallocate(elpa)
  call elpa_uninit()

#else
! Use ScaLAPACK function if ELPA not available

#ifdef R_TCOMPLEX

  call pzheev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), &
    work = rttmp, lwork = -1, rwork = ftmp, lrwork = -1, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pzheev workspace query failure, error code = ", info
    call messages_fatal(1, namespace=namespace)
  end if

  lwork = nint(abs(rttmp))
  lrwork = nint(real(ftmp, 8))

  SAFE_ALLOCATE(work(1:lwork))
  SAFE_ALLOCATE(rwork(1:lrwork))

  call pzheev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), &
    work = work(1), lwork = lwork, rwork = rwork(1), lrwork = lrwork, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pzheev call failure, error code = ", info
    call messages_fatal(1, namespace=namespace)
  end if

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(rwork)

#else

  call pdsyev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), work = rttmp, lwork = -1, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pdsyev workspace query failure, error code = ", info
    call messages_fatal(1, namespace=namespace)
  end if

  lwork = nint(abs(rttmp))
  SAFE_ALLOCATE(work(1:lwork))

  call pdsyev(jobz = 'V', uplo = 'U', n = st%nst, a = hs(1, 1) , ia = 1, ja = 1, desca = hs_desc(1), &
    w = eigenval(1), z = evectors(1, 1), iz = 1, jz = 1, descz = hs_desc(1), work = work(1), lwork = lwork, info = info)

  if(info /= 0) then
    write(message(1),'(a,i6)') "ScaLAPACK pdsyev call failure, error code = ", info
    call messages_fatal(1, namespace=namespace)
  end if
  
  SAFE_DEALLOCATE_A(work)
#endif
  
#endif
!(HAVE_ELPA)

  call profiling_out(prof_diag)

  SAFE_DEALLOCATE_A(hs)

  !$omp parallel private(ist, idim, ip)
  do ist = st%st_start, st%st_end
    do idim = 1, st%d%dim
      !$omp do
      do ip = 1, mesh%np
        hpsi(ip, idim, ist) = psi(ip, idim, ist)
      end do
      !$omp end do nowait
    end do
  end do
  !$omp end parallel

  call profiling_in(prof_gemm2, "SCALAPACK_GEMM2")
  call pblas_gemm('n', 'n', total_np, st%nst, st%nst, &
    R_TOTYPE(M_ONE), hpsi(1, 1, st%st_start), 1, 1, psi_desc(1), &
    evectors(1, 1), 1, 1, hs_desc(1), &
    R_TOTYPE(M_ZERO), psi(1, 1, st%st_start), 1, 1, psi_desc(1))
  call profiling_out(prof_gemm2)

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do ist = st%st_start, st%st_end
      call X(hamiltonian_elec_apply)(hm, namespace, mesh, psi(:, :, ist) , hpsi(:, :, st%st_start), ist, ik)
      diff(ist) = X(states_elec_residue)(mesh, st%d%dim, hpsi(:, :, st%st_start), eigenval(ist), psi(:, :, ist))
    end do
  end if
  
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(hs)
  
  POP_SUB(X(subspace_diag_scalapack))

#endif /* SCALAPACK */  
end subroutine X(subspace_diag_scalapack)

! ------------------------------------------------------

!> This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag_hamiltonian)(namespace, mesh, st, hm, ik, hmss)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(states_elec_t), target, intent(inout) :: st
  type(hamiltonian_elec_t),    intent(in)    :: hm
  integer,                     intent(in)    :: ik
  R_TYPE,                      intent(out)   :: hmss(:, :)

  integer       :: ib, ip
  R_TYPE, allocatable :: psi(:, :, :), hpsi(:, :, :)
  type(batch_t), allocatable :: hpsib(:)
  integer :: sp, size, block_size
  type(accel_mem_t) :: psi_buffer, hpsi_buffer, hmss_buffer

  PUSH_SUB(X(subspace_diag_hamiltonian))
  call profiling_in(hamiltonian_elec_prof, "SUBSPACE_HAMILTONIAN")

  SAFE_ALLOCATE(hpsib(st%group%block_start:st%group%block_end))
  
  do ib = st%group%block_start, st%group%block_end
    call st%group%psib(ib, ik)%copy_to(hpsib(ib))
    call X(hamiltonian_elec_apply_batch)(hm, namespace, mesh, st%group%psib(ib, ik), hpsib(ib), ik)
  end do
  
  if(st%are_packed() .and. accel_is_enabled()) then

    ASSERT(ubound(hmss, dim = 1) == st%nst)

    call accel_create_buffer(hmss_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*st%nst)
    call accel_set_buffer_to_zero(hmss_buffer, R_TYPE_VAL, st%nst*st%nst)

    if(.not. st%parallel_in_states .and. st%group%block_start == st%group%block_end) then
      ! all the states are stored in one block
      ! we can use blas directly

      call X(accel_gemm)(transA = ACCEL_BLAS_N, transB = ACCEL_BLAS_C, &
        M = int(st%nst, 8), N = int(st%nst, 8), K = int(mesh%np, 8), &
        alpha = R_TOTYPE(mesh%volume_element), &
        A = st%group%psib(st%group%block_start, ik)%pack%buffer, offA = 0_8, &
        lda = int(st%group%psib(st%group%block_start, ik)%pack%size(1), 8), &
        B = hpsib(st%group%block_start)%pack%buffer, offB = 0_8, &
        ldb = int(hpsib(st%group%block_start)%pack%size(1), 8), &
        beta = R_TOTYPE(CNST(0.0)), &
        C = hmss_buffer, offC = 0_8, ldc = int(st%nst, 8))

    else

      ASSERT(.not. st%parallel_in_states)
      
      ! we have to copy the blocks to a temporary array
      block_size = batch_points_block_size(st%group%psib(st%group%block_start, ik))

      call accel_create_buffer(psi_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*block_size)
      call accel_create_buffer(hpsi_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*block_size)

      do sp = 1, mesh%np, block_size
        size = min(block_size, mesh%np - sp + 1)

        do ib = st%group%block_start, st%group%block_end
          ASSERT(R_TYPE_VAL == st%group%psib(ib, ik)%type())
          call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psi_buffer, st%nst)
          call batch_get_points(hpsib(ib), sp, sp + size - 1, hpsi_buffer, st%nst)
        end do

        call X(accel_gemm)(transA = ACCEL_BLAS_N, transB = ACCEL_BLAS_C, &
          M = int(st%nst, 8), N = int(st%nst, 8), K = int(size, 8), &
          alpha = R_TOTYPE(mesh%volume_element), &
          A = psi_buffer, offA = 0_8, lda = int(st%nst, 8), &
          B = hpsi_buffer, offB = 0_8, ldb = int(st%nst, 8), beta = R_TOTYPE(CNST(1.0)), & 
          C = hmss_buffer, offC = 0_8, ldc = int(st%nst, 8))
        
        call accel_finish()

      end do


      call accel_release_buffer(psi_buffer)
      call accel_release_buffer(hpsi_buffer)
      
    end if

    call accel_read_buffer(hmss_buffer, st%nst*st%nst, hmss)
    call accel_release_buffer(hmss_buffer)

  else

#ifdef R_TREAL  
    block_size = max(40, hardware%l2%size/(2*8*st%nst))
#else
    block_size = max(20, hardware%l2%size/(2*16*st%nst))
#endif

    hmss(1:st%nst, 1:st%nst) = CNST(0.0)
    
    SAFE_ALLOCATE(psi(1:st%nst, 1:st%d%dim, 1:block_size))
    SAFE_ALLOCATE(hpsi(1:st%nst, 1:st%d%dim, 1:block_size))

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)
      
      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psi)
        call batch_get_points(hpsib(ib), sp, sp + size - 1, hpsi)
      end do

      if(st%parallel_in_states) then
        call states_elec_parallel_gather(st, (/st%d%dim, size/), psi)
        call states_elec_parallel_gather(st, (/st%d%dim, size/), hpsi)
      end if
      
      if (mesh%use_curvilinear) then
        do ip = 1, size
          psi(1:st%nst, 1:st%d%dim, ip) = psi(1:st%nst, 1:st%d%dim, ip)*mesh%vol_pp(sp + ip - 1)
        end do
      end if

      call blas_gemm(transa = 'n', transb = 'c',        &
        m = st%nst, n = st%nst, k = size*st%d%dim,      &
        alpha = R_TOTYPE(mesh%volume_element),      &
        a = hpsi(1, 1, 1), lda = ubound(hpsi, dim = 1),   &
        b = psi(1, 1, 1), ldb = ubound(psi, dim = 1), &
        beta = R_TOTYPE(CNST(1.0)),                     & 
        c = hmss(1, 1), ldc = ubound(hmss, dim = 1))
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

  end if
  
  call profiling_count_operations((R_ADD + R_MUL)*st%nst*(st%nst - CNST(1.0))*mesh%np)
  
  do ib = st%group%block_start, st%group%block_end
    call hpsib(ib)%end()
  end do
  
  SAFE_DEALLOCATE_A(hpsib)
    
  if (mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, hmss, dim = (/st%nst, st%nst/))
  
  call profiling_out(hamiltonian_elec_prof)
  POP_SUB(X(subspace_diag_hamiltonian))

end subroutine X(subspace_diag_hamiltonian)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
