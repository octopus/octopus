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
!> Orthonormalizes nst orbitals in mesh (honours state parallelization).
subroutine X(states_elec_orthogonalization_full)(st, namespace, mesh, ik)
  type(states_elec_t),    intent(inout) :: st
  type(namespace_t),      intent(in)    :: namespace
  type(mesh_t),           intent(in)    :: mesh
  integer,                intent(in)    :: ik

  R_TYPE, allocatable :: ss(:, :)
  type(profile_t), save :: prof
  integer :: nst

#ifdef HAVE_SCALAPACK
!pgi$r opt=0
!This is a pragma for the PGI compiler, forcing optimization -O0 for this subroutine
!With PGI 10.9 and ScaLAPACK, at -O2 and higher optimization levels, the test finite_systems_3d/10-fullerene fails in
!states_elec_orthogonalization_full.par_gs with error message
!glibc detected *** octopus_mpi: malloc(): memory corruption
#endif

  call profiling_in(prof, TOSTRING(X(ORTHOGONALIZATION_FULL)))
  PUSH_SUB(X(states_elec_orthogonalization_full))

  nst = st%nst

  select case(st%d%orth_method)
  case(OPTION__STATESORTHOGONALIZATION__CHOLESKY_SERIAL)
    call cholesky_serial()

  case(OPTION__STATESORTHOGONALIZATION__CHOLESKY_PARALLEL)
    call cholesky_parallel()

  case(OPTION__STATESORTHOGONALIZATION__CGS, OPTION__STATESORTHOGONALIZATION__MGS, &
       OPTION__STATESORTHOGONALIZATION__DRCGS)
    call mgs()

  case default
    write(message(1),'(a,i6)') "Internal error from states_elec_orthogonalization_full: orth_method has illegal value ", &
      st%d%orth_method
    call messages_fatal(1, namespace=namespace)
  end select

  call profiling_out(prof)
  POP_SUB(X(states_elec_orthogonalization_full))

contains
  
  subroutine cholesky_serial()

    integer :: ierr
    logical :: bof

    PUSH_SUB(X(states_elec_orthogonalization_full).cholesky_serial)

    SAFE_ALLOCATE(ss(1:nst, 1:nst))

    ss = M_ZERO

    call X(states_elec_calc_overlap)(st, mesh, ik, ss)

    bof = .false.
    ! calculate the Cholesky decomposition
    call lalg_cholesky(nst, ss, bof = bof, err_code = ierr)

    if(bof) then
      write(message(1),'(a,i6)') "The cholesky_serial orthogonalization failed with error code ", ierr
      message(2) = "There may be a linear dependence, a zero vector, or maybe a library problem."
      message(3) = "Using the Gram-Schimdt orthogonalization instead."
      call messages_warning(3, namespace=namespace)
    end if
  
    if(.not. bof) then
      call X(states_elec_trsm)(st, namespace, mesh, ik, ss)
    else
      call mgs()
    end if

    SAFE_DEALLOCATE_A(ss)

    POP_SUB(X(states_elec_orthogonalization_full).cholesky_serial)
  end subroutine cholesky_serial


  ! -----------------------------------------------------------------------------------------------
  subroutine cholesky_parallel()

    R_TYPE, allocatable :: psi(:, :, :)
    integer             :: psi_block(1:2), total_np
#ifdef HAVE_SCALAPACK
    integer             :: info, nbl, nrow, ncol
    integer             :: psi_desc(BLACS_DLEN), ss_desc(BLACS_DLEN)
    type(profile_t), save :: prof_cholesky, prof_trsm, prof_herk
#endif

    PUSH_SUB(X(states_elec_orthogonalization_full).cholesky_parallel)

! some checks
#ifndef HAVE_MPI
    message(1) = 'The cholesky_parallel orthogonalizer can only be used in parallel.'
    call messages_fatal(1, namespace=namespace)
#else
#ifndef HAVE_SCALAPACK
    message(1) = 'The cholesky_parallel orthogonalizer requires ScaLAPACK.'
    call messages_fatal(1, only_root_writes = .true., namespace=namespace)
#endif
    if(st%dom_st_mpi_grp%size == 1) then
      message(1) = 'The cholesky_parallel orthogonalizer is designed to be used with domain or state parallelization.'
      call messages_warning(1, namespace=namespace)
    end if
#endif

    call states_elec_parallel_blacs_blocksize(st, namespace, mesh, psi_block, total_np)

    SAFE_ALLOCATE(psi(1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

    call states_elec_get_state(st, mesh, ik, psi)
    
    ! We need to set to zero some extra parts of the array
    if(st%d%dim == 1) then
     psi(mesh%np + 1:psi_block(1), 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    else
     psi(mesh%np + 1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end) = M_ZERO
    end if

#ifdef HAVE_SCALAPACK
    
    call descinit(psi_desc(1), total_np, st%nst, psi_block(1), psi_block(2), 0, 0, st%dom_st_proc_grid%context, &
      st%d%dim*ubound(psi, dim = 1), info)

    if(info /= 0) then
      write(message(1),'(3a,i6)') "descinit for psi failed in ", TOSTRING(X(states_elec_orthogonalization_full)), &
        ".cholesky_parallel with error ", info
      call messages_fatal(1, namespace=namespace)
    end if

    nbl = min(32, st%nst)
    nrow = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
    ncol = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

    SAFE_ALLOCATE(ss(1:nrow, 1:ncol))

    call descinit(ss_desc(1), st%nst, st%nst, nbl, nbl, 0, 0, st%dom_st_proc_grid%context, ubound(ss, dim = 1), info)

    if(info /= 0) then
      write(message(1),'(3a,i6)') "descinit for ss failed in ", TOSTRING(X(states_elec_orthogonalization_full)), &
        ".cholesky_parallel with error ", info
      call messages_fatal(1, namespace=namespace)
    end if

    ss = M_ZERO

    call profiling_in(prof_herk, TOSTRING(X(SCALAPACK_HERK)))
    call pblas_herk(uplo = 'U', trans = 'C', n = st%nst, k = total_np, &
      alpha = R_TOTYPE(mesh%vol_pp(1)), a = psi(1, 1, st%st_start), ia = 1, ja = 1, desca = psi_desc(1), &
      beta = R_TOTYPE(M_ZERO), c = ss(1, 1), ic = 1, jc = 1, descc = ss_desc(1))
    call profiling_count_operations(TOFLOAT(mesh%np*nst)**2*(R_ADD + R_MUL))
    call profiling_out(prof_herk)

    call profiling_in(prof_cholesky, TOSTRING(X(SCALAPACK_CHOLESKY)))
    ! calculate the Cholesky decomposition
    call scalapack_potrf(uplo = 'U', n = st%nst, a = ss(1, 1), ia = 1, ja = 1, desca = ss_desc(1), info = info)
    call profiling_out(prof_cholesky)

    if(info /= 0) then
      write(message(1),'(3a,i6)') "cholesky_parallel orthogonalization with ", TOSTRING(pX(potrf)), &
        " failed with error ", info
      call messages_fatal(1, namespace=namespace)
    end if

    call profiling_in(prof_trsm, TOSTRING(X(SCALAPACK_TRSM)))
    call pblas_trsm(side = 'R', uplo = 'U', transa = 'N', diag = 'N', m = total_np, n = st%nst, &
      alpha = R_TOTYPE(M_ONE), a = ss(1, 1), ia = 1, ja = 1, desca = ss_desc(1), &
      b = psi(1, 1, st%st_start), ib = 1, jb = 1, descb = psi_desc(1))
    call profiling_out(prof_trsm)
#endif

    SAFE_DEALLOCATE_A(ss)

    call states_elec_set_state(st, mesh, ik, psi)
    
    POP_SUB(X(states_elec_orthogonalization_full).cholesky_parallel)
  end subroutine cholesky_parallel

  ! ----------------------------------------------------------------------------------

  subroutine mgs()

    integer :: ist, jst, idim, is
    FLOAT   :: cc
    R_TYPE, allocatable :: aa(:), psii(:, :), psij(:, :)
    R_TYPE, allocatable :: psii0(:, :)
    integer :: method
    
    PUSH_SUB(X(states_elec_orthogonalization_full).mgs)

    ! set method to MGS in case this method is called as bof from cholesky
    method = st%d%orth_method
    if(method == OPTION__STATESORTHOGONALIZATION__CHOLESKY_SERIAL) then
      method = OPTION__STATESORTHOGONALIZATION__MGS
    end if

    if(st%parallel_in_states .and. &
       method /= OPTION__STATESORTHOGONALIZATION__MGS) then
      message(1) = 'The mgs orthogonalization method cannot work with state-parallelization.'
      call messages_fatal(1, only_root_writes = .true., namespace=namespace)
    end if

    SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psii0(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psij(1:mesh%np, 1:st%d%dim))

    SAFE_ALLOCATE(aa(1:nst))

    do ist = 1, nst

      !The different algorithms are given in Giraud et al., 
      !Computers and Mathematics with Applications 50, 1069 (2005).
      select case(method)
      case(OPTION__STATESORTHOGONALIZATION__MGS)

        if(ist >= st%st_start .and. ist <= st%st_end) then
          call states_elec_get_state(st, mesh, ist, ik, psii)

          ! renormalize
          cc = TOFLOAT(X(mf_dotp)(mesh, st%d%dim, psii, psii))
          call lalg_scal(mesh%np, st%d%dim, M_ONE/sqrt(cc), psii)
          call states_elec_set_state(st, mesh, ist, ik, psii)
        end if

        if(st%parallel_in_states) then
#ifdef HAVE_MPI
          call MPI_Bcast(psii(1, 1), mesh%np*st%d%dim, R_MPITYPE, st%node(ist), st%mpi_grp%comm, mpi_err)
#endif
        end if

        aa = M_ZERO

        ! calculate the projections
        do jst = ist + 1, nst
          if(jst < st%st_start .or. jst > st%st_end) cycle
          call states_elec_get_state(st, mesh, jst, ik, psij)
          aa(jst) = X(mf_dotp)(mesh, st%d%dim, psii, psij, reduce = .false.)
        end do
        if(mesh%parallel_in_domains) call mesh%allreduce(aa, dim = nst)
 
        ! subtract the projections
        do jst = ist + 1, nst
          if(jst < st%st_start .or. jst > st%st_end) cycle
          call states_elec_get_state(st, mesh, jst, ik, psij)
          do idim = 1, st%d%dim
            call lalg_axpy(mesh%np, -aa(jst), psii(:, idim), psij(:, idim))
          end do
          call states_elec_set_state(st, mesh, jst, ik, psij)
        end do

      case(OPTION__STATESORTHOGONALIZATION__CGS)

        call states_elec_get_state(st, mesh, ist, ik, psii)

        ! calculate the projections first with the same vector
        do jst = 1, ist - 1
          call states_elec_get_state(st, mesh, jst, ik, psij)
          aa(jst) = X(mf_dotp)(mesh, st%d%dim, psij, psii, reduce = .false.)
        end do

        if(mesh%parallel_in_domains .and. ist > 1) call mesh%allreduce(aa, dim = ist - 1)
        ! then subtract the projections
        do jst = 1, ist - 1
          call states_elec_get_state(st, mesh, jst, ik, psij)
          do idim = 1, st%d%dim
            call lalg_axpy(mesh%np, -aa(jst), psij(:, idim), psii(:, idim))
          end do
        end do

      case(OPTION__STATESORTHOGONALIZATION__DRCGS)

        call states_elec_get_state(st, mesh, ist, ik, psii)

        !double step reorthogonalization
        do is = 1, 2
          if(ist>1) then
            call states_elec_get_state(st, mesh, ist-1, ik, psii0)
            aa(1) = X(mf_dotp)(mesh, st%d%dim, psii0, psii)
            do idim = 1, st%d%dim
              call lalg_axpy(mesh%np, -aa(1), psii0(:, idim), psii(:, idim))
            end do
          end if
      
          ! calculate the projections
          do jst = 1, ist - 1
            call states_elec_get_state(st, mesh, jst, ik, psij)
            aa(jst) = X(mf_dotp)(mesh, st%d%dim, psij, psii, reduce = .false.)
          end do
          if(mesh%parallel_in_domains .and. ist > 1) call mesh%allreduce(aa, dim = ist - 1)

          ! subtract the projections
          do jst = 1, ist - 1
            call states_elec_get_state(st, mesh, jst, ik, psij)
            do idim = 1, st%d%dim
              call lalg_axpy(mesh%np, -aa(jst), psij(:, idim), psii(:, idim))
            end do
          end do
        end do

      end select

      !In case of modified Gram-Schmidt, this was done before.
      if(method == OPTION__STATESORTHOGONALIZATION__CGS .or. &
         method == OPTION__STATESORTHOGONALIZATION__DRCGS) then
        ! renormalize
        cc = TOFLOAT(X(mf_dotp)(mesh, st%d%dim, psii, psii))

        call lalg_scal(mesh%np, st%d%dim, M_ONE/sqrt(cc), psii)

        call states_elec_set_state(st, mesh, ist, ik, psii)
      end if
      
    end do

    SAFE_DEALLOCATE_A(psii)
    SAFE_DEALLOCATE_A(psii0)
    SAFE_DEALLOCATE_A(psij)
    SAFE_DEALLOCATE_A(aa)

    POP_SUB(X(states_elec_orthogonalization_full).mgs)
  end subroutine mgs

end subroutine X(states_elec_orthogonalization_full)


! ---------------------------------------------------------

subroutine X(states_elec_trsm)(st, namespace, mesh, ik, ss)
  type(states_elec_t),    intent(inout) :: st
  type(namespace_t),      intent(in)    :: namespace
  type(mesh_t),           intent(in)    :: mesh
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(in)    :: ss(:, :)

  integer :: idim, block_size, ib, size, sp
  R_TYPE, allocatable :: psicopy(:, :, :)
  type(accel_mem_t) :: psicopy_buffer, ss_buffer
  type(profile_t), save :: prof_copy
  type(profile_t), save :: prof

  PUSH_SUB(X(states_elec_trsm))
  call profiling_in(prof, TOSTRING(X(STATES_TRSM)))

  if(.not. (st%are_packed() .and. accel_is_enabled())) then

#ifdef R_TREAL  
    block_size = max(40, hardware%l2%size/(2*8*st%nst))
#else
    block_size = max(20, hardware%l2%size/(2*16*st%nst))
#endif

    SAFE_ALLOCATE(psicopy(1:st%nst, 1:st%d%dim, 1:block_size))

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy)
      end do

      if(st%parallel_in_states) call states_elec_parallel_gather(st, (/st%d%dim, size/), psicopy)      
      
      do idim = 1, st%d%dim
        
        call blas_trsm(side = 'L', uplo = 'U', transa = 'T', diag = 'N', &
          m = st%nst, n = size, &
          alpha = R_TOTYPE(M_ONE), a = ss(1, 1), lda = ubound(ss, dim = 1), &
          b = psicopy(1, idim, 1), ldb = ubound(psicopy, dim = 1)*st%d%dim)

      end do
      
      do ib = st%group%block_start, st%group%block_end
        call batch_set_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy)
      end do

    end do 
    
    SAFE_DEALLOCATE_A(psicopy)

  else

    if(st%d%dim > 1) call messages_not_implemented('Opencl states_elec_trsm for spinors', namespace=namespace)

    block_size = batch_points_block_size()

    call accel_create_buffer(psicopy_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*block_size)

    call accel_create_buffer(ss_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, product(ubound(ss)))

    call profiling_in(prof_copy, TOSTRING(X(STATES_TRSM_COPY)))
    call accel_write_buffer(ss_buffer, product(ubound(ss)), ss)
    call profiling_count_transfers(product(ubound(ss)), ss(1, 1))

    call profiling_out(prof_copy)

    if(st%parallel_in_states) then
      SAFE_ALLOCATE(psicopy(1:st%nst, 1:st%d%dim, 1:block_size))
    end if

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      do ib = st%group%block_start, st%group%block_end
        ASSERT(R_TYPE_VAL == st%group%psib(ib, ik)%type())
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy_buffer, st%nst)
      end do

      if(st%parallel_in_states) then
        call accel_read_buffer(psicopy_buffer, st%nst*st%d%dim*block_size, psicopy)
        call states_elec_parallel_gather(st, (/st%d%dim, size/), psicopy)
        call accel_write_buffer(psicopy_buffer, st%nst*st%d%dim*block_size, psicopy)
      end if

      call X(accel_trsm)(side = ACCEL_BLAS_LEFT, uplo = ACCEL_BLAS_UPPER, &
        trans = ACCEL_BLAS_T, diag = ACCEL_BLAS_DIAG_NON_UNIT, &
        M = int(st%nst, 8), N = int(size, 8), alpha = R_TOTYPE(M_ONE), &
        A = ss_buffer, offA = 0_8, lda = int(ubound(ss, dim = 1), 8), &
        B = psicopy_buffer, offB = 0_8, ldb = int(st%nst, 8))
      
      do ib = st%group%block_start, st%group%block_end
        call batch_set_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy_buffer, st%nst)
      end do
    end do

    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_A(psicopy)
    end if

    call accel_release_buffer(ss_buffer)
    call accel_release_buffer(psicopy_buffer)

  end if

  call profiling_count_operations(mesh%np*TOFLOAT(st%nst*(st%nst + 1))*st%d%dim*CNST(0.5)*(R_ADD + R_MUL))


  call profiling_out(prof)
  POP_SUB(X(states_elec_trsm))
end subroutine X(states_elec_trsm)

! ---------------------------------------------------------
subroutine X(states_elec_orthogonalize_single)(st, mesh, nst, iqn, phi, normalize, mask, overlap, norm, Theta_fi, beta_ij, &
  against_all)
  type(states_elec_t), target, intent(in)    :: st
  type(mesh_t),                intent(in)    :: mesh
  integer,                     intent(in)    :: nst
  integer,                     intent(in)    :: iqn
  R_TYPE,                      intent(inout) :: phi(:,:)     !< phi(mesh%np_part, dim)
  logical,           optional, intent(in)    :: normalize
  logical,           optional, intent(inout) :: mask(:)      !< mask(nst)
  R_TYPE,            optional, intent(out)   :: overlap(:) 
  FLOAT,             optional, intent(out)   :: norm
  FLOAT,             optional, intent(in)    :: theta_fi
  R_TYPE,            optional, intent(in)    :: beta_ij(:)   !< beta_ij(nst)
  logical,           optional, intent(in)    :: against_all

  integer :: ist, idim, length_ss, ibind
  FLOAT   :: nrm2
  R_TYPE, allocatable  :: ss(:), psi(:, :)
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof
  logical :: against_all_
  type(wfs_elec_t), pointer :: batch
  
  call profiling_in(prof, TOSTRING(X(GRAM_SCHMIDT)))
  PUSH_SUB(X(states_elec_orthogonalize_single))

  ASSERT(nst <= st%nst)
  ASSERT(.not. st%parallel_in_states)
  ! if against_all is set to true, phi is orthogonalized to all other states except nst+1
  ! (nst + 1 is chosen because this routine is usually called with nst=ist-1 in a loop)
  against_all_ = optional_default(against_all, .false.)

  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  length_ss = nst
  if(against_all_) then
    length_ss = st%nst
  end if
  SAFE_ALLOCATE(ss(1:length_ss))
  ! Check length of optional arguments
  if(present(mask)) then
    ASSERT(ubound(mask, dim=1) >= length_ss)
  end if
  if(present(overlap)) then
    ASSERT(ubound(overlap, dim=1) >= length_ss)
  end if
  if(present(beta_ij)) then
    ASSERT(ubound(beta_ij, dim=1) >= length_ss)
  end if

  ss = M_ZERO

  do ist = 1, st%nst
    if(skip_this_iteration(ist, nst, against_all_)) cycle
    if(present(mask)) then
      if(mask(ist)) cycle
    end if
 
    !To understand this, one should look at states_elec_get_states and batch_get_states routines 
    batch => st%group%psib(st%group%iblock(ist, iqn), iqn)
    select case(batch%status())
    case(BATCH_NOT_PACKED)
      ss(ist) = R_TOTYPE(M_ZERO)
      do idim = 1, st%d%dim
        ibind = batch%inv_index((/ist, idim/)) 
        ss(ist) = ss(ist) + X(mf_dotp)(mesh, batch%X(ff_linear)(:, ibind), phi(:,idim), reduce = .false.)
      end do
    case(BATCH_PACKED, BATCH_DEVICE_PACKED)
      !Not properly implemented
      !We need to reorder the operations is these two cases
      call states_elec_get_state(st, mesh, ist, iqn, psi)
      ss(ist) = X(mf_dotp)(mesh, st%d%dim, psi, phi, reduce = .false.)
    end select
  end do
    
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(GRAM_SCHMIDT_REDUCE)))
    call mesh%allreduce(ss, dim = length_ss)
    call profiling_out(reduce_prof)
  end if

  if(present(mask)) then
    do ist = 1, st%nst
      if(skip_this_iteration(ist, nst, against_all_)) cycle
      mask(ist) = (abs(ss(ist)) <= M_EPSILON)
    end do
  end if

  if(present(beta_ij)) ss(1:length_ss) = ss(1:length_ss)*beta_ij(1:length_ss)
  
  if(present(theta_fi)) then
    if(theta_fi /= M_ONE) phi(1:mesh%np, 1:st%d%dim) = theta_fi*phi(1:mesh%np, 1:st%d%dim)
  end if

  do ist = 1, st%nst
    if(skip_this_iteration(ist, nst, against_all_)) cycle
    if(present(mask)) then
      if(mask(ist)) cycle
    end if
    
    batch => st%group%psib(st%group%iblock(ist, iqn), iqn)
    select case(batch%status())
    case(BATCH_NOT_PACKED)
      do idim = 1, st%d%dim
        ibind = batch%inv_index((/ist, idim/))
        call blas_axpy(mesh%np, -ss(ist), batch%X(ff_linear)(1, ibind), 1, phi(1, idim), 1)
      end do
    case(BATCH_PACKED, BATCH_DEVICE_PACKED)
      !Not properly implemented
      !We need to reorder the operations is these two cases
      call states_elec_get_state(st, mesh, ist, iqn, psi)
      do idim = 1, st%d%dim
        call blas_axpy(mesh%np, -ss(ist), psi(1, idim), 1, phi(1, idim), 1)
      end do
    end select
  end do

  if(optional_default(normalize, .false.)) then
    call X(mf_normalize)(mesh, st%d%dim, phi, nrm2)
  end if

  if(present(overlap)) then
    overlap(1:length_ss) = ss(1:length_ss)
  end if

  if(present(norm)) then
    ASSERT(present(normalize))
    ASSERT(normalize)
    norm = nrm2
  end if

  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(states_elec_orthogonalize_single))
  call profiling_out(prof)

  contains

    logical function skip_this_iteration(ist, nst, against_all_states)
      integer, intent(in) :: ist, nst
      logical, intent(in) :: against_all_states

      skip_this_iteration = .false.
      if(.not.against_all_states) then
        ! orthogonalize against previous states only
        if(ist > nst) skip_this_iteration = .true.
      else
        ! orthogonalize against all other states besides nst + 1
        if(ist == nst + 1) skip_this_iteration = .true.
      end if
    end function skip_this_iteration
end subroutine X(states_elec_orthogonalize_single)

! ---------------------------------------------------------
subroutine X(states_elec_orthogonalize_single_batch)(st, mesh, nst, iqn, phi, normalize, mask, overlap, norm, Theta_fi, beta_ij, &
  against_all)
  type(states_elec_t), intent(in)    :: st
  type(mesh_t),        intent(in)    :: mesh
  integer,             intent(in)    :: nst
  integer,             intent(in)    :: iqn
  R_TYPE,              intent(inout) :: phi(:,:)     !< phi(mesh%np_part, dim)
  logical, optional,   intent(in)    :: normalize
  logical, optional,   intent(inout) :: mask(:)      !< mask(nst)
  R_TYPE,  optional,   intent(out)   :: overlap(:) 
  FLOAT,   optional,   intent(out)   :: norm
  FLOAT,   optional,   intent(in)    :: theta_fi
  R_TYPE,  optional,   intent(in)    :: beta_ij(:)   !< beta_ij(nst)
  logical, optional,   intent(in)    :: against_all

  integer :: ib, minst, maxst, ist, length_ss
  FLOAT   :: nrm2
  R_TYPE, allocatable  :: ss(:)
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof
  logical :: against_all_
  
  call profiling_in(prof, TOSTRING(X(GRAM_SCHMIDT_BATCH)))
  PUSH_SUB(X(states_elec_orthogonalize_single_batch))

  ASSERT(nst <= st%nst)
  ASSERT(.not. st%parallel_in_states)
  ! if against_all is set to true, phi is orthogonalized to all other states except nst+1
  ! (nst + 1 is chosen because this routine is usually called with nst=ist-1 in a loop)
  against_all_ = optional_default(against_all, .false.)

  length_ss = nst
  if(against_all_) then
    length_ss = st%nst
  end if
  SAFE_ALLOCATE(ss(1:length_ss))
  ! Check length of optional arguments
  if(present(mask)) then
    ASSERT(ubound(mask, dim=1) >= length_ss)
  end if
  if(present(overlap)) then
    ASSERT(ubound(overlap, dim=1) >= length_ss)
  end if
  if(present(beta_ij)) then
    ASSERT(ubound(beta_ij, dim=1) >= length_ss)
  end if

  ss = R_TOTYPE(M_ZERO)

  do ib = st%group%block_start, st%group%block_end
    minst = states_elec_block_min(st, ib)
    maxst = min(states_elec_block_max(st, ib), length_ss)
    if(minst > length_ss) cycle

    if(skip_this_batch(minst, maxst, nst, against_all_)) cycle
    if(present(mask)) then
      if(all(mask(minst:maxst))) cycle
    end if
 
    call X(mesh_batch_mf_dotp)(mesh, st%group%psib(ib, iqn), phi, ss(minst:maxst), reduce = .false., nst = maxst-minst+1)

    !In case some of the states in the batche need to be skipped
    do ist = minst, maxst
      if(ist > length_ss) cycle
      if(skip_this_iteration(ist, nst, against_all_)) ss(ist) = R_TOTYPE(M_ZERO)
      if(present(mask)) then
        if(mask(ist)) ss(ist) = R_TOTYPE(M_ZERO)
      end if
    end do    
  end do
    
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(GRAM_SCHMIDT_REDUCE)))
    call mesh%allreduce(ss, dim = length_ss)
    call profiling_out(reduce_prof)
  end if

  if(present(mask)) then
    do ist = 1, st%nst
      if(skip_this_iteration(ist, nst, against_all_)) cycle
      mask(ist) = (abs(ss(ist)) <= M_EPSILON)
    end do
  end if

  if(present(beta_ij)) ss(1:length_ss) = ss(1:length_ss)*beta_ij(1:length_ss)
  
  if(present(theta_fi)) then
    if(theta_fi /= M_ONE) phi(1:mesh%np, 1:st%d%dim) = theta_fi*phi(1:mesh%np, 1:st%d%dim)
  end if

  do ib = st%group%block_start, st%group%block_end
    minst = states_elec_block_min(st, ib)
    maxst = min(states_elec_block_max(st, ib), length_ss)
    if(minst > length_ss) cycle

    if(skip_this_batch(minst, maxst, nst, against_all_)) cycle
    if(present(mask)) then
      if(all(mask(minst:maxst))) cycle
    end if

    call X(batch_axpy_function)(mesh%np, -ss(minst:maxst), st%group%psib(ib, iqn), phi, nst = maxst-minst+1) 

  end do

  if(optional_default(normalize, .false.)) then
    call X(mf_normalize)(mesh, st%d%dim, phi, nrm2)
  end if

  if(present(overlap)) then
    overlap(1:length_ss) = ss(1:length_ss)
  end if

  if(present(norm)) then
    ASSERT(present(normalize))
    ASSERT(normalize)
    norm = nrm2
  end if

  SAFE_DEALLOCATE_A(ss)

  POP_SUB(X(states_elec_orthogonalize_single_batch))
  call profiling_out(prof)

  contains
   logical function skip_this_iteration(ist, nst, against_all_states)
      integer, intent(in) :: ist, nst
      logical, intent(in) :: against_all_states

      skip_this_iteration = .false.
      if(.not.against_all_states) then
        ! orthogonalize against previous states only
        if(ist > nst) skip_this_iteration = .true.
      else
        ! orthogonalize against all other states besides nst + 1
        if(ist == nst + 1) skip_this_iteration = .true.
      end if
    end function skip_this_iteration

    logical function skip_this_batch(minst, maxst, nst, against_all_states)
      integer, intent(in) :: minst, maxst, nst
      logical, intent(in) :: against_all_states

      skip_this_batch = .false.
      if(.not.against_all_states) then
        ! orthogonalize against previous states only
        if(minst > nst) skip_this_batch = .true.
      else
        ! orthogonalize against all other states besides nst + 1
        if(minst == nst + 1 .and. maxst == nst + 1) skip_this_batch = .true.
      end if
    end function skip_this_batch
end subroutine X(states_elec_orthogonalize_single_batch)

! ---------------------------------------------------------
!> Orthonormalizes phi to the nst orbitals psi.
!! It also permits doing only the orthogonalization (no normalization).
!! And one can pass an extra optional argument, mask, which:
!!  - on input, if mask(p) = .true., the p-orbital is not used.
!!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
!! If Theta_Fi and beta_ij are present, it performs the generalized orthogonalization
!!   (Theta_Fi - sum_j beta_ij |j><j|Phi> as in De Gironcoli PRB 51, 6774 (1995).
!! This is used in response for metals
subroutine X(states_elec_orthogonalization)(mesh, nst, dim, psi, phi,  &
  normalize, mask, overlap, norm, Theta_fi, beta_ij, gs_scheme)
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: nst
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:,:,:)   !< psi(mesh%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     !< phi(mesh%np_part, dim)
  logical, optional, intent(in)    :: normalize
  logical, optional, intent(inout) :: mask(:)      !< mask(nst)
  R_TYPE,  optional, intent(out)   :: overlap(:) 
  R_TYPE,  optional, intent(out)   :: norm
  FLOAT,   optional, intent(in)    :: Theta_Fi
  R_TYPE,  optional, intent(in)    :: beta_ij(:)   !< beta_ij(nst)
  integer, optional, intent(in)    :: gs_scheme

  logical :: normalize_
  integer :: ist, idim, is
  FLOAT   :: nrm2
  R_TYPE, allocatable  :: ss(:), ss_full(:)
  integer :: block_size, size, sp, ep
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof
  logical :: drcgs
  integer :: nsteps

  call profiling_in(prof, TOSTRING(X(GRAM_SCHMIDT)))
  PUSH_SUB(X(states_elec_orthogonalization))

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)

  SAFE_ALLOCATE(ss(1:nst))

  ss = R_TOTYPE(M_ZERO)

  drcgs = .false.
  nsteps = 1
  if(present(gs_scheme)) then
    if(gs_scheme == OPTION__ARNOLDIORTHOGONALIZATION__DRCGS) then
      drcgs = .true.
      nsteps = 2
      SAFE_ALLOCATE(ss_full(1:nst))
      ss_full = R_TOTYPE(M_ZERO)
    end if
  end if

  do is = 1, nsteps
    if(nst>=1 .and. drcgs) then
      ss(1) = X(mf_dotp)(mesh, dim, psi(:, :, nst), phi)
      do idim = 1, dim
        call lalg_axpy(mesh%np, -ss(1), psi(:, idim, nst), phi(:, idim))
      end do
      call profiling_count_operations((R_ADD + R_MUL) * mesh%np * dim * 2)
      if(present(overlap)) ss_full(nst) = ss_full(nst) + ss(1)
    end if
    ss = R_TOTYPE(M_ZERO)

    if(.not. mesh%use_curvilinear) then

      do sp = 1, mesh%np, block_size
        size = min(block_size, mesh%np - sp + 1)
        do ist = 1, nst
          do idim = 1, dim
        
            if(present(mask)) then
              if(mask(ist)) cycle
            end if

            ss(ist) = ss(ist) + blas_dot(size, psi(sp, idim, ist), 1, phi(sp, idim), 1)
          end do
        end do
      end do

      ss = ss * mesh%vol_pp(1)

      call profiling_count_operations((R_ADD + R_MUL) * mesh%np * dim * nst)

    else

      do sp = 1, mesh%np, block_size
        size = min(block_size, mesh%np - sp + 1)
        ep = sp - 1 + size
        do ist = 1, nst
          do idim = 1, dim

            if(present(mask)) then
              if(mask(ist)) cycle
            end if

            ss(ist) = ss(ist) + sum(mesh%vol_pp(sp:ep)*R_CONJ(psi(sp:ep, idim, ist))*phi(sp:ep, idim))

          end do
        end do
      end do

      call profiling_count_operations((R_ADD + 2 * R_MUL) * mesh%np * dim * nst)

    end if

    if(mesh%parallel_in_domains) then
      call profiling_in(reduce_prof, TOSTRING(X(GRAM_SCHMIDT_REDUCE)))
      call mesh%allreduce(ss, dim = nst)
      call profiling_out(reduce_prof)
    end if

    if(present(mask)) then
      do ist = 1, nst
        mask(ist) = (abs(ss(ist)) <= M_EPSILON)
      end do
    end if

    if(present(beta_ij))  &
      ss(:) = ss(:) * beta_ij(:)

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      if(present(Theta_Fi)) then
        if(Theta_Fi /= M_ONE) &
          call blas_scal(size, R_TOTYPE(Theta_Fi), phi(sp, idim), 1)
      end if

      do idim = 1, dim

        do ist = 1, nst

          if(present(mask)) then
            if(mask(ist)) cycle
          end if

          call blas_axpy(size, -ss(ist), psi(sp, idim, ist), 1, phi(sp, idim), 1)

        end do
      end do
    end do

    call profiling_count_operations((R_ADD + R_MUL) * mesh%np * dim * nst)

    !We accumulate the overlap
    if(drcgs .and. present(overlap)) then
      do ist = 1, nst
        ss_full(ist) = ss_full(ist) + ss(ist)
      end do 
    end if
  end do

  ! the following ifs cannot be given as a single line (without the
  ! then) to avoid a bug in xlf 10.1

  normalize_ = .false.
  if(present(normalize)) then
    normalize_ = normalize
  end if

  if(normalize_) then
    call X(mf_normalize)(mesh, dim, phi, nrm2)
  end if

  if(present(overlap)) then
    if(drcgs) then
      overlap(1:nst) = ss_full(1:nst)
    else
      overlap(1:nst) = ss(1:nst)
    end if
  end if

  if(present(norm)) then
    if(normalize_) then
      norm = nrm2
    else
      norm = R_REAL(X(mf_nrm2)(mesh, dim, phi))
    end if
  end if

  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(ss_full)

  POP_SUB(X(states_elec_orthogonalization))
  call profiling_out(prof)
end subroutine X(states_elec_orthogonalization)


! ---------------------------------------------------------
FLOAT function X(states_elec_residue)(mesh, dim, hf, ee, ff) result(rr)
  type(mesh_t),      intent(in)  :: mesh
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:)
  FLOAT,             intent(in)  :: ee
  R_TYPE,            intent(in)  :: ff(:,:)

  R_TYPE, allocatable :: res(:,:)
  type(profile_t), save :: prof
  integer :: ip, idim

  PUSH_SUB(X(states_elec_residue))

  call profiling_in(prof, TOSTRING(X(RESIDUE)))

  SAFE_ALLOCATE(res(1:mesh%np, 1:dim))

  do idim = 1, dim
    do ip = 1, mesh%np
      res(ip, idim) = hf(ip, idim) - ee*ff(ip, idim)
    end do
  end do

  call profiling_count_operations(dim*mesh%np*(R_ADD + R_MUL))

  rr = X(mf_nrm2)(mesh, dim, res)
  SAFE_DEALLOCATE_A(res)

  call profiling_out(prof)

  POP_SUB(X(states_elec_residue))

end function X(states_elec_residue)


! ---------------------------------------------------------
!> The routine calculates the expectation value of the momentum 
!! operator
!! <p> = < phi*(ist, k) | -i \nabla | phi(ist, ik) >
!!
! ---------------------------------------------------------
subroutine X(states_elec_calc_momentum)(st, space, der, kpoints, momentum)
  type(states_elec_t),  intent(in)  :: st
  type(space_t),        intent(in)  :: space
  type(derivatives_t),  intent(in)  :: der
  type(kpoints_t),      intent(in)  :: kpoints
  FLOAT,                intent(out) :: momentum(:,:,:)

  integer             :: idim, ist, ik, idir
  CMPLX               :: expect_val_p
  R_TYPE, allocatable :: psi(:, :), grad(:,:,:)
  FLOAT               :: kpoint(1:MAX_DIM)  
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT, allocatable  :: lmomentum(:), gmomentum(:)
  FLOAT, allocatable  :: lmom(:, :, :)
  integer             :: kstart, kend, kn, ndim
#endif

  PUSH_SUB(X(states_elec_calc_momentum))

  SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:st%d%dim, 1:space%dim))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call states_elec_get_state(st, der%mesh, ist, ik, psi)

      do idim = 1, st%d%dim
        call X(derivatives_grad)(der, psi(:, idim), grad(:, idim, 1:space%dim))
      end do

      do idir = 1, space%dim
        ! since the expectation value of the momentum operator is real
        ! for square integrable wfns this integral should be purely imaginary 
        ! for complex wfns but real for real wfns (see case distinction below)
        expect_val_p = X(mf_dotp)(der%mesh, st%d%dim, psi, grad(:, :, idir))

        ! In the case of real wavefunctions we do not include the 
        ! -i prefactor of p = -i \nabla
        if (states_are_real(st)) then
          momentum(idir, ist, ik) = real( expect_val_p )
        else
          momentum(idir, ist, ik) = real( -M_zI*expect_val_p )
        end if
      end do

      ! have to add the momentum vector in the case of periodic systems, 
      ! since psi contains only u_k
      kpoint = M_ZERO
      kpoint(1:space%dim) = kpoints%get_point(st%d%get_kpoint_index(ik))
      do idir = 1, space%periodic_dim
        momentum(idir, ist, ik) = momentum(idir, ist, ik) + kpoint(idir)
      end do

    end do

    ! Exchange momenta in the parallel case.
#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      SAFE_ALLOCATE(lmomentum(1:st%lnst))
      SAFE_ALLOCATE(gmomentum(1:st%nst))

      do idir = 1, der%dim
        lmomentum(1:st%lnst) = momentum(idir, st%st_start:st%st_end, ik)
        call lmpi_gen_allgatherv(st%lnst, lmomentum, tmp, gmomentum, st%mpi_grp)
        momentum(idir, 1:st%nst, ik) = gmomentum(1:st%nst)
      end do

      SAFE_DEALLOCATE_A(lmomentum)
      SAFE_DEALLOCATE_A(gmomentum)
    end if
#endif
  end do

#if defined(HAVE_MPI)
  if(st%d%kpt%parallel) then
    kstart = st%d%kpt%start
    kend = st%d%kpt%end
    kn = st%d%kpt%nlocal
    ndim = ubound(momentum, dim = 1)
    
    ASSERT(.not. st%parallel_in_states)
    
    SAFE_ALLOCATE(lmom(1:ndim, 1:st%nst, 1:kn))
    
    lmom(1:ndim, 1:st%nst, 1:kn) = momentum(1:ndim, 1:st%nst, kstart:kend)
    
    call MPI_Allgatherv(lmom(1, 1, 1), ndim*st%nst*kn, MPI_FLOAT, &
      momentum, st%d%kpt%num(:)*st%nst*ndim, (st%d%kpt%range(1, :) - 1)*st%nst*ndim, MPI_FLOAT, &
      st%d%kpt%mpi_grp%comm, mpi_err)
    
    SAFE_DEALLOCATE_A(lmom)
  end if
#endif  

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad)

  POP_SUB(X(states_elec_calc_momentum))
end subroutine X(states_elec_calc_momentum)


! ---------------------------------------------------------
!> It calculates the expectation value of the angular
!! momentum of the states. If l2 is passed, it also
!! calculates the expectation value of the square of the
!! angular momentum of the state phi.
! ---------------------------------------------------------
subroutine X(states_elec_angular_momentum)(st, gr, ll, l2)
  type(states_elec_t),  intent(in)     :: st
  type(grid_t),    intent(in)     :: gr
  FLOAT,           intent(out)    :: ll(:, :, :) !< (st%nst, st%d%nik, 1 or 3)
  FLOAT, optional, intent(out)    :: l2(:, :)    !< (st%nst, st%d%nik)

  integer :: idim, ist, ik
  R_TYPE, allocatable :: psi(:), lpsi(:, :)

  PUSH_SUB(X(states_elec_angular_momemtum))

  ASSERT(gr%sb%dim /= 1)

  SAFE_ALLOCATE(psi(1:gr%mesh%np_part))

  select case(gr%sb%dim)
  case(3)
    SAFE_ALLOCATE(lpsi(1:gr%mesh%np_part, 1:3))
  case(2)
    SAFE_ALLOCATE(lpsi(1:gr%mesh%np_part, 1:1))
  end select

  ll = M_ZERO
  if(present(l2)) l2 = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim
        call states_elec_get_state(st, gr%mesh, idim, ist, ik, psi)

#if defined(R_TREAL)
        ll = M_ZERO
#else
        call X(physics_op_L)(gr%der, psi, lpsi)

        ll(ist, ik, 1) = ll(ist, ik, 1) + TOFLOAT(X(mf_dotp)(gr%mesh, psi, lpsi(:, 1), reduce = .false.))
        if(gr%sb%dim == 3) then
          ll(ist, ik, 2) = ll(ist, ik, 2) + TOFLOAT(X(mf_dotp)(gr%mesh, psi, lpsi(:, 2), reduce = .false.))
          ll(ist, ik, 3) = ll(ist, ik, 3) + TOFLOAT(X(mf_dotp)(gr%mesh, psi, lpsi(:, 3), reduce = .false.))
        end if
#endif
        if(present(l2)) then
          call X(physics_op_L2)(gr%der, psi(:), lpsi(:, 1))
          l2(ist, ik) = l2(ist, ik) + TOFLOAT(X(mf_dotp)(gr%mesh, psi(:), lpsi(:, 1), reduce = .false.))
        end if
      end do
    end do
  end do

  if(gr%mesh%parallel_in_domains) then
#if !defined(R_TREAL)
    call gr%mesh%allreduce(ll)
#endif
    if(present(l2)) then
      call gr%mesh%allreduce(l2)
    end if
  end if

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(lpsi)
  POP_SUB(X(states_elec_angular_momemtum))
end subroutine X(states_elec_angular_momentum)


! ---------------------------------------------------------
subroutine X(states_elec_matrix)(st1, st2, mesh, aa)
  type(states_elec_t), intent(in)  :: st1, st2
  type(mesh_t),        intent(in)  :: mesh
  R_TYPE,              intent(out) :: aa(:, :, :)

  integer :: ii, jj, dim, ik
  R_TYPE, allocatable :: psi1(:, :), psi2(:, :)
#if defined(HAVE_MPI)
  R_TYPE, allocatable :: phi2(:, :)
  integer :: kk, ll, ist
  integer :: status(MPI_STATUS_SIZE)
  integer :: request
#endif

  PUSH_SUB(X(states_elec_matrix))

  dim = st1%d%dim

  aa (:,:,:) = R_TOTYPE(M_ZERO)

  SAFE_ALLOCATE(psi1(1:mesh%np, 1:st1%d%dim))
  SAFE_ALLOCATE(psi2(1:mesh%np, 1:st1%d%dim))

  aa(:, :, :) = M_ZERO

  do ik = st1%d%kpt%start, st1%d%kpt%end

    if(st1%parallel_in_states) then

#if defined(HAVE_MPI)
      call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
      ! Each process sends the states in st2 to the rest of the processes.
      do ist = st1%st_start, st1%st_end
        call states_elec_get_state(st2, mesh, ist, ik, psi2)
        do jj = 0, st1%mpi_grp%size - 1
          if(st1%mpi_grp%rank /= jj) then
            call MPI_Isend(psi2(1, 1), st1%d%dim*mesh%np, R_MPITYPE, jj, ist, st1%mpi_grp%comm, request, mpi_err)
          end if
        end do
      end do

      ! Processes are received, and then the matrix elements are calculated.
      SAFE_ALLOCATE(phi2(1:mesh%np, 1:st1%d%dim))
      do jj = 1, st2%nst

        ll = st1%node(jj)

        if(ll /= st1%mpi_grp%rank) then
          call MPI_Irecv(phi2(1, 1), st1%d%dim*mesh%np, R_MPITYPE, ll, jj, st1%mpi_grp%comm, request, mpi_err)
          call MPI_Wait(request, status, mpi_err)
        else
          call states_elec_get_state(st2, mesh, jj, ik, phi2)
        end if

        do ist = st1%st_start, st1%st_end
          call states_elec_get_state(st1, mesh, ist, ik, psi1)
          aa(ist, jj, ik) = X(mf_dotp)(mesh, dim, psi1, phi2, reduce = .false.)
        end do

      end do
      SAFE_DEALLOCATE_A(phi2)

      if(mesh%parallel_in_domains) call mesh%allreduce(aa(:,:,ik))

      ! Each process holds some lines of the matrix. So it is broadcasted (All processes
      ! should get the whole matrix)
      call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
      do ii = 1, st1%nst
        kk = st1%node(ii)
        do jj = 1, st2%nst
          call MPI_Bcast(aa(ii, jj, ik), 1, R_MPITYPE, kk, st1%mpi_grp%comm, mpi_err)
        end do
      end do
#endif

    else

      do ii = st1%st_start, st1%st_end

        call states_elec_get_state(st1, mesh, ii, ik, psi1)

        do jj = st2%st_start, st2%st_end

          call states_elec_get_state(st2, mesh, jj, ik, psi2)

          aa(ii, jj, ik) = X(mf_dotp)(mesh, dim, psi1, psi2, reduce = .false.)

        end do
      end do
   
      if(mesh%parallel_in_domains) call mesh%allreduce(aa(:, :, ik))

    end if

  end do

  if(st1%d%kpt%parallel) then
    call comm_allreduce(st1%d%kpt%mpi_grp, aa)
  end if


  SAFE_DEALLOCATE_A(psi1)
  SAFE_DEALLOCATE_A(psi2)    

  POP_SUB(X(states_elec_matrix))
end subroutine X(states_elec_matrix)

! -----------------------------------------------------------

subroutine X(states_elec_calc_orth_test)(st, namespace, mesh, kpoints)
  type(states_elec_t),    intent(inout) :: st
  type(namespace_t),      intent(in)    :: namespace
  type(mesh_t),           intent(in)    :: mesh
  type(kpoints_t),        intent(in)    :: kpoints
  
  PUSH_SUB(X(states_elec_calc_orth_test))

  call states_elec_allocate_wfns(st, mesh, wfs_type = R_TYPE_VAL)

  call states_elec_generate_random(st, mesh, kpoints)

  message(1) = 'Info: Orthogonalizing random wavefunctions.'
  message(2) = ''
  call messages_info(2)

  if(st%d%pack_states) call st%pack()

  call X(states_elec_orthogonalization_full)(st, namespace, mesh, 1)

  if(st%d%pack_states) call st%unpack()

  call print_results()
  
  call states_elec_deallocate_wfns(st)

  POP_SUB(X(states_elec_calc_orth_test))

contains

  subroutine print_results()

    integer :: ist, jst
    R_TYPE :: dd
    R_TYPE, allocatable :: psi1(:, :), psi2(:, :)
    R_TYPE, allocatable :: spsi1(:, :), spsi2(:, :)
#ifdef HAVE_MPI
    integer :: req(1:4), nreq
#endif

    PUSH_SUB(X(states_elec_calc_orth_test).print_results)

    SAFE_ALLOCATE(psi1(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))
#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      SAFE_ALLOCATE(spsi1(1:mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(spsi2(1:mesh%np, 1:st%d%dim))
    end if
#endif

    message(1) = 'Residuals:'
    call messages_info(1)
    
    do ist = 1, st%nst
      do jst = ist, st%nst
        if(.not. st%parallel_in_states) then
          call states_elec_get_state(st, mesh, ist, 1, psi1)
          call states_elec_get_state(st, mesh, jst, 1, psi2)
        end if

#ifdef HAVE_MPI
        if(st%parallel_in_states) then

          ! we bring the two vectors to node 0 and calculate the dot
          ! product, this is very simple and very slow, we only do it for testing
          
          nreq = 0

          ! post the receptions
          if(st%mpi_grp%rank == 0) then
            call MPI_Irecv(psi1(1, 1), mesh%np*st%d%dim, R_MPITYPE, st%node(ist), ist, &
              st%mpi_grp%comm, req(nreq + 1), mpi_err)
            call MPI_Irecv(psi2(1, 1), mesh%np*st%d%dim, R_MPITYPE, st%node(jst), jst, &
              st%mpi_grp%comm, req(nreq + 2), mpi_err)
            nreq = nreq + 2
          end if

          ! if I have the wave function, I send it (note: a node could be sending to itself, this is by design)
          if(st%node(ist)  == st%mpi_grp%rank) then
            nreq = nreq + 1
            call states_elec_get_state(st, mesh, ist, 1, spsi1)
            call MPI_Isend(spsi1(1, 1), mesh%np*st%d%dim, R_MPITYPE, 0, ist, st%mpi_grp%comm, req(nreq), mpi_err)
          end if
          
          if(st%node(jst) == st%mpi_grp%rank) then
            nreq = nreq + 1
            call states_elec_get_state(st, mesh, jst, 1, spsi2)
            call MPI_Isend(spsi2(1, 1), mesh%np*st%d%dim, R_MPITYPE, 0, jst, st%mpi_grp%comm, req(nreq), mpi_err)
          end if

          if(nreq > 0) call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, mpi_err)

          if(st%mpi_grp%rank /= 0) cycle

        end if
#endif
        dd = X(mf_dotp)(mesh, st%d%dim, psi1, psi2)
        write (message(1), '(2i7, e16.6)') ist, jst, abs(dd)
        call messages_info(1)

      end do
    end do
    
    message(1) = ''
    call messages_info(1)

    SAFE_DEALLOCATE_A(psi1)
    SAFE_DEALLOCATE_A(psi2)
    SAFE_DEALLOCATE_A(spsi1)
    SAFE_DEALLOCATE_A(spsi2)

    POP_SUB(X(states_elec_calc_orth_test).print_results)
  end subroutine print_results

end subroutine X(states_elec_calc_orth_test)

! ---------------------------------------------------------

subroutine X(states_elec_rotate)(st, namespace, mesh, uu, ik)
  type(states_elec_t), intent(inout) :: st
  type(namespace_t),   intent(in)    :: namespace
  type(mesh_t),        intent(in)    :: mesh
  R_TYPE,              intent(in)    :: uu(:, :)
  integer,             intent(in)    :: ik
  
  integer       :: block_size, sp, idim, size, ib
  R_TYPE, allocatable :: psinew(:, :, :), psicopy(:, :, :)
  type(accel_mem_t) :: psinew_buffer, psicopy_buffer, uu_buffer
  type(profile_t), save :: prof

  PUSH_SUB(X(states_elec_rotate))

  call profiling_in(prof, TOSTRING(X(STATES_ROTATE)))

  ASSERT(R_TYPE_VAL == st%get_type())
  
  if(.not. st%are_packed() .or. .not. accel_is_enabled()) then
    
#ifdef R_TREAL  
    block_size = max(40, hardware%l2%size/(2*8*st%nst))
#else
    block_size = max(20, hardware%l2%size/(2*16*st%nst))
#endif

    SAFE_ALLOCATE(psinew(1:st%nst, 1:st%d%dim, 1:block_size))
    SAFE_ALLOCATE(psicopy(1:st%nst, 1:st%d%dim, 1:block_size))

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)
      
      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy)
      end do

      if(st%parallel_in_states) call states_elec_parallel_gather(st, (/st%d%dim, size/), psicopy)
      
      do idim = 1, st%d%dim
        
        call blas_gemm(transa = 'c', transb = 'n',        &
          m = st%nst, n = size, k = st%nst,               &
          alpha = R_TOTYPE(M_ONE),                        &
          a = uu(1, 1), lda = ubound(uu, dim = 1),        &
          b = psicopy(1, idim, 1), ldb = st%nst*st%d%dim, &
          beta = R_TOTYPE(M_ZERO),                        & 
          c = psinew(1, idim, 1), ldc = st%nst*st%d%dim)
        
      end do
      
      do ib = st%group%block_start, st%group%block_end
        call batch_set_points(st%group%psib(ib, ik), sp, sp + size - 1, psinew)
      end do

    end do

    call profiling_count_operations((R_ADD + R_MUL)*st%nst*st%d%dim*(st%nst - CNST(1.0))*mesh%np)

    SAFE_DEALLOCATE_A(psinew)
    SAFE_DEALLOCATE_A(psicopy)

  else

    if(st%d%dim > 1) call messages_not_implemented('Opencl states_elec_rotate for spinors', namespace=namespace)

    block_size = batch_points_block_size()

    call accel_create_buffer(uu_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, product(ubound(uu)))
    call accel_write_buffer(uu_buffer, product(ubound(uu)), uu)

    call accel_create_buffer(psicopy_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*st%d%dim*block_size)
    call accel_create_buffer(psinew_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*st%d%dim*block_size)
    if(st%parallel_in_states) then
      SAFE_ALLOCATE(psicopy(1:st%nst, 1:st%d%dim, 1:block_size))
    end if

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)
      
      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psicopy_buffer, st%nst)
      end do

     if(st%parallel_in_states) then
        call accel_read_buffer(psicopy_buffer, st%nst*st%d%dim*block_size, psicopy)
        call states_elec_parallel_gather(st, (/st%d%dim, size/), psicopy)
        call accel_write_buffer(psicopy_buffer, st%nst*st%d%dim*block_size, psicopy)
      end if

      call X(accel_gemm)(transA = CUBLAS_OP_C, transB = CUBLAS_OP_N, &
        M = int(st%nst, 8), N = int(size, 8), K = int(st%nst, 8), alpha = R_TOTYPE(M_ONE), &
        A = uu_buffer, offA = 0_8, lda = int(ubound(uu, dim = 1), 8), &
        B = psicopy_buffer, offB = 0_8, ldb = int(st%nst, 8), beta = R_TOTYPE(M_ZERO), &
        C = psinew_buffer, offC = 0_8, ldc = int(st%nst, 8))
      
      call accel_finish()

      do ib = st%group%block_start, st%group%block_end
        call batch_set_points(st%group%psib(ib, ik), sp, sp + size - 1, psinew_buffer, st%nst)
      end do
    end do

    call profiling_count_operations((R_ADD + R_MUL)*st%nst*(st%nst - CNST(1.0))*mesh%np)
   
    call accel_release_buffer(uu_buffer)
    call accel_release_buffer(psicopy_buffer)
    call accel_release_buffer(psinew_buffer)
    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_A(psicopy)
    end if

  end if

  call profiling_out(prof)
  POP_SUB(X(states_elec_rotate))
end subroutine X(states_elec_rotate)

! ---------------------------------------------------------

subroutine X(states_elec_calc_overlap)(st, mesh, ik, overlap)
  type(states_elec_t), intent(inout) :: st
  type(mesh_t),        intent(in)    :: mesh
  integer,             intent(in)    :: ik
  R_TYPE,              intent(out)   :: overlap(:, :)

  integer :: ip, ib, jb, block_size, sp, size
#ifndef R_TREAL
  integer :: ist, jst
#endif
  type(profile_t), save :: prof
  FLOAT :: vol
  R_TYPE, allocatable :: psi(:, :, :)
  type(accel_mem_t) :: psi_buffer, overlap_buffer

  PUSH_SUB(X(states_elec_calc_overlap))

  call profiling_in(prof, TOSTRING(X(STATES_OVERLAP)))

  if(.not. st%are_packed() .or. .not. accel_is_enabled() .or. &
     (st%parallel_in_states .and. .not. accel_is_enabled())) then

#ifdef R_TREAL  
    block_size = max(80, hardware%l2%size/(8*st%nst))
#else
    block_size = max(40, hardware%l2%size/(16*st%nst))
#endif

    SAFE_ALLOCATE(psi(1:st%nst, 1:st%d%dim, 1:block_size))

    overlap(1:st%nst, 1:st%nst) = R_TOTYPE(M_ZERO)

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psi)
      end do

      if(st%parallel_in_states) call states_elec_parallel_gather(st, (/st%d%dim, size/), psi)
      
      if(mesh%use_curvilinear) then
        do ip = sp, sp + size - 1
          vol = sqrt(mesh%vol_pp(ip))
          psi(1:st%nst, 1:st%d%dim, ip) = psi(1:st%nst, 1:st%d%dim, ip)*vol
        end do
      end if

      call blas_herk(uplo = 'u', trans = 'n',              &
        n = st%nst, k = size*st%d%dim,                     &
        alpha = mesh%volume_element,                       &
        a = psi(1, 1, 1), lda = ubound(psi, dim = 1),      &
        beta = CNST(1.0),                                  & 
        c = overlap(1, 1), ldc = ubound(overlap, dim = 1))

    end do

#ifndef R_TREAL
    do jst = 1, st%nst
      do ist = 1, jst
        overlap(ist, jst) = conjg(overlap(ist, jst))
      end do
    end do
#endif

    call profiling_count_operations((R_ADD + R_MUL)*CNST(0.5)*st%nst*st%d%dim*(st%nst - CNST(1.0))*mesh%np)

    if(mesh%parallel_in_domains) call mesh%allreduce(overlap, dim = (/st%nst, st%nst/))

    SAFE_DEALLOCATE_A(psi)


  else if(accel_is_enabled()) then

    ASSERT(ubound(overlap, dim = 1) == st%nst)
    
    call accel_create_buffer(overlap_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*st%nst)
    call accel_set_buffer_to_zero(overlap_buffer, R_TYPE_VAL, st%nst*st%nst)

    ! we need to use a temporary array

    block_size = batch_points_block_size()

    call accel_create_buffer(psi_buffer, ACCEL_MEM_READ_WRITE, R_TYPE_VAL, st%nst*st%d%dim*block_size)
    if(st%parallel_in_states) then
      SAFE_ALLOCATE(psi(1:st%nst, 1:st%d%dim, 1:block_size))
    end if

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      do ib = st%group%block_start, st%group%block_end
        ASSERT(R_TYPE_VAL == st%group%psib(ib, ik)%type())
        call batch_get_points(st%group%psib(ib, ik), sp, sp + size - 1, psi_buffer, st%nst)
      end do

      if(st%parallel_in_states) then
        call accel_read_buffer(psi_buffer, st%nst*st%d%dim*block_size, psi)
        call states_elec_parallel_gather(st, (/st%d%dim, size/), psi)
        call accel_write_buffer(psi_buffer, st%nst*st%d%dim*block_size, psi)
      end if

      call X(accel_herk)(uplo = ACCEL_BLAS_UPPER, trans = ACCEL_BLAS_N, &
        n = int(st%nst, 8), k = int(size*st%d%dim, 8), &
        alpha = mesh%volume_element, &
        A = psi_buffer, offa = 0_8, lda = int(st%nst, 8), &
        beta = M_ONE, &
        C = overlap_buffer, offc = 0_8, ldc = int(st%nst, 8))
      call accel_finish()
    end do

    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_A(psi)
    end if

    call accel_finish()

    call accel_release_buffer(psi_buffer)

    call profiling_count_operations((R_ADD + R_MUL)*CNST(0.5)*st%nst*st%d%dim*(st%nst - CNST(1.0))*mesh%np)

    call accel_read_buffer(overlap_buffer, st%nst*st%nst, overlap)

    call accel_finish()

#ifndef R_TREAL
    do jst = 1, st%nst
      do ist = 1, jst
        overlap(ist, jst) = conjg(overlap(ist, jst))
      end do
    end do
#endif

    if(mesh%parallel_in_domains) call mesh%allreduce(overlap, dim = (/st%nst, st%nst/))

    call accel_release_buffer(overlap_buffer)

  else

    overlap(1:st%nst, 1:st%nst) = R_TOTYPE(M_ZERO)

    do ib = st%group%block_start, st%group%block_end
      do jb = ib, st%group%block_end
        if(ib == jb) then
          call X(mesh_batch_dotp_self)(mesh, st%group%psib(ib, ik), overlap, reduce = .false.)
        else
          call X(mesh_batch_dotp_matrix)(mesh, st%group%psib(ib, ik), st%group%psib(jb, ik), overlap, reduce = .false.)
        end if
      end do
    end do

    if(mesh%parallel_in_domains) call mesh%allreduce(overlap, dim = (/st%nst, st%nst/))
  end if

  ! Debug output
  if(debug%info .and. mpi_grp_is_root(mpi_world)) then
    do ib = 1, st%nst
      do jb = 1, st%nst
#ifndef R_TREAL
        write(12, '(e12.6,a,e12.6,a)', advance = 'no') real(overlap(ib, jb)), ' ',  aimag(overlap(ib, jb)), ' '
#else
        write(12, '(e12.6,a)', advance = 'no') overlap(ib, jb), ' '
#endif
      end do
      write(12, *) ' ' 
    end do
  end if

  call profiling_out(prof)

  POP_SUB(X(states_elec_calc_overlap))
end subroutine X(states_elec_calc_overlap)

!> This routine computes the projection between two set of states
subroutine X(states_elec_calc_projections)(st, gs_st, namespace, mesh, ik, proj, gs_nst)
  type(states_elec_t),    intent(in)    :: st
  type(states_elec_t),    intent(in)    :: gs_st
  type(namespace_t),      intent(in)    :: namespace
  type(mesh_t),           intent(in)    :: mesh
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(out)   :: proj(:, :)
  integer, optional,      intent(in)    :: gs_nst

  integer       :: ib, ip
  R_TYPE, allocatable :: psi(:, :, :), gspsi(:, :, :)
  integer :: sp, size, block_size
  type(profile_t), save :: prof
  integer :: gs_nst_

  PUSH_SUB(X(states_elec_calc_projections))
  call profiling_in(prof, TOSTRING(X(STATES_PROJECTIONS)))

  if(st%are_packed() .and. accel_is_enabled()) then
   message(1) = "states_elec_calc_projections is not implemented with packed states or accel."
   call messages_fatal(1, namespace=namespace)
  else

#ifdef R_TREAL  
    block_size = max(40, hardware%l2%size/(2*8*st%nst))
#else
    block_size = max(20, hardware%l2%size/(2*16*st%nst))
#endif

   gs_nst_ = gs_st%nst
   if(present(gs_nst)) gs_nst_ = gs_nst

    proj(1:gs_nst_, 1:st%nst) = CNST(0.0)
    
    SAFE_ALLOCATE(psi(1:st%nst, 1:st%d%dim, 1:block_size))
    SAFE_ALLOCATE(gspsi(1:max(gs_nst_, st%nst), 1:gs_st%d%dim, 1:block_size))

    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)
      
      do ib = st%group%block_start, st%group%block_end
        call batch_get_points(st%group%psib(ib, ik),  sp, sp + size - 1, psi)
        call batch_get_points(gs_st%group%psib(ib,ik), sp, sp + size - 1, gspsi)
      end do

      if(st%parallel_in_states) then
        call states_elec_parallel_gather(st, (/st%d%dim, size/), psi)
        call states_elec_parallel_gather(gs_st, (/st%d%dim, size/), gspsi)
      end if
      
      if(mesh%use_curvilinear) then
        do ip = 1, size
          psi(1:st%nst, 1:st%d%dim, ip) = psi(1:st%nst, 1:st%d%dim, ip)*mesh%vol_pp(sp + ip - 1)
          gspsi(1:gs_nst_, 1:st%d%dim, ip) = gspsi(1:gs_nst_, 1:st%d%dim, ip)*mesh%vol_pp(sp + ip - 1)
        end do
      end if

      call blas_gemm(transa = 'n', transb = 'c',        &
        m = gs_nst_, n = st%nst, k = size*st%d%dim,      &
        alpha = R_TOTYPE(mesh%volume_element),      &
        a = gspsi(1, 1, 1), lda = ubound(gspsi, dim = 1),   &
        b = psi(1, 1, 1), ldb = ubound(psi, dim = 1), &
        beta = R_TOTYPE(CNST(1.0)),                     & 
        c = proj(1, 1), ldc = ubound(proj, dim = 1))
    end do

  end if
  
  call profiling_count_operations((R_ADD + R_MUL)*gs_nst_*st%d%dim*(st%nst - CNST(1.0))*mesh%np)
  
  if(mesh%parallel_in_domains) call mesh%allreduce(proj, dim = (/gs_nst_, st%nst/))
  
  call profiling_out(prof)
  POP_SUB(X(states_elec_calc_projections))

end subroutine X(states_elec_calc_projections)

! ---------------------------------------------------------
subroutine X(states_elec_me_one_body)(st, namespace, gr, nspin, vhxc, nint, iindex, jindex, oneint)
  type(states_elec_t), intent(inout) :: st
  type(namespace_t),   intent(in)    :: namespace
  type(grid_t),        intent(in)    :: gr
  integer,             intent(in)    :: nspin
  FLOAT,               intent(in)    :: vhxc(1:gr%mesh%np, nspin)
  integer,             intent(in)    :: nint
  integer,             intent(out)   :: iindex(1:nint)
  integer,             intent(out)   :: jindex(1:nint)
  R_TYPE,              intent(out)   :: oneint(1:nint)  
  
  integer ist, jst, np, iint
  R_TYPE :: me
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(states_elec_me_one_body))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np_part, 1:st%d%dim))

  if (st%d%ispin == SPINORS) then
    call messages_not_implemented("One-body integrals with spinors.", namespace=namespace)
  end if

  
  np = gr%mesh%np
  iint = 1

  do ist = 1, st%nst

    call states_elec_get_state(st, gr%mesh, ist, 1, psii)
    
    do jst = 1, st%nst
      if(jst > ist) cycle
      
      call states_elec_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:np, 1) = R_CONJ(psii(1:np, 1))*vhxc(1:np, 1)*psij(1:np, 1)

      me = - X(mf_integrate)(gr%mesh, psij(:, 1))

      if(ist==jst) me = me + st%eigenval(ist,1)

      iindex(iint) = ist
      jindex(iint) = jst
      oneint(iint) = me
      iint = iint + 1

    end do
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)

  POP_SUB(X(states_elec_me_one_body))
end subroutine X(states_elec_me_one_body)


! ---------------------------------------------------------
subroutine X(states_elec_me_two_body) (st, namespace, space, gr, kpoints, psolver, st_min, st_max, iindex, &
                                         jindex, kindex, lindex, twoint, phase, singularity, exc_k)
  type(states_elec_t), target,   intent(inout) :: st
  type(namespace_t),             intent(in)    :: namespace
  type(space_t),                 intent(in)    :: space
  type(grid_t),                  intent(in)    :: gr
  type(kpoints_t),               intent(in)    :: kpoints
  type(poisson_t),               intent(inout) :: psolver
  integer,                       intent(in)    :: st_min, st_max
  integer,                       intent(out)   :: iindex(:,:)
  integer,                       intent(out)   :: jindex(:,:)
  integer,                       intent(out)   :: kindex(:,:)
  integer,                       intent(out)   :: lindex(:,:)
  R_TYPE,                        intent(out)   :: twoint(:)  !
  CMPLX,               optional, intent(in)    :: phase(:,st%d%kpt%start:)
  type(singularity_t), optional, intent(in)    :: singularity
  logical,             optional, intent(in)    :: exc_k

  integer :: ist, jst, kst, lst, ijst, klst, ikpt, jkpt, kkpt, lkpt
  integer :: ist_global, jst_global, kst_global, lst_global, nst, nst_tot
  integer :: iint, ikpoint, jkpoint, ip, ibind, npath
  R_TYPE  :: me
  R_TYPE, allocatable :: nn(:), vv(:), two_body_int(:), tmp(:)
  R_TYPE, pointer :: psii(:), psij(:), psil(:)
  R_TYPE, allocatable :: psik(:, :)
  FLOAT :: qq(1:MAX_DIM)
  logical :: exc_k_
  class(wfs_elec_t), pointer :: wfs
  type(fourier_space_op_t) :: coulb

  PUSH_SUB(X(states_elec_me_two_body))

  SAFE_ALLOCATE(nn(1:gr%mesh%np))
  SAFE_ALLOCATE(vv(1:gr%mesh%np))
  SAFE_ALLOCATE(psik(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(tmp(1:gr%mesh%np))
  SAFE_ALLOCATE(two_body_int(1:gr%mesh%np))

  if (st%d%ispin == SPINORS) then
    call messages_not_implemented("Two-body integrals with spinors.", namespace=namespace)
  end if

  ASSERT(present(phase) .eqv. present(singularity))
#ifdef R_TCOMPLEX
  ASSERT(present(phase))
#endif

  npath = kpoints_nkpt_in_path(kpoints)

  if(st%are_packed()) call st%unpack()

  ijst = 0
  iint = 1

  nst_tot = (st_max-st_min+1)*st%d%nik
  nst = (st_max-st_min+1)

  exc_k_ = .false.
  if(present(exc_k)) exc_k_ = exc_k

  if(present(singularity)) then
    qq = M_ZERO
    call poisson_build_kernel(psolver, namespace, space, coulb, qq, M_ZERO)
  end if

  do ist_global = 1, nst_tot
    ist = mod(ist_global - 1, nst) + 1
    ikpt = (ist_global - ist) / nst + 1
    ikpoint = st%d%get_kpoint_index(ikpt)

    wfs => st%group%psib(st%group%iblock(ist+st_min-1, ikpt), ikpt)
    ASSERT(wfs%status() /= BATCH_DEVICE_PACKED)
    ibind = wfs%inv_index((/ist+st_min-1, 1/))
    if(wfs%status() == BATCH_NOT_PACKED) then
      psii => wfs%X(ff_linear)(:, ibind)
    else if(wfs%status() == BATCH_PACKED) then
      psii => wfs%X(ff_pack)(ibind, :)
    end if

    do jst_global = 1, nst_tot
      jst = mod(jst_global - 1, nst) + 1
      jkpt = (jst_global - jst) / nst + 1
      jkpoint = st%d%get_kpoint_index(jkpt)

      if(exc_k_ .and. ist /= jst) cycle

      if(present(singularity)) then
        qq(1:gr%der%dim) = kpoints%get_point(ikpoint, absolute_coordinates=.false.) &
                         - kpoints%get_point(jkpoint, absolute_coordinates=.false.)
        ! In case of k-points, the poisson solver must contains k-q 
        ! in the Coulomb potential, and must be changed for each q point
        call poisson_build_kernel(psolver, namespace, space, coulb, qq, M_ZERO, &
                  -(kpoints%full%npoints-npath)*kpoints%latt%rcell_volume*(singularity%Fk(jkpoint)-singularity%FF))
      end if

#ifndef R_TCOMPLEX
      if(jst_global > ist_global) cycle
#endif
      ijst=ijst+1
      
      wfs => st%group%psib(st%group%iblock(jst+st_min-1, jkpt), jkpt)
      ibind = wfs%inv_index((/jst+st_min-1, 1/))
      if(wfs%status() == BATCH_NOT_PACKED) then
        psij => wfs%X(ff_linear)(:, ibind)
      else if(wfs%status() == BATCH_PACKED) then
        psij => wfs%X(ff_pack)(ibind, :)
      end if

      nn(1:gr%mesh%np) = R_CONJ(psii(1:gr%mesh%np))*psij(1:gr%mesh%np)
      if(present(singularity)) then
        call X(poisson_solve)(psolver, vv, nn, all_nodes=.false., kernel=coulb)
      else
        call X(poisson_solve)(psolver, vv, nn, all_nodes=.false.)
      end if

      !We now put back the phase that we treated analytically using the Poisson solver
#ifdef R_TCOMPLEX
      do ip = 1, gr%mesh%np
        vv(ip) = vv(ip) * exp(M_zI*sum(qq(1:gr%der%dim)*gr%mesh%x(ip, 1:gr%der%dim)))
      end do
#endif

      klst=0
      do kst_global = 1, nst_tot
        kst = mod(kst_global - 1, nst) + 1
        kkpt = (kst_global - kst) / nst + 1

        if(exc_k_ .and. kkpt /= jkpt) cycle

        call states_elec_get_state(st, gr%mesh, kst+st_min-1, kkpt, psik)
#ifdef R_TCOMPLEX
        if(present(phase)) then
           call states_elec_set_phase(st%d, psik, phase(1:gr%mesh%np, kkpt), gr%mesh%np, .false.)
        end if
#endif

        tmp(1:gr%mesh%np) = vv(1:gr%mesh%np)*R_CONJ(psik(1:gr%mesh%np, 1))

        do lst_global = 1, nst_tot
          lst = mod(lst_global - 1, nst) + 1
          lkpt = (lst_global - lst)/nst + 1

#ifndef R_TCOMPLEX
          if(lst_global > kst_global) cycle
          klst=klst+1
          if(klst > ijst) cycle
#endif

          if(exc_k_ .and. kst /= lst) cycle
          if(exc_k_ .and. lkpt /= ikpt) cycle
          wfs => st%group%psib(st%group%iblock(lst+st_min-1, lkpt), lkpt)
          ibind = wfs%inv_index((/lst+st_min-1, 1/))
          if(wfs%status() == BATCH_NOT_PACKED) then
            psil => wfs%X(ff_linear)(:, ibind)
          else if(wfs%status() == BATCH_PACKED) then
            psil => wfs%X(ff_pack)(ibind, :)
          end if

          if(present(phase)) then
#ifdef R_TCOMPLEX
            !$omp parallel do
            do ip = 1, gr%mesh%np
              two_body_int(ip) = tmp(ip)*psil(ip)*phase(ip, lkpt)
            end do
            !$omp end parallel do
#endif 
          else
            !$omp parallel do
            do ip = 1, gr%mesh%np
              two_body_int(ip) = tmp(ip)*psil(ip)
            end do
            !$omp end parallel do
          end if

          me = X(mf_integrate)(gr%mesh, two_body_int(:), reduce = .false.)

          iindex(1,iint) =  ist+st_min-1
          iindex(2,iint) =  ikpt
          jindex(1,iint) =  jst+st_min-1
          jindex(2,iint) =  jkpt
          kindex(1,iint) =  kst+st_min-1
          kindex(2,iint) =  kkpt
          lindex(1,iint) =  lst+st_min-1
          lindex(2,iint) =  lkpt
          twoint(iint) =  me
          iint = iint + 1

        end do
      end do
    end do
  end do

  if(gr%mesh%parallel_in_domains) then
    call gr%mesh%allreduce(twoint)
  end if

  if(present(singularity)) then
    call fourier_space_op_end(coulb)
  end if

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(two_body_int)
  SAFE_DEALLOCATE_A(psik)

  POP_SUB(X(states_elec_me_two_body))
end subroutine X(states_elec_me_two_body)


!> Perform RRQR on the transpose states stored in the states object
!! and return the pivot vector 
!! This is not an all-purpose routine for RRQR, but only operates on the
!! specific set stored in st
subroutine X(states_elec_rrqr_decomposition)(st, namespace, mesh, nst, root, ik, jpvt)
  type(states_elec_t), intent(in)  :: st
  type(namespace_t),   intent(in)  :: namespace
  type(mesh_t),        intent(in)  :: mesh
  integer,             intent(in)  :: nst
  logical,             intent(in)  :: root !< this is needed for serial
  integer,             intent(in)  :: ik ! perform SCDM with this k-point
  integer,             intent(out) :: jpvt(:)

  integer :: total_np, nref, info, wsize
  R_TYPE, allocatable :: tau(:), work(:)
  R_TYPE :: tmp
#ifndef R_TREAL
  FLOAT, allocatable :: rwork(:)
#endif
  R_TYPE, allocatable ::  state_global(:), temp_state(:,:)
  R_TYPE, allocatable :: KSt(:,:)
  R_TYPE, allocatable :: psi(:, :)
  integer :: ii,ist,  count, lnst
  logical :: do_serial
  integer :: psi_block(2), blacs_info
  integer, allocatable :: ipiv(:)
#ifdef HAVE_SCALAPACK
  integer :: psi_desc(BLACS_DLEN)
#ifndef R_TREAL
  integer :: rwsize
  FLOAT :: tmp2
#endif
#endif
  integer :: sender
  type(profile_t), save :: prof

  PUSH_SUB(X(states_elec_rrqr_decomposition))
  call profiling_in(prof, TOSTRING(X(RRQR)))

  ASSERT(.not. mesh%use_curvilinear)
  ASSERT(nst == st%nst)

  lnst = st%lnst

  ! decide whether we can use ScaLAPACK
  do_serial = .false.
  if(mesh%parallel_in_domains .or. st%parallel_in_states) then
#ifndef HAVE_SCALAPACK
     message(1) = 'The RRQR is performed in serial. Try linking ScaLAPCK'
     call messages_warning(1, namespace=namespace)
     do_serial = .true.
#else
     if(.not.st%scalapack_compatible) then
        message(1) = 'The RRQR is performed in serial. Try setting ScaLAPACKCompatible = yes'
        call messages_warning(1, namespace=namespace)
        do_serial = .true.
     end if
#endif
  else
     do_serial = .true.
  endif

  if(.not.do_serial) then
    
    call states_elec_parallel_blacs_blocksize(st, namespace, mesh, psi_block, total_np)
    
    ! allocate local part of transpose state matrix
    SAFE_ALLOCATE(KSt(1:lnst,1:total_np))
    SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
    
    ! copy states into the transpose matrix
    count = 0
    do ist = st%st_start,st%st_end
      count = count + 1

      call states_elec_get_state(st, mesh, ist, ik, psi)
      
      ! We need to set to zero some extra parts of the array
      if(st%d%dim == 1) then
        psi(mesh%np + 1:psi_block(1), 1:st%d%dim) = M_ZERO
      else
        psi(mesh%np + 1:mesh%np_part, 1:st%d%dim) = M_ZERO
      end if

      KSt(count, 1:total_np) = psi(1:total_np, 1)
    end do

    SAFE_DEALLOCATE_A(psi)
     
    ! DISTRIBUTE THE MATRIX ON THE PROCESS GRID
    ! Initialize the descriptor array for the main matrices (ScaLAPACK)
#ifdef HAVE_SCALAPACK
    call descinit(psi_desc(1), nst, total_np, psi_block(2), psi_block(1), 0, 0, &
      st%dom_st_proc_grid%context, lnst, blacs_info)
#endif
     
    if(blacs_info /= 0) then
       write(message(1),'(a,i6)') 'descinit failed with error code: ', blacs_info
       call messages_fatal(1, namespace=namespace)
    end if
    
    nref = min(nst, total_np)
    SAFE_ALLOCATE(tau(1:nref))
    tau = M_ZERO

    ! calculate the QR decomposition
    SAFE_ALLOCATE(ipiv(1:total_np))
    ipiv(1:total_np) = 0

    ! Note: lapack routine has different number of arguments depending on type
#ifdef HAVE_SCALAPACK
#ifndef R_TREAL
    call pzgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), tmp, -1, tmp2, -1, blacs_info) 
#else 
    call pdgeqpf( nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), tmp, -1, blacs_info)
#endif
#endif
    
    if(blacs_info /= 0) then
      write(message(1),'(a,i6)') 'scalapack geqrf workspace query failed with error code: ', blacs_info
      call messages_fatal(1, namespace=namespace)
    end if
     
    wsize = nint(R_REAL(tmp))
    SAFE_ALLOCATE(work(1:wsize))
#ifdef HAVE_SCALAPACK
#ifndef R_TREAL
    rwsize = max(1,nint(R_REAL(tmp2)))
    SAFE_ALLOCATE(rwork(1:rwsize))
    call pzgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), work(1), wsize, rwork(1), rwsize, blacs_info)
    SAFE_DEALLOCATE_A(rwork)
#else
    call pdgeqpf(nst, total_np, KSt(1,1), 1, 1, psi_desc(1), ipiv(1), tau(1), work(1), wsize,  blacs_info)
#endif
#endif

    if(blacs_info /= 0) then
      write(message(1),'(a,i6)') 'scalapack geqrf call failed with error code: ', blacs_info
      call messages_fatal(1, namespace=namespace)
    end if
    SAFE_DEALLOCATE_A(work)
     
     ! copy the first nst global elements of ipiv into jpvt
     ! bcast is at the end of the routine
!     if(mpi_world%rank==0)  then
!        do ist =1,nst
!           write(123,*) ipiv(ist)
!        end do
!     end if
    jpvt(1:nst) =  ipiv(1:nst)
     
  else
    ! first gather states into one array on the root process
    ! build transpose of KS set on which RRQR is performed
    if(root) then
       SAFE_ALLOCATE(KSt(1:nst,1:mesh%np_global))
    end if
    
    ! gather states in case of domain parallelization
    if (mesh%parallel_in_domains.or.st%parallel_in_states) then
      SAFE_ALLOCATE(temp_state(1:mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(state_global(1:mesh%np_global))
      
      count = 0
      do ii = 1, nst
        !we are copying states like this:  KSt(i,:) = st%psi(:,dim,i,nik)
        state_global(1:mesh%np_global) = M_ZERO
        sender = 0
        if(state_is_local(st,ii)) then
          call states_elec_get_state(st, mesh, ii, ik, temp_state)
          call vec_gather(mesh%vp, 0, temp_state(1:mesh%np, 1), state_global)
          if(mesh%mpi_grp%rank ==0) sender = mpi_world%rank
        end if
        call comm_allreduce(mpi_world,sender)
#ifdef HAVE_MPI
        call MPI_Bcast(state_global,mesh%np_global , R_MPITYPE, sender, mpi_world%comm, mpi_err)
#endif
        ! keep full Kohn-Sham matrix only on root
        if (root)  KSt(ii,1:mesh%np_global)  = st%occ(ii, 1)*state_global(1:mesh%np_global)
      end do
      SAFE_DEALLOCATE_A(state_global)
      SAFE_DEALLOCATE_A(temp_state)
    else
      ! serial
      SAFE_ALLOCATE(temp_state(1:mesh%np, st%d%dim))
      do ii = 1, nst
        ! this call is necessary becasue we want to have only np not np_part
        call states_elec_get_state(st, mesh, ii, ik, temp_state)
        KSt(ii,:) = st%occ(ii,1)*temp_state(:,1)
      end do
      SAFE_DEALLOCATE_A(temp_state)
    end if

    ! now perform serial RRQR
    ! dummy call to obtain dimension of work
    ! Note: the lapack routine has different number of arguments depending on type
    if(root) then
      SAFE_ALLOCATE(work(1:1))
      SAFE_ALLOCATE(tau(1:nst))
#ifdef R_TREAL
      call dgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, -1, info)
#else
      SAFE_ALLOCATE(rwork(1:2*mesh%np_global))
      call zgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, -1, rwork, info)
#endif
      if (info /= 0) then
         write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
         call messages_fatal(1, namespace=namespace)
      end if

      wsize = int(work(1))
      SAFE_DEALLOCATE_A(work)
      SAFE_ALLOCATE(work(1:wsize))

      jpvt(:) = 0
      tau(:) = 0.
      ! actual call
#ifdef R_TREAL
         call dgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, wsize, info)
#else
         call zgeqp3(nst, mesh%np_global, kst, nst, jpvt, tau, work, wsize, rwork, info)
#endif
      if (info /= 0)then
         write(message(1),'(A28,I2)') 'Illegal argument in ZGEQP3: ', info
         call messages_fatal(1, namespace=namespace)
      end if
      SAFE_DEALLOCATE_A(work)
    endif

    SAFE_DEALLOCATE_A(temp_state)
    SAFE_DEALLOCATE_A(state_global)
    
   endif

#ifdef HAVE_MPI
    call MPI_Bcast(JPVT,nst, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
#endif

   call profiling_out(prof)
   POP_SUB(X(states_elec_rrqr_decomposition))
end subroutine X(states_elec_rrqr_decomposition)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
