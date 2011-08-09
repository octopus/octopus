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
!! $Id$


! ---------------------------------------------------------
!> Orthonormalizes nst orbitals in mesh (honours state
!! parallelization).
subroutine X(states_orthogonalization_full)(st, mesh, ik)
  type(states_t), target, intent(inout) :: st
  type(mesh_t),           intent(in)    :: mesh
  integer,                intent(in)    :: ik


  R_TYPE, allocatable :: ss(:, :), qq(:, :), psi_tmp(:, :, :)
  type(profile_t), save :: prof
  integer :: idim, ist, jst, kst, nst
  FLOAT   :: nrm2
  logical :: bof
  R_TYPE, pointer :: psi(:, :, :)

#ifdef HAVE_SCALAPACK
!pgi$r opt=0
!This is a pragma for the PGI compiler, forcing optimization -O0 for this subroutine
!With PGI 10.9 and ScaLAPACK, at -O2 and higher optimization levels, the test finite_systems_3d/10-fullerene fails in
!states_orthogonalization_full.par_gs with error message
!glibc detected *** octopus_mpi: malloc(): memory corruption
#endif

  ASSERT(st%d%orth_method /= 0)

  call profiling_in(prof, "GRAM_SCHMIDT_FULL")
  PUSH_SUB(X(states_orthogonalization_full))

  psi => st%X(psi)(:, :, :, ik)
  nst = st%nst

  select case(st%d%orth_method)
  case(ORTH_GS)

    if(st%parallel_in_states) then
      message(1) = 'The selected orthogonalization method cannot work with state-parallelization.'
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(ss(1:nst, 1:nst))

    call X(states_overlap)(st, mesh, ik, ss)

    bof = .false.
    ! calculate the Cholesky decomposition
    call lalg_cholesky(nst, ss, bof = bof)

    if(bof) then
      message(1) = "Orthogonalization failed; probably your eigenvectors are not independent."
      call messages_warning(1)
    end if

    do idim = 1, st%d%dim
      ! multiply by the inverse of ss
      call blas_trsm('R', 'U', 'N', 'N', mesh%np, nst, R_TOTYPE(M_ONE), ss(1, 1), nst, &
        psi(1, idim, 1), ubound(psi, dim = 1)*st%d%dim)
    end do

    call profiling_count_operations(dble(mesh%np)*dble(nst)**2*(R_ADD + R_MUL))

    SAFE_DEALLOCATE_A(ss)

  case(ORTH_PAR_GS)
    call par_gs()

  case(ORTH_QR)
    call qr()

  case(ORTH_MGS)
    call mgs()
  end select

  call profiling_out(prof)
  POP_SUB(X(states_orthogonalization_full))

contains
  
  subroutine par_gs()
#ifdef HAVE_SCALAPACK
    integer             :: info, nbl, nrow, ncol
    integer             :: psi_block(1:2), total_np, psi_desc(BLACS_DLEN), ss_desc(BLACS_DLEN)
#endif

    PUSH_SUB(X(states_orthogonalization_full).par_gs)

! some checks
#ifndef HAVE_MPI
    message(1) = 'The parallel Gram-Schmidt orthogonalizer can only be used in parallel.'
    call messages_fatal(1)
#else
#ifndef HAVE_SCALAPACK
    message(1) = 'The parallel Gram-Schmidt orthogonalizer requires ScaLAPACK.'
    call messages_fatal(1)
#endif
    if(st%dom_st_mpi_grp%size == 1) then
      message(1) = 'The parallel Gram-Schmidt orthogonalizer is designed to be used with domain or state parallelization.'
      call messages_warning(1)
    end if
#endif


#ifdef HAVE_SCALAPACK
    call states_blacs_blocksize(st, mesh, psi_block, total_np)

    ! We need to set to zero some extra parts of the array
    if(st%d%dim == 1) then
      psi(mesh%np + 1:psi_block(1), 1:st%d%dim, 1:st%lnst) = M_ZERO
    else
      psi(mesh%np + 1:mesh%np_part, 1:st%d%dim, 1:st%lnst) = M_ZERO
    end if

    call descinit(psi_desc(1), total_np, st%nst, psi_block(1), psi_block(2), 0, 0, st%dom_st_proc_grid%context, &
      st%d%dim*ubound(psi, dim = 1), info)

    if(info /= 0) then
      write(message(1),'(a,i6)') "descinit for psi failed in states_orthogonalization_full.par_gs with error ", info
      call messages_warning(1)
    end if

    nbl = min(32, st%nst)
    nrow = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%myrow, 0, st%dom_st_proc_grid%nprow))
    ncol = max(1, numroc(st%nst, nbl, st%dom_st_proc_grid%mycol, 0, st%dom_st_proc_grid%npcol))

    SAFE_ALLOCATE(ss(1:nrow, 1:ncol))

    call descinit(ss_desc(1), st%nst, st%nst, nbl, nbl, 0, 0, st%dom_st_proc_grid%context, ubound(ss, dim = 1), info)

    if(info /= 0) then
      write(message(1),'(a,i6)') "descinit for ss failed in states_orthogonalization_full.par_gs with error ", info
      call messages_warning(1)
    end if

    ss = M_ZERO

    call pblas_herk(uplo = 'U', trans = 'C', n = st%nst, k = total_np, &
      alpha = R_TOTYPE(mesh%vol_pp(1)), a = psi(1, 1, 1), ia = 1, ja = 1, desca = psi_desc(1), &
      beta = R_TOTYPE(M_ZERO), c = ss(1, 1), ic = 1, jc = 1, descc = ss_desc(1))

    ! calculate the Cholesky decomposition
    call scalapack_potrf(uplo = 'U', n = st%nst, a = ss(1, 1), ia = 1, ja = 1, desca = ss_desc(1), info = info)

    if(info /= 0) then
      write(message(1),'(a,i6,a)') &
        "Orthogonalization with potrf failed with error ", info, "; probably your eigenvectors are not independent."
      call messages_warning(1)
    end if

    call pblas_trsm(side = 'R', uplo = 'U', transa = 'N', diag = 'N', m = total_np, n = st%nst, &
      alpha = R_TOTYPE(M_ONE), a = ss(1, 1), ia = 1, ja = 1, desca = ss_desc(1), &
      b = psi(1, 1, 1), ib = 1, jb = 1, descb = psi_desc(1))

    call profiling_count_operations(dble(mesh%np)*dble(nst)**2*(R_ADD + R_MUL))

    SAFE_DEALLOCATE_A(ss)
#endif

    POP_SUB(X(states_orthogonalization_full).par_gs)
  end subroutine par_gs

  ! -----------------------------------------------------------------------------------------------

  subroutine qr()
    integer :: total_np, nref, info, wsize
    R_TYPE, allocatable :: tau(:), work(:)
    R_TYPE :: tmp
#ifdef HAVE_SCALAPACK
    integer :: psi_block(2), psi_desc(BLACS_DLEN), blacs_info
#endif

    PUSH_SUB(X(states_orthogonalization_full).qr)

    ASSERT(.not. mesh%use_curvilinear)

    if(mesh%parallel_in_domains .or. st%parallel_in_states) then

#ifndef HAVE_SCALAPACK
      message(1) = 'The QR orthogonalizer requires ScaLAPACK to work in parallel.'
      call messages_fatal(1)
#endif

#ifdef HAVE_SCALAPACK

      call states_blacs_blocksize(st, mesh, psi_block, total_np)

      ! We need to set to zero some extra parts of the array
      if(st%d%dim == 1) then
        psi(mesh%np + 1:psi_block(1), 1:st%d%dim, 1:st%lnst) = M_ZERO
      else
        psi(mesh%np + 1:mesh%np_part, 1:st%d%dim, 1:st%lnst) = M_ZERO
      end if

      ! DISTRIBUTE THE MATRIX ON THE PROCESS GRID
      ! Initialize the descriptor array for the main matrices (ScaLAPACK)
      call descinit(psi_desc(1), total_np, nst, psi_block(1), psi_block(2), 0, 0, &
        st%dom_st_proc_grid%context, mesh%np_part*st%d%dim, blacs_info)

      nref = min(nst, total_np)
      SAFE_ALLOCATE(tau(1:nref))
      tau = M_ZERO

      ! calculate the QR decomposition
      call scalapack_geqrf(total_np, nst, psi(1, 1, 1), 1, 1, psi_desc(1), tau(1), tmp, -1, blacs_info)
      wsize = nint(R_REAL(tmp))
      SAFE_ALLOCATE(work(1:wsize))
      call scalapack_geqrf(total_np, nst, psi(1, 1, 1), 1, 1, psi_desc(1), tau(1), work(1), wsize, blacs_info)
      SAFE_DEALLOCATE_A(work)

      ! now calculate Q
      call scalapack_orgqr(total_np, nst, nref, psi(1, 1, 1), 1, 1, psi_desc(1), tau(1), tmp, -1, blacs_info)
      wsize = nint(R_REAL(tmp))
      SAFE_ALLOCATE(work(1:wsize))
      call scalapack_orgqr(total_np, nst, nref, psi(1, 1, 1), 1, 1, psi_desc(1), tau(1), work(1), wsize, blacs_info)
      SAFE_DEALLOCATE_A(work)

      if(blacs_info /= 0) then
        write(message(1),'(a,I6)') 'ScaLAPACK execution failed. Failed code: ', blacs_info
        call messages_warning(1)
      end if
#else
      message(1) = 'The QR orthogonalization in parallel requires ScaLAPACK.'
      call messages_fatal(1)
#endif 
    else

      total_np = mesh%np + mesh%np_part*(st%d%dim - 1)
      psi(mesh%np + 1:mesh%np_part, 1:(st%d%dim - 1), 1:st%lnst) = M_ZERO

      nref = min(nst, total_np)
      SAFE_ALLOCATE(tau(1:nref))
      tau = M_ZERO

      ! get the optimal size of the work array
      call lapack_geqrf(total_np, nst, psi(1, 1, 1), mesh%np_part*st%d%dim, tau(1), tmp, -1, info)
      wsize = nint(R_REAL(tmp))

      ! calculate the QR decomposition
      SAFE_ALLOCATE(work(1:wsize))
      call lapack_geqrf(total_np, nst, psi(1, 1, 1), mesh%np_part*st%d%dim, tau(1), work(1), wsize, info)
      SAFE_DEALLOCATE_A(work)

      ! get the optimal size of the work array
      call lapack_orgqr(total_np, nst, nref, psi(1, 1, 1), mesh%np_part*st%d%dim, tau(1), tmp, -1, info)
      wsize = nint(R_REAL(tmp))

      ! now calculate Q
      SAFE_ALLOCATE(work(1:wsize))
      call lapack_orgqr(total_np, nst, nref, psi(1, 1, 1), mesh%np_part*st%d%dim, tau(1), work(1), wsize, info)
      SAFE_DEALLOCATE_A(work)
    end if

    SAFE_DEALLOCATE_A(tau)

    ! we need to scale by the volume element to get the proper normalization
    psi = psi/sqrt(mesh%vol_pp(1))

    POP_SUB(X(states_orthogonalization_full).qr)

  end subroutine qr

  ! ----------------------------------------------------------------------------------

  subroutine mgs()
    integer :: ist, jst, idim
    FLOAT   :: cc
    R_TYPE, allocatable :: aa(:)
    FLOAT,  allocatable :: bb(:)

    PUSH_SUB(X(states_orthogonalization_full).mgs)

    if(st%parallel_in_states) then
      message(1) = 'The selected orthogonalization method cannot work with state-parallelization.'
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(bb(1:nst))

    ! normalize the initial vectors
    do ist = 1, nst
      bb(ist) = X(mf_dotp)(mesh, st%d%dim, psi(:, :, ist), psi(:, :, ist), reduce = .false.)
    end do

    if(mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, bb, dim = nst)

    do ist = 1, nst      
      do idim = 1, st%d%dim
        call lalg_scal(mesh%np, M_ONE/sqrt(bb(ist)), psi(:, idim, ist))
      end do
    end do

    SAFE_DEALLOCATE_A(bb)

    SAFE_ALLOCATE(aa(1:nst))

    do ist = 1, nst
      ! calculate the projections
      do jst = 1, ist - 1
        aa(jst) = X(mf_dotp)(mesh, st%d%dim, psi(:, :, jst), psi(:, :, ist), reduce = .false.)
      end do

      if(mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, aa, dim = ist - 1)

      ! substract the projections
      do jst = 1, ist - 1
        do idim = 1, st%d%dim
          call lalg_axpy(mesh%np, -aa(jst), psi(:, idim, jst), psi(:, idim, ist))
        end do
      end do

      ! renormalize
      cc = X(mf_dotp)(mesh, st%d%dim, psi(:, :, ist), psi(:, :, ist))
      do idim = 1, st%d%dim
        call lalg_scal(mesh%np, M_ONE/sqrt(cc), psi(:, idim, ist))
      end do
    end do

    SAFE_DEALLOCATE_A(aa)
    POP_SUB(X(states_orthogonalization_full).mgs)
  end subroutine mgs

end subroutine X(states_orthogonalization_full)


! ---------------------------------------------------------
!> Orthonormalizes phi to the nst orbitals psi.
!! It also permits doing only the orthogonalization (no normalization).
!! And one can pass an extra optional argument, mask, which:
!!  - on input, if mask(p) = .true., the p-orbital is not used.
!!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
!! If Theta_Fi and beta_ij are present, it performs the generalized orthogonalization
!!   (Theta_Fi - sum_j beta_ij |j><j|Phi> as in De Gironcoli PRB 51, 6774 (1995).
!! This is used in response for metals
subroutine X(states_orthogonalization)(mesh, nst, dim, psi, phi,  &
  normalize, mask, overlap, norm, Theta_fi, beta_ij)
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
  R_TYPE,  optional, intent(in)    :: beta_ij(:)   ! beta_ij(nst)

  logical :: normalize_
  integer :: ist, idim
  FLOAT   :: nrm2
  R_TYPE, allocatable  :: ss(:)
  integer :: block_size, size, sp, ep
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof

  call profiling_in(prof, "GRAM_SCHMIDT")
  PUSH_SUB(X(states_orthogonalization))

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.
  block_size = hardware%X(block_size)

  SAFE_ALLOCATE(ss(1:nst))

  ss = M_ZERO

  if(.not. mesh%use_curvilinear) then

    do idim = 1, dim
      do sp = 1, mesh%np, block_size
        size = min(block_size, mesh%np - sp + 1)
        do ist = 1, nst

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

    do idim = 1, dim
      do sp = 1, mesh%np, block_size
        size = min(block_size, mesh%np - sp + 1)
        ep = sp - 1 + size
        do ist = 1, nst

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
    call profiling_in(reduce_prof, "GRAM_SCHMIDT_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, ss, dim = nst)
    call profiling_out(reduce_prof)
  end if

  if(present(mask)) then
    do ist = 1, nst
      mask(ist) = (abs(ss(ist)) <= M_EPSILON)
    end do
  end if

  if(present(beta_ij))  &
    ss(:) = ss(:) * beta_ij(:)

  do idim = 1, dim
    do sp = 1, mesh%np, block_size
      size = min(block_size, mesh%np - sp + 1)

      if(present(Theta_Fi)) then
        if(Theta_Fi .ne. M_ONE) &
          call blas_scal(size, R_TOTYPE(Theta_Fi), phi(sp, idim), 1)
      end if

      do ist = 1, nst

        if(present(mask)) then
          if(mask(ist)) cycle
        end if

        call blas_axpy(size, -ss(ist), psi(sp, idim, ist), 1, phi(sp, idim), 1)

      end do
    end do
  end do

  call profiling_count_operations((R_ADD + R_MUL) * mesh%np * dim * nst)

  ! the following ifs cannot be given as a single line (without the
  ! then) to avoid a bug in xlf 10.1

  normalize_ = .false.
  if(present(normalize)) then
    normalize_ = normalize
  end if

  if(normalize_) then
    nrm2 = X(mf_nrm2)(mesh, dim, phi)
    do idim = 1, dim
      call lalg_scal(mesh%np, M_ONE / nrm2, phi(:, idim))
    end do
  end if

  if(present(overlap)) then
    overlap(1:nst) = ss(1:nst)
  end if

  if(present(norm)) then
    ASSERT(normalize)
    norm = nrm2
  end if

  SAFE_DEALLOCATE_A(ss)

  POP_SUB(X(states_orthogonalization))
  call profiling_out(prof)
end subroutine X(states_orthogonalization)


! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(mesh, dim, psi)
  type(mesh_t),    intent(in)    :: mesh
  integer,         intent(in)    :: dim
  R_TYPE,          intent(inout) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim, ip

  PUSH_SUB(X(states_normalize_orbital))

  norm = X(mf_nrm2) (mesh, dim, psi)

  forall (idim = 1:dim, ip = 1:mesh%np) psi(ip, idim) = psi(ip, idim)/norm
  
  POP_SUB(X(states_normalize_orbital))
end subroutine X(states_normalize_orbital)


! ---------------------------------------------------------
FLOAT function X(states_residue)(mesh, dim, hf, ee, ff) result(rr)
  type(mesh_t),      intent(in)  :: mesh
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:), ff(:,:)
  FLOAT,             intent(in)  :: ee

  R_TYPE, allocatable :: res(:,:)
  type(profile_t), save :: prof
  integer :: ip, idim

  PUSH_SUB(X(states_residue))

  call profiling_in(prof, "RESIDUE")

  SAFE_ALLOCATE(res(1:mesh%np_part, 1:dim))

  forall (idim = 1:dim, ip = 1:mesh%np) res(ip, idim) = hf(ip, idim) - ee*ff(ip, idim)

  call profiling_count_operations(dim*mesh%np*(R_ADD + R_MUL))

  rr = X(mf_nrm2)(mesh, dim, res)
  SAFE_DEALLOCATE_A(res)

  call profiling_out(prof)

  POP_SUB(X(states_residue))

end function X(states_residue)


! ---------------------------------------------------------
!> The routine calculates the expectation value of the momentum 
!! operator
!! <p> = < phi*(ist, k) | -i \nabla | phi(ist, ik) >
!!
! ---------------------------------------------------------
subroutine X(states_calc_momentum)(st, der, momentum)
  type(states_t),      intent(inout) :: st
  type(derivatives_t), intent(inout) :: der
  FLOAT,               intent(out)   :: momentum(:,:,:)

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

  PUSH_SUB(X(states_calc_momentum))

  SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:st%d%dim, 1:der%mesh%sb%dim))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call states_get_state(st, der%mesh, ist, ik, psi)

      do idim = 1, st%d%dim
        call X(derivatives_grad)(der, psi(:, idim), grad(:, idim, 1:der%mesh%sb%dim))
      end do

      do idir = 1, der%mesh%sb%dim
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
      kpoint(1:der%mesh%sb%dim) = kpoints_get_point(der%mesh%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
      forall(idir = 1:der%mesh%sb%periodic_dim) momentum(idir, ist, ik) = momentum(idir, ist, ik) + kpoint(idir)

    end do

    ! Exchange momenta in the parallel case.
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

    if(st%parallel_in_states) then
      SAFE_ALLOCATE(lmomentum(1:st%lnst))
      SAFE_ALLOCATE(gmomentum(1:st%nst))

      do idir = 1, der%mesh%sb%dim
        lmomentum(1:st%lnst) = momentum(idir, st%st_start:st%st_end, ik)
        call lmpi_gen_allgatherv(st%lnst, lmomentum, tmp, gmomentum, st%mpi_grp)
        momentum(idir, 1:st%nst, ik) = gmomentum(1:st%nst)
      end do

      SAFE_DEALLOCATE_A(lmomentum)
      SAFE_DEALLOCATE_A(gmomentum)
    end if
#endif
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad)

  POP_SUB(X(states_calc_momentum))
end subroutine X(states_calc_momentum)


! ---------------------------------------------------------
!> It calculates the expectation value of the angular
!! momentum of the state phi. If l2 is passed, it also
!! calculates the expectation value of the square of the
!! angular momentum of the state phi.
! ---------------------------------------------------------
subroutine X(states_angular_momentum)(gr, phi, ll, l2)
  type(grid_t), intent(inout)  :: gr
  R_TYPE,       intent(inout)  :: phi(:, :)
  FLOAT,        intent(out)    :: ll(MAX_DIM)
  FLOAT, optional, intent(out) :: l2

  integer :: idim, dim
  R_TYPE, allocatable :: lpsi(:, :)

  PUSH_SUB(X(states_angular_momemtum))

  ASSERT(gr%mesh%sb%dim.ne.1)

  select case(gr%mesh%sb%dim)
  case(3)
    SAFE_ALLOCATE(lpsi(1:gr%mesh%np_part, 1:3))
  case(2)
    SAFE_ALLOCATE(lpsi(1:gr%mesh%np_part, 1:1))
  end select

  dim = size(phi, 2)

  ll = M_ZERO
  if(present(l2)) l2 = M_ZERO

  do idim = 1, dim
#if defined(R_TREAL)
    ll = M_ZERO
#else
    call X(physics_op_L)(gr%der, phi(:, idim), lpsi)
    select case(gr%mesh%sb%dim)
    case(3)
      ll(1) = ll(1) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
      ll(2) = ll(2) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 2))
      ll(3) = ll(3) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 3))
    case(2)
      ll(3) = ll(3) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
    end select
#endif
    if(present(l2)) then
      call X(physics_op_L2)(gr%der, phi(:, idim), lpsi(:, 1))
      l2 = l2 + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
    end if
  end do

  SAFE_DEALLOCATE_A(lpsi)
  POP_SUB(X(states_angular_momemtum))
end subroutine X(states_angular_momentum)


! ---------------------------------------------------------
subroutine X(states_matrix)(mesh, st1, st2, aa)
  type(mesh_t),   intent(in)  :: mesh
  type(states_t), intent(in)  :: st1, st2
  R_TYPE,         intent(out) :: aa(:, :, :)

  integer :: ii, jj, dim, n1, n2, ik
#if defined(HAVE_MPI)
  R_TYPE, allocatable :: phi2(:, :)
  integer :: kk, ll, ist
  integer :: status(MPI_STATUS_SIZE)
  integer :: request
#endif

  PUSH_SUB(X(states_matrix))

  n1 = st1%nst
  n2 = st2%nst

  dim = st1%d%dim

  do ik = 1, st1%d%nik

  if(st1%parallel_in_states) then

#if defined(HAVE_MPI)
    call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
    ! Each process sends the states in st2 to the rest of the processes.
    do ist = st1%st_start, st1%st_end
      do jj = 0, st1%mpi_grp%size - 1
        if(st1%mpi_grp%rank.ne.jj) then
          call MPI_Isend(st2%X(psi)(1, 1, ist, ik), st1%d%dim*mesh%np, R_MPITYPE, &
            jj, ist, st1%mpi_grp%comm, request, mpi_err)
        end if
      end do
    end do

    ! Processes are received, and then the matrix elements are calculated.
    SAFE_ALLOCATE(phi2(1:mesh%np, 1:st1%d%dim))
    do jj = 1, n2
      ll = st1%node(jj)
      if(ll.ne.st1%mpi_grp%rank) then
        call MPI_Irecv(phi2(1, 1), st1%d%dim*mesh%np, R_MPITYPE, ll, jj, st1%mpi_grp%comm, request, mpi_err)
        call MPI_Wait(request, status, mpi_err)
      else
        phi2(:, :) = st2%X(psi)(:, :, jj, ik)
      end if
      do ist = st1%st_start, st1%st_end
        aa(ist, jj, ik) = X(mf_dotp)(mesh, dim, st1%X(psi)(:, :, ist, ik), phi2(:, :))
      end do
    end do
    SAFE_DEALLOCATE_A(phi2)

    ! Each process holds some lines of the matrix. So it is broadcasted (All processes
    ! should get the whole matrix)
    call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
    do ii = 1, n1
      kk = st1%node(ii)
      do jj = 1, n2
        call MPI_Bcast(aa(ii, jj, ik), 1, R_MPITYPE, kk, st1%mpi_grp%comm, mpi_err)
      end do
    end do
#else
    write(message(1), '(a)') 'Internal error at Xstates_matrix'
    call messages_fatal(1)
#endif

  else
    do ii = 1, n1
      do jj = 1, n2
        aa(ii, jj, ik) = X(mf_dotp)(mesh, dim, st1%X(psi)(:, :, ii, ik), st2%X(psi)(:, :, jj, ik))
      end do
    end do
  end if

  end do

  POP_SUB(X(states_matrix))
end subroutine X(states_matrix)

! -----------------------------------------------------------

subroutine X(states_calc_orth_test)(st, mc, mesh)
  type(states_t),    intent(inout) :: st
  type(multicomm_t), intent(in)    :: mc
  type(mesh_t),      intent(in)    :: mesh
  
  PUSH_SUB(X(states_calc_orth_test))

  call states_allocate_wfns(st, mesh, wfs_type = R_TYPE_VAL)

  call states_generate_random(st, mesh)

  message(1) = 'Info: Orthogonalizing random wavefunctions.'
  message(2) = ''
  call messages_info(2)

  call X(states_orthogonalization_full)(st, mesh, 1)

  call print_results()
  
  call states_deallocate_wfns(st)

  POP_SUB(X(states_calc_orth_test))

contains
  subroutine print_results()
    integer :: ist, jst
    FLOAT :: dd
    R_TYPE, allocatable :: psi1(:, :), psi2(:, :)
    R_TYPE, allocatable :: spsi1(:, :), spsi2(:, :)
#ifdef HAVE_MPI
    integer :: req(1:4), nreq
#endif

    PUSH_SUB(X(states_calc_orth_test).print_results)

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
          call states_get_state(st, mesh, ist, 1, psi1)
          call states_get_state(st, mesh, jst, 1, psi2)
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
            INCR(nreq, 2)
          end if

          ! if I have the wave function, I send it (note: a node could be sending to itself, this is by design)
          if(st%node(ist)  == st%mpi_grp%rank) then
            INCR(nreq, 1)
            call states_get_state(st, mesh, ist, 1, spsi1)
            call MPI_Isend(spsi1(1, 1), mesh%np*st%d%dim, R_MPITYPE, 0, ist, st%mpi_grp%comm, req(nreq), mpi_err)
          end if
          
          if(st%node(jst) == st%mpi_grp%rank) then
            INCR(nreq, 1)
            call states_get_state(st, mesh, jst, 1, spsi2)
            call MPI_Isend(spsi2(1, 1), mesh%np*st%d%dim, R_MPITYPE, 0, jst, st%mpi_grp%comm, req(nreq), mpi_err)
          end if

          if(nreq > 0) call MPI_Waitall(nreq, req(1), MPI_STATUSES_IGNORE, mpi_err)

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

    POP_SUB(X(states_calc_orth_test).print_results)
  end subroutine print_results
end subroutine X(states_calc_orth_test)

! ---------------------------------------------------------

subroutine X(states_rotate_in_place)(mesh, st, uu, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(states_t),    intent(inout) :: st
  R_TYPE,            intent(in)    :: uu(:, :)
  integer,           intent(in)    :: ik
  
  type(batch_t) :: psib
  
  PUSH_SUB(X(states_rotate_in_place))
  
  ASSERT(associated(st%X(psi)))

  call batch_init(psib, st%d%dim, 1, st%nst, st%X(psi)(:, :, :, ik))
  call X(mesh_batch_rotate)(mesh, psib, uu)
  call batch_end(psib)

  POP_SUB(X(states_rotate_in_place))
end subroutine X(states_rotate_in_place)

! ---------------------------------------------------------

subroutine X(states_overlap)(st, mesh, ik, overlap)
  type(states_t),    intent(inout)    :: st
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: ik
  R_TYPE,            intent(out)   :: overlap(:, :)
  
  integer       :: ib, jb
  type(batch_t) :: psib
  
  PUSH_SUB(X(states_overlap))
  
  ASSERT(associated(st%X(psi)))

  if(states_are_packed(st)) then

    call batch_init(psib, st%d%dim, 1, st%nst, st%X(psi)(:, :, :, ik))
    call X(mesh_batch_dotp_self)(mesh, psib, overlap)
    call batch_end(psib)
    
  else
    
    do ib = st%block_start, st%block_end
      do jb = ib, st%block_end
        if(ib == jb) then
          call X(mesh_batch_dotp_self)(mesh, st%psib(ib, ik), overlap)
        else
          call X(mesh_batch_dotp_matrix)(mesh, st%psib(ib, ik), st%psib(jb, ik), overlap)
        end if
      end do
    end do
    
  end if

  POP_SUB(states_overlap)
end subroutine X(states_overlap)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
