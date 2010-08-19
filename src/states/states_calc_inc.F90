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
! Orthonormalizes nst orbitals in mesh (honours state
! parallelization).
subroutine X(states_orthogonalization_full)(st, nst, mesh, dim, psi)
  type(states_t),    intent(in)    :: st
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: mesh
  R_TYPE, target,    intent(inout) :: psi(:, :, st%st_start:)

  R_TYPE, allocatable :: ss(:, :), qq(:, :), psi_tmp(:, :, :)
  type(profile_t), save :: prof
  integer :: idim, ist, jst, kst
  FLOAT   :: nrm2

  call profiling_in(prof, "GRAM_SCHMIDT_FULL")
  PUSH_SUB(X(states_orthogonalization_full))

  if(.not. st%parallel_in_states) then
    call X(states_orthogonalization_block)(st, nst, mesh, dim, psi)
  else

    SAFE_ALLOCATE(qq(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(ss(1:nst, 1:nst))
    ss = M_ZERO

    call states_blockt_mul(mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
      psi, psi, ss, symm = .true.)
  
    qq = M_ZERO
    do ist = 1, st%nst

      if(ist > nst) then
        qq(ist, ist) = M_ONE
        cycle
      end if

      ! calculate the norm of the resulting vector, and use it to scale
      ! the coefficients so we get normalized vectors.
      nrm2 = ss(ist, ist)
      do jst = 1, ist - 1
        nrm2 = nrm2 - M_TWO*abs(ss(ist, jst))**2/ss(jst, jst)
        do kst = 1, ist - 1
          nrm2 = nrm2 + ss(jst, kst)*ss(ist, jst)*ss(kst, ist)/(ss(jst, jst)*ss(kst, kst))
        end do
      end do
      nrm2 = M_ONE/sqrt(nrm2)

      ! now generate the matrix with the linear combination of orbitals
      qq(ist, ist) = nrm2
      do jst = 1, ist - 1
        qq(jst, ist) = -ss(jst, ist)/ss(jst, jst)*nrm2
      end do

    end do

    SAFE_ALLOCATE(psi_tmp(1:mesh%np_part, 1:dim, st%st_start:st%st_end))

    do ist = st%st_start, st%st_end
      do idim = 1, dim
        call lalg_copy(mesh%np, psi(:, idim, ist), psi_tmp(:, idim, ist))
      end do
    end do

    call states_block_matr_mul(mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, psi_tmp, qq, psi)

    SAFE_DEALLOCATE_A(psi_tmp)

  end if

  SAFE_DEALLOCATE_A(ss)
  SAFE_DEALLOCATE_A(qq)

  POP_SUB(X(states_orthogonalization_full))
  call profiling_out(prof)
end subroutine X(states_orthogonalization_full)


! ---------------------------------------------------------
subroutine X(states_orthogonalization_block)(st, nst, mesh, dim, psi)
  type(states_t),    intent(in)    :: st
  integer,           intent(in)    :: nst
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: dim
  R_TYPE, target,    intent(inout) :: psi(:, :, :)

  R_TYPE, allocatable :: ss(:, :), work(:), tau(:)
  FLOAT :: tmp
  type(batch_t) :: psib
  logical :: bof
  integer :: idim, nref, wsize, info

  PUSH_SUB(X(states_orthogonalization_block))

  select case(st%d%orth_method)
  case(ORTH_GS)

    SAFE_ALLOCATE(ss(1:nst, 1:nst))
    ss = M_ZERO

    call batch_init(psib, st%d%dim, 1, nst, psi)
    call X(mesh_batch_dotp_self)(mesh, psib, ss)
    call batch_end(psib)

    bof = .false.
    ! calculate the Cholesky decomposition
    call lalg_cholesky(nst, ss, bof = bof)

    if(bof) then
      message(1) = "Warning: Orthogonalization failed; probably your eigenvectors are not independent."
      call write_warning(1)
    end if

    do idim = 1, st%d%dim
      ! multiply by the inverse of ss
      call blas_trsm('R', 'U', 'N', 'N', mesh%np, nst, R_TOTYPE(M_ONE), ss(1, 1), nst, &
        psi(1, idim, 1), ubound(psi, dim = 1)*st%d%dim)
    end do

    call profiling_count_operations(dble(mesh%np)*dble(nst)**2*(R_ADD + R_MUL))

    SAFE_DEALLOCATE_A(ss)

  case(ORTH_QR)

    ASSERT(states_are_real(st))
    ASSERT(.not. mesh%use_curvilinear)
    ASSERT(.not. mesh%parallel_in_domains)

    nref = min(nst, mesh%np)

    SAFE_ALLOCATE(tau(1:nref))

    tau = M_ZERO

    ! get the optimal size of the work array
    call dgeqrf(mesh%np, nst, psi(1, 1, 1), mesh%np_part, tau(1), tmp, -1, info)
    wsize = nint(tmp)

    ! calculate the QR decomposition
    SAFE_ALLOCATE(work(1:wsize))
    call dgeqrf(mesh%np, nst, psi(1, 1, 1), mesh%np_part, tau(1), work(1), mesh%np*20, info)
    SAFE_DEALLOCATE_A(work)

    ! get the optimal size of the work array
    call dorgqr(mesh%np, nst, nref, psi(1, 1, 1), mesh%np_part, tau(1), tmp, -1, info)
    wsize = nint(tmp)

    ! now calculate Q
    SAFE_ALLOCATE(work(1:wsize))
    call dorgqr(mesh%np, nst, nref, psi(1, 1, 1), mesh%np_part, tau(1), work(1), wsize, info)
    SAFE_DEALLOCATE_A(work)

    SAFE_DEALLOCATE_A(tau)
    
    ! we need to scale by the volume element to get the proper normalization
    psi = psi/sqrt(mesh%vol_pp(1))

  case(ORTH_MGS)
    call mgs()
  end select
  POP_SUB(X(states_orthogonalization_block))
  contains
    subroutine mgs()
      integer :: ist, jst, idim
      R_TYPE, allocatable :: aa(:)
      FLOAT,  allocatable :: bb(:)
#ifdef HAVE_MPI
      R_TYPE, allocatable :: aac(:)
      FLOAT,  allocatable :: bbc(:)
#endif      

      SAFE_ALLOCATE(bb(1:nst))

      do ist = 1, nst
        bb(ist) = X(mf_dotp)(mesh, dim, psi(:, :, ist), psi(:, :, ist), reduce = .false.)
      end do

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then
        SAFE_ALLOCATE(bbc(1:nst))
        call MPI_Allreduce(bb(1), bbc(1), nst, MPI_FLOAT, MPI_SUM, mesh%mpi_grp%comm, mpi_err) 
        bb = bbc
        SAFE_DEALLOCATE_A(bbc)
      end if
#endif

      do ist = 1, nst      
        do idim = 1, dim
          call lalg_scal(mesh%np, M_ONE/sqrt(bb(ist)), psi(:, idim, ist))
        end do
      end do

      SAFE_DEALLOCATE_A(bb)

      SAFE_ALLOCATE(aa(1:nst))

      do ist = 1, nst
        do jst = 1, ist - 1
          aa(jst) = X(mf_dotp)(mesh, dim, psi(:, :, jst), psi(:, :, ist), reduce = .false.)
        end do

#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) then
          SAFE_ALLOCATE(aac(1:nst))
          call MPI_Allreduce(aa(1), aac(1), ist - 1, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
          aa = aac
          SAFE_DEALLOCATE_A(aac)
        end if
#endif

        do jst = 1, ist - 1
          do idim = 1, dim
            call lalg_axpy(mesh%np, -aa(jst), psi(:, idim, jst), psi(:, idim, ist))
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(aa)

    end subroutine mgs
end subroutine X(states_orthogonalization_block)


! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits doing only the orthogonalization (no normalization).
! And one can pass an extra optional argument, mask, which:
!  - on input, if mask(p) = .true., the p-orbital is not used.
!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
! If Theta_Fi and beta_ij are present, it performs the generalized orthogonalization
!   (Theta_Fi - sum_j beta_ij |j><j|Phi>
! This is used in response for metals
subroutine X(states_orthogonalization)(mesh, nst, dim, psi, phi,  &
  normalize, mask, overlap, norm, Theta_fi, beta_ij)
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: nst
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:,:,:)   ! psi(mesh%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     ! phi(mesh%np_part, dim)
  logical, optional, intent(in)    :: normalize
  logical, optional, intent(inout) :: mask(:)      ! mask(nst)
  R_TYPE,  optional, intent(out)   :: overlap(:) 
  R_TYPE,  optional, intent(out)   :: norm
  FLOAT,   optional, intent(in)    :: Theta_Fi
  FLOAT,   optional, intent(in)    :: beta_ij(:)   ! beta_ij(nst)

  logical :: normalize_
  integer :: ist, idim
  FLOAT   :: nrm2
  R_TYPE, allocatable  :: ss(:)
  integer :: block_size, size, sp, ep
  type(profile_t), save :: prof
#ifdef HAVE_MPI
  R_TYPE, allocatable  :: ss_tmp(:)
  type(profile_t), save :: reduce_prof
#endif

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

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(ss_tmp(1:nst))
    call profiling_in(reduce_prof, "GRAM_SCHMIDT_REDUCE")
    call MPI_Allreduce(ss(1), ss_tmp(1), nst, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(reduce_prof)
    ss = ss_tmp
    SAFE_DEALLOCATE_A(ss_tmp)
  end if
#endif

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
! The routine calculates the expectation value of the momentum 
! operator
! <p> = < phi*(ist, k) | -i \nabla | phi(ist, ik) >
!
! Note, the blas routines cdotc, zdotc take care of complex 
! conjugation *. Therefore we pass phi directly.
! ---------------------------------------------------------
subroutine X(states_calc_momentum)(gr, st, momentum)
  type(grid_t),   intent(inout) :: gr
  type(states_t), intent(inout) :: st
  FLOAT,          intent(out)   :: momentum(:,:,:)

  integer             :: idim, ist, ik, idir
  CMPLX               :: expect_val_p
  R_TYPE, allocatable :: grad(:,:,:)
  FLOAT               :: kpoint(1:MAX_DIM)  
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT, allocatable  :: lmomentum(:), gmomentum(:)
  FLOAT, allocatable  :: lmom(:, :, :)
  integer             :: kstart, kend, kn, ndim
#endif

  PUSH_SUB(X(states_calc_momentum))

  SAFE_ALLOCATE(grad(1:gr%mesh%np, 1:st%d%dim, 1:gr%mesh%sb%dim))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      do idim = 1, st%d%dim
        ! compute gradient of st%X(psi)
        call X(derivatives_grad)(gr%der, st%X(psi)(:, idim, ist, ik), grad(:, idim, 1:gr%mesh%sb%dim))
      end do

      do idir = 1, gr%mesh%sb%dim
        ! since the expectation value of the momentum operator is real
        ! for square integrable wfns this integral should be purely imaginary 
        ! for complex wfns but real for real wfns (see case distinction below)
        expect_val_p = X(mf_dotp)(gr%mesh, st%d%dim, &
          st%X(psi)(1:gr%mesh%np, 1:st%d%dim, ist, ik), grad(1:gr%mesh%np, 1:st%d%dim, idir))

        ! In the case of real wavefunctions we do not include the 
        ! -i prefactor of p = -i \nabla
        if (states_are_real(st)) then
          momentum(idir, ist, ik) = real( expect_val_p )
        else
          momentum(idir, ist, ik) = real( -M_zI*expect_val_p )
        end if
      end do

      ! have to add the momentum vector in the case of periodic systems, 
      ! since st%X(psi) contains only u_k
      kpoint = M_ZERO
      kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
      do idir = 1, gr%sb%periodic_dim
        momentum(idir, ist, ik) = momentum(idir, ist, ik) + kpoint(idir)
      end do
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

      do idir = 1, gr%mesh%sb%dim
        lmomentum(1:st%lnst) = momentum(idir, st%st_start:st%st_end, ik)
        call lmpi_gen_allgatherv(st%lnst, lmomentum, tmp, gmomentum, st%mpi_grp)
        momentum(idir, 1:st%nst, ik) = gmomentum(1:st%nst)
      end do

      SAFE_DEALLOCATE_A(lmomentum)
      SAFE_DEALLOCATE_A(gmomentum)
    end if
#endif
  end do
  SAFE_DEALLOCATE_A(grad)

  POP_SUB(X(states_calc_momentum))
end subroutine X(states_calc_momentum)


! ---------------------------------------------------------
! It calculates the expectation value of the angular
! momentum of the state phi. If l2 is passed, it also
! calculates the expectation value of the square of the
! angular momentum of the state phi.
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
    call write_fatal(1)
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
