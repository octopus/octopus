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
! Orthonormalizes nst orbitals in mesh m (honours state
! parallelization).
subroutine X(states_gram_schmidt_full)(st, nst, m, dim, psi, start)
  type(states_t),    intent(in)    :: st
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE, target,    intent(inout) :: psi(:, :, st%st_start:)   ! psi(m%np_part, dim, nst)
  integer, optional, intent(in)    :: start

  R_TYPE, allocatable :: ss(:, :), qq(:, :), psi_tmp(:, :, :)
  type(profile_t), save :: prof
  integer :: idim, ist, jst, kst, start_
  FLOAT   :: nrm2

  call profiling_in(prof, "GRAM_SCHMIDT_FULL")
  call push_sub('states_calc_inc.Xstates_gram_schmidt_full')

  start_ = 1
  if(present(start)) start_ = start

  ALLOCATE(ss(st%nst, st%nst), st%nst*st%nst)
  ALLOCATE(qq(st%nst, st%nst), st%nst*st%nst)

  ! calculate the matrix of dot products
  call states_blockt_mul(m, st, st%st_start, st%st_end, st%st_start, st%st_end, psi, psi, ss, symm = .true.)

  qq = M_ZERO
  do ist = 1, st%nst

    if(ist < start_ .or. ist > nst) then
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
    nrm2 = sqrt(nrm2)

    ! now generate the matrix with the linear combination of orbitals
    qq(ist, ist) = M_ONE/nrm2
    do jst = 1, ist - 1
      qq(jst, ist) = -ss(jst, ist)/ss(jst, jst)/nrm2
    end do

  end do
  
  if(.not. st%parallel_in_states) then

    call X(states_linear_combination)(m, st%nst, st%d%dim, qq, psi)

  else

    ALLOCATE(psi_tmp(m%np_part, dim, st%st_start:st%st_end), m%np_part*dim*st%lnst)
    
    do ist = st%st_start, st%st_end
      do idim = 1, dim
        call lalg_copy(m%np, psi(:, idim, ist), psi_tmp(:, idim, ist))
      end do
    end do

    call states_block_matr_mul(m, st, st%st_start, st%st_end, st%st_start, st%st_end, psi_tmp, qq, psi)

    deallocate(psi_tmp)

  end if

  call pop_sub()
  call profiling_out(prof)
end subroutine X(states_gram_schmidt_full)


! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits doing only the orthogonalization (no normalization).
! And one can pass an extra optional argument, mask, which:
!  - on input, if mask(p) = .true., the p-orbital is not used.
!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
! If Theta_Fi and beta_ij are present, it performs the generalized orthogonalization
!   (Theta_Fi - sum_j beta_ij |j><j|Phi>
! This is used in response for metals
subroutine X(states_gram_schmidt)(m, nst, dim, psi, phi,  &
  normalize, mask, overlap, norm, Theta_fi, beta_ij)
  type(mesh_t),      intent(in)    :: m
  integer,           intent(in)    :: nst
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     ! phi(m%np_part, dim)
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
  call push_sub('states_calc_inc.Xstates_gram_schmidt')

  ! This routine uses blocking to optimize cache usage. One block of
  ! |phi> is loaded in cache L1 and then then we calculate the dot
  ! product of it with the corresponding blocks of |psi_k>, next we
  ! load another block and do the same. This way we only have to load
  ! |psi> from the L2 or memory.

  block_size = hardware%X(block_size)

  ALLOCATE(ss(1:nst), nst)

  ss = M_ZERO

  if(.not. m%use_curvlinear) then

    do idim = 1, dim
      do sp = 1, m%np, block_size
        size = min(block_size, m%np - sp + 1)
        do ist = 1, nst

          if(present(mask)) then
            if(mask(ist)) cycle
          end if

          ss(ist) = ss(ist) + blas_dot(size, psi(sp, idim, ist), 1, phi(sp, idim), 1)
        end do
      end do
    end do

    ss = ss * m%vol_pp(1)

    call profiling_count_operations((R_ADD + R_MUL) * m%np * dim * nst)

  else

    do idim = 1, dim
      do sp = 1, m%np, block_size
        size = min(block_size, m%np - sp + 1)
        ep = sp - 1 + size
        do ist = 1, nst

          if(present(mask)) then
            if(mask(ist)) cycle
          end if

          ss(ist) = ss(ist) + sum(m%vol_pp(sp:ep)*R_CONJ(psi(sp:ep, idim, ist))*phi(sp:ep, idim))

        end do
      end do
    end do

    call profiling_count_operations((R_ADD + 2 * R_MUL) * m%np * dim * nst)

  end if

#ifdef HAVE_MPI
  if(m%parallel_in_domains) then
    ALLOCATE(ss_tmp(1:nst), nst)
    call profiling_in(reduce_prof, "GRAM_SCHMIDT_REDUCE")
    call MPI_Allreduce(ss, ss_tmp, nst, R_MPITYPE, MPI_SUM, m%vp%comm, mpi_err)
    call profiling_out(reduce_prof)
    ss = ss_tmp
    deallocate(ss_tmp)
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
    do sp = 1, m%np, block_size
      size = min(block_size, m%np - sp + 1)

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

  call profiling_count_operations((R_ADD + R_MUL) * m%np * dim * nst)

  ! the following ifs cannot be given as a single line (without the
  ! then) to avoid a bug in xlf 10.1

  normalize_ = .false.
  if(present(normalize)) then
    normalize_ = normalize
  end if

  if(normalize_) then
    nrm2 = X(mf_nrm2)(m, dim, phi)
    do idim = 1, dim
      call lalg_scal(m%np, M_ONE / nrm2, phi(:, idim))
    end do
  end if

  if(present(overlap)) then
    overlap(1:nst) = ss(1:nst)
  end if

  if(present(norm)) then
    ASSERT(normalize)
    norm = nrm2
  end if

  deallocate(ss)

  call pop_sub()
  call profiling_out(prof)
end subroutine X(states_gram_schmidt)

! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(m, dim, psi)
  type(mesh_t),    intent(in)    :: m
  integer,         intent(in)    :: dim
  R_TYPE,          intent(inout) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim, ip

  call push_sub('states_calc_inc.Xstates_normalize_orbital')

  norm = X(mf_nrm2) (m, dim, psi)
  norm = sqrt(norm)

  forall (idim = 1:dim, ip = 1:m%np) psi(ip, idim) = psi(ip, idim)/norm
  
  call pop_sub()
end subroutine X(states_normalize_orbital)

! ---------------------------------------------------------
FLOAT function X(states_residue)(m, dim, hf, e, f) result(r)
  type(mesh_t),      intent(in)  :: m
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:), f(:,:)
  FLOAT,             intent(in)  :: e

  R_TYPE, allocatable :: res(:,:)
  type(profile_t), save :: prof
  integer :: ip, idim

  call push_sub('states_calc_inc.Xstates_residue')

  call profiling_in(prof, "RESIDUE")

  ALLOCATE(res(m%np_part, dim), m%np_part*dim)

  forall (idim = 1:dim, ip = 1:m%np) res(ip, idim) = hf(ip, idim) - e*f(ip, idim)

  call profiling_count_operations(dim*m%np*(R_ADD + R_MUL))

  r = X(mf_nrm2)(m, dim, res)
  deallocate(res)

  call profiling_out(prof)

  call pop_sub()

end function X(states_residue)


! ---------------------------------------------------------
! The routine calculates the expectation value of the momentum 
! operator
! <p> = < phi*(ist, k) | -i \nabla | phi(ist, ik) >
!
! Note, the blas routines cdotc, zdotc take care of complex 
! conjugation *. Therefore we pass phi directly.
! ---------------------------------------------------------
subroutine X(states_calc_momentum)(gr, st)
  type(grid_t),   intent(inout) :: gr
  type(states_t), intent(inout) :: st

  integer             :: idim, ist, ik, i
  CMPLX               :: expect_val_p
  R_TYPE, allocatable :: grad(:,:,:)  
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT, allocatable  :: lmomentum(:), gmomentum(:)
#endif

  call push_sub('states_calc_inc.Xstates_calc_momentum')

  ALLOCATE(grad(gr%mesh%np, st%d%dim, gr%mesh%sb%dim), gr%mesh%np*st%d%dim*gr%mesh%sb%dim)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      do idim = 1, st%d%dim
        ! compute gradient of st%X(psi)
        call X(derivatives_grad)(gr%der, st%X(psi)(:, idim, ist, ik), grad(:, idim, 1:gr%mesh%sb%dim))
      end do

      do i = 1, gr%mesh%sb%dim
        ! since the expectation value of the momentum operator is real
        ! for square integrable wfns this integral should be purely imaginary 
        ! for complex wfns but real for real wfns (see case distinction below)
        expect_val_p = X(mf_dotp)(gr%mesh, st%d%dim, &
          st%X(psi)(1:gr%mesh%np, 1:st%d%dim, ist, ik), grad(1:gr%mesh%np, 1:st%d%dim, i))

        ! In the case of real wave functions we do not include the 
        ! -i prefactor of p = -i \nabla
        if (st%wfs_type == M_REAL) then
          st%momentum(i, ist, ik) = real( expect_val_p )
        else
          st%momentum(i, ist, ik) = real( -M_zI*expect_val_p )
        end if
      end do

      ! have to add the momentum vector in the case of periodic systems, 
      ! since st%X(psi) contains only u_k
      do i = 1, gr%sb%periodic_dim
        st%momentum(i, ist, ik) = st%momentum(i, ist, ik) + st%d%kpoints(i, ik)
      end do
    end do

    ! Exchange momenta in the state parallel case.
#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      ALLOCATE(lmomentum(1:st%lnst), st%lnst)
      ALLOCATE(gmomentum(1:st%nst), st%nst)

      do i = 1, gr%mesh%sb%dim
        lmomentum(1:st%lnst) = st%momentum(i, st%st_start:st%st_end, ik)
        call lmpi_gen_allgatherv(st%lnst, lmomentum, tmp, gmomentum, st%mpi_grp)
        st%momentum(i, 1:st%nst, ik) = gmomentum(1:st%nst)
      end do

      deallocate(lmomentum)
      deallocate(gmomentum)
    end if
#endif
  end do

  deallocate(grad)

  call pop_sub()
end subroutine X(states_calc_momentum)


! ---------------------------------------------------------
! It calculates the expectation value of the angular
! momentum of the state phi. If l2 is passed, it also
! calculates the expectation value of the square of the
! angular momentum of the state phi.
! ---------------------------------------------------------
subroutine X(states_angular_momentum)(gr, phi, l, l2)
  type(grid_t), intent(inout)  :: gr
  R_TYPE,       intent(inout)  :: phi(:, :)
  FLOAT,        intent(out)    :: l(MAX_DIM)
  FLOAT, optional, intent(out) :: l2

  integer :: idim, dim
  R_TYPE, allocatable :: lpsi(:, :)

  call push_sub('states_calc_inc.Xstates_angular_momemtum')

  ASSERT(gr%mesh%sb%dim .ne.1)

  select case(gr%mesh%sb%dim)
  case(3)
    ALLOCATE(lpsi(gr%mesh%np_part, 3), gr%mesh%np_part*3)
  case(2)
    ALLOCATE(lpsi(gr%mesh%np_part, 1), gr%mesh%np_part*1)
  end select

  dim = size(phi, 2)

  l = M_ZERO
  if(present(l2)) l2 = M_ZERO

  do idim = 1, dim
#if defined(R_TREAL)
    l = M_ZERO
#else
    call X(f_angular_momentum)(gr%sb, gr%mesh, gr%der, phi(:, idim), lpsi)
    select case(gr%mesh%sb%dim)
    case(3)
      l(1) = l(1) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
      l(2) = l(2) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 2))
      l(3) = l(3) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 3))
    case(2)
      l(3) = l(3) + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
    end select
#endif
    if(present(l2)) then
      call X(f_l2)(gr%sb, gr%mesh, gr%der, phi(:, idim), lpsi(:, 1))
      l2 = l2 + X(mf_dotp)(gr%mesh, phi(:, idim), lpsi(:, 1))
    end if
  end do

  deallocate(lpsi)
  call pop_sub()
end subroutine X(states_angular_momentum)


! ---------------------------------------------------------
subroutine X(states_matrix)(m, st1, st2, a)
  type(mesh_t),   intent(in)  :: m
  type(states_t), intent(in)  :: st1, st2
  R_TYPE,         intent(out) :: a(:, :, :)

  integer :: i, j, dim, n1, n2, ik
#if defined(HAVE_MPI)
  R_TYPE, allocatable :: phi2(:, :)
  integer :: k, l
  integer :: status(MPI_STATUS_SIZE)
  integer :: request
#endif

  call push_sub('states_calc_inc.Xstates_matrix')

  n1 = st1%nst
  n2 = st2%nst

  dim = st1%d%dim

  do ik = 1, st1%d%nik

  if(st1%parallel_in_states) then

#if defined(HAVE_MPI)
    call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
    ! Each process sends the states in st2 to the rest of the processes.
    do i = st1%st_start, st1%st_end
      do j = 0, st1%mpi_grp%size - 1
        if(st1%mpi_grp%rank.ne.j) then
          call MPI_Isend(st2%X(psi)(1, 1, i, ik), st1%d%dim*m%np, R_MPITYPE, &
            j, i, st1%mpi_grp%comm, request, mpi_err)
        end if
      end do
    end do

    ! Processes are received, and then the matrix elements are calculated.
    ALLOCATE(phi2(m%np, st1%d%dim), m%np*st1%d%dim)
    do j = 1, n2
      l = st1%node(j)
      if(l.ne.st1%mpi_grp%rank) then
        call MPI_Irecv(phi2(1, 1), st1%d%dim*m%np, R_MPITYPE, l, j, st1%mpi_grp%comm, request, mpi_err)
        call MPI_Wait(request, status, mpi_err)
      else
        phi2(:, :) = st2%X(psi)(:, :, j, ik)
      end if
      do i = st1%st_start, st1%st_end
        a(i, j, ik) = X(mf_dotp)(m, dim, st1%X(psi)(:, :, i, ik), phi2(:, :))
      end do
    end do
    deallocate(phi2)

    ! Each process holds some lines of the matrix. So it is broadcasted (All processes
    ! should get the whole matrix)
    call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
    do i = 1, n1
      k = st1%node(i)
      do j = 1, n2
        call MPI_Bcast(a(i, j, ik), 1, R_MPITYPE, k, st1%mpi_grp%comm, mpi_err)
      end do
    end do
#else
    write(message(1), '(a)') 'Internal error at Xstates_matrix'
    call write_fatal(1)
#endif

  else
    do i = 1, n1
      do j = 1, n2
        a(i, j, ik) = X(mf_dotp)(m, dim, st1%X(psi)(:, :, i, ik), st2%X(psi)(:, :, j, ik))
      end do
    end do
  end if

  end do

  call pop_sub()
end subroutine X(states_matrix)

! ------------------------------------------------------------------
! This routine forms a new sets of states as linear combination of the
! previous one, the coefficients are given by the matrix tranfs. This
! combination is done in place, so it does not require a second set of
! wavefunctions to be allocated.
!
! The performance of this function could be improved using blocks and
! blas, but for the moment this is not necessary.
!
subroutine X(states_linear_combination)(mesh, nst, dim, transf, psi)
  type(mesh_t),        intent(in)    :: mesh
  integer,             intent(in)    :: nst
  integer,             intent(in)    :: dim
  R_TYPE,              intent(in)    :: transf(:, :)
  R_TYPE,              intent(inout) :: psi(:, :, :)
  
  R_TYPE, allocatable :: psinew(:, :)
  
  integer :: ist, idim, block_size, size, sp
  type(profile_t), save :: prof

  call profiling_in(prof, "STATES_ROTATE")

  block_size = hardware%X(block_size)

  ALLOCATE(psinew(block_size, 1:nst), block_size*nst)
  
  do sp = 1, mesh%np, block_size
    size = min(block_size, mesh%np - sp + 1)

    do idim = 1, dim
      
      call blas_gemm('N', 'N', &
           size, nst, nst, &
           R_TOTYPE(M_ONE), psi(sp, idim, 1), mesh%np_part*dim, &
           transf(1, 1), nst, &
           R_TOTYPE(M_ZERO), psinew(1, 1), block_size)
      
      do ist = 1, nst
        call blas_copy(size, psinew(1, ist), 1, psi(sp, idim, ist), 1)
      end do
     
    end do

  end do

  call profiling_count_operations((R_ADD + R_MUL)*dble(mesh%np)*dim*nst**2)

  call profiling_out(prof)

end subroutine X(states_linear_combination)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
