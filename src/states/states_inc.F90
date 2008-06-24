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
  R_TYPE, target,    intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  integer, optional, intent(in)    :: start

  integer :: p, q, stst, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss
#if defined(HAVE_MPI)
  R_TYPE, target, allocatable :: buf(:, :)
#endif
  R_TYPE, pointer             :: psi_q(:, :), psi_p(:, :)
  type(profile_t), save :: prof

  call profiling_in(prof, "GRAM_SCHMIDT_FULL")
  call push_sub('states_inc.Xstates_gram_schmidt_full')

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ALLOCATE(buf(m%np_part, st%d%dim), m%np_part*st%d%dim)
  end if
#endif

  do p = stst, nst
    ! Orthogonalize.
    do q = 1, p - 1
      ! Copy orbital p from other node if necessary.
      if(st%parallel_in_states) then
#if defined(HAVE_MPI)
        if(state_is_local(st, p).and.state_is_local(st, q)) then
          psi_q => psi(:, :, q-st%st_start+1)
          psi_p => psi(:, :, p-st%st_start+1)
        else if(state_is_local(st, p).and..not.state_is_local(st, q)) then
          call MPI_Recv(buf, m%np_part*st%d%dim, R_MPITYPE, st%node(q), 0, st%mpi_grp%comm, mpi_err)
          psi_q => buf
          psi_p => psi(:, :, p-st%st_start+1)
        else if(.not.state_is_local(st, p).and.state_is_local(st, q)) then
          call MPI_Send(psi(:, :, q-st%st_start+1), m%np_part*st%d%dim, R_MPITYPE, st%node(p), &
            0, st%mpi_grp%comm, mpi_err)
        end if
#else
        message(1) = 'Running parallel in states without MPI. This is a bug.'
        call write_fatal(1)
#endif
      else
        psi_q => psi(:, :, q)
        psi_p => psi(:, :, p)
      end if
      
      if(state_is_local(st, p)) then
        ss = X(states_dotp)(m, dim, psi_q(:,:), psi_p(:, :))
        do idim = 1, dim
          call lalg_axpy(m%np, -ss, psi_q(:, idim), psi_p(:, idim))
        end do
      end if
    end do

    ! Normalize.
    if(state_is_local(st, p)) then
      nrm2 = X(states_nrm2)(m, dim, psi(:, :, p-st%st_start+1))
      do idim = 1, dim
        call lalg_scal(m%np, M_ONE/nrm2, psi(:, idim, p-st%st_start+1))
      end do
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    deallocate(buf)
  end if
#endif
  
  call pop_sub()
  call profiling_out(prof)
end subroutine X(states_gram_schmidt_full)


! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits to do only the orthogonalization (no normalization).
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
  integer :: block_size, size, sp
  type(profile_t), save :: prof
#ifdef HAVE_MPI
  R_TYPE, allocatable  :: ss_tmp(:)
  type(profile_t), save :: reduce_prof
#endif

  call profiling_in(prof, "GRAM_SCHMIDT")
  call push_sub('states_inc.Xstates_gram_schmidt')

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

    ss = ss*m%vol_pp(1)

    call profiling_count_operations((R_ADD + R_MUL)*m%np*dim*nst)

  else

    do idim = 1, dim
      do sp = 1, m%np, block_size
        size = min(block_size, m%np - sp + 1)
        do ist = 1, nst

          if(present(mask)) then
            if(mask(ist)) cycle
          end if

          ss(ist) = ss(ist) + sum(m%vol_pp(sp:sp + size)*R_CONJ(psi(sp:sp + size, idim, ist))*phi(sp:sp+size, idim))
        end do
      end do
    end do

    call profiling_count_operations((R_ADD + 2*R_MUL)*m%np*dim*nst)

  end if

#ifdef HAVE_MPI
  if(m%parallel_in_domains) then
    ALLOCATE(ss_tmp(1:nst), nst)
    call profiling_in(reduce_prof, "GRAM_SCHMIDT_REDUCE")
    call MPI_Allreduce(ss, ss_tmp, nst, R_MPITYPE, MPI_SUM, m%vp%comm, mpi_err)
    call profiling_out(reduce_prof)
    ss = ss_tmp
  end if
#endif

  if(present(mask)) then
    do ist = 1, nst
      mask(ist) = (abs(ss(ist)) <= M_EPSILON)
    end do
  end if

  if(present(beta_ij))  &
    ss(:) = ss(:)*beta_ij(:)

  do idim = 1, dim
    do sp = 1, m%np, block_size
      size = min(block_size, m%np - sp + 1)

      if(present(Theta_Fi)) then
        if(Theta_Fi.ne.M_ONE) &
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

  call profiling_count_operations((R_ADD + R_MUL)*m%np*dim*nst)

  normalize_ = .false.
  if(present(normalize)) normalize_ = normalize
  if(normalize_) then
    nrm2 = X(states_nrm2)(m, dim, phi)
    do idim = 1, dim
      call lalg_scal(m%np, M_ONE/nrm2, phi(:, idim))
    end do
  end if

  if(present(overlap)) overlap(1:nst) = ss(1:nst)
  if(present(norm) .and. normalize) norm = nrm2

  deallocate(ss)

  call pop_sub()
  call profiling_out(prof)
end subroutine X(states_gram_schmidt)


! ---------------------------------------------------------
R_TYPE function X(states_dotp)(m, dim, f1, f2, reduce) result(dotp)
  type(mesh_t),      intent(in) :: m
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  call push_sub('states_inc.Xstates_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    if(present(reduce)) then
      dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim), reduce)
    else
      dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim))
    end if
  end do

  call pop_sub()

end function X(states_dotp)


! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(m, dim, psi)
  type(mesh_t),    intent(in)    :: m
  integer,         intent(in)    :: dim
  R_TYPE,          intent(inout) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim

  call push_sub('states_inc.Xstates_normalize_orbital')

  norm = X(states_nrm2) (m, dim, psi)
  norm = sqrt(norm)

  do idim = 1, dim
    !$omp parallel workshare
    psi(1:m%np, idim) = psi(1:m%np, idim)/norm
    !$omp end parallel workshare
  end do

  call pop_sub()
end subroutine X(states_normalize_orbital)


! ---------------------------------------------------------
FLOAT function X(states_nrm2)(m, dim, f, reduce) result(nrm2)
  type(mesh_t),      intent(in) :: m
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  call push_sub('states_inc.Xstates_nrm2')

  nrm2 = M_ZERO
  do idim = 1, dim
    nrm2 = nrm2 + X(mf_nrm2)(m, f(:, idim), reduce)**2
  end do
  nrm2 = sqrt(nrm2)

  call pop_sub()

end function X(states_nrm2)


! ---------------------------------------------------------
FLOAT function X(states_residue)(m, dim, hf, e, f) result(r)
  type(mesh_t),      intent(in)  :: m
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:), f(:,:)
  FLOAT,             intent(in)  :: e

  R_TYPE, allocatable :: res(:,:)
  type(profile_t), save :: prof

  call push_sub('states_inc.Xstates_residue')

  call profiling_in(prof, "RESIDUE")

  ALLOCATE(res(m%np_part, dim), m%np_part*dim)

  !$omp parallel workshare
  res(1:m%np, 1:dim) = hf(1:m%np, 1:dim) - e*f(1:m%np, 1:dim)
  !$omp end parallel workshare

  call profiling_count_operations(dim*m%np*(R_ADD + R_MUL))

  r = X(states_nrm2)(m, dim, res)
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
  FLOAT               :: lmomentum(st%lnst), gmomentum(st%nst)
#endif

  call push_sub('states_inc.Xstates_calc_momentum')

  ALLOCATE(grad(NP, st%d%dim, NDIM), NP*st%d%dim*NDIM)

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end

      do idim = 1, st%d%dim
        ! compute gradient of st%X(psi)
        call X(f_gradient)(gr%sb, gr%f_der, &
          st%X(psi)(1:NP_PART, idim, ist, ik), grad(1:NP, idim, 1:NDIM))
      end do

      do i = 1, NDIM
        ! since the expectation value of the momentum operator is real
        ! for square integrable wfns this integral should be purely imaginary 
        ! for complex wfns but real for real wfns (see case distinction below)
        expect_val_p = X(states_dotp)(gr%m, st%d%dim, &
          st%X(psi)(1:NP, 1:st%d%dim, ist, ik), grad(1:NP, 1:st%d%dim, i))

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
      do i = 1, NDIM
        lmomentum(1:st%lnst) = st%momentum(i, st%st_start:st%st_end, ik)
        call lmpi_gen_alltoallv(st%lnst, lmomentum, tmp, gmomentum, st%mpi_grp)
        st%momentum(i, 1:st%nst, ik) = gmomentum(1:st%nst)
      end do
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

  call push_sub('states_inc.Xstates_angular_momemtum')

  ASSERT(gr%m%sb%dim .ne.1)

  select case(gr%m%sb%dim)
  case(3)
    ALLOCATE(lpsi(NP_PART, 3), NP_PART*3)
  case(2)
    ALLOCATE(lpsi(NP_PART, 1), NP_PART*1)
  end select

  dim = size(phi, 2)

  l = M_ZERO
  if(present(l2)) l2 = M_ZERO

  do idim = 1, dim
#if defined(R_TREAL)
    l = M_ZERO
#else
    call X(f_angular_momentum)(gr%sb, gr%f_der, phi(:, idim), lpsi)
    select case(gr%m%sb%dim)
    case(3)
      l(1) = l(1) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
      l(2) = l(2) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 2))
      l(3) = l(3) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 3))
    case(2)
      l(3) = l(3) + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
    end select
#endif
    if(present(l2)) then
      call X(f_l2)(gr%sb, gr%f_der, phi(:, idim), lpsi(:, 1))
      l2 = l2 + X(mf_dotp)(gr%m, phi(:, idim), lpsi(:, 1))
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

  call push_sub('states_inc.Xstates_matrix')

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
        a(i, j, ik) = X(states_dotp)(m, dim, st1%X(psi)(:, :, i, ik), phi2(:, :))
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
        a(i, j, ik) = X(states_dotp)(m, dim, st1%X(psi)(:, :, i, ik), st2%X(psi)(:, :, j, ik))
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
subroutine X(states_linear_combination)(st, mesh, ik, transf)
  type(states_t),      intent(inout) :: st
  type(mesh_t),        intent(in)    :: mesh
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: transf(:, :)
  
  R_TYPE, allocatable :: psiold(:)
  
  R_TYPE :: aa
  integer :: ist, jst, ip, idim

  ALLOCATE(psiold(st%st_start:st%st_end), st%lnst)
  
  do ip = 1, mesh%np
    do idim = 1, st%d%dim
      
      psiold(st%st_start:st%st_end) = st%X(psi)(ip, idim, st%st_start:st%st_end, ik)
      
      do ist = st%st_start, st%st_end
        aa = M_ZERO
        do jst = st%st_start, st%st_end
          aa = aa +transf(jst, ist)*psiold(jst)
        end do
        st%X(psi)(ip, idim, ist, ik) = aa
      end do
      
    end do
  end do
  
  deallocate(psiold)
  
end subroutine X(states_linear_combination)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
