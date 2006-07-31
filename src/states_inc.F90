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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$


! ---------------------------------------------------------
! Orthonormalizes nst orbital in mesh m
subroutine X(states_gram_schmidt1)(nst, m, dim, psi, start)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  integer, optional, intent(in)    :: start

  integer :: p, q, stst, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT1)
  call push_sub('states_inc.Xstates_gram_schmidt1')

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = X(states_dotp)(m, dim, psi(:,:, q), psi(:,:, p))
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), psi(:, idim, p))
      end do
    end do

    nrm2 = X(states_nrm2)(m, dim, psi(:,:, p))
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, psi(:, idim, p))
    end do
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT1)
end subroutine X(states_gram_schmidt1)

! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits to do only the orthogonalization (no normalization).
! And one can pass an extra optional argument, mask, which:
!  - on input, if mask(p) = .true., the p-orbital is not used.
!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
subroutine X(states_gram_schmidt2)(nst, m, dim, psi, phi, normalize, mask)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     ! phi(m%np_part, dim)
  logical, optional, intent(in)    :: normalize
  logical, optional, intent(inout) :: mask(:)      ! nst

  logical :: normalize_
  integer :: q, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT2)
  call push_sub('states_inc.Xstates_gram_schmidt2')

  do q = 1, nst
    if(present(mask)) then
      if(mask(q)) cycle
    end if
    ss = X(states_dotp)(m, dim, psi(:,:, q), phi)
    if(abs(ss) > CNST(1.0e-13)) then
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), phi(:, idim))
      end do
    else
      if(present(mask)) mask(q) = .true.
    end if
  end do

  normalize_ = .false.
  if(present(normalize)) normalize_ = normalize
  if(normalize) then
    nrm2 = X(states_nrm2)(m, dim, phi)
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, phi(:, idim))
    end do
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT2)
end subroutine X(states_gram_schmidt2)


! ---------------------------------------------------------
R_TYPE function X(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f1(:,:), f2(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim))
  end do

  call pop_sub()

end function X(states_dotp)


! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(m, dim, psi)
  type(mesh_t),    intent(in)  :: m
  integer,         intent(in)  :: dim
  R_TYPE,          intent(out) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim

  call push_sub('states_inc.Xstates_normalize_orbital')

  norm = X(states_nrm2) (m, dim, psi)
  norm = sqrt(norm)

  do idim = 1, dim
    psi(1:m%np, idim) = psi(1:m%np, idim)/norm
  end do

  call pop_sub()
end subroutine X(states_normalize_orbital)


! ---------------------------------------------------------
FLOAT function X(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_nrm2')

  nrm2 = M_ZERO
  do idim = 1, dim
    nrm2 = nrm2 + X(mf_nrm2)(m, f(:, idim))**2
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

  call push_sub('states_inc.Xstates_residue')

  ALLOCATE(res(m%np_part, dim), m%np_part*dim)

  res(1:m%np, 1:dim) = hf(1:m%np, 1:dim) - e*f(1:m%np, 1:dim)

  r = X(states_nrm2)(m, dim, res)
  deallocate(res)

  call pop_sub()

end function X(states_residue)


! -------------------------------------------------------------
! Returns the dot product of two many-body states; the first
! one is an "excited_state".
!
! WARNING!!!!: periodic systems are not considered in these expressions.
! -------------------------------------------------------------
R_TYPE function X(states_mpdotp_x)(m, excited_state, st, mat) result(dotp)
  type(mesh_t),   intent(in) :: m
  type(states_excited_t), intent(in) :: excited_state
  type(states_t), intent(in)         :: st
  R_TYPE, optional, intent(in)       :: mat(:, :, :)

  integer :: i, a, sigma, j, ik
  R_TYPE, allocatable :: mat_local(:, :, :), row(:)

  call push_sub('states_inc.Xstates_mpdotp_x')

  dotp = M_ZERO

  ASSERT(excited_state%st%d%nik.eq.st%d%nik)

  ALLOCATE(mat_local(excited_state%st%nst, st%nst, st%d%nik), st%nst*excited_state%st%nst*st%d%nik)
  ALLOCATE(row(st%nst), st%nst)

  if(present(mat)) then
    do ik = 1, st%d%nik
      mat_local = mat
    end do
  else 
    do ik = 1, st%d%nik
      call X(calculate_matrix) (m, ik, excited_state%st, st, mat_local(:, :, ik))
    end do
  end if

  do j = 1, excited_state%n_pairs
    i     = excited_state%pair(j)%i
    a     = excited_state%pair(j)%a
    sigma = excited_state%pair(j)%sigma
    row(:) = mat_local(i, :, sigma)
    mat_local(i, :, sigma) = mat_local(a, :, sigma)
    dotp = dotp + excited_state%weight(j) * X(states_mpdotp_g)(m, excited_state%st, st, mat_local) 
    mat_local(i, :, sigma) = row(:)
  end do

  deallocate(mat_local, row)
  call pop_sub()
end function X(states_mpdotp_x)


! -------------------------------------------------------------
! Returns the dot product of two many-body states st1 and st2.
! Warning: it does not permit fractional occupation numbers.
! -------------------------------------------------------------
R_TYPE function X(states_mpdotp_g)(m, st1, st2, mat) result(dotp)
  type(mesh_t),     intent(in) :: m
  type(states_t),   intent(in) :: st1, st2
  R_TYPE, optional, intent(in) :: mat(:, :, :)

  integer :: ik, ispin, nik, nst, i1, j1, i2, j2, k1, k2, i, j
  integer, allocatable :: filled1(:), filled2(:), partially_filled1(:), partially_filled2(:), &
                          half_filled1(:), half_filled2(:)
  R_TYPE, allocatable :: a(:, :), b(:, :)
  call push_sub('states_inc.Xstates_mpdotp')

  ispin = st1%d%ispin
  ASSERT(ispin.eq.st2%d%ispin)
  nik   = st1%d%nik
  ASSERT(nik.eq.st2%d%nik)
  ! Can only consider the number of states of the state that comes with less states.
  nst = min(st1%nst, st2%nst)

  ALLOCATE(a(st1%nst, st2%nst), st1%nst*st2%nst)
  dotp = M_ONE

  ALLOCATE(filled1(nst), nst)
  ALLOCATE(filled2(nst), nst)
  ALLOCATE(partially_filled1(nst), nst)
  ALLOCATE(partially_filled2(nst), nst)
  ALLOCATE(half_filled1(nst), nst)
  ALLOCATE(half_filled2(nst), nst)

  select case(ispin)
  case(UNPOLARIZED)
    do ik = 1, nik
      if(present(mat)) then
        a(1:st1%nst, 1:st2%nst) = mat(1:st1%nst, 1:st2%nst, ik)
      else
        call X(calculate_matrix) (m, ik, st1, st2, a)
      end if

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calcualte many-body dot products with partially occupied orbitals'
        call write_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call write_fatal(1)
      end if

      ALLOCATE(b(i1+k1, i1+k1), (i1+k1)*(i1+k1))
      do i = 1, i1
        do j = 1, i1
          b(i, j) = a(filled1(i), filled2(j))
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = a(filled1(i), half_filled2(j))
        end do
      end do
      do i = i1 + 1, i1 + k1
        do j = 1, i1
          b(i, j) = a(half_filled1(i), filled2(j))
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = a(half_filled1(i), half_filled2(j))
        end do
      end do

      dotp = dotp * (lalg_determinant(i1+k1, b, invert = .false.)) ** st1%d%kweights(ik)
      if(i1 > 0) then
        dotp = dotp * (lalg_determinant(i1, b(1:i1, 1:i1), invert = .false.)) ** st1%d%kweights(ik)
      end if

    end do
  case(SPIN_POLARIZED, SPINORS)
    do ik = 1, nik

      if(present(mat)) then
        a(1:st1%nst, 1:st2%nst) = mat(1:st1%nst, 1:st2%nst, ik)
      else
        call X(calculate_matrix) (m, ik, st1, st2, a)
      end if

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calcualte many-body dot products with partially occupied orbitals'
        call write_fatal(1)
      end if
      if(i1 .ne. i2) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call write_fatal(1)
      end if

      if(i1 > 0) then
        ALLOCATE(b(i1, i1), i1*i1)
        do i = 1, i1
          do j = 1, i1
            b(i, j) = a(filled1(i), filled2(j))
          end do
        end do

        dotp = dotp * lalg_determinant(i1, b, invert = .false.) ** st1%d%kweights(ik)
        deallocate(b)
      end if

    end do
  end select

  deallocate(a)
  call pop_sub()
end function X(states_mpdotp_g)


! ---------------------------------------------------------
subroutine X(calculate_matrix)(m, ik, st1, st2, a)
  type(mesh_t),   intent(in)  :: m
  integer,        intent(in)  :: ik
  type(states_t), intent(in)  :: st1, st2
  R_TYPE,         intent(out) :: a(:, :)

  integer :: i, j, dim, n1, n2
#if defined(HAVE_MPI)
  R_TYPE, allocatable :: phi2(:, :)
  integer :: k, l
  integer :: status(MPI_STATUS_SIZE)
  integer :: request
#endif

  call push_sub('states_inc.calculate_matrix')

  n1 = st1%nst
  n2 = st2%nst

  dim = st1%d%dim
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
      a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(:, :, i, ik), phi2(:, :))
    end do
  end do
  deallocate(phi2)

  ! Each process holds some lines of the matrix. So it is broadcasted (All processes
  ! should get the whole matrix)
  call MPI_Barrier(st1%mpi_grp%comm, mpi_err)
  do i = 1, n1
    k = st1%node(i)
    do j = 1, n2
      call MPI_Bcast(a(i, j), 1, R_MPITYPE, k, st1%mpi_grp%comm, mpi_err)
    end do
  end do
#else
  do i = 1, n1
    do j = 1, n2
      a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(:, :, i, ik), &
        st2%X(psi)(:, :, j, ik))
    end do
  end do
#endif

  call pop_sub()
end subroutine X(calculate_matrix)


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
