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
!! $Id: states_inc.F90 2361 2006-08-12 20:24:24Z appel $


! -------------------------------------------------------------
! Returns the dot product of two many-body states; the first
! one is an "excited_state".
!
! WARNING!!!!: periodic systems are not considered in these expressions.
! -------------------------------------------------------------
R_TYPE function X(states_mpdotp_x)(m, excited_state, st, mat) result(dotp)
  type(mesh_t),           intent(in) :: m
  type(excited_states_t), intent(in) :: excited_state
  type(states_t),         intent(in) :: st
  R_TYPE,       optional, intent(in) :: mat(:, :, :)

  integer :: j
  R_TYPE, allocatable :: mat_local(:, :, :)

  call push_sub('states_inc.Xstates_mpdotp_x')

  dotp = M_ZERO

  ASSERT(excited_state%st%d%nik.eq.st%d%nik)

  ALLOCATE(mat_local(excited_state%st%nst, st%nst, st%d%nik), st%nst*excited_state%st%nst*st%d%nik)

  if(present(mat)) then
      mat_local = mat
  else 
    call X(states_matrix)(m, excited_state%st, st, mat_local)
  end if

  do j = 1, excited_state%n_pairs
    call X(states_matrix_swap)(mat_local, excited_state%pair(j))
    dotp = dotp + excited_state%weight(j) * X(states_mpdotp_g)(m, excited_state%st, st, mat_local) 
    call X(states_matrix_swap)(mat_local, excited_state%pair(j))
  end do

  SAFE_DEALLOCATE_A(mat_local)
  call pop_sub()
end function X(states_mpdotp_x)


! -------------------------------------------------------------
! The matrix mat should contain the dot products between two 
! states. One of them (the one operating on the left) is an
! excited state. The pair_t pair indicates the substitution
! of one the occupied spin-orbitals by one unoccupied one
! (electron-hole excitation). This routine returns then the
! matrix that would correspond to the dot products between
! this excited state (on the left) and the same ground state
! (on the right)
! -------------------------------------------------------------
subroutine X(states_matrix_swap)(mat, pair)
  R_TYPE,              intent(inout) :: mat(:, :, :)
  type(states_pair_t), intent(in)    :: pair

  integer :: i, a, ik
  R_TYPE, allocatable :: row(:)

  i  = pair%i
  a  = pair%a
  ik = pair%sigma

  ALLOCATE(row(size(mat, 2)), size(mat, 2))

  ! swap row
  row(:) = mat(i, :, ik)
  mat(i, :, ik) = mat(a, :, ik)
  mat(a, :, ik) = row(:)

  SAFE_DEALLOCATE_A(row)
end subroutine X(states_matrix_swap)


! -------------------------------------------------------------
! Returns <st1 | O | st2>, where both st1 and st2 are Slater
! determinants represented by states_t st1 and st2. O is a
! one-body operator.
!
! The auxiliary Slater determinant opst2 is formed by the
! orbitals that result of applying operator O on each of the
! spin-orbitals of st2.
!
! The routine directly applies Lowdin`s formula [P.-O. Lowdin,
! Phys. Rev. 97, 1474; Eq. 49].
! -------------------------------------------------------------
R_TYPE function X(states_mpmatrixelement_g)(m, st1, st2, opst2) result(st1opst2)
  type(mesh_t),     intent(in) :: m
  type(states_t),   intent(in) :: st1, st2, opst2

  integer :: ispin, nik, nst, ik, i1, j1, k1, i2, j2, k2, i, j
  integer, allocatable :: filled1(:), filled2(:), &
                          partially_filled1(:), partially_filled2(:), &
                          half_filled1(:), half_filled2(:)
  R_TYPE, allocatable :: overlap_mat(:, :, :), op_mat(:, :, :)
  R_TYPE, allocatable :: b(:, :), c(:, :)
  R_TYPE :: z, det

  call push_sub('excited_states_inc.Xstates_mpmatrixelement_g')

  ! This should go away whenever the subroutine is ready.
  st1opst2 = R_TOTYPE(M_ONE)

  ispin = st1%d%ispin
  ASSERT(ispin.eq.st2%d%ispin)
  nik   = st1%d%nik
  ASSERT(nik.eq.st2%d%nik)
  ! Can only consider the number of states of the state that comes with fewer states.
  nst = min(st1%nst, st2%nst)

  ALLOCATE(overlap_mat(st1%nst, st2%nst, st1%d%nik), st1%nst*st2%nst*st1%d%nik)
  ALLOCATE(op_mat(st1%nst, st2%nst, st1%d%nik), st1%nst*st2%nst*st1%d%nik)

  ALLOCATE(filled1(nst), nst)
  ALLOCATE(filled2(nst), nst)
  ALLOCATE(partially_filled1(nst), nst)
  ALLOCATE(partially_filled2(nst), nst)
  ALLOCATE(half_filled1(nst), nst)
  ALLOCATE(half_filled2(nst), nst)

  select case(ispin)
  case(UNPOLARIZED)

    call X(states_matrix)(m, st1, st2, overlap_mat)
    call X(states_matrix)(m, st1, opst2, op_mat)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call write_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call write_fatal(1)
      end if

      ALLOCATE(b(i1+k1, i1+k1), (i1+k1)*(i1+k1))
      ALLOCATE(c(i1+k1, i1+k1), (i1+k1)*(i1+k1))
      do i = 1, i1
        do j = 1, i1
          b(i, j) = op_mat(filled1(i), filled2(j), ik)
          c(i, j) = overlap_mat(filled1(i), filled2(j), ik)
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = op_mat(filled1(i), half_filled2(j), ik)
          c(i, j) = overlap_mat(filled1(i), half_filled2(j), ik)
        end do
      end do
      do i = i1 + 1, i1 + k1
        do j = 1, i1
          b(i, j) = op_mat(half_filled1(i), filled2(j), ik)
          c(i, j) = overlap_mat(half_filled1(i), filled2(j), ik)
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = op_mat(half_filled1(i), half_filled2(j), ik)
          c(i, j) = overlap_mat(half_filled1(i), half_filled2(j), ik)
        end do
      end do

      det = lalg_determinant(i1+k1, c, invert = .true.)
      c = det * transpose(c)

      ! And now, apply Lowdin`s formula.
      ! <U|O|V> = 2 D`_{UV} \sum_{jk} <phi^U_j|o|phi^V_k> D`_{UV}(j|k)
      ! where D`_{UV} is the determinant of the overlap matrix between the
      ! spatial orbitals of U and V (dimension = N/2), and D`_{UV}(j|k) is
      ! the (j,k) minor of this matrix.
      z = M_ZERO
      do i = 1, i1 + k1
        do j = 1, i1 + k1
           z = z + b(i, j) * c(i, j) * (-1)**(i+j)
        end do
      end do
      z = M_TWO * det * z

      st1opst2 = st1opst2 * z ** st1%d%kweights(ik)

    end do

  case(SPIN_POLARIZED, SPINORS)

    call X(states_matrix) (m, st1, st2, overlap_mat)
    call X(states_matrix) (m, st1, opst2, op_mat)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call write_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call write_fatal(1)
      end if

      if(i1 > 0) then
        ALLOCATE(b(i1, i1), i1*i1)
        ALLOCATE(c(i1, i1), i1*i1)
        do i = 1, i1
          do j = 1, i1
            b(i, j) = op_mat(filled1(i), filled2(j), ik)
            c(i, j) = overlap_mat(filled1(i), filled2(j), ik)
          end do
        end do

        ! Get the matrix of cofactors.
        z = lalg_determinant(i1, c, invert = .true.)
        c = z * transpose(c)

        ! And now, apply Lowdin`s formula.
        z = M_ZERO
        do i = 1, i1
          do j = 1, i1
             z = z + b(i, j) * c(i, j) * (-1)**(i+j)
          end do
        end do

        st1opst2 = st1opst2 * z ** st1%d%kweights(ik)

        SAFE_DEALLOCATE_A(b)
        SAFE_DEALLOCATE_A(c)
      end if

    end do
  end select  


  SAFE_DEALLOCATE_A(overlap_mat)
  SAFE_DEALLOCATE_A(op_mat)
  SAFE_DEALLOCATE_A(filled1)
  SAFE_DEALLOCATE_A(filled2)
  SAFE_DEALLOCATE_A(partially_filled1)
  SAFE_DEALLOCATE_A(partially_filled2)
  SAFE_DEALLOCATE_A(half_filled1)
  SAFE_DEALLOCATE_A(half_filled2)
  call pop_sub()
end function X(states_mpmatrixelement_g)


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
  R_TYPE, allocatable :: a(:, :, :), b(:, :)
  call push_sub('excited_states_inc.Xstates_mpdotp')

  ispin = st1%d%ispin
  ASSERT(ispin.eq.st2%d%ispin)
  nik   = st1%d%nik
  ASSERT(nik.eq.st2%d%nik)
  ! Can only consider the number of states of the state that comes with fewer states.
  nst = min(st1%nst, st2%nst)

  ALLOCATE(a(st1%nst, st2%nst, st1%d%nik), st1%nst*st2%nst*st1%d%nik)
  dotp = M_ONE

  ALLOCATE(filled1(nst), nst)
  ALLOCATE(filled2(nst), nst)
  ALLOCATE(partially_filled1(nst), nst)
  ALLOCATE(partially_filled2(nst), nst)
  ALLOCATE(half_filled1(nst), nst)
  ALLOCATE(half_filled2(nst), nst)

  if(present(mat)) then
    a(1:st1%nst, 1:st2%nst, 1:st1%d%nik) = mat(1:st1%nst, 1:st2%nst, 1:st1%d%nik)
  else
    call X(states_matrix) (m, st1, st2, a)
  end if

  select case(ispin)
  case(UNPOLARIZED)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call write_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call write_fatal(1)
      end if

      ALLOCATE(b(i1+k1, i1+k1), (i1+k1)*(i1+k1))
      do i = 1, i1
        do j = 1, i1
          b(i, j) = a(filled1(i), filled2(j), ik)
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = a(filled1(i), half_filled2(j), ik)
        end do
      end do
      do i = i1 + 1, i1 + k1
        do j = 1, i1
          b(i, j) = a(half_filled1(i), filled2(j), ik)
        end do
        do j = i1 + 1, i1 + k1
          b(i, j) = a(half_filled1(i), half_filled2(j), ik)
        end do
      end do

      dotp = dotp * (lalg_determinant(i1+k1, b, invert = .false.)) ** st1%d%kweights(ik)
      if(i1 > 0) then
        dotp = dotp * (lalg_determinant(i1, b(1:i1, 1:i1), invert = .false.)) ** st1%d%kweights(ik)
      end if

    end do
  case(SPIN_POLARIZED, SPINORS)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
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
            b(i, j) = a(filled1(i), filled2(j), ik)
          end do
        end do

        dotp = dotp * lalg_determinant(i1, b, invert = .false.) ** st1%d%kweights(ik)
        SAFE_DEALLOCATE_A(b)
      end if

    end do
  end select

  SAFE_DEALLOCATE_A(a)
  SAFE_DEALLOCATE_A(filled1)
  SAFE_DEALLOCATE_A(filled2)
  SAFE_DEALLOCATE_A(partially_filled1)
  SAFE_DEALLOCATE_A(partially_filled2)
  SAFE_DEALLOCATE_A(half_filled1)
  SAFE_DEALLOCATE_A(half_filled2)
  call pop_sub()
end function X(states_mpdotp_g)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
