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
!> Returns the dot product of two many-body states; the first
!! one is an "excited_state".
!!
!! WARNING!!!!: periodic systems are not considered in these expressions.
R_TYPE function X(states_mpdotp_x)(mesh, excited_state, st, mat) result(dotp)
  type(mesh_t),           intent(in) :: mesh
  type(excited_states_t), intent(in) :: excited_state
  type(states_t),         intent(in) :: st
  R_TYPE,       optional, intent(in) :: mat(:, :, :)

  integer :: jj
  R_TYPE, allocatable :: mat_local(:, :, :)

  PUSH_SUB(X(states_mpdotp_x))

  dotp = M_ZERO

  ASSERT(excited_state%st%d%nik.eq.st%d%nik)

  SAFE_ALLOCATE(mat_local(1:excited_state%st%nst, 1:st%nst, 1:st%d%nik))

  if(present(mat)) then
      mat_local = mat
  else 
    call X(states_matrix)(mesh, excited_state%st, st, mat_local)
  end if

  do jj = 1, excited_state%n_pairs
    call X(states_matrix_swap)(mat_local, excited_state%pair(jj))
    dotp = dotp + excited_state%weight(jj) * X(states_mpdotp_g)(mesh, excited_state%st, st, mat_local) 
    call X(states_matrix_swap)(mat_local, excited_state%pair(jj))
  end do

  SAFE_DEALLOCATE_A(mat_local)
  POP_SUB(X(states_mpdotp_x))
end function X(states_mpdotp_x)


! -------------------------------------------------------------
!> The matrix mat should contain the dot products between two 
!! states. One of them (the one operating on the left) is an
!! excited state. The pair_t pair indicates the substitution
!! of one the occupied spin-orbitals by one unoccupied one
!! (electron-hole excitation). This routine returns then the
!! matrix that would correspond to the dot products between
!! this excited state (on the left) and the same ground state
!! (on the right)
subroutine X(states_matrix_swap)(mat, pair)
  R_TYPE,              intent(inout) :: mat(:, :, :)
  type(states_pair_t), intent(in)    :: pair

  integer :: ii, aa, ik
  R_TYPE, allocatable :: row(:)

  PUSH_SUB(X(states_matrix_swap))

  ii = pair%i
  aa = pair%a
  ik = pair%sigma

  SAFE_ALLOCATE(row(1:size(mat, 2)))

  ! swap row
  row(:) = mat(ii, :, ik)
  mat(ii, :, ik) = mat(aa, :, ik)
  mat(aa, :, ik) = row(:)

  SAFE_DEALLOCATE_A(row)
  POP_SUB(X(states_matrix_swap))
end subroutine X(states_matrix_swap)


! -------------------------------------------------------------
!> Returns <st1 | O | st2>, where both st1 and st2 are Slater
!! determinants represented by states_t st1 and st2. O is a
!! one-body operator.
!!
!! The auxiliary Slater determinant opst2 is formed by the
!! orbitals that result of applying operator O on each of the
!! spin-orbitals of st2.
!!
!! The routine directly applies Lowdin`s formula [P.-O. Lowdin,
!! Phys. Rev. 97, 1474; Eq. 49].
R_TYPE function X(states_mpmatrixelement_g)(mesh, st1, st2, opst2) result(st1opst2)
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st1, st2, opst2

  integer :: ispin, nik, nst, ik, i1, j1, k1, i2, j2, k2, ii, jj
  integer, allocatable :: filled1(:), filled2(:), &
                          partially_filled1(:), partially_filled2(:), &
                          half_filled1(:), half_filled2(:)
  R_TYPE, allocatable :: overlap_mat(:, :, :), op_mat(:, :, :)
  R_TYPE, allocatable :: bb(:, :), cc(:, :)
  R_TYPE :: zz, det

  PUSH_SUB(X(states_mpmatrixelement_g))

  ! This should go away whenever the subroutine is ready.
  st1opst2 = R_TOTYPE(M_ONE)

  ispin = st1%d%ispin
  ASSERT(ispin .eq. st2%d%ispin)
  nik   = st1%d%nik
  ASSERT(nik .eq. st2%d%nik)
  ! Can only consider the number of states of the state that comes with fewer states.
  nst = min(st1%nst, st2%nst)

  SAFE_ALLOCATE(overlap_mat(1:st1%nst, 1:st2%nst, 1:st1%d%nik))
  SAFE_ALLOCATE(op_mat(1:st1%nst, 1:st2%nst, 1:st1%d%nik))

  SAFE_ALLOCATE(          filled1(1:nst))
  SAFE_ALLOCATE(          filled2(1:nst))
  SAFE_ALLOCATE(partially_filled1(1:nst))
  SAFE_ALLOCATE(partially_filled2(1:nst))
  SAFE_ALLOCATE(     half_filled1(1:nst))
  SAFE_ALLOCATE(     half_filled2(1:nst))

  select case(ispin)
  case(UNPOLARIZED)

    call X(states_matrix)(mesh, st1, st2, overlap_mat)
    call X(states_matrix)(mesh, st1, opst2, op_mat)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call messages_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(bb(1:i1+k1, 1:i1+k1))
      SAFE_ALLOCATE(cc(1:i1+k1, 1:i1+k1))
      do ii = 1, i1
        do jj = 1, i1
          bb(ii, jj) = op_mat(filled1(ii), filled2(jj), ik)
          cc(ii, jj) = overlap_mat(filled1(ii), filled2(jj), ik)
        end do
        do jj = i1 + 1, i1 + k1
          bb(ii, jj) = op_mat(filled1(ii), half_filled2(jj), ik)
          cc(ii, jj) = overlap_mat(filled1(ii), half_filled2(jj), ik)
        end do
      end do
      do ii = i1 + 1, i1 + k1
        do jj = 1, i1
          bb(ii, jj) = op_mat(half_filled1(ii), filled2(jj), ik)
          cc(ii, jj) = overlap_mat(half_filled1(ii), filled2(jj), ik)
        end do
        do jj = i1 + 1, i1 + k1
          bb(ii, jj) = op_mat(half_filled1(ii), half_filled2(jj), ik)
          cc(ii, jj) = overlap_mat(half_filled1(ii), half_filled2(jj), ik)
        end do
      end do

      det = lalg_determinant(i1+k1, cc, invert = .true.)
      cc = det * transpose(cc)

      ! And now, apply Lowdin`s formula.
      ! <U|O|V> = 2 D`_{UV} \sum_{jk} <phi^U_j|o|phi^V_k> D`_{UV}(j|k)
      ! where D`_{UV} is the determinant of the overlap matrix between the
      ! spatial orbitals of U and V (dimension = N/2), and D`_{UV}(j|k) is
      ! the (j,k) minor of this matrix.
      zz = M_ZERO
      do ii = 1, i1 + k1
        do jj = 1, i1 + k1
           zz = zz + bb(ii, jj) * cc(ii, jj) * (-1)**(ii+jj)
        end do
      end do
      zz = M_TWO * det * zz

      st1opst2 = st1opst2 * zz ** st1%d%kweights(ik)

      SAFE_DEALLOCATE_A(bb)
      SAFE_DEALLOCATE_A(cc)

    end do

  case(SPIN_POLARIZED, SPINORS)

    call X(states_matrix) (mesh, st1, st2, overlap_mat)
    call X(states_matrix) (mesh, st1, opst2, op_mat)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call messages_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call messages_fatal(1)
      end if

      if(i1 > 0) then
        SAFE_ALLOCATE(bb(1:i1, 1:i1))
        SAFE_ALLOCATE(cc(1:i1, 1:i1))
        do ii = 1, i1
          do jj = 1, i1
            bb(ii, jj) = op_mat(filled1(ii), filled2(jj), ik)
            cc(ii, jj) = overlap_mat(filled1(ii), filled2(jj), ik)
          end do
        end do

        ! Get the matrix of cofactors.
        zz = lalg_determinant(i1, cc, invert = .true.)
        cc = zz * transpose(cc)

        ! And now, apply Lowdin`s formula.
        zz = M_ZERO
        do ii = 1, i1
          do jj = 1, i1
             zz = zz + bb(ii, jj) * cc(ii, jj) * (-1)**(ii+jj)
          end do
        end do

        st1opst2 = st1opst2 * zz ** st1%d%kweights(ik)

        SAFE_DEALLOCATE_A(bb)
        SAFE_DEALLOCATE_A(cc)
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

  POP_SUB(X(states_mpmatrixelement_g))
end function X(states_mpmatrixelement_g)


! -------------------------------------------------------------
!> Returns the dot product of two many-body states st1 and st2.
!! \warning: it does not permit fractional occupation numbers.
! -------------------------------------------------------------
R_TYPE function X(states_mpdotp_g)(mesh, st1, st2, mat) result(dotp)
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st1, st2
  R_TYPE, optional, intent(in) :: mat(:, :, :)

  integer :: ik, ispin, nik, nst, i1, j1, i2, j2, k1, k2, ii, jj
  integer, allocatable :: filled1(:), filled2(:), partially_filled1(:), partially_filled2(:), &
                          half_filled1(:), half_filled2(:)
  R_TYPE, allocatable :: aa(:, :, :), bb(:, :)

  PUSH_SUB(X(states_mpdotp_g))

  ispin = st1%d%ispin
  ASSERT(ispin .eq. st2%d%ispin)
  nik   = st1%d%nik
  ASSERT(nik .eq. st2%d%nik)
  ! Can only consider the number of states of the state that comes with fewer states.
  nst = min(st1%nst, st2%nst)

  SAFE_ALLOCATE(aa(1:st1%nst, 1:st2%nst, 1:st1%d%nik))
  dotp = M_ONE

  SAFE_ALLOCATE(          filled1(1:nst))
  SAFE_ALLOCATE(          filled2(1:nst))
  SAFE_ALLOCATE(partially_filled1(1:nst))
  SAFE_ALLOCATE(partially_filled2(1:nst))
  SAFE_ALLOCATE(     half_filled1(1:nst))
  SAFE_ALLOCATE(     half_filled2(1:nst))

  if(present(mat)) then
    aa(1:st1%nst, 1:st2%nst, 1:st1%d%nik) = mat(1:st1%nst, 1:st2%nst, 1:st1%d%nik)
  else
    call X(states_matrix) (mesh, st1, st2, aa)
  end if

  select case(ispin)
  case(UNPOLARIZED)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call messages_fatal(1)
      end if
      if(  (i1 .ne. i2)  .or.  (k1 .ne. k2) ) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp_g'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(bb(1:i1+k1, 1:i1+k1))
      do ii = 1, i1
        do jj = 1, i1
          bb(ii, jj) = aa(filled1(ii), filled2(jj), ik)
        end do
        do jj = i1 + 1, i1 + k1
          bb(ii, jj) = aa(filled1(ii), half_filled2(jj), ik)
        end do
      end do
      do ii = i1 + 1, i1 + k1
        do jj = 1, i1
          bb(ii, jj) = aa(half_filled1(ii), filled2(jj), ik)
        end do
        do jj = i1 + 1, i1 + k1
          bb(ii, jj) = aa(half_filled1(ii), half_filled2(jj), ik)
        end do
      end do

      dotp = dotp * (lalg_determinant(i1+k1, bb, invert = .false.)) ** st1%d%kweights(ik)
      if(i1 > 0) then
        dotp = dotp * (lalg_determinant(i1, bb(1:i1, 1:i1), invert = .false.)) ** st1%d%kweights(ik)
      end if

    end do
  case(SPIN_POLARIZED, SPINORS)

    do ik = 1, nik

      call occupied_states(st1, ik, i1, j1, k1, filled1, partially_filled1, half_filled1)
      call occupied_states(st2, ik, i2, j2, k2, filled2, partially_filled2, half_filled2)
      if( (j1 > 0) .or. (j2 > 0) ) then
        message(1) = 'Cannot calculate many-body dot products with partially occupied orbitals'
        call messages_fatal(1)
      end if
      if(i1 .ne. i2) then
        message(1) = 'Internal Error: different number of occupied states in states_mpdotp'
        call messages_fatal(1)
      end if

      if(i1 > 0) then
        SAFE_ALLOCATE(bb(1:i1, 1:i1))
        do ii = 1, i1
          do jj = 1, i1
            bb(ii, jj) = aa(filled1(ii), filled2(jj), ik)
          end do
        end do

        dotp = dotp * lalg_determinant(i1, bb, invert = .false.) ** st1%d%kweights(ik)
        SAFE_DEALLOCATE_A(bb)
      end if

    end do
  end select

  SAFE_DEALLOCATE_A(aa)
  SAFE_DEALLOCATE_A(filled1)
  SAFE_DEALLOCATE_A(filled2)
  SAFE_DEALLOCATE_A(partially_filled1)
  SAFE_DEALLOCATE_A(partially_filled2)
  SAFE_DEALLOCATE_A(half_filled1)
  SAFE_DEALLOCATE_A(half_filled2)
  POP_SUB(X(states_mpdotp_g))
end function X(states_mpdotp_g)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
