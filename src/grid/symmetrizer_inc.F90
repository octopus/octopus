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

!> supply field and symmfield, and/or field_vector and symmfield_vector
subroutine X(symmetrizer_apply)(this, np, field, field_vector, symmfield, symmfield_vector, &
          suppress_warning, reduced_quantity)
  type(symmetrizer_t),         intent(in)    :: this
  integer,                     intent(in)    :: np !mesh%np or mesh%fine%np
  R_TYPE,    optional, target, intent(in)    :: field(:) !< (np)
  R_TYPE,    optional, target, intent(in)    :: field_vector(:, :)  !< (np, 3)
  R_TYPE,            optional, intent(out)   :: symmfield(:) !< (np)
  R_TYPE,            optional, intent(out)   :: symmfield_vector(:, :) !< (np, 3)
  logical,           optional, intent(in)    :: suppress_warning !< use to avoid output of discrepancy,
    !! for forces, where this routine is not used to symmetrize something already supposed to be symmetric,
    !! but rather to construct the quantity properly from reduced k-points
  logical,           optional, intent(in)    :: reduced_quantity

  integer :: ip, iop, nops, ipsrc, idir
  R_TYPE  :: acc, acc_vector(1:3)
  FLOAT   :: weight, maxabsdiff, maxabs
  R_TYPE, pointer :: field_global(:), field_global_vector(:, :)

  PUSH_SUB(X(symmetrizer_apply))

  ASSERT(present(field) .eqv. present(symmfield))
  ASSERT(present(field_vector) .eqv. present(symmfield_vector))
  ! we will do nothing if following condition is not met!
  ASSERT(present(field) .or. present(field_vector))

  if(present(field)) then
    ASSERT(ubound(field, dim = 1) >= np)
    ASSERT(ubound(symmfield, dim = 1) >= np)
  else
    ASSERT(ubound(field_vector, dim = 1) >= np)
    ASSERT(ubound(symmfield_vector, dim = 1) >= np)
  end if

  ASSERT(associated(this%mesh))

  !If the quantity to be symmetrized is uniformly zero, we don`t symmetrize
  if(present(field)) then
    maxabs = maxval(abs(field(1:np)))
    if(maxabs < M_EPSILON) then
      symmfield(1:np) = field(1:np)
      POP_SUB(X(symmetrizer_apply))
      return
    end if
  end if
  if(present(field_vector)) then
    maxabs = maxval(abs(field_vector(1:np, 1:3)))
    if(maxabs < M_EPSILON) then
      symmfield_vector(1:np, 1:3) = field_vector(1:np, 1:3)
      POP_SUB(X(symmetrizer_apply))
      return
    end if
  end if

  ! With domain parallelization, we collect all points of the
  ! 'field' array. This seems reasonable, since we will probably
  ! need almost all points anyway.
  !
  ! The symmfield array is kept locally.

  if(present(field)) then
    if(this%mesh%parallel_in_domains) then
      SAFE_ALLOCATE(field_global(1:this%mesh%np_global))
      call vec_allgather(this%mesh%vp, field_global, field)
    else
      field_global => field
    end if
  end if

  if(present(field_vector)) then
    ASSERT(ubound(field_vector, dim=2) == 3)
    if(this%mesh%parallel_in_domains) then
      SAFE_ALLOCATE(field_global_vector(1:this%mesh%np_global, 1:3))
      do idir = 1, 3
        call vec_allgather(this%mesh%vp, field_global_vector(:, idir), field_vector(:, idir))
      end do
    else
      field_global_vector => field_vector
    end if
  end if

  nops = symmetries_number(this%mesh%sb%symm)
  weight = M_ONE/nops

  do ip = 1, np
    if(present(field)) acc = M_ZERO
    if(present(field_vector)) acc_vector(1:3) = M_ZERO

    ! iterate over all points that go to this point by a symmetry operation
    do iop = 1, nops
      ipsrc = this%map(ip, iop)

      if(present(field)) then
        acc = acc + field_global(ipsrc)
      end if
      if(present(field_vector)) then
        if(.not.optional_default(reduced_quantity, .false.)) then
          acc_vector(1:3) = acc_vector(1:3) + symm_op_apply_inv_cart(this%mesh%sb%symm%ops(iop), field_global_vector(ipsrc, 1:3))
        else
          acc_vector(1:3) = acc_vector(1:3) + symm_op_apply_inv_red(this%mesh%sb%symm%ops(iop), field_global_vector(ipsrc, 1:3))
        end if
      end if
    end do
    if(present(field)) &
      symmfield(ip) = weight * acc
    if(present(field_vector)) &
      symmfield_vector(ip, 1:3) = weight * acc_vector(1:3)

  end do

  if(.not. optional_default(suppress_warning, .false.)) then
    if(present(field)) then
      maxabsdiff = maxval(abs(field(1:np) - symmfield(1:np)))
      if(maxabsdiff / maxabs > CNST(1e-6)) then
        write(message(1),'(a, es12.5)') 'Symmetrization discrepancy ratio (scalar) = ', maxabsdiff / maxabs
        call messages_warning(1)
      end if
    end if
    
    if(present(field_vector)) then
      maxabsdiff = maxval(abs(field_vector(1:np, 1:3) - symmfield_vector(1:np, 1:3)))
      if(maxabsdiff / maxabs > CNST(1e-6)) then
        write(message(1),'(a, es12.5)') 'Symmetrization discrepancy ratio (vector) = ', maxabsdiff / maxabs
        call messages_warning(1)
      end if
    end if
  end if

  if(this%mesh%parallel_in_domains) then
    if(present(field)) then
      SAFE_DEALLOCATE_P(field_global)
    end if
    if(present(field_vector)) then
      SAFE_DEALLOCATE_P(field_global_vector)
    end if
  end if

  POP_SUB(X(symmetrizer_apply))
end subroutine X(symmetrizer_apply)

!The same as for symmetrizer_apply, but a single symmetry operation
!Here iop can be negative, indicating the spatial symmetry plus time reversal symmetry
subroutine X(symmetrizer_apply_single)(this, np, iop, field, symmfield)
  type(symmetrizer_t),         intent(in)    :: this
  integer,                     intent(in)    :: np !mesh%np or mesh%fine%np
  integer,                     intent(in)    :: iop
  R_TYPE,              target, intent(in)    :: field(:) !< (np)
  R_TYPE,                      intent(out)   :: symmfield(:) !< (np)

  integer :: ip
  R_TYPE, pointer :: field_global(:)
  type(profile_t), save :: prof

  PUSH_SUB(X(symmetrizer_apply_single))

  call profiling_in(prof, TOSTRING(X(SYMMETRIZE_SINGLE)))

  ASSERT(ubound(field, dim = 1) >= np)
  ASSERT(ubound(symmfield, dim = 1) >= np)
  ASSERT(associated(this%mesh))

  ! With domain parallelization, we collect all points of the
  ! 'field' array. This seems reasonable, since we will probably
  ! need almost all points anyway.
  !
  ! The symmfield array is kept locally.

  if(this%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(field_global(1:this%mesh%np_global))
    call vec_allgather(this%mesh%vp, field_global, field)
  else
    field_global => field
  end if

  if(iop>0) then
    do ip = 1, np
       symmfield(ip) = field_global(this%map(ip,abs(iop)))
    end do
  else
    do ip = 1, np
       symmfield(ip) = R_CONJ(field_global(this%map(ip,abs(iop))))
    end do
  end if

  if(this%mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(field_global)
  end if

  call profiling_out(prof)

  POP_SUB(X(symmetrizer_apply_single))
end subroutine X(symmetrizer_apply_single)

! -------------------------------------------------------------------------------
subroutine X(symmetrize_tensor_cart)(symm, tensor)
  type(symmetries_t), intent(in)    :: symm
  R_TYPE,             intent(inout) :: tensor(:,:) !< (3, 3)
  
  integer :: iop, nops
  R_TYPE :: tensor_symm(3, 3)
  FLOAT :: tmp(3,3)
  
  PUSH_SUB(X(symmetrize_tensor_cart))

  nops = symmetries_number(symm)
  
  tensor_symm(:,:) = M_ZERO
  do iop = 1, nops
    ! The use of the tmp array is a workaround for a PGI bug
    tmp = symm_op_rotation_matrix_cart(symm%ops(iop))
    tensor_symm(1:3,1:3) = tensor_symm(1:3,1:3) + &
      matmul(matmul(transpose(symm_op_rotation_matrix_cart(symm%ops(iop))), tensor(1:3, 1:3)), tmp)
  end do

  tensor(1:3,1:3) = tensor_symm(1:3,1:3) / nops

  POP_SUB(X(symmetrize_tensor_cart))
end subroutine X(symmetrize_tensor_cart)

! -------------------------------------------------------------------------------
! The magneto-optical response 
! j_\nu = \alpha_{\nu \mu, \gamma} B_\gamma E_\mu has the symmetry of
! e_{\nu \mu \gamma} e_{\alpha \beta \gamma} x_\nu x_\mu x_\alpha x_\beta. 
! The response should not change upon changing signs of any direction(s) (x -> -x).
! However, it should change sign upon permutation in the order of axes (x,y -> y,x).
! Therefore, contribution from a rotation symmetry should be multiplied by the
! determinant of the rotation matrix ignoring signs of the rotation matrix elements.
subroutine X(symmetrize_magneto_optics_cart)(symm, tensor) 
  type(symmetries_t),   intent(in)    :: symm 
  R_TYPE,               intent(inout) :: tensor(:,:,:) !< (3, 3, 3)
  
  integer :: iop, nops
  R_TYPE  :: tensor_symm(3, 3, 3)
  integer :: rot(3, 3)
  integer :: idir1, idir2, idir3, ndir
  integer :: i1, i2, i3, det
  
  PUSH_SUB(X(symmetrize_magneto_optics_cart))
    
  ndir = 3
  
  nops = symmetries_number(symm)
  
  tensor_symm(:,:,:) = M_ZERO
  
  do iop = 1, nops
    rot = symm_op_rotation_matrix_red(symm%ops(iop))
    det = abs(rot(1,1) * rot(2,2) * rot(3,3)) + abs(rot(1,2) * rot(2,3) * rot(3,1)) &
      + abs(rot(1,3) * rot(2,1) * rot(3,2)) - abs(rot(1,1) * rot(2,3) * rot(3,2)) &
      - abs(rot(1,2) * rot(2,1) * rot(3,3)) - abs(rot(1,3) * rot(2,2) * rot(3,1))
          
    do idir1 = 1, ndir
      do idir2 = 1, ndir
        do idir3 = 1, ndir
          ! change of indices upon rotation
          i1 = abs(1 * rot(1,idir1) + 2 * rot(2,idir1) + 3 * rot(3,idir1))
          i2 = abs(1 * rot(1,idir2) + 2 * rot(2,idir2) + 3 * rot(3,idir2))
          i3 = abs(1 * rot(1,idir3) + 2 * rot(2,idir3) + 3 * rot(3,idir3))
          
          tensor_symm(i1,i2,i3) = tensor_symm(i1,i2,i3) + &
            tensor(idir1,idir2,idir3) * det
        end do
      end do
    end do
  end do
  tensor(1:3,1:3,1:3) = tensor_symm(1:3,1:3,1:3) / nops

  POP_SUB(X(symmetrize_magneto_optics_cart))
end subroutine X(symmetrize_magneto_optics_cart)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
