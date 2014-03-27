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
!! $Id$

!> supply field and symmfield, and/or field_vector and symmfield_vector
subroutine X(symmetrizer_apply)(this, field, field_vector, symmfield, symmfield_vector, suppress_warning)
  type(symmetrizer_t), target, intent(in)    :: this
  R_TYPE,    optional, target, intent(in)    :: field(:) !< (this%mesh%np)
  R_TYPE,    optional, target, intent(in)    :: field_vector(:, :)  !< (this%mesh%np, 3)
  R_TYPE,            optional, intent(out)   :: symmfield(:) !< (this%mesh%np)
  R_TYPE,            optional, intent(out)   :: symmfield_vector(:, :) !< (this%mesh%np, 3)
  logical,           optional, intent(in)    :: suppress_warning !< use to avoid output of discrepancy,
    !! for forces, where this routine is not used to symmetrize something already supposed to be symmetric,
    !! but rather to construct the quantity properly from reduced k-points

  integer :: ip, iop, nops, ipsrc, idir
  integer :: destpoint(1:3), srcpoint(1:3), lsize(1:3), offset(1:3)
  R_TYPE  :: acc, acc_vector(1:3)
  FLOAT   :: weight, maxabsdiff, maxabs
  R_TYPE, pointer :: field_global(:), field_global_vector(:, :)
  type(pv_t), pointer :: vp

  PUSH_SUB(X(symmetrizer_apply))

  ASSERT(present(field) .eqv. present(symmfield))
  ASSERT(present(field_vector) .eqv. present(symmfield_vector))
  ! we will do nothing if following condition is not met!
  ASSERT(present(field) .or. present(field_vector))

  ASSERT(associated(this%mesh))
  vp => this%mesh%vp

  ! With domain parallelization, we collect all points of the
  ! 'field' array. This seems reasonable, since we will probably
  ! need almost all points anyway.
  !
  ! The symmfield array is kept locally.

  if(present(field)) then
    if(this%mesh%parallel_in_domains) then
      SAFE_ALLOCATE(field_global(1:this%mesh%np_global))
#ifdef HAVE_MPI
      call X(vec_allgather)(vp, field_global, field)
#endif
    else
      field_global => field
    end if
  endif

  if(present(field_vector)) then
    if(this%mesh%parallel_in_domains) then
      SAFE_ALLOCATE(field_global_vector(1:this%mesh%np_global, 1:3))
      do idir = 1, 3
#ifdef HAVE_MPI
        call X(vec_allgather)(vp, field_global_vector(:, idir), field_vector(:, idir))
#endif
      enddo
    else
      field_global_vector => field_vector
    end if
  endif

  nops = symmetries_number(this%mesh%sb%symm)
  weight = M_ONE/nops

  lsize(1:3) = this%mesh%idx%ll(1:3)
  offset(1:3) = this%mesh%idx%nr(1, 1:3) + this%mesh%idx%enlarge(1:3)

  do ip = 1, this%mesh%np
    if(this%mesh%parallel_in_domains) then
      ! convert to global point
      destpoint(1:3) = this%mesh%idx%lxyz(vp%local(vp%xlocal + ip - 1), 1:3) - offset(1:3)
    else
      destpoint(1:3) = this%mesh%idx%lxyz(ip, 1:3) - offset(1:3)
    end if
    ! offset moves corner of cell to origin, in integer mesh coordinates

    ASSERT(all(destpoint >= 0))
    ASSERT(all(destpoint < lsize))

    ! move to center of cell in real coordinates
    destpoint = destpoint - lsize / 2

    if(present(field)) &
      acc = M_ZERO
    if(present(field_vector)) &
      acc_vector(1:3) = M_ZERO

    ! iterate over all points that go to this point by a symmetry operation
    do iop = 1, nops
      srcpoint = symm_op_apply_inv(this%mesh%sb%symm%ops(iop), destpoint)

      ! move back to reference to origin at corner of cell
      srcpoint = srcpoint + lsize / 2
      ! apply periodic boundary conditions in periodic directions
      do idir = 1, this%mesh%sb%periodic_dim
        if(srcpoint(idir) < 0) then
          srcpoint(idir) = srcpoint(idir) + lsize(idir)
          ASSERT(srcpoint(idir) >= 0)
        else if(srcpoint(idir) >= lsize(idir)) then
          srcpoint(idir) = mod(srcpoint(idir), lsize(idir))
        end if
      end do
      ASSERT(all(srcpoint >= 0))
      ASSERT(all(srcpoint < lsize))
      srcpoint(1:3) = srcpoint(1:3) + offset(1:3)

      ipsrc = this%mesh%idx%lxyz_inv(srcpoint(1), srcpoint(2), srcpoint(3))

      if(present(field)) then
        acc = acc + field_global(ipsrc)
      endif
      if(present(field_vector)) then
        acc_vector(1:3) = acc_vector(1:3) + symm_op_apply(this%mesh%sb%symm%ops(iop), field_global_vector(ipsrc, 1:3))
      endif
    end do
    if(present(field)) &
      symmfield(ip) = weight * acc
    if(present(field_vector)) &
      symmfield_vector(ip, 1:3) = weight * acc_vector(1:3)

  end do

  if(.not. optional_default(suppress_warning, .false.)) then
    if(present(field)) then
      maxabs = maxval(abs(field(1:this%mesh%np)))
      maxabsdiff = maxval(abs(field(1:this%mesh%np) - symmfield(1:this%mesh%np)))
      if(maxabsdiff / maxabs > CNST(1e-6)) then
        write(message(1),'(a, es12.6)') 'Symmetrization discrepancy ratio (scalar) = ', maxabsdiff / maxabs
        call messages_warning(1)
      endif
    endif
    
    if(present(field_vector)) then
      maxabs = maxval(abs(field_vector(1:this%mesh%np, 1:3)))
      maxabsdiff = maxval(abs(field_vector(1:this%mesh%np, 1:3) - symmfield_vector(1:this%mesh%np, 1:3)))
      if(maxabsdiff / maxabs > CNST(1e-6)) then
        write(message(1),'(a, es12.6)') 'Symmetrization discrepancy ratio (vector) = ', maxabsdiff / maxabs
        call messages_warning(1)
      endif
    endif
  endif

  if(this%mesh%parallel_in_domains) then
    if(present(field)) then
      SAFE_DEALLOCATE_P(field_global)
    endif
    if(present(field_vector)) then
      SAFE_DEALLOCATE_P(field_global_vector)
    endif
  end if

  POP_SUB(X(symmetrizer_apply))
end subroutine X(symmetrizer_apply)

! -------------------------------------------------------------------------------
subroutine X(symmetrize_tensor)(symm, tensor)
  type(symmetries_t), intent(in)    :: symm
  R_TYPE,             intent(inout) :: tensor(:,:) !< (3, 3)
  
  integer :: iop, nops, idir
  R_TYPE :: tensor_symm(3, 3)
  
  PUSH_SUB(X(symmetrize_tensor))

  nops = symmetries_number(symm)
  
  tensor_symm(:,:) = M_ZERO
  do iop = 1, nops
    tensor_symm(:,:) = tensor_symm(:,:) + &
      matmul(matmul(dble(transpose(symm_op_rotation_matrix(symm%ops(iop)))), tensor(1:3, 1:3)), &
      dble(symm_op_rotation_matrix(symm%ops(iop))))
  enddo

  tensor(:,:) = tensor_symm(:,:) / nops

  POP_SUB(X(symmetrize_tensor))
end subroutine X(symmetrize_tensor)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
