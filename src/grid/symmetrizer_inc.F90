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
!! $Id: symmetrizer_inc.F90 5918 2009-09-25 21:43:56Z dstrubbe $

subroutine X(symmetrizer_apply)(this, field, symmfield)
  type(symmetrizer_t), intent(in)    :: this
  R_TYPE, target,      intent(in)    :: field(:)
  R_TYPE,              intent(out)   :: symmfield(:)

  integer :: ip, iop, nops, ipsrc, idir
  integer :: destpoint(1:3), srcpoint(1:3), lsize(1:3), offset(1:3)
  R_TYPE  :: acc
  FLOAT   :: weight
  R_TYPE, pointer :: field_global(:)
  type(pv_t), pointer :: vp

  PUSH_SUB(X(symmetrizer_apply))

  vp => this%mesh%vp

  if(this%mesh%parallel_in_domains) then

    ! With domain parallelization, we collect all points of the
    ! 'field' array. This seems reasonable, since we will probably
    ! need almost all points anyway
    !
    ! The symmfield array is kept locally.

    SAFE_ALLOCATE(field_global(1:this%mesh%np_global))
#ifdef HAVE_MPI
    call X(vec_allgather)(vp, field_global, field)
#endif
  else
    field_global => field
  end if
  
  nops = symmetries_number(this%mesh%sb%symm)
  weight = M_ONE/nops

  lsize(1:3) = this%mesh%idx%ll(1:3)
  offset(1:3) = this%mesh%idx%nr(1, 1:3) + this%mesh%idx%enlarge(1:3)

  do ip = 1, this%mesh%np

    if(this%mesh%parallel_in_domains) then
      ! convert to global point
      destpoint(1:3) = this%mesh%idx%lxyz(vp%local(vp%xlocal(vp%partno) + ip - 1), 1:3) - offset(1:3)
    else
      destpoint(1:3) = this%mesh%idx%lxyz(ip, 1:3) - offset(1:3)
    end if

    ASSERT(all(destpoint < lsize))

    acc = M_ZERO

    ! iterate over all points that go to this point by a symmetry operation
    do iop = 1, nops
      srcpoint = symm_op_apply_inv(this%mesh%sb%symm%ops(iop), destpoint)
      do idir = 1, 3
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

      acc = acc + field_global(ipsrc)
    end do

    symmfield(ip) = weight*acc

  end do

  if(this%mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(field_global)
  end if

  POP_SUB(X(symmetrizer_apply))
end subroutine X(symmetrizer_apply)

! ----------------------------------------------------------------------------------

subroutine X(symmetrizer_apply_vector)(this, field, symmfield)
  type(symmetrizer_t), intent(in)    :: this
  R_TYPE, target,      intent(in)    :: field(:, :)
  R_TYPE,              intent(out)   :: symmfield(:, :)

  integer :: ip, iop, nops, ipsrc, idir
  integer :: destpoint(1:3), srcpoint(1:3), lsize(1:3), offset(1:3)
  FLOAT   :: weight
  R_TYPE, pointer :: field_global(:, :)
  type(pv_t), pointer :: vp

  PUSH_SUB(X(symmetrizer_apply_vector))

  nops = symmetries_number(this%mesh%sb%symm)
  weight = M_ONE/nops

  vp => this%mesh%vp

  if(this%mesh%parallel_in_domains) then
    ! When running with domain parallelization, we collect all points
    ! of the field array. The symmfield array is kept locally.

    SAFE_ALLOCATE(field_global(1:this%mesh%np_global, 1:3))
    do idir = 1, 3
#ifdef HAVE_MPI
      call X(vec_allgather)(vp, field_global(:, idir), field(:, idir))
#endif
    end do
  else
    field_global => field
  end if

  lsize(1:3) = this%mesh%idx%ll(1:3)
  offset(1:3) = this%mesh%idx%nr(1, 1:3) + this%mesh%idx%enlarge(1:3)

  do ip = 1, this%mesh%np
    if(this%mesh%parallel_in_domains) then
      ! convert to global point
      destpoint(1:3) = this%mesh%idx%lxyz(vp%local(vp%xlocal(vp%partno) + ip - 1), 1:3) - offset(1:3)
    else
      destpoint(1:3) = this%mesh%idx%lxyz(ip, 1:3) - offset(1:3)
    end if

    ASSERT(all(destpoint < lsize))

    symmfield(ip, 1:3) = M_ZERO

    ! iterate over all points that go to this point by a symmetry operation
    do iop = 1, nops
      srcpoint = symm_op_apply_inv(this%mesh%sb%symm%ops(iop), destpoint)
      do idir = 1, 3
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

      symmfield(ip, 1:3) = symmfield(ip, 1:3) + weight*symm_op_apply(this%mesh%sb%symm%ops(iop), field_global(ipsrc, 1:3))

    end do

  end do

  if(this%mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(field_global)
  end if

  POP_SUB(X(symmetrizer_apply_vector))
end subroutine X(symmetrizer_apply_vector)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
