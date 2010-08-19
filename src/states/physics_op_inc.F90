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
!! $Id: derivatives.F90 5812 2009-08-13 10:18:31Z marques $

! ---------------------------------------------------------
! The action of the angular momentum operator (three spatial components).
! In case of real functions, it does not include the -i prefactor
! (L = -i r ^ nabla).
! ---------------------------------------------------------
subroutine X(physics_op_L)(der, ff, lf, ghost_update, set_bc)
  type(derivatives_t), intent(inout) :: der
  R_TYPE,              intent(inout) :: ff(:)    ! ff(der%mesh%np_part)
  R_TYPE,              intent(out)   :: lf(:, :) ! lf(der%mesh%np, 3) in 3D, lf(der%mesh%np, 1) in 2D
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x1, x2, x3
  R_TYPE :: factor
  integer :: ip

  PUSH_SUB(X(physics_op_L))

  ASSERT(der%mesh%sb%dim .ne. 1)

  SAFE_ALLOCATE(gf(1:der%mesh%np, 1:der%mesh%sb%dim))

  call X(derivatives_grad)(der, ff, gf, ghost_update, set_bc)

#if defined(R_TCOMPLEX)
  factor = -M_ZI
#else
  factor = M_ONE
#endif

  select case(der%mesh%sb%dim)
  case(3)
    do ip = 1, der%mesh%np
      x1 = der%mesh%x(ip, 1)
      x2 = der%mesh%x(ip, 2)
      x3 = der%mesh%x(ip, 3)
      lf(ip, 1) = factor * (x2 * gf(ip, 3) - x3 * gf(ip, 2))
      lf(ip, 2) = factor * (x3 * gf(ip, 1) - x1 * gf(ip, 3))
      lf(ip, 3) = factor * (x1 * gf(ip, 2) - x2 * gf(ip, 1))
    end do

  case(2)
    do ip = 1, der%mesh%np
      x1 = der%mesh%x(ip, 1)
      x2 = der%mesh%x(ip, 2)
      lf(ip, 1) = factor * (x1 * gf(ip, 2) - x2 * gf(ip, 1))
    end do

  end select

  SAFE_DEALLOCATE_A(gf)
  POP_SUB(X(physics_op_L))
end subroutine X(physics_op_L)


! ---------------------------------------------------------
! Square of the angular momentum L. This has to be very much improved if
! accuracy is needed.
! ---------------------------------------------------------
subroutine X(physics_op_L2)(der, ff, l2f, ghost_update, set_bc)
  type(derivatives_t), intent(inout) :: der
  R_TYPE,              intent(inout) :: ff(:)   ! ff(1:der%mesh%np_part)
  R_TYPE,              intent(out)   :: l2f(:)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  R_TYPE, allocatable :: gf(:, :), ggf(:, :, :)
  integer :: idir

  PUSH_SUB(X(physics_op_L2))

  ASSERT(der%mesh%sb%dim == 2 .or. der%mesh%sb%dim == 3)

  l2f = R_TOTYPE(M_ZERO)

  select case(der%mesh%sb%dim)
  case(3)
    SAFE_ALLOCATE( gf(1:der%mesh%np_part, 1:3))
    SAFE_ALLOCATE(ggf(1:der%mesh%np_part, 1:3, 1:3))

    call X(physics_op_L)(der, ff, gf, ghost_update, set_bc)

    do idir = 1, 3
      call X(physics_op_L)(der, gf(:, idir), ggf(:, :, idir))
    end do

    do idir = 1, 3
      l2f(1:der%mesh%np) = l2f(1:der%mesh%np) + ggf(1:der%mesh%np, idir, idir)
    end do

  case(2)
    SAFE_ALLOCATE( gf(1:der%mesh%np_part, 1:1))
    SAFE_ALLOCATE(ggf(1:der%mesh%np_part, 1:1, 1:1))

    call X(physics_op_L)(der, ff, gf, ghost_update, set_bc)
    call X(physics_op_L)(der, gf(:, 1), ggf(:, :, 1))

    l2f(1:der%mesh%np) = ggf(1:der%mesh%np, 1, 1)
  end select


  ! In case of real functions, since the angular-momentum calculations
  ! lack a (-i) prefactor, we must add a (-1) factor
#if defined(R_TREAL)
  l2f = - l2f
#endif

  SAFE_DEALLOCATE_A(gf)
  SAFE_DEALLOCATE_A(ggf)
  POP_SUB(X(physics_op_L2))
end subroutine X(physics_op_L2)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
