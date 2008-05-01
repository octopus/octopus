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
! The next two subroutines convert a function between the normal
! mesh and the cube.
! Note that the function in the mesh should be defined
! globally, not just in a partition (when running in
! parallel in real-space domains).
! ---------------------------------------------------------
subroutine X(mf2cf) (m, mf, cf)
  type(mesh_t),  intent(in)    :: m
  R_TYPE,        intent(in)    :: mf(:)  ! mf(m%np_global)
  type(X(cf_t)), intent(inout) :: cf

  integer :: i, ix, iy, iz, c(MAX_DIM)

  ASSERT(associated(cf%RS))

  c(:)  =  cf%n(:)/2 + 1

  !$omp parallel workshare
  cf%RS =  M_ZERO
  !$omp end parallel workshare

  !$omp parallel do private(ix, iy, iz)
  do i = 1, m%np_global
    ix = m%Lxyz(i, 1) + c(1)
    iy = m%Lxyz(i, 2) + c(2)
    iz = m%Lxyz(i, 3) + c(3)

    cf%RS(ix, iy, iz) = mf(i)
  end do
  !$omp end parallel do

end subroutine X(mf2cf)


! ---------------------------------------------------------
subroutine X(cf2mf) (m, cf, mf)
  type(mesh_t),  intent(in)  :: m
  type(X(cf_t)), intent(in)  :: cf
  R_TYPE,        intent(out) :: mf(:)  ! mf(m%np_global)

  integer :: i, ix, iy, iz, c(MAX_DIM)

  ASSERT(associated(cf%RS))

  c(:) =  cf%n(:)/2 + 1

  !$omp parallel do private(ix, iy, iz)
  do i = 1, m%np_global
    ix = m%Lxyz(i, 1) + c(1)
    iy = m%Lxyz(i, 2) + c(2)
    iz = m%Lxyz(i, 3) + c(3)
    mf(i) = cf%RS(ix, iy, iz)
  end do
  !$omp end parallel do

end subroutine X(cf2mf)

! ---------------------------------------------------------
! The next two subroutines convert a function in Fourier space
! between the normal mesh and the cube
! Note that the function in the mesh should be defined
! globally, not just in a partition (when running in
! parallel in real-space domains).
! ---------------------------------------------------------
subroutine X(mf2cf_FS) (m, mf, cf)
  type(mesh_t),  intent(in)    :: m
  CMPLX,         intent(in)    :: mf(:)   ! mf(m%np_global)
  type(X(cf_t)), intent(inout) :: cf

  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  cf%FS = M_z0

  do i = 1, m%np_global
    ix = pad_feq(m%Lxyz(i, 1), cf%n(1), .false.)
    if(ix > cf%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(m%Lxyz(i, 2), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(i, 3), cf%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(i)
  end do
end subroutine X(mf2cf_FS)


! ---------------------------------------------------------
subroutine X(cf_FS2mf) (m, cf, mf)
  type(mesh_t),  intent(in)  :: m
  type(X(cf_t)), intent(in)  :: cf
  CMPLX,         intent(out) :: mf(:) ! mf(m%np_global)

  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  do i = 1, m%np_global
    ix = pad_feq(m%Lxyz(i, 1), cf%n(1), .false.)
    iy = pad_feq(m%Lxyz(i, 2), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(i, 3), cf%n(3), .false.)

#ifdef R_TREAL
    if(ix > cf%nx) then
      ix = pad_feq(-m%Lxyz(i, 1), cf%n(1), .false.)
      mf(i) = conjg(cf%FS(ix, iy, iz))
    else
      mf(i) = cf%FS(ix, iy, iz)
    end if
#else
    mf(i) = cf%FS(ix, iy, iz)
#endif
  end do

end subroutine X(cf_FS2mf)

! ---------------------------------------------------------
! Calculation of derivatives
! ---------------------------------------------------------
subroutine X(f_laplacian) (sb, f_der, f, lapl, cutoff_, ghost_update, set_bc)
  type(simul_box_t), intent(in)    :: sb
  type(f_der_t),     intent(inout) :: f_der
  R_TYPE,            intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,            intent(out)   :: lapl(:)  ! lapl(m%np_part)
  FLOAT, optional,   intent(in)    :: cutoff_
  logical, optional, intent(in)    :: ghost_update
  logical, optional, intent(in)    :: set_bc

  call push_sub('f_inc.Xf_laplacian')

  call X(derivatives_lapl) (f_der%der_discr, f, lapl, ghost_update, set_bc)

  call pop_sub()
end subroutine X(f_laplacian)


! ---------------------------------------------------------
subroutine X(f_gradient) (sb, f_der, f, grad, ghost_update)
  type(simul_box_t), intent(in)    :: sb
  type(f_der_t),     intent(inout) :: f_der
  R_TYPE,            intent(inout) :: f(:)         ! f(m%np_part)
  R_TYPE,            intent(out)   :: grad(:,:)    ! grad(m%np, m%sb%dim)
  logical, optional, intent(in)    :: ghost_update

  call push_sub('f_inc.Xf_gradient')

  call X(derivatives_grad) (f_der%der_discr, f, grad, ghost_update)

  call pop_sub()
end subroutine X(f_gradient)


! ---------------------------------------------------------
subroutine X(f_divergence) (sb, f_der, f, divf, ghost_update)
  type(simul_box_t), intent(in)    :: sb
  type(f_der_t),     intent(inout) :: f_der
  R_TYPE,            intent(inout) :: f(:,:)    ! f(m%np_part, sb%dim)
  R_TYPE,            intent(out)   :: divf(:)   ! divf(m%np)
  logical, optional, intent(in)    :: ghost_update

  call push_sub('f_inc.Xf_divergence')

  call X(derivatives_div) (f_der%der_discr, f, divf, ghost_update)

  call pop_sub()
end subroutine X(f_divergence)


! ---------------------------------------------------------
subroutine X(f_curl) (f_der, f, curlf, ghost_update)
  type(f_der_t), intent(inout)  :: f_der
  R_TYPE,        intent(inout)  :: f(:,:)     ! f(m%np_part, conf%dim)
  R_TYPE,        intent(out)    :: curlf(:,:) ! curlf(m%np, conf%dim) if dim = 3, curlf(m%np, 1) if dim = 2 
  logical, optional, intent(in) :: ghost_update

  call push_sub('f_inc.Xf_curl')

  call X(derivatives_curl) (f_der%der_discr, f, curlf)

  call pop_sub()
end subroutine X(f_curl)


! -----------------------------------------------------------------------------
! This routine calculates the multipoles of a function f
! distribution, defined in the following way:
! multipole(1) is the trace of f (defined to be positive; integral
!   of the f.
! multipole(2:4) contains the dipole: integral of f times x, y or z.
! multipole(5:9, is) contains the quadrupole, defined in the usual way using
!   the spherical harmonics: multipole(5) = Integral [ f * Y_{2,-2} ],
!   multipole(6, is) = Integral [ f * Y_{2, -1} ].
! And so on.
! -----------------------------------------------------------------------------
subroutine X(f_multipoles) (mesh, ff, lmax, multipole)
  type(mesh_t),   intent(in)  :: mesh
  R_TYPE,         intent(in)  :: ff(:)
  integer,        intent(in)  :: lmax
  R_TYPE,         intent(out) :: multipole(:) ! multipole((lmax + 1)**2)

  integer :: i, l, lm, add_lm
  FLOAT   :: x(MAX_DIM), r, ylm
  R_TYPE, allocatable :: ff2(:)

  call push_sub('states.states_calculate_multipoles')

  ALLOCATE(ff2(mesh%np), mesh%np)

  ff2(1:mesh%np) = ff(1:mesh%np)
  multipole(1) = X(mf_integrate)(mesh, ff2)

  if(lmax > 0) then
    do i = 1, 3
      ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, i)
      multipole(i+1) = X(mf_integrate)(mesh, ff2)
    end do
  end if

  if(lmax>1) then
    add_lm = 5
    do l = 2, lmax
      do lm = -l, l
        do i = 1, mesh%np
          call mesh_r(mesh, i, r, x=x)
          ylm = loct_ylm(x(1), x(2), x(3), l, lm)
          ff2(i) = ff(i) * ylm * r**l
        end do
        multipole(add_lm) = X(mf_integrate)(mesh, ff2)
        add_lm = add_lm + 1
      end do
    end do
  end if

  deallocate(ff2)
  call pop_sub()
end subroutine X(f_multipoles)


! ---------------------------------------------------------
! The action of the angular momentum operator (three spatial components).
! In case of real functions, it does not include the -i prefactor
! (L = -i r ^ nabla).
! ---------------------------------------------------------
subroutine X(f_angular_momentum)(sb, f_der, f, lf, ghost_update)
  type(simul_box_t), intent(in)    :: sb
  type(f_der_t),     intent(inout) :: f_der
  R_TYPE,            intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,            intent(out)   :: lf(:,:)  ! lf(m%np, 3) in 3D, lf(m%np, 1) in 2D
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x(sb%dim)
  integer :: i

  call push_sub('f_inc.Xf_angular_momentum')

  ASSERT(sb%dim.ne.1)

  ALLOCATE(gf(f_der%m%np, sb%dim), f_der%m%np*sb%dim)
  call X(f_gradient)(sb, f_der, f, gf, ghost_update)

  do i = 1, f_der%m%np
    x(1:sb%dim) = f_der%m%x(i,1:sb%dim)
    select case(sb%dim)
    case(3)
      lf(i, 1) = (x(2)*gf(i, 3) - x(3)*gf(i, 2))
      lf(i, 2) = (x(3)*gf(i, 1) - x(1)*gf(i, 3))
      lf(i, 3) = (x(1)*gf(i, 2) - x(2)*gf(i, 1))
    case(2)
      lf(i, 1) = (x(1)*gf(i, 2) - x(2)*gf(i, 1))
    end select
  end do
#if defined(R_TCOMPLEX)
  select case(sb%dim)
  case(3)
    do i = 1, 3
      call lalg_scal(f_der%m%np, -M_zI, lf(:, i))
    end do
  case(2)
    call lalg_scal(f_der%m%np, -M_zI, lf(:, 1))
  end select
#endif

  deallocate(gf)
  call pop_sub()
end subroutine X(f_angular_momentum)


! ---------------------------------------------------------
! Square of the angular momentum L. This has to be very much improved if
! accuracy is needed.
! ---------------------------------------------------------
subroutine X(f_l2)(sb, f_der, f, l2f, ghost_update)
  type(simul_box_t), intent(in)    :: sb
  type(f_der_t),     intent(inout) :: f_der
  R_TYPE,            intent(inout) :: f(:)   ! f(1:m%np_part)
  R_TYPE,            intent(out)   :: l2f(:)
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: gf(:, :), ggf(:, :, :)
  type(mesh_t), pointer :: m
  integer :: j

  call push_sub('f_inc.Xf_l2')

  ASSERT(sb%dim.ne.1)

  m => f_der%m
  l2f = R_TOTYPE(M_ZERO)

  select case(sb%dim)
  case(3)
    ALLOCATE(gf(m%np_part, 3), m%np_part*3)
    ALLOCATE(ggf(m%np_part, 3, 3), m%np_part*3*3)

    call X(f_angular_momentum)(sb, f_der, f, gf, ghost_update)

    do j = 1, 3
      call X(f_angular_momentum)(sb, f_der, gf(:,j), ggf(:,:,j), ghost_update)
    end do

    do j = 1, sb%dim
      l2f(1:m%np) = l2f(1:m%np) + ggf(1:m%np, j, j)
    end do

  case(2)
    ALLOCATE(gf(m%np_part, 1), m%np_part)
    ALLOCATE(ggf(m%np_part, 1, 1), m%np_part)

    call X(f_angular_momentum)(sb, f_der, f, gf, ghost_update)
    call X(f_angular_momentum)(sb, f_der, gf(:, 1), ggf(:, :, 1), ghost_update)

    l2f(1:m%np) = ggf(1:m%np, 1, 1)
  end select


  ! In case of real functions, since the angular momentum calculations
  ! lack a (-i) prefactor, we must add a (-1) factor
#if defined(R_TREAL)
  l2f = - l2f
#endif

  deallocate(gf, ggf)
  call pop_sub()
end subroutine X(f_l2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
