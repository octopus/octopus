!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
! mesh and the cube
! ---------------------------------------------------------
subroutine X(mf2cf) (m, mf, cf)
  type(mesh_type), intent(in)    :: m
  R_TYPE,          intent(in)    :: mf(:)  ! mf(f_der%m%np)
  type(X(cf)),     intent(inout) :: cf

  integer :: i, ix, iy, iz, c(3)

  ASSERT(associated(cf%RS))

  c(:)  =  cf%n(:)/2 + 1
  cf%RS =  M_ZERO

  do i = 1, m%np
     ix = m%Lxyz(i, 1) + c(1)
     iy = m%Lxyz(i, 2) + c(2)
     iz = m%Lxyz(i, 3) + c(3)

     cf%RS(ix, iy, iz) = mf(i)
  end do

end subroutine X(mf2cf)


! ---------------------------------------------------------
subroutine X(cf2mf) (m, cf, mf)
  type(mesh_type), intent(in)  :: m
  type(X(cf)),     intent(in)  :: cf
  R_TYPE,          intent(out) :: mf(:)  ! mf(f_der%m%np)

  integer :: i, ix, iy, iz, c(3)

  ASSERT(associated(cf%RS))

  c(:) =  cf%n(:)/2 + 1

  do i = 1, m%np
    ix = m%Lxyz(i, 1) + c(1)
    iy = m%Lxyz(i, 2) + c(2)
    iz = m%Lxyz(i, 3) + c(3)
    mf(i) = cf%RS(ix, iy, iz)
  end do

end subroutine X(cf2mf)

! ---------------------------------------------------------
! The next two subroutines convert a function in Fourier space
! between the normal mesh and the cube
! ---------------------------------------------------------
subroutine X(mf2cf_FS) (m, mf, cf)
  type(mesh_type), intent(in)    :: m
  CMPLX,               intent(in)    :: mf(:)   ! mf(f_der%m%np)
  type(X(cf)),         intent(inout) :: cf

  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  cf%FS = M_z0

  do i = 1, m%np
    ix = pad_feq(m%Lxyz(i, 1), cf%n(1), .false.)
    if(ix > cf%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(m%Lxyz(i, 2), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(i, 3), cf%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(i)
  end do
end subroutine X(mf2cf_FS)


! ---------------------------------------------------------
subroutine X(cf_FS2mf) (m, cf, mf)
  type(mesh_type), intent(in)  :: m
  type(X(cf)),     intent(in)  :: cf
  CMPLX,           intent(out) :: mf(:) ! mf(f_der%m%np)

  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  do i = 1, m%np
    ix = pad_feq(m%Lxyz(i, 1), cf%n(1), .false.)
    iy = pad_feq(m%Lxyz(i, 2), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(i, 3), cf%n(3), .false.)

#   ifdef R_TREAL
      if(ix > cf%nx) then
        ix = pad_feq(-m%Lxyz(i, 1), cf%n(1), .false.)
        mf(i) = conjg(cf%FS(ix, iy, iz))
      else
        mf(i) = cf%FS(ix, iy, iz)
      end if
#   else
      mf(i) = cf%FS(ix, iy, iz)
#   endif
  end do

end subroutine X(cf_FS2mf)


! ---------------------------------------------------------
! Calculation of derivatives
! ---------------------------------------------------------
subroutine X(f_laplacian) (sb, f_der, f, lapl, cutoff_)
  type(simul_box_type), intent(in) :: sb
  type(f_der_type), intent(inout)  :: f_der
  R_TYPE,           intent(inout)  :: f(:)     ! f(m%np_part)
  R_TYPE,           intent(out)    :: lapl(:)  ! lapl(m%np_part)
  FLOAT, optional,  intent(in)     :: cutoff_

  FLOAT :: cutoff

  call push_sub('f_inc.f_laplacian')

  ASSERT(f_der%space==REAL_SPACE.or.f_der%space==FOURIER_SPACE)

  select case(f_der%space)
  case(REAL_SPACE)
     call X(derivatives_lapl) (f_der%der_discr, f, lapl)

#if defined(HAVE_FFT)
  case(FOURIER_SPACE)

     ! Fixes the cutoff (negative value if optional argument cutoff was not passed)
     cutoff = -M_ONE
     if(present(cutoff_)) cutoff = cutoff_

     call X(cf_alloc_RS)(f_der%X(cf_der))             ! allocate cube in real space
     call X(cf_alloc_FS)(f_der%X(cf_der))             ! allocate cube in Fourier space

     call X(mf2cf)(f_der%m, f, f_der%X(cf_der))             ! convert to cube
     call X(cf_RS2FS)(f_der%X(cf_der))                ! Fourier transform
     call X(cf_FS_lapl)(sb, f_der%m, f_der%X(cf_der), cutoff)   ! calculate Laplacian
     call X(cf_FS2RS)(f_der%X(cf_der))                ! Fourier transform back
     call X(cf2mf)(f_der%m, f_der%X(cf_der), lapl)          ! convert back to mesh

     call X(cf_free_RS)(f_der%X(cf_der))              ! clean memory
     call X(cf_free_FS)(f_der%X(cf_der))

#endif
  end select

  call pop_sub()
end subroutine X(f_laplacian)


! ---------------------------------------------------------
subroutine X(f_gradient) (sb, f_der, f, grad)
  type(simul_box_type), intent(in)    :: sb
  type(f_der_type),     intent(inout) :: f_der
  R_TYPE, target,       intent(in)    :: f(:)       ! f(m%np_part)
  R_TYPE,               intent(out)   :: grad(:,:)  ! grad(m%np_part, m%sb%dim)

  integer :: i

  call push_sub('f_inc.f_gradient')

  ASSERT(f_der%space==REAL_SPACE.or.f_der%space==FOURIER_SPACE)

  grad = R_TOTYPE(M_ZERO)

  select case(f_der%space)
  case(REAL_SPACE)
    call X(derivatives_grad) (sb, f_der%der_discr, f, grad)

#if defined(HAVE_FFT)
  case(FOURIER_SPACE)

    call X(cf_alloc_RS)(f_der%X(cf_aux))           ! allocate cube in real space
    call X(mf2cf)      (f_der%m, f, f_der%X(cf_aux))  ! convert to cube

    call X(cf_alloc_FS)(f_der%X(cf_aux))           ! allocate cube in Fourier space
    call X(cf_RS2FS)   (f_der%X(cf_aux))           ! Fourier transform

    call X(cf_free_RS) (f_der%X(cf_aux))           ! cube in real space is no longer needed

    call X(cf_alloc_RS)(f_der%X(cf_der))           ! allocate cube in real space
    call X(cf_alloc_FS)(f_der%X(cf_der))           ! allocate cube in Fourier space

    do i = 1, sb%dim
      call lalg_copy(f_der%X(cf_aux)%nx, f_der%X(cf_aux)%n(2), f_der%X(cf_aux)%n(3), &
         f_der%X(cf_aux)%FS(:,:,:), f_der%X(cf_der)%FS(:,:,:))

      call X(cf_FS_grad)(sb, f_der%m, f_der%X(cf_der), i)      ! gradient in reciprocal space
      call X(cf_FS2RS)  (f_der%X(cf_der))                      ! Fourier transform
      call X(cf2mf)     (f_der%m, f_der%X(cf_der), grad(:,i))  ! convert to mesh
    end do

    call X(cf_free_FS)(f_der%X(cf_aux))      ! clean up
    call X(cf_free_RS)(f_der%X(cf_der))
    call X(cf_free_FS)(f_der%X(cf_der))

#endif
  end select

  call pop_sub()
end subroutine X(f_gradient)


! ---------------------------------------------------------
subroutine X(f_divergence) (sb, f_der, f, divf)
  type(simul_box_type), intent(in)    :: sb
  type(f_der_type),     intent(inout) :: f_der
  R_TYPE,               intent(in)    :: f(:,:)    ! f(sb%dim, m%np)
  R_TYPE,               intent(out)   :: divf(:)   ! divf(m%np)

  integer :: i
  R_TYPE, allocatable :: aux(:)

  call push_sub('f_inc.f_divergence')

  ASSERT(f_der%space==REAL_SPACE.or.f_der%space==FOURIER_SPACE)

  select case(f_der%space)
  case(REAL_SPACE)
    call X(derivatives_div) (sb, f_der%der_discr, f, divf)

#if defined(HAVE_FFT)
  case(FOURIER_SPACE)
    allocate(aux(f_der%m%np))
    call X(cf_alloc_RS)(f_der%X(cf_der))     ! allocate cube in real space
    call X(cf_alloc_FS)(f_der%X(cf_der))     ! allocate cube in real space

    do i = 1, sb%dim
      call X(mf2cf)     (f_der%m, f(i, :), f_der%X(cf_der)) ! convert to cube
      call X(cf_RS2FS)  (f_der%X(cf_der))                   ! Fourier transform
      call X(cf_FS_grad)(sb, f_der%m, f_der%X(cf_der), i)   ! gradient in reciprocal space
      call X(cf_FS2RS)  (f_der%X(cf_der))                   ! Fourier transform
      call X(cf2mf)     (f_der%m, f_der%X(cf_der), aux)     ! convert to mesh

      if(i == 1) then
        divf = aux
      else
        divf = divf + aux
      end if
    end do

    call X(cf_free_RS)(f_der%X(cf_der))
    call X(cf_free_FS)(f_der%X(cf_der))
    deallocate(aux)
#endif
  end select

  call pop_sub()
end subroutine X(f_divergence)

! ---------------------------------------------------------
subroutine X(f_curl) (sb, f_der, f, curlf)
  type(simul_box_type), intent(in) :: sb
  type(f_der_type), intent(inout) :: f_der
  R_TYPE,           intent(in)    :: f(:,:)     ! f(m%np, conf%dim)
  R_TYPE,           intent(out)   :: curlf(:,:) ! curlf(m%np, conf%dim))

  call push_sub('f_inc.f_curl')

  ASSERT(f_der%space==REAL_SPACE.or.f_der%space==FOURIER_SPACE)

  select case(f_der%space)
  case(REAL_SPACE)
    call X(derivatives_curl) (sb, f_der%der_discr, f, curlf)

#if defined(HAVE_FFT)
  case(FOURIER_SPACE)
    message(1) = "curl calculation in fourier space not yet implemented"
    call write_fatal(1)
#endif
  end select

  call pop_sub()
end subroutine X(f_curl)

! ---------------------------------------------------------
! The action of the angular momentum operator (three spatial components).
! In case of real functions, it does not include the -i prefactor
! (L = -i r ^ nabla).
! ---------------------------------------------------------
subroutine X(f_angular_momentum)(sb, f_der, f, lf)
  type(simul_box_type), intent(in) :: sb
  type(f_der_type), intent(inout) :: f_der
  R_TYPE,           intent(in)    :: f(1:f_der%m%np)     ! f(m%np)
  R_TYPE,           intent(out)   :: lf(:,:)  ! lf(m%np, 3)

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x(3)
  integer :: i

  call push_sub('f_inc.f_angular_momentum')

  allocate(gf(f_der%m%np, 3))
  call X(f_gradient)(sb, f_der, f, gf)

  do i = 1, f_der%m%np
    x = f_der%m%x(i,:)
    lf(i, 1) = (x(2)*gf(i, 3)-x(3)*gf(i, 2))
    lf(i, 2) = (x(3)*gf(i, 1)-x(1)*gf(i, 3))
    lf(i, 3) = (x(1)*gf(i, 2)-x(2)*gf(i, 1))
  end do
#if defined(R_TCOMPLEX)
  call lalg_scal(3, f_der%m%np, -M_zI, lf)
#endif

  deallocate(gf)
  call pop_sub()
end subroutine X(f_angular_momentum)


! ---------------------------------------------------------
! Square of the angular momentum L. This has to be very much improved if
! accuracy is needed.
! ---------------------------------------------------------
subroutine X(f_l2)(sb, f_der, f, l2f)
  type(simul_box_type), intent(in) :: sb
  type(f_der_type), intent(inout) :: f_der
  R_TYPE,           intent(in)    :: f(1:f_der%m%np)
  R_TYPE,           intent(out)   :: l2f(:)

  R_TYPE, allocatable :: gf(:, :), ggf(:, :, :)
  type(mesh_type), pointer :: m
  integer :: j

  m => f_der%m

  allocate(gf(m%np, sb%dim), ggf(m%np, sb%dim, sb%dim))

  call X(f_angular_momentum)(sb, f_der, f, gf)
  do j = 1, 3
    call X(f_angular_momentum)(sb, f_der, gf(:,j), ggf(:,:,j))
  end do

  l2f(:) = M_ZERO
  do j = 1, 3
    l2f(:) = l2f(:) + ggf(:, j, j)
  end do

  ! In case of real functions, since the angular momentum calculations
  ! lack a (-i) prefactor, we must add a (-1) factor
#if defined(R_TREAL)
  l2f = - l2f
#endif
  deallocate(gf, ggf)
end subroutine X(f_l2)

