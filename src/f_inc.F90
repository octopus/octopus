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

!!! The next two subroutines convert a function between the normal 
!!! mesh and the cube
subroutine X(mf2cf) (m, mf, cf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in) :: mf(m%np)
  type(X(cf)), intent(inout) :: cf
  
  integer :: i, ix, iy, iz, c(3)

  ASSERT(associated(cf%RS))

  c(:) = cf%n(:)/2 + 1
  cf%RS = M_ZERO
 
  do i = 1, m%np
    ix = m%Lxyz(1, i) + c(1)
    iy = m%Lxyz(2, i) + c(2)
    iz = m%Lxyz(3, i) + c(3)
    cf%RS(ix, iy, iz) = mf(i)
  end do
end subroutine X(mf2cf)

subroutine X(cf2mf) (m, cf, mf)
  type(mesh_type), intent(in) :: m
  type(X(cf)), intent(in) :: cf
  R_TYPE, intent(out) :: mf(m%np)

  integer :: i, ix, iy, iz, c(3)

  ASSERT(associated(cf%RS))

  c(:) = cf%n(:)/2 + 1

  do i = 1, m%np
    ix = m%Lxyz(1, i) + c(1)
    iy = m%Lxyz(2, i) + c(2)
    iz = m%Lxyz(3, i) + c(3)
    mf(i) = cf%RS(ix, iy, iz)
  end do

end subroutine X(cf2mf)

!!! The next two subroutines convert a function in Fourier space
!!! between the normal mesh and the cube
subroutine X(mf2cf_FS) (m, mf, cf)
  type(mesh_type), intent(IN)    :: m
  CMPLX,           intent(IN)    :: mf(m%np)
  type(X(cf)),     intent(inout) :: cf
  
  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  cf%FS = M_z0
 
  do i = 1, m%np
    ix = pad_feq(m%Lxyz(1, i), cf%n(1), .false.)
    if(ix > cf%nx) cycle ! negative frequencies are redundant
    iy = pad_feq(m%Lxyz(2, i), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(3, i), cf%n(3), .false.)

    cf%FS(ix, iy, iz) = mf(i)
  end do
end subroutine X(mf2cf_FS)

subroutine X(cf_FS2mf) (m, cf, mf)
  type(mesh_type), intent(in)  :: m
  type(X(cf)),     intent(in)  :: cf
  CMPLX,           intent(out) :: mf(m%np)

  integer :: i, ix, iy, iz

  ASSERT(associated(cf%FS))

  do i = 1, m%np
    ix = pad_feq(m%Lxyz(1, i), cf%n(1), .false.)
    iy = pad_feq(m%Lxyz(2, i), cf%n(2), .false.)
    iz = pad_feq(m%Lxyz(3, i), cf%n(3), .false.)

#   ifdef R_TREAL
      if(ix > cf%nx) then
        ix = pad_feq(-m%Lxyz(1, i), cf%n(1), .false.)
        mf(i) = conjg(cf%FS(ix, iy, iz))
      else
        mf(i) = cf%FS(ix, iy, iz)
      end if
#   else
      mf(i) = cf%FS(ix, iy, iz)
#   endif
  end do

end subroutine X(cf_FS2mf)

! Calculation of derivatives
subroutine X(f_laplacian) (m, f, lapl, cutoff_)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: lapl(m%np)
  FLOAT, intent(in), optional :: cutoff_

  FLOAT :: cutoff

  call push_sub("f_laplacian")

  ASSERT(derivatives_space==REAL_SPACE.or.derivatives_space==RECIPROCAL_SPACE)

  select case(derivatives_space)
  case(REAL_SPACE)
    call X(mf_laplacian) (m, f, lapl)

#if defined(HAVE_FFT)
  case(RECIPROCAL_SPACE)

    ! Fixes the cutoff (negative value if optional argument cutoff was not passed)
    cutoff = -M_ONE
    if(present(cutoff_)) cutoff = cutoff_
    
    call X(cf_alloc_RS)(X(cf_der))             ! allocate cube in real space
    call X(cf_alloc_FS)(X(cf_der))             ! allocate cube in Fourier space
    
    call X(mf2cf)(m, f, X(cf_der))             ! convert to cube
    call X(cf_RS2FS)(X(cf_der))                ! Fourier transform
    call X(cf_FS_lapl)(m, X(cf_der), cutoff)   ! calculate Laplacian
    call X(cf_FS2RS)(X(cf_der))                ! Fourier transform back
    call X(cf2mf)(m, X(cf_der), lapl)          ! convert back to mesh

    call X(cf_free_RS)(X(cf_der))              ! clean memory
    call X(cf_free_FS)(X(cf_der))

#endif
  end select

  call pop_sub()
end subroutine X(f_laplacian)

subroutine X(f_gradient) (m, f, grad)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: grad(conf%dim, m%np)

  integer :: i, n

  call push_sub("f_gradient")

  ASSERT(derivatives_space==REAL_SPACE.or.derivatives_space==RECIPROCAL_SPACE)

  select case(derivatives_space)
  case(REAL_SPACE)
    call X(mf_gradient) (m, f, grad)

#if defined(HAVE_FFT)
  case(RECIPROCAL_SPACE)

    call X(cf_alloc_RS)(X(cf_aux))     ! allocate cube in real space
    call X(mf2cf)(m, f, X(cf_aux))     ! convert to cube

    call X(cf_alloc_FS)(X(cf_aux))     ! allocate cube in Fourier space
    call X(cf_RS2FS)(X(cf_aux))        ! Fourier transform

    call X(cf_free_RS)(X(cf_aux))      ! cube in real space is no longer needed

    call X(cf_alloc_RS)(X(cf_der))     ! allocate cube in real space
    call X(cf_alloc_FS)(X(cf_der))     ! allocate cube in Fourier space

    n = X(cf_aux)%nx*X(cf_aux)%n(2)*X(cf_aux)%n(3)
    do i = 1, conf%dim
      call zlalg_copy(n, X(cf_aux)%FS, X(cf_der)%FS)
      call X(cf_FS_grad)(m, X(cf_der), i)        ! gradient in reciprocal space
      call X(cf_FS2RS)(X(cf_der))                ! Fourier transform
      call X(cf2mf)(m, X(cf_der), grad(i, :))    ! convert to mesh
    end do

    call X(cf_free_FS)(X(cf_aux))      ! clean up
    call X(cf_free_RS)(X(cf_der))
    call X(cf_free_FS)(X(cf_der))

#endif
  end select

  call pop_sub()
end subroutine X(f_gradient)

subroutine X(f_divergence) (m, f, divf)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(conf%dim, m%np)
  R_TYPE, intent(out) :: divf(m%np)

  integer :: i
  R_TYPE, allocatable :: aux(:)

  call push_sub("f_divergence")

  ASSERT(derivatives_space==REAL_SPACE.or.derivatives_space==RECIPROCAL_SPACE)

  select case(derivatives_space)
  case(REAL_SPACE)
    call X(mf_divergence) (m, f, divf)

#if defined(HAVE_FFT)
  case(RECIPROCAL_SPACE)
    allocate(aux(m%np))
    call X(cf_alloc_RS)(X(cf_der))     ! allocate cube in real space
    call X(cf_alloc_FS)(X(cf_der))     ! allocate cube in real space

    do i = 1, conf%dim
      call X(mf2cf)(m, f(i, :), X(cf_der))     ! convert to cube
      call X(cf_RS2FS)(X(cf_der))        ! Fourier transform
      call X(cf_FS_grad)(m, X(cf_der), i)        ! gradient in reciprocal space
      call X(cf_FS2RS)(X(cf_der))                ! Fourier transform
      call X(cf2mf)(m, X(cf_der), aux)    ! convert to mesh
      
      if(i == 1) then
        divf = aux
      else
        divf = divf + aux
      end if
    end do
    
    call X(cf_free_RS)(X(cf_der))
    call X(cf_free_FS)(X(cf_der))
    deallocate(aux)
#endif
  end select

  call pop_sub()
end subroutine X(f_divergence)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The action of the angular momentum operator (three spatial components).
! In case of real functions, it does not include the -i prefactor
! (L = -i r ^ nabla).
subroutine X(f_angular_momentum)(m, f, lf)
  type(mesh_type), intent(IN)   :: m
  R_TYPE,          intent(IN)   :: f(m%np)
  R_TYPE,           intent(out) :: lf(3, m%np)

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x(3)
  integer :: i

  allocate(gf(3, m%np))
  call X(f_gradient)(m, f, grad = gf)

  do i = 1, m%np
     call mesh_xyz(m, i, x)
     lf(1, i) = (x(2)*gf(3, i)-x(3)*gf(2, i))
     lf(2, i) = (x(3)*gf(1, i)-x(1)*gf(3, i))
     lf(3, i) = (x(1)*gf(2, i)-x(2)*gf(1 ,i))
  enddo
#if defined(R_TCOMPLEX)
  call X(lalg_scal)(3*m%np, -M_zI, lf(1, 1))
#endif
  deallocate(gf)
end subroutine X(f_angular_momentum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Square of the angular momentum L. This has to be very much improved if
! accuracy is needed.
subroutine X(f_l2)(m, f, l2f)
  type(mesh_type), intent(IN)  :: m
  R_TYPE,          intent(IN)  :: f(:)
  R_TYPE,          intent(out) :: l2f(:)

  R_TYPE, allocatable :: gf(:, :), ggf(:, :, :)
  integer :: i, j

  allocate(gf(conf%dim, m%np), ggf(conf%dim, conf%dim, m%np))
  call X(f_angular_momentum)(m, f, gf)
  do j = 1, conf%dim
     call X(f_angular_momentum)(m, gf(j, 1:m%np), ggf(j, 1:conf%dim, 1:m%np))
  enddo

  do i = 1, m%np
     l2f(i) = M_ZERO
     do j = 1, conf%dim
        l2f(i) = l2f(i) + ggf(j, j, i)
     enddo
  enddo

! In case of real functions, since the angular momentum calculations
! lack a (-i) prefactor, we must add a (-1) factor
#if defined(R_TREAL)
  l2f = - l2f
#endif
  deallocate(gf, ggf)
end subroutine X(f_l2)

