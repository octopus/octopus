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

!!!  integrates a function
!!! (could not find BLAS routine to do it ;))
function R_FUNC(mf_integrate) (m, f) result(d)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)
  R_TYPE :: d

  d = sum(f(1:m%np))*m%vol_pp
end function R_FUNC(mf_integrate)

!!! this function returns the dot product between two vectors
!!! if HAVE_BLAS is defined it uses the BLAS library
R_TYPE function R_FUNC(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(m%np), f2(m%np)

#if defined(HAVE_BLAS)
  R_TYPE, external :: R_DOT

  dotp = R_DOT(m%np, f1(1), 1,  f2(1), 1)*m%vol_pp
#else
  dotp = sum(R_CONJ(f1)*f2)*m%vol_pp
#endif

end function R_FUNC(mf_dotp)

!!! this function returns the norm of a vector
!!! if HAVE_BLAS is defined it uses the BLAS library
real(r8) function R_FUNC(mf_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)

#if defined(HAVE_BLAS)
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np, f, 1)*sqrt(m%vol_pp)
#else
  nrm2 = sqrt(sum(R_CONJ(f)*f)*m%vol_pp)
#endif

end function R_FUNC(mf_nrm2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following set of subroutines calculates derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! calculates the laplacian and the gradient of a function on the mesh.
!!! The optional argument "alpha" is a multiplicative factor (i.e. if the
!!! kinetic operator wants to be calculated, pass "alpha = -1/2" to the
!!! calculation of the laplacian.
!!!
!!! If the laplacian is calculated on Fourier space, an optional argument cutoff
!!! may be passed, such that the calculation that is performed is 
!!! -min(G^2,cutoff)*phi(G), instead of just -G^2*phi(G). 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mf_laplacian) (m, f, lapl, cutoff_)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: lapl(m%np)
  real(r8), intent(in), optional :: cutoff_

  real(r8) :: cutoff

  call push_sub("mf_laplacian")

  select case(m%d%space)
  case(REAL_SPACE)
    call rs_derivative()

#if defined(HAVE_FFT)
  case(RECIPROCAL_SPACE)
    ! Fixes the cutoff (negative value if optional argument cutoff was not passed)
    cutoff = -M_ONE
    if(present(cutoff_)) cutoff = cutoff_

    call fft_derivative()
#endif

  case default
    message(1) = "Internal error in mf_laplacian"
    call write_fatal(1)
  end select
  
  call pop_sub()

contains
  subroutine rs_derivative()
    integer :: k, nl
    type(der_lookup_type), pointer :: p
    
    do k = 1, m%np
      p => m%der_lookup(k)
      nl = p%lapl_n
      lapl(k) = sum(p%lapl_w(1:nl)*f(p%lapl_i(1:nl)))
    end do
  end subroutine rs_derivative
  
  subroutine fft_derivative()
    complex(r8), allocatable :: fw(:,:,:), fw2(:,:,:)
    integer :: i, n(3)
    
    ! names are misleading, I know: they will be changed
#ifdef T_REAL
    call fft_getdim_complex(m%dfft, n)
#else
    call fft_getdim_real   (m%dfft, n)
#endif

    allocate(fw(n(1), n(2), n(3)))
    call R_FUNC(mesh_rs2fs)(m, f, fw)
    call R_FUNC(mesh_laplq)(m, fw, n, cutoff = cutoff)
    call R_FUNC(mesh_fs2rs)(m, fw, lapl)
    deallocate(fw)

  end subroutine fft_derivative
end subroutine R_FUNC(mf_laplacian)

subroutine R_FUNC(mf_gradient) (m, f, grad)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: grad(conf%dim, m%np)

  call push_sub("mf_gradient")

  select case(m%d%space)
  case(REAL_SPACE)
    call rs_derivative()
#if defined(HAVE_FFT)
  case(RECIPROCAL_SPACE)
    call fft_derivative()
#endif
  case default
    message(1) = "Internal error in mf_gradient"
    call write_fatal(1)
  end select

  call pop_sub()

contains
  subroutine rs_derivative()
    integer :: j, k, ng(3)
    type(der_lookup_type), pointer :: p
    
    do k = 1, m%np
      p => m%der_lookup(k)
      ng = p%grad_n
      do j = 1, conf%dim
        grad(j, k) = sum(p%grad_w(1:ng(j),j)*f(p%grad_i(1:ng(j),j)))
      end do
    end do

  end subroutine rs_derivative

#if defined(HAVE_FFT)
  subroutine fft_derivative()
    complex(r8), allocatable :: fw(:,:,:), fw2(:,:,:)
    integer :: i, n(3)
    
    ! names are misleading, I know: they will be changed
#ifdef T_REAL
    call fft_getdim_complex(m%dfft, n)
#else
    call fft_getdim_real   (m%dfft, n)
#endif
   
    allocate(fw(n(1), n(2), n(3)), fw2(n(1), n(2), n(3)))
    call R_FUNC(mesh_rs2fs)(m, f, fw)
    do i = 1, conf%dim
      call R_FUNC(copy)(n(1)*n(2)*n(3), fw, 1, fw2, 1)
      call R_FUNC(mesh_gradq)(m, fw2, i, n)
      call R_FUNC(mesh_fs2rs)(m, fw2, grad(i, :))
    end do
    deallocate(fw, fw2)
   
  end subroutine fft_derivative
#endif
end subroutine R_FUNC(mf_gradient)

!!! Calculates the divergence of a vectorial function f.
!!! Currently it only does so in real space.
subroutine R_FUNC(mf_divergence)(m, f, divf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(conf%dim, m%np)
  R_TYPE, intent(out) :: divf(m%np)

  integer :: j, k, ng(3)
  type(der_lookup_type), pointer :: p

  call push_sub('mf_divergence')

  do k = 1, m%np
    p => m%der_lookup(k)
    ng = p%grad_n
    divf(k) = M_ZERO
    do j = 1, conf%dim
      divf(k) = divf(k) + sum(p%grad_w(1:ng(j),j)*f(j, p%grad_i(1:ng(j), j)))
    end do
  end do

  call pop_sub()
end subroutine R_FUNC(mf_divergence)
