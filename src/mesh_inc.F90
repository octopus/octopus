!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! Conversion subroutines (they actually add, do not forget it)
! They also work, in principle, for 1 and 2D

! WARNING: I think these two routines are *not* working
subroutine R_FUNC(mesh_to_cube) (m, f_mesh, f_cube, l)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_mesh(m%np)
  R_TYPE, intent(inout) :: f_cube(:, :, :)
  integer, intent(in) :: l(3)
  
  integer :: i, ix, iy, iz, n(3)

  n(:) = l(:)/2 + 1
 
  do i = 1, m%np
    ix = m%Lxyz(1, i) + n(1)
    iy = m%Lxyz(2, i) + n(2)
    iz = m%Lxyz(3, i) + n(3)
    f_cube(ix, iy, iz) = f_cube(ix, iy, iz) + f_mesh(i)
  end do
end subroutine R_FUNC(mesh_to_cube)

subroutine R_FUNC(cube_to_mesh) (m, f_cube, f_mesh, l)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_cube(:, :, :)
  R_TYPE, intent(inout) :: f_mesh(m%np)
  integer, intent(in)   :: l(3)

  integer :: i, ix, iy, iz, n(3)

  n(:) = l(:)/2 + 1

  do i = 1, m%np
    ix = m%Lxyz(1, i) + n(1)
    iy = m%Lxyz(2, i) + n(2)
    iz = m%Lxyz(3, i) + n(3)
    f_mesh(i) = f_mesh(i) + f_cube(ix, iy, iz) 
  end do
end subroutine R_FUNC(cube_to_mesh)

! this functions returns the dot product between two vectors
! it uses BLAS
R_TYPE function R_FUNC(mesh_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(1:m%np), f2(1:m%np)
  R_TYPE, external :: R_DOT
  
  dotp = R_DOT(m%np, f1(1), 1,  f2(1), 1)*m%vol_pp
end function R_FUNC(mesh_dotp)

!!$complex(r8) function R_FUNC(mesh_dotpq)(m, f1, f2) result(dotpq)
!!$  type(mesh_type), intent(IN) :: m
!!$  complex(r8), intent(IN), dimension(*) :: f1, f2
!!$  complex(r8), external :: zdotc
!!$
!!$#ifdef R_TREAL
!!$  complex(r8) :: z
!!$  ! Wish I had commented this?
!!$  z = zdotc(m%R_FUNC(npw), f1(1), 1, f2(1), 1)
!!$  dotpq = z + conjg(z) - zdotc(m%fft_n(2)*m%fft_n(3), f1(1), m%hfft_n, f2(1), m%hfft_n)
!!$  dotpq = dotpq*m%vol_ppw
!!$#else
!!$  dotpq = zdotc(m%R_FUNC(npw), f1(1), 1,  f2(1), 1)*m%vol_ppw
!!$#endif
!!$
!!$end function R_FUNC(mesh_dotpq)

real(r8) function R_FUNC(mesh_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np, f, 1)*sqrt(m%vol_pp)
end function R_FUNC(mesh_nrm2)

!!$real(r8) function R_FUNC(mesh_nrm2q)(m, f) result(nrm2)
!!$  type(mesh_type), intent(IN) :: m
!!$  complex(r8), intent(IN), dimension(*) :: f
!!$  real(r8), external :: dznrm2
!!$
!!$#ifdef R_TREAL
!!$  nrm2 = sqrt(R_FUNC(mesh_dotpq)(m, f, f))
!!$#else
!!$  nrm2 = dznrm2(m%R_FUNC(npw), f, 1)*sqrt(m%vol_ppw)
!!$#endif
!!$end function R_FUNC(mesh_nrm2q)

! integrates a function on the mesh (could not find BLAS routine to do it ;))
function R_FUNC(mesh_integrate) (m, f)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)
  R_TYPE :: R_FUNC(mesh_integrate)

  R_FUNC(mesh_integrate) = sum(f(1:m%np))*m%vol_pp
end function R_FUNC(mesh_integrate)

function R_FUNC(mesh_moment) (m, f, i, n)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(IN)          :: f(:)
  integer, intent(in)         :: i, n
  R_TYPE ::R_FUNC(mesh_moment)

  R_FUNC(mesh_moment) = sum(f(1:m%np)*(m%lxyz(i, 1:m%np)*m%h(i))**n)*m%vol_pp
end function R_FUNC(mesh_moment)

!!$function R_FUNC(mesh_integrateq) (m, f)
!!$  type(mesh_type), intent(IN) :: m
!!$  complex(r8), intent(IN), dimension(*) :: f
!!$  R_TYPE :: R_FUNC(mesh_integrateq)
!!$  R_FUNC(mesh_integrateq) = f(1)*m%vol_pp
!!$end function R_FUNC(mesh_integrateq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Outputs in g the FS representation of a function represented in the real-
! space mesh, f. So this function has to be taken first from the "mesh" to the
! "cube" (parallelpiped-type mesh), and then FFTed to Fourier space.
!
! The dimensions of g are different wether f is real or complex, because the
! FFT representation is different (FFTW scheme).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef HAVE_FFT
subroutine R_FUNC(mesh_rs2fs)(m, f, g)
  type(mesh_type), intent(in)  :: m
  R_TYPE, intent(in)       :: f(:)
  complex(r8), intent(out) :: g(:,:,:)

  R_TYPE, allocatable :: fr(:, :, :)
  integer :: n(3)

  call push_sub('mesh_rs2fs')

  call fft_getdim_real   (m%dfft, n)
  allocate(fr(n(1), n(2), n(3)))
  fr = R_TOTYPE(M_ZERO)
  call R_FUNC(mesh_to_cube) (m, f, fr, n)
#ifdef R_TREAL
  call rfft_forward(m%dfft, fr, g)
#else
  call zfft_forward(m%zfft, fr, g)
#endif

  deallocate(fr)
  call pop_sub()
end subroutine R_FUNC(mesh_rs2fs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The opposite of mesh_rs2fs.
! NB: If FFTOptimize = .true., these subroutines are not the exact inverse of
!     each other!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mesh_fs2rs)(m, g, f)
  type(mesh_type), intent(in) :: m
  complex(r8), intent(inout)  :: g(:,:,:)
  R_TYPE, intent(out)         :: f(:)

  R_TYPE, allocatable :: fr(:, :, :)
  integer :: n(3)

  call push_sub('mesh_fs2rs')

  call fft_getdim_real(m%dfft, n)

  allocate(fr(n(1), n(2), n(3)))   
#ifdef R_TREAL
  call rfft_backward(m%dfft, g, fr)
#else
  call zfft_backward(m%zfft, g, fr)
#endif

  f = R_TOTYPE(M_ZERO)
  call R_FUNC(cube_to_mesh) (m, fr, f, n)
  deallocate(fr)

  call pop_sub()
end subroutine R_FUNC(mesh_fs2rs)

#endif

subroutine R_FUNC(mesh_random)(m, f)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(out), dimension(*) :: f

  integer, save :: iseed = 123
  integer :: i
  real(r8) :: a(3), rnd, r

  call push_sub('mesh_random')

  call quickrnd(iseed, rnd)
  a(1) = 2.0_r8*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(2) = 2.0_r8*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(3) = 2.0_r8*(2*rnd - 1)

  do i = 1, m%np
     call mesh_r(m, i, r, a=a)
     f(i) = exp(-M_HALF*r*r)
  end do

  r = R_FUNC(mesh_nrm2)(m, f)
  call R_FUNC(scal)(m%np, R_TOTYPE(M_ONE/r), f, 1) 

  call pop_sub()
end subroutine R_FUNC(mesh_random)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the divergence of a vectorial function f.
! Currently it only does so in real space.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mesh_divergence)(m, f, divf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(3, m%np)
  R_TYPE, intent(out) :: divf(m%np)

  integer :: j, k, ng(3)
  type(der_lookup_type), pointer :: p

  call push_sub('mesh_divergence')

  do k = 1, m%np
     p => m%der_lookup(k)
     ng = p%grad_n
     divf(k) = M_ZERO
     do j = 1, conf%dim
        divf(k) = divf(k) + sum(p%grad_w(1:ng(j),j)*f(j, p%grad_i(1:ng(j), j)))
     enddo
  enddo

  call pop_sub(); return
end subroutine R_FUNC(mesh_divergence)

subroutine R_FUNC(mesh_angular_momentum)(m, f, lf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(m%np)
  R_TYPE, intent(out) :: lf(3, m%np)
  R_TYPE, allocatable :: gf(:, :)
  real(r8) :: x(3)
  integer :: i
  allocate(gf(3, m%np))
  call R_FUNC(mesh_derivatives)(m, f, grad = gf)
  do i = 1, m%np
     call mesh_xyz(m, i, x)
     lf(1, i) = (x(2)*gf(3, i)-x(3)*gf(2, i))
     lf(2, i) = (x(3)*gf(1, i)-x(1)*gf(3, i))
     lf(3, i) = (x(1)*gf(2, i)-x(2)*gf(1 ,i))
  enddo
#if defined(R_TCOMPLEX)
  call zscal(3*m%np, -M_zI, lf, 1)
#endif
  deallocate(gf)
end subroutine R_FUNC(mesh_angular_momentum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculates the laplacian and the gradient of a function on the mesh.
! The optional argument "alpha" is a multiplicative factor (i.e. if the
! kinetic operator wants to be calculated, pass "alpha = -1/2" to the
! calculation of the laplacian.
! If laplacian is calculated on Fourier space, an optional argument cutoff
! may be passed, such that the calculation that is performed is 
!-min(G^2,cutoff)*phi(G), instead of just -G^2*phi(G). 
! In the case of receiving both alpha and cutoff, the cutoff is internally
! multiplied by the absolute value of the inverse of alpha.
! (Rationale: if we want to calculate min(G^2/2, E_cutoff), this is equal
! to 1/2*min(G^2,2*E_cutoff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mesh_derivatives) (m, f, lapl, grad, alpha, cutoff)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)
  R_TYPE, intent(out), optional:: lapl(:), grad(:,:)
  R_TYPE, intent(in), optional :: alpha
  real(r8), intent(in), optional :: cutoff

  R_TYPE :: alp
  real(r8) :: lcutoff

  call push_sub('mesh_derivatives')
  
  ! Fixes the multiplicative factor.
  alp = R_TOTYPE(M_ONE)
  if(present(alpha)) alp = alpha

  select case(m%d%space)
    case(REAL_SPACE)
      call rs_derivative()

#ifdef HAVE_FFT
    case(RECIPROCAL_SPACE)
      ! Fixes the cutoff (negative value if optional argument cutoff was not passed)
      lcutoff = -M_ONE
      if(present(cutoff)) lcutoff = cutoff
      if(alp.ne.M_z0) lcutoff = lcutoff*M_ONE/abs(alp)

      call fft_derivative()
#endif
    case default
      message(1) = "Internal error in mesh_derivatives"
      call write_fatal(1)
  end select

  call pop_sub()
contains

#ifdef HAVE_FFT
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

    if(present(grad)) then
      allocate(fw2(n(1), n(2), n(3)))
      do i = 1, 3
        call R_FUNC(copy)(n(1)*n(2)*n(3), fw, 1, fw2, 1)
        call R_FUNC(mesh_gradq)(m, fw2, i, n)
        call R_FUNC(mesh_fs2rs)(m, fw2, grad(i, :))
      end do
      deallocate(fw2)
      if(present(alpha)) call R_FUNC(scal)(m%np*3, alpha, grad, 1)
    end if

    if(present(lapl)) then
      call R_FUNC(mesh_laplq)(m, fw, n, cutoff = lcutoff)
      call R_FUNC(mesh_fs2rs)(m, fw, lapl)
      if(present(alpha)) call R_FUNC(scal)(m%np, alp, lapl, 1)
    end if

    deallocate(fw)

    return
  end subroutine fft_derivative
#endif

  subroutine rs_derivative()
    integer :: j, k, nl, ng(3)
    type(der_lookup_type), pointer :: p

    do k = 1, m%np
      p => m%der_lookup(k)
      if(present(lapl)) then
        nl = p%lapl_n
        lapl(k) = alp*sum(p%lapl_w(1:nl)*f(p%lapl_i(1:nl)))
      end if
      if(present(grad)) then
        ng = p%grad_n
        do j = 1, conf%dim
          grad(j, k) = alp*sum(p%grad_w(1:ng(j),j)*f(p%grad_i(1:ng(j),j)))
        end do
      end if
    end do

  end subroutine rs_derivative

end subroutine R_FUNC(mesh_derivatives)

! WARNING: This is not currently working
subroutine R_FUNC(low_frequency) (m, f, lapl)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)  :: f(m%np)
  R_TYPE, intent(out) :: lapl(1:m%np)

  integer :: k, nl
  type(der_lookup_type), pointer :: p

  call push_sub('low_frequency')

!!$  do k = 1, m%np
!!$    p => m%der_lookup(k)
!!$    nl = p%lapl_n
!!$    lapl(k) = sum(p%lapl_w(1:nl)*f(p%lapl_i(1:nl)))
!!$  end do

  call pop_sub()
end subroutine R_FUNC(low_frequency)

! WARNING the next two functions need to be cleaned
#ifdef HAVE_FFT
subroutine R_FUNC(mesh_laplq)(m, f, n, cutoff, exponential)
  type(mesh_type), intent(in)    :: m
  complex(r8), dimension(*) :: f
  integer, intent(in) :: n(3)
  real(r8), intent(in), optional :: cutoff
  R_TYPE, intent(in), optional :: exponential

  real(r8) :: lcutoff, temp(3), g2
  integer :: i, k(3), ix, iy, iz, nx
  R_TYPE :: lexp

  call push_sub('mesh_laplq')

  lcutoff = 1e10_r8
  if(present(cutoff)) then
     if(cutoff>M_ZERO) then
       lcutoff = cutoff
     endif
  endif

  lexp = R_TOTYPE(M_ZERO)
  if(present(exponential)) lexp = exponential

  nx = n(1)
#if defined(R_TREAL)
  nx = n(1)/2+1
#endif

  temp = M_ZERO
  temp(1:conf%dim) = (2.0_r8*M_Pi)/(m%l(1:conf%dim)*m%h(1:conf%dim))

  if(present(exponential)) then
    i = 0
    do iz = 1, n(3)
       k(3) = pad_feq(iz, m%l(3), .true.)
       do iy = 1, n(2)
          k(2) = pad_feq(iy, m%l(2), .true.)
          do ix = 1, nx
             k(1) = pad_feq(ix, m%l(1), .true.)
             g2 = min(lcutoff, sum((temp(1:conf%dim)*k(1:conf%dim))**2))
             i = i + 1
             f(i) = exp(-lexp*g2)*f(i)
          enddo
       enddo
    enddo
  else
    i = 0
    do iz = 1, n(3)
       k(3) = pad_feq(iz, m%l(3), .true.)
       do iy = 1, n(2)
          k(2) = pad_feq(iy, m%l(2), .true.)
          do ix = 1, nx
             k(1) = pad_feq(ix, m%l(1), .true.)
             g2 = min(lcutoff, sum((temp(1:conf%dim)*k(1:conf%dim))**2))
             i = i + 1
             f(i) = - g2*f(i)
          enddo
       enddo
    enddo
  endif

  call pop_sub()
end subroutine R_FUNC(mesh_laplq)

subroutine R_FUNC(mesh_gradq)(m, f, j, n)
  type(mesh_type), intent(in) :: m
  complex(r8), dimension(*)   :: f
  integer, intent(in) :: j, n(3)

  real(r8) :: lcutoff, temp(3), g2
  integer :: i, k(3), ix, iy, iz, nx

  call push_sub('mesh_gradq')

  nx = n(1)
#if defined(R_TREAL)
  nx = n(1)/2+1
#endif

  temp = M_ZERO
  temp(1:conf%dim) = (2.0_r8*M_Pi)/(n(1:conf%dim)*m%h(1:conf%dim))
  i = 0
  k = 0
  do iz = 1, n(3)
    if(j == 3) k(3) = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      if(j == 2) k(2) = pad_feq(iy, n(2), .true.)
      do ix = 1, nx
        if(j == 1) k(1) = pad_feq(ix, n(1), .true.)
        g2 = sum(temp(1:conf%dim)*k(1:conf%dim))
        i = i + 1
        f(i) = g2 * M_zI * f(i)
      end do
    end do
  end do

  call pop_sub()
end subroutine R_FUNC(mesh_gradq)
#endif
