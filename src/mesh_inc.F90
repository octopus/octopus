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
subroutine R_FUNC(mesh_to_cube) (m, f_mesh, f_cube, t)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_mesh(m%np)
  R_TYPE, intent(inout) :: f_cube(:, :, :)
  integer, intent(in), optional :: t
  
  integer :: i, ix, iy, iz, n(3)

  n(1:3) = m%fft_n(1:3)/2 + 1
  if(present(t)) then
    if(t == 2) n(1:3) = m%fft_n2(1:3)/2 + 1
  end if
 
  do i = 1, m%np
    ix = m%Lxyz(1, i) + n(1)
    iy = m%Lxyz(2, i) + n(2)
    iz = m%Lxyz(3, i) + n(3)
    f_cube(ix, iy, iz) = f_cube(ix, iy, iz) + f_mesh(i)
  end do
end subroutine R_FUNC(mesh_to_cube)

subroutine R_FUNC(cube_to_mesh) (m, f_cube, f_mesh, t)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_cube(:, :, :)
  R_TYPE, intent(inout) :: f_mesh(m%np)
  integer, intent(in), optional :: t

  integer :: i, ix, iy, iz, n(3)

  n(1:3) = m%fft_n(1:3)/2 + 1
  if(present(t)) then
    if(t == 2) n(1:3) = m%fft_n2(1:3)/2 + 1
  end if

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

complex(r8) function R_FUNC(mesh_dotpq)(m, f1, f2) result(dotpq)
  type(mesh_type), intent(IN) :: m
  complex(r8), intent(IN), dimension(*) :: f1, f2
  complex(r8), external :: zdotc

#ifdef R_TREAL
  complex(r8) :: z
  ! Wish I had commented this?
  z = zdotc(m%R_FUNC(npw), f1(1), 1, f2(1), 1)
  dotpq = z + conjg(z) - zdotc(m%fft_n(2)*m%fft_n(3), f1(1), m%hfft_n, f2(1), m%hfft_n)
  dotpq = dotpq*m%vol_ppw
#else
  dotpq = zdotc(m%R_FUNC(npw), f1(1), 1,  f2(1), 1)*m%vol_ppw
#endif

end function R_FUNC(mesh_dotpq)

real(r8) function R_FUNC(mesh_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np, f, 1)*sqrt(m%vol_pp)
end function R_FUNC(mesh_nrm2)

real(r8) function R_FUNC(mesh_nrm2q)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  complex(r8), intent(IN), dimension(*) :: f
  real(r8), external :: dznrm2

#ifdef R_TREAL
  nrm2 = sqrt(R_FUNC(mesh_dotpq)(m, f, f))
#else
  nrm2 = dznrm2(m%R_FUNC(npw), f, 1)*sqrt(m%vol_ppw)
#endif
end function R_FUNC(mesh_nrm2q)

! integrates a function on the mesh (could not find BLAS routine to do it ;))
function R_FUNC(mesh_integrate) (m, f)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)
  R_TYPE :: R_FUNC(mesh_integrate)
  R_FUNC(mesh_integrate) = sum(f(1:m%np))*m%vol_pp
end function R_FUNC(mesh_integrate)

function R_FUNC(mesh_integrateq) (m, f)
  type(mesh_type), intent(IN) :: m
  complex(r8), intent(IN), dimension(*) :: f
  R_TYPE :: R_FUNC(mesh_integrateq)
  R_FUNC(mesh_integrateq) = f(1)*m%vol_pp
end function R_FUNC(mesh_integrateq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Outputs in g the FS representation of a function represented in the real-
! space mesh, f. So this function has to be taken first from the "mesh" to the
! "cube" (parallelpiped-type mesh), and then FFTed to Fourier space.
!
! The dimensions of g are different wether f is real or complex, because the
! FFT representation is different (FFTW scheme).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mesh_rs2fs)(m, f, g)
  type(mesh_type), intent(in)  :: m
  R_TYPE, intent(in), dimension(*)       :: f
  complex(r8), intent(out), dimension(*) :: g

  R_TYPE, allocatable :: fr(:, :, :)

  sub_name = 'mesh_rs2fs'; call push_sub()

  allocate(fr(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
  fr = R_TOTYPE(M_ZERO)
  call R_FUNC(mesh_to_cube) (m, f, fr)
#ifdef R_TREAL
  call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, g)
#else
  call fftwnd_f77_one(m%zplanf, fr, g)
#endif

  deallocate(fr)
  call pop_sub(); return
end subroutine R_FUNC(mesh_rs2fs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The opposite of mesh_rs2fs.
! NB: If FFTOptimize = .true., these subroutines are not the exact inverse of
!     each other!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_FUNC(mesh_fs2rs)(m, g, f, factor)
  type(mesh_type), intent(in)  :: m
  complex(r8), intent(inout), dimension(*) :: g
  R_TYPE, intent(out), dimension(*)     :: f
  R_TYPE, intent(in), optional :: factor

  R_TYPE, allocatable :: fr(:, :, :)
  R_TYPE :: fac

  complex(r8), allocatable :: g_aux(:)
  allocate(g_aux(m%R_FUNC(npw)))
  g_aux(1:m%R_FUNC(npw)) = g(1:m%R_FUNC(npw))
  sub_name = 'mesh_fs2rs'; call push_sub()

  fac = M_z1
  if(present(factor)) fac = factor
  allocate(fr(m%fft_n(1), m%fft_n(2), m%fft_n(3)))   
#ifdef R_TREAL
  call rfftwnd_f77_one_complex_to_real(m%dplanb, g, fr)
#else
  call fftwnd_f77_one(m%zplanb, g, fr)
#endif

  call R_FUNC(scal)(m%fft_n(1)*m%fft_n(2)*m%fft_n(3), &
              fac/(m%fft_n(1)*m%fft_n(2)*m%fft_n(3)), fr, 1)
  f(1:m%np) = R_TOTYPE(M_ZERO)
  call R_FUNC(cube_to_mesh) (m, fr, f)

  deallocate(fr)
  g(1:m%R_FUNC(npw)) = g_aux(1:m%R_FUNC(npw))
  deallocate(g_aux)
  call pop_sub(); return
end subroutine R_FUNC(mesh_fs2rs)

! calculates the laplacian and the gradient of a function on the mesh
subroutine R_FUNC(mesh_derivatives) (m, f, lapl, grad, alpha)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(0:m%np)
  R_TYPE, intent(out), optional:: lapl(1:m%np), grad(3, 1:m%np)
  R_TYPE, intent(in), optional :: alpha

  R_TYPE :: alp

  sub_name = 'mesh_derivatives'; call push_sub()
  
  alp = R_TOTYPE(1._r8)
  if(present(alpha)) alp = alpha

  select case(m%d%space)
    case(REAL_SPACE)
#if defined(BOUNDARIES_ZERO_DERIVATIVE)
      call rs_derivative()
#else
      if(m%iso .and. present(lapl) .and. (.not.present(grad)) .and. conf%periodic_dim==0 ) then
        call rs_lapl_fast()
      else
        call rs_derivative()
      endif
#endif
    case(RECIPROCAL_SPACE)
      call fft_derivative()
  end select

  call pop_sub(); return
contains

  subroutine fft_derivative()
    complex(r8), allocatable :: fw1(:,:,:), fw2(:,:,:)
    
    integer :: n(3), nx, i, j

    n(:) = m%fft_n(:)
#ifdef R_TREAL
    nx = m%hfft_n
#else
    nx = n(1)
#endif

    allocate(fw1(nx, n(2), n(3)), fw2(nx, n(2), n(3)))
    call R_FUNC(mesh_rs2fs)(m, f(1:), fw1)
    if(present(grad)) then
      do i = 1, 3
        call mesh_gradient_in_FS(m, nx, n, fw1, fw2, i)
        call R_FUNC(mesh_fs2rs)(m, fw2, grad(i, :), factor = alp)
      end do
    end if
    if(present(lapl)) then
      call mesh_laplacian_in_FS(m, nx, n, fw1, fw2)
      call R_FUNC(mesh_fs2rs)(m, fw2, lapl, factor = alp)
    end if
    
    deallocate(fw1, fw2)
    return
  end subroutine fft_derivative
  
  subroutine rs_lapl_fast()
    R_TYPE :: sd
    integer :: k, in, ix, iy, iz, i

    do k = 1, m%np
      sd = (m%d%dlidfj(0)*f(k))*real(conf%dim, r8)
      do in = 1, m%d%norder
        sd = sd + sum(f(m%ind(1:2*conf%dim, in, k)))*m%d%dlidfj(in)
      end do
      lapl(k) = sd
    end do
    call R_FUNC(scal) (m%np, alp/m%h(1)**2, lapl(1), 1)
    
    return
  end subroutine rs_lapl_fast

  subroutine rs_derivative()
    R_TYPE :: den1(conf%dim), den2(conf%dim)
    integer :: k, in, ind1(conf%dim), ind2(conf%dim), ix, iy, iz, i
    
    do k = 1, m%np
      if(present(lapl)) then
        lapl(k) = alp*(m%d%dlidfj(0)*f(k))*sum(1/m%h(1:conf%dim)**2)
      end if
      if(present(grad)) then
        grad(:, k) = m%d%dgidfj(0)*f(k)
      end if
      
      ix = m%Lxyz(1, k); iy = m%Lxyz(2, k); iz = m%Lxyz(3, k)
      do in = 1, m%d%norder
        ind1(1) = m%Lxyz_inv(ix-in, iy, iz)
        ind2(1) = m%Lxyz_inv(ix+in, iy, iz)
        if(conf%dim > 1) then
          ind1(2) = m%Lxyz_inv(ix, iy-in, iz)
          ind2(2) = m%Lxyz_inv(ix, iy+in, iz)
          if(conf%dim > 2) then
            ind1(3) = m%Lxyz_inv(ix, iy, iz-in)
            ind2(3) = m%Lxyz_inv(ix, iy, iz+in)
          end if
        end if
        
        ! cyclic boundary conditions in the periodic direction(s)
        if (conf%periodic_dim>0) then
          if(ind1(1) == 0) then
            ind1(1) = m%Lxyz_inv(2*m%nr(1) + ix - in, iy, iz )
          end if
          if(ind2(1) == 0) then
            ind1(1) = m%Lxyz_inv(ix + in - 2*m%nr(1), iy, iz)
          end if
          if (conf%periodic_dim>1) then
            if(ind1(2) == 0) then
              ind1(2) = m%Lxyz_inv(ix, 2*m%nr(2) + iy - in, iz )
            end if
            if(ind2(2) == 0) then
              ind1(2) = m%Lxyz_inv(ix, iy + in - 2*m%nr(2), iz)
            end if
            if (conf%periodic_dim>2) then
             if(ind1(3) == 0) then
               ind1(3) = m%Lxyz_inv(ix, iy, 2*m%nr(3) + iz - in )
             end if
             if(ind2(3) == 0) then
               ind1(3) = m%Lxyz_inv(ix , iy, iz + in - 2*m%nr(3))
             end if
            end if
          end if
        end if
        
        ! If you prefer 0 wave functions at the boundary, uncomment the following
        ! Just be careful with the LB94 xc potential, for it will probably not work!
#ifndef BOUNDARIES_ZERO_DERIVATIVE
        den1(1:conf%dim) = f(ind1(1:conf%dim))
        den2(1:conf%dim) = f(ind2(1:conf%dim))
#else
        ! This also sets zero wavefunction
        ! den1 = 0._r8; den2 = 0._r8
        
        ! This peace of code changes the boundary conditions
        ! to have 0 derivative at the boundary
        do i = conf%periodic_dim+1, conf%dim
          if(ind1(i) > 0)den1(i) = f(ind1(i))
          if(ind2(i) > 0)den2(i) = f(ind2(i))
        end do
        if (conf%periodic_dim>0)
        message(1) = "Zero boundary conditions are not allowed in the periodic direction(s)"
        write(message(2),'(6x,a,i1,a,i1)'),"Zero boundary applied only in directions from ",conf%periodic_dim+1,' to ',conf%dim
        call write_warning(2)
#endif
        
        if(present(lapl)) then
          do i = 1, conf%dim
            lapl(k) = lapl(k) + alp*m%d%dlidfj(in)*((den1(i)+den2(i))/m%h(i)**2)
          end do
        end if
        
        if(present(grad)) then
          grad(1:conf%dim, k) = grad(1:conf%dim, k) + &
               m%d%dgidfj(-in)*den1(1:conf%dim) + m%d%dgidfj(in)*den2(1:conf%dim)
        end if
        
      end do
    end do
    
    if(present(grad)) then
      do i = 1, conf%dim
        grad(i,:) = grad(i,:) * (alp/m%h(i))
      enddo
    end if
    
    return
  end subroutine rs_derivative

end subroutine R_FUNC(mesh_derivatives)
