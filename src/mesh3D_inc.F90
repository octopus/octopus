! this functions returns the dot product between two vectors
! it uses BLAS
function R_FUNC(mesh_dp)(m, f1, f2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(1:m%np), f2(1:m%np)
  R_TYPE :: R_FUNC(mesh_dp)
  R_TYPE, external :: R_DOT
  
  R_FUNC(mesh_dp) = R_DOT(m%np, f1(1:m%np), 1,  f2(1:m%np), 1)*m%vol_pp
end function R_FUNC(mesh_dp)

! Conversion subroutines
! they actually add, do not forget it
subroutine R_FUNC(mesh_to_cube) (m, f_mesh, f_cube)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_mesh(m%np)
  R_TYPE, intent(inout) :: f_cube(m%fft_n, m%fft_n, m%fft_n)
  
  integer :: i, ix, iy, iz
  
  do i = 1, m%np
    ix = m%lx(i) + m%fft_n/2 + 1
    iy = m%Ly(i) + m%fft_n/2 + 1
    iz = m%Lz(i) + m%fft_n/2 + 1
    f_cube(ix, iy, iz) = f_cube(ix, iy, iz) + f_mesh(i)
  end do
end subroutine R_FUNC(mesh_to_cube)

subroutine R_FUNC(cube_to_mesh) (m, f_cube, f_mesh)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_cube(m%fft_n, m%fft_n, m%fft_n)
  R_TYPE, intent(inout) :: f_mesh(m%np)
  
  integer :: i, ix, iy, iz

  do i=1, m%np
    ix = m%lx(i) + m%fft_n/2 + 1
    iy = m%Ly(i) + m%fft_n/2 + 1
    iz = m%Lz(i) + m%fft_n/2 + 1
    f_mesh(i) = f_mesh(i) + f_cube(ix, iy, iz) 
  end do
end subroutine R_FUNC(cube_to_mesh)

! integrates a function on the mesh (could not find BLAS routine to do it ;))
function R_FUNC(mesh_integrate) (m, f)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)

  R_TYPE :: R_FUNC(mesh_integrate)

  R_FUNC(mesh_integrate) = sum(f(1:m%np))*m%vol_pp

end function R_FUNC(mesh_integrate)

! calculates the laplacian and the gradient of a function on the mesh
subroutine R_FUNC(mesh_derivatives) (m, f, lapl, grad)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(0:m%np)
  R_TYPE, intent(out), optional:: lapl(1:m%np), grad(3, 1:m%np)
  
  select case(m%d%space)
  case(0)
    call rs_derivative()
  case(1)
    call fft_derivative()
  end select

  return
contains

  subroutine fft_derivative()
    R_TYPE, allocatable :: fr(:,:,:)
    complex(r8), allocatable :: fw1(:,:,:), fw2(:,:,:)
    
    integer :: n, nx, i, j, ix, iy, iz, kx, ky, kz
    real(r8) :: temp, g2

    n = m%fft_n
#ifdef R_TREAL
    nx = n/2+1
#else
    nx = n
#endif

    allocate(fr(n, n, n), fw1(nx, n, n))
    fr = 0.0_r8
    call R_FUNC(mesh_to_cube) (m, f, fr)
    
#ifdef R_TREAL
    call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, fw1)
#else
    call fftwnd_f77_one(m%zplanf, fr, fw1)
#endif

    if(present(grad)) then
      allocate(fw2(nx, n, n))
      
      temp = (2.0_r8*M_Pi)/(n*m%h)
      do i = 1, 3
        fw2 = (0._r8, 0._r8)
        kx = 0; ky = 0; kz = 0
        do iz = 1, n
          if(i == 3) kz = pad_feq(iz, n, .true.)
          do iy = 1, n
            if(i == 2) ky = pad_feq(iy, n, .true.)
            do ix = 1, nx
              if(i == 1) kx = pad_feq(ix, n, .true.)
              
              g2 = temp*real(kx + ky + kz, r8)
              fw2(ix, iy, iz) = g2*M_zI*fw1(ix, iy, iz)
            end do
          end do
        end do
#ifdef R_TREAL
        call rfftwnd_f77_one_complex_to_real(m%dplanb, fw2, fr)
#else
        call fftwnd_f77_one(m%zplanb, fw2, fr)
#endif
        call R_FUNC(scal)(n**3, REALORCOMPLEX(1.0_r8/real(n, r8)**3), fr, 1)
        call R_FUNC(cube_to_mesh) (m, fr, grad(:, j))
      end do
      
      deallocate(fw2)
    end if
    
    if(present(lapl)) then
      temp = ((2.0_r8*M_Pi)/(n*m%h))**2
      do iz = 1, n
        kz = pad_feq(iz, n, .true.)
        do iy = 1, n
          ky = pad_feq(iy, n, .true.)
          do ix = 1, nx
            kx = pad_feq(ix, n, .true.)
            
            g2 = temp*(kx**2 + ky**2 + kz**2) 
            fw1(ix,iy,iz) = -g2*fw1(ix,iy,iz)
          end do
        end do
      end do
#ifdef R_TREAL
      call rfftwnd_f77_one_complex_to_real(m%dplanb, fw1, fr)
#else
      call fftwnd_f77_one(m%zplanb, fw1, fr)
#endif
      call R_FUNC(scal)(n**3, REALORCOMPLEX(1.0_r8/real(n, r8)**3), fr, 1)
      call R_FUNC(cube_to_mesh) (m, fr, lapl)
    end if
    
    deallocate(fr, fw1)
    return
  end subroutine fft_derivative

  subroutine rs_derivative()
    R_TYPE :: den1(3), den2(3)
    integer :: k, in, ind1(3), ind2(3), ix, iy, iz

    do k = 1, m%np
      ! first we add the 0 contribution
      den1 = f(k)
      den2 = f(k)
      if(present(lapl)) then
        lapl(k) = 3._r8*m%d%dlidfj(0)*f(k)
      end if
      if(present(grad)) then
        grad(:, k) = m%d%dgidfj(0)*f(k)
      end if

      ix = m%Lx(k); iy = m%Ly(k); iz = m%Lz(k)
      do in = 1, m%d%norder
        ind1(1) = m%Lxyz_inv(ix-in, iy, iz)
        ind1(2) = m%Lxyz_inv(ix, iy-in, iz)
        ind1(3) = m%Lxyz_inv(ix, iy, iz-in)
    
        ind2(1) = m%Lxyz_inv(ix+in, iy, iz)
        ind2(2) = m%Lxyz_inv(ix, iy+in, iz)
        ind2(3) = m%Lxyz_inv(ix, iy, iz+in)

        ! If you prefer 0 wave functions at the boundary, uncomment the following
        ! Just be careful with the LB94 xc potential, for it will probably not work!
#ifndef BOUNDARIES_ZERO_DERIVATIVE
        den1(:) = f(ind1(:))
        den2(:) = f(ind2(:))
#else
        ! This also sets zero wavefunction
        ! den1 = 0._r8; den2 = 0._r8

        ! This peace of code changes the boundary conditions
        ! to have 0 derivative at the boundary 
        if(ind1(1) > 0)den1(1) = f(ind1(1))
        if(ind1(2) > 0)den1(2) = f(ind1(2))
        if(ind1(3) > 0)den1(3) = f(ind1(3))
        
        if(ind2(1) > 0)den2(1) = f(ind2(1))
        if(ind2(2) > 0)den2(2) = f(ind2(2))
        if(ind2(3) > 0)den2(3) = f(ind2(3))
#endif
    
        if(present(lapl)) then
          lapl(k) = lapl(k) + m%d%dlidfj(in)*sum(den1(:)+den2(:))
        end if
        
        if(present(grad)) then
          grad(:, k) = grad(:, k) + m%d%dgidfj(-in)*den1(:) + m%d%dgidfj(in)*den2(:)
        end if
        
      end do
    end do

    if(present(lapl)) then
      lapl = lapl / (m%h**2)
    end if

    if(present(grad)) then
      grad = grad / m%h
    end if

    return
  end subroutine rs_derivative

end subroutine R_FUNC(mesh_derivatives)
