! this functions returns the dot product between two vectors
! it uses BLAS
R_TYPE function R_FUNC(mesh_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(1:m%np), f2(1:m%np)
  R_TYPE, external :: R_DOT
  
  dotp = R_DOT(m%np, f1(1), 1,  f2(1), 1)*m%vol_pp
end function R_FUNC(mesh_dotp)

! this functions returns the norm of a vector
! it uses BLAS
real(r8) function R_FUNC(mesh_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np, f, 1)*sqrt(m%vol_pp)
end function R_FUNC(mesh_nrm2)

subroutine R_FUNC(mesh_to_cube) (m, f_mesh, f_cube, t)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_mesh(m%np)
  R_TYPE, intent(inout) :: f_cube(:)
  integer, intent(in), optional :: t
  
  integer :: i, ix, n

  n = m%fft_n(1)/2+1
  if(present(t)) then
    if(t == 2) n = m%fft_n2(1)/2+1
  end if
 
  do i = 1, m%np
    ix = m%lx(i) + n
    f_cube(ix) = f_cube(ix) + f_mesh(i)
  end do
end subroutine R_FUNC(mesh_to_cube)

subroutine R_FUNC(cube_to_mesh) (m, f_cube, f_mesh, t)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN)    :: f_cube(:)
  R_TYPE, intent(inout) :: f_mesh(m%np)
  integer, intent(in), optional :: t

  integer :: i, ix, iy, iz, n

  n = m%fft_n(1)/2+1
  if(present(t)) then
    if(t == 2) n = m%fft_n2(1)/2+1
  end if

  do i = 1, m%np
    ix = m%lx(i) + n
    f_mesh(i) = f_mesh(i) + f_cube(ix) 
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
  case(REAL_SPACE)
    call rs_derivative()
  case(RECIPROCAL_SPACE)
    call fft_derivative()
  end select

  return
contains

  subroutine fft_derivative()
    R_TYPE, allocatable :: fr(:)
    complex(r8), allocatable :: fw1(:), fw2(:)
    
    integer :: n(3), nx, i, j, ix, iy, iz, kx, ky, kz
    real(r8) :: temp(3), g2

    n(1) = m%fft_n(1); n(2) = 0; n(3) = 0
#ifdef R_TREAL
    nx = n(1)/2+1
#else
    nx = n(1)
#endif

    allocate(fr(n(1)), fw1(nx))
    fr = R_TOTYPE(0.0_r8); fw1 = R_TOTYPE(0.0_r8)
    call R_FUNC(mesh_to_cube) (m, f(1:), fr)
    
#ifdef R_TREAL
    call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, fw1)
#else
    call fftwnd_f77_one(m%zplanf, fr, fw1)
#endif

    if(present(grad)) then
      allocate(fw2(nx))
        call mesh_gradient_in_FS(m, nx, n, fw1, fw2)
#ifdef R_TREAL
        call rfftwnd_f77_one_complex_to_real(m%dplanb, fw2, fr)
#else
        call fftwnd_f77_one(m%zplanb, fw2, fr)
#endif
        call R_FUNC(scal)(n(1), R_TOTYPE(1.0_r8/n(1)), fr, 1)
        grad(i, :) = R_TOTYPE(0._r8)
        call R_FUNC(cube_to_mesh) (m, fr, grad(i, :))
      deallocate(fw2)
    end if
    
    if(present(lapl)) then
      allocate(fw2(nx))
      call mesh_laplacian_in_FS(m, nx, n, fw1, fw2)
#ifdef R_TREAL
      call rfftwnd_f77_one_complex_to_real(m%dplanb, fw2, fr)
#else
      call fftwnd_f77_one(m%zplanb, fw2, fr)
#endif
      call R_FUNC(scal)(n(1), R_TOTYPE(1.0_r8/n(1)), fr, 1)
      lapl = R_TOTYPE(0._r8)
      call R_FUNC(cube_to_mesh) (m, fr, lapl)
    end if
    
    deallocate(fr, fw1)
    return
  end subroutine fft_derivative

  subroutine rs_derivative()
    R_TYPE :: den1(3), den2(3)
    integer :: k, in, ind1(3), ind2(3), ix, iy, iz, i

    do k = 1, m%np
      ! first we add the 0 contribution
      den1 = f(k)
      den2 = f(k)
      if(present(lapl)) then
        lapl(k) = (m%d%dlidfj(0)*f(k))*(1/m%h(1)**2)
      end if
      if(present(grad)) then
        grad(:, k) = m%d%dgidfj(0)*f(k)
      end if

      ix = m%Lx(k); iy = m%Ly(k); iz = m%Lz(k)
      do in = 1, m%d%norder
        ind1(1) = m%Lxyz_inv(ix-in, iy, iz); ind1(2) = 0; ind1(3) = 0
    
        ind2(1) = m%Lxyz_inv(ix+in, iy, iz); ind2(2) = 0; ind2(3) = 0

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
        if(ind2(1) > 0)den2(1) = f(ind2(1))
#endif
    
        if(present(lapl)) then
             lapl(k) = lapl(k) + m%d%dlidfj(in)*((den1(1)+den2(1))/m%h(1)**2)
        end if
        
        if(present(grad)) then
          grad(:, k) = grad(:, k) + m%d%dgidfj(-in)*den1(:) + m%d%dgidfj(in)*den2(:)
        end if
        
      end do
    end do

    if(present(grad)) then
         grad(i,:) = grad(i,:) / m%h(1)
    end if

    return
  end subroutine rs_derivative

end subroutine R_FUNC(mesh_derivatives)

subroutine R_FUNC(mesh_write_function) (m, f, u, filename)
  type(mesh_type), intent(IN)  :: m
  R_TYPE, intent(IN)           :: f(m%np)
  real(r8), intent(IN)         :: u
  character(len=*), intent(in) :: filename

  integer :: iunit, i

  call io_assign(iunit)
  open(unit=iunit, file=trim(filename), form='formatted')
  do i = 1, m%np
     write(iunit, '(f12.6, es20.12)') &
          m%lx(i)*m%h(1)/units_out%length%factor, f(i)/u
  enddo
  call io_close(iunit)

end subroutine R_FUNC(mesh_write_function)
