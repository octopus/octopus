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

! Conversion subroutines
! they actually add, do not forget it
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
    ix = m%lx(i) + n(1)
    iy = m%Ly(i) + n(2)
    iz = m%Lz(i) + n(3)
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
    ix = m%lx(i) + n(1)
    iy = m%Ly(i) + n(2)
    iz = m%Lz(i) + n(3)
    f_mesh(i) = f_mesh(i) + f_cube(ix, iy, iz) 
  end do
end subroutine R_FUNC(cube_to_mesh)

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
#if defined(POLYMERS) || defined(BOUNDARIES_ZERO_DERIVATIVE)
      call rs_derivative()
#else
      if(m%iso .and. present(lapl) .and. (.not.present(grad)) ) then
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
    R_TYPE, allocatable :: fr(:,:,:)
    complex(r8), allocatable :: fw1(:,:,:), fw2(:,:,:)
    
    integer :: n(3), nx, i, j

    do i=1, 3
       n(i) = m%fft_n(i)
    enddo
#ifdef R_TREAL
    nx = n(1)/2 + 1
#else
    nx = n(1)
#endif

    allocate(fr(n(1), n(2), n(3)), fw1(nx, n(2), n(3)))
    fr = 0.0_r8; fw1 = R_TOTYPE(0._r8)
    call R_FUNC(mesh_to_cube) (m, f(1:), fr)
    
#ifdef R_TREAL
    call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, fw1)
#else
    call fftwnd_f77_one(m%zplanf, fr, fw1)
#endif

! I am not sure that implementation of the gradient is correct;
! should be checked.
    if(present(grad)) then
      allocate(fw2(nx, n(2), n(3)))
      
      do i = 1, 3
        call mesh_gradient_in_FS(m, nx, n, fw1, fw2, i)
#ifdef R_TREAL
        call rfftwnd_f77_one_complex_to_real(m%dplanb, fw2, fr)
#else
        call fftwnd_f77_one(m%zplanb, fw2, fr)
#endif

        call R_FUNC(scal)(n(1)*n(2)*n(3), alp/(n(1)*n(2)*n(3)), fr, 1)
        grad(i, :) = R_TOTYPE(0._r8)
        call R_FUNC(cube_to_mesh) (m, fr, grad(i, :))
      end do
      
      deallocate(fw2)
    end if
    
    if(present(lapl)) then
      allocate(fw2(nx, n(2), n(3)))

      call mesh_laplacian_in_FS(m, nx, n, fw1, fw2)
#ifdef R_TREAL
      call rfftwnd_f77_one_complex_to_real(m%dplanb, fw2, fr)
#else
      call fftwnd_f77_one(m%zplanb, fw2, fr)
#endif

      call R_FUNC(scal)(n(1)*n(2)*n(3), alp/(n(1)*n(2)*n(3)), fr, 1)
      lapl = R_TOTYPE(0._r8)
      call R_FUNC(cube_to_mesh) (m, fr, lapl)

      deallocate(fw2)
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
        lapl(k) = (m%d%dlidfj(0)*f(k))*sum(1/m%h(:)**2)
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

#ifdef POLYMERS
        ! cyclic boundary conditions in the z direction
        if(ind1(3) == 0) then
          ind1(3) = m%Lxyz_inv(ix, iy, 2*m%nr(3) + iz - in)
        end if
        if(ind2(3) == 0) then
          ind1(3) = m%Lxyz_inv(ix, iy, iz + in - 2*m%nr(3))
        end if
#endif

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
          do i = 1, 3
             lapl(k) = lapl(k) + alp*m%d%dlidfj(in)*((den1(i)+den2(i))/m%h(i)**2)
          enddo
        end if
        
        if(present(grad)) then
          grad(:, k) = grad(:, k) + m%d%dgidfj(-in)*den1(:) + m%d%dgidfj(in)*den2(:)
        end if
        
      end do
    end do

    if(present(grad)) then
      do i = 1, 3
         grad(i,:) = grad(i,:) * (alp/m%h(i))
      enddo
    end if

    return
  end subroutine rs_derivative

  subroutine rs_lapl_fast()
    R_TYPE :: sd
    integer :: k, in, ix, iy, iz, i

    do k = 1, m%np
      sd = (m%d%dlidfj(0)*f(k))*3.0_r8
      do in = 1, m%d%norder
        sd = sd + sum(f(m%ind(:, in, k)))*m%d%dlidfj(in)
      end do
      lapl(k) = sd
    end do
    call R_FUNC(scal) (m%np, alp/m%h(1)**2, lapl(1), 1)

    return
  end subroutine rs_lapl_fast

end subroutine R_FUNC(mesh_derivatives)
