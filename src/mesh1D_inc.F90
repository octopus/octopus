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

! calculates the laplacian and the gradient of a function on the mesh
subroutine R_FUNC(mesh_derivatives) (m, f, lapl, grad, alpha)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(0:m%np)
  R_TYPE, intent(out), optional:: lapl(1:m%np), grad(3, 1:m%np)
  R_TYPE, intent(in), optional :: alpha

  real(r8) :: alp

  alp = 1._r8
  if(present(alpha)) alp = alpha

  call rs_derivative()

  return
contains

  subroutine rs_derivative()
    R_TYPE :: den1(3), den2(3)
    integer :: k, in, ind1(3), ind2(3), ix, iy, iz, i
    
    do k = 1, m%np
      ! first we add the 0 contribution
      den1 = f(k)
      den2 = f(k)
      if(present(lapl)) then
        lapl(k) = (m%d%dlidfj(0)*f(k))
      end if
      if(present(grad)) then
        grad(1, k) = m%d%dgidfj(0)*f(k)
      end if
      
      ix = m%Lxyz(1, k)
      do in = 1, m%d%norder
        ind1(1) = m%Lxyz_inv(ix-in, 0, 0)
        ind2(1) = m%Lxyz_inv(ix+in, 0, 0)
        
        ! If you prefer 0 wave functions at the boundary, uncomment the following
        ! Just be careful with the LB94 xc potential, for it will probably not work!
#ifndef BOUNDARIES_ZERO_DERIVATIVE
        den1(1) = f(ind1(1))
        den2(1) = f(ind2(1))
#else
        ! This also sets zero wavefunction
        ! den1 = 0._r8; den2 = 0._r8
        ! This peace of code changes the boundary conditions
        ! to have 0 derivative at the boundary 
        if(ind1(1) > 0)den1(1) = f(ind1(1))
        if(ind2(1) > 0)den2(1) = f(ind2(1))
#endif
    
        if(present(lapl)) then
             lapl(k) = lapl(k) + m%d%dlidfj(in)*(den1(1)+den2(1))
        end if
        
        if(present(grad)) then
          grad(1, k) = grad(1, k) + m%d%dgidfj(-in)*den1(1) + m%d%dgidfj(in)*den2(1)
        end if
        
      end do
    end do

    if(present(grad)) then
         grad(1,:) = grad(1,:) / m%h(1)
         grad(2:3,:) = 0.0_r8
    end if
    if(present(lapl)) then
         lapl = lapl * alp / m%h(1)**2
    endif

    return
  end subroutine rs_derivative

end subroutine R_FUNC(mesh_derivatives)
