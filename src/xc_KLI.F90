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

subroutine X(xc_KLI_solve) (m, st, is, oep, oep_level)
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(in)    :: st
  type(xc_oep_type), intent(inout) :: oep
  integer,           intent(in)    :: oep_level, is

  integer :: i, j, n
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:)
  FLOAT, allocatable :: Ma(:,:), x(:,:), y(:,:)

  ! some intermediate quantities
  ! vxc contains the Slater part!
  allocate(rho_sigma(m%np))
  do i = 1, m%np
    rho_sigma(i) = max(sum(oep%socc*st%occ(:, is)*R_ABS(st%X(psi)(i, 1, :, is))**2), CNST(1e-20))
    oep%vxc(i)       = oep%socc* &
         sum(st%occ(:, is)*R_REAL(oep%lxc(i, :)*st%X(psi)(i, 1, :, is)))/rho_sigma(i)
  end do

  n = oep%eigen_n
  slater_approx: if(oep_level > 0) then

    allocate(v_bar_S(st%nst))
    do i = 1, st%nst
      if(st%occ(i, is) .gt. small) then
        v_bar_S(i) = sum(R_ABS(st%X(psi)(:, 1, i, is))**2 * oep%vxc(:))*m%vol_pp
      end if
    end do
    
    if(n > 0) then ! there is more than one state, so solve linear equation
      allocate(x(n, 1))
      x = M_ZERO
      allocate(Ma(n, n), y(n, 1))
      do i = 1, n
        do j = i, n
          Ma(i,j) = -sum(                                 &
               R_ABS(st%X(psi)(:, 1, oep%eigen_index(i), is))**2 * &
               R_ABS(st%X(psi)(:, 1, oep%eigen_index(j), is))**2 / &
               rho_sigma(:))*m%vol_pp
          Ma(j,i) = Ma(i,j)
        end do
        Ma(i,i) = 1 + Ma(i,i)
        
        y(i, 1) = v_bar_S(oep%eigen_index(i)) - oep%uxc_bar(oep%eigen_index(i))
      end do
      
      call dlinsyssolve(n, 1, Ma, y, x)
      deallocate(Ma, y)
      
      ! add contribution of low lying states
      do i = 1, n
        oep%vxc(:) = oep%vxc(:) + &
             oep%socc*st%occ(oep%eigen_index(i),is)* x(i,1) *  &
             R_ABS(st%X(psi)(:, 1, oep%eigen_index(i), is))**2 / rho_sigma(:)
      end do
      deallocate(x)
    end if
    
    deallocate(v_bar_S)
  end if slater_approx

  deallocate(rho_sigma)

end subroutine X(xc_KLI_solve)
