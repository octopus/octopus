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

! integrates a function on the mesh (could not find BLAS routine to do it ;))
function R_FUNC(mesh_integrate) (m, f)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(1:m%np)

  R_TYPE :: R_FUNC(mesh_integrate)

  R_FUNC(mesh_integrate) = sum(f(1:m%np))*m%vol_pp

end function R_FUNC(mesh_integrate)

#if defined(ONE_D)
#  include "mesh1D_inc.F90"
#elif defined(THREE_D)
#  include "mesh3D_inc.F90"
#endif
