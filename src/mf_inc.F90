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
function X(mf_integrate) (m, f) result(d)
  type(mesh_type), intent(IN) :: m
  R_TYPE,          intent(IN) :: f(:)  ! f(m%np)
  R_TYPE                      :: d

  d = sum(f(1:m%np)*m%vol_pp(1:m%np))
end function X(mf_integrate)

!!! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f1(:), f2(:) ! f(m%np)

  R_TYPE, allocatable :: l(:)

  if(m%use_curvlinear) then
    allocate(l(m%np))
    l(1:m%np) = f1(1:m%np) * m%vol_pp(1:m%np)
    dotp = lalg_dot(m%np, l(:),  f2(:))
    deallocate(l)
  else
    dotp = lalg_dot(m%np, f1(:),  f2(:))*m%vol_pp(1)
  end if

end function X(mf_dotp)

!!! this function returns the norm of a vector
FLOAT function X(mf_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE,          intent(IN) :: f(:) ! f(m%np)

  if(m%use_curvlinear) then
   nrm2 = sqrt(X(mf_dotp) (m, f, f))
  else
    nrm2 = lalg_nrm2(m%np, f(:))*sqrt(m%vol_pp(1))
  end if

end function X(mf_nrm2)

!!! This function calculates the x_i moment of the function f
function X(mf_moment) (m, f, i, n) result(r)
  type(mesh_type), intent(in) :: m
  R_TYPE,          intent(in) :: f(:)  ! f(m%np)
  integer,         intent(in) :: i, n
  R_TYPE                      :: r

  r = sum(f(:)*m%x(:,i)**n * m%vol_pp(:))
end function X(mf_moment)


!!! This subroutine generates a gaussian wave-function in a random
!!! position in space
subroutine X(mf_random)(m, f)
  type(mesh_type), intent(IN)  :: m
  R_TYPE,          intent(out) :: f(:)

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(3), rnd, r

  call push_sub('states_random')

  call quickrnd(iseed, rnd)
  a(1) = M_TWO*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(2) = M_TWO*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(3) = M_TWO*(2*rnd - 1)

  do i = 1, m%np
     call mesh_r(m, i, r, a=a)
     f(i) = exp(-M_HALF*r*r)
  end do

  r = X(mf_nrm2)(m, f)
  call lalg_scal(m%np, R_TOTYPE(M_ONE/r), f) 

  call pop_sub()
end subroutine X(mf_random)
