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
  R_TYPE, intent(IN) :: f(m%np)
  R_TYPE :: d

  d = sum(f(1:m%np))*m%vol_pp
end function X(mf_integrate)

!!! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(m%np), f2(m%np)

  dotp = X(lalg_dot)(m%np, f1(1),  f2(1))*m%vol_pp

end function X(mf_dotp)

!!! this function returns the norm of a vector
FLOAT function X(mf_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)

  nrm2 = X(lalg_nrm2)(m%np, f(1))*sqrt(m%vol_pp)

end function X(mf_nrm2)

!!! This function calculates the x_i moment of the function f
function X(mf_moment) (m, f, i, n) result(r)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)          :: f(1:m%np)
  integer, intent(in)         :: i, n
  R_TYPE                      :: r

  r = sum(f(1:m%np)*m%Lxyz(i, 1:m%np)**n) * m%h(i)**n * m%vol_pp
end function X(mf_moment)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following set of subroutines calculates derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! calculates the laplacian of a function on the mesh.
subroutine X(mf_laplacian) (m, f, lapl)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: lapl(m%np)
  integer :: k, nl
  call push_sub("mf_laplacian")
  nl = m%lapl%n
  do k = 1, m%np
     lapl(k) = sum(m%lapl%w(1:nl, k)*f(m%lapl%i(1:nl, k)))
  end do
  call pop_sub()
end subroutine X(mf_laplacian)

!!! applies a low-frequency filter to a function in a mesh.
subroutine X(mf_filter) (m, filter, f, filteredf)
  type(mesh_type), intent(in) :: m
  type(derivatives_type), intent(in) :: filter
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: filteredf(m%np)
  integer :: k, nl
  call push_sub("mf_filter")

  do k = 1, m%np
     filteredf(k) = sum(filter%w(1:filter%n, k)*f(filter%i(1:filter%n, k)))
  end do

  call pop_sub()
end subroutine X(mf_filter)

subroutine X(mf_gradient) (m, f, grad)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: grad(conf%dim, m%np)
  integer :: j, k, ng(3)
  call push_sub("mf_gradient")
  do j = 1, conf%dim
     ng(j) = m%grad(j)%n
     do k = 1, m%np
        grad(j, k) = sum(m%grad(j)%w(1:ng(j), k)*f(m%grad(j)%i(1:ng(j), k)))
     enddo
  enddo
  call pop_sub()
end subroutine X(mf_gradient)

!!! Calculates the divergence of a vectorial function f.
!!! Currently it only does so in real space.
subroutine X(mf_divergence)(m, f, divf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(conf%dim, m%np)
  R_TYPE, intent(out) :: divf(m%np)
  integer :: j, k, ng(3)
  call push_sub('mf_divergence')

  divf = R_TOTYPE(M_ZERO)
  do j = 1, conf%dim
     ng(j) = m%grad(j)%n
     do k = 1, m%np
        divf(k) = divf(k) + sum(m%grad(j)%w(1:ng(j), k)*f(j, m%grad(j)%i(1:ng(j), k)))
     end do
  end do

  call pop_sub()
end subroutine X(mf_divergence)

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
  call X(lalg_scal)(m%np, R_TOTYPE(M_ONE/r), f) 

  call pop_sub()
end subroutine X(mf_random)
