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

  d = sum(f(:)*m%vol_pp(:))
end function X(mf_integrate)

!!! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(:), f2(:) ! f(m%np)

  R_TYPE, allocatable :: l(:)

  if(m%use_curvlinear) then

    allocate(l(m%np))
    l(:) = f1(:) * m%vol_pp(:)
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

  r = sum(f(1:m%np)*m%x(i, 1:m%np)**n * m%vol_pp(1:m%np))
end function X(mf_moment)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following set of subroutines calculates derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! calculates the laplacian of a function on the mesh.
subroutine X(mf_laplacian) (m, f, lapl)
  type(mesh_type), intent(in)  :: m
  R_TYPE,          intent(in)  :: f(:)     ! f(m%np)
  R_TYPE,          intent(out) :: lapl(:)  ! lapl(m%np)

  integer :: j, k, nl
  R_TYPE, allocatable :: l1(:,:), l2(:,:), t(:)

  call push_sub("mf_laplacian")

  if(m%use_curvlinear) then ! hack
    allocate(l1(conf%dim, m%np), l2(conf%dim, m%np), t(m%np))
    call X(mf_gradient) (m, f, l1)

    lapl = R_TOTYPE(M_ZERO)
    do j = 1, conf%dim
      t(:) = l1(j, :)
      call X(mf_gradient) (m, t, l2)
      lapl = lapl + l2(j,:)

!      do k = 1, m%np
!        print '(i3,10f12.6)', k, m%x(1, k), f(k), t(k), l2(1,k)
!      end do
    end do
    deallocate(l1, l2, t)

!    stop
  else
    nl = m%lapl%n
    do k = 1, m%np
      lapl(k) = sum(m%lapl%w(1:nl, k)*f(m%lapl%i(1:nl, k)))
    end do
  end if

  call pop_sub()
end subroutine X(mf_laplacian)

!!! applies a low-frequency filter to a function in a mesh.
subroutine X(mf_filter) (m, filter, f, filteredf)
  type(mesh_type),        intent(in)  :: m
  type(derivatives_type), intent(in)  :: filter
  R_TYPE,                 intent(in)  :: f(:)          ! (m%np)
  R_TYPE,                 intent(out) :: filteredf(:)  ! (m%np)

  integer :: k, nl
  call push_sub("mf_filter")

  do k = 1, m%np
    filteredf(k) = sum(filter%w(1:filter%n, k)*f(filter%i(1:filter%n, k)))
  end do

  call pop_sub()
end subroutine X(mf_filter)

subroutine X(mf_gradient) (m, f, grad)
  type(mesh_type), intent(IN)  :: m
  R_TYPE,          intent(in)  :: f(:)       ! f(m%np)
  R_TYPE,          intent(out) :: grad(:,:)  ! grad(conf%dim, m%np)

  integer :: i, j, k, ng
  FLOAT   :: chi(conf%dim), Jac(conf%dim,conf%dim)

  call push_sub("mf_gradient")

  do j = 1, conf%dim
    ng = m%grad(j)%n
    do k = 1, m%np
      grad(j, k) = sum(m%grad(j)%w(1:ng, k)*f(m%grad(j)%i(1:ng, k)))
    end do
  end do

  if(m%use_curvlinear) then
    do k = 1, m%np
      call jacobian(m%geo, m%x(1:conf%dim,k), chi(:), Jac(:,:))
      grad(1:conf%dim, k) = matmul(grad(1:conf%dim,k), Jac(:,:))
    end do
  end if

  call pop_sub()
end subroutine X(mf_gradient)

!!! Calculates the divergence of a vectorial function f.
!!! Currently it only does so in real space.
subroutine X(mf_divergence)(m, f, divf)
  type(mesh_type), intent(in)  :: m
  R_TYPE,          intent(in)  :: f(:,:)   ! (conf%dim, m%np)
  R_TYPE,          intent(out) :: divf(:)  ! (m%np)

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
  call lalg_scal(m%np, R_TOTYPE(M_ONE/r), f) 

  call pop_sub()
end subroutine X(mf_random)
