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
!!! (could not find BLAS routine to do it ;))
function X(mf_integrate) (m, f) result(d)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)
  R_TYPE :: d

  d = sum(f(1:m%np))*m%vol_pp
end function X(mf_integrate)

!!! this function returns the dot product between two vectors
!!! if HAVE_BLAS is defined it uses the BLAS library
R_TYPE function X(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f1(m%np), f2(m%np)

#if defined(HAVE_BLAS)
  real(r8), external :: ddot

  dotp = R_DOT(m%np, f1(1), 1,  f2(1), 1)*m%vol_pp
#else
  dotp = sum(R_CONJ(f1)*f2)*m%vol_pp
#endif

end function X(mf_dotp)

!!! this function returns the norm of a vector
!!! if HAVE_BLAS is defined it uses the BLAS library
real(r8) function X(mf_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)

#if defined(HAVE_BLAS)
!  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np, f(1), 1)*sqrt(m%vol_pp)
#else
  nrm2 = sqrt(sum(R_CONJ(f)*f)*m%vol_pp)
#endif

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
  type(lookup_type), pointer :: p

  call push_sub("mf_laplacian")
    
  do k = 1, m%np
    p => m%laplacian%lookup(k)
    nl = p%n
    lapl(k) = sum(p%w(1:nl)*f(p%i(1:nl)))
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
  type(lookup_type), pointer :: p

  call push_sub("mf_filter")

  do k = 1, m%np
    p => filter%lookup(k)
    nl = p%n
    filteredf(k) = sum(p%w(1:nl)*f(p%i(1:nl)))
  end do

  call pop_sub()
end subroutine X(mf_filter)

subroutine X(mf_gradient) (m, f, grad)
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(in) :: f(m%np)
  R_TYPE, intent(out) :: grad(conf%dim, m%np)
  
  integer :: j, k, ng(3)
  type(lookup_type), pointer :: p

  call push_sub("mf_gradient")
    
  do k = 1, m%np
    do j = 1, conf%dim
       p => m%grad(j)%lookup(k)
       ng = p%n
       grad(j, k) = sum(p%w(1:ng(j))*f(p%i(1:ng(j))))
    end do
  end do

  call pop_sub()
end subroutine X(mf_gradient)

!!! Calculates the divergence of a vectorial function f.
!!! Currently it only does so in real space.
subroutine X(mf_divergence)(m, f, divf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(m%np, conf%dim)
  R_TYPE, intent(out) :: divf(m%np)

  integer :: j, k, ng(3)
  type(lookup_type), pointer :: p

  call push_sub('mf_divergence')

  divf = R_TOTYPE(M_ZERO)

  do j = 1, conf%dim
     do k = 1, m%np
        p => m%grad(j)%lookup(k)
        ng = p%n
        divf(k) = divf(k) + sum(p%w(1:ng(j))*f(p%i(1:ng(j)), j))
     end do
  end do

  call pop_sub()
end subroutine X(mf_divergence)
