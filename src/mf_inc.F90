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
!!
!! $Id$


!!!  integrates a function
R_TYPE function X(mf_integrate) (m, f) result(d)
  type(mesh_type), intent(in) :: m
  R_TYPE,          intent(in) :: f(1:m%np)  ! f(m%np)

  call push_sub('mf_inc.Xmf_integrate')

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  d = X(vec_integrate)(m%vp, f(1:m%np)*m%vol_pp(1:m%np))
#else 
  d = sum(f(1:m%np)*m%vol_pp(1:m%np))
#endif

  call pop_sub()

end function X(mf_integrate)


!!! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp)(m, f1, f2) result(dotp)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in) :: f1(1:m%np), f2(1:m%np) ! f(m%np)

  R_TYPE, allocatable :: l(:)
  R_TYPE              :: dotp_tmp

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  integer :: ierr
#endif

  call push_sub('mf_inc.Xmf_dotp')

  ! This is not implemented via vec_integrate
  ! because BLAS is much faster.
  if(m%use_curvlinear) then
    allocate(l(m%np))
    l(1:m%np) = f1(1:m%np) * m%vol_pp(1:m%np)
    dotp_tmp  = lalg_dot(m%np, l(:),  f2(:))
    deallocate(l)
  else
    dotp_tmp = lalg_dot(m%np, f1(:),  f2(:))*m%vol_pp(1)
  end if
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  call TS(MPI_Allreduce)(dotp_tmp, dotp, 1, R_MPITYPE, &
                     MPI_SUM, m%vp%comm, ierr)
#else
  dotp = dotp_tmp
#endif

  call pop_sub()

end function X(mf_dotp)


!!! this function returns the norm of a vector
FLOAT function X(mf_nrm2)(m, f) result(nrm2)
  type(mesh_type), intent(in) :: m
  R_TYPE,          intent(in) :: f(1:m%np) ! f(m%np)

  call push_sub('mf_inc.Xmf_dotp')
  
  nrm2 = sqrt(X(mf_dotp) (m, f, f))
  
  call pop_sub()

end function X(mf_nrm2)


!!! This function calculates the x_i moment of the function f
function X(mf_moment) (m, f, i, n) result(r)
  type(mesh_type), intent(in) :: m
  R_TYPE,          intent(in) :: f(1:m%np)  ! f(m%np)
  integer,         intent(in) :: i, n
  R_TYPE                      :: r

  call push_sub('mf_inc.Xmf_dotp')

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  r = X(vec_integrate)(m%vp, f(1:m%np)*m%x(1:m%np,i)**n * m%vol_pp(1:m%np))
#else
  r = sum(f(1:m%np)*m%x(1:m%np,i)**n * m%vol_pp(1:m%np))
#endif
  
  call pop_sub()

end function X(mf_moment)


!!! This subroutine generates a gaussian wave-function in a random
!!! position in space
subroutine X(mf_random)(m, f)
  type(mesh_type),      intent(in)  :: m
  R_TYPE,               intent(out) :: f(1:m%np)

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(3), rnd, r

  call push_sub('mf_inc.states_random')

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
