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


! ---------------------------------------------------------
! integrates a function
R_TYPE function X(mf_integrate) (mesh, f) result(d)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,          intent(in) :: f(:)  ! f(mesh%np)

  call profiling_in(C_PROFILING_MF_INTEGRATE)
  call push_sub('mf_inc.Xmf_integrate')

  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    d = X(vec_integrate)(mesh%vp, f(1:mesh%np)*mesh%vol_pp(1:mesh%np))
#else
    ASSERT(.false.)
#endif
  else
    d = sum(f(1:mesh%np)*mesh%vol_pp(1:mesh%np))
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_MF_INTEGRATE)

end function X(mf_integrate)


! ---------------------------------------------------------
! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp)(mesh, f1, f2) result(dotp)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,          intent(in) :: f1(:), f2(:)

  R_TYPE, allocatable :: l(:)
  R_TYPE              :: dotp_tmp

#if defined(HAVE_MPI)
  integer :: ierr
#endif

  call profiling_in(C_PROFILING_MF_DOTP)
  call push_sub('mf_inc.Xmf_dotp')

  ! This is not implemented via vec_integrate
  ! because BLAS is much faster.
  if(mesh%use_curvlinear) then
    ALLOCATE(l(mesh%np), mesh%np)
    l(1:mesh%np) = f1(1:mesh%np) * mesh%vol_pp(1:mesh%np)
    dotp_tmp  = lalg_dot(mesh%np, l(:),  f2(:))
    deallocate(l)
  else
    dotp_tmp = lalg_dot(mesh%np, f1(:),  f2(:))*mesh%vol_pp(1)
  end if

  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    call profiling_in(C_PROFILING_MF_DOTP_ALLREDUCE)
    call TS(MPI_Allreduce)(dotp_tmp, dotp, 1, R_MPITYPE, &
      MPI_SUM, mesh%vp%comm, ierr)
    call profiling_out(C_PROFILING_MF_DOTP_ALLREDUCE)
#else
    ASSERT(.false.)
#endif
  else
    dotp = dotp_tmp
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_MF_DOTP)

end function X(mf_dotp)


! ---------------------------------------------------------
! this function returns the norm of a vector
FLOAT function X(mf_nrm2)(m, f) result(nrm2)
  type(mesh_t), intent(in) :: m
  R_TYPE,          intent(in) :: f(:)

  call profiling_in(C_PROFILING_MF_NRM2)
  call push_sub('mf_inc.Xmf_nrm2')

  nrm2 = sqrt(X(mf_dotp) (m, f, f))

  call pop_sub()
  call profiling_out(C_PROFILING_MF_NRM2)

end function X(mf_nrm2)


! ---------------------------------------------------------
! This function calculates the x_i moment of the function f
function X(mf_moment) (m, f, i, n) result(r)
  type(mesh_t), intent(in) :: m
  R_TYPE,          intent(in) :: f(1:m%np)  ! f(m%np)
  integer,         intent(in) :: i, n
  R_TYPE                      :: r

  call push_sub('mf_inc.Xmf_moment')

  if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
    r = X(vec_integrate)(m%vp, f(1:m%np)*m%x(1:m%np,i)**n * m%vol_pp(1:m%np))
#else
    ASSERT(.false.)
#endif
  else
    r = sum(f(1:m%np) * m%x(1:m%np,i)**n * m%vol_pp(1:m%np))
  end if

  call pop_sub()

end function X(mf_moment)


! ---------------------------------------------------------
! This subroutine generates a gaussian wave-function in a
! random position in space
subroutine X(mf_random)(m, f)
  type(mesh_t),      intent(in)  :: m
  R_TYPE,               intent(out) :: f(1:m%np)

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(3), rnd, r

  call push_sub('mf_inc.Xmf_random')

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

! ---------------------------------------------------------
! Does a partial integration of a function, i.e. in N
! dimensions, it integrates over N-1 dimensions. The
! argument j is the dimension which is not integrated.
! The input function f is a normal function defined on
! the mesh, whereas the output function is a 1D array
! of mesh%l(j) + 2*mesh%enlarge(j) elements.
!
! In order to retrieve the coordinate v of u(n), one needs
! to do:
!  k = gr%m%nr(1, j) + n - 1
!  i = gr%m%lxyz_inv(k, 0, 0) ! Assuming j = 1
!  if(i>0) v = gr%m%x(i, j)   ! I do not understand why
!                             ! we need i>0...
!
! WARNING: It will stop if one is using curvilinear
! coordinates, or real-space domain parallelization.
subroutine X(mf_partial_integrate)(mesh, j, f, u)
  type(mesh_t), intent(in)  :: mesh
  integer,      intent(in)  :: j
  R_TYPE,       intent(in)  :: f(:)
  R_TYPE,       intent(out) :: u(:)

  integer :: i, k, m
  call push_sub('mf_inc.mf_partial_integrate')

  ASSERT(.not.(mesh%parallel_in_domains))
  ASSERT(.not.(mesh%use_curvlinear))

  k = mesh%l(j)+2*mesh%enlarge(j)

  u(1:k) = M_ZERO
  do i = 1, mesh%np
     m = mesh%lxyz(i, j) - mesh%nr(1, j) + 1
     u(m) = u(m) + f(i)
  end do

  u(1:k) = u(1:k) * mesh%vol_pp(1)/mesh%h(j)

  call pop_sub()
end subroutine X(mf_partial_integrate)
