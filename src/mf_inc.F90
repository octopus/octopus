!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$


! ---------------------------------------------------------
! integrates a function
R_TYPE function X(mf_integrate) (mesh, f) result(d)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: f(:)  ! f(mesh%np)

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
  R_TYPE,       intent(in) :: f1(:), f2(:)

  R_TYPE, allocatable :: l(:)
  R_TYPE              :: dotp_tmp

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
    call MPI_Allreduce(dotp_tmp, dotp, 1, R_MPITYPE, &
      MPI_SUM, mesh%vp%comm, mpi_err)
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
  R_TYPE,       intent(in) :: f(:)

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
  R_TYPE,       intent(in) :: f(1:m%np)  ! f(m%np)
  integer,      intent(in) :: i, n

  R_TYPE                   :: r

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
  type(mesh_t), intent(in)  :: m
  R_TYPE,       intent(out) :: f(1:m%np)

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(MAX_DIM), rnd, r

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


! ---------------------------------------------------------
subroutine X(mf_interpolate) (mesh_in, mesh_out, full_interpolation, u, f)
  type(mesh_t), intent(in)    :: mesh_in, mesh_out
  logical,      intent(in)    :: full_interpolation
  R_TYPE,       intent(inout) :: u(:)    ! u(mesh_in%np_global)
  R_TYPE,       intent(out)   :: f(:)    ! f(mesh%np)

  FLOAT :: xmin, ymin, dx, dy, rmax, px, py, pz, xyzmin(MAX_DIM), xyzdel(MAX_DIM)
  FLOAT, allocatable :: rsq(:), a(:, :)
#if !defined(R_TREAL)
  FLOAT :: ip
  FLOAT, allocatable :: aux_u(:)
#endif  
  R_TYPE, allocatable :: f_global(:)
  integer, allocatable :: lcell(:, :, :), lnext(:)
  integer :: i, j, k, nq, nw, nr, ier, npoints, ix, iy, iz

  FLOAT, external :: qs2val, qs3val

  call push_sub('mf_inc.Xmf_interpolate2d')

  if(full_interpolation) then

    npoints = mesh_in%np_global

    select case(mesh_in%sb%dim)
    case(2)
      nq = 13
      nw = 19
      nr = nint(sqrt(npoints/CNST(3.0)))
      ALLOCATE(lcell(nr, nr, 1), nr*nr)
      ALLOCATE(a(5, npoints), 5*npoints)
    case(3)
      nq = 17 ! This is the recommended value in qshep3d.f90
      nw = 16 ! The recommended value in qshep3d.f90 is 32, but this speeds up things.
      nr = nint((npoints/CNST(3.0))**(M_ONE/M_THREE))
      ALLOCATE(lcell(nr, nr, nr), nr*nr*nr)
      ALLOCATE(a(9, npoints), 9*npoints)
    case(1)
      stop 'Believe it or not, cannot do 1D interpolation, only 2D or 3D.'
    end select
    ALLOCATE(lnext(npoints), npoints)
    ALLOCATE(rsq(npoints), npoints)


    f = M_ZERO

#if defined(R_TREAL)
    select case(mesh_in%sb%dim)
    case(2)
      call qshep2 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                    u, nq, nw, nr, lcell(:, :, 1), lnext, xmin, ymin, &
                    dx, dy, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        f(i) = qs2val ( px, py, npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                        u, nr, lcell(:, :, 1), lnext, xmin, &
                        ymin, dx, dy, rmax, rsq, a )
      end do
    case(3)
      call qshep3 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                    u, nq, nw, nr, lcell, lnext, xyzmin, &
                    xyzdel, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        pz = mesh_out%x(i, 3)
        f(i) = qs3val ( px, py, pz, npoints, &
                        mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                        u, nr, lcell, lnext, xyzmin, xyzdel, rmax, rsq, a )
      end do
    case(1)
      stop 'Believe it or not, cannot do 1D interpolation, only 2D or 3D.'
    end select
#else
    ALLOCATE(aux_u(npoints), npoints)
    select case(mesh_in%sb%dim)
    case(2)
      aux_u = R_REAL(u)
      call qshep2 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                    aux_u, nq, nw, nr, lcell(:, :, 1), lnext, xmin, ymin, &
                    dx, dy, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        f(i) = qs2val ( px, py, npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                        aux_u, nr, lcell(:, :, 1), lnext, xmin, &
                        ymin, dx, dy, rmax, rsq, a )
      end do
      aux_u = R_AIMAG(u)
      call qshep2 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                    aux_u, nq, nw, nr, lcell(:, :, 1), lnext, xmin, ymin, &
                    dx, dy, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        ip = qs2val ( px, py, npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), &
                      aux_u, nr, lcell(:, :, 1), lnext, xmin, &
                      ymin, dx, dy, rmax, rsq, a )
        f(i) = f(i) + M_zI*ip
      end do

    case(3)
      aux_u = R_REAL(u)
      call qshep3 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                    aux_u, nq, nw, nr, lcell, lnext, xyzmin, &
                    xyzdel, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        pz = mesh_out%x(i, 3)
        f(i) = qs3val ( px, py, pz, npoints, &
                        mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                        aux_u, nr, lcell, lnext, xyzmin, xyzdel, rmax, rsq, a )
      end do
      aux_u = R_AIMAG(u)
      call qshep3 ( npoints, mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                    aux_u, nq, nw, nr, lcell, lnext, xyzmin, &
                    xyzdel, rmax, rsq, a, ier )
      do i = 1, mesh_out%np
        px = mesh_out%x(i, 1)
        py = mesh_out%x(i, 2)
        pz = mesh_out%x(i, 3)
        ip = qs3val ( px, py, pz, npoints, &
                      mesh_in%x(1:npoints, 1), mesh_in%x(1:npoints, 2), mesh_in%x(1:npoints, 3), &
                      aux_u, nr, lcell, lnext, xyzmin, xyzdel, rmax, rsq, a )
        f(i) = f(i) + M_zI*ip
      end do

    case(1)
      stop 'Believe it or not, cannot do 1D interpolation, only 2D or 3D.'
    end select
    deallocate(aux_u)
#endif
  
    deallocate(rsq, a, lnext, lcell)

  else

    if(mesh_in%parallel_in_domains) then
      ALLOCATE(f_global(mesh_in%np_global), mesh_in%np_global)
#if defined(HAVE_MPI)
      call X(vec_allgather)(mesh_in%vp, f_global, u)
#endif
    end if

    f = M_ZERO
    do i = 1, mesh_out%np
      k = i
      if(mesh_out%parallel_in_domains) &
        k = mesh_out%vp%local(mesh_out%vp%xlocal(mesh_out%vp%partno)+i-1)
      ix = mesh_out%lxyz(k, 1); if ( ix < mesh_in%nr(1, 1) .or. ix > mesh_in%nr(2, 1) ) cycle
      iy = mesh_out%lxyz(k, 2); if ( iy < mesh_in%nr(1, 2) .or. iy > mesh_in%nr(2, 2) ) cycle
      iz = mesh_out%lxyz(k, 3); if ( iz < mesh_in%nr(1, 3) .or. iz > mesh_in%nr(2, 3) ) cycle
      j = mesh_in%lxyz_inv(ix, iy, iz)

      if(mesh_in%parallel_in_domains) then
        if(j <= mesh_in%np_global) f(i) = f_global(j)
      else
        f(i) = u(j)
      end if
    end do

    if(mesh_in%parallel_in_domains) then
      deallocate(f_global)
    end if

  end if

  call pop_sub()
end subroutine X(mf_interpolate)
