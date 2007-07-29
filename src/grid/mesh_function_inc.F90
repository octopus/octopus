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
    !$omp parallel workshare
    d = sum(f(1:mesh%np)*mesh%vol_pp(1:mesh%np))
    !$omp end parallel workshare
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
  R_TYPE,       intent(in)    :: u(:)    ! u(mesh_in%np_global)
  R_TYPE,       intent(out)   :: f(:)    ! f(mesh%np)

  FLOAT :: p(MAX_DIM)
  integer :: ix, iy, iz
  R_TYPE, allocatable :: f_global(:)
  integer :: i, j, k
  type(qshep_t) :: interp

  call push_sub('mf_inc.Xmf_interpolate')

  if(full_interpolation) then

    select case(mesh_in%sb%dim)
    case(2)
      call init_qshep(interp, mesh_in%np_global, u, mesh_in%x(:, 1), mesh_in%x(:, 2))
      do i = 1, mesh_out%np
        p(1) = mesh_out%x(i, 1)
        p(2) = mesh_out%x(i, 2)
        f(i) = qshep_interpolate(interp, u, p(1:2))
      end do
      call kill_qshep(interp)
    case(3)
      call init_qshep(interp, mesh_in%np_global, u, &
                      mesh_in%x(:, 1), mesh_in%x(:, 2), mesh_in%x(:, 3))
      do i = 1, mesh_out%np
        p(1) = mesh_out%x(i, 1)
        p(2) = mesh_out%x(i, 2)
        p(3) = mesh_out%x(i, 3)
        f(i) = qshep_interpolate(interp, u, p)
      end do
      call kill_qshep(interp)
    case(1)
      stop 'Believe it or not, cannot do 1D interpolation, only 2D or 3D.'
    end select
  
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
        if(j > 0 .and. j <= mesh_in%np_global) f(i) = f_global(j)
      else
        if(j > 0 .and. j <= mesh_in%np_global) f(i) = u(j)
      end if
    end do

    if(mesh_in%parallel_in_domains) then
      deallocate(f_global)
    end if

  end if

  call pop_sub()
end subroutine X(mf_interpolate)


! ---------------------------------------------------------
! Given a function f defined on mesh, and a plane, it gives 
! back the values of f on the plane, by doing the suitable
! interpolation
subroutine X(mf_interpolate_on_plane)(mesh, plane, f, f_in_plane)
  type(mesh_t),       intent(in)  :: mesh
  type(mesh_plane_t), intent(in)  :: plane
  R_TYPE,             intent(in)  :: f(:)
  R_TYPE,             intent(out) :: f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv)

  integer :: i, j
  R_TYPE, allocatable :: f_global(:)
  FLOAT :: p(MAX_DIM)
  type(qshep_t) :: interp

  call push_sub('mf_inc.Xmf_interpolate_on_plane')

  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, mesh%x(1:mesh%np_global, 1), &
                  mesh%x(1:mesh%np_global, 2), mesh%x(1:mesh%np_global, 3) )

  do i = plane%nu, plane%mu
    do j = plane%nv, plane%mv
      p(1) = plane%origin(1) + i*plane%spacing * plane%u(1) + j * plane%spacing * plane%v(1)
      p(2) = plane%origin(2) + i*plane%spacing * plane%u(2) + j * plane%spacing * plane%v(2)
      p(3) = plane%origin(3) + i*plane%spacing * plane%u(3) + j * plane%spacing * plane%v(3)
      f_in_plane(i, j) = qshep_interpolate(interp, f_global, p)
    end do
  end do

  call kill_qshep(interp)

  deallocate(f_global)
  call pop_sub()
end subroutine X(mf_interpolate_on_plane)


! ---------------------------------------------------------
! Given a function f defined on mesh, and a plane, it gives 
! back the values of f on the plane, by doing the suitable
! interpolation
subroutine X(mf_interpolate_on_line)(mesh, line, f, f_in_line)
  type(mesh_t),       intent(in)  :: mesh
  type(mesh_line_t), intent(in)  :: line
  R_TYPE,             intent(in)  :: f(:)
  R_TYPE,             intent(out) :: f_in_line(line%nu:line%mu)

  integer :: i
  R_TYPE, allocatable :: f_global(:)
  FLOAT :: p(MAX_DIM)
  type(qshep_t) :: interp
  
  call push_sub('mf_inc.Xmf_interpolate_on_line')

  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, mesh%x(1:mesh%np_global, 1), mesh%x(1:mesh%np_global, 2))
  do i = line%nu, line%mu
    p(1) = line%origin(1) + i*line%spacing * line%u(1)
    p(2) = line%origin(2) + i*line%spacing * line%u(2)
    f_in_line(i) = qshep_interpolate(interp, f_global, p(1:2))
  end do
  call kill_qshep(interp)

  deallocate(f_global)

  call pop_sub()
end subroutine X(mf_interpolate_on_line)


! ---------------------------------------------------------
! This subroutine calculates the surface integral of a scalar
! function on a given plane
R_TYPE function X(mf_surface_integral_scalar) (mesh, f, plane) result(d)
  type(mesh_t), intent(in)       :: mesh
  R_TYPE,       intent(in)       :: f(:)  ! f(mesh%np)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: f_in_plane(:, :)

  call push_sub('mf_inc.mf_surface_integral')

  if(mesh%sb%dim .ne. 3) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality'
    call write_fatal(1)
  end if

  ALLOCATE(f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv), (plane%mu-plane%nu+1)*(plane%mv-plane%nv+1))

  call X(mf_interpolate_on_plane)(mesh, plane, f, f_in_plane)

  d = sum(f_in_plane(:, :)*plane%spacing**2)

  deallocate(f_in_plane)
  call pop_sub()
end function X(mf_surface_integral_scalar)


! ---------------------------------------------------------
! This subroutine calculates the surface integral of a vector
! function on a given plane
R_TYPE function X(mf_surface_integral_vector) (mesh, f, plane) result(d)
  type(mesh_t), intent(in)       :: mesh
  R_TYPE,       intent(in)       :: f(:, :)  ! f(mesh%np, MAX_DIM)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: fn(:)
  integer :: i

  call push_sub('mf_inc.mf_surface_integral')

  ALLOCATE(fn(mesh%np), mesh%np)
  do i = 1, mesh%np
    fn(i) = sum(f(i, :)*plane%n(:))
  end do

  d =  X(mf_surface_integral_scalar)(mesh, fn, plane)

  deallocate(fn)
  call pop_sub()
end function X(mf_surface_integral_vector)


! ---------------------------------------------------------
! This subroutine calculates the line integral of a scalar
! function on a given line
R_TYPE function X(mf_line_integral_scalar) (mesh, f, line) result(d)
  type(mesh_t), intent(in)       :: mesh
  R_TYPE,       intent(in)       :: f(:)  ! f(mesh%np)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: f_in_line(:)

  call push_sub('mf_inc.mf_line_integral_scalar')

  if(mesh%sb%dim .ne. 2) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality'
    call write_fatal(1)
  end if

  ALLOCATE(f_in_line(line%nu:line%mu), line%mu-line%nu+1)

  call X(mf_interpolate_on_line)(mesh, line, f, f_in_line)

  d = sum(f_in_line(:)*line%spacing)

  deallocate(f_in_line)
  call pop_sub()
end function X(mf_line_integral_scalar)


! ---------------------------------------------------------
! This subroutine calculates the line integral of a vector
! function on a given line
R_TYPE function X(mf_line_integral_vector) (mesh, f, line) result(d)
  type(mesh_t), intent(in)       :: mesh
  R_TYPE,       intent(in)       :: f(:, :)  ! f(mesh%np, MAX_DIM)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: fn(:)
  integer :: i

  call push_sub('mf_inc.mf_line_integral_vector')

  ALLOCATE(fn(mesh%np), mesh%np)
  do i = 1, mesh%np
    fn(i) = sum(f(i, :)*line%n(:))
  end do

  d =  X(mf_line_integral_scalar)(mesh, fn, line)

  deallocate(fn)
  call pop_sub()
end function X(mf_line_integral_vector)


! puts a spline that represents a radial function into a mesh
subroutine X(mf_put_radial_spline)(m, spl, center, f, add)
  type(mesh_t),        intent(in)    :: m
  type(loct_spline_t), intent(in)    :: spl
  FLOAT,               intent(in)    :: center(1:MAX_DIM)
  R_TYPE,              intent(inout) :: f(:)
  logical, optional,   intent(in)    :: add

  integer :: ip
  FLOAT :: r

  logical :: add_

  add_ = .false.

  if(present(add)) add_ = add
  
  if( add_ ) then 

    !$omp parallel do private(r)
    do ip = 1, m%np
      r = sqrt(sum((m%x(ip, 1:MAX_DIM) - center(1:MAX_DIM))**2))
      f(ip) = f(ip) + loct_splint(spl, r)
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(r)
     do ip = 1, m%np
      r = sqrt(sum((m%x(ip, 1:MAX_DIM) - center(1:MAX_DIM))**2))
      f(ip) = loct_splint(spl, r)
    end do
    !$omp end parallel do

  end if

end subroutine X(mf_put_radial_spline)





!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
