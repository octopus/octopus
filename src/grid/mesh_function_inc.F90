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
! This function returns the dot product between two vectors,
! but using the mesh_aux defined as a global object in this
! module. This way it can be called by external libraries,
! passing only the two vectors. First, one has to 
! make sure that mesh_aux is pointint to some defined
! mesh data structure, by calling mesh_init_mesh_aux.
! ---------------------------------------------------------
R_TYPE function X(mf_dotp_aux)(f1, f2) result(dotp)
  R_TYPE,            intent(in) :: f1(:), f2(:)
  dotp = X(mf_dotp)(mesh_aux, f1, f2)
end function X(mf_dotp_aux)

! ---------------------------------------------------------
! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp_1)(mesh, f1, f2, reduce) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: f1(:), f2(:)
  logical, optional, intent(in) :: reduce

  R_TYPE, allocatable :: l(:)
  R_TYPE              :: dotp_tmp
  logical             :: reduce_

  call profiling_in(C_PROFILING_MF_DOTP)
  call push_sub('mf_inc.Xmf_dotp')

  reduce_ = .true.
  if(present(reduce)) then
    reduce_ = reduce
  end if

  ! This is not implemented via vec_integrate
  ! because BLAS is much faster.
  if(mesh%use_curvlinear) then
    ALLOCATE(l(mesh%np), mesh%np)
    l(1:mesh%np) = f1(1:mesh%np) * mesh%vol_pp(1:mesh%np)
    dotp_tmp  = lalg_dot(mesh%np, l(:),  f2(:))

    call profiling_count_operations(mesh%np*(2*R_ADD + R_MUL))

    deallocate(l)
  else
    dotp_tmp = lalg_dot(mesh%np, f1(:),  f2(:))*mesh%vol_pp(1)
    
    call profiling_count_operations(mesh%np*(R_ADD + R_MUL))

  end if

  if(mesh%parallel_in_domains.and.reduce_) then
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

end function X(mf_dotp_1)


! ---------------------------------------------------------
R_TYPE function X(mf_dotp_2)(mesh, dim, f1, f2, reduce) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim
#ifdef HAVE_MPI
  logical :: reduce_
  R_TYPE :: dotp_tmp
#endif

  call push_sub('states_inc.Xmf_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp_1)(mesh, f1(:, idim), f2(:, idim), reduce = .false.)
  end do

#ifdef HAVE_MPI
  reduce_ = .true.
  if(present(reduce)) reduce_ = reduce

  if(mesh%parallel_in_domains.and.reduce_) then
    call profiling_in(C_PROFILING_MF_DOTP_ALLREDUCE)
    dotp_tmp = dotp
    call MPI_Allreduce(dotp_tmp, dotp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_DOTP_ALLREDUCE)
  end if
#endif

  call pop_sub()

end function X(mf_dotp_2)

! ---------------------------------------------------------
! this function returns the the norm of a vector
FLOAT function X(mf_nrm2)(mesh, f, reduce) result(nrm2)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: f(:)
  logical, optional, intent(in) :: reduce

  R_TYPE, allocatable :: l(:)
  FLOAT               :: nrm2_tmp
  logical             :: reduce_

  call profiling_in(C_PROFILING_MF_NRM2)
  call push_sub('mf_inc.Xmf_nrm2')

  ! This is not implemented via vec_integrate
  ! because BLAS is much faster.
  if(mesh%use_curvlinear) then
    ALLOCATE(l(mesh%np), mesh%np)
    l(1:mesh%np) = f(1:mesh%np) * sqrt(mesh%vol_pp(1:mesh%np))
    nrm2_tmp = lalg_nrm2(mesh%np, l)
    deallocate(l)
  else
    nrm2_tmp = lalg_nrm2(mesh%np, f) * sqrt(mesh%vol_pp(1))
  end if

  reduce_ = .true.
#if defined(HAVE_MPI)
  if(present(reduce)) reduce_ = reduce
#endif

  if(mesh%parallel_in_domains .and. reduce_) then
#if defined(HAVE_MPI)
    nrm2_tmp = nrm2_tmp**2
    call profiling_in(C_PROFILING_MF_DOTP_ALLREDUCE)
    call MPI_Allreduce(nrm2_tmp, nrm2, 1, MPI_FLOAT, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_DOTP_ALLREDUCE)
    nrm2 = sqrt(nrm2)
#else
    ASSERT(.false.)
#endif
  else
    nrm2 = nrm2_tmp
  end if

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
subroutine X(mf_random)(m, f, seed)
  type(mesh_t), intent(in)  :: m
  R_TYPE,       intent(out) :: f(1:m%np)
  integer, optional, intent(in) :: seed

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(MAX_DIM), rnd, r

  call push_sub('mf_inc.Xmf_random')

  if(present(seed)) then
    iseed = iseed + seed
  end if

  a = M_ZERO
  do i = 1, m%sb%dim
    call quickrnd(iseed, rnd)
    a(i) = M_TWO*(2*rnd - 1)
  end do

  !$omp parallel do private(r)
  do i = 1, m%np
    r = sum((m%x(i, 1:MAX_DIM) - a(1:MAX_DIM))**2)
    if ( r < CNST(100.0) ) then 
      f(i) = exp(-M_HALF*r)
    else
      f(i) = M_ZERO
    end if
  end do
  !$omp end parallel do

  r = X(mf_nrm2)(m, f)
  call lalg_scal(m%np, M_ONE/r, f)

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
  type(mesh_t),         intent(in)  :: mesh_in, mesh_out
  logical,              intent(in)  :: full_interpolation
  R_TYPE,               intent(in)  :: u(:)    ! u(mesh_in%np_global)
  R_TYPE,               intent(out) :: f(:)    ! f(mesh%np)

  integer :: ix, iy, iz
  R_TYPE, allocatable :: f_global(:)
  integer :: i, j, k

  call push_sub('mf_inc.Xmf_interpolate')

  if(full_interpolation) then
    call X(mf_interpolate_points) (mesh_in, u, mesh_out%np, mesh_out%x, f)
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
! this function receives a function u defined in a mesh, and returns
! in the interpolated values of the function over the npoints defined
! by x

subroutine X(mf_interpolate_points) (mesh_in, u, npoints, x, f)
  type(mesh_t),         intent(in)  :: mesh_in
  R_TYPE, target,       intent(in)  :: u(:)    ! u(mesh_in%np_global)
  integer,              intent(in)  :: npoints
  FLOAT,                intent(in)  :: x(:, :)
  R_TYPE,               intent(out) :: f(:)    ! f(mesh%np)

  real(8) :: p(MAX_DIM)
  R_DOUBLE, pointer :: ru(:)
  integer :: i
  type(qshep_t) :: interp
  real(8), pointer :: rx(:, :)
#ifndef R_TCOMPLEX
  type(spline_t) :: interp1d
#endif

  call push_sub('mf_inc.Xmf_interpolate')

#ifdef SINGLE_PRECISION
  ALLOCATE(rx(1:ubound(mesh_in%x, DIM=1), 1:MAX_DIM), ubound(mesh_in%x, DIM=1)*MAX_DIM)
  rx = mesh_in%x
  ALLOCATE(ru(1:ubound(u, DIM=1)), ubound(u, DIM=1))
  ru = u
#else
  rx => mesh_in%x
  ru => u
#endif

  select case(mesh_in%sb%dim)
  case(2)
    call init_qshep(interp, mesh_in%np_global, ru, rx(:, 1), rx(:, 2))
    do i = 1, npoints
      p(1) = x(i, 1)
      p(2) = x(i, 2)
      f(i) = qshep_interpolate(interp, ru, p(1:2))
    end do
    call kill_qshep(interp)
  case(3)
    call init_qshep(interp, mesh_in%np_global, ru, rx(:, 1), rx(:, 2), rx(:, 3))
    do i = 1, npoints
      p(1) = x(i, 1)
      p(2) = x(i, 2)
      p(3) = x(i, 3)
      f(i) = qshep_interpolate(interp, ru, p)
    end do
    call kill_qshep(interp)
  case(1)
#ifdef R_TCOMPLEX
    message(1) = 'Believe it or not, cannot do 1D complex interpolation, only 2D or 3D.'
    call write_fatal(1)
#else
    call spline_init(interp1d)
    call spline_fit(mesh_in%np_global, R_REAL(rx(:, 1)), ru, interp1d)
    do i = 1, npoints
      f(i) = spline_eval(interp1d, x(i, 1))
    end do
    call spline_end(interp1d)
#endif
  end select
#ifdef SINGLE_PRECISION
  deallocate(rx)
  deallocate(ru)
#endif

  call pop_sub()
end subroutine X(mf_interpolate_points)


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
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: p(MAX_DIM)
  type(qshep_t) :: interp
  real(8), pointer :: rx(:, :)

  call push_sub('mf_inc.Xmf_interpolate_on_plane')

#ifdef SINGLE_PRECISION
    ALLOCATE(rx(1:ubound(mesh%x, DIM=1), 1:MAX_DIM), ubound(mesh%x, DIM=1)*MAX_DIM)
    rx = mesh%x
#else
    rx => mesh%x
#endif

  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, rx(:, 1), rx(:, 2), rx(:, 3) )

  do i = plane%nu, plane%mu
    do j = plane%nv, plane%mv
      p(1) = plane%origin(1) + i*plane%spacing * plane%u(1) + j * plane%spacing * plane%v(1)
      p(2) = plane%origin(2) + i*plane%spacing * plane%u(2) + j * plane%spacing * plane%v(2)
      p(3) = plane%origin(3) + i*plane%spacing * plane%u(3) + j * plane%spacing * plane%v(3)
      f_in_plane(i, j) = qshep_interpolate(interp, f_global, p)
    end do
  end do

  call kill_qshep(interp)
#ifdef SINGLE_PRECISION
    deallocate(rx)
#endif
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
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: p(MAX_DIM)
  type(qshep_t) :: interp
  real(8), pointer :: rx(:, :)

  call push_sub('mf_inc.Xmf_interpolate_on_line')

#ifdef SINGLE_PRECISION
    ALLOCATE(rx(1:ubound(mesh%x, DIM=1), 1:MAX_DIM), ubound(mesh%x, DIM=1)*MAX_DIM)
    rx = mesh%x
#else
    rx => mesh%x
#endif

  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, rx(:, 1), rx(:, 2))
  do i = line%nu, line%mu
    p(1) = line%origin(1) + i*line%spacing * line%u(1)
    p(2) = line%origin(2) + i*line%spacing * line%u(2)
    f_in_line(i) = qshep_interpolate(interp, f_global, p(1:2))
  end do
  call kill_qshep(interp)

  deallocate(f_global)
#ifdef SINGLE_PRECISION
    deallocate(rx)
#endif
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
  type(spline_t), intent(in)    :: spl
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
      f(ip) = f(ip) + spline_eval(spl, r)
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(r)
     do ip = 1, m%np
      r = sqrt(sum((m%x(ip, 1:MAX_DIM) - center(1:MAX_DIM))**2))
      f(ip) = spline_eval(spl, r)
    end do
    !$omp end parallel do

  end if

end subroutine X(mf_put_radial_spline)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
