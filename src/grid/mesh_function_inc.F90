!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Verstraete
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

  integer :: ip
#ifdef HAVE_MPI
  R_TYPE :: d_local
#endif

  call profiling_in(C_PROFILING_MF_INTEGRATE, 'MF_INTEGRATE')
  call push_sub('mesh_function_inc.Xmf_integrate')

  d = M_ZERO
  if (mesh%use_curvilinear) then
    do ip = 1, mesh%np
      d = d + f(ip)*mesh%vol_pp(ip)
    end do
  else
    do ip = 1, mesh%np
      d = d + f(ip)
    end do
    d = d*mesh%vol_pp(1)
  end if

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    d_local = d
    call MPI_Allreduce(d_local, d, 1, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if
#endif
  
  call pop_sub()
  call profiling_out(C_PROFILING_MF_INTEGRATE)

end function X(mf_integrate)


! ---------------------------------------------------------
! This function returns the dot product between two vectors,
! but using the mesh_aux defined as a global object in this
! module. This way it can be called by external libraries,
! passing only the two vectors. First, one has to 
! make sure that mesh_aux is pointing to some defined
! mesh data structure, by calling mesh_init_mesh_aux.
! ---------------------------------------------------------
R_TYPE function X(mf_dotp_aux)(f1, f2) result(dotp)
  R_TYPE,            intent(in) :: f1(:), f2(:)
  dotp = X(mf_dotp)(mesh_aux, f1, f2)
end function X(mf_dotp_aux)

! ---------------------------------------------------------
! this function returns the dot product between two vectors
R_TYPE function X(mf_dotp_1)(mesh, f1, f2, reduce, dotu) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: f1(:), f2(:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
     ! if true, use lalg_dotu instead of lalg_dot;
     ! no complex conjugation.  Default is false.
     ! has no effect if working with real version

  R_TYPE              :: dotp_tmp
  logical             :: reduce_
#ifdef R_TCOMPLEX
  logical             :: dotu_
#endif
  integer             :: ip

  call profiling_in(C_PROFILING_MF_DOTP, "MF_DOTP")
  call push_sub('mesh_function_inc.Xmf_dotp_1')

  reduce_ = .true.
  if(present(reduce)) then
    reduce_ = reduce
  end if

#ifdef R_TCOMPLEX
  dotu_ = .false.
  if(present(dotu)) then
    dotu_ = dotu
  end if
#endif

  if(mesh%use_curvilinear) then
    dotp_tmp = M_ZERO
    ! preprocessor conditionals necessary since lalg_dotu only exists for complex input
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      do ip = 1, mesh%np
        dotp_tmp = dotp_tmp + mesh%vol_pp(ip)*R_CONJ(f1(ip))*f2(ip)
      end do
#ifdef R_TCOMPLEX
    else
      do ip = 1, mesh%np
        dotp_tmp = dotp_tmp + mesh%vol_pp(ip)*f1(ip)*f2(ip)
      end do
    endif
#endif
    call profiling_count_operations(mesh%np*(2*R_ADD + R_MUL))
  else
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      dotp_tmp = lalg_dot(mesh%np, f1(:),  f2(:))*mesh%vol_pp(1)
#ifdef R_TCOMPLEX
    else
      dotp_tmp = lalg_dotu(mesh%np, f1(:),  f2(:))*mesh%vol_pp(1)
    endif
#endif
    call profiling_count_operations(mesh%np * (R_ADD + R_MUL))

  end if

  if(mesh%parallel_in_domains.and.reduce_) then
#if defined(HAVE_MPI)
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    call MPI_Allreduce(dotp_tmp, dotp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_REDUCE)
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
R_TYPE function X(mf_dotp_2)(mesh, dim, f1, f2, reduce, dotu) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
     ! if true, use lalg_dotu instead of lalg_dot;
     ! no complex conjugation.  Default is false.

  integer :: idim
  logical :: dotu_
#ifdef HAVE_MPI
  logical :: reduce_
  R_TYPE :: dotp_tmp
#endif

  call push_sub('mesh_function_inc.Xmf_dotp_2')

  dotu_ = .false.
  if(present(dotu)) then
    dotu_ = dotu
  end if

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp_1)(mesh, f1(:, idim), f2(:, idim), reduce = .false., dotu = dotu_)
  end do

#ifdef HAVE_MPI
  reduce_ = .true.
  if(present(reduce)) then
    reduce_ = reduce
  endif

  if(mesh%parallel_in_domains.and.reduce_) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    dotp_tmp = dotp
    call MPI_Allreduce(dotp_tmp, dotp, 1, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if
#endif

  call pop_sub()

end function X(mf_dotp_2)

! ---------------------------------------------------------
! this function returns the the norm of a vector
FLOAT function X(mf_nrm2_1)(mesh, f, reduce) result(nrm2)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: f(:)
  logical, optional, intent(in) :: reduce

  FLOAT               :: nrm2_tmp
  logical             :: reduce_
  R_TYPE, allocatable :: l(:)

  call profiling_in(C_PROFILING_MF_NRM2, "MF_NRM2")
  call push_sub('mesh_function_inc.Xmf_nrm2_1')

  if(mesh%use_curvilinear) then
    ALLOCATE(l(mesh%np), mesh%np)
    l(1:mesh%np) = f(1:mesh%np)*sqrt(mesh%vol_pp(1:mesh%np))
    nrm2_tmp = lalg_nrm2(mesh%np, l)
    SAFE_DEALLOCATE_A(l)
  else
    nrm2_tmp = lalg_nrm2(mesh%np, f)*sqrt(mesh%vol_pp(1))
  end if

  reduce_ = .true.
#if defined(HAVE_MPI)
  if(present(reduce)) then
    reduce_ = reduce
  endif
#endif

  if(mesh%parallel_in_domains .and. reduce_) then
#if defined(HAVE_MPI)
    nrm2_tmp = nrm2_tmp**2
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    call MPI_Allreduce(nrm2_tmp, nrm2, 1, MPI_FLOAT, MPI_SUM, mesh%vp%comm, mpi_err)
    call profiling_out(C_PROFILING_MF_REDUCE)
    nrm2 = sqrt(nrm2)
#else
    ASSERT(.false.)
#endif
  else
    nrm2 = nrm2_tmp
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_MF_NRM2)

end function X(mf_nrm2_1)

! ---------------------------------------------------------
FLOAT function X(mf_nrm2_2)(mesh, dim, f, reduce) result(nrm2)
  type(mesh_t),      intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  call push_sub('mesh_function_inc.Xmf_nrm2_2')

  nrm2 = M_ZERO

  do idim = 1, dim
    if(present(reduce)) then
      nrm2 = hypot(nrm2, X(mf_nrm2)(mesh, f(:, idim), reduce = reduce))
    else
      nrm2 = hypot(nrm2, X(mf_nrm2)(mesh, f(:, idim)))
    end if
  end do

  call pop_sub()

end function X(mf_nrm2_2)

! ---------------------------------------------------------
! This function calculates the x_i moment of the function f
R_TYPE function X(mf_moment) (mesh, ff, idir, order) result(rr)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: ff(:)
  integer,      intent(in) :: idir
  integer,      intent(in) :: order

  R_TYPE, allocatable :: fxn(:)

  call push_sub('mesh_function_inc.Xmf_moment')

  ALLOCATE(fxn(1:mesh%np), mesh%np)

  fxn(1:mesh%np) = ff(1:mesh%np)*mesh%x(1:mesh%np, idir)**order
  rr = X(mf_integrate)(mesh, fxn)

  SAFE_DEALLOCATE_A(fxn)
  call pop_sub()

end function X(mf_moment)


! ---------------------------------------------------------
! This subroutine generates a Gaussian wave-function at a
! random position in space
subroutine X(mf_random)(mesh, f, seed)
  type(mesh_t),      intent(in)  :: mesh
  R_TYPE,            intent(out) :: f(:)
  integer, optional, intent(in)  :: seed

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(MAX_DIM), rnd, r

  call push_sub('mesh_function_inc.Xmf_random')

  if(present(seed)) then
    iseed = iseed + seed
  end if

  a = M_ZERO
  do i = 1, mesh%sb%dim
    call quickrnd(iseed, rnd)
    a(i) = M_TWO*(2*rnd - 1)
  end do

  !$omp parallel do private(r)
  do i = 1, mesh%np
    r = sum((mesh%x(i, 1:mesh%sb%dim) - a(1:mesh%sb%dim))**2)
    if ( r < CNST(100.0) ) then 
      f(i) = exp(-M_HALF*r)
    else
      f(i) = M_ZERO
    end if
  end do
  !$omp end parallel do

  r = X(mf_nrm2)(mesh, f)
  call lalg_scal(mesh%np, M_ONE/r, f)

  call pop_sub()

end subroutine X(mf_random)

! ---------------------------------------------------------
! Does a partial integration of a function, i.e. in N
! dimensions, it integrates over N-1 dimensions. The
! argument j is the dimension which is not integrated.
! The input function f is a normal function defined on
! the mesh, whereas the output function is a 1D array
! of mesh%idx%ll(j) + 2*mesh%enlarge(j) elements.
!
! In order to retrieve the coordinate v of u(n), one needs
! to do:
!  k = gr%mesh%idx%nr(1, j) + n - 1
!  i = gr%mesh%idx%Lxyz_inv(k, 0, 0) ! Assuming j = 1
!  if(i>0) v = gr%mesh%x(i, j)   ! I do not understand why
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
  call push_sub('mesh_function_inc.Xmf_partial_integrate')

  ASSERT(.not.(mesh%parallel_in_domains))
  ASSERT(.not.(mesh%use_curvilinear))

  k = mesh%idx%ll(j) + 2*mesh%idx%enlarge(j)

  u(1:k) = M_ZERO
  do i = 1, mesh%np
     m = mesh%idx%Lxyz(i, j) - mesh%idx%nr(1, j) + 1
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

  call push_sub('mesh_function_inc.Xmf_interpolate')

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
      ix = mesh_out%idx%Lxyz(k, 1); if ( ix < mesh_in%idx%nr(1, 1) .or. ix > mesh_in%idx%nr(2, 1) ) cycle
      iy = mesh_out%idx%Lxyz(k, 2); if ( iy < mesh_in%idx%nr(1, 2) .or. iy > mesh_in%idx%nr(2, 2) ) cycle
      iz = mesh_out%idx%Lxyz(k, 3); if ( iz < mesh_in%idx%nr(1, 3) .or. iz > mesh_in%idx%nr(2, 3) ) cycle
      j = mesh_in%idx%Lxyz_inv(ix, iy, iz)

      if(mesh_in%parallel_in_domains) then
        if(j > 0 .and. j <= mesh_in%np_global) f(i) = f_global(j)
      else
        if(j > 0 .and. j <= mesh_in%np_global) f(i) = u(j)
      end if
    end do

    if(mesh_in%parallel_in_domains) then
      SAFE_DEALLOCATE_A(f_global)
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

  call push_sub('mesh_function_inc.Xmf_interpolate_points')

#ifdef SINGLE_PRECISION
  ALLOCATE(rx(1:ubound(mesh_in%x, DIM=1), 1:mesh_in%sb%dim), ubound(mesh_in%x, DIM=1)*mesh_in%sb%dim)
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
  SAFE_DEALLOCATE_P(rx)
  SAFE_DEALLOCATE_P(ru)
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

  integer :: i, j, ip
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: p(MAX_DIM)
  type(qshep_t) :: interp
  real(8), allocatable :: xglobal(:, :)

  call push_sub('mesh_function_inc.Xmf_interpolate_on_plane')

  ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM), mesh%np_part_global*MAX_DIM)
  do ip = 1, mesh%np_part_global
    xglobal(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip)
  end do

  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2), xglobal(:, 3) )

  do i = plane%nu, plane%mu
    do j = plane%nv, plane%mv
      p(1) = plane%origin(1) + i*plane%spacing * plane%u(1) + j * plane%spacing * plane%v(1)
      p(2) = plane%origin(2) + i*plane%spacing * plane%u(2) + j * plane%spacing * plane%v(2)
      p(3) = plane%origin(3) + i*plane%spacing * plane%u(3) + j * plane%spacing * plane%v(3)
      f_in_plane(i, j) = qshep_interpolate(interp, f_global, p)
    end do
  end do

  call kill_qshep(interp)

  SAFE_DEALLOCATE_A(xglobal)
  SAFE_DEALLOCATE_A(f_global)
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

  integer :: i, ip
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: p(MAX_DIM)
  type(qshep_t) :: interp
  real(8), allocatable :: xglobal(:, :)

  call push_sub('mesh_function_inc.Xmf_interpolate_on_line')

  ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM), mesh%np_part_global*MAX_DIM)
  do ip = 1, mesh%np_part_global
     xglobal(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip)
  end do
  
  ALLOCATE(f_global(mesh%np_global), mesh%np_global)
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, f_global, f)
#else
  f_global(1:mesh%np_global) = f(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2))
  do i = line%nu, line%mu
    p(1) = line%origin(1) + i*line%spacing * line%u(1)
    p(2) = line%origin(2) + i*line%spacing * line%u(2)
    f_in_line(i) = qshep_interpolate(interp, f_global, p(1:2))
  end do
  call kill_qshep(interp)

  SAFE_DEALLOCATE_A(f_global)
  SAFE_DEALLOCATE_A(xglobal)

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

  call push_sub('mesh_function_inc.Xmf_surface_integral_scalar')

  if(mesh%sb%dim .ne. 3) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality'
    call write_fatal(1)
  end if

  ALLOCATE(f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv), (plane%mu-plane%nu+1)*(plane%mv-plane%nv+1))

  call X(mf_interpolate_on_plane)(mesh, plane, f, f_in_plane)

  d = sum(f_in_plane(:, :)*plane%spacing**2)

  SAFE_DEALLOCATE_A(f_in_plane)
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

  call push_sub('mesh_function_inc.Xmf_surface_integral_vector')

  ALLOCATE(fn(mesh%np), mesh%np)
  do i = 1, mesh%np
    fn(i) = sum(f(i, :)*plane%n(:))
  end do

  d =  X(mf_surface_integral_scalar)(mesh, fn, plane)

  SAFE_DEALLOCATE_A(fn)
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

  call push_sub('mesh_function_inc.Xmf_line_integral_scalar')

  if(mesh%sb%dim .ne. 2) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality'
    call write_fatal(1)
  end if

  ALLOCATE(f_in_line(line%nu:line%mu), line%mu-line%nu+1)

  call X(mf_interpolate_on_line)(mesh, line, f, f_in_line)

  d = sum(f_in_line(:)*line%spacing)

  SAFE_DEALLOCATE_A(f_in_line)
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

  call push_sub('mesh_function_inc.Xmf_line_integral_vector')

  ALLOCATE(fn(mesh%np), mesh%np)
  do i = 1, mesh%np
    fn(i) = sum(f(i, :)*line%n(:))
  end do

  d =  X(mf_line_integral_scalar)(mesh, fn, line)

  SAFE_DEALLOCATE_A(fn)
  call pop_sub()
end function X(mf_line_integral_vector)


! puts a spline that represents a radial function into a mesh
subroutine X(mf_put_radial_spline)(mesh, spl, center, f, add)
  type(mesh_t),        intent(in)    :: mesh
  type(spline_t),      intent(in)    :: spl
  FLOAT,               intent(in)    :: center(1:MAX_DIM)
  R_TYPE,              intent(inout) :: f(:)
  logical, optional,   intent(in)    :: add

  integer :: ip
  FLOAT :: r

  logical :: add_

  add_ = .false.

  if(present(add)) then
    add_ = add
  endif

  if( add_ ) then 

    do ip = 1, mesh%np
      r = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
      f(ip) = f(ip) + spline_eval(spl, r)
    end do

  else

     do ip = 1, mesh%np
      r = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
      f(ip) = spline_eval(spl, r)
    end do

  end if

end subroutine X(mf_put_radial_spline)


! -----------------------------------------------------------------------------
! This routine calculates the multipoles of a function f
! distribution, defined in the following way:
! multipole(1) is the trace of f (defined to be positive; integral
!   of the f.
! multipole(2:4) contains the dipole: integral of f times x, y or z.
! multipole(5:9, is) contains the quadrupole, defined in the usual way using
!   the spherical harmonics: multipole(5) = Integral [ f * Y_{2,-2} ],
!   multipole(6, is) = Integral [ f * Y_{2, -1} ].
! And so on.
! -----------------------------------------------------------------------------
subroutine X(mf_multipoles) (mesh, ff, lmax, multipole)
  type(mesh_t),   intent(in)  :: mesh
  R_TYPE,         intent(in)  :: ff(:)
  integer,        intent(in)  :: lmax
  R_TYPE,         intent(out) :: multipole(:) ! multipole((lmax + 1)**2)

  integer :: i, l, lm, add_lm
  FLOAT   :: x(MAX_DIM), r, ylm
  R_TYPE, allocatable :: ff2(:)

  call push_sub('mesh_function_inc.Xmf_multipoles')

  ALLOCATE(ff2(mesh%np), mesh%np)

  ff2(1:mesh%np) = ff(1:mesh%np)
  multipole(1) = X(mf_integrate)(mesh, ff2)

  if(lmax > 0) then
    do i = 1, 3
      ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, i)
      multipole(i+1) = X(mf_integrate)(mesh, ff2)
    end do
  end if

  if(lmax>1) then
    add_lm = 5
    do l = 2, lmax
      do lm = -l, l
        do i = 1, mesh%np
          call mesh_r(mesh, i, r, x=x)
          call loct_ylm(1, x(1), x(2), x(3), l, lm, ylm)
          ff2(i) = ff(i) * ylm * r**l
        end do
        multipole(add_lm) = X(mf_integrate)(mesh, ff2)
        add_lm = add_lm + 1
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(ff2)
  call pop_sub()
end subroutine X(mf_multipoles)

subroutine X(mf_dotp_batch)(mesh, aa, bb, dot, symm)
  type(mesh_t),      intent(in)    :: mesh
  type(batch_t),     intent(in)    :: aa
  type(batch_t),     intent(in)    :: bb
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: symm         !for the moment it is ignored

  integer :: ist, jst, idim, sp, block_size, ep, ip, lda, ldb
  R_TYPE :: ss
  R_TYPE, allocatable :: dd(:, :)
#ifdef HAVE_MPI
  R_TYPE, allocatable :: ddtmp(:, :)
#endif
  type(profile_t), save :: prof
  logical :: use_blas

  call push_sub('mesh_function_inc.Xmf_dotp_batch')
  call profiling_in(prof, "DOTP_BATCH")

  ASSERT(aa%dim == bb%dim)

  ALLOCATE(dd(1:aa%nst, 1:bb%nst), aa%nst*bb%nst)

  use_blas = associated(aa%X(psicont)) .and. associated(bb%X(psicont)) .and. (.not. mesh%use_curvilinear) .and. (aa%dim == 1)

  if(use_blas) then

    lda = size(aa%X(psicont), dim = 1)
    ldb = size(bb%X(psicont), dim = 1)
    call blas_gemm('c', 'n', aa%nst, bb%nst, mesh%np, &
         R_TOTYPE(mesh%vol_pp(1)), &
         aa%X(psicont)(1, 1, 1), lda, &
         bb%X(psicont)(1, 1, 1), ldb, &
         R_TOTYPE(M_ZERO), dd(1, 1), aa%nst)

  else

    dd = R_TOTYPE(M_ZERO)

    block_size = hardware%X(block_size)

    do idim = 1, aa%dim
      do sp = 1, mesh%np, block_size
        ep = min(mesh%np, sp + block_size - 1)

        if(mesh%use_curvilinear) then

          do ist = 1, aa%nst
            do jst = 1, bb%nst

              ss = M_ZERO
              do ip = sp, ep
                ss = ss + mesh%vol_pp(ip)*R_CONJ(aa%states(ist)%X(psi)(ip, idim))*bb%states(jst)%X(psi)(ip, idim)
              end do
              dd(ist, jst) = dd(ist, jst) + ss

            end do
          end do

        else

          do ist = 1, aa%nst
            do jst = 1, bb%nst
              dd(ist, jst) = dd(ist, jst) + &
                   mesh%vol_pp(1)*blas_dot(ep - sp + 1, aa%states(ist)%X(psi)(sp, idim), 1, bb%states(jst)%X(psi)(sp, idim), 1)
            end do
          end do

        end if
      end do
    end do

  end if

  if(mesh%use_curvilinear) then
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + 2*R_MUL))
  else
    call profiling_count_operations(dble(mesh%np)*aa%nst*bb%nst*(R_ADD + R_MUL))
  end if

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    ALLOCATE(ddtmp(1:aa%nst, 1:bb%nst), aa%nst*bb%nst)
    forall(ist = 1:aa%nst, jst = 1:bb%nst) ddtmp(ist, jst) = dd(ist, jst)
    call MPI_Allreduce(ddtmp, dd, aa%nst*bb%nst, R_MPITYPE, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(ddtmp)
  end if
#endif

  forall(ist = 1:aa%nst, jst = 1:bb%nst) dot(aa%states(ist)%ist, bb%states(jst)%ist) = dd(ist, jst)

  SAFE_DEALLOCATE_A(dd)

  call profiling_out(prof)
  call pop_sub()
end subroutine X(mf_dotp_batch)

!------------------------------------------------------------
! This routine calculates the one-body density matrix gamma 
! for particle ikeeppart, used in higher dimensional model
! hamiltonian calculations (MJV, NH) 
!------------------------------------------------------------
subroutine X(mf_calculate_gamma)(ikeeppart, nparticles_densmat, ndim1part, hypercube_1part, &
     npt_1part, nr_1part, enlarge_1part,   mesh, psi, gamma)
  integer, intent(in)      :: ikeeppart, npt_1part, ndim1part
  integer, intent(in)      :: nr_1part(ndim1part)
  integer, intent(in)      :: enlarge_1part
  integer, intent(in)      :: nparticles_densmat
  type(mesh_t), intent(in) :: mesh
  type(hypercube_t),intent(in) :: hypercube_1part
  R_TYPE, intent(in)       :: psi(:)
  CMPLX, intent(out)       :: gamma(:, :)

  integer :: imesh, imeshp, icoord, icoordp, jdim
  integer, allocatable :: ix(:), ix_1part(:), ixp(:)
  FLOAT :: volume_element
  R_TYPE, allocatable :: psi_global(:)

  call push_sub('mesh_function_inc.Xmf_calculate_gamma')

  write (*,*) ' enlarge_1part, npt_1part ', enlarge_1part, npt_1part

  ! In case of running parallel in domains, we need to operate psi_global, which 
  ! contains the full wave function after "gathering" all the domains.
  ALLOCATE(psi_global(mesh%np_part_global), mesh%np_part_global)
  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    call X(vec_allgather)(mesh%vp, psi_global, psi)
#endif
  else
    psi_global(1:mesh%np_part_global) = psi(1:mesh%np_part_global)
  end if

  ALLOCATE(ix(MAX_DIM), MAX_DIM)
  ALLOCATE(ixp(MAX_DIM), MAX_DIM)

  volume_element = 1.0d0
  do jdim = 1, MAX_DIM
    if (mesh%h(jdim) > 1.e-10) volume_element=volume_element*mesh%h(jdim)
  end do
  do jdim = (ikeeppart - 1)*ndim1part + 1, ikeeppart*ndim1part
    if (mesh%h(jdim) > 1.e-10) volume_element = volume_element/mesh%h(jdim)
  end do

  gamma = R_TOTYPE(M_ZERO)

  ALLOCATE(ix_1part(ndim1part), ndim1part)
  xloop: do imesh = 1, mesh%np_global
! find coordinates of present point in full MAX_DIM space
     call index_to_coords(mesh%idx, mesh%sb%dim, imesh, ix)

! find index of present coordinates for particle ikeeppart
     ix_1part = ix((ikeeppart - 1)*ndim1part + 1:ikeeppart*ndim1part)
     call hypercube_x_to_i(hypercube_1part, ndim1part, nr_1part, enlarge_1part, &
              ix_1part, icoord)
	      
! run over all possible position indices of particle prime
     ixp = ix
     xprimeloop: do icoordp = 1, npt_1part

! find equivalent position of particle prime
        call hypercube_i_to_x(hypercube_1part, ndim1part, nr_1part, enlarge_1part, icoordp, ix_1part)

! change coordinates of particle ikeeppart only 
        ixp((ikeeppart - 1)*ndim1part + 1:ikeeppart*ndim1part) = ix_1part

! find new index for general point prime
        imeshp = index_from_coords(mesh%idx, mesh%sb%dim, ixp)
        
! accumulate into density matrix
        gamma(icoord, icoordp) = gamma(icoord, icoordp) + &
            nparticles_densmat*volume_element*psi_global(imesh)*R_CONJ(psi_global (imeshp))
     end do xprimeloop
  end do xloop

  SAFE_DEALLOCATE_A(psi_global)
  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)

  call pop_sub()
end subroutine X(mf_calculate_gamma)

!------------------------------------------------------------
! This routine calculates the density rho for particle 
! ikeeppart, used in higher dimensional model
! hamiltonian calculations (MJV, NH) 
!------------------------------------------------------------
subroutine X(mf_calculate_rho)(ikeeppart, nparticles_dens, ndim1part, hypercube_1part, &
     npt_1part, nr_1part, enlarge_1part,   mesh, psi, rho)
  integer, intent(in)      :: ikeeppart, npt_1part, ndim1part
  integer, intent(in)      :: nr_1part(ndim1part)
  integer, intent(in)      :: enlarge_1part
  integer, intent(in)      :: nparticles_dens
  type(mesh_t), intent(in) :: mesh
  type(hypercube_t),intent(in) :: hypercube_1part
  R_TYPE, intent(in)       :: psi(:)
  FLOAT, intent(out)       :: rho(:)

  integer :: imesh, icoord, jdim
  integer, allocatable :: ix(:), ix_1part(:), ixp(:)
  FLOAT :: volume_element
  R_TYPE, allocatable :: psi_global(:)

  call push_sub('mesh_function_inc.Xmf_calculate_gamma')

  write (*,*) ' enlarge_1part, npt_1part ', enlarge_1part, npt_1part

  ! In case of running parallel in domains, we need to operate psi_global, which 
  ! contains the full wave function after "gathering" all the domains.
  ALLOCATE(psi_global(mesh%np_part_global), mesh%np_part_global)
  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    call X(vec_allgather)(mesh%vp, psi_global, psi)
#endif
  else
    psi_global(1:mesh%np_part_global) = psi(1:mesh%np_part_global)
  end if

  ALLOCATE(ix(MAX_DIM), MAX_DIM)
  ALLOCATE(ixp(MAX_DIM), MAX_DIM)

  volume_element = 1.0d0
  do jdim = 1, MAX_DIM
    if (mesh%h(jdim) > 1.e-10) volume_element=volume_element*mesh%h(jdim)
  end do
  do jdim = (ikeeppart - 1)*ndim1part + 1, ikeeppart*ndim1part
    if (mesh%h(jdim) > 1.e-10) volume_element = volume_element/mesh%h(jdim)
  end do

  rho = R_TOTYPE(M_ZERO)

  ALLOCATE(ix_1part(ndim1part), ndim1part)
  xloop: do imesh = 1, mesh%np_global
! find coordinates of present point in full MAX_DIM space
     call index_to_coords(mesh%idx, mesh%sb%dim, imesh, ix)

! find index of present coordinates for particle ikeeppart
     ix_1part = ix((ikeeppart - 1)*ndim1part + 1:ikeeppart*ndim1part)
     call hypercube_x_to_i(hypercube_1part, ndim1part, nr_1part, enlarge_1part, &
              ix_1part, icoord)
	      
! accumulate into density matrix
      rho(icoord) = rho(icoord) + &
            nparticles_dens*volume_element*psi_global(imesh)*R_CONJ(psi_global(imesh))
  end do xloop

  SAFE_DEALLOCATE_A(psi_global)
  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(ixp)
  SAFE_DEALLOCATE_A(ix_1part)

  call pop_sub()
end subroutine X(mf_calculate_rho)



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
