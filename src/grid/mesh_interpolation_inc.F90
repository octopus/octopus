!! Copyright (C) 2014 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!


subroutine X(mesh_interpolation_evaluate)(this, values, position, interpolated_value)
  type(mesh_interpolation_t), intent(in)    :: this
  R_TYPE,                     intent(in)    :: values(:)
  FLOAT,                      intent(in)    :: position(:)
  R_TYPE,                     intent(out)   :: interpolated_value
  
  FLOAT, allocatable :: positions(:, :)
  R_TYPE :: interpolated_values(1:1)
  
  PUSH_SUB(X(mesh_interpolation_evaluate))

  SAFE_ALLOCATE(positions(1:this%mesh%sb%dim, 1:1))

  positions(1:this%mesh%sb%dim, 1) = position(1:this%mesh%sb%dim)
  call X(mesh_interpolation_evaluate_vec)(this, 1, values, positions, interpolated_values)
  interpolated_value = interpolated_values(1)
  
  SAFE_DEALLOCATE_A(positions)
  
  POP_SUB(X(mesh_interpolation_evaluate))

end subroutine X(mesh_interpolation_evaluate)

! --------------------------------------------------------------------------------

subroutine X(mesh_interpolation_evaluate_vec)(this, npoints, values, positions, interpolated_values)
  type(mesh_interpolation_t), target, intent(in)    :: this
  integer,                            intent(in)    :: npoints
  R_TYPE,                             intent(in)    :: values(:)
  FLOAT,                              intent(in)    :: positions(:, :)
  R_TYPE,                             intent(out)   :: interpolated_values(:)

  type(mesh_t), pointer :: mesh
  integer :: nm(1:MAX_DIM), ipoint
  FLOAT :: xd(1:MAX_DIM), posrel(1:MAX_DIM)
  R_TYPE :: c00, c10, c01, c11, c0, c1
  R_TYPE :: lvalues(1:8)
  integer, parameter ::  &
    i000 = 1,            &
    i100 = 2,            &
    i010 = 3,            &
    i110 = 4,            &
    i001 = 5,            &
    i101 = 6,            &
    i011 = 7,            &
    i111 = 8
  integer :: pt(1:8), npt, ipt
  logical :: inner_point, boundary_point

  PUSH_SUB(X(mesh_interpolation_evaluate))

  mesh => this%mesh

  ASSERT(ubound(values, dim = 1) == mesh%np_part)
  
  ASSERT(mesh%sb%dim <= 3)

  nm = 0

  do ipoint = 1, npoints
    
    posrel(1:mesh%sb%dim) = positions(1:mesh%sb%dim, ipoint)/mesh%spacing(1:mesh%sb%dim)

    nm(1:mesh%sb%dim) = floor(posrel(1:mesh%sb%dim))

    xd(1:mesh%sb%dim) = posrel(1:mesh%sb%dim) - nm(1:mesh%sb%dim)

    ASSERT(all(xd(1:mesh%sb%dim) >= CNST(0.0)))
    ASSERT(all(xd(1:mesh%sb%dim) <= CNST(1.0)))

    npt = 2**mesh%sb%dim

    ! get the point indices (this could be done in a loop with bit tricks)
    
    pt(i000) = mesh_global_index_from_coords(mesh, [0 + nm(1), 0 + nm(2), 0 + nm(3)])
    pt(i100) = mesh_global_index_from_coords(mesh, [1 + nm(1), 0 + nm(2), 0 + nm(3)])

    if(mesh%sb%dim >= 2) then
      pt(i010) = mesh_global_index_from_coords(mesh, [0 + nm(1), 1 + nm(2), 0 + nm(3)])
      pt(i110) = mesh_global_index_from_coords(mesh, [1 + nm(1), 1 + nm(2), 0 + nm(3)])
    end if

    if(mesh%sb%dim >= 3) then
      pt(i001) = mesh_global_index_from_coords(mesh, [0 + nm(1), 0 + nm(2), 1 + nm(3)])
      pt(i101) = mesh_global_index_from_coords(mesh, [1 + nm(1), 0 + nm(2), 1 + nm(3)])
      pt(i011) = mesh_global_index_from_coords(mesh, [0 + nm(1), 1 + nm(2), 1 + nm(3)])
      pt(i111) = mesh_global_index_from_coords(mesh, [1 + nm(1), 1 + nm(2), 1 + nm(3)])
    end if

     if(mesh%parallel_in_domains) then
       do ipt = 1, npt
         pt(ipt) = vec_global2local(mesh%vp, pt(ipt), mesh%vp%partno)
         lvalues(ipt) = CNST(0.0)
         boundary_point = pt(ipt) > mesh%np + mesh%vp%np_ghost
         inner_point = pt(ipt) > 0 .and. pt(ipt) <= mesh%np
         if(boundary_point .or. inner_point) lvalues(ipt) = values(pt(ipt))
      end do
    else
      do ipt=1, npt
        lvalues(ipt) = values(pt(ipt))
      end do
    end if

    select case(mesh%sb%dim)
    case(3)

      ! trilinear interpolation : http://en.wikipedia.org/wiki/Trilinear_interpolation
      c00 = lvalues(i000)*(CNST(1.0) - xd(1)) + lvalues(i100)*xd(1)
      c10 = lvalues(i010)*(CNST(1.0) - xd(1)) + lvalues(i110)*xd(1)
      c01 = lvalues(i001)*(CNST(1.0) - xd(1)) + lvalues(i101)*xd(1)
      c11 = lvalues(i011)*(CNST(1.0) - xd(1)) + lvalues(i111)*xd(1)
      c0 = c00*(CNST(1.0) - xd(2)) + c10*xd(2)
      c1 = c01*(CNST(1.0) - xd(2)) + c11*xd(2)
      interpolated_values(ipoint) = c0*(CNST(1.0) - xd(3)) + c1*xd(3)

    case(2)
      
      ! bilinear interpolation: http://en.wikipedia.org/wiki/Bilinear_interpolation
      c0 = lvalues(i000)*(CNST(1.0) - xd(1)) + lvalues(i100)*xd(1)
      c1 = lvalues(i010)*(CNST(1.0) - xd(1)) + lvalues(i110)*xd(1)
      interpolated_values(ipoint) = c0*(CNST(1.0) - xd(2)) + c1*xd(2)

    case(1)
      
      ! linear interpolation
      interpolated_values(ipoint) = lvalues(i000)*(CNST(1.0) - xd(1)) + lvalues(i100)*xd(1)

    case default

      ! this line is simply to suppress warnings about the result being uninitialized
      interpolated_values(ipoint) = M_ZERO
      call messages_not_implemented("mesh interpolation in dimensions other than 1,2,3")

    end select

  end do

  call mesh%allreduce(interpolated_values, npoints)
  
  POP_SUB(X(mesh_interpolation_evaluate))

end subroutine X(mesh_interpolation_evaluate_vec)

! --------------------------------------------------------------

subroutine X(mesh_interpolation_test)(mesh)
  type(mesh_t),     intent(in) :: mesh

  integer, parameter :: ntest_points = 20
  integer :: ip, idir, itest
  R_TYPE :: coeff(1:MAX_DIM), calculated, interpolated, interpolated2(1:ntest_points)
  R_TYPE, allocatable :: ff(:)
  type(c_ptr)  :: random_gen_pointer
  type(mesh_interpolation_t) :: interp
  FLOAT :: xx(1:MAX_DIM, ntest_points)

  PUSH_SUB(X(mesh_interpolation_test))

  SAFE_ALLOCATE(ff(1:mesh%np_part))

  call loct_ran_init(random_gen_pointer)

  ! generate the field to be interpolated
  if(mesh%mpi_grp%rank == 0) then
    do idir = 1, mesh%sb%dim
#ifdef R_TCOMPLEX
      coeff(idir) = TOCMPLX(loct_ran_gaussian(random_gen_pointer, CNST(100.0)), loct_ran_gaussian(random_gen_pointer, CNST(100.0)))
#else
      coeff(idir) = loct_ran_gaussian(random_gen_pointer, CNST(100.0))
#endif
    end do
  end if

  if(mesh%parallel_in_domains) then
#ifdef HAVE_MPI
    call MPI_Bcast(coeff, mesh%sb%dim, R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif
  end if
  
  do ip = 1, mesh%np_part
    ff(ip) = sum(coeff(1:mesh%sb%dim)*mesh%x(ip, 1:mesh%sb%dim))
  end do
 
  ! generate the points
  if(mesh%mpi_grp%rank == 0) then
    do itest = 1, ntest_points
      ip = 1 + nint(loct_ran_flat(random_gen_pointer, CNST(0.0), CNST(1.0))*(mesh%np - 1))
      do idir = 1, mesh%sb%dim
        xx(idir, itest) = mesh%x(ip, idir) + mesh%spacing(idir)*loct_ran_flat(random_gen_pointer, CNST(-1.0), CNST(1.0))
      end do
      call messages_write('Point')
      call messages_write(itest)
      do idir = 1, mesh%sb%dim
        call messages_write(xx(idir, itest))
      end do
      call messages_info()
    end do
  end if

  if(mesh%parallel_in_domains) then
#ifdef HAVE_MPI
    call MPI_Bcast(xx, MAX_DIM*ntest_points, R_MPITYPE, 0, mesh%mpi_grp%comm, mpi_err)
#endif
  end if

  call loct_ran_end(random_gen_pointer)

  call mesh_interpolation_init(interp, mesh)

  ! test the scalar routine
  do itest = 1, ntest_points
    
    calculated = sum(coeff(1:mesh%sb%dim)*xx(1:mesh%sb%dim, itest))
    call mesh_interpolation_evaluate(interp, ff, xx(:, itest), interpolated)
    
    call messages_write('Random point')
#ifdef R_TREAL
    call messages_write(' real')
#else
    call messages_write(' complex')
#endif
    call messages_write(itest, fmt = '(i3)')
    call messages_write(' error:')
    call messages_write(abs(calculated - interpolated), fmt = '(e8.2)', align_left = .true.)
    call messages_info()

  end do

  ! and the vectorial one

  call mesh_interpolation_evaluate(interp, ntest_points, ff, xx, interpolated2)

  do itest = 1, ntest_points

    calculated = sum(coeff(1:mesh%sb%dim)*xx(1:mesh%sb%dim, itest))
    
    call messages_write('Random point')
#ifdef R_TREAL
    call messages_write(' real')
#else
    call messages_write(' complex')
#endif
    call messages_write(ntest_points + itest, fmt = '(i3)')
    call messages_write(' error:')
    call messages_write(abs(calculated - interpolated2(itest)), fmt = '(e8.2)', align_left = .true.)
    call messages_info()

  end do
  
  call mesh_interpolation_end(interp)
  

  SAFE_DEALLOCATE_A(ff)

  POP_SUB(X(mesh_interpolation_test))
end subroutine X(mesh_interpolation_test)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
