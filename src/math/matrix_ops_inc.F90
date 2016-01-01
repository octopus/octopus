!! Copyright (C) 2015 X. Andrade
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
!! $Id: matrix_inc.F90 14808 2015-11-21 05:43:46Z xavier $


subroutine X(to_submatrix)(this, bsize, proc, nproc, submatrix)
  type(matrix_t),           intent(in)    :: this
  integer,                  intent(in)    :: bsize(:)
  integer,                  intent(in)    :: proc(:)
  integer,                  intent(in)    :: nproc(:)
  R_TYPE,                   intent(inout) :: submatrix(:, :)

  integer :: ic1, ic2, ib1, ib2, i1, i2

  PUSH_SUB(X(to_submatrix))
  
  ic1 = 1
  do ib1 = 1 + bsize(1)*proc(1), this%dim(1), bsize(1)*nproc(1)
    do i1 = ib1, min(ib1 - 1 + bsize(1), this%dim(1))
      
      ic2 = 1
      do ib2 = 1 + bsize(2)*proc(2), this%dim(2), bsize(2)*nproc(2)
        do i2 = ib2, min(ib2 - 1 + bsize(2), this%dim(2))
          
          submatrix(ic1, ic2) = this%X(mat)(i1, i2)
          
          ic2 = ic2 + 1
        end do
      end do
      
      ic1 = ic1 + 1
    end do
  end do

  POP_SUB(X(to_submatrix))
end subroutine X(to_submatrix)

! ---------------------------------------------------------------

subroutine X(from_submatrix)(this, bsize, proc, nproc, submatrix)
  type(matrix_t),           intent(inout) :: this
  integer,                  intent(in)    :: bsize(:)
  integer,                  intent(in)    :: proc(:)
  integer,                  intent(in)    :: nproc(:)
  R_TYPE,                   intent(in)    :: submatrix(:, :)

  integer :: ic1, ic2, ib1, ib2, i1, i2

  PUSH_SUB(X(from_submatrix))
  
  ic1 = 1
  do ib1 = 1 + bsize(1)*proc(1), this%dim(1), bsize(1)*nproc(1)
    do i1 = ib1, min(ib1 - 1 + bsize(1), this%dim(1))
      
      ic2 = 1
      do ib2 = 1 + bsize(2)*proc(2), this%dim(2), bsize(2)*nproc(2)
        do i2 = ib2, min(ib2 - 1 + bsize(2), this%dim(2))
          
          this%X(mat)(i1, i2) = submatrix(ic1, ic2)
          
          ic2 = ic2 + 1
        end do
      end do
      
      ic1 = ic1 + 1
    end do
  end do

  POP_SUB(X(from_submatrix))
end subroutine X(from_submatrix)

! ---------------------------------------------------------------

subroutine X(matrix_gather)(this, bsize, proc_grid, submatrix)
  type(matrix_t),           intent(inout) :: this
  integer,                  intent(in)    :: bsize(:)
  type(blacs_proc_grid_t),  intent(in)    :: proc_grid
  R_TYPE,                   intent(in)    :: submatrix(:, :)

  integer :: ip1, ip2, nr1r, nr2r, iproc
  R_TYPE, allocatable :: submatrix_rem(:, :)

  PUSH_SUB(matrix_gather)

  call profiling_in(prof_gather, "MATRIX_GATHER")
  
  do ip1 = 1, proc_grid%nprow      
    do ip2 = 1, proc_grid%npcol

      nr1r = numroc(this%dim(1), bsize(1), ip1 - 1, 0, proc_grid%nprow)
      nr2r = numroc(this%dim(2), bsize(2), ip2 - 1, 0, proc_grid%npcol)

      if(nr1r == 0 .or. nr2r == 0) cycle

      SAFE_ALLOCATE(submatrix_rem(1:nr1r, 1:nr2r))

      if(ip1 - 1 == proc_grid%myrow .and. ip2 - 1 == proc_grid%mycol) then
        submatrix_rem(1:nr1r, 1:nr2r) = submatrix(1:nr1r, 1:nr2r)
      end if

      iproc = proc_grid%usermap(ip1, ip2)

#ifdef HAVE_MPI
      call MPI_Bcast(submatrix_rem(1, 1), nr1r*nr2r, R_MPITYPE, iproc, this%mpi_grp%comm, mpi_err)
      call MPI_Barrier(this%mpi_grp%comm, mpi_err)
#endif
      
      call X(from_submatrix)(this, bsize, (/ip1, ip2/) - 1, (/proc_grid%nprow, proc_grid%npcol/), submatrix_rem)

      SAFE_DEALLOCATE_A(submatrix_rem)

    end do
  end do

  call profiling_out(prof_gather)

  POP_SUB(matrix_gather)
end subroutine X(matrix_gather)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
