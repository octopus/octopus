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
!! $Id: matrix.F90 7610 2011-03-25 03:41:01Z xavier $

#include "global.h"

module matrix_ops_oct_m
  use blacs_proc_grid_oct_m
  use global_oct_m
  use lapack_oct_m
  use matrix_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use types_oct_m
  use utils_oct_m

  implicit none

  private
  public ::                          &
    matrix_diagonalize_hermitian

  type(profile_t), save :: prof_gather
  
contains

  ! ---------------------------------------------------------
  subroutine matrix_diagonalize_hermitian(this, eigenvalues, metric)
    type(matrix_t),           intent(inout) :: this            !< The matrix to diagonalize
    FLOAT,                    intent(out)   :: eigenvalues(:)  !< The eigenvalues of the matrix
    type(matrix_t), optional, intent(inout) :: metric          !< If present, a generalized eigenvalue problem is solved

    FLOAT :: worksize
    CMPLX :: zworksize
    FLOAT, allocatable :: work(:)
    CMPLX, allocatable :: zwork(:)
    integer :: info

    PUSH_SUB(matrix_diagonalize_hermitian)

    ASSERT(ubound(eigenvalues, dim = 1) == this%dim(1))
    ASSERT(this%dim(1) == this%dim(2))

    if(present(metric)) then
      ASSERT(all(this%dim(1:2) == metric%dim(1:2)))
    end if

#ifdef HAVE_SCALAPACK
    if(this%mpi_grp%size > 1) then
      call diagonalize_hermitian_scalapack(this, eigenvalues, metric)
      POP_SUB(matrix_diagonalize_hermitian)
      return
    end if
#endif
    
    if(this%mpi_grp%rank == 0) then

      if(present(metric)) then

        if(matrix_type(this) == TYPE_FLOAT) then

          call dsygv(itype = 1, jobz = 'V', uplo = 'U', n = this%dim(1), a = this%dmat(1, 1), lda = this%dim(1), &
            b = metric%dmat(1, 1), ldb = metric%dim(1), &
            w = eigenvalues(1), work = worksize, lwork = -1, info = info)

          SAFE_ALLOCATE(work(1:int(worksize)))

          call dsygv(itype = 1, jobz = 'V', uplo = 'U', n = this%dim(1), a = this%dmat(1, 1), lda = this%dim(1), &
            b = metric%dmat(1, 1), ldb = metric%dim(1), &
            w = eigenvalues(1), work = work(1), lwork = int(worksize), info = info)

        else

          SAFE_ALLOCATE(work(1:3*this%dim(1) -2 ))

          call zhegv(itype = 1, jobz = 'V', uplo = 'U', n = this%dim(1), a = this%zmat(1, 1), lda = this%dim(1), &
            b = metric%zmat(1, 1), ldb = metric%dim(1), &
            w = eigenvalues(1), work = zworksize, lwork = -1, rwork = work(1), info = info)

          SAFE_ALLOCATE(zwork(1:int(zworksize)))

          call zhegv(itype = 1, jobz = 'V', uplo = 'U', n = this%dim(1), a = this%zmat(1, 1), lda = this%dim(1), &
            b = metric%zmat(1, 1), ldb = metric%dim(1), &
            w = eigenvalues(1), work = zwork(1), lwork = int(zworksize), rwork = work(1), info = info)

        end if

      else

        if(matrix_type(this) == TYPE_FLOAT) then

          call dsyev(jobz = 'V', uplo = 'U', n = this%dim(1), a = this%dmat(1, 1), lda = this%dim(1), &
            w = eigenvalues(1), work = worksize, lwork = -1, info = info)

          SAFE_ALLOCATE(work(1:int(worksize)))

          call dsyev(jobz = 'V', uplo = 'U', n = this%dim(1), a = this%dmat(1, 1), lda = this%dim(1), &
            w = eigenvalues(1), work = work(1), lwork = int(worksize), info = info)

        else

          SAFE_ALLOCATE(work(1:3*this%dim(1) - 2))

          call zheev(jobz = 'V', uplo = 'U', n = this%dim(1), a = this%zmat(1, 1), lda = this%dim(1), &
            w = eigenvalues(1), work = zworksize, lwork = -1, rwork = work(1), info = info)

          SAFE_ALLOCATE(zwork(1:int(zworksize)))

          call zheev(jobz = 'V', uplo = 'U', n = this%dim(1), a = this%zmat(1, 1), lda = this%dim(1), &
            w = eigenvalues(1), work = zwork(1), lwork = int(zworksize), rwork = work(1), info = info)

        end if

      end if

      SAFE_DEALLOCATE_A(work)
      SAFE_DEALLOCATE_A(zwork)
      
    end if

    if(this%mpi_grp%size > 1) then
#ifdef HAVE_MPI
      if(matrix_type(this) == TYPE_FLOAT) then
        call MPI_Bcast(this%dmat(1, 1), product(this%dim(1:2)), MPI_FLOAT, 0, this%mpi_grp%comm, mpi_err)
        call MPI_Barrier(this%mpi_grp%comm, mpi_err)
      else
        call MPI_Bcast(this%zmat(1, 1), product(this%dim(1:2)), MPI_CMPLX, 0, this%mpi_grp%comm, mpi_err)
        call MPI_Barrier(this%mpi_grp%comm, mpi_err)
      end if
      call MPI_Bcast(eigenvalues(1), this%dim(1), MPI_FLOAT, 0, this%mpi_grp%comm, mpi_err)
      call MPI_Barrier(this%mpi_grp%comm, mpi_err)
#endif
    end if

    POP_SUB(matrix_diagonalize_hermitian)
  end subroutine matrix_diagonalize_hermitian

  ! -----------------------------------------------------------------------------------------------------------
  
  subroutine diagonalize_hermitian_scalapack(this, eigenvalues, metric)
    type(matrix_t),           intent(inout) :: this            !< The matrix to diagonalize
    FLOAT,                    intent(out)   :: eigenvalues(:)  !< The eigenvalues of the matrix
    type(matrix_t), optional, intent(inout) :: metric          !< If present, a generalized eigenvalue problem is solved

    integer :: ndiv, np1, np2, info, nb1, nb2, nr1, nr2
    integer, allocatable :: div(:)
    type(blacs_proc_grid_t) :: proc_grid
    integer :: desc(BLACS_DLEN)
    FLOAT, allocatable :: da(:, :), devec(:, :), dwork(:)
    CMPLX, allocatable :: za(:, :), zevec(:, :), zwork(:), rwork(:)
    FLOAT :: dworksize
    CMPLX :: zworksize, lrwork
    
    PUSH_SUB(matrix_diagonalize_hermitian_scalapack)

    ASSERT(.not. present(metric))
    
    ! divide the processors
    ndiv = 10
    SAFE_ALLOCATE(div(1:ndiv))
    call get_divisors(this%mpi_grp%size, ndiv, div)
    
    np1 = div((ndiv + 1)/2)
    np2 = this%mpi_grp%size/np1

    SAFE_DEALLOCATE_A(div)

    ASSERT(np1*np2 == this%mpi_grp%size)

    nb1 = 10
    nb2 = 10
    
    call blacs_proc_grid_init(proc_grid, this%mpi_grp, procdim = (/np1, np2/))

    nr1 = max(1, numroc(this%dim(1), nb1, proc_grid%myrow, 0, proc_grid%nprow))
    nr2 = max(1, numroc(this%dim(2), nb2, proc_grid%mycol, 0, proc_grid%npcol))

#ifdef HAVE_SCALAPACK
    call descinit(desc(1), this%dim(1), this%dim(2), nb1, nb2, 0, 0, proc_grid%context, nr1, info)

    if(matrix_type(this) == TYPE_FLOAT) then

      SAFE_ALLOCATE(da(1:nr1, 1:nr2))
      SAFE_ALLOCATE(devec(1:nr1, 1:nr2))

      call dto_submatrix(this, (/nb1, nb2/), (/proc_grid%myrow, proc_grid%mycol/), (/proc_grid%nprow, proc_grid%npcol/), da)

      call pdsyev(jobz = 'V', uplo = 'U', n = this%dim(1), a = da(1, 1) , ia = 1, ja = 1, desca = desc(1), &
        w = eigenvalues(1), z = devec(1, 1), iz = 1, jz = 1, descz = desc(1), work = dworksize, lwork = -1, info = info)
      
      if(info /= 0) then
        write(message(1),'(a,i6)') "ScaLAPACK pdsyev workspace query failure, error code = ", info
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(dwork(1:int(dworksize)))

      call pdsyev(jobz = 'V', uplo = 'U', n = this%dim(1), a = da(1, 1) , ia = 1, ja = 1, desca = desc(1), &
        w = eigenvalues(1), z = devec(1, 1), iz = 1, jz = 1, descz = desc(1), work = dwork(1), lwork = int(dworksize), info = info)

      if(info /= 0) then
        call messages_write('ScaLAPACK pdsyev call failure, error code = ')
        call messages_write(info)
        call messages_fatal()
      end if
      
      SAFE_DEALLOCATE_A(dwork)

      call dmatrix_gather(this, (/nb1, nb2/), proc_grid, devec)

      SAFE_DEALLOCATE_A(da)
      SAFE_DEALLOCATE_A(devec)      
      
    else if(matrix_type(this) == TYPE_CMPLX) then

      SAFE_ALLOCATE(za(1:nr1, 1:nr2))
      SAFE_ALLOCATE(zevec(1:nr1, 1:nr2))
      
      call zto_submatrix(this, (/nb1, nb2/), (/proc_grid%myrow, proc_grid%mycol/), (/proc_grid%nprow, proc_grid%npcol/), za)

      lrwork = M_z0
      
      call pzheev(jobz = 'V', uplo = 'U', n = this%dim(1), a = za(1, 1) , ia = 1, ja = 1, desca = desc(1), &
        w = eigenvalues(1), z = zevec(1, 1), iz = 1, jz = 1, descz = desc(1), work = zworksize, lwork = -1, &
        rwork = lrwork, lrwork = -1, info = info)
      
      if(info /= 0) then
        write(message(1),'(a,i6)') "ScaLAPACK pdsyev workspace query failure, error code = ", info
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(zwork(1:int(zworksize)))
      SAFE_ALLOCATE(rwork(1:int(lrwork)))

      call pzheev(jobz = 'V', uplo = 'U', n = this%dim(1), a = za(1, 1) , ia = 1, ja = 1, desca = desc(1), &
        w = eigenvalues(1), z = zevec(1, 1), iz = 1, jz = 1, descz = desc(1), work = zwork(1), lwork = int(zworksize), &
        rwork = rwork(1), lrwork = int(lrwork), info = info)

      if(info /= 0) then
        call messages_write('ScaLAPACK pzheev call failure, error code = ')
        call messages_write(info)
        call messages_fatal()
      end if
      
      SAFE_DEALLOCATE_A(zwork)
      SAFE_DEALLOCATE_A(rwork)

      call zmatrix_gather(this, (/nb1, nb2/), proc_grid, zevec)

      SAFE_DEALLOCATE_A(za)
      SAFE_DEALLOCATE_A(zevec)      

    else

      ASSERT(.false.)
      
    end if
#endif
    
    call blacs_proc_grid_end(proc_grid)

    POP_SUB(matrix_diagonalize_hermitian_scalapack)
  end subroutine diagonalize_hermitian_scalapack
    
#include "undef.F90"
#include "real.F90"

#include "matrix_ops_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "matrix_ops_inc.F90"

end module matrix_ops_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
