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

!> Computes all eigenvalues and eigenvectors of a real symmetric square matrix A.

subroutine X(eigensolve_scalapack)(n, a, eigenvalues, bof, proc_grid, err_code)
  integer, intent(in)                 :: n              !< Number of rows/columns
  R_TYPE,  intent(inout)              :: a(:, :)        !< Square matrix
  FLOAT,   intent(out)                :: eigenvalues(:) !< Selected eigenvalues
  logical, optional, intent(inout)    :: bof            !< Bomb on failure.
  type(blacs_proc_grid_t), intent(in) :: proc_grid      !< BLACS processor grid
  integer, optional, intent(out)      :: err_code

#ifdef HAVE_SCALAPACK
  logical              :: bof_
  integer              :: info, lwork, blockrow, blockcol, begining_row, begining_col, &
       psi_desc(BLACS_DLEN), blacs_info, i_loc, j_loc, liwork, eigenvectors_failing, &
       eigenvalues_size, eigenvectors_computed, nn,np0, mq0, neig
  R_TYPE, allocatable  :: work(:), orthonormal_eigenvectors(:,:)
  FLOAT                :: upper_bound, lower_bound, error, gap
  R_TYPE               :: a_loc(n,n)
  integer, allocatable :: iwork(:), eigenvectors_cluster(:)
  
#ifdef R_TCOMPLEX
  FLOAT, allocatable                  :: rwork(:)
  integer                             :: lrwork
#endif
#endif

  PUSH_SUB(X(eigensolve_scalapack))
  call profiling_in(eigensolver_prof, "DENSE_EIGENSOLVER_SCALAPACK")

#ifdef HAVE_SCALAPACK
  bof_ = .true.
  if(present(bof)) then
    bof_ = bof
  end if

  blockrow = 32 ! arbitrary. But usually the best size between 16 and 64
  blockcol = 32 ! arbitrary. But usually the best size between 16 and 64
  begining_row = proc_grid%iam + 1 + blockrow
  begining_col = proc_grid%iam + 1 + blockcol
  neig =n ! because we ask for range='V'
  nn = max(n, blockcol, 2)
  np0 = numroc(nn, blockrow, 0, 0, proc_grid%nprow)
  mq0 = numroc(max(neig, blockrow, 2),blockrow, 0, 0, proc_grid%npcol) 
  lwork = 5*n + max(5 * nn, np0 * mq0 + 2 * blockrow * blockrow) + iceil(neig, proc_grid%nprow * proc_grid%npcol) * nn 
  
  lower_bound = M_ZERO
  upper_bound = real(n)
  error = real(0.001) ! it has to be tuned

  SAFE_ALLOCATE(work(1:lwork))  
  liwork = 6 * max(n,proc_grid%nprow * proc_grid%npcol + 1, 4)
  SAFE_ALLOCATE(iwork(1:liwork))
  SAFE_ALLOCATE(eigenvectors_cluster(1:2 * proc_grid%npcol * proc_grid%nprow))
 
  ! DISTRIBUTE THE MATRIX ON THE PROCESS GRID
  ! Initialize the descriptor array for the main matrices (ScaLAPACK)
  call descinit(psi_desc(1), n, n * n, blockrow, blockcol, 0, 0, proc_grid%context, n / proc_grid%nprocs, blacs_info)
  ! distribute to local matrices. 
  call pcelset(a,begining_row,begining_col,psi_desc,a(1,1))

#ifdef R_TREAL
  call scalapack_syev('V','A','U', n, a(1,1), begining_row,begining_col, psi_desc(1), lower_bound, upper_bound,    &
       0, 0, error, eigenvalues_size, eigenvectors_computed, eigenvalues(1), M_ZERO, orthonormal_eigenvectors(1,1),&
       begining_row, begining_col, psi_desc(1), work(1), lwork, iwork(1), liwork, eigenvectors_failing,            &
       eigenvectors_cluster(1),gap,info)
#else
  lrwork = 4 * n + max(5*nn, np0 * mq0) + iceil(neig, proc_grid%nprow * proc_grid%npcol)*nn
  SAFE_ALLOCATE(rwork(1:lrwork))
  call scalapack_syev( jobz = 'V', &
       range = 'A', &
       uplo = 'U', &
       n = n, &
       a = a(1,1), &
       ia = begining_row, &
       ja = begining_col,&
       desca = psi_desc(1), &
       vl = lower_bound, &
       vu = upper_bound,    &
       il = 0, &
       iu = 0, &
       abstol = error, &
       m = eigenvalues_size, &
       nz = eigenvectors_computed, &
       w = eigenvalues(1), &
       orfac = M_ZERO, &
       z = orthonormal_eigenvectors(1,1),&
       iz = begining_row, &
       jz = begining_col, &
       descz = psi_desc(1), &
       work = work(1), &
       lwork = lwork, &
       rwork = rwork(1), &
       lrwork = lwork, &
       iwork = iwork(1), &
       liwork = liwork, &
       ifail = eigenvectors_failing,  &
       iclustr = eigenvectors_cluster(1), &
       gap = gap, &
       info = info)
  ASSERT(.false.)
#endif

  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(iwork)
  
  if(info /= 0) then
    if(bof_) then
#ifdef RT_REAL
      write(message(1),'(3a,i5)') 'In ' // TOSTRING(X(scalapack_eigensolve)) // ' ScaLAPACK pdsyevx returned error message ', info
#else
      write(message(1),'(3a,i5)') 'In ' // TOSTRING(X(scalapack_eigensolve)) // ' ScaLAPACK pzheevx returned error message ', info
#endif
      call messages_fatal(1)
    else
      if(present(bof)) then
        bof = .true.
      end if
    end if
  else
    if(present(bof)) then
      bof = .false.
    end if
  end if
  if(present(err_code)) then
    err_code = info
  end if

  ! save the result to the original matrix
  a = orthonormal_eigenvectors

#endif

  call profiling_out(eigensolver_prof)
  POP_SUB(X(eigensolve_scalapack))
end subroutine X(eigensolve_scalapack)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
