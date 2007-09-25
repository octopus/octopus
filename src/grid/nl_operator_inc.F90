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

! The next two routines are wrappers around the operator sum.
! The reason for those wrappers is that, by help of the parameter
! list, compilers can generate more efficient code because of they then
! know the arrays not to be shaped.
!
! ---------------------------------------------------------

subroutine X(nl_operator_tune)(op)
  type(nl_operator_t), intent(inout) :: op
  
  R_TYPE, allocatable :: in(:), out(:)
  real(8) :: noperations, flops(OP_MIN:OP_MAX), itime, ftime
  integer :: method, ii, reps, iunit
  character(len=2) :: marker

#ifdef HAVE_MPI
  integer :: ierr, rank
  real(8) :: global_flops(OP_MIN:OP_MAX)
#endif

  call push_sub('nl_operator_inc.Xnl_operator_tune')
  
  !count the total number of floating point operations  

#ifdef R_TCOMPLEX
  noperations =  op%m%np * op%n * M_FOUR
#else
  noperations =  op%m%np * op%n * M_TWO
#endif 
  
  !measure performance of each function
  ALLOCATE(in(1:op%m%np_part), op%m%np_part)
  in(1:op%m%np_part) = M_ZERO
  ALLOCATE(out(1:op%m%np), op%m%np)
  
  flops = M_ZERO

  do method = OP_MIN, OP_MAX

    !skip methods that are not available
#ifdef R_TCOMPLEX
    if (op_is_available(method, M_CMPLX) == 0) cycle
#else
    if (op_is_available(method, M_REAL)  == 0) cycle
#endif
    op%X(function) = method

    reps = 10

    if(in_debug_mode) then
      write(message(1), '(5a)') 'Info: Profiling non local operator: ', trim(op%label), ' - ', &
#ifdef R_TCOMPLEX
           'complex - ', &
#else
           'real - ', &
#endif           
        op_function_name(method)
      call write_info(1)
    end if

    itime = loct_clock()
    do ii = 1, reps
      call X(nl_operator_operate)(op, in, out, ghost_update = .true.)
    end do
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    ftime = loct_clock()

    flops(method) = noperations * reps / (ftime - itime) 
    
  end do

  !choose the best method
  op%X(function) = OP_MIN
  do method = OP_MIN + 1, OP_MAX
    if(flops(method) > flops(op%X(function))) op%X(function) = method
  end do

#ifdef HAVE_MPI
      call MPI_Allreduce(flops, global_flops, OP_MAX-OP_MIN+1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      flops = global_flops
#endif

#ifdef HAVE_MPI      
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if(rank == 0) then
#endif
    !print to file
    if(.not. initialized) call loct_rm('exec/nl_operator_prof')
    initialized = .true.
    
    iunit = io_open('exec/nl_operator_prof', action='write', position='append', is_tmp = .true.)
    
#ifdef R_TCOMPLEX
    write (iunit, '(3a)')   'Operator       = ', trim(op%label), " complex"
#else
    write (iunit, '(3a)')   'Operator       = ', trim(op%label), " real"
#endif
    write (iunit, '(a,i8)') 'Stencil points = ', op%n
    write (iunit, '(a,i8)') 'Grid points    = ', op%m%np
    
    do method = OP_MIN, OP_MAX
#ifdef R_TCOMPLEX
      if (op_is_available(method, M_CMPLX) == 0) cycle
#else
      if (op_is_available(method, M_REAL)  == 0) cycle
#endif
      marker = '  '
      if(method == op%X(function)) marker = '* '
      write (iunit, '(2a, f8.1, a)') marker, op_function_name(method), flops(method)/CNST(1e6), ' MFlops'
    end do
    
    write(iunit, '(a)') " "
    
    call io_close(iunit)
    
#ifdef HAVE_MPI      
  end if
#endif

  call pop_sub()

end subroutine X(nl_operator_tune)

subroutine X(operate)(np, np_part, nn, nri, w, ri, rimap_inv, fi, fo)
  integer, intent(in) :: np
  integer, intent(in) :: np_part
  integer, intent(in) :: nn
  integer, intent(in) :: nri
  FLOAT,   intent(in) :: w(1:nn)
  integer, intent(in) :: ri(1:nn, 1:nri)
  integer, intent(in) :: rimap_inv(0:nri+1)
  R_TYPE,   intent(in) :: fi(1:np_part)
  R_TYPE,   intent(out):: fo(1:np) 
  
  integer :: ll, ii

  do ll = 1, nri
    do ii = rimap_inv(ll-1)+1, rimap_inv(ll)
      fo(ii) = sum(w(1:nn)  * fi(ii + ri(1:nn, ll)) )
    end do
  end do

end subroutine X(operate)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
