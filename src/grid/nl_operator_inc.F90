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
  real(8) :: noperations, flops(OP_MIN:OP_MAX), itime, ftime, bwidth(OP_MIN:OP_MAX), dvolume
  integer :: method, ii, reps, iunit
  character(len=2) :: marker

#ifdef R_TCOMPLEX
  character(len=*), parameter :: type = 'complex'
  FLOAT, parameter            :: R_OPS = M_TWO
#else
  character(len=*), parameter :: type = 'real'
  FLOAT, parameter            :: R_OPS = M_ONE
#endif

#ifdef SINGLE_PRECISION
  FLOAT, parameter :: R_SIZE = M_FOUR
#else
  FLOAT, parameter :: R_SIZE = M_EIGHT
#endif

#ifdef HAVE_MPI
  integer :: ierr, rank
  real(8) :: global_flops(OP_MIN:OP_MAX)
#endif

  call profiling_in(nl_tuning_profile, "NL_OPERATOR_TUNE")
  call push_sub('nl_operator_inc.Xnl_operator_tune')
  
  !count the total number of floating point operations  
  noperations =  op%m%np * op%n * M_TWO * R_OPS

  !the volume of data that has to be moved in bytes
  dvolume = (op%m%np_part + op%m%np) * R_OPS * R_SIZE + (op%nri + 1)* op%n * M_FOUR

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
      write(message(1), '(6a)') 'Info: Profiling non local operator: ', trim(op%label), ' - ', &
        trim(type), ' - ', op_function_name(method)
      call write_info(1)
    end if

    itime = loct_clock()
    do ii = 1, reps
      call X(nl_operator_operate)(op, in, out, ghost_update = .true., profile = .true.)
    end do
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    ftime = loct_clock()

    flops(method) = noperations * reps / (ftime - itime) 
    bwidth(method) = dvolume * reps / (ftime - itime)
    
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
    write (iunit, '(a,i8)') 'Grid points    = ', op%m%np
    write (iunit, '(a,i8)') 'Stencils       = ', op%nri
    write (iunit, '(a,i8)') 'Stencil points = ', op%n
    
    do method = OP_MIN, OP_MAX
#ifdef R_TCOMPLEX
      if (op_is_available(method, M_CMPLX) == 0) cycle
#else
      if (op_is_available(method, M_REAL)  == 0) cycle
#endif
      marker = '  '
      if(method == op%X(function)) marker = '* '
      write (iunit, '(2a, f8.1, a, f8.1, a)') &
           marker, op_function_name(method), &
           flops(method)/CNST(1e6), ' MFlops  ',&
           bwidth(method)/CNST(1e6), ' MBytes/s'
    end do
    
    write(iunit, '(a)') " "
    
    call io_close(iunit)
    
#ifdef HAVE_MPI      
  end if
#endif

  call pop_sub()
  call profiling_out(nl_tuning_profile)

end subroutine X(nl_operator_tune)

subroutine X(operate)(nn, nri, w, ri, imin, imax, fi, fo)
  integer,          intent(in)  :: nn
  integer,          intent(in)  :: nri
  FLOAT,            intent(in)  :: w(:)
  integer,          intent(in)  :: ri(:, :)
  integer,          intent(in)  :: imin(:)
  integer,          intent(in)  :: imax(:)
  R_TYPE,           intent(in)  :: fi(:)
  R_TYPE,           intent(out) :: fo(:) 
  
  integer :: ll, ii
  
  !$omp do private(ii)
  do ll = 1, nri
    do ii = imin(ll) + 1, imax(ll)
      fo(ii) = sum(w(1:nn)*fi(ii + ri(1:nn, ll)))
    end do
  end do
  !$omp end do nowait

end subroutine X(operate)

! ---------------------------------------------------------
subroutine X(nl_operator_operate)(op, fi, fo, ghost_update, profile, points)
  R_TYPE,              intent(inout) :: fi(:)  ! fi(op%np)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)  ! fo(op%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points
  
  integer :: ii, nn, nri
#if defined(HAVE_MPI)
  logical :: update
#endif
  real(8) :: ws(100)
  logical :: profile_
  integer :: nri_loc, ini
  integer, pointer :: imin(:), imax(:), ri(:, :)

  profile_ = .true. 
  if(present(profile)) profile_ = profile
  
  if(profile_) call profiling_in(nl_operate_profile, "NL_OPERATOR")
  call push_sub('nl_operator.Xnl_operator_operate')
  
#if defined(HAVE_MPI)
  if(present(ghost_update)) then
    update = ghost_update
  else
    update = .true.
  end if

  if(op%m%parallel_in_domains.and.update) then
    call X(vec_ghost_update)(op%m%vp, fi)
  end if
#endif

  if(.not. present(points)) then
    nri  =  op%nri
    imin => op%rimap_inv(1:)
    imax => op%rimap_inv(2:)
    ri   => op%ri
  else
    select case(points)
    case(INNER)
      nri  =  op%inner%nri
      imin => op%inner%imin
      imax => op%inner%imax
      ri   => op%inner%ri
    case(OUTER)
      nri  =  op%outer%nri
      imin => op%outer%imin
      imax => op%outer%imax
      ri   => op%outer%ri
    end select
  end if

  !$omp parallel private(ini, nri_loc, ws)
  nn = op%n
#ifdef R_TCOMPLEX
  if(op%cmplx_op) then
    if(op%const_w) then
      !$omp do
      do ii = 1, op%np
        fo(ii) = sum(cmplx(op%w_re(1:nn, 1),  op%w_im(1:nn, 1))  * fi(ii + op%ri(1:nn, op%rimap(ii))) )
      end do
      !$omp end do
    else
      !$omp do
      do ii = 1, op%np
        fo(ii) = sum(cmplx(op%w_re(1:nn, ii), op%w_im(1:nn, ii)) * fi(ii + op%ri(1:nn, op%rimap(ii))) )
      end do
      !$omp end do
    end if
  else
#endif
    if(op%const_w) then
#ifdef USE_OMP
      call divide_range(nri, omp_get_thread_num(), omp_get_num_threads(), ini, nri_loc)
#else 
      ini = 1
      nri_loc = nri
#endif

      select case(op%X(function))
      case(OP_FORTRAN)
        call X(operate)(nn, nri, op%w_re(:, 1), ri, imin, imax, fi, fo)
      case(OP_C)
        call X(operate_ri)(nn, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), fi(1), fo(1))
      case(OP_VEC)
        call X(operate_ri_vec)(nn, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), fi(1), fo(1))
      case(OP_AS)
        call X(operate_as)(nn, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), fi(1), fo(1), ws(1))
      end select

    else
      
      !$omp do
      do ii = 1, op%np
        fo(ii) = sum(op%w_re(1:nn, ii) * fi(ii + op%ri(1:nn, op%rimap(ii))))
      end do
      !$omp end do
    end if
#ifdef R_TCOMPLEX
  end if
#endif
  !$omp end parallel

  call pop_sub()
  if(profile_) call profiling_out(nl_operate_profile)
end subroutine X(nl_operator_operate)


! ---------------------------------------------------------
subroutine X(nl_operator_operate_diag)(op, fo)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)

  integer :: ii, nn, jj
  
  call push_sub('nl_operator.Xnl_operator_operate_diag')
  
  nn = op%n
  
  if(op%cmplx_op) then
#ifdef R_TCOMPLEX
    if(op%const_w) then
      fo(1:op%np) = cmplx(op%w_re(op%stencil_center, 1), op%w_im(op%stencil_center, 1))
    else
      fo(1:op%np) = cmplx(op%w_re(op%stencil_center, 1:op%np), op%w_im(op%stencil_center, 1:op%np))
    end if
#endif
  else
    if(op%const_w) then
      fo(1:op%np) = op%w_re(op%stencil_center, 1)
    else
      fo(1:op%np) = op%w_re(op%stencil_center, 1:op%np)
    end if
  end if
  
  call pop_sub()
  
end subroutine X(nl_operator_operate_diag)



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
