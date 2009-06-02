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
! list, compilers can generate more efficient code because then
! know the arrays not to be shaped.
!
! ---------------------------------------------------------

subroutine X(nl_operator_tune)(op, best)
  type(nl_operator_t), intent(inout) :: op
  FLOAT,  optional,    intent(out)   :: best

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

  call push_sub('nl_operator_inc.Xnl_operator_tune')
  
  !count the total number of floating point operations  
  noperations =  op%m%np*op%stencil%size*M_TWO*R_OPS

  !the volume of data that has to be moved in bytes
  dvolume = (op%m%np_part + op%m%np) * R_OPS * R_SIZE + (op%nri + 1)*op%stencil%size*M_FOUR

  !measure performance of each function
  SAFE_ALLOCATE(in(1:op%m%np_part))
  in(1:op%m%np_part) = M_ZERO
  SAFE_ALLOCATE(out(1:op%m%np))
  
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
      write(message(1), '(6a)') 'Info: Profiling non-local operator: ', trim(op%label), ' - ', &
        trim(type), ' - ', op_function_name(method)
      call write_info(1)
    end if

    itime = loct_clock()
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    do ii = 1, reps
      call X(nl_operator_operate)(op, in, out, ghost_update = .false., profile = .false.)
    end do
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    ftime = loct_clock()

    flops(method) = noperations * reps / (ftime - itime) 
    bwidth(method) = dvolume * reps / (ftime - itime)
    
  end do

  SAFE_DEALLOCATE_A(in)
  SAFE_DEALLOCATE_A(out)

#ifdef HAVE_MPI
      call MPI_Allreduce(flops, global_flops, OP_MAX-OP_MIN+1, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      flops = global_flops
#endif

  !choose the best method
  op%X(function) = OP_MIN
  do method = OP_MIN + 1, OP_MAX
    if(flops(method) > flops(op%X(function))) op%X(function) = method
  end do

  if(present(best)) best = flops(op%X(function))/CNST(1e6)

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
    write (iunit, '(a,i8)') 'Stencil points = ', op%stencil%size
    
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

end subroutine X(nl_operator_tune)

subroutine X(operate)(nn, nri, w, ri, imin, imax, fi, fo, wim)
  integer,          intent(in)  :: nn
  integer,          intent(in)  :: nri
  FLOAT,            intent(in)  :: w(:)
  integer,          intent(in)  :: ri(:, :)
  integer,          intent(in)  :: imin(:)
  integer,          intent(in)  :: imax(:)
  R_TYPE,           intent(in)  :: fi(:)
  R_TYPE,           intent(out) :: fo(:) 
  FLOAT, optional,  intent(in)  :: wim(:)

  integer :: ll, ii
  
  if( .not. present(wim)) then
    
    !$omp parallel do private(ii)
    do ll = 1, nri
      do ii = imin(ll) + 1, imax(ll)
        fo(ii) = sum(w(1:nn)*fi(ii + ri(1:nn, ll)))
      end do
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(ii)
    do ll = 1, nri
      do ii = imin(ll) + 1, imax(ll)
        fo(ii) = sum(cmplx(w(1:nn), wim(1:nn))*fi(ii + ri(1:nn, ll)))
      end do
    end do
    !$omp end parallel do

  end if
end subroutine X(operate)


subroutine X(operate_nc)(nn, nri, w, ri, imin, imax, fi, fo, wim)
  integer,          intent(in)  :: nn
  integer,          intent(in)  :: nri
  FLOAT,            intent(in)  :: w(:, :)
  integer,          intent(in)  :: ri(:, :)
  integer,          intent(in)  :: imin(:)
  integer,          intent(in)  :: imax(:)
  R_TYPE,           intent(in)  :: fi(:)
  R_TYPE,           intent(out) :: fo(:) 
  FLOAT, optional,  intent(in)  :: wim(:, :)

  integer :: ll, ii

  if( .not. present(wim)) then
    
    !$omp parallel do private(ii)
    do ll = 1, nri
      do ii = imin(ll) + 1, imax(ll)
        fo(ii) = sum(w(1:nn, ii)*fi(ii + ri(1:nn, ll)))
      end do
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(ii)
    do ll = 1, nri
      do ii = imin(ll) + 1, imax(ll)
        fo(ii) = sum(cmplx(w(1:nn, ii), wim(1:nn, ii))*fi(ii + ri(1:nn, ll)))
      end do
    end do
    !$omp end parallel do

  end if

end subroutine X(operate_nc)

subroutine X(operate_nc_batch)(nn, nri, w, ri, imin, imax, fi, fo, wim)
  integer,          intent(in)  :: nn
  integer,          intent(in)  :: nri
  FLOAT,            intent(in)  :: w(:, :)
  integer,          intent(in)  :: ri(:, :)
  integer,          intent(in)  :: imin(:)
  integer,          intent(in)  :: imax(:)
  type(batch_t),    intent(in)  :: fi
  type(batch_t),    intent(out) :: fo
  FLOAT, optional,  intent(in)  :: wim(:, :)

  integer :: ll, ii, idim, ist
  
  if( .not. present(wim)) then
    
    do ll = 1, nri
      forall(ist = 1:fi%nst, idim = 1:fi%dim, ii = imin(ll) + 1:imax(ll))
        fo%states(ist)%X(psi)(ii, idim) = sum(w(1:nn, ii)*fi%states(ist)%X(psi)(ii + ri(1:nn, ll), idim))
      end forall
    end do

  else

    do ll = 1, nri
      forall(ist = 1:fi%nst, idim = 1:fi%dim, ii = imin(ll) + 1:imax(ll))
        fo%states(ist)%X(psi)(ii, idim) = sum(cmplx(w(1:nn, ii), wim(1:nn, ii))*fi%states(ist)%X(psi)(ii + ri(1:nn, ll), idim))
      end forall
    end do

  end if
end subroutine X(operate_nc_batch)

subroutine X(select_op)(op, points, nri, imin, imax, ri)
  type(nl_operator_t), intent(in)    :: op
  integer,             intent(in)    :: points
  integer,             intent(out)   :: nri
  integer,             pointer       :: imin(:)
  integer,             pointer       :: imax(:)
  integer,             pointer       :: ri(:, :)

  select case(points)
  case(OP_ALL)
    nri  =  op%nri
    imin => op%rimap_inv(1:)
    imax => op%rimap_inv(2:)
    ri   => op%ri
  case(OP_INNER)
    nri  =  op%inner%nri
    imin => op%inner%imin
    imax => op%inner%imax
    ri   => op%inner%ri
  case(OP_OUTER)
    nri  =  op%outer%nri
    imin => op%outer%imin
    imax => op%outer%imax
    ri   => op%outer%ri
  case default
    ASSERT(.false.)
  end select

end subroutine X(select_op)

! ---------------------------------------------------------
subroutine X(nl_operator_operate)(op, fi, fo, ghost_update, profile, points)
  R_TYPE,              intent(inout) :: fi(:)  ! fi(op%np_part)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)  ! fo(op%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points

  real(8) :: ws(400)
  logical :: profile_
  integer :: points_, nri, nri_loc, ini, nns(1:2)
  integer, pointer :: imin(:), imax(:), ri(:, :)
#if defined(HAVE_MPI)
  logical :: update
#endif

  call push_sub('nl_operator_inc.Xnl_operator_operate')

  profile_ = .true. 
  if(present(profile)) profile_ = profile

  if(profile_) call profiling_in(nl_operate_profile, "NL_OPERATOR")

  ! update ghost points if required

#if defined(HAVE_MPI)
  update = .true.
  if(present(ghost_update)) update = ghost_update
  if(op%m%parallel_in_domains .and. update) call X(vec_ghost_update)(op%m%vp, fi)
#endif

  ! select whether we want to apply the operator to all points or just
  ! to a part of them

  points_ = OP_ALL
  if(present(points)) points_ = points

  call X(select_op)(op, points_, nri, imin, imax, ri)

  ! now apply the operator

  if(nri > 0) then
    if(op%cmplx_op) then

      if(op%const_w) then
        call X(operate)(op%stencil%size, nri, op%w_re(:, 1), ri, imin, imax, fi, fo, op%w_im(:, 1))
      else
        call X(operate_nc)(op%stencil%size, nri, op%w_re, ri, imin, imax, fi, fo, op%w_im)
      end if

    else

      if(.not. op%const_w) then
        call X(operate_nc)(op%stencil%size, nri, op%w_re, ri, imin, imax, fi, fo)
      else

        !$omp parallel private(ini, nri_loc, ws)
#ifdef USE_OMP
        call multicomm_divide_range_omp(nri, ini, nri_loc)
#else 
        ini = 1
        nri_loc = nri
#endif

        select case(op%X(function))
        case(OP_FORTRAN)
          call X(operate)(op%stencil%size, nri, op%w_re(:, 1), ri, imin, imax, fi, fo)
        case(OP_C)
          call X(operate_ri)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), fi(1), fo(1))
#ifdef HAVE_VEC
        case(OP_VEC)
          call X(operate_ri_vec)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), fi(1), fo(1))
#endif
#ifdef HAVE_AS
        case(OP_AS)
          nns(1) = op%stencil%size
          nns(2) = nri_loc
          call X(operate_as)(nns, op%w_re(1, 1), ri(1, ini), imin(ini), imax(ini), fi(1), fo(1), ws(1))
#endif
        end select
        !$omp end parallel

      end if

    end if

    if(profile_) then
      if(op%cmplx_op) then
        call profiling_count_operations((imax(nri) - imin(1))*op%stencil%size*(R_ADD + R_MUL))
      else 
        call profiling_count_operations((imax(nri) - imin(1))*op%stencil%size*2*R_ADD)
      end if
    end if

  end if

  if(profile_) call profiling_out(nl_operate_profile)

  call pop_sub()
end subroutine X(nl_operator_operate)

! ---------------------------------------------------------
subroutine X(nl_operator_operate_batch)(op, fi, fo, ghost_update, points)
  type(nl_operator_t), intent(in)    :: op
  type(batch_t),       intent(inout) :: fi
  type(batch_t),       intent(inout) :: fo
  logical, optional,   intent(in)    :: ghost_update
  integer, optional,   intent(in)    :: points

  integer :: idim, ist, points_
  logical :: update
  integer :: nri, nri_loc, ini, nns(1:2)
  integer, pointer :: imin(:), imax(:), ri(:, :)
  R_TYPE,  pointer ::  pfi(:), pfo(:)
  real(8) :: ws(400)

  call profiling_in(operate_batch_prof, "NL_OPERATOR_BATCH")
  
  
  points_ = OP_ALL
  if(present(points)) points_ = points

  call X(select_op)(op, points_, nri, imin, imax, ri)

#ifdef HAVE_MPI
  update = .true.
  if(present(ghost_update)) update = ghost_update

  if(op%m%parallel_in_domains .and. update) then
    do ist = 1, fi%nst
      do idim = 1, fi%dim
        call X(vec_ghost_update)(op%m%vp, fi%states(ist)%X(psi)(:, idim))
      end do
    end do
  end if
#endif

  if(nri > 0) then
    if(op%const_w) then
      if(op%cmplx_op) then
        do ist = 1, fi%nst
          do idim = 1, fi%dim
            pfi => fi%states(ist)%X(psi)(:, idim)
            pfo => fo%states(ist)%X(psi)(:, idim)
            call X(operate)(op%stencil%size, nri, op%w_re(:, 1), ri, imin, imax, pfi, pfo, op%w_im(:, 1))
          end do
        end do
      else
        !$omp parallel private(ini, nri_loc, ws, idim, ist, pfi, pfo)
#ifdef USE_OMP
        call multicomm_divide_range_omp(nri, ini, nri_loc)
#else 
        ini = 1
        nri_loc = nri
#endif
        do idim = 1, fi%dim
          do ist = 1, fi%nst
            pfi => fi%states(ist)%X(psi)(:, idim)
            pfo => fo%states(ist)%X(psi)(:, idim)
            select case(op%X(function))
            case(OP_FORTRAN)
              call X(operate)(op%stencil%size, nri, op%w_re(:, 1), ri, imin, imax, pfi, pfo)
            case(OP_C)
              call X(operate_ri)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1))
#ifdef HAVE_VEC
            case(OP_VEC)
              call X(operate_ri_vec)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1))
#endif
#ifdef HAVE_AS
            case(OP_AS)
              nns(1) = op%stencil%size
              nns(2) = nri_loc
              call X(operate_as)(nns, op%w_re(1, 1), ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1), ws(1))
#endif
            end select
          end do
        end do
        !$omp end parallel
      end if
    else
      if(op%cmplx_op) then
        call X(operate_nc_batch)(op%stencil%size, nri, op%w_re, ri, imin, imax, fi, fo, wim = op%w_im)
      else
        call X(operate_nc_batch)(op%stencil%size, nri, op%w_re, ri, imin, imax, fi, fo)
      end if
    end if

    if(op%cmplx_op) then
      call profiling_count_operations(fi%nst*(imax(nri) - imin(1))*op%stencil%size*(R_ADD + R_MUL))
    else
      call profiling_count_operations(fi%nst*(imax(nri) - imin(1))*op%stencil%size*2*R_ADD)
    end if
  end if

  call profiling_out(operate_batch_prof)

end subroutine X(nl_operator_operate_batch)

! ---------------------------------------------------------
subroutine X(nl_operator_operate_diag)(op, fo)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)

  call push_sub('nl_operator_inc.Xnl_operator_operate_diag')
  
  if(op%cmplx_op) then
#ifdef R_TCOMPLEX
    if(op%const_w) then
      fo(1:op%np) = cmplx(op%w_re(op%stencil%center, 1), op%w_im(op%stencil%center, 1))
    else
      fo(1:op%np) = cmplx(op%w_re(op%stencil%center, 1:op%np), op%w_im(op%stencil%center, 1:op%np))
    end if
#endif
  else
    if(op%const_w) then
      fo(1:op%np) = op%w_re(op%stencil%center, 1)
    else
      fo(1:op%np) = op%w_re(op%stencil%center, 1:op%np)
    end if
  end if
  
  call pop_sub()
  
end subroutine X(nl_operator_operate_diag)



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
