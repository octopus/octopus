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

!----------------------------------------------------
subroutine X(nl_operator_tune)(op, best)
  type(nl_operator_t), intent(inout) :: op
  FLOAT,  optional,    intent(out)   :: best

  R_TYPE, allocatable :: in(:), out(:)
  real(8) :: noperations, flops(OP_MIN:OP_MAX), itime, ftime, bwidth(OP_MIN:OP_MAX), dvolume
  integer :: method, ii, reps, iunit
  character(len=2)  :: marker
  character(len=6)  :: filenum
  character(len=30) :: filename

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
  integer :: ierr
#endif

  call push_sub('nl_operator_inc.Xnl_operator_tune')
  
  !count the total number of floating point operations  
  noperations =  op%mesh%np*op%stencil%size*M_TWO*R_OPS

  !the volume of data that has to be moved in bytes
  dvolume = (op%mesh%np_part + op%mesh%np) * R_OPS * R_SIZE + (op%nri + 1)*op%stencil%size*M_FOUR

  !measure performance of each function
  SAFE_ALLOCATE(in(1:op%mesh%np_part))
  SAFE_ALLOCATE(out(1:op%mesh%np))
  in(:) = M_ZERO
  flops = M_ZERO

  do method = OP_MIN, OP_MAX

    !skip methods that are not available
#ifdef R_TCOMPLEX
    if (op_is_available(method, TYPE_CMPLX) == 0) cycle
#else
    if (op_is_available(method, TYPE_FLOAT)  == 0) cycle
#endif
    op%X(function) = method

    ! set the number of repetitions so the total time is about 0.1
    ! secs at 1 Gigaflop (we restrict the value between 1 and 30)
    reps = min(30, max(1, nint(CNST(1e8)/noperations)))

    if(in_debug_mode) then
      write(message(1), '(6a)') 'Info: Profiling non-local operator: ', trim(op%label), ' - ', &
        trim(type), ' - ', op_function_name(method)
      call write_info(1)
    end if

    itime = loct_clock()
#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, ierr)
#endif
    do ii = 1, reps
      call X(nl_operator_operate)(op, in, out, ghost_update = .false., profile = .false.)
    end do
#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, ierr)
#endif
    ftime = loct_clock()

    flops(method) = noperations * reps / (ftime - itime) 
    bwidth(method) = dvolume * reps / (ftime - itime)
    
  end do

  SAFE_DEALLOCATE_A(in)
  SAFE_DEALLOCATE_A(out)

  ! choose the best method
  op%X(function) = OP_MIN
  do method = OP_MIN + 1, OP_MAX
    if(flops(method) > flops(op%X(function))) op%X(function) = method
  end do

  if(present(best)) then
#ifdef HAVE_MPI
    call MPI_Allreduce(flops(op%X(function)), best, 1, &
      MPI_DOUBLE_PRECISION, MPI_SUM, mpi_world%comm, ierr)
    best = best/CNST(1e6)
#else
    best = flops(op%X(function))/CNST(1e6)
#endif
  end if

  write(filenum, '(i6.6)') mpi_world%rank
  filename = 'exec/nl_operator_prof.'//trim(filenum)
  
  !print to file
  if(.not. initialized) call loct_rm(trim(filename))
  initialized = .true.
  
  iunit = io_open(trim(filename), action='write', position='append', is_tmp = .true.)
    
#ifdef R_TCOMPLEX
    write (iunit, '(3a)')   'Operator       = ', trim(op%label), " complex"
#else
    write (iunit, '(3a)')   'Operator       = ', trim(op%label), " real"
#endif
    write (iunit, '(a,i8)') 'Grid points    = ', op%mesh%np
    write (iunit, '(a,i8)') 'Stencils       = ', op%nri
    write (iunit, '(a,i8)') 'Stencil points = ', op%stencil%size

#ifdef HAVE_MPI
    if(op%mesh%parallel_in_domains) then
      write (iunit, '(a,i8)') 'Inner points   = ', sum(op%inner%imax(1:op%inner%nri) - op%inner%imin(1:op%inner%nri))
      write (iunit, '(a,i8)') 'Inner stencils = ', op%inner%nri
      write (iunit, '(a,i8)') 'Outer points   = ', sum(op%outer%imax(1:op%outer%nri) - op%outer%imin(1:op%outer%nri))
      write (iunit, '(a,i8)') 'Outer stencils = ', op%outer%nri
    end if
#endif

    do method = OP_MIN, OP_MAX
#ifdef R_TCOMPLEX
      if (op_is_available(method, TYPE_CMPLX) == 0) cycle
#else
      if (op_is_available(method, TYPE_FLOAT)  == 0) cycle
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

  call pop_sub('nl_operator_inc.Xnl_operator_tune')

end subroutine X(nl_operator_tune)


! ---------------------------------------------------------
subroutine X(nl_operator_operate_batch)(op, fi, fo, ghost_update, profile, points)
  type(nl_operator_t), intent(in)    :: op
  type(batch_t),       intent(inout) :: fi
  type(batch_t),       intent(inout) :: fo
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points

  integer :: ist, points_, cop
  logical :: ghost_update_, profile_
  integer :: nri, nri_loc, ini
  integer, pointer :: imin(:), imax(:), ri(:, :)
  R_TYPE,  pointer ::  pfi(:), pfo(:)
  real(8) :: ws(400)
#ifdef HAVE_AS
  integer nns(1:2)
#endif

  call push_sub('nl_operator_inc.Xnl_operator_operate_batch')

  ASSERT(fi%nst_linear == fo%nst_linear)

  do ist = 1, fi%nst_linear
    ASSERT(ubound(fi%states_linear(ist)%X(psi)(:), dim=1) == op%mesh%np_part)
    ASSERT(ubound(fo%states_linear(ist)%X(psi)(:), dim=1) >= op%mesh%np)
  end do

  points_ = OP_ALL
  if(present(points)) points_ = points

  profile_ = .true.
  if(present(profile)) profile_ = profile
  if(profile_) call profiling_in(operate_batch_prof, "NL_OPERATOR_BATCH")

  call select_op()

  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

#ifdef HAVE_MPI
  if(op%mesh%parallel_in_domains .and. ghost_update_) then
    do ist = 1, fi%nst_linear
      call X(vec_ghost_update)(op%mesh%vp, fi%states_linear(ist)%X(psi)(:))
    end do
  end if
#endif

  if(nri > 0) then
    if(.not.op%const_w) then
      call operate_non_const_weights()
    else if(op%cmplx_op .or. op%X(function)==OP_FORTRAN) then
      call operate_const_weights()
#ifdef HAVE_OPENCL
    else if(opencl_is_enabled() .and. batch_is_in_buffer(fi) .and. batch_is_in_buffer(fo)) then
      call operate_opencl()
#endif
    else
      !$omp parallel private(ini, nri_loc, ws, ist, pfi, pfo)
#ifdef USE_OMP
      call multicomm_divide_range_omp(nri, ini, nri_loc)
#else 
      ini = 1
      nri_loc = nri
#endif
      do ist = 1, fi%nst_linear
        pfi => fi%states_linear(ist)%X(psi)(:)
        pfo => fo%states_linear(ist)%X(psi)(:)
        select case(op%X(function))
        case(OP_C)
          call X(operate_ri)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1))
#ifdef HAVE_VEC
        case(OP_VEC)
          call X(operate_ri_vec)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1))
#endif
#if defined(HAVE_BLUE_GENE) && defined(R_TCOMPLEX)
        case(OP_BG)
          call X(operate_bg)(op%stencil%size, op%w_re(1, 1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1))
#endif
#ifdef HAVE_AS
        case(OP_AS)
          nns(1) = op%stencil%size
          nns(2) = nri_loc
          call X(operate_as)(nns, op%w_re(1, 1), ri(1, ini), imin(ini), imax(ini), pfi(1), pfo(1), ws(1))
#endif
        end select

      end do
      !$omp end parallel
    end if

    ! count operations
    if(profile_) then
      if(op%cmplx_op) then
        cop = fi%nst_linear*(imax(nri) - imin(1))*op%stencil%size*(R_ADD + R_MUL)
      else
        cop = fi%nst_linear*(imax(nri) - imin(1))*op%stencil%size*2*R_ADD
      end if
      call profiling_count_operations(cop)
    end if
  end if

  if(profile_) call profiling_out(operate_batch_prof)
  call pop_sub('nl_operator_inc.Xnl_operator_operate_batch')

contains

  ! ---------------------------------------------------------
  subroutine select_op()
    select case(points_)
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
  end subroutine select_op


  ! ---------------------------------------------------------
  subroutine operate_const_weights()
    integer :: nn, ll, ii, ist

    ASSERT(.not. (batch_is_in_buffer(fi) .or. batch_is_in_buffer(fo)))

    nn = op%stencil%size

    if(op%cmplx_op) then
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = sum(cmplx(op%w_re(1:nn, 1), op%w_im(1:nn, 1)) * &
            fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = sum(op%w_re(1:nn, 1)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    end if
  end subroutine operate_const_weights


  ! ---------------------------------------------------------
  subroutine operate_non_const_weights()
    integer :: nn, ll, ii, ist

    ASSERT(.not. (batch_is_in_buffer(fi) .or. batch_is_in_buffer(fo)))

    if(op%cmplx_op) then
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = sum(cmplx(op%w_re(1:nn, ii), op%w_im(1:nn, ii)) * &
            fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = sum(op%w_re(1:nn, ii)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    end if

  end subroutine operate_non_const_weights

#ifdef HAVE_OPENCL

  ! ------------------------------------------
  subroutine operate_opencl()
    type(opencl_mem_t) :: buff_weights
    integer :: pnri, bsize

    ASSERT(.not. op%mesh%parallel_in_domains)

    call opencl_create_buffer(buff_weights, CL_MEM_READ_ONLY, TYPE_FLOAT, op%stencil%size)
    call opencl_write_buffer(buff_weights, op%stencil%size, op%w_re(:, 1))

    call opencl_set_kernel_arg(X(operate),  0, fi%nst_linear)
    call opencl_set_kernel_arg(X(operate),  1, op%stencil%size)
    call opencl_set_kernel_arg(X(operate),  2, nri)
    call opencl_set_kernel_arg(X(operate),  3, op%buff_ri)
    call opencl_set_kernel_arg(X(operate),  4, op%buff_imin)
    call opencl_set_kernel_arg(X(operate),  5, op%buff_imax)
    call opencl_set_kernel_arg(X(operate),  6, buff_weights)
    call opencl_set_kernel_arg(X(operate),  7, fi%buffer)
    call opencl_set_kernel_arg(X(operate),  8, batch_buffer_ubound(fi))
    call opencl_set_kernel_arg(X(operate),  9, fo%buffer)
    call opencl_set_kernel_arg(X(operate), 10, batch_buffer_ubound(fo))

    pnri = nri
    bsize = opencl_max_workgroup_size()
    if(mod(pnri, bsize) /= 0) pnri = pnri + bsize - mod(pnri, bsize)

    call opencl_kernel_run(X(operate), (/pnri/), (/bsize/))

    call batch_buffer_was_modified(fo)

    call opencl_release_buffer(buff_weights)
  end subroutine operate_opencl

#endif

end subroutine X(nl_operator_operate_batch)


! ---------------------------------------------------------
subroutine X(nl_operator_operate)(op, fi, fo, ghost_update, profile, points)
  R_TYPE,              intent(inout) :: fi(:)  ! fi(op%np_part)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)  ! fo(op%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points

  type(batch_t) :: batch_fi, batch_fo

  call push_sub('nl_operator_inc.Xnl_operator_operate')

  call batch_init     (batch_fi, 1)
  call batch_add_state(batch_fi, fi)

  call batch_init     (batch_fo, 1)
  call batch_add_state(batch_fo, fo)

  call X(nl_operator_operate_batch)(op, batch_fi, batch_fo, ghost_update, profile, points)

  call batch_end(batch_fi)
  call batch_end(batch_fo)

  call pop_sub('nl_operator_inc.Xnl_operator_operate')
end subroutine X(nl_operator_operate)


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
  
  call pop_sub('nl_operator_inc.Xnl_operator_operate_diag')
  
end subroutine X(nl_operator_operate_diag)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
