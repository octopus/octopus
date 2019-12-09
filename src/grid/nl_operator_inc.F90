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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------

subroutine X(nl_operator_operate_batch)(op, fi, fo, ghost_update, profile, points, factor)
  type(nl_operator_t), target, intent(in)    :: op
  type(batch_t),       target, intent(inout) :: fi
  type(batch_t),               intent(inout) :: fo !< this should be target, but old ifort 9.1 segfaults with it
  logical,         optional,   intent(in)    :: ghost_update
  logical,         optional,   intent(in)    :: profile
  integer,         optional,   intent(in)    :: points
  FLOAT,           optional,   intent(in)    :: factor

  integer :: ist, points_
  real(8) :: cop
  logical :: ghost_update_, profile_, use_opencl
  integer :: nri
  integer, pointer :: imin(:), imax(:), ri(:, :)
!  FLOAT, allocatable :: wre(:), wim(:)
  R_BASE, allocatable :: wre(:), wim(:)  
#ifdef R_TREAL
  integer, parameter :: logldf = 0
#else
  integer, parameter :: logldf = 1
#endif
#ifndef SINGLE_PRECISION
  integer :: nri_loc, ini
  R_TYPE,  pointer :: pfi(:), pfo(:)
#endif
  
  PUSH_SUB(X(nl_operator_operate_batch))

  ASSERT(batch_status(fi) == batch_status(fo))
  ASSERT(batch_type(fi) == R_TYPE_VAL)
  ASSERT(batch_type(fo) == R_TYPE_VAL)

  ASSERT(fi%nst_linear == fo%nst_linear)

  points_ = OP_ALL
  if(present(points)) points_ = points

  profile_ = .true.
  if(present(profile)) profile_ = profile
  if(profile_) call profiling_in(operate_batch_prof, "NL_OPERATOR_BATCH")

  call select_op()

  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

  if(op%mesh%parallel_in_domains .and. ghost_update_) then
    ASSERT(.not. batch_is_packed(fi))
    
    do ist = 1, fi%nst_linear
#ifdef HAVE_MPI
      call X(vec_ghost_update)(op%mesh%vp, fi%states_linear(ist)%X(psi)(:))
#endif
    end do
  end if

  if(op%const_w) then
    SAFE_ALLOCATE(wre(1:op%stencil%size))

    wre(1:op%stencil%size) = op%w(:, 1)

    if(present(factor)) then
      wre(1:op%stencil%size) = wre(1:op%stencil%size)*factor
    end if
  end if

  use_opencl = .false.
  
  if(nri > 0) then
    if(.not.op%const_w) then
      call operate_non_const_weights()
    else if(accel_is_enabled() .and. batch_is_packed(fi) .and. batch_is_packed(fo)) then
      use_opencl = .true.
      call operate_opencl()
    else if(X(function_global) == OP_FORTRAN) then
      call operate_const_weights()
    else
      
! for the moment this is not implemented
#ifndef SINGLE_PRECISION

      !$omp parallel private(ini, nri_loc, ist, pfi, pfo)
#ifdef HAVE_OPENMP
      call multicomm_divide_range_omp(nri, ini, nri_loc)
#else 
      ini = 1
      nri_loc = nri
#endif
      
      if(batch_is_packed(fi) .and. batch_is_packed(fo)) then
        
        ASSERT(ubound(fi%pack%X(psi), dim = 2) == op%mesh%np_part)
        ASSERT(ubound(fo%pack%X(psi), dim = 2) >= op%mesh%np)
        
        call X(operate_ri_vec)(op%stencil%size, wre(1), nri_loc, ri(1, ini), imin(ini), imax(ini), &
          fi%pack%X(psi)(1, 1), log2(fi%pack%size_real(1)), fo%pack%X(psi)(1, 1))
      else
        do ist = 1, fi%nst_linear

          ASSERT(ubound(fi%states_linear(ist)%X(psi), dim = 1) == op%mesh%np_part)
          ASSERT(ubound(fo%states_linear(ist)%X(psi), dim = 1) >= op%mesh%np)
          
          pfi => fi%states_linear(ist)%X(psi)(:)
          pfo => fo%states_linear(ist)%X(psi)(:)
          
          call X(operate_ri_vec)(op%stencil%size, wre(1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), logldf, pfo(1))
        end do
      end if
      !$omp end parallel
#endif
    end if

    ! count operations
    if(profile_ .and. .not. use_opencl) then
      cop = fi%nst_linear*dble(imax(nri) - imin(1))*op%stencil%size*2*R_ADD
      call profiling_count_operations(cop)
    end if
  end if

  SAFE_DEALLOCATE_A(wre)
  SAFE_DEALLOCATE_A(wim)

  if(profile_) call profiling_out(operate_batch_prof)
  POP_SUB(X(nl_operator_operate_batch))

contains

  ! ---------------------------------------------------------
  subroutine select_op()

    PUSH_SUB(X(nl_operator_operate_batch).select_op)

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

    POP_SUB(X(nl_operator_operate_batch).select_op)
  end subroutine select_op


!pgi$r novector
!This is a pragma for the PGI compiler, preventing vector optimization for this subroutine

  ! ---------------------------------------------------------
  subroutine operate_const_weights()
    integer :: nn, ll, ii, ist

    PUSH_SUB(X(nl_operator_operate_batch).operate_const_weights)

    nn = op%stencil%size

    select case(batch_status(fi))

    case(BATCH_DEVICE_PACKED)

      ASSERT(.false.)
      
    case(BATCH_NOT_PACKED)
      
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        do ii = imin(ll) + 1, imax(ll)
          do ist = 1, fi%nst_linear
            fo%states_linear(ist)%X(psi)(ii) = sum(wre(1:nn)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
          end do
        end do
      end do
      !$omp end parallel do

    case(BATCH_PACKED)

      ASSERT(allocated(fo%pack%X(psi)))
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        do ii = imin(ll) + 1, imax(ll)
          do ist = 1, fi%nst_linear
            fo%pack%X(psi)(ist, ii) = sum(wre(1:nn)*fi%pack%X(psi)(ist, ii + ri(1:nn, ll)))
          end do
        end do
      end do
      !$omp end parallel do

    end select

    POP_SUB(X(nl_operator_operate_batch).operate_const_weights)
  end subroutine operate_const_weights


  ! ---------------------------------------------------------
  subroutine operate_non_const_weights()
    integer :: nn, ll, ii, ist
    FLOAT :: factor_

    PUSH_SUB(X(nl_operator_operate_batch).operate_non_const_weights)

    factor_ = M_ONE
    if(present(factor)) factor_ = factor

    select case(batch_status(fi))

    case(BATCH_DEVICE_PACKED)

      ASSERT(.false.)
      
    case(BATCH_NOT_PACKED)
      
      !$omp parallel do private(ll, ist, ii, nn)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = factor_*sum(op%w(1:nn, ii)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do

    case(BATCH_PACKED)

      !$omp parallel do private(ll, ist, ii, nn)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%pack%X(psi)(ist, ii) = factor_*sum(op%w(1:nn, ii)*fi%pack%X(psi)(ist, ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
      
    end select
      
    POP_SUB(X(nl_operator_operate_batch).operate_non_const_weights)
  end subroutine operate_non_const_weights

  ! ------------------------------------------
  subroutine operate_opencl()
    integer    :: pnri, bsize, isize, localsize, ist, eff_size, iarg, npoints, dim2, dim3
    integer(8) :: local_mem_size
    type(accel_mem_t) :: buff_weights
    type(profile_t), save :: prof
    type(accel_kernel_t) :: kernel_operate

    PUSH_SUB(X(nl_operator_operate_batch).operate_opencl)
    call profiling_in(prof, "CL_NL_OPERATOR")

    ASSERT(accel_buffer_is_allocated(fi%pack%buffer))
    ASSERT(accel_buffer_is_allocated(fo%pack%buffer))
    
    kernel_operate = op%kernel

    call accel_create_buffer(buff_weights, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, op%stencil%size)

    call accel_write_buffer(buff_weights, op%stencil%size, wre)

    ASSERT(fi%pack%size_real(1) == fo%pack%size_real(1))

    eff_size = fi%pack%size_real(1)

    select case(function_opencl)
    case(OP_INVMAP)
      ASSERT(points_ == OP_ALL)
      ASSERT(accel_buffer_is_allocated(op%buff_ri))
      ASSERT(accel_buffer_is_allocated(op%buff_imin))
      ASSERT(accel_buffer_is_allocated(op%buff_imax))
      
      call accel_set_kernel_arg(kernel_operate, 0, op%stencil%size)
      call accel_set_kernel_arg(kernel_operate, 1, nri)
      call accel_set_kernel_arg(kernel_operate, 2, op%buff_ri)
      call accel_set_kernel_arg(kernel_operate, 3, op%buff_imin)
      call accel_set_kernel_arg(kernel_operate, 4, op%buff_imax)
      call accel_set_kernel_arg(kernel_operate, 5, buff_weights)
      call accel_set_kernel_arg(kernel_operate, 6, fi%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 7, log2(eff_size))
      call accel_set_kernel_arg(kernel_operate, 8, fo%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 9, log2(eff_size))

      bsize = accel_kernel_workgroup_size(kernel_operate)
      pnri = pad(nri, bsize)

      call accel_kernel_run(kernel_operate, (/eff_size, pnri/), (/eff_size, bsize/eff_size/))

    case(OP_MAP)
      ASSERT(accel_buffer_is_allocated(op%buff_ri))
      ASSERT(accel_buffer_is_allocated(op%buff_map))

      call accel_set_kernel_arg(kernel_operate, 0, op%mesh%np)
      call accel_set_kernel_arg(kernel_operate, 1, op%buff_ri)
      call accel_set_kernel_arg(kernel_operate, 2, op%buff_map)
      call accel_set_kernel_arg(kernel_operate, 3, buff_weights)
      call accel_set_kernel_arg(kernel_operate, 4, fi%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 5, log2(eff_size))
      call accel_set_kernel_arg(kernel_operate, 6, fo%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 7, log2(eff_size))

      iarg = 7

      npoints = op%mesh%np
      if(op%mesh%parallel_in_domains) then
        iarg = iarg + 1
        select case(points_)
        case(OP_INNER)
          npoints = op%ninner
          call accel_set_kernel_arg(kernel_operate, 0, op%ninner)
          call accel_set_kernel_arg(kernel_operate, iarg, op%buff_inner)
        case(OP_OUTER)
          npoints = op%nouter
          call accel_set_kernel_arg(kernel_operate, 0, op%nouter)
          call accel_set_kernel_arg(kernel_operate, iarg, op%buff_outer)
        case(OP_ALL)
          call accel_set_kernel_arg(kernel_operate, iarg, op%buff_all)
        case default
          ASSERT(.false.)
        end select
      end if
      
      if(accel_use_shared_mem()) then
        local_mem_size = accel_local_memory_size()
        localsize = int(dble(local_mem_size)/(op%stencil%size*types_get_size(TYPE_INTEGER)))
        localsize = localsize - mod(localsize, eff_size)
        bsize = eff_size*localsize
        bsize = min(accel_kernel_workgroup_size(kernel_operate), bsize)
      else
        bsize = accel_kernel_workgroup_size(kernel_operate)
      end if
      
      if(bsize < fi%pack%size_real(1)) then
        message(1) = "The value of StatesBlockSize is too large for this OpenCL implementation."
        call messages_fatal(1)
      end if

      localsize = bsize/eff_size

      ASSERT(localsize > 0)

      if(accel_use_shared_mem()) then
        ASSERT(localsize*op%stencil%size*types_get_size(TYPE_INTEGER) <= local_mem_size)

        iarg = iarg + 1
        call accel_set_kernel_arg(kernel_operate, iarg, TYPE_INTEGER, localsize*op%stencil%size)
      end if

      dim3 = op%mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(op%mesh%np, localsize))

      call accel_kernel_run(kernel_operate, (/eff_size, dim2, dim3/), (/eff_size, localsize, 1/))

      call profiling_count_transfers(npoints*(op%stencil%size + 2), localsize)
      call profiling_count_transfers(fi%nst_linear*npoints*(op%stencil%size + 1), R_TOTYPE(M_ONE))
      
    case(OP_NOMAP)
      ASSERT(points_ == OP_ALL)
      ASSERT(accel_buffer_is_allocated(op%buff_stencil))
      ASSERT(accel_buffer_is_allocated(op%buff_xyz_to_ip))
      ASSERT(accel_buffer_is_allocated(op%buff_ip_to_xyz))
      
      call accel_set_kernel_arg(kernel_operate, 0, op%mesh%np)
      call accel_set_kernel_arg(kernel_operate, 1, op%buff_stencil)
      call accel_set_kernel_arg(kernel_operate, 2, op%buff_xyz_to_ip)
      call accel_set_kernel_arg(kernel_operate, 3, op%buff_ip_to_xyz)
      call accel_set_kernel_arg(kernel_operate, 4, buff_weights)
      call accel_set_kernel_arg(kernel_operate, 5, fi%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 6, fo%pack%buffer)
      call accel_set_kernel_arg(kernel_operate, 7, log2(eff_size))

      if(accel_use_shared_mem()) then
        local_mem_size = accel_local_memory_size()
        isize = int(dble(local_mem_size)/(op%stencil%size*types_get_size(TYPE_INTEGER)))
        isize = isize - mod(isize, eff_size)
        bsize = eff_size*isize
        bsize = min(accel_kernel_workgroup_size(kernel_operate), bsize)
      else
        bsize = accel_kernel_workgroup_size(kernel_operate)
      end if

      if(bsize < fi%pack%size_real(1)) then
        call messages_write('The value of StatesBlockSize is too large for this OpenCL implementation.')
        call messages_fatal()
      end if

      isize = bsize/eff_size

      ASSERT(isize > 0)

      if(accel_use_shared_mem()) then
        ASSERT(isize*op%stencil%size*types_get_size(TYPE_INTEGER) <= local_mem_size)
        call accel_set_kernel_arg(kernel_operate, 8, TYPE_INTEGER, isize*op%stencil%size)
      end if
      
      call accel_kernel_run(kernel_operate, (/eff_size, pad(op%mesh%np, bsize)/), (/eff_size, isize/))

      call profiling_count_transfers(op%stencil%size*op%mesh%np + op%mesh%np, isize)

      do ist = 1, fi%nst_linear
        call profiling_count_transfers(op%mesh%np_part*op%stencil%size + op%mesh%np, R_TOTYPE(M_ONE))
      end do
    end select
    
    if(profile_) then
      select case(points_)
      case(OP_INNER)
        call profiling_count_operations(fi%nst_linear*dble(op%ninner)*op%stencil%size*2*R_ADD)
      case(OP_OUTER)
        call profiling_count_operations(fi%nst_linear*dble(op%nouter)*op%stencil%size*2*R_ADD)
      case(OP_ALL)
        call profiling_count_operations(fi%nst_linear*dble(op%mesh%np)*op%stencil%size*2*R_ADD)
      case default
        ASSERT(.false.)
      end select
    end if
    
    ! when doing the inner points the synchronization is done
    ! in X(ghost_update_batch_finish) after a Waitall call to
    ! overlap communication and computation
    if(points_ /= OP_INNER) then
      call accel_finish()
    end if

    call accel_release_buffer(buff_weights)

    call profiling_out(prof)
    POP_SUB(X(nl_operator_operate_batch).operate_opencl)
  end subroutine operate_opencl

end subroutine X(nl_operator_operate_batch)

! ---------------------------------------------------------

subroutine X(nl_operator_operate)(op, fi, fo, ghost_update, profile, points)
  R_TYPE,              intent(inout) :: fi(:)  !< fi(op%np_part)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,  target,     intent(out)   :: fo(:)  !< fo(op%np). target for batch_add_state
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points

  type(batch_t) :: batch_fi, batch_fo

  PUSH_SUB(X(nl_operator_operate))

  call batch_init     (batch_fi, 1)
  call batch_add_state(batch_fi, fi)

  call batch_init     (batch_fo, 1)
  call batch_add_state(batch_fo, fo)

  call X(nl_operator_operate_batch)(op, batch_fi, batch_fo, ghost_update, profile, points)

  call batch_end(batch_fi)
  call batch_end(batch_fo)

  POP_SUB(X(nl_operator_operate))
end subroutine X(nl_operator_operate)


! ---------------------------------------------------------
subroutine X(nl_operator_operate_diag)(op, fo)
  type(nl_operator_t), intent(in)    :: op
  R_TYPE,              intent(out)   :: fo(:)

  PUSH_SUB(X(nl_operator_operate_diag))
  
  if(op%const_w) then
    fo(1:op%np) = op%w(op%stencil%center, 1)
  else
    fo(1:op%np) = op%w(op%stencil%center, 1:op%np)
  end if
  
  POP_SUB(X(nl_operator_operate_diag))
  
end subroutine X(nl_operator_operate_diag)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
