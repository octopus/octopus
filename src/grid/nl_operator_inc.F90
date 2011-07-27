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

subroutine X(nl_operator_operate_batch)(op, fi, fo, ghost_update, profile, points, factor)
  type(nl_operator_t), intent(in)    :: op
  type(batch_t),       intent(inout) :: fi
  type(batch_t),       intent(inout) :: fo
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: profile
  integer, optional,   intent(in)    :: points
  FLOAT,   optional,   intent(in)    :: factor

  integer :: ist, points_
  real(8) :: cop
  logical :: ghost_update_, profile_
  integer :: nri, nri_loc, ini
  integer, pointer :: imin(:), imax(:), ri(:, :)
  R_TYPE,  pointer ::  pfi(:), pfo(:)
  FLOAT,   pointer :: wre(:), wim(:)

  PUSH_SUB(X(nl_operator_operate_batch))

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

  nullify(wre)
  nullify(wim)

  if(op%const_w) then
    if(present(factor)) then
      SAFE_ALLOCATE(wre(1:op%stencil%size))
      wre = op%w_re(:, 1)*factor
      if(op%cmplx_op) then
        SAFE_ALLOCATE(wim(1:op%stencil%size))
        wim = op%w_im(:, 1)*factor
      end if
    else
      wre => op%w_re(:, 1)
      if(op%cmplx_op) wim => op%w_im(:, 1)
    end if
  end if

  if(nri > 0) then
    if(.not.op%const_w) then
      call operate_non_const_weights()
    else if(op%cmplx_op .or. X(function_global) == OP_FORTRAN) then
      call operate_const_weights()
#ifdef HAVE_OPENCL
    else if(opencl_is_enabled() .and. batch_is_packed(fi) .and. batch_is_packed(fo)) then
      call operate_opencl()
#endif
    else 

      !$omp parallel private(ini, nri_loc, ist, pfi, pfo)
#ifdef HAVE_OPENMP
      call multicomm_divide_range_omp(nri, ini, nri_loc)
#else 
      ini = 1
      nri_loc = nri
#endif
      
      if(batch_is_packed(fi) .and. batch_is_packed(fo)) then
        call operate_ri_vec(op%stencil%size, wre(1), nri_loc, ri(1, ini), imin(ini), imax(ini), &
          fi%pack%X(psi)(1, 1), log2(fi%pack%size_real(1)), fo%pack%X(psi)(1, 1))
      else
        do ist = 1, fi%nst_linear
          pfi => fi%states_linear(ist)%X(psi)(:)
          pfo => fo%states_linear(ist)%X(psi)(:)
          
#ifdef R_TREAL
#define LOGLDF 0
#else
#define LOGLDF 1
#endif
          call operate_ri_vec(op%stencil%size, wre(1), nri_loc, ri(1, ini), imin(ini), imax(ini), pfi(1), LOGLDF, pfo(1))
#undef LOGLDF
        end do
      end if
      !$omp end parallel
    end if

    ! count operations
    if(profile_) then
      if(op%cmplx_op) then
        cop = fi%nst_linear*dble(imax(nri) - imin(1))*op%stencil%size*(R_ADD + R_MUL)
      else
        cop = fi%nst_linear*dble(imax(nri) - imin(1))*op%stencil%size*2*R_ADD
      end if
      call profiling_count_operations(cop)
    end if

    call batch_pack_was_modified(fo)

  end if

  if(op%const_w .and. present(factor)) then
    SAFE_DEALLOCATE_P(wre)
    SAFE_DEALLOCATE_P(wim)
  end if

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


  ! ---------------------------------------------------------
  subroutine operate_const_weights()
    integer :: nn, ll, ii, ist

    PUSH_SUB(X(nl_operator_operate_batch).operate_const_weights)

    nn = op%stencil%size

    if(op%cmplx_op) then
      ASSERT(.not. (batch_is_packed(fi) .or. batch_is_packed(fo)))

      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = sum(cmplx(wre(1:nn), wim(1:nn)) * &
            fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do

    else

      if(batch_is_packed(fi)) then
        !$omp parallel do private(ll, ist, ii)
        do ll = 1, nri
          forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
            fo%pack%X(psi)(ist, ii) = sum(wre(1:nn)*fi%pack%X(psi)(ist, ii + ri(1:nn, ll)))
          end forall
        end do
        !$omp end parallel do
      else   
        !$omp parallel do private(ll, ist, ii)
        do ll = 1, nri
          forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
            fo%states_linear(ist)%X(psi)(ii) = sum(wre(1:nn)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
          end forall
        end do
        !$omp end parallel do
      end if

    end if

    POP_SUB(X(nl_operator_operate_batch).operate_const_weights)
  end subroutine operate_const_weights


  ! ---------------------------------------------------------
  subroutine operate_non_const_weights()
    integer :: nn, ll, ii, ist
    FLOAT :: factor_

    PUSH_SUB(X(nl_operator_operate_batch).operate_non_const_weights)

    ASSERT(.not. (batch_is_packed(fi) .or. batch_is_packed(fo)))

    factor_ = M_ONE
    if(present(factor)) factor_ = factor

    if(op%cmplx_op) then
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = factor_*sum(cmplx(op%w_re(1:nn, ii), op%w_im(1:nn, ii)) * &
            fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(ll, ist, ii)
      do ll = 1, nri
        nn = op%nn(ll)
        forall(ist = 1:fi%nst_linear, ii = imin(ll) + 1:imax(ll))
          fo%states_linear(ist)%X(psi)(ii) = factor_*sum(op%w_re(1:nn, ii)*fi%states_linear(ist)%X(psi)(ii + ri(1:nn, ll)))
        end forall
      end do
      !$omp end parallel do
    end if

    POP_SUB(X(nl_operator_operate_batch).operate_non_const_weights)
  end subroutine operate_non_const_weights

#ifdef HAVE_OPENCL

  ! ------------------------------------------
  subroutine operate_opencl()
    integer :: pnri, bsize, isize, ist, local_mem_size
    type(opencl_mem_t) :: buff_weights
    type(profile_t), save :: prof

    PUSH_SUB(X(nl_operator_operate_batch).operate_opencl)
    call profiling_in(prof, "CL_NL_OPERATOR")
    
    ASSERT(points_ == OP_ALL)

    call opencl_create_buffer(buff_weights, CL_MEM_READ_ONLY, TYPE_FLOAT, op%stencil%size)
    call opencl_write_buffer(buff_weights, op%stencil%size, wre)

    select case(function_opencl)
    case(OP_INVMAP)
      call opencl_set_kernel_arg(kernel_operate, 0, op%stencil%size)
      call opencl_set_kernel_arg(kernel_operate, 1, nri)
      call opencl_set_kernel_arg(kernel_operate, 2, op%buff_ri)
      call opencl_set_kernel_arg(kernel_operate, 3, op%buff_imin)
      call opencl_set_kernel_arg(kernel_operate, 4, op%buff_imax)
      call opencl_set_kernel_arg(kernel_operate, 5, buff_weights)
      call opencl_set_kernel_arg(kernel_operate, 6, fi%pack%buffer)
      call opencl_set_kernel_arg(kernel_operate, 7, fi%pack%size_real(1))
      call opencl_set_kernel_arg(kernel_operate, 8, fo%pack%buffer)
      call opencl_set_kernel_arg(kernel_operate, 9, fo%pack%size_real(1))
      
      bsize = 128
      isize = 8
      pnri = pad(nri, bsize)
      
      call opencl_kernel_run(kernel_operate, &
        (/fi%pack%size_real(1), pnri/), (/fi%pack%size_real(1), bsize/(fi%pack%size_real(1))/))
      
    case(OP_MAP)
      call opencl_set_kernel_arg(kernel_operate, 0, op%stencil%size)
      call opencl_set_kernel_arg(kernel_operate, 1, op%mesh%np)
      call opencl_set_kernel_arg(kernel_operate, 2, op%buff_ri)
      call opencl_set_kernel_arg(kernel_operate, 3, op%buff_map)
      call opencl_set_kernel_arg(kernel_operate, 4, buff_weights)
      call opencl_set_kernel_arg(kernel_operate, 5, fi%pack%buffer)
      call opencl_set_kernel_arg(kernel_operate, 6, log2(fi%pack%size_real(1)))
      call opencl_set_kernel_arg(kernel_operate, 7, fo%pack%buffer)
      call opencl_set_kernel_arg(kernel_operate, 8, log2(fi%pack%size_real(1)))

      
      local_mem_size = f90_cl_device_local_mem_size(opencl%device)
      isize = int(dble(local_mem_size)/(op%stencil%size*types_get_size(TYPE_INTEGER)))
      isize = isize - mod(isize, fi%pack%size_real(1))
      bsize = fi%pack%size_real(1)*isize
      bsize = min(opencl_kernel_workgroup_size(kernel_operate), bsize)

      if(bsize < fi%pack%size_real(1)) then
        message(1) = "The value of StatesBlockSize is too large for this OpenCL implementation."
        call messages_fatal(1)
      end if

      isize = bsize/(fi%pack%size_real(1))

      ASSERT(isize > 0)
      ASSERT(isize*op%stencil%size*types_get_size(TYPE_INTEGER) <= local_mem_size)

      call opencl_set_kernel_arg(kernel_operate, 9, TYPE_INTEGER, isize*op%stencil%size)

      call opencl_kernel_run(kernel_operate, &
        (/fi%pack%size_real(1), pad(op%mesh%np, bsize)/), (/fi%pack%size_real(1), isize/))
      
      call profiling_count_transfers(op%stencil%size*op%mesh%np + op%mesh%np, isize)
      do ist = 1, fi%nst_linear
        call profiling_count_transfers(op%mesh%np_part*op%stencil%size + op%mesh%np, R_TOTYPE(M_ONE))
      end do
    end select
    
    call opencl_finish()

    call opencl_release_buffer(buff_weights)

    call profiling_out(prof)
    POP_SUB(X(nl_operator_operate_batch).operate_opencl)
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
  
  POP_SUB(X(nl_operator_operate_diag))
  
end subroutine X(nl_operator_operate_diag)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
