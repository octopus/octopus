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

! This module calculates the derivatives (gradients, Laplacians, etc.) 
! of a function. Note that the function whose derivative is to be calculated
! *has* to be defined (1:mesh%np_part), while the (1:mesh%np) values of the derivative
! are calculated. This was made to simplify the parallel mode, and has to be
! followed by the rest of the code.


! ---------------------------------------------------------
! Set all boundary points in ff to zero to implement zero
! boundary conditions for the derivatives, in finite system;
! or set according to periodic boundary conditions.
subroutine X(derivatives_batch_set_bc)(der, ffb)
  type(derivatives_t), intent(in)    :: der
  type(batch_t),       intent(inout) :: ffb

  integer :: pp, bndry_start, bndry_end

  PUSH_SUB(X(derivatives_batch_set_bc))
  call profiling_in(set_bc_prof, 'SET_BC')
  
  pp = der%mesh%vp%partno

  ! The boundary points are at different locations depending on the presence
  ! of ghost points due to domain parallelization.
  if(der%mesh%parallel_in_domains) then
    bndry_start = der%mesh%vp%np_local(pp) + der%mesh%vp%np_ghost(pp) + 1
    bndry_end   = der%mesh%vp%np_local(pp) + der%mesh%vp%np_ghost(pp) + der%mesh%vp%np_bndry(pp)
  else
    bndry_start = der%mesh%np + 1
    bndry_end   = der%mesh%np_part
  end if

  if(der%zero_bc)         call zero_boundaries()

  if(der%mesh%sb%mr_flag) call multiresolution()
  if(der%periodic_bc)     call periodic()
  

  call profiling_out(set_bc_prof)
  POP_SUB(X(derivatives_batch_set_bc))

contains

  ! ---------------------------------------------------------
  subroutine zero_boundaries()
    integer :: ist, ip
#ifdef HAVE_OPENCL
    integer :: localsize, globalsize
#endif

    PUSH_SUB(X(derivatives_batch_set_bc).zero_boundaries)

    select case(batch_status(ffb))
#ifdef HAVE_OPENCL
    case(BATCH_CL_PACKED)
      call opencl_set_kernel_arg(set_zero_part, 0, bndry_start - 1)
      call opencl_set_kernel_arg(set_zero_part, 1, bndry_end - 1)
      call opencl_set_kernel_arg(set_zero_part, 2, ffb%pack%buffer)
      call opencl_set_kernel_arg(set_zero_part, 3, log2(ffb%pack%size_real(1)))

      localsize = opencl_max_workgroup_size()/ffb%pack%size_real(1)
      globalsize = pad(bndry_end - bndry_start + 1, localsize)
      
      call opencl_kernel_run(set_zero_part, (/ffb%pack%size_real(1), globalsize/), (/ffb%pack%size_real(1), localsize/))
      call opencl_finish()

#endif
    case(BATCH_PACKED)
      forall(ip = bndry_start:bndry_end) 
        forall(ist = 1:ffb%nst_linear)
          ffb%pack%X(psi)(ist, ip) = R_TOTYPE(M_ZERO)
        end forall
      end forall

    case(BATCH_NOT_PACKED)
      do ist = 1, ffb%nst_linear
        forall (ip = bndry_start:bndry_end) ffb%states_linear(ist)%X(psi)(ip) = R_TOTYPE(M_ZERO)
      end do

    end select

    call batch_pack_was_modified(ffb)

    POP_SUB(X(derivatives_batch_set_bc).zero_boundaries)
  end subroutine zero_boundaries


  ! ---------------------------------------------------------
  subroutine multiresolution()
    integer :: ist, ip
    integer :: ii, jj, kk, ix, iy, iz, dx, dy, dz, i_lev
    FLOAT :: weight
    R_TYPE, pointer :: ff(:)

    PUSH_SUB(X(derivatives_batch_set_bc).multiresolution)

    do ist = 1, ffb%nst_linear
      ff => ffb%states_linear(ist)%X(psi)

      do ip = bndry_start, bndry_end
        ix = der%mesh%idx%lxyz(ip, 1)
        iy = der%mesh%idx%lxyz(ip, 2)
        iz = der%mesh%idx%lxyz(ip, 3)

        i_lev = der%mesh%resolution(ix,iy,iz)

        ! resolution is 2**num_radii for outer boundary points, but now we want inner boundary points
        if(i_lev .ne. 2**der%mesh%sb%hr_area%num_radii) then
          dx = abs(mod(ix, 2**(i_lev)))
          dy = abs(mod(iy, 2**(i_lev)))
          dz = abs(mod(iz, 2**(i_lev)))

          do ii = 1, der%mesh%sb%hr_area%interp%nn
            do jj = 1, der%mesh%sb%hr_area%interp%nn
              do kk = 1, der%mesh%sb%hr_area%interp%nn
                weight = der%mesh%sb%hr_area%interp%ww(ii) * &
                  der%mesh%sb%hr_area%interp%ww(jj) *        &
                  der%mesh%sb%hr_area%interp%ww(kk)

                ff(ip) = ff(ip) + weight * ff(der%mesh%idx%lxyz_inv(   &
                  ix + der%mesh%sb%hr_area%interp%posi(ii) * dx,       &
                  iy + der%mesh%sb%hr_area%interp%posi(jj) * dy,       &
                  iz + der%mesh%sb%hr_area%interp%posi(kk) * dz))
              end do
            end do
          end do
        end if

      end do ! ip
    end do ! ist

    POP_SUB(X(derivatives_batch_set_bc).multiresolution)
  end subroutine multiresolution


  ! ---------------------------------------------------------
  subroutine periodic()
    integer :: ip, ist
    R_TYPE, pointer :: ff(:)

#ifdef HAVE_MPI
    integer :: ipart, nreq
    integer, allocatable :: req(:), statuses(:, :)
#endif

    PUSH_SUB(X(derivatives_batch_set_bc).periodic)

#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains) then
      call profiling_in(set_bc_comm_prof, 'SET_BC_COMMUNICATION')
      
      ! get the points that are copies from other nodes
      SAFE_ALLOCATE(req(1:2 * der%mesh%vp%npart * ffb%nst_linear))

      nreq = 0

      do ist = 1, ffb%nst_linear
        do ipart = 1, der%mesh%vp%npart
          if(ipart == der%mesh%vp%partno .or. der%boundaries%nrecv(ipart) == 0) cycle

          ff => ffb%states_linear(ist)%X(psi)
          nreq = nreq + 1
          call MPI_Irecv(ff(1), 1, der%boundaries%X(recv_type)(ipart), ipart - 1, 3, &
            der%mesh%vp%comm, req(nreq), mpi_err)
        end do
      end do
      
      do ist = 1, ffb%nst_linear
        do ipart = 1, der%mesh%vp%npart
          if(ipart == der%mesh%vp%partno .or. der%boundaries%nsend(ipart) == 0) cycle

          ff => ffb%states_linear(ist)%X(psi)
          nreq = nreq + 1
          call MPI_Isend(ff(1), 1, der%boundaries%X(send_type)(ipart), ipart - 1, 3, &
            der%mesh%vp%comm, req(nreq), mpi_err)
        end do
      end do
      
      call profiling_count_transfers(  &
        sum(der%boundaries%nsend(1:der%mesh%vp%npart) + der%boundaries%nrecv(1:der%mesh%vp%npart)) * ffb%nst_linear, &
        R_TOTYPE(M_ONE))
      
      call profiling_out(set_bc_comm_prof)
    end if
#endif
    
    do ist = 1, ffb%nst_linear
      ff => ffb%states_linear(ist)%X(psi)
      forall (ip = 1:der%boundaries%nper)
        ff(der%boundaries%per_points(ip)) = ff(der%boundaries%per_map(ip))
      end forall
    end do
    
#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains .and. nreq > 0) then
      call profiling_in(set_bc_comm_prof)
      
      ! if nreq = 0, 1:nreq generates an illegal array reference
      SAFE_ALLOCATE(statuses(1:MPI_STATUS_SIZE, 1:nreq))
      call MPI_Waitall(nreq, req(1), statuses(1, 1), mpi_err)
      SAFE_DEALLOCATE_A(statuses)
      SAFE_DEALLOCATE_A(req)
      
      call profiling_out(set_bc_comm_prof)
    end if
#endif

    POP_SUB(X(derivatives_batch_set_bc).periodic)
  end subroutine periodic

end subroutine X(derivatives_batch_set_bc)


! ---------------------------------------------------------
subroutine X(derivatives_set_bc)(der, ff)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:)

  type(batch_t) :: batch_ff

  PUSH_SUB(X(derivatives_set_bc))

  call batch_init     (batch_ff, 1)
  call batch_add_state(batch_ff, ff)

  ASSERT(batch_is_ok(batch_ff))

  call X(derivatives_batch_set_bc)(der, batch_ff)

  call batch_end(batch_ff)
  POP_SUB(X(derivatives_set_bc))
end subroutine X(derivatives_set_bc)


! ---------------------------------------------------------
! These are the workhorse routines that handle the calculation of derivatives
subroutine X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc, factor)
  type(nl_operator_t),      target, intent(in)    :: op
  type(derivatives_t),      target, intent(in)    :: der
  type(batch_t),            target, intent(inout) :: ff
  type(batch_t),            target, intent(inout) :: opff
  type(derivatives_handle_batch_t), intent(out)   :: handle
  logical,                optional, intent(in)    :: ghost_update
  logical,                optional, intent(in)    :: set_bc
  R_TYPE,                 optional, intent(in)    :: factor

  PUSH_SUB(X(derivatives_batch_start))

  handle%ghost_update = optional_default(ghost_update, .true.)

  handle%op   => op
  handle%der  => der
  handle%ff   => ff
  handle%opff => opff

  if(present(factor)) then
    handle%factor_present = .true.
    handle%X(factor) = factor
  else
    handle%factor_present = .false.
  end if

  ASSERT(handle%ff%nst_linear == handle%opff%nst_linear)

  if(optional_default(set_bc, .true.)) call X(derivatives_batch_set_bc)(der, ff)

#ifdef HAVE_MPI

  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle%ghost_update) then
    call X(ghost_update_batch_start)(der%mesh%vp, ff, handle%pv_h)
  end if
#endif

  POP_SUB(X(derivatives_batch_start))
end subroutine X(derivatives_batch_start)


! ---------------------------------------------------------
subroutine X(derivatives_batch_finish)(handle)
  type(derivatives_handle_batch_t), intent(inout) :: handle

  logical :: done

  PUSH_SUB(X(derivatives_batch_finish))

  done = .false.

#ifdef HAVE_MPI
  if(derivatives_overlap(handle%der) .and. handle%der%mesh%parallel_in_domains .and. handle%ghost_update) then

    if(batch_status(handle%ff) /= BATCH_CL_PACKED) then
      done = .true.

      if(handle%factor_present) then
        call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, &
          ghost_update=.false., points=OP_INNER, factor = handle%X(factor))
      else
        call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update=.false., points=OP_INNER)
      end if
    end if

    call X(ghost_update_batch_finish)(handle%pv_h)

    if(batch_status(handle%ff) /= BATCH_CL_PACKED) then
      if(handle%factor_present) then
        call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, &
          ghost_update = .false., points = OP_OUTER, factor = handle%X(factor))
      else
        call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update = .false., points = OP_OUTER)
      end if

      ASSERT(done)
    end if

  end if
#endif

  if(.not. done) then
    if(handle%factor_present) then
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, &
        ghost_update = handle%ghost_update, factor = handle%X(factor))
    else
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update = handle%ghost_update)
    end if
  end if

  POP_SUB(X(derivatives_batch_finish))
end subroutine X(derivatives_batch_finish)


! ---------------------------------------------------------
subroutine X(derivatives_batch_perform)(op, der, ff, opff, ghost_update, set_bc)
  type(nl_operator_t), intent(in)    :: op
  type(derivatives_t), intent(in)    :: der
  type(batch_t),       intent(inout) :: ff
  type(batch_t),       intent(inout) :: opff
  logical,   optional, intent(in)    :: ghost_update
  logical,   optional, intent(in)    :: set_bc

  type(derivatives_handle_batch_t) :: handle

  PUSH_SUB(X(derivatives_batch_perform))

  call X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc)
  call X(derivatives_batch_finish)(handle)

  POP_SUB(X(derivatives_batch_perform))

end subroutine X(derivatives_batch_perform)


! ---------------------------------------------------------
! Now the simplified interfaces
subroutine X(derivatives_perform)(op, der, ff, op_ff, ghost_update, set_bc)
  type(nl_operator_t), target, intent(in)    :: op
  type(derivatives_t),         intent(in)    :: der
  R_TYPE,                      intent(inout) :: ff(:)     ! ff(der%mesh%np_part)
  R_TYPE,                      intent(out)   :: op_ff(:)  ! op_ff(der%mesh%np)
  logical, optional,           intent(in)    :: ghost_update
  logical, optional,           intent(in)    :: set_bc

  type(batch_t) :: batch_ff, batch_op_ff

  PUSH_SUB(X(derivatives_perform))

  call batch_init     (batch_ff, 1)
  call batch_add_state(batch_ff, ff)

  call batch_init     (batch_op_ff, 1)
  call batch_add_state(batch_op_ff, op_ff)

  ASSERT(batch_is_ok(batch_ff))
  ASSERT(batch_is_ok(batch_op_ff))

  call X(derivatives_batch_perform) (op, der, batch_ff, batch_op_ff, ghost_update, set_bc)

  call batch_end(batch_ff)
  call batch_end(batch_op_ff)
        
  POP_SUB(X(derivatives_perform))

end subroutine X(derivatives_perform)


! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  R_TYPE,                    intent(inout) :: ff(:)     ! f(der%mesh%np_part)
  R_TYPE,                    intent(out)   :: op_ff(:)  ! lapl(der%mesh%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  PUSH_SUB(X(derivatives_lapl))

  call X(derivatives_perform)(der%lapl, der, ff, op_ff, ghost_update, set_bc)
        
  POP_SUB(X(derivatives_lapl))
end subroutine X(derivatives_lapl)


! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:)        ! ff(der%mesh%np_part)
  R_TYPE,              intent(out)   :: op_ff(:, :)  ! op_ff(der%mesh%np, der%mesh%sb%dim)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  integer :: idir
  logical :: set_bc_, ghost_update_
  
  PUSH_SUB(X(derivatives_grad))
  call profiling_in(gradient_prof, "GRADIENT")

  ASSERT(ubound(op_ff, DIM=2) >= der%dim)

  set_bc_ = optional_default(set_bc, .true.)
  ghost_update_ = optional_default(ghost_update, .true.)
    
  do idir = 1, der%dim
    call X(derivatives_perform) (der%grad(idir), der, ff, op_ff(:, idir), ghost_update_, set_bc_)

    set_bc_       = .false. ! there is no need to update again
    ghost_update_ = .false. ! the boundary or ghost points
  end do

  call profiling_out(gradient_prof)
  POP_SUB(X(derivatives_grad))
end subroutine X(derivatives_grad)


! ---------------------------------------------------------
subroutine X(derivatives_div)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:,:)   ! ff(der%mesh%np_part, der%mesh%sb%dim)
  R_TYPE,              intent(out)   :: op_ff(:)  ! op_ff(der%mesh%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  R_TYPE, allocatable :: tmp(:)
  integer             :: idir, ii

  PUSH_SUB(X(derivatives_div))
  call profiling_in(divergence_prof, "DIVERGENCE")

  ASSERT(ubound(ff, DIM=2) >= der%dim)

  call X(derivatives_perform) (der%grad(1), der, ff(:, 1), op_ff, ghost_update, set_bc)

  SAFE_ALLOCATE(tmp(1:der%mesh%np))

  do idir = 2, der%dim
    call X(derivatives_perform) (der%grad(idir), der, ff(:, idir), tmp, ghost_update, set_bc)

    forall(ii = 1:der%mesh%np) op_ff(ii) = op_ff(ii) + tmp(ii)
  end do

  SAFE_DEALLOCATE_A(tmp)

  call profiling_out(divergence_prof)
  POP_SUB(X(derivatives_div))
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:,:)    ! ff(der%mesh%np_part, der%dim) 
  R_TYPE,              intent(out)   :: op_ff(:,:)
    ! if dim = 2, op_ff(der%mesh%np, der%dim)
    ! if dim = 1, op_ff(der%mesh%np, 1)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  integer, parameter  :: curl_dim(3) = (/-1, 1, 3/)
  R_TYPE, allocatable :: tmp(:)
  integer             :: ii, np

  PUSH_SUB(X(derivatives_curl))
  call profiling_in(curl_prof, "CURL")

  ASSERT(der%dim==2 .or. der%dim==3)
  ASSERT(ubound(ff,    DIM=2) >= der%dim)
  ASSERT(ubound(op_ff, DIM=2) >= curl_dim(der%dim))

  SAFE_ALLOCATE(tmp(1:der%mesh%np_part))

  op_ff(:,:) = R_TOTYPE(M_ZERO)
  np = der%mesh%np

  select case(der%dim)
  case(3)
    call X(derivatives_perform) (der%grad(3), der, ff(:,1), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 2) = op_ff(ii, 2) + tmp(ii)
    call X(derivatives_perform) (der%grad(2), der, ff(:,1), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 3) = op_ff(ii, 3) - tmp(ii)

    call X(derivatives_perform) (der%grad(3), der, ff(:,2), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) - tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,2), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 3) = op_ff(ii, 3) + tmp(ii)

    call X(derivatives_perform) (der%grad(2), der, ff(:,3), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) + tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,3), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 2) = op_ff(ii, 2) - tmp(ii)

  case(2)
    call X(derivatives_perform) (der%grad(2), der, ff(:,1), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) - tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,2), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) + tmp(ii)
  end select

  SAFE_DEALLOCATE_A(tmp)
  call profiling_out(curl_prof)
  POP_SUB(X(derivatives_curl))
end subroutine X(derivatives_curl)


! ----------------------------------------------------------

subroutine X(derivatives_test)(this)
  type(derivatives_t), intent(in)  :: this

  R_TYPE, allocatable :: ff(:), opff(:, :)
  R_TYPE :: aa, bb, cc
  integer :: ip, idir, ist
  type(batch_t) :: ffb, opffb
  integer :: blocksize, max_blocksize, itime, times
  logical :: packstates
  real(8) :: stime, etime
  character(len=7) :: type

  call parse_logical(datasets_check('StatesPack'), .true., packstates)

  SAFE_ALLOCATE(ff(1:this%mesh%np_part))
  SAFE_ALLOCATE(opff(1:this%mesh%np, 1:this%mesh%sb%dim))

#ifdef R_TREAL
  type = 'real'
#else
  type = 'complex'
#endif

  ! Note: here we need to use a constant function or anything that
  ! is constant at the borders, since we assume that all boundary
  ! points have equal values to optimize the application of the nl-operator.

  aa = CNST(1.0)/this%mesh%sb%lsize(1)
  bb = CNST(10.0)
  cc = CNST(100.0)

#ifdef R_TCOMPLEX
  ! we make things more "complex"
  aa = aa + M_ZI*CNST(0.01)
  bb = bb*exp(M_ZI*CNST(0.345))
  cc = cc - M_ZI*CNST(50.0)
#endif

  forall(ip = 1:this%mesh%np_part) ff(ip) = bb*exp(-aa*sum(this%mesh%x(ip, :)**2)) + cc

  call parse_integer(datasets_check('TestMinBlockSize'), 1, blocksize)
  call parse_integer(datasets_check('TestMaxBlockSize'), 128, max_blocksize)
  call parse_integer(datasets_check('TestRepetitions'), 10, times)

  do 

    call batch_init(ffb, 1, blocksize)
    call X(batch_new)(ffb, 1, blocksize, this%mesh%np_part)

    call batch_init(opffb, 1, blocksize)
    call X(batch_new)(opffb, 1, blocksize, this%mesh%np)

    forall(ist = 1:blocksize, ip = 1:this%mesh%np_part)
      ffb%states_linear(ist)%X(psi)(ip) = ff(ip)
    end forall

    if(packstates) then
      call batch_pack(ffb)
      call batch_pack(opffb, copy = .false.)
    end if

    stime = loct_clock()
    do itime = 1, times
      call X(derivatives_batch_perform)(this%lapl, this, ffb, opffb, set_bc = .false.)
    end do
    etime = (loct_clock() - stime)/dble(times)

    if(packstates) then
      call batch_unpack(ffb, copy = .false.)
      call batch_unpack(opffb)
    end if

    forall(ip = 1:this%mesh%np) 
      opffb%states_linear(blocksize)%X(psi)(ip) = opffb%states_linear(blocksize)%X(psi)(ip) - &
        (M_FOUR*aa**2*bb*sum(this%mesh%x(ip, :)**2)*exp(-aa*sum(this%mesh%x(ip, :)**2)) &
        - this%mesh%sb%dim*M_TWO*aa*bb*exp(-aa*sum(this%mesh%x(ip, :)**2)))
    end forall

    write(message(1), '(3a,i3,a,es17.10,a,f8.3)') &
      'Laplacian ', trim(type),  &
      ' bsize = ', blocksize,    &
      ' , error = ', X(mf_nrm2)(this%mesh, opffb%states_linear(blocksize)%X(psi)), &
      ' , gflops = ',  &
#ifdef R_TREAL
      blocksize*this%mesh%np*CNST(2.0)*this%lapl%stencil%size/(etime*CNST(1.0e9))
#else
      blocksize*this%mesh%np*CNST(4.0)*this%lapl%stencil%size/(etime*CNST(1.0e9))
#endif

    call messages_info(1)

    call batch_end(ffb)
    call batch_end(opffb)

    blocksize = 2*blocksize
    if(blocksize > max_blocksize) exit

  end do

  call X(derivatives_grad)(this, ff, opff, set_bc = .false.)

  forall(idir = 1:this%mesh%sb%dim, ip = 1:this%mesh%np) 
    opff(ip, idir) = opff(ip, idir) - (-M_TWO*aa*bb*this%mesh%x(ip, idir)*exp(-aa*sum(this%mesh%x(ip, :)**2)))
  end forall

  message(1) = ''
  call messages_info(1)


  write(message(1), '(3a, es17.10)') 'Gradient ', trim(type),  &
    ' err = ', X(mf_nrm2)(this%mesh, this%mesh%sb%dim, opff)
  call messages_info(1)

  message(1) = ''
  call messages_info(1)

  SAFE_DEALLOCATE_A(ff)
  SAFE_DEALLOCATE_A(opff)

end subroutine X(derivatives_test)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
