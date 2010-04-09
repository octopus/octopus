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
subroutine X(derivatives_batch_set_bc)(der, batch_ff)
  type(derivatives_t), intent(in)    :: der
  type(batch_t),       intent(inout) :: batch_ff

  integer :: pp, bndry_start, bndry_end

  call push_sub('derivatives_inc.Xderivatives_batch_set_bc')
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
  call pop_sub('derivatives_inc.Xderivatives_batch_set_bc')

contains
  ! ---------------------------------------------------------
  subroutine zero_boundaries()
    integer :: ist, ip

    call push_sub('derivatives_inc.Xderivatives_batch_set_bc.zero_boundaries')

    forall (ist=1:batch_ff%nst_linear, ip=bndry_start:bndry_end)
      batch_ff%states_linear(ist)%X(psi)(ip) = R_TOTYPE(M_ZERO)
    end forall

    call pop_sub('derivatives_inc.Xderivatives_batch_set_bc.zero_boundaries')
  end subroutine zero_boundaries


  ! ---------------------------------------------------------
  subroutine multiresolution()
    integer :: ist, ip
    integer :: ii, jj, kk, ix, iy, iz, dx, dy, dz, i_lev
    FLOAT :: weight
    R_TYPE, pointer :: ff(:)

    call push_sub('derivatives_inc.Xderivatives_batch_set_bc.multiresolution')

    do ist = 1, batch_ff%nst_linear
      ff => batch_ff%states_linear(ist)%X(psi)

      do ip = bndry_start, bndry_end
        ix = der%mesh%idx%Lxyz(ip, 1)
        iy = der%mesh%idx%Lxyz(ip, 2)
        iz = der%mesh%idx%Lxyz(ip, 3)

        i_lev = der%mesh%resolution(ix,iy,iz)

        ! resolution is 2**num_radii for outer boundary points, but now we want inner boundary points
        if(i_lev.ne.2**der%mesh%sb%hr_area%num_radii) then
          dx = abs(mod(ix, 2**(i_lev)))
          dy = abs(mod(iy, 2**(i_lev)))
          dz = abs(mod(iz, 2**(i_lev)))

          do ii = 1, der%mesh%sb%hr_area%interp%nn
            do jj = 1, der%mesh%sb%hr_area%interp%nn
              do kk = 1, der%mesh%sb%hr_area%interp%nn
                weight = der%mesh%sb%hr_area%interp%ww(ii)* &
                  der%mesh%sb%hr_area%interp%ww(jj)*        &
                  der%mesh%sb%hr_area%interp%ww(kk)

                ff(ip) = ff(ip) + weight * ff(der%mesh%idx%Lxyz_inv(   &
                  ix + der%mesh%sb%hr_area%interp%posi(ii)*dx,        &
                  iy + der%mesh%sb%hr_area%interp%posi(jj)*dy,        &
                  iz + der%mesh%sb%hr_area%interp%posi(kk)*dz))
              end do
            end do
          end do
        end if

      end do ! ip
    end do ! ist

    call pop_sub('derivatives_inc.Xderivatives_batch_set_bc.multiresolution')
  end subroutine multiresolution


  ! ---------------------------------------------------------
  subroutine periodic()
    integer :: ip, ist
    R_TYPE, pointer :: ff(:)

#ifdef HAVE_MPI
    integer :: ipart, nreq
    integer, allocatable :: req(:), statuses(:, :)
#endif

    call push_sub('derivatives_inc.Xderivatives_batch_set_bc.periodic')

#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains) then
      call profiling_in(set_bc_comm_prof, 'SET_BC_COMMUNICATION')
      
      ! get the points that are copies from other nodes
      SAFE_ALLOCATE(req(1:2*der%mesh%vp%npart*batch_ff%nst_linear))

      nreq = 0

      do ist = 1, batch_ff%nst_linear
        do ipart = 1, der%mesh%vp%npart
          if(ipart == der%mesh%vp%partno .or. der%mesh%nrecv(ipart) == 0) cycle

          ff => batch_ff%states_linear(ist)%X(psi)
          nreq = nreq + 1
          call MPI_Irecv(ff(1), 1, der%mesh%X(recv_type)(ipart), ipart - 1, 3, &
            der%mesh%vp%comm, req(nreq), mpi_err)
        end do
      end do
      
      do ist = 1, batch_ff%nst_linear
        do ipart = 1, der%mesh%vp%npart
          if(ipart == der%mesh%vp%partno .or. der%mesh%nsend(ipart) == 0) cycle

          ff => batch_ff%states_linear(ist)%X(psi)
          nreq = nreq + 1
          call MPI_Isend(ff(1), 1, der%mesh%X(send_type)(ipart), ipart - 1, 3, &
            der%mesh%vp%comm, req(nreq), mpi_err)
        end do
      end do
      
      call profiling_count_transfers(  &
        sum(der%mesh%nsend(1:der%mesh%vp%npart) + der%mesh%nrecv(1:der%mesh%vp%npart))*batch_ff%nst_linear, &
        R_TOTYPE(M_ONE))
      
      call profiling_out(set_bc_comm_prof)
    end if
#endif
    
    do ist = 1, batch_ff%nst_linear
      ff => batch_ff%states_linear(ist)%X(psi)
      forall (ip = 1:der%mesh%nper)
        ff(der%mesh%per_points(ip)) = ff(der%mesh%per_map(ip))
      end forall
    end do
    
#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains) then
      call profiling_in(set_bc_comm_prof)
      
      SAFE_ALLOCATE(statuses(1:MPI_STATUS_SIZE, 1:nreq))
      call MPI_Waitall(nreq, req(1), statuses(1, 1), mpi_err)
      SAFE_DEALLOCATE_A(statuses)
      SAFE_DEALLOCATE_A(req)
      
      call profiling_out(set_bc_comm_prof)
    end if
#endif

    call pop_sub('derivatives_inc.Xderivatives_batch_set_bc.periodic')
  end subroutine periodic

end subroutine X(derivatives_batch_set_bc)


! ---------------------------------------------------------
subroutine X(derivatives_set_bc)(der, ff)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:)

  type(batch_t) :: batch_ff

  call push_sub('derivatives_inc.Xderivatives_set_bc')

  call batch_init     (batch_ff, 1)
  call batch_add_state(batch_ff, ff)

  call X(derivatives_batch_set_bc)(der, batch_ff)

  call batch_end(batch_ff)
  call pop_sub('derivatives_inc.Xderivatives_set_bc')
end subroutine X(derivatives_set_bc)


! ---------------------------------------------------------
! This are the workhorse routines that handle the calculation of derivatives
subroutine X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc)
  type(nl_operator_t),      target, intent(in)    :: op
  type(derivatives_t),      target, intent(in)    :: der
  type(batch_t),            target, intent(inout) :: ff
  type(batch_t),            target, intent(inout) :: opff
  type(derivatives_handle_batch_t), intent(out)   :: handle
  logical,                optional, intent(in)    :: ghost_update
  logical,                optional, intent(in)    :: set_bc

  logical :: set_bc_

  call push_sub('derivatives_inc.Xderivatives_batch_start')

  handle%ghost_update = .true.
  if(present(ghost_update)) handle%ghost_update = ghost_update

  handle%op   => op
  handle%der  => der
  handle%ff   => ff
  handle%opff => opff

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(derivatives_batch_set_bc)(der, ff)

#ifdef HAVE_MPI

  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle%ghost_update) then
    call X(ghost_update_batch_start)(der%mesh%vp, ff, der%comm_method, handle%pv_h)
  end if
#endif

  call pop_sub('derivatives_inc.Xderivatives_batch_start')
end subroutine X(derivatives_batch_start)


! ---------------------------------------------------------
subroutine X(derivatives_batch_finish)(handle)
  type(derivatives_handle_batch_t), intent(inout) :: handle

  call push_sub('derivatives_inc.Xderivatives_batch_finish')

#ifdef HAVE_MPI
  if(derivatives_overlap(handle%der) .and. handle%der%mesh%parallel_in_domains .and. handle%ghost_update) then
    call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update=.false., points=OP_INNER)

    call X(ghost_update_batch_finish)(handle%pv_h)

    call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update=.false., points=OP_OUTER)
  else
#endif
    
    call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update=handle%ghost_update)
    
#ifdef HAVE_MPI
  end if
#endif

  call pop_sub('derivatives_inc.Xderivatives_batch_finish')
end subroutine X(derivatives_batch_finish)


! ---------------------------------------------------------
subroutine X(derivatives_batch_perform)(op, der, ff, opff, ghost_update, set_bc)
  type(nl_operator_t),      target, intent(in)    :: op
  type(derivatives_t),      target, intent(in)    :: der
  type(batch_t),            target, intent(inout) :: ff
  type(batch_t),            target, intent(inout) :: opff
  logical,                optional, intent(in)    :: ghost_update
  logical,                optional, intent(in)    :: set_bc

  type(derivatives_handle_batch_t) :: handle

  call push_sub('derivatives_inc.Xderivatives_batch_perform')

  call X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc)
  call X(derivatives_batch_finish)(handle)

  call pop_sub('derivatives_inc.Xderivatives_batch_perform')

end subroutine X(derivatives_batch_perform)


! ---------------------------------------------------------
! Now the simplified interfaces
subroutine X(derivatives_perform)(op, der, ff, op_ff, ghost_update, set_bc)
  type(nl_operator_t), target, intent(in)    :: op
  type(derivatives_t),         intent(in)    :: der
  R_TYPE,                      intent(inout) :: ff(:)     ! ff(der%mesh%np_part)
  R_TYPE,                      intent(out)   :: op_ff(:)  ! lapl(der%mesh%np)
  logical, optional,           intent(in)    :: ghost_update
  logical, optional,           intent(in)    :: set_bc

  type(batch_t) :: batch_ff, batch_op_ff

  call push_sub('derivatives_inc.Xderivatives_perform')

  call batch_init     (batch_ff, 1)
  call batch_add_state(batch_ff, ff)

  call batch_init     (batch_op_ff, 1)
  call batch_add_state(batch_op_ff, op_ff)

  call X(derivatives_batch_perform) (op, der, batch_ff, batch_op_ff, ghost_update, set_bc)

  call batch_end(batch_ff)
  call batch_end(batch_op_ff)
        
  call pop_sub('derivatives_inc.Xderivatives_perform')

end subroutine X(derivatives_perform)


! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  R_TYPE,                    intent(inout) :: ff(:)     ! f(der%mesh%np_part)
  R_TYPE,                    intent(out)   :: op_ff(:)  ! lapl(der%mesh%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  call push_sub('derivatives_inc.Xderivatives_lapl')

  call X(derivatives_perform)(der%lapl, der, ff, op_ff, ghost_update, set_bc)
        
  call pop_sub('derivatives_inc.Xderivatives_lapl')
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
  
  call push_sub('derivatives_inc.Xderivatives_grad')
  call profiling_in(gradient_prof, "GRADIENT")

  ASSERT(ubound(op_ff, DIM=2) >= der%dim)

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc

  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update
    
  do idir = 1, der%dim
    call X(derivatives_perform) (der%grad(idir), der, ff, op_ff(:, idir), ghost_update_, set_bc_)

    set_bc_       = .false. ! there is no need to update again
    ghost_update_ = .false. ! the boundary or ghost points
  end do

  call profiling_out(gradient_prof)
  call pop_sub('derivatives_inc.Xderivatives_grad')
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

  call push_sub('derivatives_inc.Xderivatives_div')
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
  call pop_sub('derivatives_inc.Xderivatives_div')
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

  call push_sub('derivatives_inc.Xderivatives_curl')
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
  call pop_sub('derivatives_inc.Xderivatives_curl')
end subroutine X(derivatives_curl)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
