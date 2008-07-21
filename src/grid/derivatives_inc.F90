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

! This module calculates the derivatives (gradients, laplacians, etc.) 
! of a function. Note that the function whose derivative is to be calculated
! *has* to be defined (1:m%np_part), while the (1:m%np) values of the derivative
! are calculated. This was made to simplify the parallel mode, and has to be
! followed by the rest of the code.

! ---------------------------------------------------------
! The transpose of the Laplacian.
subroutine X(derivatives_laplt)(der, f, lapl)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,            intent(out)   :: lapl(:)  ! lapl(m%np)

  call push_sub('derivatives_inc.Xderivatives_laplt')

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  call X(set_bc)(der, f)

  call X(nl_operator_operate) (der%laplt, f, lapl)

  call pop_sub()
end subroutine X(derivatives_laplt)

subroutine X(init_f)(der, handle, f)
  type(der_discr_t),         intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle
  R_TYPE,  target,           intent(in)    :: f(:)

  if(ubound(f, DIM = 1) == der%m%np_part) then
    handle%X(f) => f
    handle%dealloc_f = .false.
  else
    ASSERT(ubound(f, DIM = 1) == der%m%np)
    ALLOCATE(handle%X(f)(1:der%m%np_part), der%m%np_part)
    call lalg_copy(der%m%np, f, handle%X(f))
    handle%dealloc_f = .true.
  end if

end subroutine X(init_f)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_start)(der, handle, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,  target,           intent(inout) :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_

  call push_sub('derivatives_inc.Xderivatives_lapl_start')

  ASSERT(ubound(lapl, DIM=1) >= der%m%np)
  handle%X(lapl) => lapl

  call X(init_f)(der, handle, f)

  handle%ghost_update = .true.
  if(present(ghost_update)) handle%ghost_update = ghost_update

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(set_bc)(der, handle%X(f))

#ifdef HAVE_MPI
  if(der%overlap .and. der%m%parallel_in_domains .and. handle%ghost_update) then
    call X(vec_ighost_update)(der%m%vp,  handle%X(f), handle%pv_h)
    call X(nl_operator_operate)(der%lapl,  handle%X(f), handle%X(lapl), ghost_update = .false., points = OP_INNER)
  end if
#endif

  call pop_sub()
end subroutine X(derivatives_lapl_start)

! --------------------------------------------------------
! Call NBC_Test to give MPI a chance to process pending messages.
! This improves the obverlap but is only necessary due to deficiencies
! in MPI implementations.
! ---------------------------------------------------------
subroutine X(derivatives_lapl_keep_going)(der, handle)
  type(der_discr_t),         intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl_keep_going')

#ifdef HAVE_MPI
  if(der%overlap .and. der%m%parallel_in_domains .and. handle%ghost_update) then
    call pv_handle_test(handle%pv_h)
  end if
#endif

  call pop_sub()
end subroutine X(derivatives_lapl_keep_going)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_finish)(der, handle)
  type(der_discr_t),         intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl_finish')

#ifdef HAVE_MPI
  if(der%overlap .and. der%m%parallel_in_domains .and. handle%ghost_update) then

    call pv_handle_wait(handle%pv_h)
    call X(nl_operator_operate) (der%lapl, handle%X(f), handle%X(lapl), ghost_update = .false., points = OP_OUTER)
    
    call pop_sub()
    return
  end if
#endif

  call X(nl_operator_operate) (der%lapl, handle%X(f), handle%X(lapl), ghost_update = handle%ghost_update)
  
  if(handle%dealloc_f) deallocate(handle%X(f))

  call pop_sub()
end subroutine X(derivatives_lapl_finish)

! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  type(der_handle_t) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl')

  call der_handle_init(handle, der%m)
  call X(derivatives_lapl_start) (der, handle, f, lapl, ghost_update, set_bc)
  call X(derivatives_lapl_finish)(der, handle)
  call der_handle_end(handle)
        
  call pop_sub()
end subroutine X(derivatives_lapl)

! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, f, grad, ghost_update, set_bc)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,            intent(out)   :: grad(:,:)  ! grad(m%np, m%sb%dim)
  logical, optional, intent(in)    :: ghost_update
  logical, optional, intent(in)    :: set_bc

  type(der_handle_t) :: handle
  integer :: idir
  logical :: set_bc_
  
#ifdef HAVE_MPI
  logical :: ghost_update_
#endif

  call push_sub('derivatives_inc.Xderivatives_grad')

  call der_handle_init(handle, der%m)
  call X(init_f)(der, handle, f)

  ASSERT(ubound(grad, DIM=1) >= der%m%np)
  ASSERT(ubound(grad, DIM=2) >= der%m%sb%dim)

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(set_bc)(der, handle%X(f))

#ifdef HAVE_MPI
  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

  if(der%overlap .and. der%m%parallel_in_domains .and. ghost_update_) then

    call X(vec_ighost_update)(der%m%vp, handle%X(f), handle%pv_h)

    do idir = 1, der%m%sb%dim
      call X(nl_operator_operate)(der%grad(idir), handle%X(f), grad(:, idir), ghost_update = .false., points = OP_INNER)
    end do

    call pv_handle_wait(handle%pv_h)

    do idir = 1, der%m%sb%dim
      call X(nl_operator_operate)(der%grad(idir), handle%X(f), grad(:, idir), ghost_update = .false., points = OP_OUTER)
    end do

  else
#endif
    
    do idir = 1, der%m%sb%dim
      call X(nl_operator_operate) (der%grad(idir), handle%X(f), grad(:, idir), ghost_update = ghost_update)
    end do
    
#ifdef HAVE_MPI
  end if
#endif

  if(handle%dealloc_f) deallocate(handle%X(f))
  call der_handle_end(handle)

  call pop_sub()
end subroutine X(derivatives_grad)

! ---------------------------------------------------------
subroutine X(derivatives_oper)(op, der, f, opf, ghost_update, set_bc)
  type(nl_operator_t), intent(in)    :: op
  type(der_discr_t),   intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,              intent(out)   :: opf(:)     ! opf(m%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  type(der_handle_t) :: handle
#ifdef HAVE_MPI
  logical :: ghost_update_
#endif

  call push_sub('derivatives_inc.Xderivatives_oper')
  
  ASSERT(ubound(opf, DIM=1) >= der%m%np)

  call der_handle_init(handle, der%m)
  call X(init_f)(der, handle, f)

  if (.not. present(set_bc) .or. set_bc) call X(set_bc)(der, handle%X(f))

#ifdef HAVE_MPI
  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

  if(der%overlap .and. der%m%parallel_in_domains .and. ghost_update_) then

    call X(vec_ighost_update)(der%m%vp, handle%X(f), handle%pv_h)

    call X(nl_operator_operate)(op, handle%X(f), opf, ghost_update = .false., points = OP_INNER)

    call pv_handle_wait(handle%pv_h)

    call X(nl_operator_operate)(op, handle%X(f), opf, ghost_update = .false., points = OP_OUTER)

  else
#endif

    call X(nl_operator_operate) (op, handle%X(f), opf, ghost_update = ghost_update)

#ifdef HAVE_MPI
  end if
#endif

  if(handle%dealloc_f) deallocate(handle%X(f))
  call der_handle_end(handle)

  call pop_sub()
end subroutine X(derivatives_oper)


! ---------------------------------------------------------
subroutine X(derivatives_div)(der, f, div, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:,:)   ! f(m%np_part, m%sb%dim)
  R_TYPE,            intent(out)   :: div(:)   ! div(m%np)
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(ubound(f,   DIM=1) == der%m%np_part)
  ASSERT(ubound(div, DIM=1) >= der%m%np)

  do i = 1, der%m%sb%dim
    call X(set_bc)(der, f(:, i))
  end do

  ALLOCATE(tmp(der%m%np), der%m%np)

  div(:) = R_TOTYPE(M_ZERO)
  do i = 1, der%m%sb%dim
    call X(nl_operator_operate) (der%grad(i), f(:,i), tmp, ghost_update=ghost_update)
    div(:) = div(:) + tmp(:)
  end do

  deallocate(tmp)

  call pop_sub()
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, f, curl, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:,:)    ! f(m%np_part, der%m%sb%dim) 
  R_TYPE,            intent(out)   :: curl(:,:) ! curl(m%np, der%m%sb%dim) if dim = 2, curl(m%np, 1) if dim = 1.
  logical, optional, intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i, np

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(der%m%sb%dim == 3 .or. der%m%sb%dim == 2)
  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(curl, DIM=1) >= der%m%np)
  select case(der%m%sb%dim)
    case(3)
      ASSERT(ubound(curl, DIM=2) == der%m%sb%dim)
    case(2)
      ASSERT(ubound(curl, DIM=2) == 1)
    case(1)
      write(message(1),'(a)') 'INTERNAL ERROR at Xderivatives_curl: 1D not allowed'
      call write_fatal(1)
  end select

  do i = 1, der%m%sb%dim
    call X(set_bc)(der, f(:, i))
  end do
  
  ALLOCATE(tmp(der%m%np_part), der%m%np_part)

  curl(:,:) = R_TOTYPE(M_ZERO)
  np = der%m%np

  select case(der%m%sb%dim)
  case(3)
    call X(nl_operator_operate) (der%grad(3), f(:,1), tmp, ghost_update=ghost_update)
    curl(1:np,2) = curl(1:np,2) + tmp(1:np)
    call X(nl_operator_operate) (der%grad(2), f(:,1), tmp, ghost_update=ghost_update)
    curl(1:np,3) = curl(1:np,3) - tmp(1:np)

    call X(nl_operator_operate) (der%grad(3), f(:,2), tmp, ghost_update=ghost_update)
    curl(1:np,1) = curl(1:np,1) - tmp(1:np)
    call X(nl_operator_operate) (der%grad(1), f(:,2), tmp, ghost_update=ghost_update)
    curl(1:np,3) = curl(1:np,3) + tmp(1:np)

    call X(nl_operator_operate) (der%grad(2), f(:,3), tmp, ghost_update=ghost_update)
    curl(1:np,1) = curl(1:np,1) + tmp(1:np)
    call X(nl_operator_operate) (der%grad(1), f(:,3), tmp, ghost_update=ghost_update)
    curl(1:np,2) = curl(1:np,2) - tmp(1:np)
  case(2)
    call X(nl_operator_operate) (der%grad(2), f(:,1), tmp, ghost_update=ghost_update)
    curl(1:np,1) = curl(1:np,1) - tmp(1:np)
    call X(nl_operator_operate) (der%grad(1), f(:,2), tmp, ghost_update=ghost_update)
    curl(1:np,1) = curl(1:np,1) + tmp(1:np)
  end select

  deallocate(tmp)
  call pop_sub()
end subroutine X(derivatives_curl)


! ---------------------------------------------------------
! Set all boundary points in f to zero to implement zero
! boundary conditions for the derivatives.
subroutine X(set_bc)(der, f)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:)

  integer :: bndry_start, bndry_end
  integer :: p
  integer :: iper
  type(profile_t), save :: set_bc_prof
#ifdef HAVE_MPI
  integer :: ipart, nreq, status(MPI_STATUS_SIZE)
  integer, allocatable :: req(:), statuses(:, :)
#endif

  call push_sub('derivatives_inc.Xset_bc')
  call profiling_in(set_bc_prof, 'SET_BC')
  
  ASSERT(ubound(f, DIM=1) == der%m%np_part)

  p = der%m%vp%partno
   
  if(der%zero_bc) then
    
    ! The boundary points are at different locations depending on the presence
    ! of ghost points due to domain parallelization.
    if(der%m%parallel_in_domains) then
      bndry_start = der%m%vp%np_local(p) + der%m%vp%np_ghost(p) + 1
      bndry_end   = der%m%vp%np_local(p) + der%m%vp%np_ghost(p) + der%m%vp%np_bndry(p)
    else
      bndry_start = der%m%np+1
      bndry_end   = der%m%np_part
    end if
    
    !$omp parallel workshare
    f(bndry_start:bndry_end) = R_TOTYPE(M_ZERO)
    !$omp end parallel workshare
    
  end if

  if(der%periodic_bc) then

#ifdef HAVE_MPI
    if(der%m%parallel_in_domains) then
      ! get the points from other nodes
      ALLOCATE(req(2*der%m%vp%p), 2*der%m%vp%p)

      nreq = 0

      do ipart = 1, der%m%vp%p
        if(ipart == p .or. der%m%nsend(ipart) == 0) cycle
        nreq = nreq + 1
        call MPI_Isend(f, 1, der%m%X(send_type)(ipart), ipart - 1, 3, der%m%vp%comm, req(nreq), mpi_err)
      end do

      do ipart = 1, der%m%vp%p
        if(ipart == p .or. der%m%nrecv(ipart) == 0) cycle
        call MPI_Recv(f(der%m%np + 1:), 1, der%m%X(recv_type)(ipart), ipart - 1, 3, der%m%vp%comm, status, mpi_err)
      end do

    end if
#endif

    !$omp parallel do
    do iper = 1, der%m%nper
      f(der%m%per_points(iper)) = f(der%m%per_map(iper))
    end do
    !$omp end parallel do

#ifdef HAVE_MPI
    if(der%m%parallel_in_domains) then
      ALLOCATE(statuses(MPI_STATUS_SIZE, nreq), MPI_STATUS_SIZE*nreq)
      call MPI_Waitall(nreq, req, statuses, mpi_err)
    end if
#endif

  end if

  call profiling_out(set_bc_prof)
  call pop_sub()

end subroutine X(set_bc)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
