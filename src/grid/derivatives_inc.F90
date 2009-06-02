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
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,              intent(out)   :: lapl(:)  ! lapl(m%np)

  call push_sub('derivatives_inc.Xderivatives_laplt')

  ASSERT(ubound(f,    DIM=1) == der%mesh%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%mesh%np)

  call X(set_bc)(der, f)

  call X(nl_operator_operate) (der%laplt, f, lapl)

  call pop_sub()
end subroutine X(derivatives_laplt)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_start)(der, handle, f, lapl, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle
  R_TYPE,  target,           intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,  target,           intent(inout) :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_

  call push_sub('derivatives_inc.Xderivatives_lapl_start')

  ASSERT(ubound(f, DIM=1) >= der%mesh%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%mesh%np)

  handle%X(f) => f
  handle%X(lapl) => lapl
  
  handle%ghost_update = .true.
  if(present(ghost_update)) handle%ghost_update = ghost_update

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(set_bc)(der, handle%X(f))

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle%ghost_update) then
    call X(vec_ighost_update)(der%mesh%vp,  handle%X(f), handle%pv_h)
  end if
#endif

  call pop_sub()
end subroutine X(derivatives_lapl_start)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_finish)(der, handle)
  type(derivatives_t),       intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl_finish')

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle%ghost_update) then

    call X(nl_operator_operate)(der%lapl,  handle%X(f), handle%X(lapl), ghost_update = .false., points = OP_INNER)
    call pv_handle_wait(handle%pv_h)
    call X(nl_operator_operate) (der%lapl, handle%X(f), handle%X(lapl), ghost_update = .false., points = OP_OUTER)
    
    call pop_sub(); return
  end if
#endif

  call X(nl_operator_operate) (der%lapl, handle%X(f), handle%X(lapl), ghost_update = handle%ghost_update)
  
  call pop_sub()
end subroutine X(derivatives_lapl_finish)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_batch_start)(der, handle, ff, lapl, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle(:, :)
  type(batch_t),             intent(inout) :: ff
  type(batch_t),             intent(inout) :: lapl
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_
#ifdef HAVE_MPI
  integer :: ist, idim
#endif
  call push_sub('derivatives_inc.Xderivatives_lapl_start')


  handle(1, 1)%ghost_update = .true.
  if(present(ghost_update)) handle(1, 1)%ghost_update = ghost_update

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc

  if(set_bc_) call X(set_bc_batch)(der, ff)

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle(1, 1)%ghost_update) then
    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call X(vec_ighost_update)(der%mesh%vp, ff%states(ist)%X(psi)(:, idim), handle(idim, ist)%pv_h)
      end do
    end do
  end if
#endif

  call pop_sub()
end subroutine X(derivatives_lapl_batch_start)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_batch_finish)(der, handle, ff, lapl)
  type(derivatives_t),       intent(in)    :: der
  type(der_handle_t),        intent(inout) :: handle(:, :)
  type(batch_t),             intent(inout) :: ff
  type(batch_t),             intent(inout) :: lapl

#ifdef HAVE_MPI
  integer :: ist, idim
#endif

  call push_sub('derivatives_inc.Xderivatives_lapl_finish')

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle(1, 1)%ghost_update) then

    call X(nl_operator_operate_batch)(der%lapl, ff, lapl, ghost_update = .false., points = OP_INNER)
    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call pv_handle_wait(handle(idim, ist)%pv_h)
      end do
    end do
    call X(nl_operator_operate_batch)(der%lapl, ff, lapl, ghost_update = .false., points = OP_OUTER)
    
    call pop_sub(); return
  end if
#endif

  call X(nl_operator_operate_batch)(der%lapl, ff, lapl, ghost_update = handle(1, 1)%ghost_update)
  
  call pop_sub()
end subroutine X(derivatives_lapl_batch_finish)

! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  type(der_handle_t) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl')

  call der_handle_init(handle, der)
  call X(derivatives_lapl_start) (der, handle, f, lapl, ghost_update, set_bc)
  call X(derivatives_lapl_finish)(der, handle)
  call der_handle_end(handle)
        
  call pop_sub()
end subroutine X(derivatives_lapl)

! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, f, grad, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:)        ! f(mesh%np_part)
  R_TYPE,              intent(out)   :: grad(:, :)  ! grad(mesh%np, mesh%sb%dim)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  type(der_handle_t) :: handle
  integer :: idir
  logical :: set_bc_
  
#ifdef HAVE_MPI
  logical :: ghost_update_
#endif

  call push_sub('derivatives_inc.Xderivatives_grad')

  call der_handle_init(handle, der)

  ASSERT(ubound(grad, DIM=1) >= der%mesh%np)
  ASSERT(ubound(grad, DIM=2) >= der%dim)

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(set_bc)(der, f)

#ifdef HAVE_MPI
  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. ghost_update_) then

    call X(vec_ighost_update)(der%mesh%vp, f, handle%pv_h)

    do idir = 1, der%dim
      call X(nl_operator_operate)(der%grad(idir), f, grad(:, idir), ghost_update = .false., points = OP_INNER)
    end do

    call pv_handle_wait(handle%pv_h)

    do idir = 1, der%dim
      call X(nl_operator_operate)(der%grad(idir), f, grad(:, idir), ghost_update = .false., points = OP_OUTER)
    end do

  else
#endif
    
    do idir = 1, der%dim
      call X(nl_operator_operate) (der%grad(idir), f, grad(:, idir), ghost_update = ghost_update)
    end do
    
#ifdef HAVE_MPI
  end if
#endif

  call der_handle_end(handle)

  call pop_sub()
end subroutine X(derivatives_grad)

! ---------------------------------------------------------
subroutine X(derivatives_oper)(op, der, f, opf, ghost_update, set_bc)
  type(nl_operator_t),   intent(in)    :: op
  type(derivatives_t),   intent(in)    :: der
  R_TYPE,                intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,                intent(out)   :: opf(:)     ! opf(m%np)
  logical, optional,     intent(in)    :: ghost_update
  logical, optional,     intent(in)    :: set_bc

  type(der_handle_t) :: handle
  logical :: ghost_update_, set_bc_

  call push_sub('derivatives_inc.Xderivatives_oper')
  
  ASSERT(ubound(opf, DIM=1) >= der%mesh%np)

  call der_handle_init(handle, der)

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc

  if(set_bc_) call X(set_bc)(der, f)

  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. ghost_update_) then

    call X(vec_ighost_update)(der%mesh%vp, f, handle%pv_h)

    call X(nl_operator_operate)(op, f, opf, ghost_update = .false., points = OP_INNER)

    call pv_handle_wait(handle%pv_h)

    call X(nl_operator_operate)(op, f, opf, ghost_update = .false., points = OP_OUTER)

  else
#endif

    call X(nl_operator_operate) (op, f, opf, ghost_update = ghost_update_)

#ifdef HAVE_MPI
  end if
#endif

  call der_handle_end(handle)

  call pop_sub()
end subroutine X(derivatives_oper)

! ---------------------------------------------------------
subroutine X(derivatives_oper_batch)(op, der, ff, opff, ghost_update, set_bc)
  type(nl_operator_t),   intent(in)    :: op
  type(derivatives_t),   intent(in)    :: der
  type(batch_t),         intent(inout) :: ff
  type(batch_t),         intent(inout) :: opff
  logical, optional,     intent(in)    :: ghost_update
  logical, optional,     intent(in)    :: set_bc

  logical :: set_bc_
  logical :: ghost_update_
#ifdef HAVE_MPI
  type(der_handle_t), allocatable :: handle(:, :)
  integer :: ist, idim
#endif

  call push_sub('derivatives_inc.Xderivatives_oper')

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc

  if(set_bc_) call X(set_bc_batch)(der, ff)

  ghost_update_ = .true.
  if(present(ghost_update)) ghost_update_ = ghost_update

#ifdef HAVE_MPI
  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. ghost_update_) then
    SAFE_ALLOCATE(handle(1:ff%nst, 1:ff%dim))
    
    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call der_handle_init(handle(ist, idim), der)
      end do
    end do
    
    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call X(vec_ighost_update)(der%mesh%vp, ff%states(ist)%X(psi)(:, idim), handle(ist, idim)%pv_h)
      end do
    end do
    
    call X(nl_operator_operate_batch)(op, ff, opff, ghost_update = .false., points = OP_INNER)

    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call pv_handle_wait(handle(ist, idim)%pv_h)
      end do
    end do

    do ist = 1, ff%nst
      do idim = 1, ff%dim
        call der_handle_end(handle(ist, idim))
      end do
    end do
    
    call X(nl_operator_operate_batch)(op, ff, opff, ghost_update = .false., points = OP_OUTER)

  else
#endif

    call X(nl_operator_operate_batch) (op, ff, opff, ghost_update = ghost_update_)

#ifdef HAVE_MPI
  end if

  SAFE_DEALLOCATE_A(handle)
#endif

  call pop_sub()
end subroutine X(derivatives_oper_batch)

! ---------------------------------------------------------
subroutine X(derivatives_div)(der, f, div, ghost_update)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:,:)   ! f(m%np_part, m%sb%dim)
  R_TYPE,              intent(out)   :: div(:)   ! div(m%np)
  logical, optional,   intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(ubound(f,   DIM=1) == der%mesh%np_part)
  ASSERT(ubound(div, DIM=1) >= der%mesh%np)

  do i = 1, der%dim
    call X(set_bc)(der, f(:, i))
  end do

  SAFE_ALLOCATE(tmp(1:der%mesh%np))

  div(:) = R_TOTYPE(M_ZERO)
  do i = 1, der%dim
    call X(nl_operator_operate) (der%grad(i), f(:,i), tmp, ghost_update=ghost_update)
    div(:) = div(:) + tmp(:)
  end do

  SAFE_DEALLOCATE_A(tmp)

  call pop_sub()
end subroutine X(derivatives_div)


! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, f, curl, ghost_update)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:,:)    ! f(m%np_part, der%dim) 
  R_TYPE,              intent(out)   :: curl(:,:) ! curl(m%np, der%dim) if dim = 2, curl(m%np, 1) if dim = 1.
  logical, optional,   intent(in)    :: ghost_update

  R_TYPE, allocatable :: tmp(:)
  integer             :: i, np

  call push_sub('derivatives_inc.Xderivatives_div')

  ASSERT(der%dim == 3 .or. der%dim == 2)
  ASSERT(ubound(f,    DIM=1) == der%mesh%np_part)
  ASSERT(ubound(curl, DIM=1) >= der%mesh%np)
  select case(der%dim)
    case(3)
      ASSERT(ubound(curl, DIM=2) == der%dim)
    case(2)
      ASSERT(ubound(curl, DIM=2) == 1)
    case(1)
      write(message(1),'(a)') 'INTERNAL ERROR at Xderivatives_curl: 1D not allowed'
      call write_fatal(1)
  end select

  do i = 1, der%dim
    call X(set_bc)(der, f(:, i))
  end do
  
  SAFE_ALLOCATE(tmp(1:der%mesh%np_part))

  curl(:,:) = R_TOTYPE(M_ZERO)
  np = der%mesh%np

  select case(der%dim)
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

  SAFE_DEALLOCATE_A(tmp)
  call pop_sub()
end subroutine X(derivatives_curl)


! ---------------------------------------------------------
! Set all boundary points in f to zero to implement zero
! boundary conditions for the derivatives, in finite system;
! or set according to periodic boundary conditions.
subroutine X(set_bc)(der, f)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:)

  integer :: bndry_start, bndry_end
  integer :: p
  integer :: iper, ip, ix, iy, iz, dx, dy, dz
#ifdef HAVE_MPI
  integer :: ipart, nreq
  integer, allocatable :: req(:), statuses(:, :)
#endif

  call push_sub('derivatives_inc.Xset_bc')
  call profiling_in(set_bc_prof, 'SET_BC')
  
  ASSERT(ubound(f, DIM=1) >= der%mesh%np_part)

  p = der%mesh%vp%partno
  
  ! The boundary points are at different locations depending on the presence
  ! of ghost points due to domain parallelization.
  if(der%mesh%parallel_in_domains) then
    bndry_start = der%mesh%vp%np_local(p) + der%mesh%vp%np_ghost(p) + 1
    bndry_end   = der%mesh%vp%np_local(p) + der%mesh%vp%np_ghost(p) + der%mesh%vp%np_bndry(p)
  else
    bndry_start = der%mesh%np+1
    bndry_end   = der%mesh%np_part
  end if
    
  if(der%zero_bc) then
    forall(ip = bndry_start:bndry_end) f(ip) = R_TOTYPE(M_ZERO)
  end if

  if(der%mesh%sb%mr_flag) call multires()

  if(der%periodic_bc) then

#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains) then
      call profiling_in(set_bc_comm_prof, 'SET_BC_COMMUNICATION')
      ! get the points from other nodes
      SAFE_ALLOCATE(req(1:2*der%mesh%vp%npart))

      nreq = 0

      do ipart = 1, der%mesh%vp%npart
        if(ipart == p .or. der%mesh%nsend(ipart) == 0) cycle
        nreq = nreq + 1
        call MPI_Isend(f, 1, der%mesh%X(send_type)(ipart), ipart - 1, 3, der%mesh%vp%comm, req(nreq), mpi_err)
      end do

      do ipart = 1, der%mesh%vp%npart
        if(ipart == p .or. der%mesh%nrecv(ipart) == 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(f, 1, der%mesh%X(recv_type)(ipart), ipart - 1, 3, der%mesh%vp%comm, req(nreq), mpi_err)
      end do

      call profiling_count_transfers(sum(der%mesh%nrecv(1:der%mesh%vp%npart) + der%mesh%nrecv(1:der%mesh%vp%npart)), f(1))

      call profiling_out(set_bc_comm_prof)
    end if
#endif

    do iper = 1, der%mesh%nper
      f(der%mesh%per_points(iper)) = f(der%mesh%per_map(iper))
    end do

#ifdef HAVE_MPI
    if(der%mesh%parallel_in_domains) then
      call profiling_in(set_bc_comm_prof)

      SAFE_ALLOCATE(statuses(1:MPI_STATUS_SIZE, 1:nreq))
      call MPI_Waitall(nreq, req, statuses, mpi_err)
      SAFE_DEALLOCATE_A(statuses)
      SAFE_DEALLOCATE_A(req)

      call profiling_out(set_bc_comm_prof)
    end if
#endif

  end if

  call profiling_out(set_bc_prof)
  call pop_sub()

contains 

  subroutine multires()
    integer :: order, nn, ii, jj, kk, i_lev
    FLOAT, allocatable :: pos(:), ww(:)
    integer, allocatable :: posi(:)
    
    call push_sub('derivatives_inc.Xset_bc.multires')

    order = der%mesh%sb%mr_iorder !interpolation order
    
    nn = 2*order

    SAFE_ALLOCATE(ww(1:nn))
    SAFE_ALLOCATE(pos(1:nn))
    SAFE_ALLOCATE(posi(1:nn))

    do ii = 1, order
      posi(ii) = 1 + 2*(ii - 1)
      posi(order + ii) = -posi(ii)
      pos(ii) =  posi(ii)
      pos(order + ii) = -pos(ii)
    end do
    
    call interpolation_coefficients(nn, pos, M_ZERO, ww)

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

        do ii = 1, nn
          do jj = 1, nn
            do kk = 1, nn
              f(ip) = f(ip) + ww(ii)*ww(jj)*ww(kk)*&
                              f(der%mesh%idx%Lxyz_inv(ix + posi(ii)*dx, &
                                                      iy + posi(jj)*dy, &
                                                      iz + posi(kk)*dz))
            end do
          end do
        end do
        
      end if
      
    end do

    SAFE_DEALLOCATE_A(ww)
    SAFE_DEALLOCATE_A(pos)
    SAFE_DEALLOCATE_A(posi)

    call pop_sub()
  end subroutine multires

end subroutine X(set_bc)

! ---------------------------------------------------------
! Set all boundary points in f to zero to implement zero
! boundary conditions for the derivatives, in finite system;
! or set according to periodic boundary conditions.
subroutine X(set_bc_batch)(der, fb)
  type(derivatives_t), intent(in)    :: der
  type(batch_t),       intent(inout) :: fb

  integer :: ip, idim, ist
  integer :: pp, bndry_start, bndry_end
#ifdef HAVE_MPI
  integer :: ipart, nreq
  integer, allocatable :: req(:), statuses(:, :)
#endif

  call push_sub('derivatives_inc.Xset_bc_batch')

  if(der%mesh%sb%mr_flag) then

    do ist = 1, fb%nst
      do idim = 1, fb%dim
        call X(set_bc)(der, fb%states(ist)%X(psi)(:, idim))
      end do
    end do

  else

    call profiling_in(set_bc_prof, 'SET_BC')

    pp = der%mesh%vp%partno

    ! The boundary points are at different locations depending on the presence
    ! of ghost points due to domain parallelization.
    if(der%mesh%parallel_in_domains) then
      bndry_start = der%mesh%vp%np_local(pp) + der%mesh%vp%np_ghost(pp) + 1
      bndry_end   = der%mesh%vp%np_local(pp) + der%mesh%vp%np_ghost(pp) + der%mesh%vp%np_bndry(pp)
    else
      bndry_start = der%mesh%np+1
      bndry_end   = der%mesh%np_part
    end if

    if(der%zero_bc) then
      forall (ist = 1:fb%nst, idim = 1:fb%dim, ip = bndry_start:bndry_end)
        fb%states(ist)%X(psi)(ip, idim) = R_TOTYPE(M_ZERO)
      end forall
    end if

    if(der%periodic_bc) then

#ifdef HAVE_MPI
      if(der%mesh%parallel_in_domains) then
        call profiling_in(set_bc_comm_prof, 'SET_BC_COMMUNICATION')

        ! get the points that are copies from other nodes
        SAFE_ALLOCATE(req(1:2*der%mesh%vp%npart*fb%dim*fb%nst))

        nreq = 0

        do ist = 1, fb%nst
          do idim = 1, fb%dim
            do ipart = 1, der%mesh%vp%npart
              if(ipart == der%mesh%vp%partno .or. der%mesh%nrecv(ipart) == 0) cycle
              nreq = nreq + 1
              call MPI_Irecv(fb%states(ist)%X(psi)(:, idim), 1, der%mesh%X(recv_type)(ipart), ipart - 1, 3, &
                   der%mesh%vp%comm, req(nreq), mpi_err)
            end do
          end do
        end do

        do ist = 1, fb%nst
          do idim = 1, fb%dim
            do ipart = 1, der%mesh%vp%npart
              if(ipart == der%mesh%vp%partno .or. der%mesh%nsend(ipart) == 0) cycle
              nreq = nreq + 1
              call MPI_Isend(fb%states(ist)%X(psi)(:, idim), 1, der%mesh%X(send_type)(ipart), ipart - 1, 3, &
                   der%mesh%vp%comm, req(nreq), mpi_err)
            end do
          end do
        end do
        
        call profiling_count_transfers(&
        sum(der%mesh%nsend(1:der%mesh%vp%npart) + der%mesh%nrecv(1:der%mesh%vp%npart))*fb%dim*fb%nst, &
             R_TOTYPE(M_ONE))

        call profiling_out(set_bc_comm_prof)
      end if
#endif

      forall (ist = 1:fb%nst, idim = 1:fb%dim, ip = 1:der%mesh%nper)
        fb%states(ist)%X(psi)(der%mesh%per_points(ip), idim) = fb%states(ist)%X(psi)(der%mesh%per_map(ip), idim)
      end forall

#ifdef HAVE_MPI
      if(der%mesh%parallel_in_domains) then
        call profiling_in(set_bc_comm_prof)

        SAFE_ALLOCATE(statuses(1:MPI_STATUS_SIZE, 1:nreq))
        call MPI_Waitall(nreq, req, statuses, mpi_err)
        SAFE_DEALLOCATE_A(statuses)
        SAFE_DEALLOCATE_A(req)

        call profiling_out(set_bc_comm_prof)
      end if
#endif

    end if

    call profiling_out(set_bc_prof)

  end if

  call pop_sub()

end subroutine X(set_bc_batch)

! ---------------------------------------------------------
! The action of the angular momentum operator (three spatial components).
! In case of real functions, it does not include the -i prefactor
! (L = -i r ^ nabla).
! ---------------------------------------------------------
subroutine X(f_angular_momentum)(sb, mesh, der, f, lf, ghost_update, set_bc)
  type(simul_box_t),   intent(in)    :: sb
  type(mesh_t),        intent(in)    :: mesh
  type(derivatives_t), intent(inout) :: der
  R_TYPE,              intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,              intent(out)   :: lf(:, :)  ! lf(m%np, 3) in 3D, lf(m%np, 1) in 2D
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x1, x2, x3
  R_TYPE :: factor
  integer :: ip

  call push_sub('derivatives_inc.Xf_angular_momentum')

  ASSERT(sb%dim.ne.1)

  SAFE_ALLOCATE(gf(1:mesh%np, 1:sb%dim))

  if (present(ghost_update) .and. present(set_bc)) &
    call X(derivatives_grad)(der, f, gf, ghost_update, set_bc)
  if (present(ghost_update) .and. .not. present(set_bc)) &
    call X(derivatives_grad)(der, f, gf, ghost_update = ghost_update)
  if (.not. present(ghost_update) .and. present(set_bc)) &
    call X(derivatives_grad)(der, f, gf, set_bc = set_bc)
  if (.not. present(ghost_update) .and. .not. present(set_bc)) &
    call X(derivatives_grad)(der, f, gf)

#if defined(R_TCOMPLEX)
  factor = -M_ZI
#else
  factor = M_ONE
#endif

  select case(sb%dim)
  case(3)
    do ip = 1, mesh%np
      x1 = mesh%x(ip, 1)
      x2 = mesh%x(ip, 2)
      x3 = mesh%x(ip, 3)
      lf(ip, 1) = factor*(x2*gf(ip, 3) - x3*gf(ip, 2))
      lf(ip, 2) = factor*(x3*gf(ip, 1) - x1*gf(ip, 3))
      lf(ip, 3) = factor*(x1*gf(ip, 2) - x2*gf(ip, 1))
    end do

  case(2)
    do ip = 1, mesh%np
      x1 = mesh%x(ip, 1)
      x2 = mesh%x(ip, 2)
      lf(ip, 1) = factor*(x1*gf(ip, 2) - x2*gf(ip, 1))
    end do

  end select

  SAFE_DEALLOCATE_A(gf)
  call pop_sub()
end subroutine X(f_angular_momentum)

! ---------------------------------------------------------
! Square of the angular momentum L. This has to be very much improved if
! accuracy is needed.
! ---------------------------------------------------------
subroutine X(f_l2)(sb, m, der, f, l2f, ghost_update)
  type(simul_box_t),   intent(in)    :: sb
  type(mesh_t),        intent(in)    :: m
  type(derivatives_t), intent(inout) :: der
  R_TYPE,              intent(inout) :: f(:)   ! f(1:m%np_part)
  R_TYPE,              intent(out)   :: l2f(:)
  logical, optional,   intent(in)    :: ghost_update

  R_TYPE, allocatable :: gf(:, :), ggf(:, :, :)
  integer :: j

  call push_sub('derivatives_inc.Xf_l2')

  ASSERT(sb%dim == 2 .or. sb%dim == 3)

  l2f = R_TOTYPE(M_ZERO)

  select case(sb%dim)
  case(3)
    SAFE_ALLOCATE( gf(1:m%np_part, 1:3))
    SAFE_ALLOCATE(ggf(1:m%np_part, 1:3, 1:3))

    if (present(ghost_update)) then
       call X(f_angular_momentum)(sb, m, der, f, gf, ghost_update)
    else
       call X(f_angular_momentum)(sb, m, der, f, gf)
    endif

    do j = 1, 3
      if (present(ghost_update)) then
        call X(f_angular_momentum)(sb, m, der, gf(:,j), ggf(:,:,j), ghost_update)
      else
        call X(f_angular_momentum)(sb, m, der, gf(:,j), ggf(:,:,j))
      endif

    end do

    do j = 1, sb%dim
      l2f(1:m%np) = l2f(1:m%np) + ggf(1:m%np, j, j)
    end do

  case(2)
    SAFE_ALLOCATE( gf(1:m%np_part, 1:1))
    SAFE_ALLOCATE(ggf(1:m%np_part, 1:1, 1:1))

    if (present(ghost_update)) then
      call X(f_angular_momentum)(sb, m, der, f, gf, ghost_update)
      call X(f_angular_momentum)(sb, m, der, gf(:, 1), ggf(:, :, 1), ghost_update)
    else
      call X(f_angular_momentum)(sb, m, der, f, gf)
      call X(f_angular_momentum)(sb, m, der, gf(:, 1), ggf(:, :, 1))
    endif

    l2f(1:m%np) = ggf(1:m%np, 1, 1)
  end select


  ! In case of real functions, since the angular momentum calculations
  ! lack a (-i) prefactor, we must add a (-1) factor
#if defined(R_TREAL)
  l2f = - l2f
#endif

  SAFE_DEALLOCATE_A(gf)
  SAFE_DEALLOCATE_A(ggf)
  call pop_sub()
end subroutine X(f_l2)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
