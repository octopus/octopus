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

! ---------------------------------------------------------
subroutine X(derivatives_lapl_start)(der, handle, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  type(pv_handle_t),         intent(inout) :: handle
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_, update

  call push_sub('derivatives_inc.Xderivatives_lapl_start')

  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(lapl, DIM=1) >= der%m%np)

  set_bc_ = .true.
  if(present(set_bc)) set_bc_ = set_bc
  if(set_bc_) call X(set_bc)(der, f)

#ifdef HAVE_LIBNBC
  
  update = .true.
  if(present(ghost_update)) update = ghost_update
  
  if(der%overlap .and. der%m%parallel_in_domains .and. update) then
    call X(vec_ighost_update)(der%m%vp, f, handle)
    call X(nl_operator_operate)(der%lapl, f, lapl, ghost_update = .false., points = OP_INNER)
  end if

#endif

  call pop_sub()
end subroutine X(derivatives_lapl_start)

! --------------------------------------------------------
! Call NBC_Test to give MPI a chance to process pending messages.
! This improves the obverlap but is only necessary due to deficiencies
! in MPI implementations.
! ---------------------------------------------------------
subroutine X(derivatives_lapl_keep_going)(der, handle, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  type(pv_handle_t),         intent(inout) :: handle
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(inout) :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_, update

  call push_sub('derivatives_inc.Xderivatives_lapl_keep_going')

#ifdef HAVE_LIBNBC
  
  update = .true.
  if(present(ghost_update)) update = ghost_update
  
  if(der%overlap .and. der%m%parallel_in_domains .and. update) then
    call pv_handle_test(handle)
  end if

#endif

  call pop_sub()
end subroutine X(derivatives_lapl_keep_going)

! ---------------------------------------------------------
subroutine X(derivatives_lapl_finish)(der, handle, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  type(pv_handle_t),         intent(inout) :: handle
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(inout) :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_, update

  call push_sub('derivatives_inc.Xderivatives_lapl_finish')

#ifdef HAVE_LIBNBC
  update = .true.
  if(present(ghost_update)) update = ghost_update

  if(der%overlap .and. der%m%parallel_in_domains .and. update) then

    call pv_handle_wait(handle)
    call X(nl_operator_operate) (der%lapl, f, lapl, ghost_update = .false., points = OP_OUTER)
    
    call pop_sub()
    return
  end if
#endif

  call X(nl_operator_operate) (der%lapl, f, lapl, ghost_update = ghost_update)
  call pop_sub()
end subroutine X(derivatives_lapl_finish)

! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, f, lapl, ghost_update, set_bc)
  type(der_discr_t),         intent(in)    :: der
  R_TYPE,                    intent(inout) :: f(:)     ! f(m%np_part)
  R_TYPE,                    intent(out)   :: lapl(:)  ! lapl(m%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  type(pv_handle_t) :: handle

  call push_sub('derivatives_inc.Xderivatives_lapl')

  call pv_handle_init(handle)
  call X(derivatives_lapl_start) (der, handle, f, lapl, ghost_update, set_bc)
  call X(derivatives_lapl_finish)(der, handle, f, lapl, ghost_update, set_bc)
  call pv_handle_end(handle)
        
  call pop_sub()
end subroutine X(derivatives_lapl)

! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, f, grad, ghost_update)
  type(der_discr_t), intent(in)    :: der
  R_TYPE,            intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,            intent(out)   :: grad(:,:)  ! grad(m%np, m%sb%dim)
  logical, optional, intent(in)    :: ghost_update

  integer :: i

  call push_sub('derivatives_inc.Xderivatives_grad')
  
  ASSERT(ubound(f,    DIM=1) == der%m%np_part)
  ASSERT(ubound(grad, DIM=1) >= der%m%np)
  ASSERT(ubound(grad, DIM=2) >= der%m%sb%dim)

  call X(set_bc)(der, f)

  do i = 1, der%m%sb%dim
    call X(nl_operator_operate) (der%grad(i), f, grad(:,i), ghost_update=ghost_update)
  end do

  call pop_sub()
end subroutine X(derivatives_grad)

! ---------------------------------------------------------
subroutine X(derivatives_oper)(op, der, f, opf, ghost_update)
  type(nl_operator_t), intent(in)    :: op
  type(der_discr_t),   intent(in)    :: der
  R_TYPE,              intent(inout) :: f(:)       ! f(m%np_part)
  R_TYPE,              intent(out)   :: opf(:)     ! opf(m%np)
  logical, optional,   intent(in)    :: ghost_update

  call push_sub('derivatives_inc.Xderivatives_grad')
  
  ASSERT(ubound(f,   DIM=1) == der%m%np_part)
  ASSERT(ubound(opf, DIM=1) >= der%m%np)

  call X(set_bc)(der, f)
  call X(nl_operator_operate) (op, f, opf, ghost_update=ghost_update)

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

  call push_sub('derivatives_inc.Xset_bc')
  call profiling_in(set_bc_prof, 'SET_BC')
  
  if(der%zero_bc) then
    
    p = der%m%vp%partno
    
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
    !$omp parallel do
    do iper = 1, der%m%nper
      f(der%m%per_points(iper)) = f(der%m%per_map(iper))
    end do
    !$omp end parallel do
  end if

  call profiling_out(set_bc_prof)
  call pop_sub()

end subroutine X(set_bc)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
