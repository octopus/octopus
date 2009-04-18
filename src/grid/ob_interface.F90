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
!! $Id: td_transport.F90 3030 2007-06-25 16:45:05Z marques $

! Management of the mesh points in the interface regions
! for open boundary calculations.

#include "global.h"

module ob_interface_m
  use derivatives_m
  use global_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m

  implicit none

  private
  public ::               &
    interface_apply_sym_op, &
    interface_apply_op,     &
    interface_t,          &
    interface_init,       &
    interface_end,        &
    interface_index,      &
    interface_write_info, &
    get_intf_wf,          &
    put_intf_wf,          &
    member_of_interface

  ! Describes the points belonging to the left and right interface in an
  ! open boundaries calculation.
  ! (Important: the point numbers of the interface are dense and
  ! sorted in increasing order.)
  type interface_t
    integer          :: extent
    integer          :: np      ! number of interface points
    integer          :: np_uc   ! number of points in the unit cell
    integer, pointer :: index(:)       ! (np)
    integer          :: index_range(2)
    logical          :: offdiag_invertible
  end type interface_t

contains

  ! ---------------------------------------------------------
  ! Calculate the member points of the interface region.
  subroutine interface_init(m, sb, der_discr, intf, il)
    type(mesh_t),        intent(in)  :: m
    type(simul_box_t),   intent(in)  :: sb
    type(derivatives_t), intent(in)  :: der_discr
    type(interface_t),   intent(out) :: intf
    integer,             intent(in)  :: il

    logical :: ok
    integer :: i, from(MAX_DIM), to(MAX_DIM), dir, lr

    call push_sub('ob_interface.interface_init')

    intf%extent = maxval((/stencil_extent(der_discr, TRANS_DIR), lead_unit_cell_extent(sb, il)/))
    ! FIXME: not sure if it holds for lead-potentials with dependence in transport direction
    intf%extent = intf%extent - 1
    if(intf%extent.eq.stencil_extent(der_discr, TRANS_DIR)) then
      intf%offdiag_invertible = .true.
    else
      intf%offdiag_invertible = .false.
    end if

    intf%np = intf%extent*m%idx%ll(2)*m%idx%ll(3)
    intf%np_uc = (intf%extent + 1)*m%idx%ll(2)*m%idx%ll(3)

    ALLOCATE(intf%index(intf%np), intf%np)

    ! Extract submeshes for the interface regions.
    ! 
    ! In 2D for the left lead.
    !   +------to------------ ...
    !   |       |
    !   |       |
    ! from------+------------ ...
    if(il.eq.LEFT) then
      dir = 1
      lr  = 1
    else
      lr  = 2
      dir = -1
    end if
    from = m%idx%nr(lr, :) + dir*m%idx%enlarge
    to(1:3) = from(1:3) + dir*(/intf%extent, m%idx%ll(2), m%idx%ll(3)/) - dir

    call mesh_subset_indices(m, from, to, intf%index)

    ! Extract begin and end point numbers of interface region.
    ! Sort the array first to make translations between interface point
    ! number and global point number faster.
    call sort(intf%index)
    intf%index_range(1) = intf%index(1)
    intf%index_range(2) = intf%index(intf%np)

    ! In debug mode, check if x_intface is dense, i. e. no index
    ! between min(x_index) and max(x_index) is omitted.
    ! As transport is along x-axis, which is the index running slowest,
    ! this should be the case.
    if(in_debug_mode) then
      ok = .true.
      do i = intf%index_range(1), intf%index_range(2)
        if(intf%index(i-intf%index_range(1)+1).ne.i) then
          ok = .false.
        end if
      end do
      if(.not.ok) then
        message(1) = 'Failed assertion:'
        message(2) = 'ob_interface.interface_init: the interface region'
        message(3) = 'has holes.'
        call write_fatal(3)
      end if
    end if

    call pop_sub()
  end subroutine interface_init


  ! ---------------------------------------------------------
  ! Return the interface point number of the global point idx.
  ! For this to work, the intface%index(:, il) array has to
  ! be sorted in increasing order.
  integer function interface_index(idx, intf)
    integer,           intent(in) :: idx
    type(interface_t), intent(in) :: intf

    call push_sub('ob_interface.interface_index')

    interface_index = idx - intf%index_range(1) + 1

    call pop_sub()
  end function interface_index
  

  ! ---------------------------------------------------------
  ! Checks if point number idx is an interface point of interface
  ! number il.
  logical function member_of_interface(idx, intf)
    integer,           intent(in) :: idx
    type(interface_t), intent(in) :: intf

    call push_sub('ob_interface.member_of_interface')

    member_of_interface =              &
      idx.ge.intf%index_range(1) .and. &
      idx.le.intf%index_range(2)

    call pop_sub()
  end function member_of_interface


  ! ---------------------------------------------------------
  ! Get wave function points of interface intf.
  ! intf_wf must be of dimension intf%np.
  subroutine get_intf_wf(intf, zpsi, intf_wf)
    type(interface_t), intent(in)  :: intf
    CMPLX,             intent(in)  :: zpsi(:)
    CMPLX,             intent(out) :: intf_wf(:)

    call push_sub('ob_interface.get_intf_wf')

    intf_wf(1:intf%np) = zpsi(intf%index_range(1):intf%index_range(2))

    call pop_sub()
  end subroutine get_intf_wf


  ! ---------------------------------------------------------
  ! Put wave function points of interface intf.
  ! intf_wf must be of dimension intf%np.
  subroutine put_intf_wf(intf, intf_wf, zpsi)
    type(interface_t), intent(in)     :: intf
    CMPLX,             intent(in)     :: intf_wf(:)
    CMPLX,             intent(inout)  :: zpsi(:)

    call push_sub('ob_interface.put_intf_wf')

    zpsi(intf%index_range(1):intf%index_range(2)) = intf_wf(1:intf%np)

    call pop_sub()
  end subroutine put_intf_wf


  ! ---------------------------------------------------------
  ! Apply an np x np operator op to the interface region of the wavefunction.
  ! np is the number of points in the interface: wf <- wf + alpha op wf
  subroutine interface_apply_op(intf, alpha, op, wf)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: alpha
    CMPLX,             intent(in)    :: op(:, :)
    CMPLX,             intent(inout) :: wf(:)
    
    CMPLX, allocatable :: intf_wf(:), op_intf_wf(:)

    call push_sub('ob_interface.interface_apply_op')

    ALLOCATE(intf_wf(intf%np), intf%np)
    ALLOCATE(op_intf_wf(intf%np), intf%np)

    call get_intf_wf(intf, wf, intf_wf)
    op_intf_wf(1:intf%np) = intf_wf(1:intf%np)
    call lalg_gemv(intf%np, intf%np, alpha, op, intf_wf, M_z1, op_intf_wf)
    call put_intf_wf(intf, op_intf_wf, wf)

    SAFE_DEALLOCATE_A(intf_wf)
    SAFE_DEALLOCATE_A(op_intf_wf)

    call pop_sub()
  end subroutine interface_apply_op


  ! ---------------------------------------------------------
  ! Apply an np x np symmetric operator op to the interface region of the wavefunction.
  ! np is the number of points in the interface: res <- beta res + alpha op wf
  subroutine interface_apply_sym_op(intf, alpha, op, wf, beta, res)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: alpha
    CMPLX,             intent(in)    :: op(:, :)
    CMPLX,             intent(in)    :: wf(:)
    CMPLX,             intent(in)    :: beta
    CMPLX,             intent(inout) :: res(:)
    
    CMPLX, allocatable :: intf_wf(:), op_intf_wf(:)

    call push_sub('ob_interface.interface_apply_sym_op')

    ALLOCATE(intf_wf(intf%np), intf%np)
    ALLOCATE(op_intf_wf(intf%np), intf%np)

    call get_intf_wf(intf, wf, intf_wf)
    call get_intf_wf(intf, res, op_intf_wf)
    call lalg_symv(intf%np, alpha, op, intf_wf, beta, op_intf_wf)
    call put_intf_wf(intf, op_intf_wf, res)

    SAFE_DEALLOCATE_A(intf_wf)
    SAFE_DEALLOCATE_A(op_intf_wf)

    call pop_sub()
  end subroutine interface_apply_sym_op


  ! ---------------------------------------------------------
  ! Write number of interface points.
  subroutine interface_write_info(intf, il, iunit)
    type(interface_t), intent(in) :: intf
    integer,           intent(in) :: il
    integer,           intent(in) ::iunit

    call push_sub('ob_interface.interface_write_info')

    write(message(1), '(a,i6)') 'Number of points in '//LEAD_NAME(il)// &
      ' interface: ', intf%np
    call write_info(1, iunit)

    call pop_sub()
  end subroutine interface_write_info


  ! ---------------------------------------------------------
  ! Free memory.
  subroutine interface_end(intf)
    type(interface_t), intent(inout) :: intf

    call push_sub('ob_interface.intface_end')

    SAFE_DEALLOCATE_P(intf%index)

    call pop_sub()
  end subroutine interface_end
end module ob_interface_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
