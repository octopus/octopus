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
  public ::                 &
    lead_t,                 &
    interface_apply_op,     &
    interface_t,            &
    interface_init,         &
    interface_end,          &
    interface_messages_info,   &
    get_intf_wf,            &
    put_intf_wf,            &
    member_of_interface

  ! Describes the points belonging to the left and right interface in an
  ! open boundaries calculation.
  ! (Important: the point numbers of the interface are dense and
  ! sorted in increasing order.)
  type interface_t
    integer          :: extent_uc  ! extent of the lead unit cell
    integer          :: np_intf    ! number of interface points
    integer          :: np_uc      ! number of unit cell points
    integer          :: np_part_uc ! number of unit cell points including all points
    integer          :: nblocks    ! number of interfaces in a unit cell
    integer          :: il         ! which lead (1..NLEADS)
    integer, pointer :: index(:)   ! (np_uc)
    integer          :: index_range(2) ! lowest and highest index
    logical          :: reducible  ! is the lead unit cell a integer multiple of the interface
  end type interface_t

  type lead_t
    integer        :: np              ! number of points in the lead unit cell
    integer        :: np_part         ! as np including ghost and boundary points
    CMPLX, pointer :: h_diag(:, :, :) ! Diagonal block of the lead Hamiltonian.
    CMPLX, pointer :: h_offdiag(:, :) ! Offdiagonal block of the lead Hamiltonian.
    FLOAT, pointer :: vks(:, :)       ! (np_uc, nspin) Kohn-Sham potential of the leads.
    FLOAT, pointer :: vhartree(:)     ! (np_uc) Hartree potential of the leads.
!    FLOAT, pointer :: A_static(:,:)   ! static vector potential
!    CMPLX, pointer :: grad_diag(:, :, :)    ! diagonal block of the gradient operator
!    CMPLX, pointer :: grad_offdiag(:, :, :) ! offdiagonal block of the gradient operator
  end type lead_t

contains

  ! ---------------------------------------------------------
  ! Calculate the member points of the interface region.
  subroutine interface_init(der, intf, il, lsize, extent_uc)
    type(derivatives_t), intent(in)  :: der
    type(interface_t),   intent(out) :: intf
    integer,             intent(in)  :: il
    FLOAT,               intent(in)  :: lsize(:)
    integer, optional,   intent(in)  :: extent_uc ! new reduced extent of the unit cell

    integer :: from(MAX_DIM), to(MAX_DIM), ll(MAX_DIM), dir, tdir

    PUSH_SUB(interface_init)

    ASSERT(associated(der%mesh))

    intf%il = il
    tdir = (il+1)/2

    intf%np_part_uc = 0
    if (present(extent_uc)) then
      ! if the lead potential has no dependence in transport direction
      ! then the size of unit cell extent will be reduced to interface extent
      ! this happens within hamiltonian->init_lead_h which calls this subroutine
      intf%extent_uc = extent_uc
    else
      intf%extent_uc = nint(2*lsize(tdir)/der%mesh%spacing(tdir))
    endif

    intf%reducible = mod(intf%extent_uc, derivatives_stencil_extent(der, tdir)).eq.0
    if(.not.intf%reducible) then
      intf%nblocks = 1
    else
      intf%nblocks = intf%extent_uc / derivatives_stencil_extent(der, tdir)
    end if
    ! Extract submeshes for the interface regions.
    ! 
    ! In 2D for the left lead.
    !   +------to------------ ...
    !   |       |
    !   |       |
    ! from------+------------ ...
    ! valid for every lead
    ! directions: LEFT -> RIGHT, BOTTOM -> TOP, REAR -> FRONT
    
    ! direction: +1 if from L->R and -1 if R->L (other leads accordingly)
    dir = (-1)**(il+1)
    ll(:) = der%mesh%idx%ll(:)
    ! the interface region has only intf%extent points in its normal direction
    ll(tdir) =  derivatives_stencil_extent(der, tdir)
    
    ! if bottom, top, front or back lead then reduce x-extention by 2 points
    ! (covered already by left and right lead)
    if(tdir > 1) ll(1) = ll(1)-2
    ! if front or back lead then reduce also y-extention by 2 points
    ! (covered already by bottom and top lead)
    if(tdir > 2) ll(2) = ll(2)-2

    intf%np_intf = product(ll(:))

    ll(tdir) =  intf%extent_uc
    intf%np_uc = product(ll(:))

    SAFE_ALLOCATE(intf%index(1:intf%np_uc))

    ! the point where we start
    from(:) = der%mesh%idx%nr( mod(il+1,2)+1, :) + dir*der%mesh%idx%enlarge

    ! shift the starting point by 1 so that it is centralized again
    if(tdir > 1) from(1) = from(1) + dir
    if(tdir > 2) from(2) = from(2) + dir
    ! the point were we end
    ! we are 1 point too far in every direction, so go back
    to(:)   = from(:) + dir*(ll(:) - 1)

    call mesh_subset_indices(der%mesh, from, to, intf%index)

    ! Extract begin and end point numbers of interface region.
    ! Sort the array first to make translations between interface point
    ! number and global point number faster.
    call sort(intf%index)
    intf%index_range(1) = intf%index(1)
    intf%index_range(2) = intf%index(intf%np_uc)

    POP_SUB(interface_init)
  end subroutine interface_init


  ! ---------------------------------------------------------
  ! Checks if point number idx is an interface point of interface.
  logical function member_of_interface(idx, intf, index)
    integer,           intent(in)  :: idx
    type(interface_t), intent(in)  :: intf
    integer,           intent(out) :: index ! index in interface

    integer :: ii

    PUSH_SUB(member_of_interface)

    ! first test if index is between min and max
    member_of_interface = (idx >= intf%index_range(1) .and. idx <= intf%index_range(2))

    index = 0

    ! if so check exactly
    ! (not the fastest way, but works if we have non-connected points)
    if (member_of_interface) then
      member_of_interface = .false.
      do ii=1, intf%np_uc
        if (intf%index(ii).eq.idx) then
          member_of_interface = .true.
          index = ii
          exit
        end if
      end do
    end if

    POP_SUB(member_of_interface)
  end function member_of_interface


  ! ---------------------------------------------------------
  ! Get wavefunction points of interface intf.
  ! intf_wf must be of dimension intf%np_uc.
  subroutine get_intf_wf(intf, zpsi, intf_wf)
    type(interface_t), intent(in)  :: intf
    CMPLX,             intent(in)  :: zpsi(:)
    CMPLX,             intent(out) :: intf_wf(:)

    PUSH_SUB(get_intf_wf)

    intf_wf(1:intf%np_uc) = zpsi(intf%index(1:intf%np_uc))

    POP_SUB(get_intf_wf)
  end subroutine get_intf_wf


  ! ---------------------------------------------------------
  ! Put wavefunction points of interface intf.
  ! intf_wf must be of dimension intf%np_uc.
  subroutine put_intf_wf(intf, intf_wf, zpsi)
    type(interface_t), intent(in)     :: intf
    CMPLX,             intent(in)     :: intf_wf(:)
    CMPLX,             intent(inout)  :: zpsi(:)

    PUSH_SUB(put_intf_wf)
    
    ! FIXME: this will probably fail for MeshBlockSize > 1
    zpsi(intf%index(1:intf%np_uc)) = intf_wf(1:intf%np_uc)

    POP_SUB(put_intf_wf)
  end subroutine put_intf_wf


  ! ---------------------------------------------------------
  ! Apply an np x np operator op to the interface region of the wavefunction.
  ! np is the number of points in the interface: wf <- wf + alpha op wf
  subroutine interface_apply_op(intf, alpha, op, wf, res)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: alpha
    CMPLX,             intent(in)    :: op(:, :)
    CMPLX,             intent(in)    :: wf(:)
    CMPLX,             intent(inout) :: res(:)
    
    CMPLX, allocatable :: intf_wf(:), op_intf_wf(:)

    PUSH_SUB(interface_apply_op)

    SAFE_ALLOCATE(intf_wf(1:intf%np_uc))
    SAFE_ALLOCATE(op_intf_wf(1:intf%np_uc))

    call get_intf_wf(intf, wf, intf_wf)
    call get_intf_wf(intf, res, op_intf_wf)
    call lalg_gemv(intf%np_uc, intf%np_uc, alpha, op, intf_wf, M_z1, op_intf_wf)
    call put_intf_wf(intf, op_intf_wf, res)

    SAFE_DEALLOCATE_A(intf_wf)
    SAFE_DEALLOCATE_A(op_intf_wf)

    POP_SUB(interface_apply_op)
  end subroutine interface_apply_op


  ! ---------------------------------------------------------
  ! Write number of interface points.
  subroutine interface_messages_info(intf, iunit)
    type(interface_t), intent(in) :: intf
    integer,           intent(in) :: iunit

    PUSH_SUB(interface_messages_info)

    write(message(1), '(a,i6)') 'Number of points in '//LEAD_NAME(intf%il)// &
      ' interface: ', intf%np_uc
    call messages_info(1, iunit)

    POP_SUB(interface_messages_info)
  end subroutine interface_messages_info


  ! ---------------------------------------------------------
  ! Free memory.
  subroutine interface_end(intf)
    type(interface_t), intent(inout) :: intf

    PUSH_SUB(intface_end)

    SAFE_DEALLOCATE_P(intf%index)

    POP_SUB(intface_end)
  end subroutine interface_end

end module ob_interface_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
