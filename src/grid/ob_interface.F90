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
    interface_apply_sym_op, &
    interface_apply_op,     &
    interface_t,            &
    interface_init,         &
    interface_end,          &
    interface_write_info,   &
    get_intf_wf,            &
    put_intf_wf,            &
    member_of_interface

  ! Describes the points belonging to the left and right interface in an
  ! open boundaries calculation.
  ! (Important: the point numbers of the interface are dense and
  ! sorted in increasing order.)
  type interface_t
    integer          :: extent
    integer          :: np      ! number of interface points
    integer          :: np_uc   ! number of unit cell points
    integer          :: nblocks ! nblocks = np_uc / np number of interfaces in a unit cell
    integer          :: il      ! which lead (1..NLEADS)
    integer, pointer :: index(:)       ! (np)
    integer          :: index_range(2)
    logical          :: offdiag_invertible
  end type interface_t

  type lead_t
    CMPLX, pointer :: h_diag(:, :, :)      ! Diagonal block of the lead Hamiltonian.
    CMPLX, pointer :: h_offdiag(:, :)      ! Offdiagonal block of the lead Hamiltonian.
    FLOAT, pointer :: vks(:, :)            ! (np, nspin) Kohn-Sham potential of the leads.
    FLOAT, pointer :: vhartree(:)          ! (np) Hartree potential of the leads.
  end type lead_t

contains

  ! ---------------------------------------------------------
  ! Calculate the member points of the interface region.
  subroutine interface_init(m, sb, der_discr, intf, il, extent)
    type(mesh_t),        intent(in)  :: m
    type(simul_box_t),   intent(in)  :: sb
    type(derivatives_t), intent(in)  :: der_discr
    type(interface_t),   intent(out) :: intf
    integer,             intent(in)  :: il
    integer, optional,   intent(in)  :: extent

    integer :: from(MAX_DIM), to(MAX_DIM), ll(MAX_DIM), dir, tdir

    call push_sub('ob_interface.interface_init')

    intf%il = il

    tdir = (il+1)/2
    
    if (present(extent)) then
      write(*,*) 'stencil_extent(tdir)', tdir, derivatives_stencil_extent(der_discr, tdir)
      write(*,*) 'extent', extent
      ASSERT(derivatives_stencil_extent(der_discr, tdir).le.extent)
      intf%extent = extent
    else
      write(*,*) 'stencil_extent(tdir)', tdir, derivatives_stencil_extent(der_discr, tdir)
      write(*,*) 'lead_unit_cell_extent(sb, il)', il, lead_unit_cell_extent(sb, il)
      ASSERT(derivatives_stencil_extent(der_discr, tdir).le.lead_unit_cell_extent(sb, il))
      intf%extent = maxval((/derivatives_stencil_extent(der_discr, tdir), lead_unit_cell_extent(sb, il)/))
      ! if the lead potential has no dependence in transport direction
      ! then reduce size of unit cell extent to interface extent
      !if(stencil_extent(der_discr, tdir).eq.1) then
      !  intf%extent = intf%extent - 1
      !end if
    endif
    if(intf%extent.eq.derivatives_stencil_extent(der_discr, tdir)) then
      intf%offdiag_invertible = .true.
    else
      intf%offdiag_invertible = .false.
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
    ll(:) = m%idx%ll(:)
    ! the interface region has only intf%extent points in its normal direction
    ll(tdir) = intf%extent

    intf%np = product(ll(:))
    SAFE_ALLOCATE(intf%index(1:intf%np))

    ! the point where we start
    from(:) = m%idx%nr( mod(il+1,2)+1, :) + dir*m%idx%enlarge
    ! the point were we end
    ! we are 1 point too far in every direction, so go back
    to(:)   = from(:) + dir*(ll(:) - 1)


    call mesh_subset_indices(m, from, to, intf%index)

    ! Extract begin and end point numbers of interface region.
    ! Sort the array first to make translations between interface point
    ! number and global point number faster.
    call sort(intf%index)
    intf%index_range(1) = intf%index(1)
    intf%index_range(2) = intf%index(intf%np)

    call pop_sub()
  end subroutine interface_init


  ! ---------------------------------------------------------
  ! Checks if point number idx is an interface point of interface.
  logical function member_of_interface(idx, intf, index)
    integer,           intent(in) :: idx
    type(interface_t), intent(in) :: intf
    integer,           intent(out) :: index ! index in interface

    integer :: ii

    call push_sub('ob_interface.member_of_interface')

    ! first test if index is between min and max
    member_of_interface =              &
      idx.ge.intf%index_range(1) .and. &
      idx.le.intf%index_range(2)

    index = 0

    ! if so check exactly
    ! (not the fastest way, but works if we have non-connected points)
    if (member_of_interface) then
      member_of_interface = .false.
      do ii=1, intf%np
        if (intf%index(ii).eq.idx) then
          member_of_interface = .true.
          index = ii
          exit
        end if
      end do
    end if

    call pop_sub()
  end function member_of_interface


  ! ---------------------------------------------------------
  ! Get wave function points of interface intf.
  ! intf_wf must be of dimension intf%np.
  subroutine get_intf_wf(intf, zpsi, intf_wf)
    type(interface_t), intent(in)  :: intf
    CMPLX,             intent(in)  :: zpsi(:)
    CMPLX,             intent(out) :: intf_wf(:)

    integer :: ii

    call push_sub('ob_interface.get_intf_wf')

    do ii=1, intf%np
      intf_wf(ii) = zpsi(intf%index(ii))
    end do

    call pop_sub()
  end subroutine get_intf_wf


  ! ---------------------------------------------------------
  ! Put wave function points of interface intf.
  ! intf_wf must be of dimension intf%np.
  subroutine put_intf_wf(intf, intf_wf, zpsi)
    type(interface_t), intent(in)     :: intf
    CMPLX,             intent(in)     :: intf_wf(:)
    CMPLX,             intent(inout)  :: zpsi(:)

    integer :: ii

    call push_sub('ob_interface.put_intf_wf')

    do ii=1, intf%np
      zpsi(intf%index(ii)) = intf_wf(ii)
    end do

    call pop_sub()
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

    call push_sub('ob_interface.interface_apply_op')

    SAFE_ALLOCATE(intf_wf(1:intf%np))
    SAFE_ALLOCATE(op_intf_wf(1:intf%np))

    call get_intf_wf(intf, wf, intf_wf)
    call get_intf_wf(intf, res, op_intf_wf)
    call lalg_gemv(intf%np, intf%np, alpha, op, intf_wf, M_z1, op_intf_wf)
    call put_intf_wf(intf, op_intf_wf, res)

    SAFE_DEALLOCATE_A(intf_wf)
    SAFE_DEALLOCATE_A(op_intf_wf)

    call pop_sub()
  end subroutine interface_apply_op


  ! ---------------------------------------------------------
  ! Apply an np x np symmetric operator op to the interface region of the wavefunction.
  ! np is the number of points in the interface: res <- res + alpha op wf
  subroutine interface_apply_sym_op(intf, alpha, op, wf, res)
    type(interface_t), intent(in)    :: intf
    CMPLX,             intent(in)    :: alpha
    CMPLX,             intent(in)    :: op(:, :)
    CMPLX,             intent(in)    :: wf(:)
    CMPLX,             intent(inout) :: res(:)
    
    CMPLX, allocatable :: intf_wf(:), op_intf_wf(:)

    call push_sub('ob_interface.interface_apply_sym_op')

    SAFE_ALLOCATE(intf_wf(1:intf%np))
    SAFE_ALLOCATE(op_intf_wf(1:intf%np))

    call get_intf_wf(intf, wf, intf_wf)
    call get_intf_wf(intf, res, op_intf_wf)
    call lalg_symv(intf%np, alpha, op, intf_wf, M_z1, op_intf_wf)
    call put_intf_wf(intf, op_intf_wf, res)

    SAFE_DEALLOCATE_A(intf_wf)
    SAFE_DEALLOCATE_A(op_intf_wf)

    call pop_sub()
  end subroutine interface_apply_sym_op


  ! ---------------------------------------------------------
  ! Write number of interface points.
  subroutine interface_write_info(intf, iunit)
    type(interface_t), intent(in) :: intf
    integer,           intent(in) :: iunit

    call push_sub('ob_interface.interface_write_info')

    write(message(1), '(a,i6)') 'Number of points in '//LEAD_NAME(intf%il)// &
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
