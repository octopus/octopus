!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 M. Oliveira
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

#include "global.h"

module multiresolution_oct_m
  use global_oct_m
  use io_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                     &
    interp_t,                   &
    multiresolution_t,          &
    multiresolution_init,       &
    multiresolution_end,        &
    multiresolution_copy,       &
    multiresolution_use,        &
    multiresolution_dump,       &
    multiresolution_load

  type :: interp_t
    ! Components are public by default
    integer              :: nn, order  !< interpolation points and order
    FLOAT,   allocatable :: ww(:)      !< weights
    integer, allocatable :: posi(:)    !< positions
  end type interp_t

  type :: multiresolution_t
    ! Components are public by default
    type(interp_t)     :: interp          !< interpolation points
    integer,   private :: num_areas = 0   !< number of multiresolution areas
    integer            :: num_radii       !< number of radii (resolution borders)
    FLOAT, allocatable :: radius(:)       !< radius of the high-resolution area
    FLOAT, allocatable :: center(:)       !< central point
  end type multiresolution_t

contains

  !--------------------------------------------------------------
  subroutine multiresolution_init(mr, namespace, space)
    type(multiresolution_t),             intent(inout) :: mr
    type(namespace_t),                   intent(in)    :: namespace
    type(space_t),                       intent(in)    :: space

    integer              :: idir, irad, order
    type(block_t)        :: blk

    PUSH_SUB(multiresolution_init)

    !%Variable MultiResolutionArea
    !%Type block
    !%Section Mesh
    !%Description
    !% (Experimental) Multiresolution regions are set with this
    !% parameter. The first three numbers define the central
    !% point of the region, and the following ones set
    !% the radii where resolution changes (measured from the
    !% central point).
    !% NOTE: currently, only one area can be set up, and only works in 3D, and in serial.
    !%End
    if(parse_block(namespace, 'MultiResolutionArea', blk) == 0) then

      call messages_experimental('Multi-resolution')

      ! number of areas
      mr%num_areas = parse_block_n(blk)

      ! number of radii
      mr%num_radii = parse_block_cols(blk, 0) - space%dim

      ! the central point
      SAFE_ALLOCATE(mr%center(1:space%dim))
      do idir = 1, space%dim
        call parse_block_float(blk, 0, idir - 1, mr%center(idir))
      end do

      if (mr%num_areas /= 1) call messages_input_error(namespace, 'MultiResolutionArea')

      ! the radii
      SAFE_ALLOCATE(mr%radius(1:mr%num_radii))
      do irad = 1, mr%num_radii
        call parse_block_float(blk, 0, space%dim + irad - 1, mr%radius(irad))
        mr%radius(irad) = units_to_atomic(units_inp%length, mr%radius(irad))
      end do
      call parse_block_end(blk)

      ! Create interpolation points (posi) and weights (ww)

      !%Variable MultiResolutionInterpolationOrder
      !%Type integer
      !%Default 5
      !%Section Mesh
      !%Description
      !% The interpolation order in the multiresolution approach (with <tt>MultiResolutionArea</tt>).
      !%End
      call messages_obsolete_variable(namespace, 'MR_InterpolationOrder', 'MultiResolutionInterpolationOrder')
      call parse_variable(namespace, 'MultiResolutionInterpolationOrder', 5, order)
      if (order <= 0) then
        message(1) = "The value for MultiResolutionInterpolationOrder must be > 0."
        call messages_fatal(1, namespace=namespace)
      end if
      call interp_init(mr%interp, order)
    end if

    POP_SUB(multiresolution_init)
  end subroutine multiresolution_init

  ! ------------------------------------------------------------
  subroutine interp_init(this, order)
    type(interp_t),          intent(inout) :: this
    integer,                 intent(in)    :: order

    FLOAT   :: pos(1:2*order)
    integer :: ii

    PUSH_SUB(interp_init)

    this%order = order
    this%nn = 2*this%order
    SAFE_ALLOCATE(this%ww(1:this%nn))
    SAFE_ALLOCATE(this%posi(1:this%nn))

    do ii = 1, this%order
      this%posi(ii) = 1 + 2*(ii - 1)
      this%posi(this%order + ii) = -this%posi(ii)
      pos(ii) = this%posi(ii)
      pos(this%order + ii) = -pos(ii)
    end do

    call interpolation_coefficients(this%nn, pos, M_ZERO, this%ww)

    POP_SUB(interp_init)
  end subroutine interp_init

  !--------------------------------------------------------------
  subroutine multiresolution_end(this)
    type(multiresolution_t), intent(inout) :: this

    PUSH_SUB(multiresolution_end)

    SAFE_DEALLOCATE_A(this%center)
    SAFE_DEALLOCATE_A(this%radius)
    SAFE_DEALLOCATE_A(this%interp%ww)
    SAFE_DEALLOCATE_A(this%interp%posi)

    POP_SUB(multiresolution_end)
  end subroutine multiresolution_end

  ! --------------------------------------------------------------
  subroutine multiresolution_copy(mr_out, mr_in)
    type(multiresolution_t), intent(out) :: mr_out
    type(multiresolution_t), intent(in)  :: mr_in

    PUSH_SUB(multiresolution_copy)

    mr_out%num_areas = mr_in%num_areas
    mr_out%num_radii  = mr_in%num_radii
    SAFE_ALLOCATE_SOURCE_A(mr_out%center, mr_in%center)
    SAFE_ALLOCATE_SOURCE_A(mr_out%radius, mr_in%radius)

    POP_SUB(multiresolution_copy)
  end subroutine multiresolution_copy

  ! --------------------------------------------------------------
  pure logical function multiresolution_use(this)
    type(multiresolution_t), intent(in)  :: this

    multiresolution_use = this%num_areas > 0

  end function multiresolution_use

  ! --------------------------------------------------------------
  subroutine multiresolution_dump(this, iunit)
    type(multiresolution_t), intent(in) :: this
    integer,                 intent(in) :: iunit

    integer :: ir, idir

    PUSH_SUB(multiresolution_dump)

    write(iunit, '(a20,i4)') 'num_areas=         ', this%num_areas
    write(iunit, '(a20,i4)') 'num_radii=         ', this%num_radii
    do ir = 1, this%num_radii
      write(iunit, '(a10,i2.2,a9,e22.14)') 'mr_radius_', ir, '=        ', this%radius(ir)
    end do
    do idir = 1, ubound(this%center, dim=1)
      write(iunit, '(a7,i1,a13,e22.14)')   'center(', idir, ')=           ', this%center(idir)
    end do

    POP_SUB(multiresolution_dump)
  end subroutine multiresolution_dump

  ! --------------------------------------------------------------
  subroutine multiresolution_load(this, dim, iunit, mpi_grp, ierr)
    type(multiresolution_t), intent(inout) :: this
    integer,                 intent(in)    :: dim
    integer,                 intent(in)    :: iunit
    type(mpi_grp_t),         intent(in)  :: mpi_grp
    integer,                 intent(inout) :: ierr

    integer :: ir, idir, err
    character(len=20)  :: str
    character(len=100), allocatable :: lines(:)

    PUSH_SUB(multiresolution_load)

    SAFE_ALLOCATE(lines(1:2))
    call iopar_read(mpi_grp, iunit, lines, 2, err)
    if (err /= 0) then
      ierr = ierr + 2**8
    else
      read(lines(1),*) str, this%num_areas
      read(lines(2),*) str, this%num_radii
    end if
    SAFE_DEALLOCATE_A(lines)

    SAFE_ALLOCATE(this%radius(1:this%num_radii))
    SAFE_ALLOCATE(lines(1:this%num_radii))
    this%num_radii = 0
    call iopar_read(mpi_grp, iunit, lines, this%num_radii, err)
    if (err /= 0) then
      ierr = ierr + 2**9
    else
      do ir = 1, this%num_radii
        read(lines(1),*) str, this%radius(ir)
      end do
    end if
    SAFE_DEALLOCATE_A(lines)

    SAFE_ALLOCATE(lines(1:dim))
    SAFE_ALLOCATE(this%center(1:dim))
    call iopar_read(mpi_grp, iunit, lines, dim, err)
    if (err /= 0) then
      ierr = ierr + 2**10
    else
      do idir = 1, dim
        read(lines(1), *) str, this%center(idir)
      end do
    end if
    SAFE_DEALLOCATE_A(lines)

    POP_SUB(multiresolution_load)
  end subroutine multiresolution_load

end module multiresolution_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
