!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module index_oct_m
  use global_oct_m
  use iihash_oct_m
  use io_oct_m
  use io_binary_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m

  implicit none

  private

  public ::                 &
    index_t,                &
    index_from_coords,      &
    index_from_coords_vec,  &
    index_to_coords,        &
    index_hilbert_to_point, &
    index_point_to_hilbert, &
    index_dump,             &
    index_load

  type index_t
    ! Components are public by default
    integer              :: dim              !< the dimension
    integer              :: nr(2, MAX_DIM)   !< dimensions of the box where the points are contained
    integer              :: ll(MAX_DIM)      !< literally nr(2,:) - nr(1,:) + 1 - 2*enlarge(:)
    integer, allocatable :: lxyz(:,:)        !< return x, y and z for each point
    integer, allocatable :: lxyz_inv(:,:,:)  !< return points # for each xyz
    integer              :: enlarge(MAX_DIM) !< number of points to add for boundary conditions
    integer(8)           :: checksum
    integer(8), allocatable :: grid_to_hilbert(:) !< map: local grid index -> Hilbert index
    type(lihash_t)       :: hilbert_to_grid    !< inverse map: Hilbert index -> local grid index
    integer(8), allocatable :: grid_to_hilbert_global(:) !< map: global grid index -> Hilbert index
    type(lihash_t)       :: hilbert_to_grid_global       !< inverse map: Hilbert index -> global grid index
    integer              :: bits               !< bits per dimension for Hilbert index
    integer              :: offset(MAX_DIM)    !< offset for getting the indices from the Hilbert index
  end type index_t


  character(len=18), parameter :: dump_tag = '*** index_dump ***'

  interface
    subroutine hilbert_index_to_point_int(dim, nbits, index, point)
      implicit none

      integer,    intent(in)       :: dim
      integer,    intent(in)       :: nbits
      integer,    intent(in)       :: index
      integer,    intent(out)      :: point !< (1:3)
    end subroutine hilbert_index_to_point_int

    subroutine hilbert_point_to_index_int(dim, nbits, index, point)
      implicit none

      integer,    intent(in)       :: dim
      integer,    intent(in)       :: nbits
      integer,    intent(out)      :: index
      integer,    intent(in)       :: point !< (1:3)
    end subroutine hilbert_point_to_index_int
  end interface

contains

  !> This function takes care of the boundary conditions for a given
  !! vector of integer coordinates it returns the true _global_ index
  !! of the point.
  integer function index_from_coords(idx, ix) result(index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ix(:)

    integer :: ix2(MAX_DIM), idir
    integer(8) :: ihilbert
    logical :: found

    ! No PUSH SUB, called too often

    do idir = 1, idx%dim
      ix2(idir) = ix(idir)
    end do
    do idir=idx%dim + 1, MAX_DIM
      ix2(idir) = 0
    end do

    call index_point_to_hilbert(idx, idx%dim, ihilbert, ix2)
    index = lihash_lookup(idx%hilbert_to_grid_global, ihilbert, found)
    if(.not. found) index = 0
    
  end function index_from_coords

  ! -----------------------------------------------------------------

  subroutine index_from_coords_vec(idx, npoints, ix, index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: npoints
    integer,       intent(in)    :: ix(:, :)
    integer,       intent(out)   :: index(:)

    integer :: ix2(MAX_DIM), idir, ip
    integer(8) :: ihilbert
    logical :: found

    ! No PUSH SUB, called too often
    ix2 = 0

    do ip = 1, npoints
      do idir = 1, idx%dim
        ix2(idir) = ix(idir, ip)
      end do
      !index(ip) = idx%lxyz_inv(ix2(1), ix2(2), ix2(3))
      call index_point_to_hilbert(idx, idx%dim, ihilbert, ix2)
      index(ip) = lihash_lookup(idx%hilbert_to_grid_global, ihilbert, found)
      if(.not. found) index(ip) = 0
    end do
    
  end subroutine index_from_coords_vec


  !> Given a _global_ point index, this function returns the set of
  !! integer coordinates of the point.
  subroutine index_to_coords(idx, ip, ix)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ip
    integer,       intent(out)   :: ix(:)

    integer :: idir 

    ! We set all ix to zero first (otherwise the non-existent dimensions would be 
    ! undefined on exit).
    ix = 0
    call index_hilbert_to_point(idx, idx%dim, idx%grid_to_hilbert_global(ip), ix)
  end subroutine index_to_coords

  ! --------------------------------------------------------------
  subroutine index_dump(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),    intent(in)  :: idx
    integer,          intent(in)  :: np
    character(len=*), intent(in)  :: dir
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: ierr

    integer :: err

    PUSH_SUB(index_dump)

    ierr = 0

    if (int(np, 8)*idx%dim > huge(0)) then
      message(1) = "Too many global mesh points to write restart file for lxyz.obf."
      call messages_fatal(1)
    end if
    if (mpi_grp_is_root(mpi_grp)) then
      ! the index array is a global function and only root will write
      ASSERT(allocated(idx%grid_to_hilbert_global))
      call io_binary_write(trim(io_workpath(dir, namespace))//"/indices.obf", np, &
        idx%grid_to_hilbert_global, err)
      if (err /= 0) then
        ierr = ierr + 1
        message(1) = "Unable to write index function to '"//trim(dir)//"/indices.obf'."
        call messages_warning(1) 
      end if
    end if

#if defined(HAVE_MPI)
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, mpi_grp%comm, mpi_err)
#endif

    POP_SUB(index_dump)
  end subroutine index_dump


  ! --------------------------------------------------------------
  !> Load the index arrays from a file
  subroutine index_load(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),     intent(inout) :: idx
    integer,           intent(in)    :: np
    character(len=*),  intent(in)    :: dir
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(out)   :: ierr

    integer :: ip, idir, ix(MAX_DIM), err
    logical :: exists

    PUSH_SUB(index_load)

    ierr = 0

    ASSERT(allocated(idx%grid_to_hilbert_global))

    if (mpi_grp_is_root(mpi_grp)) then
      ! check for existence of lxyz.obf, print error message if found
      inquire(file=trim(trim(io_workpath(dir, namespace))//"/lxyz.obf"), exist=exists)
      if (exists) then
        message(1) = "Found lxyz.obf file. This means you created the restart files with an old version of the code."
        message(2) = "Please generate the restart files again with the current version of the code"
        message(3) = "because the internal format has changed."
        call messages_fatal(3)
      end if
      ! the index array is a global function and only root will write
      call io_binary_read(trim(io_workpath(dir, namespace))//"/indices.obf", np, &
        idx%grid_to_hilbert_global, err)
      if (err /= 0) then
        ierr = ierr + 1
        message(1) = "Unable to read index function from '"//trim(dir)//"/indices.obf'."
        call messages_warning(1)
      end if
    end if

#if defined(HAVE_MPI)
    ! Broadcast the results and synchronize
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, mpi_grp%comm, mpi_err)
    if (ierr == 0) then
      call MPI_Bcast(idx%grid_to_hilbert(1), np, MPI_LONG_LONG, 0, mpi_grp%comm, mpi_err)
    end if
#endif

    ! fill global hash map
    call lihash_init(idx%hilbert_to_grid_global)
    do ip = 1, np
      call lihash_insert(idx%hilbert_to_grid_global, idx%grid_to_hilbert_global(ip), ip)
    end do

    POP_SUB(index_load)
  end subroutine index_load

  subroutine index_hilbert_to_point(idx, dim, ihilbert, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(8),    intent(in)  :: ihilbert
    integer,       intent(out) :: point(:)

    ! no push_sub/pop_sub, called too often
    call hilbert_index_to_point(dim, idx%bits, ihilbert, point(1))
    point(1:dim) = point(1:dim) - idx%offset(1:dim)
  end subroutine index_hilbert_to_point

  subroutine index_point_to_hilbert(idx, dim, ihilbert, point)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: dim
    integer(8),    intent(out) :: ihilbert
    integer,       intent(in)  :: point(:)

    integer :: point_copy(1:MAX_DIM)

    ! no push_sub/pop_sub, called too often
    point_copy(1:dim) = point(1:dim) + idx%offset(1:dim)
    call hilbert_point_to_index(dim, idx%bits, ihilbert, point_copy(1))
  end subroutine index_point_to_hilbert
end module index_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
