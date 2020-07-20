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
  use hypercube_oct_m
  use io_oct_m
  use io_binary_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m

  implicit none

  private
  public ::                &
    index_t,               &
    index_from_coords,     &
    index_from_coords_vec, &
    index_to_coords,       &
    index_dump,            &
    index_load,            &
    index_dump_lxyz,       &
    index_load_lxyz

  type index_t
    ! Components are public by default
    type(hypercube_t)    :: hypercube
    logical              :: is_hypercube     !< true if the box shape is an hypercube
    integer              :: dim              !< the dimension
    integer              :: nr(2, MAX_DIM)   !< dimensions of the box where the points are contained
    integer              :: ll(MAX_DIM)      !< literally nr(2,:) - nr(1,:) + 1 - 2*enlarge(:)
    integer, allocatable :: lxyz(:,:)        !< return x, y and z for each point
    integer, allocatable :: lxyz_inv(:,:,:)  !< return points # for each xyz
    integer              :: enlarge(MAX_DIM) !< number of points to add for boundary conditions
    integer(8)           :: checksum
  end type index_t


  character(len=18), parameter :: dump_tag = '*** index_dump ***'

contains

  !> This function takes care of the boundary conditions for a given
  !! vector of integer coordinates it returns the true _global_ index
  !! of the point.
  integer function index_from_coords(idx, ix) result(index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ix(:)

    integer :: ix2(MAX_DIM), idir

    ! No PUSH SUB, called too often

    do idir = 1, idx%dim
      ix2(idir) = ix(idir)
    end do
    do idir=idx%dim + 1, MAX_DIM
      ix2(idir) = 0
    end do

    if(.not. idx%is_hypercube) then
      index = idx%lxyz_inv(ix2(1), ix2(2), ix2(3))
    else
      call hypercube_x_to_i(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ix, index)
    end if
    
  end function index_from_coords

  ! -----------------------------------------------------------------

  subroutine index_from_coords_vec(idx, npoints, ix, index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: npoints
    integer,       intent(in)    :: ix(:, :)
    integer,       intent(out)   :: index(:)

    integer :: ix2(MAX_DIM), idir, ip

    ! No PUSH SUB, called too often
    ix2 = 0

    if(.not. idx%is_hypercube) then
      do ip = 1, npoints
        do idir = 1, idx%dim
          ix2(idir) = ix(idir, ip)
        end do
        index(ip) = idx%lxyz_inv(ix2(1), ix2(2), ix2(3))
      end do
    else
      do ip = 1, npoints
        call hypercube_x_to_i(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ix(:, ip), index(ip))
      end do
    end if
    
  end subroutine index_from_coords_vec


  !> Given a _global_ point index, this function returns the set of
  !! integer coordinates of the point.
  pure subroutine index_to_coords(idx, ip, ix)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ip
    integer,       intent(out)   :: ix(:)

    integer :: idir 

    ! We set all ix to zero first (otherwise the non-existent dimensions would be 
    ! undefined on exit).
    ix = 0
    if(.not. idx%is_hypercube) then
      do idir = 1, idx%dim
        ix(idir) = idx%lxyz(ip, idir)
      end do
    else
      call hypercube_i_to_x(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ip, ix)
    end if
  end subroutine index_to_coords


  ! --------------------------------------------------------------
  subroutine index_dump(idx, dir, filename, mpi_grp, namespace, ierr)
    type(index_t),     intent(in)  :: idx 
    character(len=*),  intent(in)  :: dir
    character(len=*),  intent(in)  :: filename
    type(mpi_grp_t),   intent(in)  :: mpi_grp
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(out) :: ierr

    integer :: iunit, idir

    PUSH_SUB(index_dump)

    ierr = 0

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action='write', &
      position="append", die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      ierr = ierr + 1
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1)
    else
      !Only root writes to the file
      if (mpi_grp_is_root(mpi_grp)) then
        write(iunit, '(a)') dump_tag
        write(iunit, '(a20,l1)')  'is_hypercube=       ', idx%is_hypercube
        write(iunit, '(a20,i21)') 'dim=                ', idx%dim
        if (.not. idx%is_hypercube) then
          write(iunit, '(a20,7i8)') 'nr(1, :)=           ', (idx%nr(1, idir), idir = 1, idx%dim)
          write(iunit, '(a20,7i8)') 'nr(2, :)=           ', (idx%nr(2, idir), idir = 1, idx%dim)
          write(iunit, '(a20,7i8)') 'l(:)=               ', idx%ll(1:idx%dim)
          write(iunit, '(a20,7i8)') 'enlarge(:)=         ', idx%enlarge(1:idx%dim)
          ! The next two lines should always come last
          write(iunit, '(a20,i21)') 'algorithm=          ', 1
          write(iunit, '(a20,i21)') 'checksum=           ', idx%checksum
        end if
      end if

      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(index_dump)
  end subroutine index_dump


  ! --------------------------------------------------------------
  subroutine index_load(idx, dir, filename, mpi_grp, namespace, ierr)
    type(index_t),     intent(inout) :: idx
    character(len=*),  intent(in)    :: dir
    character(len=*),  intent(in)    :: filename
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(out)   :: ierr

    integer :: iunit, idir, err
    character(len=100) :: lines(6)
    character(len=20)  :: str

    PUSH_SUB(index_load)

    ierr = 0

    idx%nr = 0
    idx%ll = 0
    idx%enlarge = 0
    idx%is_hypercube = .false.

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action="read", &
      status="old", die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      ierr = ierr + 1
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1)
    else
      ! Find the dump tag.
      call iopar_find_line(mpi_grp, iunit, dump_tag, err)
      if (err /= 0) ierr = ierr + 2

      if (ierr == 0) then
        idx%is_hypercube = .false.
        call iopar_read(mpi_grp, iunit, lines, 2, err)
        if (err /= 0) then
          ierr = ierr + 4
        else
          read(lines(1), '(a20,l1)')  str, idx%is_hypercube
          read(lines(2), '(a20,i21)') str, idx%dim
        end if

        if (.not. idx%is_hypercube) then
          call iopar_read(mpi_grp, iunit, lines, 6, err)
          if (err /= 0) then
            ierr = ierr + 8            
          else
            read(lines(1), '(a20,7i8)') str, (idx%nr(1, idir), idir = 1,idx%dim)
            read(lines(2), '(a20,7i8)') str, (idx%nr(2, idir), idir = 1,idx%dim)
            read(lines(3), '(a20,7i8)') str, idx%ll(1:idx%dim)
            read(lines(4), '(a20,7i8)') str, idx%enlarge(1:idx%dim)
            !For the moment we do not read the algorithm, as it is always set to 1
            read(lines(6), '(a20,i21)') str, idx%checksum
          end if
        end if
      end if
      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(index_load)
  end subroutine index_load


  ! --------------------------------------------------------------
  subroutine index_dump_lxyz(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),    intent(in)  :: idx
    integer,          intent(in)  :: np
    character(len=*), intent(in)  :: dir
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: ierr

    integer :: err

    PUSH_SUB(index_dump_lxyz)

    ierr = 0

    if (.not. idx%is_hypercube) then
      if (mpi_grp_is_root(mpi_grp)) then
        ! lxyz is a global function and only root will write
        ASSERT(allocated(idx%lxyz))
        call io_binary_write(trim(io_workpath(dir, namespace))//"/lxyz.obf", np*idx%dim, idx%lxyz, err)
        if (err /= 0) then
          ierr = ierr + 1
          message(1) = "Unable to write index function to '"//trim(dir)//"/lxyz.obf'."
          call messages_warning(1) 
        end if
      end if

#if defined(HAVE_MPI)
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, mpi_grp%comm, mpi_err)
      call MPI_Barrier(mpi_grp%comm, mpi_err)
#endif
    end if

    POP_SUB(index_dump_lxyz)
  end subroutine index_dump_lxyz


  ! --------------------------------------------------------------
  !> Fill the lxyz and lxyz_inv arrays from a file
  subroutine index_load_lxyz(idx, np, dir, mpi_grp, namespace, ierr)
    type(index_t),     intent(inout) :: idx
    integer,           intent(in)    :: np
    character(len=*),  intent(in)    :: dir
    type(mpi_grp_t),   intent(in)    :: mpi_grp
    type(namespace_t), intent(in)    :: namespace
    integer,           intent(out)   :: ierr

    integer :: ip, idir, ix(MAX_DIM), err

    PUSH_SUB(index_load_lxyz)

    ierr = 0

    if (.not. idx%is_hypercube) then
      ASSERT(allocated(idx%lxyz))

      if (mpi_grp_is_root(mpi_grp)) then
        ! lxyz is a global function and only root will write
        call io_binary_read(trim(io_workpath(dir, namespace))//"/lxyz.obf", np*idx%dim, idx%lxyz, err)
        if (err /= 0) then
          ierr = ierr + 1
          message(1) = "Unable to read index function from '"//trim(dir)//"/lxyz.obf'."
          call messages_warning(1)
        end if
      end if

#if defined(HAVE_MPI)
      ! Broadcast the results and synchronize
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, mpi_grp%comm, mpi_err)
      if (ierr == 0) then
        call MPI_Bcast(idx%lxyz(1,1), np*idx%dim, MPI_INTEGER, 0, mpi_grp%comm, mpi_err)
      end if
      call MPI_Barrier(mpi_grp%comm, mpi_err)
#endif

      ! Compute lxyz_inv from lxyz
      if (ierr == 0) then
        do ip = 1, np
          do idir = 1, idx%dim
            ix(idir) = idx%lxyz(ip, idir)
          end do
          do idir=idx%dim + 1, MAX_DIM
            ix(idir) = 0
          end do
          idx%lxyz_inv(ix(1), ix(2), ix(3)) = ip
        end do
      end if

    end if

    POP_SUB(index_load_lxyz)
  end subroutine index_load_lxyz

end module index_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
