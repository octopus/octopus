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
!! $Id$

#include "global.h"

module index_m
  use global_m
  use hypercube_m
  use io_m
  use io_binary_m
  use messages_m
  use simul_box_m

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
    index_load_lxyz,       &
    index_subset_indices

  type index_t
    type(hypercube_t) :: hypercube
    logical           :: is_hypercube     !< true if the box shape is an hypercube
    integer           :: dim              !< the dimension
    integer           :: nr(2, MAX_DIM)   !< dimensions of the box where the points are contained
    integer           :: ll(MAX_DIM)      !< literally nr(2,:) - nr(1,:) + 1 - 2*enlarge(:)
    integer, pointer  :: lxyz(:,:)        !< return x, y and z for each point
    integer, pointer  :: lxyz_inv(:,:,:)  !< return points # for each xyz
    integer           :: enlarge(MAX_DIM) !< number of points to add for boundary conditions
    integer(8)        :: checksum
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

    forall (idir = 1:idx%dim) ix2(idir) = ix(idir)
    forall (idir = idx%dim + 1:MAX_DIM) ix2(idir) = 0

    if(.not. idx%is_hypercube) then
      index = idx%lxyz_inv(ix2(1), ix2(2), ix2(3))
    else
      call hypercube_x_to_i(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ix, index)
    end if
    
  end function index_from_coords


  subroutine index_from_coords_vec(idx, npoints, ix, index)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: npoints
    integer,       intent(in)    :: ix(:, :)
    integer,       intent(out)   :: index(:)

    integer :: ix2(MAX_DIM), idir, ip

    PUSH_SUB(index_from_coords_vec)
    ix2 = 0

    if(.not. idx%is_hypercube) then
      do ip = 1, npoints
        forall (idir = 1:idx%dim) ix2(idir) = ix(idir, ip)
        index(ip) = idx%lxyz_inv(ix2(1), ix2(2), ix2(3))
      end do
    else
      do ip = 1, npoints
        call hypercube_x_to_i(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ix, index(ip))
      end do
    end if
    
    POP_SUB(index_from_coords_vec)
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
      forall (idir = 1:idx%dim) ix(idir) = idx%lxyz(ip, idir)
    else
      call hypercube_i_to_x(idx%hypercube, idx%dim, idx%nr, idx%enlarge(1), ip, ix)
    end if
  end subroutine index_to_coords


  !> This function returns .true. if the coordinates are inside the
  !! inner cube (without the enlargement).
  logical pure function coords_in_inner_cube(idx, ix) result(inside)
    type(index_t), intent(in)    :: idx
    integer,       intent(in)    :: ix(:)
    
    integer :: idir

    inside = .true.
    do idir = 1, idx%dim
      inside = inside &
           .and. (ix(idir) > idx%nr(1, idir) + idx%enlarge(idir)) &
           .and. (ix(idir) < idx%nr(2, idir) - idx%enlarge(idir))
    end do

  end function coords_in_inner_cube


  ! --------------------------------------------------------------
  subroutine index_dump(idx, dir, filename)
    type(index_t),    intent(in)    :: idx 
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: filename

    integer :: iunit, idir

    PUSH_SUB(index_dump)

    iunit = io_open(trim(dir)//trim(filename), action='write', position="append", is_tmp=.true.)

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

    call io_close(iunit)

    POP_SUB(index_dump)
  end subroutine index_dump


  ! --------------------------------------------------------------
  subroutine index_load(idx, dir, filename)
    type(index_t),    intent(inout) :: idx
    character(len=*), intent(in)    :: dir
    character(len=*), intent(in)    :: filename

    integer :: iunit, idir
    character(len=100) :: line
    character(len=20)  :: str

    PUSH_SUB(index_load)

    iunit = io_open(trim(dir)//trim(filename), action='read', status="old", die=.false., is_tmp=.true.)

    ! Find (and throw away) the dump tag.
    do
      read(iunit, '(a)') line
      if(trim(line) == dump_tag) exit
    end do

    idx%nr = 0
    idx%ll = 0
    idx%enlarge = 0

    read(iunit, '(a20,l1)')  str, idx%is_hypercube
    read(iunit, '(a20,i21)') str, idx%dim
    if (.not. idx%is_hypercube) then
      read(iunit, '(a20,7i8)')  str, (idx%nr(1, idir), idir = 1,idx%dim)
      read(iunit, '(a20,7i8)')  str, (idx%nr(2, idir), idir = 1,idx%dim)
      read(iunit, '(a20,7i8)')  str, idx%ll(1:idx%dim)
      read(iunit, '(a20,7i8)')  str, idx%enlarge(1:idx%dim)
    end if

    call io_close(iunit)

    POP_SUB(index_load)
  end subroutine index_load


  ! --------------------------------------------------------------
  subroutine index_dump_lxyz(idx, np, dir, ierr)
    type(index_t),    intent(in)  :: idx
    integer,          intent(in)  :: np
    character(len=*), intent(in)  :: dir
    integer,          intent(out) :: ierr

    PUSH_SUB(index_dump_lxyz)

    if(.not. idx%is_hypercube) then
      ASSERT(associated(idx%lxyz))
      call io_binary_write(trim(dir)//'/lxyz.obf', np*idx%dim, idx%lxyz, ierr)
    end if

    POP_SUB(index_dump_lxyz)
  end subroutine index_dump_lxyz


  ! --------------------------------------------------------------
  !> Fill the lxyz and lxyz_inv arrays from a file
  subroutine index_load_lxyz(idx, np, filename)
    type(index_t),    intent(inout) :: idx
    integer,          intent(in)    :: np
    character(len=*), intent(in)    :: filename

    integer :: ip, idir, ix(MAX_DIM), ierr

    PUSH_SUB(index_load_lxyz)

    ! FIXME 4D
    ASSERT(.not. idx%is_hypercube)

    call io_binary_read(trim(filename)//'.obf', np*idx%dim, idx%lxyz, ierr)

    if(ierr > 0) then
      message(1) = "Failed to read file "//trim(filename)//'.obf'
      call messages_fatal(1)
    end if

    do ip = 1, np
      forall (idir = 1:idx%dim) ix(idir) = idx%lxyz(ip, idir)
      forall (idir = idx%dim + 1:MAX_DIM) ix(idir) = 0
      idx%lxyz_inv(ix(1), ix(2), ix(3)) = ip
    end do

    POP_SUB(index_load_lxyz)
  end subroutine index_load_lxyz


  ! --------------------------------------------------------------
  !> Checks if the (x, y, z) indices of point are valid, i.e.
  !! inside the dimensions of the simulation box.
  logical function index_valid(idx, point)
    type(index_t), intent(in) :: idx
    integer,       intent(in) :: point(MAX_DIM)

    integer :: idir
    logical :: valid

    PUSH_SUB(index_valid)

    valid = .true.
    do idir = 1, idx%dim
      if(point(idir)  <  idx%nr(1, idir) .or. point(idir) > idx%nr(2, idir)) then
        valid = .false.
      end if
    end do

    index_valid = valid

    POP_SUB(index_valid)
  end function index_valid


  ! --------------------------------------------------------------
  !> Extracts the point numbers of a rectangular subset spanned
  !! by the two corner points from and to.
  subroutine index_subset_indices(idx, from, to, indices)
    type(index_t), intent(in)  :: idx
    integer,       intent(in)  :: from(MAX_DIM)
    integer,       intent(in)  :: to(MAX_DIM)
    integer,       intent(out) :: indices(:)

    integer :: lb(MAX_DIM) !< Lower bound of indices.
    integer :: ub(MAX_DIM) !< Upper bound of indices.

    integer :: ix, iy, iz, ii

    PUSH_SUB(index_subset_indices)

    ! In debug mode, check for valid indices in from, to first.
    if(in_debug_mode) then
      if(.not.index_valid(idx, from).or..not.index_valid(idx, to)) then
        message(1) = 'Failed assertion:'
        message(2) = 'mesh.mesh_subset_indices has been passed points outside the box:'
        message(3) = ''
        write(message (4), '(a, i6, a, i6, a, i6, a)') &
          '  from = (', from(1), ', ', from(2), ', ', from(3), ')'
        write(message(5), '(a, i6, a, i6, a, i6, a)') & 
          '  to   = (', to(1), ', ', to(2), ', ', to(3), ')'
        call messages_fatal(5)
      end if
    end if

    lb = min(from, to)
    ub = max(from, to)

    ii = 1
    do ix = lb(1), ub(1)
      do iy = lb(2), ub(2)
        do iz = lb(3), ub(3)
          indices(ii) = idx%lxyz_inv(ix, iy, iz)
          ii         = ii + 1
        end do
      end do
    end do

    POP_SUB(index_subset_indices)
  end subroutine index_subset_indices

end module index_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
