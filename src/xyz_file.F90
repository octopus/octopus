!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module xyz_file
  use global
  use messages
  use datasets_mod
  use string
  use units
  use lib_oct_parser
  use io

  implicit none

  private
  public ::                     &
    xyz_file_atom,              &
    xyz_file_info,              &
    xyz_file_init,              &
    xyz_file_end,               &
    xyz_file_read

  integer, public, parameter :: &
    XYZ_FILE_ERR      = 0,      &
    XYZ_FILE_PDB      = 1,      &
    XYZ_FILE_XYZ      = 2,      &
    XYZ_FILE_INP      = 3

  integer, public, parameter :: &
    XYZ_FLAGS_RESIDUE = 1,      &
    XYZ_FLAGS_CHARGE  = 2,      &
    XYZ_FLAGS_MOVE    = 4

  type xyz_file_atom
    character(len=15) :: label  ! stuff that is always known
    FLOAT             :: x(3)

    FLOAT             :: charge ! stuff specific to PDB files
    character(len=3)  :: residue
    logical           :: move   ! stuff specific to the inp file
  end type xyz_file_atom

  type xyz_file_info
    integer :: file_type
    integer :: flags

    integer :: n                ! number of atoms in file
    type(xyz_file_atom), pointer :: atom(:)
  end type xyz_file_info


contains

  ! ---------------------------------------------------------
  subroutine xyz_file_init(gf)
    type(xyz_file_info), intent(out) :: gf

    gf%file_type = XYZ_FILE_ERR
    gf%flags     = 0
    gf%n         = 0
    nullify(gf%atom)
  end subroutine xyz_file_init


  ! ---------------------------------------------------------
  subroutine xyz_file_end(gf)
    type(xyz_file_info), intent(inout) :: gf

    if(associated(gf%atom)) then
      deallocate(gf%atom)
    end if
    call xyz_file_init(gf)
  end subroutine xyz_file_end


  ! ---------------------------------------------------------
  subroutine xyz_file_read(what, gf)
    character(len=*),    intent(in)    :: what
    type(xyz_file_info), intent(inout) :: gf

    integer :: i, j, iunit
    integer(POINTER_SIZE) :: blk
    character(len=80) :: str

    call push_sub('xyz_file.xyz_file_read')

    if(loct_parse_isdef(check_inp('PDB'//trim(what))).ne.0) then
      gf%file_type = XYZ_FILE_PDB
      gf%flags = ior(gf%flags, XYZ_FLAGS_RESIDUE)
      gf%flags = ior(gf%flags, XYZ_FLAGS_CHARGE)

      call loct_parse_string(check_inp('PDB'//trim(what)), 'coords.pdb', str)

      iunit = io_open(str, action='read')
      call xyz_file_read_PDB(iunit, gf)
      call io_close(iunit)

    else if(loct_parse_isdef(check_inp('XYZ'//trim(what))).ne.0) then ! read a xyz file
      gf%file_type = XYZ_FILE_XYZ
      call loct_parse_string(check_inp('XYZ'//trim(what)), 'coords.xyz', str)

      iunit = io_open(str, action='read')
      read(iunit, *) gf%n
      read(iunit, *) ! skip comment line

      ALLOCATE(gf%atom(gf%n), gf%n)

      do i = 1, gf%n
        read(iunit,*) gf%atom(i)%label, gf%atom(i)%x(:)
      end do

      call io_close(iunit)

    else if(loct_parse_block(check_inp(trim(what)), blk) == 0) then
      gf%n = loct_parse_block_n(blk)

      gf%file_type = XYZ_FILE_INP
      gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

      ALLOCATE(gf%atom(gf%n), gf%n)

      do i = 1, gf%n
        j = loct_parse_block_cols(blk, i-1)
        if((j.ne.4).and.(j.ne.5)) then
          write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', i
          call write_fatal(1)
        end if
        call loct_parse_block_string (blk, i-1, 0, gf%atom(i)%label)
        call loct_parse_block_float  (blk, i-1, 1, gf%atom(i)%x(1))
        call loct_parse_block_float  (blk, i-1, 2, gf%atom(i)%x(2))
        call loct_parse_block_float  (blk, i-1, 3, gf%atom(i)%x(3))
        if(j == 5) then
          call loct_parse_block_logical(blk, i-1, 4, gf%atom(i)%move)
        else
          gf%atom(i)%move = .true.
        end if
      end do
      call loct_parse_block_end(blk)

    end if

    ! adjust units
    do i = 1, gf%n
      gf%atom(i)%x = gf%atom(i)%x * units_inp%length%factor
    end do

    call pop_sub()
  end subroutine xyz_file_read


  ! ---------------------------------------------------------
  subroutine xyz_file_read_PDB(iunit, gf)
    integer,             intent(in)    :: iunit
    type(xyz_file_info), intent(inout) :: gf

    character(len=80) :: record
    character(len=6)  :: record_name
    integer :: na

    ! First count number of atoms
    rewind(iunit)
    do
      read(iunit, '(a80)', err=990, end=990) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        gf%n = gf%n + 1
      end if
    end do
990 continue

    ALLOCATE(gf%atom(gf%n), gf%n)

    ! read in the data
    rewind(iunit)
    na = 1
    do
      read(iunit, '(a80)', err=991, end=991) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') gf%atom(na)%label, gf%atom(na)%residue
        call str_trim(gf%atom(na)%label)
        gf%atom(na)%label = gf%atom(na)%label(1:1)
        call str_trim(gf%atom(na)%residue)

        read(record, '(30x,3f8.3)') gf%atom(na)%x
        na = na + 1
      end if
    end do
991 continue

  end subroutine xyz_file_read_PDB

end module xyz_file
