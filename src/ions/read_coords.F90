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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module read_coords_m
  use datasets_m
  use global_m
  use io_m
  use parser_m
  use messages_m
  use profiling_m
  use space_m
  use string_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                     &
    read_coords_atom,              &
    read_coords_info,              &
    read_coords_init,              &
    read_coords_end,               &
    read_coords_read

  !> for read_coords_info::file_type
  integer, public, parameter :: &
    READ_COORDS_ERR      = 0,      &
    READ_COORDS_PDB      = 1,      &
    READ_COORDS_XYZ      = 2,      &
    READ_COORDS_INP      = 3,      &
    READ_COORDS_REDUCED  = 4

  !> for read_coords_info::flags
  integer, public, parameter :: &
    XYZ_FLAGS_RESIDUE = 1,      &
    XYZ_FLAGS_CHARGE  = 2,      &
    XYZ_FLAGS_MOVE    = 4

  type read_coords_atom
    character(len=15) :: label  !< stuff that is always known
    FLOAT             :: x(MAX_DIM)

    FLOAT             :: charge !< stuff specific to PDB files
    character(len=3)  :: residue
    logical           :: move   !< stuff specific to the inp file
  end type read_coords_atom

  type read_coords_info
    integer :: source
    integer :: flags

    integer :: n                !< number of atoms in file
    type(read_coords_atom), pointer :: atom(:)
  end type read_coords_info


contains

  ! ---------------------------------------------------------
  subroutine read_coords_init(gf)
    type(read_coords_info), intent(out) :: gf

    PUSH_SUB(read_coords_init)

    gf%source = READ_COORDS_ERR
    gf%flags     = 0
    gf%n         = 0
    nullify(gf%atom)

    POP_SUB(read_coords_init)
  end subroutine read_coords_init


  ! ---------------------------------------------------------
  subroutine read_coords_end(gf)
    type(read_coords_info), intent(inout) :: gf

    PUSH_SUB(read_coords_end)

    if(associated(gf%atom)) then
      SAFE_DEALLOCATE_P(gf%atom)
    end if
    call read_coords_init(gf)

    POP_SUB(read_coords_end)
  end subroutine read_coords_end


  ! ---------------------------------------------------------
  subroutine read_coords_read(what, gf, space)
    character(len=*),    intent(in)    :: what
    type(read_coords_info), intent(inout) :: gf
    type(space_t),       intent(in)    :: space

    integer :: ia, ncol, iunit, jdir
    type(block_t) :: blk
    character(len=80) :: str
    logical :: done

    PUSH_SUB(read_coords_read)

    done = .false.

    if(parse_isdef(datasets_check('PDB'//trim(what))) /= 0) then
      call check_duplicated(done)

      gf%source = READ_COORDS_PDB
      gf%flags = ior(gf%flags, XYZ_FLAGS_RESIDUE)
      gf%flags = ior(gf%flags, XYZ_FLAGS_CHARGE)

      ! no default, since we do not do this unless the input tag is present
      call parse_string(datasets_check('PDB'//trim(what)), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1)

      iunit = io_open(str, action='read')
      call read_coords_read_PDB(what, iunit, gf)
      call io_close(iunit)
    end if

    if(parse_isdef(datasets_check('XYZ'//trim(what))) /= 0) then ! read an xyz file
      call check_duplicated(done)

      gf%source = READ_COORDS_XYZ
      ! no default, since we do not do this unless the input tag is present
      call parse_string(datasets_check('XYZ'//trim(what)), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1)

      iunit = io_open(str, status='old', action='read', is_tmp=.true.)
      read(iunit, *) gf%n
      read(iunit, *) ! skip comment line

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if
    
    if(parse_block(datasets_check(trim(what)), blk) == 0) then
      call check_duplicated(done)

      gf%n = parse_block_n(blk)

      gf%source = READ_COORDS_INP
      gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

      message(1) = "Reading " // trim(what) // " from " // trim(what) // " block"
      call messages_info(1)

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        ncol = parse_block_cols(blk, ia - 1)
        if((ncol  <  space%dim + 1) .or. (ncol > space%dim + 2)) then
          write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', ia
          call messages_fatal(1)
        end if
        call parse_block_string (blk, ia - 1, 0, gf%atom(ia)%label)
        do jdir = 1, space%dim
          call parse_block_float  (blk, ia - 1, jdir, gf%atom(ia)%x(jdir))
        end do
        if(ncol == space%dim + 2) then
          call parse_block_logical(blk, ia - 1, space%dim + 1, gf%atom(ia)%move)
        else
          gf%atom(ia)%move = .true.
        end if
      end do

      call parse_block_end(blk)
    end if

    ! This is valid only for Coordinates.
    if(trim(what) == 'Coordinates' .and. parse_block(datasets_check('Reduced'//trim(what)), blk) == 0) then
      call check_duplicated(done)

      gf%n = parse_block_n(blk)

      gf%source = READ_COORDS_REDUCED
      gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

      message(1) = "Reading " // trim(what) // " from Reduced" // trim(what) // " block"
      call messages_info(1)

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        ncol = parse_block_cols(blk, ia - 1)
        if((ncol  <  space%dim + 1) .or. (ncol > space%dim + 2)) then
          write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', ia
          call messages_fatal(1)
        end if
        call parse_block_string (blk, ia - 1, 0, gf%atom(ia)%label)
        do jdir = 1, space%dim
          call parse_block_float  (blk, ia - 1, jdir, gf%atom(ia)%x(jdir))
        end do
        if(ncol == space%dim + 2) then
          call parse_block_logical(blk, ia - 1, space%dim + 1, gf%atom(ia)%move)
        else
          gf%atom(ia)%move = .true.
        end if
      end do

      call parse_block_end(blk)
    end if

    ! adjust units
    do ia = 1, gf%n
      do jdir = space%dim + 1, MAX_DIM
        gf%atom(ia)%x(jdir) = M_ZERO
      end do
      gf%atom(ia)%x = units_to_atomic(units_inp%length, gf%atom(ia)%x)
    end do

    POP_SUB(read_coords_read)

  contains
    
    subroutine check_duplicated(done)
      logical, intent(inout) :: done
      
      PUSH_SUB(read_coords_read.check_duplicated)

      if(.not. done) then
        done = .true.
      else
        message(1) = 'Multiple definitions of '//trim(what)//' in the input file.'
        call messages_fatal(1)
      end if

      POP_SUB(read_coords_read.check_duplicated)
    end subroutine check_duplicated

  end subroutine read_coords_read


  ! ---------------------------------------------------------
  subroutine read_coords_read_PDB(what, iunit, gf)
    character(len=*),    intent(in)    :: what
    integer,             intent(in)    :: iunit
    type(read_coords_info), intent(inout) :: gf

    character(len=80) :: record
    character(len=6)  :: record_name
    integer :: na

    PUSH_SUB(read_coords_read_PDB)

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

    SAFE_ALLOCATE(gf%atom(1:gf%n))

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

        if(trim(what) == 'Classical') then
          read(record, '(30x,3f8.3,6x,f5.2)') gf%atom(na)%x(1:3), gf%atom(na)%charge
        else
          read(record, '(30x,3f8.3)') gf%atom(na)%x(1:3)
        endif

        na = na + 1
      end if
    end do
991 continue

    POP_SUB(read_coords_read_PDB)
  end subroutine read_coords_read_PDB

end module read_coords_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
