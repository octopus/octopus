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

#include "global.h"

module xyz_file_m
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
    FLOAT             :: x(MAX_DIM)

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

    PUSH_SUB(xyz_file_init)

    gf%file_type = XYZ_FILE_ERR
    gf%flags     = 0
    gf%n         = 0
    nullify(gf%atom)

    POP_SUB(xyz_file_init)
  end subroutine xyz_file_init


  ! ---------------------------------------------------------
  subroutine xyz_file_end(gf)
    type(xyz_file_info), intent(inout) :: gf

    PUSH_SUB(xyz_file_end)

    if(associated(gf%atom)) then
      SAFE_DEALLOCATE_P(gf%atom)
    end if
    call xyz_file_init(gf)

    POP_SUB(xyz_file_end)
  end subroutine xyz_file_end


  ! ---------------------------------------------------------
  subroutine xyz_file_read(what, gf, space)
    character(len=*),    intent(in)    :: what
    type(xyz_file_info), intent(inout) :: gf
    type(space_t),       intent(in)    :: space

    integer :: ia, ncol, iunit, jdir
    type(block_t) :: blk
    character(len=80) :: str
    logical :: done

    PUSH_SUB(xyz_file_read)

    done = .false.

    if(parse_isdef(datasets_check('PDB'//trim(what))) .ne. 0) then
      call check_duplicated(done)

      gf%file_type = XYZ_FILE_PDB
      gf%flags = ior(gf%flags, XYZ_FLAGS_RESIDUE)
      gf%flags = ior(gf%flags, XYZ_FLAGS_CHARGE)

      call parse_string(datasets_check('PDB'//trim(what)), 'coords.pdb', str)

      iunit = io_open(str, action='read')
      call xyz_file_read_PDB(what, iunit, gf)
      call io_close(iunit)
    end if

    if(parse_isdef(datasets_check('XYZ'//trim(what))) .ne. 0) then ! read a xyz file
      call check_duplicated(done)

      gf%file_type = XYZ_FILE_XYZ
      call parse_string(datasets_check('XYZ'//trim(what)), 'coords.xyz', str)

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

      gf%file_type = XYZ_FILE_INP
      gf%flags = ior(gf%flags, XYZ_FLAGS_MOVE)

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        ncol = parse_block_cols(blk, ia - 1)
        if((ncol .lt. space%dim + 1) .or. (ncol .gt. space%dim + 2)) then
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

    POP_SUB(xyz_file_read)

  contains
    
    subroutine check_duplicated(done)
      logical, intent(inout) :: done
      
      if(.not. done) then
        done = .true.
      else
        message(1) = 'Multiple definitions of '//trim(what)//' in the input file.'
        call messages_fatal(1)
      end if
    end subroutine check_duplicated

  end subroutine xyz_file_read


  ! ---------------------------------------------------------
  subroutine xyz_file_read_PDB(what, iunit, gf)
    character(len=*),    intent(in)    :: what
    integer,             intent(in)    :: iunit
    type(xyz_file_info), intent(inout) :: gf

    character(len=80) :: record
    character(len=6)  :: record_name
    integer :: na

    PUSH_SUB(xyz_file_read_PDB)

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
          read(record, '(30x,3f8.3,6x,f5.2)') gf%atom(na)%x, gf%atom(na)%charge
        else
          read(record, '(30x,3f8.3)') gf%atom(na)%x
        endif

        na = na + 1
      end if
    end do
991 continue

    POP_SUB(xyz_file_read_PDB)
  end subroutine xyz_file_read_PDB

end module xyz_file_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
