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

#include "global.h"

module read_coords_oct_m
  use global_oct_m
  use io_oct_m
  use parser_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m

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

    !> variables for passing info from XSF input to simul_box_init
    integer :: periodic_dim
    FLOAT :: lsize(MAX_DIM)
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

    gf%periodic_dim = -1
    gf%lsize(:) = -M_ONE

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
    character(len=*),       intent(in)    :: what
    type(read_coords_info), intent(inout) :: gf
    type(space_t),          intent(in)    :: space

    integer :: ia, ncol, iunit, jdir
    type(block_t) :: blk
    character(len=256) :: str
    logical :: done

    PUSH_SUB(read_coords_read)

    done = .false.

    !%Variable XYZCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If <tt>PDBCoordinates</tt> is not present, the program reads the atomic coordinates from
    !% the XYZ file specified by the variable <tt>XYZCoordinates</tt> -- in case this variable
    !% is present. The XYZ format is very simple: The first line of the file has an integer
    !% indicating the number of atoms. The second can contain comments that are simply ignored by
    !% <tt>Octopus</tt>. Then there follows one line per atom, containing the chemical species and
    !% the Cartesian coordinates of the atom.
    !%
    !% WARNING: By default the coordinates are treated in the units
    !% specified by <tt>Units</tt> and/or <tt>UnitsInput</tt>, which
    !% means Octopus might expect xyz files to be in atomic units. If
    !% you want the XYZ file to be read in Angstrom, as most codes do,
    !% you can set the variable <tt>UnitsXYZFiles</tt> to
    !% <tt>angstrom</tt>.
    !%End

    if(parse_is_defined('XYZ'//trim(what))) then ! read an xyz file
      call check_duplicated(done)

      gf%source = READ_COORDS_XYZ
      ! no default, since we do not do this unless the input tag is present
      call parse_variable('XYZ'//trim(what), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1)

      iunit = io_open(str, status='old', action='read')
      read(iunit, *) gf%n

      if(gf%n <= 0) then
        write(message(1),'(a,i6)') "Invalid number of atoms ", gf%n
        call messages_fatal(1)
      end if

      read(iunit, *) ! skip comment line

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if
    
    !%Variable Coordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% If <tt>XYZCoordinates</tt>, <tt>PDBCoordinates</tt>, and <tt>XSFCoordinates</tt> were not found,
    !% <tt>Octopus</tt> tries to read the coordinates for the atoms from the block <tt>Coordinates</tt>. The
    !% format is quite straightforward:
    !%
    !% <tt>%Coordinates
    !% <br>&nbsp;&nbsp;'C' |      -0.56415 | 0.0 | 0.0 | no
    !% <br>&nbsp;&nbsp;'O' | &nbsp;0.56415 | 0.0 | 0.0 | no
    !% <br>%</tt>
    !%
    !% The first line defines a carbon atom at coordinates (-0.56415, 0.0, 0.0),
    !% that is <b>not</b> allowed to move during dynamical simulations. The second line has
    !% a similar meaning. This block obviously defines a carbon monoxide molecule, if the
    !% input units are <tt>eV_Angstrom</tt>. The number of coordinates for each species
    !% must be equal to the dimension of your space (generally 3).
    !% Note that in this way it is possible to fix some of the atoms (this
    !% is not possible when specifying the coordinates through a <tt>PDBCoordinates</tt> or
    !% <tt>XYZCoordinates</tt> file). The last column is optional, and the default is yes.
    !% It is always possible to fix <b>all</b> atoms using the <tt>MoveIons</tt> directive.
    !%End

    if(parse_block(trim(what), blk) == 0) then
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

    !%Variable ReducedCoordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% This block gives the atomic coordinates relative to the real
    !% space unit cell. The format is the same as the
    !% <tt>Coordinates</tt> block.
    !%
    !% Note that in Octopus the origin of coordinates is in the center
    !% of the cell, so the coordinates inside the cell are in the
    !% range [-0.5, 0.5).
    !%
    !% This block cannot be used with the <tt>minimum</tt> box shape.
    !%End

    ! This is valid only for Coordinates.
    if(trim(what) == 'Coordinates') then
      if(parse_block('Reduced'//trim(what), blk) == 0) then
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
          do jdir = space%dim + 1, MAX_DIM
            gf%atom(ia)%x(jdir) = M_ZERO
          end do
          if(ncol == space%dim + 2) then
            call parse_block_logical(blk, ia - 1, space%dim + 1, gf%atom(ia)%move)
          else
            gf%atom(ia)%move = .true.
          end if
        end do

        call parse_block_end(blk)
      end if
    endif

    if(gf%source /= READ_COORDS_REDUCED) then
      ! adjust units
      do ia = 1, gf%n
        do jdir = space%dim + 1, MAX_DIM
          gf%atom(ia)%x(jdir) = M_ZERO
        end do

        if(gf%source == READ_COORDS_XYZ) then
          gf%atom(ia)%x = units_to_atomic(units_inp%length_xyz_file, gf%atom(ia)%x)
        else
          gf%atom(ia)%x = units_to_atomic(units_inp%length, gf%atom(ia)%x)
        end if
        
      end do
    end if

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

end module read_coords_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
