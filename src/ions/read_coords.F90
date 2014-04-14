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
    READ_COORDS_REDUCED  = 4,      &
    READ_COORDS_XSF      = 5

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

    integer :: ia, ncol, iunit, jdir, int_one
    type(block_t) :: blk
    character(len=256) :: str
    logical :: done
    FLOAT :: latvec(MAX_DIM, MAX_DIM)

    PUSH_SUB(read_coords_read)

    done = .false.

    !%Variable PDBCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If this variable is present, the program tries to read the atomic coordinates
    !% from the file specified by its value. The PDB (Protein Data Bank,
    !% <tt>http://www.rcsb.org/pdb/</tt>) format is quite complicated, and it goes 
    !% well beyond the scope of this manual. You can find a comprehensive
    !% description <a href='http://www.wwpdb.org/docs.html'>here</a>.
    !% From the plethora of instructions defined in the PDB standard, <tt>Octopus</tt>
    !% only reads two, <tt>ATOM</tt> and <tt>HETATOM</tt>. From these fields, it reads:
    !% <ul>
    !% <li> columns 13-16: The species; in fact <tt>Octopus</tt> only cares about the
    !% first letter - "CA" and "CB" will both refer to carbon - so elements whose
    !% chemical symbol has more than one letter cannot be represented in this way.
    !% So, if you want to run mercury (Hg), please use one of the other two methods
    !% to input the coordinates: <tt>XYZCoordinates</tt> or <tt>Coordinates</tt>.</li>
    !% <li> columns 18-21: The residue. If residue is <tt>QM</tt>, the atom is treated by quantum
    !% mechanics; otherwise it is simply treated as an external classical point charge.
    !% Its charge will be given by columns 61-65.</li>
    !% <li> columns 31-54: The Cartesian coordinates. The Fortran format is <tt>(3f8.3)</tt>.</li>
    !% <li> columns 61-65: Classical charge of the atom. The Fortran format is <tt>(f6.2)</tt>.</li>
    !% </ul>
    !%End

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
    !% NOTE: The coordinates are treated in the units specified by <tt>Units</tt> and/or <tt>UnitsInput</tt>.
    !%End

    if(parse_isdef(datasets_check('XYZ'//trim(what))) /= 0) then ! read an xyz file
      call check_duplicated(done)

      gf%source = READ_COORDS_XYZ
      ! no default, since we do not do this unless the input tag is present
      call parse_string(datasets_check('XYZ'//trim(what)), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1)

      iunit = io_open(str, status='old', action='read', is_tmp=.true.)
      read(iunit, *) gf%n

      if(gf%n <= 0) then
        write(message(1),'(a,i6)') "Invalid number of atoms ", gf%n
        call messages_fatal(1)
      endif

      read(iunit, *) ! skip comment line

      SAFE_ALLOCATE(gf%atom(1:gf%n))

      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if

    !%Variable XSFCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% Another option besides PDB and XYZ coordinates formats is XSF, as defined by the XCrySDen visualization
    !% program (http://www.xcrysden.org/doc/XSF.html). Specify the filename with this variable.
    !% is present. The XYZ format is very simple: The first line of the file has an integer
    !% indicating the number of atoms. The second can contain comments that are simply ignored by
    !% <tt>Octopus</tt>. Then there follows one line per atom, containing the chemical species and
    !% the Cartesian coordinates of the atom.
    !% NOTE: The coordinates are treated in the units specified by <tt>Units</tt> and/or <tt>UnitsInput</tt>.
    !%End

    if(parse_isdef(datasets_check('XSF'//trim(what))) /= 0) then ! read an xsf file
      call check_duplicated(done)

      gf%source = READ_COORDS_XSF
      ! no default, since we do not do this unless the input tag is present
      call parse_string(datasets_check('XSF'//trim(what)), '', str)

      message(1) = "Reading " // trim(what) // " from " // trim(str)
      call messages_info(1)

      iunit = io_open(str, status='old', action='read', is_tmp=.true.)

      read(iunit, *) str ! periodicity = 'CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE'; FIXME: or just 'ATOMS'
      select case(trim(str))
      case('CRYSTAL')
        gf%periodic_dim = 3
      case('SLAB')
        gf%periodic_dim = 2
      case('POLYMER')
        gf%periodic_dim = 1
      case('MOLECULE')
        gf%periodic_dim = 0
      case default
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of CRYSTAL/SLAB/POLYMER/MOLECULE.'
        call messages_fatal(1)
      end select

      read(iunit, *) str
      if(trim(str) /= 'PRIMVEC') then
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of "PRIMVEC".'
        call messages_warning(1)
      endif

      latvec(:,:) = M_ZERO
      do jdir = 1, space%dim
        read(iunit, *) latvec(1:space%dim, jdir)
        gf%lsize(jdir) = M_HALF * latvec(jdir, jdir)
        latvec(jdir, jdir) = M_ZERO
      enddo
      if(any(abs(latvec(1:space%dim, 1:space%dim)) > M_EPSILON)) then
        message(1) = 'XSF file has non-orthogonal lattice vectors. Only orthogonal is supported.'
        call messages_fatal(1)
      endif
      if(any(gf%lsize(1:space%dim) < M_EPSILON)) then
        message(1) = "XSF file must have positive lattice vectors."
        call messages_fatal(1)
      endif

      read(iunit, *) str
      if(trim(str) /= 'PRIMCOORD') then
        write(message(1),'(3a)') 'Line in file was "', trim(str), '" instead of "PRIMCOORD".'
        call messages_warning(1)
      endif

      read(iunit, *) gf%n, int_one
      if(gf%n <= 0) then
        write(message(1),'(a,i6)') "Invalid number of atoms ", gf%n
        call messages_fatal(1)
      endif
      if(int_one /= 1) then
        write(message(1),'(a,i6,a)') 'Number in file was ', int_one, ' instead of 1.'
        call messages_warning(1)
      endif
      SAFE_ALLOCATE(gf%atom(1:gf%n))

      ! TODO: add support for velocities as vectors here?
      do ia = 1, gf%n
        read(iunit,*) gf%atom(ia)%label, gf%atom(ia)%x(1:space%dim)
      end do

      call io_close(iunit)
    end if
    
    !%Variable Coordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% If neither <tt>XYZCoordinates</tt> nor <tt>PDBCoordinates</tt> was found, <tt>Octopus</tt>
    !% tries to read the coordinates for the atoms from the block <tt>Coordinates</tt>. The
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
    !% 
    !%End

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
