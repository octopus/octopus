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

module geometry_m
  use c_pointer_m
  use datasets_m
  use global_m
  use io_m
  use loct_math_m
  use loct_parser_m
  use messages_m
  use multicomm_m
  use mpi_m
  use specie_m
  use string_m
  use units_m
  use varinfo_m
  use xyz_file_m

  implicit none

  private
  public ::                &
    atom_t,                &
    atom_classical_t,      &
    geometry_t,            &
    geometry_init,         &
    geometry_init_xyz,     &
    geometry_init_species, &
    geometry_partition,    &
    geometry_nvnl,         &
    geometry_end,          &
    geometry_dipole,       &
    geometry_min_distance, &
    cm_pos,                &
    cm_vel,                &
    atom_write_xyz,        &
    loadPDB,               &
    geometry_val_charge,   &
    geometry_grid_defaults,&
    assignment(=)

  type atom_t
    character(len=15) :: label
    type(specie_t), pointer :: spec              ! pointer to specie
    FLOAT :: x(MAX_DIM), v(MAX_DIM), f(MAX_DIM)  ! position/velocity/force of atom in real space
    logical :: move                              ! should I move this atom in the optimization mode
  end type atom_t

  type atom_classical_t
    character(len=15) :: label

    FLOAT :: x(MAX_DIM), v(MAX_DIM), f(MAX_DIM)
    FLOAT :: charge
  end type atom_classical_t

  type geometry_t
    integer :: natoms
    type(atom_t), pointer :: atom(:)

    integer :: ncatoms              ! For QM+MM calculations
    type(atom_classical_t), pointer :: catom(:)

    integer :: nspecies
    type(specie_t), pointer :: specie(:)

    logical :: only_user_def        ! Do we want to treat only user defined species?

    FLOAT :: kinetic_energy    ! the ion-ion energy

    logical :: nlpp                 ! is any species having non-local pp
    logical :: nlcc                 ! is any species having non-local core corrections?

    logical :: parallel_in_atoms

    type(mpi_grp_t) :: mpi_grp

    integer :: atoms_start
    integer :: atoms_end
    integer :: nlatoms
    
    integer, pointer :: atoms_range(:, :)
    integer, pointer :: atoms_num(:)
    integer, pointer :: atoms_node(:)
  end type geometry_t

  interface assignment (=)
    module procedure atom_copy
  end interface

contains

  ! ---------------------------------------------------------
  subroutine geometry_init(geo)
    type(geometry_t),            intent(inout) :: geo

    call push_sub('geometry.geometry_init')

    ! initialize geometry
    call geometry_init_xyz(geo)
    call geometry_init_species(geo)

    geo%parallel_in_atoms = .false.

    geo%atoms_start = 1
    geo%atoms_end   = geo%natoms
    geo%nlatoms     = geo%natoms

    nullify(geo%atoms_range, geo%atoms_num)

    call mpi_grp_init(geo%mpi_grp, -1)

    call pop_sub()
  end subroutine geometry_init


  ! ---------------------------------------------------------------
  ! initializes the xyz positions of the atoms in the structure geo
  subroutine geometry_init_xyz(geo)
    type(geometry_t), intent(inout) :: geo

    integer :: i
    type(xyz_file_info) :: xyz

    call push_sub('geometry.geometry_init_xyz')

    ! load positions of the atoms
    call xyz_file_init(xyz)

    !%Variable PDBCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If this variable is present, the program tries to read the atomic coordinates
    !% from the file specified by its value. The PDB (Protein Data Bank
    !% (http://www.rcsb.org/pdb/)) format is quite complicated, and it goes 
    !% well beyond the scope of this manual. You can find a comprehensive
    !% description in <a href='http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html'>here</a>.
    !% From the plethora of instructions defined in the PDB standard, octopus
    !% only reads two, "ATOM" and "HETATOM". From these fields, it reads:
    !% <ul>
    !% <li> columns 13-16: The specie; in fact "octopus" only cares about the
    !% first letter - "CA" and "CB" will both refer to Carbon - so elements whose
    !% chemical symbol has more than one letter can not be represented in this way.
    !% So, if you want to run mercury ("Hg") please use one of the other two methods_m
    !% to input the coordinates, "XYZCoordinates" or "Coordinates".</li>
    !% <li> columns 18-21: The residue. If residue is "QM", the atom is treated in Quantum
    !% Mechanics, otherwise it is simply treated as an external classical point charge.
    !% Its charge will be given by columns 61-65.</li>
    !% <li> columns 31-54: The Cartesian coordinates. The Fortran format is "(3f8.3)".</li>
    !% <li> columns 61-65: Classical charge of the atom. The Fortran format is "(f6.2)".</li>
    !% </ul>
    !%End

    !%Variable XYZCoordinates
    !%Type string
    !%Section System::Coordinates
    !%Description
    !% If "PDBCoordinates" is not present, the program reads the atomic coordinates from
    !% the XYZ file specified by the variable "XYZCoordinates" -- in case this variable
    !% is present. The XYZ format is very simple:  The first line of the file has an integer
    !% indicating the number of atoms. The second can contain comments that are simply ignored by
    !% "octopus". Then there follows one line per each atom, containing the chemical species and
    !% the Cartesian coordinates of the atom.
    !%End

    !%Variable Coordinates
    !%Type block
    !%Section System::Coordinates
    !%Description
    !% If neither a "XYZCoordinates" nor a "PDBCoordinates" was found, octopus
    !% tries to read the coordinates for the atoms from the block "Coordinates". The
    !% format is quite straightforward:
    !%
    !% <tt>%Coordinates
    !% <br>&nbsp;&nbsp;'C' |      -0.56415 | 0.0 | 0.0 | no
    !% <br>&nbsp;&nbsp;'O' | &nbsp;0.56415 | 0.0 | 0.0 | no
    !% <br>%</tt>
    !%
    !% The first line defines a Carbon atom at coordinates ("-0.56415", "0.0", "0.0"),
    !% that is _not_ allowed to move during dynamical simulations. The second line has
    !% a similar meaning. This block obviously defines a Carbon monoxide molecule, if the
    !% input units are AA. Note that in this way it is possible to fix some of the atoms (this
    !% is not possible when specifying the coordinates through a "PDBCoordinates" or
    !% "XYZCoordinates" file). It is always possible to fix _all_ atoms using the "MoveIons" directive.
    !%End
    call xyz_file_read('Coordinates', xyz)

    ! copy information from xyz to geo
    geo%natoms = xyz%n
    nullify(geo%atom)
    ALLOCATE(geo%atom(geo%natoms), geo%natoms)
    do i = 1, geo%natoms
      geo%atom(i)%label = xyz%atom(i)%label
      geo%atom(i)%x     = xyz%atom(i)%x
      geo%atom(i)%f     = M_ZERO
      if(iand(xyz%flags, XYZ_FLAGS_MOVE).ne.0) then
        geo%atom(i)%move = xyz%atom(i)%move
      else
        geo%atom(i)%move = .true.
      end if
    end do
    call xyz_file_end(xyz)

    ! load positions of the classical atoms, if any
    call xyz_file_init(xyz)
    nullify(geo%catom)
    geo%ncatoms = 0
    call xyz_file_read('Classical', xyz)
    if(xyz%file_type.ne.XYZ_FILE_ERR) then ! found classical atoms
      if(.not.iand(xyz%flags, XYZ_FLAGS_CHARGE).ne.0) then
        message(1) = "Need to know charge for the Classical atoms"
        message(2) = "Please use a .pdb"
        call write_fatal(2)
      end if
      geo%ncatoms = xyz%n
      write(message(1), '(a,i8)') 'Info: Number of classical atoms = ', geo%ncatoms
      call write_info(1)

      ALLOCATE(geo%catom(geo%ncatoms), geo%ncatoms)
      do i = 1, geo%ncatoms
        geo%catom(i)%label  = xyz%atom(i)%label
        geo%catom(i)%x      = xyz%atom(i)%x
        geo%catom(i)%v      = M_ZERO
        geo%catom(i)%f      = M_ZERO
        geo%catom(i)%charge = xyz%atom(i)%charge
      end do
      call xyz_file_end(xyz)
    end if

    if(geometry_atoms_are_too_close(geo)) then
      write(message(1), '(a)') "Some of the atoms seem to sit too close to each other."
      write(message(2), '(a)') "Please review your input files."
      call write_fatal(2)
    end if

    call pop_sub()
  end subroutine geometry_init_xyz

  ! ---------------------------------------------------------
  subroutine geometry_init_species(geo)
    type(geometry_t), intent(inout) :: geo

    integer :: i, j, k, ispin

    call push_sub('geometry.geometry_init_species')

    ! First, count the species
    geo%nspecies = 0
    atoms1:  do i = 1, geo%natoms
      do j = 1, i - 1
        if(trim(geo%atom(j)%label) == trim(geo%atom(i)%label)) cycle atoms1
      end do
      geo%nspecies = geo%nspecies + 1
    end do atoms1

    ! Allocate the species structure.
    ALLOCATE(geo%specie(geo%nspecies), geo%nspecies)

    ! Now, read the data.
    k = 0
    geo%only_user_def = .true.
    atoms2: do i = 1, geo%natoms
      do j = 1, i - 1
        if(trim(geo%atom(j)%label) == trim(geo%atom(i)%label)) cycle atoms2
      end do
      k = k + 1
      geo%specie(k)%label = geo%atom(j)%label
      geo%specie(k)%index = k
      call specie_read(geo%specie(k), trim(geo%specie(k)%label))
      geo%only_user_def = (geo%only_user_def .and. (geo%specie(k)%type==SPEC_USDEF))
    end do atoms2

    ! Reads the spin components. This is read here, as well as in states_init,
    ! to be able to pass it to the pseudopotential initializations subroutine.
    call loct_parse_int(check_inp('SpinComponents'), 1, ispin)
    if(.not.varinfo_valid_option('SpinComponents', ispin)) call input_error('SpinComponents')
    ispin = min(2, ispin)

    call messages_print_stress(stdout, "Species")
    do i = 1, geo%nspecies
      call specie_init(geo%specie(i), ispin)
    end do
    call messages_print_stress(stdout)

    !  assign species
    do i = 1, geo%natoms
      do j = 1, geo%nspecies
        if(trim(geo%atom(i)%label) == trim(geo%specie(j)%label)) then
          geo%atom(i)%spec => geo%specie(j)
          exit
        end if
      end do
    end do

    ! find out if we need non-local core corrections
    geo%nlcc = .false.
    geo%nlpp = .false.
    do i = 1, geo%nspecies
      geo%nlcc = (geo%nlcc.or.geo%specie(i)%nlcc)
      geo%nlpp = (geo%nlpp .or. specie_is_ps(geo%specie(i)))
    end do

    call pop_sub()
  end subroutine geometry_init_species

  subroutine geometry_partition(geo, mc)
    type(geometry_t),            intent(inout) :: geo
    type(multicomm_t),           intent(in)    :: mc

#ifdef HAVE_MPI
    integer :: size, rank, kk

    call push_sub('geometry.geometry_partition')

    ! for the atoms we use the parallelization in states

    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then

      geo%parallel_in_atoms = .true.

      call mpi_grp_init(geo%mpi_grp, mc%group_comm(P_STRATEGY_STATES))

      size = mc%group_sizes(P_STRATEGY_STATES)
      rank = mc%who_am_i(P_STRATEGY_STATES)

      ALLOCATE(geo%atoms_range(2, 0:size - 1), 2*size)
      ALLOCATE(geo%atoms_num(0:size - 1), size)
      ALLOCATE(geo%atoms_node(1:geo%natoms), geo%natoms)

      call multicomm_divide_range(geo%natoms, size, geo%atoms_range(1, :), geo%atoms_range(2, :), geo%atoms_num)

      message(1) = 'Info: Parallelization in atoms:'
        call write_info(1)
      do kk = 1, size
        write(message(1),'(a,i4,a,i4,a,i4)') 'Info: Node in states-group ', kk - 1, &
          ' will manage atoms', geo%atoms_range(1, kk - 1), " - ", geo%atoms_range(2, kk - 1)
        call write_info(1)
        if(rank .eq. kk - 1) then
          geo%atoms_start = geo%atoms_range(1, kk - 1)
          geo%atoms_end   = geo%atoms_range(2, kk - 1)
          geo%nlatoms     = geo%atoms_num(kk - 1)
        endif
        
        geo%atoms_node(geo%atoms_range(1, kk - 1):geo%atoms_range(2, kk - 1)) = kk -1

      end do
      
    end if

    call pop_sub()
#endif
  end subroutine geometry_partition

  ! ---------------------------------------------------------
  subroutine loadPDB(iunit, geo)
    integer,          intent(in)    :: iunit
    type(geometry_t), intent(inout) :: geo

    character(len=80) :: record
    character(len=6)  :: record_name
    character(len=4)  :: atm
    character(len=3)  :: res
    integer :: na, nca

    ! First count number of atoms
    rewind(iunit)
    geo%natoms = 0
    geo%ncatoms = 0
    do
      read(iunit, '(a80)', err=990, end=990) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(17x,a3)') res
        if(trim(res) == 'QM') then
          geo%natoms = geo%natoms + 1
        else
          geo%ncatoms = geo%ncatoms + 1
        end if
      end if
    end do
990 continue

    ALLOCATE(geo%atom(geo%natoms), geo%natoms)
    ALLOCATE(geo%catom(geo%ncatoms), geo%ncatoms)

    ! read in the data
    rewind(iunit)
    na = 1; nca = 1
    do
      read(iunit, '(a80)', err=991, end=991) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') atm, res
        call str_trim(atm)
        if(trim(res) == 'QM') then
          read(record, '(30x,3f8.3)') geo%atom(na)%x
          geo%atom(na)%label = atm(1:1)
          na = na + 1
        else
          geo%catom(nca)%label = atm
          read(record, '(30x,3f8.3,6x,f6.2)') geo%catom(nca)%x, geo%catom(nca)%charge
          nca = nca + 1
        end if
      end if
    end do
991 continue

  end subroutine loadPDB


  ! ---------------------------------------------------------
  subroutine geometry_end(geo)
    type(geometry_t), intent(inout) :: geo

    call push_sub('geometry.geometry_end')

    if(associated(geo%atoms_range)) then
      deallocate(geo%atoms_range, geo%atoms_num, geo%atoms_node)
      nullify(geo%atoms_range, geo%atoms_num, geo%atoms_node)
    end if

    if(associated(geo%atom)) then ! sanity check
      deallocate(geo%atom); nullify(geo%atom)
    end if

    if(geo%ncatoms > 0 .and. associated(geo%catom)) then
      deallocate(geo%catom); nullify(geo%catom)
    end if

    call specie_end(geo%nspecies, geo%specie)

    call pop_sub()
  end subroutine geometry_end


  ! ---------------------------------------------------------
  ! Returns the number of non-local operator that should be defined.
  function geometry_nvnl(geo) result(res)
    type(geometry_t), intent(in) :: geo
    integer                      :: res

    type(specie_t), pointer :: s
    integer :: ia, l

    call push_sub('geometry.geometry_nvnl')

    res = 0
    do ia = 1, geo%natoms
      s => geo%atom(ia)%spec
      if(specie_is_ps(s)) then
        do l = 0, s%ps%l_max
          if(l == s%ps%l_loc) cycle
          res = res + 2*l + 1
        end do
      end if
    end do

    call pop_sub()
  end function geometry_nvnl

  ! ---------------------------------------------------------
  ! This function returns .true. if two atoms are too close.
  logical function geometry_atoms_are_too_close(geo) result(l)
    type(geometry_t), intent(in) :: geo
    FLOAT :: r

    call geometry_min_distance(geo, r)
    l = (r < CNST(1.0e-5) .and. geo%natoms > 1)

  end function geometry_atoms_are_too_close

  ! ---------------------------------------------------------
  subroutine geometry_dipole(geo, dipole)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: dipole(MAX_DIM)

    integer :: i

    dipole = M_ZERO
    do i = 1, geo%natoms
      dipole(:) = dipole(:) + geo%atom(i)%spec%z_val*geo%atom(i)%x(:)
    end do

  end subroutine geometry_dipole


  ! ---------------------------------------------------------
  subroutine geometry_min_distance(geo, rmin)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: rmin

    integer :: i, j
    FLOAT :: r

    rmin = huge(REAL_PRECISION)
    do i = 1, geo%natoms
      do j = i+1, geo%natoms
        r = sqrt(sum((geo%atom(i)%x-geo%atom(j)%x)**2))
        if(r < rmin) then
          rmin = r
        end if
      end do
    end do

  end subroutine geometry_min_distance


  ! ---------------------------------------------------------
  subroutine cm_pos(geo, pos)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: pos(MAX_DIM)

    FLOAT :: m
    integer :: i

    pos = M_ZERO; m = M_ZERO
    do i = 1, geo%natoms
      pos = pos + geo%atom(i)%spec%weight*geo%atom(i)%x
      m = m + geo%atom(i)%spec%weight
    end do
    pos = pos/m
  end subroutine cm_pos


  ! ---------------------------------------------------------
  subroutine cm_vel(geo, vel)
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: vel(MAX_DIM)

    FLOAT :: m
    integer :: i

    vel = M_ZERO; m = M_ZERO
    do i = 1, geo%natoms
      vel = vel + geo%atom(i)%spec%weight*geo%atom(i)%v
      m = m + geo%atom(i)%spec%weight
    end do
    vel = vel/m
  end subroutine cm_vel


  ! ---------------------------------------------------------
  subroutine atom_write_xyz(dir, fname, geo, append, comment)
    character(len=*),    intent(in) :: dir, fname
    type(geometry_t),    intent(in) :: geo
    logical,             intent(in), optional :: append
    character(len=*),    intent(in), optional :: comment

    integer i, iunit
    character(len=6) position

    if( .not. mpi_grp_is_root(mpi_world)) return

    call io_mkdir(dir)
    position = 'asis'
    if(present(append)) then
      if(append) position = 'append'
    end if
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xyz', action='write', position=position)

    write(iunit, '(i4)') geo%natoms
    if (present(comment)) then
      write(iunit, '(1x,a)') comment
    else
      write(iunit, '(1x)')
    endif
    do i = 1, geo%natoms
      write(iunit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)

    if(geo%ncatoms > 0) then
      iunit = io_open(trim(dir)//'/'//trim(fname)//'_classical.xyz', action='write', position=position)
      write(iunit, '(i4)') geo%ncatoms
      write(iunit, '(1x)')
      do i = 1, geo%ncatoms
        write(iunit, '(6x,a1,2x,3f12.6,a,f12.6)') &
           geo%catom(i)%label(1:1), geo%catom(i)%x(:)/units_out%length%factor, &
           " # ", geo%catom(i)%charge
      end do
      call io_close(iunit)
    end if

  end subroutine atom_write_xyz


  ! ---------------------------------------------------------
  subroutine geometry_val_charge(geo, val_charge)
    type(geometry_t), intent(in) :: geo
    FLOAT,           intent(out) :: val_charge

    integer :: i

    call push_sub('geometry.geometry_val_charge')

    val_charge = M_ZERO
    do i = 1, geo%natoms
      val_charge = val_charge - geo%atom(i)%spec%Z_val
    end do

    call pop_sub()
  end subroutine geometry_val_charge


  ! ---------------------------------------------------------
  subroutine geometry_grid_defaults(geo, def_h, def_rsize)
    type(geometry_t), intent(in) :: geo
    FLOAT,           intent(out) :: def_h, def_rsize

    integer :: i

    call push_sub('geometry.geometry_grid_defaults')

    def_h     =  huge(REAL_PRECISION)
    def_rsize = -huge(REAL_PRECISION)
    do i = 1, geo%nspecies
      def_h     = min(def_h,     geo%specie(i)%def_h)
      def_rsize = max(def_rsize, geo%specie(i)%def_rsize)
    end do

    call pop_sub()
  end subroutine geometry_grid_defaults


  !--------------------------------------------------------------
  subroutine atom_copy(aout, ain)
    type(atom_t), intent(out) :: aout
    type(atom_t), intent(in)  :: ain
    aout%label = ain%label
    aout%spec  = ain%spec
    aout%x     = ain%x
    aout%v     = ain%v
    aout%f     = ain%f
    aout%move  = ain%move
  end subroutine atom_copy

end module geometry_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
