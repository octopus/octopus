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

module geometry_m
  use global_m
  use mpi_m
  use varinfo_m
  use messages_m
  use datasets_m
  use string_m
  use units_m
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use specie_m
  use xyz_file_m

  implicit none

  private
  public ::                &
    atom_t,                &
    atom_classical_t,      &
    geometry_t,            &
    geometry_init,         &
    geometry_init_xyz,     &
    geometry_init_vel,     &
    geometry_filter,       &
    geometry_init_species, &
    geometry_debug,        &
    geometry_nvnl,         &
    geometry_end,          &
    ion_ion_energy,        &
    kinetic_energy,        &
    geometry_dipole,       &
    geometry_min_distance, &
    cm_pos,                &
    cm_vel,                &
    atom_write_xyz,        &
    loadPDB,               &
    geometry_val_charge,   &
    geometry_grid_defaults

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
    character(len=20) :: sysname    ! the name of the system we are running

    integer :: natoms
    type(atom_t), pointer :: atom(:)

    integer :: ncatoms              ! For QM+MM calculations
    type(atom_classical_t), pointer :: catom(:)

    integer :: nspecies
    type(specie_t), pointer :: specie(:)

    logical :: only_user_def        ! Do we want to treat only user defined species?

    FLOAT :: eii, kinetic_energy    ! the ion-ion energy

    logical :: nlpp                 ! is any species having non-local pp
    logical :: nlcc                 ! is any species having non-local core corrections?

  end type geometry_t

contains

  ! ---------------------------------------------------------
  subroutine geometry_init(geo)
    type(geometry_t), intent(inout) :: geo

    call push_sub('geometry.geometry_init')

    ! initialize geometry
    call geometry_init_xyz(geo)
    call geometry_init_species(geo)
    call geometry_init_vel(geo)

    call pop_sub()
  end subroutine geometry_init


  ! ---------------------------------------------------------------
  ! initializes the xyz positions of the atoms in the structure geo
  subroutine geometry_init_xyz(geo)
    type(geometry_t), intent(inout) :: geo

    integer :: i
    type(xyz_file_info) :: xyz

    call push_sub('geometry.geometry_init_xyz')

    !%Variable SystemName
    !%Type string
    !%Default "system"
    !%Section System
    !%Description
    !% A string that identifies the current run. This parameter is seldomly used, but
    !% it is sometimes useful to have in the input file.
    !%End
    call loct_parse_string(check_inp('SystemName'), 'system', geo%sysname)

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


  !-----------------------------------------------------------------
  ! initializes the velocities of the atoms in the structure geo
  subroutine geometry_init_vel(geo)
    type(geometry_t), intent(inout) :: geo

    integer :: i, j
    FLOAT   :: x(MAX_DIM), temperature, sigma, kin1, kin2
    integer(POINTER_SIZE) :: random_gen_pointer
    type(xyz_file_info) :: xyz

    call push_sub('geometry.geometry_init_vel')

    !%Variable RandomVelocityTemp
    !%Type string
    !%Section System::Velocities
    !%Description
    !% If this variable is present, octopus will assign random velocities to the atoms 
    !% following a Bolzmann distribution with temperature given by RandomVelocityTemp.
    !%End

    ! we now load the velocities, either from the temperature, from the input, or from a file
    if(loct_parse_isdef(check_inp('RandomVelocityTemp')).ne.0) then
      call loct_ran_init(random_gen_pointer)
      call loct_parse_float(check_inp('RandomVelocityTemp'), M_ZERO, temperature)
      do i = 1, geo%natoms
        sigma = sqrt( P_Kb*temperature / geo%atom(i)%spec%weight )
        do j = 1, 3
          geo%atom(i)%v(j) = loct_ran_gaussian(random_gen_pointer, sigma)
        end do
      end do
      call loct_ran_end(random_gen_pointer)

      kin1 = kinetic_energy(geo)
      call cm_vel(geo, x)
      do i = 1, geo%natoms
        geo%atom(i)%v = geo%atom(i)%v - x
      end do
      kin2 = kinetic_energy(geo)
      do i = 1, geo%natoms
        geo%atom(i)%v(:) =  sqrt(kin1/kin2)*geo%atom(i)%v(:)
      end do

      write(message(1),'(a,f10.4,1x,a)') 'Info: Initial velocities ramdomly distributed with T =', &
        temperature, 'K'
      write(message(2),'(2x,a,f8.4,1x,a)') '<K>       =', &
        (kinetic_energy(geo)/geo%natoms)/units_out%energy%factor, &
        units_out%energy%abbrev
      write(message(3),'(2x,a,f8.4,1x,a)') '3/2 k_B T =', &
        (M_THREE/M_TWO)*P_Kb*temperature/units_out%energy%factor, &
        units_out%energy%abbrev
      call write_info(3)

    else
      !%Variable XYZVelocities
      !%Type string
      !%Section System::Velocities
      !%Description
      !% octopus will try to read the starting velocities of the atoms from the XYZ file 
      !% specified by the variable XYZVelocities.
      !% Note that you do not need to specify initial velocities if you are not going
      !% to perform ion dynamics; if you are going to allow the ions to move but the velocities
      !% are not specified, they are considered to be null.
      !%End

      !%Variable Velocities
      !%Type block
      !%Section System::Velocities
      !%Description
      !% If XYZVelocities is not present, octopus will try to fetch the initial 
      !% atomic velocities from this block. If this block is not present, octopus
      !% will reset the initial velocities to zero. The format of this block can be
      !% illustrated by this example:
      !%
      !% <tt>%Velocities
      !% <br>&nbsp;&nbsp;'C'  |      -1.7 | 0.0 | 0.0
      !% <br>&nbsp;&nbsp;'O'  | &nbsp;1.7 | 0.0 | 0.0
      !% <br>%</tt>
      !%
      !% It describes one Carbon and one Oxygen moving at the relative
      !% velocity of 3.4, velocity units.
      !%
      !% Note: It is important for the velocities to maintain the ordering 
      !% in which the species were defined in the coordinates specifications.
      !%End

      call xyz_file_init(xyz)
      call xyz_file_read('Velocities', xyz)
      if(xyz%file_type.ne.XYZ_FILE_ERR) then
        if(geo%natoms.ne.xyz%n) then
          write(message(1), '(a,i4,a,i4)') 'I need exactly ', geo%natoms, ' velocities, but I found ', xyz%n
          call write_fatal(1)
        end if

        ! copy information and adjust units
        do i = 1, geo%natoms
          geo%atom(i)%v = xyz%atom(i)%x * (units_inp%velocity%factor / units_inp%length%factor)
        end do
        call xyz_file_end(xyz)

      else
        do i = 1, geo%natoms
          geo%atom(i)%v = M_ZERO
        end do
      end if
    end if

    geo%kinetic_energy = kinetic_energy(geo)

    call pop_sub()
  end subroutine geometry_init_vel


  ! ---------------------------------------------------------
  subroutine geometry_filter(geo, gmax)
    type(geometry_t), intent(inout) :: geo
    FLOAT, intent(in) :: gmax
    integer :: i

    message(1) = 'Info: filtering the potentials.'
    call write_info(1)
    do i = 1, geo%nspecies
      if(.not.geo%specie(i)%local) call specie_filter(geo%specie(i), gmax)
    end do

  end subroutine geometry_filter


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
      geo%nlpp = (geo%nlpp.or.(.not.geo%specie(i)%local))
    end do

    call pop_sub()
  end subroutine geometry_init_species


  ! ---------------------------------------------------------
  subroutine geometry_debug(geo, dir)
    type(geometry_t), intent(in) :: geo
    character(len=*), intent(in) :: dir

    character(len=256) :: dirname
    integer :: i

    call push_sub('geometry.specie_debug')

    write(dirname, '(2a)') trim(dir), '/geometry'
    call io_mkdir(dirname)
    do i = 1, geo%nspecies
      call specie_debug(trim(dirname), geo%specie(i))
    end do

    call pop_sub()
  end subroutine geometry_debug


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

    if(associated(geo%atom)) then ! sanity check
      deallocate(geo%atom); nullify(geo%atom)
    end if

    if(geo%ncatoms > 0 .and. associated(geo%catom)) then
      deallocate(geo%catom); nullify(geo%catom)
    end if

    call specie_end(geo%nspecies, geo%specie)

  end subroutine geometry_end


  ! ---------------------------------------------------------
  ! Returns the number of non-local operator that should be defined.
  function geometry_nvnl(geo) result(res)
    type(geometry_t), intent(in) :: geo
    integer                      :: res

    type(specie_t), pointer :: s
    integer :: ia, l

    call push_sub('geometry.atom_nvnl')
    res = 0
    do ia = 1, geo%natoms
      s => geo%atom(ia)%spec
      if(s%local) cycle
      do l = 0, s%ps%l_max
        if(l == s%ps%l_loc) cycle
        res = res + 2*l + 1
      end do
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
  FLOAT function ion_ion_energy(geo)
    type(geometry_t), target, intent(in) :: geo

    type(specie_t), pointer :: s
    FLOAT :: r
    integer :: i, j

    ! Note that a possible jellium-jellium interaction (the case where more
    ! than one jellium species is present) is not properly calculated.
    ! But if only one jellium sphere is present, it correctly calculates its
    ! electrostatic energy. I do not know right now if there is a closed
    ! analytical expression for the electrostatic energy of a system of two
    ! uniformly charged spheres.
    ion_ion_energy = M_ZERO
    do i = 1, geo%natoms
      s => geo%atom(i)%spec
      if(s%type .eq. SPEC_JELLI) then
        ion_ion_energy = ion_ion_energy + (M_THREE/M_FIVE)*s%z_val**2/s%jradius
      end if
      do j = 1, i - 1
        r = sqrt(sum((geo%atom(i)%x - geo%atom(j)%x)**2))
        ion_ion_energy = ion_ion_energy + s%z_val*geo%atom(j)%spec%z_val/r
      end do
    end do

  end function ion_ion_energy


  ! ---------------------------------------------------------
  FLOAT function kinetic_energy(geo)
    type(geometry_t), intent(in) :: geo

    integer :: i

    kinetic_energy = M_ZERO
    do i = 1, geo%natoms
      kinetic_energy = kinetic_energy + &
        M_HALF*geo%atom(i)%spec%weight*sum(geo%atom(i)%v(:)**2)
    end do

  end function kinetic_energy


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

    rmin = huge(PRECISION)
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
  subroutine atom_write_xyz(dir, fname, geo)
    character(len=*),    intent(in) :: dir, fname
    type(geometry_t),    intent(in) :: geo

    integer i, iunit

    if( .not. mpi_grp_is_root(mpi_world)) return

    call io_mkdir(dir)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xyz', action='write')

    write(iunit, '(i4)') geo%natoms
    write(iunit, '(1x)')
    do i = 1, geo%natoms
      write(iunit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)

    if(geo%ncatoms > 0) then
      iunit = io_open(trim(dir)//'/'//trim(fname), action='write')
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

    def_h     =  huge(PRECISION)
    def_rsize = -huge(PRECISION)
    do i = 1, geo%nspecies
      def_h     = min(def_h,     geo%specie(i)%def_h)
      def_rsize = max(def_rsize, geo%specie(i)%def_rsize)
    end do

    call pop_sub()
  end subroutine geometry_grid_defaults

end module geometry_m
