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

#include "global.h"

module geometry
  use global
  use lib_oct_parser
  use lib_oct
  use io
  use specie
  use xyz_file

  implicit none

  type atom_type
    character(len=10) :: label
    type(specie_type), pointer :: spec ! pointer to specie

    FLOAT :: x(3), v(3), f(3) ! position/velocity/force of atom in real space

    logical :: move              ! should I move this atom in the optimization mode
  end type atom_type

  type atom_classical_type
    FLOAT :: x(3), v(3), f(3)
    FLOAT :: charge

    character(len=4) :: label
  end type atom_classical_type

  type geometry_type
    character(len=20) :: sysname    ! the name of the system we are running

    integer :: natoms
    type(atom_type), pointer :: atom(:)

    integer :: ncatoms  ! For QM+MM calculations
    type(atom_classical_type), pointer :: catom(:)

    integer :: nspecies
    type(specie_type), pointer :: specie(:)

    FLOAT :: eii, kinetic_energy ! the ion-ion energy

    logical :: nlpp    ! is any species having non-local pp
    logical :: nlcc    ! is any species having non-local core corrections?

  end type geometry_type

contains

subroutine geometry_init_xyz(geo)
  type(geometry_type), intent(inout) :: geo

  integer :: iunit, i, j
  integer(POINTER_SIZE) :: random_gen_pointer
  character(len=80) :: str, label
  FLOAT :: temperature, sigma, x(3), kin1, kin2
  type(xyz_file_info) :: xyz

  call push_sub('geometry_init')

  ! get the name of the system
  call loct_parse_string('SystemName', 'system', geo%sysname)

  ! load positions of the atoms
  call xyz_file_init(xyz)
  call xyz_file_read('Coordinates', xyz)

  ! copy information from xyz to geo
  geo%natoms = xyz%n
  allocate(geo%atom(geo%natoms))
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
  call xyz_file_read('Classical', xyz)
  if(xyz%file_type.ne.XYZ_FILE_ERR) then ! found classical atoms
    if(.not.iand(xyz%flags, XYZ_FLAGS_CHARGE).ne.0) then
      message(1) = "Need to know charge for the Classical atoms"
      message(2) = "Please use a .pdb"
      call write_fatal(2)
    end if
    geo%ncatoms = xyz%n
    allocate(geo%catom(geo%ncatoms))
    do i = 1, geo%ncatoms
      geo%catom(i)%label  = xyz%atom(i)%label
      geo%catom(i)%x      = xyz%atom(i)%x
      geo%catom(i)%v      = M_ZERO
      geo%catom(i)%f      = M_ZERO
      geo%catom(i)%charge = xyz%atom(i)%charge
    end do
    call xyz_file_end(xyz)
  end if
    
  ! we now load the velocities, either from the temperature, from the input, or from a file
  if(loct_parse_isdef("RandomVelocityTemp").ne.0) then
    call loct_ran_init(random_gen_pointer)
    call loct_parse_float("RandomVelocityTemp", M_ZERO, temperature)
    do i = 1, geo%natoms
       sigma = sqrt( P_Kb*temperature / geo%atom(i)%spec%weight )
       do j = 1, 3
          geo%atom(i)%v(j) = loct_ran_gaussian(random_gen_pointer, sigma)
       enddo
    enddo
    call loct_ran_end(random_gen_pointer)
    kin1 = kinetic_energy(geo)
    call cm_vel(geo, x)
    do i = 1, geo%natoms
       geo%atom(i)%v = geo%atom(i)%v - x
    enddo
    kin2 = kinetic_energy(geo)
    do i = 1, geo%natoms
       geo%atom(i)%v(:) =  sqrt(kin1/kin2)*geo%atom(i)%v(:)
    enddo

    write(message(1),'(a,f10.4,1x,a)') 'Info: Initial velocities ramdomly distributed with T =', & 
                                   temperature, 'K'
    write(message(2),'(a,f8.4,1x,a)') 'Info: <K>       =', &
                                   (kinetic_energy(geo)/geo%natoms)/units_out%energy%factor, &
                                   units_out%energy%abbrev
    write(message(3),'(a,f8.4,1x,a)') 'Info: 3/2 k_B T =', &
                                   (M_THREE/M_TWO)*P_Kb*temperature/units_out%energy%factor, &
                                   units_out%energy%abbrev
    write(message(4),'(a)')
    call write_info(4)

  else
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

  call pop_sub()

end subroutine geometry_init_xyz

subroutine geometry_init_species(geo, val_charge_, def_h_, def_rsize_)
  type(geometry_type), intent(inout) :: geo
  FLOAT, optional, intent(out)   :: val_charge_ ! the valence charge
  FLOAT, optional, intent(out)   :: def_h_      ! the default mesh spacing
  FLOAT, optional, intent(out)   :: def_rsize_  ! the default size of the minimum box

  type specie_list_type
    character(len=10) :: label
    integer :: location ! 1 <= Species block; 2 <= defaults file
    integer :: line     ! in which line of the file the specie is defined
    logical :: used
  end type specie_list_type

  FLOAT :: val_charge, def_h, def_rsize
  type(specie_list_type), allocatable :: spec_list(:)
  integer :: i, j, n_spec, ispin
  logical :: ok

  call push_sub('geometry_init_species')

  call read_species_label() ! which species do we know about

  ! let us see which species that we are going to use
  do i = 1, geo%natoms
    ok = .false.
    do j = 1, n_spec
      if(trim(geo%atom(i)%label) == trim(spec_list(j)%label)) then
        spec_list(j)%used = .true.
        ok = .true.
        exit
      end if
    end do
    if(.not.ok) then
      write(message(1),'(aaa)') "Specie '", geo%atom(i)%label, "' was not found"
      call write_fatal(1)
    end if
  end do

  ! how many are they
  geo%nspecies = 0
  do i = 1, n_spec
    if(spec_list(i)%used) geo%nspecies = geo%nspecies + 1
  end do

  ! we allocate the necessary memory
  allocate(geo%specie(geo%nspecies))

  ! Reads the spin components. This is read here, as well as in states_init,
  ! to be able to pass it to the pseudopotential initializations subroutine.
  call loct_parse_int('SpinComponents', 1, ispin)
  if (ispin < 1 .or. ispin > 3) then
    write(message(1),'(a,i4,a)') "Input: '", ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if
  ispin = min(2, ispin)
  
  ! we now load the individual species
  def_h     = huge(PRECISION)
  def_rsize = huge(PRECISION)
  j = 0
  do i = 1, n_spec
    if(.not.spec_list(i)%used) cycle

    j = j + 1
    geo%specie(j)%label = spec_list(i)%label
    geo%specie(j)%index = j
    call specie_init(geo%specie(j), spec_list(i)%location, spec_list(i)%line, ispin)

    def_h     = min(def_h,     geo%specie(j)%def_h)
    def_rsize = min(def_rsize, geo%specie(j)%def_rsize)
  end do

  !  assign species and find total charge of the system
  val_charge = M_ZERO
  do i = 1, geo%natoms
    do j = 1, geo%nspecies
      if(trim(geo%atom(i)%label) == trim(geo%specie(j)%label)) then
        geo%atom(i)%spec => geo%specie(j)
        exit
      end if
    end do
    val_charge = val_charge - geo%atom(i)%spec%Z_val
  end do
  
  ! find out if we need non-local core corrections
  geo%nlcc = .false.
  geo%nlpp = .false.
  do i = 1, geo%nspecies
    geo%nlcc = (geo%nlcc.or.geo%specie(i)%nlcc)
    geo%nlpp = (geo%nlpp.or.(.not.geo%specie(i)%local))
  end do

  ! return values
  if(present(val_charge_)) val_charge_ = val_charge
  if(present(def_h_))      def_h_      = def_h
  if(present(def_rsize_))  def_rsize_  = def_rsize

  ! clean up
  deallocate(spec_list)
  call pop_sub()

contains
  subroutine read_species_label()
    integer :: i, n_spec_def, n_spec_block
    character(len=256) :: fname
    integer :: iunit

    write(fname, '(2a)') trim(conf%share), "/PP/defaults"

    ! how many species do we have defined in defaults file
    n_spec_def = loct_number_of_lines(fname)
    if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

    ! is the block Species defined
    n_spec_block = loct_parse_block_n("Species")
    
    ! total number of species that we have available
    n_spec = n_spec_def + n_spec_block
    if(n_spec < 0) then
      message(1) = "No species defined either in Specie block or in 'default' file"
      call write_fatal(1)
    end if
    
    ! now we load the species
    allocate(spec_list(n_spec))
    do i = 1, n_spec_block
      call loct_parse_block_string("Species", i-1, 0, spec_list(i)%label)
      spec_list(i)%location = 1
      spec_list(i)%line = i
      spec_list(i)%used = .false.
    end do

    if(n_spec_def < 1) return ! nothing else to do
    call io_assign(iunit)
    open(unit=iunit, file=fname, action='read')
    read(iunit,*)
    do i = n_spec_block+1, n_spec
      read(iunit,*) spec_list(i)%label
      spec_list(i)%location = 2
      spec_list(i)%line = i - n_spec_block
      spec_list(i)%used = .false.
    end do
    call io_close(iunit)

  end subroutine read_species_label

  integer function get_specie(label) result(n)
    character(len=*) :: label

    integer :: j

    n = -1
    do j = 1, geo%nspecies
      if( trim(label) == trim(geo%specie(j)%label) ) then
        n = j
        return
      end if
    end do
    
  end function get_specie
  
end subroutine geometry_init_species

subroutine loadPDB(iunit, geo)
  integer,             intent(in)    :: iunit
  type(geometry_type), intent(inout) :: geo
 
  character(len=80) :: record
  character(len=6) :: record_name
  character(len=4) :: atm
  character(len=3) :: res
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

  allocate(geo%atom(geo%natoms), geo%catom(geo%ncatoms))

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
          !a(na)%spec => s(get_specie(atm(1:1)))
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

subroutine geometry_end(geo)
  type(geometry_type), intent(inout) :: geo

  if(associated(geo%atom)) then ! sanity check
    deallocate(geo%atom); nullify(geo%atom)
  end if

  if(geo%ncatoms > 0 .and. associated(geo%catom)) then
    deallocate(geo%catom); nullify(geo%catom)
  end if

  call specie_end(geo%nspecies, geo%specie)

end subroutine geometry_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the number of non-local operator that should be defined.
function geometry_nvnl(geo) result(res)
  type(geometry_type), intent(IN) :: geo
  integer                         :: res

  type(specie_type), pointer :: s
  integer :: ia

  call push_sub('atom_nvnl')
  res = 0
  do ia = 1, geo%natoms
     s => geo%atom(ia)%spec
     if(s%local) cycle
     res = res + (s%ps%l_max+1)**2
  enddo

  call pop_sub()
end function geometry_nvnl

FLOAT function ion_ion_energy(geo)
  type(geometry_type), intent(IN) :: geo
  
  FLOAT :: r
  integer :: i, j

  ! calculate the ion-ion energy
  ion_ion_energy = M_ZERO
  do i = 2, geo%natoms
    do j = 1, i - 1
      r = sqrt(sum((geo%atom(i)%x - geo%atom(j)%x)**2))
      ion_ion_energy = ion_ion_energy + geo%atom(i)%spec%Z_val*geo%atom(j)%spec%Z_val/r
    end do
  end do

end function ion_ion_energy

FLOAT function kinetic_energy(geo)
  type(geometry_type), intent(IN) :: geo

  integer :: i

  kinetic_energy = M_ZERO
  do i = 1, geo%natoms
    kinetic_energy = kinetic_energy + &
                     M_HALF*geo%atom(i)%spec%weight*sum(geo%atom(i)%v(:)**2)
  end do

end function kinetic_energy

subroutine geometry_dipole(geo, dipole)
  type(geometry_type), intent(IN)  :: geo
  FLOAT,               intent(out) :: dipole(3)

  integer :: i

  dipole = M_ZERO
  do i = 1, geo%natoms
     dipole(1:conf%dim) = dipole(1:conf%dim) + geo%atom(i)%spec%z_val*geo%atom(i)%x(1:conf%dim)
  end do

end subroutine geometry_dipole

subroutine cm_pos(geo, pos)
  type(geometry_type), intent(IN)  :: geo
  FLOAT,               intent(out) :: pos(3)

  FLOAT :: m
  integer :: i

  pos = M_ZERO; m = M_ZERO
  do i = 1, geo%natoms
     pos = pos + geo%atom(i)%spec%weight*geo%atom(i)%x
     m = m + geo%atom(i)%spec%weight
  enddo
  pos = pos/m
end subroutine cm_pos

subroutine cm_vel(geo, vel)
  type(geometry_type), intent(IN)  :: geo
  FLOAT,               intent(out) :: vel(3)

  FLOAT :: m
  integer :: i

  vel = M_ZERO; m = M_ZERO
  do i = 1, geo%natoms
     vel = vel + geo%atom(i)%spec%weight*geo%atom(i)%v
     m = m + geo%atom(i)%spec%weight
  enddo
  vel = vel/m
end subroutine cm_vel

subroutine atom_write_xyz(dir, fname, geo)
  character(len=*),    intent(in) :: dir, fname
  type(geometry_type), intent(IN) :: geo
  
  integer i, iunit
  
#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    call loct_mkdir(trim(dir))

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // '.xyz', status='unknown')
    write(iunit, '(i4)') geo%natoms
    write(iunit, '(1x)')
    do i = 1, geo%natoms
      write(iunit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)
    
    if(geo%ncatoms > 0) then
      call io_assign(iunit)
      open(iunit, file=trim(dir) // "/" // trim(fname) // '_classical.xyz', status='unknown')
      write(iunit, '(i4)') geo%ncatoms
      write(iunit, '(1x)')
      do i = 1, geo%ncatoms
        write(iunit, '(6x,a1,2x,3f12.6,a,f12.6)') &
             geo%catom(i)%label(1:1), geo%catom(i)%x(:)/units_out%length%factor, &
             " # ", geo%catom(i)%charge
      end do
      call io_close(iunit)
    end if

#ifdef HAVE_MPI
  end if
#endif
  
end subroutine atom_write_xyz

end module geometry
