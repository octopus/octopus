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
  use mesh
  use mesh_function
  use atom
  use specie

  implicit none

  type geometry_type
    character(len=20)          :: sysname    ! the name of the system we are running

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

subroutine geometry_init(geo, val_charge, no_species_init)
  type(geometry_type), intent(inout) :: geo
  FLOAT, intent(out), optional :: val_charge
  logical, intent(in), optional :: no_species_init

  integer :: iunit, i, j
  integer(POINTER_SIZE) :: random_gen_pointer
  character(len=80) :: str, label
  logical :: l
  FLOAT :: temperature, sigma, x(3), kin1, kin2

  call push_sub('geometry_init')

  ! get the name of the system
  call loct_parse_string('SystemName', 'system', geo%sysname)

  ! this is not very nice, but it is needed for the xyzanim
  if (present(no_species_init)) then
    if (.not.no_species_init) geo%nspecies = specie_init(geo%specie)
  else
    geo%nspecies = specie_init(geo%specie)
  end if

  if(conf%dim == 3.and.loct_parse_isdef("PDBCoordinates").ne.0) then
    call loct_parse_string('PDBCoordinates', 'coords.pdb', label)

    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    call loadPDB(iunit, geo)
    do i = 1, geo%natoms
      geo%atom(i)%spec => geo%specie(get_specie(geo%atom(i)%label))
    enddo
    call io_close(iunit)

  else
    ! we now load the positions, either from the input, or from a file
    if(conf%dim == 3.and.loct_parse_isdef("XYZCoordinates").ne.0) then ! read a xyz file
      call loct_parse_string('XYZCoordinates', 'coords.xyz', label)

      call io_assign(iunit)
      open(iunit, status='unknown', file=trim(label))
      read(iunit, *) geo%natoms
      read(iunit, *) ! skip comment line
      allocate(geo%atom(geo%natoms))
      nullify(geo%catom); geo%ncatoms = 0;

      do i = 1, geo%natoms
        read(iunit,*) geo%atom(i)%label, geo%atom(i)%x(:)
        geo%atom(i)%move = .true.
        geo%atom(i)%spec => geo%specie(get_specie(geo%atom(i)%label))
      end do

      call io_close(iunit)
    else
      str = "Coordinates"
      geo%natoms = loct_parse_block_n(str)
      if(geo%natoms <= 0) then
        message(1) = "Input: Coordinates block not specified"
        message(2) = '% Coordinates'
        message(3) = '  specie  x  y  z  move'
        message(4) = '%'
        call write_fatal(4)
      end if
      
      allocate(geo%atom(geo%natoms))
      nullify(geo%catom); geo%ncatoms = 0
      do i = 1, geo%natoms
        call loct_parse_block_string (str, i-1, 0, geo%atom(i)%label)
        geo%atom(i)%spec => geo%specie(get_specie(geo%atom(i)%label))
        call loct_parse_block_float  (str, i-1, 1, geo%atom(i)%x(1))
        call loct_parse_block_float  (str, i-1, 2, geo%atom(i)%x(2))
        call loct_parse_block_float  (str, i-1, 3, geo%atom(i)%x(3))
        call loct_parse_block_logical(str, i-1, 4, geo%atom(i)%move)
      end do
    end if
  end if

  ! units conversion
  do i = 1, geo%natoms
    geo%atom(i)%x = geo%atom(i)%x * units_inp%length%factor
  end do
  do i = 1, geo%ncatoms
    geo%catom(i)%x = geo%catom(i)%x * units_inp%length%factor
  end do

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

  elseif(loct_parse_isdef("XYZVelocities").ne.0 .and. conf%dim==3) then ! read a xyz file
    call io_assign(iunit)
    call loct_parse_string('XYZVelocities', 'velocities.xyz', label)
    open(iunit, status='unknown', file=trim(label))
      
    read(iunit, *)
    read(iunit, *) ! skip comment line
      
    do i = 1, geo%natoms
      read(iunit,*) label, geo%atom(i)%v
      geo%atom(i)%v = units_inp%velocity%factor * geo%atom(i)%v !units conversion
    end do

    call io_close(iunit)
  else 
    str = "Velocities"
    if(loct_parse_isdef(str).ne.0) then
      do i = 1, geo%natoms
        call loct_parse_block_float(str, i-1, 1, geo%atom(i)%v(1))
        call loct_parse_block_float(str, i-1, 2, geo%atom(i)%v(2))
        call loct_parse_block_float(str, i-1, 3, geo%atom(i)%v(3))
        geo%atom(i)%v = geo%atom(i)%v * units_inp%velocity%factor
      end do
    else
      geo%atom%v(1) = M_ZERO
      geo%atom%v(2) = M_ZERO
      geo%atom%v(3) = M_ZERO
    end if
  end if

  !  find total charge of the system
  val_charge = M_ZERO
  do i = 1, geo%natoms
    val_charge = val_charge - geo%atom(i)%spec%Z_val
  enddo

  ! find out if we need non-local core corrections
  geo%nlcc = .false.
  geo%nlpp = .false.
  do i = 1, geo%nspecies
    geo%nlcc = (geo%nlcc.or.geo%specie(i)%nlcc)
    geo%nlpp = (geo%nlcc.or.(.not.geo%specie(i)%local))
  end do

  call pop_sub()

contains

  integer function get_specie(label)
    character(len=*) :: label

    integer :: j
    logical :: l

    l = .false.
    do j = 1, geo%nspecies
      if( trim(label) == trim(geo%specie(j)%label) ) then
        l = .true.
        get_specie = j
        return
      end if
    end do
    
    if(.not. l) then
      message(1) = "Specie '"+trim(label)+"' not found"
      call write_fatal(1)
    end if
    
  end function get_specie

end subroutine geometry_init

subroutine loadPDB(iunit, geo)
  integer, intent(in) :: iunit
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
  integer :: res
  type(geometry_type), intent(in) :: geo

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
  type(geometry_type), intent(in) :: geo
  
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
  type(geometry_type), intent(in) :: geo

  integer :: i

  kinetic_energy = M_ZERO
  do i = 1, geo%natoms
    kinetic_energy = kinetic_energy + &
                     M_HALF*geo%atom(i)%spec%weight*sum(geo%atom(i)%v(:)**2)
  end do

end function kinetic_energy

subroutine geometry_dipole(geo, dipole)
  type(geometry_type), intent(in) :: geo
  FLOAT, intent(out) :: dipole(3)

  integer :: i

  dipole = M_ZERO
  do i = 1, geo%natoms
     dipole(1:conf%dim) = dipole(1:conf%dim) + geo%atom(i)%spec%z_val*geo%atom(i)%x(1:conf%dim)
  end do

end subroutine geometry_dipole

subroutine cm_pos(geo, pos)
  type(geometry_type), intent(in) :: geo
  FLOAT, intent(out) :: pos(3)

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
  type(geometry_type), intent(in) :: geo
  FLOAT, intent(out) :: vel(3)

  FLOAT :: m
  integer :: i

  vel = M_ZERO; m = M_ZERO
  do i = 1, geo%natoms
     vel = vel + geo%atom(i)%spec%weight*geo%atom(i)%v
     m = m + geo%atom(i)%spec%weight
  enddo
  vel = vel/m
end subroutine cm_vel

! builds a density which is the sum of the atomic densities
subroutine lcao_dens(m, geo, qtot, nspin, spin_channels, rho)
  type(mesh_type),     intent(in)  :: m
  type(geometry_type), intent(in)  :: geo
  FLOAT,               intent(in)  :: qtot  ! the total charge of the system
  integer,             intent(in)  :: nspin, spin_channels
  FLOAT,               intent(out) :: rho(m%np, nspin)

  integer :: ia, is
  FLOAT :: r
  type(atom_type),   pointer :: a

  call push_sub('lcao_dens')

  rho = M_ZERO
  do ia = 1, geo%natoms
    rho(1:m%np, 1:spin_channels) = rho(1:m%np, 1:spin_channels) &
                                   + atom_density(m, geo%atom(ia), spin_channels)
  end do

  ! we now renormalize the density (necessary if we have a charged system)
  r = M_ZERO
  do is = 1, spin_channels
    r = r + dmf_integrate(m, rho(:, is))
  end do
  write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
  call write_info(1)
  r = qtot/r
  rho = r*rho
  r = M_ZERO
  do is = 1, spin_channels
    r = r + dmf_integrate(m, rho(:, is))
  end do
  write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', r
  call write_info(1)

  call pop_sub()

end subroutine lcao_dens

subroutine atom_write_xyz(dir, fname, geo)
  character(len=*), intent(in)  :: dir, fname
  type(geometry_type), intent(in) :: geo
  
  integer i, iunit
  
#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    call loct_mkdir(trim(dir))

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+'.xyz', status='unknown')
    write(iunit, '(i4)') geo%natoms
    write(iunit, '(1x)')
    do i = 1, geo%natoms
      write(iunit, '(6x,a,2x,3f12.6)') geo%atom(i)%label, geo%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)
    
    if(geo%ncatoms > 0) then
      call io_assign(iunit)
      open(iunit, file=trim(dir)+"/"+trim(fname)+'_classical.xyz', status='unknown')
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

subroutine geometry_adjust(geo)
  type(geometry_type), intent(inout) :: geo

  FLOAT :: center(3), from(3), from2(3), to(3)
  character(len=80) :: str

  ! is there something to do
  if(geo%natoms <= 1) return

  ! recenter
  call find_center(center)
  call translate(-center)

  ! get to axis
  str = "MainAxis"
  if(loct_parse_isdef(str) .ne. 0) then
    call loct_parse_block_float(str, 0, 0, to(1))
    call loct_parse_block_float(str, 0, 1, to(2))
    call loct_parse_block_float(str, 0, 2, to(3))
  else
    to(1) = M_ZERO; to(2) = M_ZERO; to(3) = M_ONE
  end if
  to = to / sqrt(sum(to**2))

  ! rotate to main axis
  call find_axis(from, from2)
  call rotate(from, from2, to)

  ! recenter
  call find_center(center)
  call translate(-center)
  
contains

  subroutine find_center(x)
    FLOAT, intent(out) :: x(3)

    FLOAT :: xmin(3), xmax(3)
    integer  :: i, j

    xmin =  CNST(1e10)
    xmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, 3
        if(geo%atom(i)%x(j) > xmax(j)) xmax(j) = geo%atom(i)%x(j)
        if(geo%atom(i)%x(j) < xmin(j)) xmin(j) = geo%atom(i)%x(j)
      end do
    end do

    x = (xmax + xmin)/M_TWO
  end subroutine find_center

  subroutine translate(x)
    FLOAT, intent(in) :: x(3)

    integer  :: i

    do i = 1, geo%natoms
      geo%atom(i)%x = geo%atom(i)%x + x
    end do
    do i = 1, geo%ncatoms
      geo%catom(i)%x = geo%catom(i)%x + x
    end do
  end subroutine translate

  subroutine find_axis(x, x2)
    FLOAT, intent(out) :: x(3), x2(3)

    integer  :: i, j
    FLOAT :: rmax, r, r2

    ! first get the further apart atoms
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      do j = 1, geo%natoms/2 + 1
        r = sqrt(sum((geo%atom(i)%x-geo%atom(j)%x)**2))
        if(r > rmax) then
          rmax = r
          x = geo%atom(i)%x - geo%atom(j)%x
        end if
      end do
    end do
    x  = x /sqrt(sum(x**2))

    ! now let us find out what is the second most important axis
    rmax = -CNST(1e10)
    do i = 1, geo%natoms
      r2 = sum(x * geo%atom(i)%x)
      r = sqrt(sum((geo%atom(i)%x - r2*x)**2))
      if(r > rmax) then
        rmax = r
        x2 = geo%atom(i)%x - r2*x
      end if
    end do
    if(sum(x2**2) == M_ZERO) then ! linear molecule
      if(x(1) == M_ZERO) then
        x2(1) = x(1); x2(2) = -x(3); x2(3) = x(2)
      else if(x(1) == M_ZERO) then
        x2(2) = x(2); x2(1) = -x(3); x2(3) = x(1)
      else
        x2(3) = x(3); x2(1) = -x(2); x2(2) = x(1)
      end if
    end if
    x2 = x2/sqrt(sum(x2**2))

  end subroutine find_axis

  subroutine rotate(from, from2, to)
    FLOAT, intent(in) :: from(3), from2(3), to(3) ! assumed to be normalize

    integer :: i
    FLOAT :: m1(3,3), m2(3,3), m3(3,3), f2(3), per(3)
    FLOAT :: alpha, r

    ! initialize matrices
    m1 = M_ZERO; m1(1,1) = M_ONE; m1(2,2) = M_ONE; m1(3,3) = M_ONE

    ! rotate the to axis to the z axis
    if(to(2).ne.M_ZERO) then
      alpha = atan2(to(2), to(1))
      call rotate_z(m1, alpha)
    end if
    alpha = atan2(sqrt(to(1)**2 + to(2)**2), to(3))
    call rotate_y(m1, -alpha)

    ! get perpendicular to z and from
    f2 = matmul(m1, from)
    per(1) = -f2(2)
    per(2) =  f2(1)
    per(3) = M_ZERO
    r = sqrt(sum(per**2))
    if(r > M_ZERO) then
      per = per/r
    else
      per(2) = M_ONE
    end if

    ! rotate perpendicular axis to the y axis
    m2 = M_ZERO; m2(1,1) = M_ONE; m2(2,2) = M_ONE; m2(3,3) = M_ONE
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha)

    ! rotate from => to (around the y axis)
    m3 = M_ZERO; m3(1,1) = M_ONE; m3(2,2) = M_ONE; m3(3,3) = M_ONE
    alpha = acos(sum(from*to))
    call rotate_y(m3, -alpha)
    
    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do i = 1, geo%natoms
      f2 = geo%atom(i)%x
      geo%atom(i)%x = matmul(m1, f2)
    end do

    do i = 1, geo%ncatoms
      f2 = geo%catom(i)%x
      geo%catom(i)%x = matmul(m1, f2)
    end do

  end subroutine rotate

  subroutine rotate_x(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in) :: angle

    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(1, 1) = M_ONE
    aux(2, 2) = ca
    aux(3, 3) = ca
    aux(2, 3) = sa
    aux(3, 2) = -sa

    m = matmul(aux, m)
  end subroutine rotate_x

  subroutine rotate_y(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in) :: angle

    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(2, 2) = M_ONE
    aux(1, 1) = ca
    aux(3, 3) = ca
    aux(1, 3) = sa
    aux(3, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_y

  subroutine rotate_z(m, angle)
    FLOAT, intent(inout) :: m(3,3)
    FLOAT, intent(in) :: angle

    FLOAT :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = M_ZERO
    aux(3, 3) = M_ONE
    aux(1, 1) = ca
    aux(2, 2) = ca
    aux(1, 2) = sa
    aux(2, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_z

end subroutine geometry_adjust

end module geometry
