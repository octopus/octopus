module atom
use global
use units
use fdf
use specie

implicit none

type atom_type
  type(specie_type), pointer :: spec ! pointer to specie
#ifndef ONE_D
  real(r8) :: x(3), v(3), f(3) ! position/velocity/force of atom in real space
#else
  real(r8) :: x(1), v(1), f(1)
#endif
  logical :: move              ! should I move this atom in the optimization mode

  ! the mesh around a given atom...
  integer :: Mps
  integer, pointer :: Jxyz(:)
  real(r8), pointer :: pnts_ps, &  ! # points in ps sphere
       uV(:,:), uVu(:),         &  ! the Kleinman Bylander projectors
       duV(:,:,:)                  ! the gradient of the projectors
end type atom_type

type atom_classical_type
  real(r8) :: x(3), v(3), f(3)
  real(r8) :: charge

  character(len=4) :: label
end type atom_classical_type

contains

subroutine atom_init(natoms, a, ncatoms, ca, ns, s)
  integer, intent(out) :: natoms, ncatoms
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: iunit, i
  logical :: l, xyz
  character(len=20) :: label

  sub_name = 'atom_init'; call push_sub()

  if(fdf_defined("PDBCoordinates")) then
    call io_assign(iunit)
    label = fdf_string('PDBCoordinates', 'coords.pdb')
    open(iunit, status='unknown', file=trim(label))

    call loadPDB(iunit)
    
    close(iunit)
  else
    ! we now load the positions, either from the input, or from a file
    if(fdf_defined("XYZCoordinates")) then ! read a xyz file
      call io_assign(iunit)
      label = fdf_string('XYZCoordinates', 'coords.xyz')
      open(iunit, status='unknown', file=trim(label))
      
      read(iunit, *) natoms
      read(iunit) ! skip comment line
      allocate(a(natoms))
      
      call loadXYZ(iunit, .true.)
      close(iunit)
    else
      natoms = fdf_integer('NumberAtoms', 0)
      if(.not. fdf_block("Coordinates", iunit)) then
        message(1) = "Input: Coordinates block not specified"
        message(2) = '%block Coordinates'
        message(3) = '  specie  x  y  z  move'
        message(4) = '%endblock Coordinates'
        call write_fatal(4)
      end if
      
      if (natoms < 1) then
        write(message(1), '(a,i4,a)') "Input: '", natoms, "' is not a valid NumberAtoms"
        message(2) = '(1 <= NumberAtoms)'
        call write_fatal(2)
      end if
      
      call loadXYZ(iunit, .false.)
    end if
  end if

  ! we now load the velocities, either from the input, or from a file
  l = .true.
  xyz = .false.
  if(fdf_defined("XYZVelocities")) then ! read a xyz file
    call io_assign(iunit)
    label = fdf_string('XYZVelocities', 'vel.xyz')
    open(iunit, status='unknown', file=label)

    read(iunit) ! skip number atoms
    read(iunit) ! skip comment line
    xyz = .true.
  else if(.not. fdf_block("Velocities", iunit)) then
    a(:)%v(1) = 0._r8
    a(:)%v(2) = 0._r8
    a(:)%v(3) = 0._r8

    l = .false.
  end if
  
  if(l) then
    do i = 1, natoms
      read(iunit,*) label, a(i)%v(:)
      a(i)%v(:) = units_inp%velocity%factor * a(i)%v(:) !units conversion
    end do

    if(xyz) call io_close(iunit)
  end if

  call pop_sub()

contains
  subroutine loadXYZ(iunit, xyz)
    integer, intent(in) :: iunit
    logical, intent(in) :: xyz

    character(len=20) :: label
    integer :: i

    do i = 1, natoms
      if(xyz) then
        read(iunit,*) label, a(i)%x(:)
        a(i)%move = .true.
      else
        read(iunit,*) label, a(i)%x(:), a(i)%move
      end if
      a(i)%x(:) = units_inp%length%factor * a(i)%x(:) !units conversion
      a(i)%spec => s(get_specie(label))
    end do

  end subroutine loadXYZ

  subroutine loadPDB(iunit)
    integer, intent(in) :: iunit

    character(len=80) :: record
    character(len=6) :: record_name
    character(len=4) :: atm
    character(len=3) :: res
    integer :: na, nca

    ! First count number of atoms
    rewind(iunit)
    natoms = 0
    ncatoms = 0
    do
      read(iunit, '(a80)', err=990) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(17x,a3)') res
        if(trim(res) == 'QM') then
          natoms = natoms + 1
        else
          ncatoms = ncatoms + 1
        end if
      end if      
    end do
990 continue

    allocate(a(natoms), ca(ncatoms))

    ! read in the data
    rewind(iunit)
    na = 1; nca = 1
    do
      read(iunit, '(a80)', err=991) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') atm, res
        call str_trim(atm)
        if(trim(res) == 'QM') then
          read(record, '(30x,3f8.3)') a(na)%x
          a(na)%spec => s(get_specie(atm(1:1)))
          na = na + 1
        else
          ca(nca)%label = atm
          read(record, '(30x,3f8.3,6x,f6.2)') ca(nca)%x, ca(nca)%charge
          nca = nca + 1
        end if
      end if      
    end do
991 continue

  end subroutine loadPDB

  integer function get_specie(label)
    character(len=*) :: label

    integer :: j
    logical :: l

    l = .false.
    do j = 1, ns
      if( trim(label) == trim(s(j)%label) ) then
        l = .true.
        get_specie = j
        return
      end if
    end do
    
    if(.not. l) then
      message(1) = "Specie '"//trim(label)//"' not found"
      call write_fatal(1)
    end if
    
  end function get_specie

end subroutine atom_init

subroutine atom_end(na, a, nca, ca)
  integer, intent(in) :: na, nca
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)

  call atom_dealloc(na, a)

  if(associated(a)) then ! sanity check
    deallocate(a); nullify(a)
  end if

  if(nca > 0 .and. associated(ca)) then
    deallocate(ca); nullify(ca)
  end if
end subroutine atom_end

subroutine atom_dealloc(na, a)
  integer, intent(in) :: na
  type(atom_type), pointer :: a(:)

  integer :: ia

  do ia = 1, na
    if(associated(a(ia)%Jxyz)) then
      deallocate(a(ia)%Jxyz, a(ia)%uV, a(ia)%uVu, a(ia)%duV)
      nullify(a(ia)%Jxyz, a(ia)%uV, a(ia)%uVu, a(ia)%duV)
    end if
  end do

end subroutine atom_dealloc

end module atom
