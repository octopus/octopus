module atom
use global
use units
use fdf
use specie

implicit none

type atom_type
  type(specie_type), pointer :: spec ! pointer to specie
  real(r8) :: x(3), v(3), f(3) ! position/velocity/force of atom in real space
  logical :: move              ! should I move this atom in the optimization mode
end type atom_type

contains

function atom_init(a, ns, s)
  integer :: atom_init
  type(atom_type), pointer :: a(:)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: natoms, iunit, i, j
  logical :: l, xyz
  character(len=20) :: label

  sub_name = 'atom_init'; call push_sub()

  ! we now load the positions, either from the input, or from a file
  if(fdf_defined("CoordinatesFile")) then ! read a xyz file
    call io_assign(iunit)
    label = fdf_string('CoordinatesFile', 'coords.xyz')
    open(iunit, status='unknown', file=label)

    read(iunit, *) natoms
    read(iunit) ! skip comment line
    xyz = .true.
  else
    natoms = fdf_integer('NumberAtoms', 0)
    xyz = .false.
    if(.not. fdf_block("AtomCoordinates", iunit)) then
      message(1) = "Input: AtomCoordinates block not specified"
      message(2) = '%block AtomCoordinates'
      message(3) = '  specie  x  y  z  move'
      message(4) = '%endblock AtomCoordinates'
      call write_fatal(4)
    endif
  end if

  if (natoms < 1) then
    write(message(1), '(a,i4,a)') "Input: '", natoms, "' is not a valid NumberAtoms"
    message(2) = '(1 <= NumberAtoms)'
    call write_fatal(2)
  end if
  allocate(a(natoms))

  do i = 1, natoms
    if(xyz) then
      read(iunit,*) label, a(i)%x(:)
      a(i)%move = .true.
    else
      read(iunit,*) label, a(i)%x(:), a(i)%move
    end if
    a(i)%x(:) = units_inp%length%factor * a(i)%x(:) !units conversion
      
    l = .false.
    do j = 1, ns
      if( label == s(j)%label ) then
        l = .true.
        a(i)%spec => s(j)
        exit
      end if
    end do
      
    if(.not. l) then
      message(1) = "Specie '"//trim(label)//"' not found"
      call write_fatal(1)
    end if
  end do

  if(xyz) call io_close(iunit)

  ! we now load the velocities, either from the input, or from a file
  l = .true.
  xyz = .false.
  if(fdf_defined("VelocitiesFile")) then ! read a xyz file
    call io_assign(iunit)
    label = fdf_string('VelocitiesFile', 'vel.xyz')
    open(iunit, status='unknown', file=label)

    read(iunit) ! skip number atoms
    read(iunit) ! skip comment line
    xyz = .true.
  else if(.not. fdf_block("AtomVelocities", iunit)) then
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

  atom_init = natoms
  call pop_sub()

end function atom_init

subroutine atom_end(na, a)
  integer, intent(in) :: na
  type(atom_type), pointer :: a(:)

  if(associated(a)) then ! sanity check
    deallocate(a); nullify(a)
  end if

end subroutine atom_end

end module atom
