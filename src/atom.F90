module atom
use global
use liboct
use io
use units
use specie

implicit none

type atom_type
  type(specie_type), pointer :: spec ! pointer to specie

  real(r8) :: x(3), v(3), f(3) ! position/velocity/force of atom in real space

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
  character(len=40) :: str, label

  sub_name = 'atom_init'; call push_sub()

  if(conf%dim == 3.and.oct_parse_isdef(C_string("PDBCoordinates")).ne.0) then
    call oct_parse_str('PDBCoordinates', 'coords.pdb', label)

    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    call loadPDB(iunit)
    close(iunit)

  else
    ! we now load the positions, either from the input, or from a file
    if(conf%dim == 3.and.oct_parse_isdef(C_string("XYZCoordinates")).ne.0) then ! read a xyz file
      call oct_parse_str('XYZCoordinates', 'coords.xyz', label)

      call io_assign(iunit)
      open(iunit, status='unknown', file=trim(label))
      read(iunit, *) natoms
      read(iunit) ! skip comment line
      allocate(a(natoms))
      
      do i = 1, natoms
        read(iunit,*) label, a(i)%x(:)
        a(i)%move = .true.
        a(i)%x(:) = units_inp%length%factor * a(i)%x(:) !units conversion
        a(i)%spec => s(get_specie(label))
      end do

      close(iunit)
    else
      str = C_string("Coordinates")
      natoms = oct_parse_block_n(str)
      if(natoms <= 0) then
        message(1) = "Input: Coordinates block not specified"
        message(2) = '% Coordinates'
        message(3) = '  specie  x  y  z  move'
        message(4) = '%'
        call write_fatal(4)
      end if
      
      allocate(a(natoms))
      do i = 1, natoms
        call oct_parse_block_str    (str, i-1, 0, label)
        a(i)%spec => s(get_specie(label))
        call oct_parse_block_double (str, i-1, 1, a(i)%x(1))
        call oct_parse_block_double (str, i-1, 2, a(i)%x(2))
        call oct_parse_block_double (str, i-1, 3, a(i)%x(3))
        a(i)%x = a(i)%x * units_inp%length%factor
        call oct_parse_block_logical(str, i-1, 4, a(i)%move)
      end do
    end if
  end if

  ! we now load the velocities, either from the input, or from a file
  if(oct_parse_isdef(C_string("XYZVelocities")).ne.0 .and. conf%dim==3) then ! read a xyz file
    call io_assign(iunit)
    call oct_parse_str('XYZVelocities', 'coords.xyz', label)
    open(iunit, status='unknown', file=trim(label))
      
    read(iunit)
    read(iunit) ! skip comment line
      
    do i = 1, natoms
      read(iunit,*) label, a(i)%v(:)
      a(i)%v(:) = units_inp%velocity%factor * a(i)%v(:) !units conversion
    end do

    call io_close(iunit)
  else 
    str = C_string("Velocities")
    if(oct_parse_isdef(str).ne.0) then
      do i = 1, natoms
        call oct_parse_block_double (str, i-1, 1, a(i)%v(1))
        call oct_parse_block_double (str, i-1, 2, a(i)%v(2))
        call oct_parse_block_double (str, i-1, 3, a(i)%v(3))
        a(i)%v = a(i)%v * units_inp%velocity%factor
      end do
    else
      a(:)%v(1) = 0._r8
      a(:)%v(2) = 0._r8
      a(:)%v(3) = 0._r8
    end if
  end if
  
  call pop_sub()

contains

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
