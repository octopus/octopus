program centergeom
  use global
  use units
  use lib_oct_parser
  use atom

  integer :: natoms, ncatoms, iunit, i
  character(len=80) :: label, str
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)

  call global_init()
  call units_init()

  if(loct_parse_isdef("PDBCoordinates").ne.0) then
    call loct_parse_string('PDBCoordinates', 'coords.pdb', label)
    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    call loadPDB(iunit, natoms, a, ncatoms, ca)
    call io_close(iunit)
  elseif(loct_parse_isdef("XYZCoordinates").ne.0) then
    call loct_parse_string('XYZCoordinates', 'coords.xyz', label)
    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    read(iunit, *) natoms
    read(iunit, *) ! skip comment line
    allocate(a(natoms))
    nullify(ca); ncatoms = 0;
    do i = 1, natoms
       read(iunit,*) a(i)%label, a(i)%x(:)
       a(i)%move = .true.
    end do
    call io_close(iunit)
  else
    str = "Coordinates"
    natoms = loct_parse_block_n(str)
    if(natoms <= 0) then
      message(1) = "Input: Coordinates block not specified"
      message(2) = '% Coordinates'
      message(3) = '  specie  x  y  z  move'
      message(4) = '%'
      call write_fatal(4)
    end if
      
    allocate(a(natoms))
    nullify(ca); ncatoms = 0
    do i = 1, natoms
       call loct_parse_block_string(str, i-1, 0, a(i)%label)
       call loct_parse_block_float  (str, i-1, 1, a(i)%x(1))
       call loct_parse_block_float  (str, i-1, 2, a(i)%x(2))
       call loct_parse_block_float  (str, i-1, 3, a(i)%x(3))
       call loct_parse_block_logical(str, i-1, 4, a(i)%move)
    end do
  endif

  ! units conversion
  do i = 1, natoms
    a(i)%x = a(i)%x * units_inp%length%factor
  end do
  do i = 1, ncatoms
    ca(i)%x = ca(i)%x * units_inp%length%factor
  end do
  
  call atom_adjust(natoms, a , ncatoms, ca)
  call atom_write_xyz(".", "adjusted", natoms, a, ncatoms, ca)

end program centergeom
