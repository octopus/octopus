program centergeom
  use global
  use units
  use lib_oct_parser
  use atom
  use geometry

  implicit none

  integer :: iunit, i
  character(len=80) :: label, str
  type(geometry_type) :: geo

  call global_init()
  call units_init()

  if(loct_parse_isdef("PDBCoordinates").ne.0) then
    call loct_parse_string('PDBCoordinates', 'coords.pdb', label)
    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    call loadPDB(iunit, geo)
    call io_close(iunit)
  elseif(loct_parse_isdef("XYZCoordinates").ne.0) then
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
       call loct_parse_block_string(str, i-1, 0, geo%atom(i)%label)
       call loct_parse_block_float  (str, i-1, 1, geo%atom(i)%x(1))
       call loct_parse_block_float  (str, i-1, 2, geo%atom(i)%x(2))
       call loct_parse_block_float  (str, i-1, 3, geo%atom(i)%x(3))
       call loct_parse_block_logical(str, i-1, 4, geo%atom(i)%move)
    end do
  endif

  ! units conversion
  do i = 1, geo%natoms
    geo%atom(i)%x = geo%atom(i)%x * units_inp%length%factor
  end do
  do i = 1, geo%ncatoms
    geo%catom(i)%x = geo%catom(i)%x * units_inp%length%factor
  end do
  
  call geometry_adjust(geo)
  call atom_write_xyz(".", "adjusted", geo)

end program centergeom
