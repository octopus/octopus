!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creates the mesh
! assumes that the box_shape, h, rsize/zsize of m are filled
! if norder is present, then set nk, Kx, Ky and Kz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh1D_create(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)
  
  ! some local stuff
  character(len=40) :: err
  integer :: next = 0
  integer :: il, ik, ix, iy, iz
  integer, pointer :: Lxyz(:)
  logical :: b

  sub_name = 'mesh1D_create'; call push_sub()

  call oct_parse_int(C_string('Box_shape'), SPHERE, m%box_shape)
  if (m%box_shape>1 .or. m%box_shape<1) then
    write(err, *) m%box_shape
    message(1) = "Input: '"//trim(err)//"' is not a valid Box_Shape"
    message(2) = '(1 <= Box_Shape <= 1)'
    call write_fatal(2)
  end if

  call oct_parse_double(C_string('spacing'), 0.6_r8/units_inp%length%factor, m%h)
  m%h = m%h * units_inp%length%factor
  if (m%h<0.01_r8 .or. m%h>2.0_r8) then
    write(err, *) m%h
    message(1) = "Input: '"//trim(err)//"' is not a valid spacing"
    message(2) = '(0.01 <= Spacing [b] <= 2)'
    call write_fatal(2)
  end if
  m%vol_pp = m%h

  call oct_parse_double(C_string('radius'), 20.0_r8/units_inp%length%factor, m%rsize)
  m%rsize = m%rsize * units_inp%length%factor
  if (m%rsize<1.0_r8 .or. m%rsize>500.0_r8) then
    write(err, *) m%rsize
    message(1) = "Input: '"//trim(err)//"' is not a valid radius"
    message(2) = '(1 <= radius [b] <= 500)'
    call write_fatal(2)
  end if

  ! set nr and nx
  m%nr = int(m%rsize/m%h)
  m%nx = m%nr + m%d%norder

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(-m%Nx:m%Nx))
  allocate(Lxyz(-m%Nx:m%Nx))
  m%Lxyz_inv = 0
  Lxyz = 0

  ! We label 1 the points inside the inner mesh
  ! and 2 for the points in the outer mesh
  do ix = -m%nr, m%nr
        b = in_sphere(m, ix)
        if(b) then
          Lxyz(ix-m%d%norder:ix+m%d%norder) = 2
        endif
  enddo

  do ix = -m%nr, m%nr
        b = in_sphere(m, ix)
        if(b) then
          Lxyz(ix) = 1
        endif
  enddo

  ! now, we count the number of internal and external points
  il=0
  ik=0
  do ix = -m%Nx, m%Nx
        if(Lxyz(ix) == 1) il = il + 1
        if(Lxyz(ix) == 2) ik = ik + 1
  enddo
  m%np = il
  m%nk = ik

  allocate(m%Lx(il))
  allocate(m%Kx(ik))

  il=0
  ik=0
  do ix = -m%Nx, m%Nx
        if(Lxyz(ix) == 1) then
          il = il + 1
          m%Lx(il) = ix
          m%Lxyz_inv(ix) = il
        endif
        if(Lxyz(ix) == 2) then
          ik = ik + 1
          m%Kx(ik) = ix
        endif
  enddo

  deallocate(Lxyz); nullify(Lxyz)

  call mesh_write_info(m, stdout)

  call pop_sub()
  return

contains

  ! The three next functions are self-explanatory!
  ! They check if a point is within the sphere, cilinder, 
  ! or atomic spheres
  logical function in_sphere(m, ix)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix
    
    in_sphere = (real(ix,r8)*m%H <= m%rsize)
    return
  end function in_sphere


end subroutine mesh1D_create
