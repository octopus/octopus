!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creates the mesh
! assumes that the box_shape, h, rsize/zsize of m are filled
! if norder is present, then set nk, Kx, Ky and Kz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh3D_create(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)
  
  ! some local stuff
  character(len=40) :: err
  integer :: next = 0
  integer :: il, ik, ix, iy, iz
  integer, pointer :: Lxyz(:,:,:)
  logical :: b

  sub_name = 'mesh3D_create'; call push_sub()

  ! Read box shape.
  call oct_parse_int(C_string('BoxShape'), SPHERE, m%box_shape)
  if (m%box_shape>4 .or. m%box_shape<1) then
    write(err, *) m%box_shape
    message(1) = "Input: '"//trim(err)//"' is not a valid BoxShape"
    message(2) = '(1 <= Box_Shape <= 4)'
    call write_fatal(2)
  end if
  if(m%box_shape == MINIMUM .and. (.not.present(natoms) .or. .not.present(atom))) then
    message(1) = "Internal Error"
    call write_fatal(1)
  endif
#ifdef POLYMERS
  if((m%box_shape.ne.CYLINDER).or.(m%box_shape.ne.PARALLELEPIPED)) then
    message(1) = "When running polymer calculations, you have to use"
    message(2) = "BoxShape = cylinder | parallelepiped"
    call write_fatal(2)
  end if
#endif

  ! Read the grid spacing.
  select case(m%box_shape)
  case(SPHERE, CYLINDER, MINIMUM)
    call oct_parse_double(C_string('spacing'), 0.6_r8/units_inp%length%factor, m%h(1))
    m%h(1) = m%h(1)*units_inp%length%factor 
    m%h(2) = m%h(1); m%h(3) = m%h(1)
  case(PARALLELEPIPED)
    call oct_parse_block_double(C_string('spacing'), 0, 0, m%h(1))
    call oct_parse_block_double(C_string('spacing'), 0, 1, m%h(2))
    call oct_parse_block_double(C_string('spacing'), 0, 2, m%h(3))
    m%h = m%h*units_inp%length%factor
  end select
  m%vol_pp = m%h(1)*m%h(2)*m%h(3)

  ! Read the box size.
  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER .or. m%box_shape == MINIMUM) then
    call oct_parse_double(C_string('radius'), 20.0_r8/units_inp%length%factor, m%rsize)
    m%rsize = m%rsize * units_inp%length%factor
  end if
  if(m%box_shape == CYLINDER) then
    call oct_parse_double(C_string('zlength'), 1.0_r8/units_inp%length%factor, m%zsize)
    m%zsize = m%zsize * units_inp%length%factor
  endif
  if(m%box_shape == PARALLELEPIPED) then
    call oct_parse_block_double(C_string('lsize'), 0, 0, m%lsize(1))
    call oct_parse_block_double(C_string('lsize'), 0, 1, m%lsize(2))
    call oct_parse_block_double(C_string('lsize'), 0, 2, m%lsize(3))
    m%lsize = m%lsize*units_inp%length%factor
  endif

  ! set nr and nx
  select case(m%box_shape)
  case(SPHERE)
    m%nr(1:3) = int(m%rsize/m%h(1))
  case(CYLINDER)
    m%nr(1:3) = max(int(m%rsize/m%h(1)), int(m%zsize/m%h(1)))
  case(MINIMUM)
    message(1) = 'MINIMUM cage not correctly implemented. Sorry.'
    call write_fatal(1)
  case(PARALLELEPIPED)
    m%nr(:) = int((m%lsize(:)/2.0_r8)/m%h(:))
  end select
  m%nx(1:3) = m%nr(1:3) + m%d%norder

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(-m%Nx(1):m%Nx(1),-m%Nx(2):m%Nx(2),-m%Nx(3):m%Nx(3)))
  allocate(Lxyz(-m%Nx(1):m%Nx(1),-m%Nx(2):m%Nx(2),-m%Nx(3):m%Nx(3)))
  m%Lxyz_inv = 0
  Lxyz = 0

  ! We label 1 the points inside the inner mesh
  ! and 2 for the points in the outer mesh
  do ix = -m%nr(1), m%nr(1)
    do iy = -m%nr(2), m%nr(2)
      do iz = -m%nr(3), m%nr(3)
        select case(m%box_shape)
        case(SPHERE)
          b = in_sphere(m, ix, iy, iz)
        case(CYLINDER)
          b = in_cylinder(m, ix, iy, iz)
        case(MINIMUM)
          b = in_atom(m, ix, iy, iz, natoms, atom)
        case(PARALLELEPIPED)
          b = .true.
        end select
        
        if(b) then
          Lxyz(ix-m%d%norder:ix+m%d%norder, iy, iz) = 2
          Lxyz(ix, iy-m%d%norder:iy+m%d%norder, iz) = 2
          Lxyz(ix, iy, iz-m%d%norder:iz+m%d%norder) = 2
        endif
      enddo
    enddo
  enddo

#ifdef POLYMERS
  ! in this case the leftmost point is equivalent to the rightmost
  ! so we take it out
  il = m%nr(3) - 1
#else
  il = m%nr(3)
#endif

  do ix = -m%nr(1), m%nr(1)
    do iy = -m%nr(2), m%nr(2)
      do iz = -m%nr(3), il
        select case(m%box_shape)
        case(SPHERE)
          b = in_sphere(m, ix, iy, iz)
        case(CYLINDER)
          b = in_cylinder(m, ix, iy, iz)
        case(MINIMUM)
          b = in_atom(m, ix, iy, iz, natoms, atom)
        case(PARALLELEPIPED)
          b = .true.
        end select
        
        if(b) then
          Lxyz(ix, iy, iz) = 1
        endif
      enddo
    enddo
  enddo

  ! now, we count the number of internal and external points
  il=0
  ik=0
  do ix = -m%Nx(1), m%Nx(1)
    do iy = -m%Nx(2), m%Nx(2)
      do iz = -m%Nx(3), m%Nx(3)
        if(Lxyz(ix, iy, iz) == 1) il = il + 1
        if(Lxyz(ix, iy, iz) == 2) ik = ik + 1
      enddo
    enddo
  enddo
  m%np = il
  m%nk = ik

  allocate(m%Lx(il), m%Ly(il), m%Lz(il))
  allocate(m%Kx(ik), m%Ky(ik), m%Kz(ik))

  il=0
  ik=0
  do ix = -m%Nx(1), m%Nx(1)
    do iy = -m%Nx(2), m%Nx(2)
      do iz = -m%Nx(3), m%Nx(3)
        if(Lxyz(ix,iy,iz) == 1) then
          il = il + 1
          m%Lx(il) = ix
          m%Ly(il) = iy
          m%Lz(il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        endif
        if(Lxyz(ix,iy,iz) == 2) then
          ik = ik + 1
          m%Kx(ik) = ix
          m%Ky(ik) = iy
          m%Kz(ik) = iz
        endif
      enddo
    enddo
  enddo

  deallocate(Lxyz); nullify(Lxyz)

  call mesh_write_info(m, stdout)

  call pop_sub()
  return

contains

  ! The three next functions are self-explanatory!
  ! They check if a point is within the sphere, cylinder, 
  ! or atomic spheres
  logical function in_sphere(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    in_sphere = (sqrt(real(ix**2+iy**2+iz**2,r8))*m%H(1) <= m%rsize)
    return
  end function in_sphere
  
  logical function in_cylinder(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    real(r8) r, z
    r = sqrt(real(ix**2+iy**2, r8))*m%h(1)
    z = iz*m%H(3)
    
    in_cylinder = (r<=m%rsize .and. abs(z)<=m%zsize)
    return
  end function in_cylinder
  
  logical function in_atom(m, ix, iy, iz, natoms, atom)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    integer, intent(in), optional :: natoms
    type(atom_type),  pointer, optional :: atom(:)
    
    integer :: i
    real(r8) :: r
    
    in_atom = .false.
    do i = 1, natoms
      r = sqrt((ix*m%H(1) - atom(i)%x(1))**2 + (iy*m%H(1) - atom(i)%x(2))**2&
           + (iz*m%H(1) - atom(i)%x(3))**2)
      if(r <= m%rsize) then
        in_atom = .true.
        exit
      end if
    end do

    return
  end function in_atom

end subroutine mesh3D_create
