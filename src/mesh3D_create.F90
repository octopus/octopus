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

  m%box_shape = fdf_integer('Box_shape', SPHERE)
  if (m%box_shape>3 .or. m%box_shape<1) then
    write(err, *) m%box_shape
    message(1) = "Input: '"//trim(err)//"' is not a valid Box_Shape"
    message(2) = '(1 <= Box_Shape <= 3)'
    call write_fatal(2)
  end if
  if(m%box_shape == MINIMUM .and. (.not.present(natoms) .or. .not.present(atom))) then
    message(1) = "Internal Error"
    call write_fatal(1)
  endif

  m%h = fdf_double('spacing', 0.6_r8/units_inp%length%factor)*units_inp%length%factor
  if (m%h<0.01_r8 .or. m%h>2.0_r8) then
    write(err, *) m%h
    message(1) = "Input: '"//trim(err)//"' is not a valid spacing"
    message(2) = '(0.01 <= Spacing [b] <= 2)'
    call write_fatal(2)
  end if
  m%rsize = fdf_double('radius', 20.0_r8/units_inp%length%factor)*units_inp%length%factor
  if (m%rsize<1.0_r8 .or. m%rsize>500.0_r8) then
    write(err, *) m%rsize
    message(1) = "Input: '"//trim(err)//"' is not a valid radius"
    message(2) = '(1 <= radius [b] <= 500)'
    call write_fatal(2)
  end if
  ! this is sometimes required, so we always read....
  m%zsize = fdf_double('zlength', 1.0_r8/units_inp%length%factor)*units_inp%length%factor
  if(m%zsize<1.0_r8 .or. m%zsize>500.0_r8) then
    write(err, *) m%zsize
    message(1) = "Input: '"//trim(err)//"' is not a valid zlength"
    message(2) = '(1 <= zlength [b] <= 500)'
    call write_fatal(2)
  end if

  ! set nr and nx
  m%nr = max(int(m%rsize/m%h), int(m%zsize/m%h))
  m%nx = m%nr + m%d%norder

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(-m%Nx:m%Nx,-m%Nx:m%Nx,-m%Nx:m%Nx))
  allocate(Lxyz(-m%Nx:m%Nx,-m%Nx:m%Nx,-m%Nx:m%Nx))
  m%Lxyz_inv = 0
  Lxyz = 0

  ! We label 1 the points inside the inner mesh
  ! and 2 for the points in the outer mesh
  do ix = -m%nr, m%nr
    do iy = -m%nr, m%nr
      do iz = -m%nr, m%nr
        select case(m%box_shape)
        case(SPHERE)
          b = in_sphere(m, ix, iy, iz)
        case(CILINDER)
          b = in_cilinder(m, ix, iy, iz)
        case(MINIMUM)
          b = in_atom(m, ix, iy, iz, natoms, atom)
        end select
        
        if(b) then
          Lxyz(ix-m%d%norder:ix+m%d%norder, iy, iz) = 2
          Lxyz(ix, iy-m%d%norder:iy+m%d%norder, iz) = 2
          Lxyz(ix, iy, iz-m%d%norder:iz+m%d%norder) = 2
        endif
      enddo
    enddo
  enddo

  do ix = -m%nr, m%nr
    do iy = -m%nr, m%nr
      do iz = -m%nr, m%nr
        select case(m%box_shape)
        case(SPHERE)
          b = in_sphere(m, ix, iy, iz)
        case(CILINDER)
          b = in_cilinder(m, ix, iy, iz)
        case(MINIMUM)
          b = in_atom(m, ix, iy, iz, natoms, atom)
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
  do ix = -m%Nx, m%Nx
    do iy = -m%Nx, m%Nx
      do iz = -m%Nx, m%Nx
        if(Lxyz(ix,iy,iz).eq.1) il = il + 1
        if(Lxyz(ix,iy,iz).eq.2) ik = ik + 1
      enddo
    enddo
  enddo
  m%np = il
  m%nk = ik

  allocate(m%Lx(il), m%Ly(il), m%Lz(il))
  allocate(m%Kx(ik), m%Ky(ik), m%Kz(ik))

  il=0
  ik=0
  do ix = -m%Nx, m%Nx
    do iy = -m%Nx, m%Nx
      do iz = -m%Nx, m%Nx
        if(Lxyz(ix,iy,iz).eq.1) then
          il = il + 1
          m%Lx(il) = ix
          m%Ly(il) = iy
          m%Lz(il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        endif
        if(Lxyz(ix,iy,iz).eq.2) then
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
  ! They check if a point is within the sphere, cilinder, 
  ! or atomic spheres
  logical function in_sphere(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    in_sphere = (sqrt(real(ix**2+iy**2+iz**2,r8))*m%H <= m%rsize)
    return
  end function in_sphere
  
  logical function in_cilinder(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    real(r8) r, z
    r = sqrt(real(ix**2+iy**2, r8))*m%h
    z = iz*m%H
    
    in_cilinder = (r<=m%rsize .and. abs(z)<=m%zsize)
    return
  end function in_cilinder
  
  logical function in_atom(m, ix, iy, iz, natoms, atom)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    integer, intent(in), optional :: natoms
    type(atom_type),  pointer, optional :: atom(:)
    
    integer :: i
    real(r8) :: r
    
    in_atom = .false.
    do i = 1, natoms
      r = sqrt((ix*m%H - atom(i)%x(1))**2 + (iy*m%H - atom(i)%x(2))**2&
           + (iz*m%H - atom(i)%x(3))**2)
      if(r <= m%rsize) then
        in_atom = .true.
        exit
      end if
    end do

    return
  end function in_atom
end subroutine mesh3D_create
