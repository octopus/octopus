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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creates the mesh
! assumes that the box_shape, h, rsize/xsize of m are filled
! if norder is present, then set nk, Kx, Ky and Kz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh_create(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)
  
  ! some local stuff
  integer :: i, il, ik, ix, iy, iz
  integer :: ii(3), nr_tmp(3)
  integer, pointer :: Lxyz_tmp(:,:,:)
  logical :: b
  
  call push_sub('mesh_create')
  
  ! Read box shape.
  call oct_parse_int('BoxShape', SPHERE, m%box_shape)
  if(m%box_shape<1 .or. m%box_shape>4) then
    write(message(1), '(a,i2,a)') "Input: '", m%box_shape, "' is not a valid BoxShape"
    message(2) = '(1 <= Box_Shape <= 4)'
    call write_fatal(2)
  end if
  if(m%box_shape == MINIMUM .and. (.not.present(natoms) .or. .not.present(atom))) then
    message(1) = "Internal Error"
    call write_fatal(1)
  endif
  
  select case(m%box_shape)
  case(SPHERE)
    if(conf%dim>1 .and. conf%periodic_dim>0) then
      message(1) = 'Spherical mesh is not allowed for periodic systems'
      call write_fatal(1)
    end if
  case(CYLINDER)
    if (conf%dim>2 .and. &
    ((conf%dim-conf%periodic_dim == 0) .or. (conf%dim-conf%periodic_dim == 1))) then
        message(1) = 'Cylindrical mesh is not allowed for systems'
        message(2) = 'that are periodic in more than one dimension'
        call write_fatal(2)
    end if
  case(MINIMUM)
    select case(conf%dim)
    case(1)
      message(1) = 'Minimum mesh cannot be used, switching to spherical'
      call write_warning(1)
    case default
      if (conf%periodic_dim>0) then
        message(1)= 'Minimum mesh is not allowed for periodic systems'
       call write_fatal(1)
      end if
    end select
  end select

  ! ignore box_shape in 1D
  if(conf%dim==1) m%box_shape=SPHERE

  ! Read the box size.
  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER .or. m%box_shape == MINIMUM) then
    call oct_parse_double('radius', 20.0_r8/units_inp%length%factor, m%rsize)
    m%rsize = m%rsize * units_inp%length%factor
    m%lsize(1) = m%rsize
  end if
  if(m%box_shape == CYLINDER) then
    call oct_parse_double('xlength', 1.0_r8/units_inp%length%factor, m%xsize)
    m%xsize = m%xsize * units_inp%length%factor
    m%lsize(1) = m%xsize
  end if
  if(m%box_shape == PARALLELEPIPED) then
    if(oct_parse_block_n('lsize')<1) then
      message(1) = 'Block "lsize" not found in input file.'
      call write_fatal(1)
    endif
    do i = 1, conf%dim
      call oct_parse_block_double('lsize', 0, i-1, m%lsize(i))
    end do
    m%lsize = m%lsize*units_inp%length%factor
  end if
  
  ! build primitive vectors (only simple cubic, tetra, or orthororhombic )
  m%rlat = M_ZERO
  m%klat = M_ZERO
  do i = 1, conf%periodic_dim
    m%rlat(i,i) = 2*m%lsize(i)
    m%klat(i,i) = M_PI/m%lsize(i)
  end do
 
  !build shifts to nearest neighbour primitive cells
  ii=(/0,-1,1/)
  m%shift=M_ZERO
  i=1
  do iz = 1, 3
    do iy = 1, 3
      do ix = 1, 3
        m%shift(i,1)=ii(ix)*m%rlat(1,1)
        m%shift(i,2)=ii(iy)*m%rlat(2,2)
        m%shift(i,3)=ii(iz)*m%rlat(3,3)
        i = i + 1
      end do
    end do
  end do
  
  ! Read the grid spacing.
  m%h = M_ZERO
  select case(m%box_shape)
  case(SPHERE,CYLINDER,MINIMUM)
    call oct_parse_double('spacing', 0.6_r8/units_inp%length%factor, m%h(1))
    m%h(1:conf%dim) = m%h(1)
    m%iso = .true.
  case(PARALLELEPIPED)
    if(oct_parse_block_n('spacing')<1) then
      m%h(1:3) = 0.6_r8/units_inp%length%factor
    else
      do i = 1, conf%dim
         call oct_parse_block_double('spacing', 0, i-1, m%h(i))
      end do
    endif
    m%iso = .true.
    do i = 2, conf%dim
      if(m%h(i).ne.m%h(1)) m%iso = .false.
    end do
  end select
  m%h(1:conf%dim) = m%h(1:conf%dim)*units_inp%length%factor
  
  ! set nr and nx, and adjust the mesh so that:
  ! 1) the new grid exactly fills the box;
  ! 2) the new mesh is not larger than the user defined mesh.
  m%nr = 0; m%nx = 0
  select case(m%box_shape)
  case(SPHERE)
    m%nr(1:conf%dim) = int(m%rsize/m%h(1))
  case(CYLINDER)
    m%nr(1:conf%dim) = max(int(m%rsize/m%h(1)), int(m%xsize/m%h(1)))
  case(MINIMUM)
    message(1) = 'MINIMUM cage not correctly implemented. Sorry.'
    call write_fatal(1)
  case(PARALLELEPIPED)
    m%nr(1:conf%dim) = int(m%lsize(1:conf%dim)/m%h(1:conf%dim))
    do i = 1, conf%dim
      if (mod(m%lsize(i),m%h(i)) /= M_ZERO) then
        m%nr(i) = m%nr(i) + 1
        m%h(i) = m%lsize(i)/real(m%nr(i))
      end if
    end do 
  end select
  m%nx(1:conf%dim) = m%nr(1:conf%dim) + m%d%norder

  m%vol_pp = M_ONE
  do i = 1, conf%dim
    m%vol_pp = m%vol_pp*m%h(i)
  end do

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(-m%nx(1):m%nx(1),-m%nx(2):m%nx(2),-m%nx(3):m%nx(3)))
  allocate(Lxyz_tmp(-m%nx(1):m%nx(1),-m%nx(2):m%nx(2),-m%nx(3):m%nx(3)))

  ! initialize to zero
  m%Lxyz_inv(-m%nx(1):m%nx(1),-m%nx(2):m%nx(2),-m%nx(3):m%nx(3)) = 0
  Lxyz_tmp(-m%nx(1):m%nx(1),-m%nx(2):m%nx(2),-m%nx(3):m%nx(3)) = 0

  ! We label 1 the points inside the inner mesh and 2 for the points in the outer mesh
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
          Lxyz_tmp(ix-m%d%norder:ix+m%d%norder, iy, iz) = 2
          if(conf%dim > 1) Lxyz_tmp(ix, iy-m%d%norder:iy+m%d%norder, iz) = 2
          if(conf%dim > 2) Lxyz_tmp(ix, iy, iz-m%d%norder:iz+m%d%norder) = 2
        endif
      enddo
    enddo
  enddo

  ! the last point is equivalent to the first one in periodic directions
  nr_tmp=m%nr
  do i = 1, conf%periodic_dim
    nr_tmp(i)=nr_tmp(i)-1
  end do

  do ix = -m%nr(1), nr_tmp(1)
    do iy = -m%nr(2), nr_tmp(2)
      do iz = -m%nr(3), nr_tmp(3)
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
          Lxyz_tmp(ix, iy, iz) = 1
        end if
      end do
    end do
  end do
  
  ! now, we count the number of internal and external points
  il=0
  ik=0
  do ix = -m%nx(1), m%nx(1)
    do iy = -m%nx(2), m%nx(2)
      do iz = -m%nx(3), m%nx(3)
        if(Lxyz_tmp(ix, iy, iz) == 1) il = il + 1
        if(Lxyz_tmp(ix, iy, iz) == 2) ik = ik + 1
      end do
    end do
  end do
  m%np = il
  m%nk = ik

  allocate(m%Lxyz(3, il), m%Kxyz(3, ik))
  il = 0; ik = 0
  do ix = -m%nx(1), m%nx(1)
    do iy = -m%nx(2), m%nx(2)
      do iz = -m%nx(3), m%nx(3)
        if(Lxyz_tmp(ix, iy, iz) == 1) then
          il = il + 1
          m%Lxyz(1, il) = ix
          m%Lxyz(2, il) = iy
          m%Lxyz(3, il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        end if
        if(Lxyz_tmp(ix, iy, iz) == 2) then
          ik = ik + 1
          m%Kxyz(1, ik) = ix
          m%Kxyz(2, ik) = iy
          m%Kxyz(3, ik) = iz
        end if
      enddo
    enddo
  enddo

  allocate(m%ind(2*conf%dim, m%d%norder, m%np))
  m%ind = 0
  do il = 1, m%np
    ix = m%Lxyz(1, il); iy = m%Lxyz(2, il); iz = m%Lxyz(3, il)
    do ik = 1, m%d%norder

      m%ind(1, ik, il) = m%Lxyz_inv(ix-ik, iy, iz)
      m%ind(2, ik, il) = m%Lxyz_inv(ix+ik, iy, iz)
      if(conf%dim > 1) then
        m%ind(3, ik, il) = m%Lxyz_inv(ix, iy-ik, iz)
        m%ind(4, ik, il) = m%Lxyz_inv(ix, iy+ik, iz)
        if(conf%dim > 2) then
          if(conf%dim > 2) m%ind(5, ik, il) = m%Lxyz_inv(ix, iy, iz-ik)
          if(conf%dim > 2) m%ind(6, ik, il) = m%Lxyz_inv(ix, iy, iz+ik)
        end if
      end if

    end do
  end do


  deallocate(Lxyz_tmp); nullify(Lxyz_tmp)

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
    
    in_sphere = (sqrt(real(ix**2+iy**2+iz**2,r8))*m%h(1) <= m%rsize)
    return
  end function in_sphere
  
  logical function in_cylinder(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    real(r8) r, x
    r = sqrt(real(iy**2+iz**2, r8))*m%h(1)
    x = ix*m%H(1)
    
    in_cylinder = (r<=m%rsize .and. abs(x)<=m%xsize)
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
      r = sqrt((ix*m%H(1) - atom(i)%x(1))**2 + (iy*m%H(1) - atom(i)%x(2))**2 &
           + (iz*m%H(1) - atom(i)%x(3))**2)
      if(r <= m%rsize) then
        in_atom = .true.
        exit
      end if
    end do

    return
  end function in_atom

end subroutine mesh_create
