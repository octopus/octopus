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

#define DELTA_R 1e-12_r8

subroutine mesh_create(m, natoms, atom, enlarge_)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)
  integer, intent(in), optional :: enlarge_

  ! some local stuff
  integer :: i, ix, iy, iz, ii(3), enlarge
  
  call push_sub('mesh_create')
  
  enlarge = 0
  if(present(enlarge_)) enlarge = enlarge_

  ! Read box shape.
  call oct_parse_int('BoxShape', SPHERE, m%box_shape)
  if(m%box_shape<1 .or. m%box_shape>4) then
    write(message(1), '(a,i2,a)') "Input: '", m%box_shape, "' is not a valid BoxShape"
    message(2) = '(1 <= Box_Shape <= 4)'
    call write_fatal(2)
  end if
  
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
  end select

  ! ignore box_shape in 1D
  if(conf%dim==1) m%box_shape=SPHERE

  ! Read the box size.
  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER) then
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
 
  ! build shifts to nearest neighbour primitive cells
  ii = (/0,-1,1/)
  m%shift=M_ZERO
  i = 1
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
  case(SPHERE,CYLINDER)
    call oct_parse_double('spacing', 0.6_r8/units_inp%length%factor, m%h(1))
    m%h(1:conf%dim) = m%h(1)
  case(PARALLELEPIPED)
    if(oct_parse_block_n('spacing')<1) then
      m%h(1:3) = 0.6_r8/units_inp%length%factor
    else
      do i = 1, conf%dim
         call oct_parse_block_double('spacing', 0, i-1, m%h(i))
      end do
    endif
  end select
  m%h(1:conf%dim) = m%h(1:conf%dim)*units_inp%length%factor

  ! set nr and adjust the mesh so that:
  ! 1) the new grid exactly fills the box;
  ! 2) the new mesh is not larger than the user defined mesh.
  m%nr = 0
  select case(m%box_shape)
  case(SPHERE)
    m%nr(2, 1:conf%dim) = int((m%rsize+DELTA_R)/m%h(1))
  case(CYLINDER)
    m%nr(2, 1:conf%dim) = max(int((m%rsize+DELTA_R)/m%h(1)), int((m%xsize+DELTA_R)/m%h(1)))
  case(PARALLELEPIPED)
    m%nr(2, 1:conf%dim) = int((m%lsize(1:conf%dim)+DELTA_R)/m%h(1:conf%dim))
    do i = 1, conf%dim
      if (mod(m%lsize(i),m%h(i)) /= M_ZERO) then
        m%nr(2, i) = m%nr(2, i) + 1
        m%h(i) = m%lsize(i)/real(m%nr(2, i))
      end if
    end do 
  end select
  m%nr(1,:) = -m%nr(2,:)

  ! the last point is equivalent to the first one in periodic directions
  do i = 1, conf%periodic_dim
    m%nr(2, i) = m%nr(2, i) - 1
  end do

  m%vol_pp = M_ONE
  do i = 1, conf%dim
    m%vol_pp = m%vol_pp*m%h(i)
  end do

  call mesh_create_xyz(m, enlarge)
  call mesh_write_info(m, stdout)

  call pop_sub()
  return
end subroutine mesh_create

subroutine mesh_create_xyz(m, enlarge)
  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: enlarge

  integer :: i, il, ik, ix, iy, iz
  integer, pointer :: Lxyz_tmp(:,:,:)
  logical :: b

! enlarge mesh in the non-periodic dimensions
  do i = conf%periodic_dim+1, conf%dim
    m%nr(1,i) = m%nr(1,i) - enlarge
    m%nr(2,i) = m%nr(2,i) + enlarge
  end do

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(  Lxyz_tmp(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  m%Lxyz_inv(:,:,:) = 0
  Lxyz_tmp(:,:,:)   = 0
    

  ! We label 2 the points inside the mesh + enlargement
  if(enlarge > 0) then
    do ix = m%nr(1,1), m%nr(2,1)
      do iy = m%nr(1,2), m%nr(2,2)
        do iz = m%nr(1,3), m%nr(2,3)
          if(in_mesh(m, ix, iy, iz)) then
            do i = -enlarge, 0 ! first include the points on the left
              if(ix+i>=m%nr(1,1))                Lxyz_tmp(ix+i, iy, iz) = 2
              if(iy+i>=m%nr(1,2).and.conf%dim>1) Lxyz_tmp(ix, iy+i, iz) = 2
              if(iz+i>=m%nr(1,3).and.conf%dim>2) Lxyz_tmp(ix, iy, iz+i) = 2
            end do
            do i = 1, enlarge ! and now on the right
              if(ix+i<=m%nr(2,1))                Lxyz_tmp(ix+i, iy, iz) = 2
              if(iy+i<=m%nr(2,2).and.conf%dim>1) Lxyz_tmp(ix, iy+i, iz) = 2
              if(iz+i<=m%nr(2,3).and.conf%dim>2) Lxyz_tmp(ix, iy, iz+i) = 2
            end do
          endif
        end do
      end do
    end do
  end if

  ! we label 1 the points inside the mesh
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(in_mesh(m, ix, iy, iz)) Lxyz_tmp(ix, iy, iz) = 1
      end do
    end do
  end do

  ! now, we count the number of internal points
  il = 0
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) > 0) il = il + 1
      end do
    end do
  end do
  m%np = il

  allocate(m%Lxyz(3, il))
  il = 0

  ! first we fill the poins in the inner mesh
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 1) then
          il = il + 1
          m%Lxyz(1, il) = ix
          m%Lxyz(2, il) = iy
          m%Lxyz(3, il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        end if
      enddo
    enddo
  enddo
  m%np_in = il

  ! and now the points from the enlargement
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 2) then
          il = il + 1
          m%Lxyz(1, il) = ix
          m%Lxyz(2, il) = iy
          m%Lxyz(3, il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        end if
      end do
    end do
  end do

  ! now we calculate a lookup table for the derivatives
  call lookup_table(m, m%d)

  deallocate(Lxyz_tmp); nullify(Lxyz_tmp)
end subroutine mesh_create_xyz

logical function in_mesh(m, ix, iy, iz)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: ix, iy, iz
    

  select case(m%box_shape)
  case(SPHERE)
    in_mesh = (sqrt(real(ix**2 + iy**2 + iz**2, r8))*m%h(1) <= m%rsize+DELTA_R)
  case(CYLINDER)
    in_mesh = in_cylinder()
  case(PARALLELEPIPED)
    in_mesh = (ix >= m%nr(1,1).and.ix <= m%nr(2,1)).and. &
         (iy >= m%nr(1,2).and.iy <= m%nr(2,2)).and.(iz >= m%nr(1,3).and.iz <= m%nr(2,3))
  end select

contains

  logical function in_cylinder()
    real(r8) r, x
    r = sqrt(real(iy**2 + iz**2, r8))*m%h(1)
    x = ix*m%h(1)
    
    in_cylinder = (r<=m%rsize+DELTA_R .and. abs(x)<=m%xsize+DELTA_R)
  end function in_cylinder

end function in_mesh

subroutine lookup_table(m, d)
  type(mesh_type), intent(inout) :: m
  type(mesh_derivatives_type), intent(in) :: d

  integer :: i, j, ix(3), ik, in, nl, ng(3)
  
  allocate(m%der_lookup(m%np))

  do i = 1, m%np
    ix(:) = m%Lxyz(:, i)
    
    ! first we count the points (zero is triple counted, but we do not care for now ;)
    ! nl counts the points for the laplacian and ng for the gradient
    nl = 0; ng = 0
    do j = 1, conf%dim
      do ik = -d%norder, d%norder
        ix(j) = ix(j) + ik
        if(index(ix, sign(j, ik)) > 0) then
          nl    = nl    + 1
          ng(j) = ng(j) + 1
        end if
        ix(j) = ix(j) - ik
      end do
    end do

    ! allocate
    m%der_lookup(i)%lapl_n = nl
    m%der_lookup(i)%grad_n = ng
    ng(1) = maxval(ng)
    allocate(m%der_lookup(i)%lapl_i(nl), m%der_lookup(i)%grad_i(ng(1), 1:conf%dim))
    allocate(m%der_lookup(i)%lapl_w(nl), m%der_lookup(i)%grad_w(ng(1), 1:conf%dim))
    
    ! fill in the table
    nl = 0; ng = 0
    do j = 1, conf%dim
      do ik = -d%norder, d%norder
        ix(j) = ix(j) + ik
        
        in = index(ix, sign(j, ik))
        if(in > 0) then
          nl    = nl    + 1
          ng(j) = ng(j) + 1
          
          m%der_lookup(i)%lapl_i(nl) = in
          m%der_lookup(i)%lapl_w(nl) = d%dlidfj(abs(ik))/m%h(j)**2
          
          m%der_lookup(i)%grad_i(ng(j), j) = in
          m%der_lookup(i)%grad_w(ng(j), j) = d%dgidfj(ik)/m%h(j)
          
        end if
        
        ix(j) = ix(j) - ik
      end do
    end do
    
  end do
  
contains

  ! this function takes care of the boundary conditions
  ! for a given x,y,z it returns the true index of the point
  function index(ix_, dir)
    integer :: index
    integer, intent(in) :: ix_(3), dir

    integer :: i, ix(3)

    ix = ix_ ! make a local copy that we can change
    index = 1
    do i = 1, conf%dim
      if(ix(i) < m%nr(1, i)) then       ! first look left
        if(i <= conf%periodic_dim) then ! fold point
          ix(i) = ix(i) - 2*m%nr(1, i)
        else
          ix(i) = m%nr(1, i)
          index = 0
        end if
      else if(ix(i) > m%nr(2, i)) then  ! the same, but on the right
        if(i <= conf%periodic_dim) then
          ix(i) = ix(i) + 2*m%nr(1, i)
        else
          ix(i) = m%nr(2, i)
          index = 0
        end if
      end if
    end do

    if(index.ne.0) index = m%Lxyz_inv(ix(1), ix(2), ix(3))

    if(index==0.and.conf%boundary_zero_derivative) then
      do
        index = m%Lxyz_inv(ix(1), ix(2), ix(3))
        if(index.ne.0) exit
        ix(abs(dir)) = ix(abs(dir)) - sign(1, dir)
      end do
    end if

  end function index

end subroutine lookup_table
