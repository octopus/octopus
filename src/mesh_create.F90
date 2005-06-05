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
!!
!! $Id$

#define DELTA_R CNST(1e-12)

subroutine mesh_init(m, geo)
  type(mesh_type),     intent(inout) :: m
  type(geometry_type), pointer       :: geo

  ! some local stuff
  FLOAT :: def_h, def_rsize
  integer :: i, ix, iy, iz, ii(3)
  
  call push_sub('mesh_init')

  call geometry_grid_defaults(geo, def_h, def_rsize)
 
  m%geo => geo            ! keep a pointer to the geometry

  call read_misc()          ! miscellany stuff
  call read_box()           ! parameters defining the simulation box
  call read_spacing ()      ! parameters defining the (canonical) spacing
  call read_box_offset()    ! parameters defining the offset of the origin
  call build_lattice()      ! build lattice vectors

  ! initialize curvlinear coordinates
  m%use_curvlinear = curvlinear_init(m%lsize(:), m%cv)

  call adjust_nr()          ! find out the extension of the simulation box

  call pop_sub()
  
contains
  subroutine read_misc()

    call loct_parse_float(check_inp('DoubleFFTParameter'), M_TWO, m%fft_alpha)
    if (m%fft_alpha < M_ONE .or. m%fft_alpha > M_THREE ) then
      write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
      message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
      call write_fatal(2)
    end if
    
  end subroutine read_misc

  subroutine read_box()
    integer(POINTER_SIZE) :: blk

    ! Read box shape.
    call loct_parse_int(check_inp('BoxShape'), MINIMUM, m%box_shape)
    if(m%box_shape<1 .or. m%box_shape>4) then
      write(message(1), '(a,i2,a)') "Input: '", m%box_shape, "' is not a valid BoxShape"
      message(2) = '(1 <= Box_Shape <= 4)'
      call write_fatal(2)
    end if
    
    select case(m%box_shape)
    case(SPHERE,MINIMUM)
      if(conf%dim>1 .and. conf%periodic_dim>0) then
        message(1) = 'Spherical or minimum mesh is not allowed for periodic systems'
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
    if(conf%dim==1.and.m%box_shape /= PARALLELEPIPED) m%box_shape=SPHERE

    m%rsize = -M_ONE
    if(m%box_shape == MINIMUM.and.def_rsize>M_ZERO) m%rsize = def_rsize/units_inp%length%factor

    if(m%box_shape == SPHERE.or.m%box_shape == CYLINDER.or.m%box_shape == MINIMUM) then
      call loct_parse_float(check_inp('radius'), m%rsize, m%rsize)
      if(m%rsize < 0) then
        message(1) = "Either:"
        message(2) = "   *) variable 'radius' is not defined"
        message(3) = "   *) your input for 'radius' is negative"
        message(4) = "   *) I can't find a suitable default"
        call write_fatal(4)
      end if
      m%rsize = m%rsize * units_inp%length%factor
      if(def_rsize>M_ZERO) call check_def(def_rsize, m%rsize, 'radius')
    end if

    if(m%box_shape == CYLINDER) then
      call loct_parse_float(check_inp('xlength'), M_ONE/units_inp%length%factor, m%xsize)
      m%xsize = m%xsize * units_inp%length%factor
      m%lsize(1) = m%xsize
      if(def_rsize>M_ZERO.and.conf%periodic_dim==0) call check_def(def_rsize, m%xsize, 'xlength')
    end if

    m%lsize = -M_ONE
    if(m%box_shape == PARALLELEPIPED) then
      
      if(loct_parse_block(check_inp('lsize'), blk) == 0) then
        if(loct_parse_block_cols(blk,0) < conf%dim) then
            message(1) = 'Size of Block "lsize" does not match number of dimensions'
            call write_fatal(1)
        endif
      else
        message(1) = 'Block "lsize" not found in input file.'
        call write_fatal(1)
      endif
     
      do i = 1, conf%dim
        call loct_parse_block_float(blk, 0, i-1, m%lsize(i))
        if(def_rsize>M_ZERO.and.conf%periodic_dim<i) call check_def(def_rsize, m%lsize(i), 'lsize')
      end do
      m%lsize = m%lsize*units_inp%length%factor

      call loct_parse_block_end(blk)
    end if
    
    ! fill in lsize structure
    select case(m%box_shape)
    case(SPHERE)
      m%lsize(1:conf%dim) = m%rsize
    case(CYLINDER)
      m%lsize(1)          = m%xsize
      m%lsize(2:conf%dim) = m%rsize
    case(MINIMUM)
      do i = 1, conf%dim
        m%lsize(i)        = maxval(geo%atom(:)%x(i)) + m%rsize
      end do
    end select
    
  end subroutine read_box

  subroutine read_spacing()
    integer :: i
    integer(POINTER_SIZE) :: blk

    m%h = -M_ONE
    select case(m%box_shape)
    case(SPHERE,CYLINDER,MINIMUM)
      call loct_parse_float(check_inp('spacing'), m%h(1), m%h(1))
      m%h(1:conf%dim) = m%h(1)

    case(PARALLELEPIPED)
      if(loct_parse_block(check_inp('spacing'), blk) == 0) then
        do i = 1, conf%dim
          call loct_parse_block_float(blk, 0, i-1, m%h(i))
        end do
        call loct_parse_block_end(blk)
      else
        message(1) = '"Spacing" is a block if BoxShape == parallelepiped.'
        call write_fatal(1)
      endif
    end select

    do i = 1, conf%dim
      m%h(i) = m%h(i)*units_inp%length%factor
      if(m%h(i) < M_ZERO) then
        if(def_h > M_ZERO) then
          m%h(i) = def_h
          write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", i, &
              ") [", trim(units_out%length%abbrev), "] = ",                 &
              m%h(i)/units_out%length%factor
          call write_info(1)
        else
          message(1) = 'Either:'
          message(2) = "   *) variable 'spacing' is not defined"
          message(3) = "   *) your input for 'spacing' is negative"
          message(4) = "   *) I can't find a suitable default"
          call write_fatal(4)
        end if
      end if
      if(def_rsize>M_ZERO) call check_def(m%h(i), def_rsize, 'spacing')
    end do

  end subroutine read_spacing

  subroutine read_box_offset()
    integer :: i
    integer(POINTER_SIZE) :: blk

    m%box_offset = M_ZERO
    select case(m%box_shape)
    case(PARALLELEPIPED)
      if(loct_parse_block(check_inp('BoxOffset'), blk) == 0) then
        do i = 1, conf%dim
          call loct_parse_block_float(blk, 0, i-1, m%box_offset(i))
        end do
        call loct_parse_block_end(blk)
      else
        message(1) = 'Block "BoxOffset" not properly defined in input file.'
        message(2) = 'Assuming zero offset in all directions.'
        call write_warning(2)
        m%box_offset = M_ZERO
      endif
    end select
        
  end subroutine read_box_offset

  subroutine check_def(var, def, text)
    FLOAT, intent(in) :: var, def
    character(len=*), intent(in) :: text

    if(var > def) then
      write(message(1), '(3a)') "The value for '", text, "' does not match the recommended value"
      write(message(2), '(f8.3,a,f8.3)') var, ' > ', def
      call write_warning(2)
    end if
  end subroutine check_def

  subroutine build_lattice()
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

  end subroutine build_lattice

  ! set nr and adjust the mesh so that:
  ! 1) the new grid exactly fills the box;
  ! 2) the new mesh is not larger than the user defined mesh.
  subroutine adjust_nr()
    integer :: i, j
    FLOAT   :: x(conf%dim), chi(conf%dim)
    logical :: out

    m%nr = 0
    do i = 1, conf%dim
      chi(:) = M_ZERO; j = 0
      out = .false.
      do while(.not.out)
        j      = j + 1
        chi(i) = j*m%h(i)
        call curvlinear_chi2x(m%cv, m%geo, chi(:), x(:))
        out = (x(i) > m%lsize(i))
      end do
      m%nr(2, i) = j - 1
    end do

    ! we have a symmetric mesh (for now)
    m%nr(1,:) = -m%nr(2,:)
    
    ! we have to ajust a couple of things for the periodic directions
    do i = 1, conf%periodic_dim
      m%h(i)     = m%lsize(i)/real(m%nr(2, i))
      m%nr(2, i) = m%nr(2, i) - 1
    end do

    m%l(:) = m%nr(2, :) - m%nr(1, :) + 1
  end subroutine adjust_nr

end subroutine mesh_init

#if defined(HAVE_METIS)
subroutine mesh_partition(m, Lxyz_tmp)
  type(mesh_type), intent(inout) :: m
  integer, pointer :: Lxyz_tmp(:,:,:)
  
  integer :: i, ix, iy, iz, ne
  integer :: etype, edgecut
  integer, allocatable :: elmnts(:), epart(:), npart(:)

  ASSERT(conf%dim==2.or.conf%dim==3)

  ! let us count how many elements we have
  ne = 0
  do ix = m%nr(1,1), m%nr(2,1) - 1
    do iy = m%nr(1,2), m%nr(2,2) - 1
      do iz = m%nr(1,3), max(0, m%nr(2,3) - 1) ! max is necessary when running in 2 dimensions
        if(.not.point_OK(ix, iy, iz)) cycle
        ne = ne + 1
      end do
    end do
  end do

  ! allocate space necessary for the elements
  allocate(elmnts(4*(conf%dim-1)*ne))

  ! create elements array
  i = 1
  do ix = m%nr(1,1), m%nr(2,1) - 1
    do iy = m%nr(1,2), m%nr(2,2) - 1
      do iz = m%nr(1,3), max(0, m%nr(2,3) - 1) ! max is necessary when running in 2 dimensions
        if(.not.point_OK(ix, iy, iz)) cycle
        elmnts(i+0) = m%Lxyz_inv(ix,   iy,   iz)
        elmnts(i+1) = m%Lxyz_inv(ix+1, iy,   iz)
        elmnts(i+2) = m%Lxyz_inv(ix+1, iy+1, iz)
        elmnts(i+3) = m%Lxyz_inv(ix,   iy+1, iz)
        i = i + 4

        if(conf%dim <= 2) cycle
        elmnts(i+0) = m%Lxyz_inv(ix,   iy,   iz+1)
        elmnts(i+1) = m%Lxyz_inv(ix+1, iy,   iz+1)
        elmnts(i+2) = m%Lxyz_inv(ix+1, iy+1, iz+1)
        elmnts(i+3) = m%Lxyz_inv(ix,   iy+1, iz+1)
        i = i + 4
      end do
    end do
  end do

  allocate(epart(ne), npart(m%np))
  if(conf%dim == 2) then
    etype = 4 ! quadrilaterals
  else
    etype = 3 ! hexahedra
  end if
  call oct_METIS_partition(ne, m%np, elmnts, etype, 1, mpiv%numprocs, edgecut, epart, npart)

  if(mpiv%node == 0) then
    do i = 1, m%np
      write(12,*) m%x(i,1:2), npart(i)
    end do
  end if
  close(12)
  !stop

  deallocate(epart, npart)

  contains
    logical function point_OK(ix, iy, iz) result(OK)
      integer, intent(in) :: ix, iy, iz

      ! check if all 8 vertices are inside the mesh
      ! the labelling is consistent with METIS requirements
      OK = &
         Lxyz_tmp(ix,   iy,   iz  ).eq.1.and. &   ! 1
         Lxyz_tmp(ix+1, iy,   iz  ).eq.1.and. &   ! 2
         Lxyz_tmp(ix+1, iy+1, iz  ).eq.1.and. &   ! 3
         Lxyz_tmp(ix,   iy+1, iz  ).eq.1          ! 4

      if(conf%dim <= 2) return

      OK = OK.and. &
         Lxyz_tmp(ix,   iy,   iz+1).eq.1.and. &   ! 5
         Lxyz_tmp(ix+1, iy,   iz+1).eq.1.and. &   ! 6
         Lxyz_tmp(ix+1, iy+1, iz+1).eq.1.and. &   ! 7
         Lxyz_tmp(ix,   iy+1, iz+1).eq.1          ! 8
      
    end function point_OK
end subroutine mesh_partition
#endif

subroutine mesh_create_xyz(m, enlarge)
  type(mesh_type),     intent(inout) :: m
  integer,             intent(in)    :: enlarge

  integer :: i, j, k, il, ix, iy, iz, ey, ez
  integer, pointer :: Lxyz_tmp(:,:,:)
  FLOAT :: chi(3)
  FLOAT, pointer :: x(:,:,:,:)

  call push_sub('mesh_create_xyz')

! enlarge mesh in the non-periodic dimensions
  do i = conf%periodic_dim+1, conf%dim
    m%nr(1,i) = m%nr(1,i) - enlarge
    m%nr(2,i) = m%nr(2,i) + enlarge
  end do

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(  Lxyz_tmp(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(       x(3,m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  
  m%Lxyz_inv(:,:,:) = 0
  Lxyz_tmp(:,:,:)   = 0
  x(:,:,:,:)        = M_ZERO

  ! We label 2 the points inside the mesh + enlargement
  ey = 0; ez = 0
  if(conf%dim > 1) ey = enlarge
  if(conf%dim > 2) ez = enlarge

  do ix = m%nr(1,1), m%nr(2,1)
    chi(1) = real(ix, PRECISION) * m%h(1) + m%box_offset(1)

    do iy = m%nr(1,2), m%nr(2,2)
      chi(2) = real(iy, PRECISION) * m%h(2) + m%box_offset(2)

      do iz = m%nr(1,3), m%nr(2,3)
        chi(3) = real(iz, PRECISION) * m%h(3) + m%box_offset(3)

        call curvlinear_chi2x(m%cv, m%geo, chi(:), x(:, ix, iy, iz))

        if(in_mesh(m, x(:, ix, iy, iz))) then
          do i = -enlarge, enlarge
            do j = -ey, ey
              do k = -ez, ez
                if(  &
                   ix+i>=m%nr(1,1).and.ix+i<=m%nr(2,1).and. &
                   iy+j>=m%nr(1,2).and.iy+j<=m%nr(2,2).and. &
                   iz+k>=m%nr(1,3).and.iz+k<=m%nr(2,3)) Lxyz_tmp(ix+i, iy+j, iz+k) = 2
              end do
            end do
          end do
        end if

      end do
    end do
  end do

  ! we label 1 the points inside the mesh, and we count the points
  il = 0
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(in_mesh(m, x(:, ix, iy, iz))) Lxyz_tmp(ix, iy, iz) = 1
        if(Lxyz_tmp(ix, iy, iz) > 0) il = il + 1
      end do
    end do
  end do
  m%np_tot = il

  allocate(m%Lxyz(m%np_tot, 3), m%x(m%np_tot, 3))

  ! first we fill the points in the inner mesh
  il = 0
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 1) then
          il = il + 1
          m%Lxyz(il, 1) = ix
          m%Lxyz(il, 2) = iy
          m%Lxyz(il, 3) = iz
          m%Lxyz_inv(ix,iy,iz) = il
          m%x(il,:) = x(:,ix,iy,iz)
        end if
      enddo
    enddo
  enddo
  m%np = il

  ! and now the points from the enlargement
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 2) then
          il = il + 1
          m%Lxyz(il, 1) = ix
          m%Lxyz(il, 2) = iy
          m%Lxyz(il, 3) = iz
          m%Lxyz_inv(ix,iy,iz) = il
          m%x(il,:) = x(:,ix,iy,iz)
        end if
      end do
    end do
  end do

  call get_vol_pp()

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  call mesh_partition(m, Lxyz_tmp)
#endif

  deallocate(Lxyz_tmp); nullify(Lxyz_tmp)
  deallocate(x)

  call pop_sub()
contains

  ! calculate the volume of integration
  subroutine get_vol_pp()
    integer :: i
    FLOAT :: f, chi(conf%dim)

    f = M_ONE
    do i = 1, conf%dim
      f = f*m%h(i)
    end do

    allocate(m%vol_pp(m%np))
    m%vol_pp(:) = f

    do i = 1, m%np
      chi(1:conf%dim) = m%Lxyz(i, 1:conf%dim)*m%h(1:conf%dim)
      m%vol_pp(i) = m%vol_pp(i)*curvlinear_det_Jac(m%cv, m%geo, m%x(i,:), chi)
    end do

  end subroutine get_vol_pp

end subroutine mesh_create_xyz


logical function in_mesh(m, x)
  type(mesh_type),     intent(IN) :: m
  FLOAT,               intent(in) :: x(:) ! x(3)
    
  select case(m%box_shape)
  case(SPHERE)
    in_mesh = (sqrt(sum(x**2)) <= m%rsize+DELTA_R)
  case(CYLINDER)
    in_mesh = in_cylinder()
  case(MINIMUM)
    in_mesh = in_minimum()
  case(PARALLELEPIPED)
    in_mesh =  &
       (x(1) >= m%nr(1,1)*m%h(1)+m%box_offset(1).and.x(1) <= m%nr(2,1)*m%h(1)+m%box_offset(1)).and. &
       (x(2) >= m%nr(1,2)*m%h(2)+m%box_offset(2).and.x(2) <= m%nr(2,2)*m%h(2)+m%box_offset(2)).and. &
       (x(3) >= m%nr(1,3)*m%h(3)+m%box_offset(3).and.x(3) <= m%nr(2,3)*m%h(3)+m%box_offset(3))
  end select

contains

  logical function in_cylinder()
    FLOAT :: r
    r = sqrt(x(2)**2 + x(3)**2)
    
    in_cylinder = (r<=m%rsize+DELTA_R .and. abs(x(1))<=m%xsize+DELTA_R)
  end function in_cylinder

  logical function in_minimum()
    integer :: i
    FLOAT   :: r
    
    in_minimum = .false.
    do i = 1, m%geo%natoms
      r = sqrt(sum((x(:) - m%geo%atom(i)%x(:))**2))
      if(r<=m%rsize+DELTA_R) then
        in_minimum = .true.
        exit
      end if
    end do
    
  end function in_minimum
end function in_mesh


! this function takes care of the boundary conditions
! for a given x,y,z it returns the true index of the point

! WARNING: have to get rid of dir, otherwise will not work
function mesh_index(m, ix_, dir) result(index)
    type(mesh_type), intent(in) :: m
    integer :: index
    integer, intent(in) :: ix_(:), dir

    integer :: i, ix(3)  ! ix has to go until 3, not conf%dim

    ix = 0
    ix(1:conf%dim) = ix_(1:conf%dim) ! make a local copy that we can change

    index = 1
    do i = 1, conf%dim
      if(ix(i) < m%nr(1, i)) then       ! first look left
        if(i <= conf%periodic_dim) then ! fold point
          ix(i) = ix(i) + abs(m%nr(2,i) - m%nr(1,i))
        else
          ix(i) = m%nr(1, i)
          index = 0
        end if
      else if(ix(i) > m%nr(2, i)) then  ! the same, but on the right
        if(i <= conf%periodic_dim) then
          ix(i) = ix(i) - abs(m%nr(2,i) - m%nr(1,i))
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
    
end function mesh_index
