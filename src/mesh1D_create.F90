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
! assumes that the box_shape, h, rsize/zsize of m are filled
! if norder is present, then set nk, Kx, Ky and Kz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh1D_create(m)!, natoms, atom)
  type(mesh_type), intent(inout) :: m
  !integer, intent(in), optional :: natoms
  !type(atom_type), pointer, optional :: atom(:)
  
  ! some local stuff
  character(len=40) :: err
  integer :: next = 0
  integer :: il, ik, ix, iy, iz
  integer, pointer :: Lxyz(:,:,:)
  logical :: b

  sub_name = 'mesh3D_create'; call push_sub()

  ! Always a "sphere".
  m%box_shape = 1

  call oct_parse_double(C_string('spacing'),0.6_r8/units_inp%length%factor, m%h(1))
  m%h(1) = m%h(1)*units_inp%length%factor
  if (m%h(1)<0.01_r8 .or. m%h(1)>2.0_r8) then
      write(err, '(a,f10.5,a,f10.5,a,f10.5,a)') '(',m%h(1),',',m%h(2),',',m%h(3),')'
      message(1) = "Input: '"//trim(err)//"' is not a valid spacing"
      message(2) = '(0.01 <= Spacing [b] <= 2)'
      call write_fatal(2)
  end if    
  m%h(2) = 0.0_r8; m%h(3) = 0.0_r8
  m%iso = .true. ! This is obvious in 1D...
  m%vol_pp = m%h(1) 

  call oct_parse_double(C_string('radius'), 20.0_r8/units_inp%length%factor, m%rsize)
  m%rsize = m%rsize * units_inp%length%factor
    if (m%rsize<1.0_r8 .or. m%rsize>500.0_r8) then
      write(err, *) m%rsize
      message(1) = "Input: '"//trim(err)//"' is not a valid radius"
      message(2) = '(1 <= radius [b] <= 500)'
      call write_fatal(2)
    end if

  m%nr(1) = int(m%rsize/m%h(1)); m%nr(2) = 0; m%nr(3) = 0
  m%nx(1) = m%nr(1)+m%d%norder; m%nx(2) = 0; m%nx(3) = 0

  ! From now on, all this may be simplified in 1D (some day).

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
          b = in_sphere(m, ix, iy, iz)
          if(b) Lxyz(ix-m%d%norder:ix+m%d%norder, iy, iz) = 2
      enddo
    enddo
  enddo

  do ix = -m%nr(1), m%nr(1)
    do iy = -m%nr(2), m%nr(2)
      do iz = -m%nr(3), m%nr(3)
          b = in_sphere(m, ix, iy, iz)
        if(b) Lxyz(ix, iy, iz) = 1
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

  allocate(m%Lx(il), m%Ly(il), m%Lz(il), m%Kx(ik), m%Ky(ik), m%Kz(ik))

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

  logical function in_sphere(m, ix, iy, iz)
    type(mesh_type), intent(IN) :: m
    integer, intent(in) :: ix, iy, iz
    
    !in_sphere = (sqrt(real(ix**2+iy**2+iz**2,r8))*m%H(1) <= m%rsize)
    in_sphere = ( abs(ix*m%h(1)) <= m%rsize)
    return
  end function in_sphere
  
end subroutine mesh1D_create
