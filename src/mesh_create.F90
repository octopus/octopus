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

#if defined(HAVE_METIS)
subroutine mesh_partition(m, Lxyz_tmp)
  type(mesh_type),      intent(inout) :: m
  integer,              intent(in)    :: Lxyz_tmp(:,:,:)

  integer :: i, ix, iy, iz, ne
  integer :: etype, edgecut
  integer, allocatable :: elmnts(:), epart(:), npart(:)

  ASSERT(m%sb%dim==2.or.m%sb%dim==3)

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
  allocate(elmnts(4*(m%sb%dim-1)*ne))

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

           if(m%sb%dim <= 2) cycle
           elmnts(i+0) = m%Lxyz_inv(ix,   iy,   iz+1)
           elmnts(i+1) = m%Lxyz_inv(ix+1, iy,   iz+1)
           elmnts(i+2) = m%Lxyz_inv(ix+1, iy+1, iz+1)
           elmnts(i+3) = m%Lxyz_inv(ix,   iy+1, iz+1)
           i = i + 4
        end do
     end do
  end do

  allocate(epart(ne), npart(m%np))
  if(m%sb%dim == 2) then
     etype = 4 ! quadrilaterals
  else
     etype = 3 ! hexahedra
  end if
  call oct_METIS_Part_Mesh_Nodal(ne, m%np, elmnts, etype, 1, mpiv%numprocs, edgecut, epart, npart)

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

    if(m%sb%dim <= 2) return

    OK = OK.and. &
         Lxyz_tmp(ix,   iy,   iz+1).eq.1.and. &   ! 5
         Lxyz_tmp(ix+1, iy,   iz+1).eq.1.and. &   ! 6
         Lxyz_tmp(ix+1, iy+1, iz+1).eq.1.and. &   ! 7
         Lxyz_tmp(ix,   iy+1, iz+1).eq.1          ! 8

  end function point_OK
end subroutine mesh_partition
#endif

subroutine mesh_create_xyz(sb, m, cv, geo)
  type(simul_box_type),  intent(in)    :: sb
  type(mesh_type),       intent(inout) :: m
  type(curvlinear_type), intent(in)    :: cv
  type(geometry_type),   intent(in)    :: geo

  integer :: i, j, k, il, ix, iy, iz
  integer, allocatable :: Lxyz_tmp(:,:,:)
  FLOAT,   allocatable :: x(:,:,:,:)
  FLOAT :: chi(3)

  call push_sub('mesh_create.mesh_create_xyz')

  ! enlarge mesh in the non-periodic dimensions
    m%nr(1,:) = m%nr(1,:) - m%enlarge(:)
    m%nr(2,:) = m%nr(2,:) + m%enlarge(:)

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(  Lxyz_tmp(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(       x(3,m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))

  m%Lxyz_inv(:,:,:) = 0
  Lxyz_tmp(:,:,:)   = 0
  x(:,:,:,:)        = M_ZERO

  ! We label 2 the points inside the mesh + enlargement

  do ix = m%nr(1,1), m%nr(2,1)
    chi(1) = real(ix, PRECISION) * m%h(1) + sb%box_offset(1)

    do iy = m%nr(1,2), m%nr(2,2)
      chi(2) = real(iy, PRECISION) * m%h(2) + sb%box_offset(2)

      do iz = m%nr(1,3), m%nr(2,3)
        chi(3) = real(iz, PRECISION) * m%h(3) + sb%box_offset(3)

        call curvlinear_chi2x(sb, geo, cv, chi(:), x(:, ix, iy, iz))

        if(simul_box_in_box(sb, geo, x(:, ix, iy, iz))) then
          do i = -m%enlarge(1), m%enlarge(1)
            do j = -m%enlarge(2), m%enlarge(2)
              do k = -m%enlarge(3), m%enlarge(3)
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
        if(simul_box_in_box(sb, geo, x(:, ix, iy, iz))) Lxyz_tmp(ix, iy, iz) = 1
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

  call mesh_get_vol_pp(sb, geo, cv, m)

!#if defined(HAVE_MPI) && defined(HAVE_METIS)
!  call mesh_partition(m, Lxyz_tmp)
!#endif

  deallocate(Lxyz_tmp)
  deallocate(x)

  call pop_sub()

end subroutine mesh_create_xyz


! calculate the volume of integration
subroutine mesh_get_vol_pp(sb, geo, cv, mesh)
  type(simul_box_type),  intent(in)    :: sb
  type(geometry_type),   intent(in)    :: geo
  type(curvlinear_type), intent(in)    :: cv
  type(mesh_type),       intent(inout) :: mesh

  integer :: i
  FLOAT :: f, chi(sb%dim)

  f = M_ONE
  do i = 1, sb%dim
    f = f*mesh%h(i)
  end do

  allocate(mesh%vol_pp(mesh%np_tot))
  mesh%vol_pp(:) = f

  do i = 1, mesh%np_tot
    chi(1:sb%dim) = mesh%Lxyz(i, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i,:), chi)
  end do

end subroutine mesh_get_vol_pp


! this function takes care of the boundary conditions
! for a given x,y,z it returns the true index of the point

! WARNING: have to get rid of dir, otherwise will not work
integer function mesh_index(m, ix_, dir) result(index)
  type(mesh_type),      intent(in) :: m
  integer,              intent(in) :: ix_(:), dir

  integer :: i, ix(3)  ! ix has to go until 3, not sb%dim

  ix = 0
  ix(1:m%sb%dim) = ix_(1:m%sb%dim) ! make a local copy that we can change

  index = 1
  do i = 1, m%sb%dim
    if(ix(i) < m%nr(1, i)) then       ! first look left
      if(i <= m%sb%periodic_dim) then ! fold point
        ix(i) = ix(i) + abs(m%nr(2,i) - m%nr(1,i) + 1)
      else
        ix(i) = m%nr(1, i)
        index = 0
      end if
    else if(ix(i) > m%nr(2, i)) then  ! the same, but on the right
      if(i <= m%sb%periodic_dim) then
        ix(i) = ix(i) - abs(m%nr(2,i) - m%nr(1,i) + 1)
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
