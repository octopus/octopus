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

#include "global.h"

module mesh
  use global
  use messages
  use lib_oct_parser
  use units
  use math
  use geometry
  use curvlinear
  use simul_box
  use lib_adv_alg
  
  implicit none
  
  private
  public ::  &
     mesh_type, &
     mesh_init, &
     mesh_end, &
     mesh_double_box, &
     mesh_index, &
     mesh_inborder, &
     mesh_r, &
     mesh_gcutoff, &
     mesh_write_info

  type mesh_type
    type(simul_box_type), pointer :: sb ! The simulation box
    logical :: use_curvlinear
    
    FLOAT :: h(3)           ! the (constant) spacing between the points
    
    integer  :: np         ! number of points in mesh
    integer  :: np_tot     ! total number of points including ghost points
    
    integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
    integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz
    
    ! some other vars
    integer :: nr(2,3)  ! dimensions of the box where the points are contained
    integer :: l(3)     ! literally n(2,:) - n(1,:) + 1
    
    FLOAT, pointer :: x(:,:)    ! the points
    FLOAT, pointer :: vol_pp(:) ! element of volume for integrations
    
  end type mesh_type

contains

  subroutine mesh_init(m, sb, cv, geo, enlarge)
    type(mesh_type),       intent(inout) :: m
    type(simul_box_type),  target        :: sb
    type(curvlinear_type), intent(in)    :: cv
    type(geometry_type),   intent(in)    :: geo
    integer,               intent(in)    :: enlarge
    
    call push_sub('mesh_init')

    m%sb => sb
    m%use_curvlinear = cv%method.ne.CURV_METHOD_UNIFORM
    m%h = sb%h ! this number can change in the following

    call adjust_nr()          ! find out the extension of the simulation box
    call mesh_create_xyz(m, sb, cv, geo, enlarge)

    call pop_sub()
  contains
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
          call curvlinear_chi2x(cv, geo, chi(:), x(:))
          out = (x(i) > sb%lsize(i))
        end do
        m%nr(2, i) = j - 1
      end do
      
      ! we have a symmetric mesh (for now)
      m%nr(1,:) = -m%nr(2,:)
    
      ! we have to ajust a couple of things for the periodic directions
      do i = 1, sb%periodic_dim
        m%h(i)     = sb%lsize(i)/real(m%nr(2, i))
        m%nr(2, i) = m%nr(2, i) - 1
      end do

      m%l(:) = m%nr(2, :) - m%nr(1, :) + 1
    end subroutine adjust_nr

  end subroutine mesh_init


  ! finds the dimension of a box doubled in the non-periodic dimensions
  subroutine mesh_double_box(m, sb, db)
    type(mesh_type),      intent(in)  :: m
    type(simul_box_type), intent(in) :: sb
    integer,              intent(out) :: db(3)
    
    integer :: i
    
    db = 1

    ! OLD: I let it here because maybe I revert to this method later
    ! double mesh with 2n + 1 points
    !  do i = conf%periodic_dim + 1, conf%dim
    !    db(i) = nint(m%fft_alpha*(m%nr(2,i)-m%nr(1,i))) + 1
    !  end do
    !  do i = 1, conf%periodic_dim ! for periodic systems
    !    db(i) = m%nr(2,i) - m%nr(1,i) + 1
    !  end do
    
    ! NEW
    ! double mesh with 2n points
    do i = 1, sb%periodic_dim 
      db(i) = m%l(i)
    end do
    do i = sb%periodic_dim + 1, conf%dim
      db(i) = nint(m%sb%fft_alpha*(m%l(i)-1)) + 1
    end do
    
  end subroutine mesh_double_box
  
  subroutine mesh_write_info(m, unit)
    type(mesh_type), intent(IN) :: m
    integer,         intent(in) :: unit
    
#ifdef HAVE_MPI
    if(mpiv%node .ne. 0) return
#endif
    if(unit==stdout.and.conf%verbose<VERBOSE_NORMAL) return
    
    call push_sub('mesh_write_info')
    
    write(unit,'(a)') 'Main mesh:'
    
    write(unit,'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)')     &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',           &
       '(', m%h(1)/units_out%length%factor, ',',                       &
            m%h(2)/units_out%length%factor, ',',                       &
            m%h(3)/units_out%length%factor, ')',                       &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ',   &
       m%vol_pp(1)/units_out%length%factor**3

    write(unit,'(a, i6)') '  # inner mesh = ', m%np

    write(unit,'(3a,f9.3,a)') '  Grid Cutoff [',trim(units_out%energy%abbrev),'] = ', &
       (M_PI**2/(M_TWO*maxval(m%h)**2))/units_out%energy%factor
    
    call pop_sub()
  end subroutine mesh_write_info


  subroutine mesh_r(m, i, r, a, x)
    type(mesh_type), intent(IN)  :: m
    integer,         intent(in)  :: i
    FLOAT,           intent(out) :: r
    FLOAT,           intent(in),  optional :: a(:) ! a(conf%dim)
    FLOAT,           intent(out), optional :: x(:) ! x(conf%dim)
    
    FLOAT :: xx(conf%dim)
    
    xx(1:conf%dim) = m%x(i,1:conf%dim)
    if(present(a)) xx(1:conf%dim) = xx(1:conf%dim) - a(1:conf%dim)
    r = sqrt(dot_product(xx, xx))

    if(present(x)) x(1:conf%dim) = xx(1:conf%dim)
  end subroutine mesh_r

  !/*---------------------------------------------------------------------------------
  ! Finds out if a given point of a mesh belongs to the "border" of the mesh.
  ! A point belongs to the border of the mesh if it is too close to any of the
  ! walls of the mesh. The criterium is set by input parameter "width".
  !
  ! m     : the mesh.
  ! i     : the point in the mesh.
  ! n     : on output, the number (0<=n<=3) of "walls" of the mesh that the point is
  !         too close to, in order to consider it belonging to a mesh.
  ! d     : the distances of the point to the walls, for each of the walls that the
  !         point is too close to.
  ! width : the width of the border.
  !
  ! So, if n>0, the point is in the border.
  !---------------------------------------------------------------------------------*/
  subroutine mesh_inborder(m, i, n, d, width)
    type(mesh_type), intent(in)  :: m
    integer,         intent(in)  :: i
    FLOAT,           intent(in)  :: width
    integer,         intent(out) :: n
    FLOAT,           intent(out) :: d(3)
    
    integer :: j
    FLOAT :: x(3), r, dd
    
    call mesh_r(m, i, r, x=x)
    n = 0
    select case(m%sb%box_shape)
    case(SPHERE)
      dd = r - (m%sb%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      end if
    case(CYLINDER)
      dd = sqrt(x(2)**2 + x(3)**2) - (m%sb%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      endif
      dd = abs(x(1)) - (m%sb%xsize - width)
      if(dd.gt.M_ZERO) then
        n = n + 1; d(n) = dd
      endif
    case(MINIMUM)
      message(1) = "Absorbing boundaries are not yet implemented for the 'minimum' box"
      call write_fatal(1)
    case(PARALLELEPIPED)
      do j = 1, conf%dim
        dd = abs(x(j)) - (m%sb%lsize(j) - width)
        if(dd.gt.M_ZERO) then
           n = n + 1; d(n) = dd
         end if
       end do
     end select
     
   end subroutine mesh_inborder


   !/*---------------------------------------------------------------------------------
   ! mesgh_gcutoff returns the "natural" band limitation of the grid m, in terms
   ! of the maximum G vector. For a cubic regular grid of spacing h is M_PI/h.
   !---------------------------------------------------------------------------------*/
   FLOAT function mesh_gcutoff(m) result(gmax)
     type(mesh_type), intent(in) :: m
     gmax = M_PI/(maxval(m%h))
   end function mesh_gcutoff
   

   subroutine mesh_end(m)
     type(mesh_type), intent(inout) :: m
     
     call push_sub('mesh_end')
     
     if(associated(m%Lxyz)) then
       deallocate(m%Lxyz, m%Lxyz_inv, m%x, m%vol_pp)
       nullify(m%Lxyz, m%Lxyz_inv, m%x, m%vol_pp)
     end if
     
     call pop_sub()
   end subroutine mesh_end

#include "mesh_create.F90"

end module mesh
