!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

module multigrid
  use curvlinear
  use geometry
  use mesh
  use functions

  implicit none

  private
  public ::                     &
    multigrid_level_type,       &
    multigrid_type,             &
    multigrid_init,             &
    multigrid_end,              &
    multigrid_fine2coarse,      &
    multigrid_coarse2fine

  integer, parameter, public :: &
    INJECTION  = 1,             &
    FULLWEIGHT = 2

  type multigrid_level_type
    type(mesh_type),  pointer  :: m
    type(f_der_type), pointer  :: f_der

    integer          ::  n_coarse
    integer, pointer :: to_coarse(:)

    integer          ::  n_fine
    integer          ::  n_fine1,       n_fine2,       n_fine4,       n_fine8
    integer, pointer :: to_fine1(:,:), to_fine2(:,:), to_fine4(:,:), to_fine8(:,:)
    integer, pointer :: fine_i(:)
  end type multigrid_level_type

  type multigrid_type
    integer                             :: n_levels
    type(multigrid_level_type), pointer :: level(:)
  end type multigrid_type


contains


  ! ---------------------------------------------------------
  subroutine multigrid_init(geo, cv, m, f_der, mgrid)
    type(geometry_type),      intent(in)  :: geo
    type(curvlinear_type),    intent(in)  :: cv
    type(mesh_type),  target, intent(in)  :: m
    type(f_der_type), target, intent(in)  :: f_der
    type(multigrid_type),     intent(out) :: mgrid

    integer :: i, n_levels, np

    !%Variable MultigridLevels
    !%Type integer
    !%Default 0 
    !%Section Hamiltonian
    !%Description
    !% Number of levels in the grid hierachy used for multigrid. Positive
    !% numbers indicate an absolute numbers of levels, negative
    !% numbers are substracted to maximum number of levels posible for
    !% the grid been used. Default is the maximum number of levels for
    !% the grid.
    !%Option max_levels 0
    !% Calculate the optimous number of levels for the grid.
    !%End

    call loct_parse_int(check_inp('MultigridLevels'), 0, n_levels)

    if ( n_levels <= 0 )then
      n_levels=n_levels-2
      np=m%np

      do while ( np > 1 )
        np=np/8
        n_levels=n_levels+1
      end do
    else
      n_levels=n_levels-1
    end if

    mgrid%n_levels = n_levels

    allocate(mgrid%level(0:n_levels))

    mgrid%level(0)%m     => m
    mgrid%level(0)%f_der => f_der

    mgrid%level(0)%n_fine = m%np
    allocate(mgrid%level(0)%fine_i(m%np))

    write(message(1), '(a,i3)') "Multigrid levels:", n_levels+1
    call write_info(1)

    do i = 1, mgrid%n_levels
      allocate(mgrid%level(i)%m, mgrid%level(i)%f_der)
      call multigrid_mesh_half(geo, cv, mgrid%level(i-1)%m, mgrid%level(i)%m)
      call get_transfer_tables(mgrid%level(i), mgrid%level(i-1)%m, mgrid%level(i)%m)

      call mesh_write_info(mgrid%level(i)%m, stdout)

      ! setup derivatives
      call f_der_init(mgrid%level(i)%f_der, m%sb, cv%method.ne.CURV_METHOD_UNIFORM)
      call f_der_build(m%sb, mgrid%level(i)%m, mgrid%level(i)%f_der)
    end do

  contains

    !/*---------------------------------------------------------------------------------
    ! Creates a mesh that has twice the spacing betwen the points than the in mesh.
    ! This is used in the multi-grid routines
    !---------------------------------------------------------------------------------*/
    subroutine multigrid_mesh_half(geo, cv, mesh_in, mesh_out)
      type(geometry_type),   intent(in)  :: geo
      type(curvlinear_type), intent(in)  :: cv
      type(mesh_type),       intent(in)  :: mesh_in
      type(mesh_type),       intent(inout) :: mesh_out

      mesh_out%sb             => mesh_in%sb
      mesh_out%use_curvlinear =  mesh_in%use_curvlinear

      mesh_out%h(:)    = 2*mesh_in%h(:)
      mesh_out%nr(:,:) = mesh_in%nr(:,:)/2
      mesh_out%l(:)    = mesh_out%nr(2, :) - mesh_out%nr(1, :) + 1

      mesh_out%enlarge = mesh_in%enlarge

      call mesh_init_stage_2(mesh_out%sb, mesh_out, geo, cv)
      call mesh_init_stage_3(mesh_out, geo, cv)

    end subroutine multigrid_mesh_half


    ! creates the lookup tables to go between the coarse and fine meshes
    subroutine get_transfer_tables(tt, fine, coarse)
      type(multigrid_level_type), intent(inout) :: tt
      type(mesh_type),            intent(in)  :: fine, coarse

      integer :: i, i1, i2, i4, i8
      integer :: x(3), mod2(3)

      tt%n_coarse = coarse%np
      allocate(tt%to_coarse(tt%n_coarse))

      do i = 1, tt%n_coarse
        tt%to_coarse(i) = fine%Lxyz_inv(2*coarse%Lxyz(i, 1), 2*coarse%Lxyz(i, 2), 2*coarse%Lxyz(i, 3))
      end do

      ! count
      tt%n_fine = fine%np
      allocate(tt%fine_i(tt%n_fine))

      tt%n_fine1 = 0
      tt%n_fine2 = 0
      tt%n_fine4 = 0
      tt%n_fine8 = 0
      do i = 1, tt%n_fine
        mod2 = mod(fine%Lxyz(i,:), 2)
        if(all(mod2.eq.0)) then
          tt%n_fine1 = tt%n_fine1 + 1
          tt%fine_i(i) = 1
        else if( &
          ((mod2(1).eq.0).and.(mod2(2).eq.0).and.(mod2(3).ne.0)).or. &
          ((mod2(3).eq.0).and.(mod2(1).eq.0).and.(mod2(2).ne.0)).or. &
          ((mod2(2).eq.0).and.(mod2(3).eq.0).and.(mod2(1).ne.0))) then
          tt%n_fine2 = tt%n_fine2 + 1
          tt%fine_i(i) = 2
        else if( &
          ((mod2(1).eq.0).and.(mod2(2).ne.0).and.(mod2(3).ne.0)).or. &
          ((mod2(3).eq.0).and.(mod2(1).ne.0).and.(mod2(2).ne.0)).or. &
          ((mod2(2).eq.0).and.(mod2(3).ne.0).and.(mod2(1).ne.0))) then
          tt%n_fine4 = tt%n_fine4 + 1
          tt%fine_i(i) = 4
        else if(all(mod2.ne.0)) then
          tt%n_fine8 = tt%n_fine8 + 1
          tt%fine_i(i) = 8
        end if
      end do

      ASSERT(tt%n_fine1+tt%n_fine2+tt%n_fine4+tt%n_fine8 == tt%n_fine)

      allocate(tt%to_fine1(1, tt%n_fine1), tt%to_fine2(2, tt%n_fine2), &
        tt%to_fine4(4, tt%n_fine4), tt%to_fine8(8, tt%n_fine8))

      ! and now build the tables
      i1 = 0;  i2 = 0;  i4 = 0;  i8 = 0
      do i = 1, fine%np
        x(1:3)    = fine%Lxyz(i, 1:3)/2
        mod2(1:3) = mod(fine%Lxyz(i, 1:3), 2)

        if((mod2(1).eq.0).and.(mod2(2).eq.0).and.(mod2(3).eq.0)) then
          i1 = i1 + 1
          tt%to_fine1(1, i1) = coarse%Lxyz_inv(x(1), x(2), x(3))

        else if((mod2(1).eq.0).and.(mod2(2).eq.0).and.(mod2(3).ne.0)) then
          i2 = i2 + 1
          tt%to_fine2(1, i2) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine2(2, i2) = coarse%Lxyz_inv(x(1), x(2), x(3) + mod2(3))
        else if((mod2(3).eq.0).and.(mod2(1).eq.0).and.(mod2(2).ne.0)) then
          i2 = i2 + 1
          tt%to_fine2(1, i2) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine2(2, i2) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3))
        else if((mod2(2).eq.0).and.(mod2(3).eq.0).and.(mod2(1).ne.0)) then
          i2 = i2 + 1
          tt%to_fine2(1, i2) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine2(2, i2) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3))

        else if((mod2(1).eq.0).and.(mod2(2).ne.0).and.(mod2(3).ne.0)) then
          i4 = i4 + 1
          tt%to_fine4(1, i4) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine4(3, i4) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3))
          tt%to_fine4(2, i4) = coarse%Lxyz_inv(x(1), x(2), x(3) + mod2(3))
          tt%to_fine4(4, i4) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3) + mod2(3))
        else if((mod2(3).eq.0).and.(mod2(1).ne.0).and.(mod2(2).ne.0)) then
          i4 = i4 + 1
          tt%to_fine4(1, i4) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine4(2, i4) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3))
          tt%to_fine4(3, i4) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3))
          tt%to_fine4(4, i4) = coarse%Lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3))
        else if((mod2(2).eq.0).and.(mod2(3).ne.0).and.(mod2(1).ne.0)) then
          i4 = i4 + 1
          tt%to_fine4(1, i4) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine4(2, i4) = coarse%Lxyz_inv(x(1), x(2), x(3) + mod2(3))
          tt%to_fine4(3, i4) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3))
          tt%to_fine4(4, i4) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3) + mod2(3))

        else if((mod2(1).ne.0).and.(mod2(2).ne.0).and.(mod2(3).ne.0)) then
          i8 = i8 + 1
          tt%to_fine8(1, i8) = coarse%Lxyz_inv(x(1), x(2), x(3))
          tt%to_fine8(2, i8) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3))
          tt%to_fine8(3, i8) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3))
          tt%to_fine8(4, i8) = coarse%Lxyz_inv(x(1), x(2), x(3) + mod2(3))
          tt%to_fine8(5, i8) = coarse%Lxyz_inv(x(1), x(2) + mod2(2), x(3) + mod2(3))
          tt%to_fine8(6, i8) = coarse%Lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3))
          tt%to_fine8(7, i8) = coarse%Lxyz_inv(x(1) + mod2(1), x(2), x(3) + mod2(3))
          tt%to_fine8(8, i8) = coarse%Lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3))
        end if

      end do

      ASSERT(i1==tt%n_fine1.and.i2==tt%n_fine2.and.i4==tt%n_fine4.and.i8==tt%n_fine8)

    end subroutine get_transfer_tables

  end subroutine multigrid_init


  ! ---------------------------------------------------------
  subroutine multigrid_end(mgrid)
    type(multigrid_type), intent(inout) :: mgrid

    integer :: i
    type(multigrid_level_type), pointer :: level

    do i = 1, mgrid%n_levels
      level => mgrid%level(i)

      call f_der_end(level%f_der)
      call mesh_end(level%m)
      deallocate(level%m, level%f_der)
      nullify   (level%m, level%f_der)

      deallocate(level%to_coarse, level%to_fine1, level%to_fine2, &
        level%to_fine4, level%to_fine8, level%fine_i)
      nullify   (level%to_coarse, level%to_fine1, level%to_fine2, &
        level%to_fine4, level%to_fine8, level%fine_i)
    end do

    deallocate(mgrid%level)

  end subroutine multigrid_end


  ! ---------------------------------------------------------
  subroutine multigrid_fine2coarse(mgrid, ilevel, f_fine, f_coarse, method_p)
    type(multigrid_type), intent(in)  :: mgrid
    integer,              intent(in)  :: ilevel
    FLOAT,                intent(in)  :: f_fine(:)
    FLOAT,                intent(out) :: f_coarse(:)
    integer, optional,    intent(in)  :: method_p
    integer :: method

    if(present(method_p)) then
      method=method_p
    else
      method=FULLWEIGHT
    end if

    select case(method)
    case(FULLWEIGHT)
      call multigrid_restriction(mgrid, ilevel, f_fine, f_coarse)
    case(INJECTION)
      call multigrid_injection(mgrid, ilevel, f_fine, f_coarse)
    case default
      write(message(1), '(a,i2,a)') 'Multigrid: Restriction method  = ', method, ' is not valid.'
      call write_fatal(1)
    end select

  end subroutine multigrid_fine2coarse


  ! ---------------------------------------------------------
  subroutine multigrid_injection(mgrid, ilevel, f_fine, f_coarse)
    type(multigrid_type), intent(in)  :: mgrid
    integer,              intent(in)  :: ilevel
    FLOAT,                intent(in)  :: f_fine(:)
    FLOAT,                intent(out) :: f_coarse(:)

    integer :: i
    type(multigrid_level_type), pointer :: level

    ASSERT(ilevel>0.and.ilevel<=mgrid%n_levels)

    level => mgrid%level(ilevel)
    do i = 1, level%n_coarse
      f_coarse(i) = f_fine(level%to_coarse(i))
    end do

  end subroutine multigrid_injection


  ! ---------------------------------------------------------
  subroutine multigrid_restriction(mgrid, ilevel, f_fine, f_coarse)
    type(multigrid_type), intent(in)  :: mgrid
    integer,              intent(in)  :: ilevel
    FLOAT,                intent(in)  :: f_fine(:)
    FLOAT,                intent(out) :: f_coarse(:)

    FLOAT :: weight(-1:1,-1:1,-1:1)

    integer :: n,fn,di,dj,dk,d,fi(3)
    type(multigrid_level_type), pointer :: level
    type(mesh_type), pointer :: fine_mesh,coarse_mesh


    ASSERT(ilevel>0.and.ilevel<=mgrid%n_levels)

    do di = -1, 1
      do dj = -1, 1
        do dk = -1, 1
          d = abs(di) + abs(dj) + abs(dk)
          weight(di, dj, dk) = CNST(0.5)**d
        end do
      end do
    end do

    level => mgrid%level(ilevel)

    fine_mesh => mgrid%level(ilevel-1)%m
    coarse_mesh => mgrid%level(ilevel)%m

    do n = 1, level%n_coarse
      fn   = level%to_coarse(n)
      fi(:)= fine_mesh%Lxyz(fn,:)

      f_coarse(n) = M_ZERO

      do di = -1, 1
        do dj = -1, 1
          do dk = -1, 1
            fn = fine_mesh%Lxyz_inv(fi(1)+di,fi(2)+dj,fi(3)+dk)
            if(fn <= fine_mesh%np ) then
              f_coarse(n) = f_coarse(n)+weight(di,dj,dk)*f_fine(fn)*fine_mesh%vol_pp(fn)
            end if
          end do
        end do
      end do

      f_coarse(n) = f_coarse(n)/coarse_mesh%vol_pp(n)
    end do

  end subroutine multigrid_restriction


  ! ---------------------------------------------------------
  subroutine multigrid_coarse2fine(mgrid, ilevel, f_coarse, f_fine)
    type(multigrid_type), intent(in)  :: mgrid
    integer,              intent(in)  :: ilevel
    FLOAT,                intent(in)  :: f_coarse(:)
    FLOAT,                intent(out) :: f_fine(:)

    type(multigrid_level_type), pointer :: level
    FLOAT, pointer :: vol_pp(:)

    integer :: i, i1, i2, i4, i8
    integer :: jj, j(8)
    FLOAT   :: vol_total

    ASSERT(ilevel>0.and.ilevel<=mgrid%n_levels)

    level  => mgrid%level(ilevel)
    vol_pp => level%m%vol_pp

    i1 = 0;  i2 = 0;  i4 = 0;  i8 = 0;
    do i = 1, level%n_fine
      select case(level%fine_i(i))
      case(1)
        i1 = i1 + 1
        j(1:1) = level%to_fine1(1:1, i1)
      case(2)
        i2 = i2 + 1
        j(1:2) = level%to_fine2(1:2, i2)
      case(4)
        i4 = i4 + 1
        j(1:4) = level%to_fine4(1:4, i4)
      case(8)
        i8 = i8 + 1
        j(1:8) = level%to_fine8(1:8, i8)
      end select

      f_fine(i) = M_ZERO
      vol_total = M_ZERO
      do jj = 1, level%fine_i(i)
        f_fine(i) = f_fine(i) + vol_pp(j(jj))*f_coarse(j(jj))
        vol_total = vol_total + vol_pp(j(jj))
      end do
      f_fine(i) = f_fine(i)/vol_total

    end do
  end subroutine multigrid_coarse2fine

end module multigrid
