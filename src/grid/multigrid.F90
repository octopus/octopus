!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

module multigrid_m
  use boundaries_m
  use curvilinear_m
  use derivatives_m
  use datasets_m
  use geometry_m
  use global_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_init_m
  use messages_m
  use par_vec_m
  use stencil_m
  use transfer_table_m
  use profiling_m

  implicit none

  private
  public ::                         &
    multigrid_level_t,              &
    multigrid_t,                    &
    multigrid_init,                 &
    multigrid_end,                  &
    multigrid_mesh_half,            &
    multigrid_mesh_double,          &
    multigrid_number_of_levels,     &
    dmultigrid_fine2coarse,         &
    zmultigrid_fine2coarse,         &
    dmultigrid_coarse2fine,         &
    zmultigrid_coarse2fine,         &
    multigrid_get_transfer_tables

  integer, parameter, public :: &
    INJECTION  = 1,             &
    FULLWEIGHT = 2

  type multigrid_level_t
    type(transfer_table_t)          :: tt
    type(mesh_t),          pointer  :: mesh
    type(derivatives_t),   pointer  :: der
  end type multigrid_level_t

  type multigrid_t
    integer                          :: n_levels
    type(multigrid_level_t), pointer :: level(:)

    integer          :: tp
    integer, pointer :: sp(:)
    integer, pointer :: ep(:)
    integer, pointer :: ep_part(:)
  end type multigrid_t

  type(profile_t), save :: interp_prof, injection_prof, restrict_prof

contains

  ! ---------------------------------------------------------
  subroutine multigrid_init(mgrid, geo, cv, mesh, der, stencil)
    type(multigrid_t),     target, intent(out) :: mgrid
    type(geometry_t),              intent(in)  :: geo
    type(curvilinear_t),           intent(in)  :: cv
    type(mesh_t),          target, intent(in)  :: mesh
    type(derivatives_t),   target, intent(in)  :: der
    type(stencil_t),               intent(in)  :: stencil

    integer :: i, n_levels, np

    PUSH_SUB(multigrid_init)

    !%Variable MultigridLevels
    !%Type integer
    !%Default max_levels
    !%Section Mesh
    !%Description
    !% Number of levels in the grid hierarchy used for multigrid. Positive
    !% numbers indicate an absolute number of levels, negative
    !% numbers are subtracted from the maximum number of levels possible.
    !%Option max_levels 0
    !% Calculate the optimal number of levels for the grid.
    !%End

    call parse_integer(datasets_check('MultigridLevels'), 0, n_levels)

    if ( n_levels <= 0 )then
      n_levels=n_levels-3
      np=mesh%np

      do while ( np > 1 )
        np=np/8
        n_levels=n_levels+1
      end do
    else
      n_levels=n_levels-1
    end if

    mgrid%n_levels = n_levels

    SAFE_ALLOCATE(mgrid%level(0:n_levels))

    mgrid%level(0)%mesh => mesh
    mgrid%level(0)%der => der

    mgrid%level(0)%tt%n_fine = mesh%np
    SAFE_ALLOCATE(mgrid%level(0)%tt%fine_i(1:mesh%np))

    write(message(1), '(a,i3)') "Multigrid levels:", n_levels + 1
    call messages_info(1)

    do i = 1, mgrid%n_levels
      SAFE_ALLOCATE(mgrid%level(i)%mesh)
      SAFE_ALLOCATE(mgrid%level(i)%der)
      
      call multigrid_mesh_half(geo, cv, mgrid%level(i-1)%mesh, mgrid%level(i)%mesh, stencil)

      call derivatives_init(mgrid%level(i)%der, mesh%sb, cv%method.ne.CURV_METHOD_UNIFORM)

      if(mesh%parallel_in_domains) then
        call mesh_init_stage_3(mgrid%level(i)%mesh, stencil, mesh%mpi_grp, parent = mgrid%level(i - 1)%mesh)
      else
        call mesh_init_stage_3(mgrid%level(i)%mesh)
      end if

      call multigrid_get_transfer_tables(mgrid%level(i)%tt, mgrid%level(i-1)%mesh, mgrid%level(i)%mesh)

      call derivatives_build(mgrid%level(i)%der, mgrid%level(i)%mesh)

      call mesh_write_info(mgrid%level(i)%mesh, stdout)
      
      mgrid%level(i)%der%finer => mgrid%level(i - 1)%der
      mgrid%level(i - 1)%der%coarser => mgrid%level(i)%der
      mgrid%level(i)%der%to_finer => mgrid%level(i)%tt
      mgrid%level(i - 1)%der%to_coarser => mgrid%level(i)%tt
    end do
    
    SAFE_ALLOCATE(mgrid%sp(0:mgrid%n_levels))
    SAFE_ALLOCATE(mgrid%ep(0:mgrid%n_levels))
    SAFE_ALLOCATE(mgrid%ep_part(0:mgrid%n_levels))

    mgrid%tp = 0
    do i = 0, mgrid%n_levels    
      mgrid%sp(i) = 1 + mgrid%tp
      mgrid%ep(i) = mgrid%tp + mgrid%level(i)%mesh%np
      mgrid%tp = mgrid%tp + mgrid%level(i)%mesh%np_part
      mgrid%ep_part(i) = mgrid%tp      
    end do

    POP_SUB(multigrid_init)
  end subroutine multigrid_init

  ! ---------------------------------------------------------
  !> creates the lookup tables to go between the coarse and fine meshes
  subroutine multigrid_get_transfer_tables(tt, fine, coarse)
    type(transfer_table_t), intent(inout) :: tt
    type(mesh_t),           intent(in)    :: fine, coarse

    integer :: i, i1, i2, i4, i8, pt, ig
#ifdef HAVE_MPI
    integer :: ii, jj
#endif
    integer :: x(MAX_DIM), mod2(MAX_DIM)

    PUSH_SUB(multigrid_get_transfer_tables)

    tt%n_coarse = coarse%np
    SAFE_ALLOCATE(tt%to_coarse(1:tt%n_coarse))

    ! GENERATE THE TABLE TO MAP FROM THE FINE TO THE COARSE GRID
    do i = 1, tt%n_coarse
      ig = i
#ifdef HAVE_MPI
      ! translate to a global index of the coarse grid
      if(coarse%parallel_in_domains) ig = coarse%vp%local(ig - 1 + coarse%vp%xlocal(coarse%vp%partno))
#endif
      ! locate the equivalent global fine grid point
      ig = fine%idx%lxyz_inv(2*coarse%idx%lxyz(ig, 1), 2*coarse%idx%lxyz(ig, 2), 2*coarse%idx%lxyz(ig, 3))
#ifdef HAVE_MPI
      ! translate to a local number of the fine grid
      if(fine%parallel_in_domains) ig = vec_global2local(fine%vp, ig, fine%vp%partno)
#endif
      tt%to_coarse(i) = ig
    end do

    ! count
    tt%n_fine = fine%np
    SAFE_ALLOCATE(tt%fine_i(1:tt%n_fine))

    tt%n_fine1 = 0
    tt%n_fine2 = 0
    tt%n_fine4 = 0
    tt%n_fine8 = 0
    do i = 1, tt%n_fine
      ig = i
#ifdef HAVE_MPI
      ! translate to a global index
      if(fine%parallel_in_domains) ig = fine%vp%local(ig - 1 + fine%vp%xlocal(fine%vp%partno))
#endif
      mod2 = mod(fine%idx%lxyz(ig, :), 2)
      
      pt = sum(abs(mod2(1:3)))
      
      select case(pt)
      case(0)
        tt%n_fine1 = tt%n_fine1 + 1
        tt%fine_i(i) = 1
      case(1)
        tt%n_fine2 = tt%n_fine2 + 1
        tt%fine_i(i) = 2
      case(2)
        tt%n_fine4 = tt%n_fine4 + 1
        tt%fine_i(i) = 4
      case(3)
        tt%n_fine8 = tt%n_fine8 + 1
        tt%fine_i(i) = 8
      end select
    end do
    
    ASSERT(tt%n_fine1 + tt%n_fine2 + tt%n_fine4 + tt%n_fine8 == tt%n_fine)

    SAFE_ALLOCATE(tt%to_fine1(1:1, 1:tt%n_fine1))
    SAFE_ALLOCATE(tt%to_fine2(1:2, 1:tt%n_fine2))
    SAFE_ALLOCATE(tt%to_fine4(1:4, 1:tt%n_fine4))
    SAFE_ALLOCATE(tt%to_fine8(1:8, 1:tt%n_fine8))

    ! and now build the tables
    i1 = 0;  i2 = 0;  i4 = 0;  i8 = 0
    do i = 1, fine%np
      ig = i
#ifdef HAVE_MPI
      ! translate to a global index
      if(fine%parallel_in_domains) ig = fine%vp%local(ig - 1 + fine%vp%xlocal(fine%vp%partno))
#endif
      x(1:3)    = fine%idx%lxyz(ig, 1:3)/2
      mod2(1:3) = mod(fine%idx%lxyz(ig, 1:3), 2)

      pt = sum(abs(mod2(1:3)))

      select case(pt)
      case(0)
        i1 = i1 + 1
        tt%to_fine1(1, i1) = coarse%idx%lxyz_inv(x(1), x(2), x(3))
        
      case(1)
        i2 = i2 + 1
        tt%to_fine2(1, i2) = coarse%idx%lxyz_inv(x(1)          , x(2)          , x(3)          )
        tt%to_fine2(2, i2) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3))
        
      case(2)
        i4 = i4 + 1
        tt%to_fine4(1, i4) = coarse%idx%lxyz_inv(x(1)          , x(2) + mod2(2), x(3) + mod2(3))
        tt%to_fine4(2, i4) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2)          , x(3) + mod2(3))
        tt%to_fine4(3, i4) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3)          )
        tt%to_fine4(4, i4) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3))
        
      case(3)
        i8 = i8 + 1
        tt%to_fine8(1, i8) = coarse%idx%lxyz_inv(x(1)          , x(2)          , x(3)          )
        tt%to_fine8(2, i8) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2)          , x(3)          )
        tt%to_fine8(3, i8) = coarse%idx%lxyz_inv(x(1)          , x(2) + mod2(2), x(3)          )
        tt%to_fine8(4, i8) = coarse%idx%lxyz_inv(x(1)          , x(2)          , x(3) + mod2(3))
        tt%to_fine8(5, i8) = coarse%idx%lxyz_inv(x(1)          , x(2) + mod2(2), x(3) + mod2(3))
        tt%to_fine8(6, i8) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2)          , x(3) + mod2(3))
        tt%to_fine8(7, i8) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3)          )
        tt%to_fine8(8, i8) = coarse%idx%lxyz_inv(x(1) + mod2(1), x(2) + mod2(2), x(3) + mod2(3))
        
      end select
      
    end do

    ASSERT(i1 == tt%n_fine1 .and. i2 == tt%n_fine2 .and. i4 == tt%n_fine4 .and. i8 == tt%n_fine8)

    ! translate to local points.
#ifdef HAVE_MPI
    if (coarse%parallel_in_domains) then

      do ii = 1, tt%n_fine1
        tt%to_fine1(1, ii) = vec_global2local(coarse%vp, tt%to_fine1(1, ii), coarse%vp%partno)
      end do

      do ii = 1, tt%n_fine2
        do jj = 1, 2
          tt%to_fine2(jj, ii) = vec_global2local(coarse%vp, tt%to_fine2(jj, ii), coarse%vp%partno)
        end do
      end do

      do ii = 1, tt%n_fine4
        do jj = 1, 4
          tt%to_fine4(jj, ii) = vec_global2local(coarse%vp, tt%to_fine4(jj, ii), coarse%vp%partno)
        end do
      end do

      do ii = 1, tt%n_fine8
        do jj = 1, 8
          tt%to_fine8(jj, ii) = vec_global2local(coarse%vp, tt%to_fine8(jj, ii), coarse%vp%partno)
        end do
      end do

    end if
#endif

    POP_SUB(multigrid_get_transfer_tables)
  end subroutine multigrid_get_transfer_tables

  !---------------------------------------------------------------------------------
  !> Creates a mesh that has twice the spacing betwen the points than the in mesh.
  !! This is used in the multi-grid routines
  !---------------------------------------------------------------------------------
  subroutine multigrid_mesh_half(geo, cv, mesh_in, mesh_out, stencil)
    type(geometry_t),           intent(in)    :: geo
    type(curvilinear_t),        intent(in)    :: cv
    type(mesh_t),       target, intent(in)    :: mesh_in
    type(mesh_t),               intent(inout) :: mesh_out
    type(stencil_t),            intent(in)    :: stencil

    PUSH_SUB(multigrid_mesh_half)

    mesh_out%sb              => mesh_in%sb
    mesh_out%idx%sb          => mesh_in%idx%sb
    mesh_out%use_curvilinear =  mesh_in%use_curvilinear
    mesh_out%cv              => mesh_in%cv

    mesh_out%spacing(:)  = 2*mesh_in%spacing(:)
    mesh_out%idx%nr(:,:) = mesh_in%idx%nr(:,:)/2
    mesh_out%idx%ll(:)   = mesh_out%idx%nr(2, :) - mesh_out%idx%nr(1, :) + 1

    mesh_out%idx%enlarge = mesh_in%idx%enlarge
    
    call mesh_init_stage_2(mesh_out, mesh_out%sb, geo, cv, stencil)

    POP_SUB(multigrid_mesh_half)
  end subroutine multigrid_mesh_half

  !---------------------------------------------------------------------------------
  subroutine multigrid_mesh_double(geo, cv, mesh_in, mesh_out, stencil)    
    type(geometry_t),           intent(in)    :: geo
    type(curvilinear_t),        intent(in)    :: cv
    type(mesh_t),       target, intent(in)    :: mesh_in
    type(mesh_t),               intent(inout) :: mesh_out
    type(stencil_t),            intent(in)    :: stencil

    PUSH_SUB(multigrid_mesh_double)

    mesh_out%sb             => mesh_in%sb
    mesh_out%idx%sb         => mesh_in%idx%sb
    mesh_out%use_curvilinear =  mesh_in%use_curvilinear
    mesh_out%cv             => mesh_in%cv

    mesh_out%spacing(:)  = M_HALF*mesh_in%spacing(:)
    mesh_out%idx%nr(:,:) = mesh_in%idx%nr(:,:)*2
    mesh_out%idx%ll(:)   = mesh_out%idx%nr(2, :) - mesh_out%idx%nr(1, :) + 1
    
    mesh_out%idx%enlarge = mesh_in%idx%enlarge
    
    call mesh_init_stage_2(mesh_out, mesh_out%sb, geo, cv, stencil)

    POP_SUB(multigrid_mesh_double)
  end subroutine multigrid_mesh_double

  ! ---------------------------------------------------------
  subroutine multigrid_end(mgrid)
    type(multigrid_t), target, intent(inout) :: mgrid

    integer :: i
    type(multigrid_level_t), pointer :: level

    PUSH_SUB(multigrid_end)

    SAFE_DEALLOCATE_P(mgrid%sp)
    SAFE_DEALLOCATE_P(mgrid%ep)
    SAFE_DEALLOCATE_P(mgrid%ep_part)

    SAFE_DEALLOCATE_P(mgrid%level(0)%tt%fine_i)

    do i = 1, mgrid%n_levels
      level => mgrid%level(i)

      call derivatives_end(level%der)
      call mesh_end(level%mesh)
      SAFE_DEALLOCATE_P(level%mesh)
      SAFE_DEALLOCATE_P(level%der)

      SAFE_DEALLOCATE_P(level%tt%to_coarse)
      SAFE_DEALLOCATE_P(level%tt%to_fine1)
      SAFE_DEALLOCATE_P(level%tt%to_fine2)
      SAFE_DEALLOCATE_P(level%tt%to_fine4)
      SAFE_DEALLOCATE_P(level%tt%to_fine8)
      SAFE_DEALLOCATE_P(level%tt%fine_i)
    end do

    SAFE_DEALLOCATE_P(mgrid%level)

    POP_SUB(multigrid_end)
  end subroutine multigrid_end

  !---------------------------------------------------------------------------------
  integer function multigrid_number_of_levels(base_der) result(number)
    type(derivatives_t), target, intent(in)  :: base_der

    type(derivatives_t), pointer :: next_der

    next_der => base_der%coarser

    number = 0
    do 
      number = number + 1
      next_der => next_der%coarser
      if(.not. associated(next_der)) exit
    end do

  end function multigrid_number_of_levels

#include "undef.F90"
#include "real.F90"
#include "multigrid_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "multigrid_inc.F90"

end module multigrid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
