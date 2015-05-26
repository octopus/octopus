!! Copyright (C) 2015 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.

#include "global.h"

!> This module should implement absorbing boundaries under the form of 
!! mask-function, suitable only for the TD Schroedinger equation, and
!! complex absorbing potential (CAPs), suitable also for the static case

module boundaries_abs_m
  use io_function_m
  use io_m
  use cube_function_m
  use global_m
  use geometry_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use simul_box_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::               &
    ab_t,                 &  
    ab_mask_t,            &
    ab_cap_t,             &
    ab_init,              &
    ab_end,               &
    ab_write_info

  type ab_mask_t
    CMPLX, pointer          :: mf(:)     !< The mask-function on the mesh
    type(cube_function_t)   :: cf        !< The mask-function on the cube
  end type ab_mask_t

  type ab_cap_t
    CMPLX, pointer          :: mf(:)     !< The CAP on the mesh
    type(cube_function_t)   :: cf        !< The CAP on the cube
    FLOAT                   :: height
  end type ab_cap_t
  
  type ab_t
    logical         :: use_mask
    type(ab_mask_t) :: mask

    logical         :: use_cap
    type(ab_cap_t)  :: cap

    logical         :: ab_user_def
    FLOAT, pointer  :: ab_ufn(:)
  end type ab_t

  integer, parameter :: &
    AB_MASK = 0, &
    AB_CAP  = 2

contains

  ! ---------------------------------------------------------
  subroutine ab_init(this, mesh, sb, geo)
    type(ab_t),               intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(geometry_t),         intent(in)  :: geo

    integer             :: ip
    FLOAT               :: bounds(1:2)
    integer             :: cols_abshape_block, idim

    FLOAT               :: xx(1:mesh%sb%dim), rr
    FLOAT               :: ufn_re, ufn_im
    character(len=1024) :: user_def_expr

    FLOAT, allocatable  :: mf(:)
    integer             :: ab_flags
    type(block_t)       :: blk

    PUSH_SUB(ab_init)

    !%Variable ABType
    !%Type flag
    !%Default ab_mask
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !%Sets the type of the absorbing boundaries.
    !%Option ab_mask 0
    !% Absorbing boundaries with mask function.
    !%Option ab_cap 2
    !% Absorbing boundaries with complex absorbing potential.
    !%End
    call parse_variable('ABType', AB_MASK, ab_flags)
    if(.not.varinfo_valid_option('ABType', ab_flags, is_flag = .true.)) then
      call messages_input_error('ABType')
    end if

    !%Variable ABShape
    !%Type block
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Set the shape of the absorbing boundaries. Here you can set the inner 
    !% and outer bounds by setting the block as follows:
    !%
    !% <tt>%ABShape
    !% <br>&nbsp;&nbsp; inner | outer | "user-defined"
    !% <br>%</tt>
    !%
    !% The optional 3rd column is a user-defined expression for the absorbing 
    !% boundaries. For example, <math>r</math> creates a spherical absorbing zone for 
    !% coordinates with <math>{\tt inner} < r < {\tt outer}</math>, and <math>z</math> creates an absorbing plane. 
    !% Note, values <tt>outer</tt> larger than the box size may lead in these cases to 
    !% unexpected reflection behaviours.
    !% If no expression is given, the absorbing zone follows the edges of the 
    !% box (not valid for user-defined box).
    !%End
    call messages_obsolete_variable('ABWidth', 'ABShape')

    cols_abshape_block = 0
    if(parse_block('ABShape', blk) < 0) then
      cols_abshape_block = parse_block_cols(blk, 0)
    else
      message(1) = "Input: ABShape not specified. Using default values for absorbing boundaries."
      call messages_info(1)

      if (mesh%sb%box_shape == SPHERE) then
        bounds(1)=mesh%sb%rsize/M_TWO
        bounds(2)=mesh%sb%rsize
      else if(mesh%sb%box_shape == PARALLELEPIPED) then
        bounds(1)=mesh%sb%lsize(1)/M_TWO
        bounds(2)=mesh%sb%lsize(1)
      else
        bounds(1)=M_ZERO
        bounds(2)=CNST(0.4)
      end if
    end if
 
    select case(cols_abshape_block)
    case(0)
      message(1) = "Input: ABShape block must have at least 2 columns."
      call messages_fatal(1)
    case(1)
      message(1) = "Input: ABShape block must have at least 2 columns."
      call messages_fatal(1)
    case(2)
      call parse_block_float(blk, 0, 0, bounds(1))
      call parse_block_float(blk, 0, 1, bounds(2))
 
      if (mesh%sb%box_shape == SPHERE) then
        if(bounds(2) > mesh%sb%rsize)  bounds(2) = mesh%sb%rsize 
        message(1) = "Info: using spherical absorbing boundaries."
      else if (mesh%sb%box_shape == PARALLELEPIPED) then
        if(bounds(2) > mesh%sb%lsize(1))  bounds(2) = mesh%sb%lsize(1) 
        message(1) = "Info: using cubic absorbing boundaries."
      end if    
      call messages_info(1)
 
    case(3)
      this%ab_user_def = .true.
      SAFE_ALLOCATE(this%ab_ufn(1:mesh%np))
      this%ab_ufn = M_ZERO
      call parse_block_float( blk, 0, 0, bounds(1))
      call parse_block_float( blk, 0, 1, bounds(2))
      call parse_block_string(blk, 0, 2, user_def_expr)
      do ip = 1, mesh%np
        xx = M_ZERO
        xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
        rr = units_from_atomic(units_inp%length, sqrt(sum(xx(1:mesh%sb%dim)**2)))
        forall(idim = 1:mesh%sb%dim) xx(idim) = units_from_atomic(units_inp%length, xx(idim))
        call parse_expression(ufn_re, ufn_im, mesh%sb%dim, xx, rr, M_ZERO, user_def_expr)
        this%ab_ufn(ip) = ufn_re
      end do
      message(1) = "Input: using user-defined function from expression:"
      write(message(2),'(a,a)') '   F(x,y,z) = ', trim(user_def_expr) 
      call messages_info(2)
    end select
 
    write(message(1),'(a,es10.3,3a)') & 
      "  val1 = ", units_from_atomic(units_inp%length, bounds(1) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    write(message(2),'(a,es10.3,3a)') & 
      "  val2 = ", units_from_atomic(units_inp%length, bounds(2) ),&
      ' [', trim(units_abbrev(units_inp%length)), ']'
    call messages_info(2)

    ! generate boundary function
    SAFE_ALLOCATE(mf(1:mesh%np))
    call ab_generate_mf(this, mesh, geo, bounds, mf)

    ! mask or cap
    if(this%use_mask) then
      SAFE_ALLOCATE(this%mask%mf(1:mesh%np))
      this%mask%mf = mf
    end if
  
    if(this%use_cap) then

      call messages_obsolete_variable('ABHeight', 'ABCapHeight')
      !%Variable ABCapHeight
      !%Type float
      !%Default -0.2 a.u.
      !%Section Time-Dependent::Absorbing Boundaries
      !%Description
      !% When <tt>AbsorbingBoundaries = sin2</tt>, this is the height of the imaginary potential.
      !%End
      call parse_variable('ABCapHeight', -CNST(0.2), this%cap%height, units_inp%energy)

      SAFE_ALLOCATE(this%cap%mf(1:mesh%np))   
      this%cap%mf(:) = this%cap%height * mf(:)

      call messages_not_implemented('Complex absorbing potential')
    end if

    SAFE_DEALLOCATE_A(mf)

    POP_SUB(ab_init)
  end subroutine ab_init

  ! ---------------------------------------------------------
  subroutine ab_generate_mf(this, mesh, geo, bounds, mf)
    type(ab_t),               intent(out)   :: this
    type(mesh_t),             intent(in)    :: mesh
    type(geometry_t),         intent(in)    :: geo
    FLOAT,                    intent(in)    :: bounds(1:2)
    FLOAT,                    intent(inout)   :: mf(:)

    integer :: ip, dir, ierr
    FLOAT   :: width
    FLOAT   :: xx(1:MAX_DIM), rr, dd, ddv(1:MAX_DIM), tmp(1:MAX_DIM)

    PUSH_SUB(ab_generate_mf)

    ! generate the boundaries on the mesh 

    width = bounds(2) - bounds(1)
    xx = M_ZERO
 
    do ip = 1, mesh%np
      xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
      rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
 
      if(this%ab_user_def) then
        dd = this%ab_ufn(ip) - bounds(1)
        if(dd > M_ZERO) then
          if(this%ab_ufn(ip) < bounds(2) ) then
            mf(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
          else
            mf(ip) = M_ONE
          end if
        end if
 
      else ! this%ab_user_def == .false.
 
        if(mesh%sb%box_shape == SPHERE) then
        
          dd = rr -  bounds(1) 
          if(dd > M_ZERO ) then 
            if (dd  <  width) then
              mf(ip) = M_ONE * sin(dd * M_PI / (M_TWO * (width) ))**2
            else 
              mf(ip) = M_ONE 
            end if
          end if
 
        else if (mesh%sb%box_shape == PARALLELEPIPED) then

          ! We are filling from the center opposite to the spherical case
          tmp = M_ONE
          mf(ip) = M_ONE
          ddv(:) = abs(xx(:)) -  bounds(1) 
          do dir=1, mesh%sb%dim
            if(ddv(dir) > M_ZERO ) then 
              if (ddv(dir)  <  width) then
                tmp(dir) = M_ONE - sin(ddv(dir) * M_PI / (M_TWO * (width) ))**2
              else 
                tmp(dir) = M_ZERO
              end if
            end if        
          mf(ip) = mf(ip)*tmp(dir)
          end do
          mf(ip) = M_ONE - mf(ip)

        else

          if(mesh_inborder(mesh, geo, ip, dd, width)) &
            mf(ip) = M_ONE - sin(dd * M_PI / (M_TWO * width))**2

        end if
      end if
    end do

    if(this%ab_user_def) then 
      SAFE_DEALLOCATE_P(this%ab_ufn)
    end if

    POP_SUB(ab_generate_mf)
  end subroutine ab_generate_mf

  ! ---------------------------------------------------------
  subroutine ab_end(this)
    type(ab_t),   intent(inout) :: this
    PUSH_SUB(ab_end)

    if(this%use_mask) call ab_mask_end(this%mask)
    if(this%use_cap)  call ab_cap_end(this%cap)

    POP_SUB(ab_end)
  end subroutine ab_end

  ! ---------------------------------------------------------
  subroutine ab_mask_end(this)
    type(ab_mask_t),  intent(inout) :: this
    PUSH_SUB(ab_mask_end)

    SAFE_DEALLOCATE_P(this%mf)

    POP_SUB(ab_mask_end)
  end subroutine ab_mask_end

  ! ---------------------------------------------------------
  subroutine ab_cap_end(this)
    type(ab_cap_t),  intent(inout) :: this
    PUSH_SUB(ab_cap_end)

    POP_SUB(ab_cap_end)
  end subroutine ab_cap_end

  ! ---------------------------------------------------------
  subroutine ab_write_info(this)
    type(ab_t),   intent(in) :: this
    PUSH_SUB(ab_write_info)


    POP_SUB(ab_write_info)
  end subroutine ab_write_info

end module boundaries_abs_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
