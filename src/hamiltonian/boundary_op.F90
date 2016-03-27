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

module boundary_op_oct_m
  use io_oct_m
  use io_function_oct_m
  use io_oct_m
  use cube_function_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use simul_box_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::               &
    bc_init,              &
    bc_end,               &
    bc_write_info,        &
    bc_t

  type bc_t
    integer                 :: abtype
    FLOAT, pointer          :: mf(:)     !< The mask-function on the mesh
    type(cube_function_t)   :: cf        !< The mask-function on the cube
    logical                 :: ab_user_def
    FLOAT, pointer          :: ab_ufn(:)
  end type bc_t    

  integer, public, parameter :: &
    NOT_ABSORBING       = 0, &
    MASK_ABSORBING      = 1, &
    IMAGINARY_ABSORBING = 2, &
    EXTERIOR            = 3

contains

  ! ---------------------------------------------------------
  subroutine bc_init(this, mesh, sb, geo)
    type(bc_t),               intent(out) :: this
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(geometry_t),         intent(in)  :: geo

    integer             :: ip
    FLOAT               :: bounds(1:2)
    integer             :: cols_abshape_block, imdim

    FLOAT               :: xx(1:MAX_DIM), rr
    FLOAT               :: ufn_re, ufn_im
    character(len=1024) :: user_def_expr

    FLOAT, allocatable  :: mf(:)
    FLOAT               :: abheight, abwidth
    type(block_t)       :: blk

    character(len=50)   :: str

    PUSH_SUB(bc_init)

    this%ab_user_def = .false.
    
    !%Variable AbsorbingBoundaries
    !%Type flag
    !%Default not_absorbing
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of
    !% standing density waves, one can make the boundaries of the simulation
    !% box absorbing and use exterior complex scaling.
    !%Option not_absorbing 0
    !% Reflecting boundaries.
    !%Option mask 1
    !% Absorbing boundaries with a mask function.
    !%Option cap 2
    !% Absorbing boundaries with a complex absorbing potential.
    !%Option exterior 3
    !% Exterior complex scaling (not yet implemented).
    !%End
    call parse_variable('AbsorbingBoundaries', NOT_ABSORBING, this%abtype)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', this%abtype, is_flag = .true.)) then
      call messages_input_error('AbsorbingBoundaries')
    end if

    if(this%abtype == EXTERIOR) &
      call messages_not_implemented('Exterior complex scaling')

    if(this%abtype /= NOT_ABSORBING) then
      write(str, '(a,i5)') 'Absorbing Boundaries'
      call messages_print_stress(stdout, trim(str))

      !%Variable ABCapHeight
      !%Type float
      !%Default -0.2 a.u.
      !%Section Time-Dependent::Absorbing Boundaries
      !%Description
      !% When <tt>AbsorbingBoundaries = cap</tt>, this is the height of the imaginary potential.
      !%End
      if(this%abtype == IMAGINARY_ABSORBING) then
        call parse_variable('ABCapHeight', -CNST(0.2), abheight, units_inp%energy)
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
      !% coordinates with <math>{\tt inner} < r < {\tt outer}</math>, and <math>z</math> creates an 
      !% absorbing plane. 
      !% Note, values <tt>outer</tt> larger than the box size may lead in these cases to 
      !% unexpected reflection behaviours.
      !% If no expression is given, the absorbing zone follows the edges of the 
      !% box (not valid for user-defined box).
      !%End

      cols_abshape_block = 0
      if(parse_block('ABShape', blk) < 0) then
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
      else
        cols_abshape_block = parse_block_cols(blk, 0)
      
        select case(cols_abshape_block)
        case(2)
          call parse_block_float(blk, 0, 0, bounds(1), units_inp%length)
          call parse_block_float(blk, 0, 1, bounds(2), units_inp%length)
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
          call parse_block_float( blk, 0, 0, bounds(1), units_inp%length)
          call parse_block_float( blk, 0, 1, bounds(2), units_inp%length)
          call parse_block_string(blk, 0, 2, user_def_expr)
          do ip = 1, mesh%np
            xx = M_ZERO
            xx(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM)
            rr = units_from_atomic(units_inp%length, sqrt(sum(xx(1:mesh%sb%dim)**2)))
            forall(imdim = 1:mesh%sb%dim) xx(imdim) = units_from_atomic(units_inp%length, xx(imdim))
            call parse_expression(ufn_re, ufn_im, mesh%sb%dim, xx, rr, M_ZERO, user_def_expr)
            this%ab_ufn(ip) = ufn_re
          end do
          message(1) = "Input: using user-defined function from expression:"
          write(message(2),'(a,a)') '   F(x,y,z) = ', trim(user_def_expr) 
          call messages_info(2)
        case default
          message(1) = "Input: ABShape block must have at least 2 columns."
          call messages_fatal(1)
        end select

        call parse_block_end(blk)
      end if
      
      !%Variable ABWidth
      !%Type float
      !%Section Time-Dependent::Absorbing Boundaries
      !%Description
      !% Specifies the boundary width. For a finer control over the absorbing boundary 
      !% shape use ABShape. 
      !%End
!       call messages_obsolete_variable('ABWidth', 'ABShape')
      abwidth = bounds(2)-bounds(1)
      call parse_variable('ABWidth', abwidth, abwidth, units_inp%length)
      bounds(1) = bounds(2) - abwidth
      
      write(message(1),'(a,es10.3,3a)') & 
        "  Lower bound = ", units_from_atomic(units_inp%length, bounds(1) ),&
        ' [', trim(units_abbrev(units_inp%length)), ']'
      write(message(2),'(a,es10.3,3a)') & 
        "  Upper bound = ", units_from_atomic(units_inp%length, bounds(2) ),&
        ' [', trim(units_abbrev(units_inp%length)), ']'
      call messages_info(2)
      
      ! generate boundary function
      SAFE_ALLOCATE(mf(1:mesh%np))
      call bc_generate_mf(this, mesh, geo, bounds, mf)
      
      ! mask or cap
      SAFE_ALLOCATE(this%mf(1:mesh%np))

      if(this%abtype == MASK_ABSORBING) then
        this%mf = M_ONE - mf
      else if(this%abtype == IMAGINARY_ABSORBING) then
        this%mf(:) = abheight * mf(:)
      end if
      
      if(debug%info) call bc_write_info(this, mesh)
      
      SAFE_DEALLOCATE_A(mf)
    end if

    POP_SUB(bc_init)
  end subroutine bc_init

  ! ---------------------------------------------------------
  subroutine bc_end(this)
    type(bc_t),   intent(inout) :: this
    PUSH_SUB(bc_end)

    if(this%abtype /= NOT_ABSORBING) then
      SAFE_DEALLOCATE_P(this%mf)
    end if

    POP_SUB(bc_end)
  end subroutine bc_end

  ! ---------------------------------------------------------
  subroutine bc_write_info(this, mesh)
    type(bc_t),               intent(in) :: this
    type(mesh_t),             intent(in) :: mesh

    integer :: err

    PUSH_SUB(bc_write_info)

    if(this%abtype == MASK_ABSORBING) then
      call dio_function_output(io_function_fill_how("VTK"), "./td.general", "mask", mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", "mask", mesh, &
        this%mf(1:mesh%np), unit_one, err)
    else if(this%abtype == IMAGINARY_ABSORBING) then
      call dio_function_output(io_function_fill_how("VTK"), "./td.general", "cap", mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", "cap", mesh, &
        this%mf(1:mesh%np), unit_one, err)
    end if

    POP_SUB(bc_write_info)
  end subroutine bc_write_info

  ! ---------------------------------------------------------
  subroutine bc_generate_mf(this, mesh, geo, bounds, mf)
    type(bc_t),               intent(out)   :: this
    type(mesh_t),             intent(in)    :: mesh
    type(geometry_t),         intent(in)    :: geo
    FLOAT,                    intent(in)    :: bounds(1:2)
    FLOAT,                    intent(inout) :: mf(:)

    integer :: ip, dir, ierr
    FLOAT   :: width
    FLOAT   :: xx(1:MAX_DIM), rr, dd, ddv(1:MAX_DIM), tmp(1:MAX_DIM)

    PUSH_SUB(bc_generate_mf)

    ! generate the boundaries on the mesh 

    mf = M_ZERO

    width = bounds(2) - bounds(1)
    xx = M_ZERO

    do ip = 1, mesh%np
      xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
      rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
 
      if(this%ab_user_def) then
        dd = this%ab_ufn(ip) - bounds(1)
        if(dd > M_ZERO) then
          if(this%ab_ufn(ip) < bounds(2) ) then
            mf(ip) = sin(dd * M_PI / (M_TWO * (width) ))**2
          else
            mf(ip) = M_ONE
          end if
        end if
 
      else ! this%ab_user_def == .false.
 
        if(mesh%sb%box_shape == SPHERE) then
        
          dd = rr -  bounds(1) 
          if(dd > M_ZERO ) then 
            if (dd  <  width) then
              mf(ip) = sin(dd * M_PI / (M_TWO * (width) ))**2
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
          mf(ip) = mf(ip) * tmp(dir)
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

    POP_SUB(bc_generate_mf)
  end subroutine bc_generate_mf

end module boundary_op_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
