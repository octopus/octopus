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
  use io_function_oct_m
  use cube_function_oct_m
  use geometry_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
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
    private
    integer, public         :: abtype
    FLOAT, pointer, public  :: mf(:)     !< The mask-function on the mesh
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
  subroutine bc_init(this, namespace, mesh, sb, geo)
    type(bc_t),               intent(out) :: this
    type(namespace_t),        intent(in)  :: namespace
    type(mesh_t),             intent(in)  :: mesh
    type(simul_box_t),        intent(in)  :: sb
    type(geometry_t),         intent(in)  :: geo

    integer             :: ip
    FLOAT               :: bounds(1:MAX_DIM,1:2)
    integer             :: cols_abshape_block, imdim, maxdim

    FLOAT               :: xx(1:MAX_DIM), rr
    FLOAT               :: ufn_re, ufn_im
    character(len=1024) :: user_def_expr

    FLOAT, allocatable  :: mf(:)
    FLOAT               :: abheight, abwidth, abwidth_def
    type(block_t)       :: blk

    character(len=50)   :: str

    PUSH_SUB(bc_init)
    
    bounds = M_ZERO

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
    call parse_variable(namespace, 'AbsorbingBoundaries', NOT_ABSORBING, this%abtype)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', this%abtype, is_flag = .true.)) then
      call messages_input_error(namespace, 'AbsorbingBoundaries')
    end if

    if(this%abtype == EXTERIOR) &
      call messages_not_implemented('Exterior complex scaling', namespace=namespace)

    if(this%abtype /= NOT_ABSORBING) then
      write(str, '(a,i5)') 'Absorbing Boundaries'
      call messages_print_stress(stdout, trim(str), namespace=namespace)

      !%Variable ABCapHeight
      !%Type float
      !%Default -0.2 a.u.
      !%Section Time-Dependent::Absorbing Boundaries
      !%Description
      !% When <tt>AbsorbingBoundaries = cap</tt>, this is the height of the imaginary potential.
      !%End
      if(this%abtype == IMAGINARY_ABSORBING) then
        call parse_variable(namespace, 'ABCapHeight', -CNST(0.2), abheight, units_inp%energy)
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
      if(parse_block(namespace, 'ABShape', blk) < 0) then
        message(1) = "Input: ABShape not specified. Using default values for absorbing boundaries."
        call messages_info(1)
      
        select case (sb%box_shape)
        case (SPHERE)
          bounds(1,1) = sb%rsize/M_TWO
          bounds(1,2) = sb%rsize
        case (PARALLELEPIPED)
          bounds(1:sb%dim,1)=sb%lsize(1:sb%dim)/M_TWO
          bounds(1:sb%dim,2)=sb%lsize(1:sb%dim)
        case (CYLINDER)
          bounds(1,2) = sb%xsize
          bounds(2,2) = sb%rsize
          bounds(1,1) = bounds(1,1)/M_TWO
          bounds(2,1) = bounds(2,2)/M_TWO
        case default
          bounds(:,1) = M_ZERO
          bounds(:,2) = CNST(0.4)
        end select
      else
        cols_abshape_block = parse_block_cols(blk, 0)
      
        select case(cols_abshape_block)
        case(2)
          call parse_block_float(blk, 0, 0, bounds(1,1), units_inp%length)
          call parse_block_float(blk, 0, 1, bounds(1,2), units_inp%length)
          if (sb%box_shape == SPHERE) then
            if(bounds(1,2) > sb%rsize)  bounds(1,2) = sb%rsize 
            message(1) = "Info: using spherical absorbing boundaries."
          else if (sb%box_shape == PARALLELEPIPED) then
            do imdim = 1, sb%dim
              if(bounds(imdim,2) > sb%lsize(imdim))  bounds(imdim,2) = sb%lsize(imdim) 
            end do
            message(1) = "Info: using cubic absorbing boundaries."
          else if (sb%box_shape == CYLINDER) then
            if(bounds(1,2) > sb%xsize)  bounds(1,2) = sb%xsize
            if(bounds(1,2) > sb%rsize)  bounds(1,2) = sb%rsize            
            message(1) = "Info: using cylindrical absorbing boundaries."
          end if    
          call messages_info(1)
        case(3)
          this%ab_user_def = .true.
          SAFE_ALLOCATE(this%ab_ufn(1:mesh%np))
          this%ab_ufn = M_ZERO
          call parse_block_float( blk, 0, 0, bounds(1,1), units_inp%length)
          call parse_block_float( blk, 0, 1, bounds(1,2), units_inp%length)
          call parse_block_string(blk, 0, 2, user_def_expr)
          do ip = 1, mesh%np
            xx = M_ZERO
            xx(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM)
            rr = units_from_atomic(units_inp%length, sqrt(sum(xx(1:sb%dim)**2)))
            do imdim = 1, sb%dim
              xx(imdim) = units_from_atomic(units_inp%length, xx(imdim))
            end do
            call parse_expression(ufn_re, ufn_im, sb%dim, xx, rr, M_ZERO, user_def_expr)
            this%ab_ufn(ip) = ufn_re
          end do
          message(1) = "Input: using user-defined function from expression:"
          write(message(2),'(a,a)') '   F(x,y,z) = ', trim(user_def_expr) 
          call messages_info(2)
        case default
          message(1) = "Input: ABShape block must have at least 2 columns."
          call messages_fatal(1, namespace=namespace)
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
      abwidth_def = bounds(1,2)-bounds(1,1)
      call parse_variable(namespace, 'ABWidth', abwidth_def, abwidth, units_inp%length)
      bounds(1:sb%dim,1) = bounds(1:sb%dim,2) - abwidth

      maxdim = sb%dim
      if(sb%box_shape == SPHERE .or. this%ab_user_def) maxdim = 1
      if(sb%box_shape == CYLINDER) maxdim = 2
      write(message(1),'(a,2a,9f12.6)') & 
          "  Lower bound [", trim(units_abbrev(units_inp%length)), '] =', & 
          (units_from_atomic(units_inp%length, bounds(imdim,1)), imdim=1,maxdim)
       write(message(2),'(a,2a,9f12.6)') & 
          "  Upper bound [", trim(units_abbrev(units_inp%length)), '] =', & 
          (units_from_atomic(units_inp%length, bounds(imdim,2)), imdim=1,maxdim)
      call messages_info(2)
      
      ! generate boundary function
      SAFE_ALLOCATE(mf(1:mesh%np))
      call bc_generate_mf(this, mesh, sb, geo, bounds, mf)
      
      ! mask or cap
      SAFE_ALLOCATE(this%mf(1:mesh%np))

      if(this%abtype == MASK_ABSORBING) then
        this%mf = M_ONE - mf
      else if(this%abtype == IMAGINARY_ABSORBING) then
        this%mf(:) = abheight * mf(:)
      end if
      
      if(debug%info) call bc_write_info(this, mesh, namespace)
      
      call messages_print_stress(stdout)
      
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
  subroutine bc_write_info(this, mesh, namespace)
    type(bc_t),               intent(in) :: this
    type(mesh_t),             intent(in) :: mesh
    type(namespace_t),        intent(in) :: namespace

    integer :: err

    PUSH_SUB(bc_write_info)

    if(this%abtype == MASK_ABSORBING) then
      call dio_function_output(io_function_fill_how("VTK"), "./td.general", "mask", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", "mask", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
    else if(this%abtype == IMAGINARY_ABSORBING) then
      call dio_function_output(io_function_fill_how("VTK"), "./td.general", "cap", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", "cap", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneX"), "./td.general", "cap", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneY"), "./td.general", "cap", namespace, mesh, &
        this%mf(1:mesh%np), unit_one, err)
    end if

    POP_SUB(bc_write_info)
  end subroutine bc_write_info

  ! ---------------------------------------------------------
  subroutine bc_generate_mf(this, mesh, sb, geo, bounds, mf)
    type(bc_t),               intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    type(simul_box_t),        intent(in)    :: sb
    type(geometry_t),         intent(in)    :: geo
    FLOAT,                    intent(in)    :: bounds(1:MAX_DIM, 1:2)
    FLOAT,                    intent(inout) :: mf(:)

    integer :: ip, dir
    FLOAT   :: width(MAX_DIM)
    FLOAT   :: xx(1:MAX_DIM), rr, dd, ddv(1:MAX_DIM), tmp(1:MAX_DIM)

    PUSH_SUB(bc_generate_mf)

    ! generate the boundaries on the mesh 

    mf = M_ZERO

    width(1:sb%dim) = bounds(1:sb%dim,2) - bounds(1:sb%dim,1)
    xx = M_ZERO

    do ip = 1, mesh%np
      xx(1:sb%dim) = mesh%x(ip, 1:sb%dim)
      rr = sqrt(dot_product(xx(1:sb%dim), xx(1:sb%dim)))
 
      if(this%ab_user_def) then
        dd = this%ab_ufn(ip) - bounds(1,1)
        if(dd > M_ZERO) then
          if(this%ab_ufn(ip) < bounds(1,2) ) then
            mf(ip) = sin(dd * M_PI / (M_TWO * (width(1)) ))**2
          else
            mf(ip) = M_ONE
          end if
        end if
 
      else ! this%ab_user_def == .false.
 
        select case (sb%box_shape)
        case (SPHERE)
        
          dd = rr -  bounds(1,1) 
          if(dd > M_ZERO ) then 
            if (dd  <  width(1)) then
              mf(ip) = sin(dd * M_PI / (M_TWO * (width(1)) ))**2
            else 
              mf(ip) = M_ONE
            end if
          end if
 
        case (PARALLELEPIPED)

          ! We are filling from the center opposite to the spherical case
          tmp = M_ONE
          mf(ip) = M_ONE
          ddv(:) = abs(xx(:)) -  bounds(:,1)
          do dir=1, sb%dim
            if(ddv(dir) > M_ZERO ) then 
              if (ddv(dir)  <  width(dir)) then
                tmp(dir) = M_ONE - sin(ddv(dir) * M_PI / (M_TWO * (width(dir)) ))**2
              else 
                tmp(dir) = M_ZERO
              end if
            end if        
          mf(ip) = mf(ip) * tmp(dir)
          end do
          mf(ip) = M_ONE - mf(ip)
          
        case (CYLINDER)
          
          rr = sqrt(dot_product(xx(2:sb%dim), xx(2:sb%dim)))
          tmp = M_ONE
          mf(ip) = M_ONE
          ddv(1) = abs(xx(1)) - bounds(1,1) 
          ddv(2) = rr         - bounds(2,1) 
          do dir=1, 2
            if(ddv(dir) > M_ZERO ) then 
              if (ddv(dir)  <  width(dir)) then
                tmp(dir) = M_ONE - sin(ddv(dir) * M_PI / (M_TWO * (width(dir)) ))**2
              else 
                tmp(dir) = M_ZERO
              end if
            end if        
          mf(ip) = mf(ip) * tmp(dir)
          end do
          mf(ip) = M_ONE - mf(ip)

        case default

          if (mesh_inborder(mesh, geo, ip, dd, width(1))) then
            mf(ip) = M_ONE - sin(dd * M_PI / (M_TWO * width(1)))**2
          end if

        end select
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
