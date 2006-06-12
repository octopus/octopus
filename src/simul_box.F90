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

module simul_box_m
  use global_m
  use messages_m
  use string_m
  use varinfo_m
  use units_m
  use datasets_m
  use lib_oct_parser_m
  use lib_oct_m
  use geometry_m
  use math_m

  implicit none

  private
  public ::                     &
    simul_box_t,                &
    simul_box_init,             &
    simul_box_write_info,       &
    simul_box_is_periodic,      &
    simul_box_in_box,           &
    simul_box_dump,             &
    simul_box_init_from_file,   &
    operator(.eq.)

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4,         &
    BOX_IMAGE      = 5,         &
    BOX_USDEF      = 123

  type simul_box_t
    integer  :: box_shape   ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                            ! 4->parallelpiped (orthonormal, up to now).

    FLOAT :: h(3)           ! the (canonical) spacing between the points
    FLOAT :: box_offset(3)  ! shifts of the origin in the respective direction

    FLOAT :: rsize          ! the radius of the sphere or of the cylinder
    FLOAT :: xsize          ! the length of the cylinder in the x direction
    FLOAT :: lsize(3)       ! half of the length of the parallelepiped in each direction.

    integer(POINTER_SIZE) :: image    ! for the box defined through an image
    character(len=1024)   :: user_def ! for the user defined box

    FLOAT :: rlat(3,3)      ! lattice primitive vectors
    FLOAT :: klat(3,3)      ! reciprocal lattice primitive vectors
    FLOAT :: shift(27,3)    ! shift to equivalent positions in nearest neighbour primitive cells

    FLOAT :: fft_alpha      ! enlargement factor for double box

    integer :: dim
    integer :: periodic_dim

  end type simul_box_t

  interface operator(.eq.)
    module procedure simul_box_is_eq
  end interface

contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, geo)
    type(simul_box_t), intent(inout) :: sb
    type(geometry_t),  intent(in)    :: geo

    ! some local stuff
    FLOAT :: def_h, def_rsize
    integer :: i, ix, iy, iz, ii(MAX_DIM)

    call push_sub('simul_box.simul_box_init')

    call geometry_grid_defaults(geo, def_h, def_rsize)

    call read_misc()          ! miscellany stuff
    call read_box()           ! parameters defining the simulation box
    call read_spacing ()      ! parameters defining the (canonical) spacing
    call read_box_offset()    ! parameters defining the offset of the origin
    call build_lattice()      ! build lattice vectors

    call pop_sub()

  contains

    !--------------------------------------------------------------
    subroutine read_misc()

      call push_sub('simul_box.read_misc')

      !%Variable DoubleFFTParameter
      !%Type float
      !%Default 2.0
      !%Section Mesh::FFTs
      !%Description
      !% For solving Poisson equation in Fourier space, and for applying the local potential
      !% in Fourier space, an auxiliary cubic mesh is built. This mesh will be larger than
      !% the circumscribed cube to the usual mesh by a factor <tt>DoubleFFTParameter</tt>. See
      !% the section that refers to Poisson equation, and to the local potential for details
      !% [The default value of two is typically good].
      !%End
      call loct_parse_float(check_inp('DoubleFFTParameter'), M_TWO, sb%fft_alpha)
      if (sb%fft_alpha < M_ONE .or. sb%fft_alpha > M_THREE ) then
        write(message(1), '(a,f12.5,a)') "Input: '", sb%fft_alpha, &
          "' is not a valid DoubleFFTParameter"
        message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
        call write_fatal(2)
      end if

      sb%dim = calc_dim

      !%Variable PeriodicDimensions
      !%Type integer
      !%Default 0
      !%Section Mesh::Simulation Box
      !%Description
      !% Define which directions are to be considered periodic. Of course, it has to be a number
      !% from 0 to three, and it cannot be larger than Dimensions.
      !%Option 0
      !% No direction is periodic (molecule)
      !%Option 1
      !% The x direction is periodic (wire)
      !%Option 2
      !% The x and y directions are periodic (slab)
      !%Option 3
      !% The x, y, and z directions are periodic (bulk)
      !%End
      call loct_parse_int(check_inp('PeriodicDimensions'), 0, sb%periodic_dim)
      if ((sb%periodic_dim < 0) .or. (sb%periodic_dim > 3) .or. (sb%periodic_dim > sb%dim)) &
        call input_error('PeriodicDimensions')

      call pop_sub()
    end subroutine read_misc


    !--------------------------------------------------------------
    subroutine read_box()
      integer(POINTER_SIZE) :: blk
      character(len=200) :: filename
      FLOAT :: default
      
      call push_sub('simul_box.read_box')
      ! Read box shape.
      ! need to find out calc_mode already here since some of the variables here (e.g.
      ! periodic dimensions) can be different for the subsystems

      !%Variable BoxShape
      !%Type integer
      !%Default minimum
      !%Section Mesh::Simulation Box
      !%Description
      !% This variable decides the shape of the simulation box.
      !% Note that some incompatibilities apply:
      !% <ul>
      !% <li>Spherical or minimum mesh is not allowed for periodic systems.</li>
      !% <li>Cylindrical mesh is not allowed for systems that are periodic in more than one dimension.</li>
      !% >li>Box_image is only allowed in 2D.</li>
      !% </ul>
      !%Option sphere 1
      !% The simulation box will be a sphere of radius Radius
      !%Option cylinder 2
      !% The simulation box will be a cylinder with radius Radius and height two times
      !% Xlength
      !%Option minimum 3
      !% The simulation box will be constructed by adding spheres created around each
      !% atom (or user defined potential), of radius Radius.
      !%Option parallelepiped 4
      !% The simulation box will be a parallelpiped whose dimensions are taken from
      !% the variable lsize.
      !%Option box_image 5
      !% The simulation box will be defined through an image. White means that the point
      !% is contained in the simulation box, while any other color means that the point is out.
      !%Option user_defined 123
      !% The shape of the simulation box will be read from the variable <tt>BoxShapeUsDef</tt>
      !%End
      call loct_parse_int(check_inp('BoxShape'), MINIMUM, sb%box_shape)
      if(.not.varinfo_valid_option('BoxShape', sb%box_shape)) call input_error('BoxShape')
      select case(sb%box_shape)
      case(SPHERE,MINIMUM,BOX_IMAGE,BOX_USDEF)
        if(sb%dim>1 .and. simul_box_is_periodic(sb)) call input_error('BoxShape')
      case(CYLINDER)
        if (sb%dim>2 .and. &
          ((sb%dim - sb%periodic_dim == 0) .or. (sb%dim - sb%periodic_dim == 1))) call input_error('BoxShape')
      end select

      ! ignore box_shape in 1D
      if(sb%dim==1.and.sb%box_shape /= PARALLELEPIPED) sb%box_shape=SPHERE

      ! Can not use images in 1D or 3D
      if(sb%dim/=2.and.sb%box_shape == BOX_IMAGE) call input_error('BoxShape')

      sb%rsize = -M_ONE
      !%Variable Radius
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% If BoxShape is not "parallelepiped" defines the radius of the spheres or of the cylinder.
      !% It has to be a positive number. If it is not defined in the input file, then the program
      !% will attempt to fine a suitable default, but this is not always possible, in which case
      !% the code will stop issuing this error message.
      !%End
      select case(sb%box_shape)
      case(SPHERE, CYLINDER)
        call loct_parse_float(check_inp('radius'), def_rsize/units_inp%length%factor, sb%rsize)
        if(sb%rsize < CNST(0.0)) call input_error('radius')
        sb%rsize = sb%rsize * units_inp%length%factor
        if(def_rsize>M_ZERO) call check_def(def_rsize, sb%rsize, 'radius')
      case(MINIMUM)
        default=sb%rsize
        call loct_parse_float(check_inp('radius'), default, sb%rsize)
        sb%rsize = sb%rsize * units_inp%length%factor
        if(sb%rsize < M_ZERO .and. def_rsize < M_ZERO) call input_error('Radius')
      end select

      if(sb%box_shape == CYLINDER) then
        !%Variable Xlength
        !%Type float
        !%Section Mesh::Simulation Box
        !%Description
        !% If BoxShape is "cylinder", it is half the total length of the cylinder.
        !%End
        call loct_parse_float(check_inp('xlength'), M_ONE/units_inp%length%factor, sb%xsize)
        sb%xsize = sb%xsize * units_inp%length%factor
        sb%lsize(1) = sb%xsize
        if(def_rsize>M_ZERO.and.sb%periodic_dim==0) call check_def(def_rsize, sb%xsize, 'xlength')
      end if

      sb%lsize = M_ZERO
      if(sb%box_shape == PARALLELEPIPED .or. sb%box_shape == BOX_IMAGE .or. sb%box_shape == BOX_USDEF) then

        !%Variable Lsize
        !%Type block
        !%Section Mesh::Simulation Box
        !%Description
        !% In case BoxShape is "parallelepiped", "box_image", or "user_defined", this is assumed to be a block of the form:
        !%
        !% <tt>%Lsize
        !% <br>&nbsp;&nbsp;sizex | sizey | sizez
        !% <br>%</tt>
        !%
        !% where the "size*" are half the lengths of the box in each direction.
        !%
        !% If BoxShape is "parallelpiped", this block has to be defined in the input file. The
        !% number of columns must match the dimensionality of the calculation.
        !%End
        if(loct_parse_block(check_inp('lsize'), blk) == 0) then
          if(loct_parse_block_cols(blk,0) < sb%dim) call input_error('lsize')
          do i = 1, sb%dim
            call loct_parse_block_float(blk, 0, i-1, sb%lsize(i))
          end do
        else
          call loct_parse_float(check_inp('lsize'), -M_ONE, sb%lsize(1))
          sb%lsize(1:sb%dim) = sb%lsize(1)
        end if
        sb%lsize = sb%lsize*units_inp%length%factor

        do i = 1, sb%dim
          if(def_rsize>M_ZERO.and.sb%periodic_dim<i) call check_def(def_rsize, sb%lsize(i), 'lsize')
        end do

        call loct_parse_block_end(blk)
      end if

      ! read in image for box_image
      if(sb%box_shape == BOX_IMAGE) then

        !%Variable BoxShapeImage
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Name of the file that contains the image that defines the simulation box.
        !%End
#if defined(HAVE_GDLIB)        
        call loct_parse_string(check_inp("BoxShapeImage"), "", filename)
        sb%image = loct_gdimage_create_from(filename)
        if(sb%image == 0) then
          message(1) = "Could not open file '" // filename // "'"
          call write_fatal(1)
        end if
#else
        message(1) = "To use 'BoxShape = box_image' you have to compile octopus"
        message(2) = "with GD library support"
        call write_fatal(2)
#endif
      end if     

      ! read in box shape for user defined boxes
      if(sb%box_shape == BOX_USDEF) then

        !%Variable BoxShapeUsDef
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Boolean expression that defines the interior of the simulation box. For example,
        !% <tt>BoxShapeUsDef = "(sqrt(x^2+y^2) <= 4) && z>-2 && z<2"</tt> defines a cylinder
        !% with axis parallel to the z axis
        !%End
        
        call loct_parse_string(check_inp("BoxShapeUsDef"), "x^2+y^2+z^2 < 4", sb%user_def)
        call conv_to_C_string(sb%user_def)
      end if

      ! fill in lsize structure
      select case(sb%box_shape)
      case(SPHERE)
        sb%lsize(1:sb%dim) = sb%rsize
      case(CYLINDER)
        sb%lsize(1)        = sb%xsize
        sb%lsize(2:sb%dim) = sb%rsize
      case(MINIMUM)
        do i = 1, sb%dim
          if(sb%rsize > M_ZERO) then
            sb%lsize(i) = maxval(abs(geo%atom(:)%x(i))) + sb%rsize
          else
            sb%lsize(i) = maxval(abs(geo%atom(:)%x(i))) + def_rsize
          end if
        end do
      end select

      call pop_sub()
    end subroutine read_box


    !--------------------------------------------------------------
    subroutine read_spacing()
      integer :: i, sx, sy
      integer(POINTER_SIZE) :: blk


      call push_sub('simul_box.read_spacing')

      ! initialize to -1
      sb%h = -M_ONE

#if defined(HAVE_GDLIB)
      if(sb%box_shape == BOX_IMAGE) then 
        ! spacing is determined from lsize and the size of the image
        sx = loct_gdImage_SX(sb%image)
        sy = loct_gdImage_SY(sb%image)

        sb%h(1) = M_TWO*sb%lsize(1)/real(sx, PRECISION)
        sb%h(2) = M_TWO*sb%lsize(2)/real(sy, PRECISION)
        return
      end if
#endif

      !%Variable Spacing
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% The spacing between the points in the mesh. In case of using curvilinear
      !% coordinates, this is a canonical spacing that will be changed locally by the
      !% transformation.
      !%
      !% It is possible to have a different spacing in each one of the cartesian directions
      !% if we define <tt>Spacing</tt> as block of the form
      !%
      !% <tt>%Spacing
      !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
      !% <br>%</tt>
      !%End

      if(loct_parse_block(check_inp('Spacing'), blk) == 0) then
        if(loct_parse_block_cols(blk,0) < sb%dim) call input_error('Spacing')
        do i = 1, sb%dim
          call loct_parse_block_float(blk, 0, i-1, sb%h(i))
        end do
        call loct_parse_block_end(blk)
      else
        call loct_parse_float(check_inp('Spacing'), sb%h(1), sb%h(1))
        sb%h(1:sb%dim) = sb%h(1)
      end if

      do i = 1, sb%dim
        sb%h(i) = sb%h(i)*units_inp%length%factor
        if(sb%h(i) < M_ZERO) then
          if(def_h > M_ZERO) then
            sb%h(i) = def_h
            write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", i, &
              ") [", trim(units_out%length%abbrev), "] = ",                 &
              sb%h(i)/units_out%length%factor
            call write_info(1)
          else
            message(1) = 'Either:'
            message(2) = "   *) variable 'Spacing' is not defined"
            message(3) = "   *) your input for 'Spacing' is negative"
            message(4) = "   *) I can't find a suitable default"
            call write_fatal(4)
          end if
        end if
        if(def_rsize>M_ZERO) call check_def(sb%h(i), def_rsize, 'Spacing')
      end do

      call pop_sub()
    end subroutine read_spacing


    !--------------------------------------------------------------
    subroutine read_box_offset()
      integer :: i
      integer(POINTER_SIZE) :: blk

      call push_sub('simul_box.read_box_offset')
      !%Variable BoxOffset
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% The zero of the simulation box. It can be either a float, or a block containing
      !% the (x,y,z) value of the zero.
      !%End
      sb%box_offset = M_ZERO
      if(loct_parse_block(check_inp('BoxOffset'), blk) == 0) then
        do i = 1, sb%dim
          call loct_parse_block_float(blk, 0, i-1, sb%box_offset(i))
        end do
        call loct_parse_block_end(blk)
      else
        call loct_parse_float(check_inp('BoxOffset'), M_ZERO, sb%box_offset(1))
        sb%box_offset(1:sb%dim) = sb%box_offset(1)
      end if
      sb%box_offset(:) = sb%box_offset(:)*units_inp%length%factor

      call pop_sub()
    end subroutine read_box_offset


    !--------------------------------------------------------------
    subroutine check_def(var, def, text)
      FLOAT, intent(in) :: var, def
      character(len=*), intent(in) :: text

      if(var > def) then
        write(message(1), '(3a)') "The value for '", text, "' does not match the recommended value"
        write(message(2), '(f8.3,a,f8.3)') var, ' > ', def
        call write_warning(2)
      end if
    end subroutine check_def


    !--------------------------------------------------------------
    subroutine build_lattice()

      call push_sub('simul_box.build_lattice')
      ! build primitive vectors (only simple cubic, tetra, or orthororhombic )
      sb%rlat = M_ZERO
      sb%klat = M_ZERO
      do i = 1, sb%dim
        sb%rlat(i,i) = 2*sb%lsize(i)
        sb%klat(i,i) = M_PI/sb%lsize(i)
      end do

      ! build shifts to nearest neighbour primitive cells
      ii = (/0,-1,1/)
      sb%shift=M_ZERO
      i = 1
      do iz = 1, 3
        do iy = 1, 3
          do ix = 1, 3
            sb%shift(i,1)=ii(ix)*sb%rlat(1,1)
            sb%shift(i,2)=ii(iy)*sb%rlat(2,2)
            sb%shift(i,3)=ii(iz)*sb%rlat(3,3)
            i = i + 1
          end do
        end do
      end do

      call pop_sub()
    end subroutine build_lattice

  end subroutine simul_box_init


  !--------------------------------------------------------------
  subroutine simul_box_end()

  end subroutine simul_box_end


  !--------------------------------------------------------------
  subroutine simul_box_write_info(sb, iunit)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: iunit

    character(len=15), parameter :: bs(5) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'around nuclei ', &
      'parallelepiped', &
      'image defined '/)

    call push_sub('simul_box.simul_box_write_info')

    write(message(1),'(a)') 'Simulation Box:'
    if(sb%box_shape.eq.BOX_USDEF) then
      write(message(2), '(a)') '  Type = user defined'
    else
      write(message(2), '(a,a,1x)') '  Type = ', bs(sb%box_shape)
    end if
    call write_info(2, iunit)

    if(sb%box_shape == SPHERE.or.sb%box_shape == CYLINDER.or.sb%box_shape == MINIMUM) then
      write(message(1), '(3a,f7.3)') '  Radius  [', trim(units_out%length%abbrev), '] = ', &
        sb%rsize/units_out%length%factor
      call write_info(1, iunit)
    end if
    if(sb%box_shape == CYLINDER) then
      write(message(1), '(3a,f7.3)') '  Xlength [', trim(units_out%length%abbrev), '] = ', &
        sb%xsize/units_out%length%factor
      call write_info(1, iunit)
    end if
    if(sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(3a, a, f8.3, a, f8.3, a, f8.3, a)')     &
        '  Lengths [', trim(units_out%length%abbrev), '] = ',    &
        '(', sb%lsize(1)/units_out%length%factor, ',',           &
        sb%lsize(2)/units_out%length%factor, ',',                &
        sb%lsize(3)/units_out%length%factor, ')'
    end if

    write(message(1), '(a,i1,a)') 'The octopus will run in ', sb%dim, ' dimension(s).'
    write(message(2), '(a,i1,a)') 'The octopus will treat the system as periodic in ', &
      sb%periodic_dim, ' dimension(s).'
    call write_info(2, iunit)

    if(sb%periodic_dim > 0) then
      write(message(1),'(1x)')
      write(message(2),'(a,3a,a)') 'Lattice Primitive Vectors [', trim(units_out%length%abbrev), ']'
      write(message(3),'(a,f8.3)') '    x axis ', sb%rlat(1,1)/units_out%length%factor
      write(message(4),'(a,f8.3)') '    y axis ', sb%rlat(2,2)/units_out%length%factor
      write(message(5),'(a,f8.3)') '    z axis ', sb%rlat(3,3) /units_out%length%factor
      write(message(6),'(a,3a,a)') 'Reciprocal Lattice Primitive Vectors [', trim(units_out%length%abbrev), '^-1]'
      write(message(7),'(a,f8.3)') '  k_x axis ', sb%klat(1,1)*units_out%length%factor
      write(message(8),'(a,f8.3)') '  k_y axis ', sb%klat(2,2)*units_out%length%factor
      write(message(9),'(a,f8.3)') '  k_z axis ', sb%klat(3,3)*units_out%length%factor
      call write_info(9, iunit)
    end if

    call pop_sub()
  end subroutine simul_box_write_info


  !--------------------------------------------------------------
  logical function simul_box_in_box(sb, geo, x) result(in_box)
    type(simul_box_t),  intent(in) :: sb
    type(geometry_t),   intent(in) :: geo
    FLOAT,              intent(in) :: x(:) ! x(3)

    FLOAT, parameter :: DELTA_R = CNST(1e-12)
    FLOAT :: r, re, im, xx(MAX_DIM)
    integer :: red, green, blue, ix, iy

    xx(:) = x(:) - sb%box_offset(:)

    select case(sb%box_shape)
    case(SPHERE)
      in_box = (sqrt(sum(xx(:)**2)) <= sb%rsize+DELTA_R)

    case(CYLINDER)
      r = sqrt(xx(2)**2 + xx(3)**2)
      in_box = (r<=sb%rsize+DELTA_R .and. abs(xx(1))<=sb%xsize+DELTA_R)

    case(MINIMUM)
      in_box = in_minimum()

    case(PARALLELEPIPED)
      in_box =  &
        (xx(1) >= -sb%lsize(1).and.xx(1) <= sb%lsize(1)).and. &
        (xx(2) >= -sb%lsize(2).and.xx(2) <= sb%lsize(2)).and. &
        (xx(3) >= -sb%lsize(3).and.xx(3) <= sb%lsize(3))

#if defined(HAVE_GDLIB)
    case(BOX_IMAGE)
      ix = int((xx(1) + sb%lsize(1))/sb%h(1))
      iy = int((xx(2) + sb%lsize(2))/sb%h(2))
      call loct_gdimage_get_pixel_rgb(sb%image, ix, iy, red, green, blue)
      in_box = (red == 255).and.(green == 255).and.(blue == 255)
#endif

    case(BOX_USDEF)
      ! is it inside the user given boundaries
      in_box =  &
        (xx(1) >= -sb%lsize(1).and.xx(1) <= sb%lsize(1)).and. &
        (xx(2) >= -sb%lsize(2).and.xx(2) <= sb%lsize(2)).and. &
        (xx(3) >= -sb%lsize(3).and.xx(3) <= sb%lsize(3))

      ! and inside the simulation box
      xx(:) = xx(:)/units_inp%length%factor ! convert from a.u. to input units
      r = sqrt(sum(xx(:)**2))
      call loct_parse_expression(re, im, xx(1), xx(2), xx(3), r, sb%user_def)
      in_box = in_box .and. (re .ne. M_ZERO)
    end select

  contains

    !--------------------------------------------------------------
    logical function in_minimum()
      integer :: i
      FLOAT :: radius

      in_minimum = .false.
      do i = 1, geo%natoms
        r = sqrt(sum((xx(:) - geo%atom(i)%x(:))**2))
        if(sb%rsize > M_ZERO) then
          radius = sb%rsize
        else
          radius = geo%atom(i)%spec%def_rsize
        endif
        if(r <= radius + DELTA_R) then
          in_minimum = .true.
          exit
        end if
      end do

    end function in_minimum

  end function simul_box_in_box


  !--------------------------------------------------------------
  logical function simul_box_is_periodic(sb)
    type(simul_box_t), intent(in) :: sb

    simul_box_is_periodic = sb%periodic_dim > 0
  end function simul_box_is_periodic


  !--------------------------------------------------------------
  subroutine simul_box_dump(sb, iunit)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: iunit
    write(iunit, '(a20,i4)')        'box_shape=          ',sb%box_shape
    write(iunit, '(a20,i4)')        'dim=                ',sb%dim
    write(iunit, '(a20,i4)')        'periodic_dim=       ',sb%periodic_dim
    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,3e22.14)') 'lsize=              ', sb%lsize(1:3)
    case(CYLINDER)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,e22.14)')  'xlength=            ', sb%xsize
      write(iunit, '(a20,3e22.14)') 'lsize=              ', sb%lsize(1:3)
    case(PARALLELEPIPED)
      write(iunit, '(a20,3e22.14)') 'lsize=              ', sb%lsize(1:3)
    case(BOX_USDEF)
      write(iunit, '(a20,3e22.14)') 'lsize=              ', sb%lsize(1:3)
      write(iunit, '(a20,a1024)')   'user_def=           ', sb%user_def
    end select
    write(iunit, '(a20,e22.14)')    'fft_alpha=          ', sb%fft_alpha
    write(iunit, '(a20,3e22.14)')   'h=                  ', sb%h(1:3)
    write(iunit, '(a20,3e22.14)')   'box_offset=         ', sb%box_offset(1:3)
  end subroutine simul_box_dump


  !--------------------------------------------------------------
  subroutine simul_box_init_from_file(sb, iunit)
    type(simul_box_t), intent(inout) :: sb
    integer,           intent(in)    :: iunit
    character(len=20) :: str
    read(iunit, *) str, sb%box_shape
    read(iunit, *) str, sb%dim
    read(iunit, *) str, sb%periodic_dim
    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      read(iunit, *) str, sb%rsize
      read(iunit, *) str, sb%lsize(1:3)
    case(CYLINDER)
      read(iunit, *) str, sb%rsize
      read(iunit, *) str, sb%xsize
      read(iunit, *) str, sb%lsize(1:3)
    case(PARALLELEPIPED)
      read(iunit, *) str, sb%lsize(1:3)
    case(BOX_USDEF)
      read(iunit, *) str, sb%lsize(1:3)
      read(iunit, *) str, sb%user_def
    end select
    read(iunit, *) str, sb%fft_alpha
    read(iunit, *) str, sb%h(1:3)
    read(iunit, *) str, sb%box_offset(1:3)
  end subroutine simul_box_init_from_file


  !--------------------------------------------------------------
  logical function simul_box_is_eq(sb1, sb2) result(res)
    type(simul_box_t), intent(in) :: sb1, sb2

    res = .false.
    if(sb1%box_shape .ne. sb2%box_shape)             return
    if(sb1%dim .ne. sb2%dim)                         return
    if(sb1%periodic_dim .ne. sb2%periodic_dim)       return
    select case(sb1%box_shape)
    case(SPHERE, MINIMUM)
      if(.not.(sb1%rsize .app. sb2%rsize))           return
      if(.not.(sb1%lsize .app. sb2%lsize))           return
    case(CYLINDER)
      if(.not.(sb1%rsize .app. sb2%rsize))           return
      if(.not.(sb1%xsize .app. sb2%xsize))           return
      if(.not.(sb1%lsize .app. sb2%lsize))           return
    case(PARALLELEPIPED)
      if(.not.(sb1%lsize .app. sb2%lsize))           return
    case(BOX_USDEF)
      if(.not.(sb1%lsize .app. sb2%lsize))           return
      if(trim(sb1%user_def) .ne. trim(sb2%user_def)) return
    end select
    if(.not.(sb1%fft_alpha .app. sb2%fft_alpha))     return
    if(.not.(sb1%h .app. sb2%h))                     return
    if(.not.(sb1%box_offset .app. sb2%box_offset))   return
    res = .true.

  end function simul_box_is_eq

end module simul_box_m
