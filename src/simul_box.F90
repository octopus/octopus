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

module simul_box
  use global
  use messages
  use varinfo
  use units
  use syslabels
  use lib_oct_parser
  use geometry

  implicit none

  private
  public ::                     &
    simul_box_type,             &
    simul_box_init,             &
    simul_box_write_info,       &
    simul_box_is_periodic,      &
    simul_box_in_box

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4

  type simul_box_type
    integer  :: box_shape   ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                            ! 4->parallelpiped (orthonormal, up to now).

    FLOAT :: h(3)           ! the (canonical) spacing between the points
    FLOAT :: box_offset(3)  ! shifts of the origin in the respective direction

    FLOAT :: rsize          ! the radius of the sphere or of the cylinder
    FLOAT :: xsize          ! the length of the cylinder in the x direction
    FLOAT :: lsize(3)       ! half of the length of the parallelepiped in each direction.

    FLOAT :: rlat(3,3)      ! lattice primitive vectors
    FLOAT :: klat(3,3)      ! reciprocal lattice primitive vectors
    FLOAT :: shift(27,3)    ! shift to equivalent positions in nearest neighbour primitive cells

    FLOAT :: fft_alpha      ! enlargement factor for double box

    integer :: dim
    integer :: periodic_dim

  end type simul_box_type


contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, geo)
    type(simul_box_type), intent(inout) :: sb
    type(geometry_type),  intent(in)    :: geo

    ! some local stuff
    FLOAT :: def_h, def_rsize
    integer :: i, ix, iy, iz, ii(3)

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

    end subroutine read_misc


    !--------------------------------------------------------------
    subroutine read_box()
      integer(POINTER_SIZE) :: blk

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
      !%End
      call loct_parse_int(check_inp('BoxShape'), MINIMUM, sb%box_shape)
      if(.not.varinfo_valid_option('BoxShape', sb%box_shape)) call input_error('BoxShape')
      select case(sb%box_shape)
      case(SPHERE,MINIMUM)
        if(sb%dim>1 .and. simul_box_is_periodic(sb)) call input_error('BoxShape')
      case(CYLINDER)
        if (sb%dim>2 .and. &
          ((sb%dim - sb%periodic_dim == 0) .or. (sb%dim - sb%periodic_dim == 1))) call input_error('BoxShape')
      end select

      ! ignore box_shape in 1D
      if(sb%dim==1.and.sb%box_shape /= PARALLELEPIPED) sb%box_shape=SPHERE

      sb%rsize = -M_ONE
      if(sb%box_shape == MINIMUM.and.def_rsize>M_ZERO) sb%rsize = def_rsize/units_inp%length%factor

      if(sb%box_shape == SPHERE.or.sb%box_shape == CYLINDER.or.sb%box_shape == MINIMUM) then
        !%Variable Radius
        !%Type float
        !%Section Mesh::Simulation Box
        !%Description
        !% If BoxShape is not "parallelepiped" defines the radius of the spheres or of the cylinder.
        !% It has to be a positive number. If it is not defined in the input file, then the program
        !% will attempt to fine a suitable default, but this is not always possible, in which case
        !% the code will stop issuing this error message.
        !%End
        call loct_parse_float(check_inp('radius'), sb%rsize, sb%rsize)
        if(sb%rsize < CNST(0.0)) call input_error('radius')
        sb%rsize = sb%rsize * units_inp%length%factor
        if(def_rsize>M_ZERO) call check_def(def_rsize, sb%rsize, 'radius')
      end if

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
      if(sb%box_shape == PARALLELEPIPED) then

        !%Variable Lsize
        !%Type block
        !%Section Mesh::Simulation Box
        !%Description
        !% In case BoxShape is "parallelepiped", this is assumed to be a block of the form:
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
        else
          call input_error('lsize')
        end if

        do i = 1, sb%dim
          call loct_parse_block_float(blk, 0, i-1, sb%lsize(i))
          if(def_rsize>M_ZERO.and.sb%periodic_dim<i) call check_def(def_rsize, sb%lsize(i), 'lsize')
        end do
        sb%lsize = sb%lsize*units_inp%length%factor

        call loct_parse_block_end(blk)
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
          sb%lsize(i)      = maxval(geo%atom(:)%x(i)) + sb%rsize
        end do
      end select

    end subroutine read_box


    !--------------------------------------------------------------
    subroutine read_spacing()
      integer :: i
      integer(POINTER_SIZE) :: blk

      sb%h = -M_ONE
      select case(sb%box_shape)
      case(SPHERE,CYLINDER,MINIMUM)
        call loct_parse_float(check_inp('spacing'), sb%h(1), sb%h(1))
        sb%h(1:sb%dim) = sb%h(1)

      case(PARALLELEPIPED)
        if(loct_parse_block(check_inp('spacing'), blk) == 0) then
          do i = 1, sb%dim
            call loct_parse_block_float(blk, 0, i-1, sb%h(i))
          end do
          call loct_parse_block_end(blk)
        else
          message(1) = '"Spacing" is a block if BoxShape == parallelepiped.'
          call write_fatal(1)
        end if
      end select

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
            message(2) = "   *) variable 'spacing' is not defined"
            message(3) = "   *) your input for 'spacing' is negative"
            message(4) = "   *) I can't find a suitable default"
            call write_fatal(4)
          end if
        end if
        if(def_rsize>M_ZERO) call check_def(sb%h(i), def_rsize, 'spacing')
      end do

    end subroutine read_spacing


    !--------------------------------------------------------------
    subroutine read_box_offset()
      integer :: i
      integer(POINTER_SIZE) :: blk

      sb%box_offset = M_ZERO
      select case(sb%box_shape)
      case(PARALLELEPIPED)
        if(loct_parse_block(check_inp('BoxOffset'), blk) == 0) then
          do i = 1, sb%dim
            call loct_parse_block_float(blk, 0, i-1, sb%box_offset(i))
          end do
          call loct_parse_block_end(blk)
        end if
      end select

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

    end subroutine build_lattice

  end subroutine simul_box_init


  !--------------------------------------------------------------
  subroutine simul_box_end()

  end subroutine simul_box_end


  !--------------------------------------------------------------
  subroutine simul_box_write_info(sb, iunit)
    type(simul_box_type), intent(in) :: sb
    integer,              intent(in) :: iunit

    character(len=15), parameter :: bs(4) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'around nuclei ', &
      'parallelepiped'/)

    call push_sub('simul_box.simul_box_write_info')

    write(message(1),'(a)') 'Simulation Box:'
    write(message(2), '(a,a,1x)') '  Type = ', bs(sb%box_shape)
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

    write(message(1),'(3a, a, f6.3, a, f6.3, a, f6.3, a)')     &
      '  Lengths [', trim(units_out%length%abbrev), '] = ',    &
      '(', sb%lsize(1)/units_out%length%factor, ',',           &
      sb%lsize(2)/units_out%length%factor, ',',                &
      sb%lsize(3)/units_out%length%factor, ')'


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
    type(simul_box_type),  intent(in) :: sb
    type(geometry_type),   intent(in) :: geo
    FLOAT,                 intent(in) :: x(3) ! x(3)

    FLOAT, parameter :: DELTA_R = CNST(1e-12)
    FLOAT :: r

    select case(sb%box_shape)
    case(SPHERE)
      in_box = (sqrt(sum(x**2)) <= sb%rsize+DELTA_R)
    case(CYLINDER)
      r = sqrt(x(2)**2 + x(3)**2)
      in_box = (r<=sb%rsize+DELTA_R .and. abs(x(1))<=sb%xsize+DELTA_R)
    case(MINIMUM)
      in_box = in_minimum()
    case(PARALLELEPIPED)
      in_box =  &
        (x(1) >= -sb%lsize(1) + sb%box_offset(1).and.x(1) <= sb%lsize(1) + sb%box_offset(1)).and. &
        (x(2) >= -sb%lsize(2) + sb%box_offset(2).and.x(2) <= sb%lsize(2) + sb%box_offset(2)).and. &
        (x(3) >= -sb%lsize(3) + sb%box_offset(3).and.x(3) <= sb%lsize(2) + sb%box_offset(3))
    end select

  contains

    !--------------------------------------------------------------
    logical function in_minimum()
      integer :: i

      in_minimum = .false.
      do i = 1, geo%natoms
        r = sqrt(sum((x(:) - geo%atom(i)%x(:))**2))
        if(r <= sb%rsize+DELTA_R) then
          in_minimum = .true.
          exit
        end if
      end do

    end function in_minimum

  end function simul_box_in_box


  !--------------------------------------------------------------
  logical function simul_box_is_periodic(sb)
    type(simul_box_type), intent(in) :: sb

    simul_box_is_periodic = sb%periodic_dim > 0
  end function simul_box_is_periodic

end module simul_box
