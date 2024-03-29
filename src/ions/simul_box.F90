!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2013,2021 M. Oliveira
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
!!

#include "global.h"

module simul_box_oct_m
  use box_oct_m
  use box_cylinder_oct_m
  use box_hypercube_oct_m
  use box_image_oct_m
  use box_minimum_oct_m
  use box_parallelepiped_oct_m
  use box_sphere_oct_m
  use box_user_defined_oct_m
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use iso_c_binding
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use lookup_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                     &
    simul_box_t,                &
    simul_box_init,             &
    simul_box_end,              &
    simul_box_copy

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4,         &
    BOX_IMAGE      = 5,         &
    HYPERCUBE      = 6,         &
    BOX_USDEF      = 77
  !< BOX_USDEF shares a number with other 'user_defined' input file options.

  type, extends(box_t) :: simul_box_t
    ! Components are public by default

    class(box_t), pointer :: box => NULL() !< This is a temporary placeholder
                                           !! for storing the boxes making up
                                           !! the simulation box, until we can
                                           !! make simul_box_t an extension of a
                                           !! multibox_t (and thus store this in
                                           !! a list).

    FLOAT :: lsize(MAX_DIM) !< half of the length of the parallelepiped in each direction.

    type(lattice_vectors_t), pointer :: latt => NULL()
    
    FLOAT :: stress_tensor(MAX_DIM,MAX_DIM)   !< reciprocal-lattice primitive vectors

    integer :: periodic_dim
  contains
    procedure :: contains_points => simul_box_contains_points
    procedure :: write_info => simul_box_write_info
    procedure :: short_info => simul_box_short_info
  end type simul_box_t

contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, namespace, ions, space)
    type(simul_box_t),                   intent(inout) :: sb
    type(namespace_t),                   intent(in)    :: namespace
    type(ions_t),      target,           intent(inout) :: ions
    type(space_t),                       intent(in)    :: space

    ! some local stuff
    integer :: box_shape
    FLOAT :: def_h, def_rsize, center(space%dim), rsize, xsize
    integer :: n_site_types, n_sites
    integer, allocatable :: site_type(:)
    FLOAT,   allocatable :: site_type_radius(:)
    character(len=LABEL_LEN), allocatable :: site_type_label(:)
    character(len=200) :: filename
    character(len=1024) :: user_def

    PUSH_SUB(simul_box_init)

    call ions%grid_defaults(def_h, def_rsize)

    sb%dim = space%dim
    sb%periodic_dim = space%periodic_dim

    sb%latt => ions%latt

    call read_box()                        ! Parameters defining the simulation box.

    center = M_ZERO ! Currently all the boxes have to be centered at the origin.
    select case (box_shape)
    case (SPHERE)
      sb%box => box_sphere_t(space%dim, center, rsize)
    case (CYLINDER)
      sb%box => box_cylinder_t(space%dim, center, rsize, 1, M_TWO*xsize, namespace, periodic_boundaries=(space%periodic_dim > 0))
    case (PARALLELEPIPED)
      sb%box => box_parallelepiped_t(space%dim, center, M_TWO*sb%lsize(1:space%dim), n_periodic_boundaries=space%periodic_dim)
    case (HYPERCUBE)
      sb%box => box_hypercube_t(space%dim, center, M_TWO*sb%lsize(1:space%dim), n_periodic_boundaries=space%periodic_dim)
    case (BOX_USDEF)
      sb%box => box_user_defined_t(space%dim, center, user_def, M_TWO*sb%lsize(1:space%dim))

    case (MINIMUM)
      sb%box => box_minimum_t(space%dim, n_site_types, site_type_label, site_type_radius, n_sites, site_type, ions%pos)
      SAFE_DEALLOCATE_A(site_type_label)
      SAFE_DEALLOCATE_A(site_type_radius)
      SAFE_DEALLOCATE_A(site_type)

    case (BOX_IMAGE)
      sb%box => box_image_t(center, sb%lsize, filename, space%periodic_dim, namespace)
    end select

    POP_SUB(simul_box_init)
  contains

    !--------------------------------------------------------------
    subroutine read_box()
      type(block_t) :: blk

      FLOAT :: default, factor(space%dim)
      integer :: default_boxshape, idir, iatom, ispec

      PUSH_SUB(simul_box_init.read_box)
      ! Read box shape.

      !%Variable BoxShape
      !%Type integer
      !%Section Mesh::Simulation Box
      !%Description
      !% This variable decides the shape of the simulation box.
      !% The default is <tt>minimum</tt> for finite systems and <tt>parallelepiped</tt> for periodic systems.
      !% Note that some incompatibilities apply:
      !% <ul><li>Spherical or minimum mesh is not allowed for periodic systems.
      !% <li>Cylindrical mesh is not allowed for systems that are periodic in more than one dimension.
      !% <li><tt>box_image</tt> is only allowed in 2D.</ul>
      !%Option sphere 1
      !% The simulation box will be a sphere of radius <tt>Radius</tt>. (In 2D, this is a circle.)
      !%Option cylinder 2
      !% The simulation box will be a cylinder with radius <tt>Radius</tt> and height (in the <i>x</i>-direction)
      !% of 2 <tt>Xlength</tt>.
      !%Option minimum 3
      !% The simulation box will be constructed by adding spheres created around each
      !% atom (or user-defined potential), of radius <tt>Radius</tt>.
      !%Option parallelepiped 4
      !% The simulation box will be a parallelepiped whose dimensions are taken from
      !% the variable <tt>Lsize</tt>.
      !%Option box_image 5
      !% The simulation box will be defined through an image, specified with <tt>BoxShapeImage</tt>.
      !% White (RGB = 255,255,255) means that the point
      !% is contained in the simulation box, while any other color means that the point is out.
      !% The image will be scaled to fit <tt>Lsize</tt>, while its resolution will define the default <tt>Spacing</tt>.
      !% The actual box may be slightly larger than <tt>Lsize</tt> to ensure one grid point = one pixel for
      !% default <tt>Spacing</tt>.
      !%Option user_defined 77
      !% The shape of the simulation box will be read from the variable <tt>BoxShapeUsDef</tt>.
      !%Option hypercube 6
      !% (experimental) The simulation box will be a hypercube or
      !% hyperparallelepiped. This is equivalent to the
      !% <tt>parallelepiped</tt> box but it can work with an arbitrary
      !% number of dimensions.
      !%End

      if (space%is_periodic()) then
        default_boxshape = PARALLELEPIPED
      else
        default_boxshape = MINIMUM
      end if
      call parse_variable(namespace, 'BoxShape', default_boxshape, box_shape)
      if(.not.varinfo_valid_option('BoxShape', box_shape)) call messages_input_error(namespace, 'BoxShape')
      select case(box_shape)
      case(SPHERE, MINIMUM, BOX_USDEF)
        if(sb%dim > 1 .and. space%is_periodic()) call messages_input_error(namespace, 'BoxShape')
      case(CYLINDER)
        if(sb%dim == 2) then
          message(1) = "BoxShape = cylinder is not meaningful in 2D. Use sphere if you want a circle."
          call messages_fatal(1, namespace=namespace)
        end if
        if(space%periodic_dim > 1) call messages_input_error(namespace, 'BoxShape')
      end select

      ! ignore box_shape in 1D
      if (sb%dim == 1 .and. box_shape /= PARALLELEPIPED .and. box_shape /= HYPERCUBE) then
        box_shape = SPHERE
      end if

      ! Cannot use images in 1D or 3D
      if(sb%dim /= 2 .and. box_shape == BOX_IMAGE) call messages_input_error(namespace, 'BoxShape')

      if(sb%dim > 3 .and. box_shape /= HYPERCUBE) then
        message(1) = "For more than 3 dimensions, you can only use the hypercubic box."
        call messages_fatal(1, namespace=namespace)
        ! FIXME: why not a hypersphere as another option?
        ! Also, hypercube should be unified with parallepiped.
      end if

      rsize = -M_ONE
      !%Variable Radius
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% Defines the radius for <tt>BoxShape</tt> = <tt>sphere</tt>,
      !% <tt>cylinder</tt>, or <tt>minimum</tt>.  Must be a positive
      !% number. If not specified, the code will look for values in
      !% the <tt>Species</tt> block, or, from the default
      !% pseudopotential parameters.  In these cases, for
      !% <tt>minimum</tt>, a different radius is used for each
      !% species, while for other shapes, the maximum radius is used.
      !%End
      select case(box_shape)
      case(SPHERE, CYLINDER)
        call parse_variable(namespace, 'Radius', def_rsize, rsize, units_inp%length)
        if(rsize < M_ZERO) call messages_input_error(namespace, 'radius')
        if(def_rsize>M_ZERO) call messages_check_def(rsize, .false., def_rsize, 'radius', units_out%length)
      case(MINIMUM)
        default=rsize
        call parse_variable(namespace, 'radius', default, rsize, units_inp%length)
        if(rsize < M_ZERO .and. def_rsize < M_ZERO) call messages_input_error(namespace, 'Radius')

        n_site_types = ions%nspecies
        SAFE_ALLOCATE(site_type_label(1:ions%nspecies))
        SAFE_ALLOCATE(site_type_radius(1:ions%nspecies))

        do ispec = 1, ions%nspecies
          site_type_label(ispec) = species_label(ions%species(ispec))
          if (rsize > M_ZERO) then
            site_type_radius(ispec) = rsize
          else
            if (species_def_rsize(ions%species(ispec)) < -M_EPSILON) then
              write(message(1),'(a,a,a)') 'Using default radii for minimum box, but radius for ', &
                trim(species_label(ions%species(ispec))), ' is negative or undefined.'
              message(2) = "Define it properly in the Species block or set the Radius variable explicitly."
              call messages_fatal(2, namespace=namespace)
            else
              site_type_radius(ispec) = species_def_rsize(ions%species(ispec))
            end if
          end if
        end do

        n_sites = ions%natoms
        SAFE_ALLOCATE(site_type(1:ions%natoms))
        do iatom = 1, ions%natoms
          do ispec = 1, ions%nspecies
            if (ions%atom(iatom)%label == site_type_label(ispec)) then
              site_type(iatom) = ispec
            end if
          end do
        end do

      end select

      if(box_shape == CYLINDER) then
        !%Variable Xlength
        !%Default <tt>Radius</tt>
        !%Type float
        !%Section Mesh::Simulation Box
        !%Description
        !% If <tt>BoxShape</tt> is <tt>cylinder</tt>, the total length of the cylinder is twice <tt>Xlength</tt>.
        !% Note that when PeriodicDimensions = 1, then the length of the cylinder is determined from the lattice vectors.
        !%End
        if(rsize > M_ZERO) then
          default = rsize
        else
          default = def_rsize
        end if

        if (space%is_periodic()) then
          xsize = sqrt(sum(sb%latt%rlattice(1:space%periodic_dim, 1)**2))/M_TWO
        else
          call parse_variable(namespace, 'Xlength', default, xsize, units_inp%length)
          if (def_rsize > M_ZERO .and. space%periodic_dim == 0) then
            call messages_check_def(xsize, .false., def_rsize, 'xlength', units_out%length)
          end if
        end if
      end if

      sb%lsize = M_ZERO
      if(box_shape == PARALLELEPIPED .or. box_shape == HYPERCUBE .or. &
         box_shape == BOX_IMAGE .or. box_shape == BOX_USDEF) then

        !%Variable Lsize
        !%Type block
        !%Section Mesh::Simulation Box
        !%Description
        !% If <tt>BoxShape</tt> is <tt>parallelepiped</tt>, <tt>hypercube</tt>,
        !% <tt>box_image</tt>, or <tt>user_defined</tt>, this is a
        !% block of the form:
        !%
        !% <tt>%Lsize
        !% <br>&nbsp;&nbsp;sizex | sizey | sizez | ...
        !% <br>%</tt>
        !%
        !% where the <tt>size*</tt> are half the lengths of the box in each direction.
        !%
        !% The number of columns must match the dimensionality of the
        !% calculation. If you want a cube you can also set <tt>Lsize</tt> as a
        !% single variable.
        !%End

        ! lsize along the periodic dimensions must always be set from the norm of the lattice vectors
        do idir = 1, space%periodic_dim
          sb%lsize(idir) = sqrt(sum(sb%latt%rlattice(1:space%dim, idir)**2))/M_TWO
        end do

        if (space%is_periodic()) then
          ! For mixed-periodicity, lsize along the non-periodic dimensions is
          ! by default set from the lattice parameters (this can still be
          ! overriden by setting Lsize, see bellow).
          do idir = space%periodic_dim + 1, space%dim
            sb%lsize(idir) = sqrt(sum(sb%latt%rlattice(1:space%dim, idir)**2))/M_TWO
          end do

          ! Now we renormalize the lattice parameters along the non-periodic dimensions
          do idir = 1, space%periodic_dim
            factor(idir) = M_ONE
          end do
          do idir = space%periodic_dim + 1, space%dim
            factor(idir) = M_ONE/sqrt(sum(sb%latt%rlattice(1:space%dim, idir)**2))
          end do
          call sb%latt%scale(factor)

        else
          ! Lsize must be set for finite systems, as in that case we do not have the lattice parameters
          if (.not. parse_is_defined(namespace, 'Lsize')) then
            call messages_input_error(namespace, 'Lsize', 'Lsize is required for finite systems')
          end if
        end if

        ! Note that for cases with mixed-periodicidy, the user still has the
        ! option to set Lsize to override the size of the box along the
        ! non-periodic dimensions given by the lattice parameters. This
        ! requires the user to also set Lsize for the periodic dimensions,
        ! which at the moment must match exactly the corresponding values
        ! given by the lattice vectors.
        if (parse_block(namespace, 'Lsize', blk) == 0) then
          ! Lsize is specified as a block
          if (parse_block_cols(blk, 0) < space%dim) then
            call messages_input_error(namespace, 'Lsize')
          end if

          do idir = 1, space%dim
            call parse_block_float(blk, 0, idir - 1, sb%lsize(idir), units_inp%length)
            if (def_rsize > M_ZERO) then
              call messages_check_def(sb%lsize(idir), .false., def_rsize, 'Lsize', units_out%length)
            end if
          end do
          call parse_block_end(blk)

        else if (parse_is_defined(namespace, 'Lsize')) then
          ! Lsize is specified as a scalar
          call parse_variable(namespace, 'Lsize', -M_ONE, sb%lsize(1), units_inp%length)
          if (abs(sb%lsize(1) + M_ONE) <= M_EPSILON) then
            call messages_input_error(namespace, 'Lsize')
          end if
          if (def_rsize > M_ZERO) then
            call messages_check_def(sb%lsize(1), .false., def_rsize, 'Lsize', units_out%length)
          end if
          do idir = 2, space%dim
            sb%lsize(idir) = sb%lsize(1)
          end do

        end if

        ! Check that lsize is consistent with the lattice vectors along the periodic dimensions
        do idir = 1, space%periodic_dim
          if (abs(M_TWO*sb%lsize(idir) - sqrt(sum(sb%latt%rlattice(1:space%dim, idir)**2))) > M_EPSILON) then
            call messages_input_error(namespace, 'Lsize', &
              'Lsize must be exactly half the length of the lattice vectors along periodic dimensions')
          end if
        end do
      end if

      ! read in image for box_image
      if(box_shape == BOX_IMAGE) then

        !%Variable BoxShapeImage
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Name of the file that contains the image that defines the simulation box
        !% when <tt>BoxShape = box_image</tt>. No default. Will search in current
        !% directory and <tt>OCTOPUS-HOME/share/</tt>.
        !%End
#if defined(HAVE_GDLIB)
        call parse_variable(namespace, 'BoxShapeImage', '', filename)
        if(trim(filename) == "") then
          message(1) = "Must specify BoxShapeImage if BoxShape = box_image."
          call messages_fatal(1, namespace=namespace)
        end if
#else
        message(1) = "To use 'BoxShape = box_image', you have to compile Octopus"
        message(2) = "with GD library support."
        call messages_fatal(2, namespace=namespace)
#endif
      end if

      ! read in box shape for user-defined boxes
      if(box_shape == BOX_USDEF) then

        !%Variable BoxShapeUsDef
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Boolean expression that defines the interior of the simulation box. For example,
        !% <tt>BoxShapeUsDef = "(sqrt(x^2+y^2) <= 4) && z>-2 && z<2"</tt> defines a cylinder
        !% with axis parallel to the <i>z</i>-axis.
        !%End

        call parse_variable(namespace, 'BoxShapeUsDef', 'x^2+y^2+z^2 < 4', user_def)
      end if

      ! fill in lsize structure
      select case(box_shape)
      case(SPHERE)
        sb%lsize(1:sb%dim) = rsize
      case(CYLINDER)
        sb%lsize(1)        = xsize
        sb%lsize(2:sb%dim) = rsize
      case(MINIMUM)
        do idir = 1, sb%dim
          if(rsize > M_ZERO) then
            sb%lsize(idir) = maxval(abs(ions%pos(idir, :))) + rsize
          else
            sb%lsize(idir) = maxval(abs(ions%pos(idir, :))) + def_rsize
          end if
        end do
      end select

      call messages_obsolete_variable(namespace, 'BoxOffset')
      
      POP_SUB(simul_box_init.read_box)
    end subroutine read_box

  end subroutine simul_box_init

  !--------------------------------------------------------------
  subroutine simul_box_end(sb)
    type(simul_box_t), intent(inout) :: sb    

    class(box_t), pointer :: box

    PUSH_SUB(simul_box_end)

    nullify(sb%latt)

    ! We first need to bet a pointer to the box to deallocated it because of a
    ! bug in gfortran.
    box => sb%box
    SAFE_DEALLOCATE_P(box)

    POP_SUB(simul_box_end)
  end subroutine simul_box_end


  !--------------------------------------------------------------
  recursive subroutine simul_box_write_info(this, iunit)
    class(simul_box_t), intent(in) :: this
    integer,            intent(in) :: iunit

    PUSH_SUB(simul_box_write_info)

    call this%box%write_info(iunit)

    POP_SUB(simul_box_write_info)
  end subroutine simul_box_write_info

  character(len=BOX_INFO_LEN) function simul_box_short_info(this, unit_length) result(info)
    class(simul_box_t), intent(in) :: this
    type(unit_t),       intent(in) :: unit_length

    PUSH_SUB(simul_box_short_info)

    info = this%box%short_info(unit_length)

    POP_SUB(simul_box_short_info)
  end function simul_box_short_info

  !--------------------------------------------------------------
  !> Checks if a group of mesh points belong to the actual mesh.
  function simul_box_contains_points(this, nn, xx) result(contained)
    class(simul_box_t), intent(in)  :: this
    integer,            intent(in)  :: nn
    FLOAT,              intent(in)  :: xx(:, :) !< (1:, 1:this%dim)
    logical :: contained(1:nn)

    FLOAT, allocatable :: xx_red(:, :)

    ! no push_sub because this function is called very frequently
    SAFE_ALLOCATE(xx_red(1:nn, 1:this%dim))
    xx_red = M_ZERO
    
    !convert from Cartesian to reduced lattice coord 
    if (nn == 1) then
      xx_red(1, 1:this%dim) = matmul(xx(1, 1:this%dim), this%latt%klattice_primitive(1:this%dim, 1:this%dim))
    else
      call lalg_gemm(nn, this%dim, this%dim, M_ONE, xx, this%latt%klattice_primitive, M_ZERO, xx_red)
    end if

    contained = this%box%contains_points(nn, xx_red)

    SAFE_DEALLOCATE_A(xx_red)

  end function simul_box_contains_points

  ! --------------------------------------------------------------
  recursive subroutine simul_box_copy(sbout, sbin)
    type(simul_box_t), intent(out) :: sbout
    type(simul_box_t), intent(in)  :: sbin

    PUSH_SUB(simul_box_copy)

    sbout%lsize          = sbin%lsize
    sbout%latt          => sbin%latt
    sbout%dim            = sbin%dim
    sbout%periodic_dim   = sbin%periodic_dim

    POP_SUB(simul_box_copy)
  end subroutine simul_box_copy

end module simul_box_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
