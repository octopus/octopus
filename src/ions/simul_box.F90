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
  use atom_oct_m
  use box_oct_m
  use iso_c_binding
  use gdlib_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lookup_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m
  use string_oct_m
  use symm_op_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                     &
    simul_box_t,                &
    simul_box_init,             &
    simul_box_lookup_init,      &
    simul_box_end,              &
    simul_box_write_info,       &
    simul_box_write_short_info, &
    simul_box_is_periodic,      &
    simul_box_has_zero_bc,      &
    simul_box_atoms_in_box,     &
    simul_box_copy,             &
    simul_box_periodic_atom_in_box, &
    simul_box_symmetry_check,   &
    reciprocal_lattice

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
    type(symmetries_t) :: symm
    !> 1->sphere, 2->cylinder, 3->sphere around each atom,
    !! 4->parallelepiped (orthonormal, up to now).
    integer  :: box_shape   

    FLOAT :: rsize          !< the radius of the sphere or of the cylinder
    FLOAT :: xsize          !< the length of the cylinder in the x-direction
    FLOAT :: lsize(MAX_DIM) !< half of the length of the parallelepiped in each direction.

    type(lookup_t), private :: atom_lookup

    type(geometry_t), pointer, private :: geo

    character(len=1024), private :: user_def !< for the user-defined box

    FLOAT :: rlattice_primitive(MAX_DIM,MAX_DIM)   !< lattice primitive vectors
    FLOAT :: rlattice          (MAX_DIM,MAX_DIM)   !< lattice vectors
    FLOAT :: klattice_primitive(MAX_DIM,MAX_DIM)   !< reciprocal-lattice primitive vectors
    FLOAT :: klattice          (MAX_DIM,MAX_DIM)   !< reciprocal-lattice vectors
    FLOAT, private :: volume_element               !< the volume element in real space
    FLOAT :: surface_element   (MAX_DIM)         !< surface element in real space
    FLOAT :: rcell_volume                        !< the volume of the cell in real space
    FLOAT :: alpha, beta, gamma                  !< the angles defining the cell
    FLOAT :: rmetric            (MAX_DIM,MAX_DIM) !< metric for the real space lattice vectors
    FLOAT :: stress_tensor(MAX_DIM,MAX_DIM)   !< reciprocal-lattice primitive vectors
    logical :: nonorthogonal
    
    type(kpoints_t) :: kpoints                   !< the k-points

    integer :: periodic_dim

    !> for the box defined through an image
    integer             :: image_size(1:2)
    type(c_ptr), private         :: image
    character(len=200), private  :: filename

  contains
    procedure :: contains_points => simul_box_contains_points
  end type simul_box_t

contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, namespace, geo, space)
    type(simul_box_t),                   intent(inout) :: sb
    type(namespace_t),                   intent(in)    :: namespace
    type(geometry_t), target,            intent(inout) :: geo
    type(space_t),                       intent(in)    :: space

    ! some local stuff
    FLOAT :: def_h, def_rsize
    logical :: only_gamma_kpoint

    PUSH_SUB(simul_box_init)

    call geometry_grid_defaults(geo, def_h, def_rsize)

    sb%geo => geo

    call read_misc()                       ! Miscellaneous stuff.
    call read_box()                        ! Parameters defining the simulation box.
    call simul_box_lookup_init(sb, geo)
    call simul_box_build_lattice(sb, namespace)       ! Build lattice vectors.
    call simul_box_atoms_in_box(sb, geo, namespace, .true.)   ! Put all the atoms inside the box.

    call simul_box_check_atoms_are_too_close(geo, sb, namespace)

    call symmetries_init(sb%symm, namespace, geo, sb%dim, sb%periodic_dim, sb%rlattice, sb%klattice)

    ! we need k-points for periodic systems
    only_gamma_kpoint = (sb%periodic_dim == 0)
    call kpoints_init(sb%kpoints, namespace, sb%symm, sb%dim, sb%rlattice, sb%klattice, only_gamma_kpoint)

    call simul_box_symmetry_check(sb, geo, sb%dim, namespace)

    POP_SUB(simul_box_init)

  contains


    !--------------------------------------------------------------
    subroutine read_misc()
      PUSH_SUB(simul_box_init.read_misc)

      sb%dim = space%dim

      !%Variable PeriodicDimensions
      !%Type integer
      !%Default 0
      !%Section System
      !%Description
      !% Define how many directions are to be considered periodic. It has to be a number
      !% between zero and <tt>Dimensions</tt>.
      !% (WARNING: For systems that are periodic in 1D and  2D, interaction between ions is assumed to be periodic in 3D.
      !% This affects the calculation of total energy and forces.)
      !%Option 0
      !% No direction is periodic (molecule).
      !%Option 1
      !% The <i>x</i> direction is periodic (wire, polymer).
      !%Option 2
      !% The <i>x</i> and <i>y</i> directions are periodic (slab).
      !%Option 3
      !% The <i>x</i>, <i>y</i>, and <i>z</i> directions are periodic (bulk).
      !%End

      if(geo%periodic_dim == -1) then
        call parse_variable(namespace, 'PeriodicDimensions', 0, sb%periodic_dim)
      else
        sb%periodic_dim = geo%periodic_dim
      end if
      if ((sb%periodic_dim < 0) .or. (sb%periodic_dim > MAX_DIM) .or. (sb%periodic_dim > sb%dim)) then
        call messages_input_error(namespace, 'PeriodicDimensions')
      end if

      if(sb%periodic_dim > 0 .and. sb%periodic_dim < sb%dim) then
        call messages_experimental('Support for mixed periodicity systems')
      end if

      if(sb%periodic_dim == 1) then
        call messages_write('For systems that  are periodic in 1D, interaction between', new_line = .true.)
        call messages_write('ions is assumed to be periodic in 3D. This affects the calculation', new_line = .true.)
        call messages_write('of total energy and forces.')
        call messages_warning(namespace=namespace)
      end if

      POP_SUB(simul_box_init.read_misc)
    end subroutine read_misc


    !--------------------------------------------------------------
    subroutine read_box()
      type(block_t) :: blk

      FLOAT :: default
      integer :: default_boxshape, idir, iatom
#if defined(HAVE_GDLIB)
      logical :: found
      integer :: box_npts
#endif

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

      if(simul_box_is_periodic(sb)) then
        default_boxshape = PARALLELEPIPED
      else
        default_boxshape = MINIMUM
      end if
      call parse_variable(namespace, 'BoxShape', default_boxshape, sb%box_shape)
      if(.not.varinfo_valid_option('BoxShape', sb%box_shape)) call messages_input_error(namespace, 'BoxShape')
      select case(sb%box_shape)
      case(SPHERE, MINIMUM, BOX_USDEF)
        if(sb%dim > 1 .and. simul_box_is_periodic(sb)) call messages_input_error(namespace, 'BoxShape')
      case(CYLINDER)
        if(sb%dim == 2) then
          message(1) = "BoxShape = cylinder is not meaningful in 2D. Use sphere if you want a circle."
          call messages_fatal(1, namespace=namespace)
        end if
        if(sb%periodic_dim > 1) call messages_input_error(namespace, 'BoxShape')
      end select

      ! ignore box_shape in 1D
      if(sb%dim == 1 .and. sb%box_shape /= PARALLELEPIPED .and. sb%box_shape /= HYPERCUBE) &
        sb%box_shape = SPHERE

      ! Cannot use images in 1D or 3D
      if(sb%dim /= 2 .and. sb%box_shape == BOX_IMAGE) call messages_input_error(namespace, 'BoxShape')

      if(sb%dim > 3 .and. sb%box_shape /= HYPERCUBE) then
        message(1) = "For more than 3 dimensions, you can only use the hypercubic box."
        call messages_fatal(1, namespace=namespace)
        ! FIXME: why not a hypersphere as another option?
        ! Also, hypercube should be unified with parallepiped.
      end if

      sb%rsize = -M_ONE
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
      select case(sb%box_shape)
      case(SPHERE, CYLINDER)
        call parse_variable(namespace, 'Radius', def_rsize, sb%rsize, units_inp%length)
        if(sb%rsize < M_ZERO) call messages_input_error(namespace, 'radius')
        if(def_rsize>M_ZERO) call messages_check_def(sb%rsize, .false., def_rsize, 'radius', units_out%length)
      case(MINIMUM)

        if(geo%reduced_coordinates) then
          message(1) = "The 'minimum' box shape cannot be used if atomic positions"
          message(2) = "are given as reduced coordinates."
          call messages_fatal(2, namespace=namespace)
        end if

        default=sb%rsize
        call parse_variable(namespace, 'radius', default, sb%rsize, units_inp%length)
        if(sb%rsize < M_ZERO .and. def_rsize < M_ZERO) call messages_input_error(namespace, 'Radius')

        if (sb%rsize <= M_ZERO) then
          do iatom = 1, sb%geo%natoms
            if (species_def_rsize(sb%geo%atom(iatom)%species) < -M_EPSILON) then
              write(message(1),'(a,a,a)') 'Using default radii for minimum box, but radius for ', &
                trim(species_label(sb%geo%atom(iatom)%species)), ' is negative or undefined.'
              message(2) = "Define it properly in the Species block or set the Radius variable explicitly."
              call messages_fatal(2, namespace=namespace)
            end if
          end do
        end if

      end select

      if(sb%box_shape == CYLINDER) then
        !%Variable Xlength
        !%Default <tt>Radius</tt>
        !%Type float
        !%Section Mesh::Simulation Box
        !%Description
        !% If <tt>BoxShape</tt> is <tt>cylinder</tt>, the total length of the cylinder is twice <tt>Xlength</tt>.
        !%End
        if(sb%rsize > M_ZERO) then
          default = sb%rsize
        else
          default = def_rsize
        end if

        call parse_variable(namespace, 'Xlength', default, sb%xsize, units_inp%length)
        if(def_rsize > M_ZERO .and. sb%periodic_dim == 0) &
          call messages_check_def(sb%xsize, .false., def_rsize, 'xlength', units_out%length)
      end if

      sb%lsize = M_ZERO
      if(sb%box_shape == PARALLELEPIPED .or. sb%box_shape == HYPERCUBE .or. &
         sb%box_shape == BOX_IMAGE .or. sb%box_shape == BOX_USDEF) then

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

        if(all(geo%lsize(1:sb%dim) > M_ZERO)) then
          ! use value read from XSF lattice vectors
          sb%lsize(:) = geo%lsize(:)
        else if(parse_block(namespace, 'Lsize', blk) == 0) then
          if(parse_block_cols(blk,0) < sb%dim .and. .not. parse_is_defined(namespace, 'LatticeVectors')) &
              call messages_input_error(namespace, 'Lsize')
          do idir = 1, sb%dim
            call parse_block_float(blk, 0, idir - 1, sb%lsize(idir), units_inp%length)
            if(def_rsize > M_ZERO .and. sb%periodic_dim < idir) &
              call messages_check_def(sb%lsize(idir), .false., def_rsize, 'Lsize', units_out%length)
          end do
          call parse_block_end(blk)
        else if ((parse_is_defined(namespace, 'Lsize'))) then
          call parse_variable(namespace, 'Lsize', -M_ONE, sb%lsize(1), units_inp%length)
          if(abs(sb%lsize(1)+M_ONE)  <=  M_EPSILON) then
            call messages_input_error(namespace, 'Lsize')
          end if
          if(def_rsize > M_ZERO .and. sb%periodic_dim < sb%dim) &
            call messages_check_def(sb%lsize(1), .false., def_rsize, 'Lsize', units_out%length)
          sb%lsize(1:sb%dim) = sb%lsize(1)
        else
          message(1) = "Lsize was not found in input file. Continuing anyway."
          call messages_warning(1, namespace=namespace)
        end if
      else
        ! if not a compatible box-shape
        if(all(geo%lsize(1:sb%dim) > M_ZERO)) then
          message(1) = "Ignoring lattice vectors from XSF file."
          call messages_warning(1, namespace=namespace)
        end if
      end if

      ! read in image for box_image
      if(sb%box_shape == BOX_IMAGE) then

        !%Variable BoxShapeImage
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Name of the file that contains the image that defines the simulation box
        !% when <tt>BoxShape = box_image</tt>. No default. Will search in current
        !% directory and <tt>OCTOPUS-HOME/share/</tt>.
        !%End
#if defined(HAVE_GDLIB)
        call parse_variable(namespace, 'BoxShapeImage', '', sb%filename)
        if(trim(sb%filename) == "") then
          message(1) = "Must specify BoxShapeImage if BoxShape = box_image."
          call messages_fatal(1, namespace=namespace)
        end if

        ! Find out the file and read it.
        inquire(file=trim(sb%filename), exist=found)
        if(.not. found) then
          message(1) = "Could not find file '" // trim(sb%filename) // "' for BoxShape = box_image."

          sb%filename = trim(conf%share) // '/' // trim(sb%filename)
          inquire(file=trim(sb%filename), exist=found)
          
          if(.not. found) call messages_fatal(1, namespace=namespace)
        end if

        sb%image = gdlib_image_create_from(sb%filename)
        if(.not.c_associated(sb%image)) then
          message(1) = "Could not open file '" // trim(sb%filename) // "' for BoxShape = box_image."
          call messages_fatal(1, namespace=namespace)
        end if
        sb%image_size(1) = gdlib_image_sx(sb%image)
        sb%image_size(2) = gdlib_image_sy(sb%image)

        ! adjust Lsize if necessary to ensure that one grid point = one pixel
        do idir = 1, 2
          box_npts = sb%image_size(idir)
          if((idir >  sb%periodic_dim .and. even(sb%image_size(idir))) .or. &
             (idir <= sb%periodic_dim .and.  odd(sb%image_size(idir)))) then
            box_npts = box_npts + 1
            sb%lsize(idir) = sb%lsize(idir) * box_npts / sb%image_size(idir)
          end if
        end do
#else
        message(1) = "To use 'BoxShape = box_image', you have to compile Octopus"
        message(2) = "with GD library support."
        call messages_fatal(2, namespace=namespace)
#endif
      end if

      ! read in box shape for user-defined boxes
      if(sb%box_shape == BOX_USDEF) then

        !%Variable BoxShapeUsDef
        !%Type string
        !%Section Mesh::Simulation Box
        !%Description
        !% Boolean expression that defines the interior of the simulation box. For example,
        !% <tt>BoxShapeUsDef = "(sqrt(x^2+y^2) <= 4) && z>-2 && z<2"</tt> defines a cylinder
        !% with axis parallel to the <i>z</i>-axis.
        !%End

        call parse_variable(namespace, 'BoxShapeUsDef', 'x^2+y^2+z^2 < 4', sb%user_def)
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
        do idir = 1, sb%dim
          if(sb%rsize > M_ZERO) then
            sb%lsize(idir) = maxval(abs(geo%atom(:)%x(idir))) + sb%rsize
          else
            sb%lsize(idir) = maxval(abs(geo%atom(:)%x(idir))) + def_rsize
          end if
        end do
      end select

      call messages_obsolete_variable(namespace, 'BoxOffset')
      
      POP_SUB(simul_box_init.read_box)
    end subroutine read_box

  end subroutine simul_box_init

  ! ------------------------------------------------------------
  subroutine simul_box_lookup_init(this, geo)
    type(simul_box_t), intent(inout) :: this
    type(geometry_t),  intent(in)    :: geo
    !
    FLOAT, allocatable :: pos(:, :)
    integer            :: iatom
    !
    PUSH_SUB(simul_box_lookup_init)

    SAFE_ALLOCATE(pos(1:this%dim,1:geo%natoms))

    do iatom = 1, geo%natoms
      pos(:,iatom) = geo%atom(iatom)%x(1:this%dim)
    end do

    call lookup_init(this%atom_lookup, this%dim, geo%natoms, pos)

    SAFE_DEALLOCATE_A(pos)
    POP_SUB(simul_box_lookup_init)
    return
  end subroutine simul_box_lookup_init

  !--------------------------------------------------------------
  subroutine simul_box_build_lattice(sb, namespace, rlattice_primitive)
    type(simul_box_t), intent(inout) :: sb
    type(namespace_t), intent(in)    :: namespace
    FLOAT,   optional, intent(in)    :: rlattice_primitive(:,:)

    type(block_t) :: blk
    FLOAT :: norm, lparams(3)
    integer :: idim, jdim, ncols
    logical :: has_angles
    FLOAT :: angles(1:MAX_DIM), cosang, a2, aa, cc
    FLOAT, parameter :: tol_angle = CNST(1.0e-6)  

    PUSH_SUB(simul_box_build_lattice)

    sb%alpha = CNST(90.0)
    sb%beta  = CNST(90.0)
    sb%gamma = CNST(90.0)

    if(present(rlattice_primitive)) then
      sb%rlattice_primitive(1:sb%dim, 1:sb%dim) = rlattice_primitive(1:sb%dim, 1:sb%dim)
      sb%nonorthogonal = .false.
    else
      
      
      !%Variable LatticeParameters
      !%Type block
      !%Default 1 | 1 | 1
      !%Section Mesh::Simulation Box
      !%Description
      !% The lattice parameters (a, b, c). 
      !% This option is incompatible with Lsize and either one of the 
      !% two must be specified in the input file for periodic systems.
      !% A second optional line can be used tu define the angles between the lattice vectors
      !%End
      lparams(:) = M_ONE
      has_angles = .false.
      angles = CNST(90.0)

      if (parse_block(namespace, 'LatticeParameters', blk) == 0) then
        do idim = 1, sb%dim
          call parse_block_float(blk, 0, idim - 1, lparams(idim))
        end do

        if(parse_block_n(blk) > 1) then ! we have a shift, or even more
          ncols = parse_block_cols(blk, 1)
          if(ncols /= sb%dim) then
            write(message(1),'(a,i3,a,i3)') 'LatticeParameters angle has ', ncols, ' columns but must have ', sb%dim
            call messages_fatal(1, namespace=namespace)
          end if
          do idim = 1, sb%dim
            call parse_block_float(blk, 1, idim - 1, angles(idim))
          end do
          has_angles = .true.
        end if
        call parse_block_end(blk)

        if (parse_is_defined(namespace, 'Lsize')) then
          message(1) = 'LatticeParameters is incompatible with Lsize'
          call messages_print_var_info(stdout, "LatticeParameters")
          call messages_fatal(1, namespace=namespace)
        end if 

      end if

      if( has_angles ) then
        sb%alpha = angles(1)
        sb%beta  = angles(2)
        sb%gamma = angles(3)

        !Converting the angles to LatticeVectors
        !See 57_iovars/ingeo.F90 in Abinit for details
        if( abs(angles(1)-angles(2))< tol_angle .and. abs(angles(2)-angles(3))< tol_angle .and.  &
                 (abs(angles(1)-CNST(90.0))+abs(angles(2)-CNST(90.0))+abs(angles(3)-CNST(90.0)))> tol_angle ) then

          cosang=cos(M_PI*angles(1)/CNST(180.0));
          a2=M_TWO/M_THREE*(M_ONE-cosang);
          aa=sqrt(a2);
          cc=sqrt(M_ONE-a2);
          sb%rlattice_primitive(1,1) = aa
          sb%rlattice_primitive(2,1) = M_ZERO
          sb%rlattice_primitive(3,1) = cc
          sb%rlattice_primitive(1,2) =-M_HALF*aa
          sb%rlattice_primitive(2,2) = M_HALF*sqrt(M_THREE)*aa
          sb%rlattice_primitive(3,2) = cc
          sb%rlattice_primitive(1,3) =-M_HALF*aa
          sb%rlattice_primitive(2,3) =-M_HALF*sqrt(M_THREE)*aa
          sb%rlattice_primitive(3,3) = cc
        else
          sb%rlattice_primitive(1,1) = M_ONE
          sb%rlattice_primitive(2,1) = M_ZERO
          sb%rlattice_primitive(3,1) = M_ZERO
          sb%rlattice_primitive(1,2) = cos(M_PI*angles(3)/CNST(180.0))
          sb%rlattice_primitive(2,2) = sin(M_PI*angles(3)/CNST(180.0))
          sb%rlattice_primitive(3,2) = M_ZERO
          sb%rlattice_primitive(1,3) = cos(M_PI*angles(2)/CNST(180.0))
          sb%rlattice_primitive(2,3) = (cos(M_PI*angles(1)/CNST(180.0))-sb%rlattice_primitive(1,2)* sb%rlattice_primitive(1,3))&
                                         /sb%rlattice_primitive(2,2) 
          sb%rlattice_primitive(3,3) = sqrt(M_ONE-sb%rlattice_primitive(1,3)**2-sb%rlattice_primitive(2,3)**2)
        end if

        if (parse_is_defined(namespace, 'LatticeVectors')) then
          message(1) = 'LatticeParameters with angles is incompatible with LatticeVectors'
          call messages_print_var_info(stdout, "LatticeParameters")
          call messages_fatal(1, namespace=namespace)
        end if

        if(any(abs(angles-CNST(90.0)) > M_EPSILON )) then
          sb%nonorthogonal = .true.
        else
          sb%nonorthogonal = .false.
        end if
        
      else  

        !%Variable LatticeVectors
        !%Type block
        !%Default simple cubic
        !%Section Mesh::Simulation Box
        !%Description
        !% Primitive lattice vectors. Vectors are stored in rows.
        !% Default:
        !% <br><br><tt>%LatticeVectors
        !% <br>&nbsp;&nbsp;1.0 | 0.0 | 0.0
        !% <br>&nbsp;&nbsp;0.0 | 1.0 | 0.0
        !% <br>&nbsp;&nbsp;0.0 | 0.0 | 1.0
        !% <br>%<br></tt>
        !%End
        sb%rlattice_primitive = M_ZERO
        sb%nonorthogonal = .false.
        do idim = 1, sb%dim
          sb%rlattice_primitive(idim, idim) = M_ONE
        end do

        if (parse_block(namespace, 'LatticeVectors', blk) == 0) then 
          do idim = 1, sb%dim
            do jdim = 1, sb%dim
              call parse_block_float(blk, idim - 1,  jdim - 1, sb%rlattice_primitive(jdim, idim))
              if(idim /= jdim .and. abs(sb%rlattice_primitive(jdim, idim)) > M_EPSILON) sb%nonorthogonal = .true.
            enddo
          end do
          call parse_block_end(blk)

        end if
      end if

      ! Always need Lsize for periodic systems even if LatticeVectors block is not present
      if (.not. parse_is_defined(namespace, 'Lsize') .and. sb%periodic_dim > 0) then
        do idim = 1, sb%dim
          if (sb%lsize(idim) == M_ZERO) then
            sb%lsize(idim) = lparams(idim)*M_HALF
          end if
        end do
      end if

    end if

    sb%rlattice = M_ZERO
    do idim = 1, sb%dim
      norm = sqrt(sum(sb%rlattice_primitive(1:sb%dim, idim)**2))
      sb%lsize(idim) = sb%lsize(idim) * norm
      do jdim = 1, sb%dim
        sb%rlattice_primitive(jdim, idim) = sb%rlattice_primitive(jdim, idim) / norm
        sb%rlattice(jdim, idim) = sb%rlattice_primitive(jdim, idim) * M_TWO*sb%lsize(idim)
      end do
    end do
    
    call reciprocal_lattice(sb%rlattice, sb%klattice, sb%rcell_volume, sb%dim, namespace)
    sb%klattice = sb%klattice * M_TWO*M_PI

    call reciprocal_lattice(sb%rlattice_primitive, sb%klattice_primitive, sb%volume_element, sb%dim, namespace)

    if(sb%dim == 3) then
      sb%surface_element(1) = sqrt(abs(sum(dcross_product(sb%rlattice_primitive(1:3, 2), sb%rlattice_primitive(1:3, 3))**2)))
      sb%surface_element(2) = sqrt(abs(sum(dcross_product(sb%rlattice_primitive(1:3, 3), sb%rlattice_primitive(1:3, 1))**2)))
      sb%surface_element(3) = sqrt(abs(sum(dcross_product(sb%rlattice_primitive(1:3, 1), sb%rlattice_primitive(1:3, 2))**2)))
    end if

    ! rlattice_primitive is the A matrix from Chelikowski PRB 78 075109 (2008)
    ! klattice_primitive is the transpose (!) of the B matrix, with no 2 pi factor included
    ! klattice is the proper reciprocal lattice vectors, with 2 pi factor, and in units of 1/bohr
    ! The F matrix of Chelikowski is matmul(transpose(sb%klattice_primitive), sb%klattice_primitive)
    sb%rmetric = matmul(transpose(sb%rlattice_primitive), sb%rlattice_primitive)
    if(.not. has_angles .and. sb%dim == 3) then
      !We compute the angles from the lattice vectors
      sb%alpha=acos(sb%rmetric(2,3)/sqrt(sb%rmetric(2,2)*sb%rmetric(3,3)))/M_PI*CNST(180.0)
      sb%beta =acos(sb%rmetric(1,3)/sqrt(sb%rmetric(1,1)*sb%rmetric(3,3)))/M_PI*CNST(180.0)
      sb%gamma=acos(sb%rmetric(1,2)/sqrt(sb%rmetric(1,1)*sb%rmetric(2,2)))/M_PI*CNST(180.0)
    end if

    POP_SUB(simul_box_build_lattice)
  end subroutine simul_box_build_lattice


  !> This function adjusts the coordinates defined in the geometry
  !! object. If coordinates were given in reduced coordinates it
  !! converts them to real coordinates and it checks that the atoms
  !! are inside the box.
  !!
  !! If atoms are not in the box: if the system is periodic, the atoms
  !! are moved inside the box, if the system is finite, nothing
  !! happens or a warning is written, depending on the argument
  !! warn_if_not.
  ! ---------------------------------------------------------
  subroutine simul_box_atoms_in_box(sb, geo, namespace, warn_if_not, die_if_not)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(inout) :: geo
    type(namespace_t), intent(in)    :: namespace
    logical,           intent(in)    :: warn_if_not
    logical, optional, intent(in)    :: die_if_not

    integer :: iatom, pd
    logical :: die_if_not_

    PUSH_SUB(simul_box_atoms_in_box)

    die_if_not_ = optional_default(die_if_not, .false.)
    pd = sb%periodic_dim

    do iatom = 1, geo%natoms

      call simul_box_periodic_atom_in_box(sb, geo, geo%atom(iatom)%x(:))

      if(geo%reduced_coordinates) then
        geo%atom(iatom)%x(pd + 1:sb%dim) = M_TWO*sb%lsize(pd + 1:sb%dim)*geo%atom(iatom)%x(pd + 1:sb%dim)
      end if

      if( .not. sb%contains_point(geo%atom(iatom)%x)) then
        write(message(1), '(a,i5,a)') "Atom ", iatom, " is outside the box." 
        if (sb%periodic_dim /= sb%dim) then
          ! FIXME: This could fail for partial periodicity systems
          ! because contains_point is too strict with atoms close to
          ! the upper boundary to the cell.
          if(warn_if_not) call messages_warning(1, namespace=namespace)
          if(die_if_not_) call messages_fatal(1, namespace=namespace)
        end if
      end if

    end do

    ! done with the conversion to real coordinates
    geo%reduced_coordinates =  .false.

    POP_SUB(simul_box_atoms_in_box)
  end subroutine simul_box_atoms_in_box

  ! --------------------------------------------------------
  
  subroutine simul_box_periodic_atom_in_box(sb, geo, ratom)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(in)    :: geo
    FLOAT,             intent(inout) :: ratom(:)

    FLOAT :: xx(1:MAX_DIM)
    integer :: pd, idir

    pd = sb%periodic_dim

    if (simul_box_is_periodic(sb)) then
      if(.not. geo%reduced_coordinates) then
        !convert the position to reduced coordinates
         xx(1:pd) = matmul(ratom(1:pd), sb%klattice(1:pd, 1:pd))/(M_TWO*M_PI)
      else
        ! in this case coordinates are already in reduced space
        xx(1:pd) = ratom(1:pd)
      end if

      xx(1:pd) = xx(1:pd) + M_HALF
      do idir = 1, pd
        xx(idir) = xx(idir) - anint(xx(idir))
        if(xx(idir) < -CNST(1.0e-6)) &
          xx(idir) = xx(idir) + M_ONE
      end do
      ASSERT(all(xx(1:pd) >= -CNST(1.0e-6)))
      ASSERT(all(xx(1:pd) < CNST(1.0)))

      xx(1:pd) = (xx(1:pd) - M_HALF)
      ratom(1:pd) = matmul(sb%rlattice(1:pd, 1:pd), xx(1:pd))

    end if
    

  end subroutine simul_box_periodic_atom_in_box

  !--------------------------------------------------------------
  subroutine reciprocal_lattice(rv, kv, volume, dim, namespace)
    FLOAT,             intent(in)  :: rv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: kv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: volume
    integer,           intent(in)  :: dim
    type(namespace_t), intent(in)  :: namespace

    integer :: ii
    FLOAT :: cross(1:3), rv3(1:3, 1:3)

    PUSH_SUB(reciprocal_lattice)

    kv(:,:) = M_ZERO

    select case(dim)
    case(3)
      cross(1:3) = dcross_product(rv(1:3, 2), rv(1:3, 3)) 
      volume = dot_product(rv(1:3, 1), cross(1:3))

      kv(1:3, 1) = dcross_product(rv(:, 2), rv(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv(:, 3), rv(:, 1))/volume
      kv(1:3, 3) = dcross_product(rv(:, 1), rv(:, 2))/volume    
    case(2)
      rv3(1:3, 1:3) = M_ZERO
      rv3(1:2, 1:2) = rv(1:2, 1:2)
      rv3(3, 3) = M_ONE
      cross(1:3) = dcross_product(rv3(1:3, 1), rv3(1:3, 2)) 
      volume = dot_product(rv3(1:3, 3), cross(1:3))

      kv(1:3, 1) = dcross_product(rv3(:, 2), rv3(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv3(:, 3), rv3(:, 1))/volume
    case(1)
      volume = rv(1, 1)
      kv(1, 1) = M_ONE / rv(1, 1)
    case default ! dim > 3
      message(1) = "Reciprocal lattice for dim > 3 assumes no periodicity."
      call messages_warning(1, namespace=namespace)
      volume = M_ONE
      do ii = 1, dim
        kv(ii, ii) = M_ONE/rv(ii,ii)
        !  At least initialize the thing
        volume = volume * sqrt(sum(rv(:, ii)**2))
      end do
    end select

    if ( volume < M_ZERO ) then 
      message(1) = "Your lattice vectors form a left-handed system."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(reciprocal_lattice)
  end subroutine reciprocal_lattice

  !--------------------------------------------------------------
  subroutine simul_box_end(sb)
    type(simul_box_t), intent(inout) :: sb    

    PUSH_SUB(simul_box_end)

    call symmetries_end(sb%symm)

    call lookup_end(sb%atom_lookup)
    call kpoints_end(sb%kpoints)

#ifdef HAVE_GDLIB
    if(sb%box_shape == BOX_IMAGE) &
      call gdlib_imagedestroy(sb%image)
#endif

    POP_SUB(simul_box_end)
  end subroutine simul_box_end


  !--------------------------------------------------------------
  recursive subroutine simul_box_write_info(sb, iunit)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: iunit

    character(len=15), parameter :: bs(6) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'minimum       ', &
      'parallelepiped', &
      'image-defined ', &
      'hypercube     '/)

    integer :: idir, idir2, ispec

    PUSH_SUB(simul_box_write_info)

    write(message(1),'(a)') 'Simulation Box:'
    if(sb%box_shape  ==  BOX_USDEF) then
      write(message(2), '(a)') '  Type = user-defined'
    else if(sb%box_shape == BOX_IMAGE) then
      write(message(2), '(3a,i6,a,i6)') '  Type = defined by image "', trim(sb%filename), '", ', &
        sb%image_size(1), ' x ', sb%image_size(2)
    else
      write(message(2), '(2a)') '  Type = ', trim(bs(sb%box_shape))
    end if
    call messages_info(2, iunit)

    if(sb%box_shape == SPHERE .or. sb%box_shape == CYLINDER &
      .or. (sb%box_shape == MINIMUM .and. sb%rsize > M_ZERO)) then
      write(message(1), '(3a,f7.3)') '  Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%rsize)
      call messages_info(1, iunit)
    end if

    if (sb%box_shape == MINIMUM .and. sb%rsize <= M_ZERO) then
      do ispec = 1, sb%geo%nspecies
        write(message(1), '(a,a5,5x,a,f7.3,2a)') '  Species = ', trim(species_label(sb%geo%species(ispec))), 'Radius = ', &
          units_from_atomic(units_out%length, species_def_rsize(sb%geo%species(ispec))), ' ', trim(units_abbrev(units_out%length))
        call messages_info(1, iunit)
      end do
    end if

    if(sb%box_shape == CYLINDER) then
      write(message(1), '(3a,f7.3)') '  Xlength [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%xsize)
      call messages_info(1, iunit)
    end if

    if(sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(3a, 99(a, f8.3), a)')     &
        '  Lengths [', trim(units_abbrev(units_out%length)), '] = ',    &
        '(', (units_from_atomic(units_out%length, M_TWO*sb%lsize(idir)), ',', idir = 1, sb%dim - 1),  &
        units_from_atomic(units_out%length, M_TWO*sb%lsize(sb%dim)), ')'
      call messages_info(1, iunit)
    end if

    write(message(1), '(a,i1,a)') '  Octopus will run in ', sb%dim, ' dimension(s).'
    write(message(2), '(a,i1,a)') '  Octopus will treat the system as periodic in ', &
      sb%periodic_dim, ' dimension(s).'
    call messages_info(2, iunit)

    if(sb%periodic_dim > 0 .or. sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(1x)')
      write(message(2),'(a,3a,a)') '  Lattice Vectors [', trim(units_abbrev(units_out%length)), ']'
      do idir = 1, sb%dim
        write(message(2+idir),'(9f12.6)') (units_from_atomic(units_out%length, sb%rlattice(idir2, idir)), &
          idir2 = 1, sb%dim) 
      end do
      call messages_info(2+sb%dim, iunit)

      write(message(1),'(a,f18.4,3a,i1.1,a)') &
        '  Cell volume = ', units_from_atomic(units_out%length**sb%dim, sb%rcell_volume), &
        ' [', trim(units_abbrev(units_out%length**sb%dim)), ']'
      call messages_info(1, iunit)

      write(message(1),'(a,3a,a)') '  Reciprocal-Lattice Vectors [', trim(units_abbrev(units_out%length**(-1))), ']'
      do idir = 1, sb%dim
        write(message(1+idir),'(3f12.6)') (units_from_atomic(unit_one / units_out%length, sb%klattice(idir2, idir)), &
          idir2 = 1, sb%dim)
      end do
      call messages_info(1+sb%dim, iunit)

      if(sb%dim == 3) then
        write(message(1),'(a)') '  Cell angles [degree]'
        write(message(2),'(a, f8.3)') '    alpha = ', sb%alpha
        write(message(3),'(a, f8.3)') '    beta  = ', sb%beta
        write(message(4),'(a, f8.3)') '    gamma = ', sb%gamma
        call messages_info(4, iunit)
      end if
    end if

    POP_SUB(simul_box_write_info)
  end subroutine simul_box_write_info

  subroutine simul_box_write_short_info(sb, iunit)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: iunit

    integer :: idir1, idir2

    PUSH_SUB(simul_box_write_short_info)

    write(iunit, '(a,i1,a)', advance='no') 'Dimensions = ', sb%dim, '; '
    write(iunit, '(a,i1,a)', advance='no') 'PeriodicDimensions = ', sb%periodic_dim, '; '

    write(iunit, '(a)', advance='no') 'BoxShape = '
    select case(sb%box_shape)
    case(SPHERE)
      write(iunit, '(a,f11.6,a)', advance='no') 'sphere; Radius =', units_from_atomic(unit_angstrom, sb%rsize), ' Ang'
    case(CYLINDER)
      write(iunit, '(a,f11.6,a,f11.6,a)', advance='no') 'cylinder, Radius =', units_from_atomic(unit_angstrom, sb%rsize), &
        ' Ang; Xlength =', units_from_atomic(unit_angstrom, sb%xsize), ' Ang'
    case(MINIMUM)
      write(iunit, '(a,f11.6,a)', advance='no') 'minimum; Radius =', units_from_atomic(unit_angstrom, sb%rsize), ' Ang'
    case(PARALLELEPIPED)
      write(iunit, '(a)', advance='no') 'parallelepiped; LatticeVectors[Ang] = '
      do idir1 = 1, sb%dim
        write(iunit, '(a)', advance='no') '['
        do idir2 = 1, sb%dim
          write(iunit, '(x,f11.6)', advance='no') units_from_atomic(unit_angstrom, sb%rlattice(idir2, idir1))
        end do
        write(iunit, '(a)', advance='no') ']'
        if(idir1 < sb%dim) then
          write(iunit, '(a)', advance='no') ', '
        end if
      end do
    case(BOX_IMAGE)
      write(iunit, '(a)', advance='no') 'box_image; BoxShapeImage = '//trim(sb%filename)
    case(HYPERCUBE)
      write(iunit, '(a)', advance='no') 'hypercube'  ! add parameters?
    case(BOX_USDEF)
      write(iunit, '(a)', advance='no') 'user_defined; BoxShapeUsDef = "'//trim(sb%user_def)//'"'
    end select

    write(iunit, '()')
    POP_SUB(simul_box_write_short_info)

  end subroutine simul_box_write_short_info

  !--------------------------------------------------------------
  !> Checks if a group of mesh points belong to the actual mesh.
  function simul_box_contains_points(this, nn, xx) result(contained)
    class(simul_box_t), intent(in)  :: this
    integer,            intent(in)  :: nn
    FLOAT,              intent(in)  :: xx(:, :) !< (1:, 1:this%dim)
    logical :: contained(1:nn)

    FLOAT, parameter :: DELTA = CNST(1e-12)
    FLOAT :: rr, re, im, dist2, radius
    FLOAT :: llimit(MAX_DIM), ulimit(MAX_DIM)
    FLOAT, allocatable :: xx_red(:, :)
    integer :: ip, idir, iatom, ilist
    integer, allocatable :: nlist(:)
    integer, pointer :: list(:, :)

#if defined(HAVE_GDLIB)
    integer :: red, green, blue, ix, iy
#endif

    ! no push_sub because this function is called very frequently
    SAFE_ALLOCATE(xx_red(1:nn, 1:this%dim))
    xx_red = M_ZERO
    
    !convert from Cartesian to reduced lattice coord 
    if (nn == 1) then
      xx_red(1, 1:this%dim) = matmul(xx(1, 1:this%dim), this%klattice_primitive(1:this%dim, 1:this%dim))
    else
      call lalg_gemm(nn, this%dim, this%dim, M_ONE, xx, this%klattice_primitive, M_ZERO, xx_red)
    end if

    select case(this%box_shape)
    case(SPHERE)
      do ip = 1, nn
        contained(ip) = sum(xx_red(ip, 1:this%dim)**2) <= (this%rsize + DELTA)**2
      end do

    case(CYLINDER)
      do ip = 1, nn
        rr = sqrt(sum(xx_red(ip, 2:this%dim)**2))
        contained(ip) = rr <= this%rsize + DELTA
        if(this%periodic_dim >= 1) then
          contained(ip) = contained(ip) .and. xx_red(ip, 1) >= -this%xsize - DELTA
          contained(ip) = contained(ip) .and. xx_red(ip, 1) <=  this%xsize - DELTA
        else
          contained(ip) = contained(ip) .and. abs(xx_red(ip, 1)) <= this%xsize + DELTA
        end if
      end do

    case(MINIMUM)

      if(this%rsize > M_ZERO) then
        radius = this%rsize
      else
        radius = M_ZERO
        do iatom = 1, this%geo%natoms
          radius = max(radius, species_def_rsize(this%geo%atom(iatom)%species))
        end do
      end if

      radius = radius + DELTA

      SAFE_ALLOCATE(nlist(1:nn))

      if(this%rsize > M_ZERO) then
        nullify(list)
        call lookup_get_list(this%atom_lookup, nn, xx_red, radius, nlist)
      else
        call lookup_get_list(this%atom_lookup, nn, xx_red, radius, nlist, list = list)
      end if

      if(this%rsize > M_ZERO) then
        do ip = 1, nn
          contained(ip) = (nlist(ip) /= 0)
        end do
      else
        do ip = 1, nn
          contained(ip) = .false.
          do ilist = 1, nlist(ip)
            iatom = list(ilist, ip)
            dist2 = sum((xx_red(ip, 1:this%dim) - this%geo%atom(iatom)%x(1:this%dim))**2)
            if(dist2 < species_def_rsize(this%geo%atom(iatom)%species)**2) then
              contained(ip) = .true.
              exit
            end if
          end do
        end do
      end if

      SAFE_DEALLOCATE_A(nlist)
      SAFE_DEALLOCATE_P(list)

    case(PARALLELEPIPED, HYPERCUBE) 
      llimit(1:this%dim) = -this%lsize(1:this%dim) - DELTA
      ulimit(1:this%dim) =  this%lsize(1:this%dim) + DELTA
      ulimit(1:this%periodic_dim)  = this%lsize(1:this%periodic_dim) - DELTA

      do ip = 1, nn
        contained(ip) = all(xx_red(ip, 1:this%dim) >= llimit(1:this%dim) .and. xx_red(ip, 1:this%dim) <= ulimit(1:this%dim))
      end do

#if defined(HAVE_GDLIB)
! Why the minus sign for y? Explanation: http://biolinx.bios.niu.edu/bios546/gd_mod.htm
! For reasons that probably made sense to someone at some time, computer graphic coordinates are not the same
! as in standard graphing. ... The top left corner of the screen is (0,0).

    case(BOX_IMAGE)
      do ip = 1, nn
        ix = nint(( xx_red(ip, 1) + this%lsize(1)) * this%image_size(1) / (M_TWO * this%lsize(1)))
        iy = nint((-xx_red(ip, 2) + this%lsize(2)) * this%image_size(2) / (M_TWO * this%lsize(2)))
        call gdlib_image_get_pixel_rgb(this%image, ix, iy, red, green, blue)
        contained(ip) = (red == 255) .and. (green == 255) .and. (blue == 255)
      end do
#endif

    case(BOX_USDEF)
      ! is it inside the user-given boundaries?
      do ip = 1, nn
        contained(ip) =  all(xx_red(ip, 1:this%dim) >= -this%lsize(1:this%dim) - DELTA) &
          .and. all(xx(ip, 1:this%dim) <= this%lsize(1:this%dim) + DELTA)

        ! and inside the simulation box?
        do idir = 1, this%dim
          xx_red(ip, idir) = units_from_atomic(units_inp%length, xx_red(ip, idir))
        end do
        rr = sqrt(sum(xx_red(ip, 1:this%dim)**2))
        call parse_expression(re, im, this%dim, xx_red(ip, :), rr, M_ZERO, this%user_def)
        contained(ip) = contained(ip) .and. (re /= M_ZERO)
      end do
    end select

    SAFE_DEALLOCATE_A(xx_red)

  end function simul_box_contains_points


  !--------------------------------------------------------------
  logical pure function simul_box_is_periodic(sb)
    type(simul_box_t), intent(in) :: sb

    simul_box_is_periodic = sb%periodic_dim > 0

  end function simul_box_is_periodic


  !--------------------------------------------------------------
  logical pure function simul_box_has_zero_bc(sb)
    type(simul_box_t), intent(in) :: sb

    simul_box_has_zero_bc = .not. simul_box_is_periodic(sb)

  end function simul_box_has_zero_bc

  ! --------------------------------------------------------------
  recursive subroutine simul_box_copy(sbout, sbin)
    type(simul_box_t), intent(out) :: sbout
    type(simul_box_t), intent(in)  :: sbin

    PUSH_SUB(simul_box_copy)

    sbout%box_shape               = sbin%box_shape
    sbout%rsize                   = sbin%rsize
    sbout%xsize                   = sbin%xsize
    sbout%lsize                   = sbin%lsize
    sbout%image                   = sbin%image
    sbout%user_def                = sbin%user_def
    sbout%rlattice                = sbin%rlattice
    sbout%rlattice_primitive      = sbin%rlattice_primitive
    sbout%klattice                = sbin%klattice
    sbout%klattice_primitive      = sbin%klattice_primitive
    sbout%volume_element          = sbin%volume_element
    sbout%dim                     = sbin%dim
    sbout%periodic_dim            = sbin%periodic_dim

    call kpoints_copy(sbin%kpoints, sbout%kpoints)

    call lookup_copy(sbin%atom_lookup, sbout%atom_lookup)

    if(simul_box_is_periodic(sbin)) call symmetries_copy(sbin%symm, sbout%symm)

    POP_SUB(simul_box_copy)
  end subroutine simul_box_copy

  ! -----------------------------------------------------

  subroutine simul_box_check_atoms_are_too_close(geo, sb, namespace)
    type(geometry_t),  intent(in) :: geo
    type(simul_box_t), intent(in) :: sb
    type(namespace_t), intent(in) :: namespace

    FLOAT :: mindist
    FLOAT, parameter :: threshold = CNST(1e-5)

    PUSH_SUB(simul_box_check_atoms_are_too_close)

    if(geo%natoms == 1) then
      POP_SUB(simul_box_check_atoms_are_too_close)
      return
    end if

    mindist = simul_box_min_distance(geo, sb, real_atoms_only = .false.)
    if(mindist < threshold) then
      write(message(1), '(a)') "Some of the atoms seem to sit too close to each other."
      write(message(2), '(a)') "Please review your input files and the output geometry (in 'static/')."
      write(message(3), '(a, f12.6, 1x, a)') "Minimum distance = ", &
        units_from_atomic(units_out%length, mindist), trim(units_abbrev(units_out%length))
      call messages_warning(3, namespace=namespace)

      ! then write out the geometry, whether asked for or not in Output variable
      call io_mkdir(STATIC_DIR, namespace)
      call geometry_write_xyz(geo, trim(STATIC_DIR)//'/geometry', namespace)
    end if

    if(simul_box_min_distance(geo, sb, real_atoms_only = .true.) < threshold) then
      message(1) = "It cannot be correct to run with physical atoms so close."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(simul_box_check_atoms_are_too_close)
  end subroutine simul_box_check_atoms_are_too_close

  ! ---------------------------------------------------------
  FLOAT function simul_box_min_distance(geo, sb, real_atoms_only) result(rmin)
    type(geometry_t),  intent(in) :: geo
    type(simul_box_t), intent(in) :: sb
    logical, optional, intent(in) :: real_atoms_only

    integer :: iatom, jatom, idir
    FLOAT   :: xx(MAX_DIM)
    logical :: real_atoms_only_
    type(species_t), pointer :: species

    PUSH_SUB(simul_box_min_distance)

    real_atoms_only_ = optional_default(real_atoms_only, .false.)

    rmin = huge(rmin)
    do iatom = 1, geo%natoms
      call atom_get_species(geo%atom(iatom), species)
      if(real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
      do jatom = iatom + 1, geo%natoms
        call atom_get_species(geo%atom(iatom), species)
        if(real_atoms_only_ .and. .not. species_represents_real_atom(species)) cycle
        xx(:) = abs(geo%atom(iatom)%x(:) - geo%atom(jatom)%x(:))
        do idir = 1, sb%periodic_dim
          xx(idir) = xx(idir) - M_TWO * sb%lsize(idir) * floor(xx(idir)/(M_TWO * sb%lsize(idir)) + M_HALF)
        end do
        rmin = min(sqrt(sum(xx**2)), rmin)
      end do
    end do

    if(.not. (geo%only_user_def .and. real_atoms_only_)) then
      ! what if the nearest neighbors are periodic images?
      do idir = 1, sb%periodic_dim
        rmin = min(rmin, abs(sb%lsize(idir)))
      end do
    end if

    POP_SUB(simul_box_min_distance)
  end function simul_box_min_distance


    ! ---------------------------------------------------------
  subroutine simul_box_symmetry_check(this, geo, dim, namespace)
    type(simul_box_t),  intent(in) :: this
    type(geometry_t),   intent(in) :: geo
    integer,            intent(in) :: dim
    type(namespace_t),  intent(in) :: namespace

    integer :: iop, iatom, iatom_symm
    FLOAT :: ratom(1:MAX_DIM)

    PUSH_SUB(simul_box_symmetry_check)

    ! We want to use for instance that
    !
    ! \int dr f(Rr) V_iatom(r) \nabla f(R(v)) = R\int dr f(r) V_iatom(R*r) f(r)
    !
    ! and that the operator R should map the position of atom
    ! iatom to the position of some other atom iatom_symm, so that
    !
    ! V_iatom(R*r) = V_iatom_symm(r)
    !
    do iop = 1, symmetries_number(this%symm)
      if(iop == symmetries_identity_index(this%symm)) cycle

      do iatom = 1, geo%natoms
        ratom = M_ZERO
        if(geo%reduced_coordinates) then
          ratom(1:this%dim) = symm_op_apply_red(this%symm%ops(iop), geo%atom(iatom)%x)
        else
          ratom(1:this%dim) = symm_op_apply_cart(this%symm%ops(iop), geo%atom(iatom)%x)
        end if
     
        call simul_box_periodic_atom_in_box(this, geo, ratom)

        ! find iatom_symm
        do iatom_symm = 1, geo%natoms
          if(all(abs(ratom(1:dim) - geo%atom(iatom_symm)%x(1:dim)) < CNST(1.0e-5))) exit
        end do

        if(iatom_symm > geo%natoms) then
          write(message(1),'(a,i6)') 'Internal error: could not find symmetric partner for atom number', iatom
          write(message(2),'(a,i3,a)') 'with symmetry operation number ', iop, '.'
          call messages_fatal(2, namespace=namespace)
        end if

      end do
    end do

    POP_SUB(simul_box_symmetry_check)
  end subroutine simul_box_symmetry_check

end module simul_box_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
