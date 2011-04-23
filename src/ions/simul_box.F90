!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  use c_pointer_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use loct_m
  use lookup_m
  use math_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use space_m
  use species_m
  use string_m
  use symmetries_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    simul_box_t,                &
    simul_box_ob_info_t,        &
    simul_box_init,             &
    simul_box_end,              &
    simul_box_messages_info,       &
    simul_box_is_periodic,      &
    simul_box_has_zero_bc,      &
    simul_box_in_box,           &
    simul_box_in_box_vec,       &
    simul_box_dump,             &
    simul_box_atoms_in_box,     &
    simul_box_copy,             &
    simul_box_complex_boundaries

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4,         &
    BOX_IMAGE      = 5,         &
    HYPERCUBE      = 6,         &
    BOX_USDEF      = 77

  integer, parameter, public :: &
    LEFT      = 1,              & ! Lead indices,
    RIGHT     = 2,              & ! L=1, R=2.
    BOTTOM    = 3,              & ! for 2D open system
    TOP       = 4,              & ! 
    REAR      = 5,              & ! for 3D open system
    FRONT     = 6,              & !
    BEFORE    = 7,              & ! for 4D open system
    AFTER     = 8,              & !
    TRANS_DIR = 1                 ! Transport is in x-direction.

  integer, public :: NLEADS  ! Number of leads.

  ! the lead-names of the open boundaries, maximum 4D
  character(len=6), dimension(2*4), parameter, public :: LEAD_NAME = &
    (/'left  ', 'right ', 'bottom', 'top   ', 'rear  ', 'front ', 'before', 'after '/)

  ! open boundaries stuff
  type simul_box_ob_info_t
    integer             :: ucells         ! Number of additional unit cells.
    character(len=32)   :: dataset        ! Dataset name of the periodic lead calculation.
    character(len=32)   :: restart_dir    ! Directory where to find the lead restart files.
    character(len=32)   :: static_dir     ! Static directory of the lead ground state.
  end type simul_box_ob_info_t

  type, public :: interp_t
    integer          :: nn, order  ! interpolation points and order
    FLOAT,   pointer :: ww(:)      ! weights
    integer, pointer :: posi(:)    ! positions
  end type interp_t


  type, public :: multiresolution_t
    type(interp_t) :: interp          ! interpolation points
    integer        :: num_areas       ! number of multiresolution areas
    integer        :: num_radii       ! number of radii (resolution borders)
    FLOAT, pointer :: radius(:)       ! radius of the high-resolution area
    FLOAT          :: center(MAX_DIM) ! central point
  end type multiresolution_t

  type simul_box_t
    type(symmetries_t) :: symm
    integer  :: box_shape   ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                            ! 4->parallelepiped (orthonormal, up to now).

    FLOAT :: box_offset(MAX_DIM)  ! shifts of the origin in the respective direction

    FLOAT :: rsize          ! the radius of the sphere or of the cylinder
    FLOAT :: xsize          ! the length of the cylinder in the x-direction
    FLOAT :: lsize(MAX_DIM) ! half of the length of the parallelepiped in each direction.

    type(lookup_t)        :: atom_lookup

    type(c_ptr)         :: image    ! for the box defined through an image
    character(len=1024) :: user_def ! for the user-defined box

    logical :: mr_flag                 ! .true. when using multiresolution
    type(multiresolution_t) :: hr_area ! high-resolution areas

    FLOAT :: rlattice_primitive(MAX_DIM,MAX_DIM)   ! lattice primitive vectors
    FLOAT :: rlattice          (MAX_DIM,MAX_DIM)   ! lattice vectors
    FLOAT :: klattice_primitive(MAX_DIM,MAX_DIM)   ! reciprocal-lattice primitive vectors
    FLOAT :: klattice          (MAX_DIM,MAX_DIM)   ! reciprocal-lattice vectors
    FLOAT :: volume_element                      ! the volume element in real space
    FLOAT :: rcell_volume                        ! the volume of the cell in real space

    type(kpoints_t) :: kpoints                   ! the k-points

    FLOAT :: fft_alpha      ! enlargement factor for double box

    integer :: dim
    integer :: periodic_dim
    integer :: transport_dim
#ifdef HAVE_GDLIB
    integer :: image_size(1:2)
#endif
    logical :: complex_boundaries
  end type simul_box_t

  character(len=22), parameter :: dump_tag = '*** simul_box_dump ***'

contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, geo, space, transport_mode, lead_sb, lead_info)
    type(simul_box_t),                   intent(inout) :: sb
    type(geometry_t),                    intent(inout) :: geo
    type(space_t),                       intent(in)    :: space
    logical,                   optional, intent(in)    :: transport_mode
    type(simul_box_t),         optional, intent(inout) :: lead_sb(:)
    type(simul_box_ob_info_t), optional, intent(in)    :: lead_info(:)

    ! some local stuff
    FLOAT :: def_h, def_rsize
    integer :: idir
    logical :: only_gamma_kpoint

    PUSH_SUB(simul_box_init)

    sb%transport_dim = 0
    
    call geometry_grid_defaults(geo, def_h, def_rsize)

    call read_misc()                       ! Miscellaneous stuff.
    call read_box()                        ! Parameters defining the simulation box.
    call sb_lookup_init()
    call read_box_offset()                 ! Parameters defining the offset of the origin.
    if(present(transport_mode)) then
      ASSERT(present(lead_sb) .and. present(lead_info))
      call ob_simul_box_init(sb, transport_mode, lead_sb, space, lead_info, geo)
    end if
    call simul_box_build_lattice(sb)       ! Build lattice vectors.
    call simul_box_atoms_in_box(sb, geo, .true.)   ! Put all the atoms inside the box.

    call symmetries_init(sb%symm, geo, sb%dim, sb%periodic_dim, sb%rlattice, sb%lsize)

    ! we need k-points for periodic systems or for open boundaries
    only_gamma_kpoint = sb%periodic_dim == 0 .and. .not. present(transport_mode)
    call kpoints_init(sb%kpoints, sb%symm, sb%dim, sb%rlattice, sb%klattice, geo, only_gamma_kpoint)

    POP_SUB(simul_box_init)

  contains


    !--------------------------------------------------------------
    subroutine read_misc()

      integer              :: idir, irad, ii
      type(block_t)        :: blk
      FLOAT,   allocatable :: pos(:)

      PUSH_SUB(simul_box_init.read_misc)

      !%Variable DoubleFFTParameter
      !%Type float
      !%Default 2.0
      !%Section Mesh::FFTs
      !%Description
      !% For solving the Poisson equation in Fourier space, and for applying the local potential
      !% in Fourier space, an auxiliary cubic mesh is built. This mesh will be larger than
      !% the circumscribed cube of the usual mesh by a factor <tt>DoubleFFTParameter</tt>. See
      !% the section that refers to Poisson equation, and to the local potential for details
      !% [the default value of two is typically good].
      !%End
      call parse_float(datasets_check('DoubleFFTParameter'), M_TWO, sb%fft_alpha)
      if (sb%fft_alpha < M_ONE .or. sb%fft_alpha > M_THREE ) then
        write(message(1), '(a,f12.5,a)') "Input: '", sb%fft_alpha, &
          "' is not a valid DoubleFFTParameter"
        message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
        call messages_fatal(2)
      end if

      sb%dim = space%dim

      !%Variable PeriodicDimensions
      !%Type integer
      !%Default 0
      !%Section System
      !%Description
      !% Define how many directions are to be considered periodic. It has to be a number
      !% between zero and <tt>Dimensions</tt>.
      !%Option 0
      !% No direction is periodic (molecule).
      !%Option 1
      !% The <i>x</i> direction is periodic (wire, polymer).
      !%Option 2
      !% The <i>x</i> and <i>y</i> directions are periodic (slab).
      !%Option 3
      !% The <i>x</i>, <i>y</i>, and <i>z</i> directions are periodic (bulk).
      !%End
      call parse_integer(datasets_check('PeriodicDimensions'), 0, sb%periodic_dim)
      if ((sb%periodic_dim < 0) .or. (sb%periodic_dim > MAX_DIM) .or. (sb%periodic_dim > sb%dim)) &
        call input_error('PeriodicDimensions')

      !%Variable ComplexBoundaries
      !%Type logical
      !%Default no
      !%Section System
      !%Description
      !% (Experimental) If enabled the system will have complex
      !% boundaries defined by an electrostatic potential. Must be
      !% used with the SETE poisson solver.
      !%End
      call parse_logical(datasets_check('ComplexBoundaries'), .false., sb%complex_boundaries)
      if(sb%complex_boundaries) call messages_experimental("Complex boundaries")

      !%Variable MultiResolutionArea
      !%Type block
      !%Section Mesh
      !%Description
      !% (Experimental) Multiresolution regions are set with this
      !% parameter. The first three numbers define the central
      !% point of the region, and the following ones set
      !% the radii where resolution changes (measured from the
      !% central point).
      !% NOTE: currently, only one area can be set up, and only works in 3D.
      !%End

      if(parse_block(datasets_check('MultiResolutionArea'), blk) == 0) then
        
        call messages_experimental('Multi-resolution')

        if(sb%dim /= 3) call messages_not_implemented('multi-resolution for dim != 3')
        
        ! number of areas
        sb%hr_area%num_areas = parse_block_n(blk)

        ! number of radii
        sb%hr_area%num_radii = parse_block_cols(blk, 0) - sb%dim

        sb%hr_area%center = M_ZERO

        ! the central point
        do idir = 1, sb%dim
           call parse_block_float(blk, 0, idir - 1, sb%hr_area%center(idir))
        end do

        if (sb%hr_area%num_areas /= 1) call input_error('MultiResolutionArea')

        ! the radii
        SAFE_ALLOCATE(sb%hr_area%radius(1:sb%hr_area%num_radii))
        do irad = 1, sb%hr_area%num_radii
          call parse_block_float(blk, 0, sb%dim + irad - 1, sb%hr_area%radius(irad))
          sb%hr_area%radius(irad) = units_to_atomic(units_inp%length, sb%hr_area%radius(irad))
        end do

        ! Create interpolation points (posi) and weights (ww)

        !%Variable MR_InterpolationOrder
        !%Type integer
        !%Default 5
        !%Section Mesh
        !%Description
        !% The interpolation order in multiresolution approach.
        !%End
        call parse_integer(datasets_check('MR_InterpolationOrder'), 5, sb%hr_area%interp%order)
        if(sb%hr_area%interp%order .le. 0) then
          message(1) = "The value for MR_InterpolationOrder must be > 0."
          call messages_fatal(1)
        end if

        sb%hr_area%interp%nn = 2 * sb%hr_area%interp%order
        SAFE_ALLOCATE(pos(1:sb%hr_area%interp%nn))
        SAFE_ALLOCATE(sb%hr_area%interp%ww(1:sb%hr_area%interp%nn))
        SAFE_ALLOCATE(sb%hr_area%interp%posi(1:sb%hr_area%interp%nn))
        do ii = 1, sb%hr_area%interp%order
          sb%hr_area%interp%posi(ii) = 1 + 2 * (ii - 1)
          sb%hr_area%interp%posi(sb%hr_area%interp%order + ii) = -sb%hr_area%interp%posi(ii)
          pos(ii) = sb%hr_area%interp%posi(ii)
          pos(sb%hr_area%interp%order + ii) = -pos(ii)
        end do
        call interpolation_coefficients(sb%hr_area%interp%nn, pos, M_ZERO,sb%hr_area%interp%ww)
        SAFE_DEALLOCATE_A(pos)

        sb%mr_flag = .true.
      else
        nullify(sb%hr_area%radius)
        nullify(sb%hr_area%interp%posi)
        nullify(sb%hr_area%interp%ww)
        sb%mr_flag = .false.
      end if

      POP_SUB(simul_box_init.read_misc)
    end subroutine read_misc


    !--------------------------------------------------------------
    subroutine read_box()
      type(block_t) :: blk

      FLOAT :: default
#if defined(HAVE_GDLIB)        
      character(len=200) :: filename
#endif      

      PUSH_SUB(simul_box_init.read_box)
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
      !% <li><tt>Box_image</tt> is only allowed in 2D.</li>
      !% </ul>
      !%Option sphere 1
      !% The simulation box will be a sphere of radius <tt>Radius</tt>.
      !%Option cylinder 2
      !% The simulation box will be a cylinder with radius <tt>Radius</tt> and height two times
      !% <tt>Xlength</tt>.
      !%Option minimum 3
      !% The simulation box will be constructed by adding spheres created around each
      !% atom (or user-defined potential), of radius <tt>Radius</tt>.
      !%Option parallelepiped 4
      !% The simulation box will be a parallelepiped whose dimensions are taken from
      !% the variable <tt>Lsize</tt>.
      !%Option box_image 5
      !% The simulation box will be defined through an image. White means that the point
      !% is contained in the simulation box, while any other color means that the point is out.
      !%Option user_defined 77
      !% The shape of the simulation box will be read from the variable <tt>BoxShapeUsDef</tt>.
      !%Option hypercube 6
      !% (experimental) The simulation box will be a hypercube or
      !% hyperparallelepiped. This is equivalent to the
      !% <tt>parallelepiped</tt> box but it can work with an arbitrary
      !% number of dimensions.
      !%End

      call parse_integer(datasets_check('BoxShape'), MINIMUM, sb%box_shape)
      if(.not.varinfo_valid_option('BoxShape', sb%box_shape)) call input_error('BoxShape')
      select case(sb%box_shape)
      case(SPHERE, MINIMUM, BOX_IMAGE, BOX_USDEF)
        if(sb%dim > 1 .and. simul_box_is_periodic(sb)) call input_error('BoxShape')
      case(CYLINDER)
        if (sb%dim > 2 .and. &
          ((sb%dim - sb%periodic_dim == 0) .or. (sb%dim - sb%periodic_dim == 1))) call input_error('BoxShape')
      end select

      ! ignore box_shape in 1D
      if(sb%dim == 1 .and. sb%box_shape /= PARALLELEPIPED .and. sb%box_shape /= HYPERCUBE) &
        sb%box_shape = SPHERE

      ! Cannot use images in 1D or 3D
      if(sb%dim /= 2 .and. sb%box_shape == BOX_IMAGE) call input_error('BoxShape')

      if(sb%dim > 3 .and. sb%box_shape /= HYPERCUBE) then
        message(1) = "For more than 3 dimensions, you can only use the hypercubic box."
        call messages_fatal(1)
      end if

      sb%rsize = -M_ONE
      !%Variable Radius
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% Defines the radius for <tt>BoxShape</tt> = <tt>sphere</tt>, <tt>cylinder</tt>, or <tt>minimum</tt>.
      !% Must be a positive number. If not specified, the code will look for a default value in 
      !% the <tt>Species</tt> block, or, if default pseudopotentials are used, the <tt>rsize</tt> column of
      !% <tt>share/PP/defaults</tt>. For <tt>minimum</tt>, a default radius is chosen separately for each species.
      !%End
      select case(sb%box_shape)
      case(SPHERE, CYLINDER)
        call parse_float(datasets_check('Radius'), def_rsize, sb%rsize, units_inp%length)
        if(sb%rsize < M_ZERO) call input_error('radius')
        if(def_rsize>M_ZERO) call messages_check_def(def_rsize, sb%rsize, 'radius')
      case(MINIMUM)
        default=sb%rsize
        call parse_float(datasets_check('radius'), default, sb%rsize, units_inp%length)
        if(sb%rsize < M_ZERO .and. def_rsize < M_ZERO) call input_error('Radius')
      end select

      if(sb%box_shape == CYLINDER) then
        !%Variable Xlength
        !%Type float
        !%Section Mesh::Simulation Box
        !%Description
        !% If <tt>BoxShape</tt> is <tt>cylinder</tt>, the total length of the cylinder is twice <tt>Xlength</tt>.
        !% The default is <tt>Radius</tt>.
        !%End
        if(sb%rsize > M_ZERO) then
          default = sb%rsize
        else
          default = def_rsize
        endif

        call parse_float(datasets_check('xlength'), default, sb%xsize, units_inp%length)
        sb%lsize(1) = sb%xsize
        if(def_rsize > M_ZERO .and. sb%periodic_dim == 0) &
          call messages_check_def(def_rsize, sb%xsize, 'xlength')
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

        if(parse_block(datasets_check('Lsize'), blk) == 0) then
          if(parse_block_cols(blk,0) < sb%dim) call input_error('Lsize')
          do idir = 1, sb%dim
            call parse_block_float(blk, 0, idir - 1, sb%lsize(idir))
          end do
          call parse_block_end(blk)
        else
          call parse_float(datasets_check('Lsize'), -M_ONE, sb%lsize(1))
          if(sb%lsize(1) .eq. -M_ONE) then
            call input_error('Lsize')
          end if
          sb%lsize(1:sb%dim) = sb%lsize(1)
        end if
        sb%lsize = units_to_atomic(units_inp%length, sb%lsize)

        do idir = 1, sb%dim
          if(def_rsize > M_ZERO .and. sb%periodic_dim < idir) &
            call messages_check_def(def_rsize, sb%lsize(idir), 'Lsize')
        end do
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
        call parse_string(datasets_check("BoxShapeImage"), "", filename)
        sb%image = loct_gdimage_create_from(filename)
        if(.not.c_associated(sb%image)) then
          message(1) = "Could not open file '" // filename // "'"
          call messages_fatal(1)
        end if
        sb%image_size(1) = loct_gdImage_SX(sb%image)
        sb%image_size(2) = loct_gdImage_SY(sb%image)
#else
        message(1) = "To use 'BoxShape = box_image', you have to compile Octopus"
        message(2) = "with GD library support."
        call messages_fatal(2)
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
        
        call parse_string(datasets_check("BoxShapeUsDef"), "x^2+y^2+z^2 < 4", sb%user_def)
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

      POP_SUB(simul_box_init.read_box)
    end subroutine read_box


    !--------------------------------------------------------------
    subroutine read_box_offset()
      integer :: idir
      type(block_t) :: blk

      PUSH_SUB(simul_box_init.read_box_offset)
      !%Variable BoxOffset
      !%Type float
      !%Default 0.0
      !%Section Mesh::Simulation Box
      !%Description
      !% Shifts the zero of the simulation box, relative to the atomic coordinates, by a specified vector.
      !% It can be either a float, interpreted as (x,x,x), or a block containing the (x,y,z) value of the zero.
      !% WARNING: This variable does not seem to work correctly!
      !%End
      sb%box_offset = M_ZERO
      if(parse_block(datasets_check('BoxOffset'), blk) == 0) then
        do idir = 1, sb%dim
          call parse_block_float(blk, 0, idir - 1, sb%box_offset(idir), units_inp%length)
        end do
        call parse_block_end(blk)
      else
        call parse_float(datasets_check('BoxOffset'), M_ZERO, sb%box_offset(1), units_inp%length)
        sb%box_offset(1:sb%dim) = sb%box_offset(1)
      end if

      POP_SUB(simul_box_init.read_box_offset)
    end subroutine read_box_offset


    ! ------------------------------------------------------------
    subroutine sb_lookup_init()

      FLOAT, allocatable :: pos(:, :)
      integer :: iatom

      PUSH_SUB(simul_box_init.sb_lookup_init)

      SAFE_ALLOCATE(pos(1:sb%dim, 1:geo%natoms))
     
      do iatom = 1, geo%natoms
        pos(1:sb%dim, iatom) = geo%atom(iatom)%x(1:sb%dim)
      end do
      
      call lookup_init(sb%atom_lookup, sb%dim, geo%natoms, pos)
      
      SAFE_DEALLOCATE_A(pos)
      POP_SUB(simul_box_init.sb_lookup_init)
    end subroutine sb_lookup_init

  end subroutine simul_box_init


  !--------------------------------------------------------------
  subroutine simul_box_build_lattice(sb, rlattice_primitive)
    type(simul_box_t), intent(inout) :: sb
    FLOAT,   optional, intent(in)    :: rlattice_primitive(:,:)

    type(block_t) :: blk
    FLOAT :: norm, cross(1:3)
    integer :: idim, jdim

    PUSH_SUB(simul_box_build_lattice)
    
    if(present(rlattice_primitive)) then
      sb%rlattice_primitive(1:sb%dim, 1:sb%dim) = rlattice_primitive(1:sb%dim, 1:sb%dim)
    else
      !%Variable LatticeVectors
      !%Type block
      !%Default simple cubic
      !%Section Mesh::Simulation Box
      !%Description
      !% Primitive lattice vectors. Vectors are stored in rows.
      !% Note that these vectors will be normalized to 1 after being read.
      !% Default:
      !% <tt>%LatticeVectors
      !% <br>&nbsp;&nbsp;1.0 | 0.0 | 0.0
      !% <br>&nbsp;&nbsp;0.0 | 1.0 | 0.0
      !% <br>&nbsp;&nbsp;0.0 | 0.0 | 1.0
      !% <br>%</tt>
      !%End
      
      sb%rlattice_primitive = M_ZERO
      forall(idim = 1:sb%dim) sb%rlattice_primitive(idim, idim) = M_ONE
      
      if (parse_block(datasets_check('LatticeVectors'), blk) == 0) then 
        do idim = 1, sb%dim
          do jdim = 1, sb%dim
            call parse_block_float(blk, idim - 1,  jdim - 1, sb%rlattice_primitive(jdim, idim))
          end do
        end do
      end if
    end if

    sb%rlattice = M_ZERO
    do idim = 1, sb%dim
      norm = sqrt(sum(sb%rlattice_primitive(1:sb%dim, idim)**2))
      forall(jdim = 1:sb%dim)
        sb%rlattice_primitive(jdim, idim) = sb%rlattice_primitive(jdim, idim) / norm
        sb%rlattice(jdim, idim) = sb%rlattice_primitive(jdim, idim) * M_TWO*sb%lsize(idim)
      end forall
    end do
    
    ! this has to be updated for non-orthogonal grids
    select case(sb%dim)
    case(3)
      cross = dcross_product(sb%rlattice(:,2), sb%rlattice(:,3))
      sb%rcell_volume = sum(sb%rlattice(1:3, 1)*cross(1:3))
    case(2)
      sb%rcell_volume = abs(sb%rlattice(1, 1)*sb%rlattice(2, 2) - sb%rlattice(1, 2)*sb%rlattice(2, 1))
    case(1)
      sb%rcell_volume = abs(sb%rlattice(2, 1) - sb%rlattice(1, 1))
    case default
      sb%rcell_volume = M_ONE
      do idim = 1, sb%dim
        sb%rcell_volume = sb%rcell_volume*abs(sb%rlattice(idim, idim))
      end do
    end select

    call reciprocal_lattice(sb%rlattice_primitive, sb%klattice_primitive, sb%volume_element, sb%dim)

    sb%klattice = M_ZERO
    forall(idim = 1:sb%dim, jdim = 1:sb%dim)
      sb%klattice(jdim, idim) = sb%klattice_primitive(jdim, idim)*M_TWO*M_PI/(M_TWO*sb%lsize(idim))
    end forall

    POP_SUB(simul_box_build_lattice)
  end subroutine simul_box_build_lattice


  !> This function checks that the atoms are inside the box. If not:
  !! if the system is periodic, the atoms are moved inside the box.
  !! if the system is finite, nothing happens or a warning is written,
  !! depending on the argument warn_if_not.
  ! ---------------------------------------------------------
  subroutine simul_box_atoms_in_box(sb, geo, warn_if_not)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(inout) :: geo
    logical,           intent(in)    :: warn_if_not

    integer :: iatom, pd, idir
    FLOAT :: xx(1:MAX_DIM)

    PUSH_SUB(simul_box_atoms_in_box)

    pd = sb%periodic_dim

    do iatom = 1, geo%natoms
      if (simul_box_is_periodic(sb)) then
        !convert the position to the orthogonal space
        xx(1:pd) = matmul(geo%atom(iatom)%x(1:pd) - sb%box_offset(1:pd), sb%klattice_primitive(1:pd, 1:pd))

        xx(1:pd) = xx(1:pd)/(M_TWO*sb%lsize(1:pd))
        xx(1:pd) = xx(1:pd) + M_HALF
        do idir = 1, pd
          if(xx(idir) >= M_ZERO) then
            xx(idir) = xx(idir) - aint(xx(idir))
          else
            xx(idir) = xx(idir) - aint(xx(idir)) + M_ONE
          end if
        end do
        ASSERT(all(xx(1:pd) >= M_ZERO))
        xx(1:pd) = (xx(1:pd) - M_HALF)*M_TWO*sb%lsize(1:pd) 

        geo%atom(iatom)%x(1:pd) = matmul(sb%klattice_primitive(1:pd, 1:pd), xx(1:pd) + sb%box_offset(1:pd))

      end if

      if( .not. simul_box_in_box(sb, geo, geo%atom(iatom)%x) ) then 
        write(message(1), '(a,i5,a)') "Atom ", iatom, " is outside the box." 
        if (sb%periodic_dim == sb%dim) then 
          message(2) = "This is a bug." 
          call messages_fatal(2) 
        else 
          if(warn_if_not) call messages_warning(1) 
        end if
      end if 

    end do

    POP_SUB(simul_box_atoms_in_box)
  end subroutine simul_box_atoms_in_box


  !--------------------------------------------------------------
  subroutine reciprocal_lattice(rv, kv, volume, dim)
    FLOAT,   intent(in)  :: rv(1:MAX_DIM, 1:MAX_DIM)
    FLOAT,   intent(out) :: kv(1:MAX_DIM, 1:MAX_DIM)
    FLOAT,   intent(out) :: volume
    integer, intent(in)  :: dim

    integer :: ii
    FLOAT :: cross(1:3), rv3(1:3, 1:3)

    PUSH_SUB(reciprocal_lattice)
    
    kv(1:MAX_DIM, 1:MAX_DIM) = M_ZERO

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
      message(1) = "Reciprocal lattice is not correct for dim > 3."
      call messages_warning(1)
      volume = M_ONE
      do ii = 1, dim
        kv(ii, ii) = M_ONE
        !  At least initialize the thing
        volume = volume * sqrt(sum(rv(:, ii)**2))
      end do
    end select

    if ( volume < M_ZERO ) then 
      message(1) = "Your lattice vectors form a left-handed system."
      call messages_fatal(1)
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

    SAFE_DEALLOCATE_P(sb%hr_area%radius)
    SAFE_DEALLOCATE_P(sb%hr_area%interp%ww)
    SAFE_DEALLOCATE_P(sb%hr_area%interp%posi)

    POP_SUB(simul_box_end)
  end subroutine simul_box_end


  !--------------------------------------------------------------
  recursive subroutine simul_box_messages_info(sb, geo, iunit)
    type(simul_box_t), intent(in) :: sb
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iunit

    character(len=15), parameter :: bs(6) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'minimum       ', &
      'parallelepiped', &
      'image defined ', &
      'hypercube     '/)

    integer :: idir, idir2, ispec

    PUSH_SUB(simul_box_messages_info)

    write(message(1),'(a)') 'Simulation Box:'
    if(sb%box_shape .eq. BOX_USDEF) then
      write(message(2), '(a)') '  Type = user-defined'
    else
      write(message(2), '(a,a,1x)') '  Type = ', bs(sb%box_shape)
    end if
    call messages_info(2, iunit)

    if(sb%box_shape == SPHERE .or. sb%box_shape == CYLINDER &
       .or. (sb%box_shape == MINIMUM .and. sb%rsize > M_ZERO)) then
      write(message(1), '(3a,f7.3)') '  Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%rsize)
      call messages_info(1, iunit)
    endif

    if (sb%box_shape == MINIMUM .and. sb%rsize <= M_ZERO) then
      do ispec = 1, geo%nspecies     
        write(message(1), '(a,a5,5x,a,f7.3,2a)') '  Species = ', trim(species_label(geo%species(ispec))), 'Radius = ', &
          units_from_atomic(units_out%length, species_def_rsize(geo%species(ispec))), ' ', trim(units_abbrev(units_out%length))
        call messages_info(1, iunit)
      enddo
    end if

    if(sb%box_shape == CYLINDER) then
      write(message(1), '(3a,f7.3)') '  Xlength [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%xsize)
      call messages_info(1, iunit)
    end if

    if(sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(3a, 99(a, f8.3), a)')     &
        '  Lengths [', trim(units_abbrev(units_out%length)), '] = ',    &
        '(', (units_from_atomic(units_out%length, sb%lsize(idir)), ',', idir = 1, sb%dim - 1),  &
        units_from_atomic(units_out%length, sb%lsize(sb%dim)), ')'
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
    end if

    POP_SUB(simul_box_messages_info)
  end subroutine simul_box_messages_info


  !--------------------------------------------------------------
  logical function simul_box_in_box(sb, geo, yy) result(in_box)
    type(simul_box_t),  intent(in) :: sb
    type(geometry_t),   intent(in) :: geo
    FLOAT,              intent(in) :: yy(:)

    real(8), parameter :: DELTA = CNST(1e-12)
    FLOAT :: xx(1:MAX_DIM, 1)
    logical :: in_box2(1)

    ! no push_sub because this function is called very frequently

    xx(1:sb%dim, 1) = yy(1:sb%dim)

    call simul_box_in_box_vec(sb, geo, 1, xx, in_box2)
    in_box = in_box2(1)
    
  end function simul_box_in_box


  !--------------------------------------------------------------
  subroutine simul_box_in_box_vec(sb, geo, npoints, point, in_box)
    type(simul_box_t),  intent(in)  :: sb
    type(geometry_t),   intent(in)  :: geo
    integer,            intent(in)  :: npoints
    FLOAT,              intent(in)  :: point(:, :)
    logical,            intent(out) :: in_box(:)

    real(8), parameter :: DELTA = CNST(1e-12)
    FLOAT :: rr, re, im, dist2, radius
    real(8) :: llimit(MAX_DIM), ulimit(MAX_DIM)
    FLOAT, allocatable :: xx(:, :)
    integer :: ip, idir, iatom, ilist
    integer, allocatable :: nlist(:)
    integer, pointer :: list(:, :)

#if defined(HAVE_GDLIB)
    integer :: red, green, blue, ix, iy
#endif

    ! no push_sub because this function is called very frequently
    SAFE_ALLOCATE(xx(1:sb%dim, 1:npoints))

    forall(idir = 1:sb%dim, ip = 1:npoints)  xx(idir, ip) = point(idir, ip) - sb%box_offset(idir)

    !convert to the orthogonal space
    forall(ip = 1:npoints)
      xx(1:sb%dim, ip) = matmul(xx(1:sb%dim, ip), sb%klattice_primitive(1:sb%dim, 1:sb%dim))
    end forall

    select case(sb%box_shape)
      case(SPHERE)
        forall(ip = 1:npoints)
          in_box(ip) = sum(xx(1:sb%dim, ip)**2) <= (sb%rsize+DELTA)**2
        end forall

      case(CYLINDER)
        do ip = 1, npoints
          rr = sqrt(xx(2, ip)**2 + xx(3, ip)**2)
          in_box(ip) = (rr <= sb%rsize + DELTA .and. abs(xx(1, ip)) <= sb%xsize + DELTA)
        end do

      case(MINIMUM)

        if(sb%rsize > M_ZERO) then
          radius = sb%rsize
        else
          radius = M_ZERO
          do iatom = 1, geo%natoms
            radius = max(radius, species_def_rsize(geo%atom(iatom)%spec))
          end do
        end if

        radius = radius + DELTA

        SAFE_ALLOCATE(nlist(1:npoints))

        if(sb%rsize > M_ZERO) then
          nullify(list)
          call lookup_get_list(sb%atom_lookup, npoints, xx, radius, nlist)
        else
          call lookup_get_list(sb%atom_lookup, npoints, xx, radius, nlist, list = list)
        end if

        if(sb%rsize > M_ZERO) then
          do ip = 1, npoints
            in_box(ip) = (nlist(ip) /= 0)
          end do
        else
          do ip = 1, npoints
            in_box(ip) = .false.

            do ilist = 1, nlist(ip)

              iatom = list(ilist, ip)

              dist2 = sum((xx(1:sb%dim, ip) - geo%atom(iatom)%x(1:sb%dim))**2)

              if(dist2 < species_def_rsize(geo%atom(iatom)%spec)**2) then
                in_box(ip) = .true.
                exit
              end if

            end do
          end do
        end if

        SAFE_DEALLOCATE_A(nlist)
        SAFE_DEALLOCATE_P(list)

      case(PARALLELEPIPED, HYPERCUBE) 
        llimit(1:sb%dim) = -sb%lsize(1:sb%dim) - DELTA
        ulimit(1:sb%dim) =  sb%lsize(1:sb%dim) + DELTA
        ulimit(1:sb%periodic_dim)  = sb%lsize(1:sb%periodic_dim) - DELTA
        ulimit(1:sb%transport_dim) = sb%lsize(1:sb%transport_dim) - DELTA

        forall(ip = 1:npoints)
          in_box(ip) = all(xx(1:sb%dim, ip) >= llimit(1:sb%dim) .and. xx(1:sb%dim, ip) <= ulimit(1:sb%dim))
        end forall

#if defined(HAVE_GDLIB)
      case(BOX_IMAGE)
        do ip = 1, npoints
          ix = nint((xx(1, ip) + sb%lsize(1)) * sb%image_size(1) / (M_TWO * sb%lsize(1)))
          iy = nint((xx(2, ip) + sb%lsize(2)) * sb%image_size(2) / (M_TWO * sb%lsize(2)))
          call loct_gdimage_get_pixel_rgb(sb%image, ix, iy, red, green, blue)
          in_box(ip) = (red == 255) .and. (green == 255) .and. (blue == 255)
        end do
#endif

      case(BOX_USDEF)
        ! is it inside the user-given boundaries?
        do ip = 1, npoints
          in_box(ip) =  all(xx(1:sb%dim, ip) >= -sb%lsize(1:sb%dim) - DELTA) &
            .and. all(xx(1:sb%dim, ip) <= sb%lsize(1:sb%dim) + DELTA)

          ! and inside the simulation box?
          do idir = 1, sb%dim
            xx(idir, ip) = units_from_atomic(units_inp%length, xx(idir, ip))
          enddo
          rr = sqrt(sum(xx(1:sb%dim, ip)**2))
          call parse_expression(re, im, sb%dim, xx(:, ip), rr, M_ZERO, sb%user_def)
          in_box(ip) = in_box(ip) .and. (re .ne. M_ZERO)
        end do
    end select

    SAFE_DEALLOCATE_A(xx)

  end subroutine simul_box_in_box_vec


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

  
  !--------------------------------------------------------------
  subroutine simul_box_dump(sb, iunit)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: iunit

    integer :: idir

    PUSH_SUB(simul_box_dump)

    write(iunit, '(a)')             dump_tag
    write(iunit, '(a20,i4)')        'box_shape=          ', sb%box_shape
    write(iunit, '(a20,i4)')        'dim=                ', sb%dim
    write(iunit, '(a20,i4)')        'periodic_dim=       ', sb%periodic_dim
    write(iunit, '(a20,i4)')        'transport_dim=      ', sb%transport_dim
    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,99e22.14)') 'lsize=              ', sb%lsize(1:sb%dim)
    case(CYLINDER)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,e22.14)')  'xlength=            ', sb%xsize
      write(iunit, '(a20,99e22.14)') 'lsize=              ', sb%lsize(1:sb%dim)
    case(PARALLELEPIPED)
      write(iunit, '(a20,99e22.14)') 'lsize=              ', sb%lsize(1:sb%dim)
    case(BOX_USDEF)
      write(iunit, '(a20,99e22.14)') 'lsize=              ', sb%lsize(1:sb%dim)
      write(iunit, '(a20,a1024)')   'user_def=           ', sb%user_def
    end select
    write(iunit, '(a20,e22.14)')    'fft_alpha=          ', sb%fft_alpha
    write(iunit, '(a20,99e22.14)')   'box_offset=         ', sb%box_offset(1:sb%dim)
    write(iunit, '(a20,l7)')        'mr_flag=            ', sb%mr_flag
    if(sb%mr_flag) then
      write(iunit, '(a20,i4)')        'num_areas=         ',sb%hr_area%num_areas
      write(iunit, '(a20,i4)')        'num_radii=         ',sb%hr_area%num_radii
      do idir = 1, sb%hr_area%num_radii
        write(iunit, '(a10,i2.2,a9,e22.14)') 'mr_radius_', idir, '=        ',sb%hr_area%radius(idir)
      end do
      do idir = 1, sb%dim
        write(iunit, '(a7,i1,a13,e22.14)')   'center(', idir, ')=           ',sb%hr_area%center(idir)
      end do
    end if
    do idir = 1, sb%dim
      write(iunit, '(a9,i1,a11,99e22.14)')    'rlattice(', idir, ')=         ', &
        sb%rlattice_primitive(1:sb%dim, idir)
    end do

    POP_SUB(simul_box_dump)
  end subroutine simul_box_dump


  ! --------------------------------------------------------------
  subroutine simul_box_init_from_file(sb, iunit)
    type(simul_box_t), intent(inout) :: sb
    integer,           intent(in)    :: iunit

    character(len=20)  :: str
    character(len=300) :: line
    integer            :: idim, il, ierr
    FLOAT              :: rlattice_primitive(1:MAX_DIM, 1:MAX_DIM)

    PUSH_SUB(simul_box_init_from_file)

    ! Find (and throw away) the dump tag.
    do
      call iopar_read(mpi_world, iunit, line, ierr)
      if(trim(line) .eq. dump_tag) exit
    end do

    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%box_shape
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%dim
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%periodic_dim
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%transport_dim

    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%rsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:sb%dim)       
    case(CYLINDER)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%rsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%xsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:sb%dim)
    case(PARALLELEPIPED)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:sb%dim)
    case(BOX_USDEF)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:sb%dim)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%user_def
    end select
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%fft_alpha
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%box_offset(1:sb%dim)
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line,*) str, sb%mr_flag
    if(sb%mr_flag) then
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line,*) str, sb%hr_area%num_areas
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line,*) str, sb%hr_area%num_radii
      SAFE_ALLOCATE(sb%hr_area%radius(1:sb%hr_area%num_radii))
      do il = 1, sb%hr_area%num_radii
        call iopar_read(mpi_world, iunit, line, ierr)
        read(line,*) str, sb%hr_area%radius(il)
      end do
      do idim = 1, sb%dim
        call iopar_read(mpi_world, iunit, line, ierr)
        read(line, *) str, sb%hr_area%center(idim)
      end do
    end if
    do idim = 1, sb%dim
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, rlattice_primitive(1:sb%dim, idim)
    end do

    call simul_box_build_lattice(sb, rlattice_primitive)

    call iopar_read(mpi_world, iunit, line, ierr)

    POP_SUB(simul_box_init_from_file)
  end subroutine simul_box_init_from_file


  ! --------------------------------------------------------------
  recursive subroutine simul_box_copy(sbout, sbin)
    type(simul_box_t), intent(out) :: sbout
    type(simul_box_t), intent(in)  :: sbin

    PUSH_SUB(simul_box_copy)

    sbout%box_shape               = sbin%box_shape
    sbout%box_offset              = sbin%box_offset
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
    sbout%fft_alpha               = sbin%fft_alpha
    sbout%dim                     = sbin%dim
    sbout%periodic_dim            = sbin%periodic_dim
    sbout%transport_dim           = sbin%transport_dim
    sbout%mr_flag                 = sbin%mr_flag
    sbout%hr_area%num_areas       = sbin%hr_area%num_areas
    sbout%hr_area%num_radii       = sbin%hr_area%num_radii
    sbout%hr_area%center(1:sbin%dim)=sbin%hr_area%center(1:sbin%dim)
    
    call kpoints_copy(sbin%kpoints, sbout%kpoints)

    if(sbout%mr_flag) then
      SAFE_ALLOCATE(sbout%hr_area%radius(1:sbout%hr_area%num_radii))
      sbout%hr_area%radius(1:sbout%hr_area%num_radii) = sbin%hr_area%radius(1:sbout%hr_area%num_radii)
    end if

    call lookup_copy(sbin%atom_lookup, sbout%atom_lookup)

    if(simul_box_is_periodic(sbin)) call symmetries_copy(sbin%symm, sbout%symm)

    POP_SUB(simul_box_copy)
  end subroutine simul_box_copy


  !--------------------------------------------------------------
  subroutine ob_simul_box_init(sb, transport_mode, lead_sb, space, lead_info, geo)
    type(simul_box_t), intent(inout) :: sb
    logical,           intent(in)    :: transport_mode
    type(simul_box_t), intent(inout) :: lead_sb(:)
    type(space_t),     intent(in)    :: space
    type(simul_box_ob_info_t), intent(in) :: lead_info(:)
    type(geometry_t),  intent(inout) :: geo

    ! some local stuff
    integer :: il

    PUSH_SUB(ob_simul_box_init)

    ! Open boundaries are only possible for rectangular simulation boxes.
    if(sb%box_shape.ne.PARALLELEPIPED) then
      message(1) = 'Open boundaries are only possible with a parallelepiped'
      message(2) = 'simulation box.'
      call messages_fatal(2)
    end if
    ! Simulation box must not be periodic in transport direction.
    if(sb%periodic_dim.eq.1) then
      message(1) = 'When using open boundaries, you cannot use periodic boundary'
      message(2) = 'conditions in the x-direction.'
      call messages_fatal(2)
    end if

    if(transport_mode) then
      ! lowest index must be transport direction
      ASSERT(TRANS_DIR.eq.1)
      sb%transport_dim = TRANS_DIR
      lead_sb(:)%transport_dim = TRANS_DIR
    else ! just open boundaries
      sb%transport_dim = 0
      lead_sb(:)%transport_dim = 0
    end if

    call ob_read_lead_unit_cells(sb, lead_sb, lead_info(:)%restart_dir)
    ! Adjust the size of the simulation box by adding the proper number
    ! of unit cells to the simulation region.
    do il = LEFT, RIGHT
      sb%lsize(TRANS_DIR) = sb%lsize(TRANS_DIR) + lead_info(il)%ucells*lead_sb(il)%lsize(TRANS_DIR)
    end do
    ! Add the atoms of the lead unit cells that are included in the simulation box to geo.
    call ob_simul_box_add_lead_atoms(sb, lead_sb, space, lead_info(:)%ucells, lead_info(:)%dataset, geo)

    POP_SUB(ob_simul_box_init)

  end subroutine ob_simul_box_init


  !--------------------------------------------------------------
  ! Read the simulation boxes of the leads
  subroutine ob_read_lead_unit_cells(sb, lead_sb, dir)
    type(simul_box_t), intent(inout) :: sb
    type(simul_box_t), intent(inout) :: lead_sb(:)
    character(len=*),  intent(in)    :: dir(:)

    integer :: iunit, il

    PUSH_SUB(ob_read_lead_unit_cells)

    do il = 1, NLEADS
      iunit = io_open(trim(dir(il))//'/'//GS_DIR//'mesh', action = 'read', is_tmp = .true., grp = mpi_world)
      call simul_box_init_from_file(lead_sb(il), iunit)
      call io_close(iunit)

      ! Check whether
      ! * simulation box is a parallelepiped,
      ! * the extensions in y-, z-directions fit the central box,
      ! * the central simulation box x-length is an integer multiple of
      !   the unit cell x-length,
      ! * periodic in one dimension, and
      ! * of the same dimensionality as the central system.

      if(lead_sb(il)%box_shape .ne. PARALLELEPIPED) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' lead is not a parallelepiped.'
        call messages_fatal(1)
      end if

      if(any(sb%lsize(2:sb%dim) .ne. lead_sb(il)%lsize(2:sb%dim))) then
        message(1) = 'The size in non-transport-directions of the ' // LEAD_NAME(il) // ' lead'
        message(2) = 'does not fit the size of the non-transport-directions of the central system.'
        call messages_fatal(2)
      end if

      if(.not. is_integer_multiple(sb%lsize(1), lead_sb(il)%lsize(1))) then
        message(1) = 'The length in x-direction of the central simulation'
        message(2) = 'box is not an integer multiple of the x-length of'
        message(3) = 'the ' // trim(LEAD_NAME(il)) // ' lead.'
        call messages_fatal(3)
      end if

      if(lead_sb(il)%periodic_dim .ne. 1) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' lead is not periodic in x-direction.'
        message(2) = 'For now we assume the first unit cell to be the periodic representative.'
        call messages_warning(2)
      end if
      if(lead_sb(il)%dim .ne. sb%dim) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' has a different dimension than'
        message(2) = 'the central system.'
        call messages_fatal(2)
      end if
    end do

    POP_SUB(ob_read_lead_unit_cells)
  end subroutine ob_read_lead_unit_cells


  !--------------------------------------------------------------
  ! Read the coordinates of the leads atoms and add them to the
  ! simulation box
  subroutine ob_simul_box_add_lead_atoms(sb, lead_sb, space, ucells, lead_dataset, geo)
    type(simul_box_t), intent(inout) :: sb
    type(simul_box_t), intent(inout) :: lead_sb(:)
    type(space_t),     intent(in)    :: space
    integer,           intent(in)    :: ucells(:)
    character(len=32), intent(in)    :: lead_dataset(:)
    type(geometry_t),  intent(inout) :: geo

    type(geometry_t)  :: central_geo
    type(geometry_t), allocatable  :: lead_geo(:)
    character(len=32) :: label_bak
    integer           :: il, icell, iatom, jatom, icatom, dir

    PUSH_SUB(ob_simul_box_add_lead_atoms)

    SAFE_ALLOCATE(lead_geo(1:NLEADS))
    do il = 1, NLEADS
      ! We temporarily change the current label to read the
      ! coordinates of another dataset, namely the lead dataset.
      label_bak     = current_label
      current_label = lead_dataset(il)
      call geometry_init(lead_geo(il), space, print_info=.false.)
      current_label = label_bak
      call simul_box_atoms_in_box(lead_sb(il), lead_geo(il), .true.)
    end do

    ! Merge the geometries of the lead and of the central region.
    call geometry_copy(central_geo, geo)

    ! Set the number of atoms and classical atoms to the number
    ! of atoms coming from left and right lead and central part.
    if(geo%natoms .gt. 0) then
      SAFE_DEALLOCATE_P(geo%atom)
    end if
    geo%natoms = central_geo%natoms + ucells(LEFT)*lead_geo(LEFT)%natoms + ucells(RIGHT)*lead_geo(RIGHT)%natoms
    SAFE_ALLOCATE(geo%atom(1:geo%natoms))
    if(geo%ncatoms .gt. 0) then
      SAFE_DEALLOCATE_P(geo%catom)
    end if
    geo%ncatoms = central_geo%ncatoms + ucells(LEFT)*lead_geo(LEFT)%ncatoms + ucells(RIGHT)*lead_geo(RIGHT)%ncatoms
    SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))

    geo%only_user_def = central_geo%only_user_def .and. all(lead_geo(:)%only_user_def)
    geo%nlpp          = central_geo%nlpp .or. any(lead_geo(:)%nlpp)
    geo%nlcc          = central_geo%nlcc .or. any(lead_geo(:)%nlcc)

    ! FIXME:
    ! Please do not do these atrocities... even if the object has
    ! public components it does not mean you can just initialize it by
    ! hand. XA
    geo%atoms_dist%start   = 1
    geo%atoms_dist%end     = geo%natoms
    geo%atoms_dist%nlocal  = geo%natoms

    ! 1. Put the atoms of the central region into geo.
    geo%atom(1:central_geo%natoms)   = central_geo%atom
    geo%catom(1:central_geo%ncatoms) = central_geo%catom

    ! 2. Put the atoms of the leads into geo and adjust their x-coordinates.
    iatom  = central_geo%natoms + 1
    icatom = central_geo%ncatoms + 1

    do il = 1, NLEADS
      dir = (-1)**il
      ! We start from the "outer" unit cells of the lead.
      do icell = 1, ucells(il)
        do jatom = 1, lead_geo(il)%natoms
          geo%atom(iatom) = lead_geo(il)%atom(jatom)
          geo%atom(iatom)%x(TRANS_DIR) = geo%atom(iatom)%x(TRANS_DIR) + &
            dir * (sb%lsize(TRANS_DIR) - (2*(icell - 1) + 1) * lead_sb(il)%lsize(TRANS_DIR))
          iatom = iatom + 1
        end do

        do jatom = 1, lead_geo(il)%ncatoms
          geo%catom(icatom) = lead_geo(il)%catom(jatom)
          geo%catom(icatom)%x(TRANS_DIR) = geo%catom(icatom)%x(TRANS_DIR) + &
            dir * (sb%lsize(TRANS_DIR) - (2 * (icell - 1) + 1) * lead_sb(il)%lsize(TRANS_DIR))
        end do
      end do
    end do

    ! Initialize the species of the "extended" central system.
    if(geo%nspecies .gt. 0) then
      SAFE_DEALLOCATE_P(geo%species)
    end if
    call geometry_init_species(geo, print_info=.false.)

    do il = 1, NLEADS
      call geometry_end(lead_geo(il))
    end do

    call geometry_end(central_geo)
    SAFE_DEALLOCATE_A(lead_geo)

    POP_SUB(ob_simul_box_add_lead_atoms)
  end subroutine ob_simul_box_add_lead_atoms

  ! -----------------------------------------------------

  logical pure function simul_box_complex_boundaries(sb) result(cb)
    type(simul_box_t),  intent(in) :: sb

    cb = sb%complex_boundaries
  end function simul_box_complex_boundaries
end module simul_box_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
