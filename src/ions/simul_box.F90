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
  use calc_mode_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use lalg_basic_m
  use loct_m
  use lookup_m
  use parser_m
  use kpoints_m
  use math_m
  use messages_m
  use mpi_m
  use profiling_m
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
    simul_box_init,             &
    simul_box_end,              &
    simul_box_write_info,       &
    simul_box_is_periodic,      &
    simul_box_has_zero_bc,      &
    simul_box_in_box,           &
    simul_box_in_box_vec,       &
    simul_box_dump,             &
    simul_box_atoms_in_box,     &
    lead_unit_cell_extent,      &
    simul_box_copy

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
    OUTER     = 1,              & ! Block indices of interface wavefunctions
    INNER     = 2,              & ! for source term.
    TRANS_DIR = 1                 ! Transport is in x-direction.

  integer, public :: NLEADS  ! Number of leads.

  ! the lead-names of the open boundaries, maximum 4D
  character(len=6), dimension(2*4), parameter, public :: LEAD_NAME = &
    (/'left  ', 'right ', 'bottom', 'top   ', 'rear  ', 'front ', 'before', 'after '/)


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

    FLOAT :: h(MAX_DIM)     ! the (canonical) spacing between the points
    FLOAT :: box_offset(MAX_DIM)  ! shifts of the origin in the respective direction

    FLOAT :: rsize          ! the radius of the sphere or of the cylinder
    FLOAT :: xsize          ! the length of the cylinder in the x-direction
    FLOAT :: lsize(MAX_DIM) ! half of the length of the parallelepiped in each direction.

    FLOAT                 :: cell_length(1:3)
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

    ! For open boundaries, we need reference to the lead`s unit cell.
    ! This unit cell is itself a simulation box.
    logical           :: open_boundaries             ! Use open boundaries?
    integer           :: n_ucells                    ! Number of unit cells that fit in central region.
    integer           :: add_unit_cells(2*MAX_DIM)   ! Number of additional unit cells.
    character(len=32) :: lead_dataset(2*MAX_DIM)     ! Dataset name of the periodic lead calculation.
    character(len=32) :: lead_restart_dir(2*MAX_DIM) ! Directory where to find the lead restart files.
    character(len=32) :: lead_static_dir(2*MAX_DIM)  ! Static directory of the lead ground state.
    type(simul_box_t), pointer :: lead_unit_cell(:)  ! Simulation box of the unit cells.
    ! The next one does not really belong here but since the parsing happens in the simul_box_init
    ! it makes things a bit easier.
    character(len=20000) :: lead_td_pot_formula(2*MAX_DIM) ! Td-potential of lead.
  end type simul_box_t

  character(len=22), parameter :: dump_tag = '*** simul_box_dump ***'

contains

  !--------------------------------------------------------------
  subroutine simul_box_init(sb, geo)
    type(simul_box_t), intent(inout) :: sb
    type(geometry_t),  intent(inout) :: geo

    ! some local stuff
    FLOAT :: def_h, def_rsize
    integer :: i

    call push_sub('simul_box.simul_box_init')

    call geometry_grid_defaults(geo, def_h, def_rsize)

    call read_misc()                       ! Miscellaneous stuff.
    call read_box()                        ! Parameters defining the simulation box.
    call sb_lookup_init()
    call read_spacing ()                   ! Parameters defining the (canonical) spacing.
    call read_box_offset()                 ! Parameters defining the offset of the origin.
    call simul_box_build_lattice(sb)       ! Build lattice vectors.
    call read_open_boundaries()            ! Parameters defining open boundaries.
    call simul_box_add_lead_atoms(sb, geo) ! Add the atoms of the lead unit cells that are
                                           ! included in the simulation box to geo.
    call simul_box_atoms_in_box(sb, geo)   ! Put all the atoms inside the box.

    if(simul_box_is_periodic(sb)) call symmetries_init(sb%symm, geo, sb%dim, sb%rlattice, sb%lsize)

    call kpoints_init(sb%kpoints, sb%dim, sb%periodic_dim, sb%rlattice, sb%klattice, geo)

    call pop_sub()

  contains

    !--------------------------------------------------------------
    subroutine read_open_boundaries()
      integer            :: nr, tag, nrows, ncols
      type(block_t)      :: blk

      integer, parameter ::   &
        LEAD_DATASET     = 1, &
        LEAD_RESTART_DIR = 2, &
        LEAD_STATIC_DIR  = 3, &
        ADD_UNIT_CELLS   = 4, &
        TD_POT_FORMULA   = 5

      integer :: il

      call push_sub('simul_box.simul_box_init.read_open_boundaries')

      !%Variable OpenBoundaries
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% Enables open boundaries in the <i>x</i>-direction and defines the character
      !% of the leads attached to the left and right of the finite central system.
      !%
      !% The format is as follows:
      !%
      !% <pre>
      !% %OpenBoundaries
      !%  lead_dataset     | "dataset"   | "dataset"
      !%  lead_restart_dir | "directory" | "directory"
      !%  lead_static_dir  | "directory" | "directory"
      !%  add_unit_cells   | nl          | nr
      !%  td_pot_formula   | "formula"   | "formula"
      !% %
      !% </pre>
      !%
      !% The left column specifies characteristics of the left lead and
      !% and the right column characteristics of the right lead analogously.
      !% If only one column is given, the value specified is used for both leads.
      !%
      !% All entries except <tt>lead_dataset</tt> are optional.
      !%
      !% Currently available only in development version.
      !%
      !%Option lead_dataset 1
      !% Gives the name of the dataset used for the periodic calculation of the
      !% ground states of the leads. It is used, <i>e.g.</i>, to read in the coordinates of the
      !% atoms of the lead. Both entries for left and right have to be equal.
      !%Option lead_restart_dir 2
      !% <tt>lead_restart_dir</tt> gives the name of restart directory of the periodic
      !% ground-state calculation for the lead unit cell. Though
      !% one may give different datasets for the left and right lead, they have to
      !% be identical due to the algorithm used to obtain extended eigenstates.
      !% The default is <tt>&lt;lead_dataset&gt;restart</tt>.
      !%Option lead_static_dir 3
      !% The same as <tt>lead_restart_dir</tt> for the <tt>static</tt> directory.
      !% <tt>Octopus</tt> needs the Kohn-Sham potential of the leads. Therefore, the periodic
      !% run must include <tt>Output = potential</tt> in the input file. The default
      !% of this entry is <tt>&lt;lead_dataset&gt;static</tt>.
      !%Option add_unit_cells 4
      !% <tt>add_unit_cells</tt> specifies how many lead unit cells should
      !% be included in the computational domain. Suitable values highly depend
      !% on the system under study but the numbers <tt>nl</tt> and <tt>nr</tt> should
      !% be taken large enough for the potential to equilibrate because we assume
      !% instaneous metallic screening in the leads. Furthermore, note that in
      !% a ground-state calculation, one additional unit cell is added automatically
      !% for each lead to the computational domain because the propagation
      !% algorithm needs the knowledge of the initial state for the first unit cell
      !% outside the simulation box. If omitted, no unit cells are included in the
      !% simulation region (apart from the one which is automatically added in
      !% ground-state calculations).
      !%Option td_pot_formula 5
      !% Defines a spatially local time-dependent potential in the leads as an
      !% analytic expression.
      !%End
      if(parse_block(datasets_check('OpenBoundaries'), blk).eq.0) then
        
        call messages_devel_version("Open boundaries")

        ! Open boundaries are only possible for rectangular simulation boxes.
        if(sb%box_shape.ne.PARALLELEPIPED) then
          message(1) = 'Open boundaries are only possible with a parallelepiped'
          message(2) = 'simulation box.'
          call write_fatal(2)
        end if
        ! Simulation box must not be periodic in transport direction.
        if(sb%periodic_dim.eq.1) then
          message(1) = 'When using open boundaries you cannot use periodic boundary'
          message(2) = 'conditions in x-direction.'
          call write_fatal(2)
        end if

        sb%lead_dataset     = ''
        sb%lead_restart_dir = ''
        sb%lead_static_dir  = ''
        sb%lead_td_pot_formula = '0'
        sb%add_unit_cells   = 0
        nrows = parse_block_n(blk)
        do nr = 0, nrows-1
          call parse_block_integer(blk, nr, 0, tag)
          ncols = parse_block_cols(blk, nr)
          if(ncols.gt.3.or.ncols.lt.2) then
            call input_error('OpenBoundaries')
          end if

          select case(tag)
          case(LEAD_DATASET)
            call parse_block_string(blk, nr, 1, sb%lead_dataset(LEFT))
            if(ncols.eq.3) then
              call parse_block_string(blk, nr, 2, sb%lead_dataset(RIGHT))
              if(trim(sb%lead_dataset(LEFT)).ne.trim(sb%lead_dataset(RIGHT))) then
                message(1) = 'Datasets for left and right lead unit cells must'
                message(2) = 'be equal, i.e. only symmetric leads are possible.'
                call write_fatal(2)
              end if
            else
              sb%lead_dataset(RIGHT) = sb%lead_dataset(LEFT)
            end if
          case(LEAD_RESTART_DIR)
            call parse_block_string(blk, nr, 1, sb%lead_restart_dir(LEFT))
            if(ncols.eq.3) then
              call parse_block_string(blk, nr, 2, sb%lead_restart_dir(RIGHT))
              if(trim(sb%lead_restart_dir(LEFT)).ne.trim(sb%lead_restart_dir(RIGHT))) then
                message(1) = 'Restart directories for left and right lead'
                message(2) = 'unit cells must be equal, i.e. only symmetric'
                message(3) = 'leads are possible.'
                call write_fatal(3)
              end if
            else
              sb%lead_restart_dir(RIGHT) = sb%lead_restart_dir(LEFT)
            end if
          case(LEAD_STATIC_DIR)
            call parse_block_string(blk, nr, 1, sb%lead_static_dir(LEFT))
            if(ncols.eq.3) then
              call parse_block_string(blk, nr, 2, sb%lead_static_dir(RIGHT))
              if(trim(sb%lead_static_dir(LEFT)).ne.trim(sb%lead_static_dir(RIGHT))) then
                message(1) = 'Static directories for left and right lead'
                message(2) = 'unit cells must be equal, i.e. only symmetric'
                message(3) = 'leads are possible.'
                call write_fatal(3)
              end if
            else
              sb%lead_static_dir(RIGHT) = sb%lead_static_dir(LEFT)
            end if
          case(ADD_UNIT_CELLS)
            call parse_block_integer(blk, nr, 1, sb%add_unit_cells(LEFT))
            if(ncols.eq.3) then
              call parse_block_integer(blk, nr, 2, sb%add_unit_cells(RIGHT))
            else
              sb%add_unit_cells(RIGHT) = sb%add_unit_cells(LEFT)
            end if
            if(any(sb%add_unit_cells(1:NLEADS).lt.0)) then
              message(1) = 'add_unit_cells in the OpenBoundaries block must not be negative.'
              call write_fatal(1)
            end if
            ! If we are doing a ground-state calculation add one more unit
            ! cell at both ends to calculate the extended eigenstates.
            if(calc_mode_is(CM_GS)) then
              sb%add_unit_cells(1:NLEADS) = sb%add_unit_cells(1:NLEADS) + 1
            end if
          case(TD_POT_FORMULA)
            call parse_block_string(blk, nr, 1, sb%lead_td_pot_formula(LEFT))
            if(ncols.eq.3) then
              call parse_block_string(blk, nr, 2, sb%lead_td_pot_formula(RIGHT))
            else
              sb%lead_td_pot_formula(RIGHT) = sb%lead_td_pot_formula(LEFT)
            end if
          case default
          end select
        end do
        ! Check if necessary lead_dataset line has been provided.
        if(all(sb%lead_dataset(1:NLEADS).eq.'')) then
          call input_error('OpenBoundaries')
        end if
        ! Set default restart directory.
        if(all(sb%lead_restart_dir(1:NLEADS).eq.'')) then
          do il = 1, NLEADS
            sb%lead_restart_dir(il) = trim(sb%lead_dataset(il))//'restart'
          end do
        end if
        ! Set default static directory.
        if(all(sb%lead_static_dir(1:NLEADS).eq.'')) then
          do il = 1, NLEADS
            sb%lead_static_dir(il) = trim(sb%lead_dataset(il))//'static'
          end do
        end if
        
        sb%open_boundaries = .true.
        ! Read and check the simulation boxes of the lead unit cells.
        ! it has to be allocated directly to avoid problems with xlf 
        allocate(sb%lead_unit_cell(1:NLEADS))
        do il = 1, NLEADS
          call read_lead_unit_cell(sb, il)
        end do
        ! Adjust the size of the simulation box by adding the proper number
        ! of unit cells to the simulation region.
        do il = LEFT, RIGHT
          sb%lsize(TRANS_DIR) = sb%lsize(TRANS_DIR) + sb%add_unit_cells(il)*sb%lead_unit_cell(il)%lsize(TRANS_DIR)
        end do
        sb%n_ucells = nint(sb%lsize(TRANS_DIR)/sb%lead_unit_cell(LEFT)%lsize(TRANS_DIR))
      else
        sb%open_boundaries  = .false.
        sb%add_unit_cells   = 0
        sb%lead_restart_dir = ''
        nullify(sb%lead_unit_cell)
      end if
      call pop_sub()
    end subroutine read_open_boundaries


    !--------------------------------------------------------------
    subroutine read_misc()

      integer              :: i
      type(block_t)        :: blk
      FLOAT,   allocatable :: pos(:)

      call push_sub('simul_box.simul_box_init.read_misc')

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
        call write_fatal(2)
      end if

      sb%dim = calc_dim

      !%Variable PeriodicDimensions
      !%Type integer
      !%Default 0
      !%Section System
      !%Description
      !% Define how many directions are to be considered periodic. Of course, it has to be a number
      !% from zero to three, and it cannot be larger than <tt>Dimensions</tt>.
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
      if ((sb%periodic_dim < 0) .or. (sb%periodic_dim > 3) .or. (sb%periodic_dim > sb%dim)) &
        call input_error('PeriodicDimensions')

      !%Variable MultiResolutionArea
      !%Type block
      !%Section Mesh
      !%Description
      !% (Experimental) Multiresolution regions are set with this
      !% parameter. The first three numbers define the central
      !% point of the region, and the following ones set
      !% the radii where resolution changes (measured from the
      !% central point).
      !% NOTE: up to now, only one area can be set up
      !%End

      if(parse_block(datasets_check('MultiResolutionArea'), blk) == 0) then
        
        call messages_devel_version('Multi-resolution')

        ! number of areas
        sb%hr_area%num_areas = parse_block_n(blk)

        ! number of radii
        sb%hr_area%num_radii = parse_block_cols(blk, 0) - sb%dim

        sb%hr_area%center = M_ZERO

        ! the central point
        do i = 1, sb%dim
           call parse_block_float(blk, 0, i - 1, sb%hr_area%center(i))
        end do

        if (sb%hr_area%num_areas /= 1) call input_error('MultiResolutionArea')

        ! the radii
        SAFE_ALLOCATE(sb%hr_area%radius(1:sb%hr_area%num_radii))
        do i = 1, sb%hr_area%num_radii
          call parse_block_float(blk, 0, sb%dim + i - 1, sb%hr_area%radius(i))
          sb%hr_area%radius(i) = units_to_atomic(units_inp%length, sb%hr_area%radius(i))
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
          call write_fatal(1)
        end if

        sb%hr_area%interp%nn = 2*sb%hr_area%interp%order
        SAFE_ALLOCATE(pos(1:sb%hr_area%interp%nn))
        SAFE_ALLOCATE(sb%hr_area%interp%ww(1:sb%hr_area%interp%nn))
        SAFE_ALLOCATE(sb%hr_area%interp%posi(1:sb%hr_area%interp%nn))
        do i = 1, sb%hr_area%interp%order
          sb%hr_area%interp%posi(i) = 1 + 2*(i - 1)
          sb%hr_area%interp%posi(sb%hr_area%interp%order + i) = -sb%hr_area%interp%posi(i)
          pos(i) =  sb%hr_area%interp%posi(i)
          pos(sb%hr_area%interp%order + i) = -pos(i)
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

      !%Variable OpenBoundariesNLeads
      !%Type integer
      !%Default 2
      !%Section Open Boundaries
      !%Description
      !% The number of leads connected to the central region. Defines the number
      !% of open boundaries for a parallelepiped simulation box shape.
      !%End
      call parse_integer(datasets_check('OpenBoundariesNLeads'), 2, NLEADS)
      if ((NLEADS < 0) .or. (NLEADS > 2*MAX_DIM)) &
        call input_error('OpenBoundariesNLeads')


        call pop_sub()
    end subroutine read_misc


    !--------------------------------------------------------------
    subroutine read_box()
      type(block_t) :: blk
      FLOAT :: default
#if defined(HAVE_GDLIB)        
      character(len=200) :: filename
#endif      

      call push_sub('simul_box.simul_box_init.read_box')
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
      case(SPHERE,MINIMUM,BOX_IMAGE,BOX_USDEF)
        if(sb%dim>1 .and. simul_box_is_periodic(sb)) call input_error('BoxShape')
      case(CYLINDER)
        if (sb%dim>2 .and. &
          ((sb%dim - sb%periodic_dim == 0) .or. (sb%dim - sb%periodic_dim == 1))) call input_error('BoxShape')
      end select

      ! ignore box_shape in 1D
      if(sb%dim==1.and.sb%box_shape /= PARALLELEPIPED.and.sb%box_shape /= HYPERCUBE) sb%box_shape=SPHERE

      ! Cannot use images in 1D or 3D
      if(sb%dim/=2.and.sb%box_shape == BOX_IMAGE) call input_error('BoxShape')

      if(sb%dim > 3 .and. sb%box_shape /= HYPERCUBE) then
        message(1) = "For more than 3 dimensions you can only use the hypercubic box."
        call write_fatal(1)
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
        call parse_float(datasets_check('radius'), units_from_atomic(units_inp%length, def_rsize), sb%rsize)
        if(sb%rsize < M_ZERO) call input_error('radius')
        sb%rsize = units_to_atomic(units_inp%length, sb%rsize)
        if(def_rsize>M_ZERO) call check_def(def_rsize, sb%rsize, 'radius')
      case(MINIMUM)
        default=sb%rsize
        call parse_float(datasets_check('radius'), default, sb%rsize)
        sb%rsize = units_to_atomic(units_inp%length, sb%rsize)
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

        call parse_float(datasets_check('xlength'), units_from_atomic(units_inp%length, default), sb%xsize)
        sb%xsize = units_to_atomic(units_inp%length, sb%xsize)
        sb%lsize(1) = sb%xsize
        if(def_rsize>M_ZERO.and.sb%periodic_dim==0) call check_def(def_rsize, sb%xsize, 'xlength')
      end if

      sb%lsize = M_ZERO
      if(sb%box_shape == PARALLELEPIPED .or.sb%box_shape == HYPERCUBE .or. &
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
          do i = 1, sb%dim
            call parse_block_float(blk, 0, i-1, sb%lsize(i))
          end do
          call parse_block_end(blk)
        else
          call parse_float(datasets_check('Lsize'), -M_ONE, sb%lsize(1))
          if(sb%lsize(1).eq.-M_ONE) then
            call input_error('Lsize')
          end if
          sb%lsize(1:sb%dim) = sb%lsize(1)
        end if
        sb%lsize = units_to_atomic(units_inp%length, sb%lsize)

        do i = 1, sb%dim
          if(def_rsize>M_ZERO.and.sb%periodic_dim<i) call check_def(def_rsize, sb%lsize(i), 'Lsize')
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
          call write_fatal(1)
        end if
#else
        message(1) = "To use 'BoxShape = box_image' you have to compile octopus"
        message(2) = "with GD library support."
        call write_fatal(2)
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
      type(block_t) :: blk
      integer :: i
#if defined(HAVE_GDLIB)
      integer :: sx, sy
#endif

      call push_sub('simul_box.simul_box_init.read_spacing')

      ! initialize to -1
      sb%h = -M_ONE

#if defined(HAVE_GDLIB)
      if(sb%box_shape == BOX_IMAGE) then 
        ! spacing is determined from lsize and the size of the image
        sx = loct_gdImage_SX(sb%image)
        sy = loct_gdImage_SY(sb%image)

        sb%h(1) = M_TWO*sb%lsize(1)/real(sx, REAL_PRECISION)
        sb%h(2) = M_TWO*sb%lsize(2)/real(sy, REAL_PRECISION)
        call pop_sub(); return
      end if
#endif

      !%Variable Spacing
      !%Type float
      !%Section Mesh::Simulation Box
      !%Description
      !% The spacing between the points in the mesh. If using curvilinear
      !% coordinates, this is a canonical spacing that will be changed locally by the
      !% transformation.
      !%
      !% It is possible to have a different spacing in each one of the Cartesian directions
      !% if we define <tt>Spacing</tt> as block of the form
      !%
      !% <tt>%Spacing
      !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
      !% <br>%</tt>
      !%End

      if(parse_block(datasets_check('Spacing'), blk) == 0) then
        if(parse_block_cols(blk,0) < sb%dim) call input_error('Spacing')
        do i = 1, sb%dim
          call parse_block_float(blk, 0, i - 1, sb%h(i), units_inp%length)
        end do
        call parse_block_end(blk)
      else
        call parse_float(datasets_check('Spacing'), sb%h(1), sb%h(1), units_inp%length)
        sb%h(1:sb%dim) = sb%h(1)
      end if

      do i = 1, sb%dim
        if(sb%h(i) < M_ZERO) then
          if(def_h > M_ZERO.and.def_h < huge(def_h)) then
            sb%h(i) = def_h
            write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", i, &
              ") [", trim(units_abbrev(units_out%length)), "] = ",                        &
              units_from_atomic(units_out%length, sb%h(i))
            call write_info(1)
          else
            message(1) = 'Either:'
            message(2) = "   *) variable 'Spacing' is not defined and"
            message(4) = "      I can't find a suitable default"
            message(3) = "   *) your input for 'Spacing' is negative"
            call write_fatal(4)
          end if
        end if
        if(def_rsize>M_ZERO) call check_def(sb%h(i), def_rsize, 'Spacing')
      end do

      call pop_sub()
    end subroutine read_spacing


    !--------------------------------------------------------------
    subroutine read_box_offset()
      integer :: idir
      type(block_t) :: blk

      call push_sub('simul_box.simul_box_init.read_box_offset')
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

      call pop_sub()
    end subroutine read_box_offset


    !--------------------------------------------------------------
    subroutine check_def(var, def, text)
      FLOAT, intent(in) :: var, def
      character(len=*), intent(in) :: text

      call push_sub('simul_box.simul_box_init.check_def')

      if(var > def) then
        write(message(1), '(3a)') "The value for '", text, "' does not match the recommended value"
        write(message(2), '(f8.3,a,f8.3)') var, ' > ', def
        call write_warning(2)
      end if

      call pop_sub()
    end subroutine check_def

    ! ------------------------------------------------------------
    subroutine sb_lookup_init()
      FLOAT, allocatable :: pos(:, :)
      integer :: iatom

      SAFE_ALLOCATE(pos(1:sb%dim, 1:geo%natoms))
     
      do iatom = 1, geo%natoms
        pos(1:sb%dim, iatom) = geo%atom(iatom)%x(1:sb%dim)
      end do
      
      call lookup_init(sb%atom_lookup, sb%dim, geo%natoms, pos)
      
      SAFE_DEALLOCATE_A(pos)

    end subroutine sb_lookup_init

  end subroutine simul_box_init

  !--------------------------------------------------------------

  subroutine simul_box_build_lattice(sb, rlattice_primitive)
    type(simul_box_t), intent(inout) :: sb
    FLOAT,   optional, intent(in)    :: rlattice_primitive(:,:)

    type(block_t) :: blk
    FLOAT :: norm, cross(1:3)
    integer :: idim, jdim

    call push_sub('simul_box.simul_box_build_lattice')
    
    if(present(rlattice_primitive)) then
      sb%rlattice_primitive(1:MAX_DIM,1:MAX_DIM) = rlattice_primitive(1:MAX_DIM,1:MAX_DIM)
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
      forall(idim = 1:MAX_DIM) sb%rlattice_primitive(idim, idim) = M_ONE
      
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
    end select

    call reciprocal_lattice(sb%rlattice_primitive, sb%klattice_primitive, sb%volume_element)
    
    sb%klattice = M_ZERO
    forall(idim=1:sb%periodic_dim, jdim=1:MAX_DIM)
      sb%klattice(jdim, idim) = sb%klattice_primitive(jdim, idim)*M_TWO*M_PI/(M_TWO*sb%lsize(idim))
    end forall
    
    call pop_sub()
  end subroutine simul_box_build_lattice

  
  !--------------------------------------------------------------
  ! Read the coordinates of the leads atoms and add them to the
  ! simulation box (for open boundaries only, of course).
  subroutine simul_box_add_lead_atoms(sb, geo)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(inout) :: geo

    type(geometry_t)  :: central_geo
    type(geometry_t), allocatable  :: lead_geo(:)
    character(len=32) :: label_bak
    integer           :: il, j, n, iatom, icatom, dir

    call push_sub('simul_box.simul_box_add_lead_atoms')

    if(sb%open_boundaries) then
      SAFE_ALLOCATE(lead_geo(1:NLEADS))
      do il = 1, NLEADS
        ! We temporarily change the current label to read the
        ! coordinates of another dataset, namely the lead dataset.
        label_bak     = current_label
        current_label = sb%lead_dataset(il)
        call geometry_init(lead_geo(il), print_info=.false.)
        current_label = label_bak
        call simul_box_atoms_in_box(sb%lead_unit_cell(il), lead_geo(il))
      end do

      ! Merge the geometries of the lead and of the central region.
      call geometry_copy(central_geo, geo)

      ! Set the number of atoms and classical atoms to the number
      ! of atoms coming from left and right lead and central part.
      if(geo%natoms.gt.0) then
        SAFE_DEALLOCATE_P(geo%atom)
      end if
      geo%natoms = central_geo%natoms +                 &
        sb%add_unit_cells(LEFT)*lead_geo(LEFT)%natoms + &
        sb%add_unit_cells(RIGHT)*lead_geo(RIGHT)%natoms
      SAFE_ALLOCATE(geo%atom(1:geo%natoms))
      if(geo%ncatoms.gt.0) then
        SAFE_DEALLOCATE_P(geo%catom)
      end if
      geo%ncatoms = central_geo%ncatoms +                &
        sb%add_unit_cells(LEFT)*lead_geo(LEFT)%ncatoms + &
        sb%add_unit_cells(RIGHT)*lead_geo(RIGHT)%ncatoms
      SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))

      geo%only_user_def = central_geo%only_user_def.and.all(lead_geo(:)%only_user_def)
      geo%nlpp          = central_geo%nlpp.or.any(lead_geo(:)%nlpp)
      geo%nlcc          = central_geo%nlcc.or.any(lead_geo(:)%nlcc)
      geo%atoms%start   = 1
      geo%atoms%end     = geo%natoms
      geo%atoms%nlocal  = geo%natoms
      
      ! 1. Put the atoms of the central region into geo.
      geo%atom(1:central_geo%natoms)   = central_geo%atom
      geo%catom(1:central_geo%ncatoms) = central_geo%catom

      ! 2. Put the atoms of the leads into geo and adjust their x coordinates.
      iatom  = central_geo%natoms+1
      icatom = central_geo%ncatoms+1
      do il = 1, NLEADS
        if(il.eq.LEFT) then
          dir = -1
        else
          dir = 1
        end if
        ! We start from the "outer" unit cells of the lead.
        do j = 1, sb%add_unit_cells(il)
          do n = 1, lead_geo(il)%natoms
            geo%atom(iatom) = lead_geo(il)%atom(n)
            geo%atom(iatom)%x(TRANS_DIR) = geo%atom(iatom)%x(TRANS_DIR) + &
              dir*(sb%lsize(TRANS_DIR) - (2*(j-1)+1)*sb%lead_unit_cell(il)%lsize(TRANS_DIR))
            iatom = iatom + 1
          end do
          do n = 1, lead_geo(il)%ncatoms
            geo%catom(icatom) = lead_geo(il)%catom(n)
            geo%catom(icatom)%x(TRANS_DIR) = geo%catom(icatom)%x(TRANS_DIR) + &
              dir*(sb%lsize(TRANS_DIR) - (2*(j-1)+1)*sb%lead_unit_cell(il)%lsize(TRANS_DIR))
          end do
        end do
      end do
      ! Initialize the species of the "extended" central system.
      if(geo%nspecies.gt.0) then
        SAFE_DEALLOCATE_P(geo%species)
      end if
      call geometry_init_species(geo, print_info=.false.)

      do il = 1, NLEADS
        call geometry_end(lead_geo(il))
      end do
      call geometry_end(central_geo)
      SAFE_DEALLOCATE_A(lead_geo)

    end if
    call pop_sub()
  end subroutine simul_box_add_lead_atoms


  !--------------------------------------------------------------
  subroutine simul_box_atoms_in_box(sb, geo)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(inout) :: geo

    ! This function checks that the atoms are inside the box. If not:
    ! if the system is periodic, the atoms are moved inside the box;
    ! if the system is finite, a warning is written.

    integer :: iatom, pd, idir
    FLOAT :: xx(1:MAX_DIM)

    call push_sub('simul_box.simul_box_atoms_in_box')

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

      if ( .not. simul_box_in_box(sb, geo, geo%atom(iatom)%x) ) then
        write(message(1), '(a,i5,a)') "Atom ", iatom, " is outside the box."
        if (sb%periodic_dim == sb%dim) then
          message(2) = "This is a bug."
          call write_fatal(2)
        else
          call write_warning(1)
        end if
      end if

    end do

    call pop_sub()
  end subroutine simul_box_atoms_in_box


  !--------------------------------------------------------------
  subroutine reciprocal_lattice(rv, kv, volume)
    FLOAT, intent(in)    :: rv(MAX_DIM, MAX_DIM)
    FLOAT, intent(out)   :: kv(MAX_DIM, MAX_DIM)
    FLOAT, intent(out)   :: volume

    FLOAT :: tmp(1:MAX_DIM)

    call push_sub('simul_box.reciprocal_lattice')
    
    kv = M_ZERO

    tmp(1:3) = dcross_product(rv(1:3, 2), rv(1:3, 3)) 
    volume = dot_product(rv(1:3, 1), tmp(1:3))

    if ( volume < M_ZERO ) then 
      message(1) = "Your lattice vectors form a left-handed system."
      call write_fatal(1)
    end if

    kv(1:3, 1) = dcross_product(rv(:, 2), rv(:, 3))/volume
    kv(1:3, 2) = dcross_product(rv(:, 3), rv(:, 1))/volume
    kv(1:3, 3) = dcross_product(rv(:, 1), rv(:, 2))/volume    

    call pop_sub()
  end subroutine reciprocal_lattice


  !--------------------------------------------------------------
  subroutine simul_box_end(sb)
    type(simul_box_t), intent(inout) :: sb    

    call push_sub('simul_box.simul_box_end')

    if(simul_box_is_periodic(sb)) call symmetries_end(sb%symm)

    call lookup_end(sb%atom_lookup)
    call kpoints_end(sb%kpoints)

    if(sb%open_boundaries) then
      ! deallocated directly to avoid problems with xlf
      if(associated(sb%lead_unit_cell)) deallocate(sb%lead_unit_cell)
      nullify(sb%lead_unit_cell)
    end if

    SAFE_DEALLOCATE_P(sb%hr_area%radius)
    SAFE_DEALLOCATE_P(sb%hr_area%interp%ww)
    SAFE_DEALLOCATE_P(sb%hr_area%interp%posi)

    call pop_sub()
  end subroutine simul_box_end


  !--------------------------------------------------------------
  recursive subroutine simul_box_write_info(sb, geo, iunit)
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

    call push_sub('simul_box.simul_box_write_info')

    write(message(1),'(a)') 'Simulation Box:'
    if(sb%box_shape.eq.BOX_USDEF) then
      write(message(2), '(a)') '  Type = user-defined'
    else
      write(message(2), '(a,a,1x)') '  Type = ', bs(sb%box_shape)
    end if
    call write_info(2, iunit)

    if(sb%box_shape == SPHERE .or. sb%box_shape == CYLINDER &
       .or. (sb%box_shape == MINIMUM .and. sb%rsize > M_ZERO)) then
      write(message(1), '(3a,f7.3)') '  Radius  [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%rsize)
      call write_info(1, iunit)
    endif

    if (sb%box_shape == MINIMUM .and. sb%rsize <= M_ZERO) then
      do ispec = 1, geo%nspecies     
        write(message(1), '(a,a5,5x,a,f7.3,2a)') '  Species = ', trim(species_label(geo%species(ispec))), 'Radius = ', &
          units_from_atomic(units_out%length, species_def_rsize(geo%species(ispec))), ' ', trim(units_abbrev(units_out%length))
        call write_info(1, iunit)
      enddo
    end if

    if(sb%box_shape == CYLINDER) then
      write(message(1), '(3a,f7.3)') '  Xlength [', trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, sb%xsize)
      call write_info(1, iunit)
    end if

    if(sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(3a, a, f8.3, a, f8.3, a, f8.3, a)')     &
        '  Lengths [', trim(units_abbrev(units_out%length)), '] = ',    &
        '(', units_from_atomic(units_out%length, sb%lsize(1)), ',',           &
        units_from_atomic(units_out%length, sb%lsize(2)), ',',                &
        units_from_atomic(units_out%length, sb%lsize(3)), ')'
      call write_info(1, iunit)
    end if

    write(message(1), '(a,i1,a)') '  Octopus will run in ', sb%dim, ' dimension(s).'
    write(message(2), '(a,i1,a)') '  Octopus will treat the system as periodic in ', &
      sb%periodic_dim, ' dimension(s).'
    call write_info(2, iunit)

    if(sb%periodic_dim > 0 .or. sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(1x)')
      write(message(2),'(a,3a,a)') '  Lattice Vectors [', trim(units_abbrev(units_out%length)), ']'
      do idir = 1, sb%dim
        write(message(2+idir),'(9f12.6)') (units_from_atomic(units_out%length, sb%rlattice(idir2, idir)), &
                                         idir2 = 1, sb%dim) 
      end do
      call write_info(2+sb%dim, iunit)

      write(message(1),'(a,f18.4,3a,i1.1,a)') &
        '  Cell volume = ', units_from_atomic(units_out%length**sb%dim, sb%rcell_volume), &
        ' [', trim(units_abbrev(units_out%length**sb%dim)), ']'
      call write_info(1, iunit)

      write(message(1),'(a,3a,a)') '  Reciprocal-Lattice Vectors [', trim(units_abbrev(units_out%length**(-1))), ']'
      do idir = 1, sb%dim
        write(message(1+idir),'(3f12.6)') (units_from_atomic(unit_one / units_out%length, sb%klattice(idir2, idir)), &
                                           idir2 = 1, sb%dim)
      end do
      call write_info(1+sb%dim, iunit)
    end if

    if(sb%open_boundaries) then
      write(message(1), '(a)')       'Open boundaries in x-direction:'
      write(message(2), '(a,2i4)') '  Additional unit cells left:    ', &
        sb%add_unit_cells(LEFT)
      write(message(3), '(a,2i4)') '  Additional unit cells right:   ', &
        sb%add_unit_cells(RIGHT)
      write(message(4), '(a)')     '  Left lead read from directory:  '// &
        trim(sb%lead_restart_dir(LEFT))
      write(message(5), '(a)')     '  Right lead read from directory: '// &
        trim(sb%lead_restart_dir(RIGHT))
      call write_info(5, iunit)
    end if

    call pop_sub()
  end subroutine simul_box_write_info


  !--------------------------------------------------------------
  logical function simul_box_in_box(sb, geo, x) result(in_box)
    type(simul_box_t),  intent(in) :: sb
    type(geometry_t),   intent(in) :: geo
    FLOAT,              intent(in) :: x(:)

    real(8), parameter :: DELTA = CNST(1e-12)
    FLOAT :: xx(1:MAX_DIM, 1)
    logical :: in_box2(1)

    ! no push_sub because this function is called very frequently

    xx = M_ZERO
    xx(1:sb%dim, 1) = x(1:sb%dim)

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
    FLOAT :: r, re, im, dist2, radius
    real(8) :: llimit(MAX_DIM), ulimit(MAX_DIM)
    FLOAT, allocatable :: xx(:, :)
    integer :: ip, idir, iatom, ilist
    integer, allocatable :: nlist(:)
    integer, pointer :: list(:, :)

#if defined(HAVE_GDLIB)
    integer :: red, green, blue, ix, iy
#endif

    ! no push_sub because this function is called very frequently
    SAFE_ALLOCATE(xx(1:MAX_DIM, 1:npoints))

    if(MAX_DIM /= sb%dim) xx = M_ZERO

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
          r = sqrt(xx(2, ip)**2 + xx(3, ip)**2)
          in_box(ip) = (r <= sb%rsize + DELTA .and. abs(xx(1, ip)) <= sb%xsize+DELTA)
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

              dist2 = sum((xx(1:MAX_DIM, ip) - geo%atom(iatom)%x(1:MAX_DIM))**2)

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
        ulimit(1:sb%periodic_dim) = sb%lsize(1:sb%periodic_dim) - DELTA

        if(sb%open_boundaries) then
          ulimit(TRANS_DIR) = sb%lsize(TRANS_DIR) - DELTA
        end if

        forall(ip = 1:npoints)
          in_box(ip) = all(xx(1:sb%dim, ip) >= llimit(1:sb%dim) .and. xx(1:sb%dim, ip) <= ulimit(1:sb%dim))
        end forall

#if defined(HAVE_GDLIB)
      case(BOX_IMAGE)
        do ip = 1, npoints
          ix = int((xx(1, ip) + sb%lsize(1))/sb%h(1))
          iy = int((xx(2, ip) + sb%lsize(2))/sb%h(2))
          call loct_gdimage_get_pixel_rgb(sb%image, ix, iy, red, green, blue)
          in_box(ip) = (red == 255).and.(green == 255).and.(blue == 255)
        end do
#endif

      case(BOX_USDEF)
        ! is it inside the user-given boundaries
        do ip = 1, npoints
          in_box(ip) =  &
            (xx(1, ip) >= -sb%lsize(1)-DELTA.and.xx(1, ip) <= sb%lsize(1)+DELTA).and. &
            (xx(2, ip) >= -sb%lsize(2)-DELTA.and.xx(2, ip) <= sb%lsize(2)+DELTA).and. &
            (xx(3, ip) >= -sb%lsize(3)-DELTA.and.xx(3, ip) <= sb%lsize(3)+DELTA)

          ! and inside the simulation box
          do idir = 1, sb%dim
            xx(idir, ip) = units_from_atomic(units_inp%length, xx(idir, ip))
          enddo
          r = sqrt(sum(xx(:, ip)**2))
          call parse_expression(re, im, sb%dim, xx(:, ip), r, M_ZERO, sb%user_def)
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

    integer :: i

    call push_sub('simul_box.simul_box_dump')

    write(iunit, '(a)')             dump_tag
    write(iunit, '(a20,i4)')        'box_shape=          ', sb%box_shape
    write(iunit, '(a20,i4)')        'dim=                ', sb%dim
    write(iunit, '(a20,i4)')        'periodic_dim=       ', sb%periodic_dim
    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,9e22.14)') 'lsize=              ', sb%lsize(1:MAX_DIM)
    case(CYLINDER)
      write(iunit, '(a20,e22.14)')  'rsize=              ', sb%rsize
      write(iunit, '(a20,e22.14)')  'xlength=            ', sb%xsize
      write(iunit, '(a20,9e22.14)') 'lsize=              ', sb%lsize(1:MAX_DIM)
    case(PARALLELEPIPED)
      write(iunit, '(a20,9e22.14)') 'lsize=              ', sb%lsize(1:MAX_DIM)
    case(BOX_USDEF)
      write(iunit, '(a20,9e22.14)') 'lsize=              ', sb%lsize(1:MAX_DIM)
      write(iunit, '(a20,a1024)')   'user_def=           ', sb%user_def
    end select
    write(iunit, '(a20,e22.14)')    'fft_alpha=          ', sb%fft_alpha
    write(iunit, '(a20,9e22.14)')   'h=                  ', sb%h(1:MAX_DIM)
    write(iunit, '(a20,9e22.14)')   'box_offset=         ', sb%box_offset(1:MAX_DIM)
    write(iunit, '(a20,l7)')        'mr_flag=            ', sb%mr_flag
    if(sb%mr_flag) then
      write(iunit, '(a20,i4)')        'num_areas=         ',sb%hr_area%num_areas
      write(iunit, '(a20,i4)')        'num_radii=         ',sb%hr_area%num_radii
      do i = 1, sb%hr_area%num_radii
        write(iunit, '(a10,i2.2,a9,e22.14)') 'mr_radius_', i, '=        ',sb%hr_area%radius(i)
      end do
      do i = 1, MAX_DIM
        write(iunit, '(a7,i1,a13,e22.14)')   'center(', i, ')=           ',sb%hr_area%center(i)
      end do
    end if
    do i = 1, MAX_DIM
      write(iunit, '(a9,i1,a11,9e22.14)')    'rlattice(', i, ')=         ', &
        sb%rlattice_primitive(1:MAX_DIM, i)
    end do
    write(iunit, '(a20,l7)')        'open_boundaries=    ', sb%open_boundaries
    if(sb%open_boundaries) then
      write(iunit, '(a20,2i4)')     'add_unit_cells=     ', sb%add_unit_cells(1:NLEADS) ! FIXME for NLEADS>2
      write(iunit, '(a20,a32)')     'lead_restart_dir(L)=', sb%lead_restart_dir(LEFT)
      write(iunit, '(a20,a32)')     'lead_restart_dir(R)=', sb%lead_restart_dir(RIGHT)
    end if

    call pop_sub()
  end subroutine simul_box_dump


  !--------------------------------------------------------------
  ! Read the simulation box for lead il from sb%lead_restart_dir(il) into
  ! sb%lead_unit_cell(il).
  subroutine read_lead_unit_cell(sb, il)
    type(simul_box_t), intent(inout) :: sb
    integer,           intent(in)    :: il

    integer :: iunit

    call push_sub('simul_box.read_lead_unit_cell')

    iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/mesh', action='read', is_tmp=.true., grp = mpi_world)
    call simul_box_init_from_file(sb%lead_unit_cell(il), iunit)
    call io_close(iunit)
    ! Check, if
    ! * simulation box is a parallelepiped,
    ! * lead simulation box has the same spacing than central region,
    ! * the extensions in y-, z-directions fit the central box,
    ! * the central simulation box x-length is an integer multiple of
    !   the unit cell x-length,
    ! * periodic in one dimension, and
    ! * of the same dimensionality as the central system.
    if(sb%lead_unit_cell(il)%box_shape.ne.PARALLELEPIPED) then
      message(1) = 'Simulation box of '//LEAD_NAME(il)//' lead is not a parallelepiped.'
      call write_fatal(1)
    end if
    if(any(sb%h.ne.sb%lead_unit_cell(il)%h)) then
      message(1) = 'Simulation box of '//LEAD_NAME(il)//' has a different spacing than'
      message(2) = 'the central system.'
      call write_fatal(2)
    end if
    if(any(sb%lsize(2:3).ne.sb%lead_unit_cell(il)%lsize(2:3))) then
      message(1) = 'The size in y-, z-directions of the '//LEAD_NAME(LEFT)//' lead'
      message(2) = 'does not fit the size of the y-, z-directions of the central system.'
      call write_fatal(2)
    end if
    if(.not.is_integer_multiple(sb%lsize(1), sb%lead_unit_cell(il)%lsize(1))) then
      message(1) = 'The length in x-direction of the central simulation'
      message(2) = 'box is not an integer multiple of the x-length of'
      message(3) = 'the '//trim(LEAD_NAME(il))//' lead.'
      call write_fatal(3)
    end if
    if(sb%lead_unit_cell(il)%periodic_dim.ne.1) then
      message(1) = 'Simulation box of '//LEAD_NAME(il)//' lead is not periodic in x-direction.'
!      call write_fatal(1)
      call write_warning(1)
    end if
    if(sb%lead_unit_cell(il)%dim.ne.calc_dim) then
      message(1) = 'Simulation box of '//LEAD_NAME(il)//' has a different dimension than'
      message(2) = 'the central system.'
      call write_fatal(2)
    end if

    call pop_sub()
  end subroutine read_lead_unit_cell


  ! --------------------------------------------------------------
  ! Returns the extent of the lead unit cell of lead il in transport
  ! direction. Returns -1 if open boundaries are not used.
  integer function lead_unit_cell_extent(sb, il)
    type(simul_box_t), intent(in) :: sb
    integer,           intent(in) :: il

    integer :: tdir ! transport direction
    
    call push_sub('simul_box.lead_unit_cell_extent')
    
    if(sb%open_boundaries) then
      tdir = (il+1)/2
      lead_unit_cell_extent = nint(2*sb%lead_unit_cell(il)%lsize(tdir)/sb%h(tdir))
    else
      lead_unit_cell_extent = -1
    end if

    call pop_sub()
  end function lead_unit_cell_extent


  ! --------------------------------------------------------------
  subroutine simul_box_init_from_file(sb, iunit)
    type(simul_box_t), intent(inout) :: sb
    integer,           intent(in)    :: iunit

    character(len=20)  :: str
    character(len=300) :: line
    integer            :: idim, il, ierr
    FLOAT              :: rlattice_primitive(1:MAX_DIM, 1:MAX_DIM)

    call push_sub('simul_box.simul_box_init_from_file')

    ! Find (and throw away) the dump tag.
    do
      call iopar_read(mpi_world, iunit, line, ierr)
      if(trim(line).eq.dump_tag) exit
    end do

    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%box_shape
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%dim
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%periodic_dim

    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%rsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:MAX_DIM)       
    case(CYLINDER)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%rsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%xsize
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:MAX_DIM)
    case(PARALLELEPIPED)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:MAX_DIM)
    case(BOX_USDEF)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lsize(1:MAX_DIM)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%user_def
    end select
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%fft_alpha
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%h(1:MAX_DIM)
    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%box_offset(1:MAX_DIM)
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
      do idim=1, MAX_DIM
        call iopar_read(mpi_world, iunit, line, ierr)
        read(line, *) str, sb%hr_area%center(idim)
      end do
    end if
    do idim=1, MAX_DIM
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, rlattice_primitive(1:MAX_DIM, idim)
    end do

    call simul_box_build_lattice(sb, rlattice_primitive)

    call iopar_read(mpi_world, iunit, line, ierr)
    read(line, *) str, sb%open_boundaries
    if(sb%open_boundaries) then
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%add_unit_cells(LEFT), sb%add_unit_cells(RIGHT)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lead_restart_dir(LEFT)
      call iopar_read(mpi_world, iunit, line, ierr)
      read(line, *) str, sb%lead_restart_dir(RIGHT)
      SAFE_ALLOCATE(sb%lead_unit_cell(1:NLEADS))
      do il = 1, NLEADS
        call read_lead_unit_cell(sb, il)
      end do
    end if

    call pop_sub()
  end subroutine simul_box_init_from_file

  ! --------------------------------------------------------------
  recursive subroutine simul_box_copy(sbout, sbin)
    type(simul_box_t), intent(out) :: sbout
    type(simul_box_t), intent(in)  :: sbin

    integer :: il

    call push_sub('simul_box.simul_box_copy')

    sbout%box_shape               = sbin%box_shape
    sbout%h                       = sbin%h
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
    sbout%open_boundaries         = sbin%open_boundaries
    sbout%add_unit_cells          = sbin%add_unit_cells
    sbout%lead_restart_dir        = sbin%lead_restart_dir
    sbout%mr_flag                 = sbin%mr_flag
    sbout%hr_area%num_areas       = sbin%hr_area%num_areas
    sbout%hr_area%num_radii       = sbin%hr_area%num_radii
    sbout%hr_area%center(1:MAX_DIM)=sbin%hr_area%center(1:MAX_DIM)
    
    call kpoints_copy(sbout%kpoints, sbin%kpoints)

    if(sbout%mr_flag) then
      SAFE_ALLOCATE(sbout%hr_area%radius(1:sbout%hr_area%num_radii))
      sbout%hr_area%radius(1:sbout%hr_area%num_radii) = sbin%hr_area%radius(1:sbout%hr_area%num_radii)
    end if

    if(associated(sbin%lead_unit_cell)) then
      SAFE_ALLOCATE(sbout%lead_unit_cell(1:NLEADS))
      do il = 1, NLEADS
        call simul_box_copy(sbout%lead_unit_cell(il), sbin%lead_unit_cell(il))
      end do
    end if

    call lookup_copy(sbin%atom_lookup, sbout%atom_lookup)

    if(simul_box_is_periodic(sbin)) call symmetries_copy(sbin%symm, sbout%symm)

    call pop_sub()
  end subroutine simul_box_copy

end module simul_box_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
