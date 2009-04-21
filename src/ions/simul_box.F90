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
  use loct_parser_m
  use math_m
  use messages_m
  use profiling_m
  use string_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    simul_box_t,                &
    simul_box_init,             &
    simul_box_end,              &
    simul_box_write_info,       &
    simul_box_is_periodic,      &
    simul_box_is_eq,            &
    simul_box_in_box,           &
    simul_box_dump,             &
    simul_box_init_from_file,   &
    simul_box_atoms_in_box,     &
    lead_unit_cell_extent,      &
    simul_box_multires,         &
    simul_box_copy

  integer, parameter, public :: &
    SPHERE         = 1,         &
    CYLINDER       = 2,         &
    MINIMUM        = 3,         &
    PARALLELEPIPED = 4,         &
    BOX_IMAGE      = 5,         &
    HYPERCUBE      = 6,         &
    BOX_USDEF      = 123

  integer, parameter, public :: &
    LEFT      = 1,              & ! Lead indices,
    RIGHT     = 2,              & ! L=1, R=2.
    OUTER     = 1,              & ! Block indices of interface wave functions
    INNER     = 2,              & ! for source term.
    NLEADS    = 2,              & ! Number of leads.
    TRANS_DIR = 1                 ! Transport is in x-direction.

  character(len=5), dimension(NLEADS), parameter, public :: LEAD_NAME = &
    (/'left ', 'right'/)

  type simul_box_t
    integer  :: box_shape   ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                            ! 4->parallelepiped (orthonormal, up to now).

    FLOAT :: h(MAX_DIM)     ! the (canonical) spacing between the points
    FLOAT :: box_offset(MAX_DIM)  ! shifts of the origin in the respective direction

    FLOAT :: rsize          ! the radius of the sphere or of the cylinder
    FLOAT :: xsize          ! the length of the cylinder in the x direction
    FLOAT :: lsize(MAX_DIM) ! half of the length of the parallelepiped in each direction.

    type(c_ptr)   :: image    ! for the box defined through an image
    character(len=1024) :: user_def ! for the user-defined box

    FLOAT :: inner_size

    FLOAT :: rlattice(MAX_DIM,MAX_DIM)      ! lattice primitive vectors
    FLOAT :: klattice_unitary(MAX_DIM,MAX_DIM)      ! reciprocal lattice primitive vectors
    FLOAT :: klattice(MAX_DIM,MAX_DIM)      ! reciprocal lattice primitive vectors
    FLOAT :: volume_element     ! the volume element in real space
    FLOAT :: rcell_volume       ! the volume of the cell in real space
    FLOAT :: fft_alpha      ! enlargement factor for double box

    integer :: dim
    integer :: periodic_dim

    ! For open boundaries, we need reference to the lead`s unit cell.
    ! This unit cell is itself a simulation box.
    logical                    :: open_boundaries          ! Use open boundaries?
    integer                    :: n_ucells                 ! Number of unit cells that fit in central region.
    integer                    :: add_unit_cells(NLEADS)   ! Number of additonal unit cells.
    character(len=32)          :: lead_dataset(NLEADS)     ! Dataset name of the periodic lead calculation.
    character(len=32)          :: lead_restart_dir(NLEADS) ! Directory where to find the lead restart files.
    character(len=32)          :: lead_static_dir(NLEADS)  ! Static directory of the lead ground state.
    type(simul_box_t), pointer :: lead_unit_cell(:)        ! Simulation box of the unit cells.
    ! The next one does not really belong here but since the parsing happens in the simul_box_init
    ! it makes things a bit easier.
    character(len=1000)        :: lead_td_pot_formula(NLEADS) ! Td-potential of lead.
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
    call read_spacing ()                   ! Parameters defining the (canonical) spacing.
    call read_box_offset()                 ! Parameters defining the offset of the origin.
    call read_open_boundaries()            ! Parameters defining open boundaries.
    call build_lattice()                   ! Build lattice vectors.
    call simul_box_add_lead_atoms(sb, geo) ! Add the atoms of the lead unit cells that are
                                           ! included in the simulation box to geo.
    call simul_box_atoms_in_box(sb, geo)   ! Put all the atoms inside the box.

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
      !% Enables open boundaries in x-direction and defines the character
      !% of the leads attached left and right ot the finite central system.
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
      !% leads ground state. It is used, e. g., to read in the coordinates of the
      !% atoms of the lead. Both entries for left and right have to be equal.
      !%Option lead_restart_dir 2
      !% <tt>lead_restart_dir</tt> gives the name of restart directory of the periodic
      !% ground-state calculation for the lead unit cell. Though
      !% one may give different datasets for the left and right lead, they have to
      !% be identical due to the algorithm used to obtain extended eigenstates.
      !% The default is <tt>&lt;lead_dataset&gt;restart</tt>.
      !%Option lead_static_dir 3
      !% The same as <tt>lead_restart_dir</tt> for the <tt>static</tt> directory.
      !% Octopus needs the Kohn-Sham potential of the leads. Therefore, the periodic
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
      if(loct_parse_block(datasets_check('OpenBoundaries'), blk).eq.0) then
        
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
        nrows = loct_parse_block_n(blk)
        do nr = 0, nrows-1
          call loct_parse_block_int(blk, nr, 0, tag)
          ncols = loct_parse_block_cols(blk, nr)
          if(ncols.gt.3.or.ncols.lt.2) then
            call input_error('OpenBoundaries')
          end if

          select case(tag)
          case(LEAD_DATASET)
            call loct_parse_block_string(blk, nr, 1, sb%lead_dataset(LEFT))
            if(ncols.eq.3) then
              call loct_parse_block_string(blk, nr, 2, sb%lead_dataset(RIGHT))
              if(trim(sb%lead_dataset(LEFT)).ne.trim(sb%lead_dataset(RIGHT))) then
                message(1) = 'Datasets for left and right lead unit cells must'
                message(2) = 'be equal, i. e. only symmetric leads are possible.'
                call write_fatal(2)
              end if
            else
              sb%lead_dataset(RIGHT) = sb%lead_dataset(LEFT)
            end if
          case(LEAD_RESTART_DIR)
            call loct_parse_block_string(blk, nr, 1, sb%lead_restart_dir(LEFT))
            if(ncols.eq.3) then
              call loct_parse_block_string(blk, nr, 2, sb%lead_restart_dir(RIGHT))
              if(trim(sb%lead_restart_dir(LEFT)).ne.trim(sb%lead_restart_dir(RIGHT))) then
                message(1) = 'Restart directories for left and right lead'
                message(2) = 'unit cells must be equal, i. e. only symmetric'
                message(3) = 'leads are possible.'
                call write_fatal(3)
              end if
            else
              sb%lead_restart_dir(RIGHT) = sb%lead_restart_dir(LEFT)
            end if
          case(LEAD_STATIC_DIR)
            call loct_parse_block_string(blk, nr, 1, sb%lead_static_dir(LEFT))
            if(ncols.eq.3) then
              call loct_parse_block_string(blk, nr, 2, sb%lead_static_dir(RIGHT))
              if(trim(sb%lead_static_dir(LEFT)).ne.trim(sb%lead_static_dir(RIGHT))) then
                message(1) = 'Static directories for left and right lead'
                message(2) = 'unit cells must be equal, i. e. only symmetric'
                message(3) = 'leads are possible.'
                call write_fatal(3)
              end if
            else
              sb%lead_static_dir(RIGHT) = sb%lead_static_dir(LEFT)
            end if
          case(ADD_UNIT_CELLS)
            call loct_parse_block_int(blk, nr, 1, sb%add_unit_cells(LEFT))
            if(ncols.eq.3) then
              call loct_parse_block_int(blk, nr, 2, sb%add_unit_cells(RIGHT))
            else
              sb%add_unit_cells(RIGHT) = sb%add_unit_cells(LEFT)
            end if
            if(any(sb%add_unit_cells.lt.0)) then
              message(1) = 'add_unit_cells in the OpenBoundaries block must not be negative.'
              call write_fatal(1)
            end if
            ! If we are doing a ground-state calculation add one more unit
            ! cell at both ends to calculate the extended eigenstates.
            if(calc_mode_is(CM_GS)) then
              sb%add_unit_cells = sb%add_unit_cells + 1
            end if
          case(TD_POT_FORMULA)
            call loct_parse_block_string(blk, nr, 1, sb%lead_td_pot_formula(LEFT))
            if(ncols.eq.3) then
              call loct_parse_block_string(blk, nr, 2, sb%lead_td_pot_formula(RIGHT))
            else
              sb%lead_td_pot_formula(RIGHT) = sb%lead_td_pot_formula(LEFT)
            end if
          case default
          end select
        end do
        ! Check if necessary lead_dataset line has been provided.
        if(all(sb%lead_dataset.eq.'')) then
          call input_error('OpenBoundaries')
        end if
        ! Set default restart directory.
        if(all(sb%lead_restart_dir.eq.'')) then
          do il = 1, NLEADS
            sb%lead_restart_dir(il) = trim(sb%lead_dataset(il))//'restart'
          end do
        end if
        ! Set default static directory.
        if(all(sb%lead_static_dir.eq.'')) then
          do il = 1, NLEADS
            sb%lead_static_dir(il) = trim(sb%lead_dataset(il))//'static'
          end do
        end if
        
        sb%open_boundaries = .true.
        ! Read and check the simulation boxes of the lead unit cells.
        ALLOCATE(sb%lead_unit_cell(NLEADS), NLEADS)
        do il = 1, NLEADS
          call read_lead_unit_cell(sb, il)
        end do
        ! Adjust the size of the simulation box by adding the proper number
        ! of unit cells to the simulation region.
        do il = 1, NLEADS
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

      call push_sub('simul_box.simul_box_init.read_misc')

      !%Variable DoubleFFTParameter
      !%Type float
      !%Default 2.0
      !%Section Mesh::FFTs
      !%Description
      !% For solving Poisson equation in Fourier space, and for applying the local potential
      !% in Fourier space, an auxiliary cubic mesh is built. This mesh will be larger than
      !% the circumscribed cube to the usual mesh by a factor <tt>DoubleFFTParameter</tt>. See
      !% the section that refers to Poisson equation, and to the local potential for details
      !% [the default value of two is typically good].
      !%End
      call loct_parse_float(datasets_check('DoubleFFTParameter'), M_TWO, sb%fft_alpha)
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
      call loct_parse_int(datasets_check('PeriodicDimensions'), 0, sb%periodic_dim)
      if ((sb%periodic_dim < 0) .or. (sb%periodic_dim > 3) .or. (sb%periodic_dim > sb%dim)) &
        call input_error('PeriodicDimensions')

      !%Variable HighResolutionArea
      !%Type float
      !%Default 1.0
      !%Section Mesh
      !%Description
      !% (Experimental) When this variable is set, octopus will use a
      !% multiresolution grid with a central region with the specified
      !% spacing and an outer region with double the spacing. This
      !% variable controls the size of the inner region as a fraction
      !% (in each direction) of the simulation box. So for instance if
      !% your box is a sphere of radius 10, and you set
      !% HighResolutionArea to 0.5, you will get a sphere of radius 5
      !% with higher resolution. By default this variable is set to 1,
      !% so this feature is disabled.
      !%End
      call loct_parse_float(datasets_check('HighResolutionArea'), M_ONE, sb%inner_size)
      sb%inner_size = min(sb%inner_size, M_ONE)

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
      !% <li>Box_image is only allowed in 2D.</li>
      !% </ul>
      !%Option sphere 1
      !% The simulation box will be a sphere of radius Radius
      !%Option cylinder 2
      !% The simulation box will be a cylinder with radius Radius and height two times
      !% Xlength
      !%Option minimum 3
      !% The simulation box will be constructed by adding spheres created around each
      !% atom (or user-defined potential), of radius Radius.
      !%Option parallelepiped 4
      !% The simulation box will be a parallelepiped whose dimensions are taken from
      !% the variable lsize.
      !%Option box_image 5
      !% The simulation box will be defined through an image. White means that the point
      !% is contained in the simulation box, while any other color means that the point is out.
      !%Option user_defined 123
      !% The shape of the simulation box will be read from the variable <tt>BoxShapeUsDef</tt>
      !%Option hypercube 6
      !% (experimental) The simulation box will an hypercube or
      !% hyperparallelepiped, this is equivalent to the
      !% <tt>parallelepiped</tt> box but it can work with an arbitrary
      !% number of dimensions.
      !%End

      call loct_parse_int(datasets_check('BoxShape'), MINIMUM, sb%box_shape)
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
      !% If BoxShape is not "parallelepiped", defines the radius of the spheres or of the cylinder.
      !% It has to be a positive number. If it is not defined in the input file, then the program
      !% will attempt to find a suitable default. However, this is not always possible, in which case
      !% the code will stop, issuing this error message.
      !% For the minimum option, if the radius is not specified, a radius is chosen separately for each species.
      !% If default pseudopotentials are used, the radii are read from the "rsize" column of share/PP/defaults.
      !% Otherwise the radii are read from the Species block.
      !%End
      select case(sb%box_shape)
      case(SPHERE, CYLINDER)
        call loct_parse_float(datasets_check('radius'), def_rsize/units_inp%length%factor, sb%rsize)
        if(sb%rsize < CNST(0.0)) call input_error('radius')
        sb%rsize = sb%rsize * units_inp%length%factor
        if(def_rsize>M_ZERO) call check_def(def_rsize, sb%rsize, 'radius')
      case(MINIMUM)
        default=sb%rsize
        call loct_parse_float(datasets_check('radius'), default, sb%rsize)
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
        call loct_parse_float(datasets_check('xlength'), M_ONE/units_inp%length%factor, sb%xsize)
        sb%xsize = sb%xsize * units_inp%length%factor
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
        !% In case BoxShape is "parallelepiped", "hypercube",
        !% "box_image", or "user_defined", this is assumed to be a
        !% block of the form:
        !%
        !% <tt>%Lsize
        !% <br>&nbsp;&nbsp;sizex | sizey | sizez | ...
        !% <br>%</tt>
        !%
        !% where the "size*" are half the lengths of the box in each direction.
        !%
        !% The number of columns must match the dimensionality of the
        !% calculation. If you want a cube you can also set Lsize as a
        !% single variable.
        !%End

        if(loct_parse_block(datasets_check('Lsize'), blk) == 0) then
          if(loct_parse_block_cols(blk,0) < sb%dim) call input_error('Lsize')
          do i = 1, sb%dim
            call loct_parse_block_float(blk, 0, i-1, sb%lsize(i))
          end do
          call loct_parse_block_end(blk)
        else
          call loct_parse_float(datasets_check('Lsize'), -M_ONE, sb%lsize(1))
          if(sb%lsize(1).eq.-M_ONE) then
            call input_error('Lsize')
          end if
          sb%lsize(1:sb%dim) = sb%lsize(1)
        end if
        sb%lsize = sb%lsize*units_inp%length%factor

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
        call loct_parse_string(datasets_check("BoxShapeImage"), "", filename)
        sb%image = loct_gdimage_create_from(filename)
        if(.not.c_associated(sb%image)) then
          message(1) = "Could not open file '" // filename // "'"
          call write_fatal(1)
        end if
#else
        message(1) = "To use 'BoxShape = box_image' you have to compile octopus"
        message(2) = "with GD library support"
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
        !% with axis parallel to the z axis
        !%End
        
        call loct_parse_string(datasets_check("BoxShapeUsDef"), "x^2+y^2+z^2 < 4", sb%user_def)
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

      if(loct_parse_block(datasets_check('Spacing'), blk) == 0) then
        if(loct_parse_block_cols(blk,0) < sb%dim) call input_error('Spacing')
        do i = 1, sb%dim
          call loct_parse_block_float(blk, 0, i-1, sb%h(i))
        end do
        call loct_parse_block_end(blk)
      else
        call loct_parse_float(datasets_check('Spacing'), sb%h(1), sb%h(1))
        sb%h(1:sb%dim) = sb%h(1)
      end if

      do i = 1, sb%dim
        sb%h(i) = sb%h(i)*units_inp%length%factor
        if(sb%h(i) < M_ZERO) then
          if(def_h > M_ZERO.and.def_h < huge(REAL_PRECISION)) then
            sb%h(i) = def_h
            write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", i, &
              ") [", trim(units_out%length%abbrev), "] = ",                        &
              sb%h(i)/units_out%length%factor
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
      integer :: i
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
      if(loct_parse_block(datasets_check('BoxOffset'), blk) == 0) then
        do i = 1, sb%dim
          call loct_parse_block_float(blk, 0, i-1, sb%box_offset(i))
        end do
        call loct_parse_block_end(blk)
      else
        call loct_parse_float(datasets_check('BoxOffset'), M_ZERO, sb%box_offset(1))
        sb%box_offset(1:sb%dim) = sb%box_offset(1)
      end if
      sb%box_offset(:) = sb%box_offset(:)*units_inp%length%factor

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


    !--------------------------------------------------------------
    subroutine build_lattice()
      type(block_t) :: blk
      FLOAT :: norm
      integer :: idim, jdim

      call push_sub('simul_box.simul_box_init.build_lattice')

      !%Variable LatticeVectors
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% <tt>%Spacing
      !% <br>&nbsp;&nbsp;1.0 | 0.0 | 0.0
      !% <br>&nbsp;&nbsp;0.0 | 1.0 | 0.0
      !% <br>&nbsp;&nbsp;0.0 | 0.0 | 1.0
      !% <br>%</tt>
      !%End
      
      ! this has to be updated for non-orthogonal grids
      sb%rcell_volume = product(M_TWO*sb%lsize(1:sb%periodic_dim))

      sb%rlattice = M_ZERO
      do idim = 1, MAX_DIM
        sb%rlattice(idim, idim) = M_ONE
      end do

      if (loct_parse_block(datasets_check('LatticeVectors'), blk) == 0) then 
        do idim = 1, sb%dim
          do jdim = 1, sb%dim
            call loct_parse_block_float(blk, idim - 1,  jdim - 1, sb%rlattice(jdim, idim))
          end do
        end do
      end if

      do idim = 1, sb%dim
        norm = sqrt(sum(sb%rlattice(1:sb%dim, idim)**2))
        do jdim = 1, sb%dim
          sb%rlattice(jdim, idim) = sb%rlattice(jdim, idim) / norm
        end do
      end do

      call reciprocal_lattice(sb%rlattice, sb%klattice_unitary, sb%volume_element)
      
      sb%klattice = M_ZERO
      do idim = 1, sb%periodic_dim
        sb%klattice(1:MAX_DIM, idim) = sb%klattice_unitary(1:MAX_DIM, idim)*M_TWO*M_PI/(M_TWO*sb%lsize(idim))
      end do

      call pop_sub()
    end subroutine build_lattice
  end subroutine simul_box_init

  
  !--------------------------------------------------------------
  ! Read the coordinates of the leads atoms and add the to the
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
      ALLOCATE(lead_geo(NLEADS), NLEADS)
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
      ALLOCATE(geo%atom(geo%natoms), geo%natoms)
      if(geo%ncatoms.gt.0) then
        SAFE_DEALLOCATE_P(geo%catom)
      end if
      geo%ncatoms = central_geo%ncatoms +                &
        sb%add_unit_cells(LEFT)*lead_geo(LEFT)%ncatoms + &
        sb%add_unit_cells(RIGHT)*lead_geo(RIGHT)%ncatoms
      ALLOCATE(geo%catom(geo%ncatoms), geo%ncatoms)

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
      deallocate(lead_geo)

    end if
    call pop_sub()
  end subroutine simul_box_add_lead_atoms


  !--------------------------------------------------------------
  subroutine simul_box_atoms_in_box(sb, geo)
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(inout) :: geo

    ! this function checks that the atoms are inside the box. if not:
    ! if the system is periodic, the atoms are moved inside the box;
    ! if the system is finite, a warning is emitted.

    integer :: iatom, pd, idir
    FLOAT :: xx(1:MAX_DIM)

    call push_sub('simul_box.simul_box_atoms_in_box')

    pd = sb%periodic_dim

    do iatom = 1, geo%natoms

      if (simul_box_is_periodic(sb)) then

        !convert the position to the orthogonal space
        xx(1:pd) = matmul(geo%atom(iatom)%x(1:pd) - sb%box_offset(1:pd), sb%klattice_unitary(1:pd, 1:pd))

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

        geo%atom(iatom)%x(1:pd) = matmul(sb%klattice_unitary(1:pd, 1:pd), xx(1:pd) + sb%box_offset(1:pd))

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
      message(1) = "Your lattice vectors form a left-handed system"
      call write_fatal(1)
    end if

    kv(1:3, 1) = dcross_product(rv(:, 2), rv(:, 3))/volume
    kv(1:3, 2) = dcross_product(rv(:, 3), rv(:, 1))/volume
    kv(1:3, 3) = dcross_product(rv(:, 1), rv(:, 2))/volume    

    call pop_sub()
  end subroutine reciprocal_lattice


  !--------------------------------------------------------------
  recursive subroutine simul_box_end(sb)
    type(simul_box_t), intent(inout) :: sb    

    integer :: il

    call push_sub('simul_box.simul_box_end')

    if(sb%open_boundaries) then
      do il = 1, NLEADS
        call simul_box_end(sb%lead_unit_cell(il))
      end do
      SAFE_DEALLOCATE_P(sb%lead_unit_cell)
      nullify(sb%lead_unit_cell)
    end if

    call pop_sub()
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
      'image defined ', &
      'hypercube     '/)

    integer :: ii

    call push_sub('simul_box.simul_box_write_info')

    write(message(1),'(a)') 'Simulation Box:'
    if(sb%box_shape.eq.BOX_USDEF) then
      write(message(2), '(a)') '  Type = user-defined'
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

    write(message(1), '(a,i1,a)') '  Octopus will run in ', sb%dim, ' dimension(s).'
    write(message(2), '(a,i1,a)') '  Octopus will treat the system as periodic in ', &
      sb%periodic_dim, ' dimension(s).'
    call write_info(2, iunit)

    if(sb%periodic_dim > 0 .or. sb%box_shape == PARALLELEPIPED) then
      write(message(1),'(1x)')
      write(message(2),'(a,3a,a)') '  Lattice Vectors [', trim(units_out%length%abbrev), ']'
      do ii = 1, sb%dim
        write(message(2+ii),'(9f12.6)') sb%rlattice(1:sb%dim,ii)*M_TWO*sb%lsize(ii)/units_out%length%factor
      end do
      call write_info(2+sb%dim, iunit)

      write(message(1),'(a,3a,a)') '  Reciprocal Lattice Vectors [', trim(units_out%length%abbrev), '^-1]'
      do ii = 1, sb%dim
        write(message(1+ii),'(3f12.6)')   sb%klattice(1:sb%dim,ii)*units_out%length%factor
      end do
      call write_info(1+sb%dim, iunit)

      write(message(1),'(a,f12.6,3a)') &
        '  Volume element = ', sb%volume_element/units_out%length%factor**sb%dim, &
        ' [', trim(units_out%length%abbrev), ']'
      call write_info(1, iunit)
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
  logical function simul_box_in_box(sb, geo, x, inner_box) result(in_box)
    type(simul_box_t),  intent(in) :: sb
    type(geometry_t),   intent(in) :: geo
    FLOAT,              intent(in) :: x(:)
    logical, optional,  intent(in) :: inner_box

    real(8), parameter :: DELTA = CNST(1e-12)
    FLOAT :: r, re, im, xx(MAX_DIM)
    real(8) :: llimit(MAX_DIM), ulimit(MAX_DIM)
    real(8) :: factor

#if defined(HAVE_GDLIB)
    integer :: red, green, blue, ix, iy
#endif

    in_box = .true.

    factor = M_ONE
    if(present(inner_box)) then
      if(inner_box) factor = sb%inner_size
    end if

    xx(1:sb%dim) = x(1:sb%dim) - sb%box_offset(1:sb%dim)

    !convert to the orthogonal space
    xx(1:sb%dim) = matmul(xx(1:sb%dim), sb%klattice_unitary(1:sb%dim, 1:sb%dim))

    xx(1:sb%dim) = xx(1:sb%dim)/factor

    select case(sb%box_shape)
    case(SPHERE)
      in_box = (sqrt(sum(xx(1:sb%dim)**2)) <= sb%rsize+DELTA)

    case(CYLINDER)
      r = sqrt(xx(2)**2 + xx(3)**2)
      in_box = (r<=sb%rsize+DELTA .and. abs(xx(1)) <= sb%xsize+DELTA)

    case(MINIMUM)
      in_box = in_minimum()

    case(PARALLELEPIPED, HYPERCUBE) 
      llimit(1:sb%dim) = -sb%lsize(1:sb%dim) - DELTA
      ulimit(1:sb%dim) =  sb%lsize(1:sb%dim) + DELTA
      ulimit(1:sb%periodic_dim) = sb%lsize(1:sb%periodic_dim) - DELTA
      if(sb%open_boundaries) then
        ulimit(TRANS_DIR) = sb%lsize(TRANS_DIR) - DELTA
      end if

      in_box = all(xx(1:sb%dim) >= llimit(1:sb%dim) .and. xx(1:sb%dim) <= ulimit(1:sb%dim))

#if defined(HAVE_GDLIB)
    case(BOX_IMAGE)
      ix = int((xx(1) + sb%lsize(1))/sb%h(1))
      iy = int((xx(2) + sb%lsize(2))/sb%h(2))
      call loct_gdimage_get_pixel_rgb(sb%image, ix, iy, red, green, blue)
      in_box = (red == 255).and.(green == 255).and.(blue == 255)
#endif

    case(BOX_USDEF)
      ! is it inside the user-given boundaries
      in_box =  &
        (xx(1) >= -sb%lsize(1)-DELTA.and.xx(1) <= sb%lsize(1)+DELTA).and. &
        (xx(2) >= -sb%lsize(2)-DELTA.and.xx(2) <= sb%lsize(2)+DELTA).and. &
        (xx(3) >= -sb%lsize(3)-DELTA.and.xx(3) <= sb%lsize(3)+DELTA)
      
      ! and inside the simulation box
      xx(:) = xx(:)/units_inp%length%factor ! convert from a.u. to input units
      r = sqrt(sum(xx(:)**2))
      call loct_parse_expression(re, im, sb%dim, xx, r, M_ZERO, sb%user_def)
      in_box = in_box .and. (re .ne. M_ZERO)
    end select

  contains

    !--------------------------------------------------------------
    logical function in_minimum()
      integer :: i
      FLOAT :: radius

      in_minimum = .false.
      do i = 1, geo%natoms
        r = sqrt(sum((xx(1:sb%dim) - geo%atom(i)%x(1:sb%dim))**2))
        if(sb%rsize > M_ZERO) then
          radius = sb%rsize
        else
          radius = geo%atom(i)%spec%def_rsize
        endif
        if(r <= radius + DELTA) then
          in_minimum = .true.
          exit
        end if
      end do

    end function in_minimum
  end function simul_box_in_box


  !--------------------------------------------------------------
  logical pure function simul_box_is_periodic(sb)
    type(simul_box_t), intent(in) :: sb

    simul_box_is_periodic = sb%periodic_dim > 0

  end function simul_box_is_periodic


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
    write(iunit, '(a20,e22.14)')    'inner_size=         ', sb%inner_size
    do i = 1, MAX_DIM
      write(iunit, '(a9,i1,a11,9e22.14)')   'rlattice(', i, ')=         ', sb%rlattice(1:MAX_DIM, i)
    end do
    write(iunit, '(a20,l7)')        'open_boundaries=    ', sb%open_boundaries
    if(sb%open_boundaries) then
      write(iunit, '(a20,2i4)')     'add_unit_cells=     ', sb%add_unit_cells
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

    iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/mesh', action='read', is_tmp=.true.)
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
      call write_fatal(1)
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
    integer                       :: il
    
    call push_sub('simul_box.lead_unit_cell_extent')
    
    if(sb%open_boundaries) then
      lead_unit_cell_extent = int(2*sb%lead_unit_cell(il)%lsize(TRANS_DIR)/sb%h(TRANS_DIR))
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
    character(len=100) :: line
    integer            :: idim, jdim, il
    FLOAT              :: norm

    call push_sub('simul_box.simul_box_init_from_file')

    ! Find (and throw away) the dump tag.
    do
      read(iunit, '(a)') line
      if(trim(line).eq.dump_tag) exit
    end do

    read(iunit, *) str, sb%box_shape
    read(iunit, *) str, sb%dim
    read(iunit, *) str, sb%periodic_dim
    select case(sb%box_shape)
    case(SPHERE, MINIMUM)
      read(iunit, *) str, sb%rsize
      read(iunit, *) str, sb%lsize(1:MAX_DIM)
    case(CYLINDER)
      read(iunit, *) str, sb%rsize
      read(iunit, *) str, sb%xsize
      read(iunit, *) str, sb%lsize(1:MAX_DIM)
    case(PARALLELEPIPED)
      read(iunit, *) str, sb%lsize(1:MAX_DIM)
    case(BOX_USDEF)
      read(iunit, *) str, sb%lsize(1:MAX_DIM)
      read(iunit, *) str, sb%user_def
    end select
    read(iunit, *) str, sb%fft_alpha
    read(iunit, *) str, sb%h(1:MAX_DIM)
    read(iunit, *) str, sb%box_offset(1:MAX_DIM)
    read(iunit, *) str, sb%inner_size
    do idim=1, MAX_DIM
      read(iunit, *) str, sb%rlattice(1:MAX_DIM, idim)
    end do

    do idim = 1, sb%dim
      norm = lalg_nrm2(sb%dim, sb%rlattice(1:sb%dim, idim))
      do jdim = 1, sb%dim
        sb%rlattice(jdim, idim) = sb%rlattice(jdim, idim) / norm
      end do
    end do

    call reciprocal_lattice(sb%rlattice, sb%klattice_unitary, sb%volume_element)

    sb%klattice = M_ZERO
    do idim = 1, sb%periodic_dim
      sb%klattice(:, idim) = sb%klattice_unitary(:, idim)*M_TWO*M_PI/(M_TWO*sb%lsize(idim))
    end do

    read(iunit, *) str, sb%open_boundaries
    if(sb%open_boundaries) then
      read(iunit, *) str, sb%add_unit_cells(LEFT), sb%add_unit_cells(RIGHT)
      read(iunit, *) str, sb%lead_restart_dir(LEFT)
      read(iunit, *) str, sb%lead_restart_dir(RIGHT)
      ALLOCATE(sb%lead_unit_cell(NLEADS), NLEADS)
      do il = 1, NLEADS
        call read_lead_unit_cell(sb, il)
      end do
    end if

    call pop_sub()
  end subroutine simul_box_init_from_file


  !--------------------------------------------------------------
  recursive logical function simul_box_is_eq(sb1, sb2) result(res)
    type(simul_box_t), intent(in) :: sb1, sb2

    integer :: il

    call push_sub('simul_box.simul_box_is_eq')

    res = .false.
    if(sb1%box_shape .ne. sb2%box_shape)             go to 1
    if(sb1%dim .ne. sb2%dim)                         go to 1
    if(sb1%periodic_dim .ne. sb2%periodic_dim)       go to 1

    select case(sb1%box_shape)
    case(SPHERE, MINIMUM)
      if(.not.(sb1%rsize .app. sb2%rsize))           go to 1
      if(.not.(sb1%lsize .app. sb2%lsize))           go to 1
    case(CYLINDER)
      if(.not.(sb1%rsize .app. sb2%rsize))           go to 1
      if(.not.(sb1%xsize .app. sb2%xsize))           go to 1
      if(.not.(sb1%lsize .app. sb2%lsize))           go to 1
    case(PARALLELEPIPED)
      if(.not.(sb1%lsize .app. sb2%lsize))           go to 1
    case(BOX_USDEF)
      if(.not.(sb1%lsize .app. sb2%lsize))           go to 1
      if(trim(sb1%user_def) .ne. trim(sb2%user_def)) go to 1
    end select

    if(.not.(sb1%fft_alpha .app. sb2%fft_alpha))     go to 1
    if(.not.(sb1%h .app. sb2%h))                     go to 1
    if(.not.(sb1%box_offset .app. sb2%box_offset))   go to 1
    if(sb1%open_boundaries.neqv.sb2%open_boundaries)   go to 1

    if (sb1%open_boundaries) then
      if(any(sb1%add_unit_cells.ne.sb2%add_unit_cells))     go to 1
      if(any(sb1%lead_restart_dir.ne.sb2%lead_restart_dir)) go to 1
      if(associated(sb1%lead_unit_cell).and.associated(sb2%lead_unit_cell)) then
        do il = 1, NLEADS
          if(.not.simul_box_is_eq(sb1%lead_unit_cell(il), sb2%lead_unit_cell(il))) go to 1
        end do
      end if
    endif
    res = .true.

1   call pop_sub()
  end function simul_box_is_eq


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
    sbout%klattice                = sbin%klattice
    sbout%klattice_unitary        = sbin%klattice_unitary
    sbout%volume_element          = sbin%volume_element
    sbout%fft_alpha               = sbin%fft_alpha
    sbout%dim                     = sbin%dim
    sbout%periodic_dim            = sbin%periodic_dim
    sbout%open_boundaries         = sbin%open_boundaries
    sbout%add_unit_cells          = sbin%add_unit_cells
    sbout%lead_restart_dir        = sbin%lead_restart_dir
    sbout%inner_size              = sbin%inner_size
    if(associated(sbin%lead_unit_cell)) then
      ALLOCATE(sbout%lead_unit_cell(NLEADS), NLEADS)
      do il = 1, NLEADS
        call simul_box_copy(sbout%lead_unit_cell(il), sbin%lead_unit_cell(il))
      end do
    end if

    call pop_sub()
  end subroutine simul_box_copy

  logical pure function simul_box_multires(this) result(mr)
    type(simul_box_t), intent(in) :: this

    mr = abs(this%inner_size - M_ONE) > M_EPSILON
  end function simul_box_multires

end module simul_box_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
