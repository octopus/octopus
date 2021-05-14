!! Copyright (C) 2020 N. Tancogne-Dejean
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

module external_potential_oct_m
  use clock_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use iihash_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use interaction_partner_oct_m
  use io_function_oct_m
  use linked_list_oct_m
  use lorentz_force_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  implicit none

  private
  public ::               &
    external_potential_t, &
    load_external_potentials

  type, extends(interaction_partner_t) :: external_potential_t
    private

    integer, public :: type   !< Type of external potential

    character(len=1024) :: potential_formula !< for the user-defined potential
    character(len=200)  :: density_formula   !< If we have a charge distribution creating the potential
    character(len=MAX_PATH_LEN) :: filename  !< for the potential read from a file.
    FLOAT :: omega

    FLOAT, allocatable, public :: pot(:)

    FLOAT, allocatable, public :: B_field(:)           !< static magnetic field
    integer                    :: gauge_2D
    FLOAT, allocatable, public :: A_static(:,:)        !< static vector potential
    FLOAT, allocatable, public :: E_field(:)           !< static electric field
    !Auxiliary arrays for the electrons only
    !TODO: Suppress once electrons fully use the new framework
    FLOAT, allocatable, public :: v_ext(:)             !< static scalar potential - 1:gr%mesh%np_part

  contains
    procedure :: calculate => external_potential_calculate
    procedure :: allocate_memory => external_potential_allocate
    procedure :: deallocate_memory => external_potential_deallocate
    procedure :: update_exposed_quantities => external_potential_update_exposed_quantities
    procedure :: update_exposed_quantity => external_potential_update_exposed_quantity
    procedure :: init_interaction_as_partner => external_potential_init_interaction_as_partner
    procedure :: copy_quantities_to_interaction => external_potential_copy_quantities_to_interaction
    final :: external_potential_finalize
  end type external_potential_t

  integer, public, parameter ::  &
    EXTERNAL_POT_USDEF          = 201,           & !< user-defined function for local potential
    EXTERNAL_POT_FROM_FILE      = 202,           &
    EXTERNAL_POT_CHARGE_DENSITY = 203,           & !< user-defined function for charge density
    EXTERNAL_POT_STATIC_BFIELD  = 204,           &  !< Static magnetic field
    EXTERNAL_POT_STATIC_EFIELD  = 205
 

  interface external_potential_t
    module procedure external_potential_init
  end interface external_potential_t

contains

  function external_potential_init(namespace) result(this)
    class(external_potential_t), pointer :: this
    type(namespace_t), intent(in) :: namespace

    PUSH_SUB(external_potential_init)

    SAFE_ALLOCATE(this)

    this%namespace = namespace_t("ExternalPotential", parent=namespace)
    call space_init(this%space, this%namespace)

    this%type = -1

    ! Initialize clock without a time-step, as the potential will not be propagated
    this%clock = clock_t()

    POP_SUB(external_potential_init)
  end function external_potential_init

  ! ---------------------------------------------------------
  subroutine external_potential_finalize(this)
    type(external_potential_t), intent(inout) :: this

    PUSH_SUB(external_potential_finalize)

    call this%deallocate_memory()

    POP_SUB(external_potential_finalize)
  end subroutine external_potential_finalize

  ! ---------------------------------------------------------
  subroutine external_potential_allocate(this, mesh)
    class(external_potential_t), intent(inout) :: this
    type(mesh_t),                intent(in)    :: mesh

    PUSH_SUB(external_potential_allocate)

    select case(this%type)
    case(EXTERNAL_POT_USDEF, EXTERNAL_POT_FROM_FILE, EXTERNAL_POT_CHARGE_DENSITY)
      SAFE_ALLOCATE(this%pot(1:mesh%np))
    case(EXTERNAL_POT_STATIC_BFIELD)
      SAFE_ALLOCATE(this%A_static(1:mesh%np, 1:this%space%dim))
    case(EXTERNAL_POT_STATIC_EFIELD)
      if (this%space%periodic_dim < this%space%dim) then
        SAFE_ALLOCATE(this%pot(1:mesh%np))
        SAFE_ALLOCATE(this%v_ext(1:mesh%np_part))
      end if
    end select

    POP_SUB(external_potential_allocate)
  end subroutine external_potential_allocate

  ! ---------------------------------------------------------
  subroutine external_potential_deallocate(this)
    class(external_potential_t), intent(inout) :: this

    PUSH_SUB(external_potential_deallocate)

    SAFE_DEALLOCATE_A(this%pot)
    SAFE_DEALLOCATE_A(this%B_field)
    SAFE_DEALLOCATE_A(this%A_static)
    SAFE_DEALLOCATE_A(this%E_field)
    SAFE_DEALLOCATE_A(this%v_ext)

    POP_SUB(external_potential_deallocate)
  end subroutine external_potential_deallocate

  ! ---------------------------------------------------------
  logical function external_potential_update_exposed_quantities(partner, requested_time, interaction) &
    result(allowed_to_update)
    class(external_potential_t), intent(inout) :: partner
    type(clock_t),               intent(in)    :: requested_time
    class(interaction_t),        intent(inout) :: interaction

    PUSH_SUB(external_potential_update_exposed_quantities)

    ! Always allowed to update, as the external potentials are not propagated
    allowed_to_update = .true.

    call partner%clock%set_time(requested_time)

    select type (interaction)
    type is (ghost_interaction_t)
      ! Nothing to copy. We still need to check that we are at the right
      ! time for the update though!
    class default
      call partner%copy_quantities_to_interaction(interaction)
    end select

    POP_SUB(external_potential_update_exposed_quantities)

  end function external_potential_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine external_potential_update_exposed_quantity(partner, iq)
    class(external_potential_t),      intent(inout) :: partner
    integer,                          intent(in)    :: iq

    PUSH_SUB(external_potential_update_exposed_quantities)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(external_potential_update_exposed_quantities)

  end subroutine external_potential_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine external_potential_init_interaction_as_partner(partner, interaction)
    class(external_potential_t),     intent(in)    :: partner
    class(interaction_t),            intent(inout) :: interaction

    PUSH_SUB(external_potential_init_interaction_as_partner)

    select type (interaction)
    type is (lorentz_force_t)
      ! Nothing to be initialized for the Lorentz force.
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(external_potential_init_interaction_as_partner)
  end subroutine external_potential_init_interaction_as_partner

    ! ---------------------------------------------------------
  subroutine external_potential_copy_quantities_to_interaction(partner, interaction)
    class(external_potential_t),     intent(inout) :: partner
    class(interaction_t),            intent(inout) :: interaction

    integer :: ip

    PUSH_SUB(external_potential_copy_quantities_to_interaction)

    select type (interaction)
    type is (lorentz_force_t)
      if (partner%type == EXTERNAL_POT_STATIC_EFIELD) then
        do ip = 1, interaction%system_np
          interaction%partner_E_field(:, ip) = partner%E_field
          interaction%partner_B_field(:, ip) = M_ZERO
        end do
      else if (partner%type == EXTERNAL_POT_STATIC_BFIELD) then 
        do ip = 1, interaction%system_np
          interaction%partner_E_field(:, ip) = M_ZERO
          interaction%partner_B_field(:, ip) = partner%B_field
        end do
      else
        ASSERT(.false.) !This should never occur.
      end if 

    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(external_potential_copy_quantities_to_interaction)
  end subroutine external_potential_copy_quantities_to_interaction


  ! ---------------------------------------------------------
  subroutine external_potential_calculate(this, namespace, mesh, poisson)
    class(external_potential_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace
    type(mesh_t),                intent(in)    :: mesh
    type(poisson_t),             intent(in)    :: poisson

    FLOAT :: pot_re, pot_im, r, xx(1:MAX_DIM)
    FLOAT, allocatable :: den(:), grx(:)
    integer :: ip, err

    PUSH_SUB(external_potential_calculate)

    select case(this%type)

    case(EXTERNAL_POT_USDEF)
      ASSERT(allocated(this%pot))

      do ip = 1, mesh%np
        call mesh_r(mesh, ip, r, coords = xx)
        call parse_expression(pot_re, pot_im, this%space%dim, xx, r, M_ZERO, this%potential_formula)
        this%pot(ip) = pot_re
      end do 

    case(EXTERNAL_POT_FROM_FILE)
      ASSERT(allocated(this%pot))

      call dio_function_input(trim(this%filename), namespace, this%space, mesh, this%pot, err)
      if(err /= 0) then
        write(message(1), '(a)')    'Error loading file '//trim(this%filename)//'.'
        write(message(2), '(a,i4)') 'Error code returned = ', err      
        call messages_fatal(2, namespace=namespace)
      end if

    case(EXTERNAL_POT_CHARGE_DENSITY)
      ASSERT(allocated(this%pot))

      SAFE_ALLOCATE(den(1:mesh%np))

      do ip = 1, mesh%np
        call mesh_r(mesh, ip, r, coords = xx)
        call parse_expression(pot_re, pot_im, this%space%dim, xx, r, M_ZERO, this%potential_formula)
        den(ip) = pot_re
      end do

      if(poisson_solver_is_iterative(poisson)) then
        ! pot has to be initialized before entering routine
        ! and our best guess for the potential is zero
        this%pot(1:mesh%np) = M_ZERO
      end if
      call dpoisson_solve(poisson, this%pot, den, all_nodes = .false.)

      SAFE_DEALLOCATE_A(den)

    case(EXTERNAL_POT_STATIC_BFIELD)
      ASSERT(allocated(this%B_field))

      ! Compute the vector potential
      SAFE_ALLOCATE(grx(1:this%space%dim))

      select case(this%space%dim)
      case(2)
        select case(this%gauge_2d)
        case(0) ! linear_xy
          if (this%space%periodic_dim == 1) then
            message(1) = "For 2D system, 1D-periodic, StaticMagneticField can only be "
            message(2) = "applied for StaticMagneticField2DGauge = linear_y."
            call messages_fatal(2, namespace=namespace)
          end if
          do ip = 1, mesh%np
            grx(1:this%space%dim) = mesh%x(ip, 1:this%space%dim)
            this%A_static(ip, :) = M_HALF/P_C*(/grx(2), -grx(1)/) * this%B_field(3)
          end do
        case(1) ! linear y
          do ip = 1, mesh%np
            grx(1:this%space%dim) = mesh%x(ip, 1:this%space%dim)
            this%A_static(ip, :) = M_ONE/P_C*(/grx(2), M_ZERO/) * this%B_field(3)
          end do
      end select
    case(3)
      do ip = 1, mesh%np
        grx(1:this%space%dim) = mesh%x(ip, 1:this%space%dim)
        this%A_static(ip, :) = M_HALF/P_C*(/grx(2) * this%B_field(3) - grx(3) * this%B_field(2), &
                               grx(3) * this%B_field(1) - grx(1) * this%B_field(3), &
                               grx(1) * this%B_field(2) - grx(2) * this%B_field(1)/)
      end do
    end select

    SAFE_DEALLOCATE_A(grx)

    case(EXTERNAL_POT_STATIC_EFIELD)
      ASSERT(allocated(this%E_field))

      if (this%space%periodic_dim < this%space%dim) then
        ! Compute the scalar potential
        !
        ! Note that the -1 sign is missing. This is because we
        ! consider the electrons with +1 charge. The electric field
        ! however retains the sign because we also consider protons to
        ! have +1 charge when calculating the force.
        !
        ! NTD: This comment is very confusing and prone to error
        ! TODO: Fix this to have physically sound quantities and interactions
        do ip = 1, mesh%np
          this%pot(ip) = sum(mesh%x(ip, this%space%periodic_dim + 1:this%space%dim) &
                                    * this%E_field(this%space%periodic_dim + 1:this%space%dim))
        end do
        ! The following is needed to make interpolations.
        ! It is used by PCM.
        this%v_ext(1:mesh%np) = this%pot(1:mesh%np)
        do ip = mesh%np+1, mesh%np_part
          this%v_ext(ip) = sum(mesh%x(ip, this%space%periodic_dim + 1:this%space%dim) &
                                 * this%E_field(this%space%periodic_dim + 1:this%space%dim))
        end do
      end if

    end select

    POP_SUB(external_potential_calculate)
  end subroutine external_potential_calculate

  subroutine load_external_potentials(external_potentials, namespace)
    class(partner_list_t), intent(inout)  :: external_potentials
    type(namespace_t),    intent(in)     :: namespace

    integer :: n_pot_block, row, read_data
    type(block_t) :: blk
    class(external_potential_t), pointer :: pot

    integer :: dim, periodic_dim, idir

    PUSH_SUB(load_external_potentials)

    !%Variable StaticExternalPotentials
    !%Type block
    !%Section System
    !%Description
    !% An static external potential is a model potential added to the local potential of the Hamiltonian
    !%
    !% The format of this block is the following: 
    !% The first field defines the type of species (the valid options are detailed
    !% below).
    !%
    !% Then a list of parameters follows. The parameters are specified
    !% by a first field with the parameter name and the field that
    !% follows with the value of the parameter. Some parameters are
    !% specific to a certain species while others are accepted by all
    !% species. These are <tt>mass</tt>, <tt>max_spacing</tt>, and <tt>min_radius</tt>.
    !%
    !% These are examples of possible species:
    !%
    !% <tt>%ExternalPotential
    !% <br>&nbsp;&nbsp; potential_user_defined | potential_formula | "1/2*r^2"
    !% <br>%</tt>
    !%Option potential_from_file  -202
    !% The potential is read from a file. Accepted file formats, detected by extension: obf, ncdf and csv.
    !%Option potential_user_defined -201
    !% Species with user-defined potential. The potential for the
    !% species is defined by the formula given by the <tt>potential_formula</tt>
    !% parameter.
    !%Option potential_charge_density -203
    !% The potential for this species is created from the distribution
    !% of charge given by the <tt>density_formula</tt> parameter.
    !%Option file -10010
    !% The path for the file that describes the species.
    !%Option potential_formula -10012
    !% Mathematical expression that defines the potential for <tt>species_user_defined</tt>. You can use
    !% any of the <i>x</i>, <i>y</i>, <i>z</i> or <i>r</i> variables.
    !%Option density_formula -10013
    !% Mathematical expression that defines the charge density for <tt>species_charge_density</tt>. You can use
    !% any of the <i>x</i>, <i>y</i>, <i>z</i> or <i>r</i> variables.
    !%End

    ! First, find out if there is a Species block.
    n_pot_block = 0
    if(parse_block(namespace, 'StaticExternalPotentials', blk) == 0) then
      n_pot_block = parse_block_n(blk)

      do row = 0, n_pot_block-1
        !Create a potential
        pot => external_potential_t(namespace) 
        !Parse the information from the block
        call read_from_block(pot, namespace, blk, row, read_data)
        ASSERT(read_data > 0)
        !Add this to the list
        call external_potentials%add(pot)  
      end do 
      call parse_block_end(blk)
    end if

    
    !Here I am parsing the variables Dimensions et PeriodicDimensions because we do not have access 
    !to this information here.
    !TODO: This needs to be removed and replaced by something better
    call parse_variable(namespace, 'Dimensions', 3, dim) 
    call parse_variable(namespace, 'PeriodicDimensions', 0, periodic_dim)


    !%Variable StaticMagneticField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% A static constant magnetic field may be added to the usual Hamiltonian,
    !% by setting the block <tt>StaticMagneticField</tt>.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the magnetic field vector. Note that
    !% if you are running the code in 1D mode, this will not work, and if you
    !% are running the code in 2D mode the magnetic field will have to be in
    !% the <i>z</i>-direction, so that the first two columns should be zero.
    !% Possible in periodic system only in these cases: 2D system, 1D periodic,
    !% with <tt>StaticMagneticField2DGauge = linear_y</tt>;
    !% 3D system, 1D periodic, field is zero in <i>y</i>- and <i>z</i>-directions (given
    !% currently implemented gauges).
    !%
    !% The magnetic field should always be entered in atomic units, regardless
    !% of the <tt>Units</tt> variable. Note that we use the "Gaussian" system
    !% meaning 1 au[B] = <math>1.7152553 \times 10^7</math> Gauss, which corresponds to
    !% <math>1.7152553 \times 10^3</math> Tesla.
    !%End
    if(parse_block(namespace, 'StaticMagneticField', blk) == 0) then
      !Create a potential
      pot => external_potential_t(namespace)
      pot%type = EXTERNAL_POT_STATIC_BFIELD
      call pot%supported_interactions_as_partner%add(LORENTZ_FORCE)

      !%Variable StaticMagneticField2DGauge
      !%Type integer
      !%Default linear_xy
      !%Section Hamiltonian
      !%Description
      !% The gauge of the static vector potential <math>A</math> when a magnetic field
      !% <math>B = \left( 0, 0, B_z \right)</math> is applied to a 2D-system.
      !%Option linear_xy 0
      !% Linear gauge with <math>A = \frac{1}{2c} \left( -y, x \right) B_z</math>. (Cannot be used for periodic systems.)
      !%Option linear_y 1
      !% Linear gauge with <math>A = \frac{1}{c} \left( -y, 0 \right) B_z</math>. Can be used for <tt>PeriodicDimensions = 1</tt>
      !% but not <tt>PeriodicDimensions = 2</tt>.
      !%End
      call parse_variable(namespace, 'StaticMagneticField2DGauge', 0, pot%gauge_2d)
      if(.not.varinfo_valid_option('StaticMagneticField2DGauge', pot%gauge_2d)) &
        call messages_input_error(namespace, 'StaticMagneticField2DGauge')

      SAFE_ALLOCATE(pot%B_field(1:3))
      do idir = 1, 3
        call parse_block_float(blk, 0, idir - 1, pot%B_field(idir))
      end do
      select case(dim)
      case(1)
        call messages_input_error(namespace, 'StaticMagneticField')
      case(2)
        if(periodic_dim == 2) then
          message(1) = "StaticMagneticField cannot be applied in a 2D, 2D-periodic system."
          call messages_fatal(1, namespace=namespace)
        end if
        if(pot%B_field(1)**2 + pot%B_field(2)**2 > M_ZERO) then
          call messages_input_error(namespace, 'StaticMagneticField')
        end if
      case(3)
        ! Consider cross-product below: if grx(1:this%space%periodic_dim) is used, it is not ok.
        ! Therefore, if idir is periodic, B_field for all other directions must be zero.
        ! 1D-periodic: only Bx. 2D-periodic or 3D-periodic: not allowed. Other gauges could allow 2D-periodic case.
        if(periodic_dim >= 2) then
          message(1) = "In 3D, StaticMagneticField cannot be applied when the system is 2D- or 3D-periodic."
          call messages_fatal(1, namespace=namespace)
        else if(periodic_dim == 1 .and. any(abs(pot%B_field(2:3)) > M_ZERO)) then
          message(1) = "In 3D, 1D-periodic, StaticMagneticField must be zero in the y- and z-directions."
          call messages_fatal(1, namespace=namespace)
        end if
      end select
      call parse_block_end(blk)

      if(dim > 3) call messages_not_implemented('Magnetic field for dim > 3', namespace=namespace)

      !Add this to the list
      call external_potentials%add(pot)      

      !The corresponding A field on the mesh is computed in the routine external_potential_calculate 

    end if

    !%Variable StaticElectricField
    !%Type block
    !%Default 0
    !%Section Hamiltonian
    !%Description
    !% A static constant electric field may be added to the usual Hamiltonian,
    !% by setting the block <tt>StaticElectricField</tt>.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the electric field vector.
    !% It can be applied in a periodic direction of a large supercell via
    !% the single-point Berry phase.
    !%End
    if(parse_block(namespace, 'StaticElectricField', blk)==0) then
      !Create a potential
      pot => external_potential_t(namespace) 
      pot%type = EXTERNAL_POT_STATIC_EFIELD
      call pot%supported_interactions_as_partner%add(LORENTZ_FORCE)
   
      SAFE_ALLOCATE(pot%E_field(1:dim))
      do idir = 1, dim
        call parse_block_float(blk, 0, idir - 1, pot%E_field(idir), units_inp%energy / units_inp%length)
  
        !Electron-specific checks (k-points) are done in the hamiltonian_elec.F90 file
        if(idir <= periodic_dim .and. abs(pot%E_field(idir)) > M_EPSILON) then
          message(1) = "Applying StaticElectricField in a periodic direction is only accurate for large supercells."
          call messages_warning(1, namespace=namespace)
        end if
      end do
      call parse_block_end(blk)
  
      !Add this to the list
      call external_potentials%add(pot)
  
      !The corresponding A field on the mesh is computed in the routine external_potential_calculate
    end if


    POP_SUB(load_external_potentials)
  end subroutine load_external_potentials

  ! ---------------------------------------------------------
  subroutine read_from_block(pot, namespace, blk, row, read_data)
    type(external_potential_t), intent(inout) :: pot
    type(namespace_t),          intent(in)    :: namespace
    type(block_t),              intent(in)    :: blk
    integer,                    intent(in)    :: row
    integer,                    intent(out)   :: read_data

    integer :: ncols, icol, flag
    type(iihash_t) :: read_parameters


    PUSH_SUB(read_from_block)

    ncols = parse_block_cols(blk, row)
    read_data = 0

    call parse_block_integer(blk, row, 0, pot%type)

    ! To detect the old species block format, options are represented
    ! as negative values. If we get a non-negative value we know we
    ! are reading a mass.
    if(pot%type >= 0) then
      message(1) = 'Error in reading the ExternalPotentials block'
      call messages_fatal(1, namespace=namespace)
    end if

    ! now we convert back to positive
    pot%type = -pot%type

    read_data = 1

    if(pot%type /= EXTERNAL_POT_CHARGE_DENSITY .and. pot%type /= EXTERNAL_POT_USDEF .and. pot%type /= EXTERNAL_POT_FROM_FILE) then
      call messages_input_error(namespace, 'ExternalPotentials', "Unknown type of external potential")
    end if
    
    call iihash_init(read_parameters)
    
    icol = read_data
    do
      if(icol >= ncols) exit

      call parse_block_integer(blk, row, icol, flag)
      
      select case(flag)

      case(OPTION__STATICEXTERNALPOTENTIALS__FILE)
        call check_duplication(OPTION__STATICEXTERNALPOTENTIALS__FILE)
        call parse_block_string(blk, row, icol + 1, pot%filename)

      case(OPTION__STATICEXTERNALPOTENTIALS__POTENTIAL_FORMULA)
        call check_duplication(OPTION__STATICEXTERNALPOTENTIALS__POTENTIAL_FORMULA)
        call parse_block_string(blk, row, icol + 1, pot%potential_formula)
        call conv_to_C_string(pot%potential_formula)

        if(pot%type /= EXTERNAL_POT_USDEF) then
          call messages_input_error(namespace, 'ExternalPotentials', 'potential_formula can only be used with user_defined')
        end if

      case(OPTION__STATICEXTERNALPOTENTIALS__DENSITY_FORMULA)
        call check_duplication(OPTION__STATICEXTERNALPOTENTIALS__DENSITY_FORMULA)
        call parse_block_string(blk, row, icol + 1, pot%density_formula)
        call conv_to_C_string(pot%density_formula)
              
        if(pot%type /= EXTERNAL_POT_CHARGE_DENSITY) then
          call messages_input_error(namespace, 'ExternalPotentials', 'density_formula can only be used with charge_density')
        end if

      case default
        call messages_input_error(namespace, 'ExternalPotentials', "Unknown parameter ")
        
      end select

      icol = icol + 2        
    end do
    ! CHECK THAT WHAT WE PARSED MAKES SENSE
    

    if(pot%type == EXTERNAL_POT_USDEF .and. .not. parameter_defined(OPTION__STATICEXTERNALPOTENTIALS__POTENTIAL_FORMULA)) then
      call messages_input_error(namespace, 'ExternalPotentials', "The 'potential_formula' parameter is missing.")
    end if

    if(pot%type == EXTERNAL_POT_CHARGE_DENSITY .and. .not. parameter_defined(OPTION__STATICEXTERNALPOTENTIALS__DENSITY_FORMULA)) then
      call messages_input_error(namespace, 'ExternalPotentials', "The 'density_formula' parameter is missing.")
    end if
    
    if(pot%type == EXTERNAL_POT_FROM_FILE .and. .not. (parameter_defined(OPTION__STATICEXTERNALPOTENTIALS__FILE))) then
      call messages_input_error(namespace, 'ExternalPotentials', "The 'file' parameter is missing.")
    end if

    call iihash_end(read_parameters)

    POP_SUB(read_from_block)

  contains

    logical function parameter_defined(param) result(defined)
      integer(8), intent(in) :: param

      integer :: tmp
      
      PUSH_SUB(read_from_block.parameter_defined)

      tmp = iihash_lookup(read_parameters, int(-param), defined)
      
      POP_SUB(read_from_block.parameter_defined)
    end function parameter_defined

    !------------------------------------------------------
    
    subroutine check_duplication(param)
      integer(8), intent(in) :: param

      PUSH_SUB(read_from_block.check_duplication)

      if(parameter_defined(param)) then
        call messages_input_error(namespace, 'ExternalPotentials', "Duplicated parameter in external potential.")
      end if

      call iihash_insert(read_parameters, int(-param), 1)

      POP_SUB(read_from_block.check_duplication)
    end subroutine check_duplication
    
  end subroutine read_from_block
  ! ---------------------------------------------------------


end module external_potential_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
