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
  use global_oct_m
  use iihash_oct_m
  use interaction_abst_oct_m
  use interaction_partner_oct_m
  use io_function_oct_m
  use linked_list_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use string_oct_m
  implicit none

  private
  public ::               &
    external_potential_t, &
    load_external_potentials

  type, extends(interaction_partner_t) :: external_potential_t
    private

    integer :: type   !< Type of external potential
    character(len=1024) :: potential_formula !< for the user-defined potential
    character(len=200) :: density_formula !< If we have a charge distribution creating the potential
    character(len=MAX_PATH_LEN) :: filename !< for the potential read from a file.
    FLOAT :: omega

    type(mesh_t),    pointer :: mesh  
    type(poisson_t), pointer :: poisson

    FLOAT, allocatable, public :: pot(:)

  contains
    procedure :: calculate => external_potential_calculate
    procedure :: allocate_memory => external_potential_allocate
    procedure :: deallocate_memory => external_potential_deallocate
    procedure :: update_exposed_quantities => external_potential_update_exposed_quantities
    procedure :: update_exposed_quantity => external_potential_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => external_potential_copy_quantities_to_interaction
    final :: external_potential_finalize
  end type external_potential_t

  integer, public, parameter ::  &
    EXTERNAL_POT_USDEF          = 201,           & !< user-defined function for local potential
    EXTERNAL_POT_FROM_FILE      = 202,           &
    EXTERNAL_POT_CHARGE_DENSITY = 203              !< user-defined function for charge density

  interface external_potential_t
    module procedure external_potential_init
  end interface external_potential_t

contains

  function external_potential_init(mesh) result(this)
    class(external_potential_t), pointer :: this
    type(mesh_t), target,     intent(in) :: mesh

    PUSH_SUB(external_potential_init)

    SAFE_ALLOCATE(this)

    this%mesh => mesh

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
  subroutine external_potential_allocate(this)
    class(external_potential_t), intent(inout) :: this

    PUSH_SUB(external_potential_allocate)

    SAFE_ALLOCATE(this%pot(1:this%mesh%np)) 

    POP_SUB(external_potential_allocate)

  end subroutine external_potential_allocate

  ! ---------------------------------------------------------
  subroutine external_potential_deallocate(this)
    class(external_potential_t), intent(inout) :: this

    PUSH_SUB(external_potential_deallocate)

    SAFE_DEALLOCATE_A(this%pot)

    POP_SUB(external_potential_deallocate)

  end subroutine external_potential_deallocate

  ! ---------------------------------------------------------
  logical function external_potential_update_exposed_quantities(this, requested_time, interaction)
    class(external_potential_t), intent(inout) :: this
    type(clock_t),                intent(in)    :: requested_time
    class(interaction_abst_t),    intent(inout) :: interaction

    PUSH_SUB(external_potential_update_exposed_quantities)

    POP_SUB(external_potential_update_exposed_quantities)

  end function external_potential_update_exposed_quantities

  ! ---------------------------------------------------------
  subroutine external_potential_update_exposed_quantity(this, iq, requested_time)
    class(external_potential_t),      intent(inout) :: this
    integer,                           intent(in)    :: iq
    class(clock_t),                    intent(in)    :: requested_time

    PUSH_SUB(external_potential_update_exposed_quantities)

    POP_SUB(external_potential_update_exposed_quantities)

  end subroutine external_potential_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine external_potential_copy_quantities_to_interaction(this, interaction)
    class(external_potential_t),     intent(inout) :: this
    class(interaction_abst_t),        intent(inout) :: interaction

    PUSH_SUB(external_potential_copy_quantities_to_interaction)

    POP_SUB(external_potential_copy_quantities_to_interaction)

  end subroutine external_potential_copy_quantities_to_interaction


  ! ---------------------------------------------------------
  subroutine external_potential_calculate(this, namespace)
    class(external_potential_t), intent(inout) :: this
    type(namespace_t),           intent(in)    :: namespace

    FLOAT :: pot_re, pot_im, r, xx(1:MAX_DIM)
    FLOAT, allocatable :: den(:)
    integer :: ip, err

    PUSH_SUB(external_potential_calculate)

    ASSERT(allocated(this%pot))

    select case(this%type)

    case(EXTERNAL_POT_USDEF)

      do ip = 1, this%mesh%np
        call mesh_r(this%mesh, ip, r, coords = xx)
        call parse_expression(pot_re, pot_im, this%mesh%sb%dim, xx, r, M_ZERO, this%potential_formula)
        this%pot(ip) = pot_re
      end do 

    case(EXTERNAL_POT_FROM_FILE)

      call dio_function_input(trim(this%filename), namespace, this%mesh, this%pot, err)
      if(err /= 0) then
        write(message(1), '(a)')    'Error loading file '//trim(this%filename)//'.'
        write(message(2), '(a,i4)') 'Error code returned = ', err      
        call messages_fatal(2, namespace=namespace)
      end if

    case(EXTERNAL_POT_CHARGE_DENSITY)

      SAFE_ALLOCATE(den(1:this%mesh%np))

      do ip = 1, this%mesh%np
        call mesh_r(this%mesh, ip, r, coords = xx)
        call parse_expression(pot_re, pot_im, this%mesh%sb%dim, xx, r, M_ZERO, this%potential_formula)
        den(ip) = pot_re
      end do

      if(poisson_solver_is_iterative(this%poisson)) then
        ! pot has to be initialized before entering routine
        ! and our best guess for the potential is zero
        this%pot(1:this%mesh%np) = M_ZERO
      end if
      call dpoisson_solve(this%poisson, this%pot, den, all_nodes = .false.)

      SAFE_DEALLOCATE_A(den)

    end select

    POP_SUB(external_potential_calculate)
  end subroutine external_potential_calculate

  subroutine load_external_potentials(external_potentials, namespace, mesh)
    type(linked_list_t),  intent(inout)  :: external_potentials
    type(namespace_t),    intent(in)     :: namespace
    type(mesh_t),         intent(in)     :: mesh

    integer :: n_pot_block, row, read_data
    type(block_t) :: blk
    type(external_potential_t), pointer :: pot

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
        pot => external_potential_t(mesh) 
        !Parse the information from the block
        call read_from_block(pot, namespace, blk, row, read_data)
        ASSERT(read_data > 0)
        !Add this to the list
        call external_potentials%add(pot)  
      end do 
      call parse_block_end(blk)
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

    integer :: ncols, icol, flag, set_read_data, ierr
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
