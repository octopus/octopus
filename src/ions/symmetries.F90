!! Copyright (C) 2009 X. Andrade
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

module symmetries_oct_m
  use iso_c_binding
  use geometry_oct_m
  use global_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use species_oct_m
  use spglib_f08
  use symm_op_oct_m

  implicit none

  private
  
  public ::                        &
    symmetries_t,                  &
    symmetries_init,               &
    symmetries_copy,               &
    symmetries_end,                &
    symmetries_number,             &
    symmetries_apply_kpoint_red,   &
    symmetries_space_group_number, &
    symmetries_have_break_dir,     &
    symmetries_identity_index,     &
    symmetries_write_info

  type symmetries_t
    type(symm_op_t), allocatable :: ops(:)
    integer                  :: nops
    FLOAT                    :: breakdir(1:3)
    integer                  :: space_group
    logical                  :: any_non_spherical
    logical                  :: symmetries_compute
    character(len=6)         :: group_name
    character(len=30)        :: group_elements
    character(len=10)        :: symbol
    character(len=6)         :: schoenflies
  end type symmetries_t

  real(8), parameter, public :: SYMPREC = CNST(1e-5)

  !> NOTE: unfortunately, these routines use global variables shared among them
  interface
    subroutine symmetries_finite_init(natoms, types, positions, verbosity, point_group)
      integer, intent(in)  :: natoms
      integer, intent(in)  :: types !< (natoms)
      real(8), intent(in)  :: positions !< (3, natoms)
      integer, intent(in)  :: verbosity
      integer, intent(out) :: point_group
    end subroutine symmetries_finite_init

    subroutine symmetries_finite_get_group_name(point_group, name)
      integer,          intent(in)  :: point_group
      character(len=*), intent(out) :: name
    end subroutine symmetries_finite_get_group_name

    subroutine symmetries_finite_get_group_elements(point_group, elements)
      integer,          intent(in)  :: point_group
      character(len=*), intent(out) :: elements
    end subroutine symmetries_finite_get_group_elements

    subroutine symmetries_finite_end()
    end subroutine symmetries_finite_end
  end interface

contains

  subroutine symmetries_init(this, geo, dim, periodic_dim, rlattice, klattice)
    type(symmetries_t),  intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    integer,             intent(in)  :: dim
    integer,             intent(in)  :: periodic_dim
    FLOAT,               intent(in)  :: rlattice(:, :)
    FLOAT,               intent(in)  :: klattice(:, :)

    integer :: max_size, dim4syms
    integer :: idir, iatom, iop, verbosity, point_group
    real(8) :: lattice(1:3, 1:3)
    real(8), allocatable :: position(:, :)
    integer, allocatable :: typs(:)
    type(block_t) :: blk
    type(symm_op_t) :: tmpop
    integer :: natoms, identity(3,3)
    logical :: found_identity, is_supercell
    integer                  :: fullnops
    integer, allocatable     :: rotation(:, :, :)
    real(8), allocatable     :: translation(:, :)
    character(kind=c_char) :: c_symbol(11), c_schoenflies(7) 
    logical :: def_sym_comp
    
    PUSH_SUB(symmetries_init)

    ! if someone cares, they could try to analyze the symmetry point group of the individual species too
    this%any_non_spherical = .false.
    do iatom = 1, geo%natoms
      this%any_non_spherical = this%any_non_spherical                   .or. &
        species_type(geo%atom(iatom)%species) == SPECIES_USDEF          .or. &
        species_type(geo%atom(iatom)%species) == SPECIES_JELLIUM_SLAB   .or. &
        species_type(geo%atom(iatom)%species) == SPECIES_CHARGE_DENSITY .or. &
        species_type(geo%atom(iatom)%species) == SPECIES_FROM_FILE      .or. &
        species_type(geo%atom(iatom)%species) == SPECIES_FROZEN
      if(this%any_non_spherical)exit
    end do

    dim4syms = min(3,dim)

    def_sym_comp = (geo%natoms < 100)
    def_sym_comp = def_sym_comp .and. dim == 3
    
    !%Variable SymmetriesCompute
    !%Type logical
    !%Section Execution::Symmetries
    !%Description
    !% If disabled, <tt>Octopus</tt> will not compute
    !% nor print the symmetries.
    !%
    !% By default, symmetries are computed when running in 3
    !% dimensions for systems with less than 100 atoms.
    !%End
    call parse_variable('SymmetriesCompute', def_sym_comp, this%symmetries_compute)

    if(this%symmetries_compute .and. dim /= 3) then
      call messages_experimental('symmetries for non 3D systems')
    end if
    
    if(this%any_non_spherical .or. .not. this%symmetries_compute) then
      call init_identity()

      POP_SUB(symmetries_init)
      return
    end if

    ! In all cases, we must check that the grid respects the symmetries. --DAS

    if (periodic_dim == 0) then

      call init_identity()

      ! for the moment symmetries are only used for information, so we compute them only on one node.
      if(mpi_grp_is_root(mpi_world)) then
        natoms = max(1,geo%natoms)
        verbosity = -1

        SAFE_ALLOCATE(position(1:3, natoms))
        SAFE_ALLOCATE(typs(1:natoms))

        forall(iatom = 1:geo%natoms)
          position(1:3, iatom) = M_ZERO
          position(1:dim4syms, iatom) = geo%atom(iatom)%x(1:dim4syms)
          typs(iatom) = species_index(geo%atom(iatom)%species)
        end forall

        if (this%symmetries_compute) then
          call symmetries_finite_init(geo%natoms, typs(1), position(1, 1), verbosity, point_group)
          call symmetries_finite_get_group_name(point_group, this%group_name)
          call symmetries_finite_get_group_elements(point_group, this%group_elements)
          call symmetries_finite_end()
        end if
        SAFE_DEALLOCATE_A(position)
        SAFE_DEALLOCATE_A(typs)
      end if

    else

      SAFE_ALLOCATE(position(1:3,1:geo%natoms))  ! transpose!!
      SAFE_ALLOCATE(typs(1:geo%natoms))

      do iatom = 1, geo%natoms
        position(1:3,iatom) = M_ZERO

        if(.not. geo%reduced_coordinates) then
          ! Transform atomic positions to reduced coordinates
          position(1:dim4syms,iatom) = matmul(geo%atom(iatom)%x(1:dim4syms),klattice(1:dim4syms,1:dim4syms))/(M_TWO*M_PI) 
        else
          position(1:dim4syms,iatom) = geo%atom(iatom)%x(1:dim4syms)
        end if
        position(1:dim4syms,iatom) = position(1:dim4syms,iatom)- M_HALF
        do idir = 1, dim4syms
          position(idir,iatom) = position(idir,iatom) - anint(position(idir,iatom))
        end do
        position(1:dim4syms,iatom) = position(1:dim4syms,iatom) + M_HALF

        typs(iatom) = species_index(geo%atom(iatom)%species)
      end do

      lattice = M_ZERO
      !NOTE: Why "inverse matrix" ? (NTD)
      ! get inverse matrix to extract reduced coordinates for spglib
      lattice(1:dim, 1:dim) = rlattice(1:dim, 1:dim)
      ! transpose the lattice vectors for use in spglib as row-major matrix
      lattice(:,:) = transpose(lattice(:,:))
      ! fix things for low-dimensional systems: higher dimension lattice constants set to 1
      do idir = dim + 1, 3
        lattice(idir, idir) = M_ONE
      end do

      this%space_group = spg_get_international(c_symbol, lattice(1,1), &
                 position(1,1), typs(1), geo%natoms, SYMPREC)
      this%space_group = spg_get_schoenflies(c_schoenflies, lattice(1, 1), &
                 position(1, 1), typs(1), geo%natoms, SYMPREC)
      
      if(this%space_group == 0) then
        message(1) = "Symmetry analysis failed in spglib. Disabling symmetries."
        call messages_warning(1)

        do iatom = 1, geo%natoms
          write(message(1),'(a,i6,a,3f12.6,a,3f12.6)') 'type ', typs(iatom), &
            ' reduced coords ', position(:, iatom), ' cartesian coords ', geo%atom(iatom)%x(:)
          call messages_info(1)
        end do

        call init_identity()
        POP_SUB(symmetries_init)
        return
      end if

      call c_to_f_string(c_symbol, this%symbol)
      call c_to_f_string(c_schoenflies, this%schoenflies)
      
      max_size = spg_get_multiplicity(lattice(1, 1), position(1, 1), &
                                      typs(1), geo%natoms, SYMPREC)

      ! spglib returns row-major not column-major matrices!!! --DAS
      SAFE_ALLOCATE(rotation(1:3, 1:3, 1:max_size))
      SAFE_ALLOCATE(translation(1:3, 1:max_size))

      fullnops = spg_get_symmetry(rotation(1, 1, 1), translation(1, 1), &
        max_size, lattice(1, 1), position(1, 1), typs(1), geo%natoms, SYMPREC)

      do iop = 1, fullnops
        ! transpose due to array order difference between C and fortran
        rotation(:,:,iop) = transpose(rotation(:,:,iop))
        ! sometimes spglib may return lattice vectors as 'fractional' translations        
        translation(:, iop) = translation(:, iop) - anint(translation(:, iop) + M_HALF*SYMPREC)
      end do

      ! we need to check that it is not a supercell, as in the QE routine (sgam_at)
      ! they disable fractional translations if the identity has one, because the sym ops might not form a group.
      ! spglib may return duplicate operations in this case!
      is_supercell = (fullnops > 48)
      found_identity = .false.
      identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
      do iop = 1, fullnops
        if(all(rotation(1:3, 1:3, iop) == identity(1:3, 1:3))) then
          found_identity = .true.
          if(any(abs(translation(1:3, iop)) > real(SYMPREC, REAL_PRECISION))) then
            is_supercell = .true.
            write(message(1),'(a,3f12.6)') 'Identity has a fractional translation ', translation(1:3, iop)
            call messages_info(1)
          end if
        end if
      end do
      if(.not. found_identity) then
        message(1) = "Symmetries internal error: Identity is missing from symmetry operations."
        call messages_fatal(1)
      end if
    
      if(is_supercell) then
        message(1) = "Disabling fractional translations. System appears to be a supercell."
        call messages_info(1)
      end if
      ! actually, we do not use fractional translations regardless currently

      ! this is a hack to get things working, this variable should be
      ! eliminated and the direction calculated automatically from the
      ! perturbations.

      !%Variable SymmetryBreakDir
      !%Type block
      !%Section Mesh::Simulation Box
      !%Description
      !% This variable specifies a direction in which the symmetry of
      !% the system will be broken. This is useful for generating <i>k</i>-point
      !% grids when an external perturbation is applied.
      !%End

      this%breakdir(1:3) = M_ZERO

      if(parse_block('SymmetryBreakDir', blk) == 0) then

        do idir = 1, dim4syms
          call parse_block_float(blk, 0, idir - 1, this%breakdir(idir))
        end do

        call parse_block_end(blk)

      end if

      SAFE_ALLOCATE(this%ops(1:fullnops))


      ! check all operations and leave those that kept the symmetry-breaking
      ! direction invariant and (for the moment) that do not have a translation
      this%nops = 0
      do iop = 1, fullnops
        call symm_op_init(tmpop, rotation(1:3, 1:3, iop), rlattice(1:dim4syms,1:dim4syms), &
                              klattice(1:dim4syms,1:dim4syms), dim4syms, &
                              real(translation(1:3, iop), REAL_PRECISION))

        if(symm_op_invariant_cart(tmpop, this%breakdir, real(SYMPREC, REAL_PRECISION)) &
         .and. .not. symm_op_has_translation(tmpop, real(SYMPREC, REAL_PRECISION))) then
          this%nops = this%nops + 1
          call symm_op_copy(tmpop, this%ops(this%nops))
        end if
        call symm_op_end(tmpop)
      end do

      SAFE_DEALLOCATE_A(position)
      SAFE_DEALLOCATE_A(typs)
      SAFE_DEALLOCATE_A(rotation)
      SAFE_DEALLOCATE_A(translation)

    end if

    call symmetries_write_info(this, dim, periodic_dim, stdout)

    POP_SUB(symmetries_init)
    
  contains

    subroutine init_identity()
      
      PUSH_SUB(symmetries_init.init_identity)
      
      SAFE_ALLOCATE(this%ops(1:1))
      this%nops = 1
      call symm_op_init(this%ops(1), reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/)), & 
                  rlattice, klattice, dim4syms)
      this%breakdir = M_ZERO
      this%space_group = 1
      
      POP_SUB(symmetries_init.init_identity)
      
    end subroutine init_identity

    subroutine c_to_f_string(c_string, f_string)
      character(kind=c_char,len=1), intent(in)  :: c_string(*)
      character(len=*),             intent(out) :: f_string

      integer :: i

      i = 1
      do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = c_string(i)
        i = i + 1
      end do

      if (i <= len(f_string)) f_string(i:) = ' '

    end subroutine c_to_f_string
    
  end subroutine symmetries_init

  

  ! -------------------------------------------------------------------------------
  subroutine symmetries_copy(inp, outp)
    type(symmetries_t),  intent(in)  :: inp
    type(symmetries_t),  intent(out) :: outp

    integer :: iop

    PUSH_SUB(symmetries_copy)

    outp%nops = inp%nops
    outp%breakdir = inp%breakdir
    
    SAFE_ALLOCATE(outp%ops(1:outp%nops))

    do iop = 1, outp%nops
      call symm_op_copy(inp%ops(iop), outp%ops(iop))
    end do

    outp%group_name = inp%group_name
    outp%group_elements = inp%group_elements
    outp%symbol = inp%symbol
    outp%schoenflies = inp%schoenflies

    POP_SUB(symmetries_copy)
  end subroutine symmetries_copy

  ! -------------------------------------------------------------------------------

  subroutine symmetries_end(this)
    type(symmetries_t),  intent(inout) :: this

    integer :: iop

    PUSH_SUB(symmetries_end)

    if(allocated(this%ops)) then
      do iop = 1, this%nops
        call symm_op_end(this%ops(iop))
      end do

      SAFE_DEALLOCATE_A(this%ops)
    end if

    POP_SUB(symmetries_end)
  end subroutine symmetries_end

  ! -------------------------------------------------------------------------------
  
  integer pure function symmetries_number(this) result(number)
    type(symmetries_t),  intent(in) :: this
    
    number = this%nops
  end function symmetries_number

  ! -------------------------------------------------------------------------------

  subroutine symmetries_apply_kpoint_red(this, iop, aa, bb)
    type(symmetries_t),  intent(in)  :: this
    integer,             intent(in)  :: iop
    FLOAT,               intent(in)  :: aa(1:3)
    FLOAT,               intent(out) :: bb(1:3)

    PUSH_SUB(symmetries_apply_kpoint_red)

    ASSERT(0 < iop .and. iop <= this%nops)

    bb(1:3) = symm_op_apply_transpose_red(this%ops(iop), aa(1:3))

    POP_SUB(symmetries_apply_kpoint_red)
  end subroutine symmetries_apply_kpoint_red

  ! -------------------------------------------------------------------------------

  integer pure function symmetries_space_group_number(this) result(number)
    type(symmetries_t),  intent(in) :: this
    
    number = this%space_group
  end function symmetries_space_group_number

  ! -------------------------------------------------------------------------------

  logical pure function symmetries_have_break_dir(this) result(have)
    type(symmetries_t),  intent(in)  :: this

    have = any(abs(this%breakdir(1:3)) > M_EPSILON)
  end function symmetries_have_break_dir

  ! -------------------------------------------------------------------------------

  integer pure function symmetries_identity_index(this) result(index)
    type(symmetries_t),  intent(in)  :: this

    integer :: iop
    
    do iop = 1, this%nops
      if(symm_op_is_identity(this%ops(iop))) then
        index = iop
        cycle
      end if
    end do

  end function symmetries_identity_index

  ! ---------------------------------------------------------
  subroutine symmetries_write_info(this, dim, periodic_dim, iunit)
    type(symmetries_t),    intent(in) :: this
    integer,               intent(in) :: dim, periodic_dim
    integer,               intent(in) :: iunit
    
    integer :: iop
 
    PUSH_SUB(symmetries_write_info)
    
    call messages_print_stress(iunit, 'Symmetries')

    if(this%any_non_spherical) then
      message(1) = "Symmetries are disabled since non-spherically symmetric species may be present."
      call messages_info(1,iunit = iunit)
      call messages_print_stress(iunit)
    end if

    if(.not. this%symmetries_compute) then
      message(1) = "Symmetries have been disabled by SymmetriesCompute = false."
      call messages_info(1,iunit = iunit)
      call messages_print_stress(iunit)
      POP_SUB(symmetries_write_info)
      return
    end if

    if (periodic_dim == 0) then
      ! At the moment only the root node has information about symetries of finite systems.
      if(mpi_grp_is_root(mpi_world)) then
        if (this%symmetries_compute) then
          call messages_write('Symmetry elements : '//trim(this%group_elements), new_line = .true.)
          call messages_write('Symmetry group    : '//trim(this%group_name))
          call messages_info(iunit = iunit)
        end if
      end if
    else
      write(message(1),'(a, i4)') 'Space group No. ', this%space_group
      write(message(2),'(2a)') 'International: ', trim(this%symbol)
      write(message(3),'(2a)') 'Schoenflies: ', trim(this%schoenflies)
      call messages_info(3,iunit = iunit)

      write(message(1),'(a7,a31,12x,a33)') 'Index', 'Rotation matrix', 'Fractional translations'
      call messages_info(1,iunit = iunit)
      do iop = 1, this%nops
        ! list all operations and leave those that kept the symmetry-breaking
        ! direction invariant and (for the moment) that do not have a translation
        if(dim == 1) &
        write(message(1),'(i5,1x,a,2x,1(1i4,2x),1f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
                                                                    symm_op_translation_vector_red(this%ops(iop))
        if(dim == 2) &
        write(message(1),'(i5,1x,a,2x,2(2i4,2x),2f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
                                                                    symm_op_translation_vector_red(this%ops(iop))
        if(dim == 3) &
        write(message(1),'(i5,1x,a,2x,3(3i4,2x),3f12.6)') iop, ':', symm_op_rotation_matrix_red(this%ops(iop)), &
                                                                    symm_op_translation_vector_red(this%ops(iop))
        call messages_info(1,iunit = iunit)
      end do
      write(message(1), '(a,i5,a)') 'Info: The system has ', this%nops, ' symmetries that can be used.'
      call messages_info(iunit = iunit)
    end if
    call messages_print_stress(iunit)

    POP_SUB(symmetries_write_info)
  end subroutine symmetries_write_info


end module symmetries_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
