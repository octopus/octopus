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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: symmetry.F90 3479 2007-11-09 08:36:10Z xavier $

#include "global.h"

module symmetries_m
  use datasets_m
  use global_m
  use geometry_m
  use parser_m
  use messages_m
  use profiling_m
  use species_m
  use symm_op_m

  implicit none

  private
  
  public ::                      &
       symmetries_init,          &
       symmetries_copy,          &
       symmetries_end,           &
       symmetries_number,        &
       symmetries_apply_kpoint,  &
       symmetries_t

  type symmetries_t
    type(symm_op_t), pointer :: ops(:)
    integer                  :: nops
    FLOAT                    :: breakdir(1:3)
  end type symmetries_t

  real(8), parameter :: symprec = CNST(1e-5)

  interface
    ! these functions are defined in spglib_f.c

    integer function spglib_get_max_multiplicity(lattice, position, types, num_atom, symprec)
      real(8), intent(in) :: lattice
      real(8), intent(in) :: position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
    end function spglib_get_max_multiplicity

    integer function spglib_get_symmetry(rotation, translation, max_size, lattice, position, types, num_atom, symprec)
      integer, intent(out) :: rotation
      real(8), intent(out) :: translation
      integer, intent(in)  :: max_size
      real(8), intent(in)  :: lattice
      real(8), intent(in)  :: position
      integer, intent(in)  :: types
      integer, intent(in)  :: num_atom
      real(8), intent(in)  :: symprec
    end function spglib_get_symmetry

    subroutine spglib_show_symmetry(lattice, position, types, num_atom, symprec)
      real(8), intent(in) :: lattice
      real(8), intent(in) :: position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
    end subroutine spglib_show_symmetry

  end interface
  
contains
  
  subroutine symmetries_init(this, geo, dim, rlattice, lsize)
    type(symmetries_t),  intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: rlattice(:, :)
    FLOAT,               intent(in)  :: lsize(:)
    
    integer :: max_size, fullnops
    integer :: idir, iatom, iop
    real(8) :: lattice(1:3, 1:3)
    real(8), allocatable :: position(:, :)
    integer, allocatable :: typs(:)
    type(block_t) :: blk
    integer, allocatable :: rotation(:, :, :)
    FLOAT,   allocatable :: translation(:, :)
    type(symm_op_t) :: tmpop

    lattice(1:3, 1:3) = rlattice(1:3, 1:3)
    SAFE_ALLOCATE(position(1:3, 1:geo%natoms))
    SAFE_ALLOCATE(typs(1:geo%natoms))

    ! we have to fix things for low dimensionality systems
    do idir = dim + 1, 3
      lattice(idir, idir) = M_ONE
    end do
    
    forall(iatom = 1:geo%natoms)
      !this has to be fixed for non-orthogonal cells
      position(1:3, iatom) = M_HALF
      position(1:dim, iatom) = geo%atom(iatom)%x(1:dim)/(M_TWO*lsize(1:dim)) + M_HALF
      typs(iatom) = anint(species_z(geo%atom(iatom)%spec))
    end forall

    ! This outputs information about the symmetries, I will disable it
    ! for the moment as it causes some problems.
    !    call spglib_show_symmetry(lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)

    max_size = spglib_get_max_multiplicity(lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)

    SAFE_ALLOCATE(rotation(1:3, 1:3, 1:max_size))
    SAFE_ALLOCATE(translation(1:3, 1:max_size))

    fullnops = spglib_get_symmetry(rotation(1, 1, 1), translation(1, 1), &
         max_size, lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)

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

    if(parse_block(datasets_check('SymmetryBreakDir'), blk) == 0) then
      
      do idir = 1, dim
        call parse_block_float(blk, 0, idir - 1, this%breakdir(idir))
      end do
      
      call parse_block_end(blk)
      
    end if

    SAFE_ALLOCATE(this%ops(1:fullnops))

    ! check all operation and leave those that kept the symmetry
    ! breaking direction invariant and (for the moment) that do not have a translation
    this%nops = 0
    do iop = 1, fullnops
      call symm_op_init(tmpop, rotation(:, :, iop), translation(:, iop))

      if(symm_op_invariant(tmpop, this%breakdir, symprec) &
        .and. .not. symm_op_has_translation(tmpop, symprec)) then
        this%nops = this%nops + 1
        call symm_op_copy(tmpop, this%ops(this%nops))
      end if
      call symm_op_end(tmpop)
    end do

    write(message(1), '(a,i5,a)') 'Info: The system has ', this%nops, ' symmetries that can be used.'
    call write_info(1)

    SAFE_DEALLOCATE_A(rotation)
    SAFE_DEALLOCATE_A(translation)

  end subroutine symmetries_init

  ! -------------------------------------------------------------------------------

  subroutine symmetries_copy(inp, outp)
    type(symmetries_t),  intent(in)  :: inp
    type(symmetries_t),  intent(out) :: outp

    integer :: iop

    outp%nops = inp%nops
    outp%breakdir = inp%breakdir
    
    SAFE_ALLOCATE(outp%ops(1:outp%nops))

    do iop = 1, outp%nops
      call symm_op_copy(inp%ops(iop), outp%ops(iop))
    end do

  end subroutine symmetries_copy

  ! -------------------------------------------------------------------------------

  subroutine symmetries_end(this)
    type(symmetries_t),  intent(inout) :: this

    integer :: iop

    do iop = 1, this%nops
      call symm_op_end(this%ops(iop))
    end do

    SAFE_DEALLOCATE_P(this%ops)
  end subroutine symmetries_end

  ! -------------------------------------------------------------------------------
  
  integer pure function symmetries_number(this) result(number)
    type(symmetries_t),  intent(in) :: this
    
    number = this%nops
  end function symmetries_number

  ! -------------------------------------------------------------------------------

  subroutine symmetries_apply_kpoint(this, iop, aa, bb)
    type(symmetries_t),  intent(in)  :: this
    integer,             intent(in)  :: iop
    FLOAT,               intent(in)  :: aa(1:3)
    FLOAT,               intent(out) :: bb(1:3)

    ASSERT(0 < iop .and. iop <= this%nops)

    bb(1:3) = symm_op_apply_inv(this%ops(iop), aa(1:3))
  end subroutine symmetries_apply_kpoint

end module symmetries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
