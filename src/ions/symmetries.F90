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
  use loct_parser_m
  use messages_m
  use profiling_m
  use simul_box_m
  use species_m

  implicit none

  private
  
  public ::                   &
       symmetries_init,       &
       symmetries_end,        &
       symmetries_number,     &
       symmetries_apply,      &
       symmetries_t

  type symmetries_t
    private
    integer, pointer :: rotation(:, :, :)
    real(8), pointer :: translation(:, :)
    integer          :: nops
    FLOAT            :: breakdir(1:MAX_DIM)
  end type symmetries_t

  real(8), parameter :: symprec = CNST(1e-5)

  interface
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
  
  subroutine symmetries_init(this, geo, sb)
    type(symmetries_t),  intent(out) :: this
    type(geometry_t),  intent(in)  :: geo
    type(simul_box_t), intent(in)  :: sb
    
    integer :: max_size
    integer :: idir, iatom
    real(8) :: lattice(1:3, 1:3)
    real(8), allocatable :: position(:, :)
    integer, allocatable :: typs(:)
    type(block_t) :: blk

    lattice(1:3, 1:3) = sb%rlattice(1:3, 1:3)
    SAFE_ALLOCATE(position(1:3, 1:geo%natoms))
    SAFE_ALLOCATE(typs(1:geo%natoms))
    
    forall(iatom = 1:geo%natoms)
      !this has to be fixed for non-orthogonal cells
      position(1:3, iatom) = geo%atom(iatom)%x(1:3)/(M_TWO*sb%lsize(1:3)) + M_HALF
      typs(iatom) = anint(species_z(geo%atom(iatom)%spec))
    end forall

    ! This outputs information about the symmetries, I will disable it
    ! for the moment as it causes some problems.
    !    call spglib_show_symmetry(lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)

    max_size = spglib_get_max_multiplicity(lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)

    SAFE_ALLOCATE(this%rotation(1:3, 1:3, 1:max_size))
    SAFE_ALLOCATE(this%translation(1:3, 1:max_size))

    this%nops = spglib_get_symmetry(this%rotation(1, 1, 1), this%translation(1, 1), &
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

    this%breakdir(1:MAX_DIM) = M_ZERO

    if(loct_parse_block(datasets_check('SymmetryBreakDir'), blk) == 0) then
      
      do idir = 1, sb%dim
        call loct_parse_block_float(blk, 0, idir - 1, this%breakdir(idir))
      end do
      
      call loct_parse_block_end(blk)
      
    end if

  end subroutine symmetries_init
  
  subroutine symmetries_end(this)
    type(symmetries_t),  intent(inout) :: this
    
    SAFE_DEALLOCATE_P(this%rotation)
    SAFE_DEALLOCATE_P(this%translation)
  end subroutine symmetries_end
  
  integer pure function symmetries_number(this) result(number)
    type(symmetries_t),  intent(in) :: this
    
    number = this%nops
  end function symmetries_number
  
  subroutine symmetries_apply(this, iop, aa, bb)
    type(symmetries_t),  intent(in)  :: this
    integer,             intent(in)  :: iop
    FLOAT,               intent(in)  :: aa(1:MAX_DIM)
    FLOAT,               intent(out) :: bb(1:MAX_DIM)

    FLOAT :: cc(1:MAX_DIM)

    ASSERT(0 < iop .and. iop <= this%nops)

    bb(1:MAX_DIM) = aa(1:MAX_DIM)
    
    ! if the operation leaves the vector invariant
    cc(1:3) = matmul(this%breakdir(1:3), dble(this%rotation(1:3, 1:3, iop)))
    if(all(abs(cc(1:3) - this%breakdir(1:3)) < symprec)) then
      ! we can use it
      bb(1:3) = matmul(aa(1:3), dble(this%rotation(1:3, 1:3, iop)))
    end if

  end subroutine symmetries_apply

end module symmetries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
