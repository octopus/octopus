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
  use global_m
  use geometry_m
  use messages_m
  use profiling_m
  use simul_box_m
  use species_m

  implicit none

  private
  
  public ::                   &
       symmetries_init,       &
       symmetries_end,        &
       symmetries_t

  type symmetries_t
    private
    real(8), pointer :: rotation(:, :, :)
    real(8), pointer :: translation(:, :)
    integer          :: nops
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
      real(8), intent(out) :: rotation
      real(8), intent(out) :: translation
      integer, intent(in)  :: max_size
      real(8), intent(in)  :: lattice
      real(8), intent(in)  :: position
      integer, intent(in)  :: types
      integer, intent(in)  :: num_atom
      real(8), intent(in)  :: symprec
    end function spglib_get_symmetry
  end interface
  
contains
  
  subroutine symmetries_init(this, geo, sb)
    type(symmetries_t),  intent(out) :: this
    type(geometry_t),  intent(in)  :: geo
    type(simul_box_t), intent(in)  :: sb
    
    integer :: max_size
    integer :: idir, jdir, iatom
    real(8) :: lattice(1:3, 1:3)
    real(8), allocatable :: position(:, :)
    integer, allocatable :: typs(:)
    
    forall(idir = 1:sb%dim, jdir = 1:sb%dim) lattice(jdir, idir) = sb%rlattice(idir, jdir)*sb%lsize(idir)
    ALLOCATE(position(1:3, geo%natoms), 3*geo%natoms)
    ALLOCATE(typs(geo%natoms), geo%natoms)
    
    forall(iatom = 1:geo%natoms)
      position(1:3, iatom) = geo%atom(iatom)%x(1:3)
      typs(iatom) = anint(geo%atom(iatom)%spec%z)
    end forall
      
    max_size = spglib_get_max_multiplicity(lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)
    
    ALLOCATE(this%rotation(1:3, 1:3, max_size), 9*max_size)
    ALLOCATE(this%translation(1:3, max_size), 3*max_size)
    
    this%nops = spglib_get_symmetry(this%rotation(1, 1, 1), this%translation(1, 1), &
         max_size, lattice(1, 1), position(1, 1), typs(1), geo%natoms, symprec)
  end subroutine symmetries_init
  
  subroutine symmetries_end(this)
    type(symmetries_t),  intent(inout) :: this
    
    deallocate(this%rotation)
    deallocate(this%translation)
  end subroutine symmetries_end
  
  integer pure function symmetries_number(this) result(number)
    type(symmetries_t),  intent(in) :: this
    
    number = this%nops
  end function symmetries_number
  
end module symmetries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
