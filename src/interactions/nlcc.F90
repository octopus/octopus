!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2021 Nicolas Tancogne-Dejean
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

module nlcc_oct_m
  use atom_oct_m
  use comm_oct_m
  use density_interaction_oct_m
  use distributed_oct_m
  use epot_oct_m
  use global_oct_m
  use interaction_with_partner_oct_m
  use interaction_partner_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use ps_oct_m
  use quantity_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use splines_oct_m
  use submesh_oct_m

  implicit none

  private
  public ::                &
    nlcc_t

  type, extends(density_interaction_t) :: nlcc_t
    private

    type(mesh_t),    pointer :: mesh
    type(space_t),   pointer :: space

    ! This information should be copied, and not obtained thanks to pointer
    ! This is a temporary change here
    type(distributed_t), pointer, public :: atoms_dist
    type(atom_t), pointer, public :: atom(:)
    FLOAT, pointer, public :: pos(:,:)
    type(lattice_vectors_t), pointer :: latt

  contains
    procedure :: init => nlcc_init
    procedure :: calculate => nlcc_calculate
    procedure :: end => nlcc_end
    final :: nlcc_finalize
  end type nlcc_t

  interface nlcc_t
    module procedure nlcc_constructor
  end interface nlcc_t

contains

  ! ---------------------------------------------------------
  function nlcc_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(nlcc_t),               pointer       :: this

    PUSH_SUB(nlcc_constructor)

    SAFE_ALLOCATE(this)

    this%label = "nlcc"

    this%partner => partner

    nullify(this%mesh)

    !We do not need any system quantity here
    this%n_system_quantities = 0

    POP_SUB(nlcc_constructor)
  end function nlcc_constructor

  ! ---------------------------------------------------------
  subroutine nlcc_init(this, mesh, ions)
    class(nlcc_t),         intent(inout) :: this
    type(mesh_t),     target, intent(in) :: mesh
    type(ions_t),     target, intent(in) :: ions


    PUSH_SUB(nlcc_init)

    this%mesh => mesh

    SAFE_ALLOCATE(this%density(1:mesh%np,1:1))

    this%atoms_dist => ions%atoms_dist
    this%atom => ions%atom
    this%space => ions%space
    this%pos => ions%pos
    this%latt => ions%latt

    POP_SUB(nlcc_init)
  end subroutine nlcc_init

  ! ---------------------------------------------------------
  subroutine nlcc_calculate(this)
    class(nlcc_t),             intent(inout) :: this

    integer :: ia
    type(profile_t), save :: prof

    PUSH_SUB(nlcc_calculate)

    call profiling_in(prof, "NLCC_INTER")

    this%density(:,:) = M_ZERO

    do ia = this%atoms_dist%start, this%atoms_dist%end
      if(species_has_nlcc(this%atom(ia)%species) .and. species_is_ps(this%atom(ia)%species)) then
        call species_get_nlcc(this%atom(ia)%species, this%space, this%latt, this%pos(:, ia), this%mesh, &
              this%density(:,1), accumulate=.true.)
      endif
    end do

    ! reduce over atoms if required
    if(this%atoms_dist%parallel) then
      call comm_allreduce(this%atoms_dist%mpi_grp, this%density(:,1))
    end if
 
    call profiling_out(prof)

    POP_SUB(nlcc_calculate)
  end subroutine nlcc_calculate

  ! ---------------------------------------------------------
  subroutine nlcc_end(this)
    class(nlcc_t), intent(inout) :: this

    PUSH_SUB(nlcc_end)

    SAFE_DEALLOCATE_A(this%density)
    nullify(this%mesh)

    call interaction_with_partner_end(this)

    POP_SUB(nlcc_end)
  end subroutine nlcc_end


  ! ---------------------------------------------------------
  subroutine nlcc_finalize(this)
    type(nlcc_t), intent(inout) :: this

    PUSH_SUB(nlcc_finalize)

    call nlcc_end(this)

    POP_SUB(nlcc_finalize)
  end subroutine nlcc_finalize

end module nlcc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
