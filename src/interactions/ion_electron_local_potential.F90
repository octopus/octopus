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

module ion_electron_local_potential_oct_m
  use atom_oct_m
  use comm_oct_m
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
  use potential_interaction_oct_m
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
    ion_electron_local_potential_t

  type, extends(potential_interaction_t) :: ion_electron_local_potential_t
    private

    type(poisson_t), pointer :: psolver
    type(mesh_t),    pointer :: mesh
    type(space_t),   pointer :: space

    ! This information should be copied, and not obtained thanks to pointer
    ! This is a temporary change here
    type(distributed_t), pointer, public :: atoms_dist
    type(atom_t), pointer, public :: atom(:)
    FLOAT, pointer, public :: pos(:,:)
    type(lattice_vectors_t), pointer :: latt

    ! Temporary pointer to namespace
    type(namespace_t), pointer :: namespace

    logical :: have_density

  contains
    procedure :: init => ion_electron_local_potential_init
    procedure :: calculate => ion_electron_local_potential_calculate
    procedure :: end => ion_electron_local_potential_end
    final :: ion_electron_local_potential_finalize
  end type ion_electron_local_potential_t

  interface ion_electron_local_potential_t
    module procedure ion_electron_local_potential_constructor
  end interface ion_electron_local_potential_t

contains

  ! ---------------------------------------------------------
  function ion_electron_local_potential_constructor(partner) result(this)
    class(interaction_partner_t), target, intent(inout) :: partner
    class(ion_electron_local_potential_t),               pointer       :: this

    PUSH_SUB(ion_electron_local_potential_constructor)

    SAFE_ALLOCATE(this)

    this%label = "ion-electron local"

    this%partner => partner

    nullify(this%psolver)
    nullify(this%mesh)

    !We do not need any system quantity here
    this%n_system_quantities = 0

    POP_SUB(ion_electron_local_potential_constructor)
  end function ion_electron_local_potential_constructor

  ! ---------------------------------------------------------
  subroutine ion_electron_local_potential_init(this, mesh, psolver, ions, namespace)
    class(ion_electron_local_potential_t),    intent(inout) :: this
    type(mesh_t),                     target, intent(in)    :: mesh
    type(poisson_t),                  target, intent(in)    :: psolver
    type(ions_t),                     target, intent(in)    :: ions
    type(namespace_t),                target, intent(in)    :: namespace

    integer :: ia

    PUSH_SUB(ion_electron_local_potential_init)

    this%mesh => mesh
    this%psolver => psolver

    SAFE_ALLOCATE(this%potential(1:mesh%np, 1:1))

    this%have_density = .false.
    do ia = 1, ions%nspecies
      if(local_potential_has_density(ions%space, ions%species(ia))) then
        this%have_density = .true.
        exit
      end if
    end do

    this%atoms_dist => ions%atoms_dist
    this%atom => ions%atom
    this%space => ions%space
    this%pos => ions%pos
    this%latt => ions%latt

    this%namespace => namespace

    POP_SUB(ion_electron_local_potential_init)
  end subroutine ion_electron_local_potential_init

  ! ---------------------------------------------------------
  ! Note: this code is adapted from epot_local_potential. There are code duplication at the moment
  subroutine ion_electron_local_potential_calculate(this)
    class(ion_electron_local_potential_t),             intent(inout) :: this

    FLOAT, allocatable :: density(:), rho(:), vl(:)
    type(submesh_t) :: sphere
    type(ps_t), pointer :: ps
    integer :: ia, ip
    FLOAT :: radius
    type(profile_t), save :: prof

    PUSH_SUB(ion_electron_local_potential_calculate)

    call profiling_in(prof, "ION_ELEC_LOC_INTER")

    this%potential(:,1) = M_ZERO

    if(this%have_density) then
      SAFE_ALLOCATE(density(1:this%mesh%np))
      density = M_ZERO
    end if

    do ia = this%atoms_dist%start, this%atoms_dist%end
      ! Local potential, we can get it by solving the Poisson equation
      ! (for all-electron species or pseudopotentials in periodic
      ! systems) or by applying it directly to the grid

      if(local_potential_has_density(this%space, this%atom(ia)%species)) then
        
        SAFE_ALLOCATE(rho(1:this%mesh%np))
        call species_get_long_range_density(this%atom(ia)%species, this%namespace, this%space, this%latt, &
          this%pos(:, ia), this%mesh, rho)
        call lalg_axpy(this%mesh%np, M_ONE, rho, density)
        SAFE_DEALLOCATE_A(rho)

      else

        SAFE_ALLOCATE(vl(1:this%mesh%np))
        call species_get_local(this%atom(ia)%species, this%namespace, this%space, this%latt, &
          this%pos(:, ia), this%mesh, vl)
        call lalg_axpy(this%mesh%np, M_ONE, vl, this%potential(:,1))
        SAFE_DEALLOCATE_A(vl)

      end if

      !then we add the localized part
      if(species_is_ps(this%atom(ia)%species)) then

        ps => species_ps(this%atom(ia)%species)

        radius = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold) + this%mesh%spacing(1)

        call submesh_init(sphere, this%space, this%mesh, this%latt, this%pos(:, ia), radius)
        SAFE_ALLOCATE(vl(1:sphere%np))

        do ip = 1, sphere%np
          vl(ip) = spline_eval(ps%vl, sphere%r(ip))
        end do

        call submesh_add_to_mesh(sphere, vl, this%potential(:,1))

        SAFE_DEALLOCATE_A(vl)
        call submesh_end(sphere)
        nullify(ps)
      end if
    end do

    ! reduce over atoms if required
    if(this%atoms_dist%parallel) then
      call comm_allreduce(this%atoms_dist%mpi_grp, this%potential(:,1), dim = this%mesh%np)
      if(this%have_density) then
        call comm_allreduce(this%atoms_dist%mpi_grp, density, dim = this%mesh%np)
      end if
    end if

    ! now we solve the poisson equation with the density of all nodes
    if(this%have_density) then

      SAFE_ALLOCATE(vl(1:this%mesh%np_part))
      if(poisson_solver_is_iterative(this%psolver)) then
        vl(1:this%mesh%np) = M_ZERO
      end if
      call dpoisson_solve(this%psolver, vl, density)
      call lalg_axpy(this%mesh%np, M_ONE, vl, this%potential(:,1))
      SAFE_DEALLOCATE_A(vl)

    end if
    SAFE_DEALLOCATE_A(density)
 
    call profiling_out(prof)

    POP_SUB(ion_electron_local_potential_calculate)
  end subroutine ion_electron_local_potential_calculate

  ! ---------------------------------------------------------
  subroutine ion_electron_local_potential_end(this)
    class(ion_electron_local_potential_t), intent(inout) :: this

    PUSH_SUB(ion_electron_local_potential_end)

    SAFE_DEALLOCATE_A(this%potential)
    nullify(this%mesh)
    nullify(this%psolver)

    call interaction_with_partner_end(this)

    POP_SUB(ion_electron_local_potential_end)
  end subroutine ion_electron_local_potential_end


  ! ---------------------------------------------------------
  subroutine ion_electron_local_potential_finalize(this)
    type(ion_electron_local_potential_t), intent(inout) :: this

    PUSH_SUB(ion_electron_local_potential_finalize)

    call ion_electron_local_potential_end(this)

    POP_SUB(ion_electron_local_potential_finalize)
  end subroutine ion_electron_local_potential_finalize

end module ion_electron_local_potential_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
