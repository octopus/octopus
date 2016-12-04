!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$

#include "global.h"

module lda_u_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use species_oct_m
  use species_pot_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use types_oct_m  
 
  implicit none

  private

  public ::                  &
       lda_u_t,              &
       lda_u_init,           &
       dhubbard_apply,       &
       zhubbard_apply,       &
       dcorrect_energy_dc,   &
       zcorrect_energy_dc,   &
       dupdate_occ_matrices, &
       zupdate_occ_matrices, &
       lda_u_end

  type orbital_t
    type(submesh_t)     :: sphere      !The submesh of the orbital
    FLOAT, pointer      :: dorbital_sphere(:) !The orbital, if real, on the submesh
    CMPLX, pointer      :: zorbital_sphere(:) !The orbital, if complex, on the submesh
    FLOAT, pointer      :: dorbital_mesh(:) !The orbital, if real, on the full mesh
    CMPLX, pointer      :: zorbital_mesh(:) !The orbital, if complex, on the full mesh
  end type orbital_t


  type lda_u_t
    logical                  :: apply
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:)
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    type(orbital_t), pointer :: orbitals(:,:,:)    

    integer             :: natoms        !> Number of atoms (copied from geometry_t)
    integer             :: nspins
    integer, pointer    :: norbs(:)      !> Number of orbitals
    integer             :: maxnorbs      !> Maximal number of orbitals for all the atoms
    integer             :: projection    !> The method used to perform the projection
    logical             :: truncate
  end type lda_u_t

contains

 subroutine lda_u_init(this, gr, geo, st)
  implicit none

  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st

  integer :: maxorbs

  PUSH_SUB(lda_u_init)

  ASSERT(.not. this%apply)

  call messages_print_stress(stdout, "LDA+U")
 
  if(st%parallel_in_states) call messages_not_implemented("lda+u parallel in states")
  if(gr%mesh%parallel_in_domains) call messages_not_implemented("lda+u parallel in domains")
  if(st%d%ispin == SPINORS) call messages_not_implemented("lda+u with spinors") 

  this%apply = .true.

  !%Variable OrbitalsTruncateToNLRadius
  !%Type logical
  !%Default no
  !%Section Hamiltonian::LDA+U
  !%Description
  !% If set to yes, Octopus will truncate the orbitals 
  !% to the radius of the nonlocal part of the pseudopotential.
  !% This makes the orbitals basis to b non-overlapping between different atoms
  !%End
  call parse_variable('OrbitalsTruncateToNLRadius', .false., this%truncate)

  !%Variable OrbitalsProjectionMethod
  !%Type integer
  !%Default projection_fullmesh
  !%Section Hamiltonian::LDA+U
  !%Description
  !% This variable controls how the projection on the orbitals is done.
  !%Option fullmesh 0
  !% The projection is done on the full mesh. This is the default  value.
  !%Option sphere 1
  !% The projection is done using the submesh
  !%End
  call parse_variable('OrbitalsProjectionMethod', OPTION__ORBITALSPROJECTIONMETHOD__FULLMESH, this%projection)
  if(this%projection==OPTION__ORBITALSPROJECTIONMETHOD__SPHERE) &
    call messages_not_implemented("OrbitalProjectionMethod=sphere") 


  nullify(this%dn)
  nullify(this%zn)
  nullify(this%dV)
  nullify(this%zV)

  this%natoms = geo%natoms
  this%nspins = st%d%nspin

  if (states_are_real(st)) then
    call dconstruct_orbital_basis(this, geo, gr%mesh, st)
    maxorbs = maxval(this%norbs)
    SAFE_ALLOCATE(this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = M_ZERO
    SAFE_ALLOCATE(this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = M_ZERO
  else
    call zconstruct_orbital_basis(this, geo, gr%mesh, st)
    maxorbs = maxval(this%norbs)
    SAFE_ALLOCATE(this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = cmplx(M_ZERO,M_ZERO)
    SAFE_ALLOCATE(this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = cmplx(M_ZERO,M_ZERO)
  end if

  call messages_print_stress(stdout)

  POP_SUB(lda_u_init)
 end subroutine lda_u_init


 subroutine lda_u_end(this)
   implicit none
   type(lda_u_t), intent(inout) :: this

   integer :: iat, ispin, iorb
  
   PUSH_SUB(lda_u_end)  

   this%apply = .false.

   SAFE_DEALLOCATE_P(this%dn)
   SAFE_DEALLOCATE_P(this%zn)
   SAFE_DEALLOCATE_P(this%dV)
   SAFE_DEALLOCATE_P(this%zV) 
 
   do iat = 1, this%natoms
     do ispin = 1, this%nspins 
       do iorb = 1, this%norbs(iat)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%dorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%zorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%dorbital_mesh)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%zorbital_mesh)
         call submesh_end(this%orbitals(iorb,ispin,iat)%sphere)
       end do
     end do
   end do
  
   SAFE_DEALLOCATE_P(this%norbs)
   SAFE_DEALLOCATE_P(this%orbitals)

   POP_SUB(lda_u_end)
 end subroutine lda_u_end

#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
