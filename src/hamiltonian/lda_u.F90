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
  use derivatives_oct_m
  use distributed_oct_m
  use energy_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_base_oct_m 
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use mpi_oct_m
  use orbital_oct_m
  use orbital_set_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use species_pot_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use types_oct_m  
  use unit_oct_m
  use unit_system_oct_m
 
  implicit none

  private

  public ::                             &
       lda_u_t,                         &
       lda_u_nullify,                   &
       lda_u_init,                      &
       dlda_u_apply,                    &
       zlda_u_apply,                    &
       lda_u_update_basis,              &
       lda_u_update_occ_matrices,       &
       lda_u_end,                       &
       lda_u_build_phase_correction,    &
       lda_u_freeze_occ,                &
       lda_u_freeze_u,                  &
       dlda_u_set_occupations,          &
       zlda_u_set_occupations,          &
       dlda_u_get_occupations,          &
       zlda_u_get_occupations,          &
       dlda_u_update_potential,         &
       zlda_u_update_potential,         &
       lda_u_get_effectiveU,            &
       lda_u_set_effectiveU,            &
       dlda_u_commute_r,                &
       zlda_u_commute_r,                &
       dlda_u_force,                    &
       zlda_u_force


  type lda_u_t
    logical                  :: apply
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:) !> Occupation matrices and potentials 
                                                         !> for the standard scheme
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    FLOAT, pointer           :: Vloc1(:,:,:), dVloc2(:,:,:,:) !> For the corrected ACBN0 functional
    CMPLX, pointer           :: zVloc2(:,:,:,:)
    FLOAT, pointer           :: dn_alt(:,:,:,:) !> Stores the renomalized occ. matrices
    CMPLX, pointer           :: zn_alt(:,:,:,:) !> if the ACBN0 functional is used  
  
    FLOAT, pointer           :: drenorm_occ(:,:,:,:,:) !> On-site occupations (for the ACBN0 functional)  
    CMPLX, pointer           :: zrenorm_occ(:,:,:,:,:)
 
    FLOAT, pointer           :: coulomb(:,:,:,:,:) !>Coulomb integrals for all the system
                                                   !> (for the ACBN0 functional) 
 
    type(orbital_set_t), pointer :: orbsets(:)   !> All the orbital setss of the system

    integer             :: norbsets           !> Number of orbital sets 
    integer             :: nspins
    integer             :: nspecies        
    integer             :: st_end
    integer             :: maxnorbs           !> Maximal number of orbitals for all the atoms
    integer             :: max_np             !> Maximum number of points in all orbitals submesh spheres 
 
    integer             :: truncation         !> Truncation method for the orbitals
    FLOAT               :: orbitals_threshold !> Threshold for orbital truncation
    logical             :: useACBN0           !> Do we use the ACBN0 functional
    logical             :: ACBN0_corrected    !> Do we take into account missing terms from the ACBN0 original paper
    logical             :: useAllOrbitals     !> Do we use all atomic orbitals possible
    logical             :: skipSOrbitals      !> Not using s orbitals
    logical             :: IncludeOverlap     !> Do we compute and use overlap or not
    logical             :: freeze_occ         !> Occupation matrices are not recomputed during TD evolution
    logical             :: freeze_u           !> U is not recomputed during TD evolution
    logical             :: normalizeOrbitals  !> Do we normalize the orbitals 
 
    type(distributed_t) :: orbs_dist
  end type lda_u_t

contains

 subroutine lda_u_nullify(this)
  type(lda_u_t),             intent(inout) :: this

  PUSH_SUB(lda_u_nullify)

  this%apply = .false.
  this%norbsets = 0
  this%max_np = 0
  this%maxnorbs = 0
  this%nspins = 0
  this%nspecies = 0
  this%freeze_occ = .false.
  this%freeze_u = .false.
  this%IncludeOverlap = .false.

  nullify(this%dn)
  nullify(this%zn)
  nullify(this%dn_alt)
  nullify(this%zn_alt)
  nullify(this%dV)
  nullify(this%zV)
  nullify(this%coulomb)
  nullify(this%drenorm_occ)
  nullify(this%zrenorm_occ)
  nullify(this%Vloc1)
  nullify(this%dVloc2)
  nullify(this%zVloc2)

  call distributed_nullify(this%orbs_dist, 0)

  POP_SUB(lda_u_nullify)

 end subroutine lda_u_nullify

 subroutine lda_u_init(this, gr, geo, st, mc)
  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st
  type(multicomm_t),         intent(in)    :: mc

  integer :: maxorbs, iat, ispin, iorb
  real(8) :: mem, coef

  PUSH_SUB(lda_u_init)

  ASSERT(.not. this%apply)

  call messages_print_stress(stdout, "DFT+U")
 
!  if(st%parallel_in_states) call messages_not_implemented("lda+u parallel in states")
  if(gr%mesh%parallel_in_domains) call messages_not_implemented("dft+u parallel in domains")
  if(st%d%ispin == SPINORS) call messages_not_implemented("dft+u with spinors") 

  this%apply = .true.

  !%Variable OrbitalsTruncationMethod
  !%Type flag
  !%Default full
  !%Section Hamiltonian::DFT+U
  !%Description
  !% This option determine how Octopus will truncate the orbitals used for LDA+U.
  !% Except for the full method, the other options are only there to get a quick idea.
  !%Option full bit(0)
  !% The full size of the orbitals used. The radius is controled by variable OrbitalThreshold_LDAU
  !%Option box bit(1)
  !% The radius of the orbitals are restricted to the size of the simulation box. 
  !% This reduces the number of points used to discretize the orbitals.
  !%Option NLradius bit(2)
  !% The radius of the orbitals are restricted to the radius of the non-local part of the pseudopotential 
  !% of the corresponding atom.
  !%End
  call parse_variable('OrbitalsTruncationMethod', OPTION__ORBITALSTRUNCATIONMETHOD__FULL, this%truncation)

  !%Variable OrbitalsThreshold_LDAU
  !%Type float
  !%Default 0.01
  !%Section Hamiltonian::DFT+U
  !%Description
  !% Determine the threshold used to compute the radius of the atomic orbitals for LDA+U.
  !% This radius is computed by making sure that the 
  !% absolute value of the radial part of the atomic orbital is below the specified threshold.
  !% This value should be converged to be sure that results do not depend on this value. 
  !% However increasing this value increases the number of grid points covered by the orbitals and directly affect performances.
  !%End
  call parse_variable('OrbitalsThreshold_LDAU', CNST(0.01), this%orbitals_threshold)
  if(this%orbitals_threshold <= M_ZERO) call messages_input_error('OrbitalsThreshold_LDAU')

  !%Variable UseACBN0Functional
  !%Type logical
  !%Default no
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will determine the effective U term using the 
  !% ACBN0 functional as defined in PRX 5, 011006 (2015) 
  !%End
  call parse_variable('UseACBN0Functional', .false., this%useACBN0)

  !%Variable DFTUNormalizeOrbitals
  !%Type logical
  !%Default no
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will normalize the atomic orbitals
  !%End
  call parse_variable('DFTUNormalizeOrbitals', .false., this%normalizeOrbitals)

  if( this%useACBN0) then
    !%Variable UseAllAtomicOrbitals
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the effective U for all atomic orbitals
    !% from the peusopotential. Only available with ACBN0 functional.
    !%End
    call parse_variable('UseAllAtomicOrbitals', .false., this%useAllOrbitals)
    if(this%useAllOrbitals) call messages_experimental("UseAllAtomicOrbitals")

    !%Variable SkipSOrbitals
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the effective U for all atomic orbitals
    !% from the peusopotential but s orbitals. Only available with ACBN0 functional.
    !%End
    call parse_variable('SkipSOrbitals', .true., this%skipSOrbitals)   
    if(.not.this%SkipSOrbitals) call messages_experimental("SkipSOrbitals")

    !%Variable ACBN0_corrected
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the effective U term using the 
    !% ACBN0 functional, including the missing term from the derivative of the effective U.
    !%End
    call parse_variable('ACBN0_corrected', .false., this%ACBN0_corrected)
    if(this%ACBN0_corrected) call messages_experimental("ACBN0_corrected")

    !%Variable IncludeOverlap
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the overlap between orbitals on different atomic sites
    !% and use it for the ACBN0 functional
    !%End
    call parse_variable('IncludeOverlap', .false., this%IncludeOverlap)
    if(this%IncludeOverlap) call messages_experimental("IncludeOverlap")
  end if

  !We first need to load the basis
  if (states_are_real(st)) then
    call dconstruct_orbital_basis(this, geo, gr%mesh, st)
  else
    call zconstruct_orbital_basis(this, geo, gr%mesh, st)  
  end if
  maxorbs = this%maxnorbs
  this%nspins = st%d%nspin
  this%nspecies = geo%nspecies
  this%st_end = st%st_end

  !We analyse the memeory and we print the requiered memory
  !Thus, if there is not enough memory, the user knows with the code crashes
  mem = 0.0_8
  coef = 2.0_8
  if (states_are_real(st)) coef = 1.0_8
  mem = mem + coef*REAL_PRECISION*dble(maxorbs**2*st%d%nspin*this%norbsets*2) !Occupation matrices and potentials
  mem = mem + coef*REAL_PRECISION*dble(maxorbs*st%d%nspin*this%norbsets)    !Orbital occupations
  if(this%useACBN0) then
    mem = mem + REAL_PRECISION*dble(maxorbs**4*st%d%nspin*this%norbsets) !Coulomb intergrals
    mem = mem + REAL_PRECISION*dble(18*(st%d%kpt%end-st%d%kpt%start+1)*(st%st_end-st%st_start+1)) !On-site occupations
  end if
  call messages_new_line()
  call messages_write('    Approximate memory requirement for LDA+U (for each task)   :')
  call messages_write(mem, units = unit_megabytes)
  call messages_new_line()
  call messages_info()

  !We allocate the necessary ressources
  if (states_are_real(st)) then
    SAFE_ALLOCATE(this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
    this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO
    SAFE_ALLOCATE(this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
    this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%useACBN0) then
      SAFE_ALLOCATE(this%dn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
      this%dn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO
      SAFE_ALLOCATE(this%drenorm_occ(geo%nspecies,0:5,0:3,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end))
      this%drenorm_occ(geo%nspecies,0:5,0:3,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end) = M_ZERO
      if(this%ACBN0_corrected) then
        SAFE_ALLOCATE(this%Vloc1(1:maxorbs,1:st%d%nspin,1:this%norbsets))
        this%Vloc1(1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO
        SAFE_ALLOCATE(this%dVloc2(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
        this%dVloc2(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO
      end if
    end if
  else
    SAFE_ALLOCATE(this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
    this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = cmplx(M_ZERO,M_ZERO)
    SAFE_ALLOCATE(this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
    this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = cmplx(M_ZERO,M_ZERO)

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%useACBN0) then
      SAFE_ALLOCATE(this%zn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
      this%zn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = cmplx(M_ZERO,M_ZERO)
      SAFE_ALLOCATE(this%zrenorm_occ(geo%nspecies,0:5,0:3,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end))
      this%zrenorm_occ(geo%nspecies,0:5,0:3,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end) = M_ZERO
      if(this%ACBN0_corrected) then
        SAFE_ALLOCATE(this%Vloc1(1:maxorbs,1:st%d%nspin,1:this%norbsets))
        this%Vloc1(1:maxorbs,1:st%d%nspin,1:this%norbsets) = M_ZERO
        SAFE_ALLOCATE(this%zVloc2(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets))
        this%zVloc2(1:maxorbs,1:maxorbs,1:st%d%nspin,1:this%norbsets) = cmplx(M_ZERO,M_ZERO)
      end if
    end if
  end if

  call distributed_nullify(this%orbs_dist, this%norbsets)
  call distributed_init(this%orbs_dist, this%norbsets, MPI_COMM_WORLD, "orbsets")
 
  if(this%useACBN0) then
    write(message(1),'(a)')    'Computing the Coulomb integrals localized orbital basis.'
    call messages_info(1) 
    if (states_are_real(st)) then
      call dcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
    else
      call zcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
    end if
  end if


  call messages_print_stress(stdout)

  POP_SUB(lda_u_init)
 end subroutine lda_u_init


 subroutine lda_u_end(this)
   implicit none
   type(lda_u_t), intent(inout) :: this

   integer :: ios, ispin, iorb
  
   PUSH_SUB(lda_u_end)  

   this%apply = .false.

   SAFE_DEALLOCATE_P(this%dn)
   SAFE_DEALLOCATE_P(this%zn)
   SAFE_DEALLOCATE_P(this%dn_alt)
   SAFE_DEALLOCATE_P(this%zn_alt)
   SAFE_DEALLOCATE_P(this%dV)
   SAFE_DEALLOCATE_P(this%zV) 
   SAFE_DEALLOCATE_P(this%coulomb)
   SAFE_DEALLOCATE_P(this%drenorm_occ)
   SAFE_DEALLOCATE_P(this%zrenorm_occ)
   SAFE_DEALLOCATE_P(this%Vloc1)
   SAFE_DEALLOCATE_P(this%dVloc2)
   SAFE_DEALLOCATE_P(this%zVloc2)

   do ios = 1, this%norbsets
     call orbital_set_end(this%orbsets(ios))
   end do

   SAFE_DEALLOCATE_P(this%orbsets)

   this%max_np = 0
   this%norbsets = 0

   call distributed_end(this%orbs_dist)  

   POP_SUB(lda_u_end)
 end subroutine lda_u_end

 ! When moving the ions, the basis must be reconstructed
 subroutine lda_u_update_basis(this, gr, geo, st)
  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st

  integer :: iorbset

  if(.not. this%apply) return

  PUSH_SUB(lda_u_update_basis)

  call messages_print_stress(stdout, "Updating DFT+U basis")

  !We clean the orbital basis, to be able to reconstruct it
  do iorbset = 1, this%norbsets
    call orbital_set_end(this%orbsets(iorbset))
  end do
  SAFE_DEALLOCATE_P(this%orbsets)

  !We now reconstruct the basis
  if (states_are_real(st)) then
    call dconstruct_orbital_basis(this, geo, gr%mesh, st)
  else
    call zconstruct_orbital_basis(this, geo, gr%mesh, st)
  end if 

  !if(this%useACBN0) then
  !  write(message(1),'(a)')    'Computing the Coulomb integrals localized orbital basis.'
  !  call messages_info(1)
  !  if (states_are_real(st)) then
  !    call dcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
  !  else
  !    call zcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
  !  end if
  !end if


  call messages_print_stress(stdout)

  POP_SUB(lda_u_update_basis)

 end subroutine lda_u_update_basis

 ! Interface for the X(update_occ_matrices) routines
 subroutine lda_u_update_occ_matrices(this, mesh, st, hm_base, energy )
   type(lda_u_t),             intent(inout) :: this
   type(mesh_t),              intent(in)    :: mesh 
   type(states_t),            intent(in)    :: st
   type(hamiltonian_base_t),  intent(in)    :: hm_base 
   type(energy_t),            intent(inout) :: energy

   if(.not. this%apply .or. this%freeze_occ) return
   PUSH_SUB(lda_u_update_occ_matrices)

   if (states_are_real(st)) then
     call dupdate_occ_matrices(this, mesh, st, energy%lda_u_energy)
   else
     if(associated(hm_base%phase)) then
       call zupdate_occ_matrices(this, mesh, st, energy%lda_u_energy,&
                                hm_base%phase)
     else
       call zupdate_occ_matrices(this, mesh, st, energy%lda_u_energy)
     end if
   end if

   POP_SUB(lda_u_update_occ_matrices)
 end subroutine lda_u_update_occ_matrices


 !> Build the phase correction to the global phase for all orbitals
 subroutine lda_u_build_phase_correction(this, sb, std, vec_pot, vec_pot_var)
   type(lda_u_t),                 intent(inout) :: this
   type(simul_box_t),             intent(in)    :: sb 
   type(states_dim_t),            intent(in)    :: std
   FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
   FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)

   integer :: ios
 
   PUSH_SUB(lda_u_build_phase_correction)

   do ios = 1, this%norbsets
     call  orbital_set_update_phase(this%orbsets(ios), sb, std, vec_pot, vec_pot_var)
   end do
  
   POP_SUB(lda_u_build_phase_correction)

 end subroutine lda_u_build_phase_correction

  subroutine lda_u_freeze_occ(this) 
    type(lda_u_t),     intent(inout) :: this

    this%freeze_occ = .true.
  end subroutine lda_u_freeze_occ

  subroutine lda_u_freeze_u(this)            
    type(lda_u_t),     intent(inout) :: this

    this%freeze_u = .true.
  end subroutine lda_u_freeze_u

    ! ---------------------------------------------------------
  subroutine lda_u_set_effectiveU(this, Ueff)
    type(lda_u_t),  intent(inout) :: this
    FLOAT,          intent(in)    :: Ueff(:) !< (this%norbsets)

    integer :: ios

    PUSH_SUB(lda_u_set_effectiveU)

    do ios = 1,this%norbsets
      this%orbsets(ios)%Ueff = Ueff(ios)
    end do
    POP_SUB(lda_u_set_effectiveU)
  end subroutine lda_u_set_effectiveU

  ! ---------------------------------------------------------
  subroutine lda_u_get_effectiveU(this, Ueff)
    type(lda_u_t),  intent(in)    :: this
    FLOAT,          intent(inout) :: Ueff(:) !< (this%norbsets)

    integer :: ios

    PUSH_SUB(lda_u_get_effectiveU)

    do ios = 1,this%norbsets
      Ueff(ios) = this%orbsets(ios)%Ueff
    end do
    POP_SUB(lda_u_get_effectiveU)
  end subroutine lda_u_get_effectiveU


#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
