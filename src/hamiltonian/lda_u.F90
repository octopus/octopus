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
       zlda_u_force,                    &
       l_notation

  character(len=1), parameter :: &
    l_notation(0:3) = (/ 's', 'p', 'd', 'f' /)

  type lda_u_t
    logical                  :: apply
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:) !> Occupation matrices and potentials 
                                                         !> for the standard scheme
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    FLOAT, pointer           :: Vloc1(:,:,:), dVloc2(:,:,:,:) !> For the corrected ACBN0 functional
    CMPLX, pointer           :: zVloc2(:,:,:,:)
    FLOAT, pointer           :: dn_alt(:,:,:,:) !> Stores the renomalized occ. matrices
    CMPLX, pointer           :: zn_alt(:,:,:,:) !> if the ACBN0 functional is used  
  
    FLOAT, pointer           :: drenorm_occ(:,:,:,:,:,:) !> On-site occupations (for the ACBN0 functional)  
    CMPLX, pointer           :: zrenorm_occ(:,:,:,:,:,:)
 
    FLOAT, pointer           :: coulomb(:,:,:,:,:) !>Coulomb integrals for all the system
                                                       !> (for the ACBN0 functional) 
    CMPLX, pointer           :: zcoulomb(:,:,:,:,:,:,:) !>Coulomb integrals for all the system
                                                        !> (for the ACBN0 functional with spinors) 

    type(orbital_set_t), pointer :: orbsets(:)   !> All the orbital setss of the system

    integer             :: norbsets           !> Number of orbital sets 
    integer             :: nspins
    integer             :: spin_channels
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
    logical             :: minimalAtomicSphere!> Use the smallest atomic orbital radius for all of them
 
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
  this%spin_channels = 0
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
  nullify(this%zcoulomb)
  nullify(this%drenorm_occ)
  nullify(this%zrenorm_occ)
  nullify(this%Vloc1)
  nullify(this%dVloc2)
  nullify(this%zVloc2)
  nullify(this%orbsets)

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

  PUSH_SUB(lda_u_init)

  ASSERT(.not. this%apply)

  call messages_print_stress(stdout, "DFT+U")
 
  if(gr%mesh%parallel_in_domains) call messages_not_implemented("dft+u parallel in domains")

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
  call messages_print_var_option(stdout, 'OrbitalsTruncationMethod', this%truncation)

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
  call messages_print_var_value(stdout, 'OrbitalsThreshold_LDAU', this%orbitals_threshold)

  !%Variable UseACBN0Functional
  !%Type logical
  !%Default no
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will determine the effective U term using the 
  !% ACBN0 functional as defined in PRX 5, 011006 (2015) 
  !%End
  call parse_variable('UseACBN0Functional', .false., this%useACBN0)
  call messages_print_var_value(stdout,  'UseACBN0Functional', this%useACBN0)

  !%Variable DFTUNormalizeOrbitals
  !%Type logical
  !%Default yes
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will normalize the atomic orbitals
  !%End
  call parse_variable('DFTUNormalizeOrbitals', .true., this%normalizeOrbitals)
  call messages_print_var_value(stdout, 'DFTUNormalizeOrbitals', this%normalizeOrbitals)

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

    !%Variable DFTUIncludeOverlap
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the overlap between orbitals on different atomic sites
    !% and use it for the ACBN0 functional
    !%End
    call parse_variable('DFTUIncludeOverlap', .false., this%IncludeOverlap)
    if(this%IncludeOverlap) call messages_experimental("DFTUIncludeOverlap")

    if(this%useAllOrbitals) then
      !%Variable DFTUMinimalAtomicSphere
      !%Type logical
      !%Default yes
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to yes, Octupus will set the radius of the orbitals to the smallest orbital
      !% present in the pseudopotential file. Only with the ACBN0 functional and the UseAllAtomicOrbitals 
      !% options activated.
      !%End
      call parse_variable('DFTUMinimalAtomicSphere', .false., this%minimalAtomicSphere)
      if(this%minimalAtomicSphere) call messages_experimental("DFTUMinimalAtomicSphere")
    end if
  end if

  !We first need to load the basis
  if (states_are_real(st)) then
    call dconstruct_orbital_basis(this, geo, gr%mesh, st)
  else
    call zconstruct_orbital_basis(this, geo, gr%mesh, st)  
  end if
  maxorbs = this%maxnorbs
  this%nspins = st%d%nspin
  this%spin_channels = st%d%spin_channels
  this%nspecies = geo%nspecies
  this%st_end = st%st_end

  !We allocate the necessary ressources
  if (states_are_real(st)) then
    call dlda_u_allocate(this, st)
  else
    call zlda_u_allocate(this, st)
  end if

  call distributed_nullify(this%orbs_dist, this%norbsets)
 #ifdef HAVE_MPI
  call distributed_init(this%orbs_dist, this%norbsets, MPI_COMM_WORLD, "orbsets")
 #endif 

  if(this%useACBN0) then
    write(message(1),'(a)')    'Computing the Coulomb integrals localized orbital basis.'
    call messages_info(1)
    if(st%d%dim == 1) then 
      if (states_are_real(st)) then
        call dcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
      else
        call zcompute_coulomb_integrals(this, gr%mesh, gr%der, st)
      end if
    else
      ASSERT(.not.states_are_real(st))
      call compute_complex_coulomb_integrals(this, gr%mesh, gr%der, st)
    end if
  end if


  call messages_print_stress(stdout)

  POP_SUB(lda_u_init)
 end subroutine lda_u_init


 subroutine lda_u_end(this)
   implicit none
   type(lda_u_t), intent(inout) :: this

   integer :: ios
  
   PUSH_SUB(lda_u_end)  

   this%apply = .false.

   SAFE_DEALLOCATE_P(this%dn)
   SAFE_DEALLOCATE_P(this%zn)
   SAFE_DEALLOCATE_P(this%dn_alt)
   SAFE_DEALLOCATE_P(this%zn_alt)
   SAFE_DEALLOCATE_P(this%dV)
   SAFE_DEALLOCATE_P(this%zV) 
   SAFE_DEALLOCATE_P(this%coulomb)
   SAFE_DEALLOCATE_P(this%zcoulomb)
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

 #ifdef HAVE_MPI
   call distributed_end(this%orbs_dist)  
 #endif

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
     call dupdate_occ_matrices(this, mesh, st, energy%dft_u)
   else
     if(associated(hm_base%phase)) then
       call zupdate_occ_matrices(this, mesh, st, energy%dft_u,&
                                hm_base%phase)
     else
       call zupdate_occ_matrices(this, mesh, st, energy%dft_u)
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

  ! ---------------------------------------------------------
  FLOAT function get_orbial_radius(geo, mesh, ia, iorb, ispin, truncation, threshold) result(radius)
    type(geometry_t), target, intent(in)   :: geo
    type(mesh_t),             intent(in)   :: mesh
    integer,                  intent(in)   :: ia, iorb, ispin
    integer,                  intent(in)   :: truncation
    FLOAT,                    intent(in)   :: threshold

    type(species_t), pointer :: spec
    integer :: ii, ll, mm

    PUSH_SUB(get_orbital_radius)

    spec => geo%atom(ia)%species
    call species_iwf_ilm(spec, iorb, ispin, ii, ll, mm)

    if(truncation == OPTION__ORBITALSTRUNCATIONMETHOD__FULL) then
      radius = species_get_iwf_radius(spec, ii, ispin, threshold)
    else
      radius = species_get_iwf_radius(spec, ii, ispin)

      if(truncation == OPTION__ORBITALSTRUNCATIONMETHOD__BOX) then
        ! if the orbital is larger than the size of the box, we restrict it to this size, 
        ! otherwise the orbital will overlap more than one time with the simulation box.
        ! This would induces phase problem if the complete mesh is used instead of the sphere
        radius = min(radius, minval(mesh%sb%lsize(1:mesh%sb%dim)-mesh%spacing(1:mesh%sb%dim)*CNST(1.01)))
      else
        !If asked, we truncate the orbital to the radius on the projector spheres 
        !of the NL part of the pseudopotential.
        !This is a way to garanty no overlap between orbitals of different atoms.
        if(species_is_ps(spec)) &
          radius = min(radius,species_get_ps_radius(spec))
        end if
        
      end if
      ! make sure that if the spacing is too large, the orbitals fit in a few points at least
      radius = max(radius, CNST(2.0)*maxval(mesh%spacing(1:mesh%sb%dim)))

    POP_SUB(get_orbital_radius)
  end function get_orbial_radius

  ! ---------------------------------------------------------
  subroutine find_minimal_atomic_spheres(geo, mesh, minradii, truncation, threshold)
    type(geometry_t), target, intent(in)    :: geo
    type(mesh_t),             intent(in)    :: mesh
    FLOAT,                  intent(inout)   :: minradii(:)
    integer,                  intent(in)    :: truncation
    FLOAT,                    intent(in)    :: threshold

    integer :: ia, iorb

    PUSH_SUB(find_minimal_atomic_spheres)

    do ia = 1, geo%natoms

     if(species_niwfs(geo%atom(ia)%species) < 1) cycle

      minradii(ia) = get_orbial_radius(geo, mesh, ia, 1, 1, truncation, threshold)
      do iorb = 2, species_niwfs(geo%atom(ia)%species)
        minradii(ia) = min(minradii(ia), get_orbial_radius(geo, mesh, ia, iorb, 1, truncation, threshold))
      end do !iorb
    end do !ia

    POP_SUB(find_minimal_atomic_spheres) 
  end subroutine find_minimal_atomic_spheres

! ---------------------------------------------------------
! TODO: Merge this with the two_body routine in system/output_me_inc.F90
subroutine compute_complex_coulomb_integrals (this, mesh, der, st)
  type(lda_u_t),   intent(inout)  :: this
  type(mesh_t),       intent(in)  :: mesh
  type(derivatives_t), intent(in) :: der
  type(states_t),     intent(in)  :: st

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: is1, is2
  integer :: norbs, np_sphere, ios, ip
  integer :: idone, ntodo
  CMPLX, allocatable :: tmp(:), vv(:,:), nn(:,:)
  type(orbital_set_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, "DFTU_COMPEX_COULOMB_INTEGRALS")

  PUSH_SUB(compute_complex_coulomb_integrals_complex)

  ASSERT(.not. st%parallel_in_states)

  SAFE_ALLOCATE(nn(1:this%max_np,st%d%dim))
  SAFE_ALLOCATE(vv(1:this%max_np,st%d%dim))
  SAFE_ALLOCATE(tmp(1:this%max_np))

  SAFE_ALLOCATE(this%zcoulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:st%d%dim, 1:st%d%dim, 1:this%norbsets))
  this%zcoulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:st%d%dim,1:st%d%dim,1:this%norbsets) = M_ZERO

  !Lets counts the number of orbital to treat, to display a progress bar
  ntodo = 0
  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    ntodo= ntodo + norbs**4*4!((norbs+1)*norbs/2)*((norbs+1)*norbs/2+1)/2
  end do
  idone = 0
  if(mpi_world%rank == 0) call loct_progress_bar(-1, ntodo)


  do ios = this%orbs_dist%start, this%orbs_dist%end
    os => this%orbsets(ios)
    norbs = os%norbs
    np_sphere = os%sphere%np

    call poisson_init_sm(os%poisson, psolver, der, os%sphere)

    ijst=0
    do ist = 1, norbs

      do jst = 1, norbs
       ! if(jst > ist) cycle
        ijst=ijst+1

        do is1 = 1, st%d%dim
          !$omp parallel do
          do ip=1,np_sphere
            nn(ip,is1)  = os%zorb(ip,is1,ist)*os%zorb(ip,is1,jst)
          end do
          !$omp end parallel do    

          !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
          call zpoisson_solve_sm(os%poisson, os%sphere, vv(1:np_sphere,is1), nn(1:np_sphere,is1))
        end do !is1  

        klst=0
        do kst = 1, norbs
          do lst = 1, norbs
       !     if(lst > kst) cycle
            klst=klst+1
       !     if(klst > ijst) cycle

            do is1 = 1, st%d%dim
              do is2 = 1, st%d%dim

                !$omp parallel do
                do ip=1,np_sphere
                 tmp(ip) = vv(ip,is1)*os%zorb(ip,is2,lst)*os%zorb(ip,is2,kst)
                end do
                !$omp end parallel do

                this%zcoulomb(ist,jst,kst,lst,is1,is2,ios) = zsm_integrate(mesh, os%sphere, tmp(1:np_sphere))
              end do !is2
            end do !is1

            do is1 = 1, st%d%dim
              do is2 = 1, st%d%dim
                if(abs(this%zcoulomb(ist,jst,kst,lst,is1,is2,ios))<CNST(1.0e-12)) then
                  this%zcoulomb(ist,jst,kst,lst,is1,is2,ios) = M_ZERO
           !     else
           !       this%zcoulomb(kst,lst,ist,jst,is2,is1,ios) = this%zcoulomb(ist,jst,kst,lst,is1,is2,ios)
           !       this%zcoulomb(jst,ist,lst,kst,is2,is1,ios) = conjg(this%zcoulomb(ist,jst,kst,lst,is1,is2,ios))
           !       this%zcoulomb(lst,kst,jst,ist,is2,is1,ios) = conjg(this%zcoulomb(ist,jst,kst,lst,is1,is2,ios))
                end if

              end do !is2
            end do !is1

            idone = idone + 1
            if(mpi_world%rank == 0) call loct_progress_bar(idone, ntodo)
          end do !lst
        end do !kst
      end do !jst
    end do !ist
    call poisson_end(os%poisson)
  end do !iorb

  if(this%orbs_dist%parallel) then
    do ios = 1, this%norbsets
      do is2 = 1, st%d%dim
        do is1 = 1, st%d%dim
          call comm_allreduce(this%orbs_dist%mpi_grp%comm, this%zcoulomb(:,:,:,:,is1,is2,ios))
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(compute_complex_coulomb_integrals)
  call profiling_out(prof)
end subroutine compute_complex_coulomb_integrals


#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
