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

#include "global.h"

module lda_u_oct_m
  use atomic_orbital_oct_m
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
  use loewdin_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use orbitalset_utils_oct_m
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
       lda_u_write_info


  type lda_u_t
    integer                  :: level
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:) !> Occupation matrices and potentials 
                                                         !> for the standard scheme
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    FLOAT, pointer           :: dn_alt(:,:,:,:) !> Stores the renomalized occ. matrices
    CMPLX, pointer           :: zn_alt(:,:,:,:) !> if the ACBN0 functional is used  
  
    FLOAT, pointer           :: drenorm_occ(:,:,:,:,:) !> On-site occupations (for the ACBN0 functional)  
    CMPLX, pointer           :: zrenorm_occ(:,:,:,:,:)
 
    FLOAT, pointer           :: coulomb(:,:,:,:,:) !>Coulomb integrals for all the system
                                                       !> (for the ACBN0 functional) 
    CMPLX, pointer           :: zcoulomb(:,:,:,:,:,:,:) !>Coulomb integrals for all the system
                                                        !> (for the ACBN0 functional with spinors) 

    type(orbitalbasis_t) :: basis               !> The full basis of localized orbitals
    type(orbitalset_t), pointer :: orbsets(:)   !> All the orbital setss of the system
    integer             :: norbsets

    integer             :: nspins
    integer             :: spin_channels
    integer             :: nspecies        
    integer             :: maxnorbs           !> Maximal number of orbitals for all the atoms
    integer             :: max_np             !> Maximum number of points in all orbitals submesh spheres 
 
    logical             :: useAllOrbitals     !> Do we use all atomic orbitals possible
    logical             :: skipSOrbitals      !> Not using s orbitals
    logical             :: freeze_occ         !> Occupation matrices are not recomputed during TD evolution
    logical             :: freeze_u           !> U is not recomputed during TD evolution

    type(distributed_t) :: orbs_dist
  end type lda_u_t

  integer, public, parameter ::        &
    DFT_U_NONE                    = 0, &
    DFT_U_EMPIRICAL               = 1, &
    DFT_U_ACBN0                   = 2

contains

 subroutine lda_u_nullify(this)
  type(lda_u_t),             intent(inout) :: this

  PUSH_SUB(lda_u_nullify)

  this%level = DFT_U_NONE
  this%norbsets = 0
  this%max_np = 0
  this%maxnorbs = 0
  this%nspins = 0
  this%spin_channels = 0
  this%nspecies = 0
  this%useAllOrbitals = .false.
  this%skipSOrbitals = .true.
  this%freeze_occ = .false.
  this%freeze_u = .false.

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
  nullify(this%orbsets)

  call distributed_nullify(this%orbs_dist, 0)

  call orbitalbasis_nullify(this%basis)

  POP_SUB(lda_u_nullify)

 end subroutine lda_u_nullify

 subroutine lda_u_init(this, level, gr, geo, st)
  type(lda_u_t),             intent(inout) :: this
  integer,                   intent(in)    :: level
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st

  logical :: complex_coulomb_integrals
  integer :: ios

  PUSH_SUB(lda_u_init)

  ASSERT(.not. (level == DFT_U_NONE))

  call messages_print_stress(stdout, "DFT+U")
  if(gr%mesh%parallel_in_domains) call messages_experimental("dft+u parallel in domains")
  this%level = level
  
  call lda_u_write_info(this, stdout)

  call orbitalbasis_init(this%basis)

  if( this%level == DFT_U_ACBN0 ) then
    !%Variable UseAllAtomicOrbitals
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will determine the effective U for all atomic orbitals
    !% from the peusopotential. Only available with ACBN0 functional.
    !% It is strongly recommended to set AOLoewdin=yes when using the option.
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
  end if

  if (states_are_real(st)) then
    call dorbitalbasis_build(this%basis, geo, gr%mesh, st%d%kpt, st%d%dim, &
                             this%skipSOrbitals, this%useAllOrbitals)
  else
    call zorbitalbasis_build(this%basis, geo, gr%mesh, st%d%kpt, st%d%dim, &
                             this%skipSOrbitals, this%useAllOrbitals)
  end if
  this%orbsets => this%basis%orbsets
  this%norbsets = this%basis%norbsets
  this%maxnorbs = this%basis%maxnorbs
  this%max_np = this%basis%max_np
  this%nspins = st%d%nspin
  this%spin_channels = st%d%spin_channels
  this%nspecies = geo%nspecies

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


  if( this%level == DFT_U_ACBN0 ) then

    complex_coulomb_integrals = .false.
    do ios = 1, this%norbsets
        if(this%orbsets(ios)%ndim  > 1) complex_coulomb_integrals = .true.
    end do

    call messages_info(1)
    if(.not. complex_coulomb_integrals) then 
      write(message(1),'(a)')    'Computing the Coulomb integrals of the localized basis.'
      if (states_are_real(st)) then
        call dcompute_coulomb_integrals(this, gr%mesh, gr%der)
      else
        call zcompute_coulomb_integrals(this, gr%mesh, gr%der)
      end if
    else
      ASSERT(.not.states_are_real(st))
      write(message(1),'(a)')    'Computing complex Coulomb integrals of the localized basis.'
      call compute_complex_coulomb_integrals(this, gr%mesh, gr%der, st)
    end if
  end if


  call messages_print_stress(stdout)

  POP_SUB(lda_u_init)
 end subroutine lda_u_init


 subroutine lda_u_end(this)
   implicit none
   type(lda_u_t), intent(inout) :: this

   PUSH_SUB(lda_u_end)  

   this%level = DFT_U_NONE

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

   nullify(this%orbsets)
   call orbitalbasis_end(this%basis)

   this%max_np = 0

 #ifdef HAVE_MPI
   call distributed_end(this%orbs_dist)  
 #endif

   POP_SUB(lda_u_end)
 end subroutine lda_u_end

 ! When moving the ions, the basis must be reconstructed
 subroutine lda_u_update_basis(this, gr, geo, st, has_phase)
  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st
  logical,                   intent(in)    :: has_phase

  if(this%level == DFT_U_NONE) return

  PUSH_SUB(lda_u_update_basis)

  !We clean the orbital basis, to be able to reconstruct it
  call orbitalbasis_end(this%basis)
  nullify(this%orbsets)

  !We now reconstruct the basis
  if (states_are_real(st)) then
    call dorbitalbasis_build(this%basis, geo, gr%mesh, st%d%kpt, st%d%dim, &
                             this%skipSOrbitals, this%useAllOrbitals, verbose = .false.)
  else
    call zorbitalbasis_build(this%basis, geo, gr%mesh, st%d%kpt, st%d%dim, &
                             this%skipSOrbitals, this%useAllOrbitals, verbose = .false.)
  end if
  this%orbsets => this%basis%orbsets

  ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
  ! In case of a laser field, the phase is recomputed in hamiltonian_update
  if(has_phase) then
    call lda_u_build_phase_correction(this, gr%mesh%sb, st%d)
  end if

  POP_SUB(lda_u_update_basis)

 end subroutine lda_u_update_basis

 ! Interface for the X(update_occ_matrices) routines
 subroutine lda_u_update_occ_matrices(this, mesh, st, hm_base, energy )
   type(lda_u_t),             intent(inout) :: this
   type(mesh_t),              intent(in)    :: mesh 
   type(states_t),            intent(in)    :: st
   type(hamiltonian_base_t),  intent(in)    :: hm_base 
   type(energy_t),            intent(inout) :: energy

   if(this%level == DFT_U_NONE .or. this%freeze_occ) return
   PUSH_SUB(lda_u_update_occ_matrices)

   if (states_are_real(st)) then
     call dupdate_occ_matrices(this, mesh, st, energy%dft_u)
   else
     if(associated(hm_base%phase)) then
       call zupdate_occ_matrices(this, mesh, st, energy%dft_u, hm_base%phase)
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
     call orbitalset_update_phase(this%orbsets(ios), sb, std%kpt, (std%ispin==SPIN_POLARIZED), &
                                        vec_pot, vec_pot_var)
   end do

   if(this%basis%orthogonalization) then
     call zloewdin_orthogonalize(this%basis, std%kpt)
   else
     if(debug%info) call zloewdin_info(this%basis, std%kpt)
   end if
  
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


  subroutine lda_u_write_info(this, iunit)
    type(lda_u_t),  intent(in)    :: this
    integer,        intent(in)    :: iunit

    PUSH_SUB(lda_u_write_info)

    write(message(1), '(1x)')
    call messages_info(1, iunit)
    if(this%level == DFT_U_EMPIRICAL) then
      write(message(1), '(a)') "Method:"
      write(message(2), '(a)') "  [1] Dudarev et al., Phys. Rev. B 57, 1505 (1998)"
      call messages_info(2, iunit)
    else
      write(message(1), '(a)') "Method:"
      write(message(2), '(a)') "  [1] Agapito et al., Phys. Rev. X 5, 011006 (2015)"
      call messages_info(2, iunit)
    end if
    write(message(1), '(a)') "Implementation:"
    write(message(2), '(a)') "  [1] Tancogne-Dejean, Oliveira, and Rubio, Phys. Rev. B 69, 245133 (2017)"
    write(message(3), '(1x)')
    call messages_info(3, iunit)

    POP_SUB(lda_u_write_info)

  end subroutine lda_u_write_info

#include "dft_u_noncollinear_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
