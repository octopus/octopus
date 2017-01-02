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
  use energy_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_base_oct_m 
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
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
       lda_u_init,                      &
       dlda_u_apply,                    &
       zlda_u_apply,                    &
       lda_u_update_occ_matrices,       &
       lda_u_end,                       &
       lda_u_build_phase_correction,    &
       lda_u_update_U,                  &
       lda_u_freeze_occ,                &
       lda_u_freeze_u

  type orbital_t
    type(submesh_t)     :: sphere             !> The submesh of the orbital
    FLOAT, pointer      :: dorbital_sphere(:) !> The orbital, if real, on the submesh
    CMPLX, pointer      :: zorbital_sphere(:) !> The orbital, if complex, on the submesh
    CMPLX, pointer      :: phase(:,:)         !> Correction to the global phase 
                                              !> if the sphere cross the border of the box
   integer              :: ll                 !> Angular momentum of the orbital
  end type orbital_t


  type lda_u_t
    logical                  :: apply
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:) !> Occupation matrices and potentials 
                                                         !> for the standard scheme
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    FLOAT, pointer           :: dn_alt(:,:,:,:) !> Stores the renomalized occ. matrices
    CMPLX, pointer           :: zn_alt(:,:,:,:) !> if the ACBN0 functional is used  
  
    FLOAT, pointer           :: renorm_occ(:,:,:) !> On-site occupations (for the ACBN0 functional)  
 
    FLOAT, pointer           :: coulomb(:,:,:,:,:) !>Coulomb integrals for all the system
                                                   !> (for the ACBN0 functional) 
 
    type(orbital_t), pointer :: orbitals(:,:) !>An array containing all the orbitals of the system
    FLOAT, pointer           :: Ueff(:)       !> The effective U of the simplified rotational invariant form

    integer             :: natoms             !> Number of atoms (copied from geometry_t)
    integer             :: nspins
    integer, pointer    :: norbs(:)           !> Number of orbitals
    integer             :: maxnorbs           !> Maximal number of orbitals for all the atoms
    integer             :: max_np             !> Maximum number of points in all orbitals submesh spheres 
 
    integer             :: truncation         !> Truncation method for the orbitals
    FLOAT               :: orbitals_threshold !> Threshold for orbital truncation
    logical             :: useACBN0           !> Do we use the ACBN0 functional
    logical             :: freeze_occ         !> Occupation matrices are not recomputed during TD evolution
    logical             :: freeze_u           !> U is not recomputed during TD evolution
  end type lda_u_t

contains

 subroutine lda_u_init(this, gr, geo, st)
  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st

  integer :: maxorbs, iat, ispin, iorb
  real(8) :: mem, coef

  PUSH_SUB(lda_u_init)

  ASSERT(.not. this%apply)

  call messages_print_stress(stdout, "LDA+U")
 
  if(st%parallel_in_states) call messages_not_implemented("lda+u parallel in states")
  if(st%d%ispin == SPINORS) call messages_not_implemented("lda+u with spinors") 

  this%apply = .true.

  !%Variable OrbitalsTruncationMethod
  !%Type flag
  !%Default full
  !%Section Hamiltonian::LDA+U
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
  !%Section Hamiltonian::LDA+U
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
  !%Section Hamiltonian::LDA+U
  !%Description
  !% If set to yes, Octopus will determine the effective U term using the 
  !% ACBN0 functional as defined in PRX 5, 011006 (2015) 
  !%End
  call parse_variable('UseACBN0Functional', .false., this%useACBN0)

  nullify(this%dn)
  nullify(this%zn)
  nullify(this%dn_alt)
  nullify(this%zn_alt)
  nullify(this%dV)
  nullify(this%zV)
  nullify(this%Ueff)
  nullify(this%coulomb) 
  nullify(this%renorm_occ)

  this%natoms = geo%natoms
  this%nspins = st%d%nspin
  this%max_np = 0
 
  !We first need to load the basis
  if (states_are_real(st)) then
    call dconstruct_orbital_basis(this, geo, gr%mesh, st)
  else
    call zconstruct_orbital_basis(this, geo, gr%mesh, st)  
  end if
  maxorbs = maxval(this%norbs)

  !We analyse the memeory and we print the requiered memory
  !Thus, if there is not enough memory, the user knows with the code crashes
  mem = 0.0_8
  coef = 2.0_8
  if (states_are_real(st)) coef = 1.0_8
  mem = mem + coef*REAL_PRECISION*dble(maxorbs**2*st%d%nspin*geo%natoms*2) !Occupation matrices and potentials
  mem = mem + coef*REAL_PRECISION*dble(maxorbs*st%d%nspin*geo%natoms)    !Orbital occupations
  if(this%useACBN0) then
    mem = mem + REAL_PRECISION*dble(maxorbs**4*st%d%nspin*geo%natoms) !Coulomb intergrals
    mem = mem + REAL_PRECISION*dble(10*(st%d%kpt%end-st%d%kpt%start+1)*(st%st_end-st%st_start+1)) !On-site occupations
  end if
  call messages_new_line()
  call messages_write('    Approximate memory requirement for LDA+U (for each task)   :')
  call messages_write(mem, units = unit_megabytes)
  call messages_new_line()
  call messages_info()

  !We allocate the necessary ressources
  if (states_are_real(st)) then
    SAFE_ALLOCATE(this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%dn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = M_ZERO
    SAFE_ALLOCATE(this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%dV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = M_ZERO

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%useACBN0) then
      SAFE_ALLOCATE(this%dn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
      this%dn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = M_ZERO
    end if
  else
    SAFE_ALLOCATE(this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%zn(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = cmplx(M_ZERO,M_ZERO)
    SAFE_ALLOCATE(this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
    this%zV(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = cmplx(M_ZERO,M_ZERO)

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%useACBN0) then
      SAFE_ALLOCATE(this%zn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms))
      this%zn_alt(1:maxorbs,1:maxorbs,1:st%d%nspin,1:geo%natoms) = cmplx(M_ZERO,M_ZERO)
    end if
  end if
  SAFE_ALLOCATE(this%renorm_occ(10,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end))
  this%renorm_occ(10,st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end) = M_ZERO 
  !In case we use the ab-initio scheme, we need to allocate extra resources
  if(this%useACBN0) then
    SAFE_ALLOCATE(this%coulomb(1:maxorbs,1:maxorbs,1:maxorbs,1:maxorbs,1:geo%natoms))
    this%coulomb(1:maxorbs,1:maxorbs,1:maxorbs,1:maxorbs,1:geo%natoms) = M_ZERO
  end if

 
  if(this%useACBN0) then
    write(message(1),'(a)')    'Computing the Coulomb integrals localized orbital basis.'
    call messages_info(1) 
    if (states_are_real(st)) then
      call dcompute_coulomb_integrals(this, gr%mesh, st)
    else
      call zcompute_coulomb_integrals(this, gr%mesh, st)
    end if
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
   SAFE_DEALLOCATE_P(this%dn_alt)
   SAFE_DEALLOCATE_P(this%zn_alt)
   SAFE_DEALLOCATE_P(this%dV)
   SAFE_DEALLOCATE_P(this%zV) 
   SAFE_DEALLOCATE_P(this%Ueff)
   SAFE_DEALLOCATE_P(this%coulomb)
   SAFE_DEALLOCATE_P(this%renorm_occ)

   do iat = 1, this%natoms
     do iorb = 1, this%norbs(iat)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,iat)%dorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,iat)%zorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,iat)%phase)
         call submesh_end(this%orbitals(iorb,iat)%sphere)
     end do
   end do 
 
   SAFE_DEALLOCATE_P(this%norbs)
   SAFE_DEALLOCATE_P(this%orbitals)

   POP_SUB(lda_u_end)
 end subroutine lda_u_end

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

 ! Interface for the X(compute_ACBNO_U) routines
 subroutine lda_u_update_U(this, st)
   type(lda_u_t),             intent(inout) :: this
   type(states_t),            intent(in)    :: st

   if(.not. this%apply.or..not. this%useACBN0 .or.this%freeze_u) return
   PUSH_SUB(lda_u_update_U)

   if (states_are_real(st)) then
     call dcompute_ACBNO_U(this, st) 
   else
     call zcompute_ACBNO_U(this, st)
   end if

   POP_SUB(lda_u_update_U)
 end subroutine lda_u_update_U
  

 !> Build the phase correction to the global phase for all orbitals
 subroutine lda_u_build_phase_correction(this, sb, std, vec_pot, vec_pot_var)
   type(lda_u_t),                 intent(inout) :: this
   type(simul_box_t),             intent(in)    :: sb 
   type(states_dim_t),            intent(in)    :: std
   FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
   FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)

   integer :: iat, ispin, iorb
 
   PUSH_SUB(lda_u_build_phase_correction)

   do iat = 1, this%natoms
     do iorb = 1, this%norbs(iat)
       call  orbital_update_phase_correction(this%orbitals(iorb,iat), sb, std, vec_pot, vec_pot_var)
     end do
   end do
  
   POP_SUB(lda_u_build_phase_correction)

 end subroutine lda_u_build_phase_correction

  !TODO: merge with the routine of projector.F90
  !> Build the phase correction to the global phase in case the orbital crosses the border of the simulaton box
  subroutine orbital_update_phase_correction(this, sb, std, vec_pot, vec_pot_var)
    type(orbital_t),               intent(inout) :: this
    type(simul_box_t),             intent(in)    :: sb
    type(states_dim_t),            intent(in)    :: std
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)

    integer :: ns, iq, is, ikpoint
    FLOAT   :: kr, kpoint(1:MAX_DIM)
    integer :: ndim

    PUSH_SUB(orbital_update_phase_correction)

    ns = this%sphere%np
    ndim = sb%dim

    do iq = std%kpt%start, std%kpt%end
      ikpoint = states_dim_get_kpoint_index(std, iq)

      ! if this fails, it probably means that sb is not compatible with std
      ASSERT(ikpoint <= kpoints_number(sb%kpoints))

      kpoint = M_ZERO
      kpoint(1:ndim) = kpoints_get_point(sb%kpoints, ikpoint)

      do is = 1, ns
        ! this is only the correction to the global phase, that can
        ! appear if the sphere crossed the boundary of the cell.

        kr = sum(kpoint(1:ndim)*(this%sphere%x(is, 1:ndim) - this%sphere%mesh%x(this%sphere%map(is), 1:ndim)))

        if(present(vec_pot)) then
          if(allocated(vec_pot)) kr = kr + &
            sum(vec_pot(1:ndim)*(this%sphere%x(is, 1:ndim)- this%sphere%mesh%x(this%sphere%map(is), 1:ndim)))
        end if

        if(present(vec_pot_var)) then
          if(allocated(vec_pot_var)) kr = kr + sum(vec_pot_var(1:ndim, this%sphere%map(is))*this%sphere%x(is, 1:ndim))
        end if

        this%phase(is, iq) = exp(-M_zI*kr)
      end do

    end do

    POP_SUB(orbital_update_phase_correction)

  end subroutine orbital_update_phase_correction

 subroutine lda_u_freeze_occ(this) 
   type(lda_u_t),     intent(inout) :: this

   this%freeze_occ = .true.
 end subroutine lda_u_freeze_occ

 subroutine lda_u_freeze_u(this)            
   type(lda_u_t),     intent(inout) :: this

   this%freeze_u = .true.
 end subroutine lda_u_freeze_u

#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
