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
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use simul_box_oct_m
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
       dupdate_occ_matrices, &
       zupdate_occ_matrices, &
       lda_u_end,            &
       lda_u_build_phase_correction,   &
       lda_u_write_occupation_matrices

  type orbital_t
    type(submesh_t)     :: sphere             !The submesh of the orbital
    FLOAT, pointer      :: dorbital_sphere(:) !The orbital, if real, on the submesh
    CMPLX, pointer      :: zorbital_sphere(:) !The orbital, if complex, on the submesh
    FLOAT, pointer      :: dorbital_mesh(:)   !The orbital, if real, on the full mesh
    CMPLX, pointer      :: zorbital_mesh(:)   !The orbital, if complex, on the full mesh
    CMPLX, pointer      :: phase(:,:)         !Correction to the global phase if the sphere cross the border of the box
  end type orbital_t


  type lda_u_t
    logical                  :: apply
    FLOAT, pointer           :: dn(:,:,:,:), dV(:,:,:,:)
    CMPLX, pointer           :: zn(:,:,:,:), zV(:,:,:,:)
    type(orbital_t), pointer :: orbitals(:,:,:)    
    FLOAT, pointer           :: Ueff(:) !> The effective U of the simplified rotational invariant form

    integer             :: natoms          !> Number of atoms (copied from geometry_t)
    integer             :: nspins
    integer, pointer    :: norbs(:)        !> Number of orbitals
    integer             :: maxnorbs        !> Maximal number of orbitals for all the atoms
    integer             :: projection      !> The method used to perform the projection
  
    logical             :: truncate        !> Do we truncate the orbitals to the radius 
                                           !> of the NL part of the pseudo
    logical             :: selfconsistentU !> Do we compute U self-consistenly
  end type lda_u_t

contains

 subroutine lda_u_init(this, gr, geo, st)
  implicit none

  type(lda_u_t),             intent(inout) :: this
  type(grid_t),              intent(in)    :: gr
  type(geometry_t), target,  intent(in)    :: geo
  type(states_t),            intent(in)    :: st

  integer :: maxorbs, iat, ispin, iorb

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

  !%Variable UseACBNOFunctional
  !%Type logical
  !%Default no
  !%Section Hamiltonian::LDA+U
  !%Description
  !% If set to yes, Octopus will determine the effective U term using the 
  !% ACBN0 functional as defined in PRX 5, 011006 (2015) 
  !%End
  call parse_variable('UseACBNOFunctional', .false., this%selfconsistentU)

  nullify(this%dn)
  nullify(this%zn)
  nullify(this%dV)
  nullify(this%zV)
  nullify(this%Ueff)

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
   SAFE_DEALLOCATE_P(this%Ueff)
 
   do iat = 1, this%natoms
     do ispin = 1, this%nspins 
       do iorb = 1, this%norbs(iat)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%dorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%zorbital_sphere)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%dorbital_mesh)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%zorbital_mesh)
         SAFE_DEALLOCATE_P(this%orbitals(iorb,ispin,iat)%phase)
         call submesh_end(this%orbitals(iorb,ispin,iat)%sphere)
       end do
     end do
   end do
  
   SAFE_DEALLOCATE_P(this%norbs)
   SAFE_DEALLOCATE_P(this%orbitals)

   POP_SUB(lda_u_end)
 end subroutine lda_u_end

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
     do ispin = 1, this%nspins
       do iorb = 1, this%norbs(iat)
         call  orbital_update_phase_correction(this%orbitals(iorb,ispin,iat), sb, std, vec_pot, vec_pot_var)
       end do
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

 !> Prints the occupation matrices at the end of the scf calculation.
 subroutine lda_u_write_occupation_matrices(dir, this, geo, st)
   implicit none
   type(lda_u_t),     intent(inout) :: this
   character(len=*),  intent(in)    :: dir
   type(geometry_t),  intent(in)    :: geo
   type(states_t),    intent(in)    :: st

   integer :: iunit, ia, ispin, im, imp
   FLOAT :: hubbardl
 
   PUSH_SUB(lda_u_write_occupation_matrices)

   if(.not. mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
  
   call io_mkdir(dir)
   iunit = io_open(trim(dir) // "/occ_matrices", action='write')
   write(iunit,'(a)') ' Occupation matrices '

   do ia = 1, this%natoms
     hubbardl = species_hubbard_l(geo%atom(ia)%species)
     if( hubbardl .eq. M_ZERO ) cycle

     do ispin = 1,st%d%nspin 
        write(iunit,'(a, i4, a, i4)') ' Ion ', ia, ' spin ', ispin
        do im = 1, this%norbs(ia)
          write(iunit,'(1x)',advance='no') 

          if (states_are_real(st)) then
            do imp = 1, this%norbs(ia)-1
              write(iunit,'(f14.8)',advance='no') this%dn(im,imp,ispin,ia)  
            end do
            write(iunit,'(f14.8)') this%dn(im,this%norbs(ia),ispin,ia)
          else
            do imp = 1, this%norbs(ia)-1
              write(iunit,'(f14.8,f14.8)',advance='no') this%zn(im,imp,ispin,ia)
            end do
            write(iunit,'(f14.8,f14.8)') this%zn(im,this%norbs(ia),ispin,ia) 
          end if
        end do
     end do !ispin
   end do !iatom
   call io_close(iunit)

   end if

   POP_SUB(lda_u_write_occupation_matrices)
 end subroutine lda_u_write_occupation_matrices

#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m