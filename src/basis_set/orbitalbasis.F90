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

module orbitalbasis_oct_m
  use atomic_orbital_oct_m
  use distributed_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use mesh_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use orbitalset_utils_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use submesh_oct_m
 
  implicit none

  private

  public ::                            &
       orbitalbasis_t,                 &
       orbitalbasis_nullify,           &
       orbitalbasis_init,              &
       orbitalbasis_end,               &
       dorbitalbasis_build,            &
       zorbitalbasis_build,            &
       dorbitalbasis_build_empty,      &
       zorbitalbasis_build_empty

  type orbitalbasis_t
    private
    type(orbitalset_t), pointer, public :: orbsets(:) !> All the orbital sets of the system

    integer,            public :: norbsets           !> Number of orbital sets
    integer,            public :: maxnorbs           !> Maximal number of orbitals for all the atoms
    integer,            public :: max_np             !> Max. number of points in all orbitals submesh spheres
    integer,            public :: size               !> Size of the full basis
    integer, pointer,   public :: global2os(:,:)     !> Mapping functions
    integer, pointer           :: os2global(:,:)

    integer(8)                 :: truncation         !> Truncation method for the orbitals
    FLOAT                      :: threshold          !> Threshold for orbital truncation

    logical                    :: normalize          !> Do we normalize the orbitals
    logical,            public :: submeshforperiodic !> Do we use or not submeshes for the orbitals
    logical,            public :: orthogonalization  !> Orthogonalization of the basis

    character(len=256), public :: debugdir !> For debug
  end type orbitalbasis_t

contains

 subroutine orbitalbasis_nullify(this)
  type(orbitalbasis_t),    intent(inout) :: this

  PUSH_SUB(orbitalbasis_nullify)

  this%normalize = .true.
  this%submeshforperiodic = .false.
  this%orthogonalization = .false.

  nullify(this%orbsets)
  nullify(this%global2os)
  nullify(this%os2global)

  this%norbsets = 0
  this%maxnorbs = 0
  this%max_np = 0
  this%size = 0

  this%threshold = CNST(0.01)

  POP_SUB(orbitalbasis_nullify)

 end subroutine orbitalbasis_nullify

 subroutine orbitalbasis_init(this, namespace)
   type(orbitalbasis_t),    intent(inout) :: this
   type(namespace_t),       intent(in)    :: namespace

  PUSH_SUB(orbitalbasis_init)

  call orbitalbasis_nullify(this)

  !%Variable AOTruncation
  !%Type flag
  !%Default ao_full
  !%Section Hamiltonian::DFT+U
  !%Description
  !% This option determines how Octopus will truncate the orbitals used for LDA+U.
  !% Except for the full method, the other options are only there to get a quick idea.
  !%Option ao_full bit(0)
  !% The full size of the orbitals used. The radius is controled by variable OrbitalThreshold_LDAU
  !%Option ao_box bit(1)
  !% The radius of the orbitals are restricted to the size of the simulation box. 
  !% This reduces the number of points used to discretize the orbitals.
  !%Option ao_nlradius bit(2)
  !% The radius of the orbitals are restricted to the radius of the non-local part of the pseudopotential 
  !% of the corresponding atom.
  !%End
  call parse_variable(namespace, 'AOTruncation', OPTION__AOTRUNCATION__AO_FULL, this%truncation)
  call messages_print_var_option(stdout, 'AOTruncation', this%truncation)

  !%Variable AOThreshold
  !%Type float
  !%Default 0.01
  !%Section Hamiltonian::DFT+U
  !%Description
  !% Determines the threshold used to compute the radius of the atomic orbitals for LDA+U and for Wannier90.
  !% This radius is computed by making sure that the 
  !% absolute value of the radial part of the atomic orbital is below the specified threshold.
  !% This value should be converged to be sure that results do not depend on this value. 
  !% However increasing this value increases the number of grid points covered by the orbitals and directly affect performances.
  !%End
  call parse_variable(namespace, 'AOThreshold', CNST(0.01), this%threshold)
  if(this%threshold <= M_ZERO) call messages_input_error('AOThreshold')
  call messages_print_var_value(stdout, 'AOThreshold', this%threshold)

  !%Variable AONormalize
  !%Type logical
  !%Default yes
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will normalize the atomic orbitals
  !%End
  call parse_variable(namespace, 'AONormalize', .true., this%normalize)
  call messages_print_var_value(stdout, 'AONormalize', this%normalize)

  !%Variable AOSubmeshForPeriodic
  !%Type logical
  !%Default no
  !%Section Hamiltonian::DFT+U
  !%Description
  !% If set to yes, Octopus will use submeshes to internally store the orbitals with
  !% their phase instead of storing them on the mesh. This is usually slower for small
  !% periodic systems, but becomes advantageous for large supercells.
  !% At the moment this option is not compatible with Loewdin orthogonalization
  !%End
  call parse_variable(namespace, 'AOSubmeshForPeriodic', .false., this%submeshforperiodic)
  call messages_print_var_value(stdout, 'AOSubmeshForPeriodic', this%submeshforperiodic)

  !%Variable AOLoewdin
  !%Type logical
  !%Default no
  !%Section Hamiltonian::DFT+U
  !%Description
  !% This option determines if the atomic orbital basis is orthogonalized or not.
  !% This is done for using the Loewdin orthogonalization scheme.
  !% The default is set to no for the moment as this option is
  !% not yet implemented for isolated systems, and seems to lead to important egg-box effect
  !%End
  call parse_variable(namespace, 'AOLoewdin', .false., this%orthogonalization)
  call messages_print_var_value(stdout, 'AOLoewdin', this%orthogonalization)
  if(this%orthogonalization) call messages_experimental("AOLoewdin")

  if(this%orthogonalization .and. this%submeshforperiodic) &
    call messages_not_implemented("AOLoewdin=yes with AOSubmeshForPeriodic=yes.") 

  if(debug%info) then
    write(this%debugdir, '(a)') 'debug/ao_basis'
    call io_mkdir(this%debugdir, namespace)
  end if

  POP_SUB(orbitalbasis_init)
 end subroutine orbitalbasis_init


 subroutine orbitalbasis_end(this)
   implicit none
   type(orbitalbasis_t), intent(inout) :: this

   integer :: ios

   PUSH_SUB(orbitalbasis_end)  

   do ios = 1, this%norbsets
     call orbitalset_end(this%orbsets(ios))
   end do
   SAFE_DEALLOCATE_P(this%orbsets)
  
   SAFE_DEALLOCATE_P(this%global2os)
   SAFE_DEALLOCATE_P(this%os2global)

   POP_SUB(orbitalbasis_end)
 end subroutine orbitalbasis_end

#include "undef.F90"
#include "real.F90"
#include "orbitalbasis_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbitalbasis_inc.F90"

end module orbitalbasis_oct_m
