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
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use lalg_basic_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use orbitalset_utils_oct_m
  use parser_oct_m
  use profiling_oct_m
  use species_oct_m
  use submesh_oct_m
 
  implicit none

  private

  public ::                            &
       orbitalbasis_t,                 &
       orbitalbasis_init,              &
       orbitalbasis_end,               &
       dorbitalbasis_build,            &
       zorbitalbasis_build,            &
       dorbitalbasis_build_empty,      &
       zorbitalbasis_build_empty

  type orbitalbasis_t
    private
    type(orbitalset_t), allocatable, public :: orbsets(:) !> All the orbital sets of the system

    integer,            public :: norbsets = 0       !> Number of orbital sets
    integer,            public :: maxnorbs = 0       !> Maximal number of orbitals for all the atoms
    integer,            public :: max_np   = 0       !> Max. number of points in all orbitals submesh spheres
    integer,            public :: size     = 0       !> Size of the full basis
    integer, allocatable, public :: global2os(:,:)     !> Mapping functions
    integer, allocatable         :: os2global(:,:)

    integer(8)                 :: truncation             !> Truncation method for the orbitals
    FLOAT                      :: threshold = CNST(0.01) !> Threshold for orbital truncation

    logical                    :: normalize = .true.  !> Do we normalize the orbitals
    logical,            public :: submesh   = .false. !> Do we use or not submeshes for the orbitals
    logical,            public :: orthogonalization = .false. !> Orthogonalization of the basis

    character(len=256), public :: debugdir !> For debug
  end type orbitalbasis_t

contains

subroutine orbitalbasis_init(this, namespace, periodic_dim)
  type(orbitalbasis_t),    intent(inout) :: this
  type(namespace_t),       intent(in)    :: namespace
  integer,                 intent(in)    :: periodic_dim

  PUSH_SUB(orbitalbasis_init)

  !%Variable AOTruncation
  !%Type flag
  !%Default ao_full
  !%Section Atomic Orbitals
  !%Description
  !% This option determines how Octopus will truncate the orbitals used for LDA+U.
  !% Except for the full method, the other options are only there to get a quick idea.
  !%Option ao_full bit(0)
  !% The full size of the orbitals used. The radius is controled by variable AOThreshold.
  !%Option ao_box bit(1)
  !% The radius of the orbitals are restricted to the size of the simulation box. 
  !% This reduces the number of points used to discretize the orbitals.
  !% This is mostly a debug option, and one should be aware that changing the size of the simulation box
  !% will affect the result of the calculation. It is recommended to use ao_nlradius instead.
  !%Option ao_nlradius bit(2)
  !% The radius of the orbitals are restricted to the radius of the non-local part of the pseudopotential 
  !% of the corresponding atom.
  !%End
  call parse_variable(namespace, 'AOTruncation', OPTION__AOTRUNCATION__AO_FULL, this%truncation)
  call messages_print_var_option(stdout, 'AOTruncation', this%truncation)

  !%Variable AOThreshold
  !%Type float
  !%Default 0.01
  !%Section Atomic Orbitals
  !%Description
  !% Determines the threshold used to compute the radius of the atomic orbitals for LDA+U and for Wannier90.
  !% This radius is computed by making sure that the 
  !% absolute value of the radial part of the atomic orbital is below the specified threshold.
  !% This value should be converged to be sure that results do not depend on this value. 
  !% However increasing this value increases the number of grid points covered by the orbitals and directly affect performances.
  !%End
  call parse_variable(namespace, 'AOThreshold', CNST(0.01), this%threshold)
  if(this%threshold <= M_ZERO) call messages_input_error(namespace, 'AOThreshold')
  call messages_print_var_value(stdout, 'AOThreshold', this%threshold)

  !%Variable AONormalize
  !%Type logical
  !%Default yes
  !%Section Atomic Orbitals
  !%Description
  !% If set to yes, Octopus will normalize the atomic orbitals individually.
  !% This variable is ignored is <tt>AOLoewdin<\tt> is set to yes.
  !%End
  call parse_variable(namespace, 'AONormalize', .true., this%normalize)
  call messages_print_var_value(stdout, 'AONormalize', this%normalize)

  !%Variable AOSubmesh
  !%Type logical
  !%Section Atomic Orbitals
  !%Description
  !% If set to yes, Octopus will use submeshes to internally store the orbitals with
  !% their phase instead of storing them on the mesh. This is usually slower for small
  !% periodic systems, but becomes advantageous for large supercells.
  !% Submeshes are not compatible with Loewdin orthogonalization.
  !% For periodic systems, the default is set to ye, whereas it is set to no for isolated systems.
  !%End
  if(periodic_dim > 0) then
    call parse_variable(namespace, 'AOSubmesh', .false., this%submesh)
  else
    call parse_variable(namespace, 'AOSubmesh', .true., this%submesh)
  end if
  call messages_print_var_value(stdout, 'AOSubmesh', this%submesh)

  !%Variable AOLoewdin
  !%Type logical
  !%Default no
  !%Section Atomic Orbitals
  !%Description
  !% This option determines if the atomic orbital basis is orthogonalized or not.
  !% This is done for using the Loewdin orthogonalization scheme.
  !% The default is set to no for the moment as this option is
  !% not yet implemented for isolated systems, and seems to lead to important egg-box effect
  !%End
  call parse_variable(namespace, 'AOLoewdin', .false., this%orthogonalization)
  call messages_print_var_value(stdout, 'AOLoewdin', this%orthogonalization)
  if(this%orthogonalization) call messages_experimental("AOLoewdin")
  if(this%orthogonalization) this%normalize = .false.

  if(this%orthogonalization .and. this%submesh) &
    call messages_not_implemented("AOLoewdin=yes with AOSubmesh=yes.") 

  if(debug%info) then
    write(this%debugdir, '(a)') 'debug/ao_basis'
    call io_mkdir(this%debugdir, namespace)
  end if

  POP_SUB(orbitalbasis_init)
 end subroutine orbitalbasis_init


 subroutine orbitalbasis_end(this)
   type(orbitalbasis_t), intent(inout) :: this

   integer :: ios

   PUSH_SUB(orbitalbasis_end)  

   do ios = 1, this%norbsets
     call orbitalset_end(this%orbsets(ios))
   end do
   SAFE_DEALLOCATE_A(this%orbsets)
  
   SAFE_DEALLOCATE_A(this%global2os)
   SAFE_DEALLOCATE_A(this%os2global)

   POP_SUB(orbitalbasis_end)
 end subroutine orbitalbasis_end

#include "undef.F90"
#include "real.F90"
#include "orbitalbasis_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbitalbasis_inc.F90"

end module orbitalbasis_oct_m
