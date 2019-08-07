!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module subspace_oct_m
  use accel_oct_m
  use accel_blas_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use blacs_oct_m
  use blacs_proc_grid_oct_m
  use comm_oct_m
  use derivatives_oct_m
#ifdef HAVE_ELPA
  use elpa
#endif
  use global_oct_m
  use hamiltonian_oct_m
  use hardware_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pblas_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_parallel_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none
  
  private

  public ::               &
    subspace_t,           &
    subspace_init,        &
    dsubspace_diag,       &
    zsubspace_diag

  type subspace_t
    private
    integer :: method
  end type subspace_t

  type(profile_t),     save    :: diagon_prof, hamiltonian_prof
  
contains

  subroutine subspace_init(this, namespace, st, no_sd)
    type(subspace_t),    intent(out) :: this
    type(namespace_t),   intent(in)  :: namespace
    type(states_elec_t), intent(in)  :: st
    logical,             intent(in)  :: no_sd

    integer :: default

    PUSH_SUB(subspace_init)

    if(no_sd) then

      this%method = OPTION__SUBSPACEDIAGONALIZATION__NONE

    else

      !%Variable SubspaceDiagonalization
      !%Type integer
      !%Default standard
      !%Section SCF::Eigensolver
      !%Description
      !% Selects the method to perform subspace diagonalization. The
      !% default is <tt>standard</tt>, unless states parallelization is used,
      !% when the default is <tt>scalapack</tt>.
      !%Option none 0
      !% No subspace diagonalization. WARNING: this will generally give incorrect results.
      !%Option standard 1
      !% The standard routine. Can be used with domain parallelization but not
      !% state parallelization.
      !%Option scalapack 3
      !% State-parallelized version using ScaLAPACK. (Requires that
      !% Octopus was compiled with ScaLAPACK support.)
      !%End

      default = OPTION__SUBSPACEDIAGONALIZATION__STANDARD

#ifdef HAVE_SCALAPACK
      if(st%parallel_in_states) default = OPTION__SUBSPACEDIAGONALIZATION__SCALAPACK
#endif

      call parse_variable(namespace, 'SubspaceDiagonalization', default, this%method)

      if(.not.varinfo_valid_option('SubspaceDiagonalization', this%method)) call messages_input_error('SubspaceDiagonalization')
    end if

    call messages_print_var_option(stdout, 'SubspaceDiagonalization', this%method)

    ! some checks for ingenious users
    if(this%method == OPTION__SUBSPACEDIAGONALIZATION__SCALAPACK) then
#ifndef HAVE_MPI
      message(1) = 'The scalapack subspace diagonalization can only be used in parallel.'
      call messages_fatal(1, only_root_writes = .true.)
#else
#ifndef HAVE_SCALAPACK
      message(1) = 'The scalapack subspace diagonalization requires scalapack.'
      call messages_fatal(1, only_root_writes = .true.)
#endif
      if(st%dom_st_mpi_grp%size == 1) then
        message(1) = 'The scalapack subspace diagonalization is designed to be used with domain or state parallelization.'
        call messages_warning(1)
      end if

      if(st%d%kpt%parallel) then
        message(1) = 'Currently the scalapack subspace diagonalization does not use k-point parallelization.'
        call messages_warning(1)
      end if
#endif
    end if

    POP_SUB(subspace_init)
  end subroutine subspace_init

#include "undef.F90"
#include "real.F90"
#include "subspace_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "subspace_inc.F90"

end module subspace_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
