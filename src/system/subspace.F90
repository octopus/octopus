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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: eigen.F90 4287 2008-06-15 22:20:10Z xavier $

#include "global.h"

module subspace_m
  use batch_m
  use batch_ops_m
  use blas_m
  use blacs_proc_grid_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use mesh_function_m
  use mesh_batch_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use pblas_m
  use preconditioners_m
  use profiling_m
  use scalapack_m
  use states_m
  use states_block_m
  use states_calc_m
  use types_m
  use varinfo_m

  implicit none
  
  private

  public ::             &
    subspace_t,         &
    subspace_init,      &
    subspace_end,       &
    dsubspace_diag,     &
    zsubspace_diag

  integer, parameter ::        &
    SD_STANDARD  = 1,          &
    SD_SCALAPACK = 3
  
  type subspace_t
    integer :: method
  end type subspace_t


contains

  subroutine subspace_init(this, st)
    type(subspace_t), intent(out) :: this
    type(states_t),   intent(in)  :: st

    integer :: default

    !%Variable SubspaceDiagonalization
    !%Type integer
    !%Default standard
    !%Section SCF::Eigensolver
    !%Description
    !% Selects the method to perform subspace diagonalization. The
    !% default is <tt>standard</tt>, unless states parallelization is used,
    !% when the default is <tt>scalapack</tt>.
    !%Option standard 1
    !% The standard routine. Can be used with domain parallelization but not
    !% state parallelization.
    !%Option scalapack 3
    !% State-parallelized version using ScaLAPACK. (Requires that
    !% Octopus was compiled with ScaLAPACK support.)
    !%End

    default = SD_STANDARD

    if(st%parallel_in_states) then
#ifdef HAVE_SCALAPACK
      default = SD_SCALAPACK
#else
      message(1) = 'Parallelization in states of the ground state requires scalapack.'
      call messages_fatal(1)
#endif
    end if

    call parse_integer(datasets_check('SubspaceDiagonalization'), default, this%method)

    if(.not.varinfo_valid_option('SubspaceDiagonalization', this%method)) call input_error('SubspaceDiagonalization')

    call messages_print_var_option(stdout, 'SubspaceDiagonalization', this%method)

    ! some checks for ingenious users
    if(this%method == SD_SCALAPACK) then
#ifndef HAVE_MPI
      message(1) = 'The scalapack subspace diagonalization can only be used in parallel.'
      call messages_fatal(1)
#else
#ifndef HAVE_SCALAPACK
      message(1) = 'The scalapack subspace diagonalization requires scalapack.'
      call messages_fatal(1)
#endif
      if(st%dom_st_mpi_grp%size == 1) then
        message(1) = 'The scalapack subspace diagonalization is designed to be used with domain or state parallelization.'
        call messages_warning(1)
      end if

      if(st%d%kpt%parallel) then
        message(1) = 'Currently the scalapack subspace diagonalization cannot work with subspace diagonalization.'
        call messages_warning(1)
      end if
#endif
    end if

#ifdef HAVE_MPI
    if(this%method == SD_STANDARD .and. st%parallel_in_states) then
      message(1) = 'The standard subspace diagonalization cannot work with state parallelization.'
      call messages_fatal(1)
    end if
#endif
    
  end subroutine subspace_init

  ! ---------------------------------------------------------

  subroutine subspace_end(this)
    type(subspace_t), intent(inout) :: this
    
  end subroutine subspace_end

#include "undef.F90"
#include "real.F90"
#include "subspace_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "subspace_inc.F90"

end module subspace_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
