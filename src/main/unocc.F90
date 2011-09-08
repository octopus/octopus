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
!! $Id$

#include "global.h"

module unocc_m
  use datasets_m
  use density_m
  use derivatives_m
  use eigensolver_m
  use energy_m
  use epot_m
  use global_m
  use grid_m
  use geometry_m
  use output_m
  use hamiltonian_m
  use io_m
  use lcao_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
  use profiling_m
  use projector_m
  use restart_m
  use simul_box_m
  use states_m
  use states_io_m
  use states_calc_m
  use states_dim_m
  use system_m
  use v_ks_m
  use varinfo_m

  implicit none

  private
  public :: &
    unocc_run


contains

  ! ---------------------------------------------------------
  subroutine unocc_run(sys, hm, fromscratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromscratch

    type(eigensolver_t) :: eigens
    integer :: iunit, ierr, occupied_states, total_states, iter
    logical :: converged
    integer :: max_iter, nst_calculated
    integer, allocatable :: states_read(:, :)

    PUSH_SUB(unocc_run)

    !%Variable UnoccMaximumIter
    !%Type integer
    !%Default 50
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% Maximum number of eigensolver iterations. The code will stop even if convergence
    !% has not been achieved. -1 means unlimited.
    !%End
    call parse_integer(datasets_check('UnoccMaximumIter'), 50, max_iter)
    if(max_iter < 0) max_iter = huge(max_iter)
    call messages_obsolete_variable('MaximumIter', 'UnoccMaximumIter')

    occupied_states = sys%st%nst
    call init_(sys%gr%mesh, sys%st)
    total_states = sys%st%nst

    ASSERT(total_states >= occupied_states)

    SAFE_ALLOCATE(states_read(1:sys%st%d%dim, 1:sys%st%d%kpt%nlocal))

    call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, sys%geo, ierr, number_read = states_read)

    if(any(states_read < occupied_states)) then
      message(1) = "Not all the occupied KS orbitals could be read from '"//trim(restart_dir)//GS_DIR//"'"
      message(2) = "Please run a ground-state calculation first!"
      call messages_fatal(2)
    end if

    if(ierr .ne. 0) then
      message(1) = "Info:  Could not load all wavefunctions from '"//trim(restart_dir)//GS_DIR//"'"
      call messages_info(1)
    end if

    call density_calc(sys%st, sys%gr, sys%st%rho)

    if(fromScratch .or. ierr /= 0) then
      if(ierr > 0 .and. .not. fromScratch) then
        nst_calculated = minval(states_read)
      else
        nst_calculated = occupied_states
      end if
      call lcao_run(sys, hm, st_start = nst_calculated + 1)
    end if
    
    SAFE_DEALLOCATE_A(states_read)

    message(1) = "Info:  Starting calculation of unoccupied states."
    call messages_info(1)

    ! reset this variable, so that the eigensolver passes through all states
    eigens%converged(:) = 0

    do iter = 1, max_iter
      write(message(1), '(a,i3)') "Info: Unoccupied states iteration ", iter
      call messages_info(1)
      call eigensolver_run(eigens, sys%gr, sys%st, hm, 1, converged, verbose = .true.)

      if(converged .or. clean_stop()) exit
    end do

    ! write restart information.
    call restart_write (trim(tmpdir)//GS_DIR, sys%st, sys%gr, ierr)
    if(ierr .ne. 0) then
      message(1) = 'Unsuccessful write of "'//trim(tmpdir)//GS_DIR//'"'
      call messages_fatal(1)
    end if

    ! write output file
    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(STATIC_DIR)
      iunit = io_open(STATIC_DIR//'/eigenvalues', action='write')
      
      if(converged) then
        write(iunit,'(a)') 'All unoccupied states converged.'
      else
        write(iunit,'(a)') 'Some of the unoccupied states are not fully converged!'
      end if
      write(iunit,'(a, e17.6)') 'Criterion = ', eigens%final_tol
      write(iunit,'(1x)')
      call states_write_eigenvalues(iunit, sys%st%nst, sys%st, sys%gr%sb, eigens%diff)
      call io_close(iunit)
    end if

    if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik > sys%st%d%nspin) then
      call states_write_bands(STATIC_DIR, sys%st%nst, sys%st, sys%gr%sb)
      call states_write_fermi_energy(STATIC_DIR, sys%st, sys%gr%mesh, sys%gr%sb)
    end if

    ! output wavefunctions
    call output_states(sys%st, sys%gr, sys%geo, STATIC_DIR, sys%outp)

    call end_()
    POP_SUB(unocc_run)

  contains

    ! ---------------------------------------------------------
    subroutine init_(mesh, st)
      type(mesh_t),   intent(in)    :: mesh
      type(states_t), intent(inout) :: st

      integer :: nus

      PUSH_SUB(unocc_run.init_)

      !%Variable NumberUnoccStates
      !%Type integer
      !%Default 5
      !%Section Calculation Modes::Unoccupied States
      !%Description
      !% How many unoccupied states to compute.
      !%End
      call parse_integer(datasets_check('NumberUnoccStates'), 5, nus)
      if(nus <= 0) then
        message(1) = "Input: NumberUnoccStates must be > 0"
        call messages_fatal(1)
      end if

      ! fix states: THIS IS NOT OK
      st%nst    = st%nst + nus
      st%st_end = st%nst

      SAFE_DEALLOCATE_P(st%eigenval)
      SAFE_DEALLOCATE_P(st%occ)
      call states_allocate_wfns(st, mesh)
      SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
      SAFE_ALLOCATE(st%occ(1:st%nst, 1:st%d%nik))
      if(st%d%ispin == SPINORS) then
        SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
        st%spin = M_ZERO
      end if
      st%eigenval = huge(st%eigenval)
      st%occ      = M_ZERO

      ! now the eigensolver stuff
      call eigensolver_init(eigens, sys%gr, st)

      ! Having initial and final tolerance does not make sense in this case:
      eigens%init_tol       = eigens%final_tol
      eigens%final_tol_iter = 2
      eigens%converged(1:st%d%nik) = st%nst - nus

      POP_SUB(unocc_run.init_)
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(unocc_run.end_)

      call states_deallocate_wfns(sys%st)
      call eigensolver_end(eigens)

      POP_SUB(unocc_run.end_)
    end subroutine end_

  end subroutine unocc_run


end module unocc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
