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
!! $Id$

#include "global.h"

module unocc_m
  use datasets_m
  use density_m
  use eigensolver_m
  use global_m
  use output_m
  use hamiltonian_m
  use io_m
  use lcao_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_io_m
  use states_dim_m
  use system_m
  use v_ks_m
  use xc_m

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
    integer :: iunit, ierr, iter, ierr_rho, ik
    logical :: converged, forced_finish, showoccstates, is_orbital_dependent
    integer :: max_iter, nst_calculated, showstart
    integer :: n_filled, n_partially_filled, n_half_filled
    integer, allocatable :: states_read(:, :), occ_states(:)
    character(len=50) :: str

    PUSH_SUB(unocc_run)

    !%Variable MaximumIter
    !%Type integer
    !%Default 50
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% Maximum number of eigensolver iterations. The code will stop even if convergence
    !% has not been achieved. -1 means unlimited. 0 means just do LCAO or read from
    !% restart, and stop.
    !%End
    call parse_integer(datasets_check('MaximumIter'), 50, max_iter)
    call messages_obsolete_variable('UnoccMaximumIter', 'MaximumIter')
    if(max_iter < 0) max_iter = huge(max_iter)

    !%Variable UnoccShowOccStates
    !%Type logical
    !%Default false
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% If true, the convergence for the occupied states will be shown too in the output.
    !% This is useful for testing, or if the occupied states fail to converge.
    !%End
    call parse_logical(datasets_check('UnoccShowOccStates'), .false., showoccstates)

    SAFE_ALLOCATE(occ_states(1:sys%st%d%nik))
    do ik = 1, sys%st%d%nik
      call occupied_states(sys%st, ik, n_filled, n_partially_filled, n_half_filled)
      occ_states(ik) = n_filled + n_partially_filled + n_half_filled
    enddo

    call init_(sys%gr%mesh, sys%st)
    converged = .false.

    call restart_read_rho(trim(restart_dir)//GS_DIR, sys%st, sys%gr, ierr_rho)

    SAFE_ALLOCATE(states_read(1:sys%st%d%dim, 1:sys%st%d%nik))

    call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, ierr, number_read = states_read)
      
    if(ierr_rho /= 0) then
      message(1) = "Building density from wavefunctions."
      call messages_info(1)
    endif

    is_orbital_dependent = (sys%ks%theory_level == HARTREE .or. sys%ks%theory_level == HARTREE_FOCK .or. &
      xc_is_orbital_dependent(sys%ks%xc))

    if(is_orbital_dependent) then
      message(1) = "CalculationMode = unocc is not well-defined for orbital-dependent functionals,"
      message(2) = "since merely freezing the density does not guarantee consistency with the gs run."
      message(3) = "Density and occupied orbitals from the same k-points are required to restart."
      call messages_warning(3)
    endif

    if(ierr_rho /= 0 .or. is_orbital_dependent) then
      ! the array needs to hold all states and k-points, but each node is responsible for checking its own states
      do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
        if(any(states_read(1:sys%st%d%dim, ik) < occ_states(ik))) then
          message(1) = "Not all the occupied KS orbitals could be read from '"//trim(restart_dir)//GS_DIR//"'"
          message(2) = "Please run a ground-state calculation first!"
          call messages_fatal(2)
        end if
      enddo

      if(.not. hm%cmplxscl%space) then
        call density_calc(sys%st, sys%gr, sys%st%rho)
      else
        call density_calc(sys%st, sys%gr, sys%st%zrho%Re, sys%st%zrho%Im)
      endif
    end if

    if(ierr /= 0) then
      message(1) = "Info: Could not load all wavefunctions from '"//trim(restart_dir)//GS_DIR//"'"
      call messages_info(1)
    end if

    if(sys%st%d%ispin == SPINORS) then
      message(1) = "Try gs with ExtraStates instead of unocc mode for spinors."
      call messages_warning(1)
      call messages_experimental("unocc for spinors")
    endif

    if (states_are_real(sys%st)) then
      message(1) = 'Info: Using real wavefunctions.'
    else
      message(1) = 'Info: Using complex wavefunctions.'
    end if
    call messages_info(1)

    if(fromScratch .or. ierr /= 0) then
      if(ierr > 0 .and. .not. fromScratch) then
        nst_calculated = minval(states_read)
      else
        ! if not all occupied states read, must recalculate
        nst_calculated = min(maxval(occ_states), minval(states_read))
      end if
      call lcao_run(sys, hm, st_start = nst_calculated + 1)
      showstart = nst_calculated + 1
    else
      call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval = .false.)
      showstart = minval(occ_states(:)) + 1
    end if
    
    SAFE_DEALLOCATE_A(states_read)

    if(showoccstates) showstart = 1

    message(1) = "Info: Starting calculation of unoccupied states."
    call messages_info(1)

    ! reset this variable, so that the eigensolver passes through all states
    eigens%converged(:) = 0

    do iter = 1, max_iter
      call eigensolver_run(eigens, sys%gr, sys%st, hm, 1, converged)

      write(str, '(a,i5)') 'Unoccupied states iteration #', iter
      call messages_print_stress(stdout, trim(str))
      call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb, eigens%diff, st_start = showstart)
      call messages_print_stress(stdout)

      forced_finish = clean_stop(sys%mc%master_comm)
      
      ! write restart information.
      if(converged .or. (modulo(iter, sys%outp%restart_write_interval) == 0) .or. iter == max_iter .or. forced_finish) then
        call restart_write(trim(tmpdir)//GS_DIR, sys%st, sys%gr, ierr, iter=iter)
        if(ierr /= 0) then
          message(1) = 'Unsuccessful write of "'//trim(tmpdir)//GS_DIR//'"'
          call messages_fatal(1)
        end if
      end if

      if(converged .or. forced_finish) exit
    end do

    if(any(eigens%converged(:) < occ_states(:))) then
      write(message(1),'(a)') 'Some of the occupied states are not fully converged!'
      call messages_warning(1)
    endif

    SAFE_DEALLOCATE_A(occ_states)

    if(.not. converged) then
      write(message(1),'(a)') 'Some of the unoccupied states are not fully converged!'
      call messages_warning(1)
    endif

    ! write output file
    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(STATIC_DIR)
      iunit = io_open(STATIC_DIR//'/eigenvalues', action='write')
      
      if(converged) then
        write(iunit,'(a)') 'All states converged.'
      else
        write(iunit,'(a)') 'Some of the states are not fully converged!'
      end if
      write(iunit,'(a, e17.6)') 'Criterion = ', eigens%tolerance
      write(iunit,'(1x)')
      call states_write_eigenvalues(iunit, sys%st%nst, sys%st, sys%gr%sb, eigens%diff)
      call io_close(iunit)
    end if

    if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik > sys%st%d%nspin) then
      call states_write_bands(STATIC_DIR, sys%st%nst, sys%st, sys%gr%sb)
      call states_write_fermi_energy(STATIC_DIR, sys%st, sys%gr%mesh, sys%gr%sb)
    end if

    call output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, sys%ks%xc, STATIC_DIR)

    call end_()
    POP_SUB(unocc_run)

  contains

    ! ---------------------------------------------------------
    subroutine init_(mesh, st)
      type(mesh_t),   intent(in)    :: mesh
      type(states_t), intent(inout) :: st

      PUSH_SUB(unocc_run.init_)

      call messages_obsolete_variable("NumberUnoccStates", "ExtraStates")

      if(st%d%ispin == SPINORS) then
        SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
        st%spin = M_ZERO
      end if
      call states_allocate_wfns(st, mesh)

      ! now the eigensolver stuff
      call eigensolver_init(eigens, sys%gr, st)

      POP_SUB(unocc_run.init_)
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(unocc_run.end_)

      call eigensolver_end(eigens)

      POP_SUB(unocc_run.end_)
    end subroutine end_

  end subroutine unocc_run


end module unocc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
