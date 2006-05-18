!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module ground_state_m
  use global_m
  use messages_m
  use varinfo_m
  use datasets_m
  use lib_oct_parser_m
  use system_m
  use hamiltonian_m
  use v_ks_m
  use lcao_m
  use states_m
  use restart_m
  use scf_m
  use grid_m
  use simul_box_m
  use mesh_m
  use mpi_m
  use mpi_debug_m

  implicit none

  integer, parameter ::     &
    LCAO_START_NONE    = 0, &
    LCAO_START_STATES  = 1, &
    LCAO_START_FULL    = 2

  private
  public ::                 &
    ground_state_run


contains

  ! ---------------------------------------------------------
  subroutine ground_state_run(sys, h, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    logical,             intent(inout) :: fromScratch

    integer :: lcao_start, lcao_start_default = LCAO_START_FULL
    type(lcao_t) :: lcao_data
    type(scf_t) :: scfv
    integer     :: ierr

    call push_sub('gs.ground_state_run')

    call states_allocate_wfns(sys%st, sys%gr%m)

    if(.not.fromScratch) then
      ! load wave-functions
      call restart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
      if(ierr.ne.0) then
        message(1) = "Could not load wave-functions from '"//trim(tmpdir)//"restart_gs'"
        message(2) = "Starting from scratch!"
        call write_warning(1)
        fromScratch = .true.
      end if
    end if

    ! set barrier before the first communication takes place
    ! this ensures proper debug timing of MPI calls
#if defined(HAVE_MPI)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif

    if(fromScratch) then
      ! Randomly generate the initial wave-functions
      call states_generate_random(sys%st, sys%gr%m)

      ! We do not compute the density from the random wave-functions. 
      ! Instead, we try to get a better guess for the density
      call system_guess_density(sys%gr%m, sys%gr%sb, sys%gr%geo, sys%st%qtot, sys%st%d%nspin, &
           sys%st%d%spin_channels, sys%st%rho)

      ! setup Hamiltonian (we do not call system_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call v_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
      call states_fermi(sys%st, sys%gr%m)                                ! occupations
      call hamiltonian_energy(h, sys%gr, sys%st, -1)             ! total energy

      ! The initial LCAO calculation is done by default if we have pseudopotentials.
      ! Otherwise, it is not the default value and has to be enforced in the input file.
      if(sys%gr%geo%only_user_def) lcao_start_default = LCAO_START_NONE

      !%Variable LCAOStart
      !%Type integer
      !%Default lcao_states
      !%Section SCF
      !%Description
      !% Before starting a SCF calculation, <tt>octopus</tt> can perform
      !% a LCAO calculation. These can provide <tt>octopus</tt> with a good set
      !% of initial wave-functions and with a new guess for the density.
      !% (Up to current version, only a minimal basis set used.)
      !%Option lcao_none 0
      !% Do not perform a LCAO calculation before the SCF cycle.
      !%Option lcao_states 1
      !% Do a LCAO calculation before the SCF cycle and use the resulting wave-functions as 
      !% initail wave-functions without changing the guess density.
      !% This will speed-up the convergence of the eigensolver during the first SCF iterations.
      !%Option lcao_full 2
      !% Do a LCAO calculation before the SCF cycle and use the LCAO wave-functions to build a new
      !% guess density and a new KS potential.
      !% Using the LCAO density as a new guess density may improve the convergence, but can
      !% also slow it down or yield wrong results (especially for spin-polarized calculations).
      !%End
      call loct_parse_int(check_inp('LCAOStart'), lcao_start_default, lcao_start)
      if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
      call messages_print_var_option(stdout, 'LCAOStart', lcao_start)
      if (lcao_start > LCAO_START_NONE) then
          
        lcao_data%state = 0 ! Uninitialized here.
        call lcao_init(sys%gr, lcao_data, sys%st, h)
        if(lcao_data%state == 1) then
          write(message(1),'(a,i4,a)') 'Info: Performing initial LCAO calculation with ', &
               lcao_data%st%nst,' orbitals.'
          call write_info(1)
          
          call lcao_wf(lcao_data, sys%st, sys%gr%m, h)
          call lcao_end(lcao_data, sys%st%nst)

          !Just populate again the states, so that the eigenvalues are properly written
          call states_fermi(sys%st, sys%gr%m)
          call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb)

          if (lcao_start == LCAO_START_FULL) then
            ! Update the density and the Hamiltonian
            call system_h_setup(sys, h)
          end if

        end if
      end if

    else
      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call system_h_setup(sys, h)

    end if

    ! run self consistency
    if (sys%st%d%wfs_type == M_REAL) then
      message(1) = 'Info: SCF using real wavefunctions.'
    else
      message(1) = 'Info: SCF using complex wavefunctions.'
    end if
    call write_info(1)

    call scf_init(sys%gr, scfv, sys%st, h)
    call scf_run(scfv, sys%gr, sys%st, sys%ks, h, sys%outp)
    call scf_end(scfv)

    ! clean up
    call states_deallocate_wfns(sys%st)

    call pop_sub()
  end subroutine ground_state_run

end module ground_state_m
