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

  private
  public ::           &
    ground_state_run, &
    ground_state_init


contains

  ! ---------------------------------------------------------
  subroutine ground_state_run(sys, h, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(scf_t) :: scfv
    integer        :: ierr

    call push_sub('gs.ground_state_run')

    ! allocate wfs
    ALLOCATE(sys%st%X(psi)(sys%NP_PART, sys%st%d%dim, sys%st%nst, sys%st%d%nik), sys%NP_PART*sys%st%d%dim*sys%st%nst*sys%st%d%nik)

    if(.not.fromScratch) then

      ! load wave-functions
      call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
      if(ierr.ne.0) then
        message(1) = "Could not load wave-functions from '"//trim(tmpdir)//"restart_gs'"
        message(2) = "Starting from scratch!"
        call write_warning(1)
        fromScratch = .true.
      end if
    end if

    if(fromScratch) then
      call ground_state_init(sys%gr, sys%st, sys%ks, h)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    ! run self consistency
#ifdef COMPLEX_WFNS
    message(1) = 'Info: SCF using complex wavefunctions.'
#else
    message(1) = 'Info: SCF using real wavefunctions.'
#endif
    call write_info(1)

    call scf_init(sys%gr, scfv, sys%st, h)
    call scf_run(scfv, sys%gr, sys%st, sys%ks, h, sys%outp)
    call scf_end(scfv)

    ! clean up
    deallocate(sys%st%X(psi)); nullify(sys%st%X(psi))
    call pop_sub()
  end subroutine ground_state_run


  ! ---------------------------------------------------------
  subroutine ground_state_init(gr, st, ks,  h)
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h

    integer :: err
    logical :: lcao_start, lcao_start_default = .true.
    type(lcao_t) :: lcao_data
#if defined(HAVE_MPI)
    integer :: mpi_err
#endif

    call push_sub('gs.ground_state_init')

    ! set barrier before the first communication takes place
    ! this ensures proper debug timing of MPI calls
#if defined(HAVE_MPI)
    call TS(MPI_Barrier)(MPI_COMM_WORLD, mpi_err)
#endif

    message(1) = 'Info: Random generating starting wavefunctions.'
    call write_info(1)

    ! wave functions are simply random gaussians
    call states_generate_random(st, gr%m)

    ! this is certainly a better density
    call system_guess_density(gr%m, gr%sb, gr%geo, st%qtot, st%d%nspin, &
      st%d%spin_channels, st%rho)

    ! The initial LCAO calculation is done by default if we have pseudopotentials.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    if(gr%geo%only_user_def) lcao_start_default = .false.

    !%Variable LCAOStart
    !%Type logical
    !%Default yes
    !%Section SCF
    !%Description
    !% Before starting a SCF calculation, performs
    !% a LCAO calculation. These should provide <tt>octopus</tt> with a good set
    !% of initial wave-functions, and help the convergence of the SCF cycle.
    !% (Up to current version, only a minimal basis set used.)
    !%End
    call loct_parse_logical(check_inp('LCAOStart'), lcao_start_default, lcao_start)
    if(lcao_start) then
      call X(v_ks_calc)(gr, ks, h, st, calc_eigenval=.true.)

      lcao_data%state = 0 ! Uninitialized here.
      call lcao_init(gr, lcao_data, st, h)
      if(lcao_data%state == 1) then
        write(message(1),'(a,i4,a)') 'Info: Performing initial LCAO calculation with ', &
          lcao_data%st%nst,' orbitals.'
        call write_info(1)

        call lcao_wf(lcao_data, st, gr%m, gr%sb, h)
        call lcao_end(lcao_data, st%nst)

        call states_fermi(st, gr%m)                         ! occupations
        call states_write_eigenvalues(stdout, st%nst, st, gr%sb)
      end if
    end if

    ! write wave-functions to disk
    call X(restart_write)(trim(tmpdir)//'restart_gs', st, gr, err)
    if(err.ne.0) then
      message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_gs"'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine ground_state_init

end module ground_state_m
