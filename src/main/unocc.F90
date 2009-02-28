!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any %later version.
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
  use derivatives_m
  use eigensolver_m
  use energy_m
  use external_pot_m
  use global_m
  use grid_m
  use geometry_m
  use h_sys_output_m
  use hamiltonian_m
  use io_m
  use lcao_m
  use loct_m
  use loct_parser_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use poisson_m
  use projector_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use states_calc_m
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
    logical :: converged, l
    integer :: lcao_start, lcao_start_default, max_iter
    type(lcao_t) :: lcao

    ! read the maximum number of eigensolver iterations
    call loct_parse_int(datasets_check('MaximumIter'), 20, max_iter)

    occupied_states = sys%st%nst
    call init_(sys%gr%mesh, sys%st)
    total_states = sys%st%nst

    call restart_read (trim(restart_dir)//'gs', sys%st, sys%gr, sys%geo, ierr)

    if( (ierr .ne. 0)  .and.  (ierr < occupied_states) ) then
      message(1) = "Not all the occupied KS orbitals could be read from '"//trim(restart_dir)//"gs'"
      message(2) = "Please run a ground-state calculation first!"
      call write_fatal(2)
    end if

    if(ierr.ne.0) then
      message(1) = "Info:  Could not load all wave-functions from '"//trim(restart_dir)//"gs'"
      call write_info(1)
    end if

    if(fromScratch) then ! reset unoccupied states
      ierr = occupied_states
    end if

    ! Setup Hamiltonian
    message(1) = 'Info:  Setting up Hamiltonian.'
    call write_info(1)

    call states_calc_dens(sys%st, sys%gr%mesh%np)
    call v_ks_calc(sys%gr, sys%ks, hm, sys%st, calc_eigenval=.true.) ! get potentials
    call total_energy(hm, sys%gr, sys%st, -1)             ! total energy

    ! The initial LCAO calculation is done by default if we have pseudopotentials.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    lcao_start_default = LCAO_START_FULL
    if(sys%geo%only_user_def) lcao_start_default = LCAO_START_NONE

    call loct_parse_int(datasets_check('LCAOStart'), lcao_start_default, lcao_start)
    if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
    call messages_print_var_option(stdout, 'LCAOStart', lcao_start)

    if (lcao_start > LCAO_START_NONE) then
      if( (ierr.ne.0) .and. (ierr >= occupied_states)) then
        message(1) = "Info:  I will perform a LCAO calculation to get reasonable starting points."
        call write_info(1)
        call lcao_init(lcao, sys%gr, sys%geo, sys%st)
        if(lcao_is_available(lcao)) then
          call lcao_wf(lcao, sys%st, sys%gr, sys%geo, hm, start = ierr+1)
          call lcao_end(lcao)
        end if
      end if
    end if
    
    message(1) = "Info:  Starting calculation of unoccupied states"
    call write_info(1)

    ! reset this variable, so that the eigensolver passes through all states
    eigens%converged(:) = 0

    do iter = 1, max_iter
      write(message(1), '(a,i3)') "Info: Unoccupied states iteration ", iter
      call write_info(1)
      call eigensolver_run(eigens, sys%gr, sys%st, hm, 1, converged, verbose = .true.)

      if(converged) exit
    end do

    ! write restart information.
    call restart_write (trim(tmpdir)//'gs', sys%st, sys%gr, ierr)
    if(ierr.ne.0) then
      message(1) = 'Unsuccesful write of "'//trim(tmpdir)//'gs"'
      call write_fatal(1)
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

    ! calculate momentum of KS states
    if (sys%st%wfs_type == M_REAL) then
      call dstates_calc_momentum(sys%gr, sys%st)
    else
      call zstates_calc_momentum(sys%gr, sys%st)
    end if
    
    if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik>sys%st%d%nspin) then
      call states_write_bands(STATIC_DIR, sys%st%nst, sys%st, sys%gr%sb)
      call states_write_fermi_energy(STATIC_DIR, sys%st, sys%gr%mesh, sys%gr%sb)
      call states_degeneracy_matrix(sys%st)
    end if

    !%Variable WriteMatrixElements
    !%Type logical
    !%Default no
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% If true outputs the following matrix elements:
    !% <ul>
    !% <li><math>&lt;i|T + V_{ext}|j&gt;</math></li>
    !% <li><math>&lt;ij| 1/|r_1-r_2| |kl&gt;</math></li>
    !% </ul>
    !% in the directory matrix_elements
    !%End
    call loct_parse_logical(datasets_check('WriteMatrixElements'), .false., l)
    if(l) call write_matrix_elements(sys, hm)

    ! output wave-functions
    call h_sys_output_states(sys%st, sys%gr, STATIC_DIR, sys%outp)

    call end_()


  contains

    ! ---------------------------------------------------------
    subroutine init_(m, st)
      type(mesh_t),   intent(in)    :: m
      type(states_t), intent(inout) :: st

      integer :: nus

      call push_sub('unocc.unocc_run')

      !%Variable NumberUnoccStates
      !%Type integer
      !%Default 5
      !%Section Calculation Modes::Unoccupied States
      !%Description
      !% How many unoccupied states to compute.
      !%End
      call loct_parse_int(datasets_check('NumberUnoccStates'), 5, nus)
      if(nus <= 0) then
        message(1) = "Input: NumberUnoccStates must be > 0"
        call write_fatal(1)
      end if

      ! fix states: THIS IS NOT OK
      st%nst    = st%nst + nus
      st%st_end = st%nst

      deallocate(st%eigenval, st%momentum, st%occ)
      call states_allocate_wfns(st, m)
      ALLOCATE(st%eigenval(st%nst, st%d%nik), st%nst*st%d%nik)
      ALLOCATE(st%momentum(3, st%nst, st%d%nik), st%nst*st%d%nik)
      ALLOCATE(st%occ(st%nst, st%d%nik), st%nst*st%d%nik)
      if(st%d%ispin == SPINORS) then
        ALLOCATE(st%spin(3, st%nst, st%d%nik), st%nst*st%d%nik*2)
        st%spin = M_ZERO
      end if
      st%eigenval = huge(REAL_PRECISION)
      st%occ      = M_ZERO

      ! now the eigensolver stuff
      call eigensolver_init(sys%gr, eigens, st)

      ! Having initial and final tolerance does not make sense in this case:
      eigens%init_tol       = eigens%final_tol
      eigens%final_tol_iter = 2
      eigens%converged(1:st%d%nik) = st%nst - nus

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)
      call eigensolver_end(eigens)

      call pop_sub()
    end subroutine end_

  end subroutine unocc_run


  ! ---------------------------------------------------------
  ! warning: only works for spin-unpolarized and 1 k-point
  subroutine write_matrix_elements(sys, hm)
    type(system_t),         intent(inout) :: sys
    type(hamiltonian_t),    intent(in)    :: hm

    call io_mkdir("matrix_elements")

    message(1) = "Computing Matrix Elements"
    call write_info(1)

    message(1) = "  :: one-body"
    call write_info(1)
    if (sys%st%wfs_type == M_REAL) then
      call done_body(sys%gr, sys%geo, sys%st, hm)
    else
      call zone_body(sys%gr, sys%geo, sys%st, hm)
    end if

    message(1) = "  :: two-body"
    call write_info(1)
    if (sys%st%wfs_type == M_REAL) then
      call dtwo_body(sys%gr, sys%st)
    else
      call ztwo_body(sys%gr, sys%st)
    end if

  end subroutine write_matrix_elements


#include "undef.F90"
#include "real.F90"
#include "unocc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "unocc_inc.F90"

end module unocc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
