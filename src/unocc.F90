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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module unocc_m
  use global_m
  use messages_m
  use varinfo_m
  use datasets_m
  use lib_oct_parser_m
  use lib_oct_m
  use mesh_m
  use mesh_function_m
  use grid_m
  use states_m
  use system_m
  use restart_m
  use poisson_m
  use v_ks_m
  use hamiltonian_m
  use eigen_solver_m
  use io_m
  use simul_box_m
  use lcao_m

  implicit none

  private
  public :: &
    unocc_run


contains

  ! ---------------------------------------------------------
  subroutine unocc_run(sys, h, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    logical,             intent(inout) :: fromScratch

    type(eigen_solver_t) :: eigens
    integer :: iunit, ierr, ik, p, occupied_states
    FLOAT, allocatable :: dh_psi(:,:)
    CMPLX, allocatable :: zh_psi(:,:)
    logical :: converged, l
    integer :: lcao_start, lcao_start_default
    type(lcao_t) :: lcao_data

    occupied_states = sys%st%nst
    call init_(sys%gr%m, sys%st)

    call restart_read (trim(tmpdir)//'restart_gs', sys%st, sys%gr, sys%geo, ierr)
    if( (ierr .ne. 0)  .and.  (ierr < occupied_states) ) then
      message(1) = "Not all the occupied KS orbitals could be read from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a ground-state calculation first!"
      call write_fatal(2)
    end if

    if(ierr.ne.0) then
      message(1) = "Info:  Could not load all wave-functions from '"//trim(tmpdir)//"restart_gs'"
      call write_info(1)
    end if

    if(fromScratch) then ! reset unoccupied states
      ierr = occupied_states
    end if

    ! Setup Hamiltonian
    message(1) = 'Info:  Setting up Hamiltonian.'
    call write_info(1)

    call states_calc_dens(sys%st, sys%gr%m%np, sys%st%rho)
    call v_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
    call hamiltonian_energy(h, sys%gr, sys%geo, sys%st, -1)             ! total energy

    ! The initial LCAO calculation is done by default if we have pseudopotentials.
    ! Otherwise, it is not the default value and has to be enforced in the input file.
    lcao_start_default = LCAO_START_FULL
    if(sys%geo%only_user_def) lcao_start_default = LCAO_START_NONE

    call loct_parse_int(check_inp('LCAOStart'), lcao_start_default, lcao_start)
    if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
    call messages_print_var_option(stdout, 'LCAOStart', lcao_start)
    if (lcao_start > LCAO_START_NONE) then
      if( (ierr.ne.0) .and. (ierr >= occupied_states)) then
        message(1) = "Info:  I will perform a LCAO calculation to get reasonable starting points."
        call write_info(1)
        lcao_data%state = 0
        call lcao_init(sys%gr, sys%geo, lcao_data, sys%st, h)
        if(lcao_data%state .eq. 1) then
          call lcao_wf(lcao_data, sys%st, sys%gr%m, h, start = ierr+1)
          call lcao_end(lcao_data, sys%st%nst)
        end if
      end if
    end if

    message(1) = "Info:  Starting calculation of unoccupied states"
    call write_info(1)

    ! First, get the residues of the occupied states.
    ! These are assumed to be converged; otherwise one should do a SCF calculation.
    if (sys%st%d%wfs_type == M_REAL) then
      ALLOCATE(dh_psi(sys%gr%m%np, h%d%dim), sys%gr%m%np*h%d%dim)
    else
      ALLOCATE(zh_psi(sys%gr%m%np, h%d%dim), sys%gr%m%np*h%d%dim)
    end if
    do ik = 1, sys%st%d%nik
      do p = 1, eigens%converged
        if (sys%st%d%wfs_type == M_REAL) then
          call dHpsi(h, sys%gr, sys%st%dpsi(:,:, p, ik) , dh_psi, ik)
          eigens%diff(p, ik) = dstates_residue(sys%gr%m, sys%st%d%dim, dh_psi, sys%st%eigenval(p, ik), &
               sys%st%dpsi(:, :, p, ik))
        else
          call zHpsi(h, sys%gr, sys%st%zpsi(:,:, p, ik) , zh_psi, ik)
          eigens%diff(p, ik) = zstates_residue(sys%gr%m, sys%st%d%dim, zh_psi, sys%st%eigenval(p, ik), &
               sys%st%zpsi(:, :, p, ik))
        end if
      end do
    end do
    if (sys%st%d%wfs_type == M_REAL) then
      deallocate(dh_psi)
    else
      deallocate(zh_psi)
    end if

    call eigen_solver_run(eigens, sys%gr, sys%st, h, 1, converged, verbose = .true.)

    ! write restart information.
    call restart_write (trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    if(ierr.ne.0) then
      message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_gs"'
      call write_fatal(1)
    end if

    ! write output file
    call io_mkdir('static')
    iunit = io_open('static/eigenvalues', action='write')

    if(converged) then
      write(iunit,'(a)') 'All unoccupied states converged.'
    else
      write(iunit,'(a)') 'Some of the unoccupied states are not fully converged!'
    end if
    write(iunit,'(a, e17.6)') 'Criterium = ', eigens%final_tol
    write(iunit,'(1x)')
    call states_write_eigenvalues(iunit, sys%st%nst, sys%st, sys%gr%sb, eigens%diff)
    call io_close(iunit)

    ! calculate momentum of KS states
    if (sys%st%d%wfs_type == M_REAL) then
      call dstates_calc_momentum(sys%gr, sys%st)
    else
      call zstates_calc_momentum(sys%gr, sys%st)
    end if
    
    if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik>sys%st%d%nspin) then
      call states_write_bands('static', sys%st%nst, sys%st, sys%gr%sb)
      call states_write_dos  ('static', sys%st)
      call states_write_fermi_energy('static', sys%st, sys%gr%m, sys%gr%sb)
      call states_degeneracy_matrix(sys%st)
    end if

    !%Variable WriteMatrixElements
    !%Type logical
    !%Default no
    !%Section Unoccupied States
    !%Description
    !% If true outputs the following matrix elements:
    !% <ul>
    !% <li><math>&lt;i|T + V_{ext}|j&gt;</math></li>
    !% <li><math>&lt;ij| 1/|r_1-r_2| |kl&gt;</math></li>
    !% </ul>
    !% in the directory ME
    !%End
    call loct_parse_logical(check_inp('WriteMatrixElements'), .false., l)
    if(l) call write_matrix_elements(sys, h)

    ! output wave-functions
    call states_output(sys%st, sys%gr, "static", sys%outp)

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
      !%Section Unoccupied States
      !%Description
      !% How many unoccupied states to compute.
      !%End
      call loct_parse_int(check_inp('NumberUnoccStates'), 5, nus)
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
        ALLOCATE(st%mag(st%nst, st%d%nik, 2), st%nst*st%d%nik*2)
        st%mag = M_ZERO
      end if
      st%eigenval = huge(PRECISION)
      st%occ      = M_ZERO

      ! now the eigen solver stuff
      call eigen_solver_init(sys%gr, eigens, st, 50)

      ! Having initial and final tolerance does not make sense in this case:
      eigens%init_tol       = eigens%final_tol
      eigens%final_tol_iter = 2
      eigens%converged      = st%nst - nus

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)
      call eigen_solver_end(eigens)

      call pop_sub()
    end subroutine end_

  end subroutine unocc_run


  ! ---------------------------------------------------------
  ! warning: only works for spin-unpolarized and 1 k-point
  subroutine write_matrix_elements(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(in)    :: h

    call io_mkdir("ME")

    message(1) = "Computing Matrix Elements"
    call write_info(1)

    message(1) = "  :: one-body"
    call write_info(1)
    if (sys%st%d%wfs_type == M_REAL) then
      call done_body(sys%gr%m, sys%st, h)
    else
      call zone_body(sys%gr%m, sys%st, h)
    end if

    message(1) = "  :: two-body"
    call write_info(1)
    if (sys%st%d%wfs_type == M_REAL) then
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
