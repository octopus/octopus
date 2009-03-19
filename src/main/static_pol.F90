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

module static_pol_m
  use datasets_m
  use elf_m
  use em_resp_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_function_m
  use io_m
  use loct_m
  use loct_parser_m
  use mpi_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use h_sys_output_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m

  implicit none

  private
  public :: &
    static_pol_run


contains

  ! ---------------------------------------------------------
  subroutine static_pol_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch
    type(scf_t)             :: scfv
    type(grid_t),   pointer :: gr    ! shortcuts
    type(states_t), pointer :: st

    integer :: iunit, ios, i_start, i, j, is, k, ierr
    FLOAT :: e_field
    FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
    FLOAT, allocatable :: elf(:,:), lr_elf(:,:), elfd(:,:), lr_elfd(:,:)
    FLOAT, allocatable :: lr_rho(:,:)
    FLOAT :: center_dipole(1:MAX_DIM)
    logical :: out_pol
    character(len=80) :: fname

    out_pol = .false.

    call init_()

    ! load wave-functions
    call restart_read(trim(restart_dir)//'gs', sys%st, gr, sys%geo, ierr)

    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(restart_dir)//"gs'"
      message(2) = "Please run a ground-state calculation first and/or"
      message(3) = "Give the correct RestartDataset in the input file."
      call write_fatal(3)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup (sys, hm)

    ! Allocate the dipole...
    ALLOCATE(dipole(gr%mesh%sb%dim, gr%mesh%sb%dim, 2), gr%mesh%sb%dim*gr%mesh%sb%dim*2)
    dipole = M_ZERO

    if(.not.fromScratch) then
      iunit = io_open(trim(tmpdir)//'restart.pol', action='read', status='old', die=.false.)
      if(iunit > 0) then
        ! Finds out how many dipoles have already been written.
        rewind(iunit)
        i_start = 1
        do i = 1, 3
          read(iunit, fmt=*, iostat = ios) ((dipole(i, j, k), j = 1, gr%mesh%sb%dim), k = 1, 2)
          if(ios.ne.0) exit
          i_start = i_start + 1
        end do
        call io_close(iunit)
      else
        fromScratch = .true.
      end if
    end if

    if(fromScratch) then
      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(trim(tmpdir)//'restart.pol', action='write')
        call io_close(iunit)
      end if
      i_start = 1
    end if
    if(i_start > gr%mesh%sb%dim) out_pol = .true.

    ! Save local pseudopotential
    ALLOCATE(Vpsl_save(NP), NP)
    Vpsl_save = hm%ep%Vpsl

    ! Allocate the trrho to the contain the trace of the density.
    ALLOCATE(trrho(NP), NP)
    trrho = M_ZERO

    call output_init_()

    call scf_init(scfv, gr, sys%geo, st, hm)
    do i = i_start, gr%mesh%sb%dim
      do k = 1, 2
        write(message(1), '(a)')
        write(message(2), '(a,i1,a,i1)')'Info: Calculating dipole moment for field ', i, ', #',k
        call write_info(2)

        hm%ep%vpsl(1:NP) = vpsl_save(1:NP) + (-1)**k*gr%mesh%x(1:NP, i)*e_field
        call hamiltonian_update_potential(hm, gr%mesh)

        call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)

        trrho = M_ZERO
        do is = 1, st%d%spin_channels
          trrho(1:NP) = trrho(1:NP) + st%rho(1:NP, is)
        end do

        ! calculate dipole
        do j = 1, gr%mesh%sb%dim
          dipole(i, j, k) = dmf_moment(gr%mesh, trrho, j, 1)
        end do

        call output_cycle_()

      end do

      ! Writes down the dipole to the file
      if(mpi_grp_is_root(mpi_world)) then 
        iunit = io_open(trim(tmpdir)//'restart.pol', action='write', status='old', position='append')
        write(iunit, fmt='(6e20.12)') ((dipole(i, j, k), j = 1, gr%mesh%sb%dim), k = 1, 2)
        call io_close(iunit)
      end if

      out_pol = (i == gr%mesh%sb%dim)

    end do
    
    ! now calculate the dipole without field

    hm%ep%vpsl(1:NP) = vpsl_save(1:NP)
    call hamiltonian_update_potential(hm, gr%mesh)

    call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)
    
    trrho = M_ZERO
    do is = 1, st%d%spin_channels
      trrho(1:NP) = trrho(1:NP) + st%rho(1:NP, is)
    end do

        ! calculate dipole
    do j = 1, gr%mesh%sb%dim
       center_dipole(j) = dmf_moment(gr%mesh, trrho, j, 1)
    end do

    call scf_end(scfv)
    call output_end_()

    deallocate(Vpsl_save, trrho, dipole)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol.static_pol_run')

      ! shortcuts
      gr  => sys%gr
      st => sys%st

      call states_allocate_wfns(st, gr%mesh)

      !%Variable POLStaticField
      !%Type float
      !%Default 0.01 a.u.
      !%Section Linear Response::Static Polarization
      !%Description
      !% Magnitude of the static field used to calculate the static polarizability
      !% (<tt>CalculationMode = pol</tt>)
      !%End
      call loct_parse_float(datasets_check('POLStaticField'), &
         CNST(0.01)/units_inp%energy%factor*units_inp%length%factor, e_field)
      e_field = e_field * units_inp%energy%factor / units_inp%length%factor
      if (e_field <= M_ZERO) then
        write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid POLStaticField"
        message(2) = '(0 < POLStaticField)'
        call write_fatal(2)
      end if

    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)

      call pop_sub()
    end subroutine end_
  

    !-------------------------------------------------------------
    subroutine output_init_()

      !allocate memory for what we want to output
      if(iand(sys%outp%what, output_density).ne.0) then 
        ALLOCATE(lr_rho(1:NP, 1:st%d%nspin), NP*st%d%nspin)
      end if
      
      if(iand(sys%outp%what, output_elf).ne.0) then 
        ALLOCATE(    elf(1:NP, 1:st%d%nspin), NP*st%d%nspin)
        ALLOCATE( lr_elf(1:NP, 1:st%d%nspin), NP*st%d%nspin)
        ALLOCATE(   elfd(1:NP, 1:st%d%nspin), NP*st%d%nspin)
        ALLOCATE(lr_elfd(1:NP, 1:st%d%nspin), NP*st%d%nspin)
      end if
      
    end subroutine output_init_


    !-------------------------------------------------------------
    subroutine output_cycle_()
      
      !DENSITY      
      if(iand(sys%outp%what, output_density).ne.0) then 
         
        if(k == 1) then 
          lr_rho(1:NP, 1:st%d%nspin) = st%rho(1:NP, 1:st%d%nspin)
        else
          lr_rho(1:NP, 1:st%d%nspin) = &
               (st%rho(1:NP, 1:st%d%nspin) - lr_rho(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,a,i1,a,i1)') 'fd_density', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%mesh, gr%sb, lr_rho(:, is), M_ONE, ierr, geo = sys%geo)
          end do

        end if
      end if

      !ELF
      if(iand(sys%outp%what, output_elf).ne.0) then 
         
        if(k == 1) then 
          call elf_calc(st, gr, elf, elfd)
        else
          call elf_calc(st, gr, lr_elf, lr_elfd)
          
          !numerical derivative
           lr_elf(1:NP, 1:st%d%nspin) = &
                ( lr_elf(1:NP, 1:st%d%nspin) -  elf(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)
          lr_elfd(1:NP, 1:st%d%nspin) = &
               (lr_elfd(1:NP, 1:st%d%nspin) - elfd(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,a,i1,a,i1)') 'fd_elf', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%mesh, gr%sb, lr_elf(:, is), M_ONE, ierr, geo = sys%geo)
            write(fname, '(a,a,i1,a,i1)') 'fd_D', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%mesh, gr%sb, lr_elfd(:, is), M_ONE, ierr, geo = sys%geo)
          end do
        end if

      end if

    end subroutine output_cycle_


    !-------------------------------------------------------------
    subroutine output_end_()
      FLOAT :: alpha(MAX_DIM, MAX_DIM)
      CMPLX :: beta(MAX_DIM, MAX_DIM, MAX_DIM)
      integer :: iunit, idir

      call io_mkdir('linear')
      if(out_pol  .and.  mpi_grp_is_root(mpi_world)) then ! output pol file
        iunit = io_open('linear/polarizability_fd', action='write')
        write(iunit, '(2a)', advance='no') '# Polarizability tensor [', &
          trim(units_out%length%abbrev)
        if(gr%mesh%sb%dim.ne.1) write(iunit, '(a,i1)', advance='no') '^', gr%mesh%sb%dim
        write(iunit, '(a)') ']'

        alpha(1:gr%mesh%sb%dim,1:gr%mesh%sb%dim) = (dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 1) - &
             dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 2))/(M_TWO*e_field)

        beta = M_ZERO

        do idir = 1, gr%mesh%sb%dim
          beta(1:gr%mesh%sb%dim, idir, idir) = -(dipole(idir, 1:gr%mesh%sb%dim, 1) + dipole(idir, 1:gr%mesh%sb%dim, 2) - &
            M_TWO*center_dipole(1:gr%mesh%sb%dim))/e_field**2
          beta(idir, 1:gr%mesh%sb%dim, idir) = beta(1:gr%mesh%sb%dim, idir, idir) 
          beta(idir, idir, 1:gr%mesh%sb%dim) = beta(1:gr%mesh%sb%dim, idir, idir)
        end do

        call io_output_tensor(iunit, alpha, gr%mesh%sb%dim, units_out%length%factor**gr%mesh%sb%dim)
        call io_close(iunit)
        
        call out_hyperpolarizability(gr%sb, beta, converged = .true., dirname = "linear/")

      end if

      if(iand(sys%outp%what, output_density).ne.0) then 
        deallocate(lr_rho)
      end if
      
      if(iand(sys%outp%what, output_elf).ne.0) then 
        deallocate(lr_elf)
        deallocate(elf)
        deallocate(lr_elfd)
        deallocate(elfd)
      end if
    end subroutine output_end_

  end subroutine static_pol_run

end module static_pol_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
