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

module static_pol_m
  use global_m
  use messages_m
  use datasets_m
  use units_m
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use mesh_function_m
  use mesh_m
  use grid_m
  use system_m
  use hamiltonian_m
  use states_m
  use geometry_m
  use restart_m
  use scf_m
  use output_m

  implicit none

  private
  public :: &
    static_pol_run


contains

  ! ---------------------------------------------------------
  subroutine static_pol_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(scf_t)             :: scfv
    type(grid_t),   pointer :: gr    ! shortcuts
    type(states_t), pointer :: st

    integer :: iunit, ios, i_start, i, j, is, k, ierr
    FLOAT :: e_field
    FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
    FLOAT, allocatable :: elf(:,:), lr_elf(:,:), elfd(:,:), lr_elfd(:,:)
    FLOAT, allocatable :: lr_rho(:,:)
    logical :: out_pol
    character(len=80) :: fname

    call init_()

    ! load wave-functions
    call restart_read(trim(tmpdir)//'restart_gs', sys%st, gr, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a ground-state calculation first!"
      call write_fatal(2)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup (sys, h)

    ! Allocate the dipole...
    ALLOCATE(dipole(NDIM, NDIM, 2), NDIM*NDIM*2)
    dipole = M_ZERO

    if(.not.fromScratch) then
      iunit = io_open(trim(tmpdir)//'restart.pol', action='read', status='old', die=.false.)
      if(iunit > 0) then
        ! Finds out how many dipoles have already been written.
        rewind(iunit)
        i_start = 1
        do i = 1, 3
          read(iunit, fmt=*, iostat = ios) ((dipole(i, j, k), j = 1, NDIM), k = 1, 2)
          if(ios.ne.0) exit
          i_start = i_start + 1
        end do
        call io_close(iunit)
      else
        fromScratch = .true.
      end if
    end if

    if(fromScratch) then
      iunit = io_open(trim(tmpdir)//'restart.pol', action='write')
      call io_close(iunit)
      i_start = 1
    end if

    ! Save local pseudopotential
    ALLOCATE(Vpsl_save(NP), NP)
    Vpsl_save = h%ep%Vpsl

    ! Allocate the trrho to the contain the trace of the density.
    ALLOCATE(trrho(NP), NP)
    trrho = M_ZERO

    call output_init_()

    call scf_init(gr, scfv, st, h)
    do i = i_start, NDIM
      do k = 1, 2
        write(message(1), '(a)')
        write(message(2), '(a,i1,a,i1)')'Info: Calculating dipole moment for field ', i, ', #',k
        call write_info(2)

        h%ep%vpsl(1:NP) = vpsl_save(1:NP) + (-1)**k*gr%m%x(1:NP, i)*e_field

        call scf_run(scfv, sys%gr, st, sys%ks, h, sys%outp, gs_run=.false.)

        trrho = M_ZERO
        do is = 1, st%d%spin_channels
          trrho(:) = trrho(:) + st%rho(:, is)
        end do

        ! calculate dipole
        do j = 1, NDIM
          dipole(i, j, k) = dmf_moment(gr%m, trrho, j, 1)
        end do

        call output_()

      end do

      ! Writes down the dipole to the file
      iunit = io_open(trim(tmpdir)//'restart.pol', action='write', status='old', position='append')
      write(iunit, fmt='(6e20.12)') ((dipole(i, j, k), j = 1, NDIM), k = 1, 2)
      call io_close(iunit)

      if(i == NDIM) then
        out_pol = .true.
      end if
    end do
    call scf_end(scfv)

    if(out_pol) then ! output pol file
      call io_mkdir('linear')
      iunit = io_open('linear/polarizability', action='write')
      write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
        trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM
      write(iunit, '(a)') ']'

      do j = 1, NDIM
        write(iunit, '(3f12.6)') (dipole(j, 1:NDIM, 1) - dipole(j, 1:NDIM, 2))/(M_TWO*e_field) &
          / units_out%length%factor**NDIM
      end do
      call io_close(iunit)
    end if

    deallocate(Vpsl_save, trrho, dipole)
    call output_end_()
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol.static_pol_run')

      ! shortcuts
      gr  => sys%gr
      st => sys%st

      call states_allocate_wfns(st, gr%m)

      !%Variable POLStaticField
      !%Type float
      !%Default 0.01 a.u.
      !%Section Linear Response::Static Polarization
      !%Description
      !% Magnitude of the static field used to calculate the static polarizability
      !% (<tt>CalculationMode = pol</tt>)
      !%End
      call loct_parse_float(check_inp('POLStaticField'), &
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

    subroutine output_()
      
      !DENSITY      
      if(iand(sys%outp%what, output_density).ne.0) then 
         
        if(k == 1) then 
          lr_rho(1:NP, 1:st%d%nspin) = st%rho(1:NP, 1:st%d%nspin)
        else
          lr_rho(1:NP, 1:st%d%nspin) = (st%rho(1:NP, 1:st%d%nspin) - lr_rho(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,a,i1,a,i1)') 'fd_density', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%m, gr%sb, lr_rho(:, is), M_ONE, ierr)
          end do

        end if
      end if

      !ELF
      if(iand(sys%outp%what, output_elf).ne.0) then 
         
        if(k == 1) then 
          call states_calc_elf(st, gr, elf, elfd)
        else
          call states_calc_elf(st, gr, lr_elf, lr_elfd)
          
          !numerical derivative
           lr_elf(1:NP, 1:st%d%nspin) = ( lr_elf(1:NP, 1:st%d%nspin) -  elf(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)
          lr_elfd(1:NP, 1:st%d%nspin) = (lr_elfd(1:NP, 1:st%d%nspin) - elfd(1:NP, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,a,i1,a,i1)') 'fd_elf', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%m, gr%sb, lr_elf(:, is), M_ONE, ierr)
            write(fname, '(a,a,i1,a,i1)') 'fd_D', '-', is, '-', i
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%m, gr%sb, lr_elfd(:, is), M_ONE, ierr)
          end do

        end if

      end if

    end subroutine output_

    subroutine output_end_()

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
