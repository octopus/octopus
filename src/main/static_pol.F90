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
  use simul_box_m
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

    integer :: iunit, ios, i_start, ii, jj, kk, is, isign, ierr
    FLOAT :: e_field
    FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
    FLOAT, allocatable :: elf(:,:), lr_elf(:,:), elfd(:,:), lr_elfd(:,:)
    FLOAT, allocatable :: lr_rho(:,:), lr_rho2(:,:), gs_rho(:,:)
    FLOAT :: center_dipole(1:MAX_DIM)
    logical :: out_pol
    character(len=80) :: fname

    out_pol = .false.

    call init_()

    ! load wave-functions
    call restart_read(trim(restart_dir)//'gs', sys%st, gr, sys%geo, ierr)

    if(simul_box_is_periodic(gr%sb)) then
      message(1) = "Electric field cannot be applied to a periodic system (currently)."
      call write_fatal(1)
    endif

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
    SAFE_ALLOCATE(dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 1:2))
    dipole = M_ZERO

    if(.not.fromScratch) then
      iunit = io_open(trim(tmpdir)//'restart.pol', action='read', status='old', die=.false.)
      if(iunit > 0) then
        ! Finds out how many dipoles have already been written.
        rewind(iunit)
        i_start = 1
        do ii = 1, 3
          read(iunit, fmt=*, iostat = ios) ((dipole(ii, jj, isign), jj = 1, gr%mesh%sb%dim), isign = 1, 2)
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
    SAFE_ALLOCATE(Vpsl_save(1:gr%mesh%np))
    Vpsl_save = hm%ep%Vpsl

    ! Allocate the trrho to the contain the trace of the density.
    SAFE_ALLOCATE(trrho(1:gr%mesh%np))
    SAFE_ALLOCATE(gs_rho(1:gr%mesh%np, 1:st%d%nspin))
    trrho = M_ZERO
    gs_rho = M_ZERO

    call output_init_()

    call scf_init(scfv, gr, sys%geo, st, hm)

    ! now calculate the dipole without field

    hm%ep%vpsl(1:gr%mesh%np) = vpsl_save(1:gr%mesh%np)
    call hamiltonian_update_potential(hm, gr%mesh)

    write(message(1), '(a)')
    write(message(2), '(a)')'Info: Calculating dipole moment for zero field.'
    call write_info(2)
    call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)
    
    gs_rho(1:gr%mesh%np, 1:st%d%nspin) = st%rho(1:gr%mesh%np, 1:st%d%nspin)
    trrho = M_ZERO
    do is = 1, st%d%spin_channels
      trrho(1:gr%mesh%np) = trrho(1:gr%mesh%np) + gs_rho(1:gr%mesh%np, is)
    end do

    ! calculate dipole
    do jj = 1, gr%mesh%sb%dim
       center_dipole(jj) = dmf_moment(gr%mesh, trrho, jj, 1)
    end do

    do ii = i_start, gr%mesh%sb%dim
      do isign = 1, 2
        write(message(1), '(a)')
        write(message(2), '(a,f6.4,a,a,a,a,a,i1)')'Info: Calculating dipole moment for field ', &
          units_from_atomic(units_inp%energy * units_inp%length, -(-1)**isign *e_field), ' ', &
          trim(units_inp%energy%abbrev), '/', trim(units_inp%length%abbrev), ' in direction ', ii
        call write_info(2)
        ! there is an extra factor of -1 in here that is for the electronic charge

        hm%ep%vpsl(1:gr%mesh%np) = vpsl_save(1:gr%mesh%np) + (-1)**isign*gr%mesh%x(1:gr%mesh%np, ii)*e_field
        call hamiltonian_update_potential(hm, gr%mesh)

        call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)

        trrho = M_ZERO
        do is = 1, st%d%spin_channels
          trrho(1:gr%mesh%np) = trrho(1:gr%mesh%np) + st%rho(1:gr%mesh%np, is)
        end do

        ! calculate dipole
        do jj = 1, gr%mesh%sb%dim
          dipole(ii, jj, isign) = dmf_moment(gr%mesh, trrho, jj, 1)
        end do

        call output_cycle_()

      end do

      ! Writes the dipole to file
      if(mpi_grp_is_root(mpi_world)) then 
        iunit = io_open(trim(tmpdir)//'restart.pol', action='write', status='old', position='append')
        write(iunit, fmt='(6e20.12)') ((dipole(ii, jj, isign), jj = 1, gr%mesh%sb%dim), isign = 1, 2)
        call io_close(iunit)
      end if

      out_pol = (ii == gr%mesh%sb%dim)

    end do
    
    call scf_end(scfv)
    call output_end_()

    SAFE_DEALLOCATE_A(Vpsl_save)
    SAFE_DEALLOCATE_A(trrho)
    SAFE_DEALLOCATE_A(gs_rho)
    SAFE_DEALLOCATE_A(dipole)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol.static_pol_run')

      ! shortcuts
      gr  => sys%gr
      st => sys%st

      call states_allocate_wfns(st, gr%mesh)

      !%Variable EMStaticField
      !%Type float
      !%Default 0.01 a.u.
      !%Section Linear Response::Static Polarization
      !%Description
      !% Magnitude of the static field used to calculate the static polarizability,
      !% if ResponseMethod = finite_differences.
      !%End
      call loct_parse_float(datasets_check('EMLStaticField'), &
         units_from_atomic(units_inp%energy / units_inp%length, CNST(0.01)), e_field)
      e_field = units_to_atomic(units_inp%energy / units_inp%length, e_field)
      if (e_field <= M_ZERO) then
        write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid EMStaticField"
        message(2) = '(0 < EMStaticField)'
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
      if(iand(sys%outp%what, output_density).ne.0 .or. &
         iand(sys%outp%what, output_pol_density).ne.0 ) then 
        SAFE_ALLOCATE(lr_rho (1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(lr_rho2(1:gr%mesh%np, 1:st%d%nspin))
      end if
      
      if(iand(sys%outp%what, output_elf).ne.0) then 
        SAFE_ALLOCATE(    elf(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE( lr_elf(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(   elfd(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(lr_elfd(1:gr%mesh%np, 1:st%d%nspin))
      end if
      
    end subroutine output_init_


    !-------------------------------------------------------------
    subroutine output_cycle_()
      
      !DENSITY AND POLARIZABILITY DENSITY   
      if(iand(sys%outp%what, output_density).ne.0 .or. &
         iand(sys%outp%what, output_pol_density).ne.0) then 
         
        if(isign == 1) then 
          lr_rho(1:gr%mesh%np, 1:st%d%nspin) = st%rho(1:gr%mesh%np, 1:st%d%nspin)
        else
          lr_rho2(1:gr%mesh%np, 1:st%d%nspin) = &
            -(st%rho(1:gr%mesh%np, 1:st%d%nspin) + lr_rho(1:gr%mesh%np, 1:st%d%nspin) - 2 * gs_rho(1:gr%mesh%np, 1:st%d%nspin)) &
            / e_field**2

          lr_rho(1:gr%mesh%np, 1:st%d%nspin) = &
               (st%rho(1:gr%mesh%np, 1:st%d%nspin) - lr_rho(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            if(iand(sys%outp%what, output_density).ne.0) then
              write(fname, '(a,a,i1,a,i1)') 'fd_density', '-', is, '-', ii
              call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%mesh, gr%sb, lr_rho(:, is), M_ONE, ierr, geo = sys%geo)

              ! save the trouble of writing many copies of each density, since i,j = j,i
              do jj = ii, gr%mesh%sb%dim
                write(fname, '(a,a,i1,a,i1,a,i1)') 'fd2_density', '-', is, '-', ii, '-', jj
                call doutput_function(sys%outp%how, "linear", trim(fname),&
                  gr%mesh, gr%sb, lr_rho2(:, is), M_ONE, ierr, geo = sys%geo)
              enddo
            endif

            if(iand(sys%outp%what, output_pol_density).ne.0) then
              do jj = ii, gr%mesh%sb%dim
                write(fname, '(a,a,i1,a,i1,a,i1)') 'fd_pol_density', '-', is, '-', ii, '-', jj
                call doutput_function(sys%outp%how, "linear", trim(fname),&
                  gr%mesh, gr%sb, gr%mesh%x(:, jj) * lr_rho(:, is), M_ONE, ierr, geo = sys%geo)

                write(fname, '(a,a,i1,a,i1,a,i1,a,i1)') 'fd2_pol_density', '-', is, '-', ii, '-', ii, '-', jj
                call doutput_function(sys%outp%how, "linear", trim(fname),&
                  gr%mesh, gr%sb, gr%mesh%x(:, jj) * lr_rho2(:, is), M_ONE, ierr, geo = sys%geo)
              enddo
            endif
          end do

        end if
      end if

      !ELF
      if(iand(sys%outp%what, output_elf).ne.0) then 
         
        if(isign == 1) then 
          call elf_calc(st, gr, elf, elfd)
        else
          call elf_calc(st, gr, lr_elf, lr_elfd)
          
          !numerical derivative
          lr_elf(1:gr%mesh%np, 1:st%d%nspin) = &
                ( lr_elf(1:gr%mesh%np, 1:st%d%nspin) -  elf(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO*e_field)
          lr_elfd(1:gr%mesh%np, 1:st%d%nspin) = &
               (lr_elfd(1:gr%mesh%np, 1:st%d%nspin) - elfd(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,a,i1,a,i1)') 'fd_elf', '-', is, '-', ii
            call doutput_function(sys%outp%how, "linear", trim(fname),&
                gr%mesh, gr%sb, lr_elf(:, is), M_ONE, ierr, geo = sys%geo)
            write(fname, '(a,a,i1,a,i1)') 'fd_D', '-', is, '-', ii
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
        ! WARNING: terms beta(ijk) with i,j,k all different are not calculated here, and just left as zero.

        call io_output_tensor(iunit, alpha, gr%mesh%sb%dim, units_out%length%factor**gr%mesh%sb%dim)
        call io_close(iunit)
        
        call out_hyperpolarizability(gr%sb, beta, converged = .true., dirname = "linear/")

      end if

      if(iand(sys%outp%what, output_density).ne.0 .or. &
         iand(sys%outp%what, output_pol_density).ne.0) then 
        SAFE_DEALLOCATE_A(lr_rho)
      end if
      
      if(iand(sys%outp%what, output_elf).ne.0) then 
        SAFE_DEALLOCATE_A(lr_elf)
        SAFE_DEALLOCATE_A(elf)
        SAFE_DEALLOCATE_A(lr_elfd)
        SAFE_DEALLOCATE_A(elfd)
      end if
    end subroutine output_end_

  end subroutine static_pol_run

end module static_pol_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
