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
#define RESTART_FILE 'dipoles'

module static_pol_m
  use born_charges_m
  use datasets_m
  use elf_m
  use em_resp_m
  use geometry_m
  use global_m
  use grid_m
  use h_sys_output_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lcao_m
  use loct_m
  use parser_m
  use mpi_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use restart_m
  use scf_m
  use simul_box_m
  use species_m
  use states_m
  use system_m
  use unit_m
  use unit_system_m

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

    integer :: iunit, ios, i_start, ii, jj, is, isign, ierr
    FLOAT :: e_field
    FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
    FLOAT, allocatable :: elf(:,:), lr_elf(:,:), elfd(:,:), lr_elfd(:,:)
    FLOAT, allocatable :: lr_rho(:,:), lr_rho2(:,:), gs_rho(:,:), tmp_rho(:,:)
    FLOAT :: center_dipole(1:MAX_DIM), diag_dipole(1:MAX_DIM), ionic_dipole(1:MAX_DIM), print_dipole(1:MAX_DIM)
    type(born_charges_t) :: born_charges
    logical :: diagonal_done, calc_Born, start_density_is_zero_field, center_written, calc_diagonal
    character(len=80) :: fname

    call push_sub('static_pol.static_pol_run')

    call init_()

    ! load wavefunctions
    call restart_read(trim(restart_dir)//GS_DIR, sys%st, gr, sys%geo, ierr)

    if(simul_box_is_periodic(gr%sb)) then
      message(1) = "Electric field cannot be applied to a periodic system (currently)."
      call write_fatal(1)
    endif

    if(ierr .ne. 0) then
      message(1) = "Could not read KS orbitals from '"//trim(restart_dir)//GS_DIR//"'"
      message(2) = "Please run a ground-state calculation first and/or"
      message(3) = "give the correct RestartDataset in the input file."
      call write_fatal(3)
    end if

    ! set up Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup (sys, hm)

    call io_mkdir(trim(tmpdir)//EM_RESP_FD_DIR) ! restart

    ! Allocate the dipole
    SAFE_ALLOCATE(dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 1:2))
    dipole = M_ZERO

    if(.not. fromScratch) then
      iunit = io_open(trim(tmpdir)//EM_RESP_FD_DIR//RESTART_FILE, action='read', status='old', die=.false.)

      if(iunit > 0) then
        ! Finds out how many dipoles have already been written.
        rewind(iunit)

        read(iunit, fmt=*, iostat = ios) (center_dipole(jj), jj = 1, gr%mesh%sb%dim)
        center_written = (ios .eq. 0)

        i_start = 1
        do ii = 1, 3
          read(iunit, fmt=*, iostat = ios) ((dipole(ii, jj, isign), jj = 1, gr%mesh%sb%dim), isign = 1, 2)
          if(ios .ne. 0) exit
          i_start = i_start + 1
        end do

        read(iunit, fmt=*, iostat = ios) (diag_dipole(jj), jj = 1, gr%mesh%sb%dim)
        diagonal_done = (ios .eq. 0)

        call io_close(iunit)
      else
        fromScratch = .true.
      end if
    end if

    if(iand(sys%outp%what, output_density) .ne. 0 .or. &
       iand(sys%outp%what, output_pol_density) .ne. 0 .or. calc_Born) then
       fromScratch = .true.
       diagonal_done = .false.
       ! since only dipoles, not densities or forces, are stored and reloaded, we need to
       ! recalculate to get the densities and forces again
    endif

    if(fromScratch) then
      if(mpi_grp_is_root(mpi_world)) then
        ! open new file, and erase old data
        iunit = io_open(trim(tmpdir)//EM_RESP_FD_DIR//RESTART_FILE, action='write', status='replace')
        call io_close(iunit)
      end if
      center_written = .false.
      i_start = 1
    end if

    ! Save local potential
    SAFE_ALLOCATE(Vpsl_save(1:gr%mesh%np))
    Vpsl_save = hm%ep%Vpsl

    ! Allocate the trrho to contain the trace of the density.
    SAFE_ALLOCATE(trrho(1:gr%mesh%np))
    SAFE_ALLOCATE(gs_rho(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(tmp_rho(1:gr%mesh%np, 1:st%d%nspin))
    trrho = M_ZERO
    gs_rho = M_ZERO

    call output_init_()
    call scf_init(scfv, gr, sys%geo, st, hm)
    call Born_charges_init(Born_charges, sys%geo, st, gr%mesh%sb%dim)

    ! now calculate the dipole without field

    hm%ep%vpsl(1:gr%mesh%np) = vpsl_save(1:gr%mesh%np)
    call hamiltonian_update_potential(hm, gr%mesh)

    ! set barrier before the first communication takes place
    ! this ensures proper debug timing of MPI calls
#if defined(HAVE_MPI)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif

    write(message(1), '(a)')
    write(message(2), '(a)') 'Info: Calculating dipole moment for zero field.'
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

    ! Writes the dipole to file
    if(mpi_grp_is_root(mpi_world) .and. .not. center_written) then 
      iunit = io_open(trim(tmpdir)//EM_RESP_FD_DIR//RESTART_FILE, action='write', status='old', position='append')
      write(iunit, fmt='(6e20.12)') (center_dipole(jj), jj = 1, gr%mesh%sb%dim)
      call io_close(iunit)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call geometry_dipole(sys%geo, ionic_dipole)
      print_dipole(1:gr%mesh%sb%dim) = center_dipole(1:gr%mesh%sb%dim) + ionic_dipole(1:gr%mesh%sb%dim)
      call io_output_dipole(stdout, print_dipole, gr%mesh%sb%dim)
    endif

    do ii = i_start, gr%mesh%sb%dim
      do isign = 1, 2
        write(message(1), '(a)')
        write(message(2), '(a,f6.4,5a)') 'Info: Calculating dipole moment for field ', &
          units_from_atomic(units_out%force, -(-1)**isign * e_field), ' ', &
          trim(units_abbrev(units_out%force)), ' in the ', index2axis(ii), '-direction.'
        call write_info(2)
        ! there is an extra factor of -1 in here that is for the electronic charge

        hm%ep%vpsl(1:gr%mesh%np) = vpsl_save(1:gr%mesh%np) + (-1)**isign * gr%mesh%x(1:gr%mesh%np, ii) * e_field
        call hamiltonian_update_potential(hm, gr%mesh)

        if(start_density_is_zero_field) then
          st%rho(1:gr%mesh%np, 1:st%d%nspin) = gs_rho(1:gr%mesh%np, 1:st%d%nspin)
        else
          call lcao_run(sys, hm)
        endif

        call scf_mix_clear(scfv)
        call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)

        trrho = M_ZERO
        do is = 1, st%d%spin_channels
          trrho(1:gr%mesh%np) = trrho(1:gr%mesh%np) + st%rho(1:gr%mesh%np, is)
        end do

        ! calculate dipole
        do jj = 1, gr%mesh%sb%dim
          dipole(ii, jj, isign) = dmf_moment(gr%mesh, trrho, jj, 1)
        end do

        if(mpi_grp_is_root(mpi_world)) then
          print_dipole(1:gr%mesh%sb%dim) = dipole(ii, 1:gr%mesh%sb%dim, isign) + ionic_dipole(1:gr%mesh%sb%dim)
          call io_output_dipole(stdout, print_dipole, gr%mesh%sb%dim)
        endif

        call output_cycle_()

      end do

      ! Writes the dipole to file
      if(mpi_grp_is_root(mpi_world)) then 
        iunit = io_open(trim(tmpdir)//EM_RESP_FD_DIR//RESTART_FILE, action='write', status='old', position='append')
        write(iunit, fmt='(6e20.12)') ((dipole(ii, jj, isign), jj = 1, gr%mesh%sb%dim), isign = 1, 2)
        call io_close(iunit)
      end if
    end do
    
    if(.not. diagonal_done .and. calc_diagonal) then
      write(message(1), '(a)')
      write(message(2), '(a,f6.4,3a, f6.4, 3a)') 'Info: Calculating dipole moment for field ', &
         units_from_atomic(units_out%force, e_field), ' ', &
         trim(units_abbrev(units_out%force)), ' in the '//index2axis(2)//'-direction plus ', &
         units_from_atomic(units_out%force, e_field), ' ', &
         trim(units_abbrev(units_out%force)), ' in the '//index2axis(3)//'-direction.'
      call write_info(2)
  
      hm%ep%vpsl(1:gr%mesh%np) = vpsl_save(1:gr%mesh%np) &
        - (gr%mesh%x(1:gr%mesh%np, 2) + gr%mesh%x(1:gr%mesh%np, 3)) * e_field
      call hamiltonian_update_potential(hm, gr%mesh)
  
      if(start_density_is_zero_field) then
        st%rho(1:gr%mesh%np, 1:st%d%nspin) = gs_rho(1:gr%mesh%np, 1:st%d%nspin)
      else
        call lcao_run(sys, hm)
      endif

      call scf_mix_clear(scfv)
      call scf_run(scfv, sys%gr, sys%geo, st, sys%ks, hm, sys%outp, gs_run=.false., verbosity = VERB_COMPACT)
  
      trrho = M_ZERO
      do is = 1, st%d%spin_channels
        trrho(1:gr%mesh%np) = trrho(1:gr%mesh%np) + st%rho(1:gr%mesh%np, is)
      end do
  
      ! calculate dipole
      do jj = 1, gr%mesh%sb%dim
        diag_dipole(jj) = dmf_moment(gr%mesh, trrho, jj, 1)
      end do

      if(mpi_grp_is_root(mpi_world)) then
        print_dipole(1:gr%mesh%sb%dim) = diag_dipole(1:gr%mesh%sb%dim) + ionic_dipole(1:gr%mesh%sb%dim)
        call io_output_dipole(stdout, print_dipole, gr%mesh%sb%dim)
      endif
  
      ! Writes the dipole to file
      if(mpi_grp_is_root(mpi_world)) then 
        iunit = io_open(trim(tmpdir)//EM_RESP_FD_DIR//RESTART_FILE, action='write', status='old', position='append')
        write(iunit, fmt='(3e20.12)') (diag_dipole(jj), jj = 1, gr%mesh%sb%dim)
        call io_close(iunit)
      end if
    endif

    call scf_end(scfv)
    call output_end_()
    call Born_charges_end(Born_charges)

    SAFE_DEALLOCATE_A(Vpsl_save)
    SAFE_DEALLOCATE_A(trrho)
    SAFE_DEALLOCATE_A(gs_rho)
    SAFE_DEALLOCATE_A(tmp_rho)
    SAFE_DEALLOCATE_A(dipole)
    call end_()
    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol.static_pol_run.init_')

      ! shortcuts
      gr => sys%gr
      st => sys%st

      call states_allocate_wfns(st, gr%mesh)

      call messages_obsolete_variable("EMStaticField", "EMStaticElectricField")
      !%Variable EMStaticElectricField
      !%Type float
      !%Default 0.01 a.u.
      !%Section Linear Response::Static Polarization
      !%Description
      !% Magnitude of the static electric field used to calculate the static polarizability,
      !% if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_float(datasets_check('EMStaticElectricField'), &
         CNST(0.01), e_field, units_inp%force)
      if (e_field <= M_ZERO) then
        write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid EMElectricStaticField."
        message(2) = '(0 < EMStaticElectricField)'
        call write_fatal(2)
      end if

      ! variable defined in em_resp
      call parse_logical(datasets_check('EMCalcBornCharges'), .false., calc_Born)
      if (calc_Born) call messages_devel_version("Calculation of Born effective charges")

      !%Variable EMStartDensityIsZeroField
      !%Type logical
      !%Default true
      !%Section Linear Response::Static Polarization
      !%Description
      !% Use the charge density from the zero-field calculation as the starting density for
      !% SCF calculations with applied fields. For small fields, this will be fastest.
      !% If there is trouble converging with larger fields, set to false,
      !% to initialize the calculation for each field from scratch, as specified by the LCAO variables. 
      !% Only applies if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_logical(datasets_check('EMStartDensityIsZeroField'), .true., start_density_is_zero_field)

      !%Variable EMCalcDiagonalField
      !%Type logical
      !%Default true
      !%Section Linear Response::Static Polarization
      !%Description
      !% Calculate <i>yz</i>-field for beta_<i>xyz</i> hyperpolarizability, which is sometimes harder to converge.
      !% Only applies if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_logical(datasets_check('EMCalcDiagonalField'), .true., calc_diagonal)

      call pop_sub()
    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()

      call push_sub('static_pol.end_')

      call states_deallocate_wfns(sys%st)

      call pop_sub()
    end subroutine end_
  

    !-------------------------------------------------------------
    subroutine output_init_()
      call push_sub('static_pol.output_init_')

      !allocate memory for what we want to output
      if(iand(sys%outp%what, output_density) .ne. 0 .or. &
         iand(sys%outp%what, output_pol_density) .ne. 0 ) then 
        SAFE_ALLOCATE(lr_rho (1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(lr_rho2(1:gr%mesh%np, 1:st%d%nspin))
      end if
      
      if(iand(sys%outp%what, output_elf) .ne. 0) then 
        SAFE_ALLOCATE(    elf(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE( lr_elf(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(   elfd(1:gr%mesh%np, 1:st%d%nspin))
        SAFE_ALLOCATE(lr_elfd(1:gr%mesh%np, 1:st%d%nspin))
      end if
      
      call pop_sub()
    end subroutine output_init_


    !-------------------------------------------------------------
    subroutine output_cycle_()
      integer iatom
      type(unit_t) fn_unit
      
      call push_sub('static_pol.output_cycle_')

      ! BORN CHARGES
      if(calc_Born) then
        do iatom = 1, sys%geo%natoms
          if(isign == 1) then 
          ! temporary assignment for use in next cycle when isign == 2
            Born_charges%charge(ii, 1:gr%mesh%sb%dim, iatom) = sys%geo%atom(iatom)%f(1:gr%mesh%sb%dim)
          else
            Born_charges%charge(ii, 1:gr%mesh%sb%dim, iatom) = &
              (sys%geo%atom(iatom)%f(1:gr%mesh%sb%dim) - Born_charges%charge(ii, 1:gr%mesh%sb%dim, iatom)) / (M_TWO*e_field)
            Born_charges%charge(ii, ii, iatom) = Born_charges%charge(ii, ii, iatom) + species_zval(sys%geo%atom(iatom)%spec)
            ! since the efield is applied in the SCF calculation by just altering the external potential felt by the electrons,
            ! the ionic force due to the efield is not included in the forces returned by the SCF run, and so the ionic
            ! contribution to the Born charge must be added by hand here
          endif
        enddo
      endif

      !DENSITY AND POLARIZABILITY DENSITY   
      if(iand(sys%outp%what, output_density) .ne. 0 .or. &
         iand(sys%outp%what, output_pol_density) .ne. 0) then 
         
        if(isign == 1 .and. ii == 2) then
          tmp_rho(1:gr%mesh%np, 1:st%d%nspin) = st%rho(1:gr%mesh%np, 1:st%d%nspin)
          ! for use in off-diagonal non-linear densities
        endif

        if(isign == 1) then 
          ! temporary assignment for use in next cycle when isign == 2
          lr_rho(1:gr%mesh%np, 1:st%d%nspin) = st%rho(1:gr%mesh%np, 1:st%d%nspin)
        else
          lr_rho2(1:gr%mesh%np, 1:st%d%nspin) = &
            -(st%rho(1:gr%mesh%np, 1:st%d%nspin) + lr_rho(1:gr%mesh%np, 1:st%d%nspin) - &
              2 * gs_rho(1:gr%mesh%np, 1:st%d%nspin)) / e_field**2

          lr_rho(1:gr%mesh%np, 1:st%d%nspin) = &
               (st%rho(1:gr%mesh%np, 1:st%d%nspin) - lr_rho(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, st%d%nspin
            if(iand(sys%outp%what, output_density) .ne. 0) then
              fn_unit = units_out%length**(1-gr%sb%dim) / units_out%energy
              write(fname, '(a,i1,2a)') 'fd_density-sp', is, '-', index2axis(ii)
              call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
                gr%mesh, lr_rho(:, is), fn_unit, ierr, geo = sys%geo)

              ! save the trouble of writing many copies of each density, since ii,jj = jj,ii
              fn_unit = units_out%length**(2-gr%sb%dim) / units_out%energy**2
              do jj = ii, gr%mesh%sb%dim
                write(fname, '(a,i1,4a)') 'fd2_density-sp', is, '-', index2axis(ii), '-', index2axis(jj)
                call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
                  gr%mesh, lr_rho2(:, is), fn_unit, ierr, geo = sys%geo)
              enddo
            endif

            if(iand(sys%outp%what, output_pol_density) .ne. 0) then
              do jj = ii, gr%mesh%sb%dim
                fn_unit = units_out%length**(2-gr%sb%dim) / units_out%energy
                write(fname, '(a,i1,4a)') 'alpha_density-sp', is, '-', index2axis(ii), '-', index2axis(jj)
                call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname), &
                  gr%mesh, -gr%mesh%x(:, jj) * lr_rho(:, is), fn_unit, ierr, geo = sys%geo)

                fn_unit = units_out%length**(3-gr%sb%dim) / units_out%energy**2
                write(fname, '(a,i1,6a)') 'beta_density-sp', is, '-', index2axis(ii), &
                  '-', index2axis(ii), '-', index2axis(jj)
                call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname), &
                  gr%mesh, -gr%mesh%x(:, jj) * lr_rho2(:, is), fn_unit, ierr, geo = sys%geo)
              enddo
            endif
          end do

        end if
      end if

      !ELF
      if(iand(sys%outp%what, output_elf) .ne. 0) then 
         
        if(isign == 1) then 
          call elf_calc(st, gr, elf, elfd)
        else
          call elf_calc(st, gr, lr_elf, lr_elfd)
          
          !numerical derivative
          lr_elf(1:gr%mesh%np, 1:st%d%nspin) = &
               ( lr_elf(1:gr%mesh%np, 1:st%d%nspin) -  elf(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO * e_field)
          lr_elfd(1:gr%mesh%np, 1:st%d%nspin) = &
               (lr_elfd(1:gr%mesh%np, 1:st%d%nspin) - elfd(1:gr%mesh%np, 1:st%d%nspin)) / (M_TWO * e_field)

          !write
          do is = 1, st%d%nspin
            write(fname, '(a,i1,2a)') 'lr_elf-sp', is, '-', index2axis(ii)
            call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
                gr%mesh, lr_elf(:, is), unit_one, ierr, geo = sys%geo)
            write(fname, '(a,i1,2a)') 'lr_elf_D-sp', is, '-', index2axis(ii)
            call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
                gr%mesh, lr_elfd(:, is), unit_one, ierr, geo = sys%geo)
          end do
        end if

      end if

      call pop_sub()
    end subroutine output_cycle_


    !-------------------------------------------------------------
    subroutine output_end_()
      FLOAT :: alpha(MAX_DIM, MAX_DIM)
      CMPLX :: beta(MAX_DIM, MAX_DIM, MAX_DIM)
      integer :: iunit, idir
      type(unit_t) fn_unit

      call push_sub('static_pol.output_end_')

      call io_mkdir(EM_RESP_FD_DIR)

      if((iand(sys%outp%what, output_density) .ne. 0 .or. &
         iand(sys%outp%what, output_pol_density) .ne. 0) .and. calc_diagonal) then 
        lr_rho2(1:gr%mesh%np, 1:st%d%nspin) = &
          -(st%rho(1:gr%mesh%np, 1:st%d%nspin) - lr_rho(1:gr%mesh%np, 1:st%d%nspin) &
          - tmp_rho(1:gr%mesh%np, 1:st%d%nspin) + gs_rho(1:gr%mesh%np, 1:st%d%nspin)) / e_field**2
  
        do is = 1, st%d%nspin
          if(iand(sys%outp%what, output_density) .ne. 0) then
            fn_unit = units_out%length**(2-gr%sb%dim) / units_out%energy**2
            write(fname, '(a,i1,a)') 'fd2_density-sp', is, '-y-z'
            call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
              gr%mesh, lr_rho2(:, is), unit_one, ierr, geo = sys%geo)
          endif
  
          if(iand(sys%outp%what, output_pol_density) .ne. 0) then
            fn_unit = units_out%length**(3-gr%sb%dim) / units_out%energy**2
            write(fname, '(a,i1,a)') 'beta_density-sp', is, '-x-y-z'
            call doutput_function(sys%outp%how, EM_RESP_FD_DIR, trim(fname),&
              gr%mesh, -gr%mesh%x(:, 1) * lr_rho2(:, is), unit_one, ierr, geo = sys%geo)
          endif
        end do
      endif

      if(mpi_grp_is_root(mpi_world)) then ! output pol file
        iunit = io_open(EM_RESP_FD_DIR//'alpha', action='write')
        write(iunit, '(3a)') '# Polarizability tensor [', trim(units_abbrev(units_out%polarizability)), ']'

        alpha(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim) = (dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 1) - &
             dipole(1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim, 2)) / (M_TWO * e_field)

        beta = M_ZERO

        do idir = 1, gr%mesh%sb%dim
          beta(1:gr%mesh%sb%dim, idir, idir) = -(dipole(idir, 1:gr%mesh%sb%dim, 1) + dipole(idir, 1:gr%mesh%sb%dim, 2) - &
            M_TWO * center_dipole(1:gr%mesh%sb%dim)) / e_field**2
          beta(idir, 1:gr%mesh%sb%dim, idir) = beta(1:gr%mesh%sb%dim, idir, idir) 
          beta(idir, idir, 1:gr%mesh%sb%dim) = beta(1:gr%mesh%sb%dim, idir, idir)
        end do

        if(calc_diagonal) then
          beta(1, 2, 3) = -(diag_dipole(1) - dipole(2, 1, 1) - dipole(3, 1, 1) + center_dipole(1)) / e_field**2
        else
          beta(1, 2, 3) = M_ZERO
        endif

        beta(2, 3, 1) = beta(1, 2, 3)
        beta(3, 1, 2) = beta(1, 2, 3)
        beta(3, 2, 1) = beta(1, 2, 3)
        beta(1, 3, 2) = beta(1, 2, 3)
        beta(2, 1, 3) = beta(1, 2, 3)

        call io_output_tensor(iunit, alpha, gr%mesh%sb%dim, units_out%polarizability)
        call io_close(iunit)
        
        call out_hyperpolarizability(gr%sb, beta, .true., EM_RESP_FD_DIR)

        if(calc_Born) call out_Born_charges(Born_charges, sys%geo, gr%mesh%sb%dim, EM_RESP_FD_DIR)
      end if

      if(iand(sys%outp%what, output_density) .ne. 0 .or. &
         iand(sys%outp%what, output_pol_density) .ne. 0) then 
        SAFE_DEALLOCATE_A(lr_rho)
      end if
      
      if(iand(sys%outp%what, output_elf) .ne. 0) then 
        SAFE_DEALLOCATE_A(lr_elf)
        SAFE_DEALLOCATE_A(elf)
        SAFE_DEALLOCATE_A(lr_elfd)
        SAFE_DEALLOCATE_A(elfd)
      end if

      call pop_sub()
    end subroutine output_end_

  end subroutine static_pol_run

end module static_pol_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
