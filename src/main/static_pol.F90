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

#include "global.h"
#define RESTART_FILE 'dipoles'

module static_pol_oct_m
  use born_charges_oct_m
  use elf_oct_m
  use em_resp_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lcao_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m

  implicit none

  private
  public :: &
    static_pol_run


contains

  ! ---------------------------------------------------------
  subroutine static_pol_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(static_pol_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = static_pol not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call static_pol_run_legacy(system, from_scratch)
    end select

    POP_SUB(static_pol_run)
  end subroutine static_pol_run

  ! ---------------------------------------------------------
  subroutine static_pol_run_legacy(sys, fromScratch)
    type(electrons_t),    intent(inout) :: sys
    logical,              intent(in)    :: fromScratch

    type(scf_t) :: scfv
    integer :: iunit, ios, i_start, ii, jj, is, isign, ierr, read_count, verbosity
    FLOAT :: e_field, e_field_saved
    FLOAT, allocatable :: Vpsl_save(:), trrho(:), dipole(:, :, :)
    FLOAT, allocatable :: elf(:,:), lr_elf(:,:), elfd(:,:), lr_elfd(:,:)
    FLOAT, allocatable :: lr_rho(:,:), lr_rho2(:,:), gs_rho(:,:), tmp_rho(:,:)
    FLOAT :: center_dipole(1:sys%space%dim), diag_dipole(1:sys%space%dim), ionic_dipole(1:sys%space%dim), &
      print_dipole(1:sys%space%dim)
    type(born_charges_t) :: born_charges
    logical :: calc_Born, start_density_is_zero_field, write_restart_densities, calc_diagonal, verbose
    logical :: diagonal_done, center_written, fromScratch_local, field_written
    character(len=MAX_PATH_LEN) :: fname, dir_name
    character(len=120) :: line(1)
    character :: sign_char
    type(restart_t) :: gs_restart, restart_load, restart_dump

    PUSH_SUB(static_pol_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%kpoints%use_symmetries) call messages_experimental("KPoints symmetries with CalculationMode = em_resp")

    call init_()

    ! load wavefunctions
    call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_load(gs_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
    end if
    if (ierr /= 0) then
      message(1) = "Unable to read wavefunctions."
      call messages_fatal(1)
    end if
    call restart_end(gs_restart)

    if(sys%space%is_periodic()) then
      message(1) = "Electric field cannot be applied to a periodic system (currently)."
      call messages_fatal(1)
    end if

    ! set up Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call messages_info(1)
    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, calc_eigenval = .false.) ! we read them from restart

    ! Allocate the dipole
    SAFE_ALLOCATE(dipole(1:sys%space%dim, 1:sys%space%dim, 1:2))
    dipole = M_ZERO

    i_start = 1
    center_written = .false.
    diagonal_done = .false.
    field_written = .false.

    call restart_init(restart_dump, sys%namespace, RESTART_EM_RESP_FD, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)

    if(.not. fromScratch) then
      call restart_init(restart_load, sys%namespace, RESTART_EM_RESP_FD, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)
      if(ierr == 0) then
        iunit = restart_open(restart_load, RESTART_FILE)
      else
        iunit = 0
      end if

      if(iunit > 0) then
        ! Finds out how many dipoles have already been written.
        rewind(iunit)

        read(iunit, fmt=*, iostat = ios) e_field_saved
        field_written = (ios  ==  0)

        read(iunit, fmt=*, iostat = ios) (center_dipole(jj), jj = 1, sys%space%dim)
        center_written = (ios  ==  0)

        do ii = 1, 3
          read(iunit, fmt=*, iostat = ios) ((dipole(ii, jj, isign), jj = 1, sys%space%dim), isign = 1, 2)
          if(ios /= 0) exit
          i_start = i_start + 1
        end do

        read(iunit, fmt=*, iostat = ios) (diag_dipole(jj), jj = 1, sys%space%dim)
        diagonal_done = (ios  ==  0)

        call restart_close(restart_load, iunit)
      end if

      ! if saved dipoles used a different e_field, we cannot use them
      if(field_written .and. abs(e_field_saved - e_field) > CNST(1e-15)) then
        message(1) = "Saved dipoles are from a different electric field, cannot use them."
        call messages_warning(1)
        center_written = .false.
        diagonal_done = .false.
        i_start = 1
      end if

      read_count = (i_start - 1) * 2
      if(center_written) read_count = read_count + 1
      if(diagonal_done)  read_count = read_count + 1
      write(message(1),'(a,i1,a)') "Using ", read_count, " dipole(s) from file."
      call messages_info(1)
    end if

    if (sys%outp%what(OPTION__OUTPUT__DENSITY) .or. &
       sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then
       if (i_start > 2 .and. calc_diagonal) then
          i_start = 2
          diagonal_done = .false.
          !FIXME: take derivatives between yz and z (not y) so can restart from only last (z) calc
       end if
    end if

    if(i_start  ==  1) then
      ! open new file, erase old data, write e_field
      iunit = restart_open(restart_dump, RESTART_FILE, status='replace')
      write(line(1), fmt='(e20.12)') e_field
      call restart_write(restart_dump, iunit, line, 1, ierr)
      if (ierr /= 0) then
        message(1) = "Unsuccessful write of electric field."
        call messages_warning(1)
      end if
      call restart_close(restart_dump, iunit)
      center_written = .false.
    end if

    ! Save local potential
    SAFE_ALLOCATE(Vpsl_save(1:sys%gr%mesh%np))
    Vpsl_save = sys%hm%ep%Vpsl

    ! Allocate the trrho to contain the trace of the density.
    SAFE_ALLOCATE(trrho(1:sys%gr%mesh%np))
    SAFE_ALLOCATE(gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
    SAFE_ALLOCATE(tmp_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
    trrho = M_ZERO
    gs_rho = M_ZERO

    call output_init_()
    call scf_init(scfv, sys%namespace, sys%gr, sys%ions, sys%st, sys%mc, sys%hm, sys%ks, sys%space)
    call born_charges_init(Born_charges, sys%namespace, sys%ions, sys%st, sys%space%dim)

    ! now calculate the dipole without field

    sys%hm%ep%vpsl(1:sys%gr%mesh%np) = vpsl_save(1:sys%gr%mesh%np)
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)

    write(message(1), '(a)')
    write(message(2), '(a)') 'Info: Calculating dipole moment for zero field.'
    call messages_info(2)
    call scf_run(scfv, sys%namespace, sys%space, sys%mc, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, sys%outp, &
      gs_run=.false., verbosity = verbosity)

    gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)
    trrho = M_ZERO
    do is = 1, sys%st%d%spin_channels
      trrho(1:sys%gr%mesh%np) = trrho(1:sys%gr%mesh%np) + gs_rho(1:sys%gr%mesh%np, is)
    end do

    ! calculate dipole
    do jj = 1, sys%space%dim
       center_dipole(jj) = dmf_moment(sys%gr%mesh, trrho, jj, 1)
    end do

    ! Writes the dipole to file
    if(.not. center_written) then 
      iunit = restart_open(restart_dump, RESTART_FILE, position='append')
      write(line(1), fmt='(6e20.12)') (center_dipole(jj), jj = 1, sys%space%dim)
      call restart_write(restart_dump, iunit, line, 1, ierr)
      if (ierr /= 0) then
        message(1) = "Unsuccessful write of center dipole."
        call messages_warning(1)
      end if
      call restart_close(restart_dump, iunit)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      ionic_dipole(1:sys%space%dim) = sys%ions%dipole()
      print_dipole(1:sys%space%dim) = center_dipole(1:sys%space%dim) + ionic_dipole(1:sys%space%dim)
      call output_dipole(stdout, print_dipole, sys%space%dim)
    end if

    do ii = i_start, sys%space%dim
      do isign = 1, 2
        write(message(1), '(a)')
        write(message(2), '(a,f6.4,5a)') 'Info: Calculating dipole moment for field ', &
          units_from_atomic(units_out%force, (-1)**isign * e_field), ' ', &
          trim(units_abbrev(units_out%force)), ' in the ', index2axis(ii), '-direction.'
        call messages_info(2)
        ! there would be an extra factor of -1 in here that is for the electronic charge
        ! except that we treat electrons as positive

        sys%hm%ep%vpsl(1:sys%gr%mesh%np) = vpsl_save(1:sys%gr%mesh%np) + (-1)**isign * sys%gr%mesh%x(1:sys%gr%mesh%np, ii) * e_field
        call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)

        if(isign == 1) then
          sign_char = '+'
        else
          sign_char = '-'
        end if

        write(dir_name,'(a)') "field_"//index2axis(ii)//sign_char
        fromScratch_local = fromScratch

        if(.not. fromScratch) then
          call restart_open_dir(restart_load, trim(dir_name), ierr)
          if (ierr == 0) then
            call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
          end if
          call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
          if(ierr /= 0) fromScratch_local = .true.
          call restart_close_dir(restart_load)
        end if

        if(fromScratch_local) then
          if(start_density_is_zero_field) then
            sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)
            call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
          else
            call lcao_run(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, lmm_r = scfv%lmm_r)
          end if
        end if

        call scf_mix_clear(scfv)
        call scf_run(scfv, sys%namespace, sys%space, sys%mc, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, sys%outp, &
          gs_run=.false., verbosity = verbosity)

        trrho = M_ZERO
        do is = 1, sys%st%d%spin_channels
          trrho(1:sys%gr%mesh%np) = trrho(1:sys%gr%mesh%np) + sys%st%rho(1:sys%gr%mesh%np, is)
        end do

        ! calculate dipole
        do jj = 1, sys%space%dim
          dipole(ii, jj, isign) = dmf_moment(sys%gr%mesh, trrho, jj, 1)
        end do

        if(mpi_grp_is_root(mpi_world)) then
          print_dipole(1:sys%space%dim) = dipole(ii, 1:sys%space%dim, isign) + ionic_dipole(1:sys%space%dim)
          call output_dipole(stdout, print_dipole, sys%space%dim)
        end if

        call output_cycle_()

        if(write_restart_densities) then
          call restart_open_dir(restart_dump, trim(dir_name), ierr)
          if (ierr == 0) then
            call states_elec_dump(restart_dump, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
          end if
          call restart_close_dir(restart_dump)
          if(ierr /= 0) then
            message(1) = 'Unable to write states wavefunctions.'
            call messages_warning(1)
          end if
        end if
      end do

      ! Writes the dipole to file
      iunit = restart_open(restart_dump, RESTART_FILE, position='append')
      write(line(1), '(6e20.12)') ((dipole(ii, jj, isign), jj = 1, sys%space%dim), isign = 1, 2)
      call restart_write(restart_dump, iunit, line, 1, ierr)
      if (ierr /= 0) then
        message(1) = "Unsuccessful write of dipole."
        call messages_warning(1)
      end if
      call restart_close(restart_dump, iunit)
    end do
    
    if(.not. diagonal_done .and. calc_diagonal) then
      write(message(1), '(a)')
      write(message(2), '(a,f6.4,3a, f6.4, 3a)') 'Info: Calculating dipole moment for field ', &
         units_from_atomic(units_out%force, e_field), ' ', &
         trim(units_abbrev(units_out%force)), ' in the '//index2axis(2)//'-direction plus ', &
         units_from_atomic(units_out%force, e_field), ' ', &
         trim(units_abbrev(units_out%force)), ' in the '//index2axis(3)//'-direction.'
      call messages_info(2)
  
      sys%hm%ep%vpsl(1:sys%gr%mesh%np) = vpsl_save(1:sys%gr%mesh%np) &
        - (sys%gr%mesh%x(1:sys%gr%mesh%np, 2) + sys%gr%mesh%x(1:sys%gr%mesh%np, 3)) * e_field
      call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)
  
      if(isign == 1) then
        sign_char = '+'
      else
        sign_char = '-'
      end if

      fromScratch_local = fromScratch

      if(.not. fromScratch) then
        call restart_open_dir(restart_load, "field_yz+", ierr)
        if (ierr == 0) then
          call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
        end if
        call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
        if(ierr /= 0) fromScratch_local = .true.
        call restart_close_dir(restart_load)
      end if

      if(fromScratch_local) then
        if(start_density_is_zero_field) then
          sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)
          call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
        else
          call lcao_run(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, lmm_r = scfv%lmm_r)
        end if
      end if

      call scf_mix_clear(scfv)
      call scf_run(scfv, sys%namespace, sys%space, sys%mc, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, sys%outp, &
        gs_run=.false., verbosity = verbosity)
  
      trrho = M_ZERO
      do is = 1, sys%st%d%spin_channels
        trrho(1:sys%gr%mesh%np) = trrho(1:sys%gr%mesh%np) + sys%st%rho(1:sys%gr%mesh%np, is)
      end do
  
      ! calculate dipole
      do jj = 1, sys%space%dim
        diag_dipole(jj) = dmf_moment(sys%gr%mesh, trrho, jj, 1)
      end do

      if (mpi_grp_is_root(mpi_world)) then
        print_dipole(1:sys%space%dim) = diag_dipole(1:sys%space%dim) + ionic_dipole(1:sys%space%dim)
        call output_dipole(stdout, print_dipole, sys%space%dim)
      end if
  
      ! Writes the dipole to file
      iunit = restart_open(restart_dump, RESTART_FILE, position='append')
      write(line(1), fmt='(3e20.12)') (diag_dipole(jj), jj = 1, sys%space%dim)
      call restart_write(restart_dump, iunit, line, 1, ierr)
      if (ierr /= 0) then
        message(1) = "Unsuccessful write of dipole."
        call messages_warning(1)
      end if
      call restart_close(restart_dump, iunit)

      if(write_restart_densities) then
        call restart_open_dir(restart_dump, "field_yz+", ierr)
        if (ierr == 0) then
          call states_elec_dump(restart_dump, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
        end if
        call restart_close_dir(restart_dump)
        if(ierr /= 0) then
          message(1) = 'Unable to write states wavefunctions.'
          call messages_warning(1)
        end if
      end if

    end if

    if(.not. fromScratch) call restart_end(restart_load)
    call scf_end(scfv)
    call output_end_()
    call Born_charges_end(Born_charges)

    SAFE_DEALLOCATE_A(Vpsl_save)
    SAFE_DEALLOCATE_A(trrho)
    SAFE_DEALLOCATE_A(gs_rho)
    SAFE_DEALLOCATE_A(tmp_rho)
    SAFE_DEALLOCATE_A(dipole)
    call end_()
    POP_SUB(static_pol_run_legacy)

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      PUSH_SUB(static_pol_run_legacy.init_)

      call states_elec_allocate_wfns(sys%st, sys%gr%mesh)

      call messages_obsolete_variable(sys%namespace, "EMStaticField", "EMStaticElectricField")
      !%Variable EMStaticElectricField
      !%Type float
      !%Default 0.01 a.u.
      !%Section Linear Response::Static Polarization
      !%Description
      !% Magnitude of the static electric field used to calculate the static polarizability,
      !% if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_variable(sys%namespace, 'EMStaticElectricField', CNST(0.01), e_field, units_inp%force)
      if (e_field <= M_ZERO) then
        write(message(1), '(a,e14.6,a)') "Input: '", e_field, "' is not a valid EMStaticElectricField."
        message(2) = '(Must have EMStaticElectricField > 0)'
        call messages_fatal(2)
      end if

      ! variable defined in em_resp
      call parse_variable(sys%namespace, 'EMCalcBornCharges', .false., calc_Born)
      if (calc_Born) call messages_experimental("Calculation of Born effective charges")

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
      call parse_variable(sys%namespace, 'EMStartDensityIsZeroField', .true., start_density_is_zero_field)

      !%Variable EMCalcDiagonalField
      !%Type logical
      !%Default true
      !%Section Linear Response::Static Polarization
      !%Description
      !% Calculate <i>yz</i>-field for <math>\beta_{xyz}</math> hyperpolarizability, which is sometimes harder to converge.
      !% Only applies if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_variable(sys%namespace, 'EMCalcDiagonalField', .true., calc_diagonal)

      !%Variable EMWriteRestartDensities
      !%Type logical
      !%Default true
      !%Section Linear Response::Static Polarization
      !%Description
      !% Write density after each calculation for restart, rather than just the resulting electronic dipole moment.
      !% Only applies if <tt>ResponseMethod = finite_differences</tt>. Restarting from calculations at smaller
      !% fields can be helpful if there are convergence problems.
      !%End
      call parse_variable(sys%namespace, 'EMWriteRestartDensities', .true., write_restart_densities)

      !%Variable EMVerbose
      !%Type logical
      !%Default false
      !%Section Linear Response::Static Polarization
      !%Description
      !% Write full SCF output.
      !% Only applies if <tt>ResponseMethod = finite_differences</tt>.
      !%End
      call parse_variable(sys%namespace, 'EMVerbose', .false., verbose)

      if(verbose) then
        verbosity = VERB_FULL
      else
        verbosity = VERB_COMPACT
      end if

      POP_SUB(static_pol_run_legacy.init_)
    end subroutine init_

    ! ---------------------------------------------------------
    subroutine end_()

      PUSH_SUB(end_)

      call states_elec_deallocate_wfns(sys%st)

      POP_SUB(end_)
    end subroutine end_
  

    !-------------------------------------------------------------
    subroutine output_init_()
      PUSH_SUB(output_init_)

      !allocate memory for what we want to output
      if (sys%outp%what(OPTION__OUTPUT__DENSITY) .or. &
         sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then 
        SAFE_ALLOCATE(lr_rho (1:sys%gr%mesh%np, 1:sys%st%d%nspin))
        SAFE_ALLOCATE(lr_rho2(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
      end if
      
      if (sys%outp%what(OPTION__OUTPUT__ELF)) then 
        SAFE_ALLOCATE(    elf(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
        SAFE_ALLOCATE( lr_elf(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
        SAFE_ALLOCATE(   elfd(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
        SAFE_ALLOCATE(lr_elfd(1:sys%gr%mesh%np, 1:sys%st%d%nspin))
      end if
      
      POP_SUB(output_init_)
    end subroutine output_init_


    !-------------------------------------------------------------
    subroutine output_cycle_()
      integer :: iatom
      type(unit_t) :: fn_unit
      
      PUSH_SUB(output_cycle_)

      ! BORN CHARGES
      if(calc_Born) then
        do iatom = 1, sys%ions%natoms
          if(isign == 1) then
          ! temporary assignment for use in next cycle when isign == 2
            Born_charges%charge(ii, 1:sys%space%dim, iatom) = sys%ions%tot_force(:, iatom)
          else
            Born_charges%charge(ii, 1:sys%space%dim, iatom) = &
              (sys%ions%tot_force(:, iatom) - Born_charges%charge(ii, 1:sys%space%dim, iatom)) &
              / (M_TWO*e_field)
            Born_charges%charge(ii, ii, iatom) = Born_charges%charge(ii, ii, iatom) + species_zval(sys%ions%atom(iatom)%species)
            ! since the efield is applied in the SCF calculation by just altering the external potential felt by the electrons,
            ! the ionic force due to the efield is not included in the forces returned by the SCF run, and so the ionic
            ! contribution to the Born charge must be added by hand here
          end if
        end do
      end if

      !DENSITY AND POLARIZABILITY DENSITY   
      if (sys%outp%what(OPTION__OUTPUT__DENSITY) .or. &
         sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then 
         
        if(isign == 1 .and. ii == 2) then
          tmp_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)
          ! for use in off-diagonal non-linear densities
        end if

        if(isign == 1) then 
          ! temporary assignment for use in next cycle when isign == 2
          lr_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)
        else
          lr_rho2(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = &
            -(sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) + lr_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) - &
              2 * gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)) / e_field**2

          lr_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = &
               (sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) - lr_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)) / (M_TWO*e_field)

          !write
          do is = 1, sys%st%d%nspin
            if (sys%outp%what(OPTION__OUTPUT__DENSITY)) then
              fn_unit = units_out%length**(1-sys%space%dim) / units_out%energy
              write(fname, '(a,i1,2a)') 'fd_density-sp', is, '-', index2axis(ii)
              call dio_function_output(sys%outp%how(OPTION__OUTPUT__DENSITY), EM_RESP_FD_DIR, trim(fname),&
                sys%namespace, sys%space, sys%gr%mesh, lr_rho(:, is), fn_unit, ierr, ions = sys%ions)

              ! save the trouble of writing many copies of each density, since ii,jj = jj,ii
              fn_unit = units_out%length**(2-sys%space%dim) / units_out%energy**2
              do jj = ii, sys%space%dim
                write(fname, '(a,i1,4a)') 'fd2_density-sp', is, '-', index2axis(ii), '-', index2axis(jj)
                call dio_function_output(sys%outp%how(OPTION__OUTPUT__DENSITY), EM_RESP_FD_DIR, trim(fname),&
                  sys%namespace, sys%space, sys%gr%mesh, lr_rho2(:, is), fn_unit, ierr, ions = sys%ions)
              end do
            end if

            if (sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then
              do jj = ii, sys%space%dim
                fn_unit = units_out%length**(2-sys%space%dim) / units_out%energy
                write(fname, '(a,i1,4a)') 'alpha_density-sp', is, '-', index2axis(ii), '-', index2axis(jj)
                call dio_function_output(sys%outp%how(OPTION__OUTPUT__POL_DENSITY), EM_RESP_FD_DIR, trim(fname), &
                  sys%namespace, sys%space, sys%gr%mesh, -sys%gr%mesh%x(:, jj) * lr_rho(:, is), fn_unit, ierr, ions = sys%ions)

                fn_unit = units_out%length**(3-sys%space%dim) / units_out%energy**2
                write(fname, '(a,i1,6a)') 'beta_density-sp', is, '-', index2axis(ii), &
                  '-', index2axis(ii), '-', index2axis(jj)
                call dio_function_output(sys%outp%how(OPTION__OUTPUT__POL_DENSITY), EM_RESP_FD_DIR, trim(fname), &
                  sys%namespace, sys%space, sys%gr%mesh, -sys%gr%mesh%x(:, jj) * lr_rho2(:, is), fn_unit, ierr, ions = sys%ions)
              end do
            end if
          end do

        end if
      end if

      !ELF
      if (sys%outp%what(OPTION__OUTPUT__ELF)) then 
         
        if(isign == 1) then 
          call elf_calc(sys%st, sys%gr, sys%kpoints, elf, elfd)
        else
          call elf_calc(sys%st, sys%gr, sys%kpoints, lr_elf, lr_elfd)
          
          !numerical derivative
          lr_elf(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = &
               ( lr_elf(1:sys%gr%mesh%np, 1:sys%st%d%nspin) -  elf(1:sys%gr%mesh%np, 1:sys%st%d%nspin)) / (M_TWO * e_field)
          lr_elfd(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = &
               (lr_elfd(1:sys%gr%mesh%np, 1:sys%st%d%nspin) - elfd(1:sys%gr%mesh%np, 1:sys%st%d%nspin)) / (M_TWO * e_field)

          !write
          do is = 1, sys%st%d%nspin
            write(fname, '(a,i1,2a)') 'lr_elf-sp', is, '-', index2axis(ii)
            call dio_function_output(sys%outp%how(OPTION__OUTPUT__ELF), EM_RESP_FD_DIR, trim(fname),&
                sys%namespace, sys%space, sys%gr%mesh, lr_elf(:, is), unit_one, ierr, ions = sys%ions)
            write(fname, '(a,i1,2a)') 'lr_elf_D-sp', is, '-', index2axis(ii)
            call dio_function_output(sys%outp%how(OPTION__OUTPUT__ELF), EM_RESP_FD_DIR, trim(fname),&
                sys%namespace, sys%space, sys%gr%mesh, lr_elfd(:, is), unit_one, ierr, ions = sys%ions)
          end do
        end if

      end if

      POP_SUB(output_cycle_)
    end subroutine output_cycle_


    !-------------------------------------------------------------
    subroutine output_end_()
      FLOAT :: alpha(MAX_DIM, MAX_DIM)
      CMPLX :: beta(MAX_DIM, MAX_DIM, MAX_DIM)
      integer :: iunit, idir
      type(unit_t) :: fn_unit
      FLOAT :: freq_factor(3)

      PUSH_SUB(output_end_)

      call io_mkdir(EM_RESP_FD_DIR, sys%namespace)

      if ((sys%outp%what(OPTION__OUTPUT__DENSITY) .or. &
         sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) .and. calc_diagonal) then 
        lr_rho2(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = &
          -(sys%st%rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) - lr_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) &
          - tmp_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin) + gs_rho(1:sys%gr%mesh%np, 1:sys%st%d%nspin)) / e_field**2
  
        do is = 1, sys%st%d%nspin
          if (sys%outp%what(OPTION__OUTPUT__DENSITY)) then
            fn_unit = units_out%length**(2-sys%space%dim) / units_out%energy**2
            write(fname, '(a,i1,a)') 'fd2_density-sp', is, '-y-z'
            call dio_function_output(sys%outp%how(OPTION__OUTPUT__DENSITY), EM_RESP_FD_DIR, trim(fname),&
              sys%namespace, sys%space, sys%gr%mesh, lr_rho2(:, is), fn_unit, ierr, ions = sys%ions)
          end if
  
          if (sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then
            fn_unit = units_out%length**(3-sys%space%dim) / units_out%energy**2
            write(fname, '(a,i1,a)') 'beta_density-sp', is, '-x-y-z'
            call dio_function_output(sys%outp%how(OPTION__OUTPUT__POL_DENSITY), EM_RESP_FD_DIR, trim(fname),&
              sys%namespace, sys%space, sys%gr%mesh, -sys%gr%mesh%x(:, 1) * lr_rho2(:, is), fn_unit, ierr, ions = sys%ions)
          end if
        end do
      end if

      if(mpi_grp_is_root(mpi_world)) then ! output pol file
        iunit = io_open(EM_RESP_FD_DIR//'alpha', sys%namespace, action='write')
        write(iunit, '(3a)') '# Polarizability tensor [', trim(units_abbrev(units_out%polarizability)), ']'

        alpha(1:sys%space%dim, 1:sys%space%dim) = (dipole(1:sys%space%dim, 1:sys%space%dim, 1) - &
             dipole(1:sys%space%dim, 1:sys%space%dim, 2)) / (M_TWO * e_field)

        beta = M_ZERO

        do idir = 1, sys%space%dim
          beta(1:sys%space%dim, idir, idir) = &
            -(dipole(idir, 1:sys%space%dim, 1) + dipole(idir, 1:sys%space%dim, 2) - &
            M_TWO * center_dipole(1:sys%space%dim)) / e_field**2
          beta(idir, 1:sys%space%dim, idir) = beta(1:sys%space%dim, idir, idir) 
          beta(idir, idir, 1:sys%space%dim) = beta(1:sys%space%dim, idir, idir)
        end do

        if(calc_diagonal) then
          beta(1, 2, 3) = -(diag_dipole(1) - dipole(2, 1, 1) - dipole(3, 1, 1) + center_dipole(1)) / e_field**2
        else
          beta(1, 2, 3) = M_ZERO
        end if

        beta(2, 3, 1) = beta(1, 2, 3)
        beta(3, 1, 2) = beta(1, 2, 3)
        beta(3, 2, 1) = beta(1, 2, 3)
        beta(1, 3, 2) = beta(1, 2, 3)
        beta(2, 1, 3) = beta(1, 2, 3)

        call output_tensor(iunit, alpha, sys%space%dim, units_out%polarizability)
        call io_close(iunit)
        
        freq_factor(1:3) = M_ZERO ! for compatibility with em_resp version
        call out_hyperpolarizability(sys%gr%sb, beta, freq_factor(1:3), .true., EM_RESP_FD_DIR, sys%namespace)

        if (calc_Born) then
          call out_Born_charges(Born_charges, sys%ions, sys%namespace, sys%space%dim, &
            EM_RESP_FD_DIR, states_are_real(sys%st))
        end if
      end if

      if (sys%outp%what(OPTION__OUTPUT__DENSITY) .or. &
         sys%outp%what(OPTION__OUTPUT__POL_DENSITY)) then 
        SAFE_DEALLOCATE_A(lr_rho)
        SAFE_DEALLOCATE_A(lr_rho2)
      end if
      
      if (sys%outp%what(OPTION__OUTPUT__ELF)) then 
        SAFE_DEALLOCATE_A(lr_elf)
        SAFE_DEALLOCATE_A(elf)
        SAFE_DEALLOCATE_A(lr_elfd)
        SAFE_DEALLOCATE_A(elfd)
      end if

      POP_SUB(output_end_)
    end subroutine output_end_

  end subroutine static_pol_run_legacy

end module static_pol_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
