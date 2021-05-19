!! Copyright (C) 2007-2012 Xavier Andrade, David Strubbe
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

module phonons_lr_oct_m
  use born_charges_oct_m
  use epot_oct_m
  use forces_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use kdotp_oct_m
  use kdotp_calc_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use linear_response_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use smear_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use sternheimer_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use vibrations_oct_m
  use v_ks_oct_m

  implicit none

  private
  public :: &
       phonons_lr_run,    &
       phn_nm_wfs_tag,    &
       phn_wfs_tag,       &
       phn_rho_tag,       &
       axsf_mode_output
  
contains

  ! ---------------------------------------------------------
  subroutine phonons_lr_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(phonons_lr_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = vib_modes not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call phonons_lr_run_legacy(system, from_scratch)
    end select

    POP_SUB(phonons_lr_run)
  end subroutine phonons_lr_run

  ! ---------------------------------------------------------
  subroutine phonons_lr_run_legacy(sys, fromscratch)
    type(electrons_t), target, intent(inout) :: sys
    logical,                   intent(in)    :: fromscratch

    type(sternheimer_t) :: sh
    type(lr_t)          :: lr(1:1), kdotp_lr(MAX_DIM)
    type(vibrations_t)  :: vib
    type(pert_t)        :: ionic_pert

    type(ions_t),     pointer :: ions
    type(states_elec_t),   pointer :: st
    type(grid_t),     pointer :: gr

    integer :: natoms, ndim, iatom, idir, jatom, jdir, imat, jmat, iunit_restart, ierr, start_mode
    CMPLX, allocatable :: force_deriv(:,:)
    character(len=80) :: str_tmp
    character(len=300) :: line(1)
    type(Born_charges_t) :: born
    logical :: normal_mode_wfs, do_infrared, symmetrize
    type(restart_t) :: restart_load, restart_dump, kdotp_restart, gs_restart

    PUSH_SUB(phonons_lr_run_legacy)

    !some shortcuts

    ions => sys%ions
    st  => sys%st
    gr  => sys%gr

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if(sys%space%is_periodic()) then
      call messages_not_implemented('linear-response vib_modes for periodic systems')
    end if

    if(sys%hm%ep%nlcc) then
      call messages_not_implemented('linear-response vib_modes with non-linear core corrections')
    end if

    !%Variable CalcNormalModeWfs
    !%Type logical
    !%Default false
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, the response wavefunctions for each normal mode will be calculated
    !% and written in directory <tt>restart/vib_modes/phn_nm_wfs_XXXXX</tt>.
    !% This part is time-consuming and not parallel, but not needed for most purposes.
    !%End
    call parse_variable(sys%namespace, 'CalcNormalModeWfs', .false., normal_mode_wfs)

    !%Variable CalcInfrared
    !%Type logical
    !%Default true
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, infrared intensities (and Born charges) will be calculated
    !% and written in <tt>vib_modes/infrared</tt>.
    !%End
    call parse_variable(sys%namespace, 'CalcInfrared', .true., do_infrared)

    !%Variable SymmetrizeDynamicalMatrix
    !%Type logical
    !%Default true
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, all entries of the dynamical matrix will be calculated and then
    !% the matrix will be symmetrized to enforce <math>D_{ij} = D_{ji}</math>. If set to false,
    !% only the upper half of the matrix will be calculated.
    !%End
    call parse_variable(sys%namespace, 'SymmetrizeDynamicalMatrix', .true., symmetrize)

    ! replaced by properly saving and reading the dynamical matrix
    call messages_obsolete_variable(sys%namespace, 'UseRestartDontSolve')

    natoms = ions%natoms
    ndim = sys%space%dim

    call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)
    if(ierr == 0) then
      call states_elec_look_and_load(gs_restart, sys%namespace, sys%space, st, sys%gr%mesh, sys%kpoints)
      call restart_end(gs_restart)
    else
      message(1) = "Previous gs calculation is required."
      call messages_fatal(1)
    end if

    ! read kdotp wavefunctions if necessary (for IR intensities)
    if (sys%space%is_periodic() .and. do_infrared) then
      message(1) = "Reading kdotp wavefunctions for periodic directions."
      call messages_info(1)

      call restart_init(kdotp_restart, sys%namespace, RESTART_KDOTP, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)
      if(ierr /= 0) then
        message(1) = "Unable to read kdotp wavefunctions."
        message(2) = "Previous kdotp calculation required."
        call messages_fatal(2)
      end if

      do idir = 1, sys%space%periodic_dim
        call lr_init(kdotp_lr(idir))
        call lr_allocate(kdotp_lr(idir), sys%st, sys%gr%mesh)

        ! load wavefunctions
        str_tmp = trim(kdotp_wfs_tag(idir))
        call restart_open_dir(kdotp_restart, wfs_tag_sigma(str_tmp, 1), ierr)
        if (ierr == 0) then
          call states_elec_load(kdotp_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr, lr=kdotp_lr(idir))
        end if
        call restart_close_dir(kdotp_restart)

        if(ierr /= 0) then
          message(1) = "Unable to read kdotp wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'."
          message(2) = "Previous kdotp calculation required."
          call messages_fatal(2)
        end if
      end do
      call restart_end(kdotp_restart)
    end if

    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call messages_info(1)

    call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
    call sternheimer_init(sh, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, &
        wfs_are_cplx = states_are_complex(st))
    call vibrations_init(vib, ions, "lr", sys%namespace)

    call epot_precalc_local_potential(sys%hm%ep, sys%namespace, sys%gr, sys%ions)

    if(do_infrared) then
      call born_charges_init(born, sys%namespace, ions, st, ndim)
    end if
    SAFE_ALLOCATE(force_deriv(1:ndim, 1:natoms))

    !CALCULATE

    !the ionic contribution
    call build_ionic_dyn_matrix()

    !the  <phi0 | v2 | phi0> term
    if(states_are_real(st)) then
      call dionic_pert_matrix_elements_2(sys%gr, sys%namespace, sys%ions, sys%hm, 1, st, vib, CNST(-1.0), vib%dyn_matrix)
    else
      call zionic_pert_matrix_elements_2(sys%gr, sys%namespace, sys%ions, sys%hm, 1, st, vib, CNST(-1.0), vib%dyn_matrix)
    end if

    call pert_init(ionic_pert, sys%namespace, PERTURBATION_IONIC, gr, ions)

    call lr_init(lr(1))
    call lr_allocate(lr(1), st, gr%mesh)

    call restart_init(restart_dump, sys%namespace, RESTART_VIB_MODES, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=gr%mesh)
    call restart_init(restart_load, sys%namespace, RESTART_VIB_MODES, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh)

    if (fromScratch) then
      start_mode = 1
    else
      call phonons_load(restart_load, vib, start_mode)
    end if

    ! Delete, if fromScratch, or trying to open it failed and there is something wrong with it.
    if(start_mode == 1) call restart_rm(restart_dump, 'restart')

    do imat = 1, start_mode - 1
      call vibrations_out_dyn_matrix_row(vib, imat)
    end do

    do imat = start_mode, vib%num_modes
      iatom = vibrations_get_atom(vib, imat)
      idir  = vibrations_get_dir (vib, imat)
      
      write(message(1),'(a,i5,a,a1,a)') &
        "Calculating response to displacement of atom ", iatom, " in ", index2axis(idir), "-direction."
      call messages_info(1)

      ! the converged wfns for the previous mode are probably not a good starting point
      call lr_zero(lr(1), st)

      if (.not. fromscratch) then
        message(1) = "Loading restart wavefunctions for linear response."
        call messages_info(1)
        call restart_open_dir(restart_load, wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1), ierr)
        if (ierr == 0) then
          call states_elec_load(restart_load, sys%namespace, sys%space, st, sys%gr%mesh, sys%kpoints, ierr, lr = lr(1))
        end if
        if (ierr /= 0) then
          message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))//"'."
          call messages_warning(1)
        end if
        call restart_close_dir(restart_load)
      end if
      
      call pert_setup_atom(ionic_pert, iatom)
      call pert_setup_dir(ionic_pert, idir)
      
      if(states_are_real(st)) then
        call dsternheimer_solve(sh, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, sys%ions, &
          lr, 1, M_ZERO, ionic_pert, restart_dump, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))
      else
        call zsternheimer_solve(sh, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, sys%ions, &
          lr, 1, M_z0, ionic_pert, restart_dump, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))
      end if
      
      if(states_are_real(st)) then
        call dforces_derivative(gr, sys%namespace, sys%space, ions, sys%hm%ep, st, sys%kpoints, lr(1), lr(1), force_deriv, &
          sys%hm%lda_u_level)
      else
        call zforces_derivative(gr, sys%namespace, sys%space, ions, sys%hm%ep, st, sys%kpoints, lr(1), lr(1), force_deriv, &
          sys%hm%lda_u_level)
      end if

      do jmat = 1, vib%num_modes
        if(.not. symmetrize .and. jmat < imat) then
          vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(imat, jmat)
          cycle
        end if

        jatom = vibrations_get_atom(vib, jmat)
        jdir  = vibrations_get_dir (vib, jmat)

        vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(jmat, imat) + TOFLOAT(force_deriv(jdir, jatom))
        vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(jmat, imat) * vibrations_norm_factor(vib, ions, iatom, jatom)
      end do
      call vibrations_out_dyn_matrix_row(vib, imat)
      
      if(do_infrared) then
        if(states_are_real(st)) then
          call dphonons_lr_infrared(gr%mesh, ions, st, lr(1), kdotp_lr, imat, iatom, idir, vib%infrared)
        else
          call zphonons_lr_infrared(gr%mesh, ions, st, lr(1), kdotp_lr, imat, iatom, idir, vib%infrared)
        end if
      end if

      iunit_restart = restart_open(restart_dump, 'restart', position='append')
      ! open and close makes sure output is not buffered
      do jmat = 1, vib%num_modes
        write(line(1), *) jmat, imat, vib%dyn_matrix(jmat, imat)
        call restart_write(restart_dump, iunit_restart, line, 1, ierr)
        if (ierr /= 0) then
          message(1) = "Could not write restart information."
          call messages_warning(1)
        end if
      end do
      write(line(1), *) imat, (vib%infrared(imat, idir), idir = 1, ndim)
      call restart_write(restart_dump, iunit_restart, line, 1, ierr)
      if (ierr /= 0) then
        message(1) = "Could not write restart information."
        call messages_warning(1)
      end if
      call restart_close(restart_dump, iunit_restart)

      message(1) = ""
      call messages_info(1)
    end do 

    call pert_end(ionic_pert)

    if(symmetrize) call vibrations_symmetrize_dyn_matrix(vib)
    call vibrations_diag_dyn_matrix(vib)
    call vibrations_output(vib)
    call axsf_mode_output(vib, ions, gr%mesh, sys%namespace)

    if(do_infrared) then
      if(sys%space%is_periodic() .and. .not. smear_is_semiconducting(st%smear)) then
        message(1) = "Cannot calculate infrared intensities for periodic system with smearing (i.e. without a gap)."
        call messages_info(1)
      else
        call born_from_infrared(vib, born)
        call out_Born_charges(born, ions, sys%namespace, ndim, VIB_MODES_DIR, write_real = .true.)
        call calc_infrared()
      end if

      call Born_charges_end(born)
    end if

    if(normal_mode_wfs) then
      message(1) = "Calculating response wavefunctions for normal modes."
      call messages_info(1)
      if(states_are_real(st)) then
        call dphonons_lr_wavefunctions(lr(1), sys%namespace, sys%space, st, sys%gr%mesh, sys%kpoints, vib, restart_load, &
          restart_dump)
      else
        call zphonons_lr_wavefunctions(lr(1), sys%namespace, sys%space, st, sys%gr%mesh, sys%kpoints, vib, restart_load, &
          restart_dump)
      end if
    end if

    !DESTRUCT

    SAFE_DEALLOCATE_A(force_deriv)
    call lr_dealloc(lr(1))
    call vibrations_end(vib)
    call sternheimer_end(sh)
    call states_elec_deallocate_wfns(st)
    if (sys%space%is_periodic() .and. do_infrared) then
      do idir = 1, sys%space%periodic_dim
        call lr_dealloc(kdotp_lr(idir))
      end do
    end if
    call restart_end(restart_load)
    call restart_end(restart_dump)

    POP_SUB(phonons_lr_run_legacy)

  contains

    ! ---------------------------------------------------------
    !> this formulation is only valid for finite systems, or an Ewald sum is required
    !! as in Baroni et al. RMP 2001, Appendix B.
    subroutine build_ionic_dyn_matrix()

      FLOAT :: term, xi(1:MAX_DIM), xj(1:MAX_DIM), r2

      PUSH_SUB(phonons_lr_run_legacy.build_ionic_dyn_matrix)

      vib%dyn_matrix(:,:) = M_ZERO

      do iatom = 1, natoms
        xi(1:ndim) = ions%pos(:, iatom)

        do idir = 1, ndim

          do jatom = 1, natoms
            if(iatom == jatom) cycle

            do jdir = 1, ndim         

              xj(1:ndim) = ions%pos(:, jatom)
              r2 = sum((xi(1:ndim) - xj(1:ndim))**2)

              term = species_zval(ions%atom(iatom)%species) * species_zval(ions%atom(jatom)%species) &
                /(r2**CNST(1.5))*(ddelta(idir, jdir) - (M_THREE*(xi(idir)-xj(idir))*(xi(jdir)-xj(jdir)))/r2)

              ! note: this accomplishes the sum over k for diagonal terms, using the j loop
              vib%dyn_matrix(vibrations_get_index(vib, iatom, jdir), vibrations_get_index(vib, iatom, idir)) = &
                vib%dyn_matrix(vibrations_get_index(vib, iatom, jdir), vibrations_get_index(vib, iatom, idir)) + term

              vib%dyn_matrix(vibrations_get_index(vib, jatom, jdir), vibrations_get_index(vib, iatom, idir)) = &
                vib%dyn_matrix(vibrations_get_index(vib, jatom, jdir), vibrations_get_index(vib, iatom, idir)) - term
            end do
          end do
        end do
      end do
      POP_SUB(phonons_lr_run_legacy.build_ionic_dyn_matrix)

    end subroutine build_ionic_dyn_matrix

    ! ---------------------------------------------------------
    !> calculate infrared intensities
    subroutine calc_infrared()

      integer :: iunit_ir
      FLOAT :: lir(1:MAX_DIM+1)

      PUSH_SUB(phonons_lr_run_legacy.calc_infrared)

      iunit_ir = io_open(VIB_MODES_DIR//'infrared', sys%namespace, action='write')

      write(iunit_ir, '(a)', advance = 'no') '#   freq ['//trim(units_abbrev(unit_invcm))//']'
      do idir = 1, ndim
        write(iunit_ir, '(a14)', advance = 'no') '<' // index2axis(idir) // '> [' // trim(units_abbrev(units_out%length)) // ']'
      end do
      write(iunit_ir, '(a14)') 'average [' // trim(units_abbrev(units_out%length)) // ']'

      do iatom = 1, natoms
        do idir = 1, ndim

          imat = vibrations_get_index(vib, iatom, idir)

          write(iunit_ir, '(f17.8)', advance = 'no') units_from_atomic(unit_invcm, vib%freq(imat))
          do jdir = 1, ndim
            lir(jdir) = dot_product(vib%infrared(:, jdir), vib%normal_mode(:, imat))
            write(iunit_ir, '(f14.5)', advance = 'no') units_from_atomic(units_out%length, lir(jdir))
          end do

          lir(ndim+1) = sqrt(sum(lir(1:ndim)**2)/ndim)
          write(iunit_ir, '(f17.8)') units_from_atomic(units_out%length, lir(ndim + 1))
        end do
      end do

      call io_close(iunit_ir)
      POP_SUB(phonons_lr_run_legacy.calc_infrared)
    end subroutine calc_infrared

  end subroutine phonons_lr_run_legacy


  ! ---------------------------------------------------------
  subroutine born_from_infrared(vib, born)
    type(vibrations_t),   intent(in)    :: vib
    type(Born_charges_t), intent(inout) :: born

    integer :: imat, idir, iatom

    PUSH_SUB(born_from_infrared)

    do imat = 1, vib%num_modes
      idir = vibrations_get_dir(vib, imat)
      iatom = vibrations_get_atom(vib, imat)
      born%charge(1:vib%ndim, idir, iatom) = -vib%infrared(imat, 1:vib%ndim)
    end do

    POP_SUB(born_from_infrared)
  end subroutine born_from_infrared


  ! ---------------------------------------------------------
  character(len=100) function phn_rho_tag(iatom, dir) result(str)
    integer, intent(in) :: iatom, dir
    
    PUSH_SUB(phn_rho_tag)
    
    write(str, '(a,i4.4,a,i1)') 'phn_rho_', iatom, '_',  dir

    POP_SUB(phn_rho_tag)

  end function phn_rho_tag
  

  ! ---------------------------------------------------------
  character(len=100) function phn_wfs_tag(iatom, dir) result(str)
    integer, intent(in) :: iatom, dir

    PUSH_SUB(phn_wfs_tag)

    write(str, '(a,i4.4,a,a)') "phn_wfs_", iatom, "_", index2axis(dir)

    POP_SUB(phn_wfs_tag)
    
  end function phn_wfs_tag


  ! ---------------------------------------------------------
  character(len=100) function phn_nm_wfs_tag(inm) result(str)
    integer, intent(in) :: inm

    PUSH_SUB(phn_nm_wfs_tag)

    write(str, '(a,i5.5)') "phn_nm_wfs_", inm

    POP_SUB(phn_nm_wfs_tag)
    
  end function phn_nm_wfs_tag


  ! ---------------------------------------------------------
  !> output eigenvectors as animated XSF file, one per frame, displacements as forces
  subroutine axsf_mode_output(this, ions, mesh, namespace)
    type(vibrations_t), intent(in) :: this
    type(ions_t),       intent(in) :: ions
    type(mesh_t),       intent(in) :: mesh
    type(namespace_t),  intent(in) :: namespace
    
    integer :: iunit, iatom, idir, imat, jmat
    FLOAT, allocatable :: forces(:,:)
    character(len=2) :: suffix

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(axsf_mode_output)

    ! for some reason, direct usage of this%suffix gives an odd result
    suffix = vibrations_get_suffix(this)
    iunit = io_open(VIB_MODES_DIR//'normal_modes_'//suffix//'.axsf', namespace, action='write')

    write(iunit, '(a,i6)') 'ANIMSTEPS ', this%num_modes
    SAFE_ALLOCATE(forces(1:ions%space%dim, 1:ions%natoms))
    do imat = 1, this%num_modes
      do jmat = 1, this%num_modes
        iatom = vibrations_get_atom(this, jmat)
        idir  = vibrations_get_dir (this, jmat)
        forces(idir, iatom) = this%normal_mode(jmat, imat)
      end do
      call write_xsf_geometry(iunit, ions, mesh, forces = forces, index = imat)
    end do
    SAFE_DEALLOCATE_A(forces)
    call io_close(iunit)

    POP_SUB(axsf_mode_output)
  end subroutine axsf_mode_output

  ! ---------------------------------------------------------
  subroutine phonons_load(restart, vib, start_mode)
    type(restart_t),    intent(in)    :: restart
    type(vibrations_t), intent(inout) :: vib
    integer,            intent(out)   :: start_mode

    integer :: iunit, ierr, imode, jmode, imode_read, jmode_read
    character(len=120) :: line(1)

    PUSH_SUB(phonons_load)

    iunit = restart_open(restart, 'restart')
    if (iunit > 0) then
      imode_loop: do imode = 1, vib%num_modes
        do jmode = 1, vib%num_modes
          call restart_read(restart, iunit, line, 1, ierr)
          if(ierr /= 0) exit imode_loop
          read(line(1), fmt=*, iostat=ierr) jmode_read, imode_read, vib%dyn_matrix(jmode, imode)
          if(imode_read /= imode) then
            write(message(1),'(a,i9,a,i9)') "Corruption of restart data: row ", imode, " is labeled as ", imode_read
            call messages_fatal(1)
          end if
          if(jmode_read /= jmode) then
            write(message(1),'(a,i9,a,i9)') "Corruption of restart data: column ", jmode, " is labeled as ", jmode_read
            call messages_fatal(1)
          end if
        end do

        call restart_read(restart, iunit, line, 1, ierr)
        if(ierr /= 0) exit

        start_mode = imode + 1

        read(line(1), fmt=*, iostat=ierr) imode_read, vib%infrared(imode, 1:vib%ndim)
        if(imode_read /= imode) then
          write(message(1),'(a,i9,a,i9)') "Corruption of restart data: infrared row ", imode, " is labeled as ", imode_read
          call messages_fatal(1)
        end if
      end do imode_loop
        
      write(message(1),'(a,i9,a,i9)') 'Info: Read saved dynamical-matrix rows for ', &
        start_mode - 1, ' modes out of ', vib%num_modes
      call messages_info(1)
      
      call restart_close(restart, iunit)
    else
      start_mode = 1

      message(1) = "Could not open restart file 'restart'. Starting from scratch."
      call messages_warning(1)
    end if

    POP_SUB(phonons_load)
  end subroutine phonons_load

#include "complex.F90"
#include "phonons_lr_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "phonons_lr_inc.F90"

end module phonons_lr_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
