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
!! $Id$

#include "global.h"

module phonons_lr_m
  use born_charges_m
  use datasets_m
  use epot_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kdotp_m
  use kdotp_calc_m
  use lalg_basic_m
  use linear_response_m
  use loct_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use output_m
  use pert_m
  use profiling_m
  use projector_m
  use restart_m
  use simul_box_m
  use smear_m
  use species_m
  use states_m
  use states_dim_m
  use states_restart_m
  use sternheimer_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m
  use vibrations_m

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
  subroutine phonons_lr_run(sys, hm, fromscratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(in)    :: fromscratch

    type(sternheimer_t) :: sh
    type(lr_t)          :: lr(1:1), kdotp_lr(MAX_DIM)
    type(vibrations_t)  :: vib
    type(pert_t)        :: ionic_pert

    type(geometry_t), pointer :: geo
    type(states_t),   pointer :: st
    type(grid_t),     pointer :: gr

    integer :: natoms, ndim, iatom, idir, jatom, jdir, imat, jmat, iunit_restart, ierr, start_mode
    CMPLX, allocatable :: force_deriv(:,:)
    character(len=80) :: str_tmp
    character(len=300) :: line(1)
    type(Born_charges_t) :: born
    logical :: normal_mode_wfs, do_infrared, symmetrize
    type(restart_t) :: restart_load, restart_dump, kdotp_restart, gs_restart

    PUSH_SUB(phonons_lr_run)

    !some shortcuts

    geo => sys%geo
    st  => sys%st
    gr  => sys%gr

    if(simul_box_is_periodic(gr%sb)) then
      call messages_not_implemented('linear-response vib_modes for periodic systems')
    endif

    if(geo%nlcc) then
      call messages_not_implemented('linear-response vib_modes with non-linear core corrections')
    endif

    !%Variable CalcNormalModeWfs
    !%Type logical
    !%Default false
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, the response wavefunctions for each normal mode will be calculated
    !% and written in directory <tt>restart/vib_modes/phn_nm_wfs_XXXXX</tt>.
    !% This part is time-consuming and not parallel, but not needed for most purposes.
    !%End
    call parse_logical(datasets_check('CalcNormalModeWfs'), .false., normal_mode_wfs)

    !%Variable CalcInfrared
    !%Type logical
    !%Default true
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, infrared intensities (and Born charges) will be calculated
    !% and written in <tt>vib_modes/infrared</tt>.
    !%End
    call parse_logical(datasets_check('CalcInfrared'), .true., do_infrared)

    !%Variable SymmetrizeDynamicalMatrix
    !%Type logical
    !%Default true
    !%Section Linear Response::Vibrational Modes
    !%Description
    !% If set to true, all entries of the dynamical matrix will be calculated and then
    !% the matrix will be symmetrized to enforce mat(i,j) = mat(j,i). If set to false,
    !% only the upper half of the matrix will be calculated.
    !%End
    call parse_logical(datasets_check('SymmetrizeDynamicalMatrix'), .true., symmetrize)

    ! replaced by properly saving and reading the dynamical matrix
    call messages_obsolete_variable('UseRestartDontSolve')

    natoms = geo%natoms
    ndim = gr%mesh%sb%dim

    call restart_init(gs_restart, RESTART_GS, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, mesh=gr%mesh, sb=gr%sb, exact=.true.)
    call states_look_and_load(gs_restart, st, gr)
    call restart_end(gs_restart)

    ! read kdotp wavefunctions if necessary (for IR intensities)
    if (simul_box_is_periodic(gr%sb) .and. do_infrared) then
      message(1) = "Reading kdotp wavefunctions for periodic directions."
      call messages_info(1)

      call restart_init(kdotp_restart, RESTART_KDOTP, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, mesh=gr%mesh, sb=gr%sb)
      do idir = 1, gr%sb%periodic_dim
        call lr_init(kdotp_lr(idir))
        call lr_allocate(kdotp_lr(idir), sys%st, sys%gr%mesh)

        ! load wavefunctions
        str_tmp = trim(kdotp_wfs_tag(idir))
        call restart_cd(kdotp_restart, dirname=wfs_tag_sigma(str_tmp, 1))
        call states_load(kdotp_restart, sys%st, sys%gr, ierr, lr=kdotp_lr(idir))
        call restart_cd(kdotp_restart)

        if(ierr /= 0) then
          message(1) = "Unable to read kdotp wavefunctions from '"//trim(wfs_tag_sigma(str_tmp, 1))//"'."
          message(2) = "Previous kdotp calculation required."
          call messages_fatal(2)
        end if
      end do
      call restart_end(kdotp_restart)
    endif

    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call messages_info(1)

    call system_h_setup(sys, hm)
    call sternheimer_init(sh, sys, hm, "VM", wfs_are_cplx = states_are_complex(st))
    call vibrations_init(vib, geo, gr%sb, "lr")

    call epot_precalc_local_potential(hm%ep, sys%gr, sys%geo)

    if(do_infrared) then
      call Born_charges_init(born, geo, st, ndim)
    endif
    SAFE_ALLOCATE(force_deriv(1:ndim, 1:natoms))

    !CALCULATE

    !the ionic contribution
    call build_ionic_dyn_matrix()

    !the  <phi0 | v2 | phi0> term
    if(states_are_real(st)) then
      call dionic_pert_matrix_elements_2(sys%gr, sys%geo, hm, 1, st, st%dpsi(:, :, :, 1), vib, CNST(-1.0), vib%dyn_matrix)
    else
      call zionic_pert_matrix_elements_2(sys%gr, sys%geo, hm, 1, st, st%zpsi(:, :, :, 1), vib, CNST(-1.0), vib%dyn_matrix)
    endif

    call pert_init(ionic_pert, PERTURBATION_IONIC, gr, geo)

    call lr_init(lr(1))
    call lr_allocate(lr(1), st, gr%mesh)

    call restart_init(restart_dump, RESTART_VIB_MODES, RESTART_TYPE_DUMP, st%dom_st_kpt_mpi_grp, mesh=gr%mesh, sb=gr%sb)
    call restart_init(restart_load, RESTART_VIB_MODES, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, mesh=gr%mesh, sb=gr%sb)

    if (fromScratch) then
      start_mode = 1
      call restart_rm(restart_dump, 'restart')
    else
      call phonons_load(restart_load, fromScratch, vib, start_mode)
    end if

    do imat = 1, start_mode - 1
      call vibrations_out_dyn_matrix_row(vib, imat)
    enddo

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
        call restart_cd(restart_load, dirname=wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))
        call states_load(restart_load, st, gr, ierr, lr = lr(1))
        if (ierr /= 0) then
          message(1) = "Unable to read response wavefunctions from '"//trim(wfs_tag_sigma(phn_wfs_tag(iatom, idir), 1))//"'."
          call messages_fatal(1)
        end if
        call restart_cd(restart_load)
      end if
      
      call pert_setup_atom(ionic_pert, iatom)
      call pert_setup_dir(ionic_pert, idir)
      
      if(states_are_real(st)) then
        call dsternheimer_solve(sh, sys, hm, lr, 1, M_ZERO, ionic_pert, &
          restart_dump, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))
      else
        call zsternheimer_solve(sh, sys, hm, lr, 1, M_z0, ionic_pert, &
          restart_dump, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))
      endif
      
      if(states_are_real(st)) then
        call dforces_derivative(gr, geo, hm%ep, st, lr(1), lr(1), force_deriv)
      else
        call zforces_derivative(gr, geo, hm%ep, st, lr(1), lr(1), force_deriv)
      endif

      do jmat = 1, vib%num_modes
        if(.not. symmetrize .and. jmat < imat) then
          vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(imat, jmat)
          cycle
        endif

        jatom = vibrations_get_atom(vib, jmat)
        jdir  = vibrations_get_dir (vib, jmat)

        vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(jmat, imat) + TOFLOAT(force_deriv(jdir, jatom))
        vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(jmat, imat) * vibrations_norm_factor(vib, geo, iatom, jatom)
      end do
      call vibrations_out_dyn_matrix_row(vib, imat)
      
      if(do_infrared) then
        if(states_are_real(st)) then
          call dphonons_lr_infrared(gr, geo, st, lr(1), kdotp_lr, imat, iatom, idir, vib%infrared)
        else
          call zphonons_lr_infrared(gr, geo, st, lr(1), kdotp_lr, imat, iatom, idir, vib%infrared)
        endif
      endif

      iunit_restart = restart_open(restart_dump, 'restart', position='append')
      ! open and close makes sure output is not buffered
      write(line(1), *) imat, vib%dyn_matrix(:, imat), (vib%infrared(imat, idir), idir = 1, ndim)
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
    call axsf_mode_output(vib, geo, gr%mesh)

    if(do_infrared) then
      if(simul_box_is_periodic(gr%sb) .and. .not. smear_is_semiconducting(st%smear)) then
        message(1) = "Cannot calculate infrared intensities for periodic system with smearing (i.e. without a gap)."
        call messages_info(1)
      else
        call born_from_infrared(vib, born)
        call out_Born_charges(born, geo, ndim, VIB_MODES_DIR, write_real = .true.)
        call calc_infrared()
      endif

      call Born_charges_end(born)
    endif

    if(normal_mode_wfs) then
      message(1) = "Calculating response wavefunctions for normal modes."
      call messages_info(1)
      if(states_are_real(st)) then
        call dphonons_lr_wavefunctions(lr(1), st, gr, vib, restart_load, restart_dump)
      else
        call zphonons_lr_wavefunctions(lr(1), st, gr, vib, restart_load, restart_dump)
      endif
    endif

    !DESTRUCT

    SAFE_DEALLOCATE_A(force_deriv)
    call lr_dealloc(lr(1))
    call vibrations_end(vib)
    call sternheimer_end(sh)
    call states_deallocate_wfns(st)
    if (simul_box_is_periodic(gr%sb) .and. do_infrared) then
      do idir = 1, gr%sb%periodic_dim
        call lr_dealloc(kdotp_lr(idir))
      enddo
    endif
    call restart_end(restart_load)
    call restart_end(restart_dump)

    POP_SUB(phonons_lr_run)

  contains

    ! ---------------------------------------------------------
    !> this formulation is only valid for finite systems, or an Ewald sum is required
    !! as in Baroni et al. RMP 2001, Appendix B.
    subroutine build_ionic_dyn_matrix()

      FLOAT :: ac, xi(1:MAX_DIM), xj(1:MAX_DIM), xk(1:MAX_DIM), r2
      integer :: katom

      PUSH_SUB(phonons_lr_run.build_ionic_dyn_matrix)

      do iatom = 1, natoms
        do idir = 1, ndim

          do jatom = 1, natoms
            do jdir = 1, ndim         

              xi(1:ndim) = geo%atom(iatom)%x(1:ndim)

              !ion - ion
              if( iatom == jatom) then 

                ac = M_ZERO
                do katom = 1, natoms
                  if ( katom == iatom ) cycle

                  xk(1:ndim) = geo%atom(katom)%x(1:ndim)
                  r2 = sum((xi(1:ndim) - xk(1:ndim))**2)

                  ac = ac + species_zval(geo%atom(iatom)%spec) * &
                            species_zval(geo%atom(katom)%spec) &
                       /(r2**CNST(1.5)) *(&
                       -ddelta(idir, jdir) + &
                       (M_THREE*(xi(idir)-xk(idir))*(xi(jdir)-xk(jdir)))/r2 &
                       )

                end do

              else ! iatom /= jatom

                xj(1:ndim) = geo%atom(jatom)%x(1:ndim)

                r2 = sum((xi(1:ndim) - xj(1:ndim))**2)
                ac = species_zval(geo%atom(iatom)%spec) * species_zval(geo%atom(jatom)%spec) &
                     /(r2**CNST(1.5))*(&
                     ddelta(idir, jdir) - (M_THREE*(xi(idir)-xj(idir))*(xi(jdir)-xj(jdir)))/r2)

              end if


              vib%dyn_matrix(vibrations_get_index(vib, jatom, jdir), vibrations_get_index(vib, iatom, idir)) = -ac

            end do
          end do
        end do
      end do
      POP_SUB(phonons_lr_run.build_ionic_dyn_matrix)

    end subroutine build_ionic_dyn_matrix

    ! ---------------------------------------------------------
    !> calculate infrared intensities
    subroutine calc_infrared()

      integer :: iunit_ir
      FLOAT :: lir(1:MAX_DIM+1)

      PUSH_SUB(phonons_lr_run.calc_infrared)

      iunit_ir = io_open(VIB_MODES_DIR//'infrared', action='write')

      write(iunit_ir, '(a)', advance = 'no') '#   freq ['//trim(units_abbrev(unit_invcm))//']'
      do idir = 1, ndim
        write(iunit_ir, '(a14)', advance = 'no') '<' // index2axis(idir) // '> [' // trim(units_abbrev(units_out%length)) // ']'
      enddo
      write(iunit_ir, '(a14)') 'average [' // trim(units_abbrev(units_out%length)) // ']'

      do iatom = 1, natoms
        do idir = 1, ndim

          imat = vibrations_get_index(vib, iatom, idir)

          write(iunit_ir, '(f14.5)', advance = 'no') units_from_atomic(unit_invcm, vib%freq(imat))
          do jdir = 1, ndim
            lir(jdir) = dot_product(vib%infrared(:, jdir), vib%normal_mode(:, imat))
            write(iunit_ir, '(f14.5)', advance = 'no') units_from_atomic(units_out%length, lir(jdir))
          end do

          lir(ndim+1) = sqrt(sum(lir(1:ndim)**2)/M_THREE)
          write(iunit_ir, '(f14.5)') units_from_atomic(units_out%length, lir(ndim + 1))
        end do
      end do

      call io_close(iunit_ir)
      POP_SUB(phonons_lr_run.calc_infrared)
    end subroutine calc_infrared

  end subroutine phonons_lr_run


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
    enddo

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
  subroutine axsf_mode_output(this, geo, mesh)
    type(vibrations_t), intent(in) :: this
    type(geometry_t),   intent(in) :: geo
    type(mesh_t),       intent(in) :: mesh
    
    integer :: iunit, iatom, idir, imat, jmat
    FLOAT, allocatable :: forces(:,:)
    character(len=2) :: suffix

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(axsf_mode_output)

    ! for some reason, direct usage of this%suffix gives an odd result
    suffix = vibrations_get_suffix(this)
    iunit = io_open(VIB_MODES_DIR//'normal_modes_'//suffix//'.axsf', action='write')

    write(iunit, '(a,i6)') 'ANIMSTEPS ', this%num_modes
    SAFE_ALLOCATE(forces(1:geo%natoms, 1:mesh%sb%dim))
    do imat = 1, this%num_modes
      do jmat = 1, this%num_modes
        iatom = vibrations_get_atom(this, jmat)
        idir  = vibrations_get_dir (this, jmat)
        forces(iatom, idir) = this%normal_mode(jmat, imat)
      enddo
      call write_xsf_geometry(iunit, geo, mesh, forces = forces, index = imat)
    end do
    SAFE_DEALLOCATE_A(forces)
    call io_close(iunit)

    POP_SUB(axsf_mode_output)
  end subroutine axsf_mode_output

  ! ---------------------------------------------------------
  subroutine phonons_load(restart, fromScratch, vib, start_mode)
    type(restart_t),    intent(in)    :: restart
    logical,            intent(in)    :: fromScratch
    type(vibrations_t), intent(inout) :: vib
    integer,            intent(out)   :: start_mode

    integer :: iunit, ierr, imode, number
    FLOAT, allocatable :: dyn_row(:)
    FLOAT :: infrared(MAX_DIM)
    character(len=120) :: line(1)

    PUSH_SUB(phonons_load)

    iunit = restart_open(restart, 'restart')
    if (iunit > 0) then
      SAFE_ALLOCATE(dyn_row(1:vib%num_modes))

      do imode = 1, vib%num_modes
        call restart_read(restart, iunit, line, 1, ierr)
        if(ierr /= 0) exit
        read(line(1), fmt=*, iostat=ierr) number, dyn_row(:), infrared(1:vib%ndim)
        if(number /= imode) then
          write(message(1),'(a,i9,a,i9)') "Corruption of restart data: line ", imode, " is labeled as ", number
          call messages_fatal(1)
        endif
        start_mode = imode + 1
        vib%dyn_matrix(:, imode) = dyn_row(:)
        vib%infrared(imode, 1:vib%ndim) = infrared(1:vib%ndim)
      end do

      write(message(1),'(a,i9,a,i9)') 'Info: Read saved dynamical-matrix rows for ', &
           start_mode - 1, ' modes out of ', vib%num_modes
      call messages_info(1)

      SAFE_DEALLOCATE_A(dyn_row)
    else
      start_mode = 1
      
      message(1) = "Could not find restart file 'restart'. Starting from scratch."
      call messages_warning(1)
    endif
    call restart_close(restart, iunit)

#ifdef HAVE_MPI
    call MPI_Bcast(start_mode, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
#endif

    POP_SUB(phonons_load)
  end subroutine phonons_load

#include "complex.F90"
#include "phonons_lr_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "phonons_lr_inc.F90"

end module phonons_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
