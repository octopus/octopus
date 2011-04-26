!! Copyright (C) 2007 Xavier Andrade
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
!! $Id: em_resp.F90 2686 2007-02-03 22:10:51Z xavier $

#include "global.h"

module phonons_lr_m
  use datasets_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use linear_response_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use pert_m
  use profiling_m
  use projector_m
  use restart_m
  use simul_box_m
  use species_m
  use states_m
  use states_dim_m
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
       phn_rho_tag
  
contains

  ! ---------------------------------------------------------
  subroutine phonons_lr_run(sys, hm, fromscratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(in)    :: fromscratch

    type(sternheimer_t) :: sh
    type(lr_t)          :: lr(1:1)
    type(vibrations_t)  :: vib
    type(pert_t)        :: ionic_pert, electric_pert

    type(geometry_t), pointer :: geo
    type(states_t),   pointer :: st
    type(grid_t),     pointer :: gr

    integer :: natoms, ndim, iatom, idir, jatom, jdir, imat, jmat, iunit, ierr
    FLOAT, allocatable   :: infrared(:,:)

    PUSH_SUB(phonons_lr_run)

    !some shortcuts

    geo => sys%geo
    st  => sys%st
    gr  => sys%gr

    if(simul_box_is_periodic(gr%sb)) then
      call messages_not_implemented('linear-response vib_modes for periodic systems')
    endif

    natoms = geo%natoms
    ndim = gr%mesh%sb%dim

    !CONSTRUCT

    call restart_look_and_read(st, gr, geo)

    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call messages_info(1)

    call system_h_setup(sys, hm)
    call sternheimer_init(sh, sys, hm, "VM", .false.)
    call vibrations_init(vib, geo, gr%sb)

    call epot_precalc_local_potential(hm%ep, sys%gr, sys%geo, time = M_ZERO)

    SAFE_ALLOCATE(infrared(1:natoms*ndim, 1:ndim))

    !CALCULATE

    !the ionic contribution
    call build_ionic_dyn_matrix()

    !the  <phi0 | v2 | phi0> term
    call dionic_pert_matrix_elements_2(ionic_pert, sys%gr, sys%geo, hm, 1, st, st%dpsi(:, :, :, 1), vib, CNST(-1.0), vib%dyn_matrix)

    call pert_init(ionic_pert, PERTURBATION_IONIC, gr, geo)
    call pert_init(electric_pert, PERTURBATION_ELECTRIC, gr, geo)

    call lr_init(lr(1))
    call lr_allocate(lr(1), st, gr%mesh)

    do imat = 1, vib%num_modes
      iatom = vibrations_get_atom(vib, imat)
      idir  = vibrations_get_dir (vib, imat)
      
      write(message(1),'(a,i5,a,a1,a)') &
        "Calculating response to displacement of atom ", iatom, " in ", index2axis(idir), "-direction."
      call messages_info(1)

      if (.not. fromscratch) then
        message(1) = "Loading restart wavefunctions for linear response."
        call messages_info(1)
        call restart_read(trim(restart_dir)//VIB_MODES_DIR//trim(phn_wfs_tag(iatom, idir)), &
          st, gr, geo, ierr, lr = lr(1))
      end if
      
      call pert_setup_atom(ionic_pert, iatom)
      call pert_setup_dir(ionic_pert, idir)
      
      call dsternheimer_solve(sh, sys, hm, lr, 1, M_ZERO, ionic_pert, &
        VIB_MODES_DIR, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))
      
      do jmat = imat, vib%num_modes 
        jatom = vibrations_get_atom(vib, jmat)
        jdir  = vibrations_get_dir (vib, jmat)

        call pert_setup_atom(ionic_pert, jatom, iatom)
        call pert_setup_dir(ionic_pert, jdir, idir)
        
        vib%dyn_matrix(imat, jmat) = vib%dyn_matrix(imat, jmat) &
          -M_TWO*dpert_expectation_value(ionic_pert, gr, geo, hm, st, st%dpsi, lr(1)%ddl_psi)
        
        vib%dyn_matrix(jmat, imat) = vib%dyn_matrix(imat, jmat)

      end do
      
      ! this part is not correct for periodic systems, a different electric perturbation is required with d/dk
      do jdir = 1, ndim
        call pert_setup_dir(electric_pert, jdir)
        infrared(imat, jdir) = &
          dpert_expectation_value(electric_pert, gr, geo, hm, st, lr(1)%ddl_psi, st%dpsi) + &
          dpert_expectation_value(electric_pert, gr, geo, hm, st, st%dpsi, lr(1)%ddl_psi)
      end do

      message(1) = ""
      call messages_info(1)
    end do 

    call pert_end(ionic_pert)
    call pert_end(electric_pert)

    call vibrations_normalize_dyn_matrix(vib, geo)
    call vibrations_diag_dyn_matrix(vib)
    call vibrations_output(vib, "_lr")
    call calc_infrared

    message(1) = "Calculating response wavefunctions for normal modes."
    call messages_info(1)
    call vib_modes_wavefunctions

    !DESTRUCT

    SAFE_DEALLOCATE_A(infrared)
    call lr_dealloc(lr(1))
    call vibrations_end(vib)
    call sternheimer_end(sh)
    call states_deallocate_wfns(st)

    POP_SUB(phonons_lr_run)

  contains

    ! ---------------------------------------------------------
    ! this formulation is only valid for finite systems, or an Ewald sum is required
    ! as in Baroni et al. RMP 2001, Appendix B.
    subroutine build_ionic_dyn_matrix()
      FLOAT :: ac, xi(1:MAX_DIM), xj(1:MAX_DIM), xk(1:MAX_DIM), r2
      integer :: katom

      PUSH_SUB(phonons_lr_run.build_ionic_dyn_matrix)

      do iatom = 1, natoms
        do idir = 1, gr%mesh%sb%dim

          do jatom = 1, natoms
            do jdir = 1, gr%mesh%sb%dim         

              xi(1:MAX_DIM) = geo%atom(iatom)%x(1:MAX_DIM)

              !ion - ion
              if( iatom == jatom) then 

                ac = M_ZERO
                do katom = 1, natoms
                  if ( katom == iatom ) cycle

                  xk(1:MAX_DIM) = geo%atom(katom)%x(1:MAX_DIM)
                  r2 = sum((xi(1:gr%mesh%sb%dim) - xk(1:gr%mesh%sb%dim))**2)

                  ac = ac + species_zval(geo%atom(iatom)%spec) * &
                            species_zval(geo%atom(katom)%spec) &
                       /(r2**CNST(1.5)) *(&
                       -ddelta(idir, jdir) + &
                       (M_THREE*(xi(idir)-xk(idir))*(xi(jdir)-xk(jdir)))/r2 &
                       )

                end do

              else ! iatom /= jatom

                xj(1:MAX_DIM) = geo%atom(jatom)%x(1:MAX_DIM)

                r2 = sum((xi(1:gr%mesh%sb%dim) - xj(1:gr%mesh%sb%dim))**2)
                ac = species_zval(geo%atom(iatom)%spec) * species_zval(geo%atom(jatom)%spec) &
                     /(r2**CNST(1.5))*(&
                     ddelta(idir, jdir) - (M_THREE*(xi(idir)-xj(idir))*(xi(jdir)-xj(jdir)))/r2)

              end if


              vib%dyn_matrix(vibrations_get_index(vib, iatom, idir), vibrations_get_index(vib, jatom, jdir)) = -ac

            end do
          end do
        end do
      end do
      POP_SUB(phonons_lr_run.build_ionic_dyn_matrix)

    end subroutine build_ionic_dyn_matrix

    ! ---------------------------------------------------------
    subroutine calc_infrared
      FLOAT :: lir(1:MAX_DIM+1)

      PUSH_SUB(phonons_lr_run.calc_infrared)
      !calculate infrared intensities

      iunit = io_open(VIB_MODES_DIR//'infrared', action='write')

      write(iunit, '(a)') '#   freq ['//trim(units_abbrev(unit_invcm))//']     <x>           <y>           <z>           average'

      do iatom = 1, natoms
        do idir = 1, ndim

          imat = vibrations_get_index(vib, iatom, idir)

          do jdir = 1, ndim
            lir(jdir) = dot_product(infrared(:, jdir), vib%normal_mode(:, imat))
          end do
          lir(ndim+1) = sqrt(sum(lir(1:ndim)**2))

          write(iunit, '(5f14.5)') units_from_atomic(unit_invcm, sqrt(abs(vib%freq(imat)))), lir(1:ndim+1)
        end do
      end do

      call io_close(iunit)
      POP_SUB(phonons_lr_run.calc_infrared)
    end subroutine calc_infrared

    ! ---------------------------------------------------------
    subroutine vib_modes_wavefunctions
      ! now calculate the wavefunction associated with each normal mode
      type(lr_t) :: lrtmp
      integer :: ik, ist, idim, inm

      PUSH_SUB(phonons_lr_run.vib_modes_wavefunctions)

      call lr_init(lrtmp)
      call lr_allocate(lrtmp, st, gr%mesh)

      lr(1)%ddl_psi = M_ZERO

      do inm = 1, vib%num_modes

        do iatom = 1, natoms
          do idir = 1, ndim

            imat = vibrations_get_index(vib, iatom, idir)

            call restart_read(trim(restart_dir)//VIB_MODES_DIR//trim(phn_wfs_tag(iatom, idir)), &
                 st, gr, geo, ierr, lr = lrtmp)
            
            do ik = 1, st%d%nik
              do ist = st%st_start, st%st_end
                do idim = 1, st%d%dim

                  call lalg_axpy(gr%mesh%np, vib%normal_mode(imat, inm), &
                    lrtmp%ddl_psi(:, idim, ist, ik), lr(1)%ddl_psi(:, idim, ist, ik))
                  
                end do
              end do
            end do

          end do
        end do
        
        call restart_write(io_workpath(trim(tmpdir)//VIB_MODES_DIR//trim(phn_nm_wfs_tag(inm))), &
          st, gr, ierr, lr = lr(1))
        
      end do

      call lr_dealloc(lrtmp)
      POP_SUB(phonons_lr_run.vib_modes_wavefunctions)
    end subroutine vib_modes_wavefunctions

  end subroutine phonons_lr_run


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
    str = wfs_tag_sigma(str, 1)

    POP_SUB(phn_wfs_tag)
    
  end function phn_wfs_tag


  ! ---------------------------------------------------------
  character(len=100) function phn_nm_wfs_tag(inm) result(str)
    integer, intent(in) :: inm

    PUSH_SUB(phn_wfs_tag)

    write(str, '(a,i5.5)') "phn_nm_wfs_", inm
    str = wfs_tag_sigma(str, 1)

    POP_SUB(phn_wfs_tag)
    
  end function phn_nm_wfs_tag

end module phonons_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
