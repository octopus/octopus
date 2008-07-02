!!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: kdotp.F90 4145 2008-05-02 23:29:41Z xavier $

#include "global.h"
#define RESTART_DIR "kdotp/"

module kdotp_lr_m
  use datasets_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kdotp_calc_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use h_sys_output_m
  use pert_m
  use pol_lr_m
  use restart_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       kdotp_lr_run,       &
       int2str      

  type kdotp_t
    type(pert_t) :: perturbation

    FLOAT, pointer :: eff_mass_inv(:, :, :, :)  ! inverse effective mass tensor
                                                ! (ik, ist, idir1, idir2)

    type(lr_t), pointer :: lr(:,:) ! linear response for (NDIM,1)
                                   ! second index is dummy; should only be 1
                                   ! for compatibility with em_resp routines

    logical :: ok                   ! is converged?
    integer :: occ_solution_method  ! how to get occupied components of response
    FLOAT   :: degen_thres          ! maximum energy difference to be considered
                                    ! degenerate
  end type kdotp_t

contains

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(kdotp_t)           :: kdotp_vars
    type(sternheimer_t)     :: sh

    integer :: idir, ierr
    character(len=100) :: dirname, str_tmp

    call push_sub('kdotp.static_kdotp_lr_run')

    gr => sys%gr
!    ndim = sys%gr%sb%dim

    ALLOCATE(kdotp_vars%eff_mass_inv(sys%st%d%nik, sys%st%nst, NDIM, NDIM), sys%st%d%nik * sys%st%nst * NDIM * NDIM)
    kdotp_vars%eff_mass_inv(:,:,:,:)=0  

    call pert_init(kdotp_vars%perturbation, PERTURBATION_KDOTP, sys%gr, sys%geo)

    ALLOCATE(kdotp_vars%lr(1:NDIM, 1), NDIM)

    call parse_input()

    call read_wfs(sys%st, sys%gr, sys%geo, .true.)
    ! even if wfs are real, the response must be allowed to be complex

    kdotp_vars%lr(1:NDIM,:)%nst = sys%st%nst

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    call system_h_setup(sys, h)
    
    call sternheimer_init(sh, sys, h, "KdotP", hermitian = wfs_are_real(sys%st), &
         set_ham_var = 0, set_occ_response = (kdotp_vars%occ_solution_method == 0))
    ! ham_var_set = 0 results in HamiltonianVariation = V_ext_only

    do idir = 1, NDIM
      call lr_init(kdotp_vars%lr(idir, 1))
      call lr_allocate(kdotp_vars%lr(idir, 1), sys%st, sys%gr%m)

      ! load wave-functions
      if(.not.fromScratch) then
         str_tmp =  kdotp_wfs_tag(idir)
         write(dirname,'(3a, i1)') RESTART_DIR, trim(str_tmp), '_1'
         ! 1 is the sigma index which is used in em_resp
         call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
               ierr, lr=kdotp_vars%lr(idir, 1))
          
          if(ierr.ne.0) then
             message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//dirname
             call write_warning(1)
          end if
          
       end if

    end do

    call io_mkdir(trim(tmpdir)//RESTART_DIR)

    call info()

    message(1) = "Info: Calculating effective masses."
    call write_info(1)

    call io_mkdir('kdotp/')

    kdotp_vars%ok = .true.

    do idir = 1, NDIM
      write(message(1), '(a,i3)') 'Info: Calculating response for direction ', idir
      call write_info(1)
      call pert_setup_dir(kdotp_vars%perturbation, idir)
      call zsternheimer_solve(sh, sys, h, kdotp_vars%lr(idir,:), 1, M_Z0, &
        kdotp_vars%perturbation, RESTART_DIR, &
        kdotp_rho_tag(idir), kdotp_wfs_tag(idir), have_restart_rho=(ierr==0))
      kdotp_vars%ok = kdotp_vars%ok .and. sternheimer_has_converged(sh)         
    end do ! idir

    call zlr_calc_eff_mass_inv(sys, h, kdotp_vars%lr, kdotp_vars%perturbation, &
         kdotp_vars%eff_mass_inv, kdotp_vars%occ_solution_method, kdotp_vars%degen_thres)

    call kdotp_output(sys%st, sys%gr, kdotp_vars)

    do idir = 1, NDIM
      call lr_dealloc(kdotp_vars%lr(idir, 1))
    end do

    call sternheimer_end(sh)
    call pert_end(kdotp_vars%perturbation)

    deallocate(kdotp_vars%lr)
    call states_deallocate_wfns(sys%st)
    deallocate(kdotp_vars%eff_mass_inv)

    call pop_sub()

  contains

    ! ---------------------------------------------------------

    subroutine parse_input()

      !%Variable KdotP_OccupiedSolutionMethod
      !%Type integer
      !%Default sternheimer
      !%Section Linear Response::KdotP
      !%Description
      !% Method of calculating the contribution of the projection of the
      !%  linear-response wavefunctions in the occupied subspace.
      !%Option sternheimer_eqn 0
      !% The Sternheimer equation is solved including the occupied subspace,
      !% to get the full linear-response wavefunctions.
      !%Option sum_over_states 1
      !% The Sternheimer equation is solved only in the unoccupied subspace,
      !% and a sum-over-states perturbation-theory expression is used to
      !% evaluate the contributions in the occupied subspace.
      !%End      

      call loct_parse_int(check_inp('KdotP_OccupiedSolutionMethod'), &
           0, kdotp_vars%occ_solution_method)

      call loct_parse_float(check_inp('DegeneracyThreshold'), &
           CNST(1e-5), kdotp_vars%degen_thres)
      ! Note: this variable is defined in src/states.F90, in states_degeneracy_matrix

    end subroutine parse_input

    ! ---------------------------------------------------------
    subroutine info()

      call pert_info(kdotp_vars%perturbation, stdout)

      write(message(1),'(a)') 'Effective masses'
      call messages_print_stress(stdout, trim(message(1)))

      if (kdotp_vars%occ_solution_method == 0) then
          message(1) = 'Occupied solution method: Sternheimer equation'
      else
          message(1) = 'Occupied solution method: sum over states'
      endif

      call write_info(1)

      call messages_print_stress(stdout)

    end subroutine info

  end subroutine kdotp_lr_run

  ! ---------------------------------------------------------
  subroutine kdotp_output(st, gr, kdotp_vars)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(kdotp_t),        intent(inout) :: kdotp_vars

    character(len=80) :: filename
    integer :: iunit, ik, ist
    FLOAT :: determinant

    do ik = 1, st%d%nik
        write(filename, '(a, a)') 'kdotp/kpoint_', int2str(ik)
        iunit = io_open(trim(filename), action='write')
        write(iunit,'(a, a)') '# k-point index = ', int2str(ik)
        write(iunit,'(a, 3f12.8)') '# k-point coordinates = ', st%d%kpoints(1:MAX_DIM, ik)
        if (.not.kdotp_vars%ok) write(iunit, '(a)') "# WARNING: not converged"      

        write(iunit,'(a)')
        write(iunit,'(a)') '# Inverse effective mass tensors'
        do ist = 1, st%nst
          write(iunit,'(a)')
          write(iunit,'(a, a, a, f12.8, a, a)') 'State #', int2str(ist), ', Energy = ', &
            st%eigenval(ist, ik)/units_out%energy%factor, ' ', units_out%energy%abbrev
          call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), NDIM, M_ONE)
        enddo

        write(iunit,'(a)')
        write(iunit,'(a)') '# Effective mass tensors'
        do ist = 1, st%nst
          write(iunit,'(a)')
          write(iunit,'(a, a, a, f12.8, a, a)') 'State #', int2str(ist), ', Energy = ', &
            st%eigenval(ist, ik)/units_out%energy%factor, ' ', units_out%energy%abbrev
          determinant = lalg_inverter(gr%sb%dim, kdotp_vars%eff_mass_inv(ik, ist, :, :), .true.)
          call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), NDIM, M_ONE)
        enddo

     end do

  end subroutine kdotp_output

  character(len=12) function int2str(i) result(str)
    integer, intent(in) :: i
    
    write(str, '(i11)') i
    str = trim(adjustl(str))
    
  end function int2str
            
end module kdotp_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
