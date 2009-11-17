!!! Copyright (C) 2008-2009 David Strubbe
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

module kdotp_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use h_sys_output_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kdotp_calc_m
  use lalg_adv_m
  use lalg_basic_m
  use linear_response_m
  use linear_solver_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use pert_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_calc_m
  use states_dim_m
  use sternheimer_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       kdotp_lr_run,       &
       int2str      

  type kdotp_t
    type(pert_t) :: perturbation

    FLOAT, pointer :: eff_mass_inv(:, :, :, :)  ! inverse effective-mass tensor
                                                ! (ik, ist, idir1, idir2)

    type(lr_t), pointer :: lr(:,:) ! linear response for (gr%mesh%sb%dim,1)
                                   ! second index is dummy; should only be 1
                                   ! for compatibility with em_resp routines

    logical :: ok                   ! is converged?
    integer :: occ_solution_method  ! how to get occupied components of response
    FLOAT   :: degen_thres          ! maximum energy difference to be considered
                                    ! degenerate
    FLOAT   :: eta                  ! imaginary freq. added to Sternheimer eqn.
  end type kdotp_t

contains

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(kdotp_t)           :: kdotp_vars
    type(sternheimer_t)     :: sh
    logical                 :: calc_eff_mass, complex_response

    integer              :: idir, ierr, is, ist, ik
    character(len=100)   :: dirname, str_tmp
    logical, allocatable :: orth_mask(:)

    call push_sub('kdotp.kdotp_lr_run')

    call messages_devel_version("k.p perturbation and calculation of effective masses")

    gr => sys%gr

    if(.not. simul_box_is_periodic(gr%sb)) then
       message(1) = "k.p perturbation cannot be used for a finite system."
       call write_fatal(1)
    endif

    SAFE_ALLOCATE(kdotp_vars%eff_mass_inv(1:sys%st%d%nik, 1:sys%st%nst, 1:gr%mesh%sb%dim, 1:gr%mesh%sb%dim))

    kdotp_vars%eff_mass_inv(:,:,:,:) = 0  

    call pert_init(kdotp_vars%perturbation, PERTURBATION_KDOTP, sys%gr, sys%geo)

    SAFE_ALLOCATE(kdotp_vars%lr(1:gr%mesh%sb%dim, 1:1))

    call parse_input()

    complex_response = (kdotp_vars%eta /= M_ZERO ) .or. states_are_complex(sys%st)
    call restart_look_and_read(sys%st, sys%gr, sys%geo, is_complex = complex_response)

    kdotp_vars%lr(1:gr%mesh%sb%dim, 1:1)%nst = sys%st%nst

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call write_info(1)
    call system_h_setup(sys, hm)
    
    if (states_are_real(sys%st)) then
      message(1) = 'Info: SCF using real wavefunctions.'
    else
      message(1) = 'Info: SCF using complex wavefunctions.'
    end if
    call write_info(1)

    call sternheimer_init(sh, sys, hm, "KdotP_", &
         set_ham_var = 0, set_occ_response = (kdotp_vars%occ_solution_method == 0), default_solver = LS_QMR_DOTP)
    ! ham_var_set = 0 results in HamiltonianVariation = V_ext_only

    do idir = 1, gr%mesh%sb%dim
      call lr_init(kdotp_vars%lr(idir, 1))
      call lr_allocate(kdotp_vars%lr(idir, 1), sys%st, sys%gr%mesh)

      ! load wavefunctions
      if(.not.fromScratch) then
        str_tmp =  kdotp_wfs_tag(idir)
        write(dirname,'(3a)') KDOTP_DIR, trim(str_tmp), '_1'
        ! 1 is the sigma index which is used in em_resp
        call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
          ierr, lr=kdotp_vars%lr(idir, 1))
          
        if(ierr.ne.0) then
          message(1) = "Could not load response wavefunctions from '"//trim(tmpdir)//trim(dirname)//"'"
          call write_warning(1)
        end if      
      end if
    end do

    call io_mkdir(trim(tmpdir)//KDOTP_DIR) ! restart
    call io_mkdir(KDOTP_DIR)               ! data output
    call info()
    message(1) = "Info: Calculating k.p linear response of ground-state wavefunctions."
    call write_info(1)
    kdotp_vars%ok = .true.

    ! solve the Sternheimer equation
    do idir = 1, gr%mesh%sb%dim
      write(message(1), '(3a)') 'Info: Calculating response for the ', index2axis(idir), &
                                '-direction.' 
      call write_info(1)
      call pert_setup_dir(kdotp_vars%perturbation, idir)

      if(states_are_real(sys%st)) then
        call dsternheimer_solve(sh, sys, hm, kdotp_vars%lr(idir,:), 1, &
          M_ZERO, kdotp_vars%perturbation, KDOTP_DIR, &
          kdotp_rho_tag(idir), kdotp_wfs_tag(idir), have_restart_rho=(ierr==0))
      else
        call zsternheimer_solve(sh, sys, hm, kdotp_vars%lr(idir,:), 1, &
          M_zI * kdotp_vars%eta, kdotp_vars%perturbation, KDOTP_DIR, &
          kdotp_rho_tag(idir), kdotp_wfs_tag(idir), have_restart_rho=(ierr==0))
      endif

      kdotp_vars%ok = kdotp_vars%ok .and. sternheimer_has_converged(sh)         
    end do ! idir

    ! calculate effective masses
    if (calc_eff_mass) then
      message(1) = "Info: Calculating effective masses."
      call write_info(1)

      if(states_are_real(sys%st)) then
        call dcalc_eff_mass_inv(sys, hm, kdotp_vars%lr, kdotp_vars%perturbation, &
          kdotp_vars%eff_mass_inv, kdotp_vars%occ_solution_method, kdotp_vars%degen_thres)
      else
        call zcalc_eff_mass_inv(sys, hm, kdotp_vars%lr, kdotp_vars%perturbation, &
          kdotp_vars%eff_mass_inv, kdotp_vars%occ_solution_method, kdotp_vars%degen_thres)
      endif
      call kdotp_output(sys%st, sys%gr, kdotp_vars)
    endif

    ! clean up some things
    do idir = 1, gr%mesh%sb%dim
      call lr_dealloc(kdotp_vars%lr(idir, 1))
    end do

    call sternheimer_end(sh)
    call pert_end(kdotp_vars%perturbation)

    SAFE_DEALLOCATE_P(kdotp_vars%lr)
    call states_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(kdotp_vars%eff_mass_inv)

    call pop_sub()

  contains

    ! ---------------------------------------------------------

    subroutine parse_input()

      call push_sub('kdotp.kdotp_lr_run.parse_input')

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

      call parse_integer(datasets_check('KdotP_OccupiedSolutionMethod'), &
        0, kdotp_vars%occ_solution_method)

      call parse_float(datasets_check('DegeneracyThreshold'), &
        units_from_atomic(units_inp%energy, CNST(1e-5)), kdotp_vars%degen_thres)
      kdotp_vars%degen_thres = units_to_atomic(units_inp%energy, kdotp_vars%degen_thres)
      ! Note: this variable is defined in src/states_calc.F90, in states_degeneracy_matrix

      !%Variable KdotP_Eta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::KdotP
      !%Description
      !% Imaginary frequency added to Sternheimer equation which may improve convergence.
      !% Not recommended.
      !%End

      call parse_float(datasets_check('KdotP_Eta'), M_ZERO, kdotp_vars%eta)
      kdotp_vars%eta = units_to_atomic(units_inp%energy, kdotp_vars%eta)

      !%Variable KdotP_CalculateEffectiveMasses
      !%Type logical
      !%Default true
      !%Section Linear Response::KdotP
      !%Description
      !% If true, uses <tt>kdotp</tt> perturbations of ground-state wavefunctions
      !% to calculate effective masses.
      !%End      

      call parse_logical(datasets_check('KdotP_CalculateEffectiveMasses'), &
        .true., calc_eff_mass)

      call pop_sub()

   end subroutine parse_input

    ! ---------------------------------------------------------
    subroutine info()

      call push_sub('kdotp.kdotp_lr_run.info')

      call pert_info(kdotp_vars%perturbation, stdout)

      write(message(1),'(a)') 'k.p perturbation theory'
      call messages_print_stress(stdout, trim(message(1)))

      if (kdotp_vars%occ_solution_method == 0) then
        message(1) = 'Occupied solution method: Sternheimer equation.'
      else
        message(1) = 'Occupied solution method: sum over states.'
      endif

      call write_info(1)

      call messages_print_stress(stdout)
      
      call pop_sub()

    end subroutine info

  end subroutine kdotp_lr_run

  ! ---------------------------------------------------------
  subroutine kdotp_output(st, gr, kdotp_vars)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(kdotp_t),        intent(inout) :: kdotp_vars

    character(len=80) :: filename, tmp
    integer :: iunit, ik, ist, ist2, ik2, ispin
    FLOAT :: determinant

    call push_sub('kdotp.kdotp_output')

    call messages_print_stress(stdout, 'Degenerate subspaces')

    do ik = 1, st%d%nik
      ispin = states_dim_get_spin_index(st%d, ik)
      ik2 = kpoint_index(st%d, ik)

      tmp = int2str(ik2)
      write(message(1), '(3a, i1)') 'k-point ', trim(tmp), ', spin ', ispin 
      call write_info(1)

      ist = 1
      do while (ist <= st%nst)
      ! test for degeneracies
         write(message(1),'(a)') '===='
         tmp = int2str(ist)
         write(message(2),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
           units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
         call write_info(2)

         ist2 = ist + 1
         do while (ist2 <= st%nst .and. &
           abs(st%eigenval(min(ist2, st%nst), ik) - st%eigenval(ist, ik)) < kdotp_vars%degen_thres)
           ! eigenvalues are supposed to be in ascending order; if they are not, it is a sign
           ! of being in a degenerate subspace?
           ! write(*,*) ist2, ist, sys%st%eigenval(ist2, ik) - sys%st%eigenval(ist, ik)
            tmp = int2str(ist2)
            write(message(1),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
              units_from_atomic(units_out%energy, st%eigenval(ist2, ik)), ' ', units_abbrev(units_out%energy)
            call write_info(1)
            ist2 = ist2 + 1
         enddo

         ist = ist2
      enddo
      write(message(1),'()')
      call write_info(1)

      tmp = int2str(ik2)
      write(filename, '(3a, i1)') KDOTP_DIR//'kpoint_', trim(tmp), '_', ispin
      iunit = io_open(trim(filename), action='write')

      write(iunit,'(a, i10)') '# spin    index = ', ispin
      write(iunit,'(a, i10)') '# k-point index = ', ik2
      write(iunit,'(a, 3f12.8)') '# k-point coordinates = ', st%d%kpoints(1:MAX_DIM, ik)
      if (.not. kdotp_vars%ok) write(iunit, '(a)') "# WARNING: not converged"      
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Inverse effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), gr%mesh%sb%dim, unit_one)
      enddo
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        determinant = lalg_inverter(gr%sb%dim, kdotp_vars%eff_mass_inv(ik, ist, :, :), .true.)
        call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), gr%mesh%sb%dim, unit_one)
      enddo

      call io_close(iunit)
    enddo

    call pop_sub()

  end subroutine kdotp_output

  character(len=12) function int2str(i) result(str)
    integer, intent(in) :: i
    
    call push_sub('kdotp.int2str')

    write(str, '(i11)') i
    str = trim(adjustl(str))

    call pop_sub()
    
  end function int2str
            
end module kdotp_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
