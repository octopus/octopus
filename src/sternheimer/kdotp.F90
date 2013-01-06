!!! Copyright (C) 2008-2010 David Strubbe
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
  use output_m
  use hamiltonian_m
  use io_m
  use kdotp_calc_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use linear_response_m
  use linear_solver_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use mpi_m
  use parser_m
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
  use utils_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       kdotp_lr_run,       &
       int2str      

  type kdotp_t
    type(pert_t) :: perturbation
    type(pert_t) :: perturbation2

    FLOAT, pointer :: eff_mass_inv(:,:,:,:)  !< inverse effective-mass tensor
                                             !! (idir1, idir2, ist, ik)
    FLOAT, pointer :: velocity(:,:,:) !< group velocity vector (idir, ist, ik)

    type(lr_t), pointer :: lr(:,:) !< linear response for (sys%gr%sb%periodic_dim,1)
                                   !! second index is dummy; should only be 1
                                   !! for compatibility with em_resp routines

    type(lr_t), pointer :: lr2(:,:,:) !< second-order response for 
                                      !! (sys%gr%sb%periodic_dim,sys%gr%sb%periodic_dim,1)

    logical :: ok                   !< is converged?
    integer :: occ_solution_method  !< how to get occupied components of response
    FLOAT   :: degen_thres          !< maximum energy difference to be considered degenerate
    FLOAT   :: eta                  !< imaginary freq. added to Sternheimer eqn.
  end type kdotp_t

contains

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run(sys, hm, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromScratch

    type(kdotp_t)           :: kdotp_vars
    type(sternheimer_t)     :: sh, sh2
    logical                 :: calc_eff_mass, calc_2nd_order, complex_response

    integer              :: idir, idir2, ierr, pdim, ispin
    character(len=100)   :: dirname, str_tmp
    real(8)              :: errornorm

    PUSH_SUB(kdotp_lr_run)

    call messages_experimental("k.p perturbation and calculation of effective masses")

    if(hm%theory_level == HARTREE_FOCK) then
      call messages_not_implemented('Commutator of Fock operator')
    endif

    pdim = sys%gr%sb%periodic_dim

    if(.not. simul_box_is_periodic(sys%gr%sb)) then
       message(1) = "k.p perturbation cannot be used for a finite system."
       call messages_fatal(1)
    endif

    SAFE_ALLOCATE(kdotp_vars%eff_mass_inv(1:pdim, 1:pdim, 1:sys%st%nst, 1:sys%st%d%nik))
    SAFE_ALLOCATE(kdotp_vars%velocity(1:pdim, 1:sys%st%nst, 1:sys%st%d%nik))
    kdotp_vars%eff_mass_inv(:,:,:,:) = 0 
    kdotp_vars%velocity(:,:,:) = 0 

    call pert_init(kdotp_vars%perturbation, PERTURBATION_KDOTP, sys%gr, sys%geo)
    SAFE_ALLOCATE(kdotp_vars%lr(1:pdim, 1:1))

    call parse_input()

    if(calc_2nd_order) then
      call pert_init(kdotp_vars%perturbation2, PERTURBATION_NONE, sys%gr, sys%geo)
      call pert_setup_dir(kdotp_vars%perturbation2, 1) ! direction is irrelevant
      SAFE_ALLOCATE(kdotp_vars%lr2(1:pdim, 1:pdim, 1:1))
    endif

    complex_response = (kdotp_vars%eta /= M_ZERO ) .or. states_are_complex(sys%st)
    call restart_look_and_read(sys%st, sys%gr, is_complex = complex_response)

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response.'
    call messages_info(1)
    call system_h_setup(sys, hm)
    
    if(states_are_real(sys%st)) then
      message(1) = 'Info: Using real wavefunctions.'
    else
      message(1) = 'Info: Using complex wavefunctions.'
    end if
    call messages_info(1)

    message(1) = 'Calculating band velocities.'
    call messages_info(1)

    if(states_are_real(sys%st)) then
      kdotp_vars%velocity(:,:,:) = M_ZERO
    else
      call zcalc_band_velocity(sys, hm, kdotp_vars%perturbation, kdotp_vars%velocity(:,:,:))
    endif

    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(trim(tmpdir)//KDOTP_DIR, is_tmp=.true.) ! restart
      call io_mkdir(KDOTP_DIR) ! data output
      call kdotp_write_band_velocity(sys%st, pdim, kdotp_vars%velocity(:,:,:))
    endif

    call sternheimer_obsolete_variables('KdotP_', 'KdotP')
    call sternheimer_init(sh, sys, hm, 'KdotP', complex_response, set_ham_var = 0, &
      set_occ_response = (kdotp_vars%occ_solution_method == 0), set_last_occ_response = (kdotp_vars%occ_solution_method == 0))
    ! ham_var_set = 0 results in HamiltonianVariation = V_ext_only
    if(calc_2nd_order) then
      call sternheimer_init(sh2, sys, hm, 'KdotP', complex_response, set_ham_var = 0, &
        set_occ_response = .false., set_last_occ_response = .false.)
    endif

    do idir = 1, pdim
      call lr_init(kdotp_vars%lr(idir, 1))
      call lr_allocate(kdotp_vars%lr(idir, 1), sys%st, sys%gr%mesh)

      if(calc_2nd_order) then
        do idir2 = 1, pdim
          call lr_init(kdotp_vars%lr2(idir, idir2, 1))
          call lr_allocate(kdotp_vars%lr2(idir, idir2, 1), sys%st, sys%gr%mesh)
        enddo
      endif

      ! load wavefunctions
      if(.not. fromScratch) then
        str_tmp = kdotp_wfs_tag(idir)
        write(dirname,'(2a)') KDOTP_DIR, trim(wfs_tag_sigma(str_tmp, 1))
        call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, ierr, lr=kdotp_vars%lr(idir, 1))
          
        if(ierr .ne. 0) then
          message(1) = "Could not load response wavefunctions from '"//trim(tmpdir)//trim(dirname)//"'"
          call messages_warning(1)
        end if

        if(calc_2nd_order) then
          do idir2 = 1, pdim
            str_tmp = kdotp_wfs_tag(idir, idir2)
            write(dirname,'(2a)') KDOTP_DIR, trim(wfs_tag_sigma(str_tmp, 1))
            call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, ierr, lr=kdotp_vars%lr2(idir, idir2, 1))
          
            if(ierr .ne. 0) then
              message(1) = "Could not load response wavefunctions from '"//trim(tmpdir)//trim(dirname)//"'"
              call messages_warning(1)
            end if
          enddo
        endif
      end if
    end do

    call info()
    message(1) = "Info: Calculating k.p linear response of ground-state wavefunctions."
    call messages_info(1)
    kdotp_vars%ok = .true.

    ! solve the Sternheimer equation
    do idir = 1, pdim
      write(message(1), '(3a)') 'Info: Calculating response for the ', index2axis(idir), &
                                '-direction.' 
      call messages_info(1)
      call pert_setup_dir(kdotp_vars%perturbation, idir)

      if(states_are_real(sys%st)) then
        call dsternheimer_solve(sh, sys, hm, kdotp_vars%lr(idir,:), 1, &
          M_ZERO, kdotp_vars%perturbation, KDOTP_DIR, &
          "", kdotp_wfs_tag(idir), have_restart_rho = .false.)
      else
        call zsternheimer_solve(sh, sys, hm, kdotp_vars%lr(idir,:), 1, &
          M_zI * kdotp_vars%eta, kdotp_vars%perturbation, KDOTP_DIR, &
          "", kdotp_wfs_tag(idir), have_restart_rho = .false.)
      endif

      kdotp_vars%ok = kdotp_vars%ok .and. sternheimer_has_converged(sh)         

      errornorm = M_ZERO
      if(states_are_real(sys%st)) then 
        call doutput_lr(sys%st, sys%gr, kdotp_vars%lr(idir, 1), KDOTP_DIR, idir, 1, sys%outp, sys%geo, units_out%force)

        do ispin = 1, sys%st%d%nspin
          errornorm = hypot(errornorm, real(dmf_nrm2(sys%gr%mesh, kdotp_vars%lr(idir, 1)%ddl_rho(:, ispin)), 8))
        end do
      else
        call zoutput_lr(sys%st, sys%gr, kdotp_vars%lr(idir, 1), KDOTP_DIR, idir, 1, sys%outp, sys%geo, units_out%force)

        do ispin = 1, sys%st%d%nspin
          errornorm = hypot(errornorm, real(zmf_nrm2(sys%gr%mesh, kdotp_vars%lr(idir, 1)%zdl_rho(:, ispin)), 8))
        end do
      endif

      write(message(1),'(a,f12.6)') "Norm of relative density variation = ", errornorm / sys%st%qtot
      call messages_info(1)

      ! TODO: use equality of mixed partials to calculate fewer of these
      if(calc_2nd_order) then
        do idir2 = 1, pdim
          write(message(1), '(3a)') 'Info: Calculating second-order response in the ', index2axis(idir2), &
            '-direction.' 
          call messages_info(1)

          if(states_are_real(sys%st)) then
            call dsternheimer_solve_order2(sh, sh, sh2, sys, hm, kdotp_vars%lr(idir, 1:1), kdotp_vars%lr(idir, 1:1), &
              1, M_ZERO, M_ZERO, kdotp_vars%perturbation, kdotp_vars%perturbation, &
              kdotp_vars%lr2(idir, idir2, 1:1), kdotp_vars%perturbation2, KDOTP_DIR, "", kdotp_wfs_tag(idir, idir2), &
              have_restart_rho = .false., have_exact_freq = .true.)
          else
            call zsternheimer_solve_order2(sh, sh, sh2, sys, hm, kdotp_vars%lr(idir, 1:1), kdotp_vars%lr(idir, 1:1), &
              1, M_zI * kdotp_vars%eta, M_zI * kdotp_vars%eta, kdotp_vars%perturbation, kdotp_vars%perturbation, &
              kdotp_vars%lr2(idir, idir2, 1:1), kdotp_vars%perturbation2, KDOTP_DIR, "", kdotp_wfs_tag(idir, idir2), &
              have_restart_rho = .false., have_exact_freq = .true.)
          endif

        enddo
        message(1) = ""
        call messages_info(1)
      endif
    end do ! idir

    ! calculate effective masses
    if (calc_eff_mass) then
      message(1) = "Info: Calculating effective masses."
      call messages_info(1)

      if(states_are_real(sys%st)) then
        call dcalc_eff_mass_inv(sys, hm, kdotp_vars%lr, kdotp_vars%perturbation, &
          kdotp_vars%eff_mass_inv, kdotp_vars%occ_solution_method, kdotp_vars%degen_thres)
      else
        call zcalc_eff_mass_inv(sys, hm, kdotp_vars%lr, kdotp_vars%perturbation, &
          kdotp_vars%eff_mass_inv, kdotp_vars%occ_solution_method, kdotp_vars%degen_thres)
      endif

      call kdotp_write_degeneracies(sys%st, kdotp_vars%degen_thres)
      call kdotp_write_eff_mass(sys%st, sys%gr, kdotp_vars)
    endif

    ! clean up some things
    do idir = 1, pdim
      call lr_dealloc(kdotp_vars%lr(idir, 1))

      if(calc_2nd_order) then
        do idir2 = 1, pdim
          call lr_dealloc(kdotp_vars%lr2(idir, idir2, 1))
        enddo
      endif
    end do

    call sternheimer_end(sh)
    call pert_end(kdotp_vars%perturbation)
    SAFE_DEALLOCATE_P(kdotp_vars%lr)

    if(calc_2nd_order) then
      call sternheimer_end(sh2)
      call pert_end(kdotp_vars%perturbation2)
      SAFE_DEALLOCATE_P(kdotp_vars%lr2)
    endif

    call states_deallocate_wfns(sys%st)
    SAFE_DEALLOCATE_P(kdotp_vars%eff_mass_inv)
    SAFE_DEALLOCATE_P(kdotp_vars%velocity)

    POP_SUB(kdotp_lr_run)

  contains

    ! ---------------------------------------------------------

    subroutine parse_input()

      PUSH_SUB(kdotp_lr_run.parse_input)

      !%Variable KdotPOccupiedSolutionMethod
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

      call messages_obsolete_variable('KdotP_OccupiedSolutionMethod', 'KdotPOccupiedSolutionMethod')

      call parse_integer(datasets_check('KdotPOccupiedSolutionMethod'), 0, kdotp_vars%occ_solution_method)

      call parse_float(datasets_check('DegeneracyThreshold'), &
        units_from_atomic(units_inp%energy, CNST(1e-5)), kdotp_vars%degen_thres)
      kdotp_vars%degen_thres = units_to_atomic(units_inp%energy, kdotp_vars%degen_thres)
      ! Note: this variable is defined in src/states_calc.F90, in states_degeneracy_matrix

      !%Variable KdotPEta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::KdotP
      !%Description
      !% Imaginary frequency added to Sternheimer equation which may improve convergence.
      !% Not recommended.
      !%End
      call messages_obsolete_variable('KdotP_Eta', 'KdotPEta')
      call parse_float(datasets_check('KdotPEta'), M_ZERO, kdotp_vars%eta)
      kdotp_vars%eta = units_to_atomic(units_inp%energy, kdotp_vars%eta)

      !%Variable KdotPCalculateEffectiveMasses
      !%Type logical
      !%Default true
      !%Section Linear Response::KdotP
      !%Description
      !% If true, uses <tt>kdotp</tt> perturbations of ground-state wavefunctions
      !% to calculate effective masses.
      !%End      
      call messages_obsolete_variable('KdotP_CalculateEffectiveMasses', 'KdotPCalculateEffectiveMasses')
      call parse_logical(datasets_check('KdotPCalculateEffectiveMasses'), .true., calc_eff_mass)

      !%Variable KdotPCalcSecondOrder
      !%Type logical
      !%Default false
      !%Section Linear Response::KdotP
      !%Description
      !% If true, calculates second-order response of wavefunctions as well as first-order response.
      !% Note that the second derivative of the Hamiltonian is NOT included in this calculation.
      !%End      
      call parse_logical(datasets_check('KdotPCalcSecondOrder'), .false., calc_2nd_order)

      POP_SUB(kdotp_lr_run.parse_input)

   end subroutine parse_input

    ! ---------------------------------------------------------
    subroutine info()

      PUSH_SUB(kdotp_lr_run.info)

      call pert_info(kdotp_vars%perturbation, stdout)

      write(message(1),'(a)') 'k.p perturbation theory'
      call messages_print_stress(stdout, trim(message(1)))

      if (kdotp_vars%occ_solution_method == 0) then
        message(1) = 'Occupied solution method: Sternheimer equation.'
      else
        message(1) = 'Occupied solution method: sum over states.'
      endif

      call messages_info(1)

      call messages_print_stress(stdout)
      
      POP_SUB(kdotp_lr_run.info)

    end subroutine info

  end subroutine kdotp_lr_run

  ! ---------------------------------------------------------
  subroutine kdotp_write_band_velocity(st, periodic_dim, velocity)
    type(states_t), intent(inout) :: st
    integer,        intent(in)    :: periodic_dim
    FLOAT,          intent(in)    :: velocity(:,:,:)

    character(len=80) :: filename, tmp
    integer :: iunit, ik, ist, ik2, ispin, idir

    PUSH_SUB(kdotp_write_band_velocity)

    write(filename, '(a)') KDOTP_DIR//'velocity'
    iunit = io_open(trim(filename), action='write')
    write(iunit,'(a)') '# Band velocities'

    do ik = 1, st%d%nik
      ispin = states_dim_get_spin_index(st%d, ik)
      ik2 = states_dim_get_kpoint_index(st%d, ik)
      tmp = int2str(ik2)

      write(iunit,'(a,i1,a,a)') '# spin = ', ispin, ', k-point = ', trim(tmp)

      write(iunit,'(a)',advance='no') '# state    energy       '
      do idir = 1, periodic_dim
        write(iunit,'(3a)',advance='no') 'vg(', trim(index2axis(idir)), ')       '
      enddo
      write(iunit,'(a)')

      write(iunit,'(3a)',advance='no')       '#           [', trim(units_abbrev(units_out%energy)), ']     '
      do idir = 1, periodic_dim
        write(iunit,'(3a)',advance='no') '[', trim(units_abbrev(units_out%velocity)), '] '
      enddo
      write(iunit,'(a)')

      do ist = 1, st%nst
        write(iunit,'(i5,f12.5,3f12.5)') ist, units_from_atomic(units_out%energy, st%eigenval(ist, ik)), &
          velocity(1:periodic_dim, ist, ik)
      enddo
    enddo

    call io_close(iunit)
    POP_SUB(kdotp_write_band_velocity)
  end subroutine kdotp_write_band_velocity

  ! ---------------------------------------------------------
  subroutine kdotp_write_eff_mass(st, gr, kdotp_vars)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(kdotp_t),        intent(inout) :: kdotp_vars

    character(len=80) :: filename, tmp
    integer :: iunit, ik, ist, ik2, ispin
    FLOAT :: determinant

    PUSH_SUB(kdotp_write_eff_mass)

    do ik = 1, st%d%nik
      ispin = states_dim_get_spin_index(st%d, ik)
      ik2 = states_dim_get_kpoint_index(st%d, ik)

      tmp = int2str(ik2)
      write(filename, '(3a, i1)') KDOTP_DIR//'kpoint_', trim(tmp), '_', ispin
      iunit = io_open(trim(filename), action='write')

      write(iunit,'(a, i10)')    '# spin    index = ', ispin
      write(iunit,'(a, i10)')    '# k-point index = ', ik2
      write(iunit,'(a, 99f12.8)') '# k-point coordinates = ', kpoints_get_point(gr%sb%kpoints, ik2)
      if (.not. kdotp_vars%ok) write(iunit, '(a)') "# WARNING: not converged"      
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Inverse effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        call output_tensor(iunit, kdotp_vars%eff_mass_inv(:, :, ist, ik), gr%sb%periodic_dim, unit_one)
      enddo
      
      write(iunit,'(a)')
      write(iunit,'(a)') '# Effective-mass tensors'
      do ist = 1, st%nst
        write(iunit,'(a)')
        tmp = int2str(ist)
        write(iunit,'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
          units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
        determinant = lalg_inverter(gr%sb%periodic_dim, kdotp_vars%eff_mass_inv(:, :, ist, ik), .true.)
        call output_tensor(iunit, kdotp_vars%eff_mass_inv(:, :, ist, ik), gr%sb%periodic_dim, unit_one)
      enddo

      call io_close(iunit)
    enddo

    POP_SUB(kdotp_write_eff_mass)
  end subroutine kdotp_write_eff_mass

  ! ---------------------------------------------------------
  subroutine kdotp_write_degeneracies(st, threshold)
    type(states_t), intent(inout) :: st
    FLOAT,          intent(in)    :: threshold

    character(len=80) :: tmp
    integer :: ik, ist, ist2, ik2, ispin

    PUSH_SUB(kdotp_write_degeneracies)

    call messages_print_stress(stdout, 'Degenerate subspaces')

    do ik = 1, st%d%nik
      ispin = states_dim_get_spin_index(st%d, ik)
      ik2 = states_dim_get_kpoint_index(st%d, ik)

      tmp = int2str(ik2)
      write(message(1), '(3a, i1)') 'k-point ', trim(tmp), ', spin ', ispin 
      call messages_info(1)

      ist = 1
      do while (ist <= st%nst)
      ! test for degeneracies
         write(message(1),'(a)') '===='
         tmp = int2str(ist)
         write(message(2),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
           units_from_atomic(units_out%energy, st%eigenval(ist, ik)), ' ', units_abbrev(units_out%energy)
         call messages_info(2)

         ist2 = ist + 1
         do while (ist2 <= st%nst .and. &
           abs(st%eigenval(min(ist2, st%nst), ik) - st%eigenval(ist, ik)) < threshold)
           tmp = int2str(ist2)
           write(message(1),'(a, a, a, f12.8, a, a)') 'State #', trim(tmp), ', Energy = ', &
             units_from_atomic(units_out%energy, st%eigenval(ist2, ik)), ' ', units_abbrev(units_out%energy)
           call messages_info(1)
           ist2 = ist2 + 1
         enddo

         ist = ist2
      enddo

      write(message(1),'()')
      call messages_info(1)
    enddo

    message(1) = "Velocities and effective masses are not correct within degenerate subspaces."
    call messages_warning(1)

    POP_SUB(kdotp_write_degeneracies)
  end subroutine kdotp_write_degeneracies

  ! ---------------------------------------------------------
  character(len=12) pure function int2str(ii) result(str)
    integer, intent(in) :: ii
    
    write(str, '(i11)') ii
    str = trim(adjustl(str))

  end function int2str
            
end module kdotp_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
