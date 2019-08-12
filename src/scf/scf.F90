!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module scf_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use berry_oct_m
  use density_oct_m
  use eigensolver_oct_m
  use energy_calc_oct_m
  use forces_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use kpoints_oct_m
  use lcao_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use lda_u_mixer_oct_m
  use loct_oct_m
  use magnetic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mix_oct_m
  use modelmb_exchange_syms_oct_m
  use mpi_oct_m
  use multigrid_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use partial_charges_oct_m
  use poisson_oct_m
  use preconditioners_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scdm_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_group_oct_m
  use states_elec_io_oct_m
  use states_elec_restart_oct_m
  use stress_oct_m
  use symmetries_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m
  use varinfo_oct_m
  use vdw_ts_oct_m
!  use xc_functl_oct_m
  use walltimer_oct_m
  use XC_F90(lib_m)
  use xc_oep_oct_m
  
  implicit none

  private
  public ::             &
    scf_t,              &
    scf_init,           &
    scf_mix_clear,      &
    scf_run,            &
    scf_end

  integer, public, parameter :: &
    VERB_NO      = 0,   &
    VERB_COMPACT = 1,   &
    VERB_FULL    = 3
  
  !> some variables used for the SCF cycle
  type scf_t
    private
    integer, public :: max_iter   !< maximum number of SCF iterations
    integer :: max_iter_berry  !< max number of electronic iterations before updating density, for Berry potential

    FLOAT, public :: lmm_r

    ! several convergence criteria
    FLOAT :: conv_abs_dens, conv_rel_dens, conv_abs_ev, conv_rel_ev, conv_abs_force
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev, abs_force
    FLOAT :: conv_energy_diff
    FLOAT :: energy_diff
    logical :: conv_eigen_error
    logical :: check_conv

    integer :: mix_field
    logical :: lcao_restricted
    logical :: calc_force
    logical :: calc_stress
    logical :: calc_dipole
    logical :: calc_partial_charges
    type(mix_t) :: smix
    type(mixfield_t), pointer :: mixfield
    type(eigensolver_t) :: eigens
    integer :: mixdim1
    logical :: forced_finish !< remember if 'touch stop' was triggered earlier.
    type(lda_u_mixer_t) :: lda_u_mix
    type(grid_t), pointer :: gr 
  end type scf_t

contains

  ! ---------------------------------------------------------
  subroutine scf_init(scf, namespace, gr, geo, st, mc, hm, ks, conv_force)
    type(scf_t),          intent(inout) :: scf
    type(namespace_t),    intent(in)    :: namespace
    type(grid_t), target, intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_elec_t),  intent(in)    :: st
    type(multicomm_t),    intent(in)    :: mc
    type(hamiltonian_elec_t),  intent(inout) :: hm
    type(v_ks_t),         intent(in)    :: ks
    FLOAT,   optional,    intent(in)    :: conv_force

    FLOAT :: rmin
    integer :: mixdefault, ierr
    type(type_t) :: mix_type

    PUSH_SUB(scf_init)

    !%Variable MaximumIter
    !%Type integer
    !%Default 200
    !%Section SCF::Convergence
    !%Description
    !% Maximum number of SCF iterations. The code will stop even if convergence
    !% has not been achieved. -1 means unlimited.
    !% 0 means just do LCAO (or read from restart), compute the eigenvalues and energy,
    !% and stop, without updating the wavefunctions or density.
    !%End
    call parse_variable(namespace, 'MaximumIter', 200, scf%max_iter)

    !%Variable MaximumIterBerry
    !%Type integer
    !%Default 10
    !%Section SCF::Convergence
    !%Description
    !% Maximum number of iterations for the Berry potential, within each SCF iteration.
    !% Only applies if a <tt>StaticElectricField</tt> is applied in a periodic direction.
    !% The code will move on to the next SCF iteration even if convergence
    !% has not been achieved. -1 means unlimited.
    !%End
    if(associated(hm%vberry)) then
      call parse_variable(namespace, 'MaximumIterBerry', 10, scf%max_iter_berry)
      if(scf%max_iter_berry < 0) scf%max_iter_berry = huge(scf%max_iter_berry)
    end if
    
    !%Variable ConvEnergy
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Stop the SCF when the magnitude of change in energy during at
    !% one SCF iteration is smaller than this value.
    !%
    !%A zero value (the default) means do not use this criterion.
    !%End
    call parse_variable(namespace, 'ConvEnergy', M_ZERO, scf%conv_energy_diff, unit = units_inp%energy)
    
    !%Variable ConvAbsDens
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the density: 
    !%
    !% <math>\varepsilon = \int {\rm d}^3r \left| \rho^{out}(\bf r) -\rho^{inp}(\bf r) \right|</math>.
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call parse_variable(namespace, 'ConvAbsDens', M_ZERO, scf%conv_abs_dens)

    !%Variable ConvRelDens
    !%Type float
    !%Default 1e-5
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the density: 
    !%
    !% <math>\varepsilon = \frac{1}{N} \mathrm{ConvAbsDens}</math>.
    !% 
    !% <i>N</i> is the total number of electrons in the problem.  A
    !% zero value means do not use this criterion.
    !%
    !% If you reduce this value, you should also reduce
    !% <tt>EigensolverTolerance</tt> to a value of roughly 1/10 of
    !% <tt>ConvRelDens</tt> to avoid convergence problems.
    !%End
    call parse_variable(namespace, 'ConvRelDens', CNST(1e-5), scf%conv_rel_dens)

    !%Variable ConvAbsEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the sum of the eigenvalues:
    !%
    !% <math> \varepsilon = \left| \sum_{j=1}^{N_{occ}} \varepsilon_j^{out} -
    !% \sum_{j=1}^{N_{occ}} \varepsilon_j^{inp} \right| </math>
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call parse_variable(namespace, 'ConvAbsEv', M_ZERO, scf%conv_abs_ev, unit = units_inp%energy)

    !%Variable ConvRelEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the sum of the eigenvalues:
    !%
    !% <math>\varepsilon = \frac{ \left| \sum_{j=1}^{N_{occ}} ( \varepsilon_j^{out} -  \varepsilon_j^{inp} ) \right|}
    !% {\left| \sum_{j=1}^{N_{occ}} \varepsilon_j^{out} \right|} </math>
    !%
    !%A zero value (the default) means do not use this criterion.
    !%End
    call parse_variable(namespace, 'ConvRelEv', M_ZERO, scf%conv_rel_ev, unit = units_inp%energy)

    call messages_obsolete_variable(namespace, 'ConvAbsForce', 'ConvForce')
    call messages_obsolete_variable(namespace, 'ConvRelForce', 'ConvForce')

    !%Variable ConvForce
    !%Type float
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the forces: maximum variation of any
    !% component of the ionic forces in consecutive iterations.  A
    !% zero value means do not use this criterion. The default is
    !% zero, except for geometry optimization, which sets a default of
    !% 1e-8 H/b.
    !%End
    call parse_variable(namespace, 'ConvForce', optional_default(conv_force, M_ZERO), scf%conv_abs_force, unit = units_inp%force)

    scf%check_conv = &
      scf%conv_energy_diff > M_ZERO .or. &
      scf%conv_abs_dens > M_ZERO .or. scf%conv_rel_dens > M_ZERO .or. &
      scf%conv_abs_ev > M_ZERO .or. scf%conv_rel_ev > M_ZERO .or. &
      scf%conv_abs_force > M_ZERO

    if(.not. scf%check_conv .and. scf%max_iter < 0) then
      call messages_write("All convergence criteria are disabled. Octopus is cowardly refusing")
      call messages_new_line()
      call messages_write("to enter an infinite loop.")
      call messages_new_line()
      call messages_new_line()
      call messages_write("Please set one of the following variables to a positive value:")
      call messages_new_line()
      call messages_new_line()
      call messages_write(" | MaximumIter | ConvEnergy | ConvAbsDens | ConvRelDens |")
      call messages_new_line()
      call messages_write(" |  ConvAbsEv  | ConvRelEv  |  ConvForce  |")
      call messages_new_line()
      call messages_fatal()
    end if

    !%Variable ConvEigenError
    !%Type logical
    !%Default false
    !%Section SCF::Convergence
    !%Description
    !% If true, the calculation will not be considered converged unless all states have
    !% individual errors less than <tt>EigensolverTolerance</tt>.
    !%End
    call parse_variable(namespace, 'ConvEigenError', .false., scf%conv_eigen_error)

    if(scf%max_iter < 0) scf%max_iter = huge(scf%max_iter)

    call messages_obsolete_variable(namespace, 'What2Mix', 'MixField')

    !%Variable MixField
    !%Type integer
    !%Section SCF::Mixing
    !%Description
    !% Selects what should be mixed during the SCF cycle.  Note that
    !% currently the exact-exchange part of hybrid functionals is not
    !% mixed at all, which would require wavefunction-mixing, not yet
    !% implemented. This may lead to instabilities in the SCF cycle,
    !% so starting from a converged LDA/GGA calculation is recommended
    !% for hybrid functionals. The default depends on the <tt>TheoryLevel</tt>
    !% and the exchange-correlation potential used.
    !%Option none 0
    !% No mixing is done. This is the default for independent
    !% particles.
    !%Option potential 1
    !% The Kohn-Sham potential is mixed. This is the default for other cases.
    !%Option density 2
    !% Mix the density.
    !%Option states 3
    !% (Experimental) Mix the states. In this case, the mixing is always linear.
    !%End

    mixdefault = OPTION__MIXFIELD__POTENTIAL
    if(hm%theory_level == INDEPENDENT_PARTICLES) mixdefault = OPTION__MIXFIELD__NONE

    call parse_variable(namespace, 'MixField', mixdefault, scf%mix_field)
    if(.not.varinfo_valid_option('MixField', scf%mix_field)) call messages_input_error('MixField')
    call messages_print_var_option(stdout, 'MixField', scf%mix_field, "what to mix during SCF cycles")

    if (scf%mix_field == OPTION__MIXFIELD__POTENTIAL .and. hm%theory_level == INDEPENDENT_PARTICLES) then
      call messages_write('Input: Cannot mix the potential for non-interacting particles.')
      call messages_fatal()
    end if

    if (scf%mix_field == OPTION__MIXFIELD__POTENTIAL .and. hm%pcm%run_pcm) then
      call messages_write('Input: You have selected to mix the potential.', new_line = .true.)
      call messages_write('       This might produce convergence problems for solvated systems.', new_line = .true.)
      call messages_write('       Mix the Density instead.')
      call messages_warning()
    end if

    if(scf%mix_field == OPTION__MIXFIELD__DENSITY &
      .and. bitand(hm%xc_family, XC_FAMILY_OEP + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0) then

      call messages_write('Input: You have selected to mix the density with OEP or MGGA XC functionals.', new_line = .true.)
      call messages_write('       This might produce convergence problems. Mix the potential instead.')
      call messages_warning()
    end if

    if(scf%mix_field == OPTION__MIXFIELD__STATES) then
      call messages_experimental('MixField = states')
    end if
    
    ! Handle mixing now...
    select case(scf%mix_field)
    case(OPTION__MIXFIELD__POTENTIAL)
      scf%mixdim1 = gr%mesh%np
    case(OPTION__MIXFIELD__DENSITY)
      scf%mixdim1 = gr%fine%mesh%np
    case(OPTION__MIXFIELD__STATES)
      ! we do not really need the mixer, except for the value of the mixing coefficient
      scf%mixdim1 = 1
    end select

    mix_type = TYPE_FLOAT

    if(scf%mix_field == OPTION__MIXFIELD__DENSITY) then
      call mix_init(scf%smix, namespace, gr%fine%der, scf%mixdim1, 1, st%d%nspin, func_type_ = mix_type)
    else if(scf%mix_field /= OPTION__MIXFIELD__NONE) then
      call mix_init(scf%smix, namespace, gr%der, scf%mixdim1, 1, st%d%nspin, func_type_ = mix_type)
    end if

    !If we use LDA+U, we also have do mix it
    if(scf%mix_field /= OPTION__MIXFIELD__STATES) then
      call lda_u_mixer_init(hm%lda_u, scf%lda_u_mix, st)
      call lda_u_mixer_init_auxmixer(hm%lda_u, scf%lda_u_mix, scf%smix, st)
    end if
    call mix_get_field(scf%smix, scf%mixfield)

    if(hm%lda_u_level /= DFT_U_NONE .and. hm%lda_u%basisfromstates) then
      call lda_u_loadbasis(hm%lda_u, namespace, st, gr%mesh, mc, ierr)
      if (ierr /= 0) then
        message(1) = "Unable to load LDA+U basis from selected states."
        call messages_fatal(1)
      end if
      call lda_u_periodic_coulomb_integrals(hm%lda_u, namespace, st, gr%der, mc, associated(hm%hm_base%phase))
    end if

    ! now the eigensolver stuff
    call eigensolver_init(scf%eigens, namespace, gr, st, ks%xc)

    if(preconditioner_is_multigrid(scf%eigens%pre)) then
      SAFE_ALLOCATE(gr%mgrid_prec)
      call multigrid_init(gr%mgrid_prec, namespace, geo, gr%cv,gr%mesh, gr%der, gr%stencil, mc, used_for_preconditioner = .true.)
    end if

    !%Variable SCFinLCAO
    !%Type logical
    !%Default no
    !%Section SCF
    !%Description
    !% Performs the SCF cycle with the calculation restricted to the LCAO subspace.
    !% This may be useful for systems with convergence problems (first do a 
    !% calculation within the LCAO subspace, then restart from that point for
    !% an unrestricted calculation).
    !%End
    call parse_variable(namespace, 'SCFinLCAO', .false., scf%lcao_restricted)
    if(scf%lcao_restricted) then
      call messages_experimental('SCFinLCAO')
      message(1) = 'Info: SCF restricted to LCAO subspace.'
      call messages_info(1)

      if(scf%conv_eigen_error) then
        message(1) = "ConvEigenError cannot be used with SCFinLCAO, since error is unknown."
        call messages_fatal(1)
      end if
    end if


    !%Variable SCFCalculateForces
    !%Type logical
    !%Section SCF
    !%Description
    !% This variable controls whether the forces on the ions are
    !% calculated at the end of a self-consistent iteration. The
    !% default is yes, unless the system only has user-defined
    !% species.
    !%End
    call parse_variable(namespace, 'SCFCalculateForces', .not. geo%only_user_def, scf%calc_force)


    !%Variable SCFCalculateStress
    !%Type logical
    !%Section SCF
    !%Description
    !% This variable controls whether the stress on the lattice is
    !% calculated at the end of a self-consistent iteration. The
    !% default is no.
    !%End
    call parse_variable(namespace, 'SCFCalculateStress', .false. , scf%calc_stress)
    
    !%Variable SCFCalculateDipole
    !%Type logical
    !%Section SCF
    !%Description
    !% This variable controls whether the dipole is calculated at the
    !% end of a self-consistent iteration. For finite systems the
    !% default is yes. For periodic systems the default is no, unless
    !% an electric field is being applied in a periodic direction.
    !% The single-point Berry`s phase approximation is used for
    !% periodic directions. Ref:
    !% E Yaschenko, L Fu, L Resca, and R Resta, <i>Phys. Rev. B</i> <b>58</b>, 1222-1229 (1998).
    !%End
    call parse_variable(namespace, 'SCFCalculateDipole', .not. simul_box_is_periodic(gr%sb), scf%calc_dipole)
    if(associated(hm%vberry)) scf%calc_dipole = .true.

    !%Variable SCFCalculatePartialCharges
    !%Type logical
    !%Default no
    !%Section SCF
    !%Description
    !% (Experimental) This variable controls whether partial charges
    !% are calculated at the end of a self-consistent iteration.
    !%End
    call parse_variable(namespace, 'SCFCalculatePartialCharges', .false., scf%calc_partial_charges)
    if(scf%calc_partial_charges) call messages_experimental('SCFCalculatePartialCharges')

    rmin = geometry_min_distance(geo)
    if(geo%natoms == 1) then
      if(simul_box_is_periodic(gr%sb)) then
        rmin = minval(gr%sb%lsize(1:gr%sb%periodic_dim))
      else
        rmin = CNST(100.0)
      end if
    end if

    !%Variable LocalMagneticMomentsSphereRadius
    !%Type float
    !%Section Output
    !%Description
    !% The local magnetic moments are calculated by integrating the
    !% magnetization density in spheres centered around each atom.
    !% This variable controls the radius of the spheres.
    !% The default is half the minimum distance between two atoms
    !% in the input coordinates, or 100 a.u. if there is only one atom (for isolated systems).
    !%End
    call parse_variable(namespace, 'LocalMagneticMomentsSphereRadius', rmin*M_HALF, scf%lmm_r, unit = units_inp%length)
    ! this variable is also used in td/td_write.F90

    scf%forced_finish = .false.

    scf%gr => gr
    
    POP_SUB(scf_init)
  end subroutine scf_init


  ! ---------------------------------------------------------
  subroutine scf_end(scf)
    type(scf_t),  intent(inout) :: scf
    
    PUSH_SUB(scf_end)

    if(preconditioner_is_multigrid(scf%eigens%pre)) then
      call multigrid_end(scf%gr%mgrid_prec)
      SAFE_DEALLOCATE_P(scf%gr%mgrid_prec)
    end if

    call eigensolver_end(scf%eigens)
    if(scf%mix_field /= OPTION__MIXFIELD__NONE) call mix_end(scf%smix)

    nullify(scf%mixfield)

    if(scf%mix_field /= OPTION__MIXFIELD__STATES) call lda_u_mixer_end(scf%lda_u_mix, scf%smix)


    POP_SUB(scf_end)
  end subroutine scf_end


  ! ---------------------------------------------------------
  subroutine scf_mix_clear(scf)
    type(scf_t), intent(inout) :: scf

    PUSH_SUB(scf_mix_clear)

    call mix_clear(scf%smix)

    if(scf%mix_field /= OPTION__MIXFIELD__STATES) call lda_u_mixer_clear(scf%lda_u_mix, scf%smix)

    POP_SUB(scf_mix_clear)
  end subroutine scf_mix_clear


  ! ---------------------------------------------------------
  subroutine scf_run(scf, namespace, mc, gr, geo, st, ks, hm, outp, gs_run, verbosity, iters_done, &
    restart_load, restart_dump)
    type(scf_t),               intent(inout) :: scf !< self consistent cycle
    type(namespace_t),         intent(in)    :: namespace
    type(multicomm_t),         intent(in)    :: mc
    type(grid_t),              intent(inout) :: gr !< grid
    type(geometry_t),          intent(inout) :: geo !< geometry
    type(states_elec_t),       intent(inout) :: st !< States
    type(v_ks_t),              intent(inout) :: ks !< Kohn-Sham
    type(hamiltonian_elec_t),  intent(inout) :: hm !< Hamiltonian
    type(output_t),            intent(in)    :: outp
    logical,         optional, intent(in)    :: gs_run
    integer,         optional, intent(in)    :: verbosity 
    integer,         optional, intent(out)   :: iters_done
    type(restart_t), optional, intent(in)    :: restart_load
    type(restart_t), optional, intent(in)    :: restart_dump

    logical :: finish, gs_run_, berry_conv, forced_finish_tmp
    integer :: iter, is, iatom, nspin, ierr, iberry, idir, verbosity_, ib, iqn
    FLOAT :: evsum_out, evsum_in, forcetmp, dipole(MAX_DIM), dipole_prev(MAX_DIM)
    real(8) :: etime, itime
    character(len=MAX_PATH_LEN) :: dirname
    type(lcao_t) :: lcao    !< Linear combination of atomic orbitals
    type(profile_t), save :: prof
    FLOAT, allocatable :: rhoout(:,:,:), rhoin(:,:,:)
    FLOAT, allocatable :: vhxc_old(:,:)
    FLOAT, allocatable :: forceout(:,:), forcein(:,:), forcediff(:), tmp(:)
    type(batch_t), allocatable :: psioutb(:, :)

    PUSH_SUB(scf_run)

    if(scf%forced_finish) then
      message(1) = "Previous clean stop, not doing SCF and quitting."
      call messages_fatal(1, only_root_writes = .true.)
    end if

    gs_run_ = .true.
    if(present(gs_run)) gs_run_ = gs_run

    verbosity_ = VERB_FULL
    if(present(verbosity)) verbosity_ = verbosity

    if(scf%lcao_restricted) then
      call lcao_init(lcao, namespace, gr, geo, st)
      if(.not. lcao_is_available(lcao)) then
        message(1) = 'LCAO is not available. Cannot do SCF in LCAO.'
        call messages_fatal(1)
      end if
    end if

    nspin = st%d%nspin

    if (present(restart_load)) then
      if (restart_has_flag(restart_load, RESTART_FLAG_RHO)) then
        ! Load density and used it to recalculated the KS potential.
        call states_elec_load_rho(restart_load, st, gr, ierr)
        if (ierr /= 0) then
          message(1) = 'Unable to read density. Density will be calculated from states.'
          call messages_warning(1)
        else
          if(bitand(ks%xc_family, XC_FAMILY_OEP) == 0) then
            call v_ks_calc(ks, namespace, hm, st, geo)
          else
            if (.not. restart_has_flag(restart_load, RESTART_FLAG_VHXC) .and. ks%oep%level /= XC_OEP_FULL) then
              call v_ks_calc(ks, namespace, hm, st, geo)
            end if
          end if
        end if
      end if

      if (restart_has_flag(restart_load, RESTART_FLAG_VHXC)) then
        call hamiltonian_elec_load_vhxc(restart_load, hm, gr%mesh, ierr)
        if (ierr /= 0) then
          message(1) = 'Unable to read Vhxc. Vhxc will be calculated from states.'
          call messages_warning(1)
        else
          call hamiltonian_elec_update(hm, gr%mesh, namespace)
          if(bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
            if (ks%oep%level == XC_OEP_FULL) then
              do is = 1, st%d%nspin
                ks%oep%vxc(1:gr%mesh%np, is) = hm%vhxc(1:gr%mesh%np, is) - hm%vhartree(1:gr%mesh%np)
              end do
              call v_ks_calc(ks, namespace, hm, st, geo)
            end if
          end if
        end if
      end if

      if (restart_has_flag(restart_load, RESTART_FLAG_MIX)) then
        select case (scf%mix_field)
        case (OPTION__MIXFIELD__DENSITY)
          call mix_load(restart_load, scf%smix, gr%fine%mesh, ierr)
        case (OPTION__MIXFIELD__POTENTIAL)
          call mix_load(restart_load, scf%smix, gr%mesh, ierr)
        end select
        if (ierr /= 0) then
          message(1) = "Unable to read mixing information. Mixing will start from scratch."
          call messages_warning(1)
        end if
      end if

      if(hm%lda_u_level /= DFT_U_NONE) then
        call lda_u_load(restart_load, hm%lda_u, st, ierr) 
        if (ierr /= 0) then
          message(1) = "Unable to read LDA+U information. LDA+U data will be calculated from states."
          call messages_warning(1)
        end if
      end if 
    end if

    SAFE_ALLOCATE(rhoout(1:gr%fine%mesh%np, 1:1, 1:nspin))
    SAFE_ALLOCATE(rhoin (1:gr%fine%mesh%np, 1:1, 1:nspin))

    rhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)
    rhoout = M_ZERO

    !We store the Hxc potential for the contribution to the forces
    if(scf%calc_force .or. scf%conv_abs_force > M_ZERO &
        .or. (outp%duringscf .and. bitand(outp%what, OPTION__OUTPUT__FORCES) /= 0)) then
      SAFE_ALLOCATE(vhxc_old(1:gr%mesh%np, 1:nspin))
      vhxc_old(1:gr%mesh%np, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
    end if
    

    select case(scf%mix_field)
    case(OPTION__MIXFIELD__POTENTIAL)
      call mixfield_set_vin(scf%mixfield, hm%vhxc)
    case(OPTION__MIXFIELD__DENSITY)
      call mixfield_set_vin(scf%mixfield, rhoin)

    case(OPTION__MIXFIELD__STATES)

      SAFE_ALLOCATE(psioutb(st%group%block_start:st%group%block_end, st%d%kpt%start:st%d%kpt%end))

      do iqn = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call batch_copy(st%group%psib(ib, iqn), psioutb(ib, iqn))
        end do
      end do
      
    end select

    call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy)
    !If we use LDA+U, we also have do mix it
    if(scf%mix_field /= OPTION__MIXFIELD__STATES) call lda_u_mixer_set_vin(hm%lda_u, scf%lda_u_mix)

    evsum_in = states_elec_eigenvalues_sum(st)

    ! allocate and compute forces only if they are used as convergence criteria
    if (scf%conv_abs_force > M_ZERO) then
      SAFE_ALLOCATE(  forcein(1:geo%natoms, 1:gr%sb%dim))
      SAFE_ALLOCATE( forceout(1:geo%natoms, 1:gr%sb%dim))
      SAFE_ALLOCATE(forcediff(1:gr%sb%dim))
      call forces_calculate(gr, namespace, geo, hm, st, ks)
      do iatom = 1, geo%natoms
        forcein(iatom, 1:gr%sb%dim) = geo%atom(iatom)%f(1:gr%sb%dim)
      end do
    end if

    call create_convergence_file(STATIC_DIR, "convergence")
    
    if ( verbosity_ /= VERB_NO ) then
      if(scf%max_iter > 0) then
        write(message(1),'(a)') 'Info: Starting SCF iteration.'
      else
        write(message(1),'(a)') 'Info: No SCF iterations will be done.'
        ! we cannot tell whether it is converged.
        finish = .false.
      end if
      call messages_info(1)
    end if

    ! SCF cycle
    itime = loct_clock()
    
    do iter = 1, scf%max_iter
      call profiling_in(prof, "SCF_CYCLE")

      ! reset scdm flag
      scdm_is_local = .false.
       
      ! this initialization seems redundant but avoids improper optimization at -O3 by PGI 7 on chum,
      ! which would cause a failure of testsuite/linear_response/04-vib_modes.03-vib_modes_fd.inp
      scf%eigens%converged = 0

      scf%energy_diff = hm%energy%total

      !Used for the contribution to the forces
      if(scf%calc_force .or. scf%conv_abs_force > M_ZERO .or. &
          (outp%duringscf .and. bitand(outp%what, OPTION__OUTPUT__FORCES) /= 0)) & 
        vhxc_old(1:gr%mesh%np, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
      
      if(scf%lcao_restricted) then
        call lcao_init_orbitals(lcao, st, gr, geo)
        call lcao_wf(lcao, st, gr, geo, hm, namespace)
      else
        if(associated(hm%vberry)) then
          ks%frozen_hxc = .true.
          do iberry = 1, scf%max_iter_berry
            scf%eigens%converged = 0
            call eigensolver_run(scf%eigens, gr, st, hm, iter)

            call v_ks_calc(ks, namespace, hm, st, geo, calc_current=outp%duringscf)

            dipole_prev = dipole
            call calc_dipole(dipole)
            write(message(1),'(a,9f12.6)') 'Dipole = ', dipole(1:gr%sb%dim)
            call messages_info(1)

            berry_conv = .true.
            do idir = 1, gr%sb%periodic_dim
              berry_conv = berry_conv .and. &
                (abs((dipole(idir) - dipole_prev(idir)) / dipole_prev(idir)) < CNST(1e-5) &
                .or. abs(dipole(idir) - dipole_prev(idir)) < CNST(1e-5))
            end do
            if(berry_conv) exit
          end do
          ks%frozen_hxc = .false.
        else
          scf%eigens%converged = 0
          call eigensolver_run(scf%eigens, gr, st, hm, iter)
        end if
      end if

      ! occupations
      call states_elec_fermi(st, gr%mesh)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )

      ! compute output density, potential (if needed) and eigenvalues sum
      call density_calc(st, gr, st%rho)

      rhoout(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)

      select case(scf%mix_field)
      case(OPTION__MIXFIELD__POTENTIAL)
        call v_ks_calc(ks, namespace, hm, st, geo, calc_current=outp%duringscf)
        call mixfield_set_vout(scf%mixfield, hm%vhxc)
      case (OPTION__MIXFIELD__DENSITY)
        call mixfield_set_vout(scf%mixfield, rhoout)
      case(OPTION__MIXFIELD__STATES)

        do iqn = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call batch_copy_data(gr%mesh%np, st%group%psib(ib, iqn), psioutb(ib, iqn))
          end do
        end do
      end select
      
      if(scf%mix_field /= OPTION__MIXFIELD__STATES) call lda_u_mixer_set_vout(hm%lda_u, scf%lda_u_mix)
 
      evsum_out = states_elec_eigenvalues_sum(st)

      ! recalculate total energy
      call energy_calc_total(hm, gr, st, iunit = 0)

      ! compute convergence criteria
      scf%energy_diff = hm%energy%total - scf%energy_diff
      scf%abs_dens = M_ZERO
      SAFE_ALLOCATE(tmp(1:gr%fine%mesh%np))
      do is = 1, nspin
        tmp = abs(rhoin(1:gr%fine%mesh%np, 1, is) - rhoout(1:gr%fine%mesh%np, 1, is))
        scf%abs_dens = scf%abs_dens + dmf_integrate(gr%fine%mesh, tmp)
      end do
      SAFE_DEALLOCATE_A(tmp)

      ! compute forces only if they are used as convergence criterion
      if (scf%conv_abs_force > M_ZERO) then
        call forces_calculate(gr, namespace, geo, hm, st, ks, vhxc_old=vhxc_old)
        scf%abs_force = M_ZERO
        do iatom = 1, geo%natoms
          forceout(iatom,1:gr%sb%dim) = geo%atom(iatom)%f(1:gr%sb%dim)
          forcediff(1:gr%sb%dim) = abs( forceout(iatom,1:gr%sb%dim) - forcein(iatom,1:gr%sb%dim) )
          forcetmp = maxval( forcediff )
          if ( forcetmp > scf%abs_force ) then
            scf%abs_force = forcetmp
          end if
        end do
      else
        if(outp%duringscf .and. bitand(outp%what, OPTION__OUTPUT__FORCES) /= 0 &
           .and. outp%output_interval /= 0 &
           .and. gs_run_ .and. mod(iter, outp%output_interval) == 0)  &
          call forces_calculate(gr, namespace, geo, hm, st, ks, vhxc_old=vhxc_old)
      end if

      if(abs(st%qtot) <= M_EPSILON) then
        scf%rel_dens = M_HUGE
      else
        scf%rel_dens = scf%abs_dens / st%qtot
      end if

      scf%abs_ev = abs(evsum_out - evsum_in)
      if(abs(evsum_out) <= M_EPSILON) then
        scf%rel_ev = M_HUGE
      else
        scf%rel_ev = scf%abs_ev / abs(evsum_out)
      end if

      scf%eigens%current_rel_dens_error = scf%rel_dens

      ! are we finished?
      finish = scf%check_conv .and. &
        (scf%conv_abs_dens  <= M_ZERO .or. scf%abs_dens  <= scf%conv_abs_dens)  .and. &
        (scf%conv_rel_dens  <= M_ZERO .or. scf%rel_dens  <= scf%conv_rel_dens)  .and. &
        (scf%conv_abs_force <= M_ZERO .or. scf%abs_force <= scf%conv_abs_force) .and. &
        (scf%conv_abs_ev    <= M_ZERO .or. scf%abs_ev    <= scf%conv_abs_ev)    .and. &
        (scf%conv_rel_ev    <= M_ZERO .or. scf%rel_ev    <= scf%conv_rel_ev)    .and. &
        (scf%conv_energy_diff <= M_ZERO .or. abs(scf%energy_diff) <= scf%conv_energy_diff) .and. &
        (.not. scf%conv_eigen_error .or. all(scf%eigens%converged == st%nst))

      etime = loct_clock() - itime
      itime = etime + itime
      call scf_write_iter()

      ! mixing
      select case (scf%mix_field)
      case (OPTION__MIXFIELD__DENSITY)
        ! mix input and output densities and compute new potential
        call mixing(scf%smix)
        call mixfield_get_vnew(scf%mixfield, st%rho)
        ! for spinors, having components 3 or 4 be negative is not unphysical
        if(minval(st%rho(1:gr%fine%mesh%np, 1:st%d%spin_channels)) < -CNST(1e-6)) then
          write(message(1),*) 'Negative density after mixing. Minimum value = ', &
            minval(st%rho(1:gr%fine%mesh%np, 1:st%d%spin_channels))
          call messages_warning(1)
        end if
        call lda_u_mixer_get_vnew(hm%lda_u, scf%lda_u_mix, st)
        call v_ks_calc(ks, namespace, hm, st, geo, calc_current=outp%duringscf)
      case (OPTION__MIXFIELD__POTENTIAL)
        ! mix input and output potentials
        call mixing(scf%smix)
        call mixfield_get_vnew(scf%mixfield, hm%vhxc)
        call lda_u_mixer_get_vnew(hm%lda_u, scf%lda_u_mix, st)
        call hamiltonian_elec_update(hm, gr%mesh, namespace)
        
      case(OPTION__MIXFIELD__STATES)

        do iqn = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call batch_scal(gr%mesh%np, CNST(1.0) - mix_coefficient(scf%smix), st%group%psib(ib, iqn))
            call batch_axpy(gr%mesh%np, mix_coefficient(scf%smix), psioutb(ib, iqn), st%group%psib(ib, iqn))
          end do
        end do

        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, namespace, hm, st, geo, calc_current=outp%duringscf)
        
      case(OPTION__MIXFIELD__NONE)
        call v_ks_calc(ks, namespace, hm, st, geo, calc_current=outp%duringscf)
      end select


      ! Are we asked to stop? (Whenever Fortran is ready for signals, this should go away)
      scf%forced_finish = clean_stop(mc%master_comm) .or. walltimer_alarm()

#ifdef HAVE_MPI
      call MPI_Allreduce(scf%forced_finish, forced_finish_tmp, 1, MPI_LOGICAL, MPI_LOR, mc%master_comm, mpi_err)
      scf%forced_finish = forced_finish_tmp
#endif      

      if (finish .and. st%modelmbparticles%nparticle > 0) then
        call modelmb_sym_all_states (gr, st)
      end if

      if (gs_run_ .and. present(restart_dump)) then 
        ! save restart information
         
        if ( (finish .or. (modulo(iter, outp%restart_write_interval) == 0) &
          .or. iter == scf%max_iter .or. scf%forced_finish) ) then

          call states_elec_dump(restart_dump, st, gr, ierr, iter=iter) 
          if (ierr /= 0) then
            message(1) = 'Unable to write states wavefunctions.'
            call messages_warning(1)
          end if

          call states_elec_dump_rho(restart_dump, st, gr, ierr, iter=iter)
          if (ierr /= 0) then
            message(1) = 'Unable to write density.'
            call messages_warning(1)
          end if

          if(hm%lda_u_level /= DFT_U_NONE) then
            call lda_u_dump(restart_dump, hm%lda_u, st, ierr, iter=iter)
            if (ierr /= 0) then
              message(1) = 'Unable to write LDA+U information.'
              call messages_warning(1)
            end if
          end if

          select case (scf%mix_field)
          case (OPTION__MIXFIELD__DENSITY)
            call mix_dump(restart_dump, scf%smix, gr%fine%mesh, ierr)
            if (ierr /= 0) then
              message(1) = 'Unable to write mixing information.'
              call messages_warning(1)
            end if
          case (OPTION__MIXFIELD__POTENTIAL)
            call hamiltonian_elec_dump_vhxc(restart_dump, hm, gr%mesh, ierr)
            if (ierr /= 0) then
              message(1) = 'Unable to write Vhxc.'
              call messages_warning(1)
            end if

            call mix_dump(restart_dump, scf%smix, gr%mesh, ierr)
            if (ierr /= 0) then
              message(1) = 'Unable to write mixing information.'
              call messages_warning(1)
            end if
          end select
        end if
      end if

      call write_convergence_file(STATIC_DIR, "convergence")
      
      if(finish) then
        if(present(iters_done)) iters_done = iter
        if(verbosity_ >= VERB_COMPACT) then
          write(message(1), '(a, i4, a)') 'Info: SCF converged in ', iter, ' iterations'
          write(message(2), '(a)')        '' 
          call messages_info(2)
        end if
        call profiling_out(prof)
        exit
      end if

      if((outp%what+outp%what_lda_u+outp%whatBZ)/=0 .and. outp%duringscf .and. outp%output_interval /= 0 &
        .and. gs_run_ .and. mod(iter, outp%output_interval) == 0) then
        write(dirname,'(a,a,i4.4)') trim(outp%iter_dir),"scf.",iter
        call output_all(outp, namespace, gr, geo, st, hm, ks, dirname)
      end if

      ! save information for the next iteration
      rhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)

      select case(scf%mix_field)
        case(OPTION__MIXFIELD__POTENTIAL)
          call mixfield_set_vin(scf%mixfield, hm%vhxc(1:gr%mesh%np, 1:nspin))
        case (OPTION__MIXFIELD__DENSITY)
          call mixfield_set_vin(scf%mixfield, rhoin)
      end select
      if(scf%mix_field /= OPTION__MIXFIELD__STATES) call lda_u_mixer_set_vin(hm%lda_u, scf%lda_u_mix)

      evsum_in = evsum_out
      if (scf%conv_abs_force > M_ZERO) then
        forcein(1:geo%natoms, 1:gr%sb%dim) = forceout(1:geo%natoms, 1:gr%sb%dim)
      end if


      if(scf%forced_finish) then
        call profiling_out(prof)
        exit
      end if

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly(namespace)

      call profiling_out(prof)
    end do !iter

    if(scf%lcao_restricted) call lcao_end(lcao)

    if((scf%max_iter > 0 .and. scf%mix_field == OPTION__MIXFIELD__POTENTIAL) .or. bitand(outp%what, OPTION__OUTPUT__CURRENT) /= 0) then
      call v_ks_calc(ks, namespace, hm, st, geo)
    end if

    select case(scf%mix_field)
    case(OPTION__MIXFIELD__STATES)

      do iqn = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call batch_end(psioutb(ib, iqn))
        end do
      end do
      
      SAFE_DEALLOCATE_A(psioutb)
    end select

    SAFE_DEALLOCATE_A(rhoout)
    SAFE_DEALLOCATE_A(rhoin)

    if(scf%max_iter > 0 .and. any(scf%eigens%converged < st%nst) .and. .not. scf%lcao_restricted) then
      write(message(1),'(a)') 'Some of the states are not fully converged!'
      call messages_warning(1)
    end if

    if(.not.finish) then
      write(message(1), '(a,i4,a)') 'SCF *not* converged after ', iter - 1, ' iterations.'
      call messages_warning(1)
    end if

    ! calculate forces
    if(scf%calc_force) then
      call forces_calculate(gr, namespace, geo, hm, st, ks, vhxc_old=vhxc_old)
    end if

    ! calculate stress
    if(scf%calc_stress) call stress_calculate(gr, hm, st, geo, ks) 
    
    if(scf%max_iter == 0) then
      call energy_calc_eigenvalues(hm, gr%der, st)
      call states_elec_fermi(st, gr%mesh)
      call states_elec_write_eigenvalues(stdout, st%nst, st, gr%sb)
    end if

    if(gs_run_) then 
      ! output final information
      call scf_write_static(STATIC_DIR, "info")
      call output_all(outp, namespace, gr, geo, st, hm, ks, STATIC_DIR)
    end if

    if(simul_box_is_periodic(gr%sb) .and. st%d%nik > st%d%nspin) then
      if(bitand(gr%sb%kpoints%method, KPOINTS_PATH) /= 0)  then
        call states_elec_write_bandstructure(STATIC_DIR, namespace, st%nst, st, gr%sb, geo, gr%mesh, &
          hm%hm_base%phase, vec_pot = hm%hm_base%uniform_vector_potential, &
          vec_pot_var = hm%hm_base%vector_potential)
      end if
    end if

    if( ks%vdw_correction == OPTION__VDWCORRECTION__VDW_TS) then
      call vdw_ts_write_c6ab(ks%vdw_ts, geo, STATIC_DIR, 'c6ab_eff', namespace)
    end if

    SAFE_DEALLOCATE_A(vhxc_old)

    POP_SUB(scf_run)

  contains


    ! ---------------------------------------------------------
    subroutine scf_write_iter()
      character(len=50) :: str
      FLOAT :: mem
#ifdef HAVE_MPI
      FLOAT :: mem_tmp
#endif

      PUSH_SUB(scf_run.scf_write_iter)

      if ( verbosity_ == VERB_FULL ) then

        write(str, '(a,i5)') 'SCF CYCLE ITER #' ,iter
        call messages_print_stress(stdout, trim(str))
        write(message(1),'(a,es15.8,2(a,es9.2))') ' etot  = ', units_from_atomic(units_out%energy, hm%energy%total), &
          ' abs_ev   = ', units_from_atomic(units_out%energy, scf%abs_ev), ' rel_ev   = ', scf%rel_ev
        write(message(2),'(a,es15.2,2(a,es9.2))') &
          ' ediff = ', scf%energy_diff, ' abs_dens = ', scf%abs_dens, ' rel_dens = ', scf%rel_dens
        ! write info about forces only if they are used as convergence criteria
        if (scf%conv_abs_force > M_ZERO) then
          write(message(3),'(23x,a,es9.2)') &
            ' force    = ', units_from_atomic(units_out%force, scf%abs_force)
          call messages_info(3)
        else
          call messages_info(2)
        end if

        if(.not.scf%lcao_restricted) then
          write(message(1),'(a,i6)') 'Matrix vector products: ', scf%eigens%matvec
          write(message(2),'(a,i6)') 'Converged eigenvectors: ', sum(scf%eigens%converged(1:st%d%nik))
          call messages_info(2)
          call states_elec_write_eigenvalues(stdout, st%nst, st, gr%sb, scf%eigens%diff, compact = .true.)
        else
          call states_elec_write_eigenvalues(stdout, st%nst, st, gr%sb, compact = .true.)
        end if

        if(associated(hm%vberry)) then
          call calc_dipole(dipole)
          call write_dipole(stdout, dipole)
        end if

        if(st%d%ispin > UNPOLARIZED) then
          call write_magnetic_moments(stdout, gr%mesh, st, geo, scf%lmm_r)
        end if

        if(hm%lda_u_level == DFT_U_ACBN0) then
          call lda_u_write_U(hm%lda_u, stdout)
          call lda_u_write_V(hm%lda_u, stdout)
        end if

        write(message(1),'(a)') ''
        write(message(2),'(a,i5,a,f14.2)') 'Elapsed time for SCF step ', iter,':', etime
        call messages_info(2)

        if(conf%report_memory) then
          mem = loct_get_memory_usage()/(CNST(1024.0)**2)
#ifdef HAVE_MPI
          call MPI_Allreduce(mem, mem_tmp, 1, MPI_FLOAT, MPI_SUM, mpi_world%comm, mpi_err)
          mem = mem_tmp
#endif
          write(message(1),'(a,f14.2)') 'Memory usage [Mbytes]     :', mem
          call messages_info(1)
        end if

        call messages_print_stress(stdout)

      end if

      if ( verbosity_ == VERB_COMPACT ) then
        ! write info about forces only if they are used as convergence criteria
        if (scf%conv_abs_force > M_ZERO) then
          write(message(1),'(a,i4,a,es15.8, 2(a,es9.2), a, f7.1, a)') &
            'iter ', iter, &
            ' : etot ', units_from_atomic(units_out%energy, hm%energy%total), &
            ' : abs_dens', scf%abs_dens, &
            ' : force ', units_from_atomic(units_out%force, scf%abs_force), &
            ' : etime ', etime, 's'
        else
          write(message(1),'(a,i4,a,es15.8, a,es9.2, a, f7.1, a)') &
            'iter ', iter, &
            ' : etot ', units_from_atomic(units_out%energy, hm%energy%total), &
            ' : abs_dens', scf%abs_dens, &
            ' : etime ', etime, 's'
        end if
        call messages_info(1)
      end if

      POP_SUB(scf_run.scf_write_iter)
    end subroutine scf_write_iter


    ! ---------------------------------------------------------
    subroutine scf_write_static(dir, fname)
      character(len=*), intent(in) :: dir, fname

      type(partial_charges_t) :: partial_charges
      integer :: iunit, iatom
      FLOAT, allocatable :: hirshfeld_charges(:)

      PUSH_SUB(scf_run.scf_write_static)

      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir, namespace)
        iunit = io_open(trim(dir) // "/" // trim(fname), namespace, action='write')

        call grid_write_info(gr, geo, iunit)
 
        call symmetries_write_info(gr%mesh%sb%symm, gr%sb%dim, gr%sb%periodic_dim, iunit)

        if(simul_box_is_periodic(gr%sb)) then
          call kpoints_write_info(gr%mesh%sb%kpoints, iunit)
          write(iunit,'(1x)')
        end if

        call v_ks_write_info(ks, iunit)

        ! scf information
        if(finish) then
          write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
        else
          write(iunit, '(a)') 'SCF *not* converged!'
        end if
        write(iunit, '(1x)')

        if(any(scf%eigens%converged < st%nst) .and. .not. scf%lcao_restricted) then
          write(iunit,'(a)') 'Some of the states are not fully converged!'
        end if

        call states_elec_write_eigenvalues(iunit, st%nst, st, gr%sb)
        write(iunit, '(1x)')

        if(simul_box_is_periodic(gr%sb)) then
          call states_elec_write_gaps(iunit, st, gr%sb)
          write(iunit, '(1x)')
        end if

        write(iunit, '(3a)') 'Energy [', trim(units_abbrev(units_out%energy)), ']:'
      else
        iunit = 0
      end if

      call energy_calc_total(hm, gr, st, iunit, full = .true.)

      if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      if(st%d%ispin > UNPOLARIZED) then
        call write_magnetic_moments(iunit, gr%mesh, st, geo, scf%lmm_r)
        if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      end if

      if(hm%lda_u_level == DFT_U_ACBN0) then
          call lda_u_write_U(hm%lda_u, iunit)
          call lda_u_write_V(hm%lda_u, iunit)
          if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
        end if 

      if(scf%calc_dipole) then
        call calc_dipole(dipole)
        call write_dipole(iunit, dipole)
      end if

      if(mpi_grp_is_root(mpi_world)) then
        if(scf%max_iter > 0) then
          write(iunit, '(a)') 'Convergence:'
          write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'abs_dens = ', scf%abs_dens, &
            ' (', scf%conv_abs_dens, ')'
          write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_dens = ', scf%rel_dens, &
            ' (', scf%conv_rel_dens, ')'
          write(iunit, '(6x, a, es15.8,a,es15.8,4a)') 'abs_ev = ', scf%abs_ev, &
            ' (', units_from_atomic(units_out%energy, scf%conv_abs_ev), ')', &
            ' [',  trim(units_abbrev(units_out%energy)), ']'
          write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_ev = ', scf%rel_ev, &
            ' (', scf%conv_rel_ev, ')'
          write(iunit,'(1x)')
        end if
        ! otherwise, these values are uninitialized, and unknown.

        if(scf%calc_force) call forces_write_info(iunit, geo, gr%sb, dir, namespace)

        if(scf%calc_stress) then
           write(iunit,'(a)') "Stress tensor [H/b^3]"
           write(iunit,'(a9,2x,3a18)')"T_{ij}","x","y","z"
           write(iunit,'(a9,2x,3es18.6)')"x", gr%sb%stress_tensor(1,1:3)
           write(iunit,'(a9,2x,3es18.6)')"y", gr%sb%stress_tensor(2,1:3)
           write(iunit,'(a9,2x,3es18.6)')"z", gr%sb%stress_tensor(3,1:3)
        end if
        
      end if

      if(scf%calc_partial_charges) then
        SAFE_ALLOCATE(hirshfeld_charges(1:geo%natoms))

        call partial_charges_init(partial_charges)
        call partial_charges_calculate(partial_charges, namespace, gr%fine%mesh, st, geo, hirshfeld_charges = hirshfeld_charges)
        call partial_charges_end(partial_charges)

        if(mpi_grp_is_root(mpi_world)) then

          write(iunit,'(a)') 'Partial ionic charges'
          write(iunit,'(a)') ' Ion                     Hirshfeld'

          do iatom = 1, geo%natoms
            write(iunit,'(i4,a10,f16.3)') iatom, trim(species_label(geo%atom(iatom)%species)), hirshfeld_charges(iatom)

          end do

        end if

        SAFE_DEALLOCATE_A(hirshfeld_charges)

      end if

      if(mpi_grp_is_root(mpi_world)) then
        call io_close(iunit)
      end if

      POP_SUB(scf_run.scf_write_static)
    end subroutine scf_write_static


    ! ---------------------------------------------------------
    subroutine calc_dipole(dipole)
      FLOAT, intent(out) :: dipole(:)

      integer :: ispin, idir
      FLOAT :: e_dip(MAX_DIM + 1, st%d%nspin), n_dip(MAX_DIM), nquantumpol

      PUSH_SUB(scf_run.calc_dipole)

      dipole(1:MAX_DIM) = M_ZERO

      do ispin = 1, st%d%nspin
        call dmf_multipoles(gr%fine%mesh, st%rho(:, ispin), 1, e_dip(:, ispin))
      end do

      call geometry_dipole(geo, n_dip)

      do idir = 1, gr%sb%dim
        ! in periodic directions use single-point Berry`s phase calculation
        if(idir  <=  gr%sb%periodic_dim) then
          dipole(idir) = -n_dip(idir) - berry_dipole(st, gr%mesh, idir)

          ! use quantum of polarization to reduce to smallest possible magnitude
          nquantumpol = NINT(dipole(idir)/(CNST(2.0)*gr%sb%lsize(idir)))
          dipole(idir) = dipole(idir) - nquantumpol * (CNST(2.0) * gr%sb%lsize(idir))

          ! in aperiodic directions use normal dipole formula
        else
          e_dip(idir + 1, 1) = sum(e_dip(idir + 1, :))
          dipole(idir) = -n_dip(idir) - e_dip(idir + 1, 1)
        end if
      end do

      POP_SUB(scf_run.calc_dipole)
    end subroutine calc_dipole


    ! ---------------------------------------------------------
    subroutine write_dipole(iunit, dipole)
      integer, intent(in) :: iunit
      FLOAT,   intent(in) :: dipole(:)

      PUSH_SUB(scf_run.write_dipole)

      if(mpi_grp_is_root(mpi_world)) then
        call output_dipole(iunit, dipole, gr%mesh%sb%dim)

        if (simul_box_is_periodic(gr%sb)) then
          write(iunit, '(a)') "Defined only up to quantum of polarization (e * lattice vector)."
          write(iunit, '(a)') "Single-point Berry's phase method only accurate for large supercells."

          if (gr%sb%kpoints%full%npoints > 1) then
            write(iunit, '(a)') &
              "WARNING: Single-point Berry's phase method for dipole should not be used when there is more than one k-point."
            write(iunit, '(a)') &
              "Instead, finite differences on k-points (not yet implemented) are needed."
          end if

          if(.not. smear_is_semiconducting(st%smear)) then
            write(iunit, '(a)') "Single-point Berry's phase dipole calculation not correct without integer occupations."
          end if
        end if

        write(iunit, *)
      end if

      POP_SUB(scf_run.write_dipole)
    end subroutine write_dipole

    ! -----------------------------------------------------
    
    subroutine create_convergence_file(dir, fname)
      character(len=*), intent(in) :: dir
      character(len=*), intent(in) :: fname

      integer :: iunit
      character(len=12) :: label
      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir, namespace)
        iunit = io_open(trim(dir) // "/" // trim(fname), namespace, action='write')
        write(iunit, '(a)', advance = 'no') '#iter energy           '
        label = 'energy_diff'
        write(iunit, '(1x,a)', advance = 'no') label
        label = 'abs_dens'
        write(iunit, '(1x,a)', advance = 'no') label
        label = 'rel_dens'
        write(iunit, '(1x,a)', advance = 'no') label
        label = 'abs_ev'
        write(iunit, '(1x,a)', advance = 'no') label
        label = 'rel_ev'
        write(iunit, '(1x,a)', advance = 'no') label
        if (scf%conv_abs_force > M_ZERO) then
          label = 'force_diff'
          write(iunit, '(1x,a)', advance = 'no') label
        end if
        if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0 .and. ks%theory_level /= HARTREE_FOCK) then
          if (ks%oep%level == XC_OEP_FULL) then
            label = 'OEP norm2ss'
            write(iunit, '(1x,a)', advance = 'no') label
          end if
        end if
        write(iunit,'(a)') ''
        call io_close(iunit)
      end if
      
    end subroutine create_convergence_file

    ! -----------------------------------------------------
    
    subroutine write_convergence_file(dir, fname)
      character(len=*), intent(in) :: dir
      character(len=*), intent(in) :: fname

      integer :: iunit
      
      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir, namespace)
        iunit = io_open(trim(dir) // "/" // trim(fname), namespace, action='write', position='append')
        write(iunit, '(i5,es18.8)', advance = 'no') iter, units_from_atomic(units_out%energy, hm%energy%total)
        write(iunit, '(es13.5)', advance = 'no') units_from_atomic(units_out%energy, scf%energy_diff)
        write(iunit, '(es13.5)', advance = 'no') scf%abs_dens
        write(iunit, '(es13.5)', advance = 'no') scf%rel_dens
        write(iunit, '(es13.5)', advance = 'no') units_from_atomic(units_out%energy, scf%abs_ev)
        write(iunit, '(es13.5)', advance = 'no') units_from_atomic(units_out%energy, scf%rel_ev)
        if (scf%conv_abs_force > M_ZERO) then
          write(iunit, '(es13.5)', advance = 'no') units_from_atomic(units_out%force, scf%abs_force)
        end if
        if (bitand(ks%xc_family, XC_FAMILY_OEP) /= 0 .and. ks%theory_level /= HARTREE_FOCK) then
          if (ks%oep%level == XC_OEP_FULL) &
            write(iunit, '(es13.5)', advance = 'no') ks%oep%norm2ss
        end if
        write(iunit,'(a)') ''
        call io_close(iunit)
      end if
      
    end subroutine write_convergence_file
    
  end subroutine scf_run

end module scf_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
