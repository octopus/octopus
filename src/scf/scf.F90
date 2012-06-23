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

module scf_m
  use berry_m
  use datasets_m
  use density_m
  use eigensolver_m
  use energy_calc_m
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kpoints_m
  use lcao_m
  use loct_m
  use magnetic_m
  use mesh_m
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use mix_m
  use mpi_m
  use mpi_lib_m
  use multigrid_m
  use ob_lippmann_schwinger_m
  use parser_m
  use preconditioners_m
  use profiling_m
  use restart_m
  use simul_box_m
  use smear_m
  use solids_m
  use species_m
  use states_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use types_m
  use unit_m
  use unit_system_m
  use utils_m
  use v_ks_m
  use varinfo_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::             &
    scf_t,              &
    scf_init,           &
    scf_mix_clear,      &
    scf_run,            &
    scf_end

  integer, parameter :: &
    MIXNONE = 0,        &
    MIXPOT  = 1,        &
    MIXDENS = 2

  integer, public, parameter :: &
    VERB_NO      = 0,   &
    VERB_COMPACT = 1,   &
    VERB_FULL    = 3
  
  !> some variables used for the SCF cycle
  type scf_t
    integer :: max_iter   !< maximum number of SCF iterations
    integer :: max_iter_berry  !< max number of electronic iterations before updating density, for Berry potential

    FLOAT :: lmm_r

    ! several convergence criteria
    FLOAT :: conv_abs_dens, conv_rel_dens, conv_abs_ev, conv_rel_ev, conv_abs_force
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev, abs_force

    integer :: mix_field
    logical :: lcao_restricted
    logical :: calc_force
    logical :: calc_dipole
    type(mix_t) :: smix
    type(eigensolver_t) :: eigens
    integer :: mixdim1
    integer :: mixdim2
  end type scf_t

contains

  ! ---------------------------------------------------------
  subroutine scf_init(scf, gr, geo, st, hm, conv_force)
    type(scf_t),         intent(inout) :: scf
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT, optional,     intent(in)    :: conv_force

    FLOAT :: rmin
    integer :: mixdefault

    PUSH_SUB(scf_init)

    !%Variable MaximumIter
    !%Type integer
    !%Default 200
    !%Section SCF::Convergence
    !%Description
    !% Maximum number of SCF iterations. The code will stop even if convergence
    !% has not been achieved. -1 means unlimited.
    !%End
    call parse_integer  (datasets_check('MaximumIter'), 200, scf%max_iter)

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
      call parse_integer  (datasets_check('MaximumIterBerry'), 10, scf%max_iter_berry)
      if(scf%max_iter_berry < 0) scf%max_iter_berry = huge(scf%max_iter_berry)
    end if

    !%Variable ConvAbsDens
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the density: 
    !%
    !% <math>\epsilon = \int {\rm d}^3r \vert \rho^{out}(\bf r) -\rho^{inp}(\bf r) \vert</math>.
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call parse_float(datasets_check('ConvAbsDens'), M_ZERO, scf%conv_abs_dens)

    !%Variable ConvRelDens
    !%Type float
    !%Default 1e-5
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the density: 
    !%
    !% <math>\epsilon = {1\over N} ConvAbsDens</math>.
    !% 
    !% <i>N</i> is the total number of electrons in the problem.  A
    !% zero value means do not use this criterion.
    !%End
    call parse_float(datasets_check('ConvRelDens'), CNST(1e-5), scf%conv_rel_dens)

    !%Variable ConvAbsEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the sum of the eigenvalues:
    !%
    !% <math> \epsilon = \vert \sum_{j=1}^{N_{occ}} \epsilon_j^{out} - \sum_{j=1}^{N_{occ}} \epsilon_j^{inp} \vert </math>
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call parse_float(datasets_check('ConvAbsEv'), M_ZERO, scf%conv_abs_ev)
    scf%conv_abs_ev = units_to_atomic(units_inp%energy, scf%conv_abs_ev)

    !%Variable ConvRelEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the sum of the eigenvalues:
    !%
    !% <math>\epsilon = \vert \sum_{j=1}^{N_{occ}} ( \epsilon_j^{out} -  \epsilon_j^{inp} ) \vert
    !% \over \vert \sum_{j=1}^{N_{occ}} \epsilon_j^{out} \vert </math>
    !%
    !%A zero value (the default) means do not use this criterion.
    !%End
    call parse_float(datasets_check('ConvRelEv'), M_ZERO, scf%conv_rel_ev)

    call messages_obsolete_variable("ConvAbsForce", "ConvForce")
    call messages_obsolete_variable("ConvRelForce", "ConvForce")

    !%Variable ConvForce
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the forces: maximum variation of any
    !% component of the ionic forces in consecutive iterations.  A
    !% zero value means do not use this criterion. The default is
    !% zero, except for geometry optimization, which sets a default of
    !% 1e-8.
    !%End
    call parse_float(datasets_check('ConvForce'), optional_default(conv_force, M_ZERO), scf%conv_abs_force)
    scf%conv_abs_force = units_to_atomic(units_inp%force, scf%conv_abs_force)

    if(scf%max_iter < 0 .and. &
      scf%conv_abs_dens <= M_ZERO .and. scf%conv_rel_dens <= M_ZERO .and. &
      scf%conv_abs_ev <= M_ZERO .and. scf%conv_rel_ev <= M_ZERO .and. &
      scf%conv_abs_force <= M_ZERO) then
      message(1) = "Input: Not all convergence criteria can be <= 0"
      message(2) = "Please set one of the following:"
      message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEv | ConvRelEv | ConvForce "
      call messages_fatal(3)
    end if

    if(scf%max_iter < 0) scf%max_iter = huge(scf%max_iter)

    call messages_obsolete_variable('What2Mix', 'MixField')

    !%Variable MixField
    !%Type integer
    !%Default density
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
    !% The Kohn-Sham potential is mixed. This is the default for OEP
    !% or MGGA calculations, or if <tt>StaticElectricField</tt> is applied in
    !% a periodic direction.
    !%Option density 2
    !% Mix the density. This is the default for other cases, including
    !% LDA/GGA calculations.
    !%End

    mixdefault = MIXDENS
    if(hm%theory_level==INDEPENDENT_PARTICLES) mixdefault = MIXNONE
    if(iand(hm%xc_family, XC_FAMILY_OEP + XC_FAMILY_MGGA) /= 0) mixdefault = MIXPOT
    if(associated(hm%vberry)) mixdefault = MIXPOT

    call parse_integer(datasets_check('MixField'), mixdefault, scf%mix_field)
    if(.not.varinfo_valid_option('MixField', scf%mix_field)) call input_error('MixField')
    call messages_print_var_option(stdout, 'MixField', scf%mix_field, "what to mix during SCF cycles")

    if (scf%mix_field == MIXPOT.and.hm%theory_level==INDEPENDENT_PARTICLES) then
      message(1) = "Input: Cannot mix the potential for non-interacting particles."
      call messages_fatal(1)
    end if

    if(scf%mix_field == MIXDENS .and. iand(hm%xc_family, XC_FAMILY_OEP + XC_FAMILY_MGGA) /= 0) then
      message(1) = "Input: You have selected to mix the density with OEP or MGGA XC functionals."
      message(2) = "       This might produce convergence problems. Mix the potential instead."
      call messages_warning(2)
    end if

    ! Handle mixing now...
    select case(scf%mix_field)
    case(MIXPOT)
      scf%mixdim1 = gr%mesh%np
    case(MIXDENS)
      scf%mixdim1 = gr%fine%mesh%np
    end select

    scf%mixdim2 = 1
    if(hm%d%cdft) scf%mixdim2 = 1 + gr%mesh%sb%dim
    if(scf%mix_field /= MIXNONE) then
      if(.not. hm%cmplxscl) then
        call mix_init(scf%smix, scf%mixdim1, scf%mixdim2, st%d%nspin)
      else
        call mix_init(scf%smix, scf%mixdim1, scf%mixdim2, st%d%nspin, func_type = TYPE_CMPLX)
      end if
    end if

    ! now the eigensolver stuff
    call eigensolver_init(scf%eigens, gr, st)

    if(scf%eigens%es_type == RS_MG .or. preconditioner_is_multigrid(scf%eigens%pre)) then
      if(.not. associated(gr%mgrid)) then
        SAFE_ALLOCATE(gr%mgrid)
        call multigrid_init(gr%mgrid, geo, gr%cv,gr%mesh, gr%der, gr%stencil)
      end if
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
    call parse_logical(datasets_check('SCFinLCAO'), .false., scf%lcao_restricted)
    if(scf%lcao_restricted) then
      message(1) = 'Info: SCF restricted to LCAO subspace.'
      call messages_info(1)
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
    call parse_logical(datasets_check('SCFCalculateForces'), .not. geo%only_user_def, scf%calc_force)

    !%Variable SCFCalculateDipole
    !%Type logical
    !%Section SCF
    !%Description
    !% This variable controls whether the dipole is calculated at the
    !% end of a self-consistent iteration. For finite systems the
    !% default is yes. For periodic systems the default is no, unless
    !% an electric field is being applied in a periodic direction.
    !% The single-point Berry`s phase approximation is used for
    !% periodic directions.
    !%End
    call parse_logical(datasets_check('SCFCalculateDipole'), .not. simul_box_is_periodic(gr%sb), scf%calc_dipole)
    if(associated(hm%vberry)) scf%calc_dipole = .true.

    rmin = geometry_min_distance(geo)
    if(geo%natoms == 1) rmin = CNST(100.0)

    !%Variable LocalMagneticMomentsSphereRadius
    !%Type float
    !%Section Output
    !%Description
    !% The local magnetic moments are calculated by integrating the
    !% magnetization density in spheres centered around each atom.
    !% This variable controls the radius of the spheres.
    !% The default is half the minimum distance between two atoms
    !% in the input coordinates.
    !%End
    call parse_float(datasets_check('LocalMagneticMomentsSphereRadius'), &
      units_from_atomic(units_inp%length, rmin * M_HALF), scf%lmm_r)
    ! this variable is also used in td/td_write.F90
    scf%lmm_r = units_to_atomic(units_inp%length, scf%lmm_r)

    POP_SUB(scf_init)
  end subroutine scf_init


  ! ---------------------------------------------------------
  subroutine scf_end(scf)
    type(scf_t), intent(inout) :: scf

    PUSH_SUB(scf_end)

    call eigensolver_end(scf%eigens)
    if(scf%mix_field /= MIXNONE) call mix_end(scf%smix)

    POP_SUB(scf_end)
  end subroutine scf_end


  ! ---------------------------------------------------------
  subroutine scf_mix_clear(scf)
    type(scf_t), intent(inout) :: scf

    PUSH_SUB(scf_mix_clear)

    call mix_clear(scf%smix)

    POP_SUB(scf_mix_clear)
  end subroutine scf_mix_clear


  ! ---------------------------------------------------------
  subroutine scf_run(scf, gr, geo, st, ks, hm, outp, gs_run, verbosity, iters_done)
    type(scf_t),          intent(inout) :: scf !< self consistent cycle
    type(grid_t),         intent(inout) :: gr !< grid
    type(geometry_t),     intent(inout) :: geo !< geometry
    type(states_t),       intent(inout) :: st !< States
    type(v_ks_t),         intent(inout) :: ks !< Kohn-Sham
    type(hamiltonian_t),  intent(inout) :: hm !< Hamiltonian
    type(output_t),       intent(in)    :: outp
    logical, optional,    intent(in)    :: gs_run
    integer, optional,    intent(in)    :: verbosity 
    integer, optional,    intent(out)   :: iters_done

    type(lcao_t) :: lcao    !< Linear combination of atomic orbitals
    type(profile_t), save :: prof

    integer :: iter, is, idim, iatom, nspin, err, iberry, idir
    FLOAT :: evsum_out, evsum_in, forcetmp, dipole(MAX_DIM), dipole_prev(MAX_DIM)
    real(8) :: etime, itime
    FLOAT, allocatable :: rhoout(:,:,:), rhoin(:,:,:), rhonew(:,:,:)
    FLOAT, allocatable :: vout(:,:,:), vin(:,:,:), vnew(:,:,:)
    FLOAT, allocatable :: forceout(:,:), forcein(:,:), forcediff(:), tmp(:)
    CMPLX, allocatable :: zrhoout(:,:,:), zrhoin(:,:,:), zrhonew(:,:,:)
    FLOAT, allocatable :: Imvout(:,:,:), Imvin(:,:,:), Imvnew(:,:,:)
    character(len=8) :: dirname
    logical :: finish, forced_finish, gs_run_, berry_conv, cmplxscl
    integer :: verbosity_

    PUSH_SUB(scf_run)
    
    cmplxscl = hm%cmplxscl
  
    gs_run_ = .true.
    if(present(gs_run)) gs_run_ = gs_run
    
    verbosity_ = VERB_FULL
    if(present(verbosity)) verbosity_ = verbosity

    if(ks%theory_level == CLASSICAL) then
      ! calculate forces
      if(scf%calc_force) call forces_calculate(gr, geo, hm%ep, st)
      
      if(gs_run_) then 
        ! output final information
        call scf_write_static(STATIC_DIR, "info")
        call output_all(outp, gr, geo, st, hm, ks%xc, STATIC_DIR)
      end if

      POP_SUB(scf_run)
      return
    end if

    if(scf%lcao_restricted) then
      call lcao_init(lcao, gr, geo, st)
      if(.not. lcao_is_available(lcao)) then
        message(1) = 'LCAO is not available. Cannot do SCF in LCAO.'
        call messages_fatal(1)
      end if
    end if

    nspin = st%d%nspin

    if(.not. cmplxscl) then
      
      SAFE_ALLOCATE(rhoout(1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))
      SAFE_ALLOCATE(rhoin (1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))

      rhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)
      rhoout = M_ZERO
    else
  
      SAFE_ALLOCATE(zrhoout(1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))
      SAFE_ALLOCATE(zrhoin (1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))

      zrhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%zrho%Im(1:gr%fine%mesh%np, 1:nspin)
      zrhoout = M_z0
    end if

    if (st%d%cdft) then
      rhoin(1:gr%fine%mesh%np, 2:scf%mixdim2, 1:nspin) = st%current(1:gr%fine%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
    end if
    
    select case(scf%mix_field)
    case(MIXPOT)
      SAFE_ALLOCATE(vout(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))
      SAFE_ALLOCATE( vin(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))
      SAFE_ALLOCATE(vnew(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))

      vin(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
      vout = M_ZERO
      if (st%d%cdft) vin(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      if(cmplxscl) then
        SAFE_ALLOCATE(Imvout(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))
        SAFE_ALLOCATE( Imvin(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))
        SAFE_ALLOCATE(Imvnew(1:gr%mesh%np, 1:scf%mixdim2, 1:nspin))

        Imvin(1:gr%mesh%np, 1, 1:nspin) = hm%Imvhxc(1:gr%mesh%np, 1:nspin)
        Imvout = M_ZERO
      end if
    case(MIXDENS)
      if(.not. cmplxscl) then
        SAFE_ALLOCATE(rhonew(1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))
      else
        SAFE_ALLOCATE(zrhonew(1:gr%fine%mesh%np, 1:scf%mixdim2, 1:nspin))
      end if
    end select

    evsum_in = states_eigenvalues_sum(st)

    ! allocate and compute forces only if they are used as convergence criteria
    if (scf%conv_abs_force > M_ZERO) then
      SAFE_ALLOCATE(  forcein(1:geo%natoms, 1:gr%sb%dim))
      SAFE_ALLOCATE( forceout(1:geo%natoms, 1:gr%sb%dim))
      SAFE_ALLOCATE(forcediff(1:gr%sb%dim))
      call forces_calculate(gr, geo, hm%ep, st)
      do iatom = 1, geo%natoms
        forcein(iatom, 1:gr%sb%dim) = geo%atom(iatom)%f(1:gr%sb%dim)
      end do
    endif

    if ( verbosity_ /= VERB_NO ) then
      if(scf%max_iter > 0) then
        write(message(1),'(a)') 'Info: Starting SCF iteration.'
      else
        write(message(1),'(a)') 'Info: No SCF iterations will be done.'
        finish = .true.
      endif
      call messages_info(1)
    end if

    ! SCF cycle
    itime = loct_clock()
    do iter = 1, scf%max_iter
      call profiling_in(prof, "SCF_CYCLE")

      if(scf%lcao_restricted) then
        call lcao_wf(lcao, st, gr, geo, hm)
      else
        ! FIXME: Currently, only the eigensolver or the
        ! Lippmann-Schwinger approach can be used (exclusively),
        ! i.e. no bound states for open boundaries.
        if(gr%ob_grid%open_boundaries) then
          call lippmann_schwinger(scf%eigens, hm, gr, st)
        else
          if(associated(hm%vberry)) then
            ks%frozen_hxc = .true.
            do iberry = 1, scf%max_iter_berry
              scf%eigens%converged = 0
              call eigensolver_run(scf%eigens, gr, st, hm, iter)

              call v_ks_calc(ks, hm, st, geo)
              call hamiltonian_update(hm, gr%mesh)

              dipole_prev = dipole
              call calc_dipole(dipole)
              write(message(1),'(a,9f12.6)') 'Dipole = ', dipole(1:gr%sb%dim)
              call messages_info(1)

              berry_conv = .true.
              do idir = 1, gr%sb%periodic_dim
                berry_conv = berry_conv .and. &
                  (abs((dipole(idir) - dipole_prev(idir)) / dipole_prev(idir)) < CNST(1e-5) &
                  .or. abs(dipole(idir) - dipole_prev(idir)) < CNST(1e-5))
              enddo
              if(berry_conv) exit
            enddo
            ks%frozen_hxc = .false.
          else
            scf%eigens%converged = 0
            call eigensolver_run(scf%eigens, gr, st, hm, iter)
          endif
        end if
      end if

      ! occupations
      call states_fermi(st, gr%mesh)

      ! compute output density, potential (if needed) and eigenvalues sum
      if(cmplxscl) then
        call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
!         print *,"Density integral", sum(st%zrho%Re(:,1) + M_zI * st%zrho%Im(:,1))*gr%mesh%volume_element
        print *,"Density integral", zmf_integrate(gr%mesh, st%zrho%Re(:,1) + M_zI * st%zrho%Im(:,1))
      else
        call density_calc(st, gr, st%rho)
      end if
      
      if(.not. cmplxscl) then
        rhoout(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)
      else
        zrhoout(1:gr%fine%mesh%np, 1, 1:nspin) = st%zrho%Re(1:gr%fine%mesh%np, 1:nspin) +&
                                            M_zI * st%zrho%Im(1:gr%fine%mesh%np, 1:nspin)
      end if
      
      if (hm%d%cdft) then
        call calc_physical_current(gr%der, st, st%current)
        rhoout(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin) = st%current(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      if (scf%mix_field == MIXPOT) then
        call v_ks_calc(ks, hm, st, geo)
        vout(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
        if(cmplxscl) Imvout(1:gr%mesh%np, 1, 1:nspin) = hm%Imvhxc(1:gr%mesh%np, 1:nspin)
        if (hm%d%cdft) vout(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      evsum_out = states_eigenvalues_sum(st)

      ! recalculate total energy
      call energy_calc_total(hm, gr, st, iunit = 0)

      ! compute convergence criteria
      scf%abs_dens = M_ZERO
      SAFE_ALLOCATE(tmp(1:gr%fine%mesh%np))
      do is = 1, nspin
        do idim = 1, scf%mixdim2
          if(.not. cmplxscl) then
            tmp = abs(rhoin(1:gr%fine%mesh%np, idim, is) - rhoout(1:gr%fine%mesh%np, idim, is))
          else
            tmp = abs(zrhoin(1:gr%fine%mesh%np, idim, is) - zrhoout(1:gr%fine%mesh%np, idim, is))
          end if
          scf%abs_dens = scf%abs_dens + dmf_integrate(gr%fine%mesh, tmp)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp)

      ! compute forces only if they are used as convergence criterion
      if (scf%conv_abs_force > M_ZERO) then
        call forces_calculate(gr, geo, hm%ep, st)
        scf%abs_force = M_ZERO
        do iatom = 1, geo%natoms
          forceout(iatom,1:gr%sb%dim) = geo%atom(iatom)%f(1:gr%sb%dim)
          forcediff(1:gr%sb%dim) = abs( forceout(iatom,1:gr%sb%dim) - forcein(iatom,1:gr%sb%dim) )
          forcetmp = maxval( forcediff )
          if ( forcetmp > scf%abs_force ) then
            scf%abs_force = forcetmp
          end if
        end do
      end if

      scf%rel_dens = scf%abs_dens / st%qtot
      scf%abs_ev = abs(evsum_out - evsum_in)
      scf%rel_ev = scf%abs_ev / abs(evsum_out)

      ! are we finished?
      finish = &
        (scf%conv_abs_dens  <= M_ZERO .or. scf%abs_dens  <= scf%conv_abs_dens)  .and. &
        (scf%conv_rel_dens  <= M_ZERO .or. scf%rel_dens  <= scf%conv_rel_dens)  .and. &
        (scf%conv_abs_force <= M_ZERO .or. scf%abs_force <= scf%conv_abs_force) .and. &
        (scf%conv_abs_ev    <= M_ZERO .or. scf%abs_ev    <= scf%conv_abs_ev)    .and. &
        (scf%conv_rel_ev    <= M_ZERO .or. scf%rel_ev    <= scf%conv_rel_ev)

      etime = loct_clock() - itime
      itime = etime + itime
      call scf_write_iter()

      ! mixing
      select case (scf%mix_field)
      case (MIXDENS)
        !set the pointer for dmf_dotp_aux
        call mesh_init_mesh_aux(gr%fine%mesh)
        ! mix input and output densities and compute new potential
        if(.not. cmplxscl) then
        call dmixing(scf%smix, iter, rhoin, rhoout, rhonew, dmf_dotp_aux)
          st%rho(1:gr%fine%mesh%np, 1:nspin) = rhonew(1:gr%fine%mesh%np, 1, 1:nspin)
        else
          call zmixing(scf%smix, iter, zrhoin, zrhoout, zrhonew, zmf_dotp_aux)
          st%zrho%Re(1:gr%fine%mesh%np, 1:nspin) =  real(zrhonew(1:gr%fine%mesh%np, 1, 1:nspin))                   
          st%zrho%Im(1:gr%fine%mesh%np, 1:nspin) = aimag(zrhonew(1:gr%fine%mesh%np, 1, 1:nspin))                    
        end if
        if (hm%d%cdft) st%current(1:gr%mesh%np,1:gr%mesh%sb%dim,1:nspin) = rhonew(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin)
        call v_ks_calc(ks, hm, st, geo)
      case (MIXPOT)
        !set the pointer for dmf_dotp_aux
        call mesh_init_mesh_aux(gr%mesh)
        ! mix input and output potentials
        call dmixing(scf%smix, iter, vin, vout, vnew, dmf_dotp_aux)
        hm%vhxc(1:gr%mesh%np, 1:nspin) = vnew(1:gr%mesh%np, 1, 1:nspin)
        if(cmplxscl) then
          call dmixing(scf%smix, iter, Imvin, Imvout, Imvnew, dmf_dotp_aux)
          hm%Imvhxc(1:gr%mesh%np, 1:nspin) = Imvnew(1:gr%mesh%np, 1, 1:nspin)
        end if
        call hamiltonian_update(hm, gr%mesh)
        if (hm%d%cdft) hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin) = vnew(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin)
      case(MIXNONE)
        call v_ks_calc(ks, hm, st, geo)
      end select

      ! Are we asked to stop? (Whenever Fortran is ready for signals, this should go away)
      forced_finish = clean_stop()

      if(gs_run_) then 
        ! save restart information
        if(finish .or. (modulo(iter, outp%iter) == 0) .or. iter == scf%max_iter .or. forced_finish) then
          call restart_write(trim(tmpdir) // GS_DIR, st, gr, geo, err, iter=iter)
          if(err .ne. 0) then
            message(1) = 'Unsuccessful write of "'//trim(tmpdir)//GS_DIR//'"'
            call messages_fatal(1)
          end if
        end if
      end if

      if(finish) then
        if(present(iters_done)) iters_done = iter
        if(verbosity_ >= VERB_COMPACT) then
          write(message(1), '(a, i4, a)') 'Info: SCF converged in ', iter, ' iterations'
          write(message(2), '(a)')        '' 
          call messages_info(2)
        end if
        if(scf%lcao_restricted) call lcao_end(lcao)
        call profiling_out(prof)
        exit
      end if

      if(outp%duringscf .and. gs_run_) then
        write(dirname,'(a,i4.4)') "scf.",iter
        call output_all(outp, gr, geo, st, hm, ks%xc, dirname)
      end if

      ! save information for the next iteration
      if(.not. cmplxscl) then
        rhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%rho(1:gr%fine%mesh%np, 1:nspin)
      else
        zrhoin(1:gr%fine%mesh%np, 1, 1:nspin) = st%zrho%Re(1:gr%fine%mesh%np, 1:nspin) +&
                                          M_zI * st%zrho%Im(1:gr%fine%mesh%np, 1:nspin)  
      end if
      if (hm%d%cdft) rhoin(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin) = st%current(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      if (scf%mix_field == MIXPOT) then
        vin(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
        if (cmplxscl) Imvin(1:gr%mesh%np, 1, 1:nspin) = hm%Imvhxc(1:gr%mesh%np, 1:nspin)
        if (hm%d%cdft) vin(1:gr%mesh%np, 2:scf%mixdim2, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      evsum_in = evsum_out
      if (scf%conv_abs_force > M_ZERO) then
        forcein(1:geo%natoms, 1:gr%sb%dim) = forceout(1:geo%natoms, 1:gr%sb%dim)
      end if

      if(forced_finish) then
        call profiling_out(prof)
        exit
      end if

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(prof)
    end do !iter
    
    select case(scf%mix_field)
    case(MIXPOT)
      call v_ks_calc(ks, hm, st, geo)
      SAFE_DEALLOCATE_A(vout)
      SAFE_DEALLOCATE_A(vin)
      SAFE_DEALLOCATE_A(vnew)
      SAFE_DEALLOCATE_A(Imvout)
      SAFE_DEALLOCATE_A(Imvin)
      SAFE_DEALLOCATE_A(Imvnew)
    case(MIXDENS)
      SAFE_DEALLOCATE_A(rhonew)
      SAFE_DEALLOCATE_A(zrhonew)
    case(MIXNONE)
!       call v_ks_calc(ks, hm, st, geo)
    end select

    SAFE_DEALLOCATE_A(rhoout)
    SAFE_DEALLOCATE_A(rhoin)
    SAFE_DEALLOCATE_A(zrhoout)
    SAFE_DEALLOCATE_A(zrhoin)


    if(.not.finish) then
      write(message(1), '(a,i4,a)') 'SCF *not* converged after ', iter - 1, ' iterations.'
      call messages_warning(1)
    end if

    ! calculate forces
    if(scf%calc_force) call forces_calculate(gr, geo, hm%ep, st)

    if(gs_run_) then 
      ! output final information
      call scf_write_static(STATIC_DIR, "info")
      call output_all(outp, gr, geo, st, hm, ks%xc, STATIC_DIR)

      ! write part of the source term s(0)
      if(gr%ob_grid%open_boundaries) call states_write_proj_lead_wf(gr%sb, 'open_boundaries/', gr%intf, st)
    end if

    if(simul_box_is_periodic(gr%sb) .and. st%d%nik > st%d%nspin) then
      call states_write_bands(STATIC_DIR, st%nst, st, gr%sb)
      call states_write_fermi_energy(STATIC_DIR, st, gr%mesh, gr%sb)
    end if

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
        if(cmplxscl) then
          write(message(1),'(a,es15.8,2(a,es9.2))') ' Re(etot) = ', units_from_atomic(units_out%energy, hm%energy%total), &
               ' abs_ev   = ', units_from_atomic(units_out%energy, scf%abs_ev), ' rel_ev   = ', scf%rel_ev
          write(message(2),'(a,es15.8,2(a,es9.2))') ' Im(etot) = ', units_from_atomic(units_out%energy, hm%energy%Imtotal), &
               ' abs_dens = ', scf%abs_dens, ' rel_dens = ', scf%rel_dens
        else
          write(message(1),'(a,es15.8,2(a,es9.2))') ' etot = ', units_from_atomic(units_out%energy, hm%energy%total), &
               ' abs_ev   = ', units_from_atomic(units_out%energy, scf%abs_ev), ' rel_ev   = ', scf%rel_ev
          write(message(2),'(23x,2(a,es9.2))') &
               ' abs_dens = ', scf%abs_dens, ' rel_dens = ', scf%rel_dens
        end if      
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
          call states_write_eigenvalues(stdout, st%nst, st, gr%sb, scf%eigens%diff)
        else
          call states_write_eigenvalues(stdout, st%nst, st, gr%sb)
        end if

        if(st%smear%method .ne. SMEAR_SEMICONDUCTOR .and. st%smear%method .ne. SMEAR_FIXED_OCC) then
          write(message(1), '(a,f12.6,a)') "Fermi energy = ", units_from_atomic(units_out%energy, st%smear%e_fermi)
          call messages_info(1)
        endif

        if(associated(hm%vberry)) then
          call calc_dipole(dipole)
          call write_dipole(stdout, dipole)
        endif

        if(st%d%ispin > UNPOLARIZED) then
          call write_magnetic_moments(stdout, gr%mesh, st)
        end if

        write(message(1),'(a)') ''
        write(message(2),'(a,i5,a,f14.2)') 'Elapsed time for SCF step ', iter,':', etime
        call messages_info(2)

        if(conf%report_memory) then
          mem = get_memory_usage()/(CNST(1024.0)**2)
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

      integer :: iunit, idir, iatom, ii

      PUSH_SUB(scf_run.scf_write_static)

      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir)
        iunit = io_open(trim(dir) // "/" // trim(fname), action='write')

        call grid_write_info(gr, geo, iunit)

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

        call states_write_eigenvalues(iunit, st%nst, st, gr%sb)
        if(st%smear%method .ne. SMEAR_SEMICONDUCTOR .and. st%smear%method .ne. SMEAR_FIXED_OCC) then
          write(message(1), '(a,f12.6,a)') "Fermi energy = ", units_from_atomic(units_out%energy, st%smear%e_fermi)
          call messages_info(1, iunit)
        endif
        write(iunit, '(1x)')

        if(cmplxscl .and. hm%energy%Imtotal < M_ZERO) then
          write(message(1), '(3a,es18.6), ')'Lifetime [',trim(units_abbrev(units_out%time)), '] = ', & 
            units_from_atomic(units_out%time, - M_ONE/(M_TWO * hm%energy%Imtotal))
          call messages_info(1, iunit)
        end if

        write(iunit, '(3a)') 'Energy [', trim(units_abbrev(units_out%energy)), ']:'
      else
        iunit = 0
      end if

      call energy_calc_total(hm, gr, st, iunit, full = .true.)

      if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      if(st%d%ispin > UNPOLARIZED) then
        call write_magnetic_moments(iunit, gr%mesh, st)
        if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      end if

      if(scf%calc_dipole) then
        call calc_dipole(dipole)
        call write_dipole(iunit, dipole)
      end if

      if(mpi_grp_is_root(mpi_world)) then
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

        if(scf%calc_force) then
          write(iunit,'(3a)') 'Forces on the ions [', trim(units_abbrev(units_out%force)), "]"
          write(iunit,'(a,10x,99(14x,a))') ' Ion', (index2axis(idir), idir = 1, gr%sb%dim)
          do iatom = 1, geo%natoms
            write(iunit,'(i4,a10,10f15.6)') iatom, trim(species_label(geo%atom(iatom)%spec)), &
              (units_from_atomic(units_out%force, geo%atom(iatom)%f(idir)), idir=1, gr%sb%dim)
          end do
          write(iunit,'(1x,100a1)') ("-", ii = 1, 13 + gr%sb%dim * 15)
          write(iunit,'(a14, 10f15.6)') " Max abs force", &
            (units_from_atomic(units_out%force, maxval(abs(geo%atom(1:geo%natoms)%f(idir)))), idir=1, gr%sb%dim)
          write(iunit,'(a14, 10f15.6)') " Total force", &
            (units_from_atomic(units_out%force, sum(geo%atom(1:geo%natoms)%f(idir))), idir=1, gr%sb%dim)
        end if

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
        if(idir .le. gr%sb%periodic_dim) then
          dipole(idir) = -n_dip(idir) - berry_dipole(st, gr%mesh, idir)

          ! use quantum of polarization to reduce to smallest possible magnitude
          nquantumpol = NINT(dipole(idir)/(CNST(2.0)*gr%sb%lsize(idir)))
          dipole(idir) = dipole(idir) - nquantumpol * (CNST(2.0) * gr%sb%lsize(idir))

        ! in aperiodic directions use normal dipole formula
        else
          e_dip(idir + 1, 1) = sum(e_dip(idir + 1, :))
          dipole(idir) = -n_dip(idir) - e_dip(idir + 1, 1)
        endif
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

          if (st%d%nik * st%smear%el_per_state .ne. 2) then
            write(iunit, '(a)') &
              "WARNING: Single-point Berry's phase method for dipole should not be used when there is more than one k-point."
            write(iunit, '(a)') &
              "Instead, finite differences on k-points (not yet implemented) are needed."
          endif

          if(.not. smear_is_semiconducting(st%smear)) then
            write(iunit, '(a)') "Single-point Berry's phase dipole calculation not correct without integer occupations."
          endif
        endif

        write(iunit, *)
      endif

      POP_SUB(scf_run.write_dipole)
    end subroutine write_dipole


    ! ---------------------------------------------------------
    subroutine write_magnetic_moments(iunit, mesh, st)
      integer,        intent(in) :: iunit
      type(mesh_t),   intent(in) :: mesh
      type(states_t), intent(in) :: st

      integer :: ia
      FLOAT :: mm(max(mesh%sb%dim, 3))
      FLOAT, allocatable :: lmm(:,:)

      PUSH_SUB(scf_run.write_magnetic_moments)

      call magnetic_moment(mesh, st, st%rho, mm)
      SAFE_ALLOCATE(lmm(1:max(mesh%sb%dim, 3), 1:geo%natoms))
      call magnetic_local_moments(mesh, st, geo, st%rho, scf%lmm_r, lmm)

      if(mpi_grp_is_root(mpi_world)) then

        write(iunit, '(a)') 'Total Magnetic Moment:'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit, '(a,f10.6)') ' mz = ', mm(3)
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit, '(1x,3(a,f10.6,3x))') 'mx = ', mm(1),'my = ', mm(2),'mz = ', mm(3)
        end if

        write(iunit, '(a,a,a,f7.3,a)') 'Local Magnetic Moments (sphere radius [', &
             trim(units_abbrev(units_out%length)),'] = ', units_from_atomic(units_out%length, scf%lmm_r), '):'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit,'(a,6x,14x,a)') ' Ion','mz'
          do ia = 1, geo%natoms
            write(iunit,'(i4,a10,f15.6)') ia, trim(species_label(geo%atom(ia)%spec)), lmm(3, ia)
          end do
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit,'(a,8x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
          do ia = 1, geo%natoms
            write(iunit,'(i4,a10,9f15.6)') ia, trim(species_label(geo%atom(ia)%spec)), lmm(1:mesh%sb%dim, ia)
          end do
        end if

      end if
      
      SAFE_DEALLOCATE_A(lmm)

      POP_SUB(scf_run.write_magnetic_moments)
    end subroutine write_magnetic_moments

  end subroutine scf_run

end module scf_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
