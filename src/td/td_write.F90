! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module td_write_oct_m
  use blas_oct_m
  use clock_oct_m
  use comm_oct_m
  use current_oct_m
  use density_oct_m
  use derivatives_oct_m
  use excited_states_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use output_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_mxll_oct_m
  use io_function_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use ions_oct_m
  use iso_c_binding
  use kick_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use lalg_adv_oct_m
  use lda_u_oct_m
  use loct_math_oct_m
  use magnetic_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use partial_charges_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use states_mxll_oct_m
  use states_mxll_restart_oct_m
  use td_calc_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::         &
    td_write_t,     &
    td_write_init,  &
    td_write_end,   &
    td_write_iter,  &
    td_write_data,  &
    td_write_kick,  &
    td_write_output, &
    td_write_mxll_end, &
    td_write_mxll_init, &
    td_write_mxll_iter, &
    td_write_mxll_free_data

  ! The following routine is currently unused, but will be used in the near future.
  ! In order noy to generate warnings about it, we declared it as public
  public :: td_dump_mxll

  integer, parameter ::   &
    OUT_MULTIPOLES  =  1, &
    OUT_ANGULAR     =  2, &
    OUT_SPIN        =  3, &
    OUT_POPULATIONS =  4, &
    OUT_COORDS      =  5, &
    OUT_ACC         =  6, & 
    OUT_LASER       =  7, &
    OUT_ENERGY      =  8, &
    OUT_PROJ        =  9, &
    OUT_MAGNETS     = 10, &
    OUT_GAUGE_FIELD = 11, &
    OUT_TEMPERATURE = 12, &
    OUT_FTCHD       = 13, &
    OUT_VEL         = 14, &
    OUT_EIGS        = 15, &
    OUT_ION_CH      = 16, &
    OUT_TOTAL_CURRENT = 17, &
    OUT_PARTIAL_CHARGES = 18, &
    OUT_KP_PROJ     = 19, &
    OUT_FLOQUET     = 20, &
    OUT_N_EX        = 21, &
    OUT_SEPARATE_COORDS  = 22, &
    OUT_SEPARATE_VELOCITY= 23, &
    OUT_SEPARATE_FORCES  = 24, &
    OUT_TOTAL_HEAT_CURRENT = 25, &
    OUT_TOT_M            = 26, &
    OUT_MAX              = 26
  
  integer, parameter ::      &
    OUT_DFTU_EFFECTIVE_U = 1, &
    OUT_DFTU_MAX         = 1

  integer, parameter ::   &
    OUT_MAXWELL_ENERGY          = 1, &
    OUT_MAXWELL_FIELDS          = 2, &
    OUT_MEAN_POYNTING           = 3, &
    OUT_E_FIELD_SURFACE_X       = 4, &
    OUT_E_FIELD_SURFACE_Y       = 5, &
    OUT_E_FIELD_SURFACE_Z       = 6, &
    OUT_B_FIELD_SURFACE_X       = 7, &
    OUT_B_FIELD_SURFACE_Y       = 8, &
    OUT_B_FIELD_SURFACE_Z       = 9

    integer, parameter, public ::   &
    OUT_MAXWELL_MAX             = 9

  type td_write_prop_t
    type(c_ptr)              :: handle
    type(c_ptr), allocatable :: mult_handles(:)
    type(mpi_grp_t)          :: mpi_grp
    integer                  :: hand_start
    integer                  :: hand_end
    logical                  :: write = .false.
    logical                  :: resolve_states = .false.    !< Whether to resolve output by state
  end type td_write_prop_t


  type td_write_t
    private
    type(td_write_prop_t) :: out(OUT_MAX)
    type(td_write_prop_t) :: out_dftu(OUT_DFTU_MAX)

    integer        :: lmax     !< maximum multipole moment to output
    FLOAT          :: lmm_r    !< radius of the sphere used to compute the local magnetic moments
    !> The states_elec_type where the ground state is stored, in order to
    !! calculate the projections(s) onto it.
    type(states_elec_t) :: gs_st    
    integer        :: n_excited_states  !< number of excited states onto which the projections are calculated.
    type(excited_states_t), allocatable :: excited_st(:) !< The excited states.
    integer :: compute_interval     !< Compute every compute_interval
  end type td_write_t

contains


  ! ---------------------------------------------------------
  subroutine td_write_kick(outp, namespace, space, mesh, kick, ions, iter)
    type(output_t),    intent(in) :: outp
    type(namespace_t), intent(in) :: namespace
    type(space_t),     intent(in) :: space
    type(mesh_t),      intent(in) :: mesh
    type(kick_t),      intent(in) :: kick
    type(ions_t),      intent(in) :: ions
    integer,           intent(in) :: iter

    character(len=256) :: filename
    PUSH_SUB(td_write_kick)

    write(filename, '(a,i7.7)') "td.", iter  ! name of directory
    call output_kick(outp, namespace, space, filename, mesh, ions, kick)

    POP_SUB(td_write_kick)
  end subroutine td_write_kick


  ! ---------------------------------------------------------
  subroutine td_write_init(writ, namespace, space, outp, gr, st, hm, ions, ks, ions_move, with_gauge_field, kick, iter, max_iter, &
    dt, mc)
    type(td_write_t), target, intent(out)   :: writ
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(output_t),           intent(inout) :: outp
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(ions_t),             intent(in)    :: ions
    type(v_ks_t),             intent(inout) :: ks
    logical,                  intent(in)    :: ions_move
    logical,                  intent(in)    :: with_gauge_field
    type(kick_t),             intent(in)    :: kick
    integer,                  intent(in)    :: iter
    integer,                  intent(in)    :: max_iter
    FLOAT,                    intent(in)    :: dt
    type(multicomm_t),        intent(in)    :: mc

    FLOAT :: rmin
    integer :: ierr, first, ii, ist, jj, flags, iout, default
    logical :: output_options(MAX_OUTPUT_TYPES)
    integer :: output_interval(0:MAX_OUTPUT_TYPES)
    integer(8) :: how(0:MAX_OUTPUT_TYPES)
    type(block_t) :: blk
    character(len=MAX_PATH_LEN) :: filename
    type(restart_t) :: restart_gs
    logical :: resolve_states

    PUSH_SUB(td_write_init)

    !%Variable TDOutput
    !%Type block
    !%Default multipoles + energy (+ others depending on other options)
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent
    !% simulation. Many of the options can increase the computational
    !% cost of the simulation, so only use the ones that you need. In
    !% most cases the default value is enough, as it is adapted to the
    !% details of the TD run. If the ions are allowed to be moved, additionally
    !% the geometry and the temperature are output. If a laser is
    !% included it will output by default.
    !%
    !% Note: the output files generated by this option are updated
    !% every <tt>RestartWriteInterval</tt> steps.
    !%
    !% Example:
    !% <br><br><tt>%TDOutput
    !% <br>&nbsp;&nbsp;multipoles
    !% <br>&nbsp;&nbsp;energy
    !% <br>%<br></tt>
    !%
    !%Option multipoles 1
    !% Outputs the (electric) multipole moments of the density to the file <tt>td.general/multipoles</tt>.
    !% This is required to, <i>e.g.</i>, calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>TDMultipoleLmax</tt>.
    !%Option angular 2
    !% Outputs the orbital angular momentum of the system to <tt>td.general/angular</tt>, which can be used to calculate circular
    !% dichroism.
    !%Option spin 3
    !% (Experimental) Outputs the expectation value of the spin, which can be used to calculate magnetic
    !% circular dichroism.
    !%Option populations 4
    !% (Experimental) Outputs the projection of the time-dependent
    !% Kohn-Sham Slater determinant onto the ground state (or
    !% approximations to the excited states) to the file
    !% <tt>td.general/populations</tt>. Note that the calculation of
    !% populations is expensive in memory and computer time, so it
    !% should only be used if it is really needed. See <tt>TDExcitedStatesToProject</tt>.
    !%Option geometry 5
    !% If set (and if the atoms are allowed to move), outputs the coordinates, velocities,
    !% and forces of the atoms to the the file <tt>td.general/coordinates</tt>. On by default if <tt>MoveIons = yes</tt>.
    !%Option dipole_acceleration 6
    !% When set, outputs the acceleration of the electronic dipole, calculated from the Ehrenfest theorem,
    !% in the file <tt>td.general/acceleration</tt>. This file can then be
    !% processed by the utility <tt>oct-harmonic-spectrum</tt> in order to obtain the harmonic spectrum.
    !%Option laser 7
    !% If set, outputs the laser field to the file <tt>td.general/laser</tt>.
    !% On by default if <tt>TDExternalFields</tt> is set.
    !%Option energy 8
    !% If set, <tt>octopus</tt> outputs the different components of the energy
    !% to the file <tt>td.general/energy</tt>. Will be zero except for every <tt>TDEnergyUpdateIter</tt> iterations.
    !%Option td_occup 9
    !% (Experimental) If set, outputs the projections of the
    !% time-dependent Kohn-Sham wavefunctions onto the static
    !% (zero-time) wavefunctions to the file
    !% <tt>td.general/projections.XXX</tt>. Only use this option if
    !% you really need it, as it might be computationally expensive. See <tt>TDProjStateStart</tt>.
    !% The output interval of this quantity is controled by the variable <tt>TDOutputComputeInterval</tt>
    !% In case of states parallelization, all the ground-state states are stored by each task.
    !%Option local_mag_moments 10
    !% If set, outputs the local magnetic moments, integrated in sphere centered around each atom.
    !% The radius of the sphere can be set with <tt>LocalMagneticMomentsSphereRadius</tt>.
    !%Option gauge_field 11
    !% If set, outputs the vector gauge field corresponding to a spatially uniform (but time-dependent) 
    !% external electrical potential. This is only useful in a time-dependent periodic run.
    !% On by default if <tt>GaugeVectorField</tt> is set.
    !%Option temperature 12
    !% If set, the ionic temperature at each step is printed. On by default if <tt>MoveIons = yes</tt>.
    !%Option ftchd 13
    !% Write Fourier transform of the electron density to the file <tt>ftchd.X</tt>,
    !% where X depends on the kick (e.g. with sin-shaped perturbation X=sin).
    !% This is needed for calculating the dynamic structure factor.
    !% In the case that the kick mode is qbessel, the written quantity is integral over
    !% density, multiplied by spherical Bessel function times real spherical harmonic.
    !% On by default if <tt>TDMomentumTransfer</tt> is set.
    !%Option dipole_velocity 14
    !% When set, outputs the dipole velocity, calculated from the Ehrenfest theorem,
    !% in the file <tt>td.general/velocity</tt>. This file can then be
    !% processed by the utility <tt>oct-harmonic-spectrum</tt> in order to obtain the harmonic spectrum.
    !%Option eigenvalues 15
    !% Write the KS eigenvalues. 
    !%Option ionization_channels 16
    !% Write the multiple-ionization channels using the KS orbital densities as proposed in  
    !% C. Ullrich, Journal of Molecular Structure: THEOCHEM 501, 315 (2000).
    !%Option total_current 17
    !% Output the total current (average of the current density over the cell).
    !%Option partial_charges 18
    !% Bader and Hirshfeld partial charges. The output file is called 'td.general/partial_charges'.
    !%Option td_kpoint_occup 19                                                                  
    !% Project propagated Kohn-Sham states to the states at t=0 given in the directory 
    !% restart_proj (see %RestartOptions). This is an alternative to the option
    !% td_occup, with a formating more suitable for k-points and works only in 
    !% k- and/or state parallelization
    !%Option td_floquet 20
    !% Compute non-interacting Floquet bandstructure according to further options: 
    !% TDFloquetFrequency, TDFloquetSample, TDFloquetDimension.
    !% This is done only once per td-run at t=0.
    !% works only in k- and/or state parallelization
    !%Option n_excited_el 21
    !% Output the number of excited electrons, based on the projections 
    !% of the time evolved wave-functions on the ground-state wave-functions. 
    !% The output interval of this quantity is controled by the variable <tt>TDOutputComputeInterval</tt>
    !%Option coordinates_sep 22
    !% Writes geometries in a separate file.
    !%Option velocities_sep 23
    !% Writes velocities in a separate file.
    !%Option forces_sep 24
    !% Writes forces in a separate file.
    !%Option total_heat_current 25
    !% Output the total heat current (average of the heat current density over the cell).
    !%Option total_magnetization 26
    !% Writes the total magnetization, where the total magnetization is calculated at the momentum
    !% defined by <tt>TDMomentumTransfer</tt>. 
    !% This is used to extract the magnon frequency in case of a magnon kick.
    !%End



    !defaults
    output_options = .false.
    output_options(OUT_MULTIPOLES) = .true.
    output_options(OUT_ENERGY) = .true.
    if (ions_move) then
      output_options(OUT_COORDS) = .true.
      output_options(OUT_TEMPERATURE) = .true.
    end if
    if (with_gauge_field) output_options(OUT_GAUGE_FIELD) = .true.
    if (hm%ext_lasers%no_lasers > 0) output_options(OUT_LASER) = .true.
    if (kick%qkick_mode /= QKICKMODE_NONE) output_options(OUT_FTCHD) = .true.

    call io_function_read_what_how_when(namespace, space, output_options, how, output_interval, &
    'TDOutput')

    do iout = 1, OUT_MAX
      writ%out(iout)%write = output_options(iout)
    end do

    ! experimental stuff
    if(writ%out(OUT_SPIN)%write) call messages_experimental('TDOutput = spin')
    if(writ%out(OUT_POPULATIONS)%write) call messages_experimental('TDOutput = populations')
    if(writ%out(OUT_PROJ)%write) call messages_experimental('TDOutput = td_occup')
    if(writ%out(OUT_ION_CH)%write) call messages_experimental('TDOutput = ionization_channels')
    if(writ%out(OUT_PARTIAL_CHARGES)%write) call messages_experimental('TDOutput = partial_charges')
    if(writ%out(OUT_KP_PROJ)%write) call messages_experimental('TDOutput = td_kpoint_occup')
    if(writ%out(OUT_FLOQUET)%write) call messages_experimental('TDOutput = td_floquet')
    if(writ%out(OUT_N_EX)%write) call messages_experimental('TDOutput = n_excited_el')

    !See comment in zstates_elec_mpdotp
    if (space%is_periodic() .and. writ%out(OUT_POPULATIONS)%write) then
      call messages_not_implemented("TDOutput populations for periodic systems", namespace=namespace)
    end if


    if(writ%out(OUT_KP_PROJ)%write.or.writ%out(OUT_FLOQUET)%write) then
      ! make sure this is not domain distributed
      if(gr%mesh%np /= gr%mesh%np_global) then
        message(1) = "TDOutput option td_kpoint_occup and td_floquet do not work with domain parallelization"
        call messages_fatal(1, namespace=namespace)
      end if
    end if
    
    if ((writ%out(OUT_SEPARATE_FORCES)%write .or. writ%out(OUT_COORDS)%write) .and. space%periodic_dim == 1) then
      call messages_input_error(namespace, 'TDOutput', &
        'Forces for systems periodic in 1D are not currently implemented and options that output the forces are not allowed.')
    end if

    !%Variable TDOutputResolveStates
    !%Type logical
    !%Default No
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines whether the output should be resolved by states. 
    !% 
    !% So far only TDOutput = multipoles is supported.
    !%
    !%End
    call parse_variable(namespace, 'TDOutputResolveStates', .false., resolve_states)
    if(.not. writ%out(OUT_MULTIPOLES)%write) then
       write(message(1),'(a)') "TDOutputResolveStates works only for TDOutput = multipoles."
       call messages_warning(1)
    end if
    
    

    !%Variable TDMultipoleLmax
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::TD Output
    !%Description
    !% Maximum electric multipole of the density output to the file <tt>td.general/multipoles</tt>
    !% during a time-dependent simulation. Must be non-negative.
    !%End
    call parse_variable(namespace, 'TDMultipoleLmax', 1, writ%lmax)
    if (writ%lmax < 0) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid TDMultipoleLmax."
      message(2) = '(Must be TDMultipoleLmax >= 0 )'
      call messages_fatal(2, namespace=namespace)
    end if
    call messages_obsolete_variable(namespace, 'TDDipoleLmax', 'TDMultipoleLmax')

    ! Compatibility test
    if ((writ%out(OUT_ACC)%write) .and. ions_move) then
      message(1) = 'If harmonic spectrum is to be calculated, atoms should not be allowed to move.'
      call messages_fatal(1, namespace=namespace)
    end if

    rmin = ions%min_distance()

    ! This variable is documented in scf/scf.F90
    call parse_variable(namespace, 'LocalMagneticMomentsSphereRadius', min(M_HALF*rmin, CNST(100.0)), writ%lmm_r, &
      unit=units_inp%length)

    if(writ%out(OUT_PROJ)%write .or. writ%out(OUT_POPULATIONS)%write &
      .or.writ%out(OUT_KP_PROJ)%write .or. writ%out(OUT_N_EX)%write) then

      if(st%parallel_in_states .and. writ%out(OUT_POPULATIONS)%write) then
        message(1) = "Options TDOutput = populations are not implemented for parallel in states."
        call messages_fatal(1, namespace=namespace)
      end if

      if (writ%out(OUT_N_EX)%write .and. st%parallel_in_states) then
        message(1) = "Options TDOutput = n_excited_el is not implemented for parallel in states."
        call messages_fatal(1, namespace=namespace)
      end if
      
      if(.not.writ%out(OUT_KP_PROJ)%write.and..not.writ%out(OUT_N_EX)%write) then
        call states_elec_copy(writ%gs_st, st, exclude_wfns = .true., exclude_eigenval = .true.)
      else
        ! we want the same layout of gs_st as st
        call states_elec_copy(writ%gs_st, st)
      end if

      ! clean up all the stuff we have to reallocate
      SAFE_DEALLOCATE_A(writ%gs_st%node)

      call restart_init(restart_gs, namespace, RESTART_PROJ, RESTART_TYPE_LOAD, mc, ierr, mesh=gr%mesh)

      if(.not.writ%out(OUT_KP_PROJ)%write.and..not.writ%out(OUT_N_EX)%write) then
        if(ierr == 0) &
          call states_elec_look(restart_gs, ii, jj, writ%gs_st%nst, ierr)
          writ%gs_st%st_end = writ%gs_st%nst
        if(ierr /= 0) then
          message(1) = "Unable to read states information."
          call messages_fatal(1, namespace=namespace)
        end if
 
        ! do this only when not calculating populations, since all states are needed then
        if(.not. writ%out(OUT_POPULATIONS)%write) then
          ! We will store the ground-state Kohn-Sham system for all processors.
          !%Variable TDProjStateStart
          !%Type integer
          !%Default 1
          !%Section Time-Dependent::TD Output
          !%Description
          !% To be used with <tt>TDOutput = td_occup</tt>. Not available if <tt>TDOutput = populations</tt>.
          !% Only output projections to states above <tt>TDProjStateStart</tt>. Usually one is only interested
          !% in particle-hole projections around the HOMO, so there is no need to calculate (and store)
          !% the projections of all TD states onto all static states. This sets a lower limit. The upper limit
          !% is set by the number of states in the propagation and the number of unoccupied states
          !% available.
          !%End
          call parse_variable(namespace, 'TDProjStateStart', 1, writ%gs_st%st_start)
        else
           writ%gs_st%st_start = 1
        end if

        call states_elec_deallocate_wfns(writ%gs_st)

        writ%gs_st%parallel_in_states = .false.

        ! allocate memory
        SAFE_ALLOCATE(writ%gs_st%occ(1:writ%gs_st%nst, 1:writ%gs_st%d%nik))
        SAFE_ALLOCATE(writ%gs_st%eigenval(1:writ%gs_st%nst, 1:writ%gs_st%d%nik))

        !We want all the task to have all the states
        !States can be distibuted for the states we propagate.
        SAFE_ALLOCATE(writ%gs_st%node(1:writ%gs_st%nst))
        writ%gs_st%node(:)  = 0

        writ%gs_st%eigenval = huge(writ%gs_st%eigenval)
        writ%gs_st%occ      = M_ZERO
        if(writ%gs_st%d%ispin == SPINORS) then
          SAFE_DEALLOCATE_A(writ%gs_st%spin)
          SAFE_ALLOCATE(writ%gs_st%spin(1:3, 1:writ%gs_st%nst, 1:writ%gs_st%d%nik))
        end if

        call states_elec_allocate_wfns(writ%gs_st, gr%mesh, TYPE_CMPLX)
        
      end if
 
      call states_elec_load(restart_gs, namespace, space, writ%gs_st, gr%mesh, hm%kpoints, ierr, label = ': gs for TDOutput')

      if(ierr /= 0 .and. ierr /= (writ%gs_st%st_end-writ%gs_st%st_start+1)*writ%gs_st%d%nik &
                                      *writ%gs_st%d%dim*writ%gs_st%mpi_grp%size) then
        message(1) = "Unable to read wavefunctions for TDOutput."
        call messages_fatal(1, namespace=namespace)
      end if
      call restart_end(restart_gs)
    end if
 
    ! Build the excited states...
    if(writ%out(OUT_POPULATIONS)%write) then
      !%Variable TDExcitedStatesToProject
      !%Type block
      !%Section Time-Dependent::TD Output
      !%Description
      !% <b>[WARNING: This is a *very* experimental feature]</b>
      !% To be used with <tt>TDOutput = populations</tt>.
      !% The population of the excited states
      !% (as defined by <Phi_I|Phi(t)> where |Phi(t)> is the many-body time-dependent state at
      !% time <i>t</i>, and |Phi_I> is the excited state of interest) can be approximated -- it is not clear 
      !% how well -- by substituting for those real many-body states the time-dependent Kohn-Sham
      !% determinant and some modification of the Kohn-Sham ground-state determinant (<i>e.g.</i>,
      !% a simple HOMO-LUMO substitution, or the Casida ansatz for excited states in linear-response
      !% theory. If you set <tt>TDOutput</tt> to contain <tt>populations</tt>, you may ask for these approximated
      !% populations for a number of excited states, which will be described in the files specified
      !% in this block: each line should be the name of a file that contains one excited state.
      !%
      !% This file structure is the one written by the casida run mode, in the files in the directory <tt>*_excitations</tt>.
      !% The file describes the "promotions" from occupied
      !% to unoccupied levels that change the initial Slater determinant
      !% structure specified in ground_state. These promotions are a set
      !% of electron-hole pairs. Each line in the file, after an optional header, has four
      !% columns:
      !%
      !% <i>i  a  <math>\sigma</math> weight</i>
      !% 
      !% where <i>i</i> should be an occupied state, <i>a</i> an unoccupied one, and <math>\sigma</math>
      !% the spin of the corresponding orbital. This pair is then associated with a
      !% creation-annihilation pair <math>a^{\dagger}_{a,\sigma} a_{i,\sigma}</math>, so that the
      !% excited state will be a linear combination in the form:
      !% 
      !% <math>\left|{\rm ExcitedState}\right> =
      !% \sum weight(i,a,\sigma) a^{\dagger}_{a,\sigma} a_{i,\sigma} \left|{\rm GroundState}\right></math>
      !%
      !% where <i>weight</i> is the number in the fourth column.
      !% These weights should be normalized to one; otherwise the routine
      !% will normalize them, and write a warning.
      !%End
      if(parse_block(namespace, 'TDExcitedStatesToProject', blk) == 0) then
        writ%n_excited_states = parse_block_n(blk)
        SAFE_ALLOCATE(writ%excited_st(1:writ%n_excited_states))
        do ist = 1, writ%n_excited_states
          call parse_block_string(blk, ist-1, 0, filename)
          call excited_states_init(writ%excited_st(ist), writ%gs_st, trim(filename), namespace) 
        end do
      else
        writ%n_excited_states = 0
      end if
    end if

    !%Variable TDOutputComputeInterval
    !%Type integer
    !%Default 50
    !%Section Time-Dependent::TD Output
    !%Description
    !% The TD output requested are computed
    !% when the iteration number is a multiple of the <tt>TDOutputComputeInterval</tt> variable.
    !% Must be >= 0. If it is 0, then no output is written. 
    !% Implemented only for projections and number of excited electrons for the moment.
    !%End
    call parse_variable(namespace, 'TDOutputComputeInterval', 50, writ%compute_interval)
    if(writ%compute_interval < 0) then
      message(1) = "TDOutputComputeInterval must be >= 0."
      call messages_fatal(1, namespace=namespace)
    end if


    if (iter == 0) then
      first = 0
    else
      first = iter + 1
    end if

    call io_mkdir('td.general', namespace)
    

    ! default
    writ%out(:)%mpi_grp        = mpi_world
    
    if(writ%out(OUT_MULTIPOLES)%write .and. resolve_states) then       
      !resolve state contribution on multipoles
      writ%out(OUT_MULTIPOLES)%hand_start       = st%st_start
      writ%out(OUT_MULTIPOLES)%hand_end         = st%st_end
      writ%out(OUT_MULTIPOLES)%resolve_states   = .true.
      writ%out(OUT_MULTIPOLES)%mpi_grp          = gr%mesh%mpi_grp

      SAFE_ALLOCATE(writ%out(OUT_MULTIPOLES)%mult_handles(st%st_start:st%st_end))
      
      if (mpi_grp_is_root(writ%out(OUT_MULTIPOLES)%mpi_grp)) then
        do ist = st%st_start, st%st_end
          write(filename, '(a,i4.4)') 'td.general/multipoles-ist', ist
          call write_iter_init(writ%out(OUT_MULTIPOLES)%mult_handles(ist), &
                               first, units_from_atomic(units_out%time, dt), &
          trim(io_workpath(filename, namespace)))                    
        end do

      end if
    end if


    if(mpi_grp_is_root(mpi_world)) then

      if(writ%out(OUT_MULTIPOLES)%write .and. .not. resolve_states) then 
        call write_iter_init(writ%out(OUT_MULTIPOLES)%handle, &
                             first, units_from_atomic(units_out%time, dt), &
        trim(io_workpath("td.general/multipoles", namespace)))
      end if 

      if(writ%out(OUT_FTCHD)%write) then
        select case(kick%qkick_mode)
          case (QKICKMODE_SIN)
            write(filename, '(a)') 'td.general/ftchd.sin'
          case (QKICKMODE_COS)
            write(filename, '(a)') 'td.general/ftchd.cos'
          case (QKICKMODE_BESSEL)
            write(filename, '(a, SP, I0.3, a, I0.3)') 'td.general/ftchd.l', kick%qbessel_l, '_m', kick%qbessel_m
          case default
            write(filename, '(a)') 'td.general/ftchd'
        end select
        call write_iter_init(writ%out(OUT_FTCHD)%handle, &
          first, units_from_atomic(units_out%time, dt), trim(io_workpath(filename, namespace)))
      end if  

      if(writ%out(OUT_ANGULAR)%write) &
        call write_iter_init(writ%out(OUT_ANGULAR)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/angular", namespace)))

      if(writ%out(OUT_SPIN)%write) &
        call write_iter_init(writ%out(OUT_SPIN)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/spin", namespace)))

      if(writ%out(OUT_MAGNETS)%write) &
        call write_iter_init(writ%out(OUT_MAGNETS)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/magnetic_moments", namespace)))

      if(writ%out(OUT_COORDS)%write) &
        call write_iter_init(writ%out(OUT_COORDS)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/coordinates", namespace)))

      if(writ%out(OUT_SEPARATE_COORDS)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_COORDS)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/onlyCoordinates", namespace)))

      if(writ%out(OUT_SEPARATE_VELOCITY)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_VELOCITY)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/onlyVelocities", namespace)))

      if(writ%out(OUT_SEPARATE_FORCES)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_FORCES)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/onlyForces", namespace)))

      if(writ%out(OUT_TEMPERATURE)%write) &
        call write_iter_init(writ%out(OUT_TEMPERATURE)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/temperature", namespace)))

      if(writ%out(OUT_POPULATIONS)%write) &
        call write_iter_init(writ%out(OUT_POPULATIONS)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/populations", namespace)))

      if(writ%out(OUT_ACC)%write) &
        call write_iter_init(writ%out(OUT_ACC)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/acceleration", namespace)))
          
      if(writ%out(OUT_VEL)%write) &
        call write_iter_init(writ%out(OUT_VEL)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/velocity", namespace)))

      if(writ%out(OUT_LASER)%write) then
        ! The laser file is written for the full propagation in one go, so that
        ! the user can check that the laser is correct and as intended before letting
        ! the code run for a possibly large period of time. This is done even after
        ! a restart, so that it takes into account any changes to max_iter.
        call io_rm("td.general/laser", namespace=namespace)
        call write_iter_init(writ%out(OUT_LASER)%handle, 0, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/laser", namespace)))
        do ii = 0, max_iter
          call td_write_laser(writ%out(OUT_LASER)%handle, gr, hm, dt, ii)
        end do
        call write_iter_end(writ%out(OUT_LASER)%handle)
      end if

      if(writ%out(OUT_ENERGY)%write) &
        call write_iter_init(writ%out(OUT_ENERGY)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/energy", namespace)))

      if(writ%out(OUT_PROJ)%write) &
        call write_iter_init(writ%out(OUT_PROJ)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/projections", namespace)))

      if(writ%out(OUT_KP_PROJ)%write) &
        call write_iter_init(writ%out(OUT_KP_PROJ)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/projections", namespace)))

      if(writ%out(OUT_FLOQUET)%write) &
        call write_iter_init(writ%out(OUT_FLOQUET)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/floquet_bands", namespace)))

      if(writ%out(OUT_GAUGE_FIELD)%write) &
        call write_iter_init(writ%out(OUT_GAUGE_FIELD)%handle, &
        first, units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/gauge_field", namespace)))
        
      if(writ%out(OUT_EIGS)%write) &
        call write_iter_init(writ%out(OUT_EIGS)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/eigenvalues", namespace)))

      if(writ%out(OUT_ION_CH)%write) &
        call write_iter_init(writ%out(OUT_ION_CH)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/ion_ch", namespace)))

      if(writ%out(OUT_TOTAL_CURRENT)%write) then
        call write_iter_init(writ%out(OUT_TOTAL_CURRENT)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/total_current", namespace)))
      end if

      if(writ%out(OUT_TOTAL_HEAT_CURRENT)%write) then
        call write_iter_init(writ%out(OUT_TOTAL_HEAT_CURRENT)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/total_heat_current", namespace)))
      end if

      if(writ%out(OUT_PARTIAL_CHARGES)%write) then
        call write_iter_init(writ%out(OUT_PARTIAL_CHARGES)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/partial_charges", namespace)))
      end if

     if(writ%out(OUT_N_EX)%write) &
        call write_iter_init(writ%out(OUT_N_EX)%handle, first, &
          units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/n_ex", namespace)))
    
     if(writ%out(OUT_TOT_M)%write) &
        call write_iter_init(writ%out(OUT_TOT_M)%handle, first, &
          units_from_atomic(units_out%time, dt), &
          trim(io_workpath("td.general/total_magnetization", namespace)))
      
    end if
    
    if(writ%out(OUT_TOTAL_CURRENT)%write .or. writ%out(OUT_TOTAL_HEAT_CURRENT)%write) then
      !TODO: we should only compute the current here, not v_ks
      call v_ks_calculate_current(ks, .true.)
      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval=.false., time = iter*dt, calc_energy = .false.)
    end if

    if(writ%out(OUT_N_EX)%write .and. writ%compute_interval > 0) then
      call io_mkdir(outp%iter_dir, namespace)
    end if

    if (all(outp%how == 0) .and. writ%out(OUT_N_EX)%write) then
      call io_function_read_what_how_when(namespace, space, outp%what, outp%how, outp%output_interval)
    end if

    !%Variable TDOutputDFTU
    !%Type flag
    !%Default none
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent
    !% simulation, related to LDA+U.
    !%
    !% Note: the output files generated by this option are updated
    !% every <tt>RestartWriteInterval</tt> steps.
    !%Option effective_u 1
    !% Writes the effective U for each orbital set as a function of time.
    !%End
    default = 0
    if(hm%lda_u_level == DFT_U_ACBN0) default = default + 2**(OUT_DFTU_EFFECTIVE_U - 1)
    call parse_variable(namespace, 'TDOutputDFTU', default, flags)

    if(.not.varinfo_valid_option('TDOutputDFTU', flags, is_flag = .true.)) then
      call messages_input_error(namespace, 'TDOutputDFTU')
    end if

    do iout = 1, OUT_DFTU_MAX
      writ%out_dftu(iout)%write = (iand(flags, 2**(iout - 1)) /= 0)
    end do

    if(mpi_grp_is_root(mpi_world)) then
      if(writ%out_dftu(OUT_DFTU_EFFECTIVE_U)%write) &
        call write_iter_init(writ%out_dftu(OUT_DFTU_EFFECTIVE_U)%handle, &
          first, units_from_atomic(units_out%time, dt),  &
          trim(io_workpath("td.general/effectiveU", namespace))) 
    end if
    
    do iout = 1, OUT_MAX
      if(iout == OUT_MULTIPOLES) cycle
      if(writ%out(iout)%write)  writ%out(iout)%mpi_grp = mpi_world
    end do
    
     
    POP_SUB(td_write_init)
  end subroutine td_write_init


  ! ---------------------------------------------------------
  subroutine td_write_end(writ)
    type(td_write_t), intent(inout) :: writ

    integer :: ist, iout

    PUSH_SUB(td_write_end)

    do iout = 1, OUT_MAX
      if(iout == OUT_LASER) cycle      
      if(writ%out(iout)%write) then  
        if (mpi_grp_is_root(writ%out(iout)%mpi_grp)) then
          if(writ%out(iout)%resolve_states) then
            do ist = writ%out(iout)%hand_start, writ%out(iout)%hand_end
              call write_iter_end(writ%out(iout)%mult_handles(ist))
            end do
            SAFE_DEALLOCATE_A(writ%out(iout)%mult_handles)
          else
            call write_iter_end(writ%out(iout)%handle)
          end if
        end if
      end if
    end do
     
    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_DFTU_MAX
        if(writ%out_dftu(iout)%write) call write_iter_end(writ%out_dftu(iout)%handle)  
      end do
      
    end if

    if( writ%out(OUT_POPULATIONS)%write ) then
      do ist = 1, writ%n_excited_states
        call excited_states_kill(writ%excited_st(ist))
      end do
      writ%n_excited_states = 0
    end if

    if(writ%out(OUT_PROJ)%write.or.writ%out(OUT_POPULATIONS)%write &
       .or. writ%out(OUT_N_EX)%write) then
      call states_elec_end(writ%gs_st)
    end if

    POP_SUB(td_write_end)
  end subroutine td_write_end


  ! ---------------------------------------------------------
  subroutine td_write_iter(writ, namespace, space, outp, gr, st, hm, ions, kick, dt, iter)
    type(td_write_t),         intent(inout) :: writ !< Write object
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(output_t),           intent(in)    :: outp
    type(grid_t),             intent(in)    :: gr   !< The grid
    type(states_elec_t),      intent(inout) :: st   !< State object
    type(hamiltonian_elec_t), intent(inout) :: hm   !< Hamiltonian object
    type(ions_t),             intent(inout) :: ions  !< Geometry object
    type(kick_t),             intent(in)    :: kick !< The kick
    FLOAT,                    intent(in)    :: dt   !< Delta T, time interval
    integer,                  intent(in)    :: iter !< Iteration number

    type(profile_t), save :: prof

    PUSH_SUB(td_write_iter)
    call profiling_in(prof, "TD_WRITE_ITER")


    if (writ%out(OUT_MULTIPOLES)%write) then
      call td_write_multipole(writ%out(OUT_MULTIPOLES), gr, ions, st, writ%lmax, kick, iter)
    end if

    if (writ%out(OUT_FTCHD)%write) then
      call td_write_ftchd(writ%out(OUT_FTCHD)%handle, gr, st, kick, iter)
    end if

    if (writ%out(OUT_ANGULAR)%write) then
      call td_write_angular(writ%out(OUT_ANGULAR)%handle, namespace, gr, ions, hm, st, kick, iter)
    end if

    if (writ%out(OUT_SPIN)%write) then
      call td_write_spin(writ%out(OUT_SPIN)%handle, gr, st, iter)
    end if

    if (writ%out(OUT_MAGNETS)%write) then
      call td_write_local_magnetic_moments(writ%out(OUT_MAGNETS)%handle, gr, st, ions, writ%lmm_r, iter)
    end if

    if (writ%out(OUT_TOT_M)%write) then
      call td_write_tot_mag(writ%out(OUT_TOT_M)%handle, gr, st, kick, iter)
    end if

    if(writ%out(OUT_PROJ)%write .and. mod(iter, writ%compute_interval) == 0) then
      if (mpi_grp_is_root(mpi_world)) call write_iter_set(writ%out(OUT_PROJ)%handle, iter)
      call td_write_proj(writ%out(OUT_PROJ)%handle, space, gr, ions, st, writ%gs_st, kick, iter)
    end if

    if (writ%out(OUT_FLOQUET)%write) then
      call td_write_floquet(namespace, space, hm, gr, st, iter)
    end if

    if (writ%out(OUT_KP_PROJ)%write) then
      call td_write_proj_kp(gr, hm%kpoints, st, writ%gs_st, namespace, iter)
    end if

    if (writ%out(OUT_COORDS)%write) then
      call td_write_coordinates(writ%out(OUT_COORDS)%handle, ions, iter)
    end if

    if (writ%out(OUT_SEPARATE_COORDS)%write) then
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_COORDS)%handle, ions, iter,1)
    end if

    if (writ%out(OUT_SEPARATE_VELOCITY)%write) then
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_VELOCITY)%handle, ions, iter,2)
    end if

    if (writ%out(OUT_SEPARATE_FORCES)%write) then
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_FORCES)%handle, ions, iter,3)
    end if

    if (writ%out(OUT_TEMPERATURE)%write) then
      call td_write_temperature(writ%out(OUT_TEMPERATURE)%handle, ions, iter)
    end if

    if (writ%out(OUT_POPULATIONS)%write) then
      call td_write_populations(writ%out(OUT_POPULATIONS)%handle, namespace, space, gr%mesh, st, writ, dt, iter)
    end if

    if (writ%out(OUT_ACC)%write) then
      call td_write_acc(writ%out(OUT_ACC)%handle, namespace, space, gr, ions, st, hm, dt, iter)
    end if
      
    if (writ%out(OUT_VEL)%write) then
      call td_write_vel(writ%out(OUT_VEL)%handle, space, gr%der, st, hm%kpoints, iter)
    end if

    ! td_write_laser no longer called here, because the whole laser is printed
    ! out at the beginning.

    if (writ%out(OUT_ENERGY)%write) then
      call td_write_energy(writ%out(OUT_ENERGY)%handle, hm, iter, ions%kinetic_energy)
    end if

    if (writ%out(OUT_GAUGE_FIELD)%write) then
      call gauge_field_output_write(hm%ep%gfield, writ%out(OUT_GAUGE_FIELD)%handle, iter)
    end if

    if (writ%out(OUT_EIGS)%write) then
      call td_write_eigs(writ%out(OUT_EIGS)%handle, st, iter)
    end if

    if (writ%out(OUT_ION_CH)%write) then
      call td_write_ionch(writ%out(OUT_ION_CH)%handle, gr, st, iter)
    end if

    if (writ%out(OUT_TOTAL_CURRENT)%write) then
      call td_write_total_current(writ%out(OUT_TOTAL_CURRENT)%handle, gr, st, iter)
    end if
    
    if (writ%out(OUT_TOTAL_HEAT_CURRENT)%write) then
      call td_write_total_heat_current(writ%out(OUT_TOTAL_HEAT_CURRENT)%handle, space, hm, gr, st, iter)
    end if
    
    if (writ%out(OUT_PARTIAL_CHARGES)%write) then
      call td_write_partial_charges(writ%out(OUT_PARTIAL_CHARGES)%handle, gr%mesh, st, &
        ions, iter)
    end if
    
    if (writ%out(OUT_N_EX)%write .and. mod(iter, writ%compute_interval) == 0) then
      if (mpi_grp_is_root(mpi_world)) call write_iter_set(writ%out(OUT_N_EX)%handle, iter)
      call td_write_n_ex(writ%out(OUT_N_EX)%handle, outp, namespace, gr, hm%kpoints, st, writ%gs_st, iter)
    end if

    !LDA+U outputs
    if (writ%out_dftu(OUT_DFTU_EFFECTIVE_U)%write) then
      call td_write_effective_u(writ%out_dftu(OUT_DFTU_EFFECTIVE_U)%handle, hm%lda_u, iter)
    end if

    call profiling_out(prof)
    POP_SUB(td_write_iter)
  end subroutine td_write_iter


  ! ---------------------------------------------------------
  subroutine td_write_data(writ)
    type(td_write_t),     intent(inout) :: writ

    integer :: iout, ii
    type(profile_t), save :: prof

    PUSH_SUB(td_write_data)
    call profiling_in(prof, "TD_WRITE_DATA")

    do iout = 1, OUT_MAX
      if(iout == OUT_LASER) cycle
      if(writ%out(iout)%write) then
        if (mpi_grp_is_root(writ%out(iout)%mpi_grp)) then
          if (writ%out(iout)%resolve_states) then
            do ii = writ%out(iout)%hand_start, writ%out(iout)%hand_end
              call write_iter_flush(writ%out(iout)%mult_handles(ii))
            end do
          else
            call write_iter_flush(writ%out(iout)%handle)
          end if
        end if
      end if
    end do

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_DFTU_MAX
        if(writ%out_dftu(iout)%write) call write_iter_flush(writ%out_dftu(iout)%handle)
      end do
    end if

    call profiling_out(prof)
    POP_SUB(td_write_data)
  end subroutine td_write_data

  ! ---------------------------------------------------------
  subroutine td_write_output(namespace, space, gr, st, hm, ks, outp, ions, iter, dt)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(v_ks_t),             intent(inout) :: ks
    type(output_t),           intent(in)    :: outp
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iter
    FLOAT, optional,          intent(in)    :: dt
    
    character(len=256) :: filename
    type(profile_t), save :: prof

    PUSH_SUB(td_write_output)
    call profiling_in(prof, "TD_WRITE_OUTPUT")

    ! now write down the rest
    write(filename, '(a,a,i7.7)') trim(outp%iter_dir),"td.", iter  ! name of directory

    call output_all(outp, namespace, space, filename, gr, ions, iter, st, hm, ks)
    if(present(dt)) then
      call output_scalar_pot(outp, namespace, space, filename, gr, ions, hm, iter*dt)
    else
      if(iter == 0) call output_scalar_pot(outp, namespace, space, filename, gr, ions, hm)
    end if
 
    call profiling_out(prof)
    POP_SUB(td_write_output)
  end subroutine td_write_output

  ! ---------------------------------------------------------
  subroutine td_write_spin(out_spin, gr, st, iter)
    type(c_ptr),         intent(inout) :: out_spin
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: iter

    character(len=130) :: aux
    FLOAT :: spin(3)

    PUSH_SUB(td_write_spin)

    ! The expectation value of the spin operator is half the total magnetic moment
    ! This has to be calculated by all nodes
    call magnetic_moment(gr%mesh, st, st%rho, spin)
    spin = M_HALF*spin

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_spin)

        !second line -> columns name
        call write_iter_header_start(out_spin)
        if (st%d%ispin == SPINORS) then
          write(aux, '(a2,18x)') 'Sx'
          call write_iter_header(out_spin, aux)
          write(aux, '(a2,18x)') 'Sy'
          call write_iter_header(out_spin, aux)
        end if
        write(aux, '(a2,18x)') 'Sz'
        call write_iter_header(out_spin, aux)
        call write_iter_nl(out_spin)

        call td_write_print_header_end(out_spin)
      end if

      call write_iter_start(out_spin)
      select case (st%d%ispin)
      case (SPIN_POLARIZED)
         call write_iter_double(out_spin, spin(3), 1)
      case (SPINORS)
        call write_iter_double(out_spin, spin(1:3), 3)
      end select
      call write_iter_nl(out_spin)

    end if

    POP_SUB(td_write_spin)
  end subroutine td_write_spin


  ! ---------------------------------------------------------
  subroutine td_write_local_magnetic_moments(out_magnets, gr, st, ions, lmm_r, iter)
    type(c_ptr),              intent(inout) :: out_magnets
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(in)    :: st
    type(ions_t),             intent(in)    :: ions
    FLOAT,                    intent(in)    :: lmm_r
    integer,                  intent(in)    :: iter

    integer :: ia
    character(len=50) :: aux
    FLOAT, allocatable :: lmm(:,:)

    PUSH_SUB(td_write_local_magnetic_moments)

    !get the atoms` magnetization. This has to be calculated by all nodes
    SAFE_ALLOCATE(lmm(1:3, 1:ions%natoms))
    call magnetic_local_moments(gr%mesh, st, ions, gr%der%boundaries, st%rho, lmm_r, lmm)

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_magnets)

        !second line -> columns name
        call write_iter_header_start(out_magnets)
        do ia = 1, ions%natoms
          if (st%d%ispin == SPINORS) then
            write(aux, '(a2,i2.2,16x)') 'mx', ia
            call write_iter_header(out_magnets, aux)
            write(aux, '(a2,i2.2,16x)') 'my', ia
            call write_iter_header(out_magnets, aux)
          end if
          write(aux, '(a2,i2.2,16x)') 'mz', ia
          call write_iter_header(out_magnets, aux)
        end do
        call write_iter_nl(out_magnets)

        call td_write_print_header_end(out_magnets)
      end if

      call write_iter_start(out_magnets)
      do ia = 1, ions%natoms
        select case (st%d%ispin)
        case (SPIN_POLARIZED)
          call write_iter_double(out_magnets, lmm(3, ia), 1)
        case (SPINORS)
          call write_iter_double(out_magnets, lmm(1:3, ia), 3)
        end select
      end do
      call write_iter_nl(out_magnets)
      SAFE_DEALLOCATE_A(lmm)
    end if

    POP_SUB(td_write_local_magnetic_moments)
  end subroutine td_write_local_magnetic_moments

  ! ---------------------------------------------------------
  subroutine td_write_tot_mag(out_magnets, gr, st, kick, iter)
    type(c_ptr),              intent(inout) :: out_magnets
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(in)    :: st
    type(kick_t),             intent(in)    :: kick
    integer,                  intent(in)    :: iter

    CMPLX, allocatable :: tm(:,:)
    integer :: ii, iq

    PUSH_SUB(td_write_tot_mag)

    SAFE_ALLOCATE(tm(1:6,1:kick%nqvec))

    do iq = 1, kick%nqvec
      call magnetic_total_magnetization(gr%mesh, st, kick%qvector(:,iq), tm(1:6,iq))
    end do

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_magnets)
        call kick_write(kick, out = out_magnets)

        !second line -> columns name
        call write_iter_header_start(out_magnets)
        call write_iter_header(out_magnets, 'Re[m_x(q)]')
        call write_iter_header(out_magnets, 'Im[m_x(q)]')
        call write_iter_header(out_magnets, 'Re[m_y(q)]')
        call write_iter_header(out_magnets, 'Im[m_y(q)]')
        call write_iter_header(out_magnets, 'Re[m_z(q)]')
        call write_iter_header(out_magnets, 'Im[m_z(q)]')
        call write_iter_header(out_magnets, 'Re[m_x(-q)]')
        call write_iter_header(out_magnets, 'Im[m_x(-q)]')
        call write_iter_header(out_magnets, 'Re[m_y(-q)]')
        call write_iter_header(out_magnets, 'Im[m_y(-q)]')
        call write_iter_header(out_magnets, 'Re[m_z(-q)]')
        call write_iter_header(out_magnets, 'Im[m_z(-q)]')
        call write_iter_nl(out_magnets)

        call td_write_print_header_end(out_magnets)
      end if

      call write_iter_start(out_magnets)
      do iq = 1, kick%nqvec
        do ii = 1, 6
          call write_iter_double(out_magnets, TOFLOAT(tm(ii, iq)), 1)
          call write_iter_double(out_magnets, aimag(tm(ii, iq)), 1)
        end do
      end do
      call write_iter_nl(out_magnets)
    end if

    SAFE_DEALLOCATE_A(tm)

    POP_SUB(td_write_tot_mag)
  end subroutine td_write_tot_mag



  ! ---------------------------------------------------------
  subroutine td_write_angular(out_angular, namespace, gr, ions, hm, st, kick, iter)
    type(c_ptr),              intent(inout) :: out_angular
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(inout) :: ions
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(states_elec_t),      intent(inout) :: st
    type(kick_t),             intent(in)    :: kick
    integer,                  intent(in)    :: iter

    integer :: idir
    character(len=130) :: aux
    FLOAT :: angular(MAX_DIM)
    type(pert_t)        :: angular_momentum

    PUSH_SUB(td_write_angular)

    call pert_init(angular_momentum, namespace, PERTURBATION_MAGNETIC, gr, ions)

    do idir = 1, 3
       call pert_setup_dir(angular_momentum, idir)
       !we have to multiply by 2, because the perturbation returns L/2
       angular(idir) = &
         M_TWO*TOFLOAT(zpert_states_elec_expectation_value(angular_momentum, namespace, gr, ions, hm, st))
    end do

    call pert_end(angular_momentum)

    if(mpi_grp_is_root(mpi_world)) then ! Only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_angular)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        call kick_write(kick, out = out_angular)

        !second line -> columns name
        call write_iter_header_start(out_angular)
        write(aux, '(a4,18x)') '<Lx>'
        call write_iter_header(out_angular, aux)
        write(aux, '(a4,18x)') '<Ly>'
        call write_iter_header(out_angular, aux)
        write(aux, '(a4,18x)') '<Lz>'
        call write_iter_header(out_angular, aux)
        call write_iter_nl(out_angular)

        !third line -> should hold the units.
        call write_iter_string(out_angular, '#[Iter n.]')
        call write_iter_header(out_angular, '[' // trim(units_abbrev(units_out%time)) // ']')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_nl(out_angular)

        call td_write_print_header_end(out_angular)
      end if

      call write_iter_start(out_angular)
      call write_iter_double(out_angular, angular(1:3), 3)
      call write_iter_nl(out_angular)

    end if

    POP_SUB(td_write_angular)
  end subroutine td_write_angular

  ! ---------------------------------------------------------
  !> resolve state interface
  subroutine td_write_multipole(out_multip, gr, ions, st, lmax, kick, iter)
    type(td_write_prop_t), intent(inout) :: out_multip
    type(grid_t),             intent(in) :: gr   !< The grid
    type(ions_t),             intent(in) :: ions  !< Geometry object
    type(states_elec_t),      intent(in) :: st   !< State object
    integer,                  intent(in) :: lmax
    type(kick_t),             intent(in) :: kick !< Kick object
    integer,                  intent(in) :: iter !< Iteration number

    integer :: ist
    FLOAT, allocatable :: rho(:,:)

    PUSH_SUB(td_write_multipole)
    
    
    if (out_multip%resolve_states) then
      SAFE_ALLOCATE(rho(1:gr%mesh%np_part, 1:st%d%nspin))
      rho(:,:)     = M_ZERO

      do ist = st%st_start, st%st_end
        call density_calc(st, gr, rho, istin = ist)      
        call td_write_multipole_r(out_multip%mult_handles(ist), gr, ions, st, lmax, kick, rho, iter, &
                                  mpi_grp = out_multip%mpi_grp)
      end do

      SAFE_DEALLOCATE_A(rho)

    else
      
      call td_write_multipole_r(out_multip%handle, gr, ions, st, lmax, kick, st%rho, iter)

    end if
    
    POP_SUB(td_write_multipole)
  end subroutine td_write_multipole

  ! ---------------------------------------------------------
  !> Subroutine to write multipoles to the corresponding file
  subroutine td_write_multipole_r(out_multip, gr, ions, st, lmax, kick, rho, iter, mpi_grp)
    type(c_ptr),      intent(inout) :: out_multip !< C pointer
    type(grid_t),        intent(in) :: gr   !< The grid
    type(ions_t),        intent(in) :: ions  !< Geometry object
    type(states_elec_t), intent(in) :: st   !< State object
    integer,             intent(in) :: lmax
    type(kick_t),        intent(in) :: kick !< Kick object
    FLOAT,               intent(in) :: rho(:,:)
    integer,             intent(in) :: iter !< Iteration number
    type(mpi_grp_t), optional, intent(in) :: mpi_grp   
    

    integer :: is, idir, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT :: ionic_dipole(ions%space%dim)
    FLOAT, allocatable :: multipole(:,:)
    type(mpi_grp_t)    :: mpi_grp_

    PUSH_SUB(td_write_multipole_r)

    ! We cannot output multipoles beyond the dipole for higher dimensions
    ASSERT(.not. (lmax > 1 .and. gr%sb%dim > 3))

    mpi_grp_ = mpi_world 
    if (present(mpi_grp)) mpi_grp_ = mpi_grp

    if(mpi_grp_is_root(mpi_grp_).and.iter == 0) then
      call td_write_print_header_init(out_multip)

      write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,i2)')      '# lmax         ', lmax
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      call kick_write(kick, out = out_multip)

      call write_iter_header_start(out_multip)

      do is = 1, st%d%nspin
        write(aux,'(a18,i1,a1)') 'Electronic charge(', is,')'; call write_iter_header(out_multip, aux)
        if (lmax > 0) then
          do idir = 1, gr%sb%dim
            write(aux, '(4a1,i1,a1)') '<', index2axis(idir), '>', '(', is,')'
            call write_iter_header(out_multip, aux)
          end do
        end if
        do ll = 2, lmax
          do mm = -ll, ll
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', ll, ', m=', mm, ' (', is,')'
            call write_iter_header(out_multip, aux)
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      ! units
      call write_iter_string(out_multip, '#[Iter n.]')
      call write_iter_header(out_multip, '[' // trim(units_abbrev(units_out%time)) // ']')

      do is = 1, st%d%nspin
        call write_iter_header(out_multip, 'Electrons')
        if (lmax > 0) then
          do idir = 1, gr%sb%dim
            call write_iter_header(out_multip, '[' // trim(units_abbrev(units_out%length)) // ']')
          end do
        end if
        do ll = 2, lmax
          do mm = -ll, ll
            write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
            call write_iter_header(out_multip, '[' // trim(aux) // ']')
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      call td_write_print_header_end(out_multip)
    end if

    if (gr%sb%dim > 3 .and. lmax == 1) then
      ! For higher dimensions we can only have up to the dipole
      SAFE_ALLOCATE(multipole(1:gr%sb%dim+1, 1:st%d%nspin))
    else
      SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin))
    end if
    multipole(:,:) = M_ZERO

    do is = 1, st%d%nspin
      call dmf_multipoles(gr%mesh, rho(:,is), lmax, multipole(:,is))
    end do

    if (lmax > 0) then
      ionic_dipole = ions%dipole()
      do is = 1, st%d%nspin
        multipole(2:gr%sb%dim+1, is) = -ionic_dipole(1:gr%sb%dim)/st%d%nspin - multipole(2:gr%sb%dim+1, is)
      end do
    end if

    if(mpi_grp_is_root(mpi_grp_)) then
      call write_iter_start(out_multip)
      do is = 1, st%d%nspin
        call write_iter_double(out_multip, units_from_atomic(units_out%length**0, multipole(1, is)), 1)
        if (lmax > 0) then
          do idir = 1, gr%sb%dim
            call write_iter_double(out_multip, units_from_atomic(units_out%length, multipole(1+idir, is)), 1)
          end do
        end if
        add_lm = gr%sb%dim + 2
        do ll = 2, lmax
          do mm = -ll, ll
            call write_iter_double(out_multip, units_from_atomic(units_out%length**ll, multipole(add_lm, is)), 1)
            add_lm = add_lm + 1
          end do
        end do
      end do
      call write_iter_nl(out_multip)
    end if

    SAFE_DEALLOCATE_A(multipole)
    POP_SUB(td_write_multipole_r)
  end subroutine td_write_multipole_r

  ! ---------------------------------------------------------
  subroutine td_write_ftchd(out_ftchd, gr, st, kick, iter)
    type(c_ptr),         intent(inout) :: out_ftchd
    type(grid_t),        intent(in) :: gr
    type(states_elec_t), intent(in) :: st
    type(kick_t),        intent(in) :: kick
    integer,             intent(in) :: iter

    integer :: is, ip, idir
    character(len=120) :: aux, aux2
    FLOAT   :: ftchd_bessel
    CMPLX   :: ftchd
    FLOAT   :: ylm
    FLOAT, allocatable :: integrand_bessel(:)
    CMPLX, allocatable :: integrand(:)

    PUSH_SUB(td_write_ftchd)

    if(mpi_grp_is_root(mpi_world).and.iter == 0) then
      call td_write_print_header_init(out_ftchd)

      write(aux,'(a15, i2)') '# qkickmode    ', kick%qkick_mode
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        write(aux,'(a15, i0.3, 1x, i0.3)') '# ll, mm       ', kick%qbessel_l, kick%qbessel_m
        call write_iter_string(out_ftchd, aux)
        call write_iter_nl(out_ftchd)
      end if

      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        write(aux, '(a15, f9.6)') '# qlength      ', kick%qlength
      else ! sin or cos
        write(aux, '(a15)')       '# qvector      '
        do idir = 1, gr%sb%dim
          write(aux2, '(f9.5)') kick%qvector(idir,1)
          aux = trim(aux) // trim(aux2)
        end do
      end if
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      call write_iter_header_start(out_ftchd)
      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        write(aux,'(a17)') 'int(j_l*Y_lm*rho)'
      else
        write(aux,'(a12)') 'Real, Imag'
      end if
      call write_iter_header(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      ! units
      call write_iter_string(out_ftchd, '#[Iter n.]')
      call write_iter_header(out_ftchd, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_nl(out_ftchd)
      call td_write_print_header_end(out_ftchd)

    end if

    ftchd = M_ZERO

    ! If kick mode is exp, sin, or cos, apply the normal Fourier transform
    if(kick%qkick_mode /= QKICKMODE_BESSEL) then
      SAFE_ALLOCATE(integrand(1:gr%mesh%np))
      integrand = M_ZERO
      do is = 1, st%d%nspin
        do ip = 1, gr%mesh%np
          integrand(ip) = integrand(ip) + st%rho(ip, is) * exp(-M_zI*sum(gr%mesh%x(ip, 1:gr%sb%dim)*kick%qvector(1:gr%sb%dim, 1)))
        end do
      end do
      ftchd = zmf_integrate(gr%mesh, integrand)
      SAFE_DEALLOCATE_A(integrand)
    else
      ftchd_bessel = M_ZERO
      SAFE_ALLOCATE(integrand_bessel(1:gr%mesh%np))
      integrand_bessel = M_ZERO
      do is = 1, st%d%nspin
        do ip = 1, gr%mesh%np
          call grylmr(gr%mesh%x(ip, 1), gr%mesh%x(ip, 2), gr%mesh%x(ip, 3), kick%qbessel_l, kick%qbessel_m, ylm)
          integrand_bessel(ip) = integrand_bessel(ip) + st%rho(ip, is) * &
                                 loct_sph_bessel(kick%qbessel_l, kick%qlength*sqrt(sum(gr%mesh%x(ip, :)**2)))*ylm
        end do
      end do
      ftchd_bessel = dmf_integrate(gr%mesh, integrand_bessel)
      SAFE_DEALLOCATE_A(integrand_bessel)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_ftchd)
      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        call write_iter_double(out_ftchd, ftchd_bessel, 1)
      else ! exp, sin, cos
        call write_iter_double(out_ftchd, real(ftchd), 1)
        call write_iter_double(out_ftchd, aimag(ftchd), 1)
      end if
      call write_iter_nl(out_ftchd)
    end if

    POP_SUB(td_write_ftchd)
  end subroutine td_write_ftchd

  ! ---------------------------------------------------------
  subroutine td_write_coordinates(out_coords, ions, iter)
    type(c_ptr),       intent(inout) :: out_coords
    type(ions_t),      intent(in)    :: ions
    integer,           intent(in)    :: iter

    integer :: iatom, idir
    character(len=50) :: aux
    FLOAT :: tmp(1:MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_coordinates)

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do iatom = 1, ions%natoms
        do idir = 1, ions%space%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', iatom, ',', idir, ')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do iatom = 1, ions%natoms
        do idir = 1, ions%space%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', iatom, ',', idir,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do iatom = 1, ions%natoms
        do idir = 1, ions%space%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', iatom, ',', idir,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(out_coords, &
        'Positions in '   // trim(units_abbrev(units_out%length))   //   &
        ', Velocities in '// trim(units_abbrev(units_out%velocity)) //   &
        ', Forces in '    // trim(units_abbrev(units_out%force)))
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    do iatom = 1, ions%natoms
      tmp(1:ions%space%dim) = units_from_atomic(units_out%length, ions%pos(:, iatom))
      call write_iter_double(out_coords, tmp, ions%space%dim)
    end do
    do iatom = 1, ions%natoms
      tmp(1:ions%space%dim) = units_from_atomic(units_out%velocity, ions%vel(:, iatom))
      call write_iter_double(out_coords, tmp, ions%space%dim)
    end do
    do iatom = 1, ions%natoms
      tmp(1:ions%space%dim) = units_from_atomic(units_out%force, ions%tot_force(:, iatom))
      call write_iter_double(out_coords, tmp, ions%space%dim)
    end do
    call write_iter_nl(out_coords)

    POP_SUB(td_write_coordinates)
  end subroutine td_write_coordinates

  ! ---------------------------------------------------------
  subroutine td_write_sep_coordinates(out_coords, ions, iter, which)
    type(c_ptr),       intent(inout) :: out_coords
    type(ions_t),      intent(in)    :: ions
    integer,           intent(in)    :: iter
    integer,           intent(in)    :: which !1=xyz, 2=velocity, 3=force

    integer, parameter :: COORDINATES=1
    integer, parameter :: VELOCITIES=2
    integer, parameter :: FORCES=3
    integer :: iatom, idir
    character(len=50) :: aux
    FLOAT :: tmp(1:MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_sep_coordinates)

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do iatom = 1, ions%natoms
        do idir = 1, ions%space%dim
          select case (which)
            case (COORDINATES)
              write(aux, '(a2,i3,a1,i3,a1)') 'x(', iatom, ',', idir, ')'
            case (VELOCITIES)
              write(aux, '(a2,i3,a1,i3,a1)') 'v(', iatom, ',', idir,')'
            case (FORCES)
              write(aux, '(a2,i3,a1,i3,a1)') 'f(', iatom, ',', idir,')'
          end select
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_abbrev(units_out%time)) // ']')
      select case (which)
        case (COORDINATES)
          call write_iter_string(out_coords, &
            'Positions in '   // trim(units_abbrev(units_out%length))) 
        case (VELOCITIES)
          call write_iter_string(out_coords, &
            'Velocities in '  // trim(units_abbrev(units_out%velocity)))
        case (FORCES)
          call write_iter_string(out_coords, &
            'Forces in '    // trim(units_abbrev(units_out%force)))
      end select
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    select case (which)
      case (COORDINATES)
        do iatom = 1, ions%natoms
          tmp(1:ions%space%dim) = units_from_atomic(units_out%length, ions%pos(:, iatom))
          call write_iter_double(out_coords, tmp, ions%space%dim)
        end do
      case (VELOCITIES)
        do iatom = 1, ions%natoms
          tmp(1:ions%space%dim) = units_from_atomic(units_out%velocity, ions%vel(:, iatom))
          call write_iter_double(out_coords, tmp, ions%space%dim)
        end do
      case (FORCES)
        do iatom = 1, ions%natoms
          tmp(1:ions%space%dim) = units_from_atomic(units_out%force, ions%tot_force(:, iatom))
          call write_iter_double(out_coords, tmp, ions%space%dim)
        end do
    end select
       
    call write_iter_nl(out_coords)

    POP_SUB(td_write_sep_coordinates)
  end subroutine td_write_sep_coordinates


  ! ---------------------------------------------------------
  subroutine td_write_temperature(out_temperature, ions, iter)
    type(c_ptr),       intent(inout) :: out_temperature
    type(ions_t),      intent(in) :: ions
    integer,           intent(in) :: iter

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_temperature)

    if(iter == 0) then
      call td_write_print_header_init(out_temperature)

      ! first line: column names
      call write_iter_header_start(out_temperature)
      call write_iter_header(out_temperature, 'Temperature')
      call write_iter_nl(out_temperature)

      ! second line: units
      call write_iter_string(out_temperature, '#[Iter n.]')
      call write_iter_header(out_temperature, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(out_temperature, '        [K]')
      call write_iter_nl(out_temperature)

      call td_write_print_header_end(out_temperature)
    end if

    call write_iter_start(out_temperature)

    call write_iter_double(out_temperature, units_from_atomic(unit_kelvin, ion_dynamics_temperature(ions)), 1)

    call write_iter_nl(out_temperature)

    POP_SUB(td_write_temperature)
  end subroutine td_write_temperature


  ! ---------------------------------------------------------
  subroutine td_write_populations(out_populations, namespace, space, mesh, st, writ, dt, iter)
    type(c_ptr),            intent(inout) :: out_populations
    type(namespace_t),      intent(in)    :: namespace
    type(space_t),          intent(in)    :: space
    type(mesh_t),           intent(in)    :: mesh
    type(states_elec_t),    intent(inout) :: st
    type(td_write_t),       intent(in)    :: writ
    FLOAT,                  intent(in)    :: dt
    integer,                intent(in)    :: iter
 
    integer :: ist
    character(len=6) :: excited_name
    CMPLX :: gsp
    CMPLX, allocatable :: excited_state_p(:)
    CMPLX, allocatable :: dotprodmatrix(:, :, :)


    PUSH_SUB(td_write_populations)

    SAFE_ALLOCATE(dotprodmatrix(1:writ%gs_st%nst, 1:st%nst, 1:st%d%nik))
    call zstates_elec_matrix(writ%gs_st, st, mesh, dotprodmatrix)


    !See comment in zstates_elec_mpdotp
    ASSERT(.not. space%is_periodic())

    ! all processors calculate the projection
    gsp = zstates_elec_mpdotp(namespace, mesh, writ%gs_st, st, dotprodmatrix)

    if(writ%n_excited_states > 0) then
      SAFE_ALLOCATE(excited_state_p(1:writ%n_excited_states))
      do ist = 1, writ%n_excited_states
        excited_state_p(ist) = zstates_elec_mpdotp(namespace, mesh, writ%excited_st(ist), st, dotprodmatrix)
      end do
    end if

    if(mpi_grp_is_root(mpi_world)) then
      if(iter == 0) then
        call td_write_print_header_init(out_populations)

        ! first line -> column names
        call write_iter_header_start(out_populations)
        call write_iter_header(out_populations, 'Re<Phi_gs|Phi(t)>')
        call write_iter_header(out_populations, 'Im<Phi_gs|Phi(t)>')
        do ist = 1, writ%n_excited_states
          write(excited_name,'(a2,i3,a1)') 'P(', ist,')'
          call write_iter_header(out_populations, 'Re<'//excited_name//'|Phi(t)>')
          call write_iter_header(out_populations, 'Im<'//excited_name//'|Phi(t)>')
        end do
        call write_iter_nl(out_populations)

        ! second line -> units
        call write_iter_string(out_populations, '#[Iter n.]')
        call write_iter_header(out_populations, '[' // trim(units_abbrev(units_out%time)) // ']')
        call write_iter_nl(out_populations)

        call td_write_print_header_end(out_populations)
      end if

      ! cannot call write_iter_start, for the step is not 1
      call write_iter_int(out_populations, iter, 1)
      call write_iter_double(out_populations, units_from_atomic(units_out%time, iter*dt),  1)
      call write_iter_double(out_populations, real(gsp),  1)
      call write_iter_double(out_populations, aimag(gsp), 1)
      do ist = 1, writ%n_excited_states
        call write_iter_double(out_populations, real(excited_state_p(ist)),  1)
        call write_iter_double(out_populations, aimag(excited_state_p(ist)), 1)
      end do
      call write_iter_nl(out_populations)
    end if

    if(writ%n_excited_states > 0) then
      SAFE_DEALLOCATE_A(excited_state_p)
    end if
    SAFE_DEALLOCATE_A(dotprodmatrix)
    POP_SUB(td_write_populations)
  end subroutine td_write_populations


  ! ---------------------------------------------------------
  subroutine td_write_acc(out_acc, namespace, space, gr, ions, st, hm, dt, iter)
    type(c_ptr),              intent(inout) :: out_acc
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                    intent(in)    :: dt
    integer,                  intent(in)    :: iter

    integer :: idim
    character(len=7) :: aux
    FLOAT :: acc(MAX_DIM)

    PUSH_SUB(td_write_acc)

    if(iter == 0 .and. mpi_grp_is_root(mpi_world)) then
      call td_write_print_header_init(out_acc)

      ! first line -> column names
      call write_iter_header_start(out_acc)
      do idim = 1, space%dim
        write(aux, '(a4,i1,a1)') 'Acc(', idim, ')'
        call write_iter_header(out_acc, aux)
      end do
      call write_iter_nl(out_acc)

      ! second line: units
      call write_iter_string(out_acc, '#[Iter n.]')
      call write_iter_header(out_acc, '[' // trim(units_abbrev(units_out%time)) // ']')
      do idim = 1, space%dim
        call write_iter_header(out_acc, '[' // trim(units_abbrev(units_out%acceleration)) // ']')
      end do
      call write_iter_nl(out_acc)
      call td_write_print_header_end(out_acc)
    end if

    call td_calc_tacc(namespace, space, gr, ions, st, hm, acc, dt*iter)

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_acc)
      acc = units_from_atomic(units_out%acceleration, acc)
      call write_iter_double(out_acc, acc, space%dim)
      call write_iter_nl(out_acc)
    end if

    POP_SUB(td_write_acc)
  end subroutine td_write_acc
  
  ! ---------------------------------------------------------
  subroutine td_write_vel(out_vel, space, der, st, kpoints, iter)
    type(c_ptr),         intent(inout) :: out_vel
    type(space_t),       intent(in)    :: space
    type(derivatives_t), intent(in)    :: der
    type(states_elec_t), intent(inout) :: st
    type(kpoints_t),     intent(in)    :: kpoints
    integer,             intent(in)    :: iter

    integer :: idim
    character(len=7) :: aux
    FLOAT :: vel(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_vel)

    if(iter == 0) then
      call td_write_print_header_init(out_vel)

      ! first line -> column names
      call write_iter_header_start(out_vel)
      do idim = 1, space%dim
        write(aux, '(a4,i1,a1)') 'Vel(', idim, ')'
        call write_iter_header(out_vel, aux)
      end do
      call write_iter_nl(out_vel)

      ! second line: units
      call write_iter_string(out_vel, '#[Iter n.]')
      call write_iter_header(out_vel, '[' // trim(units_abbrev(units_out%time)) // ']')
      do idim = 1, space%dim
        call write_iter_header(out_vel, '[' // trim(units_abbrev(units_out%velocity)) // ']')
      end do
      call write_iter_nl(out_vel)
      call td_write_print_header_end(out_vel)
    end if

    call td_calc_tvel(space, der, st, kpoints, vel)

    call write_iter_start(out_vel)
    vel = units_from_atomic(units_out%velocity, vel)
    call write_iter_double(out_vel, vel, space%dim)
    call write_iter_nl(out_vel)

    POP_SUB(td_write_vel)
  end subroutine td_write_vel


  ! ---------------------------------------------------------
  subroutine td_write_laser(out_laser, gr, hm, dt, iter)
    type(c_ptr),         intent(inout) :: out_laser
    type(grid_t),        intent(in) :: gr
    type(hamiltonian_elec_t), intent(in) :: hm
    FLOAT,               intent(in) :: dt
    integer,             intent(in) :: iter

    integer :: il, idir
    FLOAT :: field(MAX_DIM)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    ! no PUSH SUB, called too often

    if(iter == 0) then
      call td_write_print_header_init(out_laser)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, dt), &
        " [", trim(units_abbrev(units_out%time)), "]"
      call write_iter_string(out_laser, aux)
      call write_iter_nl(out_laser)

      call write_iter_header_start(out_laser)
      do il = 1, hm%ext_lasers%no_lasers
        select case(laser_kind(hm%ext_lasers%lasers(il)))
        case(E_FIELD_ELECTRIC)
          do idir = 1, gr%sb%dim
            write(aux, '(a,i1,a)') 'E(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_MAGNETIC)
          do idir = 1, gr%sb%dim
            write(aux, '(a,i1,a)') 'B(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          do idir = 1, gr%sb%dim
            write(aux, '(a,i1,a)') 'A(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_SCALAR_POTENTIAL)
          write(aux, '(a,i1,a)') 'e(t)'
          call write_iter_header(out_laser, aux)
        end select
      end do
      call write_iter_nl(out_laser)

      call write_iter_string(out_laser, '#[Iter n.]')
      call write_iter_header(out_laser, '[' // trim(units_abbrev(units_out%time)) // ']')

      ! Note that we do not print out units of E, B, or A, but rather units of e*E, e*B, e*A.
      ! (force, force, and energy, respectively). The reason is that the units of E, B or A 
      ! are ugly.
      do il = 1, hm%ext_lasers%no_lasers
        select case(laser_kind(hm%ext_lasers%lasers(il)))
        case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
          aux = '[' // trim(units_abbrev(units_out%force)) // ']'
          do idir = 1, gr%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          aux = '[' // trim(units_abbrev(units_out%energy)) // ']'
          do idir = 1, gr%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_SCALAR_POTENTIAL)
          aux = '[adim]'
          call write_iter_header(out_laser, aux)
        end select
      end do
      call write_iter_nl(out_laser)

      call td_write_print_header_end(out_laser)
    end if

    call write_iter_start(out_laser)

    do il = 1, hm%ext_lasers%no_lasers
      field = M_ZERO
      call laser_field(hm%ext_lasers%lasers(il), field(1:gr%sb%dim), iter*dt)
      select case(laser_kind(hm%ext_lasers%lasers(il)))
      case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
        field = units_from_atomic(units_out%force, field)
        call write_iter_double(out_laser, field, gr%sb%dim)
      case(E_FIELD_VECTOR_POTENTIAL)
        field = units_from_atomic(units_out%energy, field)
        call write_iter_double(out_laser, field, gr%sb%dim)
      case(E_FIELD_SCALAR_POTENTIAL)
        call write_iter_double(out_laser, field(1), 1)
      end select
    end do

    call write_iter_nl(out_laser)

  end subroutine td_write_laser


  ! ---------------------------------------------------------
  subroutine td_write_energy(out_energy, hm, iter, ke)
    type(c_ptr),         intent(inout)   :: out_energy
    type(hamiltonian_elec_t), intent(in) :: hm
    integer,             intent(in)      :: iter
    FLOAT, intent(in)                    :: ke

    integer :: ii

    integer :: n_columns

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_energy)

    n_columns = 9          

    if(iter == 0) then
      call td_write_print_header_init(out_energy)

      ! first line -> column names
      call write_iter_header_start(out_energy)
      call write_iter_header(out_energy, 'Total')
      call write_iter_header(out_energy, 'Kinetic (ions)')
      call write_iter_header(out_energy, 'Ion-Ion')
      call write_iter_header(out_energy, 'Electronic')
      call write_iter_header(out_energy, 'Eigenvalues')
      call write_iter_header(out_energy, 'Hartree')
      call write_iter_header(out_energy, 'Int[n v_xc]')
      call write_iter_header(out_energy, 'Exchange')
      call write_iter_header(out_energy, 'Correlation')

      if (hm%pcm%run_pcm) then 
          call write_iter_header(out_energy, 'E_M-solvent')
          n_columns = n_columns + 1    
      end if   

      if (hm%lda_u_level /= DFT_U_NONE) then
          call write_iter_header(out_energy, 'Hubbard')
          n_columns = n_columns + 1
      end if

      call write_iter_nl(out_energy)

      ! units

      call write_iter_string(out_energy, '#[Iter n.]')
      call write_iter_header(out_energy, '[' // trim(units_abbrev(units_out%time)) // ']')

      do ii = 1, n_columns
        call write_iter_header(out_energy, '[' // trim(units_abbrev(units_out%energy)) // ']')
      end do
      call write_iter_nl(out_energy)
      
      
      call td_write_print_header_end(out_energy)
    end if

    call write_iter_start(out_energy)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%total+ke), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, ke), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%ep%eii), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%total-hm%ep%eii), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%eigenvalues), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%hartree), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%intnvxc), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%exchange), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%correlation), 1)
    
    !adding the molecule-solvent electrostatic interaction
    if (hm%pcm%run_pcm) call write_iter_double(out_energy, &
                             units_from_atomic(units_out%energy, hm%energy%int_ee_pcm + hm%energy%int_en_pcm + &
                                                                 hm%energy%int_nn_pcm + hm%energy%int_ne_pcm), 1)

    if(hm%lda_u_level /= DFT_U_NONE) &
         call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%dft_u), 1)

    call write_iter_nl(out_energy)

    POP_SUB(td_write_energy)
  end subroutine td_write_energy

  ! ---------------------------------------------------------
  subroutine td_write_eigs(out_eigs, st, iter)
    type(c_ptr),         intent(inout) :: out_eigs
    type(states_elec_t),      intent(in) :: st
    integer,             intent(in) :: iter

    integer             :: ii, is
    character(len=68)   :: buf

    PUSH_SUB(td_write_eigs)
  
    if(.not.mpi_grp_is_root(mpi_world)) then 
      POP_SUB(td_write_eigs)
      return ! only first node outputs        
    end if


    if(iter == 0) then
      call td_write_print_header_init(out_eigs)

      write(buf, '(a15,i2)')      '# nst          ', st%nst
      call write_iter_string(out_eigs, buf)
      call write_iter_nl(out_eigs)

      write(buf, '(a15,i2)')      '# nspin        ', st%d%nspin
      call write_iter_string(out_eigs, buf)
      call write_iter_nl(out_eigs)

      ! first line -> column names
      call write_iter_header_start(out_eigs)
      do is = 1, st%d%kpt%nglobal
        do ii = 1, st%nst
          write(buf, '(a,i4)') 'Eigenvalue ',ii
          call write_iter_header(out_eigs, buf)
        end do
      end do
      call write_iter_nl(out_eigs)

      ! second line: units
      call write_iter_string(out_eigs, '#[Iter n.]')
      call write_iter_header(out_eigs, '[' // trim(units_abbrev(units_out%time)) // ']')
      do is = 1, st%d%kpt%nglobal
        do ii = 1, st%nst
          call write_iter_header(out_eigs, '[' // trim(units_abbrev(units_out%energy)) // ']')
        end do
      end do
      call write_iter_nl(out_eigs)
      call td_write_print_header_end(out_eigs)
    end if

    call write_iter_start(out_eigs)
    do is = 1, st%d%kpt%nglobal
      do ii =1 , st%nst
        call write_iter_double(out_eigs, units_from_atomic(units_out%energy, st%eigenval(ii,is)), 1)
      end do
    end do
    call write_iter_nl(out_eigs)

    POP_SUB(td_write_eigs)
  end subroutine td_write_eigs

  ! ---------------------------------------------------------
  subroutine td_write_ionch(out_ionch, gr, st, iter)
    type(c_ptr),         intent(inout) :: out_ionch
    type(grid_t),        intent(in) :: gr
    type(states_elec_t),      intent(in) :: st
    integer,             intent(in) :: iter

    integer             :: ii, ist, Nch, ik, idim
    character(len=68)   :: buf
    FLOAT, allocatable  :: ch(:), occ(:)

#if defined(HAVE_MPI) 
    FLOAT, allocatable :: occbuf(:)
#endif

    PUSH_SUB(td_write_ionch)
    

    Nch =  st%nst * st%d%kpt%nglobal * st%d%dim 
    SAFE_ALLOCATE(ch(0: Nch)) 
    SAFE_ALLOCATE(occ(0: Nch))

    occ(:) = M_ZERO
    ii = 1
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim        
          if (st%st_start <= ist .and. ist <= st%st_end .and. &
              st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
            occ(ii) = st%occ(ist, ik)
          end if
          ii = ii+1
        end do
      end do
    end do
      
      
#if defined(HAVE_MPI) 
    if(st%parallel_in_states) then
      SAFE_ALLOCATE(occbuf(0: Nch)) 
      occbuf(:) = M_ZERO
      call MPI_Allreduce(occ(0), occbuf(0), Nch+1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      occ(:) = occbuf(:)
      SAFE_DEALLOCATE_A(occbuf) 
    end if
#endif      
    
    !Calculate the channels
    call td_calc_ionch(gr, st, ch, Nch)
  
  
    if(.not.mpi_grp_is_root(mpi_world)) then 
      SAFE_DEALLOCATE_A(ch)
      POP_SUB(td_write_ionch)
      return ! only first node outputs        
    end if


    if(iter == 0) then
      call td_write_print_header_init(out_ionch)
      
      ! first line -> column names
      call write_iter_header_start(out_ionch)
      
      do ii = 0, Nch
        if(occ(ii)>M_ZERO .or. ii == 0) then
          write(buf, '(a,f4.1,a)') 'Pion(',occ(ii)*ii,'+, t)'
          call write_iter_header(out_ionch, buf)
        end if
      end do
      call write_iter_nl(out_ionch)

      ! second line: units
      call write_iter_string(out_ionch, '#[Iter n.]')
      call write_iter_header(out_ionch, '[' // trim(units_abbrev(units_out%time)) // ']')
      do ii = 0, Nch
        if(occ(ii)>M_ZERO .or. ii == 0) then
          call write_iter_header(out_ionch, '[' // trim(units_abbrev(unit_one)) // ']')
        end if
      end do
      call write_iter_nl(out_ionch)
      call td_write_print_header_end(out_ionch)
    end if

    call write_iter_start(out_ionch)
    do ii =0 , Nch
      if(occ(ii)>M_ZERO .or. ii == 0) then
        call write_iter_double(out_ionch, units_from_atomic(unit_one, ch(ii)), 1)
      end if
    end do
    call write_iter_nl(out_ionch)

    SAFE_DEALLOCATE_A(ch)
    SAFE_DEALLOCATE_A(occ)

    POP_SUB(td_write_ionch)
  end subroutine td_write_ionch

  ! ---------------------------------------------------------
  subroutine td_write_proj(out_proj, space, gr, ions, st, gs_st, kick, iter)
    type(c_ptr),         intent(inout) :: out_proj
    type(space_t),       intent(in)    :: space
    type(grid_t),        intent(in)    :: gr
    type(ions_t),        intent(in)    :: ions
    type(states_elec_t), intent(inout) :: st
    type(states_elec_t), intent(in)    :: gs_st
    type(kick_t),        intent(in)    :: kick
    integer,             intent(in)    :: iter

    CMPLX, allocatable :: projections(:,:,:)
    character(len=80) :: aux
    integer :: ik, ist, uist, idir

    PUSH_SUB(td_write_proj)

    if(iter == 0) then
      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_init(out_proj)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        call kick_write(kick, out = out_proj)

        call write_iter_string(out_proj, "#%")
        call write_iter_nl(out_proj)
 
        write(aux, '(a,i8)') "# nik  ", st%d%nik
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        write(aux, '(a,2i8)') "#  st  ", gs_st%st_start, st%nst
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        write(aux, '(a,2i8)') "# ust  ", gs_st%st_start, gs_st%st_end
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        do ik = 1, st%d%nik
          call write_iter_string(out_proj, "# w(ik)*occ(ist,ik)  ")
          do ist = gs_st%st_start, st%nst
            call write_iter_double(out_proj, st%d%kweights(ik)*st%occ(ist, ik), 1)
          end do
          call write_iter_nl(out_proj)
        end do

        call write_iter_header_start(out_proj)
        do ik = 1, st%d%nik
          do ist = gs_st%st_start, st%nst
            do uist = gs_st%st_start, gs_st%st_end
              write(aux, '(i4,a,i4)') ist, ' -> ', uist
              call write_iter_header(out_proj, 'Re {'//trim(aux)//'}')
              call write_iter_header(out_proj, 'Im {'//trim(aux)//'}')
            end do
          end do
        end do
        call write_iter_nl(out_proj)

      end if

      !The dipole matrix elements cannot be computed like that for solids
      if(.not. space%is_periodic()) then

        SAFE_ALLOCATE(projections(1:st%nst, gs_st%st_start:gs_st%st_end, 1:st%d%nik))
        do idir = 1, ions%space%dim
          projections = M_Z0

          call dipole_matrix_elements(idir)

          if(mpi_grp_is_root(mpi_world)) then
            write(aux, '(a,i1,a)') "<i|x_", idir, "|a>"
            call write_iter_string(out_proj, "# ------")
            call write_iter_header(out_proj, aux)
            do ik = 1, st%d%nik
              do ist = gs_st%st_start, st%st_end
                do uist = gs_st%st_start, gs_st%st_end
                  call write_iter_double(out_proj,  real(projections(ist, uist, ik)), 1)
                  call write_iter_double(out_proj, aimag(projections(ist, uist, ik)), 1)
                end do
              end do
            end do
            call write_iter_nl(out_proj)
          
          end if
        end do
        SAFE_DEALLOCATE_A(projections)

      end if

      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_end(out_proj)
      end if

    end if

    SAFE_ALLOCATE(projections(1:st%nst, gs_st%st_start:gs_st%st_end, 1:st%d%nik))
    projections(:,:,:) = M_Z0
    call calc_projections(gr%mesh, st, gs_st, projections)

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_proj)
      do ik = 1, st%d%nik
        do ist = gs_st%st_start, st%nst
          do uist = gs_st%st_start, gs_st%st_end
            call write_iter_double(out_proj, TOFLOAT(projections(ist, uist, ik)), 1)
            call write_iter_double(out_proj, aimag(projections(ist, uist, ik)), 1)
          end do
        end do
      end do
      call write_iter_nl(out_proj)
    end if

    SAFE_DEALLOCATE_A(projections)
    POP_SUB(td_write_proj)

  contains
    ! ---------------------------------------------------------
    subroutine dipole_matrix_elements(dir)
      integer, intent(in) :: dir

      integer :: uist, ist, ik, idim
      FLOAT   :: n_dip(space%dim)
      CMPLX, allocatable :: xpsi(:,:)
      CMPLX, allocatable :: psi(:, :), gspsi(:, :)
      
      PUSH_SUB(td_write_proj.dipole_matrix_elements)
 
      SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(gspsi(1:gr%mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(xpsi(1:gr%mesh%np, 1:st%d%dim))

      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          call states_elec_get_state(st, gr%mesh, ist, ik, psi)
          do uist = gs_st%st_start, gs_st%st_end
            call states_elec_get_state(gs_st, gr%mesh, uist, ik, gspsi)

            do idim = 1, st%d%dim
              xpsi(1:gr%mesh%np, idim) = gr%mesh%x(1:gr%mesh%np, dir)*gspsi(1:gr%mesh%np, idim)
            end do
            projections(ist, uist, ik) = -zmf_dotp(gr%mesh, st%d%dim, psi, xpsi, reduce = .false.)

          end do
        end do
      end do
      
      SAFE_DEALLOCATE_A(xpsi)
      SAFE_DEALLOCATE_A(gspsi)
      SAFE_DEALLOCATE_A(psi)

      call comm_allreduce(st%dom_st_kpt_mpi_grp, projections)

      ! n_dip is not defined for more than space%dim
      n_dip = ions%dipole()
      do ik = 1, st%d%nik
        do ist = gs_st%st_start, st%nst
          do uist = gs_st%st_start, gs_st%st_end
            projections(ist, uist, ik) = projections(ist, uist, ik) - n_dip(dir)
          end do
        end do
      end do


      POP_SUB(td_write_proj.dipole_matrix_elements)
    end subroutine dipole_matrix_elements

  end subroutine td_write_proj

  ! ---------------------------------------------------------
  !> This routine computes the total number of excited electrons
  !> based on projections on the GS orbitals
  !> The procedure is very similar to the td_write_proj
  ! ---------------------------------------------------------
  subroutine td_write_n_ex(out_nex, outp, namespace, gr, kpoints, st, gs_st, iter)
    type(c_ptr),         intent(inout) :: out_nex
    type(output_t),      intent(in)    :: outp
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(kpoints_t),     intent(in)    :: kpoints
    type(states_elec_t), intent(inout) :: st
    type(states_elec_t), intent(in)    :: gs_st
    integer,             intent(in)    :: iter

    CMPLX, allocatable :: projections(:,:)
    character(len=80) :: aux, dir
    integer :: ik, ikpt, ist, uist, err
    FLOAT :: Nex, weight
    integer :: gs_nst
    FLOAT, allocatable :: Nex_kpt(:)
    

    PUSH_SUB(td_write_n_ex)

    if(iter == 0) then
      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_init(out_nex)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_nex, aux)
        call write_iter_nl(out_nex)

        call write_iter_string(out_nex, "#%")
        call write_iter_nl(out_nex)

        write(aux, '(a,i8)') "# nik  ", st%d%nik
        call write_iter_string(out_nex, aux)
        call write_iter_nl(out_nex)

        write(aux, '(a,2i8)') "#  st  ", gs_st%st_start, st%nst
        call write_iter_string(out_nex, aux)
        call write_iter_nl(out_nex)

        write(aux, '(a,2i8)') "# ust  ", gs_st%st_start, gs_st%st_end
        call write_iter_string(out_nex, aux)
        call write_iter_nl(out_nex)

        call write_iter_header_start(out_nex)
        call write_iter_header(out_nex, '#  iter t Nex(t)')
        call write_iter_nl(out_nex)

      end if

      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_end(out_nex)
      end if

    end if

    !We only need the occupied GS states
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = 1, gs_st%nst
        if(st%occ(ist, ik)>M_EPSILON) gs_nst = ist
      end do
    end do

    SAFE_ALLOCATE(projections(1:gs_nst, 1:st%nst))
     
    SAFE_ALLOCATE(Nex_kpt(1:st%d%nik)) 
    Nex_kpt = M_ZERO 
    do ik = st%d%kpt%start, st%d%kpt%end
      ikpt = st%d%get_kpoint_index(ik)
      call zstates_elec_calc_projections(st, gs_st, namespace, gr%mesh, ik, projections, gs_nst)
      do ist = 1, gs_nst
        weight = st%d%kweights(ik) * st%occ(ist, ik)/ st%smear%el_per_state 
        do uist = st%st_start, st%st_end
          Nex_kpt(ikpt) = Nex_kpt(ikpt) - weight * st%occ(uist, ik) * abs(projections(ist, uist))**2
        end do
      end do
     if(st%d%ispin == SPIN_POLARIZED) then
       Nex_kpt(ikpt) = Nex_kpt(ikpt) + st%qtot*M_HALF*st%d%kweights(ik)
     else
       Nex_kpt(ikpt) = Nex_kpt(ikpt) + st%qtot*st%d%kweights(ik)
     end if
    end do

   if(st%parallel_in_states .or. st%d%kpt%parallel) then
     call comm_allreduce(st%st_kpt_mpi_grp, Nex_kpt)
   end if

  Nex = sum(Nex_kpt)

  if(mpi_grp_is_root(mpi_world)) then
    call write_iter_start(out_nex)
    call write_iter_double(out_nex, Nex, 1)
    call write_iter_nl(out_nex)
  end if
 
  ! now write down the k-resolved part
  write(dir, '(a,a,i7.7)') trim(outp%iter_dir),"td.", iter  ! name of directory

  call io_function_output_global_BZ(outp%how(OPTION__OUTPUT__CURRENT_KPT) &
    + outp%how(OPTION__OUTPUT__DENSITY_KPT), dir, "n_excited_el_kpt", namespace, &
    kpoints, Nex_kpt, unit_one, err) 

 
  SAFE_DEALLOCATE_A(projections)
  SAFE_DEALLOCATE_A(Nex_kpt)

  POP_SUB(td_write_n_ex)
 end subroutine td_write_n_ex

   ! ---------------------------------------------------------
   !> This subroutine calculates:
   !! \f[
   !! p(uist, ist, ik) = < \phi(ist, ik) (t) | \phi_0(uist, k) >
   !! \f]
   ! ---------------------------------------------------------
   subroutine calc_projections(mesh, st, gs_st, projections)
     implicit none 
    
     type(mesh_t),        intent(in)    :: mesh
     type(states_elec_t), intent(inout) :: st
     type(states_elec_t), intent(in)    :: gs_st
     CMPLX,               intent(inout) :: projections(1:st%nst, &
                                         gs_st%st_start:gs_st%nst, 1:st%d%nik)
 
     integer :: uist, ist, ik
     CMPLX, allocatable :: psi(:, :), gspsi(:, :)
     PUSH_SUB(calc_projections)
    
     SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
     SAFE_ALLOCATE(gspsi(1:mesh%np, 1:st%d%dim))

     projections(:,:,:) = M_ZERO
     
     do ik = st%d%kpt%start, st%d%kpt%end 
       do ist = st%st_start, st%st_end
         call states_elec_get_state(st, mesh, ist, ik, psi)
         do uist = gs_st%st_start, gs_st%nst
           call states_elec_get_state(gs_st, mesh, uist, ik, gspsi)
           projections(ist, uist, ik) = zmf_dotp(mesh, st%d%dim, psi, gspsi, reduce = .false.)
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(gspsi)

    call comm_allreduce(st%dom_st_kpt_mpi_grp, projections)


    POP_SUB(calc_projections)
  end subroutine calc_projections


  subroutine td_write_proj_kp(gr, kpoints, st, gs_st, namespace, iter)
    type(grid_t),        intent(in)    :: gr
    type(kpoints_t),     intent(in)    :: kpoints
    type(states_elec_t), intent(in)    :: st
    type(states_elec_t), intent(inout) :: gs_st
    type(namespace_t),   intent(in)    :: namespace
    integer,             intent(in)    :: iter

    CMPLX, allocatable :: proj(:,:), psi(:,:,:), gs_psi(:,:,:), temp_state(:,:)
    character(len=80) :: filename1, filename2
    integer :: ik,ist, jst, file, idim, nk_proj
    type(mesh_t) :: mesh

    PUSH_SUB(td_write_proj_kp)

    ! this is slow, so we don`t do it every step
    if(.not.mod(iter,50) == 0) then
       POP_SUB(td_write_proj_kp)
       return
    end if

    mesh = gr%mesh

    write(filename1,'(I10)') iter
    filename1 = 'td.general/projections_iter_'//trim(adjustl(filename1))
    file = 9845623
    
    SAFE_ALLOCATE(proj(1:gs_st%nst, 1:gs_st%nst))
    SAFE_ALLOCATE(psi(1:gs_st%nst,1:gs_st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(gs_psi(1:gs_st%nst,1:gs_st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(temp_state(1:mesh%np,1:gs_st%d%dim))
    
    ! Project only k-points that have a zero weight.
    ! Why? It is unlikely that one is interested in the projections 
    ! of the Monkhorst-Pack kpoints, but instead we assume that
    ! the user has specified a k-path with zero weights
    nk_proj = kpoints%nik_skip

    do ik = kpoints%reduced%npoints-nk_proj+1, kpoints%reduced%npoints
      ! reset arrays
      psi(1:gs_st%nst, 1:gs_st%d%dim, 1:mesh%np)= M_ZERO
      gs_psi(1:gs_st%nst, 1:gs_st%d%dim, 1:mesh%np)= M_ZERO
      ! open file for writing
      if(mpi_world%rank==0) then
        write(filename2,'(I10)') ik
        filename2 = trim(adjustl(filename1))//'_ik_'//trim(adjustl(filename2))
        file = io_open(filename2, namespace, action='write')
      end if
      ! get all states at ik that are locally stored (ground state and td-states)
      do ist=gs_st%st_start,gs_st%st_end
        if(state_kpt_is_local(gs_st, ist, ik)) then
          call states_elec_get_state(st, mesh, ist, ik,temp_state )
          do idim = 1,gs_st%d%dim
            psi(ist,idim,1:mesh%np) =  temp_state(1:mesh%np,idim)
          end do
          call states_elec_get_state(gs_st, mesh, ist, ik, temp_state )
          do idim = 1,gs_st%d%dim
            gs_psi(ist,idim,1:mesh%np) =  temp_state(1:mesh%np,idim)
          end do
        end if
      end do
      ! collect states at ik from all processes in one array
      call comm_allreduce(mpi_world, psi)
      call comm_allreduce(mpi_world, gs_psi)
       
      ! compute the overlaps as a matrix product
      proj(1:gs_st%nst,1:gs_st%nst) = M_ZERO
      call zgemm('n',                                  &
                 'c',                                  &
                 gs_st%nst,                            &
                 gs_st%nst,                            &
                 mesh%np_global*gs_st%d%dim,           &
                 TOCMPLX(mesh%volume_element, M_ZERO), &
                 psi(1, 1, 1),                         &
                 ubound(psi, dim = 1),                 &
                 gs_psi(1, 1, 1),                      &
                 ubound(gs_psi, dim = 1),              &
                 M_z0,                                 &
                 proj(1, 1),                           &
                 ubound(proj, dim = 1))

      ! write to file 
      if(mpi_world%rank==0) then
        do ist = 1,gs_st%nst
          do jst=1,gs_st%nst
            write(file,'(I3,1x,I3,1x,e12.6,1x,e12.6,2x)') ist, jst, proj(ist,jst)
          end do
        end do
        call io_close(file)
      end if

    end do! ik

    SAFE_DEALLOCATE_A(proj)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(gs_psi)
    SAFE_DEALLOCATE_A(temp_state)
     
    POP_SUB(td_write_proj_kp)
  end subroutine td_write_proj_kp

  !---------------------------------------
  subroutine td_write_floquet(namespace, space, hm, gr, st, iter)
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st !< at iter=0 this is the groundstate
    integer,                  intent(in)    :: iter 

    CMPLX, allocatable :: hmss(:,:), psi(:,:,:), hpsi(:,:,:), temp_state1(:,:)
    CMPLX, allocatable :: HFloquet(:,:,:), HFloq_eff(:,:), temp(:,:)
    FLOAT, allocatable :: eigenval(:), bands(:,:)
    character(len=80) :: filename
    integer :: it, nT, ik, ist, in, im, file, idim, nik, ik_count
    integer :: Forder, Fdim, m0, n0, n1, nst, ii, jj, lim_nst
    logical :: downfolding = .false.
    type(mesh_t) :: mesh
    type(states_elec_t) :: hm_st

    FLOAT :: dt, Tcycle, omega

    PUSH_SUB(td_write_floquet)

    ! this does not depend on propagation, so we do it only once 
    if(.not. iter == 0) then
       POP_SUB(td_write_floquet)
       return
    end if

    mesh = gr%mesh
    nst = st%nst

    !for now no domain distributionallowed
    ASSERT(mesh%np == mesh%np_global)

   ! this is used to initialize the hpsi (more effiecient ways?)
    call states_elec_copy(hm_st, st)

    !%Variable TDFloquetFrequency
    !%Type float
    !%Default 0
    !%Section Time-Dependent::TD Output
    !%Description 
    !% Frequency for the Floquet analysis, this should be the carrier frequency or integer multiples of it.
    !% Other options will work, but likely be nonsense.
    !% 
    !%End
    call parse_variable(namespace, 'TDFloquetFrequency', M_ZERO, omega, units_inp%energy)
    call messages_print_var_value(stdout,'Frequency used for Floquet analysis', omega)
    if(abs(omega)<=M_EPSILON) then
       message(1) = "Please give a non-zero value for TDFloquetFrequency"
       call messages_fatal(1, namespace=namespace)
    endif

    ! get time of one cycle
    Tcycle=M_TWO*M_PI/omega

    !%Variable TDFloquetSample
    !%Type integer
    !%Default 20
    !%Section Time-Dependent::TD Output
    !%Description 
    !% Number of points on which one Floquet cycle is sampled in the time-integral of the Floquet analysis.
    !%
    !%End 
    call parse_variable(namespace, 'TDFloquetSample',20 ,nt)
    call messages_print_var_value(stdout,'Number of Floquet time-sampling points', nT)
    dt = Tcycle/TOFLOAT(nT)

    !%Variable TDFloquetDimension
    !%Type integer
    !%Default -1
    !%Section Time-Dependent::TD Output
    !%Description
    !% Order of Floquet Hamiltonian. If negative number is given, downfolding is performed.
    !%End
    call parse_variable(namespace, 'TDFloquetDimension',-1,Forder)
    if(Forder.ge.0) then
       call messages_print_var_value(stdout,'Order of multiphoton Floquet-Hamiltonian', Forder)
       !Dimension of multiphoton Floquet-Hamiltonian
       Fdim = 2*Forder+1
    else
       message(1) = 'Floquet-Hamiltonian is downfolded'
       call messages_info(1)
       downfolding = .true.
       Forder = 1
       Fdim = 3
    endif

    dt = Tcycle/TOFLOAT(nT)

    ! we are only interested for k-point with zero weight
    nik = hm%kpoints%nik_skip

    SAFE_ALLOCATE(hmss(1:nst,1:nst))
    SAFE_ALLOCATE( psi(1:nst,1:st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(hpsi(1:nst,1:st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(temp_state1(1:mesh%np,1:st%d%dim))

    ! multiphoton Floquet Hamiltonian, layout:
    !     (H_{-n,-m} ...  H_{-n,0} ...  H_{-n,m}) 
    !     (    .      .      .      .      .    )
    ! H = (H_{0,-m}  ...  H_{0,0}  ...  H_{0,m} )
    !     (    .      .      .      .      .    )
    !     (H_{n,-m}  ...  H_{n,0}  ...  H_{n,m} )    
    SAFE_ALLOCATE(HFloquet(1:nik,1:nst*Fdim, 1:nst*Fdim))
    HFloquet(1:nik,1:nst*Fdim, 1:nst*Fdim) = M_ZERO

    ! perform time-integral over one cycle
    do it=1,nT
      ! get non-interacting Hamiltonian at time (offset by one cycle to allow for ramp)
      call hamiltonian_elec_update(hm, gr%mesh, namespace, space, time=Tcycle+it*dt)
      ! get hpsi
      call zhamiltonian_elec_apply_all(hm, namespace, gr%mesh, st, hm_st)

      ! project Hamiltonian into grounstates for zero weight k-points
      ik_count = 0

      do ik = hm%kpoints%reduced%npoints-nik+1, hm%kpoints%reduced%npoints
        ik_count = ik_count + 1

        psi(1:nst, 1:st%d%dim, 1:mesh%np)= M_ZERO
        hpsi(1:nst, 1:st%d%dim, 1:mesh%np)= M_ZERO

        do ist=st%st_start,st%st_end
          if(state_kpt_is_local(st, ist, ik)) then
            call states_elec_get_state(st, mesh, ist, ik,temp_state1 )
            do idim = 1, st%d%dim
              psi(ist, idim, 1:mesh%np) =  temp_state1(1:mesh%np,idim)
            end do
            call states_elec_get_state(hm_st, mesh, ist, ik,temp_state1 )
            do idim = 1, st%d%dim
              hpsi(ist, idim,1:mesh%np) = temp_state1(1:mesh%np,idim)
            end do
          end if
        end do
        call comm_allreduce(mpi_world, psi)
        call comm_allreduce(mpi_world, hpsi)
        hmss(1:nst,1:nst) = M_ZERO
        call zgemm( 'n',                                  &
                    'c',                                  &
                    nst,                                  &
                    nst,                                  &
                    mesh%np_global*st%d%dim,              &
                    TOCMPLX(mesh%volume_element, M_ZERO), &
                    hpsi(1, 1, 1),                        &
                    ubound(hpsi, dim = 1),                &
                    psi(1, 1, 1),                         &
                    ubound(psi, dim = 1),                 &
                    M_z0,                                 &
                    hmss(1, 1),                           &
                    ubound(hmss, dim = 1))

        hmss(1:nst,1:nst) = CONJG(hmss(1:nst,1:nst))

        ! accumulate the Floqeut integrals
        do in=-Forder,Forder
           do im=-Forder,Forder
              ii=(in+Forder)*nst
              jj=(im+Forder)*nst
              HFloquet(ik_count,ii+1:ii+nst,jj+1:jj+nst) =  &
                HFloquet(ik_count,ii+1:ii+nst,jj+1:jj+nst) + hmss(1:nst,1:nst)*exp(-(in-im)*M_zI*omega*it*dt)
              ! diagonal term
              if(in==im) then
                 do ist = 1,nst
                    HFloquet(ik_count,ii+ist,ii+ist) = HFloquet(ik_count,ii+ist,ii+ist) + in*omega
                 end do
              end if
           end do
        end do
      end do !ik

    end do ! it

    HFloquet(:,:,:) = M_ONE/nT*HFloquet(:,:,:)

    ! diagonalize Floquet Hamiltonian
    if(downfolding) then
       ! here perform downfolding
       SAFE_ALLOCATE(HFloq_eff(1:nst,1:nst))
       SAFE_ALLOCATE(eigenval(1:nst))
       SAFE_ALLOCATE(bands(1:nik,1:nst))

       HFloq_eff(1:nst,1:nst) = M_ZERO
       do ik=1,nik
          ! the HFloquet blocks are copied directly out of the super matrix
          m0 = nst ! the m=0 start position
          n0 = nst ! the n=0 start postion
          n1 = 2*nst ! the n=+1 start postion
          HFloq_eff(1:nst,1:nst) = HFloquet(ik,n0+1:n0+nst,m0+1:m0+nst) + &
               M_ONE/omega*(matmul(HFloquet(ik,1:nst,m0+1:m0+nst), HFloquet(ik,n1+1:n1+nst,m0+1:m0+nst))- &
                            matmul(HFloquet(ik,n1+1:n1+nst,m0+1:m0+nst), HFloquet(ik,1:nst,m0+1:m0+nst)))

          call lalg_eigensolve(nst, HFloq_eff, eigenval)
          bands(ik,1:nst) = eigenval(1:nst)
       end do
       SAFE_DEALLOCATE_A(HFloq_eff)
    else
      ! the full Floquet 
      SAFE_ALLOCATE(eigenval(1:nst*Fdim))
      SAFE_ALLOCATE(bands(1:nik,1:nst*Fdim))
      SAFE_ALLOCATE(temp(1:nst*Fdim, 1:nst*Fdim))

      do ik=1,nik
         temp(1:nst*Fdim,1:nst*Fdim) = HFloquet(ik,1:nst*Fdim,1:nst*Fdim)
         call lalg_eigensolve(nst*Fdim, temp, eigenval)
         bands(ik,1:nst*Fdim) = eigenval(1:nst*Fdim)
      end do
    end if

    !write bandstructure to file
    if(downfolding) then
      lim_nst = nst
      filename="downfolded_floquet_bands"
    else
       lim_nst = nst*Fdim
       filename="floquet_bands"
    end if
    ! write bands (full or downfolded)
    if(mpi_world%rank==0) then
      file=987254
      file = io_open(filename, namespace, action = 'write')
      do ik=1,nik
        do ist = 1,lim_nst
          write(file,'(e12.6, 1x)',advance='no') bands(ik,ist)
        end do
        write(file,'(1x)')
      end do
      call io_close(file)
    endif
    
    if(.not.downfolding) then
      ! for the full Floquet case compute also the trivially shifted
      ! Floquet bands for reference (i.e. setting H_{nm}=0 for n!=m)
      bands(1:nik,1:nst*Fdim) = M_ZERO
      do ik=1,nik
        temp(1:nst*Fdim,1:nst*Fdim) = M_ZERO
        do jj=0,Fdim-1
          ii=jj*nst
          temp(ii+1:ii+nst,ii+1:ii+nst) = HFloquet(ik,ii+1:ii+nst,ii+1:ii+nst)
        end do
        call lalg_eigensolve(nst*Fdim, temp, eigenval)
        bands(ik,1:nst*Fdim) = eigenval(1:nst*Fdim)
      end do
    
      if(mpi_world%rank==0) then
        filename='trivial_floquet_bands'
        file = io_open(filename, namespace, action = 'write')
        do ik=1,nik
          do ist = 1,lim_nst
            write(file,'(e12.6, 1x)',advance='no') bands(ik,ist)
          end do
          write(file,'(1x)')
        end do
        call io_close(file)
      endif
     end if
  
    ! reset time in Hamiltonian
    call hamiltonian_elec_update(hm, gr%mesh, namespace, space, time=M_ZERO)

    SAFE_DEALLOCATE_A(hmss)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(temp_state1)
    SAFE_DEALLOCATE_A(HFloquet)
    SAFE_DEALLOCATE_A(eigenval)
    SAFE_DEALLOCATE_A(bands)
    SAFE_DEALLOCATE_A(temp)

   POP_SUB(td_write_floquet)

  end subroutine td_write_floquet

  ! ---------------------------------------------------------
  subroutine td_write_total_current(out_total_current, gr, st, iter)
    type(c_ptr),         intent(inout) :: out_total_current
    type(grid_t),        intent(in)    :: gr
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: iter

    integer :: idir, ispin
    character(len=50) :: aux
    FLOAT :: total_current(1:MAX_DIM), abs_current(1:MAX_DIM)

    PUSH_SUB(td_write_total_current)

    if(mpi_grp_is_root(mpi_world) .and. iter == 0) then
      call td_write_print_header_init(out_total_current)

      ! first line: column names
      call write_iter_header_start(out_total_current)
      
      do idir = 1, gr%sb%dim
        write(aux, '(a2,i1,a1)') 'I(', idir, ')'
        call write_iter_header(out_total_current, aux)
      end do

      do idir = 1, gr%sb%dim
        write(aux, '(a2,i1,a1)') 'IntAbs(j)(', idir, ')'
        call write_iter_header(out_total_current, aux)
      end do
      
      do ispin = 1, st%d%nspin
        do idir = 1, gr%sb%dim
          write(aux, '(a4,i1,a1,i1,a1)') 'I-sp', ispin, '(', idir, ')'
          call write_iter_header(out_total_current, aux)
        end do
      end do      

      call write_iter_nl(out_total_current)

      call td_write_print_header_end(out_total_current)
    end if
    
    ASSERT(allocated(st%current))

    if(mpi_grp_is_root(mpi_world)) &
      call write_iter_start(out_total_current)

    total_current = CNST(0.0)
    do idir = 1, gr%sb%dim
      do ispin = 1, st%d%spin_channels
        total_current(idir) =  total_current(idir) + dmf_integrate(gr%mesh, st%current(:, idir, ispin), reduce = .false.)
      end do
      total_current(idir) = units_from_atomic(units_out%length/units_out%time, total_current(idir))
    end do
    if(gr%mesh%parallel_in_domains) then
      call gr%mesh%allreduce(total_current, dim = gr%sb%dim)
    end if

    abs_current = CNST(0.0)
    do idir = 1, gr%sb%dim
      do ispin = 1, st%d%spin_channels
        abs_current(idir) =  abs_current(idir) + dmf_integrate(gr%mesh, abs(st%current(:, idir, ispin)), reduce = .false.)
      end do
      abs_current(idir) = units_from_atomic(units_out%length/units_out%time, abs_current(idir))
    end do
    if(gr%mesh%parallel_in_domains) then
      call gr%mesh%allreduce(abs_current, dim = gr%sb%dim)
    end if

   if(mpi_grp_is_root(mpi_world)) then
      call write_iter_double(out_total_current, total_current, gr%sb%dim)
      call write_iter_double(out_total_current, abs_current, gr%sb%dim)
   end if
  
    do ispin = 1, st%d%nspin
      total_current = CNST(0.0)
      do idir = 1, gr%sb%dim
        total_current(idir) = units_from_atomic(units_out%length/units_out%time, &
                               dmf_integrate(gr%mesh, st%current(:, idir, ispin), reduce = .false.))
      end do
      if(gr%mesh%parallel_in_domains) then
        call gr%mesh%allreduce(total_current, dim = gr%sb%dim)
      end if
      if(mpi_grp_is_root(mpi_world)) &
        call write_iter_double(out_total_current, total_current, gr%sb%dim)
    end do

    if(mpi_grp_is_root(mpi_world)) &
      call write_iter_nl(out_total_current)
      
    POP_SUB(td_write_total_current)
  end subroutine td_write_total_current

  ! ---------------------------------------------------------
  
  subroutine td_write_total_heat_current(write_obj, space, hm, gr, st, iter)
    type(c_ptr),              intent(inout) :: write_obj
    type(space_t),            intent(in)    :: space
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(in)    :: st
    integer,                  intent(in)    :: iter

    integer :: idir, ispin
    character(len=50) :: aux
    FLOAT, allocatable :: heat_current(:, :, :)
    FLOAT :: total_current(1:MAX_DIM)

    PUSH_SUB(td_write_total_current)

    if(mpi_grp_is_root(mpi_world) .and. iter == 0) then
      call td_write_print_header_init(write_obj)

      ! first line: column names
      call write_iter_header_start(write_obj)
      
      do idir = 1, space%dim
        write(aux, '(a2,i1,a1)') 'Jh(', idir, ')'
        call write_iter_header(write_obj, aux)
      end do

      call write_iter_nl(write_obj)

      call td_write_print_header_end(write_obj)
    end if

    SAFE_ALLOCATE(heat_current(1:gr%mesh%np, 1:space%dim, 1:st%d%nspin))  

    call current_heat_calculate(space, gr%der, hm, st, heat_current)
    
    if(mpi_grp_is_root(mpi_world)) call write_iter_start(write_obj)

    total_current = CNST(0.0)
    do idir = 1, space%dim
      do ispin = 1, st%d%spin_channels
        total_current(idir) =  total_current(idir) + dmf_integrate(gr%mesh, heat_current(:, idir, ispin))
      end do
      total_current(idir) = units_from_atomic(units_out%energy*units_out%length/units_out%time, total_current(idir))
    end do

    SAFE_DEALLOCATE_A(heat_current)
    
    if(mpi_grp_is_root(mpi_world)) call write_iter_double(write_obj, total_current, space%dim)
  
    if(mpi_grp_is_root(mpi_world)) call write_iter_nl(write_obj)
      
    POP_SUB(td_write_total_current)
  end subroutine td_write_total_heat_current


  ! ---------------------------------------------------------
  subroutine td_write_partial_charges(out_partial_charges, mesh, st, ions, iter)
    type(c_ptr),             intent(inout) :: out_partial_charges
    type(mesh_t),            intent(in)    :: mesh
    type(states_elec_t),     intent(in)    :: st
    type(ions_t),            intent(in)    :: ions
    integer,                 intent(in)    :: iter

    integer :: idir
    character(len=50) :: aux
    FLOAT, allocatable :: hirshfeld_charges(:)

    PUSH_SUB(td_write_partial_charges)

    SAFE_ALLOCATE(hirshfeld_charges(1:ions%natoms))

    call partial_charges_calculate(mesh, st, ions, hirshfeld_charges = hirshfeld_charges)
        
    if(mpi_grp_is_root(mpi_world)) then

      if(iter == 0) then
        
        call td_write_print_header_init(out_partial_charges)
        
        ! first line: column names
        call write_iter_header_start(out_partial_charges)
        
        do idir = 1, ions%natoms
          write(aux, '(a13,i3,a1)') 'hirshfeld(atom=', idir, ')'
          call write_iter_header(out_partial_charges, aux)
        end do
        
        call write_iter_nl(out_partial_charges)
        
        call td_write_print_header_end(out_partial_charges)
      end if
      
      call write_iter_start(out_partial_charges)
      
      call write_iter_double(out_partial_charges, hirshfeld_charges, ions%natoms)
      
      call write_iter_nl(out_partial_charges)
    end if

    SAFE_DEALLOCATE_A(hirshfeld_charges)

    POP_SUB(td_write_partial_charges)
  end subroutine td_write_partial_charges


  ! ---------------------------------------------------------
  subroutine td_write_effective_u(out_coords, lda_u, iter)
    type(c_ptr),       intent(inout) :: out_coords
    type(lda_u_t),     intent(in) :: lda_u
    integer,           intent(in) :: iter

    integer :: ios, inn
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_effective_u)

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do ios = 1, lda_u%norbsets
        write(aux, '(a2,i3,a1)') 'Ueff(', ios, ')'
        call write_iter_header(out_coords, aux)
      end do
      
      do ios = 1, lda_u%norbsets
        write(aux, '(a2,i3,a1)') 'U(', ios, ')'
        call write_iter_header(out_coords, aux)
      end do

      do ios = 1, lda_u%norbsets
        write(aux, '(a2,i3,a1)') 'J(', ios, ')'
        call write_iter_header(out_coords, aux)
      end do

      if(lda_u%intersite) then
        do ios = 1, lda_u%norbsets
          do inn = 1, lda_u%orbsets(ios)%nneighbors
            write(aux, '(a2,i3,a1,i3,a1)') 'V(', ios,'-', inn, ')'
            call write_iter_header(out_coords, aux)
          end do
        end do
      end if
        

      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(out_coords, &
        'Effective U '   // trim(units_abbrev(units_out%energy))   //   &
        ', U in '// trim(units_abbrev(units_out%energy)) //   &
        ', J in '    // trim(units_abbrev(units_out%energy)))
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    do ios = 1, lda_u%norbsets
      call write_iter_double(out_coords, units_from_atomic(units_out%energy,  &
                                                lda_u%orbsets(ios)%Ueff), 1)
    end do
 
    do ios = 1, lda_u%norbsets
      call write_iter_double(out_coords, units_from_atomic(units_out%energy,  &
                                                lda_u%orbsets(ios)%Ubar), 1)
    end do

    do ios = 1, lda_u%norbsets
      call write_iter_double(out_coords, units_from_atomic(units_out%energy,  &
                                                lda_u%orbsets(ios)%Jbar), 1)
    end do

    if(lda_u%intersite) then
      do ios = 1, lda_u%norbsets
        do inn = 1, lda_u%orbsets(ios)%nneighbors
          call write_iter_double(out_coords, units_from_atomic(units_out%energy,  &
                                              lda_u%orbsets(ios)%V_ij(inn,0)), 1)
        end do
      end do
    end if   

    call write_iter_nl(out_coords)

    POP_SUB(td_write_effective_u)
  end subroutine td_write_effective_u


  ! ---------------------------------------------------------
  subroutine td_write_mxll_init(writ, namespace, iter, dt)
    type(td_write_t),         intent(out)   :: writ
    type(namespace_t),        intent(in)    :: namespace
    integer,                  intent(in)    :: iter
    FLOAT,                    intent(in)    :: dt

    integer :: default, flags, iout, first, idim
    logical :: out_flag(5)

    PUSH_SUB(td_write_mxll_init)

    !%Variable MaxwellTDOutput
    !%Type flag
    !%Default maxwell_energy
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent
    !% Maxwell simulation. Many of the options can increase the computational
    !% cost of the simulation, so only use the ones that you need. In
    !% most cases the default value is enough, as it is adapted to the
    !% details of the TD run.
    !%Option maxwell_energy 1
    !% Output of the electromagnetic field energy into the folder <tt>td.general/maxwell</tt>.
    !%Option maxwell_fields 2
    !% Output of the electromagnetic field at the origin of the simulation box into the
    !% folder <tt>td.general/fields</tt>
    !%Option mean_poynting 4
    !% Output of the mean Poynting vector
    !%Option e_field_surface 8
    !% Output of the E field sliced along the planes x=0, y=0, z=0 for each field component
    !%Option b_field_surface 16
    !% Output of the B field sliced along the planes x=0, y=0, z=0 for each field component
    !%End

    default = 2**(OUT_MAXWELL_ENERGY - 1)
    call parse_variable(namespace, 'MaxwellTDOutput', default, flags)

    if(.not.varinfo_valid_option('MaxwellTDOutput', flags, is_flag = .true.)) &
        call messages_input_error(namespace, 'MaxwellTDOutput')

    do iout = 1, 5
      out_flag(iout) = (iand(flags, 2**(iout - 1)) /= 0)
    end do

    ! TODO: Improve the way the output option labels are handled
    do iout = 1, 3
      writ%out(iout)%write = out_flag(iout)
    end do
    do iout = 4, 5
      if (iout == 4) then
        if (out_flag(4)) then
          do idim=1, 3
            writ%out(4-1+idim)%write = .true.
          end do
        end if
      else if (iout == 5) then
        if (out_flag(5)) then
          do idim=1, 3
            writ%out(5-1+idim)%write = .true.
          end do
        end if
      end if
    end do

    if (iter == 0) then
      first = 0
    else
      first = iter + 1
    end if

    call io_mkdir('td.general', namespace)

    if (writ%out(OUT_MAXWELL_ENERGY)%write) then
       call write_iter_init(writ%out(OUT_MAXWELL_ENERGY)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/maxwell_energy", namespace)))
    end if

    if (writ%out(OUT_MAXWELL_FIELDS)%write) then
       call write_iter_init(writ%out(OUT_MAXWELL_FIELDS)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/fields", namespace)))
    end if

    if (writ%out(OUT_MEAN_POYNTING)%write) then
       call write_iter_init(writ%out(OUT_MEAN_POYNTING)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/mean_poynting_vector", namespace)))
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_X)%write) then
       call write_iter_init(writ%out(OUT_E_FIELD_SURFACE_X)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/electric_field_surface-x", namespace)))
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_Y)%write) then
       call write_iter_init(writ%out(OUT_E_FIELD_SURFACE_Y)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/electric_field_surface-y", namespace)))
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_Z)%write) then
       call write_iter_init(writ%out(OUT_E_FIELD_SURFACE_Z)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/electric_field_surface-z", namespace)))
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_X)%write) then
       call write_iter_init(writ%out(OUT_B_FIELD_SURFACE_X)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/magnetic_field_surface-x", namespace)))
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_Y)%write) then
       call write_iter_init(writ%out(OUT_B_FIELD_SURFACE_Y)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/magnetic_field_surface-y", namespace)))
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_Z)%write) then
       call write_iter_init(writ%out(OUT_B_FIELD_SURFACE_Z)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/magnetic_field_surface-z", namespace)))
    end if

    POP_SUB(td_write_mxll_init)
  end subroutine td_write_mxll_init

  
  ! ---------------------------------------------------------
  subroutine td_write_mxll_end(writ)
    type(td_write_t), intent(inout) :: writ

    integer :: iout

    PUSH_SUB(td_write_mxll_end)

    if(mpi_grp_is_root(mpi_world)) then    
      do iout = 1, OUT_MAXWELL_MAX
        if(writ%out(iout)%write)  call write_iter_end(writ%out(iout)%handle)
      end do
    end if

    POP_SUB(td_write_mxll_end)
  end subroutine td_write_mxll_end
    

  ! ---------------------------------------------------------
  subroutine td_write_mxll_iter(writ, gr, st, hm, dt, iter)
    type(td_write_t),              intent(inout) :: writ
    type(grid_t),                  intent(inout) :: gr
    type(states_mxll_t),           intent(inout) :: st
    type(hamiltonian_mxll_t),      intent(inout) :: hm
    FLOAT,                         intent(in)    :: dt
    integer,                       intent(in)    :: iter

    type(profile_t), save :: prof

    PUSH_SUB(td_write_mxll_iter)

    call profiling_in(prof, "TD_WRITE_ITER_MAXWELL")

    if(writ%out(OUT_MAXWELL_ENERGY)%write) then
!      if (present(hm_elec)) then
!        call td_write_maxwell_energy(writ%out(OUT_MAXWELL_ENERGY)%handle, hm, st, iter, &
!                                             hm, ions%kinetic_energy)
!      else
        call td_write_maxwell_energy(writ%out(OUT_MAXWELL_ENERGY)%handle, hm, iter)
!      end if
    end if

    if (writ%out(OUT_MAXWELL_FIELDS)%write) then
      call td_write_fields(writ%out(OUT_MAXWELL_FIELDS)%handle, st, gr, iter, dt)
    end if

    if (writ%out(OUT_MEAN_POYNTING)%write) then
      call td_write_poynting_vector(writ%out(OUT_MEAN_POYNTING)%handle, st, gr, iter, dt, hm%plane_waves)
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_X)%write) then
      call td_write_electric_field_box_surface(writ%out(OUT_E_FIELD_SURFACE_X)%handle, st, 1, iter)
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_Y)%write) then
      call td_write_electric_field_box_surface(writ%out(OUT_E_FIELD_SURFACE_Y)%handle, st, 2, iter)
    end if

    if (writ%out(OUT_E_FIELD_SURFACE_Z)%write) then
      call td_write_electric_field_box_surface(writ%out(OUT_E_FIELD_SURFACE_Z)%handle, st, 3, iter)
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_X)%write) then
      call td_write_magnetic_field_box_surface(writ%out(OUT_B_FIELD_SURFACE_X)%handle, st, 1, iter)
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_Y)%write) then
      call td_write_magnetic_field_box_surface(writ%out(OUT_B_FIELD_SURFACE_Y)%handle, st, 2, iter)
    end if

    if (writ%out(OUT_B_FIELD_SURFACE_Z)%write) then
      call td_write_magnetic_field_box_surface(writ%out(OUT_B_FIELD_SURFACE_Z)%handle, st, 3, iter)
    end if

    call profiling_out(prof)

    POP_SUB(td_write_mxll_iter)
  end subroutine td_write_mxll_iter


  !----------------------------------------------------------
  subroutine td_dump_mxll(restart, space, mesh, st, hm, iter, ierr, bc_plane_waves)
    type(restart_t),            intent(in)  :: restart
    type(space_t),              intent(in)  :: space
    type(mesh_t),               intent(in)  :: mesh
    type(states_mxll_t),        intent(in)  :: st
    type(hamiltonian_mxll_t),   intent(in)  :: hm
    integer,                    intent(in)  :: iter
    integer,                    intent(out) :: ierr
    logical,                    intent(in)  :: bc_plane_waves

    integer :: err, zff_dim, id, id1, id2, ip_in
    logical :: pml_check
    CMPLX, allocatable :: zff(:,:)

    PUSH_SUB(td_dump_mxll)

    ierr = 0

    pml_check = any(hm%bc%bc_ab_type(1:3) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML)

    if (debug%info) then
      message(1) = "Debug: Writing td_maxwell restart."
      call messages_info(1)
    end if

    if (bc_plane_waves) then
      zff_dim = 2 * st%dim
    else
      zff_dim = 1 * st%dim
    end if
    if (pml_check) then
      zff_dim = zff_dim + 18
    end if

    SAFE_ALLOCATE(zff(1:mesh%np,1:zff_dim))
    zff = M_z0

    if (bc_plane_waves) then
      zff(1:mesh%np, 1:st%dim)   = st%rs_state(1:mesh%np, 1:st%dim)
      zff(1:mesh%np, st%dim+1:st%dim+st%dim) = st%rs_state_plane_waves(1:mesh%np, 1:st%dim)
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id+1
            do ip_in = 1, hm%bc%pml%points_number
              zff(ip_in, 2*st%dim+id) = hm%bc%pml%conv_plus(ip_in, id1, id2)
              zff(ip_in, 2*st%dim+9+id) = hm%bc%pml%conv_minus(ip_in, id1, id2)
            end do
          end do
        end do
       end if
    else
      zff(1:mesh%np, 1:st%dim) = st%rs_state(1:mesh%np, 1:st%dim)
      if (pml_check) then
        id = 0
        do id1 = 1, 3
          do id2 = 1, 3
            id = id+1
            do ip_in = 1, hm%bc%pml%points_number
              zff(ip_in, st%dim+id) = hm%bc%pml%conv_plus(ip_in, id1, id2)
              zff(ip_in, st%dim+9+id) = hm%bc%pml%conv_minus(ip_in, id1, id2)
            end do
          end do
        end do
      end if
    end if

    call states_mxll_dump(restart, st, space, mesh, zff, zff_dim, err, iter)
    if (err /= 0) ierr = ierr + 1

    if (debug%info) then
      message(1) = "Debug: Writing td_maxwell restart done."
      call messages_info(1)
    end if

    SAFE_DEALLOCATE_A(zff)

    POP_SUB(td_dump_mxll)
  end subroutine td_dump_mxll


  ! ---------------------------------------------------------
  subroutine td_write_maxwell_energy(out_maxwell_energy, hm, iter)
    type(c_ptr),                   intent(inout) :: out_maxwell_energy
    type(hamiltonian_mxll_t),      intent(in)    :: hm
    integer,                       intent(in)    :: iter

    integer :: ii

    integer :: n_columns

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_maxwell_energy)

    n_columns = 7

    if(iter == 0) then
      call td_write_print_header_init(out_maxwell_energy)

      ! first line -> column names
      call write_iter_header_start(out_maxwell_energy)
      call write_iter_header(out_maxwell_energy, 'Mx energy')
      call write_iter_header(out_maxwell_energy, 'E energy')
      call write_iter_header(out_maxwell_energy, 'B energy')
      call write_iter_header(out_maxwell_energy, 'Mx energy s. b.')
      call write_iter_header(out_maxwell_energy, 'Mx energy bdry')
      call write_iter_header(out_maxwell_energy, 'Mx energy tr. f.')
      call write_iter_header(out_maxwell_energy, 'Mx energy long. f.')
      call write_iter_header(out_maxwell_energy, 'Mx energy inc. w.')

      call write_iter_nl(out_maxwell_energy)

      ! units

      call write_iter_string(out_maxwell_energy, '#[Iter n.]')
      call write_iter_header(out_maxwell_energy, '[' // trim(units_abbrev(units_out%time)) // ']')

      do ii = 1, n_columns
        call write_iter_header(out_maxwell_energy, '[' // trim(units_abbrev(units_out%energy)) // ']')
      end do
      call write_iter_nl(out_maxwell_energy)
      
      call td_write_print_header_end(out_maxwell_energy)
    end if

    call write_iter_start(out_maxwell_energy)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%energy), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%e_energy), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%b_energy), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, &
         hm%energy%energy+hm%energy%boundaries), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%boundaries), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%energy_trans), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%energy_long), 1)
    call write_iter_double(out_maxwell_energy, units_from_atomic(units_out%energy, hm%energy%energy_plane_waves), 1)
    call write_iter_nl(out_maxwell_energy)

    POP_SUB(td_write_maxwell_energy)
  end subroutine td_write_maxwell_energy


  ! ---------------------------------------------------------
  subroutine td_write_electric_field_box_surface(out_field_surf, st, dim, iter)
    type(c_ptr),                   intent(inout) :: out_field_surf
    type(states_mxll_t),           intent(in)    :: st
    integer,                       intent(in)    :: dim
    integer,                       intent(in)    :: iter

    integer :: ii

    integer :: n_columns

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_electric_field_box_surface)

    n_columns = 12

    if(iter == 0) then
      call td_write_print_header_init(out_field_surf)

      ! first line -> column names
      call write_iter_header_start(out_field_surf)
      call write_iter_header(out_field_surf, '- x direction')
      call write_iter_header(out_field_surf, '+ x direction')
      call write_iter_header(out_field_surf, '- y direction')
      call write_iter_header(out_field_surf, '+ y direction')
      call write_iter_header(out_field_surf, '- z direction')
      call write_iter_header(out_field_surf, '+ z direction')
      call write_iter_header(out_field_surf, '- x dir. p. w.')
      call write_iter_header(out_field_surf, '+ x dir. p. w.')
      call write_iter_header(out_field_surf, '- y dir. p. w.')
      call write_iter_header(out_field_surf, '+ y dir. p. w.')
      call write_iter_header(out_field_surf, '- z dir. p. w.')
      call write_iter_header(out_field_surf, '+ z dir. p. w.')

      call write_iter_nl(out_field_surf)

      ! units
      call write_iter_string(out_field_surf, '#[Iter n.]')
      call write_iter_header(out_field_surf, '[' // trim(units_abbrev(units_out%time)) // ']')

      do ii = 1, n_columns
        call write_iter_header(out_field_surf, '[' // trim(units_abbrev(units_out%energy/units_out%length)) // ']')
      end do
      call write_iter_nl(out_field_surf)

      call td_write_print_header_end(out_field_surf)
    end if

    call write_iter_start(out_field_surf)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(1,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(2,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(1,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(2,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(1,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface(2,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(1,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(2,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(1,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(2,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(1,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(units_out%energy/units_out%length, &
         st%electric_field_box_surface_plane_waves(2,3,dim)), 1)
    call write_iter_nl(out_field_surf)

    POP_SUB(td_write_electric_field_box_surface)
  end subroutine td_write_electric_field_box_surface


  ! ---------------------------------------------------------
  subroutine td_write_magnetic_field_box_surface(out_field_surf, st, dim, iter)
    type(c_ptr),                   intent(inout) :: out_field_surf
    type(states_mxll_t),           intent(in)    :: st
    integer,                       intent(in)    :: dim
    integer,                       intent(in)    :: iter

    integer :: ii

    integer :: n_columns

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_magnetic_field_box_surface)

    n_columns = 12

    if(iter == 0) then
      call td_write_print_header_init(out_field_surf)

      ! first line -> column names
      call write_iter_header_start(out_field_surf)
      call write_iter_header(out_field_surf, '- x direction')
      call write_iter_header(out_field_surf, '+ x direction')
      call write_iter_header(out_field_surf, '- y direction')
      call write_iter_header(out_field_surf, '+ y direction')
      call write_iter_header(out_field_surf, '- z direction')
      call write_iter_header(out_field_surf, '+ z direction')
      call write_iter_header(out_field_surf, '- x dir. p. w.')
      call write_iter_header(out_field_surf, '+ x dir. p. w.')
      call write_iter_header(out_field_surf, '- y dir. p. w.')
      call write_iter_header(out_field_surf, '+ y dir. p. w.')
      call write_iter_header(out_field_surf, '- z dir. p. w.')
      call write_iter_header(out_field_surf, '+ z dir. p. w.')

      call write_iter_nl(out_field_surf)

      ! units
      call write_iter_string(out_field_surf, '#[Iter n.]')
      call write_iter_header(out_field_surf, '[' // trim(units_abbrev(units_out%time)) // ']')

      do ii = 1, n_columns
        call write_iter_header(out_field_surf, '[' // trim(units_abbrev(unit_one/units_out%length**2)) // ']')
      end do
      call write_iter_nl(out_field_surf)

      call td_write_print_header_end(out_field_surf)
    end if

    call write_iter_start(out_field_surf)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(1,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(2,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(1,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(2,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(1,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface(2,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(1,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(2,1,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(1,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(2,2,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(1,3,dim)), 1)
    call write_iter_double(out_field_surf, units_from_atomic(unit_one/units_out%length**2, &
         st%magnetic_field_box_surface_plane_waves(2,3,dim)), 1)
    call write_iter_nl(out_field_surf)

    POP_SUB(td_write_magnetic_field_box_surface)
  end subroutine td_write_magnetic_field_box_surface


 ! ---------------------------------------------------------
  subroutine td_write_poynting_vector(out_poynting, st, gr, iter, dt, plane_wave_flag)
    type(c_ptr),         intent(inout) :: out_poynting
    type(states_mxll_t), intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    logical,             intent(in)    :: plane_wave_flag

    integer            :: idir
    FLOAT              :: field(MAX_DIM), field_2(MAX_DIM)
    FLOAT, allocatable :: dtmp(:,:)
    character(len=80)  :: aux

    PUSH_SUB(td_write_poynting_vector)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    SAFE_ALLOCATE(dtmp(1:gr%mesh%np,1:st%dim))

    call get_poynting_vector(gr, st, st%rs_state, st%rs_sign, dtmp, ep_field=st%ep, mu_field=st%mu, &
      mean_value=field)

    if (plane_wave_flag) then
      call get_poynting_vector_plane_waves(gr, st, st%rs_sign, dtmp, mean_value=field)
    end if

    SAFE_DEALLOCATE_A(dtmp)

    if(iter == 0) then
      call td_write_print_header_init(out_poynting)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, dt), &
        " [", trim(units_abbrev(units_out%time)), "]"
      call write_iter_string(out_poynting, aux)
      call write_iter_nl(out_poynting)

      call write_iter_header_start(out_poynting)
      do idir = 1, gr%sb%dim
        write(aux, '(a,i1,a)') 'poynting (', idir, ')'
        call write_iter_header(out_poynting, aux)
      end do
      if (plane_wave_flag) then
        do idir = 1, gr%sb%dim
          write(aux, '(a,i1,a)') 'poynting pl. w.(', idir, ')'
          call write_iter_header(out_poynting, aux)
        end do
      end if

      call write_iter_nl(out_poynting)
      call write_iter_string(out_poynting, '#[Iter n.]')
      call write_iter_header(out_poynting, '[' // trim(units_abbrev(units_out%time)) // ']')

      aux = '[' // trim(units_abbrev(units_out%force)) // ']'
      do idir = 1, 2 * gr%sb%dim
        call write_iter_header(out_poynting, aux)
      end do
      call write_iter_nl(out_poynting)
      call td_write_print_header_end(out_poynting)
    end if

    call write_iter_start(out_poynting)

    ! Output of mean poynting vector
    field = units_from_atomic(unit_one/units_out%length**2, field)
    call write_iter_double(out_poynting, field, gr%sb%dim)

    ! Output of mean poynting vector plane wave
    if (plane_wave_flag) then
      field_2 = units_from_atomic(unit_one/units_out%length**2, field_2)
      call write_iter_double(out_poynting, field_2, gr%sb%dim)
    end if

    call write_iter_nl(out_poynting)

    POP_SUB(td_write_poynting_vector)
  end subroutine td_write_poynting_vector


 ! ---------------------------------------------------------
  subroutine td_write_fields(out_fields, st, gr, iter, dt)
    type(c_ptr),         intent(inout) :: out_fields
    type(states_mxll_t),      intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer :: idir
    FLOAT :: field(gr%sb%dim)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_fields)

    if(iter == 0) then
      call td_write_print_header_init(out_fields)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, dt), &
        " [", trim(units_abbrev(units_out%time)), "]"
      call write_iter_string(out_fields, aux)
      call write_iter_nl(out_fields)

      write(aux, '(a10)') '# position'
      call write_iter_nl(out_fields)
      call write_iter_header(out_fields, aux)

      write(aux, '(a10)') '# position'
      call write_iter_nl(out_fields)
      call write_iter_header(out_fields, aux)

      write(aux, '(a10)') '# position'
      call write_iter_nl(out_fields)
      call write_iter_header(out_fields, aux)

      call write_iter_header_start(out_fields)
      do idir = 1, gr%sb%dim
        write(aux, '(a,i1,a)') 'E(', idir, ')'
        call write_iter_header(out_fields, aux)
      end do
      do idir = 1, gr%sb%dim
        write(aux, '(a,i1,a)') 'B(', idir, ')'
        call write_iter_header(out_fields, aux)
      end do
      do idir = 1, gr%sb%dim
        write(aux, '(a,i1,a)') 'E(', idir, ')'
        call write_iter_header(out_fields, aux)
      end do
      do idir = 1, gr%sb%dim
        write(aux, '(a,i1,a)') 'B(', idir, ')'
        call write_iter_header(out_fields, aux)
      end do

      call write_iter_nl(out_fields)
      call write_iter_string(out_fields, '#[Iter n.]')
      call write_iter_header(out_fields, '[' // trim(units_abbrev(units_out%time)) // ']')

      ! Note that we do not print out units of E, B, or A, but rather units of e*E, e*B, e*A.
      ! (force, force, and energy, respectively). The reason is that the units of E, B or A
      ! are ugly.
      aux = '[' // trim(units_abbrev(units_out%force)) // ']'
      do idir = 1, 4 * gr%sb%dim
        call write_iter_header(out_fields, aux)
      end do
      call write_iter_nl(out_fields)
      call td_write_print_header_end(out_fields)
    end if

    call write_iter_start(out_fields)

    ! Output of electric field at selected point
    call get_electric_field_vector(st%selected_points_rs_state(:,1), field(1:st%dim))
    field(1:st%dim) = units_from_atomic(units_out%energy/units_out%length, field(1:st%dim))
    call write_iter_double(out_fields, field(1:st%dim), gr%sb%dim)
    ! Output of magnetic field at selected point
    call get_magnetic_field_vector(st%selected_points_rs_state(:,1), st%rs_sign, field(1:st%dim))
    field(1:st%dim) = units_from_atomic(unit_one/units_out%length**2, field(1:st%dim))
    call write_iter_double(out_fields, field(1:st%dim), gr%sb%dim)

    ! Output of transverse electric field at selected point
    call get_electric_field_vector(st%selected_points_rs_state_trans(:,1), field(1:st%dim))
    field(1:st%dim) = units_from_atomic(units_out%energy/units_out%length, field(1:st%dim))
    call write_iter_double(out_fields, field(1:st%dim), gr%sb%dim)
    ! Output of transverse magnetic field at selected point
    call get_magnetic_field_vector(st%selected_points_rs_state_trans(:,1), &
         st%rs_sign, field(1:st%dim))
    field(1:st%dim) = units_from_atomic(unit_one/units_out%length**2, field(1:st%dim))
    call write_iter_double(out_fields, field(1:st%dim), gr%sb%dim)

    call write_iter_nl(out_fields)

    POP_SUB(td_write_fields)
  end subroutine td_write_fields


  !----------------------------------------------------------
  subroutine td_write_mxll_free_data(writ, namespace, space, gr, st, hm, ions, outp, clock)
    type(td_write_t),         intent(inout) :: writ
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(inout) :: gr
    type(states_mxll_t),      intent(inout) :: st
    type(hamiltonian_mxll_t), intent(inout) :: hm
    type(ions_t),             intent(inout) :: ions
    type(output_t),           intent(in)    :: outp
    type(clock_t),            intent(in)    :: clock

    character(len=256) :: filename
    integer :: iout
    type(profile_t), save :: prof

    PUSH_SUB(td_write_maxwell_free_data)
    call profiling_in(prof, "TD_WRITE_MAXWELL_DATA")

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAXWELL_MAX
        if(writ%out(iout)%write)  call write_iter_flush(writ%out(iout)%handle)
      end do
    end if

    ! now write down the rest
    write(filename, '(a,a,i7.7)') trim(outp%iter_dir),"td.", clock%get_tick()  ! name of directory

    call output_mxll(outp, namespace, space, gr, st, hm, clock%time(), ions, filename)

    call profiling_out(prof)
    POP_SUB(td_write_maxwell_free_data)
  end subroutine td_write_mxll_free_data


  ! ---------------------------------------------------------
  subroutine td_write_print_header_init(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(td_write_print_header_init)

    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_init)
  end subroutine td_write_print_header_init


  ! ---------------------------------------------------------
  subroutine td_write_print_header_end(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(td_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_end)
  end subroutine td_write_print_header_end


end module td_write_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
