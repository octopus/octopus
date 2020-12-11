!! Copyright (C) 2014 Alain Delgado Gran, Carlo Andrea Rozzi, Stefano Corni, Gabriel Gil
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

module pcm_oct_m
  use comm_oct_m
  use global_oct_m
  use geometry_oct_m
  use grid_oct_m
  use io_oct_m
  use index_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_interpolation_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use pcm_eom_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use varinfo_oct_m
  
  ! to output debug info
  use unit_oct_m
  use unit_system_oct_m
  use io_function_oct_m
  use mesh_function_oct_m

  implicit none

  private

  public ::                 &
    pcm_t,                  &
    pcm_min_t,              &
    pcm_sphere_t,           &
    pcm_init,               &
    pcm_end,                &
    pcm_charges,            &
    pcm_get_vdw_radius,     &
    pcm_pot_rs,             &
    pcm_elect_energy,       &
    pcm_v_nuclei_cav,       &
    pcm_v_cav_li,           &
    pcm_update,             &
    pcm_calc_pot_rs,        &
    pcm_charge_density,     &
    pcm_dipole,             &
    pcm_field,              &
    pcm_eps,                &
    pcm_min_input_parsing_for_spectrum, &
    PCM_TD_EQ,              &
    PCM_TD_NEQ,             &
    PCM_TD_EOM


  !> The cavity hosting the solute molecule is built from a set of 
  !! interlocking spheres with optimized radii centered at the nuclear positions.  
  type pcm_sphere_t
    private
    FLOAT :: x !
    FLOAT :: y !< center of the sphere
    FLOAT :: z !
    FLOAT :: r !< radius of the sphere (different for each species)
  end type pcm_sphere_t

  !> The resulting cavity is discretized by a set of tesserae defined in 3d.  
  integer, parameter :: PCM_DIM_SPACE = 3

  type pcm_t
    private
    logical, public                  :: run_pcm          !< If .true., PCM calculation is enabled
    integer, public                  :: tdlevel          !< flag to use non-equilibrium TD-PCM, either inertial/dynamic split or EOM
    integer                          :: n_spheres        !< Number of spheres used to build the VdW cavity
    integer, public                  :: n_tesserae       !< Total number of tesserae
    type(pcm_sphere_t),  allocatable :: spheres(:)       !< See type pcm_sphere_t
    type(pcm_tessera_t), allocatable, public :: tess(:)          !< See type pcm_tessera_t
    FLOAT                            :: scale_r          !< scaling factor for the radii of the spheres used in PCM
    FLOAT, allocatable               :: matrix(:,:)      !< static PCM response matrix (for epsilon_0)
    FLOAT, allocatable               :: matrix_d(:,:)    !< dynamical PCM response matrix (for epsilon_infty)
    FLOAT, allocatable               :: matrix_lf(:,:)   !< static PCM response matrix (for epsilon_0)        - local field effects
    FLOAT, allocatable               :: matrix_lf_d(:,:) !< dynamical PCM response matrix (for epsilon_infty) - local field effects
    FLOAT, allocatable, public       :: q_e(:)           !< polarization charges due to the solute electrons
    FLOAT, allocatable               :: q_n(:)           !< polarization charges due to the solute nuclei
    FLOAT, allocatable, public       :: q_e_in(:)        !< inertial polarization charges due to the solute electrons
    FLOAT, allocatable               :: rho_e(:)         !< polarization density due to the solute electrons
    FLOAT, allocatable               :: rho_n(:)         !< polarization density due to the solute nuclei
    FLOAT                            :: qtot_e           !< total polarization charge due to electrons
    FLOAT                            :: qtot_n           !< total polarization charge due to nuclei
    FLOAT                            :: qtot_e_in        !< total inertial polarization charge due to electrons
    FLOAT                            :: q_e_nominal      !< total (nominal) electronic charge
    FLOAT                            :: q_n_nominal      !< total (nominal) nuclear charge
    logical                          :: renorm_charges   !< flag to renormalized polarization charges
    FLOAT                            :: q_tot_tol        !< tolerance to trigger normalization of the polarization charges
    FLOAT                            :: deltaQ_e         !< difference between the calculated and nominal electronic charge
    FLOAT                            :: deltaQ_n         !< difference between the calculated and nominal nuclear charge
    FLOAT, allocatable               :: v_e(:)           !< Hartree potential at each tessera
    FLOAT, allocatable               :: v_n(:)           !< Nuclear potential at each tessera
    FLOAT, allocatable, public       :: v_e_rs(:)        !< PCM potential in real-space produced by q_e(:)
    FLOAT, allocatable, public       :: v_n_rs(:)        !< PCM potential in real-space produced by q_n(:)
    FLOAT, allocatable               :: q_ext(:)         !< polarization charges due to an ext. pot.
    FLOAT, allocatable               :: q_ext_in(:)      !< inertial polarization charges due to an ext. pot.
    FLOAT, allocatable               :: rho_ext(:)       !< polarization density due to an ext. pot.
    FLOAT                            :: qtot_ext         !< total polarization charge due to an ext. pot.
    FLOAT                            :: qtot_ext_in      !< total inertial polarization charge due to an ext. pot.
    FLOAT, allocatable               :: v_ext(:)         !< external potential at each tessera
    FLOAT, allocatable, public       :: v_ext_rs(:)      !< PCM potential in real-space produced by q_ext(:)
    FLOAT, allocatable               :: q_kick(:)        !< polarization charges due to kick
    FLOAT, allocatable               :: rho_kick(:)      !< polarization density due to kick   
    FLOAT                            :: qtot_kick        !< total polarization charge due to kick
    FLOAT, allocatable               :: v_kick(:)        !< kick potential at each tessera
    FLOAT, allocatable, public       :: v_kick_rs(:)     !< PCM potential in real-space produced by q_kick(:)
    FLOAT, public                    :: epsilon_0        !< Static dielectric constant of the solvent
    FLOAT, public                    :: epsilon_infty    !< Infinite-frequency dielectric constant of the solvent
    integer, public                  :: which_eps        !< Dielectric function model, either Debye or Drude-Lorentz
    type(debye_param_t)              :: deb              !< Debye parameters
    type(drude_param_t)              :: drl              !< Drude-Lorentz parameters
    logical, public                  :: localf           !< Logical flag to include polarization charges due to external field
    logical, public                  :: solute           !< Logical flag to include polarization charges due to the solute
    logical                          :: kick_is_present  !< .true. if there are kicks in the calculation
                                                         !< (if localf is .false. this is irrelevant to PCM)
    logical, public                  :: kick_like        !< Logical flag to consider kick-like polarization due to kick
    integer                          :: initial_asc      !< Flag to read or not pol.charges from input file
    FLOAT                            :: gaussian_width   !< Parameter to change the width of density of polarization charges
    integer                          :: info_unit        !< unit for pcm info file
    integer, public                  :: counter          !< used to print the number of SCF or TD iterations in energy_calc
    character(len=80)                :: input_cavity     !< file name containing the geometry of the VdW cavity
    integer                          :: update_iter      !< how often the pcm potential is updated
    integer, public                  :: iter             !< update iteration counter    
    integer                          :: calc_method      !< which method should be used to obtain the pcm potential
    integer                          :: tess_nn          !< number of tessera center mesh-point nearest neighbors
    FLOAT, public                    :: dt               !< time-step of propagation
    type(namespace_t), pointer       :: namespace        !< namespace, needed for output
  end type pcm_t

  type pcm_min_t
    ! Components are public by default
    logical              :: run_pcm   !< If .true., PCM calculation is enabled
    logical              :: localf    !< flag to include polarization charges due to external field
    integer              :: tdlevel   !< flag to use non-equilibrium TD-PCM, either inertial/dynamic split or EOM
    integer              :: which_eps !< dielectric function model, either Debye or Drude-Lorentz
    type(debye_param_t)	 :: deb 	    !< Debye parameters
    type(drude_param_t)	 :: drl 	    !< Drude-Lorentz parameters
  end type pcm_min_t

  FLOAT, allocatable :: s_mat_act(:,:) !< S_I matrix 
  FLOAT, allocatable :: d_mat_act(:,:) !< D_I matrix
  FLOAT, allocatable :: Sigma(:,:)     !< S_E matrix
  FLOAT, allocatable :: Delta(:,:)     !< D_E matrix in JCP 139, 024105 (2013).

  logical            :: gamess_benchmark !< Decide to output pcm_matrix in a GAMESS format 
  FLOAT, allocatable :: mat_gamess(:,:)  !< PCM matrix formatted to be inputed to GAMESS

  integer, parameter :: &
    PCM_TD_EQ  = 0,     &
    PCM_TD_NEQ = 1,     &
    PCM_TD_EOM = 2

  integer, parameter, public :: &
    PCM_CALC_DIRECT  = 1,       &
    PCM_CALC_POISSON = 2

  integer, parameter ::    &
    PCM_VDW_OPTIMIZED = 1, &
    PCM_VDW_SPECIES   = 2

  integer, parameter :: N_TESS_SPHERE = 60 !< minimum number of tesserae per sphere
                                           !< cannot be changed without changing cav_gen subroutine

contains

  !-------------------------------------------------------------------------------------------------------
  !> Initializes the PCM calculation: reads the VdW molecular cavity and generates the PCM response matrix.
  subroutine pcm_init(pcm, namespace, geo, grid, qtot, val_charge, external_potentials_present, kick_present)
    type(pcm_t),               intent(out) :: pcm
    type(geometry_t),          intent(in)  :: geo
    type(namespace_t), target, intent(in)  :: namespace
    type(grid_t),              intent(in)  :: grid
    FLOAT,                     intent(in)  :: qtot
    FLOAT,                     intent(in)  :: val_charge
    logical,                   intent(in)  :: external_potentials_present
    logical,                   intent(in)  :: kick_present

    integer :: ia, ii, itess, jtess, pcm_vdw_type, subdivider
    integer :: cav_unit_test, iunit, pcmmat_unit
    integer :: pcmmat_gamess_unit, cav_gamess_unit
    FLOAT   :: min_distance      

    integer, parameter :: MXTS = 10000

    FLOAT :: default_value
    FLOAT :: vdw_radius

    type(pcm_tessera_t), allocatable :: dum2(:)

    logical :: band
    logical :: add_spheres_h
    logical :: changed_default_nn

    integer :: default_nn
    FLOAT   :: max_area
    
    PUSH_SUB(pcm_init)

    pcm%kick_is_present = kick_present

    pcm%localf = .false.
    pcm%iter = 0
    pcm%update_iter = 1
    pcm%kick_like = .false.

    pcm%namespace => namespace

    !%Variable PCMCalculation
    !%Type logical
    !%Default no
    !%Section Hamiltonian::PCM
    !%Description
    !% If true, the calculation is performed accounting for solvation effects
    !% by using the Integral Equation Formalism Polarizable Continuum Model IEF-PCM
    !% formulated in real-space and real-time (<i>J. Chem. Phys.</i> <b>143</b>, 144111 (2015),
    !% <i>Chem. Rev.</i> <b>105</b>, 2999 (2005), <i>J. Chem. Phys.</i> <b>139</b>, 024105 (2013)).
    !% At the moment, this option is available only for <tt>TheoryLevel = DFT</tt>.
    !% PCM is tested for CalculationMode = gs, while still experimental for other values (in particular, CalculationMode = td).
    !%End
    call parse_variable(namespace, 'PCMCalculation', .false., pcm%run_pcm)
    if (pcm%run_pcm) then
      call messages_print_stress(stdout, trim('PCM'), namespace=namespace)
      if ( (grid%sb%box_shape /= MINIMUM) .or. (grid%sb%dim /= PCM_DIM_SPACE) ) then
        message(1) = "PCM is only available for BoxShape = minimum and 3d calculations"
        call messages_fatal(1, namespace=namespace)
      end if
    else
      POP_SUB(pcm_init)
      return
    end if

    !%Variable PCMVdWRadii
    !%Type integer
    !%Default pcm_vdw_optimized
    !%Section Hamiltonian::PCM
    !%Description
    !% This variable selects which van der Waals radius will be used to generate the solvent cavity.
    !%Option pcm_vdw_optimized  1
    !% Use the van der Waals radius optimized by Stefan Grimme in J. Comput. Chem. 27: 1787-1799, 2006,
    !% except for C, N and O, reported in J. Chem. Phys. 120, 3893 (2004).
    !%Option pcm_vdw_species  2
    !% The vdW radii are set from the <tt>share/pseudopotentials/elements</tt> file. These values are obtained from
    !% Alvarez S., Dalton Trans., 2013, 42, 8617-8636. Values can be changed in the <tt>Species</tt> block.
    !%End
    call parse_variable(namespace, 'PCMVdWRadii', PCM_VDW_OPTIMIZED, pcm_vdw_type)
    call messages_print_var_option(stdout, "PCMVdWRadii", pcm_vdw_type)

    select case (pcm_vdw_type)
    case (PCM_VDW_OPTIMIZED)
      default_value = CNST(1.2)
    case (PCM_VDW_SPECIES)
      default_value = M_ONE
    end select
   
    !%Variable PCMRadiusScaling
    !%Type float
    !%Section Hamiltonian::PCM
    !%Description
    !% Scales the radii of the spheres used to build the solute cavity surface.
    !% The default value depends on the choice of <tt>PCMVdWRadii</tt>:
    !% 1.2 for <tt>pcm_vdw_optimized</tt> and 1.0 for <tt>pcm_vdw_species</tt>.
    !%End
    call parse_variable(namespace, 'PCMRadiusScaling', default_value, pcm%scale_r)
    call messages_print_var_value(stdout, "PCMRadiusScaling", pcm%scale_r)

    !%Variable PCMTDLevel
    !%Type integer
    !%Default eq
    !%Section Hamiltonian::PCM
    !%Description
    !% When CalculationMode=td, PCMTDLevel it sets the way the time-depenendent solvent polarization is propagated.
    !%Option eq 0
    !% If PCMTDLevel=eq, the solvent is always in equilibrium with the solute or the external field, i.e.,
    !% the solvent polarization follows instantaneously the changes in solute density or in the external field.
    !% PCMTDLevel=neq and PCMTDLevel=eom are both nonequilibrium runs.
    !%Option neq 1
    !% If PCMTDLevel=neq, solvent polarization charges are splitted in two terms:
    !% one that follows instantaneously the changes in the solute density or in the external field (dynamical polarization charges),
    !% and another that lag behind in the evolution w.r.t. the solute density or the external field (inertial polarization charges).
    !%Option eom 2
    !% If PCMTDLevel=eom, solvent polarization charges evolves following an equation of motion, generalizing 'neq' propagation.
    !% The equation of motion used here depends on the value of PCMEpsilonModel.
    !%End
    call parse_variable(namespace, 'PCMTDLevel', PCM_TD_EQ, pcm%tdlevel)
    call messages_print_var_value(stdout, "PCMTDLevel", pcm%tdlevel)

    if (pcm%tdlevel /= PCM_TD_EQ .and. (.not. pcm%run_pcm)) then
      call messages_write('Sorry, you have set PCMTDLevel /= eq, but PCMCalculation = no.')
      call messages_new_line()
      call messages_write('To spare you some time, Octopus will proceed as if PCMCalculation = yes.')       
      call messages_warning(namespace=namespace)
    endif
    
    !%Variable PCMStaticEpsilon
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Static dielectric constant of the solvent (<math>\varepsilon_0</math>). 1.0 indicates gas phase.
    !%End
    call parse_variable(namespace, 'PCMStaticEpsilon', M_ONE, pcm%epsilon_0)
    call messages_print_var_value(stdout, "PCMStaticEpsilon", pcm%epsilon_0)

    !%Variable PCMDynamicEpsilon
    !%Type float
    !%Default PCMStaticEpsilon
    !%Section Hamiltonian::PCM
    !%Description
    !% High-frequency dielectric constant of the solvent (<math>\varepsilon_d</math>). 
    !% <math>\varepsilon_d=\varepsilon_0</math> indicate equilibrium with solvent. 
    !%End
    call parse_variable(namespace, 'PCMDynamicEpsilon', pcm%epsilon_0, pcm%epsilon_infty)
    call messages_print_var_value(stdout, "PCMDynamicEpsilon", pcm%epsilon_infty)
    
    !%Variable PCMEpsilonModel
    !%Type integer
    !%Default pcm_debye
    !%Section Hamiltonian::PCM
    !%Description
    !% Define the dielectric function model.
    !%Option pcm_debye 1
    !% Debye model: <math>\varepsilon(\omega)=\varepsilon_d+\frac{\varepsilon_0-\varepsilon_d}{1-i\omega\tau}</math>
    !%Option pcm_drude 2
    !% Drude-Lorentz model: <math>\varepsilon(\omega)=1+\frac{A}{\omega_0^2-\omega^2+i\gamma\omega}</math>
    !%End
    call parse_variable(namespace, 'PCMEpsilonModel', PCM_DEBYE_MODEL, pcm%which_eps)
    call messages_print_var_value(stdout, "PCMEpsilonModel", pcm%which_eps)
    if (.not. varinfo_valid_option('PCMEpsilonModel', pcm%which_eps)) call messages_input_error(namespace, 'PCMEpsilonModel')

    if (pcm%tdlevel == PCM_TD_EOM .and. pcm%which_eps == PCM_DRUDE_MODEL) &
      call messages_experimental("Drude-Lorentz EOM-PCM is experimental")

    if (pcm%which_eps /= PCM_DEBYE_MODEL .and. pcm%which_eps /= PCM_DRUDE_MODEL) then
      call messages_write('Sorry, only Debye or Drude-Lorentz dielectric models are available.')
      call messages_new_line()
      call messages_write('To spare you some time, Octopus will proceed with the default choice (Debye).')
      call messages_new_line()
      call messages_write('You may change PCMEpsilonModel value for a Drude-Lorentz run.')
      call messages_warning(namespace=namespace)
    end if

    if (pcm%tdlevel /= PCM_TD_EQ .and. pcm%which_eps == PCM_DEBYE_MODEL .and. pcm%epsilon_0 == pcm%epsilon_infty) then
      call messages_write('Sorry, inertial/dynamic polarization splitting scheme for TD-PCM or Debye equation-of-motion TD-PCM')
      call messages_new_line()
      call messages_write('require both static and dynamic dielectric constants, and they must be different.')
      call messages_new_line()
      call messages_write('Octopus will run using TD-PCM version in equilibrium with solvent at each time.')        
      call messages_warning(namespace=namespace)
      pcm%tdlevel = PCM_TD_EQ
    end if

    !%Variable PCMEoMInitialCharges
    !%Type integer
    !%Default 0
    !%Section Hamiltonian::PCM
    !%Description
    !% If =0 the propagation of the solvent polarization charges starts from internally generated initial charges 
    !%  in equilibrium with the initial potential.
    !% For Debye EOM-PCM, if >0 the propagation of the solvent polarization charges starts from initial charges from input file.
    !% 										if =1, initial pol. charges due to solute electrons are read from input file.
    !% 										else if =2, initial pol. charges due to external potential are read from input file.
    !% 										else if =3, initial pol. charges due to solute electrons and external potential are read from input file.
    !% Files should be located in pcm directory and are called ASC_e.dat and ASC_ext.dat, respectively.
    !% The latter files are generated after any PCM run and contain the last values of the polarization charges.
    !%End
    call parse_variable(namespace, 'PCMEoMInitialCharges', 0, pcm%initial_asc)
    call messages_print_var_value(stdout, "PCMEoMInitialCharges", pcm%initial_asc)

    if (pcm%initial_asc /= 0) then
      call messages_experimental("The use of external initialization for the EOM-PCM charges is experimental")
      if ((pcm%tdlevel /= PCM_TD_EOM .or. (pcm%tdlevel == PCM_TD_EOM .and. pcm%which_eps /= PCM_DEBYE_MODEL))) then
        call messages_write('Sorry, initial polarization charges can only be read from input file for a Debye EOM-PCM run.')
        call messages_new_line()
        call messages_write('To spare you some time, Octopus will proceed as if PCMEoMInitialCharges = 0.')
        call messages_warning(namespace=namespace)
        pcm%initial_asc = 0
      endif
    endif

    !> packing Debye parameters for convenience
    pcm%deb%eps_0 = pcm%epsilon_0
    pcm%deb%eps_d = pcm%epsilon_infty

    !< re-parse TDTimeStep to propagate polarization charges
    call parse_variable(namespace, 'TDTimeStep', M_ZERO, pcm%dt, unit = units_inp%time)

    !%Variable PCMDebyeRelaxTime
    !%Type float
    !%Default 0.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Relaxation time of the solvent within Debye model (<math>\tau</math>). Recall Debye dieletric function: 
    !% <math>\varepsilon(\omega)=\varepsilon_d+\frac{\varepsilon_0-\varepsilon_d}{1-i\omega\tau}</math>    
    !%End
    call parse_variable(namespace, 'PCMDebyeRelaxTime', M_ZERO, pcm%deb%tau)
    call messages_print_var_value(stdout, "PCMDebyeRelaxTime", pcm%deb%tau)

    if (pcm%tdlevel == PCM_TD_EOM .and. pcm%which_eps == PCM_DEBYE_MODEL .and. &
      (abs(pcm%deb%tau) <= M_EPSILON .or. pcm%deb%eps_0 == pcm%deb%eps_d)) then
      call messages_write('Sorry, you have set PCMTDLevel = eom, but you have not included all required Debye model parameters.')
      call messages_new_line()
      call messages_write('You need PCMEpsilonStatic, PCMEpsilonDynamic and PCMDebyeRelaxTime for an EOM TD-PCM run.')
      call messages_new_line()
      call messages_write('Octopus will run using TD-PCM version in equilibrium with solvent at each time.')        
      call messages_warning(namespace=namespace)
      pcm%tdlevel = PCM_TD_EQ
    end if

    if (abs(pcm%epsilon_0 - M_ONE) <= M_EPSILON ) then
      if (pcm%tdlevel == PCM_TD_EOM .and. pcm%which_eps == PCM_DRUDE_MODEL) then
        message(1) = "PCMEpsilonStatic = 1 is incompatible with a Drude-Lorentz EOM-PCM run."
        call messages_fatal(1, namespace=namespace)
      end if
    else
      !%Variable PCMDrudeLOmega
      !%Type float
      !%Default <math>\sqrt{1/(\varepsilon_0-1)}</math>
      !%Section Hamiltonian::PCM
      !%Description
      !% Resonance frequency of the solvent within Drude-Lorentz model (<math>\omega_0</math>).
      !% Recall Drude-Lorentz dielectric function: <math>\varepsilon(\omega)=1+\frac{A}{\omega_0^2-\omega^2+i\gamma\omega}</math>   
      !% Default values of <math>\omega_0</math> guarantee to recover static dielectric constant.   
      !%End
      call parse_variable(namespace, 'PCMDrudeLOmega', sqrt(M_ONE/(pcm%epsilon_0 - M_ONE)), pcm%drl%w0)
      call messages_print_var_value(stdout, "PCMDrudeLOmega", pcm%drl%w0)
    end if    

    if (pcm%tdlevel == PCM_TD_EOM .and. pcm%which_eps == PCM_DRUDE_MODEL .and. pcm%drl%w0 == M_ZERO) then
      call messages_write('Sorry, you have set PCMDrudeLOmega = 0 but this is incompatible with a Drude-Lorentz EOM-PCM run.')
      call messages_new_line()
      if (pcm%epsilon_0 /= M_ONE) then
        call messages_write('Octopus will run using the default value of PCMDrudeLOmega.')        
        call messages_warning(namespace=namespace)
        pcm%drl%w0 = sqrt(M_ONE/(pcm%epsilon_0 - M_ONE))
      else
        message(1) = "PCMEpsilonStatic = 1 is incompatible with a Drude-Lorentz EOM-PCM run."
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    !%Variable PCMDrudeLDamping
    !%Type float
    !%Default 0.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Damping factor of the solvent charges oscillations within Drude-Lorentz model (<math>\gamma</math>).
    !% Recall Drude-Lorentz dielectric function: <math>\varepsilon(\omega)=1+\frac{A}{\omega_0^2-\omega^2+i\gamma\omega}</math>   
    !%End
    call parse_variable(namespace, 'PCMDrudeLDamping', M_ZERO, pcm%drl%gm)
    call messages_print_var_value(stdout, "PCMDrudeLDamping", pcm%drl%gm)

    !< Parameter (<math>A</math>) interpolating Drude-Lorentz dielectric function to its static value (<math>\varepsilon_0</math>).
    pcm%drl%aa = (pcm%epsilon_0 - M_ONE)*pcm%drl%w0**2

    !%Variable PCMLocalField
    !%Type logical
    !%Default no
    !%Section Hamiltonian::PCM
    !%Description
    !% This variable is a flag for including local field effects when an external field is applied. The total field interacting with
    !% the molecule (also known as cavity field) is not the bare field in the solvent (the so-called Maxwell field), but it also
    !% include a contribution due to the polarization of the solvent. The latter is calculated here within the PCM framework.
    !% See [G. Gil, et al., J. Chem. Theory Comput., 2019, 15 (4), pp 2306–2319].
    !%End
    call parse_variable(namespace, 'PCMLocalField', .false., pcm%localf)
    call messages_print_var_value(stdout, "PCMLocalField", pcm%localf)

    if (pcm%localf .and. ((.not.external_potentials_present) .and. (.not.pcm%kick_is_present))) then
      message(1) = "Sorry, you have set PCMLocalField = yes, but you have not included any external potentials."
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable PCMSolute
    !%Type logical
    !%Default yes
    !%Section Hamiltonian::PCM
    !%Description
    !% This variable is a flag for including polarization effects of the solvent due to the solute. 
    !% (Useful for analysis) When external fields are applied, turning off the solvent-molecule interaction (PCMSolute=no) and 
    !% activating the solvent polarization due to the applied field (PCMLocalField=yes) allows to include only local field effects. 
    !%End
    call parse_variable(namespace, 'PCMSolute', .true., pcm%solute)
    call messages_print_var_value(stdout, "PCMSolute", pcm%solute)

    if (pcm%run_pcm .and. (.not. pcm%solute)) then
      call messages_write('N.B. This PCM run do not consider the polarization effects due to the solute.')        
      call messages_warning(namespace=namespace)
      if (.not. pcm%localf) then
        message(1) = "You have activated a PCM run without polarization effects. Octopus will halt."
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    !%Variable PCMKick
    !%Type logical
    !%Default no
    !%Section Hamiltonian::PCM
    !%Description
    !% This variable controls the effect the kick has on the polarization of the solvent.
    !% If .true.  ONLY the FAST degrees-of-freedom of the solvent follow the kick. The potential due to polarization charges behaves 
    !%  as another kick, i.e., it is a delta-perturbation. 
    !% If .false. ALL           degrees-of-freedom of the solvent follow the kick. The potential due to polarization charges evolves 
    !%  following an equation of motion. When Debye dielectric model is used, just a part of the potential behaves as another kick.
    !%End
    call parse_variable(namespace, 'PCMKick', .false., pcm%kick_like)
    call messages_print_var_value(stdout, "PCMKick", pcm%kick_like)

    if (pcm%kick_like .and. (.not. pcm%run_pcm)) then
      message(1) = "PCMKick option can only be activated when PCMCalculation = yes. Octopus will halt."
      call messages_fatal(1, namespace=namespace)
    end if

    if (pcm%kick_like .and. (.not. pcm%localf)) then
      message(1) = "PCMKick option can only be activated when a PCMLocalField = yes. Octopus will halt."
      call messages_fatal(1, namespace=namespace)
    endif

    if (pcm%kick_like .and. (.not. pcm%kick_is_present)) then
      message(1) = "Sorry, you have set PCMKick = yes, but you have not included any kick."
      call messages_fatal(1, namespace=namespace)
    endif

    if (pcm%kick_is_present .and. pcm%run_pcm .and. (.not. pcm%localf)) then
      message(1) = "You have set up a PCM calculation with a kick without local field effects."
      message(2) = "Please, reconsider if you want the kick to be relevant for the PCM run."
      call messages_warning(2, namespace=namespace)
    end if
    
    !%Variable PCMUpdateIter
    !%Type integer
    !%Default 1
    !%Section Hamiltonian::PCM
    !%Description
    !% Defines how often the PCM potential is updated during time propagation.
    !%End
    call parse_variable(namespace, 'PCMUpdateIter', 1, pcm%update_iter)
    call messages_print_var_value(stdout, "PCMUpdateIter", pcm%update_iter)

    !%Variable PCMGamessBenchmark
    !%Type logical
    !%Default .false.
    !%Section Hamiltonian::PCM
    !%Description
    !% If PCMGamessBenchmark is set to "yes", the pcm_matrix is also written in a Gamess format.
    !% for benchamarking purposes.
    !%End
    call parse_variable(namespace, 'PCMGamessBenchmark', .false., gamess_benchmark)

    !%Variable PCMRenormCharges
    !%Type logical
    !%Default .false.
    !%Section Hamiltonian::PCM
    !%Description
    !% If .true. renormalization of the polarization charges is performed to enforce fulfillment
    !% of the Gauss law, <math>\sum_i q_i^{e/n} = -[(\epsilon-1)/\epsilon] Q_M^{e/n}</math> where 
    !% <math>q_i^{e/n}</math> are the polarization charges induced by the electrons/nuclei of the molecule
    !% and <math>Q_M^{e/n}</math> is the nominal electronic/nuclear charge of the system. This can be needed
    !% to treat molecules in weakly polar solvents.
    !%End
    call parse_variable(namespace, 'PCMRenormCharges', .false., pcm%renorm_charges)

    !%Variable PCMQtotTol
    !%Type float
    !%Default 0.5
    !%Section Hamiltonian::PCM
    !%Description
    !% If <tt>PCMRenormCharges=.true.</tt> and  <math>\delta Q = |[\sum_i q_i| - ((\epsilon-1)/\epsilon)*|Q_M]|>PCMQtotTol</math>
    !% the polarization charges will be normalized as 
    !% <math>q_i^\prime=q_i + signfunction(e, n, \delta Q) (q_i/q_{tot})*\delta Q</math>
    !% with <math>q_{tot} = \sum_i q_i</math>. For values of <math>\delta Q > 0.5</math>
    !% (printed by the code in the file pcm/pcm_info.out) even, if polarization charges are renormalized, 
    !% the calculated results might be inaccurate or erroneous.
    !%End
    call parse_variable(namespace, 'PCMQtotTol', CNST(0.5), pcm%q_tot_tol)

    if (pcm%renorm_charges) then
      message(1) = "Info: Polarization charges will be renormalized"
      message(2) = "      if |Q_tot_PCM - Q_M| > PCMQtotTol"
      call messages_info(2)
    end if

    !%Variable PCMSmearingFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Parameter used to control the width (area of each tessera) of the Gaussians used to distribute
    !% the polarization charges on each tessera (arXiv:1507.05471). If set to zero, the solvent 
    !% reaction potential in real-space is defined by using point charges.
    !%End
    call parse_variable(namespace, 'PCMSmearingFactor', M_ONE, pcm%gaussian_width)
    call messages_print_var_value(stdout, "PCMSmearingFactor", pcm%gaussian_width)

    if (abs(pcm%gaussian_width) <= M_EPSILON) then
      message(1) = "Info: PCM potential will be defined in terms of polarization point charges"
      call messages_info(1)
    else
      message(1) = "Info: PCM potential is regularized to avoid Coulomb singularity"
      call messages_info(1)        
    end if

    call io_mkdir('pcm', namespace)

    !%Variable PCMCavity
    !%Type string
    !%Section Hamiltonian::PCM
    !%Description
    !% Name of the file containing the geometry of the cavity hosting the solute molecule.
    !% The data must be in atomic units and the file must contain the following information sequentially:
    !%  T               < Number of tesserae
    !%  s_x(1:T)        < coordinates x of the tesserae 
    !%  s_y(1:T)        < coordinates y of the tesserae        
    !%  s_z(1:T)        < coordinates z of the tesserae            
    !%  A(1:T)          < areas of the tesserae
    !%  R_sph(1:T)      < Radii of the spheres to which the tesserae belong
    !%  normal(1:T,1:3) < Outgoing unitary vectors at the tesserae surfaces 
    !%End
    call parse_variable(namespace, 'PCMCavity', '', pcm%input_cavity)

    if (pcm%input_cavity == '') then

      !%Variable PCMSpheresOnH
      !%Type logical
      !%Default no
      !%Section Hamiltonian::PCM
      !%Description
      !% If true, spheres centered at the Hydrogens atoms are included to build the solute cavity surface.
      !%End
      call parse_variable(namespace, 'PCMSpheresOnH', .false., add_spheres_h)

      pcm%n_spheres = 0
      band = .false.
      do ia = 1, geo%natoms
        if ((.not. add_spheres_h) .and. geo%atom(ia)%label == 'H') cycle
        pcm%n_spheres = pcm%n_spheres + 1 !counting the number of species different from Hydrogen
      end do

      SAFE_ALLOCATE(pcm%spheres(1:pcm%n_spheres))
      pcm%spheres(:)%x = M_ZERO
      pcm%spheres(:)%y = M_ZERO
      pcm%spheres(:)%z = M_ZERO        
      pcm%spheres(:)%r = M_ZERO

      pcm%n_spheres = 0
      do ia = 1, geo%natoms
        if ((.not. add_spheres_h) .and. geo%atom(ia)%label == 'H') cycle
        pcm%n_spheres = pcm%n_spheres + 1

        !> These coordinates are already in atomic units (Bohr)
        pcm%spheres(pcm%n_spheres)%x = geo%atom(ia)%x(1)
        pcm%spheres(pcm%n_spheres)%y = geo%atom(ia)%x(2)
        pcm%spheres(pcm%n_spheres)%z = geo%atom(ia)%x(3)

        vdw_radius = pcm_get_vdw_radius(geo%atom(ia)%species, pcm_vdw_type, namespace)
        pcm%spheres(pcm%n_spheres)%r = vdw_radius*pcm%scale_r     
      end do

      if (mpi_grp_is_root(mpi_world)) then
        pcm%info_unit = io_open(PCM_DIR//'pcm_info.out', pcm%namespace, action='write')
      
        write(pcm%info_unit, '(A35)') '# Configuration: Molecule + Solvent'
        write(pcm%info_unit, '(A35)') '# ---------------------------------'
        write(pcm%info_unit, '(A21,F12.3)') '# Epsilon(Solvent) = ', pcm%epsilon_0
        write(pcm%info_unit, '(A1)')'#' 
        write(pcm%info_unit, '(A35,I4)') '# Number of interlocking spheres = ', pcm%n_spheres
        write(pcm%info_unit, '(A1)')'#'  
      
        write(pcm%info_unit, '(A8,3X,A7,8X,A26,20X,A10)') '# SPHERE', 'ELEMENT', 'CENTER  (X,Y,Z) (A)', 'RADIUS (A)'
        write(pcm%info_unit, '(A8,3X,A7,4X,A43,7X,A10)') '# ------', '-------', &
          '-------------------------------------------', '----------'  
      end if
      
      pcm%n_spheres = 0
      do ia = 1, geo%natoms
        if ((.not. add_spheres_h) .and. geo%atom(ia)%label == 'H') cycle
        pcm%n_spheres = pcm%n_spheres + 1
        if (mpi_grp_is_root(mpi_world)) & 
          write(pcm%info_unit,'(A1,2X,I3,7X,A2,3X,F14.8,2X,F14.8,2X,F14.8,4X,F14.8)')'#', pcm%n_spheres, &
          geo%atom(ia)%label,          &
          geo%atom(ia)%x(1)*P_a_B,     &
          geo%atom(ia)%x(2)*P_a_B,     &
          geo%atom(ia)%x(3)*P_a_B,     &
          pcm%spheres(pcm%n_spheres)%r*P_a_B
      end do

      !%Variable PCMTessSubdivider
      !%Type integer
      !%Default 1
      !%Section Hamiltonian::PCM
      !%Description
      !% Allows to subdivide further each tessera refining the discretization of the cavity tesselation. 
      !% Can take only two values, 1 or 4. 1 corresponds to 60 tesserae per sphere, while 4 corresponds to 240 tesserae per sphere.
      !%End
      call parse_variable(namespace, 'PCMTessSubdivider', 1, subdivider)

      SAFE_ALLOCATE(dum2(1:subdivider*N_TESS_SPHERE*pcm%n_spheres))

      !%Variable PCMTessMinDistance
      !%Type float
      !%Default 0.1
      !%Section Hamiltonian::PCM
      !%Description
      !% Minimum distance between tesserae. 
      !% Any two tesserae having smaller distance in the starting tesselation will be merged together. 
      !%End
      call parse_variable(namespace, 'PCMTessMinDistance', CNST(0.1), min_distance)
      call messages_print_var_value(stdout, "PCMTessMinDistance", min_distance)

      !> Counting the number of tesserae and generating the Van der Waals discretized surface of the solute system
      call cav_gen(subdivider, min_distance, pcm%n_spheres, pcm%spheres, pcm%n_tesserae, dum2, pcm%info_unit)

      SAFE_ALLOCATE(pcm%tess(1:pcm%n_tesserae))

      do ia=1, pcm%n_tesserae
        pcm%tess(ia)%point    = M_ZERO
        pcm%tess(ia)%area     = M_ZERO
        pcm%tess(ia)%r_sphere = M_ZERO
        pcm%tess(ia)%normal   = M_ZERO
      end do

      pcm%tess = dum2(1:pcm%n_tesserae)

      SAFE_DEALLOCATE_A(dum2)

      message(1) = "Info: van der Waals surface has been calculated"
      call messages_info(1)

    else

      !> The cavity surface will be read from a external file
      iunit = io_open(trim(pcm%input_cavity), pcm%namespace, action='read', status='old')
      read(iunit,*) pcm%n_tesserae

      if (pcm%n_tesserae > MXTS) then
        write(message(1),'(a,I5,a,I5)') "total number of tesserae", pcm%n_tesserae, ">", MXTS
        call messages_warning(1, namespace=namespace)
      end if

      SAFE_ALLOCATE(pcm%tess(1:pcm%n_tesserae))

      do ia = 1, pcm%n_tesserae
        pcm%tess(ia)%point    = M_ZERO
        pcm%tess(ia)%area     = M_ZERO
        pcm%tess(ia)%r_sphere = M_ZERO
        pcm%tess(ia)%normal   = M_ZERO
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%point(1)
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%point(2)
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%point(3)
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%area
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%r_sphere
      end do

      do ia = 1, pcm%n_tesserae
        read(iunit,*) pcm%tess(ia)%normal
      end do

      call io_close(iunit)
      message(1) = "Info: van der Waals surface has been read from " // trim(pcm%input_cavity)
      call messages_info(1)
    end if

    if (mpi_grp_is_root(mpi_world)) then
      cav_unit_test = io_open(PCM_DIR//'cavity_mol.xyz', pcm%namespace, action='write')

      write (cav_unit_test,'(2X,I4)') pcm%n_tesserae + geo%natoms
      write (cav_unit_test,'(2X)')

      do ia = 1, pcm%n_tesserae
        write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') 'H', pcm%tess(ia)%point*P_a_B
      end do

      do ia = 1, geo%natoms
        write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') geo%atom(ia)%label,      &
          geo%atom(ia)%x*P_a_B
      end do

      call io_close(cav_unit_test)

      write(pcm%info_unit,'(A1)')'#'
      if (pcm%localf) then
        write(pcm%info_unit,'(A1,4X,A4,14X,A4,21X,A4,21X,A4,21X,A4,21X,A7,18X,A7,18X,A8,17X,A5,20X,A8,17X,A5,20X,A8,17X,A5)') &
          '#','iter', 'E_ee', 'E_en', 'E_nn', 'E_ne', 'E_e_ext', 'E_n_ext', 'E_M-solv', &
          'Q_pol^e', 'deltaQ^e', 'Q_pol^n', 'deltaQ^n', 'Q_pol^ext'  
      else
        write(pcm%info_unit,'(A1,4X,A4,14X,A4,21X,A4,21X,A4,21X,A4,21X,A8,17X,A5,20X,A8,17X,A5,20X, A8)') &
          '#','iter', 'E_ee', 'E_en', 'E_nn', 'E_ne', 'E_M-solv', 'Q_pol^e', 'deltaQ^e', 'Q_pol^n', 'deltaQ^n'
      end if
    end if
    pcm%counter = 0

    !>printing out the cavity surface
    if (gamess_benchmark .and. mpi_grp_is_root(mpi_world)) then 
      cav_gamess_unit = io_open(PCM_DIR//'geom_cavity_gamess.out', pcm%namespace, action='write')

      write(cav_gamess_unit,*) pcm%n_tesserae

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%point(1)
      end do

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%point(2)
      end do

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%point(3)
      end do

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%area
      end do

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%r_sphere
      end do

      do ia = 1, pcm%n_tesserae
        write(cav_gamess_unit,*) pcm%tess(ia)%normal
      end do

      call io_close(cav_gamess_unit)
    end if

    !>Generating the dynamical PCM matrix
    if (gamess_benchmark) then
      SAFE_ALLOCATE(mat_gamess(1:pcm%n_tesserae, 1:pcm%n_tesserae) )
      mat_gamess = M_ZERO
    end if

    if (pcm%epsilon_infty /= pcm%epsilon_0) then
       
      SAFE_ALLOCATE(pcm%matrix_d(1:pcm%n_tesserae, 1:pcm%n_tesserae))
      pcm%matrix_d = M_ZERO

      call pcm_matrix(pcm%epsilon_infty, pcm%tess, pcm%n_tesserae, pcm%matrix_d)

      if (gamess_benchmark .and. mpi_grp_is_root(mpi_world)) then 
        pcmmat_gamess_unit = io_open(PCM_DIR//'pcm_matrix_gamess_dyn.out', pcm%namespace, action='write')

        do jtess = 1, pcm%n_tesserae
          do itess = 1, pcm%n_tesserae
            write(pcmmat_gamess_unit,*) mat_gamess(itess,jtess) !< for benchmarking with GAMESS
          end do
        end do

        call io_close(pcmmat_gamess_unit)
        mat_gamess = M_ZERO
      end if

      if (pcm%localf) then

        SAFE_ALLOCATE(pcm%matrix_lf_d(1:pcm%n_tesserae, 1:pcm%n_tesserae))
        pcm%matrix_lf_d = M_ZERO

        call pcm_matrix(pcm%epsilon_infty, pcm%tess, pcm%n_tesserae, pcm%matrix_lf_d, .true.)

        if (mpi_grp_is_root(mpi_world)) then
          ! only to benchmark
          pcmmat_unit = io_open(PCM_DIR//'pcm_matrix_dynamic_lf.out', pcm%namespace, action='write')
          do jtess = 1, pcm%n_tesserae
            do itess = 1, pcm%n_tesserae
              write(pcmmat_unit,*) pcm%matrix_lf_d(itess,jtess)
            end do
          end do
          call io_close(pcmmat_unit)
        end if

      end if

    end if

    SAFE_ALLOCATE(pcm%matrix(1:pcm%n_tesserae, 1:pcm%n_tesserae))
    pcm%matrix = M_ZERO

    call pcm_matrix(pcm%epsilon_0, pcm%tess, pcm%n_tesserae, pcm%matrix) 

    if (mpi_grp_is_root(mpi_world)) then
      pcmmat_unit = io_open(PCM_DIR//'pcm_matrix.out', pcm%namespace, action='write')
      if (gamess_benchmark) pcmmat_gamess_unit = io_open(PCM_DIR//'pcm_matrix_gamess.out', &
        pcm%namespace, action='write')

      do jtess = 1, pcm%n_tesserae
        do itess = 1, pcm%n_tesserae
          write(pcmmat_unit,*) pcm%matrix(itess,jtess)
          if (gamess_benchmark) write(pcmmat_gamess_unit,*) mat_gamess(itess,jtess) !< for benchmarking with GAMESS
        end do
      end do
      call io_close(pcmmat_unit)
      if (gamess_benchmark) call io_close(pcmmat_gamess_unit)
    end if
#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
    if (gamess_benchmark) then 
      SAFE_DEALLOCATE_A(mat_gamess)
    end if

    if (pcm%localf) then !< static local field case

      SAFE_ALLOCATE(pcm%matrix_lf(1:pcm%n_tesserae, 1:pcm%n_tesserae))
      pcm%matrix_lf = M_ZERO

      call pcm_matrix(pcm%epsilon_0, pcm%tess, pcm%n_tesserae, pcm%matrix_lf, .true.)

      if (mpi_grp_is_root(mpi_world)) then
        ! only to benchmark
        pcmmat_unit = io_open(PCM_DIR//'pcm_matrix_static_lf.out', pcm%namespace, action='write')
        do jtess = 1, pcm%n_tesserae
          do itess = 1, pcm%n_tesserae
            write(pcmmat_unit,*) pcm%matrix_lf(itess,jtess)
          end do
        end do
        call io_close(pcmmat_unit)
      end if

    end if

    message(1) = "Info: PCM response matrices has been evaluated"
    call messages_info(1)
    
    !%Variable PCMCalcMethod
    !%Type integer
    !%Default pcm_direct
    !%Section Hamiltonian::PCM
    !%Description
    !% Defines the method to be used to obtain the PCM potential.
    !%Option pcm_direct 1
    !% Direct sum of the potential generated by the polarization charges regularized 
    !% with a Gaussian smearing [A. Delgado, et al., J Chem Phys 143, 144111 (2015)].
    !%Option pcm_poisson 2
    !% Solving the Poisson equation for the polarization charge density.
    !%End
    call parse_variable(namespace, 'PCMCalcMethod', PCM_CALC_DIRECT, pcm%calc_method)
    call messages_print_var_option(stdout, "PCMCalcMethod", pcm%calc_method)


    if (pcm%calc_method == PCM_CALC_POISSON) then

      max_area = M_EPSILON
      do ia = 1, pcm%n_tesserae
        if (pcm%tess(ia)%area > max_area) max_area = pcm%tess(ia)%area
      end do
    
      !default is as many neighbor to contain 1 gaussian width 
      default_nn = int(max_area*pcm%gaussian_width/minval(grid%mesh%spacing(1:grid%mesh%sb%dim)))
      
      changed_default_nn = .false.

      do ii = default_nn, 1, -1
        pcm%tess_nn = ii 
        if (pcm_nn_in_mesh(pcm,grid%mesh)) then 
          exit
        else
          changed_default_nn = .true.
        end if
      end do
      if (changed_default_nn) then
        call messages_write('PCM nearest neighbors have been reduced from ')
        call messages_write(default_nn)
        call messages_write(' to ')
        call messages_write(pcm%tess_nn)
        call messages_new_line()
        call messages_write('in order to fit them into the mesh.')
        call messages_new_line()
        call messages_write('This may produce unexpected results. ')
        call messages_warning(namespace=namespace)
      end if
      default_nn = pcm%tess_nn

      !%Variable PCMChargeSmearNN
      !%Type integer
      !%Default 2 * max_area * PCMSmearingFactor
      !%Section Hamiltonian::PCM
      !%Description
      !% Defines the number of nearest neighbor mesh-points to be taken around each 
      !% cavity tessera in order to smear the charge when PCMCalcMethod = pcm_poisson.
      !% Setting PCMChargeSmearNN = 1 means first nearest neighbors, PCMChargeSmearNN = 2
      !% second nearest neighbors, and so on.
      !% The default value is such that the neighbor mesh contains points in a radius 
      !% equal to the width used for the gaussian smearing.
      !%End
      call parse_variable(namespace, 'PCMChargeSmearNN', default_nn, pcm%tess_nn)
      call messages_print_var_value(stdout, "PCMChargeSmearNN", pcm%tess_nn)
      
      call pcm_poisson_sanity_check(pcm, grid%mesh)
      
    end if
    
    if (pcm%run_pcm) call messages_print_stress(stdout, namespace=namespace)

    if (pcm%calc_method == PCM_CALC_POISSON) then
      SAFE_ALLOCATE(pcm%rho_n(1:grid%mesh%np_part))
      SAFE_ALLOCATE(pcm%rho_e(1:grid%mesh%np_part))
      if (pcm%localf) then
        SAFE_ALLOCATE(pcm%rho_ext(1:grid%mesh%np_part))
        if( pcm%kick_is_present ) SAFE_ALLOCATE(pcm%rho_kick(1:grid%mesh%np_part))
      end if
    end if 


    SAFE_ALLOCATE(pcm%v_n(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%q_n(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%v_n_rs(1:grid%mesh%np))
    pcm%v_n    = M_ZERO
    pcm%q_n    = M_ZERO
    pcm%v_n_rs = M_ZERO

    SAFE_ALLOCATE(pcm%v_e(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%q_e(1:pcm%n_tesserae))
    SAFE_ALLOCATE(pcm%v_e_rs(1:grid%mesh%np))
    pcm%v_e    = M_ZERO
    pcm%q_e    = M_ZERO
    pcm%v_e_rs = M_ZERO

    SAFE_ALLOCATE(pcm%q_e_in(1:pcm%n_tesserae))
    pcm%q_e_in = M_ZERO


    if (pcm%localf) then
      SAFE_ALLOCATE(pcm%v_ext(1:pcm%n_tesserae))
      SAFE_ALLOCATE(pcm%q_ext(1:pcm%n_tesserae))
      SAFE_ALLOCATE(pcm%v_ext_rs(1:grid%mesh%np))
      pcm%v_ext    = M_ZERO
      pcm%q_ext    = M_ZERO
      pcm%v_ext_rs = M_ZERO

      SAFE_ALLOCATE(pcm%q_ext_in(1:pcm%n_tesserae))
      pcm%q_ext_in = M_ZERO

      if (pcm%kick_is_present) then
        SAFE_ALLOCATE(pcm%v_kick(1:pcm%n_tesserae))
        SAFE_ALLOCATE(pcm%q_kick(1:pcm%n_tesserae))
        SAFE_ALLOCATE(pcm%v_kick_rs(1:grid%mesh%np))
        pcm%v_ext    = M_ZERO
        pcm%q_ext    = M_ZERO
        pcm%v_ext_rs = M_ZERO
      end if 

    end if

    pcm%q_e_nominal = qtot
    pcm%q_n_nominal = val_charge
    pcm%deltaQ_e = M_ZERO
    pcm%deltaQ_n = M_ZERO
    pcm%qtot_e = M_ZERO
    pcm%qtot_n = M_ZERO
    pcm%qtot_ext = M_ZERO

    POP_SUB(pcm_init)
  end subroutine pcm_init

  ! -----------------------------------------------------------------------------
  subroutine pcm_calc_pot_rs(pcm, mesh, psolver, geo, v_h, v_ext, kick, time_present, kick_time)
    save
    type(pcm_t),                intent(inout) :: pcm
    type(mesh_t),               intent(in)    :: mesh
    type(poisson_t),            intent(in)    :: psolver
    type(geometry_t), optional, intent(in)    :: geo
    FLOAT,            optional, intent(in)    :: v_h(:)
    FLOAT,            optional, intent(in)    :: v_ext(:)
    FLOAT,            optional, intent(in)    :: kick(:)
    logical,          optional, intent(in)    :: time_present
    logical,          optional, intent(in)    :: kick_time
    
    integer :: calc

    logical :: input_asc_e
    logical :: input_asc_ext

    ! for debuggin - it should be cleaned up
    !character*5 :: iteration 

    logical :: not_yet_called = .false.
    logical :: is_time_for_kick = .false.
    logical :: after_kick = .false.

    logical :: td_calc_mode = .false.

    integer :: asc_vs_t_unit_e

    character(len=23) :: asc_vs_t_unit_format
    character(len=16) :: asc_vs_t_unit_format_tail

    integer :: ia

    PUSH_SUB(pcm_calc_pot_rs)  
    
    ASSERT(present(v_h) .or. present(geo) .or. present(v_ext) .or. present(kick) )

    if (present(v_h))                               calc = PCM_ELECTRONS
    if (present(geo))                               calc = PCM_NUCLEI
    if (present(v_ext) .and. (.not.present(kick)))  calc = PCM_EXTERNAL_POTENTIAL
    if ((.not. present(v_ext)) .and. present(kick)) calc = PCM_KICK
    if (present(v_ext) .and. present(kick))         calc = PCM_EXTERNAL_PLUS_KICK

    if (present(time_present)) then
      if (time_present) td_calc_mode = .true.
    end if

    if (present(kick_time)) then
      is_time_for_kick = kick_time
      if (kick_time) after_kick = .true.
    end if

    select case(pcm%initial_asc)
    case(1)
      input_asc_e = .true.
      input_asc_ext = .false.
    case(2)
      input_asc_e = .false.
      input_asc_ext = .true.
    case(3)
      input_asc_e = .true.
      input_asc_ext = .true.
    case default
      input_asc_e = .false.
      input_asc_ext = .false.
    end select
     
    if (calc == PCM_NUCLEI) then
      call pcm_v_nuclei_cav(pcm%v_n, geo, pcm%tess, pcm%n_tesserae)
      call pcm_charges(pcm%q_n, pcm%qtot_n, pcm%v_n, pcm%matrix, pcm%n_tesserae, &
                       pcm%q_n_nominal, pcm%epsilon_0, pcm%renorm_charges, pcm%q_tot_tol, pcm%deltaQ_n)
      if (pcm%calc_method == PCM_CALC_POISSON) call pcm_charge_density(pcm, pcm%q_n, pcm%qtot_n, mesh, pcm%rho_n)
      call pcm_pot_rs(pcm, pcm%v_n_rs, pcm%q_n, pcm%rho_n, mesh, psolver)
    end if

    if (calc == PCM_ELECTRONS) then
      call pcm_v_cav_li(pcm%v_e, v_h, pcm, mesh)
      if (td_calc_mode .and. pcm%epsilon_infty /= pcm%epsilon_0 .and. pcm%tdlevel /= PCM_TD_EQ) then
        !< BEGIN - equation-of-motion propagation or inertial/dynamical charge splitting
        select case (pcm%tdlevel)
        case (PCM_TD_EOM) !< equation-of-motion propagation
          select case (pcm%which_eps)
          case (PCM_DRUDE_MODEL)
            call pcm_charges_propagation(pcm%q_e, pcm%v_e, pcm%dt, pcm%tess, input_asc_e, calc, &
              PCM_DRUDE_MODEL, pcm%namespace, this_drl = pcm%drl)
          case default
            call pcm_charges_propagation(pcm%q_e, pcm%v_e, pcm%dt, pcm%tess, input_asc_e, calc, &
              PCM_DEBYE_MODEL, pcm%namespace, this_deb = pcm%deb)
          end select
          if ((.not. pcm%localf) .and. not_yet_called) call pcm_eom_enough_initial(not_yet_called)
          !< total pcm charges due to solute electrons
          pcm%qtot_e = sum(pcm%q_e)
        case (PCM_TD_NEQ) !< inertial/dynamical partition
          select case (pcm%iter)
          case(1)
            !< calculating inertial polarization charges (once and for all)
            !< calculated in two steps to guarantee correct renormalization: Gauss law is recovered at each step.
            !< (first step) calculating polarization charges equilibrated with the initial Hartree potential (just once)
            call pcm_charges(pcm%q_e, pcm%qtot_e, pcm%v_e, pcm%matrix, pcm%n_tesserae, &
                             pcm%q_e_nominal, pcm%epsilon_0, pcm%renorm_charges, pcm%q_tot_tol, pcm%deltaQ_e)
	          !< (second step) calculating initial dynamic polarization charges 
	          !< dont pay attention to the use of q_e_in and qtot_e_in, whose role here is only auxiliary
            call pcm_charges(pcm%q_e_in, pcm%qtot_e_in, pcm%v_e, pcm%matrix_d, pcm%n_tesserae, &
                             pcm%q_e_nominal, pcm%epsilon_infty, pcm%renorm_charges, pcm%q_tot_tol, pcm%deltaQ_e)
            !< (finally) the inertial polarization charges
            pcm%q_e_in = pcm%q_e - pcm%q_e_in
            pcm%qtot_e_in = pcm%qtot_e - pcm%qtot_e_in
          case default
            !< calculating dynamical polarization charges (each time)
            call pcm_charges(pcm%q_e, pcm%qtot_e, pcm%v_e, pcm%matrix_d, pcm%n_tesserae, &
                             pcm%q_e_nominal, pcm%epsilon_infty, pcm%renorm_charges, pcm%q_tot_tol, pcm%deltaQ_e)

            pcm%q_e = pcm%q_e + pcm%q_e_in
            pcm%qtot_e = pcm%qtot_e + pcm%qtot_e_in
          end select
        end select !< END - equation-of-motion propagation or inertial/dynamical charge splitting
      else !< BEGIN - pcm charges propagation in equilibrium with solute

        !< calculating polarization charges to be equil. (upon self-consitency) with the ground state Hartree potential (always).
        call pcm_charges(pcm%q_e, pcm%qtot_e, pcm%v_e, pcm%matrix, pcm%n_tesserae, &
                         pcm%q_e_nominal, pcm%epsilon_0, pcm%renorm_charges, pcm%q_tot_tol, pcm%deltaQ_e)

      end if !< END - pcm charges propagation in equilibrium with solute
      if (pcm%calc_method == PCM_CALC_POISSON) call pcm_charge_density(pcm, pcm%q_e, pcm%qtot_e, mesh, pcm%rho_e)
      call pcm_pot_rs(pcm, pcm%v_e_rs, pcm%q_e, pcm%rho_e, mesh, psolver)
    end if

    if ((calc == PCM_EXTERNAL_PLUS_KICK .or. calc == PCM_KICK) .and. .not. pcm%kick_like) then
      pcm%q_ext    = M_ZERO
      pcm%v_ext_rs = M_ZERO
    end if

    if (calc == PCM_EXTERNAL_POTENTIAL .or. calc == PCM_EXTERNAL_PLUS_KICK) then
      call pcm_v_cav_li(pcm%v_ext, v_ext, pcm, mesh)
      if (td_calc_mode .and. pcm%epsilon_infty /= pcm%epsilon_0 .and. pcm%tdlevel /= PCM_TD_EQ) then
        !< BEGIN - equation-of-motion propagation or inertial/dynamical charge splitting
        select case (pcm%tdlevel)
        case (PCM_TD_EOM) !< equation-of-motion propagation
          select case (pcm%which_eps)
          case(PCM_DRUDE_MODEL)
            call pcm_charges_propagation(pcm%q_ext, pcm%v_ext, pcm%dt, pcm%tess, input_asc_ext, calc, &
                                         pcm%which_eps, pcm%namespace, this_drl=pcm%drl)
          case default
            call pcm_charges_propagation(pcm%q_ext, pcm%v_ext, pcm%dt, pcm%tess, input_asc_ext, calc, &
                                         pcm%which_eps, pcm%namespace, this_deb=pcm%deb)
          end select
          if ((.not. pcm%kick_is_present) .and. not_yet_called) call pcm_eom_enough_initial(not_yet_called)
          !< total pcm charges due to time-dependent external field
          pcm%qtot_ext = sum(pcm%q_ext)
          !< summing pcm charges due to time-dependent and static external potentials
          !< dont pay attention to the use of q_ext_in and qtot_ext_in, whose role here is only auxiliary
          pcm%q_ext = pcm%q_ext + pcm%q_ext_in
          pcm%qtot_ext = pcm%qtot_ext + pcm%qtot_ext_in
        case (PCM_TD_NEQ) !< inertial/dynamical partition -> sudden switch-on
          !< calculating dynamical polarization charges (each time)
          call pcm_charges(pcm%q_ext, pcm%qtot_ext, pcm%v_ext, pcm%matrix_lf_d, pcm%n_tesserae)

          !< summing pcm charges due to time-dependent and static external potentials
          !< dont pay attention to the use of q_ext_in and qtot_ext_in, whose role here is only auxiliary
          pcm%q_ext = pcm%q_ext + pcm%q_ext_in
          pcm%qtot_ext = pcm%qtot_ext + pcm%qtot_ext_in
        end select !< END - equation-of-motion propagation or inertial/dynamical charge splitting
      else !< BEGIN - pcm charges propagation in equilibrium with external potential

        !< calculating polarization charges equilibrated with the external potential (always).
        !< if there are (together) time-dependent and static external potentials, 
        !< the pcm charges due to the static one is calculated here.
        call pcm_charges(pcm%q_ext, pcm%qtot_ext, pcm%v_ext, pcm%matrix_lf, pcm%n_tesserae)

        !< dont pay attention to the use of q_ext_in and qtot_ext_in, whose role here is only auxiliary
        pcm%q_ext_in = pcm%q_ext
        pcm%qtot_ext_in = pcm%qtot_ext
      end if !< END - pcm charges propagation in equilibrium with external field
      if (pcm%calc_method == PCM_CALC_POISSON) call pcm_charge_density(pcm, pcm%q_ext, pcm%qtot_ext, mesh, pcm%rho_ext)
      call pcm_pot_rs(pcm, pcm%v_ext_rs, pcm%q_ext, pcm%rho_ext, mesh, psolver)

    end if

    if (calc == PCM_EXTERNAL_PLUS_KICK .or. calc == PCM_KICK) then
      pcm%q_kick = M_ZERO
      pcm%v_kick_rs = M_ZERO
      if ( is_time_for_kick ) then
        call pcm_v_cav_li(pcm%v_kick, kick, pcm, mesh)
      else
        pcm%v_kick = M_ZERO
      end if
      if (pcm%kick_like) then
        if (is_time_for_kick) then
          !< kick-like polarization charges
          if (pcm%epsilon_infty /= pcm%epsilon_0) then
            call pcm_charges(pcm%q_kick, pcm%qtot_kick, pcm%v_kick, pcm%matrix_lf_d, pcm%n_tesserae)
          else
            call pcm_charges(pcm%q_kick, pcm%qtot_kick, pcm%v_kick, pcm%matrix_lf, pcm%n_tesserae)
         end if
       else
         pcm%q_kick = M_ZERO
         pcm%qtot_kick = M_ZERO
         if (pcm%calc_method == PCM_CALC_POISSON) pcm%rho_kick = M_ZERO
         pcm%v_kick_rs = M_ZERO

         POP_SUB(pcm_calc_pot_rs)
         return
        end if
      else
        if (after_kick) then
          !< BEGIN - equation-of-motion propagation
          select case (pcm%which_eps)
          case(PCM_DRUDE_MODEL)
            call pcm_charges_propagation(pcm%q_kick, pcm%v_kick, pcm%dt, pcm%tess, input_asc_ext, calc, &
                                         pcm%which_eps, pcm%namespace, this_drl=pcm%drl)
          case default
            call pcm_charges_propagation(pcm%q_kick, pcm%v_kick, pcm%dt, pcm%tess, input_asc_ext, calc, &
                                         pcm%which_eps, pcm%namespace, this_deb=pcm%deb)
          end select
          if (not_yet_called) call pcm_eom_enough_initial(not_yet_called)
          pcm%qtot_kick  = sum(pcm%q_kick)
          !< END - equation-of-motion propagation
        else
          pcm%q_kick = M_ZERO
          pcm%qtot_kick = M_ZERO
          if (pcm%calc_method == PCM_CALC_POISSON) pcm%rho_kick = M_ZERO
          pcm%v_kick_rs = M_ZERO

          POP_SUB(pcm_calc_pot_rs)
          return
        end if
      end if

      if (pcm%calc_method == PCM_CALC_POISSON) call pcm_charge_density(pcm, pcm%q_kick, pcm%qtot_kick, mesh, pcm%rho_kick)
      call pcm_pot_rs(pcm, pcm%v_kick_rs, pcm%q_kick, pcm%rho_kick, mesh, psolver)
      if (.not. pcm%kick_like) then
        if (calc == PCM_EXTERNAL_PLUS_KICK) then
          pcm%q_ext    = pcm%q_ext    + pcm%q_kick
          pcm%v_ext_rs = pcm%v_ext_rs + pcm%v_kick_rs
        else
          pcm%q_ext    = pcm%q_kick
          pcm%v_ext_rs = pcm%v_kick_rs
        end if
      end if

    end if


    if (mpi_grp_is_root(mpi_world)) then
      write (asc_vs_t_unit_format_tail,'(I5,A11)') pcm%n_tesserae,'(1X,F14.8))'
      write (asc_vs_t_unit_format,'(A)') '(F14.8,'//trim(adjustl(asc_vs_t_unit_format_tail))
  
      if (pcm%solute .and. pcm%localf .and. td_calc_mode .and. calc == PCM_ELECTRONS ) then
        asc_vs_t_unit_e = io_open(PCM_DIR//'ASC_e_vs_t.dat', pcm%namespace, &
          action='write', position='append', form='formatted')
        write(asc_vs_t_unit_e,trim(adjustl(asc_vs_t_unit_format))) pcm%iter*pcm%dt, &
          (pcm%q_e(ia) , ia = 1,pcm%n_tesserae)
        call io_close(asc_vs_t_unit_e)
      end if
    end if

    POP_SUB(pcm_calc_pot_rs)  
  end subroutine pcm_calc_pot_rs

  ! -----------------------------------------------------------------------------
  !> Calculates the Hartree/external/kick potential at the tessera representative points by doing 
  !! a 3D linear interpolation. 
  subroutine pcm_v_cav_li(v_cav, v_mesh, pcm, mesh)
    type(pcm_t), intent(in)  :: pcm
    type(mesh_t), intent(in) :: mesh
    FLOAT, intent(in)        :: v_mesh(:) !< (1:mesh%np)
    FLOAT, intent(out)       :: v_cav(:)   !< (1:n_tess)

    integer :: ia

    type(mesh_interpolation_t)  :: mesh_interpolation

    PUSH_SUB(pcm_v_cav_li)    

    v_cav = M_ZERO
    
    call mesh_interpolation_init(mesh_interpolation, mesh)

    do ia = 1, pcm%n_tesserae
      call mesh_interpolation_evaluate(mesh_interpolation, v_mesh, pcm%tess(ia)%point, v_cav(ia))
    end do

    call mesh_interpolation_end(mesh_interpolation)

    POP_SUB(pcm_v_cav_li)
  end subroutine pcm_v_cav_li

  !--------------------------------------------------------------------------------
  
  !> Calculates the classical electrostatic potential geneated by the nuclei at the tesserae.
  !! v_n_cav(ik) = \sum_{I=1}^{natoms} Z_val / |s_{ik} - R_I|
  subroutine pcm_v_nuclei_cav(v_n_cav, geo, tess, n_tess)
    FLOAT,               intent(out) :: v_n_cav(:) !< (1:n_tess)
    type(geometry_t),    intent(in)  :: geo
    type(pcm_tessera_t), intent(in)  :: tess(:)    !< (1:n_tess)
    integer,             intent(in)  :: n_tess

    FLOAT   :: diff(1:PCM_DIM_SPACE), dist, z_ia
    integer :: ik, ia

    PUSH_SUB(pcm_v_nuclei_cav)

    v_n_cav = M_ZERO

    do ik = 1, n_tess
      do ia = 1, geo%natoms

        diff = geo%atom(ia)%x(1:PCM_DIM_SPACE) - tess(ik)%point
        dist = sqrt(dot_product(diff, diff))

        z_ia = species_zval(geo%atom(ia)%species)

        v_n_cav(ik) = v_n_cav(ik) - z_ia/dist
      end do
    end do

    POP_SUB(pcm_v_nuclei_cav)
  end subroutine pcm_v_nuclei_cav

  ! -----------------------------------------------------------------------------
  
  !> Calculates the solute-solvent electrostatic interaction energy
  !! E_M-solv = \sum{ik=1}^n_tess { [VHartree(ik) + Vnuclei(ik)]*[q_e(ik) + q_n(ik)] }
  !!            (if external potential)                                   + q_ext(ik)   
  subroutine pcm_elect_energy(geo, pcm, E_int_ee, E_int_en, E_int_ne, E_int_nn, E_int_e_ext, E_int_n_ext)
    type(geometry_t), intent(in)  :: geo
    type(pcm_t),      intent(in)  :: pcm
    FLOAT,            intent(out) :: E_int_ee 
    FLOAT,            intent(out) :: E_int_en 
    FLOAT,            intent(out) :: E_int_ne 
    FLOAT,            intent(out) :: E_int_nn
    FLOAT, optional,  intent(out) :: E_int_e_ext
    FLOAT, optional,  intent(out) :: E_int_n_ext 

    FLOAT   :: diff(1:PCM_DIM_SPACE)
    FLOAT   :: dist, z_ia
    integer :: ik, ia
    type(species_t), pointer :: spci 

    PUSH_SUB(pcm_elect_energy)

    E_int_ee = M_ZERO
    E_int_en = M_ZERO
    E_int_ne = M_ZERO
    E_int_nn = M_ZERO
    
    if (pcm%localf .and. ( (.not.present(E_int_e_ext)) .or. &
                           (.not.present(E_int_n_ext))      ) ) then
      message(1) = "pcm_elect_energy: There are lacking terms in subroutine call."
      call messages_fatal(1, namespace=pcm%namespace)
    else if (pcm%localf .and. ( present(E_int_e_ext) .and. &
                                present(E_int_n_ext)       ) ) then
      E_int_e_ext = M_ZERO
      E_int_n_ext = M_ZERO
    end if

    diff = M_ZERO

    do ik = 1, pcm%n_tesserae

      E_int_ee = E_int_ee + pcm%v_e(ik)*pcm%q_e(ik)
      E_int_en = E_int_en + pcm%v_e(ik)*pcm%q_n(ik)
      if (pcm%localf) &
        E_int_e_ext = E_int_e_ext + pcm%v_e(ik)*pcm%q_ext(ik)

      do ia = 1, geo%natoms

        diff = geo%atom(ia)%x(1:PCM_DIM_SPACE) - pcm%tess(ik)%point
        dist = dot_product(diff, diff)
        dist = sqrt(dist)

        spci => geo%atom(ia)%species
        z_ia = -species_zval(spci)

        E_int_ne = E_int_ne + z_ia*pcm%q_e(ik) / dist
        E_int_nn = E_int_nn + z_ia*pcm%q_n(ik) / dist
        if (pcm%localf) E_int_n_ext = E_int_n_ext + z_ia*pcm%q_ext(ik) / dist

      end do
    end do

    E_int_ee = M_HALF*E_int_ee
    E_int_en = M_HALF*E_int_en
    E_int_ne = M_HALF*E_int_ne
    E_int_nn = M_HALF*E_int_nn
    ! E_int_e_ext and E_int_n_ext do not need 1/2 factor

    ! print results of the iteration in pcm_info file
    if ( mpi_grp_is_root(mpi_world) ) then
      if (pcm%localf) then
        write(pcm%info_unit, &
 '(3X,I5,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X)') &
              pcm%counter, &
              units_from_atomic(units_out%energy, E_int_ee), & 
              units_from_atomic(units_out%energy, E_int_en), &
              units_from_atomic(units_out%energy, E_int_nn), &
              units_from_atomic(units_out%energy, E_int_ne), &
              units_from_atomic(units_out%energy, E_int_e_ext), &
              units_from_atomic(units_out%energy, E_int_n_ext), &
              units_from_atomic(units_out%energy, E_int_ee + E_int_en + E_int_nn + E_int_ne + E_int_e_ext + E_int_n_ext), &
              pcm%qtot_e, &
              pcm%deltaQ_e, &
              pcm%qtot_n, &
              pcm%deltaQ_n, &
              pcm%qtot_ext
      else
        write(pcm%info_unit, &
              '(3X,I5,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8,5X,F20.8)') &
              pcm%counter, &
              units_from_atomic(units_out%energy, E_int_ee), & 
              units_from_atomic(units_out%energy, E_int_en), &
              units_from_atomic(units_out%energy, E_int_nn), &
              units_from_atomic(units_out%energy, E_int_ne), &
              units_from_atomic(units_out%energy, E_int_ee +  E_int_en + E_int_nn + E_int_ne), &
              pcm%qtot_e, &
              pcm%deltaQ_e, &
              pcm%qtot_n, &
              pcm%deltaQ_n
      end if
    end if

    POP_SUB(pcm_elect_energy)
  end subroutine pcm_elect_energy

  ! -----------------------------------------------------------------------------
  
  !> Calculates the polarization charges at each tessera by using the response matrix 'pcm_mat',
  !! provided the value of the molecular electrostatic potential at 
  !! the tesserae: q_pcm(ia) = \sum_{ib}^{n_tess} pcm_mat(ia,ib)*v_cav(ib).
  subroutine pcm_charges(q_pcm, q_pcm_tot, v_cav, pcm_mat, n_tess, qtot_nominal, epsilon, renorm_charges, q_tot_tol, deltaQ)
    FLOAT,             intent(out) :: q_pcm(:)     !< (1:n_tess)
    FLOAT,             intent(out) :: q_pcm_tot
    FLOAT,             intent(in)  :: v_cav(:)     !< (1:n_tess)
    FLOAT,             intent(in)  :: pcm_mat(:,:) !< (1:n_tess, 1:n_tess)
    integer,           intent(in)  :: n_tess
    FLOAT,   optional, intent(in)  :: qtot_nominal
    FLOAT,   optional, intent(in)  :: epsilon
    logical, optional, intent(in)  :: renorm_charges
    FLOAT,   optional, intent(in)  :: q_tot_tol
    FLOAT,   optional, intent(out) :: deltaQ

    integer :: ia, ib
    type(profile_t), save :: prof_init
    FLOAT :: q_pcm_tot_norm
    FLOAT :: coeff

    PUSH_SUB(pcm_charges)
    call profiling_in(prof_init, 'PCM_CHARGES') 

    q_pcm     = M_ZERO
    q_pcm_tot = M_ZERO

    do ia = 1, n_tess
      do ib = 1, n_tess
        q_pcm(ia) = q_pcm(ia) + pcm_mat(ia,ib)*v_cav(ib) !< transpose matrix might speed up
      end do
      q_pcm_tot = q_pcm_tot + q_pcm(ia)
    end do

    if (present(qtot_nominal)   .and. present(epsilon) .and.  &
        present(renorm_charges) .and. present(q_tot_tol) .and. present(deltaQ) ) then
      !< renormalization of the polarization charges to enforce fulfillment of the Gauss law.
      deltaQ = abs(q_pcm_tot) - ((epsilon - M_ONE)/epsilon )*abs(qtot_nominal)
      if ((renorm_charges) .and. (abs(deltaQ) > q_tot_tol)) then
        q_pcm_tot_norm = M_ZERO
        coeff = sign(M_ONE, qtot_nominal)*sign(M_ONE, deltaQ)
        do ia = 1, n_tess
          q_pcm(ia) = q_pcm(ia) + coeff*q_pcm(ia)/q_pcm_tot*abs(deltaQ)
          q_pcm_tot_norm = q_pcm_tot_norm + q_pcm(ia)
        end do
        q_pcm_tot = q_pcm_tot_norm
      end if
    end if

    call profiling_out(prof_init)
    POP_SUB(pcm_charges)
  end subroutine pcm_charges

  ! -----------------------------------------------------------------------------    
  !> Check wether the nearest neighbor requested are in the mesh or not
  logical function pcm_nn_in_mesh(pcm, mesh) result(in_mesh)
    type(pcm_t),  intent(in) :: pcm 
    type(mesh_t), intent(in) :: mesh

    integer :: ia, nm(1:MAX_DIM), ipt, i1, i2, i3
    FLOAT :: posrel(1:MAX_DIM)
    integer :: pt
      
    PUSH_SUB(pcm_nn_in_mesh)
  
    in_mesh = .true.
    do ia = 1, pcm%n_tesserae
  
      posrel(1:mesh%sb%dim) = pcm%tess(ia)%point(1:mesh%sb%dim)/mesh%spacing(1:mesh%sb%dim)

      nm(1:mesh%sb%dim) = floor(posrel(1:mesh%sb%dim))
  
      ! Get the nearest neighboring points
      ipt = 0
      do i1 = -pcm%tess_nn + 1 , pcm%tess_nn
        do i2 = -pcm%tess_nn + 1 , pcm%tess_nn
          do i3 = -pcm%tess_nn + 1 , pcm%tess_nn
            ipt = ipt + 1
            pt = index_from_coords(mesh%idx, [i1 + nm(1), i2 + nm(2), i3 + nm(3)])

            if (pt <= 0 .or. pt > mesh%np_part_global) then 
              in_mesh = .false.
              POP_SUB(pcm_nn_in_mesh)
              return 
            end if
          
          end do
        end do
      end do
    end do

    POP_SUB(pcm_nn_in_mesh)
  end function pcm_nn_in_mesh
  
  ! -----------------------------------------------------------------------------  
  !> Check that all the required nearest neighbors are prensent in the mesh
  subroutine pcm_poisson_sanity_check(pcm, mesh)
    type(pcm_t),     intent(in) :: pcm 
    type(mesh_t),    intent(in) :: mesh
    
    PUSH_SUB(pcm_poisson_sanity_check)

    if (.not. pcm_nn_in_mesh(pcm, mesh)) then 
      message(1) = 'The simulation box is too small to contain all the requested'
      message(2) = 'nearest neighbors for each tessera.'
      message(3) = 'Consider using a larger box or reduce PCMChargeSmearNN.'
      call messages_warning(3, namespace=pcm%namespace)
    end if

    POP_SUB(pcm_poisson_sanity_check)
  end subroutine pcm_poisson_sanity_check
  
  ! -----------------------------------------------------------------------------  
  !> Generates the polarization charge density smearing the charge with a gaussian 
  !> distribution on the mesh nearest neighboring points   of each tessera.
  subroutine pcm_charge_density(pcm, q_pcm, q_pcm_tot, mesh, rho)
    type(pcm_t),     intent(inout) :: pcm 
    FLOAT,           intent(in)    :: q_pcm(:)     !< (1:n_tess)
    FLOAT,           intent(in)    :: q_pcm_tot
    type(mesh_t),    intent(in)    :: mesh
    FLOAT,           intent(out)   :: rho(:)

    integer :: ia
    FLOAT   :: norm, qtot, rr, xx(1:MAX_DIM), pp(1:MAX_DIM)
    
    ! nearest neighbor variables 
    integer :: nm(1:MAX_DIM)
    FLOAT :: posrel(1:MAX_DIM)
    integer :: npt, ipt
    integer :: i1, i2, i3
    integer, allocatable :: pt(:)  
    FLOAT,   allocatable :: lrho(:) ! local charge density on a tessera NN 
    logical :: inner_point, boundary_point
    
    !profiling and debug 
    type(profile_t), save :: prof_init
    integer  :: ierr

    PUSH_SUB(pcm_charge_density)
    call profiling_in(prof_init, 'PCM_CHARGE_DENSITY') 
    
    npt = (2*pcm%tess_nn)**mesh%sb%dim
    SAFE_ALLOCATE(pt(1:npt))
    SAFE_ALLOCATE(lrho(1:npt))

    pt = 0 
    rho = M_ZERO

    do ia = 1, pcm%n_tesserae
      
      pp(1:mesh%sb%dim) = pcm%tess(ia)%point(1:mesh%sb%dim)
      posrel(1:mesh%sb%dim) = pp(1:mesh%sb%dim)/mesh%spacing(1:mesh%sb%dim)

      nm(1:mesh%sb%dim) = floor(posrel(1:mesh%sb%dim))
      
      ! Get the nearest neighboring points
      ipt = 0
      do i1 = -pcm%tess_nn + 1 , pcm%tess_nn
        do i2 = -pcm%tess_nn + 1 , pcm%tess_nn
          do i3 = -pcm%tess_nn + 1 , pcm%tess_nn
            ipt = ipt + 1
            pt(ipt) = index_from_coords(mesh%idx, [i1 + nm(1), i2 + nm(2), i3 + nm(3)])
          end do
        end do
      end do
      
      
      ! Extrapolate the tessera point charge with a gaussian distritibution
      ! to the neighboring points 
      ! rho(r) = N exp(-|r-sk|^2/(alpha*Ak))
      norm = M_ZERO
      lrho = M_ZERO
      do ipt = 1, npt
        
        ! Check the point is inside the mesh skip otherwise
        if (pt(ipt) > 0 .and. pt(ipt) <= mesh%np_part_global) then
          
          if (mesh%parallel_in_domains) then
            pt(ipt) = vec_global2local(mesh%vp, pt(ipt), mesh%vp%partno)
            boundary_point = pt(ipt) > mesh%np + mesh%vp%np_ghost
            inner_point = pt(ipt) > 0 .and. pt(ipt) <= mesh%np

            if (boundary_point .or. inner_point) then
              xx(1:mesh%sb%dim) = mesh%x(pt(ipt),1:mesh%sb%dim)
            else 
              cycle
            end if
          
          else
            xx(1:mesh%sb%dim) = mesh%x(pt(ipt), 1:mesh%sb%dim)
          end if
        
          rr = sum((xx(1:mesh%sb%dim) - pp(1:mesh%sb%dim))**2)
          norm = norm + exp(-rr/(pcm%tess(ia)%area*pcm%gaussian_width))
          lrho(ipt) = lrho(ipt) + exp(-rr/(pcm%tess(ia)%area*pcm%gaussian_width))
          
        end if
      end do

      ! reduce the local density scattered among nodes
      call comm_allreduce(mesh%mpi_grp%comm, lrho, npt)
      
      ! normalize the integral to the tessera point charge q_pcm(ia)
      norm = sum(lrho(1:npt)) * mesh%volume_element
      if (norm > M_EPSILON) then
        norm = q_pcm(ia)/norm
      else   
        norm = M_ZERO
      end if
      lrho(:) = lrho(:) * norm

      ! Add up the local density to the full charge density 
      do ipt = 1, npt
        
        if (pt(ipt) > 0 .and. pt(ipt) <= mesh%np_part_global) then
        
          if (mesh%parallel_in_domains) then
            boundary_point = pt(ipt) > mesh%np + mesh%vp%np_ghost
            inner_point = pt(ipt) > 0 .and. pt(ipt) <= mesh%np
            if(boundary_point .or. inner_point) rho(pt(ipt)) = rho(pt(ipt)) + lrho(ipt)
          else
            rho(pt(ipt)) = rho(pt(ipt)) + lrho(ipt)           
          end if

        end if
      end do
            
    end do

    if (debug%info) then  
      qtot = dmf_integrate(mesh, rho)
      call messages_write(' PCM charge density integrates to q = ')
      call messages_write(qtot)
      call messages_write(' (qtot = ')
      call messages_write(q_pcm_tot)
      call messages_write(')')
      call messages_info()
      
      ! Keep this here for debug purposes.    
      call dio_function_output(io_function_fill_how("VTK"), ".", "rho_pcm", pcm%namespace, &
        mesh, rho, unit_one, ierr)
    end if  

    SAFE_DEALLOCATE_A(pt)
    SAFE_DEALLOCATE_A(lrho)
    
    call profiling_out(prof_init)
    
    POP_SUB(pcm_charge_density)
  end subroutine pcm_charge_density


  ! -----------------------------------------------------------------------------  
  !> Generates the potential 'v_pcm' in real-space.
  subroutine pcm_pot_rs(pcm, v_pcm, q_pcm, rho, mesh, psolver)
    type(pcm_t),     intent(inout) :: pcm 
    FLOAT,           intent(inout) :: v_pcm(:)!< (1:mesh%np) running serially np=np_global
    FLOAT,           intent(in)    :: q_pcm(:)!< (1:n_tess)
    FLOAT,           intent(inout) :: rho(:)
    type(mesh_t),    intent(in)    :: mesh
    type(poisson_t), intent(in)    :: psolver

    type(profile_t), save :: prof_init
    integer  :: ierr

    PUSH_SUB(pcm_pot_rs)
    call profiling_in(prof_init, 'PCM_POT_RS') 

    v_pcm = M_ZERO

    select case(pcm%calc_method)
    case(PCM_CALC_DIRECT)
      call pcm_pot_rs_direct(v_pcm, q_pcm, pcm%tess, pcm%n_tesserae, mesh, pcm%gaussian_width)

    case(PCM_CALC_POISSON)
      call pcm_pot_rs_poisson(v_pcm, psolver, rho)

    case default
      message(1) = "BAD BAD BAD"
      call messages_fatal(1,only_root_writes = .true., namespace=pcm%namespace)

    end select
    
    if (debug%info) then  
      !   Keep this here for debug purposes.    
      call dio_function_output(io_function_fill_how("VTK"), ".", "v_pcm", pcm%namespace, &
        mesh, v_pcm, unit_one, ierr)
    end if

    call profiling_out(prof_init)
    POP_SUB(pcm_pot_rs)
  end subroutine pcm_pot_rs

  
  ! -----------------------------------------------------------------------------  
  !> Generates the potential 'v_pcm' in real-space solving the poisson equation for rho
  subroutine pcm_pot_rs_poisson(v_pcm, psolver, rho)
    FLOAT,           intent(inout) :: v_pcm(:)
    type(poisson_t), intent(in)    :: psolver
    FLOAT,           intent(inout) :: rho(:)
      
    PUSH_SUB(pcm_pot_rs_poisson)
    
    call dpoisson_solve(psolver, v_pcm, rho)
    
    POP_SUB(pcm_pot_rs_poisson)
  end subroutine pcm_pot_rs_poisson


  ! -----------------------------------------------------------------------------  
  !> Generates the potential 'v_pcm' in real-space by direct sum.
  subroutine pcm_pot_rs_direct(v_pcm, q_pcm, tess, n_tess, mesh, width_factor)
    FLOAT,               intent(out) :: v_pcm(:)!< (1:mesh%np) running serially np=np_global
    FLOAT,               intent(in)  :: q_pcm(:)!< (1:n_tess)
    FLOAT,               intent(in)  :: width_factor
    integer,             intent(in)  :: n_tess  
    type(mesh_t),        intent(in)  :: mesh
    type(pcm_tessera_t), intent(in)  :: tess(:) !< (1:n_tess)

    FLOAT, parameter :: P_1 = CNST(0.119763)
    FLOAT, parameter :: P_2 = CNST(0.205117)
    FLOAT, parameter :: Q_1 = CNST(0.137546)
    FLOAT, parameter :: Q_2 = CNST(0.434344)
    FLOAT            :: arg
    FLOAT            :: term
    integer 	     :: ip
    integer          :: ia

    PUSH_SUB(pcm_pot_rs_direct)

    v_pcm = M_ZERO

    if (width_factor /= M_ZERO) then
      !< regularized PCM field
      do ia = 1, n_tess
        do ip = 1, mesh%np
          ! Computing the distances between tesserae and grid points.
          call mesh_r(mesh, ip, term, origin=tess(ia)%point)
          arg = term/sqrt( tess(ia)%area*width_factor )        
          term = (M_ONE + P_1*arg + P_2*arg**2 )/(M_ONE + Q_1*arg + Q_2*arg**2 + P_2*arg**3)
          v_pcm(ip) = v_pcm(ip) + q_pcm(ia)*term/sqrt(tess(ia)%area*width_factor)
        end do
      end do
 
      v_pcm = M_TWO*v_pcm/sqrt(M_PI)

    else
      !< standard PCM field
      do ia = 1, n_tess
        do ip = 1, mesh%np
          ! Computing the distances between tesserae and grid points.
          call mesh_r(mesh, ip, term, origin=tess(ia)%point)
          v_pcm(ip) = v_pcm(ip) + q_pcm(ia)/term         
        end do
      end do
    end if

    POP_SUB(pcm_pot_rs_direct)
  end subroutine pcm_pot_rs_direct

  ! -----------------------------------------------------------------------------
  
  !> Generates the PCM response matrix. J. Tomassi et al. Chem. Rev. 105, 2999 (2005). 
  subroutine pcm_matrix(eps, tess, n_tess, pcm_mat, localf )
    FLOAT,               intent(in)  :: eps
    type(pcm_tessera_t), intent(in)  :: tess(:)      !< (1:n_tess)
    integer,             intent(in)  :: n_tess
    FLOAT,               intent(out) :: pcm_mat(:,:) !< (1:n_tess, 1:n_tess)
    logical,   optional, intent(in)  :: localf

    integer :: i, info
    integer, allocatable :: iwork(:)
    FLOAT, allocatable :: mat_tmp(:,:)

    FLOAT :: sgn_lf

    PUSH_SUB(pcm_matrix)

    !> Conforming the S_I matrix
    SAFE_ALLOCATE(s_mat_act(1:n_tess, 1:n_tess))
    call s_i_matrix(n_tess, tess)

    !> Defining the matrix S_E=S_I/eps
    SAFE_ALLOCATE(Sigma(1:n_tess, 1:n_tess))
    Sigma = s_mat_act/eps

    !> Conforming the D_I matrix
    SAFE_ALLOCATE(d_mat_act(1:n_tess, 1:n_tess))
    call d_i_matrix(n_tess, tess)

    !> Defining the matrix D_E=D_I 
    SAFE_ALLOCATE(Delta(1:n_tess, 1:n_tess))
    Delta = d_mat_act

    sgn_lf = M_ONE
    !> 'local field' differ from 'standard' PCM response matrix in some sign changes
    if (present(localf)) then 
      if (localf) sgn_lf = -M_ONE
    end if

    !> Start conforming the PCM matrix
    pcm_mat = - sgn_lf * d_mat_act

    do i = 1, n_tess
      pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi  
    end do

    SAFE_DEALLOCATE_A(d_mat_act) 

    SAFE_ALLOCATE(iwork(1:n_tess))

    !> Solving for X = S_I^-1*(2*Pi - D_I)
    !! for local field effects ---> X = S_I^-1*(2*Pi + D_I) 
    ! FIXME: use interface, or routine in lalg_adv_lapack_inc.F90
    call dgesv(n_tess, n_tess, s_mat_act, n_tess, iwork, pcm_mat, n_tess, info)        

    SAFE_DEALLOCATE_A(iwork)

    SAFE_DEALLOCATE_A(s_mat_act)

    !> Computing -S_E*S_I^-1*(2*Pi - D_I)
    !! for local field effects ---> -S_E*S_I^-1*(2*Pi + D_I)
    pcm_mat = -matmul(Sigma, pcm_mat)

    do i = 1, n_tess
      pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi
    end do

    pcm_mat = pcm_mat - sgn_lf * Delta

    SAFE_ALLOCATE(mat_tmp(1:n_tess,1:n_tess))
    mat_tmp = M_ZERO

    SAFE_ALLOCATE(d_mat_act(1:n_tess,1:n_tess))
    call d_i_matrix(n_tess, tess)

    mat_tmp = transpose(d_mat_act)

    mat_tmp = matmul(Sigma, mat_tmp)

    mat_tmp = mat_tmp + M_TWO*M_Pi*Sigma

    SAFE_DEALLOCATE_A(d_mat_act)

    SAFE_ALLOCATE(s_mat_act(1:n_tess, 1:n_tess))
    call s_i_matrix(n_tess, tess)

    mat_tmp = mat_tmp + M_TWO*M_Pi*s_mat_act - matmul(Delta, s_mat_act)

    SAFE_DEALLOCATE_A(s_mat_act)
    SAFE_DEALLOCATE_A(Sigma)
    SAFE_DEALLOCATE_A(Delta)

    SAFE_ALLOCATE(iwork(1:n_tess))

    !> Solving for [(2*pi - D_E)*S_I + S_E*(2*Pi + D_I*)]*X = [(2*Pi - D_E) - S_E*S_I^-1*(2*Pi - D_I)]    
    !! for local field ---> [(2*pi - D_E)*S_I + S_E*(2*Pi + D_I*)]*X = [(2*Pi + D_E) - S_E*S_I^-1*(2*Pi + D_I)]
    call dgesv(n_tess, n_tess, mat_tmp, n_tess, iwork, pcm_mat, n_tess, info)		  		  

    SAFE_DEALLOCATE_A(iwork)
    SAFE_DEALLOCATE_A(mat_tmp)

    pcm_mat = - sgn_lf * pcm_mat

    ! Testing
    if (gamess_benchmark) then
      do i = 1, n_tess
        mat_gamess(i,:) = pcm_mat(i,:)/tess(i)%area
      end do
    end if

    POP_SUB(pcm_matrix)
  end subroutine pcm_matrix

  ! -----------------------------------------------------------------------------
  
  subroutine s_i_matrix(n_tess, tess)
    integer,             intent(in)    :: n_tess
    type(pcm_tessera_t), intent(in)    :: tess(:)

    integer :: ii, jj

    s_mat_act = M_ZERO 

    do ii = 1, n_tess
      do jj = ii, n_tess
        s_mat_act(ii,jj) = s_mat_elem_I(tess(ii), tess(jj))
        if (ii /= jj) s_mat_act(jj,ii) = s_mat_act(ii,jj) !symmetric matrix
      end do
    end do

  end subroutine s_i_matrix

  ! -----------------------------------------------------------------------------

  subroutine d_i_matrix(n_tess, tess)
    integer,             intent(in) :: n_tess
    type(pcm_tessera_t), intent(in) :: tess(:)

    integer :: ii, jj

    d_mat_act = M_ZERO 

    do ii = 1, n_tess
      do jj = 1, n_tess !< non-symmetric matrix
        d_mat_act(ii,jj) = d_mat_elem_I(tess(ii), tess(jj))
      end do
    end do

  end subroutine d_i_matrix

  ! -----------------------------------------------------------------------------

  !> electrostatic Green function in vacuo:
  !! G_I(r,r^\prime) = 1 / | r - r^\prime |
  FLOAT function s_mat_elem_I(tessi, tessj)
    type(pcm_tessera_t), intent(in) :: tessi
    type(pcm_tessera_t), intent(in) :: tessj

    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.1)

    FLOAT :: diff(1:PCM_DIM_SPACE)
    FLOAT :: dist
    FLOAT :: s_diag
    FLOAT :: s_off_diag

    s_diag = M_ZERO
    s_off_diag = M_ZERO

    diff = tessi%point - tessj%point

    dist = sqrt(dot_product(diff, diff))

    if (abs(dist) <= M_EPSILON) then
      s_diag = M_SD_DIAG*sqrt(M_FOUR*M_Pi/tessi%area)
      s_mat_elem_I = s_diag
    else
      if (dist > M_DIST_MIN) s_off_diag = M_ONE/dist 
      s_mat_elem_I = s_off_diag
    end if

  end function s_mat_elem_I

  ! -----------------------------------------------------------------------------

  !> Gradient of the Green function in vacuo GRAD[G_I(r,r^\prime)]
  FLOAT function d_mat_elem_I(tessi, tessj)
    type(pcm_tessera_t), intent(in) :: tessi
    type(pcm_tessera_t), intent(in) :: tessj

    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.04)

    FLOAT :: diff(1:PCM_DIM_SPACE)
    FLOAT :: dist
    FLOAT :: d_diag
    FLOAT :: d_off_diag

    d_diag = M_ZERO
    d_off_diag = M_ZERO
    diff = M_ZERO

    diff = tessi%point - tessj%point

    dist = dot_product(diff, diff)
    dist = sqrt(dist)

    if (abs(dist) <= M_EPSILON) then
      !> Diagonal matrix elements  
      d_diag = M_SD_DIAG*sqrt(M_FOUR*M_Pi*tessi%area)
      d_diag = -d_diag/(M_TWO*tessi%r_sphere)
      d_mat_elem_I = d_diag

    else
      !> off-diagonal matrix elements  
      if (dist > M_DIST_MIN) then
        d_off_diag = dot_product( diff, tessj%normal(:) )	
        d_off_diag = d_off_diag*tessj%area/dist**3
      end if

      d_mat_elem_I = d_off_diag

    end if

  end function d_mat_elem_I

  ! -----------------------------------------------------------------------------

  !> It builds the solute cavity surface and calculates the vertices,
  !! representative points and areas of the tesserae by using the 
  !! Gauss-Bonnet theorem.
  subroutine cav_gen(tess_sphere, tess_min_distance, nesf, sfe, nts, cts, unit_pcminfo)
    integer,              intent(in)    :: tess_sphere
    FLOAT  ,              intent(in)    :: tess_min_distance
    type(pcm_sphere_t),   intent(inout) :: sfe(:) !< (1:pcm%n_spheres)
    integer,              intent(in)    :: nesf
    integer,              intent(out)   :: nts
    type(pcm_tessera_t),  intent(out)   :: cts(:) !< (1:pcm%n_tesserae)
    integer,              intent(in)    :: unit_pcminfo

    integer, parameter :: DIM_ANGLES = 24
    integer, parameter :: DIM_TEN = 10
    integer, parameter :: DIM_VERTICES = 122
    integer, parameter :: MAX_VERTICES = 6
    integer, parameter :: MXTS = 10000

    FLOAT :: thev(1:DIM_ANGLES)
    FLOAT :: fiv(1:DIM_ANGLES)
    FLOAT :: fir
    FLOAT :: cv(1:DIM_VERTICES, 1:PCM_DIM_SPACE)
    FLOAT :: th
    FLOAT :: fi
    FLOAT :: cth
    FLOAT :: sth

    FLOAT :: xctst(tess_sphere*N_TESS_SPHERE)
    FLOAT :: yctst(tess_sphere*N_TESS_SPHERE)
    FLOAT :: zctst(tess_sphere*N_TESS_SPHERE)
    FLOAT :: ast(tess_sphere*N_TESS_SPHERE)
    FLOAT :: nctst(PCM_DIM_SPACE, tess_sphere*N_TESS_SPHERE)

    FLOAT :: pts(1:PCM_DIM_SPACE, 1:DIM_TEN)
    FLOAT :: pp(1:PCM_DIM_SPACE)
    FLOAT :: pp1(1:PCM_DIM_SPACE)
    FLOAT :: ccc(1:PCM_DIM_SPACE, 1:DIM_TEN)

    integer :: idum(1:N_TESS_SPHERE*MAX_VERTICES)
    integer :: jvt1(1:MAX_VERTICES,1:N_TESS_SPHERE)
    integer :: isfet(1:DIM_TEN*DIM_ANGLES)

    integer :: ii
    integer :: ia
    integer :: ja
    integer :: nn
    integer :: nsfe
    integer :: its
    integer :: n1
    integer :: n2
    integer :: n3
    integer :: nv
    integer :: i_tes

    FLOAT :: xen
    FLOAT :: yen
    FLOAT :: zen
    FLOAT :: ren
    FLOAT :: area
    FLOAT :: test2
    FLOAT :: rij
    FLOAT :: dnorm

    FLOAT :: xi
    FLOAT :: yi
    FLOAT :: zi
    FLOAT :: xj
    FLOAT :: yj
    FLOAT :: zj

    FLOAT :: vol
    FLOAT :: stot
    FLOAT :: prod
    FLOAT :: dr  

    logical :: band_iter

    PUSH_SUB(cav_gen)

    !> Angles corresponding to the vertices and centres of a polyhedron
    !! within a sphere of unitary radius and centered at the origin
    data thev/ CNST(0.6523581398) , CNST(1.107148718)  , CNST(1.382085796) , &
      CNST(1.759506858)  , CNST(2.034443936)  , CNST(2.489234514) , &
      CNST(0.3261790699) , CNST(0.5535743589), &
      CNST(0.8559571251) , CNST(0.8559571251) , CNST(1.017221968) , &
      CNST(1.229116717)  , CNST(1.229116717)  , CNST(1.433327788) , &
      CNST(1.570796327)  , CNST(1.570796327)  , CNST(1.708264866) , &
      CNST(1.912475937)  , CNST(1.912475937)  , CNST(2.124370686) , &
      CNST(2.285635528)  , CNST(2.285635528)  , CNST(2.588018295) , &
      CNST(2.815413584) /
    data fiv/                       CNST(0.6283185307) , M_ZERO            , &
      CNST(0.6283185307) , M_ZERO             , CNST(0.6283185307), &
      M_ZERO             , CNST(0.6283185307) , M_ZERO, 	     &
      CNST(0.2520539002) , CNST(1.004583161)  , CNST(0.6283185307), &
      CNST(0.3293628477) , CNST(0.9272742138) , M_ZERO, 	     &
      CNST(0.3141592654) , CNST(0.9424777961) , CNST(0.6283185307), &
      CNST(0.2989556830) , CNST(0.9576813784) , M_ZERO, 	     &
      CNST(0.3762646305) , CNST(0.8803724309) , CNST(0.6283188307), &
      M_ZERO /
    data fir / CNST(1.256637061)  /

    !> the vector idum, contained in the matrix jvt1, indicates the vertices 
    !! of the tesserae (using less than 19 continuations)
    data (idum(ii),ii = 1, 280) /                                   &
      1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,         &
      33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51, &
      42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,    &
      49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,     &
      3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,   &
      55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,      &
      43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,   &
      16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,   &
      78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,     &
      17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73,  &
      9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,     &
      61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,    &
      24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,    &
      91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,    &
      24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,    &
      86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98,  &
      24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,     &
      31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93 /
    data (idum(ii),ii = 281,360) / 				      &
      108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,   &
      26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,    &
      29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,      &
      110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,    &
      118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,   &
      114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /

    if (mpi_grp_is_root(mpi_world)) then
      if (tess_sphere == 1) then
        write(unit_pcminfo,'(A1)')  '#' 
        write(unit_pcminfo,'(A34)') '# Number of tesserae / sphere = 60'
        write(unit_pcminfo,'(A1)')  '#' 
      else
        write(unit_pcminfo,'(A1)')  '#' 
        write(unit_pcminfo,'(A35)') '# Number of tesserae / sphere = 240' 
        write(unit_pcminfo,'(A1)')  '#' 
      end if
    end if

    !> geometrical data are converted to Angstrom and back transformed
    !! to Bohr at the end of the subroutine.
    dr = CNST(0.01)
    dr = dr*P_a_B

    sfe(:)%x = sfe(:)%x*P_a_B
    sfe(:)%y = sfe(:)%y*P_a_B
    sfe(:)%z = sfe(:)%z*P_a_B
    sfe(:)%r = sfe(:)%r*P_a_B

    vol  = M_ZERO
    stot = M_ZERO
    jvt1 = reshape(idum,(/6,60/))

    !> Coordinates of vertices of tesserae in a sphere with unit radius.
    !! the matrix 'cv' and 'xc', 'yc', 'zc' conatin the vertices and 
    !! the centers of 240 tesserae. The matrix 'jvt1(i,j)' denotes the index
    !! of the i-th vertex of the j-th big tessera. On each big tessera
    !! the 6 vertices are ordered as follows:
    !!  
    !!                      1
    !!
    !!                   4     5
    !!
    !!                3     6     2

    cv(1,1)   =  M_ZERO
    cv(1,2)   =  M_ZERO
    cv(1,3)   =  M_ONE

    cv(122,1) =  M_ZERO
    cv(122,2) =  M_ZERO
    cv(122,3) = -M_ONE

    ii = 1
    do ia = 1, DIM_ANGLES
      th = thev(ia)
      fi = fiv(ia)
      cth = cos(th)
      sth = sin(th)
      do ja = 1, 5
        fi = fi + fir
        if (ja == 1) fi = fiv(ia)
        ii = ii + 1
        cv(ii,1) = sth*cos(fi)
        cv(ii,2) = sth*sin(fi)
        cv(ii,3) = cth
      end do
    end do

    !> Controls whether the tessera is covered or need to be reshaped it
    nn = 0
    do nsfe = 1, nesf
      xen = sfe(nsfe)%x
      yen = sfe(nsfe)%y
      zen = sfe(nsfe)%z
      ren = sfe(nsfe)%r

      xctst(:) = M_ZERO
      yctst(:) = M_ZERO
      zctst(:) = M_ZERO
      ast(:)   = M_ZERO

      do its = 1, N_TESS_SPHERE 
        do i_tes = 1, tess_sphere
          if (tess_sphere == 1) then
            n1 = jvt1(1,its)
            n2 = jvt1(2,its)
            n3 = jvt1(3,its)
          else
            if (i_tes == 1)      then
              n1 = jvt1(1,its)
              n2 = jvt1(5,its)
              n3 = jvt1(4,its)
            elseif (i_tes == 2)  then 
              n1 = jvt1(4,its)
              n2 = jvt1(6,its)
              n3 = jvt1(3,its)
            elseif (i_tes == 3)  then
              n1 = jvt1(4,its)
              n2 = jvt1(5,its)
              n3 = jvt1(6,its)
            elseif (i_tes == 4)  then
              n1 = jvt1(2,its)
              n2 = jvt1(6,its)
              n3 = jvt1(5,its)
            end if
          end if

          pts(1,1) = cv(n1,1)*ren + xen
          pts(2,1) = cv(n1,3)*ren + yen
          pts(3,1) = cv(n1,2)*ren + zen

          pts(1,2) = cv(n2,1)*ren + xen
          pts(2,2) = cv(n2,3)*ren + yen
          pts(3,2) = cv(n2,2)*ren + zen

          pts(1,3) = cv(n3,1)*ren + xen
          pts(2,3) = cv(n3,3)*ren + yen
          pts(3,3) = cv(n3,2)*ren + zen

          pp(:)  = M_ZERO
          pp1(:) = M_ZERO
          nv = 3

          call subtessera(sfe, nsfe, nesf, nv, pts ,ccc, pp, pp1, area)

          if (abs(area) <= M_EPSILON) cycle

          xctst(tess_sphere*(its-1) + i_tes)   = pp(1)
          yctst(tess_sphere*(its-1) + i_tes)   = pp(2)
          zctst(tess_sphere*(its-1) + i_tes)   = pp(3)
          nctst(:,tess_sphere*(its-1) + i_tes) = pp1(:)
          ast(tess_sphere*(its-1) + i_tes)     = area
          isfet(tess_sphere*(its-1) + i_tes)   = nsfe

        end do
      end do !> loop through the tesseare on the sphere 'nsfe'

      do its = 1, N_TESS_SPHERE*tess_sphere

        if (abs(ast(its)) <= M_EPSILON) cycle
        nn = nn + 1

        if (nn > MXTS) then !> check the total number of tessera
          write(message(1),'(a,I5,a,I5)') "total number of tesserae", nn, ">",MXTS
          call messages_warning(1)
        end if

          cts(nn)%point(1)  = xctst(its)
          cts(nn)%point(2)  = yctst(its)
          cts(nn)%point(3)  = zctst(its)
          cts(nn)%normal(:) = nctst(:,its)
          cts(nn)%area      = ast(its)
          cts(nn)%r_sphere  = sfe(isfet(its))%r

      end do
    end do !> loop through the spheres

    nts = nn

    !> checks if two tesseare are too close
    test2 = tess_min_distance*tess_min_distance

    band_iter = .false.
    do while (.not.(band_iter))
      band_iter = .true.

      loop_ia: do ia = 1, nts-1
        if (abs(cts(ia)%area) <= M_EPSILON) cycle
        xi = cts(ia)%point(1)
        yi = cts(ia)%point(2)
        zi = cts(ia)%point(3)

        loop_ja: do ja = ia+1, nts
          if (abs(cts(ja)%area) <= M_EPSILON) cycle
          xj = cts(ja)%point(1)
          yj = cts(ja)%point(2)
          zj = cts(ja)%point(3)

          rij = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2

          if (rij > test2) cycle

          if (mpi_grp_is_root(mpi_world)) &
            write(unit_pcminfo,'(A40,I4,A5,I4,A4,F8.4,A13,F8.4,A3)' ) &
            '# Warning: The distance between tesserae', &
            ia,' and ', ja,' is ',sqrt(rij),' A, less than', tess_min_distance,' A.'

          !> calculating the coordinates of the new tessera weighted by the areas
          xi = (xi*cts(ia)%area + xj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
          yi = (yi*cts(ia)%area + yj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
          zi = (zi*cts(ia)%area + zj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)

          cts(ia)%point(1) = xi
          cts(ia)%point(2) = yi
          cts(ia)%point(3) = zi

          !> calculating the normal vector of the new tessera weighted by the areas
          cts(ia)%normal = (cts(ia)%normal*cts(ia)%area + cts(ja)%normal*cts(ja)%area)
          dnorm = sqrt(dot_product(cts(ia)%normal, cts(ia)%normal))
          cts(ia)%normal = cts(ia)%normal/dnorm

          !> calculating the sphere radius of the new tessera weighted by the areas
          cts(ia)%r_sphere = (cts(ia)%r_sphere*cts(ia)%area + cts(ja)%r_sphere*cts(ja)%area) / &
            (cts(ia)%area + cts(ja)%area)

          !> calculating the area of the new tessera
          cts(ia)%area = cts(ia)%area + cts(ja)%area

          !> deleting tessera ja
          do ii = ja+1, nts
            cts(ii-1) = cts(ii)
          end do
          nts = nts - 1 ! updating number of tesserae
          band_iter = .false.
          exit loop_ia

        end do loop_ja
      end do loop_ia
    end do !> while loop

    !> Calculates the cavity volume: vol = \sum_{its=1}^nts A_{its} s*n/3.
    vol = M_ZERO
    do its = 1, nts
      prod = dot_product( cts(its)%point, cts(its)%normal ) 
      vol  = vol + cts(its)%area * prod / M_THREE
      stot = stot + cts(its)%area
    end do

    if ( mpi_grp_is_root(mpi_world) ) then
      write(unit_pcminfo, '(A2)')  '# '
      write(unit_pcminfo, '(A29,I4)')    '# Total number of tesserae = ', nts
      write(unit_pcminfo, '(A30,F12.6)') '# Cavity surface area (A^2) = ' , stot
      write(unit_pcminfo, '(A24,F12.6)') '# Cavity volume (A^3) = '       , vol
      write(unit_pcminfo, '(A2)')  '# '
    end if

    !> transforms results into Bohr.
    cts(:)%area     = cts(:)%area*P_Ang*P_Ang
    cts(:)%point(1) = cts(:)%point(1)*P_Ang
    cts(:)%point(2) = cts(:)%point(2)*P_Ang
    cts(:)%point(3) = cts(:)%point(3)*P_Ang
    cts(:)%r_sphere = cts(:)%r_sphere*P_Ang

    sfe(:)%x=sfe(:)%x*P_Ang
    sfe(:)%y=sfe(:)%y*P_Ang
    sfe(:)%z=sfe(:)%z*P_Ang
    sfe(:)%r=sfe(:)%r*P_Ang

    POP_SUB(cav_gen)
  end subroutine cav_gen
  
  ! -----------------------------------------------------------------------------

  !> find the uncovered region for each tessera and computes the area,
  !! the representative point (pp) and the unitary normal vector (pp1)
  subroutine subtessera(sfe, ns, nesf, nv, pts, ccc, pp, pp1, area)
    type(pcm_sphere_t), intent(in)    :: sfe(:) !< (1:nesf)
    integer,            intent(in)    :: ns 
    integer,            intent(in)    :: nesf
    integer,            intent(inout) :: nv
    FLOAT,              intent(inout) :: pts(:,:) !< (1:PCM_DIM_SPACE, 1:DIM_TEN)
    FLOAT,              intent(out)   :: ccc(:,:) !< (1:PCM_DIM_SPACE, 1:DIM_TEN)
    FLOAT,              intent(out)   :: pp(:)    !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(out)   :: pp1(:)   !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(out)   :: area

    FLOAT, parameter   :: TOL = -CNST(1.0e-10)
    integer, parameter :: DIM_TEN = 10

    integer :: intsph(1:DIM_TEN)
    integer :: nsfe1
    integer :: na
    integer :: icop
    integer :: ll
    integer :: iv1
    integer :: iv2
    integer :: ii
    integer :: icut
    integer :: jj

    FLOAT  :: p1(1:PCM_DIM_SPACE)
    FLOAT  :: p2(1:PCM_DIM_SPACE)
    FLOAT  :: p3(1:PCM_DIM_SPACE)
    FLOAT  :: p4(1:PCM_DIM_SPACE)
    FLOAT  :: point(1:PCM_DIM_SPACE)
    FLOAT  :: pscr(1:PCM_DIM_SPACE,1:DIM_TEN)
    FLOAT  :: cccp(1:PCM_DIM_SPACE,1:DIM_TEN)
    FLOAT  :: pointl(1:PCM_DIM_SPACE,1:DIM_TEN)
    FLOAT  :: diff(1:PCM_DIM_SPACE)

    integer :: ind(1:DIM_TEN)
    integer :: ltyp(1:DIM_TEN)
    integer :: intscr(1:DIM_TEN)

    FLOAT  :: delr
    FLOAT  :: delr2
    FLOAT  :: rc
    FLOAT  :: rc2
    FLOAT  :: dnorm
    FLOAT  :: dist
    FLOAT  :: de2

    PUSH_SUB(subtessera)

    p1     = M_ZERO
    p2     = M_ZERO
    p3     = M_ZERO
    p4     = M_ZERO
    point  = M_ZERO
    pscr   = M_ZERO
    cccp   = M_ZERO
    pointl = M_ZERO
    diff   = M_ZERO
    area   = M_ZERO

    do jj = 1, 3
      ccc(1,jj) = sfe(ns)%x
      ccc(2,jj) = sfe(ns)%y
      ccc(3,jj) = sfe(ns)%z
    end do

    intsph = ns
    do nsfe1 = 1, nesf 
      if (nsfe1 == ns) cycle
      do jj = 1, nv
        intscr(jj) = intsph(jj)
        pscr(:,jj) = pts(:,jj)
        cccp(:,jj) = ccc(:,jj)
      end do

      icop = 0
      ind = 0
      ltyp = 0

      do ii = 1, nv
        delr2 = ( pts(1,ii) - sfe(nsfe1)%x )**2 + ( pts(2,ii) - sfe(nsfe1)%y )**2 + (pts(3,ii) - sfe(nsfe1)%z)**2
        delr = sqrt(delr2)
        if (delr < sfe(nsfe1)%r) then
          ind(ii) = 1
          icop = icop + 1
        end if
      end do

      if (icop == nv) then 
        POP_SUB(subtessera)
        return
      end if

      do ll = 1, nv
        iv1 = ll
        iv2 = ll+1
        if (ll == nv) iv2 = 1
        if ((ind(iv1) == 1) .and. (ind(iv2) == 1)) then
          ltyp(ll) = 0
        else if ((ind(iv1) == 0) .and. (ind(iv2) == 1)) then
          ltyp(ll) = 1
        else if ((ind(iv1) == 1) .and. (ind(iv2) == 0)) then
          ltyp(ll) = 2
        else if ((ind(iv1) == 0) .and. (ind(iv2) == 0)) then
          ltyp(ll) = 4
          diff = ccc(:,ll) - pts(:,ll)
          rc2 = dot_product(diff, diff)
          rc = sqrt(rc2)

          do ii = 1, 11
            point = pts(:,iv1) + ii * (pts(:,iv2) - pts(:,iv1)) / 11
            point = point - CCC(:,ll)
            dnorm = sqrt( dot_product(point, point) )
            point = point * rc / dnorm + CCC(:,ll)

            dist = sqrt((point(1) - sfe(nsfe1)%x)**2 + (point(2) - sfe(nsfe1)%y)**2 + (point(3) - sfe(nsfe1)%z)**2 )

            if ((dist - sfe(nsfe1)%r) < TOL) then
              ltyp(ll) = 3
              pointl(:, ll) = point
              exit
            end if

          end do
        end if
      end do

      icut = 0
      do ll = 1, nv
        if ((ltyp(ll) == 1) .or. (ltyp(ll) == 2)) icut = icut + 1
        if (ltyp(ll) == 3) icut = icut + 2
      end do
      icut = icut / 2
      if (icut > 1) then 
        POP_SUB(subtessera)
        return
      end if

      na = 1
      do ll = 1, nv

        if (ltyp(ll) == 0) cycle
        iv1 = ll
        iv2 = ll + 1
        if (ll == nv) iv2 = 1

        if (ltyp(ll) == 1) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
          p1 = pscr(:,iv1)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter(sfe, p1, p2, p3, p4, nsfe1, 0)
          pts(:,na) = p4

          de2 = ( sfe(nsfe1)%x - sfe(ns)%x )**2 + ( sfe(nsfe1)%y - sfe(ns)%y )**2 + &
            ( sfe(nsfe1)%z - sfe(ns)%z )**2

          ccc(1,na) = sfe(ns)%x + ( sfe(nsfe1)%x - sfe(ns)%x)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

          ccc(2,na) = sfe(ns)%y + ( sfe(nsfe1)%y - sfe(ns)%y)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

          ccc(3,na) = sfe(ns)%z + ( sfe(nsfe1)%z - sfe(ns)%z)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

          intsph(na) = nsfe1
          na = na + 1
        end if

        if (ltyp(ll) == 2) then
          p1 = pscr(:,iv1)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter(sfe, p1, p2, p3, p4, nsfe1, 1)
          pts(:,na) = p4
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if

        if (ltyp(ll) == 3) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
          p1 = pscr(:,iv1)
          p2 = pointl(:,ll)
          p3 = cccp(:,iv1)

          call inter( sfe, p1, p2, p3, p4, nsfe1, 0 )
          pts(:,na) = p4

          de2 = (sfe(nsfe1)%x - sfe(ns)%x)**2 + (sfe(nsfe1)%y - sfe(ns)%y)**2 + (sfe(nsfe1)%z - sfe(ns)%z)**2

          ccc(1,na) = sfe(ns)%x + (sfe(nsfe1)%x - sfe(ns)%x) * (sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2) / (M_TWO*de2)

          ccc(2,na) = sfe(ns)%y + (sfe(nsfe1)%y - sfe(ns)%y) * (sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2) / (M_TWO*de2)

          ccc(3,na) = sfe(ns)%z + (sfe(nsfe1)%z - sfe(ns)%z) * (sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2) / (M_TWO*de2)

          intsph(na) = nsfe1
          na = na + 1
          p1 = pointl(:,ll)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter(sfe, p1, p2, p3, p4, nsfe1, 1)
          pts(:,na) = p4
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if

        if (ltyp(ll) == 4) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if
      end do

      nv = na - 1
      if (nv > 10) then
        message(1) = "Too many vertices on the tessera"
        call messages_fatal(1)
      end if
    end do

    call gaubon(sfe, nv, ns, pts, ccc, pp, pp1, area, intsph)
    
    POP_SUB(subtessera)
  end subroutine subtessera
  
  ! -----------------------------------------------------------------------------

  !    !> Finds the point 'p4', on the arc 'p1'-'p2' developed from 'p3',
  !    !! which is on the surface of sphere 'ns'. p4 is a linear combination 
  !! of p1 and p2 with the 'alpha' parameter optimized iteratively.
  subroutine inter(sfe, p1, p2, p3, p4, ns, ia)
    type(pcm_sphere_t), intent(in)  :: sfe(:) !< (1:nesf)
    FLOAT,              intent(in)  :: p1(1:PCM_DIM_SPACE)  !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(in)  :: p2(1:PCM_DIM_SPACE)  !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(in)  :: p3(1:PCM_DIM_SPACE)  !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(out) :: p4(1:PCM_DIM_SPACE)  !< (1:PCM_DIM_SPACE)
    integer,            intent(in)  :: ns
    integer,            intent(in)  :: ia

    FLOAT, parameter :: TOL = CNST(1.0e-08)

    integer :: m_iter
    FLOAT  :: r
    FLOAT  :: alpha
    FLOAT  :: delta
    FLOAT  :: dnorm
    FLOAT  :: diff
    FLOAT  :: diff_vec(1:PCM_DIM_SPACE)
    logical :: band_iter

    PUSH_SUB(inter)

    diff_vec = p1 - p3
    r = sqrt(dot_product(diff_vec, diff_vec))

    alpha = M_HALF
    delta = M_ZERO
    m_iter = 1

    band_iter = .false.
    do while(.not.(band_iter))
      if (m_iter > 1000) then
        message(1) = "Too many iterations inside subrotuine inter"
        call messages_fatal(1)
      end if

      band_iter = .true.

      alpha = alpha + delta

      p4 = p1 + alpha*(p2 - p1) - p3
      dnorm = sqrt(dot_product(p4, p4))
      p4 = p4*r/dnorm + p3
      diff = (p4(1) - sfe(ns)%x)**2 + (p4(2) - sfe(ns)%y)**2 + (p4(3) - sfe(ns)%z)**2
      diff = sqrt(diff) - sfe(ns)%r

      if (abs(diff) < TOL) then
        POP_SUB(inter)
        return
      end if

      if (ia == 0) then
        if (diff > M_ZERO) delta =  M_ONE/(M_TWO**(m_iter + 1))
        if (diff < M_ZERO) delta = -M_ONE/(M_TWO**(m_iter + 1))
        m_iter = m_iter + 1
        band_iter = .false.
      end if

      if (ia == 1) then
        if (diff > M_ZERO) delta = -M_ONE/(M_TWO**(m_iter + 1))
        if (diff < M_ZERO) delta =  M_ONE/(M_TWO**(m_iter + 1))
        m_iter = m_iter + 1
        band_iter = .false.
      end if
    end do
    
    POP_SUB(inter)
  end subroutine inter

  ! -----------------------------------------------------------------------------

  !> Use the Gauss-Bonnet theorem to calculate the area of the 
  !! tessera with vertices 'pts(3,nv)'. 
  !! Area = R^2 [ 2pi + S(Phi(N)cosT(N)) - S(Beta(N)) ]
  !! Phi(n): length of the arc in radians of the side 'n'. 
  !! T(n): azimuthal angle for the side 'n'
  !! Beta(n): external angle respect to vertex 'n'.
  subroutine gaubon( sfe, nv, ns, pts, ccc, pp, pp1, area, intsph )
    type(pcm_sphere_t), intent(in)    :: sfe(:)    !< (1:nesf)
    FLOAT,              intent(in)    :: pts(:,:)  !< (1:PCM_DIM_SPACE,1:DIM_TEN) 
    FLOAT,              intent(in)    :: ccc(:,:)  !< (1:PCM_DIM_SPACE,1:DIM_TEN)
    FLOAT,              intent(inout) :: pp(:)     !< (1:PCM_DIM_SPACE)
    FLOAT,              intent(inout) :: pp1(:)    !< (1:PCM_DIM_SPACE)
    integer,            intent(in)    :: intsph(:) !< (1:DIM_TEN)
    FLOAT,              intent(out)   :: area
    integer,            intent(in)    :: nv
    integer,            intent(in)    :: ns

    FLOAT :: p1(1:PCM_DIM_SPACE), p2(1:PCM_DIM_SPACE), p3(1:PCM_DIM_SPACE)
    FLOAT :: u1(1:PCM_DIM_SPACE), u2(1:PCM_DIM_SPACE)
    FLOAT :: point_1(1:PCM_DIM_SPACE), point_2(1:PCM_DIM_SPACE)
    FLOAT :: tpi, sum1, dnorm, dnorm1, dnorm2
    FLOAT :: cosphin, phin, costn, sum2, betan
    integer :: nsfe1, ia, nn, n0, n1, n2

    PUSH_SUB(gaubon)

    point_1 = M_ZERO
    point_2 = M_ZERO
    p1      = M_ZERO
    p2      = M_ZERO
    p3      = M_ZERO
    u1      = M_ZERO
    u2      = M_ZERO

    tpi = M_TWO*M_Pi
    sum1 = M_ZERO
    do nn = 1, nv
      point_1 = pts(:,nn) - ccc(:,nn)
      if (nn < nv) then
        point_2 = pts(:,nn+1) - ccc(:,nn)
      else
        point_2 = pts(:,1) - ccc(:,nn)
      end if

      dnorm1 = sqrt( dot_product(point_1, point_1) )
      dnorm2 = sqrt( dot_product(point_2, point_2) )
      cosphin = dot_product(point_1, point_2) / (dnorm1*dnorm2)

      if (cosphin >  M_ONE) cosphin =  M_ONE
      if (cosphin < -M_ONE) cosphin = -M_ONE

      phin = acos(cosphin)
      nsfe1 = intsph(nn)

      point_1(1) = sfe(nsfe1)%x - sfe(ns)%x
      point_1(2) = sfe(nsfe1)%y - sfe(ns)%y
      point_1(3) = sfe(nsfe1)%z - sfe(ns)%z

      dnorm1 = sqrt( dot_product(point_1, point_1) )

      if (abs(dnorm1) <= M_EPSILON) dnorm1 = M_ONE

      point_2(1) = pts(1,nn) - sfe(ns)%x
      point_2(2) = pts(2,nn) - sfe(ns)%y
      point_2(3) = pts(3,nn) - sfe(ns)%z

      dnorm2 = sqrt( dot_product(point_2, point_2) )

      costn  = dot_product(point_1, point_2)/(dnorm1*dnorm2)
      sum1 = sum1 + phin * costn
    end do

    sum2 = M_ZERO
    !> Loop over the vertices
    do nn = 1, nv
      p1 = M_ZERO
      p2 = M_ZERO    
      p3 = M_ZERO  

      n1 = nn
      if (nn > 1)   n0 = nn - 1
      if (nn == 1)  n0 = nv
      if (nn < nv)  n2 = nn + 1
      if (nn == nv) n2 = 1

      p1 = pts(:,n1) - ccc(:,n0)
      p2 = pts(:,n0) - ccc(:,n0)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3    

      call vecp(p1, p2, p3, dnorm)
      u1 = p3/dnorm

      p1 = pts(:,n1) - ccc(:,n1)
      p2 = pts(:,n2) - ccc(:,n1)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3

      call vecp(p1, p2, p3, dnorm)
      u2 = p3/dnorm

      betan = acos( dot_product(u1, u2) )
      sum2 = sum2 + (M_Pi - betan)
    end do

    !> computes the area of the tessera
    area = sfe(ns)%r*sfe(ns)%r*(tpi + sum1 - sum2)

    !> computes the representative point
    pp = M_ZERO

    do ia = 1, nv
      pp(1) = pp(1) + ( pts(1,ia) - sfe(ns)%x )
      pp(2) = pp(2) + ( pts(2,ia) - sfe(ns)%y )
      pp(3) = pp(3) + ( pts(3,ia) - sfe(ns)%z )
    end do

    dnorm = M_ZERO
    dnorm = sqrt( dot_product(pp,pp) )

    pp(1) = sfe(ns)%x + pp(1) * sfe(ns)%r / dnorm
    pp(2) = sfe(ns)%y + pp(2) * sfe(ns)%r / dnorm
    pp(3) = sfe(ns)%z + pp(3) * sfe(ns)%r / dnorm

    !> finds the internal normal at the representative point
    pp1(1) = (pp(1) - sfe(ns)%x) / sfe(ns)%r
    pp1(2) = (pp(2) - sfe(ns)%y) / sfe(ns)%r
    pp1(3) = (pp(3) - sfe(ns)%z) / sfe(ns)%r

    !> If the area of the tessera is negative (0^-), due to numerical errors, is discarded
    if (area < M_ZERO) area = M_ZERO

    POP_SUB(gaubon)
  end subroutine gaubon

  ! -----------------------------------------------------------------------------

  !> calculates the vectorial product p3 = p1 x p2
  subroutine vecp(p1, p2, p3, dnorm)
    FLOAT, intent(in)  :: P1(:) !< (1:PCM_DIM_SPACE)
    FLOAT, intent(in)  :: P2(:) !< (1:PCM_DIM_SPACE)
    FLOAT, intent(out) :: P3(:) !< (1:PCM_DIM_SPACE)
    FLOAT, intent(out) :: dnorm

    p3 = M_ZERO
    p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
    p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
    P3(3) = p1(1)*p2(2) - p1(2)*p2(1)

    dnorm = sqrt(dot_product(p3, p3))
 
  end subroutine vecp

  subroutine pcm_end(pcm)
    type(pcm_t), intent(inout) :: pcm

    integer :: asc_unit_test, & 
               asc_unit_test_sol, &
               asc_unit_test_e, &
               asc_unit_test_n, &
               asc_unit_test_ext
    integer :: ia

    PUSH_SUB(pcm_end)

    if (pcm%solute .and. pcm%localf) then
      asc_unit_test     = io_open(PCM_DIR//'ASC.dat', pcm%namespace, action='write')
      asc_unit_test_sol = io_open(PCM_DIR//'ASC_sol.dat', pcm%namespace, action='write')
      asc_unit_test_e   = io_open(PCM_DIR//'ASC_e.dat', pcm%namespace, action='write')
      asc_unit_test_n   = io_open(PCM_DIR//'ASC_n.dat', pcm%namespace, action='write')
      asc_unit_test_ext = io_open(PCM_DIR//'ASC_ext.dat', pcm%namespace, action='write')
      do ia = 1, pcm%n_tesserae
        write(asc_unit_test,*)     pcm%tess(ia)%point, pcm%q_e(ia) + pcm%q_n(ia) + pcm%q_ext(ia), ia
        write(asc_unit_test_sol,*) pcm%tess(ia)%point, pcm%q_e(ia) + pcm%q_n(ia), ia
        write(asc_unit_test_e,*)   pcm%tess(ia)%point, pcm%q_e(ia), ia
        write(asc_unit_test_n,*)   pcm%tess(ia)%point, pcm%q_n(ia), ia
        write(asc_unit_test_ext,*) pcm%tess(ia)%point, pcm%q_ext(ia), ia
      end do
      call io_close(asc_unit_test)
      call io_close(asc_unit_test_sol)
      call io_close(asc_unit_test_e)
      call io_close(asc_unit_test_n)
      call io_close(asc_unit_test_ext)

    else if (pcm%solute .and. .not.pcm%localf) then
      asc_unit_test_sol = io_open(PCM_DIR//'ASC_sol.dat', pcm%namespace, action='write')
      asc_unit_test_e   = io_open(PCM_DIR//'ASC_e.dat', pcm%namespace, action='write')
      asc_unit_test_n   = io_open(PCM_DIR//'ASC_n.dat', pcm%namespace, action='write')
      do ia = 1, pcm%n_tesserae
        write(asc_unit_test_sol,*) pcm%tess(ia)%point, pcm%q_e(ia) + pcm%q_n(ia), ia
        write(asc_unit_test_e,*)   pcm%tess(ia)%point, pcm%q_e(ia), ia
        write(asc_unit_test_n,*)   pcm%tess(ia)%point, pcm%q_n(ia), ia
      end do
      call io_close(asc_unit_test_sol)
      call io_close(asc_unit_test_e)
      call io_close(asc_unit_test_n)

    else if (.not. pcm%solute .and. pcm%localf) then
      asc_unit_test_ext = io_open(PCM_DIR//'ASC_ext.dat', pcm%namespace, action='write')
      do ia = 1, pcm%n_tesserae
        write(asc_unit_test_ext,*) pcm%tess(ia)%point, pcm%q_ext(ia), ia
      end do
      call io_close(asc_unit_test_ext)
    end if

    SAFE_DEALLOCATE_A(pcm%spheres)
    SAFE_DEALLOCATE_A(pcm%tess)
    SAFE_DEALLOCATE_A(pcm%matrix)
    if (pcm%epsilon_infty /= pcm%epsilon_0) then
     SAFE_DEALLOCATE_A(pcm%matrix_d)
    end if
    SAFE_DEALLOCATE_A(pcm%q_e)
    SAFE_DEALLOCATE_A(pcm%q_e_in)
    SAFE_DEALLOCATE_A(pcm%q_n) 
    SAFE_DEALLOCATE_A(pcm%v_e)
    SAFE_DEALLOCATE_A(pcm%v_n)
    SAFE_DEALLOCATE_A(pcm%v_e_rs)
    SAFE_DEALLOCATE_A(pcm%v_n_rs)
    if (pcm%localf) then
      SAFE_DEALLOCATE_A(pcm%matrix_lf)
      if (pcm%epsilon_infty /= pcm%epsilon_0) then
        SAFE_DEALLOCATE_A(pcm%matrix_lf_d)
      end if

      SAFE_DEALLOCATE_A(pcm%q_ext)
      SAFE_DEALLOCATE_A(pcm%q_ext_in)
      SAFE_DEALLOCATE_A(pcm%v_ext)
      SAFE_DEALLOCATE_A(pcm%v_ext_rs)
      if (pcm%kick_is_present) then
        SAFE_DEALLOCATE_A(pcm%q_kick)
        SAFE_DEALLOCATE_A(pcm%v_kick)
        SAFE_DEALLOCATE_A(pcm%v_kick_rs)
      end if 
    end if

    if (pcm%calc_method == PCM_CALC_POISSON) then
      SAFE_DEALLOCATE_A( pcm%rho_n)
      SAFE_DEALLOCATE_A( pcm%rho_e)
      if (pcm%localf) then
        SAFE_DEALLOCATE_A( pcm%rho_ext)
        if (pcm%kick_is_present) then
          SAFE_DEALLOCATE_A( pcm%rho_kick)
        end if
      endif
    end if

    if (pcm%tdlevel == PCM_TD_EOM) call pcm_eom_end()

    if (mpi_grp_is_root(mpi_world)) call io_close(pcm%info_unit)

    POP_SUB(pcm_end)
  end subroutine pcm_end

  ! -----------------------------------------------------------------------------
  !> Update pcm potential
  logical function pcm_update(this) result(update)
    type(pcm_t), intent(inout) :: this

    this%iter = this%iter + 1 
    update = this%iter <= 6 .or. mod(this%iter, this%update_iter) == 0

    if (debug%info .and. update) then
      call messages_write(' PCM potential updated')
      call messages_new_line()
      call messages_write(' PCM update iteration counter: ')
      call messages_write(this%iter)
      call messages_info()
    end if

  end function pcm_update

  ! -----------------------------------------------------------------------------
  !> get the vdw radius
  FLOAT function pcm_get_vdw_radius(species, pcm_vdw_type, namespace)  result(vdw_r)
    type(species_t),   intent(in) :: species
    integer,           intent(in) :: pcm_vdw_type
    type(namespace_t), intent(in) :: namespace
  
    integer            :: ia
    integer, parameter :: UPTO_XE = 54
    FLOAT              :: vdw_radii(1:UPTO_XE) !< van der Waals radii in Angstrom for elements H-Xe reported 
         !  by Stefan Grimme in J. Comput. Chem. 27: 1787-1799, 2006 
         !  except for C, N and O, reported in J. Chem. Phys. 120, 3893 (2004). 
    data (vdw_radii(ia), ia=1, UPTO_XE)                                                                                  / &
      !H                                                                                                       He 
      CNST(1.001),                                                                                            CNST(1.012), &
      !Li           Be                        B            C            N            O            F            Ne
      CNST(0.825), CNST(1.408),              CNST(1.485), CNST(2.000), CNST(1.583), CNST(1.500), CNST(1.287), CNST(1.243), & 
      !Na           Mg                        Al           Si           P            S            Cl           Ar 
      CNST(1.144), CNST(1.364),              CNST(1.639), CNST(1.716), CNST(1.705), CNST(1.683), CNST(1.639), CNST(1.595), & 
      !K            Ca 
      CNST(1.485), CNST(1.474),                                                                                            & 
      !>      Sc -- Zn       <!                                        
      CNST(1.562), CNST(1.562),                                                                    & 
      CNST(1.562), CNST(1.562),                                                                    & 
      CNST(1.562), CNST(1.562),                                                                    & 
      CNST(1.562), CNST(1.562),                                                                    & 
      CNST(1.562), CNST(1.562),                                                                  & 
      !Ga           Ge           As           Se           Br           Kr  
      CNST(1.650), CNST(1.727), CNST(1.760), CNST(1.771), CNST(1.749), CNST(1.727), & 
      !Rb           Sr           !>      Y -- Cd        <!                                        
      CNST(1.628), CNST(1.606), CNST(1.639), CNST(1.639),                                                                  & 
      CNST(1.639), CNST(1.639),                                                            & 
      CNST(1.639), CNST(1.639),                                                            & 
      CNST(1.639), CNST(1.639),                                                                  & 
      CNST(1.639), CNST(1.639),                                                                  & 
      !In                  Sn           Sb           Te           I            Xe 
      CNST(2.672), CNST(1.804), CNST(1.881), CNST(1.892), CNST(1.892), CNST(1.881)  / 

    select case (pcm_vdw_type)
    case (PCM_VDW_OPTIMIZED)
      if (species_z(species) > UPTO_XE) then
        write(message(1),'(a,a)') "The van der Waals radius is missing for element ", trim(species_label(species))
        write(message(2),'(a)') "Use PCMVdWRadii = pcm_vdw_species, for other vdw radii values" 
        call messages_fatal(2, namespace=namespace)
      end if
      ia = int(species_z(species))
      vdw_r = vdw_radii(ia)*P_Ang

    case (PCM_VDW_SPECIES)
      vdw_r = species_vdw_radius(species)
      if(vdw_r < M_ZERO) then
        call messages_write('The default vdW radius for species '//trim(species_label(species))//':')
        call messages_write(' is not defined. ')
        call messages_write(' Add a positive vdW radius value in %Species block. ')
        call messages_fatal(namespace=namespace)
      end if
    end select

  end function pcm_get_vdw_radius


  ! -----------------------------------------------------------------------------
  !> Computes the dipole moment mu_pcm due to a distribution of charges q_pcm
  subroutine pcm_dipole(mu_pcm, q_pcm, tess, n_tess)
    FLOAT,               intent(out) :: mu_pcm(:) !< (1:PCM_DIM_SPACE)
    FLOAT,               intent(in)  :: q_pcm(:)  !< (1:n_tess)
    integer,             intent(in)  :: n_tess
    type(pcm_tessera_t), intent(in)  :: tess(:) !< (1:n_tess)

    integer :: ia

    PUSH_SUB(pcm_dipole)

    mu_pcm = M_ZERO
    do ia = 1, n_tess
      mu_pcm = mu_pcm + q_pcm(ia) * tess(ia)%point
    end do

    POP_SUB(pcm_dipole)
  end subroutine pcm_dipole

  ! -----------------------------------------------------------------------------
  !> Computes the field e_pcm at the reference point ref_point due to a distribution of charges q_pcm
  subroutine pcm_field(e_pcm, q_pcm, ref_point, tess, n_tess)
    FLOAT,               intent(out) :: e_pcm(:) !< (1:PCM_DIM_SPACE)
    FLOAT,               intent(in)  :: q_pcm(:)  !< (1:n_tess)
    integer,             intent(in)  :: n_tess
    type(pcm_tessera_t), intent(in)  :: tess(:) !< (1:n_tess)
    FLOAT,               intent(in)  :: ref_point(1:3)

    FLOAT :: diff(1:3)
    FLOAT :: dist

    integer :: ia

    PUSH_SUB(pcm_field)

    e_pcm = M_ZERO
    do ia = 1, n_tess
      diff = ref_point - tess(ia)%point
      dist = sqrt(dot_product(diff, diff))
      e_pcm = e_pcm + q_pcm(ia) * diff / dist**3
    end do

    POP_SUB(pcm_field)
  end subroutine pcm_field

  ! -----------------------------------------------------------------------------
  ! Driver function to evaluate eps(omega)
  subroutine pcm_eps(pcm, eps, omega)
    type(pcm_min_t), intent(in)  :: pcm
    CMPLX,           intent(out) :: eps
    FLOAT,           intent(in)  :: omega

    PUSH_SUB(pcm_eps)

    if (pcm%tdlevel == PCM_TD_EOM) then
      if (pcm%which_eps == PCM_DEBYE_MODEL) then
        call pcm_eps_deb(eps, pcm%deb, omega)
      else if (pcm%which_eps == PCM_DRUDE_MODEL) then
        call pcm_eps_drl(eps, pcm%drl, omega)
      end if
    else if (pcm%tdlevel == PCM_TD_NEQ) then
      eps = pcm%deb%eps_d
    else if (pcm%tdlevel == PCM_TD_EQ) then
      eps = pcm%deb%eps_0
    end if

    POP_SUB(pcm_eps)
  end subroutine pcm_eps

  ! -----------------------------------------------------------------------------

  subroutine pcm_min_input_parsing_for_spectrum(pcm, namespace)
    type(pcm_min_t),   intent(out)   :: pcm
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(pcm_min_input_parsing_for_spectrum)

    ! re-parsing PCM keywords
    call parse_variable(namespace, 'PCMCalculation', .false., pcm%run_pcm)
    call messages_print_stress(stdout, trim('PCM'), namespace=namespace)
    call parse_variable(namespace, 'PCMLocalField', .false., pcm%localf)
    call messages_print_var_value(stdout, "PCMLocalField", pcm%localf)
    if ( pcm%localf ) then
      call messages_experimental("PCM local field effects in the optical spectrum")
      call messages_write('Beware of possible numerical errors in the optical spectrum due to PCM local field effects,')
      call messages_new_line()
      call messages_write('particularly, when static and high-frequency values of the dielectric functions are large')
      call messages_write(' (>~10 in units of the vacuum permittivity \epsilon_0).')
      call messages_new_line()
      call messages_write('However, PCM local field effects in the optical spectrum work well for polar or non-polar solvents')
      call messages_new_line()
      call messages_write('in the nonequilibrium or equation-of-motion TD-PCM propagation schemes.')
      call messages_warning(namespace=namespace)
    end if
    call parse_variable(namespace, 'PCMTDLevel' , PCM_TD_EQ, pcm%tdlevel)
    call messages_print_var_value(stdout, "PCMTDLevel", pcm%tdlevel)

    ! reading dielectric function model parameters
    call parse_variable(namespace, 'PCMStaticEpsilon' , M_ONE, pcm%deb%eps_0)
    call messages_print_var_value(stdout, "PCMStaticEpsilon", pcm%deb%eps_0)
    if (pcm%tdlevel == PCM_TD_EOM) then
      call parse_variable(namespace, 'PCMEpsilonModel', PCM_DEBYE_MODEL, pcm%which_eps)
      call messages_print_var_value(stdout, "PCMEpsilonModel", pcm%which_eps)
      if (pcm%which_eps == PCM_DEBYE_MODEL) then
        call parse_variable(namespace, 'PCMDynamicEpsilon', pcm%deb%eps_0, pcm%deb%eps_d)
        call messages_print_var_value(stdout, "PCMDynamicEpsilon", pcm%deb%eps_d)
        call parse_variable(namespace, 'PCMDebyeRelaxTime', M_ZERO, pcm%deb%tau)
        call messages_print_var_value(stdout, "PCMDebyeRelaxTime", pcm%deb%tau)
      else if (pcm%which_eps == PCM_DRUDE_MODEL) then
        call parse_variable(namespace, 'PCMDrudeLOmega', sqrt(M_ONE/(pcm%deb%eps_0-M_ONE)), pcm%drl%w0)
        call messages_print_var_value(stdout, "PCMDrudeLOmega", pcm%drl%w0)
        call parse_variable(namespace, 'PCMDrudeLDamping', M_ZERO, pcm%drl%gm)
        call messages_print_var_value(stdout, "PCMDrudeLDamping", pcm%drl%gm)
      end if
    else if (pcm%tdlevel == PCM_TD_NEQ) then
      call parse_variable(namespace, 'PCMDynamicEpsilon', pcm%deb%eps_0, pcm%deb%eps_d)
      call messages_print_var_value(stdout, "PCMDynamicEpsilon", pcm%deb%eps_d)
    end if

    POP_SUB(pcm_min_input_parsing_for_spectrum)
  end subroutine pcm_min_input_parsing_for_spectrum

  ! -----------------------------------------------------------------------------
  ! Debye dielectric function
  subroutine pcm_eps_deb(eps, deb, omega)
    CMPLX,               intent(out) :: eps
    type(debye_param_t), intent(in) :: deb
    FLOAT,               intent(in) :: omega

    PUSH_SUB(pcm_eps_deb)

    eps = deb%eps_d +  (deb%eps_0 - deb%eps_d)/(M_ONE + (omega*deb%tau)**2) + &
    M_zI*omega*deb%tau*(deb%eps_0 - deb%eps_d)/(M_ONE + (omega*deb%tau)**2)

    POP_SUB(pcm_eps_deb)
  end subroutine pcm_eps_deb

  ! -----------------------------------------------------------------------------
  ! Drude-Lorentz dielectric function
  subroutine pcm_eps_drl(eps, drl, omega)
    CMPLX,               intent(out) :: eps
    type(drude_param_t), intent(in)  :: drl
    FLOAT,               intent(in)  :: omega

    PUSH_SUB(pcm_eps_drl)

    eps = M_ONE + (drl%w0**2 - omega**2)*drl%aa/((drl%w0**2 - omega**2)**2 + (omega*drl%gm)**2) + &
                   M_zI*omega*drl%gm*drl%aa/((drl%w0**2 - omega**2)**2 + (omega*drl%gm)**2)

    POP_SUB(pcm_eps_drl)
  end subroutine pcm_eps_drl

end module pcm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

