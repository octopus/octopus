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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module species_oct_m
  use element_oct_m
  use global_oct_m
  use iihash_oct_m
  use io_oct_m
  use loct_pointer_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ps_oct_m
  use pseudo_oct_m
  use share_directory_oct_m
  use pseudo_set_oct_m
  use splines_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                        &
    read_from_default_file,        &
    species_t,                     &
    species_init,                  &
    species_read,                  &
    species_build,                 &
    species_init_global,           &
    species_end_global,            &
    species_read_delta,            &
    species_pot_init,              &
    species_type,                  &
    species_label,                 &
    species_index,                 &
    species_has_nlcc,              &
    species_has_density,           &
    species_ps,                    &
    species_zval,                  &
    species_z,                     &
    species_sc_alpha,              &
    species_def_rsize,             &
    species_def_h,                 &
    species_jradius,               &
    species_jthick,                &
    species_hubbard_l,             &
    species_hubbard_u,             &
    species_hubbard_j,             &
    species_hubbard_alpha,         &
    species_sigma,                 &
    species_omega,                 &
    species_mass,                  &
    species_vdw_radius,            &
    species_rho_string,            &
    species_filename,              &
    species_niwfs,                 &
    species_iwf_ilm,               &
    species_iwf_n,                 &
    species_iwf_j,                 &
    species_userdef_pot,           &
    species_is_ps,                 &
    species_is_full,               &
    species_is_local,              &
    species_represents_real_atom,  &
    species_real_nl_projector,     &
    species_nl_projector,          &
    species_get_iwf_radius,        &
    species_get_ps_radius,         &
    species_copy,                  &
    species_end,                   &
    species_x_functional,          &
    species_c_functional

  integer, public, parameter :: LABEL_LEN=15

  integer, public, parameter ::  &
    SPECIES_JELLIUM        = 3,             & !< jellium sphere.
    SPECIES_JELLIUM_SLAB   = 4,             & !< jellium slab.
    SPECIES_JELLIUM_CHARGE_DENSITY = 129,   & !< jellium volume read from file
    SPECIES_FROZEN         = 5,             & !< frozen species.
    SPECIES_PSEUDO         = 7,             & !< pseudopotential
    SPECIES_PSPIO          = 110,           & !< pseudopotential parsed by pspio library
    SPECIES_USDEF          = 123,           & !< user-defined function for local potential
    SPECIES_FULL_GAUSSIAN  = 124,           & !< full-potential atom
    SPECIES_CHARGE_DENSITY = 125,           & !< user-defined function for charge density
    SPECIES_FROM_FILE      = 126,           &
    SPECIES_FULL_DELTA     = 127,           & !< full-potential atom
    SPECIES_SOFT_COULOMB   = 128              !< soft-Coulomb potential

  type species_t
    private
    integer :: index                  !< just a counter

    character(len=LABEL_LEN) :: label !< Identifier for the species
    integer :: type                   !< what type of species
    FLOAT   :: z                      !< charge of the species
    FLOAT   :: z_val                  !< valence charge of the species -- the total charge
                                      !< minus the core charge in the case of the pseudopotentials
    FLOAT   :: mass                   !< mass, in atomic mass units (!= atomic units of mass)
    FLOAT   :: vdw_radius             !< vdw radius, in atomic length units.

    logical :: has_density            !< true if the species has an electronic density


    character(len=1024) :: potential_formula !< for the user-defined potential
    FLOAT :: omega                  !< harmonic frequency for Hermite polynomials


    character(len=MAX_PATH_LEN) :: filename !< for the potential read from a file.


    FLOAT :: jradius              !< jellium stuff
    FLOAT :: jthick               !< jellium stuff

    FLOAT :: sc_alpha                !< the soft-Coulomb parameter

    type(ps_t), pointer :: ps
    logical             :: nlcc   !< true if we have non-local core corrections


    FLOAT :: sigma                !< If we have an all-electron atom:


    character(len=200) :: density_formula !< If we have a charge distribution creating the potential:


    FLOAT :: def_rsize, def_h     !< the default values for the atomic radius and  spacing


    integer :: niwfs              !< The number of initial wavefunctions
    integer, pointer :: iwf_l(:, :), iwf_m(:, :), iwf_i(:, :), iwf_n(:, :) !< i, n, l, m as a function of iorb and ispin
    FLOAT, pointer :: iwf_j(:)    !< j as a function of iorb

    integer :: hubbard_l          !< For the LDA+U, the angular momentum for the applied U
    FLOAT   :: hubbard_U          !< For the LDA+U, the effective U
    FLOAT   :: hubbard_j          !< For the LDA+U, j (l-1/2 or l+1/2)
    FLOAT   :: hubbard_alpha      !< For the LDA+U, a potential contraining the occupations
    integer :: user_lmax          !< For the TM pseudos, user defined lmax 
    integer :: user_llocal        !< For the TM pseudos, used defined llocal
    integer :: pseudopotential_set_id !< to which set this pseudopotential belongs
    logical :: pseudopotential_set_initialized
    type(pseudo_set_t) :: pseudopotential_set
  end type species_t

  interface species_end
    module procedure species_end_species
    module procedure species_end_array
  end interface species_end

  logical :: initialized = .false.
  integer :: default_pseudopotential_set_id
  type(pseudo_set_t) :: default_pseudopotential_set
  real(8) :: energy_tolerance
  logical :: automatic
  
contains

  
  ! ---------------------------------------------------------
  subroutine species_nullify(this)
    type(species_t), intent(out) :: this

    PUSH_SUB(species_nullify)

    this%index=0
    this%label=""
    this%type=0
    this%z = CNST(-1.0)
    this%z_val = CNST(-1.0)
    this%mass = CNST(-1.0)
    this%vdw_radius = CNST(-1.0)
    this%has_density=.false.
    this%potential_formula=""
    this%omega=M_ZERO
    this%filename=""
    this%jradius=M_ZERO
    this%jthick=M_ZERO
    nullify(this%ps)
    this%nlcc=.false.
    this%sigma=M_ZERO
    this%density_formula=""
    this%def_rsize = CNST(-1.0)
    this%def_h = CNST(-1.0)
    this%niwfs=-1
    nullify(this%iwf_l)
    nullify(this%iwf_m)
    nullify(this%iwf_i)
    nullify(this%iwf_n)
    nullify(this%iwf_j)
    this%hubbard_l=-1
    this%hubbard_U=M_ZERO
    this%hubbard_j=M_ZERO
    this%hubbard_alpha = M_ZERO
    this%user_lmax   = INVALID_L
    this%user_llocal = INVALID_L
    this%pseudopotential_set_id = OPTION__PSEUDOPOTENTIALSET__NONE
    this%pseudopotential_set_initialized = .false.
    call pseudo_set_nullify(this%pseudopotential_set)
    
    POP_SUB(species_nullify)
  end subroutine species_nullify


  ! ---------------------------------------------------------
  subroutine species_init_global()
    integer :: ierr
    
    PUSH_SUB(species_init_global)

    initialized = .true.

    call share_directory_set(conf%share)    

    !%Variable PseudopotentialAutomaticParameters
    !%Type logical
    !%Default false
    !%Section System::Species
    !%Description
    !% (Experimental) This enables a new automatic method for
    !% determining the grid parameters for the pseudopotential
    !% (spacing and radius). For the moment, only the spacing can be
    !% adjusted for a few pseudopotentials.
    !%
    !% This does not affect Octopus fixed default parameters for the standard
    !% pseudopotential set.
    !%End
    call parse_variable('PseudopotentialAutomaticParameters', .false., automatic)
    
    if(automatic) call messages_experimental('PseudopotentialAutomaticParameters')
    
    !%Variable PseudopotentialEnergyTolerance
    !%Type float
    !%Default 0.005
    !%Section System::Species
    !%Description
    !% For some pseudopotentials, Octopus can select the grid
    !% spacing automatically so that the discretization error
    !% when calculating the total energy is below a certain
    !% threshold. This variable controls the value of that threshold.
    !% Note that other quantities of interest might require a
    !% different spacing to be considered converged within a similar threshold.
    !%End
    call parse_variable('PseudopotentialEnergyTolerance', CNST(0.005), energy_tolerance)
    
    !%Variable PseudopotentialSet
    !%Type integer
    !%Default standard
    !%Section System::Species
    !%Description
    !% Selects the set of pseudopotentials used by default for species
    !% not defined in the <tt>Species</tt> block.
    !%
    !% These sets of pseudopotentials come from different
    !% sources. Octopus developers have not validated them. We include
    !% them with the code for convenience of the users, but you are
    !% expected to check the quality and suitability of the
    !% pseudopotential for your application.
    !%
    !%Option none 0
    !% Do not load any pseudopotential by default. All species must be
    !% specified in the Species block.
    !%Option standard 1
    !% The standard set of Octopus that provides LDA pseudopotentials
    !% in the PSF format for some elements: H, Li, C, N, O, Na, Si, S, Ti, Se, Cd.
    !%Option sg15 2
    !% The set of Optimized Norm-Conserving Vanderbilt
    !% PBE pseudopotentials. Ref: M. Schlipf and F. Gygi, <i>Comp. Phys. Commun.</i> <b>196</b>, 36 (2015).
    !% This set provides pseudopotentials for elements up to Z = 83
    !% (Bi), excluding Lanthanides.
    !%Option hgh_lda 3
    !% The set of Hartwigsen-Goedecker-Hutter LDA pseudopotentials for elements from H to Rn.
    !% Ref: C. Hartwigsen, S. Goedecker, and J. Hutter, <i>Phys. Rev. B</i> <b>58</b>, 3641 (1998).
    !%Option hgh_lda_sc 31
    !% The semicore set of Hartwigsen-Goedecker-Hutter LDA pseudopotentials.
    !% Ref: C. Hartwigsen, S. Goedecker, and J. Hutter, <i>Phys. Rev. B</i> <b>58</b>, 3641 (1998).
    !%Option hscv_lda 4
    !% (experimental) The set of Hamann-Schlueter-Chiang-Vanderbilt (HSCV) potentials
    !% for LDA exchange and correlation downloaded from http://fpmd.ucdavis.edu/potentials/index.htm.
    !% These pseudopotentials were originally intended for the QBox
    !% code. They were generated using the method of Hamann, Schluter
    !% and Chiang. Ref: D. Vanderbilt, <i>Phys. Rev. B</i> <b>32</b>, 8412 (1985).
    !% Warning from the original site: The potentials provided in this
    !% site are distributed without warranty. In most cases,
    !% potentials were not tested. Potentials should be thoroughly
    !% tested before being used in simulations.
    !%Option hscv_pbe 5
    !% (experimental) PBE version of the HSCV pseudopotentials. Check the
    !% documentation of the option <tt>hscv_lda</tt> for details and warnings.
    !%Option pseudodojo_pbe 100
    !% (experimental) PBE version of the pseudopotentials of http://pseudo-dojo.org. Version 0.4.
    !%Option pseudodojo_pbe_stringent 102
    !% (experimental) High-accuracy PBE version of the pseudopotentials of http://pseudo-dojo.org. Version 0.4.
    !%Option pseudodojo_lda 103
    !% (experimental) LDA pseudopotentials of http://pseudo-dojo.org. Version 0.4.
    !%Option pseudodojo_lda_stringent 104
    !% (experimental) High-accuracy LDA pseudopotentials of http://pseudo-dojo.org. Version 0.4.
    !%Option pseudodojo_pbesol 105
    !% (experimental) PBEsol version of the pseudopotentials of http://pseudo-dojo.org. Version 0.3.
    !%Option pseudodojo_pbesol_stringent 106
    !% (experimental) High-accuracy PBEsol version of the pseudopotentials of http://pseudo-dojo.org. Version 0.4.
    !%End

    call parse_variable('PseudopotentialSet', OPTION__PSEUDOPOTENTIALSET__STANDARD, default_pseudopotential_set_id)
    call messages_print_var_option(stdout, 'PseudopotentialSet', default_pseudopotential_set_id)
    select case (default_pseudopotential_set_id)
    case (OPTION__PSEUDOPOTENTIALSET__NONE)
      call messages_experimental('PseudopotentialSet = none')
    case (OPTION__PSEUDOPOTENTIALSET__SG15)
      call messages_experimental('PseudopotentialSet = sg15')
    case (OPTION__PSEUDOPOTENTIALSET__HSCV_LDA)
      call messages_experimental('PseudopotentialSet = hscv_lda')
    case (OPTION__PSEUDOPOTENTIALSET__HSCV_PBE)
      call messages_experimental('PseudopotentialSet = hscv_pbe')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_LDA)
      call messages_experimental('PseudopotentialSet = pseudodojo_lda')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_LDA_STRINGENT)
      call messages_experimental('PseudopotentialSet = pseudodojo_lda_stringent')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBE)
      call messages_experimental('PseudopotentialSet = pseudodojo_pbe')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBE_STRINGENT)
      call messages_experimental('PseudopotentialSet = pseudodojo_pbe_stringent')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBESOL)
      call messages_experimental('PseudopotentialSet = pseudodojo_pbesol')
    case (OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBESOL_STRINGENT)
      call messages_experimental('PseudopotentialSet = pseudodojo_pbesol_stringent')
    end select
    if(default_pseudopotential_set_id /= OPTION__PSEUDOPOTENTIALSET__NONE) then
      call pseudo_set_init(default_pseudopotential_set, get_set_directory(default_pseudopotential_set_id), ierr, automatic)
    else
      call pseudo_set_nullify(default_pseudopotential_set)
    end if

    POP_SUB(species_init_global)
  end subroutine species_init_global

  ! ---------------------------------------------------------

  subroutine species_end_global()
    PUSH_SUB(species_end_global)

    call pseudo_set_end(default_pseudopotential_set)
    
    POP_SUB(species_end_global)
  end subroutine species_end_global

  ! ---------------------------------------------------------
  !> Initializes a species object. This should be the
  !! first routine to be called (before species_read and species_build).
  !! The label argument must match one of the labels given in the %Species block
  !! in the input file or one of the labels in the defaults file.
  ! ---------------------------------------------------------
  subroutine species_init(this, label, index)
    type(species_t),  intent(out)   :: this
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: index

    PUSH_SUB(species_init)
    
    ASSERT(initialized)

    call species_nullify(this)
    this%label = trim(label)
    this%index = index

    POP_SUB(species_init)
  end subroutine species_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Reads the information (from the input file) about a species_t variable, initializing
  !! part of it (it has to be completed later with "species_build").
  !! Note that species_read has to be called only after species_init has been called.
  ! ---------------------------------------------------------
  subroutine species_read(spec, parser)
    type(species_t), intent(inout) :: spec
    type(parser_t),  intent(in)    :: parser

    character(len=LABEL_LEN)  :: lab
    integer :: ib, row, n_spec_block, read_data
    type(block_t) :: blk

    PUSH_SUB(species_read)

    spec%has_density = .false. ! there is no density associated
    spec%nlcc      = .false.   ! without non-local core corrections
    spec%def_h     = -M_ONE    ! not defined
    spec%def_rsize = -M_ONE    ! not defined
    spec%potential_formula  = ""
    spec%pseudopotential_set_id = OPTION__PSEUDOPOTENTIALSET__NONE
    read_data   = 0

    !%Variable Species
    !%Type block
    !%Section System::Species
    !%Description
    !% A species is by definition either an "ion" (nucleus + core electrons) described
    !% through a pseudopotential, or a model potential.
    !%
    !% Note that some sets of pseudopotentials are distributed with
    !% the code. To use these pseudopotentials, you do not need to define them
    !% explicitly in the <tt>Species</tt> block, as default parameters
    !% are provided.
    !% You can select the set for default pseudopotentials using the
    !% <tt>PseudopotentialSet</tt> variable.
    !%
    !% Additional pseudopotentials can be downloaded from the <a
    !% href='http://octopus-code.org/wiki/Pseudopotentials'>
    !% octopus homepage</a> or from other sources. Supported norm-conserving pseudopotential formats are
    !% detected by the file extension: UPF (<tt>.upf</tt>), PSF (SIESTA, <tt>.psf</tt>), FHI (ABINIT 6, <tt>.fhi</tt>),
    !% CPI (Fritz-Haber, <tt>.cpi</tt>), QSO (quantum-simulation.org, for Qbox, <tt>.xml</tt>),
    !% HGH (Hartwigsen-Goedecker-Hutter, <tt>.hgh</tt>).
    !% PSPIO format can also be used via <tt>species_pspio</tt> if that library is linked.
    !% Note: pseudopotentials may only be used in 3D.
    !%
    !% The format of this block is the following: The first field is a
    !% string that defines the name of the species. The second field
    !% defines the type of species (the valid options are detailed
    !% below).
    !%
    !% Then a list of parameters follows. The parameters are specified
    !% by a first field with the parameter name and the field that
    !% follows with the value of the parameter. Some parameters are
    !% specific to a certain species while others are accepted by all
    !% species. These are <tt>mass</tt>, <tt>max_spacing</tt>, and <tt>min_radius</tt>.
    !%
    !% These are examples of possible species:
    !%
    !% <tt>%Species
    !% <br>&nbsp;&nbsp;'O'       | species_pseudo         | file | 'O.psf' | lmax |  1 | lloc | 1
    !% <br>&nbsp;&nbsp;'H'       | species_pseudo         | file | '../H.hgh'
    !% <br>&nbsp;&nbsp;'Xe'      | species_pseudo         | set | pseudojo_pbe_stringent
    !% <br>&nbsp;&nbsp;'C'       | species_pseudo         | file | "carbon.xml"
    !% <br>&nbsp;&nbsp;'jlm'     | species_jellium        | jellium_radius | 5.0
    !% <br>&nbsp;&nbsp;'rho'     | species_charge_density | density_formula | "exp(-r/a)" | mass | 17.0 | valence | 6
    !% <br>&nbsp;&nbsp;'udf'     | species_user_defined   | potential_formula | "1/2*r^2" | valence | 8
    !% <br>&nbsp;&nbsp;'He_all'  | species_full_delta
    !% <br>&nbsp;&nbsp;'H_all'   | species_full_gaussian  |  gaussian_width |  0.2
    !% <br>&nbsp;&nbsp;'Li1D'    | species_soft_coulomb   |  softening | 1.5 | valence | 3
    !% <br>%</tt>
    !%Option species_pseudo  -7
    !% The species is a pseudopotential. How to get the
    !% pseudopotential can be specified by the <tt>file</tt> or
    !% the <tt>set</tt> parameters. If both are missing, the
    !% pseudopotential will be taken from the <tt>PseudopotentialSet</tt>
    !% specified for the run, this is useful if you want to change
    !% some parameters of the pseudo, like the <tt>mass</tt>.
    !%
    !% The optional parameters for this type of species are
    !% <tt>lmax</tt>, that defines the maximum angular momentum
    !% component to be used, and <tt>lloc</tt>, that defines the
    !% angular momentum to be considered as local. When these
    !% parameters are not set, the value for lmax is the maximum
    !% angular component from the pseudopotential file. The default
    !% value for <tt>lloc</tt> is taken from the pseudopotential if
    !% available, if not, it is set to 0. Note that, depending on the
    !% type of pseudopotential, it might not be possible to select
    !% <tt>lmax</tt> and <tt>lloc</tt>, if that is the case the
    !% parameters will be ignored.
    !%
    !%Option species_pspio  -110
    !% (experimental) Alternative method to read pseudopotentials
    !% using the PSPIO library. This species uses the same parameters
    !% as <tt>species_pseudo</tt>.
    !%Option species_user_defined -123
    !% Species with user-defined potential. The potential for the
    !% species is defined by the formula given by the <tt>potential_formula</tt>
    !% parameter.
    !% The
    !% <tt>valence</tt> parameter determines the number of electrons
    !% associated with the species. By default, a valence of 0 is assumed.
    !%Option species_charge_density -125
    !% The potential for this species is created from the distribution
    !% of charge given by the <tt>density_formula</tt> parameter.
    !% The
    !% <tt>valence</tt> parameter determines the number of electrons
    !% associated with the species. By default, a valence of 0 is assumed.
    !%Option species_point  -3
    !%Option species_jellium  -3
    !% Jellium sphere.
    !% The charge associated with this species must be given by the <tt>valence</tt> parameter.
    !%Option species_jellium_slab  -4
    !% A slab of jellium that extends across the simulation box in the
    !% <i>xy</i>-plane. The dimension along the <i>z</i> direction is
    !% determined by the required parameter <tt>thickness</tt>.
    !% The charge associated with this species must be given by the <tt>valence</tt> parameter.    
    !%Option species_full_delta   -127
    !% Full atomic potential represented by a delta charge
    !% distribution. The atom will be displaced to the nearest grid
    !% point. The atomic number is determined from the name of the species.
    !%Option species_full_gaussian   -124
    !% A full-potential atom is defined by a Gaussian accumulation of
    !% positive charge (distorted if curvilinear coordinates are
    !% used), in the form:
    !%
    !% <math>q(r) = z \beta \exp[ - (\vec{r}-\vec{r_0})^2 / (\sqrt{2} \delta \sigma) ] </math>
    !%
    !% <math>\beta</math> is chosen in order to maintain proper
    !% normalization (the integral of <math>q</math> should sum up to
    !% <math>z</math>). <math>\delta</math> is the grid spacing (the
    !% grid spacing in the first dimension, to be precise).
    !% <math>\vec{r_0}</math> is calculated in such a way that the the
    !% first moment of <math>q(r)/z</math> is equal to the atomic
    !% position. For a precise description, see N. A. Modine,
    !% <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997). The width of the
    !% Gaussian is set by parameter <tt>gaussian_width</tt>. The
    !% atomic number is determined from the name of the species.
    !%Option species_from_file  -126
    !% The potential is read from a file. Accepted file formats, detected by extension: obf, ncdf and csv.
    !% The
    !% <tt>valence</tt> parameter determines the number of electrons
    !% associated with the species. By default, a valence of 0 is assumed.
    !%Option species_soft_coulomb -128
    !% The potential is a soft-Coulomb function, <i>i.e.</i> a function in the form:
    !%
    !% <math>v(r) = - z_{val} / \sqrt{a^2 + r^2}</math>
    !%
    !% The value of <i>a</i> should be given by the mandatory <tt>softening</tt> parameter.
    !% The charge associated with this species must be given by the <tt>valence</tt> parameter.
    !%Option species_jellium_charge_density -129
    !% The parameter is the name of a volume block specifying the shape of the jellium.
    !%Option min_radius -10001
    !% The minimum radius of the box that will be used for this species.
    !%Option max_spacing -10002
    !% The maximum spacing allowed for converged results with this species.
    !%Option lmax -10003
    !% The maximum angular-momentum channel that will be used for the pseudopotential.
    !%Option lloc -10004
    !% The angular-momentum channel of the pseudopotential to be considered local.
    !%Option mass -10005
    !% The mass of the species in atomic mass units, <i>i.e.</i> the mass of a proton is
    !% roughly one. It is set automatically for pseudopotentials from the
    !% <a href=http://www.nist.gov/pml/data/comp.cfm>NIST values</a>.
    !% For other species, the default is 1.0.
    !%Option valence -10006
    !% The number of electrons of the species. It is set automatically for pseudopotentials,
    !% but is mandatory for other species.
    !%Option jellium_radius -10007
    !% The radius of the sphere for <tt>species_jellium</tt>. If this value is not specified,
    !% the default of 0.5 bohr is used.
    !%Option set -10017
    !% For a <tt>species_pseudo</tt>, get the pseudopotential from a
    !% particular set. This flag must be followed with one of the
    !% valid values for the variable <tt>PseudopotentialSet</tt>.
    !%Option gaussian_width -10008
    !% The width of the Gaussian (in units of spacing) used to represent
    !% the nuclear charge for <tt>species_full_gaussian</tt>. If not present,
    !% the default is 0.25.
    !%Option softening -10009
    !% The softening parameter <i>a</i> for <tt>species_soft_coulomb</tt> in units of length.
    !%Option file -10010
    !% The path for the file that describes the species.
    !%Option db_file -10011
    !% Obsolete. Use the <tt>set</tt> option of the <tt>PseudopotentialSet</tt> variable instead.
    !%Option potential_formula -10012
    !% Mathematical expression that defines the potential for <tt>species_user_defined</tt>. You can use
    !% any of the <i>x</i>, <i>y</i>, <i>z</i> or <i>r</i> variables.
    !%Option density_formula -10013
    !% Mathematical expression that defines the charge density for <tt>species_charge_density</tt>. You can use
    !% any of the <i>x</i>, <i>y</i>, <i>z</i> or <i>r</i> variables.
    !%Option thickness -10014
    !% The thickness of the slab for species_jellium_slab. Must be positive.
    !%Option vdw_radius -10015
    !% The van der Waals radius that will be used for this species.
    !%Option volume -10016
    !% Name of a volume block
    !%Option hubbard_l -10018
    !% The angular-momentum for which the effective U will be applied.
    !%Option hubbard_u -10019
    !% The effective U that will be used for the LDA+U calculations.
    !%Option hubbard_j -10020
    !% The value of j (hubbard_l-1/2 or hubbard_l+1/2) on which the effective U is applied.
    !%Option hubbard_alpha -10021
    !% The strength of the potential constraining the occupations of the localized subspace
    !% as defined in PRB 71, 035105 (2005)
    !%End

    call messages_obsolete_variable(parser, 'SpecieAllElectronSigma', 'Species')
    call messages_obsolete_variable(parser, 'SpeciesAllElectronSigma', 'Species')

    ! First, find out if there is a Species block.
    n_spec_block = 0
    if(parse_block(parser, 'Species', blk) == 0) then
      n_spec_block = parse_block_n(blk)
    end if

    ! Find out if the sought species is in the block
    row = -1
    block: do ib = 1, n_spec_block
      call parse_block_string(blk, ib-1, 0, lab)
      if(trim(lab)==trim(spec%label)) then
        row = ib - 1
        exit block
      end if
    end do block

    ! Read whatever may be read from the block
    if(row>=0) then
      call read_from_block(blk, row, spec, read_data)
      call parse_block_end(blk)

      ASSERT(read_data > 0)

      POP_SUB(species_read)
      return
    end if

    ! We get here if there is a Species block but it does not contain
    ! the species we are looking for.
    if(n_spec_block > 0) call parse_block_end(blk)

    spec%pseudopotential_set_id = default_pseudopotential_set_id
    spec%pseudopotential_set = default_pseudopotential_set
    call read_from_set(spec, read_data)

   if(read_data == 0) then
      message(1) = 'Species '//trim(spec%label)//' not found.'
      call messages_fatal(1)
    end if

    POP_SUB(species_read)
  end subroutine species_read

  ! ---------------------------------------------------------

  subroutine read_from_set(spec, read_data)
    type(species_t), intent(inout) :: spec
    integer,         intent(out)   :: read_data

    type(element_t) :: el

    PUSH_SUB(read_from_set)
    
    call element_init(el, get_symbol(spec%label))
    
    if(spec%pseudopotential_set_id /= OPTION__PSEUDOPOTENTIALSET__NONE .and. pseudo_set_has(spec%pseudopotential_set, el)) then
      spec%type = SPECIES_PSEUDO
      spec%filename = pseudo_set_file_path(spec%pseudopotential_set, el)

      ! these might have been set before
      if(spec%z < 0) spec%z = element_atomic_number(el)
      if(spec%user_lmax == INVALID_L) spec%user_lmax = pseudo_set_lmax(spec%pseudopotential_set, el)
      if(spec%user_llocal == INVALID_L) spec%user_llocal = pseudo_set_llocal(spec%pseudopotential_set, el)
      if(spec%def_h < 0) spec%def_h = pseudo_set_spacing(spec%pseudopotential_set, el, energy_tolerance)
      if(spec%def_rsize < 0) spec%def_rsize = pseudo_set_radius(spec%pseudopotential_set, el)
      if(spec%mass < 0) spec%mass = element_mass(el)
      if(spec%vdw_radius < 0) spec%vdw_radius = element_vdw_radius(el)
      read_data = 8
    else
      read_data = 0
    end if

    call element_end(el)

    POP_SUB(read_from_set)
  end subroutine read_from_set

    ! ---------------------------------------------------------

  character(len=MAX_PATH_LEN) function get_set_directory(set_id) result(filename)
    integer,         intent(in)   :: set_id

    PUSH_SUB(get_set_directory)
    
    select case(set_id)
    case(OPTION__PSEUDOPOTENTIALSET__STANDARD)
      filename = trim(conf%share)//'/pseudopotentials/PSF'
    case(OPTION__PSEUDOPOTENTIALSET__SG15)
      filename = trim(conf%share)//'/pseudopotentials/quantum-simulation.org/sg15/'
    case(OPTION__PSEUDOPOTENTIALSET__HGH_LDA)
      filename = trim(conf%share)//'/pseudopotentials/HGH/lda/'
    case(OPTION__PSEUDOPOTENTIALSET__HGH_LDA_SC)
      filename = trim(conf%share)//'/pseudopotentials/HGH/lda_sc/'
    case(OPTION__PSEUDOPOTENTIALSET__HSCV_LDA)
      filename = trim(conf%share)//'/pseudopotentials/quantum-simulation.org/hscv/lda/'
    case(OPTION__PSEUDOPOTENTIALSET__HSCV_PBE)
      filename = trim(conf%share)//'/pseudopotentials/quantum-simulation.org/hscv/pbe/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_LDA)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pw_standard/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_LDA_STRINGENT)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pw_stringent/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBE)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBE_STRINGENT)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_stringent/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBESOL)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbesol_standard/'
    case(OPTION__PSEUDOPOTENTIALSET__PSEUDODOJO_PBESOL_STRINGENT)
      filename = trim(conf%share)//'/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbesol_stringent/'
    case(OPTION__PSEUDOPOTENTIALSET__NONE)
      filename = ''
    case default
      ASSERT(.false.)
    end select

    POP_SUB(get_set_directory)
  end function get_set_directory

  ! ---------------------------------------------------------
  subroutine species_build(spec, parser, ispin, dim, print_info)
    type(species_t),   intent(inout) :: spec
    type(parser_t),    intent(in)    :: parser
    integer,           intent(in)    :: ispin
    integer,           intent(in)    :: dim
    logical, optional, intent(in)    :: print_info

    logical :: print_info_
    integer :: i
    FLOAT   :: pot_re, pot_im, xx(MAX_DIM), rr

    PUSH_SUB(species_build)

    print_info_ = .true.
    if(present(print_info)) then
      print_info_ = print_info
    end if

    ! masses are always in amu, so convert them to a.u.
    spec%mass =  units_to_atomic(unit_amu, spec%mass)

    spec%has_density = .false.

    select case(spec%type)
    case(SPECIES_SOFT_COULOMB)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a soft-Coulomb potential.'
        call messages_info(1)
      end if
      spec%niwfs = species_closed_shell_size(2*nint(spec%z_val+M_HALF))
      spec%omega = CNST(0.1)

    case(SPECIES_PSEUDO, SPECIES_PSPIO)

      ! allocate structure
      SAFE_ALLOCATE(spec%ps)
      if(spec%type == SPECIES_PSPIO) then
        call ps_pspio_init(spec%ps, spec%label, spec%Z, spec%user_lmax, spec%user_llocal, ispin, spec%filename)
      else
        call ps_init(spec%ps, parser, spec%label, spec%Z, spec%user_lmax, spec%user_llocal, ispin, spec%filename)
      end if
      spec%z_val = spec%ps%z_val
      spec%nlcc = spec%ps%nlcc
      spec%niwfs = ps_bound_niwfs(spec%ps)

      ! invalidate these variables as they should not be used after
      spec%user_lmax = INVALID_L
      spec%user_llocal = INVALID_L

    case(SPECIES_USDEF)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a user-defined potential.'
        i = min(237, len_trim(spec%potential_formula)-1) ! I subtract 1 to avoid the non-printable C "end-of-string" character.
        write(message(2),'(a,a)')      '   Potential = ', trim(spec%potential_formula(1:i))
        if(len(trim(spec%potential_formula)) > 237) then
          message(2) = trim(message(2))//'...'
        end if
        call messages_info(2)
      end if
      spec%niwfs = int(max(2*spec%z_val, CNST(1.0)))

      xx    = M_ZERO
      xx(1) = CNST(0.01)
      rr    = sqrt(sum(xx**2))
      call parse_expression(pot_re, pot_im, MAX_DIM, xx, rr, M_ZERO, spec%potential_formula)
      spec%omega = sqrt( abs(M_TWO / CNST(1.0e-4) * pot_re )) ! why...?
      ! To avoid problems with constant potentials.
      if(spec%omega <= M_ZERO) spec%omega = CNST(0.1) 

    case(SPECIES_FROM_FILE)
      if(print_info_) then
        write(message(1),'(a)') 'Species read from file "'//trim(spec%filename)//'".'
        call messages_info(1)
      end if
      spec%niwfs = 2*nint(spec%z_val+M_HALF)
      spec%omega = CNST(0.1)

    case(SPECIES_JELLIUM)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "', trim(spec%label), &
                                       '" is a jellium sphere / approximated point particle.'
        write(message(2),'(a,f11.6)')  '   Valence charge = ', spec%z_val
        write(message(3),'(a,f11.6)')  '   Radius [a.u]   = ', spec%jradius
        write(message(4),'(a,f11.6)')  '   Rs [a.u]       = ', spec%jradius * spec%z_val ** (-M_ONE/M_THREE)
        call messages_info(4)
      end if
      spec%niwfs = species_closed_shell_size(2*nint(spec%z_val+M_HALF))
      spec%omega = CNST(0.1)

    case(SPECIES_JELLIUM_SLAB)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a jellium slab.'
        write(message(2),'(a,f11.6)')  '   Valence charge  = ', spec%z_val
        write(message(3),'(a,f11.6)')  '   Thickness [a.u] = ', spec%jthick
        !write(message(4),'(a,f11.6)')  '   Rs [a.u]       = ', ( M_THREE /( M_FOUR *M_PI ) &
        !& *spec%z_val /( *sb%lsize(1) *sb%lsize(2) ) )**(1.0/3.0) 
        call messages_info(3)
      end if
      spec%niwfs = 2*nint(spec%z_val+M_HALF)
      spec%omega = CNST(0.1)

    case(SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN)
      spec%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is an all-electron atom.'
        write(message(2),'(a,f11.6)')  '   Z = ', spec%z_val
        write(message(3),'(a)')  '   Potential will be calculated solving the Poisson equation'
        write(message(4),'(a)')  '   for a delta density distribution.'
        call messages_info(4)
      end if
      spec%niwfs = species_closed_shell_size(2*nint(spec%z_val+M_HALF))
      spec%omega = spec%z_val

    case(SPECIES_CHARGE_DENSITY, SPECIES_JELLIUM_CHARGE_DENSITY)
      spec%niwfs = int(max(2*spec%z_val, CNST(1.0)))
      spec%omega = spec%z_val
      spec%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "', trim(spec%label), '" is a distribution of charge:'
        if(spec%type == SPECIES_CHARGE_DENSITY) then
          write(message(2),'(a,a)')   '   rho = ', trim(spec%density_formula)
        else
          write(message(2),'(a,a,a)') '   rho is enclosed in volume defined by the "', &
                                      trim(spec%density_formula), '" block'
        end if
        write(message(3),'(a,f11.6)')  '   Z = ', spec%z_val
        call messages_info(3)
      end if
    case default
      call messages_input_error('Species', 'Unknown species type')
    end select

    if(.not. species_is_ps(spec)) then
      ! since there is no real cap, make sure there are at least a few available
      spec%niwfs = max(5, spec%niwfs)
    end if

    SAFE_ALLOCATE(spec%iwf_n(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_l(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_m(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_i(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_j(1:spec%niwfs))

    call species_iwf_fix_qn(spec, ispin, dim)

    if(.not. species_is_ps(spec)) then
      write(message(1),'(a,i6,a,i6)') 'Number of orbitals: ', spec%niwfs
      if(print_info_) call messages_info(1)
      nullify(spec%ps)
    end if

    POP_SUB(species_build)
  end subroutine species_build
  ! ---------------------------------------------------------

  subroutine species_read_delta(spec, zz)
    type(species_t), intent(out) :: spec
    FLOAT,           intent(in)  :: zz

    PUSH_SUB(species_read_delta)

    spec%type = SPECIES_FULL_DELTA
    spec%z     = zz
    spec%z_val = zz
    spec%sigma = CNST(0.25)

    POP_SUB(species_read_delta)
  end subroutine species_read_delta


  ! ---------------------------------------------------------


  !> find size of closed shell for hydrogenic atom with size at least min_niwfs
  integer function species_closed_shell_size(min_niwfs) result(size)
    integer, intent(in) :: min_niwfs

    integer :: nn

    size = 0
    do nn = 1, min_niwfs
      if(size >= min_niwfs) exit
      size = size + nn**2
    end do

  end function species_closed_shell_size

  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> This routine performs some operations on the pseudopotential
  !! functions (filtering, etc), some of which depend on the grid
  !! cutoff value.
  ! ---------------------------------------------------------
  subroutine species_pot_init(this, grid_cutoff, filter)
    type(species_t),     intent(inout) :: this
    FLOAT,               intent(in)    :: grid_cutoff
    integer,             intent(in)    :: filter

    character(len=256) :: dirname
    integer            :: iorb
    FLOAT :: local_radius, orbital_radius

    PUSH_SUB(species_pot_init)
    
    if(species_is_ps(this)) then
      call ps_separate(this%ps)
      
      call ps_getradius(this%ps)

      if(filter /= PS_FILTER_NONE) then 
        call ps_filter(this%ps, filter, grid_cutoff)
        call ps_getradius(this%ps) ! radius may have changed
      end if

      call ps_derivatives(this%ps)

      local_radius = spline_cutoff_radius(this%ps%vl, this%ps%projectors_sphere_threshold)

      orbital_radius = M_ZERO
      ! FIXME: should take max over spins too here.
      do iorb = 1, species_niwfs(this)
        orbital_radius = max(orbital_radius, species_get_iwf_radius(this, this%iwf_i(iorb, 1), is = 1))
      end do

      call messages_write('Info: Pseudopotential for '//trim(this%label), new_line = .true.)
      call messages_write('  Radii for localized parts:', new_line = .true.)
      call messages_write('    local part     = ')
      call messages_write(local_radius, fmt = 'f5.1', units = units_out%length, new_line = .true.)
      call messages_write('    non-local part = ')
      call messages_write(this%ps%rc_max, fmt = 'f5.1', units = units_out%length, new_line = .true.)
      call messages_write('    orbitals       = ')
      call messages_write(orbital_radius, fmt = 'f5.1', units = units_out%length, new_line = .true.)
      call messages_info()

      if(max(local_radius, this%ps%rc_max) > CNST(6.0)) then
        call messages_write("One of the radii of your pseudopotential's localized parts seems", new_line = .true.)
        call messages_write("unusually large; check that your pseudopotential is correct.")
        call messages_warning()
      end if

      if(orbital_radius > CNST(20.0)) then
        call messages_write("The radius of the atomic orbitals given by your pseudopotential seems", new_line = .true.)
        call messages_write("unusually large; check that your pseudopotential is correct.")
        call messages_warning()
      end if

      if(debug%info) then
        write(dirname, '(a)') 'debug/geometry'
        call io_mkdir(dirname)
        call species_debug(trim(dirname), this)
      end if
    end if

    POP_SUB(species_pot_init)
  end subroutine species_pot_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  integer pure function species_type(spec)
    type(species_t), intent(in) :: spec
    species_type = spec%type
  end function species_type
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=LABEL_LEN) pure function species_label(spec)
    type(species_t), intent(in) :: spec
    species_label = trim(spec%label)
  end function species_label
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function species_index(spec)
    type(species_t), intent(in) :: spec
    species_index = spec%index
  end function species_index
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_has_nlcc(spec)
    type(species_t), intent(in) :: spec
    species_has_nlcc = spec%nlcc
  end function species_has_nlcc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_has_density(spec)
    type(species_t), intent(in) :: spec
    species_has_density = spec%has_density
  end function species_has_density
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function species_ps(spec)
    type(ps_t), pointer :: species_ps
    type(species_t), intent(in) :: spec
    species_ps => spec%ps
  end function species_ps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_sc_alpha(spec)
    type(species_t), intent(in) :: spec
    species_sc_alpha = spec%sc_alpha
  end function species_sc_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_zval(spec)
    type(species_t), intent(in) :: spec
    species_zval = spec%z_val
  end function species_zval
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_z(spec)
    type(species_t), intent(in) :: spec
    species_z = spec%z
  end function species_z
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_def_rsize(spec)
    type(species_t), intent(in) :: spec
    species_def_rsize = spec%def_rsize
  end function species_def_rsize
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_def_h(spec)
    type(species_t), intent(in) :: spec
    species_def_h = spec%def_h
  end function species_def_h
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_jradius(spec)
    type(species_t), intent(in) :: spec
    species_jradius = spec%jradius
  end function species_jradius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_jthick(spec)
    type(species_t), intent(in) :: spec
    species_jthick = spec%jthick
  end function species_jthick
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_sigma(spec)
    type(species_t), intent(in) :: spec
    species_sigma = spec%sigma
  end function species_sigma
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_omega(spec)
    type(species_t), intent(in) :: spec
    species_omega = spec%omega
  end function species_omega
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_mass(spec)
    type(species_t), intent(in) :: spec
    species_mass = spec%mass
  end function species_mass
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_vdw_radius(spec)
    type(species_t), intent(in) :: spec
    species_vdw_radius = spec%vdw_radius
  end function species_vdw_radius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=200) pure function species_rho_string(spec)
    type(species_t), intent(in) :: spec
    species_rho_string = trim(spec%density_formula)
  end function species_rho_string
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=200) pure function species_filename(spec)
    type(species_t), intent(in) :: spec
    species_filename = trim(spec%filename)
  end function species_filename
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function species_niwfs(spec)
    type(species_t), intent(in) :: spec
    species_niwfs = spec%niwfs
  end function species_niwfs
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  integer pure function species_hubbard_l(spec)
    type(species_t), intent(in) :: spec
    species_hubbard_l = spec%hubbard_l
  end function species_hubbard_l
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_hubbard_u(spec)
    type(species_t), intent(in) :: spec
    species_hubbard_u = spec%hubbard_u
  end function species_hubbard_u
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  FLOAT pure function species_hubbard_j(spec)
    type(species_t), intent(in) :: spec
    species_hubbard_j = spec%hubbard_j
  end function species_hubbard_j
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  FLOAT pure function species_hubbard_alpha(spec)
    type(species_t), intent(in) :: spec
    species_hubbard_alpha = spec%hubbard_alpha
  end function species_hubbard_alpha
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  pure subroutine species_iwf_ilm(spec, j, is, i, l, m)
    type(species_t), intent(in) :: spec
    integer, intent(in)         :: j, is
    integer, intent(out)        :: i, l, m

    i = spec%iwf_i(j, is)
    l = spec%iwf_l(j, is)
    m = spec%iwf_m(j, is)
  end subroutine species_iwf_ilm
  ! ---------------------------------------------------------

   ! ---------------------------------------------------------
  pure subroutine species_iwf_n(spec, j, is, n)
    type(species_t), intent(in) :: spec
    integer, intent(in)         :: j, is
    integer, intent(out)        :: n

    n = spec%iwf_n(j, is)
  end subroutine species_iwf_n
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  pure subroutine species_iwf_j(spec, iorb, j)
    type(species_t), intent(in) :: spec
    integer, intent(in)         :: iorb
    FLOAT,   intent(out)        :: j

    j = spec%iwf_j(iorb)
  end subroutine species_iwf_j
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  CMPLX function species_userdef_pot(spec, dim, xx, r)
    type(species_t),   intent(in) :: spec
    integer,           intent(in) :: dim
    FLOAT,             intent(in) :: xx(1:MAX_DIM), r

    FLOAT :: pot_re, pot_im

    PUSH_SUB(species_userdef_pot)
    
    call parse_expression(pot_re, pot_im, dim, xx, r, M_ZERO, spec%potential_formula)
    species_userdef_pot = pot_re + M_zI * pot_im  

    POP_SUB(species_userdef_pot)
  end function species_userdef_pot
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_is_ps(spec)
    type(species_t), intent(in) :: spec
    
    species_is_ps = spec%type == SPECIES_PSEUDO .or. spec%type == SPECIES_PSPIO
 
  end function species_is_ps

  ! ---------------------------------------------------------

  integer pure function species_x_functional(spec)
    type(species_t), intent(in) :: spec

    if(species_is_ps(spec)) then
      species_x_functional = spec%ps%exchange_functional

      ! if we do not know, try the pseudpotential set
      if(species_x_functional == PSEUDO_EXCHANGE_UNKNOWN) then
        select case(spec%pseudopotential_set_id)
        case(                                     &
          OPTION__PSEUDOPOTENTIALSET__STANDARD,   &
          OPTION__PSEUDOPOTENTIALSET__HGH_LDA,    &
          OPTION__PSEUDOPOTENTIALSET__HGH_LDA_SC, &
          OPTION__PSEUDOPOTENTIALSET__HSCV_LDA)
          
          species_x_functional = OPTION__XCFUNCTIONAL__LDA_X
          
        case(OPTION__PSEUDOPOTENTIALSET__HSCV_PBE)
          species_x_functional = OPTION__XCFUNCTIONAL__GGA_X_PBE
        end select
      end if
      
    else
      species_x_functional = PSEUDO_EXCHANGE_ANY
    end if

  end function species_x_functional

  ! ---------------------------------------------------------

  integer pure function species_c_functional(spec)
    type(species_t), intent(in) :: spec

    if(species_is_ps(spec)) then
      species_c_functional = spec%ps%correlation_functional

      ! if we do not know, try the pseudpotential set
      if(species_c_functional == PSEUDO_CORRELATION_UNKNOWN) then
        select case(spec%pseudopotential_set_id)
        case(                                     &
          OPTION__PSEUDOPOTENTIALSET__STANDARD,   &
          OPTION__PSEUDOPOTENTIALSET__HGH_LDA,    &
          OPTION__PSEUDOPOTENTIALSET__HGH_LDA_SC, &
          OPTION__PSEUDOPOTENTIALSET__HSCV_LDA)
          
          species_c_functional = OPTION__XCFUNCTIONAL__LDA_C_PZ_MOD/1000
          
        case(OPTION__PSEUDOPOTENTIALSET__HSCV_PBE)
          species_c_functional = OPTION__XCFUNCTIONAL__GGA_C_PBE/1000
        end select
      end if
    else
      species_c_functional = PSEUDO_CORRELATION_ANY
    end if

  end function species_c_functional

  ! ---------------------------------------------------------

  logical elemental function species_is_full(spec)
    type(species_t), intent(in) :: spec
    
    species_is_full = &
         ( spec%type == SPECIES_FULL_GAUSSIAN) .or. &
         ( spec%type == SPECIES_FULL_DELTA)
    
  end function species_is_full

  ! ---------------------------------------------------------

  logical function species_is_local(spec)
    type(species_t), intent(in) :: spec

    PUSH_SUB(species_is_local)

    species_is_local = .true.
      
    if( species_is_ps(spec) ) then
      species_is_local = .false.
      if ( spec%ps%lmax == 0 .and. spec%ps%llocal == 0) species_is_local = .true.
    end if

    POP_SUB(species_is_local)
  end function species_is_local
  ! ---------------------------------------------------------

  logical function species_represents_real_atom(spec)
    type(species_t), intent(in) :: spec
    
    integer :: type
    species_represents_real_atom = .true.

    ! NO PUSH_SUB, called too often
    
    type = species_type(spec)
    species_represents_real_atom =                    &
         type /= SPECIES_USDEF                        &
         .and. type /= SPECIES_CHARGE_DENSITY         &
         .and. type /= SPECIES_FROM_FILE              &
         .and. type /= SPECIES_JELLIUM_CHARGE_DENSITY &
         .and. type /= SPECIES_JELLIUM                &
         .and. type /= SPECIES_JELLIUM_SLAB
    
  end function species_represents_real_atom

  ! ---------------------------------------------------------
  !> This routine returns the non-local projector and its
  !! derivative, built using real spherical harmonics
  subroutine species_real_nl_projector(spec, x, l, lm, i, uV, duV)
    type(species_t),   intent(in)  :: spec
    FLOAT,             intent(in)  :: x(:)
    integer,           intent(in)  :: l, lm, i
    FLOAT,             intent(out) :: uV
    FLOAT, optional,   intent(out) :: duV(:)

    FLOAT :: r, uVr0, duvr0, ylm, gylm(1:3)
    FLOAT, parameter :: ylmconst = CNST(0.488602511902920) !  = sqrt(3/(4*pi))

    ! no push_sub because this function is called very frequently

    ASSERT(species_is_ps(spec))

    r = sqrt(sum(x(1:3)**2))

    uVr0  = spline_eval(spec%ps%kb(l, i), r)

    if(present(duV)) then
      duVr0 = spline_eval(spec%ps%dkb(l, i), r)
      gylm = M_ZERO
      call grylmr(x(1), x(2), x(3), l, lm, ylm, grylm = gylm)
      uv = uvr0*ylm
      if(r >= r_small) then
        duv(1:3) = duvr0*ylm*x(1:3)/r + uvr0*gylm(1:3)
      else
        if(l == 1) then
          duv = M_ZERO
          if(lm == -1) then
            duv(2) = -ylmconst*duvr0
          else if(lm == 0) then
            duv(3) =  ylmconst*duvr0
          else if(lm == 1) then
            duv(1) = -ylmconst*duvr0
          end if
        else
          duv = M_ZERO
        end if
      end if   
    else
      call grylmr(x(1), x(2), x(3), l, lm, ylm)
      uv = uvr0*ylm
    end if

  end subroutine species_real_nl_projector
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> This routine returns the non-local projector, built using
  !! spherical harmonics
  subroutine species_nl_projector(spec, np, x, l, lm, i, uV)
    type(species_t),   intent(in)  :: spec
    integer,           intent(in)  :: np
    FLOAT,             intent(in)  :: x(:,0:) !< (np_part, 3)
    integer,           intent(in)  :: l, lm, i
    CMPLX,             intent(out) :: uV(:) !< (np)

    integer :: ip
    CMPLX :: ylm

    PUSH_SUB(species_nl_projector)

    if(np > 0) then
      uv(1:np) = x(1:np, 0)
      call spline_eval_vec(spec%ps%kb(l, i), np, uv)

      do ip = 1, np
        call ylmr(x(ip, 1), x(ip, 2), x(ip, 3), l, lm, ylm)
        uv(ip) = uv(ip) * ylm
      end do
    end if

    POP_SUB(species_nl_projector)
  end subroutine species_nl_projector
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Return radius outside which orbital is less than threshold value 0.001
  FLOAT function species_get_iwf_radius(spec, ii, is, threshold) result(radius)
    type(species_t),   intent(in) :: spec
    integer,           intent(in) :: ii !< principal quantum number
    integer,           intent(in) :: is !< spin component
    FLOAT, optional,   intent(in) :: threshold

    FLOAT threshold_

    PUSH_SUB(species_get_iwf_radius)

    if(species_is_ps(spec)) then
      threshold_ = optional_default(threshold, spec%ps%projectors_sphere_threshold)
    else
      threshold_ = optional_default(threshold, CNST(0.001))
    end if

    if(species_is_ps(spec)) then
      ASSERT(ii <= spec%ps%conf%p)
      radius = spline_cutoff_radius(spec%ps%ur(ii, is), threshold_)
    else if(species_represents_real_atom(spec)) then
      radius = -ii*log(threshold_)/spec%Z_val
    else
      radius = sqrt(-M_TWO*log(threshold_)/spec%omega)
    end if

    ! The values for hydrogenic and harmonic-oscillator wavefunctions
    ! come from taking the exponential part (i.e. the one that controls
    ! the asymptotic behavior at large r), and setting it equal to
    ! the threshold.

    POP_SUB(species_get_iwf_radius)
  end function species_get_iwf_radius
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Return radius of the pseudopotential if this is a pseudo, zero otherwise
  FLOAT function species_get_ps_radius(spec) result(radius)
    type(species_t),   intent(in) :: spec

    PUSH_SUB(species_get_ps_radius)

    if(species_is_ps(spec)) then
      radius = spec%ps%rc_max
    else
      radius = M_ZERO
    end if

    POP_SUB(species_get_ps_radius)
  end function species_get_ps_radius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_copy(this, that, index)
    type(species_t),         intent(inout) :: this
    type(species_t), target, intent(in)    :: that
    integer,       optional, intent(in)    :: index

    PUSH_SUB(species_copy)

    call species_end(this)
    if(present(index))then
      this%index=index
    else
      this%index=that%index
    end if
    this%label=that%label
    this%type=that%type
    this%z=that%z
    this%z_val=that%z_val
    this%mass=that%mass
    this%vdw_radius=that%vdw_radius
    this%has_density=that%has_density
    this%potential_formula=that%potential_formula
    this%omega=that%omega
    this%filename=that%filename
    this%jradius=that%jradius
    this%jthick=that%jthick
    nullify(this%ps)
    !> To be implemented.
    !> ps_t has no copy procedure.
    !> if(associated(that%ps))then
    !>   allocate(this%ps)
    !>   call ps_copy(this%ps, that%ps)
    !> end if
    this%nlcc=that%nlcc
    this%sigma=that%sigma
    this%density_formula=that%density_formula
    this%def_rsize=that%def_rsize
    this%def_h=that%def_h
    this%niwfs=that%niwfs
    nullify(this%iwf_n, this%iwf_l, this%iwf_m, this%iwf_i)
    call loct_pointer_copy(this%iwf_n, that%iwf_n)
    call loct_pointer_copy(this%iwf_l, that%iwf_l)
    call loct_pointer_copy(this%iwf_m, that%iwf_m)
    call loct_pointer_copy(this%iwf_i, that%iwf_i)
    call loct_pointer_copy(this%iwf_j, that%iwf_j)
    this%hubbard_l=that%hubbard_l
    this%hubbard_U=that%hubbard_U
    this%hubbard_alpha=that%hubbard_alpha
    this%hubbard_j=that%hubbard_j
    this%user_lmax=that%user_lmax
    this%user_llocal=that%user_llocal

    POP_SUB(species_copy)
  end subroutine species_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine species_end_species(spec)
    type(species_t), intent(inout) :: spec
    
    PUSH_SUB(species_end_species)

    if(spec%pseudopotential_set_initialized) call pseudo_set_end(spec%pseudopotential_set)
    
    if (species_is_ps(spec)) then 
      if(associated(spec%ps)) then 
        call ps_end(spec%ps)
        SAFE_DEALLOCATE_P(spec%ps)
      end if
    end if
    SAFE_DEALLOCATE_P(spec%iwf_n)
    SAFE_DEALLOCATE_P(spec%iwf_l)
    SAFE_DEALLOCATE_P(spec%iwf_m)
    SAFE_DEALLOCATE_P(spec%iwf_i)
    SAFE_DEALLOCATE_P(spec%iwf_j)

    POP_SUB(species_end_species)
  end subroutine species_end_species
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine species_end_array(ns, spec)
    integer,         intent(in) :: ns
    type(species_t), pointer    :: spec(:)

    integer :: i

    PUSH_SUB(species_end_array)

    do i = 1, ns
      call species_end_species(spec(i))
    end do

    POP_SUB(species_end_array)
  end subroutine species_end_array
  ! ---------------------------------------------------------





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private procedures

  ! ---------------------------------------------------------
  subroutine species_debug(dir, spec)
    character(len=*), intent(in) :: dir
    type(species_t),  intent(in) :: spec

    character(len=256) :: dirname
    integer :: iunit
    logical :: bool

    if(.not.mpi_grp_is_root(mpi_world)) then
      call messages_debug_newlines(4)
      return
    end if

    PUSH_SUB(species_debug)

    dirname = trim(dir)//'/'//trim(spec%label)

    call io_mkdir(dirname)

    iunit = io_open(trim(dirname)//'/info', action='write')

    write(iunit, '(a,i3)')    'Index  = ', spec%index
    write(iunit, '(2a)')      'Label  = ', trim(spec%label)
    write(iunit, '(a,i3)')    'Type   = ', spec%type
    if (spec%type /= SPECIES_USDEF ) write(iunit, '(a,f15.2)') 'z      = ', spec%z
    if (spec%type == SPECIES_FROM_FILE) then
      write(iunit,'(a)')      'Species read from file "'//trim(spec%filename)//'".'
    end if
    write(iunit, '(a,f15.2)') 'z_val  = ', spec%z_val
    write(iunit, '(a,f15.2)') 'mass = ', spec%mass
    write(iunit, '(a,f15.2)') 'vdw_radius = ', spec%vdw_radius
    bool = species_is_local(spec)
    write(iunit, '(a,l1)')    'local  = ', bool
    write(iunit, '(2a)')      'usdef  = ', trim(spec%potential_formula)
    if (spec%type == SPECIES_JELLIUM) then
      write(iunit, '(a,f15.2)') 'jradius= ', spec%jradius
    end if
    if (spec%type == SPECIES_JELLIUM_SLAB) then
      write(iunit, '(a,f15.2)') 'jthick= ', spec%jthick
    end if
    write(iunit, '(a,l1)')    'nlcc   = ', spec%nlcc
    write(iunit, '(a,f15.2)') 'def_rsize = ', spec%def_rsize
    write(iunit, '(a,f15.2)') 'def_h = ', spec%def_h
    write(iunit, '(a,i3)')    'hubbard_l = ', spec%hubbard_l
    write(iunit, '(a,f15.2)') 'hubbard_U = ', spec%hubbard_U
    write(iunit, '(a,f15.2)') 'hubbard_j = ', spec%hubbard_j
    write(iunit, '(a,f15.2)') 'hubbard_alpha = ', spec%hubbard_alpha

    if(species_is_ps(spec)) then
       if(debug%info) call ps_debug(spec%ps, trim(dirname))
    end if

    call io_close(iunit)
    POP_SUB(species_debug)
  end subroutine species_debug
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine read_from_default_file(iunit, read_data, spec)
    integer,         intent(in)    :: iunit
    integer,         intent(inout) :: read_data
    type(species_t), intent(inout) :: spec

    character(len=LABEL_LEN) :: label
    type(element_t) :: element
    integer :: lmax, llocal
    FLOAT :: zz, spacing, rsize

    PUSH_SUB(read_from_default_file)

    backspace(iunit)

    spec%type = SPECIES_PSEUDO
    
    read(iunit,*) label, spec%filename, zz, lmax, llocal, spacing, rsize

    spec%filename = trim(conf%share)//'/pseudopotentials/'//trim(spec%filename)
    
    ASSERT(trim(label) == trim(spec%label))

    read_data = 8

    ! get the mass, vdw radius and atomic number for this element
    call element_init(element, get_symbol(label))

    ASSERT(element_valid(element))

    ! these might have been set before
    if(spec%z < 0) spec%z = zz
    if(spec%z < 0) spec%z = element_atomic_number(element)
    if(spec%user_lmax == INVALID_L) spec%user_lmax = lmax
    if(spec%user_llocal == INVALID_L) spec%user_llocal = llocal
    if(spec%def_h < 0) spec%def_h = spacing
    if(spec%def_rsize < 0) spec%def_rsize = rsize
    if(spec%mass < 0) spec%mass = element_mass(element)
    if(spec%vdw_radius < 0) spec%vdw_radius = element_vdw_radius(element)
    
    call element_end(element)
    
    POP_SUB(read_from_default_file)
  end subroutine read_from_default_file
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine read_from_block(blk, row, spec, read_data)
    type(block_t),   intent(in)    :: blk
    integer,         intent(in)    :: row
    type(species_t), intent(inout) :: spec
    integer,         intent(out)   :: read_data

    integer :: ncols, icol, flag, set_read_data, ierr
    type(element_t) :: element
    type(iihash_t) :: read_parameters


    PUSH_SUB(read_from_block)

    ncols = parse_block_cols(blk, row)
    read_data = 0

    call parse_block_integer(blk, row, 1, spec%type)

    ! To detect the old species block format, options are represented
    ! as negative values. If we get a non-negative value we know we
    ! are reading a mass.
    if(spec%type >= 0) then
      call messages_write('Found  a species  with the old format.  Please update', new_line = .true.)
      call messages_write('the Species block to the new format, where the second', new_line = .true.)
      call messages_write('column indicates the type of the species.')
      call messages_fatal()
    end if

    ! now we convert back to positive
    spec%type = -spec%type

    read_data = 2

    spec%Z=M_ZERO

    select case(spec%type)

    case(SPECIES_SOFT_COULOMB)

    case(SPECIES_USDEF) ! user-defined

    case(SPECIES_FROM_FILE)

    case(SPECIES_JELLIUM)
      spec%jradius = CNST(0.5)
        
    case(SPECIES_JELLIUM_SLAB)

    case(SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN)
      spec%sigma = CNST(0.25)

    case(SPECIES_CHARGE_DENSITY, SPECIES_JELLIUM_CHARGE_DENSITY)

    case(SPECIES_PSEUDO)

    case(SPECIES_PSPIO) ! a pseudopotential file to be handled by the pspio library

    case default
      call messages_input_error('Species', "Unknown type for species '"//trim(spec%label)//"'")
    end select

    spec%mass = -CNST(1.0)
    spec%vdw_radius = -CNST(1.0)
    spec%z_val = -CNST(1.0)
    spec%sc_alpha = -CNST(1.0)
    spec%jthick = -CNST(1.0)
    
    call iihash_init(read_parameters, 11)
    
    icol = read_data
    do
      if(icol >= ncols) exit

      call parse_block_integer(blk, row, icol, flag)
      
      select case(flag)

      case(OPTION__SPECIES__MIN_RADIUS)
        call check_duplication(OPTION__SPECIES__MIN_RADIUS)
        call parse_block_float(blk, row, icol + 1, spec%def_rsize, unit = units_inp%length)

      case(OPTION__SPECIES__MAX_SPACING)
        call check_duplication(OPTION__SPECIES__MAX_SPACING)
        call parse_block_float(blk, row, icol + 1, spec%def_h, unit = units_inp%length)

      case(OPTION__SPECIES__LMAX)
        call check_duplication(OPTION__SPECIES__LMAX)
        call parse_block_integer(blk, row, icol + 1, spec%user_lmax)

        if(spec%type /= SPECIES_PSEUDO .and. spec%type /= SPECIES_PSPIO) then
          call messages_input_error('Species', &
            "The 'lmax' parameter in species "//trim(spec%label)//" can only be used with pseudopotential species")
        end if
        
        if(spec%user_lmax < 0) then
          call messages_input_error('Species', &
                                    "The 'lmax' parameter in species "//trim(spec%label)//" cannot be negative")
        end if

      case(OPTION__SPECIES__LLOC)
        call check_duplication(OPTION__SPECIES__LLOC)
        call parse_block_integer(blk, row, icol + 1, spec%user_llocal)

        if(spec%type /= SPECIES_PSEUDO .and. spec%type /= SPECIES_PSPIO) then
          call messages_input_error('Species', &
            "The 'lloc' parameter in species "//trim(spec%label)//" can only be used with pseudopotential species")
        end if

        if(spec%user_llocal < 0) then
          call messages_input_error('Species', &
                                    "The 'lloc' parameter in species "//trim(spec%label)//" cannot be negative")
        end if

      case(OPTION__SPECIES__HUBBARD_L)
        call check_duplication(OPTION__SPECIES__HUBBARD_L)
        call parse_block_integer(blk, row, icol + 1, spec%hubbard_l)

        if(spec%type /= SPECIES_PSEUDO .and. spec%type /= SPECIES_PSPIO) then
          call messages_input_error('Species', &
            "The 'hubbard_l' parameter in species "//trim(spec%label)//" can only be used with pseudopotential species")
        end if

        if(spec%hubbard_l < 0) then
          call messages_input_error('Species', &
                                    "The 'hubbard_l' parameter in species "//trim(spec%label)//" cannot be negative")
        end if

     case(OPTION__SPECIES__HUBBARD_U)
        call check_duplication(OPTION__SPECIES__HUBBARD_U)
        call parse_block_float(blk, row, icol + 1, spec%hubbard_u, unit = units_inp%energy)

     case(OPTION__SPECIES__HUBBARD_ALPHA)
        call check_duplication(OPTION__SPECIES__HUBBARD_ALPHA)
        call parse_block_float(blk, row, icol + 1, spec%hubbard_alpha, unit = units_inp%energy)

     case(OPTION__SPECIES__HUBBARD_J)
        call check_duplication(OPTION__SPECIES__HUBBARD_J)
        call parse_block_float(blk, row, icol + 1, spec%hubbard_j)

        if(abs(spec%hubbard_j-spec%hubbard_l) /= M_HALF) then
          call messages_input_error('Species', "The 'hubbard_j' parameter in species "// &
                                    trim(spec%label)//" can only be hubbard_l +/- 1/2")
        end if


      case(OPTION__SPECIES__MASS)
        call check_duplication(OPTION__SPECIES__MASS)
        call parse_block_float(blk, row, icol + 1, spec%mass, unit = units_inp%mass)

      case(OPTION__SPECIES__VALENCE)
        call check_duplication(OPTION__SPECIES__VALENCE)
        call parse_block_float(blk, row, icol + 1, spec%z_val)
        spec%z = spec%z_val

      case(OPTION__SPECIES__JELLIUM_RADIUS)
        call check_duplication(OPTION__SPECIES__JELLIUM_RADIUS)
        call parse_block_float(blk, row, icol + 1, spec%jradius)
        spec%jradius = units_to_atomic(units_inp%length, spec%jradius)
        if(spec%jradius <= M_ZERO) call messages_input_error('Species', 'jellium_radius must be positive')
        if(spec%type /= SPECIES_JELLIUM) then
          call messages_input_error('Species', 'jellium_radius can only be used with species_jellium')
        end if
        
      case(OPTION__SPECIES__GAUSSIAN_WIDTH)
        call check_duplication(OPTION__SPECIES__GAUSSIAN_WIDTH)
        call parse_block_float(blk, row, icol + 1, spec%sigma)
        if(spec%sigma <= M_ZERO) call messages_input_error('Species', 'gaussian_width must be positive')
        if(spec%type /= SPECIES_FULL_GAUSSIAN) then
          call messages_input_error('Species', 'gaussian_width can only be used with species_full_gaussian')
        end if

      case(OPTION__SPECIES__SOFTENING)
        call check_duplication(OPTION__SPECIES__SOFTENING)
        call parse_block_float(blk, row, icol + 1, spec%sc_alpha)
        spec%sc_alpha = units_to_atomic(units_inp%length, spec%sc_alpha)**2
        if(spec%type /= SPECIES_SOFT_COULOMB) then
          call messages_input_error('Species', 'softening can only be used with species_soft_coulomb')
        end if

      case(OPTION__SPECIES__FILE)
        call check_duplication(OPTION__SPECIES__FILE)
        call parse_block_string(blk, row, icol + 1, spec%filename)

      case(OPTION__SPECIES__DB_FILE)
        call messages_write("The 'db_file' option for 'Species' block is obsolete. Please use", new_line = .true.)
        call messages_write("the option 'set' or the variable 'PseudopotentialSet' instead.")
        call messages_fatal()

      case(OPTION__SPECIES__SET)
        call check_duplication(OPTION__SPECIES__SET)
        call parse_block_integer(blk, row, icol + 1, spec%pseudopotential_set_id)
        spec%pseudopotential_set_initialized = .true.
        call pseudo_set_init(spec%pseudopotential_set, get_set_directory(spec%pseudopotential_set_id), ierr, automatic)
        
      case(OPTION__SPECIES__POTENTIAL_FORMULA)
        call check_duplication(OPTION__SPECIES__POTENTIAL_FORMULA)
        call parse_block_string(blk, row, icol + 1, spec%potential_formula)
        call conv_to_C_string(spec%potential_formula)

        if(spec%type /= SPECIES_USDEF) then
          call messages_input_error('Species', 'potential_formula can only be used with species_user_defined')
        end if

      case(OPTION__SPECIES__VOLUME)
        call check_duplication(OPTION__SPECIES__VOLUME)
        call parse_block_string(blk, row, icol + 1, spec%density_formula)
        call conv_to_C_string(spec%density_formula)

        if(spec%type /= SPECIES_JELLIUM_CHARGE_DENSITY) then
          call messages_input_error('Species', 'volume can only be used with species_jellium_charge_density')
        end if

      case(OPTION__SPECIES__DENSITY_FORMULA)
        call check_duplication(OPTION__SPECIES__DENSITY_FORMULA)
        call parse_block_string(blk, row, icol + 1, spec%density_formula)
        call conv_to_C_string(spec%density_formula)
              
        if(spec%type /= SPECIES_CHARGE_DENSITY) then
          call messages_input_error('Species', 'density_formula can only be used with species_charge_density')
        end if

      case(OPTION__SPECIES__THICKNESS)
        call check_duplication(OPTION__SPECIES__THICKNESS)
        call parse_block_float(blk, row, icol + 1, spec%jthick) ! thickness of the jellium slab

        if(spec%jthick <= M_ZERO) call messages_input_error('Species', &
          'the value of the thickness parameter in species '//trim(spec%label)//' must be positive.')

        spec%jthick = units_to_atomic(units_inp%length, spec%jthick) ! units conversion

        if(spec%type /= SPECIES_JELLIUM_SLAB) then
          call messages_input_error('Species', 'thickness can only be used with species_jellium_slab')
        end if
        
      case(OPTION__SPECIES__VDW_RADIUS)
        call check_duplication(OPTION__SPECIES__VDW_RADIUS)
        call parse_block_float(blk, row, icol + 1, spec%vdw_radius, unit = units_inp%length)

      case default
        call messages_input_error('Species', "Unknown parameter in species '"//trim(spec%label)//"'")
        
      end select

      icol = icol + 2        
    end do
    ! CHECK THAT WHAT WE PARSED MAKES SENSE
    
    if(spec%type == SPECIES_SOFT_COULOMB .and. .not. parameter_defined(OPTION__SPECIES__SOFTENING)) then
      call messages_input_error('Species', &
        "The 'softening' parameter is missing for species "//trim(spec%label))
    end if

    if(spec%type == SPECIES_USDEF .and. .not. parameter_defined(OPTION__SPECIES__POTENTIAL_FORMULA)) then
      call messages_input_error('Species', &
        "The 'potential_formula' parameter is missing for species '"//trim(spec%label)//"'")
    end if

    if(spec%type == SPECIES_CHARGE_DENSITY .and. .not. parameter_defined(OPTION__SPECIES__DENSITY_FORMULA)) then
      call messages_input_error('Species', &
        "The 'density_formula' parameter is missing for species '"//trim(spec%label)//"'")
    end if
    
    if(spec%type == SPECIES_FROM_FILE &
      .and. .not. (parameter_defined(OPTION__SPECIES__FILE) .or. parameter_defined(OPTION__SPECIES__DB_FILE))) then
      call messages_input_error('Species', &
        "The 'file' or 'db_file' parameter is missing for species '"//trim(spec%label)//"'")
    end if

    if(spec%type == SPECIES_JELLIUM_SLAB .and. .not. parameter_defined(OPTION__SPECIES__THICKNESS)) then
      call messages_input_error('Species', &
        "The 'thickness' parameter is missing for species '"//trim(spec%label)//"'")
    end if

    if(spec%type == SPECIES_JELLIUM_CHARGE_DENSITY .and. .not. parameter_defined(OPTION__SPECIES__VOLUME)) then
      call messages_input_error('Species', &
        "The 'volume' parameter is missing for species '"//trim(spec%label)//"'")
    end if

    if(parameter_defined(OPTION__SPECIES__LMAX) .and. parameter_defined(OPTION__SPECIES__LLOC)) then
      if(spec%user_llocal > spec%user_lmax) then
        call messages_input_error('Species', &
          "the 'lloc' parameter cannot be larger than the 'lmax' parameter in species "//trim(spec%label))
      end if
    end if
    
    select case(spec%type)
    case(SPECIES_PSEUDO, SPECIES_PSPIO, SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN)

      if( (spec%type == SPECIES_PSEUDO .or. spec%type == SPECIES_PSPIO) &
        .and. .not. (parameter_defined(OPTION__SPECIES__FILE) .or. parameter_defined(OPTION__SPECIES__DB_FILE))) then
        ! we need to read the species from the pseudopotential set

        !if the set was not defined, use the default set
        if(.not. parameter_defined(OPTION__SPECIES__SET)) then
          spec%pseudopotential_set_id = default_pseudopotential_set_id
          spec%pseudopotential_set = default_pseudopotential_set
        end if

        call read_from_set(spec, set_read_data)

        if(set_read_data == 0) then
          call messages_write('Species '//trim(spec%label)//' is not defined in the requested pseudopotential set.')
          call messages_fatal()
        end if
        
      end if

      call element_init(element, get_symbol(spec%label))
      
      if(.not. element_valid(element)) then
        call messages_write('Cannot determine the element for species '//trim(spec%label)//'.')
        call messages_fatal()
      end if

      spec%z = element_atomic_number(element)

      if(spec%type == SPECIES_FULL_DELTA .or. spec%type == SPECIES_FULL_GAUSSIAN) then
        spec%z_val = spec%z
      end if
        
      if(spec%mass < CNST(0.0)) then
        spec%mass = element_mass(element)
        call messages_write('Info: default mass for species '//trim(spec%label)//':')
        call messages_write(spec%mass)
        call messages_write(' amu.')
        call messages_info()
      end if
        
      if(spec%vdw_radius < CNST(0.0)) then
        spec%vdw_radius = element_vdw_radius(element)
        if(spec%vdw_radius < CNST(0.0)) then
          spec%vdw_radius = CNST(0.0)
          call messages_write("The default vdW radius for species '"//trim(spec%label)//"' is not defined.", &
                              new_line = .true.)
          call messages_write("You can specify the vdW radius in %Species block.")
          call messages_warning()
        end if
        call messages_write('Info: default vdW radius for species '//trim(spec%label)//':')
        call messages_write(spec%vdw_radius)
        call messages_write(' [b]')
        call messages_info()
      end if

      call element_end(element)

    case default
      if(.not. parameter_defined(OPTION__SPECIES__MASS)) then
        spec%mass = M_ONE
        call messages_write('Info: default mass for species '//trim(spec%label)//':')
        call messages_write(spec%mass)
        call messages_write(' amu.')
        call messages_info()
      end if

      if(.not. parameter_defined(OPTION__SPECIES__VDW_RADIUS)) then
        spec%vdw_radius = M_ZERO
        call messages_write('Info: default mass for species '//trim(spec%label)//':')
        call messages_write(spec%vdw_radius)
        call messages_write(' [b]')
        call messages_info()
      end if

      if(.not. parameter_defined(OPTION__SPECIES__VALENCE)) then
        if(spec%type == SPECIES_USDEF .or. spec%type == SPECIES_CHARGE_DENSITY .or. &
          spec%type == SPECIES_FROM_FILE) then
          spec%z_val = CNST(0.0)
        else
          call messages_input_error('Species', &
            "The 'valence' parameter is missing for species '"//trim(spec%label)//"'")
        end if
      end if
      
    end select

    call iihash_end(read_parameters)

    POP_SUB(read_from_block)

  contains

    logical function parameter_defined(param) result(defined)
      integer(8), intent(in) :: param

      integer :: tmp
      
      PUSH_SUB(read_from_block.parameter_defined)

      tmp = iihash_lookup(read_parameters, int(-param), defined)
      
      POP_SUB(read_from_block.parameter_defined)
    end function parameter_defined

    !------------------------------------------------------
    
    subroutine check_duplication(param)
      integer(8), intent(in) :: param

      PUSH_SUB(read_from_block.check_duplication)

      if(parameter_defined(param)) then
        call messages_input_error('Species', "Duplicated parameter in species '"//trim(spec%label)//"'")
      end if

      call iihash_insert(read_parameters, int(-param), 1)

      POP_SUB(read_from_block.check_duplication)
    end subroutine check_duplication
    
  end subroutine read_from_block
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> set up quantum numbers of orbitals
  subroutine species_iwf_fix_qn(spec, ispin, dim)
    type(species_t), intent(inout) :: spec
    integer,         intent(in)    :: ispin
    integer,         intent(in)    :: dim

    integer :: is, n, i, l, m, n1, n2, n3

    PUSH_SUB(species_iwf_fix_qn)

    if(species_is_ps(spec)) then
      
      do is = 1, ispin
        n = 1
        do i = 1, spec%ps%conf%p
          if(n > spec%niwfs) exit          
          l = spec%ps%conf%l(i)

          if(.not. spec%ps%bound(i,is)) cycle
          
          do m = -l, l
            spec%iwf_i(n, is) = i
            spec%iwf_n(n, is) = spec%ps%conf%n(i)
            spec%iwf_l(n, is) = l
            spec%iwf_m(n, is) = m
            spec%iwf_j(n) = spec%ps%conf%j(i)
            n = n + 1
          end do
          
        end do
      end do

    else if(species_represents_real_atom(spec) .and. dim == 3) then

      do is = 1, ispin
        n = 1
        ! just up to the highest principal quantum number, actually
        do i = 1, spec%niwfs
          if(n > spec%niwfs) exit
          do l = 0, i-1
            do m = -l, l
              spec%iwf_i(n, is) = i
              spec%iwf_n(n, is) = i
              spec%iwf_l(n, is) = l
              spec%iwf_m(n, is) = m
              spec%iwf_j(n) = M_ZERO
              n = n + 1
            end do
          end do
        end do
      end do

    else

      select case(dim)
      case(1)
        do is = 1, ispin
          do i = 1, spec%niwfs
            spec%iwf_i(i, is) = i
            spec%iwf_n(i, is) = 0
            spec%iwf_l(i, is) = 0
            spec%iwf_m(i, is) = 0
            spec%iwf_j(i) = M_ZERO
          end do
        end do

      case(2)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1
          do
            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1 
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = 0
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = 0
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = 0
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1
          end do
        end do

      case(3)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1; n3 = 1
          do
            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = 0
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3+1
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = n3
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3+1
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_n(i, is) = 1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = n3+1
            spec%iwf_j(i) = M_ZERO
            i = i + 1; if(i>spec%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1; n3 = n3 + 1
          end do
        end do
      end select
    end if

    POP_SUB(species_iwf_fix_qn)
  end subroutine species_iwf_fix_qn

  ! ---------------------------------------------------------

  character(len=LABEL_LEN) function get_symbol(label) result(symbol)
    character(len=*), intent(in)    :: label
    
    integer :: ilend

    ! use only the first part of the label to determine the element
    do ilend = 1, len(label)
      if( iachar(label(ilend:ilend)) >= iachar('a') .and. iachar(label(ilend:ilend)) <= iachar('z') ) cycle
      if( iachar(label(ilend:ilend)) >= iachar('A') .and. iachar(label(ilend:ilend)) <= iachar('Z') ) cycle
      exit
    end do
    ilend = ilend - 1

    symbol = label(1:ilend)
    
  end function get_symbol
      
  
end module species_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
