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
!! $Id$

#include "global.h"

module species_m
  use element_m
  use global_m
  use io_m
  use json_m
  use loct_m
  use loct_math_m
  use loct_pointer_m
  use logrid_m
  use math_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use ps_m
  use space_m
  use splines_m
  use string_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                        &
    read_from_default_file,        &
    species_t,                     &
    species_init,                  &
    species_read,                  &
    species_build,                 &
    species_init_global,           &
    species_read_delta,            &
    species_pot_init,              &
    species_init_from_data_object, &
    species_create_data_object,    &
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
    species_sigma,                 &
    species_omega,                 &
    species_mass,                  &
    species_rho_string,            &
    species_filename,              &
    species_niwfs,                 &
    species_iwf_ilm,               &
    species_userdef_pot,           &
    species_is_ps,                 &
    species_is_full,               &
    species_is_local,              &
    species_represents_real_atom,  &
    species_real_nl_projector,     &
    species_nl_projector,          &
    species_get_iwf_radius,        &
    species_copy,                  &
    species_end

  integer, public, parameter :: LABEL_LEN=15

  integer, public, parameter ::  &
    SPECIES_JELLIUM        = 3,             & !< jellium sphere.
    SPECIES_JELLIUM_SLAB   = 4,             & !< jellium slab.
    SPECIES_FROZEN         = 5,             & !< frozen species.
    SPECIES_PSEUDO         = 7,             & !< pseudopotential
    SPECIES_PSPIO          = 110,           & !< pseudopotential parsed by pspio library
    SPECIES_USDEF          = 123,           & !< user-defined function for local potential
    SPECIES_FULL_GAUSSIAN  = 124,           & !< full-potential atom
    SPECIES_CHARGE_DENSITY = 125,           & !< user-defined function for charge density
    SPECIES_FROM_FILE      = 126,           &
    SPECIES_FULL_DELTA     = 127,           & !< full-potential atom
    SPECIES_SOFT_COULOMB   = 128              !< soft-Coulomb potential

  integer, public, parameter ::             &
    PSEUDO_SET_STANDARD = 1,                &
    PSEUDO_SET_SG15     = 2,                &
    PSEUDO_SET_HGH      = 3

  type species_t
    private
    integer :: index                  !< just a counter

    character(len=LABEL_LEN) :: label !< Identifier for the species
    integer :: type                   !< what type of species
    FLOAT   :: z                      !< charge of the species
    FLOAT   :: z_val                  !< valence charge of the species -- the total charge
                                      !< minus the core charge in the case of the pseudopotentials
    FLOAT   :: mass                 !< mass, in atomic mass units (!= atomic units of mass)

    logical :: has_density            !< true if the species has an electronic density


    character(len=1024) :: user_def !< for the user-defined potential
    FLOAT :: omega                  !< harmonic frequency for Hermite polynomials


    character(len=MAX_PATH_LEN) :: filename !< for the potential read from a file.


    FLOAT :: jradius              !< jellium stuff
    FLOAT :: jthick               !< jellium stuff

    FLOAT :: sc_alpha                !< the soft-Coulomb parameter

    type(ps_t), pointer :: ps
    logical             :: nlcc   !< true if we have non-local core corrections


    FLOAT :: sigma                !< If we have an all-electron atom:


    character(len=200) :: rho     !< If we have a charge distribution creating the potential:


    FLOAT :: def_rsize, def_h     !< the default values for the spacing and atomic radius


    integer :: niwfs              !< The number of initial wavefunctions
    integer, pointer :: iwf_l(:, :), iwf_m(:, :), iwf_i(:, :) !< i, l, m as a function of iorb and ispin

    integer :: lmax, lloc         !< For the TM pseudos, the lmax and lloc.
  end type species_t

  interface species_end
    module procedure species_end_species
    module procedure species_end_array
  end interface species_end

  logical :: initialized = .false.
  integer :: pseudo_set
  
contains

  
  ! ---------------------------------------------------------
  subroutine species_nullify(this)
    type(species_t), intent(out) :: this

    PUSH_SUB(species_nullify)

    this%index=0
    this%label=""
    this%type=0
    this%z=M_ZERO
    this%z_val=M_ZERO
    this%mass=M_ZERO
    this%has_density=.false.
    this%user_def=""
    this%omega=M_ZERO
    this%filename=""
    this%jradius=M_ZERO
    this%jthick=M_ZERO
    nullify(this%ps)
    this%nlcc=.false.
    this%sigma=M_ZERO
    this%rho=""
    this%def_rsize=M_ZERO
    this%def_h=-M_ONE
    this%niwfs=-1
    nullify(this%iwf_l)
    nullify(this%iwf_m)
    nullify(this%iwf_i)
    this%lmax=0
    this%lloc=0

    POP_SUB(species_nullify)
  end subroutine species_nullify


  ! ---------------------------------------------------------
  subroutine species_init_global()
    PUSH_SUB(species_init_global)

    initialized = .true.
    
    !%Variable PseudopotentialSet
    !%Type flag
    !%Default standard
    !%Section System::Species
    !%Description
    !% Selects the set of pseudopotentials used by default.
    !%
    !%Option standard 1
    !% The standard set of Octopus that provides LDA pseudopotentials
    !% for some elements: H, Li, C, N, O, Na, Si, S, Ti, Se, Cd.
    !%Option sg15 2
    !% (experimental) The set of Optimized Norm-Conserving Vanderbilt
    !% PBE pseudopotentials developed by Schlipf and Gygi (M. Schlipf
    !% and F. Gygi, (2015) arXiv:1502.00995). This set provides
    !% pseudopotentials for most elements.
    !%Option hgh 3
    !% (experimental) The set of Hartwigsen-Goedecker-Hutter pseudopotentials.
    !%End

    call parse_integer('PseudopotentialSet', PSEUDO_SET_STANDARD, pseudo_set)
    call messages_print_var_option(stdout, 'PseudopotentialSet', pseudo_set)
    if(pseudo_set == PSEUDO_SET_SG15) call messages_experimental('PseudopotentialSet = sg15')
    if(pseudo_set == PSEUDO_SET_HGH) call messages_experimental('PseudopotentialSet = hgh')

    POP_SUB(species_init_global)
  end subroutine species_init_global
  
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
    
    if(.not. initialized) call species_init_global()

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
  subroutine species_read(spec)
    type(species_t), intent(inout) :: spec

    character(len=256) :: fname
    character(len=LABEL_LEN)  :: lab
    integer :: ib, ispec, row, n_spec_block, n_spec_def, iunit, read_data
    type(block_t) :: blk

    PUSH_SUB(species_read)

    spec%has_density = .false. ! there is no density associated
    spec%nlcc      = .false.   ! without non-local core corrections
    spec%def_h     = -M_ONE    ! not defined
    spec%def_rsize = -M_ONE    ! not defined
    spec%user_def  = ""
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
    !%
    !% You can select the set for default pseudopotentials using the
    !% <tt>PseudopotentialSet</tt> variable.
    !%
    !% Additional pseudopotentials can be downloaded from the <a
    !% href='http://www.tddft.org/programs/octopus/wiki/index.php/Pseudopotentials'>
    !% octopus homepage</a> or from other sources. As explained below,
    !% several pseudopotential formats are supported.
    !%
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
    !% species. These are 'mass', 'max_spacing', and 'max_radius'.
    !%
    !% This is an example of possible species:
    !%
    !% <tt>%Species
    !% <br>&nbsp;&nbsp;'O'       | species_pseudo         | file | 'O.psf' | lmax |  1 | lloc | 1
    !% <br>&nbsp;&nbsp;'H'       | species_pseudo         | file | '../H.hgh'
    !% <br>&nbsp;&nbsp;'Xe'      | species_pseudo         | mass | 131.29 | file | db_file | "UPF/Xe.UPF"
    !% <br>&nbsp;&nbsp;'C'       | species_pseudo         | file | "carbon.xml"
    !% <br>&nbsp;&nbsp;'jlm'     | species_jellium        | jellium_radius | 5.0
    !% <br>&nbsp;&nbsp;'rho'     | species_charge_density | "exp(-r/a)" | mass | 17.0 | valence | 6
    !% <br>&nbsp;&nbsp;'udf'     | species_user_defined   | "1/2*r^2" | valence | 8
    !% <br>&nbsp;&nbsp;'He_all'  | species_full_delta
    !% <br>&nbsp;&nbsp;'H_all'   | species_full_gaussian  |  gaussian_width |  0.2
    !% <br>&nbsp;&nbsp;'Li1D'    | species_soft_coulomb   |  softening | 1.5 | valence | 3
    !% <br>%</tt>
    !%
    !%Option species_pseudo  -7
    !% The species is a pseudopotential. The pseudopotential file must
    !% be defined by the 'file' or 'db_file' paramaters. Optional
    !% arguments are 'lmax' and 'lloc'.
    !%Option species_user_defined -123
    !% Species with user-defined potential. In this case, the fifth
    !% field is a string with a mathematical expression that defines the
    !% potential (you can use any of the <i>x</i>, <i>y</i>, <i>z</i>
    !% or <i>r</i> variables).
    !%Option species_point  -3
    !%Option species_jellium  -3
    !% Jellium sphere: the optional fifth field is the radius of the sphere (default = 0.5 a.u.).
    !%Option species_jellium_slab  -4
    !% Jellium slab: the fifth field is the thickness of the slab.
    !% The slab extends across the simulation box in the <i>xy</i>-plane.
    !%Option species_pspio  -110
    !% (experimental) PSPIO library: the pseudopotential will be read from a file,
    !% either in the working directory or in the <tt>OCTOPUS-HOME/share/pseudopotentials/UPF</tt> 
    !% directory, using the PSPIO library.
    !% No extra columns, as the maximum <i>l</i>-component of the pseudopotential to
    !% consider in the calculation and the <i>l</i>-component to consider as
    !% local are indicated in the pseudopotential file are cannot be changed.
    !%Option species_full_delta   -127
    !% Full atomic potential represented by a delta charge
    !% distribution. The atom will be displaced to the nearest grid
    !% point. No extra columns.
    !%Option species_full_gaussian   -124
    !% A full-potential atom is defined by a Gaussian accumulation of
    !% positive charge (distorted if curvilinear coordinates are
    !% used), in the form:
    !%
    !% <math>
    !% q(r) = z * \beta * exp[ - (\vec{r}-\vec{r0})**2 / (sqrt(2) * \delta * \sigma) ]
    !% </math>
    !%
    !% <math>\beta</math> is chosen in order to maintain proper
    !% normalization (the integral of <math>q</math> should sum up to
    !% <math>z</math>). <math>\delta</math> is the grid spacing (the
    !% grid spacing in the first dimension, to be precise).
    !% <math>\vec{r0}</math> is calculated in such a way that the the
    !% first moment of <math>q(r)/z</math> is equal to the atomic
    !% position. For a precise description, see N. A. Modine,
    !% <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997).
    !%
    !% Column 5 is <math>sigma</math>, the width of the Gaussian that should be
    !% small, but you may run into numerical difficulties if it is too
    !% small (0.25 by default).
    !%Option species_charge_density -125
    !% The potential is created by a distribution of charge.
    !% Column 5 is an expression for the charge distribution.
    !%Option species_from_file  -126
    !% The potential is read from a file, whose name is given in column 5.
    !% Accepted file formats, detected by extension: obf, ncdf and csv.
    !%Option species_soft_coulomb -128
    !% The potential is a soft-Coulomb function, <i>i.e.</i> a function in the form:
    !%
    !% <math>
    !% v(r) = - z_val / sqrt(a^2 + r^2)
    !% </math>
    !%
    !% The value of a should be given by the mandatory 'softening' parameter.
    !%Option min_radius -10001
    !% The minimum radius of the box required for this species to be properly converged.
    !%Option max_spacing -10002
    !% The maximum spacing required for this species for converged results.
    !%Option lmax -10003
    !% The maximum angular momentum of the pseudopotential.
    !%Option lloc -10004
    !% The component of the pseudopotenital to be considered local.
    !%Option mass -10005
    !% The mass of the species in atomic mass units, <i>i.e.</i> the mass of a proton is
    !% roughly one.
    !%Option valence -10006
    !% The number of electrons of the species.
    !%Option jellium_radius -10007
    !% The radius of jellium sphere. If this value is not specified,
    !% the default of 0.5 bohr is used.
    !%Option gaussian_width -10008
    !% The width of the gaussian in units of spacing used to represent
    !% the nuclear charge for species_full_gaussian. If not present,
    !% the default is 0.25.
    !%Option softening -10009
    !% The softening parameter a for species_soft_coulomb in units of length.
    !%Option file -10010
    !% The path for the file that describes the species.
    !%Option db_file -10011
    !% The path for the file, in the Octopus directory of
    !% pseudopotentials, that describes the species.
    !%End

    call messages_obsolete_variable('SpecieAllElectronSigma', 'Species')
    call messages_obsolete_variable('SpeciesAllElectronSigma', 'Species')

    ! First, find out if there is a Species block.
    n_spec_block = 0
    if(parse_block('Species', blk) == 0) then
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

    ! Find out if the species is in the pseudo potential set
    select case(pseudo_set)
    case(PSEUDO_SET_STANDARD)
      fname = trim(conf%share)//'/pseudopotentials/standard.set'
    case(PSEUDO_SET_SG15)
      fname = trim(conf%share)//'/pseudopotentials/sg15.set'
    case(PSEUDO_SET_HGH)
      fname = trim(conf%share)//'/pseudopotentials/hgh.set'
    case default
      ASSERT(.false.)
    end select
      
    n_spec_def = max(0, loct_number_of_lines(fname))
    if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

    iunit = io_open(fname, action='read', status='old', die=.false.)

    if(iunit > 0) then
      read(iunit,*)

      default_file: do ispec = 1, n_spec_def
        read(iunit,*) lab
        if(trim(lab) == trim(spec%label)) then
          call read_from_default_file(iunit, read_data, spec)
          exit default_file
        end if
      end do default_file

      call io_close(iunit)

    else

      call messages_write('Cannot open the octopus internal file:', new_line = .true.)
      call messages_write(" '"//trim(fname)//"'", new_line = .true.)
      call messages_write('There is something wrong with your octopus installation.')
      call messages_fatal()
      
    end if

    if(read_data == 0) then
      message(1) = 'Species '//trim(spec%label)//' not found.'
      call messages_fatal(1)
    end if

    POP_SUB(species_read)
  end subroutine species_read
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_build(spec, ispin, dim, print_info)
    type(species_t),   intent(inout) :: spec
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
        call ps_pspio_init(spec%ps, spec%Z, spec%lmax, spec%lloc, ispin, spec%filename)
      else
        call ps_init(spec%ps, spec%label, spec%Z, spec%lmax, spec%lloc, ispin, spec%filename)
      endif
      spec%z_val = spec%ps%z_val
      spec%nlcc = spec%ps%nlcc
      spec%niwfs = ps_niwfs(spec%ps)

    case(SPECIES_USDEF)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a user-defined potential.'
        i = min(237, len_trim(spec%user_def)-1) ! I subtract 1 to avoid the non-printable C "end-of-string" character.
        write(message(2),'(a,a)')      '   Potential = ', trim(spec%user_def(1:i))
        if(len(trim(spec%user_def)) > 237) then
          message(2) = trim(message(2))//'...'
        end if
        call messages_info(2)
      end if
      spec%niwfs = int(max(2*spec%z_val, CNST(1.0)))

      xx    = M_ZERO
      xx(1) = CNST(0.01)
      rr    = sqrt(sum(xx**2))
      call parse_expression(pot_re, pot_im, MAX_DIM, xx, rr, M_ZERO, spec%user_def)
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
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a jellium sphere / approximated point particle.'
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

    case(SPECIES_CHARGE_DENSITY)
      spec%niwfs = int(max(2*spec%z_val, CNST(1.0)))
      spec%omega = spec%z_val
      spec%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a distribution of charge:'
        write(message(2),'(a,a)')      '   rho = ', trim(spec%rho)
        write(message(3),'(a,f11.6)')  '   Z = ', spec%z_val
        call messages_info(3)
      end if
    case default
      call messages_input_error('Species', 'Unknown species type')
    end select

    if(.not. species_is_ps(spec)) then
      ! since there is no real cap, make sure there are at least a few available
      spec%niwfs = max(5, spec%niwfs)
    endif

    SAFE_ALLOCATE(spec%iwf_l(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_m(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_i(1:spec%niwfs, 1:ispin))

    call species_iwf_fix_qn(spec, ispin, dim)

    if(species_is_ps(spec)) then
      write(message(1),'(a,i6,a,i6)') 'Number of orbitals: total = ', ps_niwfs(spec%ps), ', bound = ', spec%niwfs
    else
      write(message(1),'(a,i6,a,i6)') 'Number of orbitals: ', spec%niwfs
      nullify(spec%ps)
    endif
    if(print_info_) call messages_info(1)

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
    enddo

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

      if(in_debug_mode) then
        write(dirname, '(a)') 'debug/geometry'
        call io_mkdir(dirname)
        call species_debug(trim(dirname), this)
      end if
    end if

    POP_SUB(species_pot_init)
  end subroutine species_pot_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine species_init_from_data_object(this, index, json)
    type(species_t),     intent(out) :: this
    integer,             intent(in)  :: index
    type(json_object_t), intent(in)  :: json

    integer :: ierr

    PUSH_SUB(species_init_from_data_object)

    if(.not. initialized) call species_init_global()

    call species_nullify(this)
       
    this%index=index
    call json_get(json, "label", this%label, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "label" from species data object.'
      call messages_fatal(1)
      return
    end if
    call json_get(json, "type", this%type, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "type" from species data object.'
      call messages_fatal(1)
      return
    end if
    call json_get(json, "z_val", this%z_val, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "z_val" from species data object.'
      call messages_fatal(1)
      return
    end if
    call json_get(json, "mass", this%mass, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "mass" from species data object.'
      call messages_fatal(1)
      return
    end if
    this%has_density=.false.
    this%user_def=""
    nullify(this%ps)
    this%nlcc=.false.
    call json_get(json, "def_rsize", this%def_rsize, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "def_rsize" from species data object.'
      call messages_fatal(1)
      return
    end if
    this%def_h=-M_ONE
    this%niwfs=-1
    nullify(this%iwf_l)
    nullify(this%iwf_m)
    nullify(this%iwf_i)

    POP_SUB(species_init_from_data_object)
  end subroutine species_init_from_data_object
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine species_create_data_object(this, json)
    type(species_t),     intent(in)  :: this
    type(json_object_t), intent(out) :: json

    PUSH_SUB(species_create_data_object)

    call json_init(json)
    call json_set(json, "label", trim(adjustl(this%label)))
    call json_set(json, "type", this%type)
    call json_set(json, "z_val", this%z_val)
    call json_set(json, "mass", this%mass)
    call json_set(json, "def_rsize", this%def_rsize)

    POP_SUB(species_create_data_object)
  end subroutine species_create_data_object
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
  character(len=200) pure function species_rho_string(spec)
    type(species_t), intent(in) :: spec
    species_rho_string = trim(spec%rho)
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
  CMPLX function species_userdef_pot(spec, dim, xx, r)
    type(species_t),   intent(in) :: spec
    integer,           intent(in) :: dim
    FLOAT,             intent(in) :: xx(1:MAX_DIM), r

    FLOAT :: pot_re, pot_im

    PUSH_SUB(species_userdef_pot)
    
    call parse_expression(pot_re, pot_im, dim, xx, r, M_ZERO, spec%user_def)
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
      if ( spec%ps%l_max /= 0 ) species_is_local = .false. 
    end if

    POP_SUB(species_is_local)
  end function species_is_local
  ! ---------------------------------------------------------



  logical function species_represents_real_atom(spec)
    type(species_t), intent(in) :: spec
    
    integer :: type
    species_represents_real_atom = .true.

    PUSH_SUB(species_represents_real_atom)
    
    type = species_type(spec)
    species_represents_real_atom = (type /= SPECIES_USDEF .and. type /= SPECIES_CHARGE_DENSITY &
      .and. type /= SPECIES_FROM_FILE .and. type /= SPECIES_JELLIUM_SLAB)
    
    POP_SUB(species_represents_real_atom)
  end function species_represents_real_atom

  ! ---------------------------------------------------------
  !> This routine returns the non-local projector and its
  !! derivative, built using real spherical harmonics
  subroutine species_real_nl_projector(spec, x, l, lm, i, uV, duV)
    type(species_t),   intent(in)  :: spec
    FLOAT,             intent(in)  :: x(:)
    integer,           intent(in)  :: l, lm, i
    FLOAT,             intent(out) :: uV
    FLOAT,             intent(out) :: duV(:)

    FLOAT :: r, uVr0, duvr0, ylm, gylm(1:3)
    FLOAT, parameter :: ylmconst = CNST(0.488602511902920) !  = sqrt(3/(4*pi))

    ! no push_sub because this function is called very frequently

    r = sqrt(sum(x(1:3)**2))

    uVr0  = spline_eval(spec%ps%kb(l, i), r)
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
      enddo
    endif

    POP_SUB(species_nl_projector)
  end subroutine species_nl_projector
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Return radius outside which orbital is less than threshold value 0.001
  FLOAT function species_get_iwf_radius(spec, ii, is) result(radius)
    type(species_t),   intent(in) :: spec
    integer,           intent(in) :: ii !< principal quantum number
    integer,           intent(in) :: is !< spin component

    FLOAT, parameter :: threshold = CNST(0.001)

    PUSH_SUB(species_get_iwf_radius)

    if(species_is_ps(spec)) then
      ASSERT(ii <= spec%ps%conf%p)
      radius = spline_cutoff_radius(spec%ps%ur(ii, is), threshold)
    else if(species_represents_real_atom(spec)) then
      radius = -ii*log(threshold)/spec%Z_val
    else
      radius = sqrt(-M_TWO*log(threshold)/spec%omega)
    end if

    ! The values for hydrogenic and harmonic-oscillator wavefunctions
    ! come from taking the exponential part (i.e. the one that controls
    ! the asymptotic behavior at large r), and setting it equal to
    ! the threshold.

    POP_SUB(species_get_iwf_radius)
  end function species_get_iwf_radius
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
    this%has_density=that%has_density
    this%user_def=that%user_def
    this%omega=that%omega
    this%filename=that%filename
    this%jradius=that%jradius
    this%jthick=that%jthick
    !> To be implemented.
    !> ps_t has no copy procedure.
    nullify(this%ps)
    if(associated(that%ps))this%ps=>that%ps
    this%nlcc=that%nlcc
    this%sigma=that%sigma
    this%rho=that%rho
    this%def_rsize=that%def_rsize
    this%def_h=that%def_h
    this%niwfs=that%niwfs
    nullify(this%iwf_l, this%iwf_m, this%iwf_i)
    call loct_pointer_copy(this%iwf_l, that%iwf_l)
    call loct_pointer_copy(this%iwf_m, that%iwf_m)
    call loct_pointer_copy(this%iwf_i, that%iwf_i)
    this%lmax=that%lmax
    this%lloc=that%lloc

    POP_SUB(species_copy)
  end subroutine species_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine species_end_species(spec)
    type(species_t), intent(inout) :: spec
    
    PUSH_SUB(species_end_species)

    if (species_is_ps(spec)) then 
      if(associated(spec%ps)) then 
        call ps_end(spec%ps)
        SAFE_DEALLOCATE_P(spec%ps)
      end if
    end if
    SAFE_DEALLOCATE_P(spec%iwf_l)
    SAFE_DEALLOCATE_P(spec%iwf_m)
    SAFE_DEALLOCATE_P(spec%iwf_i)

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
    bool = species_is_local(spec)
    write(iunit, '(a,l1)')    'local  = ', bool
    write(iunit, '(2a)')      'usdef  = ', trim(spec%user_def)
    if (spec%type == SPECIES_JELLIUM) then
      write(iunit, '(a,f15.2)') 'jradius= ', spec%jradius
    end if
    if (spec%type == SPECIES_JELLIUM_SLAB) then
      write(iunit, '(a,f15.2)') 'jthick= ', spec%jthick
    end if
    write(iunit, '(a,l1)')    'nlcc   = ', spec%nlcc
    write(iunit, '(a,f15.2)') 'def_rsize = ', spec%def_rsize
    write(iunit, '(a,f15.2)') 'def_h = ', spec%def_h
    if (spec%type /= SPECIES_USDEF ) write(iunit, '(a,i3)')    'lmax  = ', spec%lmax
    if (spec%type /= SPECIES_USDEF ) write(iunit, '(a,i3)')    'lloc  = ', spec%lloc

    if(species_is_ps(spec)) then
       if(in_debug_mode) call ps_debug(spec%ps, trim(dirname))
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

    PUSH_SUB(read_from_default_file)

    backspace(iunit)

    spec%type = SPECIES_PSEUDO
    
    read(iunit,*) label, spec%filename, spec%z, spec%lmax, spec%lloc, spec%def_h, spec%def_rsize

    spec%filename = trim(conf%share)//'/pseudopotentials/'//trim(spec%filename)
    
    ASSERT(trim(label) == trim(spec%label))

    read_data = 8

    ! get the mass for this element
    call element_init(element, label)
    ASSERT(element_valid(element))
    spec%mass = element_mass(element)
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

    integer :: ncols, icol, flag
    type(element_t) :: element
    
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
      call parse_block_string(blk, row, read_data, spec%user_def)
      call conv_to_C_string(spec%user_def)
      read_data = read_data + 1

    case(SPECIES_FROM_FILE)

    case(SPECIES_JELLIUM)
      spec%jradius = CNST(0.5)
        
    case(SPECIES_JELLIUM_SLAB)
      call parse_block_float(blk, row, read_data, spec%jthick) ! thickness of the jellium slab
      read_data = read_data + 1
      if(spec%jthick <= M_ZERO) call messages_input_error('Species')
      spec%jthick = units_to_atomic(units_inp%length, spec%jthick) ! units conversion

    case(SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN)
      spec%sigma = CNST(0.25)

    case(SPECIES_CHARGE_DENSITY)
      call parse_block_string(blk, row, read_data + 1, spec%rho)
      read_data = read_data + 1

    case(SPECIES_PSEUDO)

    case(SPECIES_PSPIO) ! a pseudopotential file to be handled by the pspio library

    case default
      call messages_input_error('Species', "Unknown type for species '"//trim(spec%label)//"'")
    end select

    spec%mass = -CNST(1.0)
    spec%z_val = -CNST(1.0)
    spec%sc_alpha = -CNST(1.0)
    
    icol = read_data
    do
      if(icol >= ncols) exit

      call parse_block_integer(blk, row, icol, flag)
      
      select case(flag)

      case(OPTION_MIN_RADIUS)
        call parse_block_float(blk, row, icol + 1, spec%def_h)

      case(OPTION_MAX_SPACING)
        call parse_block_float(blk, row, icol + 1, spec%def_rsize)

      case(OPTION_LMAX)
        call parse_block_integer(blk, row, icol + 1, spec%lmax)

      case(OPTION_LLOC)
        call parse_block_integer(blk, row, icol + 1, spec%lloc)

      case(OPTION_MASS)
        call parse_block_float(blk, row, icol + 1, spec%mass)

      case(OPTION_VALENCE)
        call parse_block_float(blk, row, icol + 1, spec%z_val)
        spec%z = spec%z_val

      case(OPTION_JELLIUM_RADIUS)
        call parse_block_float(blk, row, icol + 1, spec%jradius)
        spec%jradius = units_to_atomic(units_inp%length, spec%jradius)
        if(spec%jradius <= M_ZERO) call messages_input_error('Species', 'jellium_radius must be positive')
        if(spec%type /= SPECIES_JELLIUM) then
          call messages_input_error('Species', 'jellium_radius can only be used with species_jellium')
        end if
        
      case(OPTION_GAUSSIAN_WIDTH)
        call parse_block_float(blk, row, icol + 1, spec%sigma)
        if(spec%sigma <= M_ZERO) call messages_input_error('Species', 'gaussian_width must be positive')
        if(spec%type /= SPECIES_FULL_GAUSSIAN) then
          call messages_input_error('Species', 'gaussian_width can only be used with species_full_gaussian')
        end if

      case(OPTION_SOFTENING)
        call parse_block_float(blk, row, icol + 1, spec%sc_alpha)
        spec%sc_alpha = units_to_atomic(units_inp%length, spec%sc_alpha)**2
        if(spec%type /= SPECIES_SOFT_COULOMB) then
          call messages_input_error('Species', 'softening can only be used with species_soft_coulomb')
        end if

      case(OPTION_FILE)
        call parse_block_string(blk, row, icol + 1, spec%filename)

      case(OPTION_DB_FILE)
        call parse_block_string(blk, row, icol + 1, spec%filename)
        spec%filename = trim(conf%share)//'/pseudopotentials/'//trim(spec%filename)
        
      case default
        call messages_input_error('Species', "Unknown parameter in species '"//trim(spec%label)//"'")
        
      end select

      icol = icol + 2        
    end do

    if(spec%type == SPECIES_SOFT_COULOMB .and. spec%sc_alpha <= CNST(0.0)) then
      call messages_write("The mandatory 'softening' parameter is missing for species "//trim(spec%label)//'.')
      call messages_fatal()
    end if

    if(spec%filename == '' .and. &
      (spec%type == SPECIES_PSEUDO .or. spec%type == SPECIES_PSPIO .or. spec%type == SPECIES_FROM_FILE)) then
      call messages_input_error('Species', "The file or db_file parameter is missing for species '"//trim(spec%label)//"'")
    end if
      
    
    select case(spec%type)
    case(SPECIES_PSEUDO, SPECIES_PSPIO, SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN)

      call element_init(element, spec%label)
      
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
        
      call element_end(element)
        
    case default
      if(spec%mass < CNST(0.0)) then
        spec%mass = 1.0
        call messages_write('Info: default mass for species '//trim(spec%label)//':')
        call messages_write(spec%mass)
        call messages_write(' amu.')
        call messages_info()
      end if

      if(spec%z_val < CNST(0.0)) then
        call messages_write('Cannot determine valence for species '//trim(spec%label)//'.')
        call messages_fatal()
      end if
      
    end select
    
    POP_SUB(read_from_block)
  end subroutine read_from_block
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> set up quantum numbers of orbitals, and reject those that are unbound (for pseudopotentials)
  subroutine species_iwf_fix_qn(spec, ispin, dim)
    type(species_t), intent(inout) :: spec
    integer,         intent(in)    :: ispin
    integer,         intent(in)    :: dim

    integer :: is, n, i, l, m, n1, n2, n3
    FLOAT   :: radius
    logical, allocatable :: bound(:)

    PUSH_SUB(species_iwf_fix_qn)

    if(species_is_ps(spec)) then
      
      SAFE_ALLOCATE(bound(1:spec%ps%conf%p))

      ! we check if the orbitals are bound by looking at the atomic radius
      do i = 1, spec%ps%conf%p
        radius = M_ZERO
        do is = 1, ispin
          radius = max(radius, spline_cutoff_radius(spec%ps%ur(i, is), threshold = CNST(0.001)))
        end do
        ! we consider as bound a state that is localized to less than half the radius of the radial grid
        bound(i) = radius < CNST(0.5)*logrid_radius(spec%ps%g)
      end do
      
      do is = 1, ispin
        n = 1
        do i = 1, spec%ps%conf%p
          if(n > spec%niwfs) exit          
          l = spec%ps%conf%l(i)
           
          if(.not. bound(i)) cycle
          
          do m = -l, l
            spec%iwf_i(n, is) = i
            spec%iwf_l(n, is) = l
            spec%iwf_m(n, is) = m
            n = n + 1
          end do
          
        end do
        ! FIXME: this is wrong when spin-polarized or spinors!
        spec%niwfs = n - 1
      end do

      SAFE_DEALLOCATE_A(bound)

    else if(species_represents_real_atom(spec) .and. dim == 3) then

      do is = 1, ispin
        n = 1
        ! just up to the highest principal quantum number, actually
        do i = 1, spec%niwfs
          if(n > spec%niwfs) exit
          do l = 0, i-1
            do m = -l, l
              spec%iwf_i(n, is) = i
              spec%iwf_l(n, is) = l
              spec%iwf_m(n, is) = m
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
            spec%iwf_l(i, is) = 0
            spec%iwf_m(i, is) = 0
          end do
        end do

      case(2)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1
          do
            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = 0
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = 0
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = 0
            i = i + 1; if(i>spec%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1
          end do
        end do

      case(3)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1; n3 = 1
          do
            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = 0
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3+1
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = n3
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1+1
            spec%iwf_l(i, is) = n2
            spec%iwf_m(i, is) = n3+1
            i = i + 1; if(i>spec%niwfs) exit

            spec%iwf_i(i, is) = n1
            spec%iwf_l(i, is) = n2+1
            spec%iwf_m(i, is) = n3+1
            i = i + 1; if(i>spec%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1; n3 = n3 + 1
          end do
        end do
      end select
    end if

    POP_SUB(species_iwf_fix_qn)
  end subroutine species_iwf_fix_qn
  ! ---------------------------------------------------------

end module species_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
