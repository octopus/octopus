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

module species_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use loct_math_m
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
  public ::                    &
    species_t,                 &
    species_set_label,         &
    species_set_index,         &
    species_read,              &
    species_init,              &
    species_pot_init,          &
    species_type,              &
    species_label,             &
    species_index,             &
    species_has_nlcc,          &
    species_has_density,       &
    species_ps,                &
    species_zval,              &
    species_z,                 &
    species_def_rsize,         &
    species_def_h,             &
    species_jradius,           &
    species_jthick,            &
    species_sigma,             &
    species_omega,             &
    species_weight,            &
    species_rho_string,        &
    species_filename,          &
    species_niwfs,             &
    species_iwf_ilm,           &
    species_userdef_pot,       &
    species_is_ps,             &
    species_is_full,           &
    species_is_local,          &
    species_real_nl_projector, &
    species_nl_projector,      &
    species_get_iwf_radius,    &
    species_end


  integer, public, parameter :: &
    SPEC_USDEF  = 123,          & !< user-defined function
    SPEC_POINT  = 2,            & !< point charge: jellium sphere of radius 0.5 a.u.
    SPEC_JELLI  = 3,            & !< jellium sphere.
    SPEC_JELLI_SLAB     = 4,    & !< jellium slab.
    SPEC_FULL_DELTA     = 127,  & !< full-potential atom
    SPEC_FULL_GAUSSIAN  = 124,  & !< full-potential atom
    SPEC_CHARGE_DENSITY = 125,  &
    SPEC_FROM_FILE = 126,       &
    SPEC_PS_PSF = PS_TYPE_PSF,  & !< SIESTA pseudopotential
    SPEC_PS_HGH = PS_TYPE_HGH,  & !< HGH pseudopotential
    SPEC_PS_CPI = PS_TYPE_CPI,  & !< FHI pseudopotential (cpi format)
    SPEC_PS_FHI = PS_TYPE_FHI,  & !< FHI pseudopotential (ABINIT format)
    SPEC_PS_UPF = PS_TYPE_UPF     !< UPF pseudopotential

  type species_t
    private
    integer :: index              !< just a counter

    character(len=15) :: label    !< Identifier for the species
    integer :: type               !< what type of species
    FLOAT   :: z                  !< charge of the species
    FLOAT   :: z_val              !< valence charge of the species -- the total charge
                                  !< minus the core charge in the case of the pseudopotentials
    FLOAT   :: weight             !< mass, in atomic mass units (!= atomic units of mass)

    logical :: has_density        !< true if the species has an electronic density


    character(len=1024) :: user_def !< for the user-defined potential
    FLOAT :: omega


    character(len=200) :: filename !< for the potential read from a file.


    FLOAT :: jradius              !< jellium stuff
    FLOAT :: jthick               !< jellium stuff


    type(ps_t), pointer :: ps
    logical             :: nlcc   !< true if we have non-local core corrections


    FLOAT :: sigma                !< If we have an all-electron atom:


    character(len=200) :: rho     !< If we have a charge distribution creating the potential:


    FLOAT :: def_rsize, def_h     !< the default values for the spacing and atomic radius


    integer :: niwfs              !< The number of initial wavefunctions
    integer, pointer :: iwf_l(:, :), iwf_m(:, :), iwf_i(:, :)


    integer :: lmax, lloc         !< For the TM pseudos, the lmax and lloc.
  end type species_t


contains


  ! ---------------------------------------------------------
  !> Assigns a label to a species_t variable. This should be the
  !! first routine to be called (before species_index, species_read and species_init).
  !! This label must match one of the labels given in the %Species block
  !! in the input file -- or else one of the labels in the defaults file.
  ! ---------------------------------------------------------
  pure subroutine species_set_label(spec, label)
    type(species_t), intent(inout) :: spec
    character(len=*), intent(in)   :: label
    spec%label = trim(label)
  end subroutine species_set_label
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Assigns an index to a species_t variable. This should be the
  !! second routine to be called (before species_read and species_init),
  !! when initializing a species_t variable.
  ! ---------------------------------------------------------
  pure subroutine species_set_index(spec, k)
    type(species_t), intent(inout) :: spec
    integer, intent(in)   :: k
    spec%index = k
  end subroutine species_set_index
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> Reads the information (from the input file) about a species_t variable, initializing
  !! part of it (it has to be completed later with "species_init").
  !! Note that species_read has to be called only after species_set_label
  !! and species_set index have been called.
  ! ---------------------------------------------------------
  subroutine species_read(spec)
    type(species_t), intent(inout) :: spec

    character(len=256) :: fname
    character(len=10)  :: lab
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
    !% through a pseudopotential, or a user-defined model potential.
    !%
    !% Note that some common pseudopotentials are distributed with the code in the
    !% directory <tt>OCTOPUS-HOME/share/PP/</tt>. To use these pseudopotentials you are
    !% not required to define them explicitly in the <tt>Species</tt> block, as defaults 
    !% are provided by the program (you can override these defaults in any case). 
    !% Additional pseudopotentials can be downloaded from the 
    !% <a href='http://www.tddft.org/programs/octopus/wiki/index.php/Pseudopotentials'>
    !% octopus homepage</a>.
    !%
    !% The format of this block is the following: The first field is
    !% the name of the species, followed by the atomic mass (in atomic mass
    !% units). The third field defines the type of species (the valid options
    !% are detailed below). Each type also needs some parameters given in
    !% the remaining fields of the row.
    !%
    !% In 3D, <i>e.g.</i>
    !%
    !% <tt>%Species
    !% <br>&nbsp;&nbsp;'O'       | 15.9994 | spec_ps_psf         | 8   | 1 | 1
    !% <br>&nbsp;&nbsp;'H'       |  1.0079 | spec_ps_hgh         | 1   | 0 | 0
    !% <br>&nbsp;&nbsp;'jlm'     | 23.2    | spec_jelli          | 8   | 5.0
    !% <br>&nbsp;&nbsp;'rho'     | 17.0    | spec_charge_density | 6   | "exp(-r/a)"
    !% <br>&nbsp;&nbsp;'pnt'     | 32.3    | spec_point          | 2.0
    !% <br>&nbsp;&nbsp;'udf'     |  0.0    | spec_user_defined   | 8   | "1/2*r^2"
    !% <br>&nbsp;&nbsp;'H_all'   |  1.0079 | spec_full_delta     | 1
    !% <br>&nbsp;&nbsp;'H_all'   |  1.0079 | spec_full_gaussian  | 1
    !% <br>%</tt>
    !%
    !% Additionally, all the pseudopotential types (PSF, HGH, FHI, UPF) can take two extra
    !% fields: default spacing, and default radius (used for minimum simulation box if the
    !% radius is not specified).
    !%Option spec_user_defined 123
    !% Species with user-defined potential. In this case, the fourth
    !% field is the valence charge and the fifth
    !% field is a string with a mathematical expression that defines the
    !% potential (you can use any of the <i>x</i>, <i>y</i>, <i>z</i>
    !% or <i>r</i> variables).
    !%Option spec_point  2
    !% Point charge: the fourth field is the value of the charge.
    !%Option spec_jelli  3
    !% Jellium sphere: the extra parameters are the charge of the jellium
    !% sphere (an equal value of valence charge is assumed) and the radius of
    !% the sphere.
    !%Option spec_jelli_slab  4
    !% Jellium slab: the extra parameters are the charge of the jellium
    !% slab (an equal value of valence charge is assumed) and the thickness of
    !% the slab. The slab extends across the simulation box in the <i>xy</i>-plane.
    !%Option spec_ps_psf  100
    !% Troullier Martins pseudopotential in <tt>SIESTA</tt> format: the pseudopotential will be
    !% read from a <tt>.psf</tt> file, either in the working
    !% directory or in the <tt>OCTOPUS-HOME/share/octopus/PP/PSF</tt> directory.
    !% Columns 4, 5, 6 are the atomic number, the maximum
    !% <i>l</i>-component of the pseudopotential to consider in the
    !% calculation, and the <i>l</i>-component to consider as local.
    !%Option spec_ps_hgh  101
    !% Hartwigsen-Goedecker-Hutter pseudopotentials: column 4 is
    !% the atomic number and columns 5 and 6 are irrelevant, since they
    !% are not necessary to define the HGH pseudopotential.
    !%Option spec_ps_cpi  102
    !% Fritz-Haber pseudopotential: the pseudopotential will be
    !% read from a <tt>.cpi</tt> file, either in the working
    !% directory or in the <tt>OCTOPUS-HOME/share/PP/CPI</tt> directory.
    !% Columns 4, 5, 6 are the atomic number, the maximum
    !% <i>l</i>-component of the pseudopotential to consider in the
    !% calculation, and the <i>l</i>-component to consider as local.
    !%Option spec_ps_fhi  103
    !% Fritz-Haber pseudopotential (<tt>ABINIT</tt> format): the pseudopotential will be
    !% read from a <tt>.fhi</tt> file, either in the working
    !% directory or in the <tt>OCTOPUS-HOME/share/PP/FHI</tt> directory.
    !% Columns 4, 5, 6 are the atomic number, the maximum
    !% <i>l</i>-component of the pseudopotential to consider in the
    !% calculation, and the <i>l</i>-component to consider as local.
    !% Note that you can use the pseudopotentials from <tt>ABINIT</tt> homepage.
    !%Option spec_ps_upf  104
    !% UPF format: the pseudopotential will be
    !% read from a <tt>.UPF</tt> file, either in the working
    !% directory or in the <tt>OCTOPUS-HOME/share/PP/UPF</tt> directory.
    !% Column 4 is the atomic number. Columns 5 and 6 are 
    !% ignored, as the maximum <i>l</i>-component of the pseudopotential to
    !% consider in the calculation and the <i>l</i>-component to consider as
    !% local are indicated in the pseudopotential file are cannot be changed.
    !%Option spec_full_delta   127
    !% Full atomic potential represented by a delta charge
    !% distribution. The atom will be displaced to the nearest grid
    !% point. Column 4 is the atomic number.
    !%Option spec_full_gaussian   124
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
    !% Column 4 is the atomic number and column 5 is
    !% <math>sigma</math>, the width of the gaussian that should be
    !% small, but you may run into numerical difficulties if it is too
    !% small (0.25 by default).
    !%Option spec_charge_density 125
    !% The potential is created by a distribution of charge. The extra parameters are the
    !% valence charge of the species, and an expression for the charge distribution.
    !%Option species_from_file  126
    !% The potential is read from a file, whose name is given in column 5.
    !%End

    call messages_obsolete_variable('SpecieAllElectronSigma', 'Species')
    call messages_obsolete_variable('SpeciesAllElectronSigma', 'Species')

    ! First, find out if there is a Species block.
    n_spec_block = 0
    if(parse_block(datasets_check('Species'), blk) == 0) then
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
    end if

    ! Find out if the species is in the defaults file.
    write(fname, '(2a)') trim(conf%share), "/PP/defaults"
    n_spec_def = max(0, loct_number_of_lines(fname))
    if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

    iunit = io_open(fname, action='read', status='old', die=.false., is_tmp=.true.)
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
    end if

    if(read_data == 0) then
      message(1) = 'Species '//trim(spec%label)//' not found.'
      call messages_fatal(1)
    end if

    POP_SUB(species_read)
  end subroutine species_read
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_init(spec, ispin, space, print_info)
    type(species_t),   intent(inout) :: spec
    integer,           intent(in)    :: ispin
    type(space_t),     intent(in)    :: space
    logical, optional, intent(in)    :: print_info

    logical :: print_info_
    integer :: i
    FLOAT   :: pot_re, pot_im, xx(MAX_DIM), rr

    PUSH_SUB(species_init)

    print_info_ = .true.
    if(present(print_info)) then
      print_info_ = print_info
    end if

    ! masses are always in amu, so convert them to a.u.
    spec%weight =  units_to_atomic(units_inp%mass, spec%weight)

    spec%has_density = .false.

    select case(spec%type)
    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
      ! allocate structure
      SAFE_ALLOCATE(spec%ps) 
      call ps_init(spec%ps, spec%label, spec%type, spec%Z, spec%lmax, spec%lloc, ispin)
      spec%z_val = spec%ps%z_val
      spec%nlcc = (spec%ps%icore /= 'nc  ' )
      spec%niwfs = ps_niwfs(spec%ps)

    case(SPEC_USDEF)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a user-defined potential.'
        i = min(237, len_trim(spec%user_def)-1) ! I subtract 1 to avoid the non-printable C "end-of-string" character.
        write(message(2),'(a,a)')      '   Potential = ', trim(spec%user_def(1:i))
        if(len(trim(spec%user_def)).gt.237) then
          message(2) = trim(message(2))//'...'
        end if
        call messages_info(2)
      end if
      spec%niwfs = int(max(2*spec%z_val, CNST(1.0)))

      xx    = M_ZERO
      xx(1) = CNST(0.01)
      rr    = sqrt(sum(xx**2))
      call parse_expression(pot_re, pot_im, MAX_DIM, xx, rr, M_ZERO, spec%user_def)
      spec%omega = sqrt( abs(M_TWO / CNST(1.0e-4) * pot_re ))
      ! To avoid problems with constant potentials.
      if(spec%omega <= M_ZERO) spec%omega = CNST(0.1) 

    case(SPEC_FROM_FILE)
      if(print_info_) then
        write(message(1),'(a)') 'Species read from file "'//trim(spec%filename)//'".'
        call messages_info(1)
      end if
      spec%niwfs = 2*spec%z_val
      spec%omega = CNST(0.1)

    case(SPEC_JELLI, SPEC_POINT)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a jellium sphere / approximated point particle.'
        write(message(2),'(a,f11.6)')  '   Valence charge = ', spec%z_val
        write(message(3),'(a,f11.6)')  '   Radius [a.u]   = ', spec%jradius
        write(message(4),'(a,f11.6)')  '   Rs [a.u]       = ', spec%jradius * spec%z_val ** (-M_ONE/M_THREE)
        call messages_info(4)
      end if
      spec%niwfs = 2*spec%z_val
      spec%omega = CNST(0.1)

    case(SPEC_JELLI_SLAB)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is a jellium slab.'
        write(message(2),'(a,f11.6)')  '   Valence charge  = ', spec%z_val
        write(message(3),'(a,f11.6)')  '   Thickness [a.u] = ', spec%jthick
        !write(message(4),'(a,f11.6)')  '   Rs [a.u]       = ', ( M_THREE /( M_FOUR *M_PI ) &
        !& *spec%z_val /( *sb%lsize(1) *sb%lsize(2) ) )**(1.0/3.0) 
        call messages_info(3)
      end if
      spec%niwfs = 2*spec%z_val
      spec%omega = CNST(0.1)

    case(SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN)
      spec%niwfs = 2*spec%z_val
      spec%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(spec%label),'" is an all-electron atom.'
        write(message(2),'(a,f11.6)')  '   Z = ', spec%z_val
        write(message(3),'(a)')  '   Potential will be calculated solving the Poisson equation'
        write(message(4),'(a)')  '   for a delta density distribution.'
        call messages_info(4)
      end if
      spec%omega = spec%z_val 

    case(SPEC_CHARGE_DENSITY)
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
      call input_error('Species')
    end select

    SAFE_ALLOCATE(spec%iwf_l(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_m(1:spec%niwfs, 1:ispin))
    SAFE_ALLOCATE(spec%iwf_i(1:spec%niwfs, 1:ispin))
    call species_iwf_fix_qn(spec, ispin, space)

    POP_SUB(species_init)
  end subroutine species_init
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

      if(filter .ne. PS_FILTER_NONE) then 
        call ps_filter(this%ps, filter, grid_cutoff)
        call ps_getradius(this%ps) ! radius may have changed
      end if

      call ps_derivatives(this%ps)

      local_radius = spline_cutoff_radius(this%ps%vl, this%ps%projectors_sphere_threshold)

      orbital_radius = M_ZERO
      do iorb = 1, species_niwfs(this)
        orbital_radius = max(orbital_radius, species_get_iwf_radius(this, iorb, is = 1))
      end do

      write(message(1), '(3a)') "Info: Pseudopotential for ", trim(this%label), ". Radii for localized parts:"
      write(message(2), '(a,f5.1, 3a)') "     local part     = ", &
        units_from_atomic(units_out%length, local_radius), " [", trim(units_abbrev(units_out%length)), "] "
      write(message(3), '(a,f5.1,3a)')  "     non-local part = ", &
        units_from_atomic(units_out%length, this%ps%rc_max), " [", trim(units_abbrev(units_out%length)), "] "
      write(message(4), '(a,f5.1,3a)')  "     orbitals       = ", &
        units_from_atomic(units_out%length, orbital_radius), " [", trim(units_abbrev(units_out%length)), "] "
      call messages_info(4)

      if(max(local_radius, this%ps%rc_max) > CNST(6.0)) then
        message(1) = "One of the radii of your pseudopotential's localized parts seems"
        message(2) = "unusually large; check that your pseudopotential is correct."
        call messages_warning(2)
      end if

      if(orbital_radius > CNST(20.0)) then
        message(1) = "The radius of the atomic orbitals given by your pseudopotential seems"
        message(2) = "unusually large; check that your pseudopotential is correct."
        call messages_warning(2)
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
  integer pure function species_type(spec)
    type(species_t), intent(in) :: spec
    species_type = spec%type
  end function species_type
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=15) pure function species_label(spec)
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
  FLOAT pure function species_weight(spec)
    type(species_t), intent(in) :: spec
    species_weight = spec%weight
  end function species_weight
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
  FLOAT function species_userdef_pot(spec, dim, xx, r)
    type(species_t), intent(in) :: spec
    integer, intent(in) :: dim
    FLOAT, intent(in) :: xx(1:MAX_DIM), r
    FLOAT :: pot_re, pot_im
    call parse_expression(pot_re, pot_im, dim, xx, r, M_ZERO, spec%user_def)
    species_userdef_pot = pot_re
  end function species_userdef_pot
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_is_ps(spec)
    type(species_t), intent(in) :: spec
    
    species_is_ps = &
         ( spec%type == SPEC_PS_PSF) .or. &
         ( spec%type == SPEC_PS_HGH) .or. &
         ( spec%type == SPEC_PS_CPI) .or. &
         ( spec%type == SPEC_PS_FHI) .or. &
         ( spec%type == SPEC_PS_UPF)
 
  end function species_is_ps

  ! ---------------------------------------------------------

  logical elemental function species_is_full(spec)
    type(species_t), intent(in) :: spec
    
    species_is_full = &
         ( spec%type == SPEC_FULL_GAUSSIAN) .or. &
         ( spec%type == SPEC_FULL_DELTA)
    
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

    call grylmr(x(1), x(2), x(3), l, lm, ylm, gylm)
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
  subroutine species_nl_projector(spec, x, l, lm, i, uV)
    type(species_t),   intent(in)  :: spec
    FLOAT,             intent(in)  :: x(1:MAX_DIM)
    integer,           intent(in)  :: l, lm, i
    CMPLX,             intent(out) :: uV

    FLOAT :: r, uVr0
    CMPLX :: ylm

    ! no push_sub because this function is called very frequently
    r = sqrt(sum(x(1:MAX_DIM)**2))

    uVr0 = spline_eval(spec%ps%kb(l, i), r)

    call ylmr(x(1), x(2), x(3), l, lm, ylm)
    uv = uvr0*ylm

  end subroutine species_nl_projector
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  FLOAT function species_get_iwf_radius(spec, j, is) result(radius)
    type(species_t),   intent(in) :: spec
    integer,           intent(in) :: j
    integer,           intent(in) :: is

    integer :: i, l, m
    FLOAT, parameter :: threshold = CNST(0.001)

    PUSH_SUB(species_get_iwf_radius)

    i = spec%iwf_i(j, is)
    l = spec%iwf_l(j, is)
    m = spec%iwf_m(j, is)

    if(species_is_ps(spec)) then
      radius = spline_cutoff_radius(spec%ps%ur(i, is), threshold)
    else
      radius = sqrt(-M_TWO*log(threshold)/spec%omega)
    end if

    POP_SUB(species_get_iwf_radius)

  end function species_get_iwf_radius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_end(ns, spec)
    integer,         intent(in) :: ns
    type(species_t), pointer    :: spec(:)

    integer :: i

    PUSH_SUB(species_end)

    do i = 1, ns
      if (species_is_ps(spec(i))) then 
        if(associated(spec(i)%ps)) then 
          call ps_end(spec(i)%ps)
          SAFE_DEALLOCATE_P(spec(i)%ps)
        end if
      end if
      SAFE_DEALLOCATE_P(spec(i)%iwf_l)
      SAFE_DEALLOCATE_P(spec(i)%iwf_m)
      SAFE_DEALLOCATE_P(spec(i)%iwf_i)
    end do

    POP_SUB(species_end)
  end subroutine species_end
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
    if (spec%type /= SPEC_USDEF ) write(iunit, '(a,f15.2)') 'z      = ', spec%z
    if (spec%type == SPEC_FROM_FILE) then
      write(iunit,'(a)')      'Species read from file "'//trim(spec%filename)//'".'
    end if
    write(iunit, '(a,f15.2)') 'z_val  = ', spec%z_val
    write(iunit, '(a,f15.2)') 'weight = ', spec%weight
    bool = species_is_local(spec)
    write(iunit, '(a,l1)')    'local  = ', bool
    write(iunit, '(2a)')      'usdef  = ', trim(spec%user_def)
    if (spec%type == SPEC_JELLI .or. spec%type == SPEC_POINT) then
      write(iunit, '(a,f15.2)') 'jradius= ', spec%jradius
    end if
    if (spec%type == SPEC_JELLI_SLAB) then
      write(iunit, '(a,f15.2)') 'jthick= ', spec%jthick
    end if
    write(iunit, '(a,l1)')    'nlcc   = ', spec%nlcc
    write(iunit, '(a,f15.2)') 'def_rsize = ', spec%def_rsize
    write(iunit, '(a,f15.2)') 'def_h = ', spec%def_h
    if (spec%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lmax  = ', spec%lmax
    if (spec%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lloc  = ', spec%lloc

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

    character(len=10) :: label
    FLOAT :: weight, z, def_h, def_rsize
    integer :: lloc, lmax
    integer :: type

    PUSH_SUB(read_from_default_file)

    backspace(iunit)

    read(iunit,*) label, weight, type, z, lmax, lloc, def_h, def_rsize

    ASSERT(trim(label) == trim(spec%label))

    if(read_data == 0) then ! The Species was not supplied in the block.
      spec%weight    = weight
      spec%type      = type
    end if
    if(read_data < 4) spec%z         = z
    if(read_data < 5) spec%lmax      = lmax
    if(read_data < 6) spec%lloc      = lloc
    if(read_data < 7) spec%def_h     = def_h
    if(read_data < 8) spec%def_rsize = def_rsize

    read_data = 8

    POP_SUB(read_from_default_file)
  end subroutine read_from_default_file
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine read_from_block(blk, row, spec, read_data)
    type(block_t),   intent(in)    :: blk
    integer,         intent(in)    :: row
    type(species_t), intent(inout) :: spec
    integer,         intent(out)   :: read_data

    integer :: nn

    PUSH_SUB(read_from_block)

    read_data = 0

    call parse_block_float (blk, row, 1, spec%weight)
    call parse_block_integer (blk, row, 2, spec%type)

    select case(spec%type)
    case(SPEC_USDEF) ! user-defined
      call parse_block_float (blk, row, 3, spec%Z_val)
      call parse_block_string(blk, row, 4, spec%user_def)
      call conv_to_C_string(spec%user_def)
      read_data = 5

    case(SPEC_FROM_FILE)
      call parse_block_float (blk, row, 3, spec%Z_val)
      call parse_block_string(blk, row, 4, spec%filename)
      read_data = 5

    case(SPEC_POINT) ! this is treated as a jellium with radius 0.5
      call parse_block_float(blk, row, 3, spec%Z)
      spec%jradius = M_HALF
      spec%Z_val = 0
      read_data = 4

    case(SPEC_JELLI)
      call parse_block_float(blk, row, 3, spec%Z)      ! charge of the jellium sphere
      call parse_block_float(blk, row, 4, spec%jradius)! radius of the jellium sphere
      spec%jradius = units_to_atomic(units_inp%length, spec%jradius) ! units conversion
      spec%Z_val = spec%Z
      read_data = 5

    case(SPEC_JELLI_SLAB)
      call parse_block_float(blk, row, 3, spec%Z)      ! charge of the jellium slab
      call parse_block_float(blk, row, 4, spec%jthick) ! thickness of the jellium slab
      spec%jthick = units_to_atomic(units_inp%length, spec%jthick) ! units conversion
      spec%Z_val = spec%Z
      read_data = 5

    case(SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN)
      call parse_block_float(blk, row, 3, spec%Z)
      spec%Z_val = spec%Z
      read_data = 4
      
      if (parse_block_cols(blk, row) <= 4) then
        spec%sigma = CNST(0.25)
      else
        call parse_block_float(blk, row, 4, spec%sigma)
        if(spec%sigma <= M_ZERO) call input_error('Species')
      end if

    case(SPEC_CHARGE_DENSITY)
      call parse_block_float(blk, row, 3, spec%Z)
      call parse_block_string(blk, row, 4, spec%rho)
      spec%Z_val = spec%Z
      read_data = 5

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF) ! a pseudopotential file
      nn = parse_block_cols(blk, row)

      call parse_block_float (blk, row, 3, spec%Z)

      if(spec%type == SPEC_PS_UPF) then 
        read_data = 4
      end if

      if(nn > 4) then
        call parse_block_integer (blk, row, 4, spec%lmax)
        read_data = 5
      end if

      if(nn > 5) then
        call parse_block_integer (blk, row, 5, spec%lloc)
        read_data = 6
      end if

      if(nn > 6) then
        call parse_block_float (blk, row, 6, spec%def_h)
        spec%def_h = units_to_atomic(units_inp%length, spec%def_h)
        read_data = 7
      end if

      if(nn > 7) then
        call parse_block_float (blk, row, 7, spec%def_rsize)
        spec%def_rsize = units_to_atomic(units_inp%length, spec%def_rsize)
        read_data = 8
      end if

    case default
      call input_error('Species')
    end select

    POP_SUB(read_from_block)
  end subroutine read_from_block
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_iwf_fix_qn(spec, ispin, space)
    type(species_t), intent(inout) :: spec
    integer,         intent(in)    :: ispin
    type(space_t),   intent(in)    :: space

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
        ! we consider as bound a state that localized to less than half the radius of the radial grid
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
        spec%niwfs = n - 1
      end do

      SAFE_DEALLOCATE_A(bound)

    else

      select case(space%dim)
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
