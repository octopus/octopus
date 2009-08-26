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
  use splines_m
  use loct_m
  use loct_math_m
  use loct_parser_m
  use math_m
  use messages_m
  use mpi_m
  use profiling_m
  use ps_m
  use string_m
  use units_m

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
    species_sigma,             &
    species_omega,             &
    species_weight,            &
    species_rho_string,        &
    species_filename,          &
    species_niwfs,             &
    species_iwf_ilm,           &
    species_userdef_pot,       &
    species_is_ps,             &
    species_is_local,          &
    species_real_nl_projector, &
    species_get_nlcc,          &
    species_nl_projector,      &
    species_get_iwf_radius,    &
    species_end


  integer, public, parameter :: &
    SPEC_USDEF  = 123,          & ! user-defined function
    SPEC_POINT  = 2,            & ! point charge: jellium sphere of radius 0.5 a.u.
    SPEC_JELLI  = 3,            & ! jellium sphere.
    SPEC_ALL_E  = 124,          & ! all-electron atom
    SPEC_CHARGE_DENSITY = 125,  &
    SPEC_FROM_FILE = 126,       &
    SPEC_PS_PSF = PS_TYPE_PSF,  & ! SIESTA pseudopotential
    SPEC_PS_HGH = PS_TYPE_HGH,  & ! HGH pseudopotential
    SPEC_PS_CPI = PS_TYPE_CPI,  & ! FHI pseudopotential (cpi format)
    SPEC_PS_FHI = PS_TYPE_FHI,  & ! FHI pseudopotential (ABINIT format)
    SPEC_PS_UPF = PS_TYPE_UPF     ! UPF pseudopotential

  type species_t
    private
    integer :: index              ! just a counter

    character(len=15) :: label    ! Identifier for the species
    integer :: type               ! what type of species
    FLOAT   :: z                  ! charge of the species
    FLOAT   :: z_val              ! valence charge of the species -- the total charge
                                  ! minus the core charge in the case of the pseudopotentials
    FLOAT   :: weight             ! mass, in atomic mass units (!= atomic units of mass)

    logical :: has_density        ! true if the species has an electronic density 

    ! for the user-defined potential
    character(len=1024) :: user_def
    FLOAT :: omega

    ! for the potential read from a file.
    character(len=200) :: filename

    ! jellium stuff
    FLOAT :: jradius

    ! For the pseudopotential
    type(ps_t), pointer :: ps
    logical             :: nlcc       ! true if we have non-local core corrections

    ! If we have an all-electron atom:
    FLOAT :: sigma

    ! If we have a charge distribution creating the potential:
    character(len=200) :: rho

    ! the default values for the spacing and atomic radius
    FLOAT :: def_rsize, def_h

    ! The number of initial wave functions
    integer :: niwfs
    integer, pointer :: iwf_l(:, :), iwf_m(:, :), iwf_i(:, :)

    ! For the TM pseudos, the lmax and lloc.
    integer :: lmax, lloc
  end type species_t


contains


  ! ---------------------------------------------------------
  ! Assigns a label to a species_t variable. This should be the
  ! first routine to be called (before species_index, species_read and species_init).
  ! This label must match one of the labels given in the %Species block
  ! in the input file -- or else one of the labels in the defaults file.
  ! ---------------------------------------------------------
  pure subroutine species_set_label(s, label)
    type(species_t), intent(inout) :: s
    character(len=*), intent(in)   :: label
    s%label = trim(label)
  end subroutine species_set_label
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Assigns an index to a species_t variable. This should be the
  ! second routine to be called (before species_read and species_init),
  ! when initializing a species_t variable.
  ! ---------------------------------------------------------
  pure subroutine species_set_index(s, k)
    type(species_t), intent(inout) :: s
    integer, intent(in)   :: k
    s%index = k
  end subroutine species_set_index
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Reads the information (from the input file) about a species_t variable, initializing
  ! part of it (it has to be completed later with "species_init").
  ! Note that species_read has to be called only after species_set_label
  ! and species_set index have been called.
  ! ---------------------------------------------------------
  subroutine species_read(s)
    type(species_t), intent(inout) :: s

    character(len=256) :: fname
    character(len=10)  :: lab
    integer :: i, row, n_spec_block, n_spec_def, iunit, read_data
    type(block_t) :: blk

    call push_sub('species.species_read')

    s%has_density = .false. ! there is no density associated
    s%nlcc      = .false.   ! without non-local core corrections
    s%def_h     = -M_ONE    ! not defined
    s%def_rsize = -M_ONE    ! not defined
    s%user_def  = ""
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
    !% not required to define them explicitly in the Species block, as defaults 
    !% are provided by the program (you can override these defaults in any case). 
    !% Additional pseudopotentials can be downloaded from the 
    !% <a href='http://www.tddft.org/programs/octopus/wiki/index.php/Pseudopotentials'>
    !% octopus homepage</a>.
    !%
    !% The format of this block is the following: The first field is
    !% the name of the species, followed by the atomic mass (in atomic mass
    !% units). The third field defines the type of species (the valid options
    !% are detailed below), each type needs some extra parameters given in
    !% the following fields of the row.
    !%
    !% In 3D, <i>e.g.</i>
    !%
    !% <tt>%Species
    !% <br>&nbsp;&nbsp;'O'       | 15.9994 | spec_ps_psf  | 8   | 1 | 1
    !% <br>&nbsp;&nbsp;'H'       |  1.0079 | spec_ps_hgh  | 1   | 0 | 0
    !% <br>&nbsp;&nbsp;'jlm'     | 23.2    | spec_jelli   | 8   | 5.0
    !% <br>&nbsp;&nbsp;'pnt'     | 32.3    | spec_point   | 2.0
    !% <br>&nbsp;&nbsp;'udf'     |  0.0    | user_defined | 8   | "1/2*r^2"
    !% <br>&nbsp;&nbsp;'H_all'   |  1.0079 | spec_all_e   | 1   
    !% <br>%</tt>
    !%
    !% Additionally, all the pseudopotential types (PSF, HGH, FHI, UPF) can take two extra
    !% fields: default spacing, and default radius (used for minimum simulation box if the
    !% radius is not specified). 
    !%
    !%Option user_defined  123
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
    !%Option spec_all_e   124
    !% Atom represented with all electrons; the extra parameter is the
    !% atomic number. See the documentation of the variable 
    !% <tt>SpeciesAllElectronSigma</tt>.
    !% WARNING: Currently you can not use LCAO with this species.
    !%Option spec_charge_density 125
    !% The potential is created by a distribution of charge.
    !%Option species_from_file  126
    !% The potential is read from a file, whose name is given in column 5.
    !%End

    call obsolete_variable('SpecieAllElectronSigma', 'SpeciesAllElectronSigma')

    !%Variable SpeciesAllElectronSigma
    !%Type float
    !%Default 0.25
    !%Section System::Species
    !%Description
    !% An all-electron atom is defined by a Gaussian accumulation of positive charge 
    !% (distorted if curvilinear coordinates are used), in the form:
    !%
    !% <math>
    !% q(r) = z * \beta * exp[ - (\vec{r}-\vec{r0})**2 / (sqrt(2) * \delta * \sigma) ]
    !% </math>
    !%
    !% <math>\beta</math> is chosen in order to maintain proper normalization (the integral of
    !% <math>q</math> should sum up to <math>z</math>). <math>\delta</math> is 
    !% the grid spacing (the grid spacing in the first dimension, to be precise). 
    !% <math>\vec{r0}</math> is calculated in such a way that the the first moment of 
    !% <math>q(r)/z</math> is equal to the atomic position. This number should be small,
    !% but you may run into numerical difficulties if it is too small.
    !%
    !% For a precise description, see N. A. Modine, <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997).
    !%End
    call loct_parse_float(datasets_check('SpeciesAllElectronSigma'), CNST(0.25), s%sigma)
    if(s%sigma <= M_ZERO) call input_error('SpeciesAllElectronSigma')

    ! First, find out if there is a Species block.
    n_spec_block = 0
    if(loct_parse_block(datasets_check('Species'), blk) == 0) then
      n_spec_block = loct_parse_block_n(blk)
    end if

    ! Find out if the sought species is in the block
    row = -1
    block: do i = 1, n_spec_block
      call loct_parse_block_string(blk, i-1, 0, lab)
      if(trim(lab)==trim(s%label)) then
        row = i - 1
        exit block
      end if
    end do block

    ! Read whatever may be read from the block
    if(row>=0) then
      call read_from_block(blk, row, s, read_data)
      call loct_parse_block_end(blk)
    end if

    ! Find out if the species is in the defaults file.
    write(fname, '(2a)') trim(conf%share), "/PP/defaults"
    n_spec_def = max(0, loct_number_of_lines(fname))
    if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

    iunit = io_open(fname, action='read', status='old', die=.false., is_tmp=.true.)
    if(iunit > 0) then
      read(iunit,*)

      default_file: do i = 1, n_spec_def
        read(iunit,*) lab
        if(trim(lab) == trim(s%label)) then
          call read_from_default_file(iunit, read_data, s)
          exit default_file
        end if
      end do default_file

      call io_close(iunit)
    end if

    if(read_data == 0) then
      message(1) = 'Species '//trim(s%label)//' not found.'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine species_read
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_init(s, ispin, print_info)
    type(species_t),   intent(inout) :: s
    integer,           intent(in)    :: ispin
    logical, optional, intent(in)    :: print_info

    logical :: print_info_
    integer :: i, l
    FLOAT   :: pot_re, pot_im, xx(1)

    call push_sub('species.species_init')

    print_info_ = .true.
    if(present(print_info)) then
      print_info_ = print_info
    end if

    ! masses are always in amu, so convert them to a.u.
    s%weight =  units_to_atomic(units_inp%mass, s%weight)

    s%has_density = .false.

    select case(s%type)
    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
      ! allocate structure
      SAFE_ALLOCATE(s%ps) 
      call ps_init(s%ps, s%label, s%type, s%Z, s%lmax, s%lloc, ispin)
      s%z_val = s%ps%z_val
      s%nlcc = (s%ps%icore /= 'nc  ' )
      s%niwfs = 0
      do i = 1, s%ps%conf%p
        l = s%ps%conf%l(i)

        ! the input pseudo determines the maximum l we use
        if(l <= s%lmax) s%niwfs = s%niwfs + (2*l+1)

        ! Another choice would be to use only the occupied atomic states,
        ! usually yielding a slightly smaller basis. To use this strategy
        ! uncomment the following line
        !if(sum(s%ps%conf%occ(i, :)).ne.M_ZERO) s%niwfs = s%niwfs + (2*l+1)
      end do

    case(SPEC_USDEF)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(s%label),'" is a user-defined potential.'
        i = min(237, len_trim(s%user_def)-1) ! I substract 1 to avoid the non-printable C "end-of-string" character.
        write(message(2),'(a,a)')      '   Potential = ', trim(s%user_def(1:i))
        if(len(trim(s%user_def)).gt.237) then
          message(2) = trim(message(2))//'...'
        end if
        call write_info(2)
      end if
      s%niwfs = int(max(2*s%z_val, CNST(1.0)))
      xx(1) = CNST(0.01)
      call loct_parse_expression(pot_re, pot_im, 1, xx, xx(1), M_ZERO, s%user_def)
      s%omega = sqrt( abs(M_TWO / CNST(1.0e-4) * pot_re ))
      ! To avoid problems with constant potentials.
      if(s%omega <= M_ZERO) s%omega = CNST(0.1) 

    case(SPEC_FROM_FILE)
      if(print_info_) then
        write(message(1),'(a)') 'Species read from file "'//trim(s%filename)//'".'
        call write_info(1)
      end if
      s%omega = CNST(0.1)

    case(SPEC_JELLI, SPEC_POINT)
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(s%label),'" is a jellium sphere / approximated point particle.'
        write(message(2),'(a,f11.6)')  '   Valence charge = ', s%z_val
        write(message(3),'(a,f11.6)')  '   Radius [a.u]   = ', s%jradius
        write(message(4),'(a,f11.6)')  '   Rs [a.u]       = ', s%jradius * s%z_val ** (-M_ONE/M_THREE)
        call write_info(4)
      end if
      s%niwfs = 2*s%z_val
      s%omega = CNST(0.1)

    case(SPEC_ALL_E)
      s%niwfs = 2*s%z_val
      s%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(s%label),'" is an all-electron atom.'
        write(message(2),'(a,f11.6)')  '   Z = ', s%z_val
        write(message(3),'(a)')  '   Potential will be calculated solving Poisson equation'
        write(message(4),'(a)')  '   for a delta density distribution.'
        call write_info(4)
      end if

    case(SPEC_CHARGE_DENSITY)
      s%niwfs = int(max(2*s%z_val, CNST(1.0)))
      s%has_density = .true.
      if(print_info_) then
        write(message(1),'(a,a,a)')    'Species "',trim(s%label),'" is a distribution of charge:'
        write(message(2),'(a,a)')      '   rho = ', trim(s%rho)
        write(message(3),'(a,f11.6)')  '   Z = ', s%z_val
        call write_info(3)
      end if
    end select

    SAFE_ALLOCATE(s%iwf_l(1:s%niwfs, 1:ispin))
    SAFE_ALLOCATE(s%iwf_m(1:s%niwfs, 1:ispin))
    SAFE_ALLOCATE(s%iwf_i(1:s%niwfs, 1:ispin))
    call species_iwf_fix_qn(s, ispin)

    call pop_sub()
  end subroutine species_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! This routine performs some operations on the pseudopotential
  ! functions (filtering, etc), some of which depend on the grid
  ! cutoff value.
  ! ---------------------------------------------------------
  subroutine species_pot_init(this, grid_cutoff, filter)
    type(species_t),     intent(inout) :: this
    FLOAT,               intent(in)    :: grid_cutoff
    integer,             intent(in)    :: filter

    character(len=256) :: dirname

    call push_sub('species.species_pot_init')
    
    if(species_is_ps(this)) then
      call ps_separate(this%ps)

      call ps_getradius(this%ps)

      if(filter .ne. PS_FILTER_NONE) then 
        call ps_filter(this%ps, filter, grid_cutoff)
        call ps_getradius(this%ps) ! radius may have changed
      end if

      call ps_derivatives(this%ps)
       
      if(in_debug_mode) then
        write(dirname, '(a)') 'debug/geometry'
        call io_mkdir(dirname)
        call species_debug(trim(dirname), this)
      end if
    end if

    call pop_sub()
  end subroutine species_pot_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function species_type(s)
    type(species_t), intent(in) :: s
    species_type = s%type
  end function species_type
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=15) pure function species_label(s)
    type(species_t), intent(in) :: s
    species_label = trim(s%label)
  end function species_label
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function species_index(s)
    type(species_t), intent(in) :: s
    species_index = s%index
  end function species_index
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_has_nlcc(s)
    type(species_t), intent(in) :: s
    species_has_nlcc = s%nlcc
  end function species_has_nlcc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical pure function species_has_density(s)
    type(species_t), intent(in) :: s
    species_has_density = s%has_density
  end function species_has_density
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  function species_ps(s)
    type(ps_t), pointer :: species_ps
    type(species_t), intent(in) :: s
    species_ps => s%ps
  end function species_ps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_zval(s)
    type(species_t), intent(in) :: s
    species_zval = s%z_val
  end function species_zval
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_z(s)
    type(species_t), intent(in) :: s
    species_z = s%z
  end function species_z
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_def_rsize(s)
    type(species_t), intent(in) :: s
    species_def_rsize = s%def_rsize
  end function species_def_rsize
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_def_h(s)
    type(species_t), intent(in) :: s
    species_def_h = s%def_h
  end function species_def_h
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_jradius(s)
    type(species_t), intent(in) :: s
    species_jradius = s%jradius
  end function species_jradius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_sigma(s)
    type(species_t), intent(in) :: s
    species_sigma = s%sigma
  end function species_sigma
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_omega(s)
    type(species_t), intent(in) :: s
    species_omega = s%omega
  end function species_omega
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT pure function species_weight(s)
    type(species_t), intent(in) :: s
    species_weight = s%weight
  end function species_weight
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=200) pure function species_rho_string(s)
    type(species_t), intent(in) :: s
    species_rho_string = trim(s%rho)
  end function species_rho_string
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  character(len=200) pure function species_filename(s)
    type(species_t), intent(in) :: s
    species_filename = trim(s%filename)
  end function species_filename
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer pure function species_niwfs(s)
    type(species_t), intent(in) :: s
    species_niwfs = s%niwfs
  end function species_niwfs
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  pure subroutine species_iwf_ilm(s, j, is, i, l, m)
    type(species_t), intent(in) :: s
    integer, intent(in)         :: j, is
    integer, intent(out)        :: i, l, m
    i = s%iwf_i(j, is)
    l = s%iwf_l(j, is)
    m = s%iwf_m(j, is)
  end subroutine species_iwf_ilm
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function species_userdef_pot(s, dim, xx, r, t)
    type(species_t), intent(in) :: s
    integer, intent(in) :: dim
    FLOAT, intent(in) :: xx(1:MAX_DIM), r, t
    FLOAT :: pot_re, pot_im
    call loct_parse_expression(                            &
               pot_re, pot_im, dim, xx, r, t, s%user_def)
    species_userdef_pot = pot_re
  end function species_userdef_pot
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical function species_is_ps(s)
    type(species_t), intent(in) :: s
    
    call push_sub('species.species_is_ps')

    species_is_ps = &
         ( s%type == SPEC_PS_PSF) .or. &
         ( s%type == SPEC_PS_HGH) .or. &
         ( s%type == SPEC_PS_CPI) .or. &
         ( s%type == SPEC_PS_FHI) .or. &
         ( s%type == SPEC_PS_UPF)
    
    call pop_sub()
  end function species_is_ps
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical function species_is_local(s)
    type(species_t), intent(in) :: s

    call push_sub('species.species_is_local')

    species_is_local = .true.
      
    if( species_is_ps(s) ) then 
      if ( s%ps%l_max /= 0 ) species_is_local = .false. 
    end if

    call pop_sub()
  end function species_is_local
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! This routine returns the non-local projector and its 
  ! derivative build using real spherical harmonics
  subroutine species_real_nl_projector(s, x, l, lm, i, uV, duV)
    type(species_t),   intent(in)  :: s
    FLOAT,             intent(in)  :: x(1:MAX_DIM)
    integer,           intent(in)  :: l, lm, i
    FLOAT,             intent(out) :: uV, duV(1:MAX_DIM)

    FLOAT :: r, uVr0, duvr0, ylm, gylm(MAX_DIM)
    FLOAT, parameter :: ylmconst = CNST(0.488602511902920) !  = sqrt(3/(4*pi))

    r = sqrt(sum(x(1:MAX_DIM)**2))

    uVr0  = spline_eval(s%ps%kb(l, i), r)
    duVr0 = spline_eval(s%ps%dkb(l, i), r)

    call grylmr(x(1), x(2), x(3), l, lm, ylm, gylm)
    uv = uvr0*ylm
    if(r >= r_small) then
      duv(:) = duvr0 * ylm * x(:)/r + uvr0 * gylm(:)
    else
      if(l == 1) then
        duv = M_ZERO
        if(lm == -1) then
          duv(2) = -ylmconst * duvr0
        else if(lm == 0) then
          duv(3) =  ylmconst * duvr0
        else if(lm == 1) then
          duv(1) = -ylmconst * duvr0
        end if
      else
        duv = M_ZERO
      end if
    end if

  end subroutine species_real_nl_projector
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! This routine returns the non-local projector build using 
  ! spherical harmonics
  subroutine species_nl_projector(s, x, l, lm, i, uV)
    type(species_t),    intent(in)  :: s
    FLOAT,             intent(in)  :: x(1:MAX_DIM)
    integer,           intent(in)  :: l, lm, i
    CMPLX,             intent(out) :: uV

    FLOAT :: r, uVr0
    CMPLX :: ylm

    r = sqrt(sum(x(1:MAX_DIM)**2))

    uVr0 = spline_eval(s%ps%kb(l, i), r)

    call ylmr(x(1), x(2), x(3), l, lm, ylm)
    uv = uvr0*ylm

  end subroutine species_nl_projector
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  FLOAT function species_get_nlcc(s, x) result(l)
    type(species_t), intent(in) :: s
    FLOAT, intent(in) :: x(MAX_DIM)
    FLOAT :: r

    ! only for 3D pseudopotentials, please
    if(species_is_ps(s)) then
      r = sqrt(sum(x(:)**2))
      l = spline_eval(s%ps%core, r)
    else
      l = M_ZERO
    end if

  end function species_get_nlcc
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  FLOAT function species_get_iwf_radius(s, j, is) result(radius)
    type(species_t),   intent(in) :: s
    integer,           intent(in) :: j
    integer,           intent(in) :: is

    integer :: i, l, m
    FLOAT, parameter :: threshold = CNST(0.001)

    call push_sub('species.species_get_iwf_radius')

    i = s%iwf_i(j, is)
    l = s%iwf_l(j, is)
    m = s%iwf_m(j, is)

    if(species_is_ps(s)) then
      radius = spline_cutoff_radius(s%ps%ur(i, is), threshold)
    else
      radius = sqrt(-M_TWO*log(threshold)/s%omega)
    end if

    call pop_sub()

  end function species_get_iwf_radius
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_end(ns, s)
    integer,        intent(in) :: ns
    type(species_t), pointer    :: s(:)

    integer :: i

    call push_sub('species.species_end')

    do i = 1, ns
      if (species_is_ps(s(i))) then 
        if(associated(s(i)%ps)) then 
          call ps_end(s(i)%ps)
          SAFE_DEALLOCATE_P(s(i)%ps)
        end if
      end if
      SAFE_DEALLOCATE_P(s(i)%iwf_l)
      SAFE_DEALLOCATE_P(s(i)%iwf_m)
      SAFE_DEALLOCATE_P(s(i)%iwf_i)
    end do

    call pop_sub()
  end subroutine species_end
  ! ---------------------------------------------------------





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private procedures

  ! ---------------------------------------------------------
  subroutine species_debug(dir, s)
    character(len=*), intent(in)  :: dir
    type(species_t),    intent(in) :: s

    character(len=256) :: dirname
    integer :: iunit
    logical :: bool

    if(.not.mpi_grp_is_root(mpi_world)) then
      call write_debug_newlines(4)
      return
    end if

    call push_sub('species.species_debug')

    dirname = trim(dir)//'/'//trim(s%label)

    call io_mkdir(dirname)

    iunit = io_open(trim(dirname)//'/info', action='write')

    write(iunit, '(a,i3)')    'Index  = ', s%index
    write(iunit, '(2a)')      'Label  = ', trim(s%label)
    write(iunit, '(a,i3)')    'Type   = ', s%type
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,f15.2)') 'z      = ', s%z
    if (s%type == SPEC_FROM_FILE) then
      write(iunit,'(a)')      'Specie read from file "'//trim(s%filename)//'".'
    end if
    write(iunit, '(a,f15.2)') 'z_val  = ', s%z_val
    write(iunit, '(a,f15.2)') 'weight = ', s%weight
    bool = species_is_local(s)
    write(iunit, '(a,l1)')    'local  = ', bool
    write(iunit, '(2a)')      'usdef  = ', trim(s%user_def)
    if (s%type == SPEC_JELLI .or. s%type == SPEC_POINT) then
      write(iunit, '(a,f15.2)') 'jradius= ', s%jradius
    end if
    write(iunit, '(a,l1)')    'nlcc   = ', s%nlcc
    write(iunit, '(a,f15.2)') 'def_rsize = ', s%def_rsize
    write(iunit, '(a,f15.2)') 'def_h = ', s%def_h
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lmax  = ', s%lmax
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lloc  = ', s%lloc

    if(species_is_ps(s)) then
       if(in_debug_mode) call ps_debug(s%ps, trim(dirname))
    end if

    call io_close(iunit)
    call pop_sub()
  end subroutine species_debug
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine read_from_default_file(iunit, read_data, s)
    integer,        intent(in)    :: iunit
    integer,        intent(inout) :: read_data
    type(species_t), intent(inout) :: s

    character(len=10) :: label
    FLOAT :: weight, z, def_h, def_rsize
    integer :: lloc, lmax
    integer :: type

    call push_sub('species.read_from_default_file')

    backspace(iunit)

    read(iunit,*) label, weight, type, z, lmax, lloc, def_h, def_rsize
    def_h     = def_h     * P_ANG  ! These units are always in Angstrom
    def_rsize = def_rsize * P_ANG
    ASSERT(trim(label) == trim(s%label))

    if(read_data == 0) then ! The Species was not supplied in the block.
      s%weight    = weight
      s%type      = type
    end if
    if(read_data < 4) s%z         = z
    if(read_data < 5) s%lmax      = lmax
    if(read_data < 6) s%lloc      = lloc
    if(read_data < 7) s%def_h     = def_h
    if(read_data < 8) s%def_rsize = def_rsize

    read_data = 8

    call pop_sub()
  end subroutine read_from_default_file
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine read_from_block(blk, row, s, read_data)
    type(block_t),      intent(in) :: blk
    integer,        intent(in) :: row
    type(species_t), intent(inout) :: s
    integer,        intent(out) :: read_data

    integer :: n

    call push_sub('species.read_from_block')

    read_data = 0

    call loct_parse_block_float (blk, row, 1, s%weight)
    call loct_parse_block_int   (blk, row, 2, s%type)

    select case(s%type)
    case(SPEC_USDEF) ! user-defined
      call loct_parse_block_float (blk, row, 3, s%Z_val)
      call loct_parse_block_string(blk, row, 4, s%user_def)
      call conv_to_C_string(s%user_def)
      read_data = 5

    case(SPEC_FROM_FILE)
      call loct_parse_block_float (blk, row, 3, s%Z_val)
      call loct_parse_block_string(blk, row, 4, s%filename)
      read_data = 5

    case(SPEC_POINT) ! this is treated as a jellium with radius 0.5
      call loct_parse_block_float(blk, row, 3, s%Z)
      s%jradius = M_HALF
      s%Z_val = 0
      read_data = 4

    case(SPEC_JELLI)
      call loct_parse_block_float(blk, row, 3, s%Z)      ! charge of the jellium sphere
      call loct_parse_block_float(blk, row, 4, s%jradius)! radius of the jellium sphere
      s%jradius = units_inp%length%factor * s%jradius    ! units conversion
      s%Z_val = s%Z
      read_data = 5

    case(SPEC_ALL_E)
      call loct_parse_block_float(blk, row, 3, s%Z)
      s%Z_val = s%Z
      read_data = 4

    case(SPEC_CHARGE_DENSITY)
      call loct_parse_block_float(blk, row, 3, s%Z)
      call loct_parse_block_string(blk, row, 4, s%rho)
      s%Z_val = s%Z
      read_data = 5

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF) ! a pseudopotential file
      n = loct_parse_block_cols(blk, row)

      call loct_parse_block_float (blk, row, 3, s%Z)

      if(s%type == SPEC_PS_UPF) read_data = 4

      if(n>4) then
        call loct_parse_block_int (blk, row, 4, s%lmax)
        read_data = 5
      end if

      if(n>5) then
        call loct_parse_block_int (blk, row, 5, s%lloc)
        read_data = 6
      end if

      if(n>6) then
        call loct_parse_block_float (blk, row, 6, s%def_h)
        s%def_h = s%def_h * units_inp%length%factor
        read_data = 7
      end if

      if(n>7) then
        call loct_parse_block_float (blk, row, 7, s%def_rsize)
        s%def_rsize = s%def_rsize * units_inp%length%factor
        read_data = 8
      end if

    case default
      write(message(1), '(a,i2,a)') "Unknown pseudopotential type: '", s%type, "'"
      call write_fatal(1)
    end select

    call pop_sub()
  end subroutine read_from_block
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_iwf_fix_qn(s, ispin)
    type(species_t), intent(inout) :: s
    integer, intent(in) :: ispin

    integer :: is, n, i, l, m, n1, n2, n3

    call push_sub('species.species_iwf_fix_qn')

    if(species_is_ps(s)) then
      do is = 1, ispin
        n = 1
        do i = 1, s%ps%conf%p
          l = s%ps%conf%l(i)

          ! the input pseudo determines the maximum l we use
          if(l > s%lmax) cycle

          do m = -l, l
            s%iwf_i(n, is) = i
            s%iwf_l(n, is) = l
            s%iwf_m(n, is) = m
            n = n + 1
          end do
        end do
      end do
    else

      select case(calc_dim)
      case(1)
        do is = 1, ispin
          do i = 1, s%niwfs
            s%iwf_i(i, is) = i
            s%iwf_l(i, is) = 0
            s%iwf_m(i, is) = 0
          end do
        end do

      case(2)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1
          do
            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = 0
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1+1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = 0
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2+1
            s%iwf_m(i, is) = 0
            i = i + 1; if(i>s%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1
          end do
        end do

      case(3)
        do is = 1, ispin
          i = 1; n1 = 1; n2 = 1; n3 = 1
          do
            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = n3
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1+1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = n3
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2+1
            s%iwf_m(i, is) = 0
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = n3+1
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1+1
            s%iwf_l(i, is) = n2+1
            s%iwf_m(i, is) = n3
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1+1
            s%iwf_l(i, is) = n2
            s%iwf_m(i, is) = n3+1
            i = i + 1; if(i>s%niwfs) exit

            s%iwf_i(i, is) = n1
            s%iwf_l(i, is) = n2+1
            s%iwf_m(i, is) = n3+1
            i = i + 1; if(i>s%niwfs) exit

            n1 = n1 + 1; n2 = n2 + 1; n3 = n3 + 1
          end do
        end do
      end select
    end if

    call pop_sub()
  end subroutine species_iwf_fix_qn
  ! ---------------------------------------------------------




end module species_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
