!! Copyright (C) 2002-2012 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module ps_oct_m
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use loct_math_oct_m
  use parser_oct_m
  use logrid_oct_m
  use messages_oct_m
  use profiling_oct_m
  use ps_cpi_oct_m
  use ps_fhi_oct_m
  use ps_hgh_oct_m
  use ps_xml_oct_m
  use ps_in_grid_oct_m
#ifdef HAVE_PSPIO
  use fpspio_m
#endif
  use ps_psf_oct_m
  use ps_upf_oct_m
  use splines_oct_m
  use spline_filter_oct_m
  implicit none

  private
  public ::                     &
    ps_t,                       &
    ps_init,                    &
    ps_pspio_init,              &
    ps_separate,                &
    ps_filter,                  &
    ps_getradius,               &
    ps_derivatives,             &
    ps_debug,                   &
    ps_niwfs,                   &
    ps_end,                     &
    ps_type,                    &
    ps_has_density,             &
    ps_density_volume
  
  integer, parameter, public :: &
    PS_TYPE_PSF = 100,          &
    PS_TYPE_HGH = 101,          &
    PS_TYPE_CPI = 102,          &
    PS_TYPE_FHI = 103,          &
    PS_TYPE_UPF = 104,          &
    PS_TYPE_XML = 105

  integer, parameter, public :: &
    PS_FILTER_NONE = 0,         &
    PS_FILTER_TS   = 2,         &
    PS_FILTER_BSB  = 3

  integer, parameter, public :: HUGE_L = 100
  
  character(len=3), parameter  :: ps_name(PS_TYPE_PSF:PS_TYPE_XML) = (/"tm2", "hgh", "cpi", "fhi", "upf", "qso"/)

  type ps_t
    character(len=10) :: label
    integer           :: flavour

    integer  :: ispin    !< Consider spin (ispin = 2) or not (ispin = 1)
    FLOAT    :: z, z_val
    type(valconf_t)   :: conf
    type(logrid_t) :: g
    type(spline_t), pointer :: ur(:, :)     !< (1:conf%p, 1:ispin) atomic wavefunctions, as a function of r
    type(spline_t), pointer :: ur_sq(:, :)  !< (1:conf%p, 1:ispin) atomic wavefunctions, as a function of r^2

    ! Kleinman-Bylander projectors stuff
    integer  :: lmax    !< maximum value of l to take
    integer  :: llocal  !< which component to take as local

    type(spline_t) :: vl         !< local part

    FLOAT :: projectors_sphere_threshold !< The projectors are localized in real
                                         !! space, and so they are contained in a
                                         !! sphere whose radius is computed by
                                         !! making sure that the projector
                                         !! functions absolute value is below this
                                         !! threshold, for points outside the
                                         !! sphere.
    FLOAT :: rc_max !< The radius of the spheres that contain the projector functions.

    integer  :: kbc      !< Number of KB components (1 or 2 for TM ps, 3 for HGH)
    FLOAT, pointer :: h(:,:,:), k(:, :, :)
    type(spline_t), pointer :: kb(:, :)     !< Kleinman-Bylander projectors
    type(spline_t), pointer :: dkb(:, :)    !< derivatives of KB projectors

    logical :: nlcc    !< .true. if the pseudo has non-linear core corrections.
    type(spline_t) :: core !< normalization \int dr 4 pi r^2 rho(r) = N


    !LONG-RANGE PART OF THE LOCAL POTENTIAL
    
    logical :: has_long_range

    type(spline_t) :: vlr         !< the long-range part of the local potential
    type(spline_t) :: vlr_sq      !< the long-range part of the
                                  !< local potential in terms of r^2, to avoid the sqrt
    type(spline_t) :: nlr         !< the charge density associated with the long-range part

    FLOAT :: sigma_erf            !< the a constant in erf(r/(sqrt(2)*sigma))/r

    logical :: has_density                     !< does the species have a density?
    type(spline_t), pointer :: density(:)      !< the atomic density for each spin
    type(spline_t), pointer :: density_der(:)  !< the radial derivative for the atomic density for each spin
    
    logical :: is_separated
    logical :: local
    logical :: hamann
  end type ps_t

  FLOAT, parameter :: eps = CNST(1.0e-8)

contains


  ! ---------------------------------------------------------
  subroutine ps_init(ps, label, z, user_lmax, user_llocal, ispin, filename)
    type(ps_t),        intent(out)   :: ps
    character(len=10), intent(in)    :: label
    integer,           intent(in)    :: user_lmax
    integer,           intent(in)    :: user_llocal
    integer,           intent(in)    :: ispin
    FLOAT,             intent(in)    :: z
    character(len=*),  intent(in)    :: filename
    
    integer :: l, ii, ll, is, ierr
    type(ps_psf_t) :: ps_psf !< SIESTA pseudopotential
    type(ps_cpi_t) :: ps_cpi !< Fritz-Haber pseudopotential
    type(ps_fhi_t) :: ps_fhi !< Fritz-Haber pseudopotential (from abinit)
    type(ps_upf_t) :: ps_upf !< In case UPF format is used
    type(ps_hgh_t) :: ps_hgh !< In case Hartwigsen-Goedecker-Hutter ps are used.
    type(ps_xml_t) :: ps_xml !< For xml based pseudopotentials

    PUSH_SUB(ps_init)

    ! Fix the threshold to calculate the radius of the projector-function localization spheres:

    call messages_obsolete_variable('SpecieProjectorSphereThreshold', 'SpeciesProjectorSphereThreshold')

    !%Variable SpeciesProjectorSphereThreshold
    !%Type float
    !%Default 0.001
    !%Section System::Species
    !%Description
    !% The pseudopotentials may be composed of a local part, and a linear combination of nonlocal
    !% operators. These nonlocal projectors have "projector" form, <math> \left| v \right> \left< v \right| </math>
    !% (or, more generally speaking, <math> \left| u \right> \left< v \right| </math>).
    !% These projectors are localized in real space -- that is, the function <math>v</math>
    !% has a finite support around the nucleus. This region where the projectors are localized should
    !% be small or else the computation time required to operate with them will be very large.
    !% 
    !% In practice, this localization is fixed by requiring the definition of the projectors to be
    !% contained in a sphere of a certain radius. This radius is computed by making sure that the 
    !% absolute value of the projector functions, at points outside the localization sphere, is 
    !% below a certain threshold. This threshold is set by <tt>SpeciesProjectorSphereThreshold</tt>.
    !%End
    call parse_variable('SpeciesProjectorSphereThreshold', CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call messages_input_error('SpeciesProjectorSphereThreshold')
   
    ! Sets the flavour, label, and number of spin channels.
    ps%flavour = ps_get_type(filename)
    ps%label   = label
    ps%ispin   = ispin
    ps%hamann  = .false.
    select case(ps%flavour)
    case(PS_TYPE_PSF, PS_TYPE_HGH)
      ps%has_density = .true.
    case default
      ps%has_density = .false.
    end select
    
    if(.not. (ps%flavour >= PS_TYPE_PSF .and. ps%flavour <= PS_TYPE_XML)) then
      call messages_write("Cannot determine the pseudopotential type for species '"//trim(label)//"' from", new_line = .true.)
      call messages_write("file '"//trim(filename)//"'.")
      call messages_fatal()
    end if

    select case(ps%flavour)
    case(PS_TYPE_PSF)
      call ps_psf_init(ps_psf, ispin, filename)

      call valconf_copy(ps%conf, ps_psf%conf)
      ps%z      = z
      ps%conf%z = nint(z) ! atomic number
      ps%kbc    = 1     ! only one projector per angular momentum
      
      ps%lmax = ps_psf%ps_grid%no_l_channels - 1

      if(user_lmax /= HUGE_L) then
        ps%lmax = min(ps%lmax, user_lmax) ! Maybe the file does not have enough components.
        if(user_lmax /= ps%lmax) then
          message(1) = "lmax in Species block for " // trim(label) // " is larger than number available in pseudopotential."
          call messages_fatal(1)
        end if
      end if

      ps%conf%p = ps_psf%ps_grid%no_l_channels
      if(ps%lmax == 0) ps%llocal = 0 ! Vanderbilt is not acceptable if ps%lmax == 0.

      ! the local part of the pseudo
      if(user_llocal == HUGE_L) then
        ps%llocal = 0
      else
        ps%llocal = user_llocal
      end if
      
      call ps_psf_process(ps_psf, ps%lmax, ps%llocal)
      call logrid_copy(ps_psf%ps_grid%g, ps%g)

    case(PS_TYPE_CPI, PS_TYPE_FHI)
      call valconf_null(ps%conf)

      if(ps%flavour == PS_TYPE_CPI) then
        call ps_cpi_init(ps_cpi, trim(filename))
        ps%conf%p      = ps_cpi%ps_grid%no_l_channels
      else
        call ps_fhi_init(ps_fhi, trim(filename))
        ps%conf%p      = ps_fhi%ps_grid%no_l_channels
      end if

      ps%conf%z      = nint(z)
      ps%conf%symbol = label(1:2)
      ps%conf%type   = 1
      do l = 1, ps%conf%p
        ps%conf%l(l) = l - 1
      end do

      ps%z      = z
      ps%kbc    = 1     ! only one projector per angular momentum

      ps%lmax  = ps%conf%p - 1

      if(user_lmax /= HUGE_L) then
        ps%lmax = min(ps%lmax, user_lmax) ! Maybe the file does not have enough components.
        if(user_lmax /= ps%lmax) then
          message(1) = "lmax in Species block for " // trim(label) // " is larger than number available in pseudopotential."
          call messages_fatal(1)
        end if
      end if

      if(ps%lmax == 0) ps%llocal = 0 ! Vanderbilt is not acceptable if ps%lmax == 0.

      ! the local part of the pseudo
      if(user_llocal == HUGE_L) then
        ps%llocal = 0
      else
        ps%llocal = user_llocal
      end if
      
      if(ps%flavour == PS_TYPE_CPI) then
        call ps_cpi_process(ps_cpi, ps%llocal)
        call logrid_copy(ps_cpi%ps_grid%g, ps%g)
      else
        call ps_fhi_process(ps_fhi, ps%lmax, ps%llocal)
        call logrid_copy(ps_fhi%ps_grid%g, ps%g)
      end if

    case(PS_TYPE_HGH)
      call hgh_init(ps_hgh, trim(filename))
      call valconf_copy(ps%conf, ps_hgh%conf)

      ps%z        = z
      ps%kbc      = 3
      ps%llocal    = -1
      ps%lmax    = ps_hgh%l_max

      call hgh_process(ps_hgh)
      call logrid_copy(ps_hgh%g, ps%g)

    case(PS_TYPE_XML, PS_TYPE_UPF)
      
      call messages_experimental('XML (QSO, UPF2, and PSML) pseudopotential support')
      
      call ps_xml_init(ps_xml, trim(filename), ierr)

      if(ierr == 0) then
        
        call valconf_null(ps%conf)
        
        ps%z      = z
        ps%conf%z = nint(z)
        
        if(ps_xml%kleinman_bylander) then
          ps%conf%p = ps_xml%nwavefunctions
        else
          ps%conf%p = ps_xml%lmax + 1
        end if
        
        do ll = 0, ps_xml%lmax
          ps%conf%l(ll + 1) = ll
        end do
        
        ps%kbc    = ps_xml%nchannels
        ps%lmax  = ps_xml%lmax

        if(ps_xml%kleinman_bylander) then
          ps%llocal = ps_xml%llocal
        else
          ! we have several options
          ps%llocal = user_llocal                     ! user supplied local component
          if(ps%llocal < 0) ps%llocal = ps_xml%llocal ! the one given in the pseudopotential file
          if(ps%llocal < 0) ps%llocal = ps_xml%lmax   ! we use the maximum l possible as local
          ASSERT(ps%llocal >= 0)
          ASSERT(ps%llocal <= ps%lmax)
        end if
        
        nullify(ps%g%drdi, ps%g%s)
        
        ! use a larger grid
        ps%g%nrval = max(ps_xml%grid_size, nint(CNST(20.0)/(ps_xml%mesh_spacing)))
        
        SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
        SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))
        
        do ii = 1, ps%g%nrval
          ps%g%rofi(ii) = (ii - 1)*ps_xml%mesh_spacing
          ps%g%r2ofi(ii) = ps%g%rofi(ii)**2
        end do

        ps%flavour = PS_TYPE_XML
        
      else !read failed, this must be a UPF 1 file
        
        call ps_upf_init(ps_upf, trim(filename))
        
        call valconf_copy(ps%conf, ps_upf%conf)
        ps%z      = z
        ps%conf%z = nint(z)
        ps%kbc    = ps_upf%kb_nc
        ps%lmax  = ps_upf%l_max
        ps%llocal  = ps_upf%l_loc
        ps%has_density = .true.
        
        nullify(ps%g%drdi, ps%g%s)
        ps%g%nrval = ps_upf%np
        SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
        SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))
        ps%g%rofi = ps_upf%r
        ps%g%r2ofi = ps%g%rofi**2
        
      end if
      
    end select

    write(message(1), '(a,i2,a)') "Info: l = ", ps%lmax, " is maximum angular momentum considered."
    call messages_info(1)

    ps%local = (ps%lmax == 0 .and. ps%llocal == 0 ) .or. (ps%lmax == -1 .and. ps%llocal == -1)
    
    ! We allocate all the stuff
    SAFE_ALLOCATE(ps%kb   (0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb  (0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%h    (0:ps%lmax, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(ps%density(1:ps%ispin))
    SAFE_ALLOCATE(ps%density_der(1:ps%ispin))

    nullify(ps%k)

    call spline_init(ps%kb)
    call spline_init(ps%dkb)
    call spline_init(ps%vl)
    call spline_init(ps%core)
    call spline_init(ps%density)
    call spline_init(ps%density_der)
    
    ! Now we load the necessary information.
    select case(ps%flavour)
    case(PS_TYPE_PSF)
      call ps_grid_load(ps, ps_psf%ps_grid)
      call ps_psf_end(ps_psf)
    case(PS_TYPE_CPI)
      call ps_grid_load(ps, ps_cpi%ps_grid)
      call ps_cpi_end(ps_cpi)
    case(PS_TYPE_FHI)
      call ps_grid_load(ps, ps_fhi%ps_grid)
      call ps_fhi_end(ps_fhi)
    case(PS_TYPE_HGH)
      SAFE_ALLOCATE(ps%k    (0:ps%lmax, 1:ps%kbc, 1:ps%kbc))
      call hgh_load(ps, ps_hgh)
      call hgh_end(ps_hgh)
    case(PS_TYPE_XML, PS_TYPE_UPF)
      if(ps_xml%initialized) then
        call ps_xml_load(ps, ps_xml)
        call ps_xml_end(ps_xml)
      else
        call ps_upf_load(ps, ps_upf)
        call ps_upf_end(ps_upf)
      end if
    end select

    if(ps_has_density(ps)) then 
      do is = 1, ps%ispin
        call spline_der(ps%density(is), ps%density_der(is))
      end do
    end if

    ps%has_long_range = .true.

    ps%is_separated = .false.

    POP_SUB(ps_init)
  end subroutine ps_init

  ! ---------------------------------------------------------
  !> separate the local potential into (soft) long-ranged and (hard) short-ranged parts
  subroutine ps_separate(ps)
    type(ps_t),        intent(inout) :: ps

    FLOAT, allocatable :: vsr(:), vlr(:), nlr(:)
    FLOAT :: r
    integer :: ii
    
    PUSH_SUB(ps_separate)

    ASSERT(ps%g%nrval > 0)

    SAFE_ALLOCATE(vsr(1:ps%g%nrval))
    SAFE_ALLOCATE(vlr(1:ps%g%nrval))
    SAFE_ALLOCATE(nlr(1:ps%g%nrval))

    ps%sigma_erf = CNST(0.625) ! This is hard-coded to a reasonable value
    
    vlr(1) = -ps%z_val*M_TWO/(sqrt(M_TWO*M_PI)*ps%sigma_erf)

    do ii = 1, ps%g%nrval
      r = ps%g%rofi(ii)
      if (ii > 1) then
        vlr(ii)  = -ps%z_val*loct_erf(r/(ps%sigma_erf*sqrt(M_TWO)))/r
      end if
      vsr(ii) = spline_eval(ps%vl, r) - vlr(ii)
      nlr(ii) = -ps%z_val*M_ONE/(ps%sigma_erf*sqrt(M_TWO*M_PI))**3*exp(-M_HALF*r**2/ps%sigma_erf**2)
    end do
    
    call spline_init(ps%vlr)
    call spline_fit(ps%g%nrval, ps%g%rofi, vlr, ps%vlr)

    call spline_init(ps%vlr_sq)
    call spline_fit(ps%g%nrval, ps%g%r2ofi, vlr, ps%vlr_sq)
    
    call spline_init(ps%nlr)
    call spline_fit(ps%g%nrval, ps%g%rofi, nlr, ps%nlr)
    
    !overwrite vl
    call spline_end(ps%vl)
    call spline_init(ps%vl)
    call spline_fit(ps%g%nrval, ps%g%rofi, vsr, ps%vl)
    
    SAFE_DEALLOCATE_A(vsr)
    SAFE_DEALLOCATE_A(vlr)
    SAFE_DEALLOCATE_A(nlr)
    
    ps%is_separated = .true.

    POP_SUB(ps_separate)
  end subroutine ps_separate
  
  
  ! ---------------------------------------------------------
  subroutine ps_getradius(ps)
    type(ps_t), intent(inout) :: ps
    integer :: l, j

    PUSH_SUB(ps_getradius)

    ps%rc_max = CNST(0.0)

    do l = 0, ps%lmax
      do j = 1, ps%kbc
        ps%rc_max = max(ps%rc_max, spline_cutoff_radius(ps%kb(l, j), ps%projectors_sphere_threshold))
      end do
    end do
    
    POP_SUB(ps_getradius)
  end subroutine ps_getradius


  ! ---------------------------------------------------------
  subroutine ps_derivatives(ps)
    type(ps_t), intent(inout) :: ps
    integer :: l, j

    PUSH_SUB(ps_derivatives)

    do l = 0, ps%lmax
      do j = 1, ps%kbc
        call spline_der(ps%kb(l, j), ps%dkb(l, j))
      end do
    end do

    POP_SUB(ps_derivatives)
  end subroutine ps_derivatives


  ! ---------------------------------------------------------
  subroutine ps_filter(ps, filter, gmax)
    type(ps_t), intent(inout) :: ps
    integer,    intent(in)    :: filter
    FLOAT,      intent(in)    :: gmax

    integer :: l, k, ispin
    type(profile_t), save:: prof

    FLOAT :: alpha, beta_fs, rmax, rcut, gamma, beta_rs

    PUSH_SUB(ps_filter)
    call profiling_in(prof, "PS_FILTER")

    select case(filter)
    case(PS_FILTER_NONE)

    case(PS_FILTER_TS)
      alpha = CNST(1.1)
      gamma = CNST(2.0)

      rmax = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold)
      call spline_filter_mask(ps%vl, max(0, ps%llocal), rmax, gmax, alpha, gamma)
      do l = 0, ps%lmax
        if(l == ps%llocal) cycle
        do k = 1, ps%kbc
          call spline_filter_mask(ps%kb(l, k), l, ps%rc_max, gmax, alpha, gamma)
        end do
      end do
      
      if(ps%nlcc) then
        rmax = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold)
        call spline_filter_mask(ps%core, 0, rmax, gmax, alpha, gamma)
      end if

      if(ps_has_density(ps)) then
        do ispin = 1, ps%ispin
          if(abs(spline_integral(ps%density(ispin))) > CNST(1.0e-12)) then
            rmax = spline_cutoff_radius(ps%density(ispin), ps%projectors_sphere_threshold)
            call spline_filter_mask(ps%density(ispin), 0, rmax, gmax, alpha, gamma)
          end if
        end do
      end if

    case(PS_FILTER_BSB)
      alpha   = CNST(0.7) ! The original was M_FOUR/CNST(7.0)
      beta_fs = CNST(18.0)
      rcut    = CNST(2.5)
      beta_rs = CNST(0.4)

      call spline_filter_bessel(ps%vl, ps%llocal, gmax, alpha, beta_fs, rcut, beta_rs)
      do l = 0, ps%lmax
        if(l == ps%llocal) cycle
        do k = 1, ps%kbc
          call spline_filter_bessel(ps%kb(l, k), l, gmax, alpha, beta_fs, rcut, beta_rs)
        end do
      end do
      
      if(ps%nlcc) then
        call spline_filter_bessel(ps%core, 0, gmax, alpha, beta_fs, rcut, beta_rs)
      end if
      
      if(ps_has_density(ps)) then
        do ispin = 1, ps%ispin
          call spline_filter_bessel(ps%density(ispin), 0, gmax, alpha, beta_fs, rcut, beta_rs)
        end do
      end if

    end select

    call profiling_out(prof)
    POP_SUB(ps_filter)
  end subroutine ps_filter


  ! ---------------------------------------------------------
  subroutine ps_debug(ps, dir)
    type(ps_t), intent(in) :: ps
    character(len=*), intent(in) :: dir

    ! We will plot also some Fourier transforms.
    type(spline_t), allocatable :: fw(:, :)
    FLOAT, parameter :: gmax = CNST(40.0)

    integer  :: iunit
    integer  :: j, k, l

    PUSH_SUB(ps_debug)

    ! A text file with some basic data.
    iunit = io_open(trim(dir)//'/pseudo-info', action='write')
    write(iunit,'(a,/)')      ps%label
    write(iunit,'(a,a,/)')    'Flavour : ', ps_name(ps%flavour)
    write(iunit,'(a,f6.3)')   'z       : ', ps%z
    write(iunit,'(a,f6.3,/)') 'zval    : ', ps%z_val
    write(iunit,'(a,i4)')     'lmax    : ', ps%lmax
    write(iunit,'(a,i4)')     'lloc    : ', ps%llocal
    write(iunit,'(a,i4,/)')   'kbc     : ', ps%kbc
    write(iunit,'(a,f9.5,/)') 'rcmax   : ', ps%rc_max
    write(iunit,'(a,/)')    'h matrix:'
    do l = 0, ps%lmax
      do k = 1, ps%kbc
        write(iunit,'(10f9.5)') (ps%h(l, k, j), j = 1, ps%kbc)
      end do
    end do
    if(associated(ps%k)) then
      write(iunit,'(/,a,/)')    'k matrix:'
      do l = 0, ps%lmax
        do k = 1, ps%kbc
          write(iunit,'(10f9.5)') (ps%k(l, k, j), j = 1, ps%kbc)
        end do
      end do
    end if
    call io_close(iunit)

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local', action='write')
    call spline_print(ps%vl, iunit)
    call io_close(iunit)

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local_long_range', action='write')
    call spline_print(ps%vlr, iunit)
    call io_close(iunit)
    
    ! Fourier transform of the local part
    iunit = io_open(trim(dir)//'/local_ft', action='write')
    SAFE_ALLOCATE(fw(1:1, 1:1))
    call spline_init(fw(1, 1))
    call spline_3dft(ps%vl, fw(1, 1), gmax = gmax)
    call spline_print(fw(1, 1), iunit)
    call spline_end(fw(1, 1))
    SAFE_DEALLOCATE_A(fw)
    call io_close(iunit)

    ! Kleinman-Bylander projectors
    iunit = io_open(trim(dir)//'/nonlocal', action='write')
    call spline_print(ps%kb, iunit)
    call io_close(iunit)

    iunit = io_open(trim(dir)//'/nonlocal_derivative', action='write')
    call spline_print(ps%dkb, iunit)
    call io_close(iunit)

    iunit = io_open(trim(dir)//'/nonlocal_ft', action='write')
    SAFE_ALLOCATE(fw(0:ps%lmax, 1:ps%kbc))
    call spline_init(fw)
    do k = 0, ps%lmax
      do j = 1, ps%kbc
        call spline_3dft(ps%kb(k, j), fw(k, j), gmax = gmax)
      end do
    end do
    call spline_print(fw, iunit)
    call spline_end(fw)
    SAFE_DEALLOCATE_A(fw)
    call io_close(iunit)

    ! Pseudo-wavefunctions
    iunit = io_open(trim(dir)//'/wavefunctions', action='write')
    call spline_print(ps%ur, iunit)
    call io_close(iunit)

    ! Density
    if (ps%has_density) then
      iunit = io_open(trim(dir)//'/density', action='write')
      call spline_print(ps%density, iunit)
      call io_close(iunit)

      iunit = io_open(trim(dir)//'/density_derivative', action='write')
      call spline_print(ps%density_der, iunit)
      call io_close(iunit)
    end if

    ! Non-linear core-corrections
    if(ps%nlcc) then
      iunit = io_open(trim(dir)//'/nlcc', action='write')
      call spline_print(ps%core, iunit)
      call io_close(iunit)
    end if

    POP_SUB(ps_debug)
  end subroutine ps_debug


  ! ---------------------------------------------------------
  subroutine ps_end(ps)
    type(ps_t), intent(inout) :: ps

    if(.not. associated(ps%kb)) return

    PUSH_SUB(ps_end)

    if(ps%is_separated) then
      call spline_end(ps%vlr)
      call spline_end(ps%vlr_sq)
      call spline_end(ps%nlr)
    end if

    call spline_end(ps%kb)
    call spline_end(ps%dkb)
    call spline_end(ps%ur)
    call spline_end(ps%ur_sq)

    call spline_end(ps%vl)
    call spline_end(ps%core)

    if(associated(ps%density)) call spline_end(ps%density)
    if(associated(ps%density_der)) call spline_end(ps%density_der)

    call logrid_end(ps%g)

    SAFE_DEALLOCATE_P(ps%kb)
    SAFE_DEALLOCATE_P(ps%dkb)
    SAFE_DEALLOCATE_P(ps%ur)
    SAFE_DEALLOCATE_P(ps%ur_sq)
    SAFE_DEALLOCATE_P(ps%h)
    SAFE_DEALLOCATE_P(ps%k)
    SAFE_DEALLOCATE_P(ps%density)
    SAFE_DEALLOCATE_P(ps%density_der)

    POP_SUB(ps_end)
  end subroutine ps_end


  ! ---------------------------------------------------------
  subroutine hgh_load(ps, ps_hgh)
    type(ps_t),     intent(inout) :: ps
    type(ps_hgh_t), intent(inout) :: ps_hgh

    integer :: l, ll
    FLOAT :: x

    PUSH_SUB(hgh_load)

    ! Fixes some components of ps
    ps%z_val = ps_hgh%z_val
    ps%nlcc = .false.
    if(ps%lmax>=0) then
      ps%rc_max = CNST(1.1) * maxval(ps_hgh%kbr(0:ps%lmax)) ! Increase a little.
    else
      ps%rc_max = M_ZERO
    end if
    ps%h(0:ps%lmax, 1:ps%kbc, 1:ps%kbc) = ps_hgh%h(0:ps%lmax, 1:ps%kbc, 1:ps%kbc)
    ps%k(0:ps%lmax, 1:ps%kbc, 1:ps%kbc) = ps_hgh%k(0:ps%lmax, 1:ps%kbc, 1:ps%kbc)

    ! Fixes the occupations
    if(ps%ispin == 2) then
      do l = 1, ps%conf%p
        ll = ps%conf%l(l)
        x = ps%conf%occ(l, 1)
        ps%conf%occ(l, 1) = min(x, real(2*ll+1, REAL_PRECISION))
        ps%conf%occ(l, 2) = x - ps%conf%occ(l, 1)
      end do
    end if

    ! now we fit the splines
    call get_splines()

    POP_SUB(hgh_load)

  contains

    ! ---------------------------------------------------------
    subroutine get_splines()
      integer :: l, is, nrc, j, ip
      FLOAT, allocatable :: hato(:), dens(:)

      PUSH_SUB(hgh_load.get_splines)

      SAFE_ALLOCATE(hato(1:ps_hgh%g%nrval))
      SAFE_ALLOCATE(dens(1:ps_hgh%g%nrval))

      ! Interpolate the KB-projection functions
      do l = 0, ps_hgh%l_max
        do j = 1, 3
          hato = M_ZERO
          nrc = nint(log(ps_hgh%kbr(l)/ps_hgh%g%b + M_ONE)/ps_hgh%g%a) + 1
          hato(1:nrc) = ps_hgh%kb(1:nrc, l, j)
          call spline_fit(ps_hgh%g%nrval, ps_hgh%g%rofi, hato, ps%kb(l, j))
        end do
      end do

      ! Now the part corresponding to the local pseudopotential
      ! where the asymptotic part is subtracted
      call spline_fit(ps_hgh%g%nrval, ps_hgh%g%rofi, ps_hgh%vlocal, ps%vl)

      ! Define the table for the pseudo-wavefunction components (using splines)
      ! with a correct normalization function
      do is = 1, ps%ispin
        dens = CNST(0.0)
        do l = 1, ps%conf%p
          hato(2:ps_hgh%g%nrval) = ps_hgh%rphi(2:ps_hgh%g%nrval, l)/ps_hgh%g%rofi(2:ps_hgh%g%nrval)
          hato(1) = hato(2)

          forall(ip = 1:ps_hgh%g%nrval) dens(ip) = dens(ip) + ps%conf%occ(l, is)*hato(ip)**2/(M_FOUR*M_PI)
          
          call spline_fit(ps_hgh%g%nrval, ps_hgh%g%rofi, hato, ps%ur(l, is))
          call spline_fit(ps_hgh%g%nrval, ps_hgh%g%r2ofi, hato, ps%ur_sq(l, is))
        end do
        call spline_fit(ps_hgh%g%nrval, ps_hgh%g%rofi, dens, ps%density(is))
      end do

      SAFE_DEALLOCATE_A(hato)
      SAFE_DEALLOCATE_A(dens)

      POP_SUB(hgh_load.get_splines)
    end subroutine get_splines
  end subroutine hgh_load


  ! ---------------------------------------------------------
  subroutine ps_grid_load(ps, ps_grid)
    type(ps_t),         intent(inout) :: ps
    type(ps_in_grid_t), intent(in)  :: ps_grid

    PUSH_SUB(ps_grid_load)

    ! Fixes some components of ps, read in ps_grid
    ps%z_val = ps_grid%zval

    ps%nlcc = ps_grid%core_corrections

    ps%h(0:ps%lmax, 1, 1) = ps_grid%dkbcos(1:ps%lmax+1)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = maxval(ps_grid%kb_radius(1:ps%lmax+1)) * CNST(1.5)

    ! now we fit the splines
    call get_splines(ps_grid%g)

    ! Passes from Rydbergs to Hartrees.
    ps%h(0:ps%lmax,:,:)    = ps%h(0:ps%lmax,:,:)    / M_TWO

    POP_SUB(ps_grid_load)

  contains

    subroutine get_splines(g)
      type(logrid_t), intent(in) :: g

      FLOAT, allocatable :: hato(:), dens(:)
      integer :: is, l, ir, nrc, ip

      PUSH_SUB(ps_grid_load.get_splines)

      SAFE_ALLOCATE(hato(1:g%nrval))
      SAFE_ALLOCATE(dens(1:g%nrval))

      ! the wavefunctions
      do is = 1, ps%ispin

        dens = CNST(0.0)

        do l = 1, ps_grid%no_l_channels
          hato(2:) = ps_grid%rphi(2:, l, 1+is)/g%rofi(2:)
          hato(1)  = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
            hato(2), hato(3))

          forall(ip = 1:g%nrval) dens(ip) = dens(ip) + ps%conf%occ(l, is)*hato(ip)**2/(M_FOUR*M_PI)
          
          call spline_fit(g%nrval, g%rofi, hato, ps%ur(l, is))
          call spline_fit(g%nrval, g%r2ofi, hato, ps%ur_sq(l, is))

        end do

        call spline_fit(g%nrval, g%rofi, dens, ps%density(is))
      end do
      

      ! the Kleinman-Bylander projectors
      do l = 1, ps%lmax+1
        nrc = logrid_index(g, ps_grid%kb_radius(l)) + 1
        hato(1:nrc)         = ps_grid%KB(1:nrc, l)
        hato(nrc+1:g%nrval) = M_ZERO

        call spline_fit(g%nrval, g%rofi, hato, ps%kb(l-1, 1))
      end do

      ! Now the part corresponding to the local pseudopotential
      ! where the asymptotic part is subtracted
      hato(:) = ps_grid%vlocal(:)/M_TWO
      call spline_fit(g%nrval, g%rofi, hato, ps%vl)
      
      if(ps_grid%core_corrections) then
        ! find cutoff radius
        hato(2:) = ps_grid%chcore(2:)/(M_FOUR*M_PI*g%rofi(2:)**2)

        do ir = g%nrval-1, 2, -1
          if(hato(ir) > eps) then
            nrc = ir + 1
            exit
          end if
        end do

        hato(nrc:g%nrval) = M_ZERO
        hato(1) = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
          hato(2), hato(3))

        call spline_fit(g%nrval, g%rofi, hato, ps%core)
      end if

      SAFE_DEALLOCATE_A(hato)
      SAFE_DEALLOCATE_A(dens)
      
      POP_SUB(ps_grid_load.get_splines)
    end subroutine get_splines
  end subroutine ps_grid_load


  ! ---------------------------------------------------------
  subroutine ps_upf_load(ps, ps_upf)
    type(ps_t),     intent(inout) :: ps
    type(ps_upf_t), intent(in)    :: ps_upf

    integer :: i, l, ll, is, nrc, ir, j, ij, ispin, ip
    FLOAT :: x
    FLOAT, allocatable :: hato(:), dens(:)

    PUSH_SUB(ps_upf_load)

    ! Fixes some components of ps, read in ps_upf
    ps%z_val = ps_upf%z_val

    ps%nlcc = ps_upf%nlcc

    ! if there are two projectors for l==0, this is a hamann
    ps%hamann = ps_upf%nchannels(0) == 2

    ! The spin-dependent pseudopotentials are not supported yet, so we need to fix the occupations
    ! if we want to have a spin-dependent atomic density.
    if(ps%ispin == 2) then
      do l = 1, ps%conf%p
        ll = ps%conf%l(l)
        x = ps%conf%occ(l, 1)
        ps%conf%occ(l, 1) = min(x, real(2*ll+1, REAL_PRECISION))
        ps%conf%occ(l, 2) = x - ps%conf%occ(l, 1)
      end do
    end if

    SAFE_ALLOCATE(hato(1:ps%g%nrval))

    ! only ps%g%rofi(1) is allowed to be zero
    if(any(abs(ps%g%rofi(2:ps%g%nrval)) < M_EPSILON)) then
      message(1) = "Illegal zero values in UPF radial grid ps%g%rofi(2:ps%g%nrval)"
      call messages_fatal(1)
    end if

    !Non-linear core-corrections
    if(ps_upf%nlcc) then
      ! find cutoff radius
      hato = ps_upf%core_density

      do ir = ps%g%nrval-1, 1, -1
        if(hato(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do

      hato(nrc:ps%g%nrval) = M_ZERO

      call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%core)
    end if

    ! Now the part corresponding to the local pseudopotential
    ! where the asymptotic part is subtracted
    hato(:) = ps_upf%v_local/M_TWO
    call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%vl)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = maxval(ps_upf%kb_radius)
    ps%rc_max = max(ps_upf%local_radius, ps%rc_max) * CNST(1.5)

    ! Interpolate the KB-projection functions
    if (ps_upf%l_loc >= 0) then
      hato = M_ZERO
      do j = 1, ps%kbc
        call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%l_loc, j))
      end do
    end if

    ps%local = ps_upf%n_proj == 0
    
    ps%h = M_ZERO
    do i = 1, ps_upf%n_proj

      if(associated(ps_upf%proj_j)) then
        
        ij = 1
        if (ps_upf%kb_nc == 2) then
          if (ps_upf%proj_j(i) == ps_upf%proj_l(i) - M_HALF) ij = 2
        end if

      else

        ij = 1 + count(ps_upf%proj_l(1:i - 1) == ps_upf%proj_l(i))
        
      end if

      ASSERT(ij <= ps%kbc)
      
      ps%h(ps_upf%proj_l(i), ij, ij) = ps_upf%e(i)

      if(.not. ps_upf%version2) then
        ! in UPF 1 this value is in Ry^-1
        ps%h(ps_upf%proj_l(i), ij, ij) = ps%h(ps_upf%proj_l(i), ij, ij)*M_TWO
      else
        ! in UPF 2 this value is in Ry
        ps%h(ps_upf%proj_l(i), ij, ij) = ps%h(ps_upf%proj_l(i), ij, ij)/M_TWO
      end if
      
      nrc = logrid_index(ps%g, ps_upf%kb_radius(i)) + 1
      hato(2:nrc) = ps_upf%proj(2:nrc, i)/ps%g%rofi(2:nrc) ! in upf the projector is given in Rydbergs and is multiplied by r
      hato(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), ps%g%rofi(3), hato(2), hato(3)) !take care of the point at zero
      hato(nrc+1:ps%g%nrval) = M_ZERO

      if(.not. ps_upf%version2) then
        ! in UPF 1 the projectors are in Ry.
        hato(1:nrc) = hato(1:nrc)/M_TWO
        ! in v2 they are in Bohr^{-1/2}, so no conversion is required
      end if

      call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), ij))

      if(.not. ps%hamann) then
        if (ps_upf%proj_l(i) == 0 .and. ps_upf%kb_nc == 2) then
          hato = M_ZERO
          call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), 2))
        end if
      end if

    end do

    if(ps%conf%p > 0) then
      
      ! Define the table for the pseudo-wavefunction components (using splines)
      ! with a correct normalization function
      do is = 1, ps%ispin
        do l = 1, ps%conf%p
          ! do not divide by zero
          if(ps%g%rofi(1) > M_EPSILON) then
            hato(1) = ps_upf%wfs(1, l)/ps%g%rofi(1)
          else
            hato(1) = M_ZERO
          end if
          ! rofi /= 0 except rofi(1) possibly
          hato(2:ps%g%nrval) = ps_upf%wfs(2:ps%g%nrval, l)/ps%g%rofi(2:ps%g%nrval)
          hato(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), ps%g%rofi(3), hato(2), hato(3)) !take care of the point at zero
          
          call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%ur(l, is))
          call spline_fit(ps%g%nrval, ps%g%r2ofi, hato, ps%ur_sq(l, is))
        end do
      end do

    end if

    SAFE_DEALLOCATE_A(hato)

    
    SAFE_ALLOCATE(dens(1:ps%g%nrval))
    
    dens(2:ps%g%nrval) = ps_upf%rho(2:ps%g%nrval)/ps%g%r2ofi(2:ps%g%nrval)/ps%ispin/CNST(4.0)/M_PI
    dens(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), ps%g%rofi(3), dens(2), dens(3)) !take care of the point at zero
      
    do is = 1, ps%ispin
      call spline_fit(ps%g%nrval, ps%g%rofi, dens, ps%density(is))
    end do

    SAFE_DEALLOCATE_A(dens)

    POP_SUB(ps_upf_load)
  end subroutine ps_upf_load

  
  ! ---------------------------------------------------------
  subroutine ps_xml_load(ps, ps_xml)
    type(ps_t),     intent(inout) :: ps
    type(ps_xml_t), intent(in)    :: ps_xml

    integer :: ll, ip, is, ic, jc, ir, nrc, ii
    FLOAT :: rr, kbcos, kbnorm, dnrm, avgv, volume_element
    FLOAT, allocatable :: vlocal(:), kbprojector(:), wavefunction(:), nlcc_density(:), dens(:)

    PUSH_SUB(ps_xml_load)

    if(ps_xml%kleinman_bylander .and. ps_xml%nchannels == 2) then
      ps%hamann = .true.
    end if
    
    ps%nlcc = ps_xml%nlcc

    ps%z_val = ps_xml%valence_charge

    ! the local potential
    SAFE_ALLOCATE(vlocal(1:ps%g%nrval))

    do ip = 1, ps%g%nrval
      rr = (ip - 1)*ps_xml%mesh_spacing
      if(ip <= ps_xml%grid_size) then
        vlocal(ip) = ps_xml%potential(ip, ps%llocal)
      else
        vlocal(ip) = -ps_xml%valence_charge/rr
      end if
    end do

    call spline_fit(ps%g%nrval, ps%g%rofi, vlocal, ps%vl)

    SAFE_DEALLOCATE_A(vlocal)

    SAFE_ALLOCATE(kbprojector(1:ps%g%nrval))
    SAFE_ALLOCATE(wavefunction(1:ps%g%nrval))

    kbprojector = CNST(0.0)
    wavefunction = CNST(0.0)

    ! the projectors and the orbitals
    if(ps_xml%kleinman_bylander) then

      do ll = 0, ps_xml%lmax
        do ic = 1, ps_xml%nchannels
          
          do ip = 1, ps%g%nrval
            if(ip <= ps_xml%grid_size) then
              kbprojector(ip) = ps_xml%projector(ip, ll, ic)
            else
              kbprojector(ip) = 0.0
            end if
          end do

          call spline_fit(ps%g%nrval, ps%g%rofi, kbprojector, ps%kb(ll, ic))

          do jc = 1, ps_xml%nchannels
            ps%h(ll, ic, jc) = ps_xml%dij(ll, ic, jc)
          end do

        end do
      end do

      do ii = 1, ps_xml%nwavefunctions
        
        do ip = 1, ps%g%nrval
          if(ip <= ps_xml%grid_size) then
            wavefunction(ip) = ps_xml%wavefunction(ip, ii)
          else
            wavefunction(ip) = CNST(0.0)
          end if
        end do
      
        do is = 1, ps%ispin
          call spline_fit(ps%g%nrval, ps%g%rofi, wavefunction, ps%ur(ii, is))
          call spline_fit(ps%g%nrval, ps%g%r2ofi, wavefunction, ps%ur_sq(ii, is))
        end do

      end do
      
    else

      do ll = 0, ps_xml%lmax
        ! we need to build the KB projectors
        ! the procedure was copied from ps_in_grid.F90 (r12967)
        dnrm = M_ZERO
        avgv = M_ZERO
        do ip = 1, ps_xml%grid_size
          rr = (ip - 1)*ps_xml%mesh_spacing
          volume_element = rr**2*ps_xml%mesh_spacing
          kbprojector(ip) = (ps_xml%potential(ip, ll) - ps_xml%potential(ip, ps%llocal))*ps_xml%wavefunction(ip, ll)
          dnrm = dnrm + kbprojector(ip)**2*volume_element
          avgv = avgv + kbprojector(ip)*ps_xml%wavefunction(ip, ll)*volume_element
        end do

        kbcos = dnrm/(avgv + CNST(1.0e-20))
        kbnorm = M_ONE/(sqrt(dnrm) + CNST(1.0e-20))

        if(ll /= ps%llocal) then
          ps%h(ll, 1, 1) = kbcos        
          kbprojector = kbprojector*kbnorm
        else
          ps%h(ll, 1, 1) = CNST(0.0)
        end if

        call spline_fit(ps%g%nrval, ps%g%rofi, kbprojector, ps%kb(ll, 1))

        ! wavefunctions, for the moment we pad them with zero
        do ip = 1, ps%g%nrval
          rr = (ip - 1)*ps_xml%mesh_spacing
          if(ip <= ps_xml%grid_size) then
            wavefunction(ip) = ps_xml%wavefunction(ip, ll)
          else
            wavefunction(ip) = CNST(0.0)
          end if
        end do

        do is = 1, ps%ispin
          call spline_fit(ps%g%nrval, ps%g%rofi, wavefunction, ps%ur(ll + 1, is))
          call spline_fit(ps%g%nrval, ps%g%r2ofi, wavefunction, ps%ur_sq(ll + 1, is))
        end do
      end do
        
    end if

    ps%has_density = ps_xml%has_density
    
    if(ps_has_density(ps)) then
      
      SAFE_ALLOCATE(dens(1:ps%g%nrval))
      
      dens(1:ps_xml%grid_size) = ps_xml%density(1:ps_xml%grid_size)/ps%ispin
      dens(ps_xml%grid_size + 1:ps%g%nrval) = CNST(0.0)
      
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%rofi, dens, ps%density(is))
      end do
      
      SAFE_DEALLOCATE_A(dens)
    end if
    
    !Non-linear core-corrections
    if(ps_xml%nlcc) then

      SAFE_ALLOCATE(nlcc_density(1:ps%g%nrval))

      nlcc_density(1:ps_xml%grid_size) = ps_xml%nlcc_density(1:ps_xml%grid_size)
      
      ! find cutoff radius
      do ir = ps_xml%grid_size - 1, 1, -1
        if(nlcc_density(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do

      nlcc_density(nrc:ps%g%nrval) = M_ZERO
      
      call spline_fit(ps%g%nrval, ps%g%rofi, nlcc_density, ps%core)

      SAFE_DEALLOCATE_A(nlcc_density)
    end if
    
    call ps_getradius(ps)

    SAFE_DEALLOCATE_A(kbprojector)
    SAFE_DEALLOCATE_A(wavefunction)

    POP_SUB(ps_xml_load)
  end subroutine ps_xml_load


  ! ---------------------------------------------------------
  !> Returns the number of atomic orbitals that can be used for LCAO calculations.
  pure integer function ps_niwfs(ps)
    type(ps_t), intent(in) :: ps

    integer :: i, l

    ps_niwfs = 0
    do i = 1, ps%conf%p
      l = ps%conf%l(i)
      ps_niwfs = ps_niwfs + (2*l+1)
    end do

  end function ps_niwfs


  ! ---------------------------------------------------------
  integer function ps_get_type(filename) result(type)
    character(len=*), intent(in) :: filename

    PUSH_SUB(ps_get_type)

    type = 0
    
    if(index(filename, ".psf ") /= 0) type = PS_TYPE_PSF
    if(index(filename, ".PSF ") /= 0) type = PS_TYPE_PSF
    if(index(filename, ".hgh ") /= 0) type = PS_TYPE_HGH
    if(index(filename, ".HGH ") /= 0) type = PS_TYPE_HGH
    if(index(filename, ".cpi ") /= 0) type = PS_TYPE_CPI
    if(index(filename, ".CPI ") /= 0) type = PS_TYPE_CPI
    if(index(filename, ".fhi ") /= 0) type = PS_TYPE_FHI
    if(index(filename, ".FHI ") /= 0) type = PS_TYPE_FHI
    if(index(filename, ".upf ") /= 0) type = PS_TYPE_UPF
    if(index(filename, ".UPF ") /= 0) type = PS_TYPE_UPF
    if(index(filename, ".xml ") /= 0) type = PS_TYPE_XML
    if(index(filename, ".XML ") /= 0) type = PS_TYPE_XML
    if(index(filename, ".psml") /= 0) type = PS_TYPE_XML
    if(index(filename, ".PSML") /= 0) type = PS_TYPE_XML
    
    POP_SUB(ps_get_type)    
  end function ps_get_type


  !---------------------------------------
  pure integer function ps_type(ps)
    type(ps_t), intent(in) :: ps

    ps_type = ps%flavour
  end function ps_type

  
  !---------------------------------------
  pure logical function ps_has_density(ps) result(has_density)
    type(ps_t), intent(in) :: ps

    has_density = ps%has_density

  end function ps_has_density

  
  !---------------------------------------
  FLOAT function ps_density_volume(ps) result(volume)
    type(ps_t), intent(in) :: ps

    integer :: ip, ispin
    FLOAT :: rr
    FLOAT, allocatable ::vol(:)
    type(spline_t) :: volspl
    
    PUSH_SUB(ps_density_volume)

    SAFE_ALLOCATE(vol(1:ps%g%nrval))
    
    do ip = 1, ps%g%nrval
      rr = ps%g%rofi(ip)
      vol(ip) = CNST(0.0)
      do ispin = 1, ps%ispin
        vol(ip) = vol(ip) + spline_eval(ps%density(ispin), rr)*CNST(4.0)*M_PI*rr**5
      end do
    end do

    call spline_init(volspl)
    call spline_fit(ps%g%nrval, ps%g%rofi, vol, volspl)
    volume = spline_integral(volspl)
    call spline_end(volspl)

    SAFE_DEALLOCATE_A(vol)
    
    POP_SUB(ps_density_volume)
  end function ps_density_volume
  
#include "ps_pspio_inc.F90"
 
end module ps_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
