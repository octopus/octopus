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

module ps_m
  use atomic_m
  use datasets_m
  use global_m
  use io_m
  use loct_math_m
  use parser_m
  use logrid_m
  use messages_m
  use profiling_m
  use ps_cpi_m
  use ps_fhi_m
  use ps_hgh_m
  use ps_in_grid_m
  use ps_psf_m
  use ps_upf_m
  use splines_m
  use spline_filter_m
  implicit none

  private
  public ::                     &
    ps_t,                       &
    ps_init,                    &
    ps_separate,                &
    ps_filter,                  &
    ps_getradius,               &
    ps_derivatives,             &
    ps_debug,                   &
    ps_niwfs,                   &
    ps_end

  integer, parameter, public :: &
    PS_TYPE_PSF = 100,          &
    PS_TYPE_HGH = 101,          &
    PS_TYPE_CPI = 102,          &
    PS_TYPE_FHI = 103,          &
    PS_TYPE_UPF = 104

  integer, parameter, public :: &
    PS_FILTER_NONE = 1,         &
    PS_FILTER_TS   = 2,         &
    PS_FILTER_BSB  = 3

  character(len=3), parameter  :: ps_name(PS_TYPE_PSF:PS_TYPE_UPF) = (/"tm2", "hgh", "cpi", "fhi", "upf"/)

  type ps_t
    character(len=10) :: label
    integer           :: flavour

    integer  :: ispin    !< Consider spin (ispin = 2) or not (ispin = 1)
    FLOAT    :: z, z_val
    type(valconf_t)   :: conf
    type(logrid_t) :: g
    type(spline_t), pointer :: ur(:, :)     !< atomic wavefunctions
    type(spline_t), pointer :: ur_sq(:, :)  !< atomic wavefunctions

    ! Kleinman-Bylander projectors stuff
    integer  :: l_max    !< maximum value of l to take
    integer  :: l_loc    !< which component to take as local

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


    character(len=4) :: icore !< nonlinear core corrections
    type(spline_t) :: core !< normalization \int dr 4 pi r^2 rho(r) = N


    !LONG-RANGE PART OF THE LOCAL POTENTIAL
    
    logical :: has_long_range

    type(spline_t) :: vlr         !< the long-range part of the local potential
    type(spline_t) :: vlr_sq      !< the long-range part of the
                                  !< local potential in terms of r^2, to avoid the sqrt
    type(spline_t) :: nlr         !< the charge density associated with the long-range part
    
    FLOAT :: sigma_erf            !< the a constant in erf(r/(sqrt(2)*sigma))/r
    FLOAT :: a_erf                !< the a constant in erf(ar)/r

    type(spline_t) :: vion        !< the potential that other ions see
    type(spline_t) :: dvion       !< the potential that other ions see
    
    logical :: is_separated
  end type ps_t

  FLOAT, parameter :: eps = CNST(1.0e-8)

contains


  ! ---------------------------------------------------------
  subroutine ps_init(ps, label, flavour, z, lmax, lloc, ispin)
    type(ps_t),        intent(out)   :: ps
    character(len=10), intent(in)    :: label
    integer,           intent(in)    :: flavour
    integer,           intent(inout) :: lmax
    integer,           intent(in)    :: lloc, ispin
    FLOAT,             intent(in)    :: z

    type(ps_psf_t) :: ps_psf !< SIESTA pseudopotential
    type(ps_cpi_t) :: ps_cpi !< Fritz-Haber pseudopotential
    type(ps_fhi_t) :: ps_fhi !< Fritz-Haber pseudopotential (from abinit)
    type(ps_upf_t) :: ps_upf !< In case UPF format is used
    type(hgh_t)    :: psp    !< In case Hartwigsen-Goedecker-Hutter ps are used.

    PUSH_SUB(ps_init)

    ! Sets the flavour, label, and number of spin channels.
    ps%flavour = flavour
    ps%label   = label
    ps%ispin   = ispin
    ! Initialization and processing.
    ASSERT(flavour >= PS_TYPE_PSF .and. flavour <= PS_TYPE_UPF)

    select case(flavour)
    case(PS_TYPE_PSF)
      call ps_psf_init(ps_psf, trim(label), ispin)

      call valconf_copy(ps%conf, ps_psf%conf)
      ps%z      = z
      ps%conf%z = z     ! atomic number
      ps%kbc    = 1     ! only one projector per angular momentum
      ps%l_loc  = lloc  ! the local part of the pseudo

      ps%l_max  = min(ps_psf%ps_grid%no_l_channels - 1, lmax) ! Maybe the file has not enough components.
      ps%conf%p = ps_psf%ps_grid%no_l_channels
      if(ps%l_max == 0) ps%l_loc = 0 ! Vanderbilt is not acceptable if ps%l_max == 0.

      call ps_psf_process(ps_psf, lmax, ps%l_loc)
      call logrid_copy(ps_psf%ps_grid%g, ps%g)

    case(PS_TYPE_CPI)
      call ps_cpi_init(ps_cpi, trim(label))

      call valconf_null(ps%conf)
      ps%conf%z      = z
      ps%conf%symbol = label(1:2)
      ps%conf%type   = 1
      ps%conf%p      = ps_cpi%ps_grid%no_l_channels

      ps%z      = z
      ps%kbc    = 1     ! only one projector per angular momentum
      ps%l_loc  = lloc  ! the local part of the pseudo

      ps%l_max  = min(ps%conf%p - 1, lmax)   ! Maybe the file has not enough components.
      if(ps%l_max == 0) ps%l_loc = 0 ! Vanderbilt is not acceptable if ps%l_max == 0.

      call ps_cpi_process(ps_cpi, ps%l_loc)
      call logrid_copy(ps_cpi%ps_grid%g, ps%g)

    case(PS_TYPE_FHI)
      call ps_fhi_init(ps_fhi, trim(label))

      call valconf_null(ps%conf)
      ps%conf%z      = z
      ps%conf%symbol = label(1:2)
      ps%conf%type   = 1
      ps%conf%p      = ps_fhi%ps_grid%no_l_channels

      ps%z      = z
      ps%kbc    = 1     ! only one projector per angular momentum
      ps%l_loc  = lloc  ! the local part of the pseudo

      ps%l_max  = min(ps%conf%p - 1, lmax)   ! Maybe the file has not enough components.
      if(ps%l_max == 0) ps%l_loc = 0 ! Vanderbilt is not acceptable if ps%l_max == 0.

      call ps_fhi_process(ps_fhi, lmax, ps%l_loc)
      call logrid_copy(ps_fhi%ps_grid%g, ps%g)

    case(PS_TYPE_HGH)
      call hgh_init(psp, trim(label))
      call valconf_copy(ps%conf, psp%conf)

      ps%z        = z
      ps%kbc      = 3
      ps%l_loc    = -1
      ps%l_max    = psp%l_max

      call hgh_process(psp)
      call logrid_copy(psp%g, ps%g)

    case(PS_TYPE_UPF)
      call ps_upf_init(ps_upf, trim(label))

      call valconf_copy(ps%conf, ps_upf%conf)
      ps%z      = z
      ps%conf%z = z
      ps%kbc    = ps_upf%kb_nc
      lmax      = ps_upf%l_max
      ps%l_max  = ps_upf%l_max
      ps%l_loc  = ps_upf%l_local

      nullify(ps%g%drdi, ps%g%s)
      ps%g%nrval = ps_upf%np
      SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
      SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))
      ps%g%rofi = ps_upf%r
      ps%g%r2ofi = ps%g%rofi**2

    end select

    write(message(1), '(a,i2,a)') "Info: l = ", ps%l_max, " is maximum angular momentum considered."
    call messages_info(1)

    ! We allocate all the stuff
    SAFE_ALLOCATE(ps%kb   (0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb  (0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%h    (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(ps%k    (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    call spline_init(ps%kb)
    call spline_init(ps%dkb)
    call spline_init(ps%vl)
    call spline_init(ps%core)

    ! Now we load the necessary information.
    select case(flavour)
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
      call hgh_load(ps, psp)
      call hgh_end(psp)
    case(PS_TYPE_UPF)
      call ps_upf_load(ps, ps_upf)
      call ps_upf_end(ps_upf)
    end select

    ! Fix the threshold to calculate the radius of the projector-function localization spheres:

    call messages_obsolete_variable('SpecieProjectorSphereThreshold', 'SpeciesProjectorSphereThreshold')

    !%Variable SpeciesProjectorSphereThreshold
    !%Type float
    !%Default 0.001
    !%Section System::Species
    !%Description
    !% The pseudopotentials may be composed of a local part, and a linear combination of nonlocal
    !% operators. These nonlocal projectors have "projector" form, |<i>v</i>&gt;&lt;<i>v</i>| (or, more generally
    !% speaking, |<i>u</i>&gt;&lt;<i>v</i>|). These projectors are localized in real space -- that is, the function <i>v</i>
    !% has a finite support around the nucleus. This region where the projectors are localized should
    !% be small or else the computation time required to operate with them will be very large.
    !% 
    !% In practice, this localization is fixed by requiring the definition of the projectors to be
    !% contained in a sphere of a certain radius. This radius is computed by making sure that the 
    !% absolute value of the projector functions, at points outside the localization sphere, is 
    !% below a certain threshold. This threshold is set by <tt>SpeciesProjectorSphereThreshold</tt>.
    !%End
    call parse_float(datasets_check('SpeciesProjectorSphereThreshold'), &
      CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call input_error('SpeciesProjectorSphereThreshold')

    ps%has_long_range = .true.

    ps%is_separated = .false.

    POP_SUB(ps_init)

  end subroutine ps_init

  subroutine ps_separate(ps)
    type(ps_t),        intent(inout) :: ps

    FLOAT, allocatable :: vsr(:), vlr(:), nlr(:), vion(:)
    FLOAT :: r
    integer :: ii
    
    PUSH_SUB(ps_separate)

    !separate the local potential into (soft) long-ranged and (hard) short-ranged parts
    
    ps%sigma_erf = CNST(0.625) ! This is hard-coded to a reasonable value.
    ps%a_erf = M_ONE/(ps%sigma_erf*sqrt(M_TWO))
    
    ASSERT(ps%g%nrval > 0)

    SAFE_ALLOCATE( vsr(1:ps%g%nrval))
    SAFE_ALLOCATE( vlr(1:ps%g%nrval))
    SAFE_ALLOCATE( nlr(1:ps%g%nrval))
    SAFE_ALLOCATE(vion(1:ps%g%nrval))
    
    vlr(1) = -ps%z_val*M_TWO/(sqrt(M_TWO*M_PI)*ps%sigma_erf)

    do ii = 1, ps%g%nrval
      r = ps%g%rofi(ii)
      if ( ii > 1) then
        vlr(ii)  = -ps%z_val*loct_erf(r/(ps%sigma_erf*sqrt(M_TWO)))/r
        vion(ii) = -ps%z_val/r - vlr(ii)
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
    
    ! The ion-ion interaction
    vion(1) = vion(2)
    
    call spline_init(ps%vion)
    call spline_fit(ps%g%nrval, ps%g%rofi, vion, ps%vion)

    call spline_init(ps%dvion)
    call spline_der(ps%vion, ps%dvion)

    SAFE_DEALLOCATE_A(vsr)
    SAFE_DEALLOCATE_A(vlr)
    SAFE_DEALLOCATE_A(nlr)
    SAFE_DEALLOCATE_A(vion)
    
    ps%is_separated = .true.

    POP_SUB(ps_separate)
  end subroutine ps_separate
  
  ! ---------------------------------------------------------
  subroutine ps_getradius(ps)
    type(ps_t), intent(inout) :: ps
    integer :: l, j

    PUSH_SUB(ps_getradius)

    ps%rc_max = CNST(0.0)

    do l = 0, ps%l_max
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

    do l = 0, ps%l_max
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
    integer :: l, k

    FLOAT :: alpha, beta_fs, rmax, rcut, gamma, beta_rs

    PUSH_SUB(ps_filter)

    select case(filter)
    case(PS_FILTER_NONE)

    case(PS_FILTER_TS)
      alpha = CNST(1.1)
      gamma = CNST(2.0)


      rmax = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold)
      call spline_filter_mask(ps%vl, max(0, ps%l_loc), rmax, gmax, alpha, gamma)
      do l = 0, ps%l_max
        if(l == ps%l_loc) cycle
        do k = 1, ps%kbc
          call spline_filter_mask(ps%kb(l, k), l, ps%rc_max, gmax, alpha, gamma)
        end do
      end do
      
      if(trim(ps%icore).ne.'nc') then
        rmax = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold)
        call spline_filter_mask(ps%core, 0, rmax, gmax, alpha, gamma)
      end if

    case(PS_FILTER_BSB)
      alpha   = CNST(0.7) ! The original was M_FOUR/M_SEVEN
      beta_fs = CNST(18.0)
      rcut    = CNST(2.5)
      beta_rs = CNST(0.4)

      call spline_filter_bessel(ps%vl, ps%l_loc, gmax, alpha, beta_fs, rcut, beta_rs)
      do l = 0, ps%l_max
        if(l == ps%l_loc) cycle
        do k = 1, ps%kbc
          call spline_filter_bessel(ps%kb(l, k), l, gmax, alpha, beta_fs, rcut, beta_rs)
        end do
      end do
      
      if(trim(ps%icore).ne.'nc') then
        call spline_filter_bessel(ps%core, 0, gmax, alpha, beta_fs, rcut, beta_rs)
      end if

    end select

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
    write(iunit,'(a,i4)')     'lmax    : ', ps%l_max
    write(iunit,'(a,i4)')     'lloc    : ', ps%l_loc
    write(iunit,'(a,i4,/)')   'kbc     : ', ps%kbc
    write(iunit,'(a,f9.5,/)') 'rcmax   : ', ps%rc_max
    write(iunit,'(a,/)')    'h matrix:'
    do l = 0, ps%l_max
      do k = 1, ps%kbc
        write(iunit,'(10f9.5)') (ps%h(l, k, j), j = 1, ps%kbc)
      end do
    end do
    write(iunit,'(/,a,/)')    'k matrix:'
    do l = 0, ps%l_max
      do k = 1, ps%kbc
        write(iunit,'(10f9.5)') (ps%k(l, k, j), j = 1, ps%kbc)
      end do
    end do
    call io_close(iunit);

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local', action='write')
    call spline_print(ps%vl, iunit)
    call io_close(iunit)

    ! Fourier transform of the local part
    iunit = io_open(trim(dir)//'/local_ft', action='write')
    SAFE_ALLOCATE(fw(1, 1))
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
    SAFE_ALLOCATE(fw(0:ps%l_max, 1:ps%kbc))
    call spline_init(fw)
    do k = 0, ps%l_max
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

    if(trim(ps%icore).ne.'nc') then
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
      call spline_end(ps%vion)
      call spline_end(ps%dvion)
    end if

    call spline_end(ps%kb)
    call spline_end(ps%dkb)
    call spline_end(ps%ur)
    call spline_end(ps%ur_sq)

    call spline_end(ps%vl)
    call spline_end(ps%core)

    call logrid_end(ps%g)

    SAFE_DEALLOCATE_P(ps%kb)
    SAFE_DEALLOCATE_P(ps%dkb)
    SAFE_DEALLOCATE_P(ps%ur)
    SAFE_DEALLOCATE_P(ps%ur_sq)
    SAFE_DEALLOCATE_P(ps%h)
    SAFE_DEALLOCATE_P(ps%k)

    POP_SUB(ps_end)
  end subroutine ps_end


  ! ---------------------------------------------------------
  subroutine hgh_load(ps, psp)
    type(ps_t),  intent(inout) :: ps
    type(hgh_t), intent(inout) :: psp

    integer :: l, ll
    FLOAT :: x

    PUSH_SUB(hgh_load)

    ! Fixes some components of ps, read in psf
    ps%z_val = psp%z_val
    ps%icore = 'nc'
    if(ps%l_max>=0) then
      ps%rc_max = CNST(1.1) * maxval(psp%kbr(0:ps%l_max)) ! Increase a little.
    else
      ps%rc_max = M_ZERO
    end if
    ps%h(0:ps%l_max, 1:ps%kbc, 1:ps%kbc) = psp%h(0:ps%l_max, 1:ps%kbc, 1:ps%kbc)
    ps%k(0:ps%l_max, 1:ps%kbc, 1:ps%kbc) = psp%k(0:ps%l_max, 1:ps%kbc, 1:ps%kbc)

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
      integer :: l, is, nrc, j
      FLOAT, allocatable :: hato(:)

      PUSH_SUB(hgh_load.get_splines)

      SAFE_ALLOCATE(hato(1:psp%g%nrval))

      ! Interpolate the KB-projection functions
      do l = 0, psp%l_max
        do j = 1, 3
          hato = M_ZERO
          nrc = nint(log(psp%kbr(l)/psp%g%b + M_ONE)/psp%g%a) + 1
          hato(1:nrc) = psp%kb(1:nrc, l, j)
          call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%kb(l, j))
        end do
      end do

      ! Now the part corresponding to the local pseudopotential
      ! where the asymptotic part is subtracted
      call spline_fit(psp%g%nrval, psp%g%rofi, psp%vlocal, ps%vl)

      ! Define the table for the pseudo-wavefunction components (using splines)
      ! with a correct normalization function
      do is = 1, ps%ispin
        do l = 1, ps%conf%p
          hato(2:psp%g%nrval) = psp%rphi(2:psp%g%nrval, l)/psp%g%rofi(2:psp%g%nrval)
          hato(1) = hato(2)
          call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%ur(l, is))
          call spline_fit(psp%g%nrval, psp%g%r2ofi, hato, ps%ur_sq(l, is))
        end do
      end do
      
      SAFE_DEALLOCATE_A(hato)

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

    ps%icore = 'nc'
    if(ps_grid%core_corrections) ps%icore=''

    ps%h(0:ps%l_max, 1, 1) = ps_grid%dkbcos(1:ps%l_max+1)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = maxval(ps_grid%kb_radius(1:ps%l_max+1)) * CNST(1.5)

    ! now we fit the splines
    call get_splines(ps_grid%g)

    ! Passes from Rydbergs to Hartrees.
    ps%h(0:ps%l_max,:,:)    = ps%h(0:ps%l_max,:,:)    / M_TWO

    POP_SUB(ps_grid_load)

  contains

    subroutine get_splines(g)
      type(logrid_t), intent(in) :: g

      FLOAT, allocatable :: hato(:)
      integer :: is, l, ir, nrc

      PUSH_SUB(ps_grid_load.get_splines)

      SAFE_ALLOCATE(hato(1:g%nrval))

      ! the wavefunctions
      do is = 1, ps%ispin
        do l = 1, ps_grid%no_l_channels
          hato(2:) = ps_grid%rphi(2:, l, 1+is)/g%rofi(2:)
          hato(1)  = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
            hato(2), hato(3))

          call spline_fit(g%nrval, g%rofi, hato, ps%ur(l, is))
          call spline_fit(g%nrval, g%r2ofi, hato, ps%ur_sq(l, is))

        end do
      end do

      ! the Kleinman-Bylander projectors
      do l = 1, ps%l_max+1
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
      POP_SUB(ps_grid_load.get_splines)
    end subroutine get_splines
  end subroutine ps_grid_load


  ! ---------------------------------------------------------
  subroutine ps_upf_load(ps, ps_upf)
    type(ps_t),     intent(inout) :: ps
    type(ps_upf_t), intent(in)    :: ps_upf

    integer :: i, l, ll, is, nrc, ir, j, ij
    FLOAT :: x
    FLOAT, allocatable :: hato(:)

    PUSH_SUB(ps_upf_load)

    ! Fixes some components of ps, read in ps_upf
    ps%z_val = ps_upf%z_val

    ps%icore = 'nc'
    if(ps_upf%nlcc) ps%icore=''

    ! The spin-dependent pseudopotentials are not suported yet, so we need to fix the occupations
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
    if (ps_upf%l_local >= 0) then
      hato = M_ZERO
      do j = 1, ps%kbc
        call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%l_local, j))
      end do
    end if

    ps%h = M_ZERO
    do i = 1, ps_upf%n_proj

      ij = 1
      if (ps_upf%kb_nc == 2) then
        if (ps_upf%proj_j(i) == ps_upf%proj_l(i) - M_HALF) ij = 2
      end if
      ps%h(ps_upf%proj_l(i), ij, ij) = ps_upf%e(i)*M_TWO

      nrc = logrid_index(ps%g, ps_upf%kb_radius(i)) + 1
      hato(2:nrc) = ps_upf%proj(2:nrc, i)/ps%g%rofi(2:nrc)/M_TWO ! in upf the projector is given in Rydbergs and is multiplied by r
      hato(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), ps%g%rofi(3), hato(2), hato(3)) !take care of the point at zero
      hato(nrc+1:ps%g%nrval) = M_ZERO

      call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), ij))
      if (ps_upf%proj_l(i) == 0 .and. ps_upf%kb_nc == 2) then
        hato = M_ZERO
        call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), 2))
      end if
    end do
 
    ! Define the table for the pseudo-wavefunction components (using splines)
    ! with a correct normalization function
    do is = 1, ps%ispin
      do l = 1, ps%conf%p
        hato = ps_upf%wfs(:, l)/ps%g%rofi
        if (ps%g%rofi(1) == M_ZERO) hato(1) = M_ZERO
        call spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%ur(l, is))
        call spline_fit(ps%g%nrval, ps%g%r2ofi, hato, ps%ur_sq(l, is))
      end do
    end do

    SAFE_DEALLOCATE_A(hato)

    POP_SUB(ps_upf_load)
  end subroutine ps_upf_load


  !> Returns the number of atomic orbitals that can be used for LCAO calculations.
  !!
  pure integer function ps_niwfs(ps)
    type(ps_t), intent(in) :: ps

    integer :: i, l

    ps_niwfs = 0
    select case(ps%flavour)
    case(PS_TYPE_PSF, PS_TYPE_CPI, PS_TYPE_FHI)
      ! In this case we do not count the orbitals whose "l" number is larger than the "l_max" of the pseudo-potential.
      do i = 1, ps%conf%p
        l = ps%conf%l(i)
        if(l <= ps%l_max) ps_niwfs = ps_niwfs + (2*l+1)
      end do

    case(PS_TYPE_HGH, PS_TYPE_UPF)
      do i = 1, ps%conf%p
        l = ps%conf%l(i)
        ps_niwfs = ps_niwfs + (2*l+1)
      end do
    end select

  end function ps_niwfs
  ! ---------------------------------------------------------

end module ps_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
