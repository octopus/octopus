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
  use datasets_m
  use atomic_m
  use global_m
  use io_m
  use loct_gsl_spline_m
  use loct_math_m
  use loct_parser_m
  use logrid_m
  use messages_m
  use ps_cpi_m
  use ps_fhi_m
  use ps_hgh_m
  use ps_in_grid_m
  use ps_psf_m
  use ps_upf_m

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
    ps_end

  integer, parameter, public :: &
    PS_TYPE_PSF = 100,          &
    PS_TYPE_HGH = 101,          &
    PS_TYPE_CPI = 102,          &
    PS_TYPE_FHI = 103,          &
    PS_TYPE_UPF = 104

  character(len=3), parameter  :: ps_name(PS_TYPE_PSF:PS_TYPE_UPF) = (/"tm2", "hgh", "cpi", "fhi", "upf"/)

  type ps_t
    character(len=10) :: label
    integer           :: flavour

    integer  :: ispin    ! Consider spin (ispin = 2) or not (ispin = 1)
    FLOAT    :: z, z_val
    type(valconf_t)   :: conf
    type(logrid_t) :: g
    type(loct_spline_t), pointer :: Ur(:, :)     ! atomic wavefunctions

    ! Kleynman and Bylander projectors stuff
    integer  :: l_max    ! maximum value of l to take
    integer  :: l_loc    ! which component to take as local

    type(loct_spline_t) :: vl         ! local part
    type(loct_spline_t) :: vlocal_f   ! local potential in Fourier space (for periodic)

    FLOAT :: projectors_sphere_threshold ! The projectors are localized in real
                                         ! space, and so they are contained in a 
                                         ! sphere whose radius is computed by
                                         ! making sure that the projector
                                         ! functions absolute value is below this
                                         ! threshold, for points outside the
                                         ! sphere.
    FLOAT :: rc_max ! The radius of the spheres that contain the projector functions.

    integer  :: kbc      ! Number of KB components (1 or 2 for TM ps, 3 for HGH)
    FLOAT, pointer :: h(:,:,:), k(:, :, :)
    type(loct_spline_t), pointer :: kb(:, :)     ! Kleynman-Bylander projectors
    type(loct_spline_t), pointer :: dkb(:, :)    ! derivatives of KB projectors

    ! NLCC
    character(len=4) :: icore
    type(loct_spline_t) :: core ! core charge


    !LONG RANGE PART OF THE LOCAL POTENTIAL
    
    logical :: has_long_range

    type(loct_spline_t) :: vlr         ! the long range part of the local potential
    type(loct_spline_t) :: nlr         ! the charge density associated to the long range part
    
    FLOAT :: sigma_erf                 ! the a constant in erf(r/(sqrt(2)*sigma))/r
    FLOAT :: a_erf                     ! the a constant in erf(ar)/r

    type(loct_spline_t) :: vion        ! the potential that other ions see
    type(loct_spline_t) :: dvion       ! the potential that other ions see
  end type ps_t

  FLOAT, parameter :: eps = CNST(1.0e-8)

contains


  ! ---------------------------------------------------------
  subroutine ps_init(ps, label, flavour, z, lmax, lloc, ispin)
    type(ps_t),        intent(out) :: ps
    character(len=10), intent(in)  :: label
    integer,           intent(in)  :: flavour
    integer,           intent(in)  :: lmax, lloc, ispin
    FLOAT,             intent(in)  :: z

    type(ps_psf_t) :: ps_psf ! SIESTA pseudopotential
    type(ps_cpi_t) :: ps_cpi ! Fritz-haber pseudopotential
    type(ps_fhi_t) :: ps_fhi ! Fritz-haber pseudopotential (from abinit)
    type(ps_upf_t) :: ps_upf ! In case UPF format is used
    type(hgh_t)    :: psp    ! In case Hartwigsen-Goedecker-Hutter ps are used.

    call push_sub('ps.ps_init')

    ! Sets the flavour, label, and number of spin channels.
    ps%flavour = flavour
    ps%label   = label
    ps%ispin   = ispin

    ! Initialization and processing.
    ASSERT(flavour>=PS_TYPE_PSF.and.flavour<=PS_TYPE_UPF)

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
      ps%kbc    = 1     ! only one projector epr angular momentum
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
      ps%l_max  = ps_upf%l_max
      ps%l_loc  = ps_upf%l_local

      nullify(ps%g%drdi, ps%g%s)
      ps%g%nrval = ps_upf%np
      ALLOCATE(ps%g%rofi(ps%g%nrval), ps%g%nrval)
      ps%g%rofi = ps_upf%r

    end select

    ! We allocate all the stuff
    ALLOCATE(ps%kb   (0:ps%l_max, ps%kbc),             (ps%l_max+1)*ps%kbc)
    ALLOCATE(ps%dkb  (0:ps%l_max, ps%kbc),             (ps%l_max+1)*ps%kbc)
    ALLOCATE(ps%ur   (ps%conf%p, ps%ispin),            ps%conf%p*ps%ispin)
    ALLOCATE(ps%h    (0:ps%l_max, 1:ps%kbc, 1:ps%kbc), (ps%l_max+1)*ps%kbc*ps%kbc)
    ALLOCATE(ps%k    (0:ps%l_max, 1:ps%kbc, 1:ps%kbc), (ps%l_max+1)*ps%kbc*ps%kbc)
    call loct_spline_init(ps%kb)
    call loct_spline_init(ps%dkb)
    call loct_spline_init(ps%vl)
    call loct_spline_init(ps%core)
    call loct_spline_init(ps%vlocal_f)

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

    ! Fix the threshold to calculate the radius of the projector function localization spheres:
    !%Variable SpecieProjectorSphereThreshold
    !%Type float
    !%Default 0.001
    !%Section System::Species
    !%Description
    !% The pseudopotentials may be composed of a local part, and a linear combination of nonlocal
    !% operators. These nonlocal projectors have "projector" form, |v><v| (or, more generally
    !% speaking, |u><v|). These projectors are localized in real space -- that is, the function v
    !% has a finite support around the nucleus. This region where the projectors are localized should
    !% be small or else the computation time required to operate with them will be very large.
    !% 
    !% In practice, this localization is fixed by requiring the definition of the projectors to be
    !% contained in a sphere of a certain radius. This radius is computed by making sure that the 
    !% absolute value of the projector functions, at points outside the localization sphere, is 
    !% below a certain threshold. This threshold is set the SpecieProjectorSphereThreshold.
    !%End
    call loct_parse_float(check_inp('SpecieProjectorSphereThreshold'), &
      CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call input_error('SpecieProjectorSphereThreshold')

    ps%has_long_range = .true.

    call pop_sub()

  end subroutine ps_init

  subroutine ps_separate(ps, spacing)
    type(ps_t),        intent(out) :: ps
    FLOAT,             intent(in)  :: spacing

    FLOAT, allocatable :: vsr(:), vlr(:), nlr(:), vion(:)
    FLOAT :: r
    integer :: ii
    
    !separate the local potential in (soft) long range and (hard) short range parts
    
    ps%sigma_erf = CNST(0.625) ! This is hard-coded to a reasonable value.
    ps%a_erf = M_ONE/(ps%sigma_erf*sqrt(M_TWO))
    
    ALLOCATE(vsr(ps%g%nrval), ps%g%nrval)
    ALLOCATE(vlr(ps%g%nrval), ps%g%nrval)
    ALLOCATE(nlr(ps%g%nrval), ps%g%nrval)
    ALLOCATE(vion(ps%g%nrval), ps%g%nrval)
    
    vlr(1) = -ps%z_val*M_TWO/(sqrt(M_TWO*M_PI)*ps%sigma_erf)

    do ii = 1, ps%g%nrval
      r = ps%g%rofi(ii)
      if ( ii > 1) then
        vlr(ii)  = -ps%z_val*loct_erf(r/(ps%sigma_erf*sqrt(M_TWO)))/r
        vion(ii) = -ps%z_val/r - vlr(ii)
      end if
      vsr(ii) = loct_splint(ps%vl, r) - vlr(ii)
      nlr(ii) = -ps%z_val*M_ONE/(ps%sigma_erf*sqrt(M_TWO*M_PI))**3*exp(-M_HALF*r**2/ps%sigma_erf**2)
    end do
    
    call loct_spline_init(ps%vlr)
    call loct_spline_fit(ps%g%nrval, ps%g%rofi, vlr, ps%vlr)
    
    call loct_spline_init(ps%nlr)
    call loct_spline_fit(ps%g%nrval, ps%g%rofi, nlr, ps%nlr)
    
    !overwrite vl
    call loct_spline_end(ps%vl)
    call loct_spline_init(ps%vl)
    call loct_spline_fit(ps%g%nrval, ps%g%rofi, vsr, ps%vl)
    
    ! And take the Fourier transform
    call loct_spline_3dft(ps%vl, ps%vlocal_f, CNST(50.0))
    call loct_spline_times(CNST(1.0)/(M_FOUR*M_PI), ps%vlocal_f)

    ! The ion-ion interaction
    vion(1) = vion(2)
    
    call loct_spline_init(ps%vion)
    call loct_spline_fit(ps%g%nrval, ps%g%rofi, vion, ps%vion)

    call loct_spline_init(ps%dvion)
    call loct_spline_der(ps%vion, ps%dvion)

    deallocate(vsr, vlr, nlr, vion)
    
  end subroutine ps_separate
  
  ! ---------------------------------------------------------
  subroutine ps_getradius(ps)
    type(ps_t), intent(inout) :: ps
    integer :: l, j

    call push_sub('ps.ps_getradius')

    ps%rc_max = CNST(0.0)

    do l = 0, ps%l_max
      do j = 1, ps%kbc
        ps%rc_max = max(ps%rc_max, loct_spline_cutoff_radius(ps%kb(l, j), ps%projectors_sphere_threshold))
      end do
    end do
    
    call pop_sub()
  end subroutine ps_getradius


  ! ---------------------------------------------------------
  subroutine ps_derivatives(ps)
    type(ps_t), intent(inout) :: ps
    integer :: l, j

    call push_sub('ps.ps_derivatives')

    do l = 0, ps%l_max
      do j = 1, ps%kbc
        call loct_spline_der(ps%kb(l, j), ps%dkb(l, j))
      end do
    end do


    call pop_sub()
  end subroutine ps_derivatives


  ! ---------------------------------------------------------
  subroutine ps_filter(ps, gmax, alpha, beta, rcut, beta2)
    type(ps_t), intent(inout) :: ps
    FLOAT, intent(in) :: gmax
    FLOAT, intent(in) :: alpha, beta, rcut, beta2
    integer :: l, k

    call push_sub('ps.ps_filter')

    call loct_spline_filter(ps%vl, fs = (/ alpha*gmax, CNST(100.0) /) )

    do l = 0, ps%l_max
      do k = 1, ps%kbc
        call loct_spline_filter(ps%kb(l, k), l, fs = (/ alpha*gmax, beta /), &
          rs = (/ rcut, beta2 /))
      end do
    end do

    if(trim(ps%icore).ne.'nc') call loct_spline_filter(ps%core, fs = (/ alpha*gmax, CNST(100.0) /) )

    call pop_sub()
  end subroutine ps_filter


  ! ---------------------------------------------------------
  subroutine ps_debug(ps, dir)
    type(ps_t), intent(in) :: ps
    character(len=*), intent(in) :: dir

    ! We will plot also some Fourier transforms.
    type(loct_spline_t), allocatable :: fw(:, :)
    FLOAT, parameter :: gmax = CNST(40.0)

    character(len=30) :: dirname
    integer  :: info_unit                            ! A text file with some basic data.
    integer  :: local_unit, dlocal_unit, localw_unit ! The local part, derivative, and FT.
    integer  :: nl_unit, dnl_unit, nlw_unit          ! Nonlocal part
    integer  :: wave_unit                            ! pseudowavefunctions
    integer  :: so_unit, dso_unit, sow_unit          ! The spin-orbit non-local terms.
    integer  :: j, k, l

    call push_sub('ps.ps_debug')

    ! Opens the files.
    dirname = trim(dir)//'/ps.'//trim(ps%label)
    call io_mkdir(dirname)
    info_unit   = io_open(trim(dirname)//'/info', action='write')
    local_unit  = io_open(trim(dirname)//'/local', action='write')
    dlocal_unit = io_open(trim(dirname)//'/local_derivative', action='write')
    localw_unit = io_open(trim(dirname)//'/local_ft', action='write')
    nl_unit     = io_open(trim(dirname)//'/nonlocal', action='write')
    dnl_unit    = io_open(trim(dirname)//'/nonlocal_derivative', action='write')
    nlw_unit    = io_open(trim(dirname)//'/nonlocal_ft', action='write')
    so_unit     = io_open(trim(dirname)//'/so', action='write')
    dso_unit    = io_open(trim(dirname)//'/so_derivative', action='write')
    sow_unit    = io_open(trim(dirname)//'/so_ft', action='write')
    wave_unit   = io_open(trim(dirname)//'/wavefunctions', action='write')

    ! Writes down the info.
    write(info_unit,'(a,/)')      ps%label
    write(info_unit,'(a,a,/)')    'Flavour : ', ps_name(ps%flavour)
    write(info_unit,'(a,f6.3)')   'z       : ', ps%z
    write(info_unit,'(a,f6.3,/)') 'zval    : ', ps%z_val
    write(info_unit,'(a,i4)')     'lmax    : ', ps%l_max
    write(info_unit,'(a,i4)')     'lloc    : ', ps%l_loc
    write(info_unit,'(a,i4,/)')   'kbc     : ', ps%kbc
    write(info_unit,'(a,f9.5,/)') 'rcmax   : ', ps%rc_max
    write(info_unit,'(/,a,/)')    'h matrix:'
    do l = 0, ps%l_max
      do k = 1, ps%kbc
        write(info_unit,'(3f9.5)') (ps%h(l, k, j), j = 1, ps%kbc)
      end do
      write(info_unit, '(a)')
    end do
    write(info_unit,'(/,a,/)')    'k matrix:'
    do l = 0, ps%l_max
      do k = 1, ps%kbc
        write(info_unit,'(3f9.5)') (ps%k(l, k, j), j = 1, ps%kbc)
      end do
      write(info_unit, '(a)')
    end do

    ! Local part.
    call loct_spline_print(ps%vl, local_unit)
    ALLOCATE(fw(1, 1), 1*1)
    call loct_spline_init(fw(1, 1))
    call loct_spline_3dft(ps%vl, fw(1, 1), gmax = gmax)
    call loct_spline_print(fw(1, 1), localw_unit)
    call loct_spline_end(fw(1, 1))
    deallocate(fw)

    ! Kleinman-Bylander projectors
    call loct_spline_print(ps%kb, nl_unit)
    call loct_spline_print(ps%dkb, dnl_unit)
    ALLOCATE(fw(0:ps%l_max, 1:ps%kbc), (ps%l_max+1)*ps%kbc)
    call loct_spline_init(fw)
    do k = 0, ps%l_max
      do j = 1, ps%kbc
        call loct_spline_3dft(ps%kb(k, j), fw(k, j), gmax = gmax)
      end do
    end do
    call loct_spline_print(fw, nlw_unit)
    call loct_spline_end(fw)
    deallocate(fw)

    ! Pseudo-wavefunctions
    call loct_spline_print(ps%ur, wave_unit)

    ! Closes files and exits
    call io_close(local_unit); call io_close(dlocal_unit); call io_close(localw_unit)
    call io_close(nl_unit)   ; call io_close(dnl_unit)   ; call io_close(nlw_unit)
    call io_close(so_unit)   ; call io_close(dso_unit)   ; call io_close(sow_unit)
    call io_close(info_unit) ; call io_close(wave_unit)

    call pop_sub()
  end subroutine ps_debug


  ! ---------------------------------------------------------
  subroutine ps_end(ps)
    type(ps_t), intent(inout) :: ps

    call push_sub('ps.ps_end')

    if(.not. associated(ps%kb)) return


    call loct_spline_end(ps%kb)
    call loct_spline_end(ps%dkb)
    call loct_spline_end(ps%ur)

    call loct_spline_end(ps%vlr)
    call loct_spline_end(ps%nlr)

    call loct_spline_end(ps%vl)
    call loct_spline_end(ps%core)
    call loct_spline_end(ps%vlocal_f)
    call loct_spline_end(ps%vion)
    call loct_spline_end(ps%dvion)

    call logrid_end(ps%g)

    deallocate(ps%kb, ps%dkb, ps%ur, ps%h, ps%k)

    call pop_sub()
  end subroutine ps_end


  ! ---------------------------------------------------------
  subroutine hgh_load(ps, psp)
    type(ps_t),  intent(inout) :: ps
    type(hgh_t), intent(inout) :: psp

    integer :: l, ll
    FLOAT :: x

    call push_sub('ps.hgh_load')

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

    call pop_sub()

  contains
    ! ---------------------------------------------------------
    subroutine get_splines()
      integer :: l, is, nrc, j
      FLOAT, allocatable :: hato(:)

      call push_sub('ps.get_splines_hgh')

      ALLOCATE(   hato(psp%g%nrval), psp%g%nrval)

      ! Interpolate the KB-projection functions
      do l = 0, psp%l_max
        do j = 1, 3
          hato = M_ZERO
          nrc = nint(log(psp%kbr(l)/psp%g%b + M_ONE)/psp%g%a) + 1
          hato(1:nrc) = psp%kb(1:nrc, l, j)
          call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%kb(l, j))
        end do
      end do

      ! Now the part corresponding to the local pseudopotential
      ! where the asymptotic part is substracted
      call loct_spline_fit(psp%g%nrval, psp%g%rofi, psp%vlocal, ps%vl)

      ! Define the table for the pseudo-wavefunction components (using splines)
      ! with a correct normalization function
      do is = 1, ps%ispin
        do l = 1, ps%conf%p
          hato(2:psp%g%nrval) = psp%rphi(2:psp%g%nrval, l)/psp%g%rofi(2:psp%g%nrval)
          hato(1) = hato(2)
          call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%Ur(l, is))
        end do
      end do
      
      call pop_sub()
    end subroutine get_splines
  end subroutine hgh_load


  ! ---------------------------------------------------------
  subroutine ps_grid_load(ps, ps_grid)
    type(ps_t),         intent(inout) :: ps
    type(ps_in_grid_t), intent(in)  :: ps_grid

    call push_sub('ps.ps_grid_load')

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

    call pop_sub()

  contains

    subroutine get_splines(g)
      type(logrid_t), intent(in) :: g

      FLOAT, allocatable :: hato(:)
      integer :: is, l, ir, nrc

      ALLOCATE(hato(g%nrval), g%nrval)

      ! the wave-functions
      do is = 1, ps%ispin
        do l = 1, ps_grid%no_l_channels
          hato(2:) = ps_grid%rphi(2:, l, 1+is)/g%rofi(2:)
          hato(1)  = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
            hato(2), hato(3))

          call loct_spline_fit(g%nrval, g%rofi, hato, ps%ur(l, is))

        end do
      end do

      ! the Kleinman-Bylander projectors
      do l = 1, ps%l_max+1
        nrc = logrid_index(g, ps_grid%kb_radius(l)) + 1
        hato(1:nrc)         = ps_grid%KB(1:nrc, l)
        hato(nrc+1:g%nrval) = M_ZERO

        call loct_spline_fit(g%nrval, g%rofi, hato, ps%kb(l-1, 1))
      end do

      ! Now the part corresponding to the local pseudopotential
      ! where the asymptotic part is substracted
      hato(:) = ps_grid%vlocal(:)/M_TWO
      call loct_spline_fit(g%nrval, g%rofi, hato, ps%vl)
      
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

        call loct_spline_fit(g%nrval, g%rofi, hato, ps%core)
      end if

      deallocate(hato)
    end subroutine get_splines
  end subroutine ps_grid_load


  ! ---------------------------------------------------------
  subroutine ps_upf_load(ps, ps_upf)
    type(ps_t),         intent(out) :: ps
    type(ps_upf_t), intent(in)  :: ps_upf

    integer :: i, l, ll, is, nrc, ir, j, ij
    FLOAT :: x
    FLOAT, allocatable :: hato(:)

    call push_sub('ps.ps_upf_load')

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

    ALLOCATE(hato(ps%g%nrval), ps%g%nrval)


    !Non-linear core-corrections
    if(ps_upf%nlcc) then
      ! find cutoff radius
      hato = ps_upf%core_density/(M_FOUR*M_PI*ps%g%rofi**2)

      do ir = ps%g%nrval-1, 1, -1
        if(hato(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do

      hato(nrc:ps%g%nrval) = M_ZERO

      call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%core)
    end if

    ! Now the part corresponding to the local pseudopotential
    ! where the asymptotic part is substracted
    hato(:) = ps_upf%v_local/M_TWO
    call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%vl)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = maxval(ps_upf%kb_radius)
    ps%rc_max = max(ps_upf%local_radius, ps%rc_max) * CNST(1.5)

    ! Interpolate the KB-projection functions
    if (ps_upf%l_local >= 0) then
      hato = M_ZERO
      do j = 1, ps%kbc
        call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%l_local, j))
      end do
    end if

    ps%h = M_ZERO
    do i = 1, ps_upf%n_proj

      ij = 1
      if (ps_upf%kb_nc == 2) then
        if (ps_upf%proj_j(i) == ps_upf%proj_l(i) - M_HALF) ij = 2
      end if
      ps%h(ps_upf%proj_l(i), ij, ij) = ps_upf%e(i)

      nrc = logrid_index(ps%g, ps_upf%kb_radius(i)) + 1
      hato(1:nrc)         = ps_upf%proj(1:nrc, i)
      hato(nrc+1:ps%g%nrval) = M_ZERO

      call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), ij))
      if (ps_upf%proj_l(i) == 0 .and. ps_upf%kb_nc == 2) then
        hato = M_ZERO
        call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%kb(ps_upf%proj_l(i), 2))
      end if
    end do

    ! Passes from Rydbergs to Hartrees.
    ps%h(0:ps%l_max,:,:) = ps%h(0:ps%l_max,:,:)/M_TWO

    ! Define the table for the pseudo-wavefunction components (using splines)
    ! with a correct normalization function
    do is = 1, ps%ispin
      do l = 1, ps%conf%p
        hato = ps_upf%wfs(:, l)/ps%g%rofi
        call loct_spline_fit(ps%g%nrval, ps%g%rofi, hato, ps%Ur(l, is))
      end do
    end do

    deallocate(hato)

    call pop_sub()
  end subroutine ps_upf_load

end module ps_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
