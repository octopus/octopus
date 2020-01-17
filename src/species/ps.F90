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
  use lalg_adv_oct_m
  use loct_math_oct_m
  use parser_oct_m
  use logrid_oct_m
  use messages_oct_m
  use namespace_oct_m
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
  use pseudo_oct_m
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
    ps_bound_niwfs,             &
    ps_end,                     &
    ps_has_density,             &
    ps_has_nlcc,                &
    ps_density_volume

  integer, parameter, public :: &
    PS_FILTER_NONE = 0,         &
    PS_FILTER_TS   = 2,         &
    PS_FILTER_BSB  = 3

  integer, public, parameter ::  &
    PROJ_NONE = 0,  &
    PROJ_HGH  = 1,  &
    PROJ_KB   = 2,  &
    PROJ_RKB  = 3

  integer, parameter, public :: INVALID_L = 333

  character(len=4), parameter  :: ps_name(PSEUDO_FORMAT_UPF1:PSEUDO_FORMAT_PSP8) = &
    (/"upf1", "upf2", "qso ", "psml", "psf ", "cpi ", "fhi ", "hgh ", "psp8"/)

  type ps_t
    ! Components are public by default
    integer :: projector_type
    character(len=10), private :: label

    integer, private  :: ispin    !< Consider spin (ispin = 2) or not (ispin = 1)
    FLOAT, private    :: z
    FLOAT    :: z_val
    type(valconf_t)   :: conf
    type(logrid_t), private :: g
    type(spline_t), pointer :: ur(:, :)     !< (1:conf%p, 1:ispin) atomic wavefunctions, as a function of r
    type(spline_t), pointer, private :: ur_sq(:, :)  !< (1:conf%p, 1:ispin) atomic wavefunctions, as a function of r^2
    logical, allocatable    :: bound(:, :)  !< (1:conf%p, 1:ispin) is the state bound or not

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
    type(spline_t) :: core_der !< derivative of the core correction


    !LONG-RANGE PART OF THE LOCAL POTENTIAL

    logical, private :: has_long_range

    type(spline_t), private :: vlr !< the long-range part of the local potential
    type(spline_t) :: vlr_sq       !< the long-range part of the
    !< local potential in terms of r^2, to avoid the sqrt
    type(spline_t) :: nlr          !< the charge density associated with the long-range part

    FLOAT :: sigma_erf             !< the a constant in erf(r/(sqrt(2)*sigma))/r

    logical,        private :: has_density     !< does the species have a density?
    type(spline_t), pointer :: density(:)      !< the atomic density for each spin
    type(spline_t), pointer :: density_der(:)  !< the radial derivative for the atomic density for each spin

    logical, private :: is_separated
    logical          :: local
    logical          :: hamann
    integer, private :: file_format
    integer, private :: pseudo_type
    integer          :: exchange_functional
    integer          :: correlation_functional
  end type ps_t

  FLOAT, parameter :: eps = CNST(1.0e-8)

contains


  ! ---------------------------------------------------------
  subroutine ps_init(ps, namespace, label, z, user_lmax, user_llocal, ispin, filename)
    type(ps_t),        intent(out)   :: ps
    type(namespace_t), intent(in)    :: namespace
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
    type(ps_hgh_t) :: ps_hgh !< In case Hartwigsen-Goedecker-Hutter ps are used.
    type(ps_xml_t) :: ps_xml !< For xml based pseudopotentials
    logical, save :: xml_warned = .false.
    FLOAT, allocatable :: eigen(:, :)  !< eigenvalues

    PUSH_SUB(ps_init)

    ps%exchange_functional = PSEUDO_EXCHANGE_UNKNOWN
    ps%correlation_functional = PSEUDO_CORRELATION_UNKNOWN

    ! Fix the threshold to calculate the radius of the projector-function localization spheres:

    call messages_obsolete_variable(namespace, 'SpecieProjectorSphereThreshold', 'SpeciesProjectorSphereThreshold')

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
    call parse_variable(namespace, 'SpeciesProjectorSphereThreshold', CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call messages_input_error('SpeciesProjectorSphereThreshold')

    ps%file_format = pseudo_detect_format(filename)

    if(ps%file_format == PSEUDO_FORMAT_FILE_NOT_FOUND) then
      call messages_write("Cannot open pseudopotential file '"//trim(filename)//"'.")
      call messages_fatal(namespace=namespace)
    end if

    if(ps%file_format == PSEUDO_FORMAT_UNKNOWN) then
      call messages_write("Cannot determine the pseudopotential type for species '"//trim(label)//"' from", &
        new_line = .true.)
      call messages_write("file '"//trim(filename)//"'.")
      call messages_fatal(namespace=namespace)
    end if

    ps%label   = label
    ps%ispin   = ispin
    ps%hamann  = .false.
    ps%projector_type = PROJ_KB

    select case(ps%file_format)
    case(PSEUDO_FORMAT_PSF, PSEUDO_FORMAT_HGH)
      ps%has_density = .true.
    case default
      ps%has_density = .false.
    end select

    select case(ps%file_format)
    case(PSEUDO_FORMAT_PSF)
      ps%pseudo_type   = PSEUDO_TYPE_SEMILOCAL

      call ps_psf_init(ps_psf, ispin, filename, namespace)

      call valconf_copy(ps%conf, ps_psf%conf)
      ps%z      = z
      ps%conf%z = nint(z) ! atomic number
      ps%kbc    = 1     ! only one projector per angular momentum

      ps%lmax = ps_psf%ps_grid%no_l_channels - 1

      if(user_lmax /= INVALID_L) then
        ps%lmax = min(ps%lmax, user_lmax) ! Maybe the file does not have enough components.
        if(user_lmax /= ps%lmax) then
          message(1) = "lmax in Species block for " // trim(label) // &
            " is larger than number available in pseudopotential."
          call messages_fatal(1, namespace=namespace)
        end if
      end if

      ps%conf%p = ps_psf%ps_grid%no_l_channels
      if(ps%lmax == 0) ps%llocal = 0 ! Vanderbilt is not acceptable if ps%lmax == 0.

      ! the local part of the pseudo
      if(user_llocal == INVALID_L) then
        ps%llocal = 0
      else
        ps%llocal = user_llocal
      end if

      call ps_psf_process(ps_psf, namespace, ps%lmax, ps%llocal)
      call logrid_copy(ps_psf%ps_grid%g, ps%g)

    case(PSEUDO_FORMAT_CPI, PSEUDO_FORMAT_FHI)
      ps%pseudo_type   = PSEUDO_TYPE_SEMILOCAL

      call valconf_null(ps%conf)

      if(ps%file_format == PSEUDO_FORMAT_CPI) then
        call ps_cpi_init(ps_cpi, trim(filename), namespace)
        ps%conf%p      = ps_cpi%ps_grid%no_l_channels
      else
        call ps_fhi_init(ps_fhi, trim(filename), namespace)
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

      if(user_lmax /= INVALID_L) then
        ps%lmax = min(ps%lmax, user_lmax) ! Maybe the file does not have enough components.
        if(user_lmax /= ps%lmax) then
          message(1) = "lmax in Species block for " // trim(label) // &
            " is larger than number available in pseudopotential."
          call messages_fatal(1, namespace=namespace)
        end if
      end if

      if(ps%lmax == 0) ps%llocal = 0 ! Vanderbilt is not acceptable if ps%lmax == 0.

      ! the local part of the pseudo
      if(user_llocal == INVALID_L) then
        ps%llocal = 0
      else
        ps%llocal = user_llocal
      end if

      if(ps%file_format == PSEUDO_FORMAT_CPI) then
        call ps_cpi_process(ps_cpi, ps%llocal, namespace)
        call logrid_copy(ps_cpi%ps_grid%g, ps%g)
      else
        call ps_fhi_process(ps_fhi, ps%lmax, ps%llocal, namespace)
        call logrid_copy(ps_fhi%ps_grid%g, ps%g)
      end if

    case(PSEUDO_FORMAT_HGH)
      ps%pseudo_type   = PSEUDO_TYPE_KLEINMAN_BYLANDER
      ps%projector_type = PROJ_HGH

      call hgh_init(ps_hgh, trim(filename), namespace)
      call valconf_copy(ps%conf, ps_hgh%conf)

      ps%z        = z
      ps%kbc      = 3
      ps%llocal    = -1
      ps%lmax    = ps_hgh%l_max

      call hgh_process(ps_hgh, namespace)
      call logrid_copy(ps_hgh%g, ps%g)

    case(PSEUDO_FORMAT_QSO, PSEUDO_FORMAT_UPF1, PSEUDO_FORMAT_UPF2, PSEUDO_FORMAT_PSML, PSEUDO_FORMAT_PSP8)

      if(.not. xml_warned) then
        call messages_experimental('XML (QSO, UPF, and PSML, PSP8) pseudopotential support')
        xml_warned = .true.
      end if

      call ps_xml_init(ps_xml, namespace, trim(filename), ps%file_format, ierr)

      ps%pseudo_type   = pseudo_type(ps_xml%pseudo)
      ps%exchange_functional = pseudo_exchange(ps_xml%pseudo)
      ps%correlation_functional = pseudo_correlation(ps_xml%pseudo)

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
        ps%llocal = 0                                     ! the default
        if(ps_xml%llocal >= 0) ps%llocal = ps_xml%llocal  ! the one given in the pseudopotential file
        if(user_llocal /= INVALID_L) ps%llocal = user_llocal ! user supplied local component
        ASSERT(ps%llocal >= 0)
        ASSERT(ps%llocal <= ps%lmax)
      end if

      nullify(ps%g%drdi, ps%g%s)

      ps%g%nrval = ps_xml%grid_size

      SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
      SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))

      do ii = 1, ps%g%nrval
        ps%g%rofi(ii) = ps_xml%grid(ii)
        ps%g%r2ofi(ii) = ps_xml%grid(ii)**2
      end do

    end select

    ps%local = (ps%lmax == 0 .and. ps%llocal == 0 ) .or. (ps%lmax == -1 .and. ps%llocal == -1)

    ! We allocate all the stuff
    SAFE_ALLOCATE(ps%kb   (0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb  (0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%bound(1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%h    (0:ps%lmax, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(ps%density(1:ps%ispin))
    SAFE_ALLOCATE(ps%density_der(1:ps%ispin))

    nullify(ps%k)

    call spline_init(ps%kb)
    call spline_init(ps%dkb)
    call spline_init(ps%vl)
    call spline_init(ps%core)
    call spline_init(ps%core_der)
    call spline_init(ps%density)
    call spline_init(ps%density_der)

    SAFE_ALLOCATE(eigen(1:ps%conf%p, 1:ps%ispin))
    eigen = M_ZERO

    ! Now we load the necessary information.
    select case(ps%file_format)
    case(PSEUDO_FORMAT_PSF)
      call ps_psf_get_eigen(ps_psf, eigen)
      call ps_grid_load(ps, ps_psf%ps_grid)
      call ps_psf_end(ps_psf)
    case(PSEUDO_FORMAT_CPI)
      call ps_grid_load(ps, ps_cpi%ps_grid)
      call ps_cpi_end(ps_cpi)
    case(PSEUDO_FORMAT_FHI)
      call ps_grid_load(ps, ps_fhi%ps_grid)
      call ps_fhi_end(ps_fhi)
    case(PSEUDO_FORMAT_HGH)
      call hgh_get_eigen(ps_hgh, eigen)
      SAFE_ALLOCATE(ps%k    (0:ps%lmax, 1:ps%kbc, 1:ps%kbc))
      call hgh_load(ps, ps_hgh)
      call hgh_end(ps_hgh)
    case(PSEUDO_FORMAT_QSO, PSEUDO_FORMAT_UPF1, PSEUDO_FORMAT_UPF2, PSEUDO_FORMAT_PSML, PSEUDO_FORMAT_PSP8)
      call ps_xml_load(ps, ps_xml)
      call ps_xml_end(ps_xml)
    end select

    if(ps_has_density(ps)) then
      do is = 1, ps%ispin
        call spline_der(ps%density(is), ps%density_der(is))
      end do
    end if

    if(ps_has_nlcc(ps)) then
      call spline_der(ps%core, ps%core_der)
    end if

    call ps_check_bound(ps, eigen)

    ps%has_long_range = .true.
    ps%is_separated = .false.

    call ps_info(ps, filename)

    SAFE_DEALLOCATE_A(eigen)

    POP_SUB(ps_init)
  end subroutine ps_init

  !------------------------------------------------------------------------

  subroutine ps_info(ps, filename)
    type(ps_t),       intent(in) :: ps
    character(len=*), intent(in) :: filename

    call messages_write("  Species '"//trim(ps%label)//"'", new_line = .true.)
    call messages_write("    type             : pseudopotential", new_line = .true.)
    call messages_write("    file             : '"//trim(filename)//"'")
    call messages_info()

    call messages_write("    file format      :")
    select case(ps%file_format)
    case(PSEUDO_FORMAT_UPF1)
      call messages_write(" UPF1")
    case(PSEUDO_FORMAT_UPF2)
      call messages_write(" UPF2")
    case(PSEUDO_FORMAT_QSO)
      call messages_write(" QSO")
    case(PSEUDO_FORMAT_PSML)
      call messages_write(" PSML")
    case(PSEUDO_FORMAT_PSP8)
      call messages_write(" PSP8")
    case(PSEUDO_FORMAT_PSF)
      call messages_write(" PSF")
    case(PSEUDO_FORMAT_CPI)
      call messages_write(" CPI")
    case(PSEUDO_FORMAT_FHI)
      call messages_write(" FHI")
    case(PSEUDO_FORMAT_HGH)
      call messages_write(" HGH")
    end select
    call messages_new_line()

    call messages_write("    valence charge   :")
    call messages_write(ps%z_val, align_left = .true., fmt = '(f4.1)')
    call messages_info()

    call messages_write("    atomic number    :")
    call messages_write(nint(ps%z), fmt = '(i4)')
    call messages_info()

    call messages_write("    form on file     :")
    select case(ps%pseudo_type)
    case(PSEUDO_TYPE_ULTRASOFT)
      call messages_write(" ultrasoft")
    case(PSEUDO_TYPE_SEMILOCAL)
      call messages_write(" semilocal")
    case(PSEUDO_TYPE_KLEINMAN_BYLANDER)
      call messages_write(" kleinman-bylander")
    case(PSEUDO_TYPE_PAW)
      call messages_write(" paw")
    end select
    call messages_info()

    if(ps%pseudo_type == PSEUDO_TYPE_SEMILOCAL) then
      call messages_write("    orbital origin   :")
      select case(ps%file_format)
      case(PSEUDO_FORMAT_PSF)
        call messages_write(" calculated");
      case default
        call messages_write(" from file");
      end select
      call messages_info()
    end if

    call messages_write("    lmax             :")
    call messages_write(ps%lmax, fmt = '(i2)')
    call messages_info()

    call messages_write("    llocal           :")
    if(ps%llocal >= 0) then
      call messages_write(ps%llocal, fmt = '(i2)')
    else
      call messages_write(ps%llocal, fmt = '(i3)')
    end if
    call messages_info()

    call messages_write("    projectors per l :")
    call messages_write(ps%kbc, fmt = '(i2)')
    call messages_info()

    call messages_write("    total projectors :")
    if(ps%llocal < 0) then
      call messages_write(ps%kbc*(ps%lmax + 1), fmt = '(i2)')
    else
      call messages_write(ps%kbc*ps%lmax, fmt = '(i2)')
    end if
    call messages_info()

    if(ps%local) then
      call messages_write("    application form : local")
    else
      call messages_write("    application form : kleinman-bylander")
    end if
    call messages_info()

    call messages_write("    orbitals         :")
    call messages_write(ps_niwfs(ps), fmt='(i3)')
    call messages_info()
    call messages_write("    bound orbitals   :")
    call messages_write(ps_bound_niwfs(ps), fmt='(i3)')
    call messages_info()

    call messages_info()

  end subroutine ps_info


  ! ---------------------------------------------------------
  !> separate the local potential into (soft) long-ranged and (hard) short-ranged parts
  subroutine ps_separate(ps)
    type(ps_t),        intent(inout) :: ps

    FLOAT, allocatable :: vsr(:), vlr(:), nlr(:)
    FLOAT :: r, exp_arg
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

      exp_arg = -M_HALF*r**2/ps%sigma_erf**2
      if(exp_arg > CNST(-100)) then
        nlr(ii) = -ps%z_val*M_ONE/(ps%sigma_erf*sqrt(M_TWO*M_PI))**3*exp(exp_arg)
      else
        nlr(ii) = M_ZERO
      end if
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
      if(l == ps%llocal) cycle
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

      if(ps_has_nlcc(ps)) then
        rmax = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold)
        call spline_filter_mask(ps%core, 0, rmax, gmax, alpha, gamma)
      end if

      if(ps_has_density(ps)) then
        do ispin = 1, ps%ispin
          if(abs(spline_integral(ps%density(ispin))) > CNST(1.0e-12)) then
            rmax = spline_cutoff_radius(ps%density(ispin), ps%projectors_sphere_threshold)
            call spline_filter_mask(ps%density(ispin), 0, rmax, gmax, alpha, gamma)
            call spline_force_pos(ps%density(ispin))
          end if

          if(abs(spline_integral(ps%density_der(ispin))) > CNST(1.0e-12)) then
            rmax = spline_cutoff_radius(ps%density_der(ispin), ps%projectors_sphere_threshold)
            call spline_filter_mask(ps%density_der(ispin), 0, rmax, gmax, alpha, gamma)
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

      if(ps_has_nlcc(ps)) then
        call spline_filter_bessel(ps%core, 0, gmax, alpha, beta_fs, rcut, beta_rs)
      end if

      if(ps_has_density(ps)) then
        do ispin = 1, ps%ispin
          call spline_filter_bessel(ps%density(ispin), 0, gmax, alpha, beta_fs, rcut, beta_rs)
          call spline_force_pos(ps%density(ispin))
          call spline_filter_bessel(ps%density_der(ispin), 0, gmax, alpha, beta_fs, rcut, beta_rs)
        end do
      end if

    end select

    call profiling_out(prof)
    POP_SUB(ps_filter)
  end subroutine ps_filter

  ! ---------------------------------------------------------
  subroutine ps_check_bound(ps, eigen)
    type(ps_t), intent(inout) :: ps
    FLOAT,      intent(in)    :: eigen(:,:)

    integer :: i, is, ir
    FLOAT :: ur1, ur2

    PUSH_SUB(ps_check_bound)

    ! Unbound states have positive eigenvalues
    where(eigen > M_ZERO)
      ps%bound = .false.
    elsewhere
      ps%bound = .true.
    end where

    ! We might not have information about the eigenvalues, so we need to check the wavefunctions
    do i = 1, ps%conf%p
      do is = 1, ps%ispin
        if (.not. ps%bound(i, is)) cycle

        do ir = ps%g%nrval, 3, -1
          ! First we look for the outmost value that is not zero
          if (abs(spline_eval(ps%ur(i, is), ps%g%rofi(ir))*ps%g%rofi(ir)) > M_ZERO) then
            ! Usually bound states have exponentially decaying wavefunctions,
            ! while unbound states have exponentially diverging
            ! wavefunctions. Therefore we check if the wavefunctions
            ! value is increasing with increasing radius. The fact
            ! that we do not use the wavefunctions outmost value that
            ! is not zero is on purpose, as some pseudopotential
            ! generators do funny things with that point.
            ur1 = spline_eval(ps%ur(i, is), ps%g%rofi(ir-2))*ps%g%rofi(ir-2)
            ur2 = spline_eval(ps%ur(i, is), ps%g%rofi(ir-1))*ps%g%rofi(ir-1)
            if ((ur1*ur2 > M_ZERO) .and. (abs(ur2) > abs(ur1))) ps%bound(i, is) = .false.
            exit
          end if
        end do
      end do
    end do

    POP_SUB(ps_check_bound)
  end subroutine ps_check_bound


  ! ---------------------------------------------------------
  subroutine ps_debug(ps, dir, namespace)
    type(ps_t), intent(in) :: ps
    character(len=*), intent(in) :: dir
    type(namespace_t), intent(in) :: namespace

    ! We will plot also some Fourier transforms.
    type(spline_t), allocatable :: fw(:, :)
    FLOAT, parameter :: gmax = CNST(40.0)

    integer  :: iunit
    integer  :: j, k, l

    PUSH_SUB(ps_debug)

    ! A text file with some basic data.
    iunit = io_open(trim(dir)//'/pseudo-info', namespace, action='write')
    write(iunit,'(a,/)')      ps%label
    write(iunit,'(a,a,/)')    'Format  : ', ps_name(ps%file_format)
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

    write(iunit,'(/,a)')    'orbitals:'
    do j = 1, ps%conf%p
      write(iunit,'(1x,a,i2,3x,a,i2,3x,a,f5.1,3x,a,l1)') 'n = ', ps%conf%n(j), 'l = ', ps%conf%l(j), &
        'j = ', ps%conf%j(j), 'bound = ', all(ps%bound(j,:))
    end do


    call io_close(iunit)

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local', namespace, action='write')
    call spline_print(ps%vl, iunit)
    call io_close(iunit)

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local_long_range', namespace, action='write')
    call spline_print(ps%vlr, iunit)
    call io_close(iunit)

    ! Local part of the pseudopotential
    iunit  = io_open(trim(dir)//'/local_long_range_density', namespace, action='write')
    call spline_print(ps%nlr, iunit)
    call io_close(iunit)

    ! Fourier transform of the local part
    iunit = io_open(trim(dir)//'/local_ft', namespace, action='write')
    SAFE_ALLOCATE(fw(1:1, 1:1))
    call spline_init(fw(1, 1))
    call spline_3dft(ps%vl, fw(1, 1), gmax = gmax)
    call spline_print(fw(1, 1), iunit)
    call spline_end(fw(1, 1))
    SAFE_DEALLOCATE_A(fw)
    call io_close(iunit)

    ! Kleinman-Bylander projectors
    iunit = io_open(trim(dir)//'/nonlocal', namespace, action='write')
    call spline_print(ps%kb, iunit)
    call io_close(iunit)

    iunit = io_open(trim(dir)//'/nonlocal_derivative', namespace, action='write')
    call spline_print(ps%dkb, iunit)
    call io_close(iunit)

    iunit = io_open(trim(dir)//'/nonlocal_ft', namespace, action='write')
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
    iunit = io_open(trim(dir)//'/wavefunctions', namespace, action='write')
    call spline_print(ps%ur, iunit)
    call io_close(iunit)

    ! Density
    if (ps%has_density) then
      iunit = io_open(trim(dir)//'/density', namespace, action='write')
      call spline_print(ps%density, iunit)
      call io_close(iunit)

      iunit = io_open(trim(dir)//'/density_derivative', namespace, action='write')
      call spline_print(ps%density_der, iunit)
      call io_close(iunit)
    end if

    ! Non-linear core-corrections
    if(ps_has_nlcc(ps)) then
      iunit = io_open(trim(dir)//'/nlcc', namespace, action='write')
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
    call spline_end(ps%core_der)

    if(associated(ps%density)) call spline_end(ps%density)
    if(associated(ps%density_der)) call spline_end(ps%density_der)

    call logrid_end(ps%g)

    SAFE_DEALLOCATE_P(ps%kb)
    SAFE_DEALLOCATE_P(ps%dkb)
    SAFE_DEALLOCATE_P(ps%ur)
    SAFE_DEALLOCATE_P(ps%ur_sq)
    SAFE_DEALLOCATE_A(ps%bound)
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
          hato(1)  = first_point_extrapolate(g%rofi, hato)

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
        hato(1) = first_point_extrapolate(g%rofi, hato)

        call spline_fit(g%nrval, g%rofi, hato, ps%core)
      end if

      SAFE_DEALLOCATE_A(hato)
      SAFE_DEALLOCATE_A(dens)

      POP_SUB(ps_grid_load.get_splines)
    end subroutine get_splines
  end subroutine ps_grid_load

  ! ---------------------------------------------------------

  subroutine ps_xml_load(ps, ps_xml)
    type(ps_t),     intent(inout) :: ps
    type(ps_xml_t), intent(in)    :: ps_xml

    integer :: ll, ip, is, ic, jc, ir, nrc, ii
    FLOAT :: rr, kbcos, kbnorm, dnrm, avgv, volume_element
    FLOAT, allocatable :: vlocal(:), kbprojector(:), wavefunction(:), nlcc_density(:), dens(:)
    integer, allocatable :: cmap(:, :)
    FLOAT, allocatable :: matrix(:, :), eigenvalues(:)

    PUSH_SUB(ps_xml_load)

    ps%hamann = (ps_xml%kleinman_bylander .and. ps_xml%nchannels == 2 .and. ps_xml%llocal == -1)

    ps%nlcc = ps_xml%nlcc

    ps%z_val = ps_xml%valence_charge

    ! the local potential
    SAFE_ALLOCATE(vlocal(1:ps%g%nrval))

    do ip = 1, ps%g%nrval
      rr = ps_xml%grid(ip)
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

      SAFE_ALLOCATE(cmap(0:ps_xml%lmax, 1:ps_xml%nchannels))

      ! the order of the channels is determined by spin orbit and the j value
      do ll = 0, ps_xml%lmax
        do ic = 1, ps_xml%nchannels
          cmap(ll, ic) = ic

          if(ll == 0) cycle
          if(ll == ps_xml%llocal) cycle
          if(.not. pseudo_has_total_angular_momentum(ps_xml%pseudo)) cycle

          ASSERT(ps_xml%nchannels == 2)
          if(pseudo_projector_2j(ps_xml%pseudo, ll, ic) == 2*ll - 1) then
            ! this is Octopus convention
            cmap(ll, ic) = 2
          else
            ASSERT(pseudo_projector_2j(ps_xml%pseudo, ll, ic) == 2*ll + 1)
            cmap(ll, ic) = 1
          end if

        end do

        ! check that all numbers are present for each l
        ASSERT(sum(cmap(ll, 1:ps_xml%nchannels)) == (ps_xml%nchannels + 1)*ps_xml%nchannels/2)
      end do

      ASSERT(all(cmap >= 0 .and. cmap <= ps_xml%nchannels))

      SAFE_ALLOCATE(matrix(1:ps_xml%nchannels, 1:ps_xml%nchannels))
      SAFE_ALLOCATE(eigenvalues(1:ps_xml%nchannels))

      ps%h = CNST(0.0)


      if(pseudo_nprojectors(ps_xml%pseudo) > 0) then
        do ll = 0, ps_xml%lmax

          if (is_diagonal(ps_xml%nchannels, ps_xml%dij(ll, :, :)) .or. &
            pseudo_has_total_angular_momentum(ps_xml%pseudo)) then
            matrix = CNST(0.0)
            forall(ic = 1:ps_xml%nchannels)
              eigenvalues(ic) = ps_xml%dij(ll, ic, ic)
              matrix(ic, ic) = CNST(1.0)
            end forall
          else
            ! diagonalize the coefficient matrix
            matrix(1:ps_xml%nchannels, 1:ps_xml%nchannels) = ps_xml%dij(ll, 1:ps_xml%nchannels, 1:ps_xml%nchannels)
            call lalg_eigensolve(ps_xml%nchannels, matrix, eigenvalues)
          end if

          do ic = 1, ps_xml%nchannels

            do ip = 1, ps%g%nrval
              kbprojector(ip) = M_ZERO
              if(ip <= ps_xml%grid_size) then
                do jc = 1, ps_xml%nchannels
                  kbprojector(ip) = kbprojector(ip) + matrix(jc, ic)*ps_xml%projector(ip, ll, jc)
                end do
              end if
            end do

            call spline_fit(ps%g%nrval, ps%g%rofi, kbprojector, ps%kb(ll, cmap(ll, ic)))

            ps%h(ll, cmap(ll, ic), cmap(ll, ic)) = eigenvalues(ic)

          end do
        end do
      end if

      SAFE_DEALLOCATE_A(matrix)
      SAFE_DEALLOCATE_A(eigenvalues)

      ps%conf%p = ps_xml%nwavefunctions

      do ii = 1, ps_xml%nwavefunctions

        ps%conf%n(ii) = ps_xml%wf_n(ii)
        ps%conf%l(ii) = ps_xml%wf_l(ii)

        if(ps%ispin == 2) then
          ps%conf%occ(ii, 1) = min(ps_xml%wf_occ(ii), CNST(2.0)*ps_xml%wf_l(ii) + CNST(1.0))
          ps%conf%occ(ii, 2) = ps_xml%wf_occ(ii) - ps%conf%occ(ii, 1)
        else
          ps%conf%occ(ii, 1) = ps_xml%wf_occ(ii)
        end if

        ps%conf%j(ii) = M_ZERO
        if(pseudo_has_total_angular_momentum(ps_xml%pseudo)) then
          ps%conf%j(ii) = M_HALF*pseudo_wavefunction_2j(ps_xml%pseudo, ii)
        end if

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

      SAFE_DEALLOCATE_A(cmap)

    else

      do ll = 0, ps_xml%lmax
        ! we need to build the KB projectors
        ! the procedure was copied from ps_in_grid.F90 (r12967)
        dnrm = M_ZERO
        avgv = M_ZERO
        do ip = 1, ps_xml%grid_size
          rr = ps_xml%grid(ip)
          volume_element = rr**2*ps_xml%weights(ip)
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

  logical function is_diagonal(dim, matrix)
    integer, intent(in)    :: dim
    FLOAT,   intent(in)    :: matrix(:, :)

    integer :: ii, jj

    is_diagonal = .true.
    do ii = 1, dim
      do jj = 1, dim
        if(ii == jj) cycle
        if(abs(matrix(ii, jj)) > CNST(1e10)) is_diagonal = .false.
      end do
    end do

  end function is_diagonal

  ! ---------------------------------------------------------
  !> Returns the number of atomic orbitals taking into account then m quantum number multiplicity
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
  !> Returns the number of bound atomic orbitals taking into account then m quantum number multiplicity
  pure integer function ps_bound_niwfs(ps)
    type(ps_t), intent(in) :: ps

    integer :: i, l

    ps_bound_niwfs = 0
    do i = 1, ps%conf%p
      l = ps%conf%l(i)
      if (any(.not. ps%bound(i,:))) cycle
      ps_bound_niwfs = ps_bound_niwfs + (2*l+1)
    end do

  end function ps_bound_niwfs

  !---------------------------------------

  pure logical function ps_has_density(ps) result(has_density)
    type(ps_t), intent(in) :: ps

    has_density = ps%has_density

  end function ps_has_density

  !---------------------------------------

  pure logical function ps_has_nlcc(ps) result(has_nlcc)
    type(ps_t), intent(in) :: ps

    has_nlcc = ps%nlcc

  end function ps_has_nlcc

  !---------------------------------------
  FLOAT function ps_density_volume(ps, namespace) result(volume)
    type(ps_t),        intent(in) :: ps
    type(namespace_t), intent(in) :: namespace

    integer :: ip, ispin
    FLOAT :: rr
    FLOAT, allocatable ::vol(:)
    type(spline_t) :: volspl

    PUSH_SUB(ps_density_volume)

    if (.not. ps_has_density(ps)) then
      message(1) = "The pseudopotential does not contain an atomic density"
      call messages_fatal(1, namespace=namespace)
    end if

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
