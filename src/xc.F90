#include "config_F90.h"

module xc
use liboct
use math
use mesh
use fft
use hartree
use states

implicit none

private
public :: xc_type, xc_write_info, xc_init, &
     xc_end, dxc_pot, zxc_pot

integer, parameter ::     &
    XC_FAMILY_START= 0,    &
    XC_FAMILY_ZER  = 0,    &
    XC_FAMILY_LDA  = 1,    &
    XC_FAMILY_GGA  = 2,    &
    XC_FAMILY_KLI  = 3,    &
    XC_FAMILY_MGGA = 4,    &
    XC_FAMILY_END  = 4,    &
    X_FUNC_START     = 10, &
    X_FUNC_ZER       = 10, &
    X_FUNC_LDA_REL   = 11, &
    X_FUNC_LDA_NREL  = 12, &
    X_FUNC_GGA_PBE   = 13, &
    X_FUNC_GGA_PBER  = 14, &
    X_FUNC_GGA_LB94  = 15, &
    X_FUNC_MGGA_PKZB = 16, &
    X_FUNC_KLI_X     = 17, &
    X_FUNC_KLI_SIC   = 18, &
    X_FUNC_END       = 18, &
    C_FUNC_START     = 48, &
    C_FUNC_ZER       = 48, &
    C_FUNC_LDA_RPA   = 49, &
    C_FUNC_LDA_PZ    = 50, &
    C_FUNC_LDA_PW92  = 51, &
    C_FUNC_GGA_PBE   = 52, &
    C_FUNC_KLI_SIC   = 53, &
    C_FUNC_MGGA_PKZB = 54, &
    C_FUNC_END       = 54

character(len=4), parameter :: name_xc(XC_FAMILY_END-XC_FAMILY_START+1) = (/ &
    'None', &
    'LDA ',  &
    'GGA ',  &
    'KLI ',  &
    'MGGA' /)
character(len=18), parameter :: name_x(X_FUNC_END-X_FUNC_START+1) = (/ &
    'none              ', &
    'relativistic      ', &
    'non-relativistic  ', &
    'PBE               ', &
    'PBE - relativistic', &
    'LB94              ', &
    'PKZB              ', &
    'exact exchange    ', &
    'SIC (LDA)         ' /)
character(len=14), parameter :: name_c(C_FUNC_END-C_FUNC_START+1) = (/ &
    'none          ', &
    'RPA           ', &
    'Perdew-Zunger ', &
    'Perdew-Wang 92', &
    'PBE           ', &
    'SIC (LDA)     ', &
    'PKZB          ' /)

type xc_type
  integer :: x_family, x_func, c_family, c_func
end type xc_type

real(r8), parameter :: small = 1e-5_r8
real(r8), parameter :: denom_eps = 1e-20_r8 ! added to denominators to avoid overflows...

contains

subroutine xc_write_info(xcs, iunit)
  type(xc_type), intent(IN) :: xcs
  integer, intent(in) :: iunit

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    write(iunit, '(6x,a,a)') 'Exchange    family    : ', &
         name_xc(xcs%x_family-XC_FAMILY_START+1)
    write(iunit, '(6x,a,a)') '            functional: ', &
         name_x(xcs%x_func-X_FUNC_START+1)
    write(iunit, '(6x,a,a)') 'Correlation family    : ', &
         name_xc(xcs%c_family-XC_FAMILY_START+1)
    write(iunit, '(6x,a,a)') '            functional: ', &
         name_c(xcs%c_func-C_FUNC_START+1)

#ifdef HAVE_MPI
  end if
#endif
  return
end subroutine xc_write_info

subroutine xc_init(xcs, m, ispin)
  type(xc_type), intent(out) :: xcs
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: ispin

  character(len=5) :: xfam, cfam, xfunc, cfunc

  sub_name = 'xc_init'; call push_sub()

  call oct_parse_str('XFamily',     'LDA',  xfam)
  call oct_parse_str('XFunctional', 'NREL', xfunc)
  call oct_parse_str('CFamily',     'LDA',  cfam)
  call oct_parse_str('CFunctional', 'PZ',   cfunc)

  ! the exchange
  select case(trim(xfam))
  case('ZER')
    xcs%x_family = XC_FAMILY_ZER
    xcs%x_func   = X_FUNC_ZER
  case('LDA')
    xcs%x_family = XC_FAMILY_LDA
    select case(trim(xfunc))
    case('NREL')
      xcs%x_func = X_FUNC_LDA_NREL
    case('REL')
      xcs%x_func = X_FUNC_LDA_REL
    case default
      write(message(1), '(a,a,a)') "'", trim(xfam), &
          "' is not a known exchange LDA functional!"
      message(2) = "(XFunc = NREL | REL)"
      call write_fatal(2)
    end select
  case('GGA')
    xcs%x_family = XC_FAMILY_GGA
    select case(trim(xfunc))
    case('PBE')
      xcs%x_func = X_FUNC_GGA_PBE
    case('PBER')
      xcs%x_func = X_FUNC_GGA_PBER
    case('LB94')
      xcs%x_func = X_FUNC_GGA_LB94
    case default
      write(message(1), '(a,a,a)') "'", trim(xfam), &
          "' is not a known exchange GGA functional!"
      message(2) = "(XFunc = PBE | LB94)"
      call write_fatal(2)
    end select
  case('MGGA')
    xcs%x_family = XC_FAMILY_MGGA
    select case(trim(xfunc))
    case('PKZB')
      xcs%x_func = X_FUNC_MGGA_PKZB
    case default
      write(message(1), '(a,a,a)') "'", trim(xfam), &
          "' is not a known exchange MGGA functional!"
      message(2) = "(XFunc = PKZB)"
      call write_fatal(2)
    end select
#ifdef HAVE_LAPACK
  case('KLI')
#if defined(HAVE_MPI) && defined(MPI_TD)
    message(1) = "KLI is not allowed with MPI_TD!"
    call write_fatal(1)
#endif
    xcs%x_family = XC_FAMILY_KLI
    select case(trim(xfunc))
    case('EXX')
      xcs%x_func = X_FUNC_KLI_X
    case('SIC')
      xcs%x_func = X_FUNC_KLI_SIC
    case default
      write(message(1), '(a,a,a)') "'", trim(xfunc), &
          "' is not a known exchange KLI functional!"
      message(2) = "(XFunc = X | SIC)"
      call write_fatal(2)
    end select
#endif
  case default
    write(message(1), '(a,a,a)') "'", trim(xfam), &
        "' is not a known exchange functional family!"
    message(2) = "(XFamily = ZER | LDA | GGA | KLI)"
    call write_fatal(2)
  end select

  ! now the correlation
  select case(trim(cfam))
  case('ZER')
    xcs%c_family = XC_FAMILY_ZER
    xcs%c_func   = C_FUNC_ZER
  case('LDA')
    xcs%c_family = XC_FAMILY_LDA
    select case(trim(cfunc))
    case('PZ')
      xcs%c_func = C_FUNC_LDA_PZ
    case('PW92')
      xcs%c_func = C_FUNC_LDA_PW92
    case default
      write(message(1), '(a,a,a)') "'", trim(xfam), &
          "' is not a known correlation LDA functional!"
      message(2) = "(CFunc = PZ | PW92)"
      call write_fatal(2)
    end select
  case('GGA')
    xcs%c_family = XC_FAMILY_GGA
    select case(trim(cfunc))
    case('PBE')
      xcs%c_func = C_FUNC_GGA_PBE
    case default
      write(message(1), '(a,a,a)') "'", trim(cfam), &
          "' is not a known correlation GGA functional!"
      message(2) = "(CFunc = PBE)"
      call write_fatal(2)
    end select
#ifdef HAVE_LAPACK
  case('MGGA')
    xcs%c_family = XC_FAMILY_MGGA
    select case(trim(cfunc))
    case('PKZB')
      xcs%c_func = C_FUNC_MGGA_PKZB
    case default
      write(message(1), '(a,a,a)') "'", trim(cfam), &
          "' is not a known exchange MGGA functional!"
      message(2) = "(CFunc = PKZB)"
      call write_fatal(2)
    end select
  case('KLI')
    xcs%c_family = XC_FAMILY_KLI
    select case(trim(xfunc))
    case('SIC')
      xcs%c_func = C_FUNC_KLI_SIC
    case default
      write(message(1), '(a,a,a)') "'", trim(cfunc), &
          "' is not a known correlation KLI functional!"
      message(2) = "(XFunc = SIC)"
      call write_fatal(2)
    end select
#endif
  case default
    write(message(1), '(a,a,a)') "'", trim(cfam), &
        "' is not a known correlation functional family!"
    message(2) = "(CFamily = ZER | LDA | GGA | KLI)"
    call write_fatal(2)
  end select

  ! can only do non-colinear spin with LDA
  if(ispin == 4 .and. &
       (.not.(xcs%x_family == XC_FAMILY_ZER.or.xcs%x_family == XC_FAMILY_LDA) .or. &
       .not.(xcs%c_family == XC_FAMILY_ZER.or.xcs%c_family == XC_FAMILY_LDA))) then
    message(1) = "Can only handle non-colinear spin within the LDA!"
  end if

  call pop_sub()
  return
end subroutine xc_init

subroutine xc_end(xcs)
  type(xc_type), intent(inout) :: xcs

  return
end subroutine xc_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A couple of auxiliary functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getSpinFactor(nspin, socc, sfact)
  integer, intent(in) :: nspin
  real(r8), intent(out) :: socc, sfact

  select case(nspin)
  case(1) ! we need to correct for the spin occupancies
     socc  = 0.5_r8
     sfact = 2.0_r8
  case(2)
     socc  = 1.0_r8
     sfact = 1.0_r8
  case default
     write(6,'(a,I2)') 'KLI: error cannot handle nspin=', nspin
  end select
end subroutine getSpinFactor

function my_sign(a)
  real(r8), intent(in) ::  a
  real(r8) :: my_sign

  if(a < 0) then
     my_sign = -1.0_r8
  else
     my_sign = 1.0_r8
  end if
  return
end function my_sign

! include the xc potentials

#include "xc_LDA.F90"
#include "xc_GGA.F90"
!#include "xc_MGGA.F90"

#include "undef.F90"
#include "real.F90"
#include "xc_pot.F90"
#include "xc_KLI.F90"
#include "xc_KLI_x.F90"
#include "xc_KLI_SIC.F90"
#include "undef.F90"
#include "complex.F90"
#include "xc_pot.F90"
#include "xc_KLI.F90"
#include "xc_KLI_x.F90"
#include "xc_KLI_SIC.F90"

end module xc
