!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module xc
use poisson
use states
use vxc

implicit none

private
public :: xc_type, xc_write_info, xc_init, &
     xc_end, dxc_pot, zxc_pot

!!! This parameters have to be updated when implementing
!!! new functionals
integer, parameter ::     &
     N_XC_FAMILIES = 4,   &
     N_X_FUNCTL    = 9,   &
     N_C_FUNCTL    = 7

!!! Families of xc functionals
integer, parameter ::     &
     XC_FAMILY_LDA  = ibset(0, 0),    &
     XC_FAMILY_GGA  = ibset(0, 1),    &
     XC_FAMILY_KLI  = ibset(0, 2),    &
     XC_FAMILY_MGGA = ibset(0, 3)

!!! Functionals
integer, parameter ::     &
     X_FUNC_LDA_REL   = ibset(0,  0), &
     X_FUNC_LDA_NREL  = ibset(0,  1), &
     X_FUNC_GGA_PBE   = ibset(0,  2), &
     X_FUNC_GGA_PBER  = ibset(0,  3), &
     X_FUNC_GGA_LB94  = ibset(0,  4), &
     X_FUNC_MGGA_PKZB = ibset(0,  5), &
     X_FUNC_KLI_X     = ibset(0,  6), &
     X_FUNC_KLI_SIC   = ibset(0,  7), &
     X_FUNC_KLI_HJU   = ibset(0,  8), &
     C_FUNC_LDA_RPA   = ibset(0,  9), &
     C_FUNC_LDA_PZ    = ibset(0, 10), &
     C_FUNC_LDA_PW92  = ibset(0, 11), &
     C_FUNC_GGA_PBE   = ibset(0, 12), &
     C_FUNC_KLI_SIC   = ibset(0, 13), &
     C_FUNC_KLI_HJU   = ibset(0, 14), &
     C_FUNC_MGGA_PKZB = ibset(0, 15)


character(len=4), parameter :: xc_name_families(N_XC_FAMILIES) =   &
     (/ 'LDA ', 'GGA ', 'KLI ', 'MGGA' /)

!!! Exchange functionals should come before correlation functionals
character(len=20), parameter :: xc_name_functionals(N_X_FUNCTL+N_C_FUNCTL) = (/ &
     'LDA  - relat.       ', &  ! X_FUNC_LDA_REL
     'LDA  - non-relat.   ', &  ! X_FUNC_LDA_NREL
     'GGA  - PBE          ', &  ! X_FUNC_GGA_PBE
     'GGA  - PBE - relat. ', &  ! X_FUNC_GGA_PBER
     'GGA  - LB94         ', &  ! X_FUNC_GGA_LB94
     'MGGA - PKZB         ', &  ! X_FUNC_MGGA_PKZB
     'KLI  - EXX          ', &  ! X_FUNC_KLI_X
     'KLI  - SIC (LDA)    ', &  ! X_FUNC_KLI_SIC
     'KLI  - SIC (HJU)    ', &  ! X_FUNC_KLI_HJU
     'LDA  - RPA          ', &  ! C_FUNC_LDA_RPA
     'LDA  - PZ81         ', &  ! C_FUNC_LDA_PZ
     'LDA  - PW92         ', &  ! C_FUNC_LDA_PW92
     'GGA  - PBE          ', &  ! C_FUNC_GGA_PBE
     'KLI  - SIC (LDA)    ', &  ! C_FUNC_KLI_SIC
     'KLI  - SIC (HJU)    ', &  ! C_FUNC_KLI_HJU
     'MGGA - PKZB         '  &  ! C_FUNC_MGGA_PKZB
     /)

type xc_type
  integer :: family, functl
  logical :: nlcc ! repeated from system

  ! For LB94 fine-tuning
  logical  :: lb94_modified
  real(r8) :: lb94_beta
  real(r8) :: lb94_threshold
end type xc_type

real(r8), parameter :: small     = 1e-5_r8
real(r8), parameter :: tiny      = 1.0e-12_r8
real(r8), parameter :: denom_eps = 1e-20_r8 ! added to denominators to avoid overflows...

contains

subroutine xc_write_info(xcs, iunit)
  type(xc_type), intent(IN) :: xcs
  integer, intent(in) :: iunit

  integer :: i

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    do i = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(btest(xcs%functl, i)) then
        if(i < N_X_FUNCTL) then
          write(iunit, '(6x,a,a)') 'Exchange    : ', xc_name_functionals(i+1)
        else
          write(iunit, '(6x,a,a)') 'Correlation : ', xc_name_functionals(i+1)
        end if
      end if
    end do

#ifdef HAVE_MPI
  end if
#endif
  return
end subroutine xc_write_info

subroutine xc_init(xcs, nlcc)
  type(xc_type), intent(out) :: xcs
  logical, intent(in) :: nlcc

  character(len=50) :: func

  call push_sub('xc_init')

  xcs%nlcc   = nlcc  ! make a copy of falg indicating non-local core corrections
  xcs%family = 0     ! initialize xc family and functional to zero
  xcs%functl = 0

  ! set the appropriate flags
  call oct_parse_string('XFunctional', 'LDA', func)
  select case(func(1:4))
  case('ZER ')
  case('LDA ')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, X_FUNC_LDA_NREL)
  case('RLDA')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, X_FUNC_LDA_REL)
  case('PBE')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_PBE)
  case('RPBE')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_PBER)
  case('LB94')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_LB94)
  case('PKZB')
    xcs%family = ior(xcs%family, XC_FAMILY_MGGA)
    xcs%functl = ior(xcs%functl, X_FUNC_MGGA_PKZB)
  case('EXX ')
    xcs%family = ior(xcs%family, XC_FAMILY_KLI)
    xcs%functl = ior(xcs%functl, X_FUNC_KLI_X)
  case('SIC ')
    xcs%family = ior(xcs%family, XC_FAMILY_KLI)
    xcs%functl = ior(xcs%functl, X_FUNC_KLI_SIC)
  case('HJU ')
    xcs%family = ior(xcs%family, XC_FAMILY_KLI)
    xcs%functl = ior(xcs%functl, X_FUNC_KLI_HJU)
  case default
    write(message(1), '(a,a,a)') "'", trim(func), "' is not a known exchange functional!"
    message(2) = "Please check the manual for a list of possible values."
    call write_fatal(2)
  end select

  ! now the correlation
  call oct_parse_string('CFunctional', 'PZ81', func)
  select case(func(1:4))
  case('ZER ')
  case('PZ81')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, C_FUNC_LDA_PZ)
  case('PW92')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, C_FUNC_LDA_PW92)
  case('PBE ')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, C_FUNC_GGA_PBE)
  case('PKZB')
    xcs%family = ior(xcs%family, XC_FAMILY_MGGA)
    xcs%functl = ior(xcs%functl, C_FUNC_MGGA_PKZB)
  case('SIC ')
    xcs%family = ior(xcs%family, XC_FAMILY_KLI)
    xcs%functl = ior(xcs%functl, C_FUNC_KLI_SIC)
  case('HJU ')
    xcs%family = ior(xcs%family, XC_FAMILY_KLI)
    xcs%functl = ior(xcs%functl, C_FUNC_KLI_HJU)
  case default
    write(message(1), '(a,a,a)') "'", trim(func), &
         "' is not a known correlation functional family!"
    message(2) = "(CFamily = ZER | LDA | GGA | KLI)"
    call write_fatal(2)
  end select

#if defined(HAVE_MPI) && defined(MPI_TD)
  if(btest(xcs%family, XC_FAMILY_KLI)) then
    message(1) = "KLI is not allowed with MPI_TD!"
    call write_fatal(1)
  end if
#endif

  ! Extra parameters to fine-tune LB94 potential.
  if(iand(xcs%functl, X_FUNC_GGA_LB94).ne.0) then
    call oct_parse_double("LB94_beta", 0.05_r8, xcs%lb94_beta)
    call oct_parse_double("LB94_threshold", 1.0e-6_r8, xcs%lb94_threshold)
    call oct_parse_logical("LB94_modified", .false., xcs%lb94_modified)
  end if

  call pop_sub()
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
  case(2, 4)
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
!#include "xc_KLI.F90"
!#include "xc_KLI_x.F90"
!#include "xc_KLI_SIC.F90"
!#include "xc_HJU.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_pot.F90"
!#include "xc_KLI.F90"
!#include "xc_KLI_x.F90"
!#include "xc_KLI_SIC.F90"
!#include "xc_HJU.F90"

#include "undef.F90"

end module xc
