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
!!
!! $Id$

#include "global.h"

module xc_functl
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use lib_xc

  implicit none

  private
  public :: xc_functl_type,             &
            xc_j_functl_init,           &
            xc_functl_init_exchange,    &
            xc_functl_init_correlation, &
            xc_functl_end,              &
            xc_functl_write_info


  ! This adds to the constants defined in lib_xc. But since in that module
  ! the OEP functionals are not included, it is better to put it here.
  integer, public, parameter :: &
     XC_OEP_X             = 401     ! Exact exchange

  type xc_functl_type
    integer :: family              ! LDA, GGA, etc.
    integer :: id                  ! identifier

    integer :: spin_channels       ! XC_UNPOLARIZED | XC_POLARIZED

    integer(POINTER_SIZE) :: conf  ! the pointer used to call the library
    integer(POINTER_SIZE) :: info  ! information about the functional
  end type xc_functl_type

contains

  ! -----------------------------------------------------------
  subroutine xc_functl_init(functl, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: spin_channels

    functl%family = 0
    functl%id     = 0
    functl%spin_channels = spin_channels

  end subroutine xc_functl_init

  ! -----------------------------------------------------------
  subroutine xc_j_functl_init(functl, cdft, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    logical,              intent(in)  :: cdft
    integer,              intent(in)  :: spin_channels

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    if (.not.cdft) return

    ! read input
    call loct_parse_int(check_inp('JFunctional'), XC_LCA_OMC, functl%id)

    ! initialize
    select case(functl%id)
    case(0)

    case(XC_LCA_OMC, XC_LCA_LCH)
      functl%family = XC_FAMILY_LCA
      call xc_lca_init(functl%conf, functl%info, functl%id, &
         spin_channels)

    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
         "' is not a known current functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select

  end subroutine xc_j_functl_init


  ! -----------------------------------------------------------
  subroutine xc_functl_init_exchange(functl, ndim, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: ndim
    integer,              intent(in)  :: spin_channels

    integer :: j
    FLOAT :: alpha

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    ! read input
    call loct_parse_int(check_inp('XFunctional'), XC_LDA_X, functl%id)

    ! initialize
    select case(functl%id)
    case(0)

    case(XC_LDA_X)
      functl%family = XC_FAMILY_LDA
      call xc_lda_init(functl%conf, functl%info, XC_LDA_X, &
         spin_channels, ndim)

    case(XC_GGA_X_PBE, XC_GGA_XC_LB)
      functl%family = XC_FAMILY_GGA

      if(functl%id == XC_GGA_XC_LB) then
        call loct_parse_int  (check_inp('LB94_modified'), 0, j)
        call loct_parse_float(check_inp('LB94_threshold'), CNST(1.0e-6), alpha)
        call xc_gga_init(functl%conf, functl%info, functl%id, &
           spin_channels, j, alpha)
      else
        call xc_gga_init(functl%conf, functl%info, functl%id, spin_channels)
      end if

    case(XC_MGGA_X_TPSS)
      functl%family = XC_FAMILY_MGGA
      call xc_mgga_init(functl%conf, functl%info, functl%id, spin_channels)

    case(XC_OEP_X)
      functl%family = XC_FAMILY_OEP

    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
         "' is not a known exchange functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select

  end subroutine xc_functl_init_exchange


  ! -----------------------------------------------------------
  subroutine xc_functl_init_correlation(functl, ndim, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: ndim
    integer,              intent(in)  :: spin_channels

    FLOAT :: alpha

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    ! read input
    call loct_parse_int(check_inp('CFunctional'), XC_LDA_C_PZ, functl%id)

    ! initialize
    select case(functl%id)
    case(0)

    case(XC_LDA_C_WIGNER, XC_LDA_C_RPA, XC_LDA_C_HL, XC_LDA_C_GL, XC_LDA_C_XALPHA, &
       XC_LDA_C_VWN, XC_LDA_C_PZ, XC_LDA_C_OB_PZ, XC_LDA_C_PW, XC_LDA_C_OB_PW,     &
       XC_LDA_C_LYP, XC_LDA_C_AMGB)

      if(functl%id==XC_LDA_C_AMGB.and.ndim.ne.2) then
        message(1) = 'Functional AMGB only allowed in 2D'
        call write_fatal(1)
      end if

      functl%family = XC_FAMILY_LDA

      if(functl%id.ne.XC_LDA_C_XALPHA) then
        call xc_lda_init(functl%conf, functl%info, functl%id, spin_channels)
      else
        call loct_parse_float(check_inp('Xalpha'), M_ONE, alpha)
        call xc_lda_init(functl%conf, functl%info, XC_LDA_C_XALPHA, &
           spin_channels, ndim, alpha)
      end if

    case(XC_GGA_C_PBE)
      functl%family = XC_FAMILY_GGA
      call xc_gga_init(functl%conf, functl%info, functl%id, spin_channels)

    case(XC_MGGA_C_TPSS)
      functl%family = XC_FAMILY_MGGA
      call xc_mgga_init(functl%conf, functl%info, functl%id, spin_channels)

    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
         "' is not a known correlation functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select

  end subroutine xc_functl_init_correlation


  ! -----------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_type), intent(inout) :: functl

    select case(functl%family)
    case(XC_FAMILY_LDA);  call xc_lda_end (functl%conf)
    case(XC_FAMILY_GGA);  call xc_gga_end (functl%conf)
    case(XC_FAMILY_MGGA); call xc_mgga_end(functl%conf)
    case(XC_FAMILY_LCA);  call xc_lca_end (functl%conf)
    end select

  end subroutine xc_functl_end


  ! -----------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_type), intent(in) :: functl
    integer,              intent(in) :: iunit

    character(len=120) :: s1, s2
    integer(POINTER_SIZE) :: i, j

    select case (functl%family)
    case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_LCA)
      ! we hapilly call the xc library

      i = xc_info_kind(functl%info)
      select case(i)
      case(int(XC_EXCHANGE, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Exchange'
      case(int(XC_CORRELATION, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Correlation'
      case(int(XC_EXCHANGE_CORRELATION, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Exchange-correlation'
      end select

      call xc_info_name  (functl%info, s1)
      call xc_info_family(functl%info, s2)
      write(iunit, '(4x,4a)') trim(s1), ' (', trim(s2), ')'

      i = 0; j = 1
      call xc_info_ref(functl%info, i, s1)
      do while(i>=0)
        write(iunit, '(4x,a,i1,2a)') '[', j, '] ', trim(s1)
        call xc_info_ref(functl%info, i, s1)
        j = j + 1
      end do

    case (XC_FAMILY_OEP)
      ! this is handled separately

      select case(functl%id)
      case(XC_OEP_X)
        write(iunit, '(2x,a)') 'Exchange'
      end select

      select case(functl%id)
      case(XC_OEP_X)
        s1 = 'Exact exchange'
      end select
      write(iunit, '(4x,2a)') trim(s1), ' (OEP)'
    end select

  end subroutine xc_functl_write_info


end module xc_functl
