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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module xc_functl_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use lib_xc_m

  implicit none

  private
  public ::                     &
    xc_functl_t,                &
    xc_j_functl_init,           &
    xc_functl_init_exchange,    &
    xc_functl_init_correlation, &
    xc_functl_end,              &
    xc_functl_write_info


  ! This adds to the constants defined in lib_xc. But since in that module
  ! the OEP functionals are not included, it is better to put it here.
  integer, public, parameter :: &
    XC_OEP_X             = 401     ! Exact exchange

  type xc_functl_t
    integer   :: family            ! LDA, GGA, etc.
    integer   :: id                ! identifier

    integer   :: spin_channels     ! XC_UNPOLARIZED | XC_POLARIZED

    C_POINTER :: conf              ! the pointer used to call the library
    C_POINTER :: info              ! information about the functional
  end type xc_functl_t

contains

  ! ---------------------------------------------------------
  subroutine xc_functl_init(functl, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: spin_channels

    functl%family = 0
    functl%id     = 0
    functl%spin_channels = spin_channels

  end subroutine xc_functl_init


  ! ---------------------------------------------------------
  subroutine xc_j_functl_init(functl, cdft, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    logical,           intent(in)  :: cdft
    integer,           intent(in)  :: spin_channels

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    if (.not.cdft) return

    !%Variable JFunctional
    !%Type integer
    !%Default lca_omc
    !%Section Hamiltonian::XC
    !%Description
    !% Defines the current functional
    !%Option lca_omc 301
    !% Orestes, Marcasso & Capelle 
    !%Option lca_lch 302
    !% Lee, Colwell & Handy
    !%End
    call loct_parse_int(check_inp('JFunctional'), XC_LCA_OMC, functl%id)

    ! initialize
    select case(functl%id)
    case(0)

    case(XC_LCA_OMC, XC_LCA_LCH)
      functl%family = XC_FAMILY_LCA
      call xc_f90_lca_init(functl%conf, functl%info, functl%id, &
        spin_channels)

    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
        "' is not a known current functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select

  end subroutine xc_j_functl_init


  ! ---------------------------------------------------------
  subroutine xc_functl_init_exchange(functl, ndim, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: spin_channels

    integer :: j
    FLOAT :: alpha

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    !%Variable XFunctional
    !%Type integer
    !%Default lda_x
    !%Section Hamiltonian::XC
    !%Description
    !% Defines the exchange functional
    !%Option lda_x 1
    !% LDA: Slater exchange
    !%Option gga_x_pbe 101
    !% GGA: Perdew, Burke & Ernzerhof
    !%Option gga_x_pbe_r 102
    !% GGA: Perdew, Burke & Ernzerhof (revised)
    !%Option gga_x_b86 103
    !% GGA: Becke 86 Xalpha,beta,gamma
    !%Option gga_x_b86_r 104
    !% GGA: Becke 86 Xalpha,beta,gamma reoptimized
    !%Option gga_x_b86_mgc 105
    !% GGA: Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
    !%Option gga_x_b88 106
    !% GGA: Becke 88
    !%Option gga_x_g96 107
    !% GGA: Gill 96
    !%Option gga_x_pw86 108
    !% GGA: Perdew & Wang 86
    !%Option gga_x_pw91 109
    !% GGA: Perdew & Wang 91
    !%Option gga_x_optx 110
    !% GGA: Handy & Cohen OPTX 01
    !%Option gga_xc_dk87_r1 111
    !% GGA: dePristo & Kress 87 version R1
    !%Option gga_xc_dk87_r2 112
    !% GGA: dePristo & Kress 87 version R2
    !%Option gga_xc_lg93 113
    !% GGA: Lacks & Gordon 93
    !%Option gga_xc_lb 160
    !% GGA: van Leeuwen & Baerends (GGA)
    !%Option mgga_x_tpss 201
    !% MGGA (not working)
    !%Option oep_x 401
    !% OEP: Exact exchange
    !%End
    call loct_parse_int(check_inp('XFunctional'), XC_LDA_X, functl%id)

    if(functl%id.ne.0) then
      ! get the family of the functional
      functl%family = xc_f90_family_from_id(functl%id)

      if(functl%family == XC_FAMILY_UNKNOWN) then
        if(functl%id == XC_OEP_X) then
          functl%family = XC_FAMILY_OEP
        else
          call input_error('XFunctional')
        end if
      end if
    end if

    ! initialize
    select case(functl%family)
    case(XC_FAMILY_LDA)
      call xc_f90_lda_init(functl%conf, functl%info, XC_LDA_X, &
         spin_channels, ndim, XC_NON_RELATIVISTIC)

    case(XC_FAMILY_GGA)
      if(functl%id == XC_GGA_XC_LB) then
        call loct_parse_int  (check_inp('LB94_modified'), 0, j)
        call loct_parse_float(check_inp('LB94_threshold'), CNST(1.0e-6), alpha)
        call xc_f90_gga_init(functl%conf, functl%info, functl%id, &
          spin_channels, j, alpha)
      else
        call xc_f90_gga_init(functl%conf, functl%info, functl%id, spin_channels)
      end if

    case(XC_FAMILY_MGGA)
      call xc_f90_mgga_init(functl%conf, functl%info, functl%id, spin_channels)

    end select

  end subroutine xc_functl_init_exchange


  ! ---------------------------------------------------------
  subroutine xc_functl_init_correlation(functl, ndim, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: spin_channels

    FLOAT :: alpha

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    !%Variable CFunctional
    !%Type integer
    !%Default lda_c_pz_mod
    !%Section Hamiltonian::XC
    !%Description
    !% Defines the correlation functional
    !%Option lda_c_wigner 2
    !% LDA: Wigner parametrization
    !%Option lda_c_rpa 3
    !% LDA: Random Phase Approximation
    !%Option lda_c_hl 4
    !% LDA: Hedin & Lundqvist
    !%Option lda_c_gl 5
    !% LDA: Gunnarson & Lundqvist
    !%Option lda_c_xalpha 6
    !% LDA: Slater s Xalpha
    !%Option lda_c_vwn 7
    !% LDA: Vosko, Wilk, & Nussair
    !%Option lda_c_vwn_rpa 8
    !% LDA: Vosko, Wilk, & Nussair (fit to the RPA correlation energy)
    !%Option lda_c_pz 9
    !% LDA: Perdew & Zunger
    !%Option lda_c_pz_mod 10
    !% LDA: Perdew & Zunger (Modified to improve the matching between
    !%      the high and the low rs region)
    !%Option lda_c_ob_pz 11
    !% LDA: Ortiz & Ballone (PZ-type parametrization)
    !%Option lda_c_pw 12
    !% LDA: Perdew & Wang
    !%Option lda_c_ob_pw 13
    !% LDA: Ortiz & Ballone (PW-type parametrization)
    !%Option lda_c_amgb 14
    !% LDA: Attacalite et al functional for the 2D electron gas
    !%Option gga_c_pbe 130
    !% GGA: Perdew, Burke & Ernzerhof correlation
    !%Option lda_c_lyp 131
    !% GGA: Lee, Yang, & Parr LDA
    !%Option lda_c_p86 132
    !% GGA: Perdew 86
    !%Option mgga_c_tpss 202
    !% MGGA (not working)
    !%End
    select case(calc_dim)
      case(3); call loct_parse_int(check_inp('CFunctional'), XC_LDA_C_PZ_MOD, functl%id)
      case(2); call loct_parse_int(check_inp('CFunctional'), XC_LDA_C_AMGB, functl%id)
      case(1); call loct_parse_int(check_inp('CFunctional'), 0, functl%id)
    end select

    if(functl%id.ne.0) then
      ! get the family of the functional
      functl%family = xc_f90_family_from_id(functl%id)

      if(functl%family == XC_FAMILY_UNKNOWN) then
        call input_error('CFunctional')
      end if
    end if

    ! initialize
    select case(functl%family)
    case(XC_FAMILY_LDA)

      if(functl%id==XC_LDA_C_AMGB.and.ndim.ne.2) then
        message(1) = 'Functional AMGB only allowed in 2D'
        call write_fatal(1)
      end if

      if(functl%id.ne.XC_LDA_C_XALPHA) then
        call xc_f90_lda_init(functl%conf, functl%info, functl%id, spin_channels)
      else
        call loct_parse_float(check_inp('Xalpha'), M_ONE, alpha)
        call xc_f90_lda_init(functl%conf, functl%info, XC_LDA_C_XALPHA, &
          spin_channels, ndim, alpha)
      end if

    case(XC_FAMILY_GGA)
      call xc_f90_gga_init(functl%conf, functl%info, functl%id, spin_channels)

    case(XC_FAMILY_MGGA)
      call xc_f90_mgga_init(functl%conf, functl%info, functl%id, spin_channels)

    end select

  end subroutine xc_functl_init_correlation


  ! ---------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_t), intent(inout) :: functl

    select case(functl%family)
    case(XC_FAMILY_LDA);  call xc_f90_lda_end (functl%conf)
    case(XC_FAMILY_GGA);  call xc_f90_gga_end (functl%conf)
    case(XC_FAMILY_MGGA); call xc_f90_mgga_end(functl%conf)
    case(XC_FAMILY_LCA);  call xc_f90_lca_end (functl%conf)
    end select

  end subroutine xc_functl_end


  ! ---------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_t), intent(in) :: functl
    integer,           intent(in) :: iunit

    character(len=120) :: s1, s2
    C_POINTER :: str
    integer :: i

    call push_sub('xc_functl.xc_functl_write_info')

    select case (functl%family)
    case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_LCA)
      ! we hapilly call the xc library

      i = xc_f90_info_kind(functl%info)
      select case(i)
      case(XC_EXCHANGE)
        write(message(1), '(2x,a)') 'Exchange'
      case(XC_CORRELATION)
        write(message(1), '(2x,a)') 'Correlation'
      case(XC_EXCHANGE_CORRELATION)
        write(message(1), '(2x,a)') 'Exchange-correlation'
      end select

      call xc_f90_info_name  (functl%info, s1)
      select case(functl%family)
        case (XC_FAMILY_LDA);  write(s2,'(a)') "LDA"
        case (XC_FAMILY_GGA);  write(s2,'(a)') "GGA"
        case (XC_FAMILY_MGGA); write(s2,'(a)') "MGGA"
        case (XC_FAMILY_LCA);  write(s2,'(a)') "LCA"
      end select
      write(message(2), '(4x,4a)') trim(s1), ' (', trim(s2), ')'
      call write_info(2, iunit)
      
      i = 1; str = 0
      call xc_f90_info_ref(functl%info, str, s1)
      do while(str >= 0)
        write(message(1), '(4x,a,i1,2a)') '[', i, '] ', trim(s1)
        call write_info(1, iunit)
        call xc_f90_info_ref(functl%info, str, s1)
        i = i + 1
      end do

    case (XC_FAMILY_OEP)
      ! this is handled separately

      select case(functl%id)
      case(XC_OEP_X)
        write(message(1), '(2x,a)') 'Exchange'
        write(message(2), '(4x,a)') 'Exact exchange (OEP)'
        call write_info(2, iunit)
      end select

    end select

    call pop_sub()
  end subroutine xc_functl_write_info

end module xc_functl_m
