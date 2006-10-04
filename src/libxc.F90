!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module lib_xc_m
  implicit none

  private
  public ::                         &
    xc_info_number,                 &
    xc_info_kind,                   &
    xc_info_name,                   &
    xc_info_family,                 &
    xc_info_refs,                   &
    xc_family_from_id,              &
    xc_lda_init,                    &
    xc_lda,                         &
    xc_lda_fxc,                     &
    xc_lda_kxc,                     &
    xc_lda_end,                     &
    xc_lca_init,                    &
    xc_lca,                         &
    xc_gga_init,                    &
    xc_gga,                         &
    xc_gga_end,                     &
    xc_gga_lb,                      &
    xc_mgga_init,                   &
    xc_mgga,                        &
    xc_mgga_end

  ! Families of xc functionals
  integer, public, parameter ::     &
    XC_FAMILY_UNKNOWN       =  -1,  &
    XC_FAMILY_LDA           =   1,  &
    XC_FAMILY_GGA           =   2,  &
    XC_FAMILY_MGGA          =   4,  &
    XC_FAMILY_LCA           =   8,  &
    XC_FAMILY_OEP           =  16

  integer, public, parameter ::     &
    XC_UNPOLARIZED          =   1,  &  ! Spin unpolarized
    XC_POLARIZED            =   2      ! Spin polarized

  integer, public, parameter ::     &
    XC_NON_RELATIVISTIC     =   0,  &  ! Functional includes or not realtivistic
    XC_RELATIVISTIC         =   1      ! corrections. Only available in some functionals.

  ! Kinds
  integer, public, parameter ::     &
    XC_EXCHANGE             =   0,  &
    XC_CORRELATION          =   1,  &
    XC_EXCHANGE_CORRELATION =   2

  ! the LDAs
  integer, public, parameter ::     &
    XC_LDA_X                =   1,  &  ! Exchange
    XC_LDA_C_WIGNER         =   2,  &  ! Wigner parametrization
    XC_LDA_C_RPA            =   3,  &  ! Random Phase Approximation
    XC_LDA_C_HL             =   4,  &  ! Hedin & Lundqvist
    XC_LDA_C_GL             =   5,  &  ! Gunnarson & Lundqvist
    XC_LDA_C_XALPHA         =   6,  &  ! Slaters Xalpha
    XC_LDA_C_VWN            =   7,  &  ! Vosko, Wilk, & Nussair
    XC_LDA_C_VWN_RPA        =   8,  &  ! Vosko, Wilk, & Nussair (RPA)
    XC_LDA_C_PZ             =   9,  &  ! Perdew & Zunger
    XC_LDA_C_PZ_MOD         =  10,  &  ! Perdew & Zunger (Modified)
    XC_LDA_C_OB_PZ          =  11,  &  ! Ortiz & Ballone (PZ)
    XC_LDA_C_PW             =  12,  &  ! Perdew & Wang
    XC_LDA_C_OB_PW          =  13,  &  ! Ortiz & Ballone (PW)
    XC_LDA_C_AMGB           =  14      ! Attacalite et al

  ! the GGAs
  integer, public, parameter ::     &
    XC_GGA_X_PBE            = 101,  &  ! Perdew, Burke & Ernzerhof exchange
    XC_GGA_X_PBE_R          = 102,  &  ! Perdew, Burke & Ernzerhof exchange (revised)
    XC_GGA_X_B86            = 103,  &  ! Becke 86 Xalpha,beta,gamma
    XC_GGA_X_B86_R          = 104,  &  ! Becke 86 Xalpha,beta,gamma reoptimized
    XC_GGA_X_B86_MGC        = 105,  &  ! Becke 88 Xalfa,beta,gamma (with mod. grad. correction)
    XC_GGA_X_B88            = 106,  &  ! Becke 88
    XC_GGA_X_G96            = 107,  &  ! Gill 96
    XC_GGA_C_PBE            = 130,  &  ! Perdew, Burke & Ernzerhof correlation
    XC_GGA_C_LYP            = 131,  &  ! Lee, Yang & Parr
    XC_GGA_XC_LB            = 160      ! van Leeuwen & Baerends


  ! the meta-GGAs
  integer, public, parameter ::     &
    XC_MGGA_X_TPSS          = 201,  &  ! Perdew, Tao, Staroverov & Scuseria exchange
    XC_MGGA_C_TPSS          = 202      ! Perdew, Tao, Staroverov & Scuseria correlation

  ! the LCAs
  integer, public, parameter ::     &
    XC_LCA_OMC              = 301,  &  ! Orestes, Marcasso & Capelle
    XC_LCA_LCH              = 302      ! Lee, Colwell & Handy

  ! info
  interface
    integer function xc_info_number(info)
      integer(POINTER_SIZE), intent(in) :: info
    end function xc_info_number

    integer function xc_info_kind(info)
      integer(POINTER_SIZE), intent(in) :: info
    end function xc_info_kind

    subroutine xc_info_name(info, s)
      integer(POINTER_SIZE), intent(in)  :: info
      character(len=*),      intent(out) :: s
    end subroutine xc_info_name

    integer function xc_info_family(info)
      integer(POINTER_SIZE), intent(in)  :: info
    end function xc_info_family

    subroutine xc_info_refs(info, n, s)
      integer(POINTER_SIZE), intent(in)    :: info
      integer(POINTER_SIZE), intent(inout) :: n
      character(len=*),      intent(out)   :: s
    end subroutine xc_info_refs
  end interface


  ! functionals
  interface
    integer function xc_family_from_id(id)
      integer, intent(in) :: id
    end function xc_family_from_id
  end interface


  ! We will use the same public interface (xc_lda_init) for the three C procedures
  interface xc_lda_init
    subroutine xc_lda_init_(p, info, functional, nspin)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_lda_init_

    subroutine xc_lda_x_init(p, info, functional, nspin, dim, irel)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,               intent(in)  :: dim    ! 2 or 3 dimensions
      integer,               intent(in)  :: irel   ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_lda_x_init

    subroutine xc_lda_c_xalpha_init(p, info, functional, nspin, dim, alpha)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,               intent(in)  :: dim    ! 2 or 3 dimensions
      FLOAT,                 intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_lda_c_xalpha_init
  end interface

  interface
    subroutine xc_lda_end(p)
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine xc_lda_end

    subroutine xc_lda(p, rho, e, v)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(out) :: e     ! the energy per unit particle
      FLOAT,                 intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_lda

    subroutine xc_lda_fxc(p, rho, fxc)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(out) :: fxc   ! v(nspin,nspin) the xc kernel
    end subroutine xc_lda_fxc

    subroutine xc_lda_kxc(p, rho, kxc)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(out) :: kxc
    end subroutine xc_lda_kxc

  end interface


  ! We will use the same public procedure for the two C procedures.
  interface xc_gga_init
    subroutine xc_gga_init_(p, info, functional, nspin)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_gga_init_

    subroutine xc_gga_lb_init(p, info, functional, nspin, modified, threshold)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
      integer,               intent(in)  :: modified
      FLOAT,                 intent(in)  :: threshold
    end subroutine xc_gga_lb_init
  end interface

  interface
    subroutine xc_gga_end(p)
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine xc_gga_end

    subroutine xc_gga(p, rho, grho, e, dedd, dedgd)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      FLOAT,                 intent(out) :: e     ! the energy per unit particle
      FLOAT,                 intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      FLOAT,                 intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_gga

    subroutine xc_gga_lb(p, rho, grho, r, ip, qtot, dedd)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      FLOAT,                 intent(in)  :: r     ! distance from center of finite system
      FLOAT,                 intent(in)  :: ip    ! ionization potential
      FLOAT,                 intent(in)  :: qtot  ! total charge
      FLOAT,                 intent(out) :: dedd
    end subroutine xc_gga_lb
  end interface


  ! the meta-GGAs
  interface
    subroutine xc_mgga_init(p, info, functional, nspin)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_mgga_init

    subroutine xc_mgga_end(p)
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine xc_mgga_end

    subroutine xc_mgga(p, rho, grho, tau, e, dedd, dedgd, dedtau)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(in)  :: grho  ! grho(3,nspin) the gradient of the density
      FLOAT,                 intent(in)  :: tau   ! tau(nspin) the kinetic energy density
      FLOAT,                 intent(out) :: e     ! the energy per unit particle
      FLOAT,                 intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      FLOAT,                 intent(out) :: dedgd ! in terms of the gradient of the density
      FLOAT,                 intent(out) :: dedtau! and in terms of tau
    end subroutine xc_mgga
  end interface

  ! the LCAs
  interface
    subroutine xc_lca_init(p, info, functional, nspin)
      integer(POINTER_SIZE), intent(out) :: p
      integer(POINTER_SIZE), intent(out) :: info
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_lca_init

    subroutine xc_lca(p, rho, v, e, dedd, dedv)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(in)  :: v     ! v(3,nspin) the vorticity
      FLOAT,                 intent(out) :: e     ! the energy per unit particle
      FLOAT,                 intent(out) :: dedd  ! dedd(nspin) the derivative of the energy
      ! in terms of the density
      FLOAT,                 intent(out) :: dedv  ! and in terms of the vorticity
    end subroutine xc_lca
  end interface

end module lib_xc_m
