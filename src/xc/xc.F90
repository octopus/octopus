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

module xc_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use io_function_m
  use lalg_basic_m
  use loct_parser_m
  use XC_F90(lib_m)
  use mesh_function_m
  use mesh_m
  use messages_m
  use profiling_m
  use states_m
  use states_dim_m
  use xc_functl_m
  use varinfo_m

  implicit none

  private
  public ::             &
    xc_t,               &
    xc_init,            &
    xc_end,             &
    xc_write_info,      &
    xc_get_vxc,         &
    xc_get_vxc_and_axc, &
    xc_get_fxc,         &
    xc_get_kxc


  type xc_t
    logical :: cdft

    integer :: family                   ! the families present
    integer :: kernel_family
    type(xc_functl_t) :: functl(2,2)    ! (1,:) => exchange,    (2,:) => correlation
                                        ! (:,1) => unpolarized, (:,2) => polarized

    type(xc_functl_t) :: kernel(2,2)
    type(xc_functl_t) :: j_functl       ! current-depent part of the functional

    FLOAT   :: exx_coef                 ! amount of EXX to add for the hybrids
    integer :: mGGA_implementation      ! how to implement the MGGAs
  end type xc_t


  FLOAT, parameter :: tiny      = CNST(1.0e-12)
  FLOAT, parameter :: denom_eps = CNST(1.0e-20) ! added to denominators to avoid overflows...

contains

  ! ---------------------------------------------------------
  subroutine xc_write_info(xcs, iunit)
    type(xc_t), intent(in) :: xcs
    integer,    intent(in) :: iunit

    integer :: i

    call push_sub('xc.xc_write_info')

    if (xcs%cdft .and. iand(xcs%family, XC_FAMILY_LCA) /= 0) then
      write(message(1), '(a)') "Current-dependent exchange and correlation:"
      call write_info(1, iunit)
      call xc_functl_write_info(xcs%j_functl, iunit)

      write(message(1), '(1x)')
      write(message(2), '(a)') "Auxiliary exchange and correlation functionals:"
      call write_info(2, iunit)
    else
      write(message(1), '(a)') "Exchange and correlation:"
      call write_info(1, iunit)
    end if

    do i = 1, 2
      call xc_functl_write_info(xcs%functl(i, 1), iunit)
    end do

    if(xcs%exx_coef.ne.M_ZERO) then
      write(message(1), '(1x)')
      write(message(2), '(a,f8.5)') "Exact exchange mixing = ", xcs%exx_coef
      call write_info(2, iunit)
    end if

    call pop_sub()
  end subroutine xc_write_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, ndim, nel, spin_channels, cdft, hartree_fock)
    type(xc_t), intent(out) :: xcs
    integer,    intent(in)  :: ndim
    FLOAT,      intent(in)  :: nel
    integer,    intent(in)  :: spin_channels
    logical,    intent(in)  :: cdft
    logical,    intent(in)  :: hartree_fock

    integer :: i, x_id, c_id, xk_id, ck_id
    logical :: ll

    call push_sub('xc.xc_init')

    xcs%cdft   = cdft  ! make a copy of flag indicating the use of current-dft_m

    ! get current-dependent functional
    call xc_j_functl_init (xcs%j_functl, cdft, spin_channels)
    xcs%family = xcs%j_functl%family
    xcs%kernel_family = 0

    if (xcs%family == XC_FAMILY_LCA .or. xcs%family == 0) then
      
      call parse()

      !we also need xc functionals that do not depend on the current
      !get both spin polarized and unpolarized
      do i = 1, 2

        call xc_functl_init_exchange   (xcs%functl(1,i), x_id, ndim, i)
        call xc_functl_init_correlation(xcs%functl(2,i), c_id, ndim, nel, i)
        
        call xc_functl_init_exchange   (xcs%kernel(1,i), xk_id, ndim, i)
        call xc_functl_init_correlation(xcs%kernel(2,i), ck_id, ndim, nel, i)

      end do

      xcs%family = ior(xcs%family, xcs%functl(1,1)%family)
      xcs%family = ior(xcs%family, xcs%functl(2,1)%family)

      xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(1,1)%family)
      xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(2,1)%family)

      ! Take care of hybrid functionals (they appear in the correlation functional)
      xcs%exx_coef = M_ZERO
      ll =  (hartree_fock) &
        .or.(xcs%functl(1,1)%id.eq.XC_OEP_X) &
        .or.(iand(xcs%functl(2,1)%family, XC_FAMILY_HYB_GGA).ne.0)
      if(ll) then
        if((xcs%functl(1,1)%id.ne.0).and.(xcs%functl(1,1)%id.ne.XC_OEP_X)) then
          message(1) = "You can not use an exchange functional when performing"
          message(2) = "an Hartree-Fock calculation or using a hybrid functional"
          call write_fatal(2)
        end if

        ! get the mixing coefficient for hybrids
        if(iand(xcs%functl(2,1)%family, XC_FAMILY_HYB_GGA).ne.0) then
          call XC_F90(hyb_gga_exx_coef)(xcs%functl(2,1)%conf, xcs%exx_coef)
        else
          ! we are doing Hartree-Fock plus possibly a correlation functional
          xcs%exx_coef = M_ONE
        end if

        ! reset certain variables
        xcs%functl(1,1)%family = XC_FAMILY_OEP
        xcs%functl(1,1)%id     = XC_OEP_X
        xcs%family             = ior(xcs%family, XC_FAMILY_OEP)
      end if

      ! Now it is time for these current functionals
      if (iand(xcs%family, XC_FAMILY_LCA).ne.0 .and. &
        iand(xcs%family, XC_FAMILY_MGGA + XC_FAMILY_OEP).ne.0) then
        message(1) = "LCA functional can only be used along with LDA or GGA functionals"
        call write_fatal(1)
      end if

      if(iand(xcs%family, XC_FAMILY_MGGA).ne.0) then
        !%Variable MGGAimplementation
        !%Type integer
        !%Default mgga_gea
        !%Section Hamiltonian::XC
        !%Description
        !% Decides how to implement the meta-GGAs (NOT WORKING).
        !%Option mgga_dphi 1
        !% Use for <math>v_xc</math> the derivative of the energy functional with respect
        !% to <math>\phi^*(r)</math>. This is the approach used in most quantum chemistry
        !% (and not only) programs
        !%Option mgga_gea 2
        !% Use gradient expansion (GEA) of the kinetic energy density
        !%Option mgga_oep 3
        !% Use the OEP equation to obtain the xc potential. This is the "correct" way
        !% to do it within DFT
        !%End
        call loct_parse_int(datasets_check('MGGAimplementation'), 1, xcs%mGGA_implementation)
        if(.not.varinfo_valid_option('MGGAimplementation', xcs%mGGA_implementation)) &
          call input_error('xcs%mGGA_implementation')
      end if

    end if

    call pop_sub()
    
  contains 
    
    subroutine parse()
      integer :: val, default

      ! the first 3 digits of the number indicate the X functional and
      ! the next 3 the C functional.

      if(hartree_fock) then
        default = 0
      else
        default = XC_LDA_X
        
        select case(calc_dim)
        case(3); default = XC_LDA_X    + XC_LDA_C_PZ_MOD *1000
        case(2); default = XC_LDA_X_2D + XC_LDA_C_2D_AMGB*1000
        end select
      end if

      !%Variable XCFunctional
      !%Type integer
      !%Default lda_x+lda_c_pz_mod
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the exchange and correlation functional to be used,
      !% they should be specified as a sum of a correlation and an
      !% exchange term.
      !%
      !% The value by default is lda_x + lda_c_pz_mod.
      !%
      !%Option lda_x                      1
      !% LDA: Slater exchange
      !%Option lda_x_2d                   19
      !% LDA: Slater exchange
      !%Option lda_c_wigner               2000
      !% LDA: Wigner parametrization
      !%Option lda_c_rpa                  3000
      !% LDA: Random Phase Approximation
      !%Option lda_c_hl                   4000
      !% LDA: Hedin & Lundqvist
      !%Option lda_c_gl                   5000
      !% LDA: Gunnarson & Lundqvist
      !%Option lda_c_xalpha               6000
      !% LDA: Slater s Xalpha
      !%Option lda_c_vwn                  7000
      !% LDA: Vosko, Wilk, & Nussair
      !%Option lda_c_vwn_rpa              8000
      !% LDA: Vosko, Wilk, & Nussair (fit to the RPA correlation energy)
      !%Option lda_c_pz                   9000
      !% LDA: Perdew & Zunger
      !%Option lda_c_pz_mod              10000
      !% LDA: Perdew & Zunger (Modified to improve the matching between
      !% the high and the low rs region)
      !%Option lda_c_ob_pz               11000
      !% LDA: Ortiz & Ballone (PZ-type parametrization)
      !%Option lda_c_pw                  12000
      !% LDA: Perdew & Wang
      !%Option lda_c_pw_mod              13000
      !% LDA: Perdew & Wang (Modified to match the original PBE routine)
      !%Option lda_c_ob_pw               14000
      !% LDA: Ortiz & Ballone (PW-type parametrization)
      !%Option lda_c_2d_amgb             15000
      !% LDA: Attacalite et al functional for the 2D electron gas
      !%Option lda_c_2d_prm08            16000
      !% LDA: Pittalis, Rasanen & Marques correlation in 2D
      !%Option lda_c_vBH                 17000
      !% LDA: von Barth & Hedin
      !%Option lda_xc_teter93            20
      !% LDA: Teter 93 pade parametrization
      !%Option gga_x_pbe                101
      !% GGA: Perdew, Burke & Ernzerhof
      !%Option gga_x_pbe_r              102
      !% GGA: Perdew, Burke & Ernzerhof (revised)
      !%Option gga_x_b86                103
      !% GGA: Becke 86 Xalpha,beta,gamma
      !%Option gga_x_b86_r              104
      !% GGA: Becke 86 Xalpha,beta,gamma reoptimized
      !%Option gga_x_b86_mgc            105
      !% GGA: Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
      !%Option gga_x_b88                106
      !% GGA: Becke 88
      !%Option gga_x_g96                107
      !% GGA: Gill 96
      !%Option gga_x_pw86               108
      !% GGA: Perdew & Wang 86
      !%Option gga_x_pw91               109
      !% GGA: Perdew & Wang 91
      !%Option gga_x_optx               110
      !% GGA: Handy & Cohen OPTX 01
      !%Option gga_x_dk87_r1            111
      !% GGA: dePristo & Kress 87 (version R1)
      !%Option gga_x_dk87_r2            112
      !% GGA: dePristo & Kress 87 (version R2)
      !%Option gga_x_lg93               113
      !% GGA: Lacks & Gordon 93
      !%Option gga_x_ft97_a             114
      !% GGA: Filatov & Thiel 97 (version A)
      !%Option gga_x_ft97_b             115
      !% GGA: Filatov & Thiel 97 (version B)
      !%Option gga_x_pbe_sol            116
      !% GGA: Perdew, Burke & Ernzerhof (solids)
      !%Option gga_x_rpbe               117
      !% GGA: Hammer, Hansen & Norskov (PBE-like)
      !%Option gga_x_wc                 118
      !% GGA:  Wu & Cohen
      !%Option gga_x_mpw91              119
      !% GGA: Modified form of PW91 by Adamo & Barone
      !%Option gga_x_am05               120
      !% GGA: Armiento & Mattsson 05 exchange
      !%Option gga_x_pbea               121
      !% GGA: Madsen 07 (PBE-like)
      !%Option gga_x_mpbe               122
      !% GGA: Adamo & Barone modification to PBE
      !%Option gga_x_xpbe               123
      !% GGA: xPBE reparametrization by Xu & Goddard
      !%Option gga_x_2d_b86_mgc         124
      !% GGA: Becke 86 MGC for 2D systems
      !%Option gga_c_pbe                130000
      !% GGA: Perdew, Burke & Ernzerhof correlation
      !%Option gga_c_lyp                131000
      !% GGA: Lee, Yang, & Parr
      !%Option gga_c_p86                132000
      !% GGA: Perdew 86
      !%Option gga_c_pbe_sol            133000
      !% GGA: Perdew, Burke & Ernzerhof correlation (solids)
      !%Option gga_c_pw91               134000
      !% GGA: Perdew & Wang 91
      !%Option gga_c_am05               135000
      !% GGA: Armiento & Mattsson 05 correlation
      !%Option gga_c_xpbe               136000
      !% GGA: xPBE reparametrization by Xu & Goddard
      !%Option gga_c_lm                 137000
      !% GGA: Langreth and Mehl correlation
      !%Option gga_xc_lb                160
      !% GGA: van Leeuwen & Baerends (xc)
      !%Option gga_xc_hcth_93           161
      !% GGA: HCTH fit to 93 molecules (xc)
      !%Option gga_xc_hcth_120          162
      !% GGA: HCTH fit to 120 molecules (xc)
      !%Option gga_xc_hcth_147          163
      !% GGA: HCTH fit to 147 molecules (xc)
      !%Option gga_xc_hcth_407          164
      !% GGA: HCTH fit to 407 molecules (xc)
      !%Option gga_xc_edf1              165
      !% GGA: Empirical functional from Adamson, Gill, & Pople
      !%Option gga_xc_xlyp              166
      !% GGA: Empirical functional from Xu & Goddard
      !%Option hyb_gga_xc_b3pw91        401000
      !% Hybrid (GGA): the first hybrid by Becke
      !%Option hyb_gga_xc_b3lyp         402000
      !% Hybrid (GGA): The (in)famous B3LYP
      !%Option hyb_gga_xc_b3p86         403000
      !% Hybrid (GGA): Perdew 86 hybrid similar to B3PW91
      !%Option hyb_gga_xc_o3lyp         404000
      !% Hybrid (GGA): Hybrid using the optx functional
      !%Option hyb_gga_xc_pbeh          406000
      !% Hybrid (GGA): aka PBE0 or PBE1PBE
      !%Option hyb_gga_xc_x3lyp         411000
      !% Hybrid (GGA): maybe the best hybrid around
      !%Option hyb_gga_xc_b1wc          412000
      !% Hybrid (GGA): Becke 1-parameter mixture of WC and EXX
      !%Option mgga_x_lta               201
      !% MGGA: Local tau approximation of Ernzerhof & Scuseria
      !%Option mgga_x_tpss              202
      !% MGGA: Perdew, Tao, Staroverov & Scuseria exchange
      !%Option mgga_x_m06l              203
      !% MGGA: Zhao, Truhlar exchange
      !%Option mgga_x_gvt4              204
      !% MGGA: GVT4 from Van Voorhis and Scuseria (exchange part of VSxc)
      !%Option mgga_x_tau_hcth          205
      !% MGGA: tau-HCTH from Boese and Handy
      !%Option mgga_c_tpss              231000
      !% MGGA: Perdew, Tao, Staroverov & Scuseria correlation
      !%Option mgga_c_vsxc              232000
      !% MGGA: VSxc from Van Voorhis and Scuseria (correlation part)
      !%Option oep_x                    901
      !% OEP: Exact exchange
      !%End

      call loct_parse_int(datasets_check('XCFunctional'), default, val)

      c_id = val / 1000
      x_id = val - c_id*1000

      call obsolete_variable('XFunctional', 'XCFunctional')
      call obsolete_variable('CFunctional', 'XCFunctional')

      !%Variable XCKernel
      !%Type integer
      !%Default xc_functional
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the exchange and correlation kernel
      !%Option xc_functional -1
      !% The same functional defined by XCFunctional
      !%End

      call loct_parse_int(datasets_check('XCKernel'), -1, val)

      if( -1 == val ) then
        ck_id = c_id
        xk_id = x_id
      else
        ck_id = val / 1000
        xk_id = val - ck_id*1000  
      end if

    end subroutine parse

  end subroutine xc_init


  ! ---------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_t), intent(inout) :: xcs

    integer :: i

    if (xcs%cdft) then
      call xc_functl_end(xcs%j_functl)
    end if
    do i = 1, 2
      call xc_functl_end(xcs%functl(1,i))
      call xc_functl_end(xcs%functl(2,i))
      call xc_functl_end(xcs%kernel(1,i))
      call xc_functl_end(xcs%kernel(2,i))
    end do
    xcs%family = 0

  end subroutine xc_end


#include "vxc_inc.F90"
#include "axc_inc.F90"
#include "fxc_inc.F90"
#include "kxc_inc.F90"

end module xc_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
