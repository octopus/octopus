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
  use functions_m
  use global_m
  use grid_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use lib_xc_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use profiling_m
  use states_m
  use xc_functl_m

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
    type(xc_functl_t) :: functl(2,2) ! (1,:) => exchange,    (2,:) => correlation
                                        ! (:,1) => unpolarized, (:,2) => polarized

    type(xc_functl_t) :: kernel(2,2)

    type(xc_functl_t) :: j_functl    ! current-depent part of the functional

    ! the meta-GGA can be implemented in two ways
    integer :: mGGA_implementation      ! 1 => as a GGA like functional
                                        ! 2 => using the OEP method
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

    call pop_sub()
  end subroutine xc_write_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, ndim, spin_channels, cdft)
    type(xc_t), intent(out) :: xcs
    integer,    intent(in)  :: ndim
    integer,    intent(in)  :: spin_channels
    logical,    intent(in)  :: cdft

    integer :: i, x_id, c_id, xk_id, ck_id

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
        call xc_functl_init_correlation(xcs%functl(2,i), c_id, ndim, i)
        
        call xc_functl_init_exchange   (xcs%kernel(1,i), xk_id, ndim, i)
        call xc_functl_init_correlation(xcs%kernel(2,i), ck_id, ndim, i)

      end do

      xcs%family = ior(xcs%family, xcs%functl(1,1)%family)
      xcs%family = ior(xcs%family, xcs%functl(2,1)%family)

      xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(1,1)%family)
      xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(2,1)%family)

      if (iand(xcs%family, XC_FAMILY_LCA).ne.0 .and. &
        iand(xcs%family, XC_FAMILY_MGGA + XC_FAMILY_OEP).ne.0) then
        message(1) = "LCA functional can only be used along with LDA or GGA functionals"
        call write_fatal(1)
      end if

      if(iand(xcs%family, XC_FAMILY_MGGA).ne.0) then
        call loct_parse_int(check_inp('MGGAimplementation'), 1, xcs%mGGA_implementation)
        if(xcs%mGGA_implementation.ne.1.and.xcs%mGGA_implementation.ne.2) then
          message(1) = 'MGGAimplementation can only assume the values:'
          message(2) = '  1 : GEA implementation'
          message(3) = '  2 : OEP implementation'
          call write_fatal(3)
        end if
      end if

    end if

    call pop_sub()
    
  contains 
    
    subroutine parse()
      integer :: val, default

      ! the first 3 digits of the number indicate the X functional and
      ! the next 3 the C functional.

      default = XC_LDA_X

      select case(calc_dim)
      case(3); default = default + XC_LDA_C_PZ_MOD*1000
      case(2); default = default + XC_LDA_C_AMGB*1000
      end select

      !%Variable XCFunctional
      !%Type integer
      !%Default lda_x+lda_c_pz_mod
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the exchange and correlation functionals
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
      !% GGA: dePristo & Kress 87 (version R1)
      !%Option gga_xc_dk87_r2 112
      !% GGA: dePristo & Kress 87 (version R2)
      !%Option gga_xc_ft97_a 114
      !% GGA: Filatov & Thiel 97 (version A)
      !%Option gga_xc_ft97_b 115
      !% GGA: Filatov & Thiel 97 (version B)
      !%Option gga_xc_lb 160
      !% GGA: van Leeuwen & Baerends (GGA)
      !%Option mgga_x_tpss 201
      !% MGGA (not working)
      !%Option oep_x 401
      !% OEP: Exact exchange
      !%Option lda_c_wigner 2000
      !% LDA: Wigner parametrization
      !%Option lda_c_rpa 3000
      !% LDA: Random Phase Approximation
      !%Option lda_c_hl 4000
      !% LDA: Hedin & Lundqvist
      !%Option lda_c_gl 5000
      !% LDA: Gunnarson & Lundqvist
      !%Option lda_c_xalpha 6000
      !% LDA: Slater s Xalpha
      !%Option lda_c_vwn 7000
      !% LDA: Vosko, Wilk, & Nussair
      !%Option lda_c_vwn_rpa 8000
      !% LDA: Vosko, Wilk, & Nussair (fit to the RPA correlation energy)
      !%Option lda_c_pz 9000
      !% LDA: Perdew & Zunger
      !%Option lda_c_pz_mod 10000
      !% LDA: Perdew & Zunger (Modified to improve the matching between
      !% the high and the low rs region)
      !%Option lda_c_ob_pz 11000
      !% LDA: Ortiz & Ballone (PZ-type parametrization)
      !%Option lda_c_pw 12000
      !% LDA: Perdew & Wang
      !%Option lda_c_pw_mod 13000
      !% LDA: Perdew & Wang (Modified to match the original PBE routine)
      !%Option lda_c_ob_pw 14000
      !% LDA: Ortiz & Ballone (PW-type parametrization)
      !%Option lda_c_amgb 15000
      !% LDA: Attacalite et al functional for the 2D electron gas
      !%Option gga_c_pbe 130000
      !% GGA: Perdew, Burke & Ernzerhof correlation
      !%Option lda_c_lyp 131000
      !% GGA: Lee, Yang, & Parr LDA
      !%Option lda_c_p86 132000
      !% GGA: Perdew 86
      !%Option mgga_c_tpss 202000
      !% MGGA (not working)
      !%End

      call loct_parse_int(check_inp('XCFunctional'), default, val)

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

      call loct_parse_int(check_inp('XCKernel'), -1, val)

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


#include "xc_vxc.F90"
#include "xc_axc.F90"
#include "xc_fxc.F90"
#include "xc_kxc.F90"

end module xc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
