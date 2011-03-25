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
  use lalg_basic_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use parser_m
  use profiling_m
  use states_m
  use states_dim_m
  use symmetrizer_m
  use varinfo_m
  use XC_F90(lib_m)
  use xc_functl_m

  implicit none

  private
  public ::             &
    xc_t,               &
    xc_init,            &
    xc_end,             &
    xc_messages_info,      &
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
    type(xc_functl_t) :: j_functl       ! current-dependent part of the functional

    FLOAT   :: exx_coef                 ! amount of EXX to add for the hybrids
    integer :: mGGA_implementation      ! how to implement the MGGAs
  end type xc_t


  FLOAT, parameter :: tiny      = CNST(1.0e-12)
  FLOAT, parameter :: denom_eps = CNST(1.0e-20) ! added to denominators to avoid overflows...

contains

  ! ---------------------------------------------------------
  subroutine xc_messages_info(xcs, iunit)
    type(xc_t), intent(in) :: xcs
    integer,    intent(in) :: iunit

    integer :: isp

    PUSH_SUB(xc_messages_info)

    if (xcs%cdft .and. iand(xcs%family, XC_FAMILY_LCA) /= 0) then
      write(message(1), '(a)') "Current-dependent exchange-correlation:"
      call messages_info(1, iunit)
      call xc_functl_messages_info(xcs%j_functl, iunit)

      write(message(1), '(1x)')
      write(message(2), '(a)') "Auxiliary exchange-correlation functionals:"
      call messages_info(2, iunit)
    else
      write(message(1), '(a)') "Exchange-correlation:"
      call messages_info(1, iunit)
    end if

    do isp = 1, 2
      call xc_functl_messages_info(xcs%functl(isp, 1), iunit)
    end do

    if(xcs%exx_coef.ne.M_ZERO) then
      write(message(1), '(1x)')
      write(message(2), '(a,f8.5)') "Exact exchange mixing = ", xcs%exx_coef
      call messages_info(2, iunit)
    end if

    POP_SUB(xc_messages_info)
  end subroutine xc_messages_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, ndim, nel, spin_channels, cdft, hartree_fock)
    type(xc_t), intent(out) :: xcs
    integer,    intent(in)  :: ndim
    FLOAT,      intent(in)  :: nel
    integer,    intent(in)  :: spin_channels
    logical,    intent(in)  :: cdft
    logical,    intent(in)  :: hartree_fock

    integer :: isp, x_id, c_id, xk_id, ck_id
    logical :: ll

    PUSH_SUB(xc_init)

    xcs%cdft   = cdft  ! make a copy of flag indicating the use of current-dft_m

    ! get current-dependent functional
    call xc_j_functl_init (xcs%j_functl, cdft, spin_channels)
    xcs%family = xcs%j_functl%family
    xcs%kernel_family = 0

    if (xcs%family == XC_FAMILY_LCA .or. xcs%family == 0) then
      
      call parse()

      !we also need XC functionals that do not depend on the current
      !get both spin-polarized and unpolarized
      do isp = 1, 2

        call xc_functl_init_functl(xcs%functl(1,isp),  x_id, ndim, nel, isp)
        call xc_functl_init_functl(xcs%functl(2,isp),  c_id, ndim, nel, isp)
        
        call xc_functl_init_functl(xcs%kernel(1,isp), xk_id, ndim, nel, isp)
        call xc_functl_init_functl(xcs%kernel(2,isp), ck_id, ndim, nel, isp)

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
          message(1) = "You cannot use an exchange functional when performing"
          message(2) = "a Hartree-Fock calculation or using a hybrid functional."
          call messages_fatal(2)
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
        message(1) = "LCA functional can only be used along with LDA or GGA functionals."
        call messages_fatal(1)
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
        !% to <math>\phi^*(r)</math>. This is the approach used in most quantum-chemistry
        !% (and other) programs.
        !%Option mgga_gea 2
        !% Use gradient expansion (GEA) of the kinetic-energy density.
        !%Option mgga_oep 3
        !% Use the OEP equation to obtain the XC potential. This is the "correct" way
        !% to do it within DFT.
        !%End
        call parse_integer(datasets_check('MGGAimplementation'), 1, xcs%mGGA_implementation)
        if(.not.varinfo_valid_option('MGGAimplementation', xcs%mGGA_implementation)) &
          call input_error('xcs%mGGA_implementation')
      end if

    end if

    POP_SUB(xc_init)
    
  contains 
    
    subroutine parse()
      integer :: val, default

      PUSH_SUB(xc_init.parse)

      ! the first 3 digits of the number indicate the X functional and
      ! the next 3 the C functional.

      default = 0
      if(.not.hartree_fock) then
        select case(ndim)
        case(3); default = XC_LDA_X    + XC_LDA_C_PZ_MOD *1000
        case(2); default = XC_LDA_X_2D + XC_LDA_C_2D_AMGB*1000
        case(1); default = XC_LDA_X_1D + XC_LDA_C_1D_CSC *1000
        end select
      end if

      ! The description of this variable can be found in file src/xc/functionals_list.F90
      call parse_integer(datasets_check('XCFunctional'), default, val)

      c_id = val / 1000
      x_id = val - c_id*1000

      call messages_obsolete_variable('XFunctional', 'XCFunctional')
      call messages_obsolete_variable('CFunctional', 'XCFunctional')

      !%Variable XCKernel
      !%Type integer
      !%Default lda_x+lda_c_pz_mod
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the exchange-correlation kernel. Only LDA kernels are available currently.
      !%Option xc_functional -1
      !% The same functional defined by <tt>XCFunctional</tt>.
      !%End

      call parse_integer(datasets_check('XCKernel'), default, val)

      if( -1 == val ) then
        ck_id = c_id
        xk_id = x_id
      else
        ck_id = val / 1000
        xk_id = val - ck_id*1000  
      end if

      POP_SUB(xc_init.parse)
    end subroutine parse

  end subroutine xc_init


  ! ---------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_t), intent(inout) :: xcs

    integer :: isp

    PUSH_SUB(xc_end)

    if (xcs%cdft) then
      call xc_functl_end(xcs%j_functl)
    end if
    do isp = 1, 2
      call xc_functl_end(xcs%functl(1,isp))
      call xc_functl_end(xcs%functl(2,isp))
      call xc_functl_end(xcs%kernel(1,isp))
      call xc_functl_end(xcs%kernel(2,isp))
    end do
    xcs%family = 0

    POP_SUB(xc_end)
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
