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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module xc_m
  use cube_m
  use cube_function_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use index_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use lalg_adv_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use par_vec_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use states_dim_m
  use symmetrizer_m
  use unit_system_m
  use varinfo_m
  use XC_F90(lib_m)
  use xc_functl_m

  implicit none

  private
  public ::             &
    xc_t,               &
    xc_init,            &
    xc_end,             &
    xc_write_info,      &
    xc_get_vxc,         &
    xc_get_vxc_cmplx,   &
    xc_get_fxc,         &
    xc_get_kxc,         &
    xc_is_orbital_dependent


  type xc_t
    integer :: family                   !< the families present
    integer :: kernel_family
    type(xc_functl_t) :: functl(2,2)    !< (FUNC_X,:) => exchange,    (FUNC_C,:) => correlation
                                        !< (:,1) => unpolarized, (:,2) => polarized

    type(xc_functl_t) :: kernel(2,2)
    FLOAT   :: kernel_lrc_alpha         !< long-range correction alpha parameter for kernel in solids

    FLOAT   :: exx_coef                 !< amount of EXX to add for the hybrids
    integer :: mGGA_implementation      !< how to implement the MGGAs
    logical :: use_gi_ked               !< should we use the gauge-independent kinetic energy density?

    integer :: xc_density_correction
    logical :: xcd_optimize_cutoff
    FLOAT   :: xcd_ncutoff
    logical :: xcd_minimum
    logical :: xcd_normalize
  end type xc_t

  FLOAT, parameter :: tiny      = CNST(1.0e-12)
  FLOAT, parameter :: denom_eps = CNST(1.0e-20) !< added to denominators to avoid overflows...

  integer, parameter :: &
    LR_NONE = 0,        &
    LR_X    = 1

  integer, parameter :: &
    FUNC_X = 1,         &
    FUNC_C = 2

contains

  ! ---------------------------------------------------------
  subroutine xc_write_info(xcs, iunit)
    type(xc_t), intent(in) :: xcs
    integer,    intent(in) :: iunit

    integer :: ifunc

    PUSH_SUB(xc_write_info)

    write(message(1), '(a)') "Exchange-correlation:"
    call messages_info(1, iunit)

    do ifunc = FUNC_X, FUNC_C
      call xc_functl_write_info(xcs%functl(ifunc, 1), iunit)
    end do
    
    if(xcs%exx_coef /= M_ZERO) then
      write(message(1), '(1x)')
      write(message(2), '(a,f8.5)') "Exact exchange mixing = ", xcs%exx_coef
      call messages_info(2, iunit)
    end if

    POP_SUB(xc_write_info)
  end subroutine xc_write_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, ndim, nel, hartree_fock)
    type(xc_t), intent(out) :: xcs
    integer,    intent(in)  :: ndim
    FLOAT,      intent(in)  :: nel
    logical,    intent(in)  :: hartree_fock

    integer :: isp, x_id, c_id, xk_id, ck_id
    logical :: ll

    PUSH_SUB(xc_init)

    xcs%family = 0
    xcs%kernel_family = 0

    call parse()

    !we also need XC functionals that do not depend on the current
    !get both spin-polarized and unpolarized
    do isp = 1, 2

      call xc_functl_init_functl(xcs%functl(FUNC_X, isp),  x_id, ndim, nel, isp)
      call xc_functl_init_functl(xcs%functl(FUNC_C, isp),  c_id, ndim, nel, isp)

      call xc_functl_init_functl(xcs%kernel(FUNC_X, isp), xk_id, ndim, nel, isp)
      call xc_functl_init_functl(xcs%kernel(FUNC_C, isp), ck_id, ndim, nel, isp)

    end do

    xcs%family = ior(xcs%family, xcs%functl(FUNC_X,1)%family)
    xcs%family = ior(xcs%family, xcs%functl(FUNC_C,1)%family)

    xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(FUNC_X,1)%family)
    xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(FUNC_C,1)%family)

    ! Take care of hybrid functionals (they appear in the correlation functional)
    xcs%exx_coef = M_ZERO
    ll =  (hartree_fock) &
      .or.(xcs%functl(FUNC_X,1)%id == XC_OEP_X) &
      .or.(iand(xcs%functl(FUNC_C,1)%family, XC_FAMILY_HYB_GGA) /= 0)
    if(ll) then
      if((xcs%functl(FUNC_X,1)%id /= 0).and.(xcs%functl(FUNC_X,1)%id /= XC_OEP_X)) then
        message(1) = "You cannot use an exchange functional when performing"
        message(2) = "a Hartree-Fock calculation or using a hybrid functional."
        call messages_fatal(2)
      end if

      ! get the mixing coefficient for hybrids
      if(iand(xcs%functl(FUNC_C,1)%family, XC_FAMILY_HYB_GGA) /= 0) then
        call XC_F90(hyb_exx_coef)(xcs%functl(FUNC_C,1)%conf, xcs%exx_coef)
      else
        ! we are doing Hartree-Fock plus possibly a correlation functional
        xcs%exx_coef = M_ONE
      end if

      ! reset certain variables
      xcs%functl(FUNC_X,1)%family = XC_FAMILY_OEP
      xcs%functl(FUNC_X,1)%id     = XC_OEP_X
      xcs%family             = ior(xcs%family, XC_FAMILY_OEP)
    end if

    ! Now it is time for these current functionals
    if (iand(xcs%family, XC_FAMILY_LCA) /= 0 .and. &
      iand(xcs%family, XC_FAMILY_MGGA + XC_FAMILY_OEP) /= 0) then
      message(1) = "LCA functional can only be used along with LDA or GGA functionals."
      call messages_fatal(1)
    end if

    if(iand(xcs%family, XC_FAMILY_MGGA) /= 0) then
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


      call messages_obsolete_variable('CurrentInTau', 'XCUseGaugeIndependentKED')

      !%Variable XCUseGaugeIndependentKED
      !%Type logical
      !%Default yes
      !%Section Hamiltonian::XC
      !%Description
      !% If true, when evaluating the XC functional, a term including the (paramagnetic or total) current
      !% is added to the kinetic-energy density such as to make it gauge-independent.
      !%End
      call parse_logical(datasets_check('XCUseGaugeIndependentKED'), .true., xcs%use_gi_ked)

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

      !%Variable XCKernelLRCAlpha
      !%Type float
      !%Default 0.0
      !%Section Hamiltonian::XC
      !%Description
      !% Set to a non-zero value to add a long-range correction for solids to the kernel.
      !% This is the alpha parameter defined in S. Botti <i>et al.</i>, <i>Phys. Rev. B</i>
      !% 69, 155112 (2004), which results in multiplying the Hartree term by
      !% <math>1 - \alpha / 4 \pi</math>. 
      !%End

      call parse_float(datasets_check('XCKernelLRCAlpha'), M_ZERO, xcs%kernel_lrc_alpha)
      if(abs(xcs%kernel_lrc_alpha) > M_EPSILON) &
        call messages_experimental("Long-range correction to kernel")

      !%Variable XCDensityCorrection
      !%Type integer
      !%Default none
      !%Section Hamiltonian::XC
      !%Description
      !% This variable controls the long range correction of the XC
      !% potential using the XC density representation
      !% (http://arxiv.org/abs/1107.4339). By default, no correction
      !% is applied.
      !%Option none 0
      !% No correction is applied.
      !%Option long_range_x 1
      !% The correction is applied to the exchange potential.
      !%End
      call parse_integer(datasets_check('XCDensityCorrection'), LR_NONE, xcs%xc_density_correction)

      if(xcs%xc_density_correction /= LR_NONE) then 
        call messages_experimental('XC density correction')

        !%Variable XCDensityCorrectionOptimize
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC
        !%Description
        !% When enabled, the default, the density cutoff will be
        !% optimized to replicate the boundary conditions of the exact
        !% XC potential. If the variable is set to no, the value of
        !% the cutoff must be given by the XCDensityCorrectionCutoff
        !% variable.
        !%End
        call parse_logical(datasets_check('XCDensityCorrectionOptimize'), .true., xcs%xcd_optimize_cutoff)

        !%Variable XCDensityCorrectionCutoff
        !%Type float
        !%Default 0.0
        !%Section Hamiltonian::XC
        !%Description
        !% The value of the cutoff applied to the XC density. The default value is 0.
        !%End
        call parse_float(datasets_check('XCDensityCorrectionCutoff'), CNST(0.0), xcs%xcd_ncutoff)

        !%Variable XCDensityCorrectionMinimum
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC
        !%Description
        !% When enabled, the default, the cutoff optimization will
        !% return the first minimum of the q_xc function if it does
        !% not find a value of -1 (See http://arxiv.org/abs/1107.4339
        !% for details). This is required for atoms or small
        !% molecules, but may cause numerical problems.
        !%End
        call parse_logical(datasets_check('XCDensityCorrectionMinimum'), .true., xcs%xcd_minimum)

        !%Variable XCDensityCorrectionNormalize
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC
        !%Description
        !% When enabled, the default, the correction will be
        !% normalized to reproduce the exact boundary conditions of
        !% the XC potential.
        !%End
        call parse_logical(datasets_check('XCDensityCorrectionNormalize'), .true., xcs%xcd_normalize)
  
      end if


      POP_SUB(xc_init.parse)
    end subroutine parse

  end subroutine xc_init


  ! ---------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_t), intent(inout) :: xcs

    integer :: isp

    PUSH_SUB(xc_end)

    do isp = 1, 2
      call xc_functl_end(xcs%functl(FUNC_X, isp))
      call xc_functl_end(xcs%functl(FUNC_C, isp))
      call xc_functl_end(xcs%kernel(FUNC_X, isp))
      call xc_functl_end(xcs%kernel(FUNC_C, isp))
    end do
    xcs%family = 0

    POP_SUB(xc_end)
  end subroutine xc_end

  ! ---------------------------------------------------------
  logical function xc_is_orbital_dependent(xcs)
    type(xc_t), intent(in) :: xcs

    PUSH_SUB(xc_is_orbital_dependent)

    xc_is_orbital_dependent = xcs%exx_coef /= M_ZERO .or. &
      iand(xcs%functl(FUNC_X,1)%family, XC_FAMILY_OEP) /= 0 .or. &
      iand(xcs%family, XC_FAMILY_MGGA) /= 0

    POP_SUB(xc_is_orbital_dependent)
  end function xc_is_orbital_dependent


#include "vxc_inc.F90"
#include "fxc_inc.F90"
#include "kxc_inc.F90"

end module xc_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
