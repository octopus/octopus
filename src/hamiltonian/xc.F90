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

#include "global.h"

module xc_oct_m
  use distributed_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use exchange_operator_oct_m
  use global_oct_m
  use io_oct_m
  use io_function_oct_m
  use iso_c_binding
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lalg_adv_oct_m
  use libvdwxc_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use unit_system_oct_m
  use xc_f03_lib_m
  use xc_functl_oct_m

  implicit none

  private
  public ::             &
    xc_t,               &
    xc_init,            &
    xc_end,             &
    xc_write_info,      &
    xc_get_vxc,         &
    xc_get_fxc,         &
    xc_get_kxc,         &
    xc_is_orbital_dependent, &
    family_is_mgga,     &
    family_is_mgga_with_exc, &
    family_is_hybrid


  type xc_t
    private
    integer,           public :: family              !< the families present
    integer,           public :: flags               !<flags of the xc functional
    integer,           public :: kernel_family
    type(xc_functl_t), public :: functional(2,2)     !< (FUNC_X,:) => exchange,    (FUNC_C,:) => correlation
                                                     !! (:,1) => unpolarized, (:,2) => polarized

    type(xc_functl_t), public :: kernel(2,2)
    FLOAT,             public   :: kernel_lrc_alpha  !< long-range correction alpha parameter for kernel in solids

    FLOAT, public   :: cam_omega                !< Cam coefficients omega, alpha, beta
    FLOAT, public   :: cam_alpha                !< amount of EXX to add for the hybrids
    FLOAT, public   :: cam_beta
    logical                     :: use_gi_ked        !< should we use the gauge-independent kinetic energy density?

    integer :: xc_density_correction
    logical :: xcd_optimize_cutoff
    FLOAT   :: xcd_ncutoff
    logical :: xcd_minimum
    logical :: xcd_normalize
    logical :: parallel

  end type xc_t

  FLOAT, parameter :: tiny      = CNST(1.0e-12)

  integer, parameter :: &
    LR_NONE = 0,        &
    LR_X    = 1

  integer, public, parameter :: &
    FUNC_X = 1,         &
    FUNC_C = 2

contains

  ! ---------------------------------------------------------
  subroutine xc_write_info(xcs, iunit, namespace)
    type(xc_t),        intent(in) :: xcs
    integer,           intent(in) :: iunit
    type(namespace_t), intent(in) :: namespace

    integer :: ifunc

    PUSH_SUB(xc_write_info)

    write(message(1), '(a)') "Exchange-correlation:"
    call messages_info(1, iunit)

    do ifunc = FUNC_X, FUNC_C
      call xc_functl_write_info(xcs%functional(ifunc, 1), iunit, namespace)
    end do
    
    if(xcs%cam_alpha + xcs%cam_beta /= M_ZERO) then
      write(message(1), '(1x)')
      write(message(2), '(a,f8.5)') "Exact exchange mixing = ", xcs%cam_alpha
      if(xcs%cam_beta > M_EPSILON) &
        write(message(2), '(a,f8.5)') "Exact exchange mixing for the short-range part = ", xcs%cam_beta
      call messages_info(2, iunit)
    end if


    POP_SUB(xc_write_info)
  end subroutine xc_write_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, namespace, ndim, periodic_dim, nel, x_id, c_id, xk_id, ck_id, hartree_fock)
    type(xc_t),        intent(out) :: xcs
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: ndim
    integer,           intent(in)  :: periodic_dim
    FLOAT,             intent(in)  :: nel
    integer,           intent(in)  :: x_id
    integer,           intent(in)  :: c_id
    integer,           intent(in)  :: xk_id
    integer,           intent(in)  :: ck_id
    logical,           intent(in)  :: hartree_fock

    integer :: isp
    logical :: ll

    PUSH_SUB(xc_init)

    xcs%family = 0
    xcs%flags  = 0
    xcs%kernel_family = 0

    call parse()

    !we also need XC functionals that do not depend on the current
    !get both spin-polarized and unpolarized
    do isp = 1, 2

      call xc_functl_init(xcs%functional(FUNC_X, isp), namespace, x_id, ndim, nel, isp)
      call xc_functl_init(xcs%functional(FUNC_C, isp), namespace, c_id, ndim, nel, isp)

      call xc_functl_init(xcs%kernel(FUNC_X, isp), namespace, xk_id, ndim, nel, isp)
      call xc_functl_init(xcs%kernel(FUNC_C, isp), namespace, ck_id, ndim, nel, isp)

    end do

    xcs%family = ior(xcs%family, xcs%functional(FUNC_X,1)%family)
    xcs%family = ior(xcs%family, xcs%functional(FUNC_C,1)%family)

    xcs%flags = ior(xcs%flags, xcs%functional(FUNC_X,1)%flags)
    xcs%flags = ior(xcs%flags, xcs%functional(FUNC_C,1)%flags)

    xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(FUNC_X,1)%family)
    xcs%kernel_family = ior(xcs%kernel_family, xcs%kernel(FUNC_C,1)%family)

    ! Take care of hybrid functionals (they appear in the correlation functional)
    xcs%cam_omega = M_ZERO
    xcs%cam_alpha = M_ZERO
    xcs%cam_beta = M_ZERO

    ll =  (hartree_fock) &
      .or.(xcs%functional(FUNC_X,1)%id == XC_OEP_X) &
      .or. family_is_hybrid(xcs)
    if(ll) then
      if((xcs%functional(FUNC_X,1)%id /= 0).and.(xcs%functional(FUNC_X,1)%id /= XC_OEP_X)) then
        message(1) = "You cannot use an exchange functional when performing"
        message(2) = "a Hartree-Fock calculation or using a hybrid functional."
        call messages_fatal(2, namespace=namespace)
      end if

      if(periodic_dim == ndim) &
        call messages_experimental("Fock operator (Hartree-Fock, OEP, hybrids) in fully periodic systems")

      ! get the mixing coefficient for hybrids
      if(family_is_hybrid(xcs)) then
        call xc_f03_hyb_cam_coef(xcs%functional(FUNC_C,1)%conf, xcs%cam_omega, &
                                     xcs%cam_alpha, xcs%cam_beta)

      else
        ! we are doing Hartree-Fock plus possibly a correlation functional
        xcs%cam_omega = M_ZERO
        xcs%cam_alpha = M_ONE
        xcs%cam_beta = M_ZERO
      end if

      ! reset certain variables
      xcs%functional(FUNC_X,1)%family = XC_FAMILY_OEP
      xcs%functional(FUNC_X,1)%id     = XC_OEP_X
      if(.not. hartree_fock) then
        xcs%family             = ior(xcs%family, XC_FAMILY_OEP)
      end if
    end if

    if (bitand(xcs%family, XC_FAMILY_LCA) /= 0) &
      call messages_not_implemented("LCA current functionals", namespace) ! not even in libxc!

    call messages_obsolete_variable(namespace, 'MGGAimplementation')
    call messages_obsolete_variable(namespace, 'CurrentInTau', 'XCUseGaugeIndependentKED')

    if (xcs%functional(FUNC_X, 1)%id == XC_MGGA_X_TB09 .and. periodic_dim /= 3) then
      message(1) = "mgga_x_tb09 functional can only be used for 3D periodic systems"
      call messages_fatal(1, namespace=namespace)
    end if


    if(family_is_mgga(xcs%family)) then
      !%Variable XCUseGaugeIndependentKED
      !%Type logical
      !%Default yes
      !%Section Hamiltonian::XC
      !%Description
      !% If true, when evaluating the XC functional, a term including the (paramagnetic or total) current
      !% is added to the kinetic-energy density such as to make it gauge-independent.
      !% Applies only to meta-GGA (and hybrid meta-GGA) functionals.
      !%End
      call parse_variable(namespace, 'XCUseGaugeIndependentKED', .true., xcs%use_gi_ked)
    end if

    POP_SUB(xc_init)

  contains 

    subroutine parse()

      PUSH_SUB(xc_init.parse)

      ! the values of x_id,  c_id, xk_id, and c_id are read outside the routine
      
      !%Variable XCKernelLRCAlpha
      !%Type float
      !%Default 0.0
      !%Section Hamiltonian::XC
      !%Description
      !% Set to a non-zero value to add a long-range correction for solids to the kernel.
      !% This is the <math>\alpha</math> parameter defined in S. Botti <i>et al.</i>, <i>Phys. Rev. B</i>
      !% 69, 155112 (2004). The <math>\Gamma = \Gamma` = 0</math> term <math>-\alpha/q^2</math> is taken 
      !% into account by introducing an additional pole to the polarizability (see R. Stubner  
      !% <i>et al.</i>, <i>Phys. Rev. B</i> 70, 245119 (2004)). The rest of the terms are included by  
      !% multiplying the Hartree term by <math>1 - \alpha / 4 \pi</math>. The use of non-zero 
      !% <math>\alpha</math> in combination with <tt>HamiltonianVariation</tt> = <tt>V_ext_only</tt>  
      !% corresponds to account of only the <math>\Gamma = \Gamma` = 0</math> term. 
      !% Applicable only to isotropic systems. (Experimental)
      !%End

      call parse_variable(namespace, 'XCKernelLRCAlpha', M_ZERO, xcs%kernel_lrc_alpha)
      if(abs(xcs%kernel_lrc_alpha) > M_EPSILON) &
        call messages_experimental("Long-range correction to kernel")

      !%Variable XCDensityCorrection
      !%Type integer
      !%Default none
      !%Section Hamiltonian::XC::DensityCorrection
      !%Description
      !% This variable controls the long-range correction of the XC
      !% potential using the <a href=http://arxiv.org/abs/1107.4339>XC density representation</a>.
      !%Option none 0
      !% No correction is applied.
      !%Option long_range_x 1
      !% The correction is applied to the exchange potential.
      !%End
      call parse_variable(namespace, 'XCDensityCorrection', LR_NONE, xcs%xc_density_correction)

      if(xcs%xc_density_correction /= LR_NONE) then 
        call messages_experimental('XC density correction')

        !%Variable XCDensityCorrectionOptimize
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC::DensityCorrection
        !%Description
        !% When enabled, the density cutoff will be
        !% optimized to replicate the boundary conditions of the exact
        !% XC potential. If the variable is set to no, the value of
        !% the cutoff must be given by the <tt>XCDensityCorrectionCutoff</tt>
        !% variable.
        !%End
        call parse_variable(namespace, 'XCDensityCorrectionOptimize', .true., xcs%xcd_optimize_cutoff)

        !%Variable XCDensityCorrectionCutoff
        !%Type float
        !%Default 0.0
        !%Section Hamiltonian::XC::DensityCorrection
        !%Description
        !% The value of the cutoff applied to the XC density.
        !%End
        call parse_variable(namespace, 'XCDensityCorrectionCutoff', CNST(0.0), xcs%xcd_ncutoff)

        !%Variable XCDensityCorrectionMinimum
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC::DensityCorrection
        !%Description
        !% When enabled, the cutoff optimization will
        !% return the first minimum of the <math>q_{xc}</math> function if it does
        !% not find a value of -1 (<a href=http://arxiv.org/abs/1107.4339>details</a>).
        !% This is required for atoms or small
        !% molecules, but may cause numerical problems.
        !%End
        call parse_variable(namespace, 'XCDensityCorrectionMinimum', .true., xcs%xcd_minimum)

        !%Variable XCDensityCorrectionNormalize
        !%Type logical
        !%Default true
        !%Section Hamiltonian::XC::DensityCorrection
        !%Description
        !% When enabled, the correction will be
        !% normalized to reproduce the exact boundary conditions of
        !% the XC potential.
        !%End
        call parse_variable(namespace, 'XCDensityCorrectionNormalize', .true., xcs%xcd_normalize)

      end if

      !%Variable ParallelXC
      !%Type logical
      !%Default true
      !%Section Execution::Parallelization
      !%Description
      !% When enabled, additional parallelization
      !% will be used for the calculation of the XC functional.
      !%End
      call messages_obsolete_variable(namespace, 'XCParallel', 'ParallelXC')
      call parse_variable(namespace, 'ParallelXC', .true., xcs%parallel)
      
      POP_SUB(xc_init.parse)
    end subroutine parse

  end subroutine xc_init


  ! ---------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_t), intent(inout) :: xcs

    integer :: isp

    PUSH_SUB(xc_end)

    do isp = 1, 2
      call xc_functl_end(xcs%functional(FUNC_X, isp))
      call xc_functl_end(xcs%functional(FUNC_C, isp))
      call xc_functl_end(xcs%kernel(FUNC_X, isp))
      call xc_functl_end(xcs%kernel(FUNC_C, isp))
    end do
    xcs%family = 0
    xcs%flags  = 0

    POP_SUB(xc_end)
  end subroutine xc_end

  ! ---------------------------------------------------------
  logical function xc_is_orbital_dependent(xcs)
    type(xc_t), intent(in) :: xcs

    PUSH_SUB(xc_is_orbital_dependent)

    xc_is_orbital_dependent = (xcs%cam_alpha + xcs%cam_beta) /= M_ZERO .or. &
      bitand(xcs%functional(FUNC_X,1)%family, XC_FAMILY_OEP) /= 0 .or. &
      bitand(xcs%family, XC_FAMILY_MGGA) /= 0

    POP_SUB(xc_is_orbital_dependent)
  end function xc_is_orbital_dependent
  
  logical function family_is_mgga_with_exc(xcs)
    type(xc_t), intent(in) :: xcs

    integer :: ixc  

    PUSH_SUB(family_is_mgga_with_exc)

    family_is_mgga_with_exc = .false.
    do ixc = 1, 2
      if ((bitand(xcs%functional(ixc, 1)%family, XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0) &
        .and. (bitand(xcs%functional(ixc, 1)%flags, XC_FLAGS_HAVE_EXC) /= 0 )) family_is_mgga_with_exc = .true.
    end do

    POP_SUB(family_is_mgga_with_exc)
  end function family_is_mgga_with_exc

  logical function family_is_hybrid(xcs)
    type(xc_t), intent(in) :: xcs

    integer :: ixc

    PUSH_SUB(family_is_hybrid)

    family_is_hybrid = .false.
    do ixc = 1, 2
      if ((bitand(xcs%functional(ixc, 1)%family, XC_FAMILY_HYB_GGA + XC_FAMILY_HYB_MGGA) /= 0)) &
        family_is_hybrid = .true.
    end do

    POP_SUB(family_is_hybrid)
  end function family_is_hybrid

#include "vxc_inc.F90"
#include "fxc_inc.F90"
#include "kxc_inc.F90"

end module xc_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
