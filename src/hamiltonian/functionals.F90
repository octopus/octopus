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

module xc_functl_oct_m
  use global_oct_m
  use libvdwxc_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use xc_f03_lib_m

  implicit none

  ! Although the following file only contain comments, we include it here to make sure it exists.
  ! Otherwise the code might compile, but not run properly, as the variables documentations
  ! will be incomplete.
#include "functionals_list.F90"
  
  private
  public ::                     &
    xc_functl_t,                &
    xc_functl_init,             &
    xc_functl_end,              &
    xc_functl_write_info


  !> This adds to the constants defined in libxc. But since in that module
  !! the OEP functionals are not included, it is better to put it here.
  integer, public, parameter ::   &
    XC_KS_INVERSION = 801,        &  !< inversion of Kohn-Sham potential
    XC_OEP_X = 901,               &  !< Exact exchange
    XC_OEP_X_SLATER = 902,        &  !< Slater approximation to the exact exchange
    XC_OEP_X_FBE    = 903,        &  !< Exchange approximation based on the force balance equation
    XC_HALF_HARTREE = 917,        &  !< half-Hartree exchange for two electrons (supports complex scaling)
    XC_VDW_C_VDWDF = 918,         &  !< vdw-df correlation from libvdwxc
    XC_VDW_C_VDWDF2 = 919,        &  !< vdw-df2 correlation from libvdwxc
    XC_VDW_C_VDWDFCX = 920,       &  !< vdw-df-cx correlation from libvdwxc
    XC_HYB_GGA_XC_MVORB_HSE06 = 921, &  !< Density-based mixing parameter of HSE06
    XC_HYB_GGA_XC_MVORB_PBEH = 922,  &  !< Density-based mixing parameter of PBE0 
    XC_MGGA_X_NC_BR = 923,           &  !< Noncollinear version of the Becke-Roussel functional
    XC_RDMFT_XC_M = 601              !< RDMFT Mueller functional

  !> declaring 'family' constants for 'functionals' not handled by libxc
  !! careful not to use a value defined in libxc for another family!
  integer, public, parameter :: &
    XC_FAMILY_KS_INVERSION = 1024, &
    XC_FAMILY_RDMFT = 2048, &
    XC_FAMILY_LIBVDWXC = 4096, &
    XC_FAMILY_NC_LDA = 8192, &
    XC_FAMILY_NC_MGGA = 16384

  type xc_functl_t
    ! Components are public by default
    integer         :: family = 0            !< LDA, GGA, etc.
    integer         :: type   = 0            !< exchange, correlation, or exchange-correlation
    integer         :: id     = 0            !< identifier

    integer         :: spin_channels = 0     !< XC_UNPOLARIZED | XC_POLARIZED
    integer         :: flags         = 0     !< XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...

    type(xc_f03_func_t)               :: conf         !< the pointer used to call the library
    type(xc_f03_func_info_t), private :: info         !< information about the functional
    type(libvdwxc_t)                  :: libvdwxc     !< libvdwxc data for van der Waals functionals
  end type xc_functl_t

contains

  ! ---------------------------------------------------------
  subroutine xc_functl_init(functl, namespace, id, ndim, nel, spin_channels)
    type(xc_functl_t),  intent(inout) :: functl
    type(namespace_t),  intent(in)    :: namespace
    integer,            intent(in)    :: id
    integer,            intent(in)    :: ndim
    FLOAT,              intent(in)    :: nel
    integer,            intent(in)    :: spin_channels

    integer :: interact_1d
    FLOAT   :: alpha, parameters(2)
    logical :: ok

    PUSH_SUB(xc_functl_init)

    ! initialize structure
    functl%id = id
    functl%spin_channels = spin_channels

    if(functl%id == 0) then
      functl%family = XC_FAMILY_NONE
    else
      ! get the family of the functional
      functl%family = xc_f03_family_from_id(functl%id)
      ! this also ensures it is actually a functional defined by the linked version of libxc

      if(functl%family == XC_FAMILY_UNKNOWN) then

        select case(functl%id)
        case(XC_OEP_X, XC_OEP_X_SLATER, XC_OEP_X_FBE)
          functl%family = XC_FAMILY_OEP

        case(XC_KS_INVERSION)
          functl%family = XC_FAMILY_KS_INVERSION

        case(XC_HALF_HARTREE)
          call messages_experimental("half-Hartree exchange")
          functl%family = XC_FAMILY_LDA ! XXX not really

        case(XC_VDW_C_VDWDF, XC_VDW_C_VDWDF2, XC_VDW_C_VDWDFCX)
          functl%family = XC_FAMILY_LIBVDWXC
          !functl%flags = functl%flags + XC_FLAGS_HAVE_VXC + XC_FLAGS_HAVE_EXC

        case(XC_RDMFT_XC_M)
          functl%family = XC_FAMILY_RDMFT

        case(XC_HYB_GGA_XC_MVORB_HSE06, XC_HYB_GGA_XC_MVORB_PBEH)
          functl%family = XC_FAMILY_HYB_GGA
 
        case(XC_MGGA_X_NC_BR) 
          functl%family = XC_FAMILY_NC_MGGA

        case default
          call messages_input_error(namespace, 'XCFunctional', 'Unknown functional')

        end select
      end if
    end if

    if(functl%family == XC_FAMILY_OEP) then
      functl%type = XC_EXCHANGE

    else if(functl%family == XC_FAMILY_KS_INVERSION .or. functl%family == XC_FAMILY_RDMFT) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%family == XC_FAMILY_LIBVDWXC) then
      call xc_f03_func_init(functl%conf, XC_LDA_C_PW, spin_channels)
      functl%info = xc_f03_func_get_info(functl%conf)
      functl%type = xc_f03_func_info_get_kind(functl%info)
      functl%flags = xc_f03_func_info_get_flags(functl%info)
      ! Convert Octopus code for functional into corresponding libvdwxc code:
      call libvdwxc_init(functl%libvdwxc, namespace, functl%id - XC_VDW_C_VDWDF + 1)

    else if(functl%id == XC_HALF_HARTREE) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%id == XC_MGGA_X_NC_BR) then
      functl%type = XC_EXCHANGE
      functl%flags = XC_FLAGS_HAVE_VXC + XC_FLAGS_HAVE_EXC

    else if(functl%family == XC_FAMILY_NONE) then
      functl%type = -1
      functl%flags = 0

    else ! handled by libxc
      ! initialize

      !For the two MVORB functionals, we initialize libxc with the non-MVORB functionals
      select case(functl%id)
      case(XC_HYB_GGA_XC_MVORB_HSE06)
        call xc_f03_func_init(functl%conf, XC_HYB_GGA_XC_HSE06, spin_channels)

      case(XC_HYB_GGA_XC_MVORB_PBEH)
        call xc_f03_func_init(functl%conf, XC_HYB_GGA_XC_PBEH, spin_channels)

      case default
        call xc_f03_func_init(functl%conf, functl%id, spin_channels)
      end select
      functl%info     = xc_f03_func_get_info(functl%conf)
      functl%type     = xc_f03_func_info_get_kind(functl%info)
      functl%flags    = xc_f03_func_info_get_flags(functl%info)

      ! FIXME: no need to say this for kernel
      if(bitand(functl%flags, XC_FLAGS_HAVE_EXC) == 0) then
        message(1) = 'Specified functional does not have total energy available.'
        message(2) = 'Corresponding component of energy will just be left as zero.'
        call messages_warning(2, namespace=namespace)
      end if

      if(bitand(functl%flags, XC_FLAGS_HAVE_VXC) == 0) then
        message(1) = 'Specified functional does not have XC potential available.'
        message(2) = 'Cannot run calculations. Choose another XCFunctional.'
        call messages_fatal(2, namespace=namespace)
      end if

      ok = bitand(functl%flags, XC_FLAGS_1D) /= 0
      if((ndim /= 1).and.ok) then
        message(1) = 'Specified functional is only allowed in 1D.'
        call messages_fatal(1, namespace=namespace)
      end if
      if(ndim==1.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 1D.'
        call messages_fatal(1, namespace=namespace)
      end if

      ok = bitand(functl%flags, XC_FLAGS_2D) /= 0
      if((ndim /= 2).and.ok) then
        message(1) = 'Specified functional is only allowed in 2D.'
        call messages_fatal(1, namespace=namespace)
      end if
      if(ndim==2.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 2D.'
        call messages_fatal(1, namespace=namespace)
      end if

      ok = bitand(functl%flags, XC_FLAGS_3D) /= 0
      if((ndim /= 3).and.ok) then
        message(1) = 'Specified functional is only allowed in 3D.'
        call messages_fatal(1, namespace=namespace)
      end if
      if(ndim==3.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 3D.'
        call messages_fatal(1, namespace=namespace)
      end if
    end if
    
!    XC_NON_RELATIVISTIC     =   0
!    XC_RELATIVISTIC         =   1

    ! FIXME: aren`t there other parameters that can or should be set?
    ! special parameters that have to be configured
    select case(functl%id)
      ! FIXME: aren`t there other Xalpha functionals?
    case(XC_LDA_C_XALPHA)

      !%Variable Xalpha
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% The parameter of the Slater X<math>\alpha</math> functional. Applies only for
      !% <tt>XCFunctional = xc_lda_c_xalpha</tt>.
      !%End
      call parse_variable(namespace, 'Xalpha', M_ONE, parameters(1))
      
      call xc_f03_func_set_ext_params(functl%conf, parameters(1))
      
      ! FIXME: doesn`t this apply to other 1D functionals?
    case(XC_LDA_X_1D, XC_LDA_C_1D_CSC)
      !%Variable Interaction1D
      !%Type integer
      !%Default interaction_soft_coulomb
      !%Section Hamiltonian::XC
      !%Description
      !% When running in 1D, one has to soften the Coulomb interaction. This softening
      !% is not unique, and several possibilities exist in the literature.
      !%Option interaction_exp_screened 0
      !% Exponentially screened Coulomb interaction.
      !% See, <i>e.g.</i>, M Casula, S Sorella, and G Senatore, <i>Phys. Rev. B</i> <b>74</b>, 245427 (2006).
      !%Option interaction_soft_coulomb 1
      !% Soft Coulomb interaction of the form <math>1/\sqrt{x^2 + \alpha^2}</math>.
      !%End
      call messages_obsolete_variable(namespace, 'SoftInteraction1D_alpha', 'Interaction1D')
      call parse_variable(namespace, 'Interaction1D', OPTION__INTERACTION1D__INTERACTION_SOFT_COULOMB, interact_1d)

      !%Variable Interaction1DScreening
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the screening parameter <math>\alpha</math> of the softened Coulomb interaction
      !% when running in 1D.
      !%End
      call messages_obsolete_variable(namespace, 'SoftInteraction1D_alpha', 'Interaction1DScreening')
      call parse_variable(namespace, 'Interaction1DScreening', M_ONE, alpha)
      parameters(1) = TOFLOAT(interact_1d)
      parameters(2) = alpha
      call xc_f03_func_set_ext_params(functl%conf, parameters(1))
      
    case(XC_LDA_C_2D_PRM)
      parameters(1) = nel
      call xc_f03_func_set_ext_params(functl%conf, parameters(1))

    case (XC_GGA_X_LB)
      if (parse_is_defined(namespace, 'LB94_modified')) then
        call messages_obsolete_variable(namespace, 'LB94_modified')
      end if

      if (parse_is_defined(namespace, 'LB94_threshold')) then
        call messages_obsolete_variable(namespace, 'LB94_threshold')
      end if
    end select

    POP_SUB(xc_functl_init)
  end subroutine xc_functl_init
  

  ! ---------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_t), intent(inout) :: functl

    PUSH_SUB(xc_functl_end)

    if (functl%family /= XC_FAMILY_NONE .and. functl%family /= XC_FAMILY_OEP .and. &
        functl%family /= XC_FAMILY_KS_INVERSION .and. functl%id /= XC_HALF_HARTREE &
        .and. functl%family /= XC_FAMILY_NC_LDA .and. functl%family /= XC_FAMILY_NC_MGGA) then
      call xc_f03_func_end(functl%conf)
    end if

    if(functl%family == XC_FAMILY_LIBVDWXC) then
      call libvdwxc_end(functl%libvdwxc)
    end if

    POP_SUB(xc_functl_end)
  end subroutine xc_functl_end


  ! ---------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit, namespace)
    type(xc_functl_t), intent(in) :: functl
    integer,           intent(in) :: iunit
    type(namespace_t), intent(in) :: namespace

    character(len=120) :: family
    integer :: ii
    type(xc_f03_func_reference_t) :: ref

    PUSH_SUB(xc_functl_write_info)

    if(functl%family == XC_FAMILY_OEP) then
      ! this is handled separately

      select case(functl%id)
      case(XC_OEP_X)
        write(message(1), '(2x,a)') 'Exchange'
        write(message(2), '(4x,a)') 'Exact exchange'
        call messages_info(2, iunit)
      end select

    else if(functl%family == XC_FAMILY_KS_INVERSION) then
      ! this is handled separately
      select case(functl%id)
      case(XC_KS_INVERSION)
        write(message(1), '(2x,a)') 'Exchange-Correlation:'
        write(message(2), '(4x,a)') '  KS Inversion'
        call messages_info(2, iunit)
      end select

    else if(functl%family == XC_FAMILY_LIBVDWXC) then
      call libvdwxc_write_info(functl%libvdwxc, iunit)

    else if(functl%id == XC_MGGA_X_NC_BR) then
      write(message(1), '(2x,a)') 'Exchange'
      write(message(2), '(4x,a)') 'Noncollinear Becke-Roussel'
      call messages_info(2, iunit)

    else if(functl%id == XC_HALF_HARTREE) then
      write(message(1), '(2x,a)') 'Exchange-Correlation:'
      write(message(2), '(4x,a)') 'Half-Hartree two-electron exchange'
      call messages_info(2, iunit)
      
    else if(functl%family /= XC_FAMILY_NONE) then ! all the other families
      select case(functl%type)
      case(XC_EXCHANGE)
        write(message(1), '(2x,a)') 'Exchange'
      case(XC_CORRELATION)
        write(message(1), '(2x,a)') 'Correlation'
      case(XC_EXCHANGE_CORRELATION)
        write(message(1), '(2x,a)') 'Exchange-correlation'
      case(XC_KINETIC)
        call messages_not_implemented("kinetic-energy functionals", namespace=namespace)
      case default
        write(message(1), '(a,i6,a,i6)') "Unknown functional type ", functl%type, ' for functional ', functl%id
        call messages_fatal(1, namespace=namespace)
      end select

      select case(functl%family)
      case (XC_FAMILY_LDA);       write(family,'(a)') "LDA"
      case (XC_FAMILY_GGA);       write(family,'(a)') "GGA"
      case (XC_FAMILY_HYB_GGA);   write(family,'(a)') "Hybrid GGA"
      case (XC_FAMILY_HYB_MGGA);  write(family,'(a)') "Hybrid MGGA"
      case (XC_FAMILY_MGGA);      write(family,'(a)') "MGGA"
      end select
      write(message(2), '(4x,4a)') trim(xc_f03_func_info_get_name(functl%info)), ' (', trim(family), ')'
      call messages_info(2, iunit)

      ii = 0
      ref = xc_f03_func_info_get_references(functl%info, ii)
      do while(ii >= 0)
        write(message(1), '(4x,a,i1,2a)') '[', ii, '] ', trim(xc_f03_func_reference_get_ref(ref))
        call messages_info(1, iunit)
        ref = xc_f03_func_info_get_references(functl%info, ii)
      end do
    end if

    POP_SUB(xc_functl_write_info)
  end subroutine xc_functl_write_info

end module xc_functl_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
