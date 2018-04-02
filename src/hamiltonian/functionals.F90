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
  use parser_oct_m
  use libvdwxc_oct_m
  use messages_oct_m
  use XC_F90(lib_m)

  implicit none

  ! Although the following file only contain comments, we include it here to make sure it exists.
  ! Otherwise the code might compile, but not run properly, as the variables documentations
  ! will be incomplete.
#include "functionals_list.F90"
  
  private
  public ::                     &
    xc_functl_t,                &
    xc_functl_init_functl,      &
    xc_functl_end,              &
    xc_functl_write_info


  !> This adds to the constants defined in libxc. But since in that module
  !! the OEP functionals are not included, it is better to put it here.
  integer, public, parameter :: &
    XC_KS_INVERSION = 801,      &  !< inversion of Kohn-Sham potential
    XC_OEP_X = 901,             &  !< Exact exchange
    XC_LDA_XC_CMPLX = 701,      &  !< complex-scaled LDA exchange and correlation
    XC_PBE_XC_CMPLX = 702,      &  !< complex-scaled PBE exchange and correlation
    XC_LB94_XC_CMPLX = 703,     &  !< complex-scaled LB94 exchange and correlation
    XC_HALF_HARTREE = 917,      &  !< half-Hartree exchange for two electrons (supports complex scaling)
    XC_VDW_C_VDWDF = 918,       &  !< vdw-df correlation from libvdwxc
    XC_VDW_C_VDWDF2 = 919,      &  !< vdw-df2 correlation from libvdwxc
    XC_VDW_C_VDWDFCX = 920,     &  !< vdw-df-cx correlation from libvdwxc
    XC_RDMFT_XC_M = 601            !< RDMFT Mueller functional

  !> declaring 'family' constants for 'functionals' not handled by libxc
  !! careful not to use a value defined in libxc for another family!
  integer, public, parameter :: &
    XC_FAMILY_KS_INVERSION = 1024, &
    XC_FAMILY_RDMFT = 2048, &
    XC_FAMILY_LIBVDWXC = 4096

#ifndef HAVE_LIBXC_HYB_MGGA
  integer, public, parameter :: XC_FAMILY_HYB_MGGA = 64
#endif

  type xc_functl_t
    integer         :: family            !< LDA, GGA, etc.
    integer         :: type              !< exchange, correlation, or exchange-correlation
    integer         :: id                !< identifier

    integer         :: spin_channels     !< XC_UNPOLARIZED | XC_POLARIZED
    integer         :: flags             !< XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...

    type(XC_F90(pointer_t)) :: conf         !< the pointer used to call the library
    type(XC_F90(pointer_t)) :: info         !< information about the functional
    type(libvdwxc_t)        :: libvdwxc     !< libvdwxc data for van der Waals functionals

    integer         :: LB94_modified     !< should I use a special version of LB94 that
    FLOAT           :: LB94_threshold    !< needs to be handled specially
  end type xc_functl_t

contains

  ! ---------------------------------------------------------
  subroutine xc_functl_init(functl, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: spin_channels

    PUSH_SUB(xc_functl_init)

    functl%family = 0
    functl%type   = 0
    functl%id     = 0
    functl%flags  = 0
    functl%spin_channels = spin_channels

    POP_SUB(xc_functl_init)
  end subroutine xc_functl_init

  ! ---------------------------------------------------------

 subroutine xc_functl_init_functl(functl, id, ndim, nel, spin_channels)
    type(xc_functl_t), intent(out) :: functl
    integer,           intent(in)  :: id
    integer,           intent(in)  :: ndim
    FLOAT,             intent(in)  :: nel
    integer,           intent(in)  :: spin_channels

    integer :: interact_1d
    FLOAT   :: alpha
#ifdef HAVE_LIBXC4
    FLOAT   :: parameters(2)
#endif
    logical :: ok, lb94_modified
    integer, parameter :: INT_EXP_SCREENED = 0, INT_SOFT_COULOMB = 1

    PUSH_SUB(xc_functl_init_functl)

    ! initialize structure
    call xc_functl_init(functl, spin_channels)

    functl%id = id

    if(functl%id == 0) then
      functl%family = XC_FAMILY_NONE
    else
      ! get the family of the functional
      functl%family = XC_F90(family_from_id)(functl%id)
      ! this also ensures it is actually a functional defined by the linked version of libxc

      if(functl%family == XC_FAMILY_UNKNOWN) then
        if(functl%id == XC_OEP_X) then
          functl%family = XC_FAMILY_OEP
        else if (functl%id == XC_KS_INVERSION) then
          functl%family = XC_FAMILY_KS_INVERSION
        else if (functl%id == XC_LDA_XC_CMPLX) then
          call messages_experimental("complex-scaled LDA exchange and correlation")
          functl%family = XC_FAMILY_LDA
        else if(functl%id == XC_HALF_HARTREE) then
          call messages_experimental("half-Hartree exchange")
          functl%family = XC_FAMILY_LDA ! XXX not really
        else if(functl%id == XC_PBE_XC_CMPLX) then
          call messages_experimental("complex-scaled PBE exchange and correlation")
          functl%family = XC_FAMILY_GGA
        else if(functl%id == XC_LB94_XC_CMPLX) then
          call messages_experimental("complex-scaled LB94 exchange and correlation")
          functl%family = XC_FAMILY_GGA
        else if(functl%id == XC_VDW_C_VDWDF .or. functl%id == XC_VDW_C_VDWDF2 .or. functl%id == XC_VDW_C_VDWDFCX) then
          call messages_experimental("van der Waals functionals from libvdwxc")
          functl%family = XC_FAMILY_LIBVDWXC
          !functl%flags = functl%flags + XC_FLAGS_HAVE_VXC + XC_FLAGS_HAVE_EXC
        else if (functl%id == XC_RDMFT_XC_M) then
          functl%family = XC_FAMILY_RDMFT
        else
          call messages_input_error('XCFunctional', 'Unknown functional')
        end if
      end if
    end if

    if(functl%family == XC_FAMILY_OEP) then
      functl%type = XC_EXCHANGE

    else if(functl%family == XC_FAMILY_KS_INVERSION .or. functl%family == XC_FAMILY_RDMFT) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%family == XC_FAMILY_LIBVDWXC) then
      call XC_F90(func_init)(functl%conf, functl%info, XC_LDA_C_PW, spin_channels)
      functl%type = XC_F90(info_kind)(functl%info)
      functl%flags = XC_F90(info_flags)(functl%info)
      ! Convert Octopus code for functional into corresponding libvdwxc code:
      call libvdwxc_init(functl%libvdwxc, functl%id - XC_VDW_C_VDWDF + 1)
    else if(functl%id == XC_LDA_XC_CMPLX &
      .or. functl%id == XC_PBE_XC_CMPLX &
      .or. functl%id == XC_LB94_XC_CMPLX) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%id == XC_HALF_HARTREE) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%family  ==  XC_FAMILY_NONE) then
      functl%type = -1
      functl%flags = 0

    else ! handled by libxc
      ! initialize
      call XC_F90(func_init)(functl%conf, functl%info, functl%id, spin_channels)
      functl%type     = XC_F90(info_kind)(functl%info)
      functl%flags    = XC_F90(info_flags)(functl%info)

      ! FIXME: no need to say this for kernel
      if(iand(functl%flags, XC_FLAGS_HAVE_EXC) == 0) then
        message(1) = 'Specified functional does not have total energy available.'
        message(2) = 'Corresponding component of energy will just be left as zero.'
        call messages_warning(2)
      end if

      if(iand(functl%flags, XC_FLAGS_HAVE_VXC) == 0) then
        message(1) = 'Specified functional does not have XC potential available.'
        message(2) = 'Cannot run calculations. Choose another XCFunctional.'
        call messages_fatal(2)
      end if

      ok = iand(functl%flags, XC_FLAGS_1D) /= 0
      if((ndim /= 1).and.ok) then
        message(1) = 'Specified functional is only allowed in 1D.'
        call messages_fatal(1)
      end if
      if(ndim==1.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 1D.'
        call messages_fatal(1)
      end if

      ok = iand(functl%flags, XC_FLAGS_2D) /= 0
      if((ndim /= 2).and.ok) then
        message(1) = 'Specified functional is only allowed in 2D.'
        call messages_fatal(1)
      end if
      if(ndim==2.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 2D.'
        call messages_fatal(1)
      end if

      ok = iand(functl%flags, XC_FLAGS_3D) /= 0
      if((ndim /= 3).and.ok) then
        message(1) = 'Specified functional is only allowed in 3D.'
        call messages_fatal(1)
      end if
      if(ndim==3.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 3D.'
        call messages_fatal(1)
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
      call parse_variable('Xalpha', M_ONE, alpha)
#ifdef HAVE_LIBXC4
      parameters(1) = alpha
      call XC_F90(func_set_ext_params)(functl%conf, parameters(1))
#else
      call XC_F90(lda_c_xalpha_set_par)(functl%conf, alpha)
#endif
      
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
      call messages_obsolete_variable('SoftInteraction1D_alpha', 'Interaction1D')
      call parse_variable('Interaction1D', INT_SOFT_COULOMB, interact_1d)

      !%Variable Interaction1DScreening
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the screening parameter <math>\alpha</math> of the softened Coulomb interaction
      !% when running in 1D.
      !%End
      call messages_obsolete_variable('SoftInteraction1D_alpha', 'Interaction1DScreening')
      call parse_variable('Interaction1DScreening', M_ONE, alpha)
#ifdef HAVE_LIBXC4
      parameters(1) = real(interact_1d, REAL_PRECISION)
      parameters(2) = alpha
      call XC_F90(func_set_ext_params)(functl%conf, parameters(1))
#else
      if(functl%id == XC_LDA_X_1D) then
        call XC_F90(lda_x_1d_set_par)(functl%conf, interact_1d, alpha)
      else
        call XC_F90(lda_c_1d_csc_set_par)(functl%conf, interact_1d, alpha)
      end if
#endif
      
    case(XC_LDA_C_2D_PRM)
#ifdef HAVE_LIBXC4
      parameters(1) = nel
      call XC_F90(func_set_ext_params)(functl%conf, parameters(1))
#else
      call XC_F90(lda_c_2d_prm_set_par)(functl%conf, nel)
#endif
      
    ! FIXME: libxc has XC_GGA_X_LBM, isn`t that the modified one?
    case(XC_GGA_X_LB)
      !%Variable LB94_modified
      !%Type logical
      !%Default no
      !%Section Hamiltonian::XC
      !%Description
      !% Whether to use a modified form of the LB94 functional (<tt>XCFunctional = xc_gga_x_lb</tt>).
      !%End
      call parse_variable('LB94_modified', .false., lb94_modified)
      if(lb94_modified) then
        functl%LB94_modified = 1
      else
        functl%LB94_modified = 0
      end if

      ! FIXME: libxc seems to have 1e-32 as a threshold, should we not use that?
      !%Variable LB94_threshold
      !%Type float
      !%Default 1.0e-6
      !%Section Hamiltonian::XC
      !%Description
      !% A threshold for the LB94 functional (<tt>XCFunctional = xc_gga_x_lb</tt>).
      !%End
      call parse_variable('LB94_threshold', CNST(1.0e-6), functl%LB94_threshold)
      
    end select
    
    POP_SUB(xc_functl_init_functl)
  end subroutine xc_functl_init_functl
  

  ! ---------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_t), intent(inout) :: functl

    PUSH_SUB(xc_functl_end)

    if(functl%family /= XC_FAMILY_NONE .and. functl%family /= XC_FAMILY_OEP .and. &
      functl%family /= XC_FAMILY_KS_INVERSION .and. &
      functl%id /= XC_LDA_XC_CMPLX .and. functl%id /= XC_HALF_HARTREE .and. &
      functl%id /= XC_PBE_XC_CMPLX .and. functl%id /= XC_LB94_XC_CMPLX) then
      call XC_F90(func_end)(functl%conf)
    end if

    if(functl%family == XC_FAMILY_LIBVDWXC) then
      call libvdwxc_end(functl%libvdwxc)
    end if

    POP_SUB(xc_functl_end)
  end subroutine xc_functl_end


  ! ---------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_t), intent(in) :: functl
    integer,           intent(in) :: iunit

    character(len=120) :: s1, s2
    integer :: ii
#ifndef HAVE_LIBXC3
    type(XC_F90(pointer_t)) :: str
#endif
    
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

    else if(functl%id == XC_LDA_XC_CMPLX) then
      ! this is handled separately for the moment
      ! we will include it in libxc when done with the tests
      write(message(1), '(2x,a)') 'Exchange-Correlation:'
      write(message(2), '(4x,a)') 'Complex-scaled LDA'
      call messages_info(2, iunit)

    else if(functl%id == XC_HALF_HARTREE) then
      write(message(1), '(2x,a)') 'Exchange-Correlation:'
      write(message(2), '(4x,a)') 'Half-Hartree two-electron exchange'
      call messages_info(2, iunit)

    else if(functl%id == XC_PBE_XC_CMPLX) then
      write(message(1), '(2x,a)') 'Exchange-Correlation:'
      write(message(2), '(4x,a)') 'Complex-scaled PBE'
      call messages_info(2, iunit)

    else if(functl%id == XC_LB94_XC_CMPLX) then
      write(message(1), '(2x,a)') 'Exchange-Correlation:'
      write(message(2), '(4x,a)') 'Complex-scaled LB94'
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
        call messages_not_implemented("kinetic-energy functionals")
      case default
        write(message(1), '(a,i6,a,i6)') "Unknown functional type ", functl%type, ' for functional ', functl%id
        call messages_fatal(1)
      end select

      call XC_F90(info_name)  (functl%info, s1)
      select case(functl%family)
        case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
        case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
        case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
        case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
        case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
      end select
      write(message(2), '(4x,4a)') trim(s1), ' (', trim(s2), ')'
      call messages_info(2, iunit)
      
      ii = 0
#if defined HAVE_LIBXC3 || defined HAVE_LIBXC4
      call XC_F90(info_refs)(functl%info, ii, s1)
#else
      call XC_F90(info_refs)(functl%info, ii, str, s1)
#endif
      do while(ii >= 0)
        write(message(1), '(4x,a,i1,2a)') '[', ii, '] ', trim(s1)
        call messages_info(1, iunit)
#if defined HAVE_LIBXC3 || defined HAVE_LIBXC4
        call XC_F90(info_refs)(functl%info, ii, s1)
#else
        call XC_F90(info_refs)(functl%info, ii, str, s1)
#endif
      end do
    end if

    POP_SUB(xc_functl_write_info)
  end subroutine xc_functl_write_info

end module xc_functl_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
