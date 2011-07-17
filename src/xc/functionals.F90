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

module xc_functl_m
  use datasets_m
  use global_m
  use parser_m
  use messages_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::                     &
    xc_functl_t,                &
    xc_functl_init_functl,      &
    xc_functl_end,              &
    xc_functl_write_info


  ! This adds to the constants defined in lib_xc. But since in that module
  ! the OEP functionals are not included, it is better to put it here.
  integer, public, parameter :: &
    XC_KS_INVERSION = 801,      &  ! inversion of Kohn-Sham potential
    XC_OEP_X = 901,             &  ! Exact exchange
    XC_FAMILY_KS_INVERSION = 64

  type xc_functl_t
    integer         :: family            ! LDA, GGA, etc.
    integer         :: type              ! exchange, correlation, or exchange-correlation
    integer         :: id                ! identifier

    integer         :: spin_channels     ! XC_UNPOLARIZED | XC_POLARIZED
    integer         :: flags             ! XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...

    type(XC_F90(pointer_t)) :: conf         ! the pointer used to call the library
    type(XC_F90(pointer_t)) :: info         ! information about the functional

    integer         :: LB94_modified     ! should I use a special version of LB94 that
    FLOAT           :: LB94_threshold    ! needs to be handled specially
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
    logical :: ok
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

      if(functl%family == XC_FAMILY_UNKNOWN) then
        if(functl%id == XC_OEP_X) then
          functl%family = XC_FAMILY_OEP
        else if (functl%id == XC_KS_INVERSION) then
          functl%family = XC_FAMILY_KS_INVERSION
        else
          call input_error('XCFunctional')
        end if
      end if
    end if

    if(functl%family == XC_FAMILY_OEP) then
      functl%type = XC_EXCHANGE

    else if(functl%family == XC_FAMILY_KS_INVERSION) then
      functl%type = XC_EXCHANGE_CORRELATION

    else if(functl%family .eq. XC_FAMILY_NONE) then
      functl%type = -1

    else ! handled by libxc
      ! initialize
      call XC_F90(func_init)(functl%conf, functl%info, functl%id, spin_channels)
      functl%type     = XC_F90(info_kind)(functl%info)
      functl%flags    = XC_F90(info_flags)(functl%info)

      ok = iand(functl%flags, XC_FLAGS_1D).ne.0
      if((ndim.ne.1).and.ok) then
        message(1) = 'Specified functional is only allowed in 1D.'
        call messages_fatal(1)
      end if
      if(ndim==1.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 1D.'
        call messages_fatal(1)
      end if

      ok = iand(functl%flags, XC_FLAGS_2D).ne.0
      if((ndim.ne.2).and.ok) then
        message(1) = 'Specified functional is only allowed in 2D.'
        call messages_fatal(1)
      end if
      if(ndim==2.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 2D.'
        call messages_fatal(1)
      end if

      ok = iand(functl%flags, XC_FLAGS_3D).ne.0
      if((ndim.ne.3).and.ok) then
        message(1) = 'Specified functional is only allowed in 3D.'
        call messages_fatal(1)
      end if
      if(ndim==3.and.(.not.ok)) then
        message(1) = 'Cannot use the specified functionals in 3D.'
        call messages_fatal(1)
      end if
    end if

    ! special parameters that have to be configured
    select case(functl%id)
    case(XC_LDA_C_XALPHA)
      call parse_float(datasets_check('Xalpha'), M_ONE, alpha)
      call XC_F90(lda_c_xalpha_set_par)(functl%conf, alpha)

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
      !% Soft Coulomb interaction of the form 1/sqrt(x^2 + alpha^2). This is the default.
      !%End
      call messages_obsolete_variable('SoftInteraction1D_alpha', 'Interaction1D')
      call parse_integer(datasets_check('Interaction1D'), INT_SOFT_COULOMB, interact_1d)

      !%Variable Interaction1DScreening
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% Defines the screening parameter, alpha, of the softened Coulomb interaction
      !% when running in 1D. The default value is 1.0.
      !%End
      call messages_obsolete_variable('SoftInteraction1D_alpha', 'Interaction1DScreening')
      call parse_float(datasets_check('Interaction1DScreening'), M_ONE, alpha)
      
      if(functl%id == XC_LDA_X_1D) then
        call XC_F90(lda_x_1d_set_par)(functl%conf, interact_1d, alpha)
      else
        call XC_F90(lda_c_1d_csc_set_par)(functl%conf, interact_1d, alpha)
      end if

    case(XC_LDA_C_2D_PRM)
      call XC_F90(lda_c_2d_prm_set_par)(functl%conf, nel)

    case(XC_GGA_X_LB)
      call parse_integer(datasets_check('LB94_modified'), 0, functl%LB94_modified)
      call parse_float(datasets_check('LB94_threshold'), CNST(1.0e-6), functl%LB94_threshold)
      
    end select
    
    POP_SUB(xc_functl_init_functl)
  end subroutine xc_functl_init_functl
  

  ! ---------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_t), intent(inout) :: functl

    PUSH_SUB(xc_functl_end)

    if(functl%family.ne.XC_FAMILY_NONE .and. functl%family.ne.XC_FAMILY_OEP .and.  &
         functl%family.ne.XC_FAMILY_KS_INVERSION ) then
      call XC_F90(func_end)(functl%conf)
    end if

    POP_SUB(xc_functl_end)
  end subroutine xc_functl_end


  ! ---------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_t), intent(in) :: functl
    integer,           intent(in) :: iunit

    character(len=120) :: s1, s2
    type(XC_F90(pointer_t)) :: str
    integer :: ii

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

    else if(functl%family .ne. XC_FAMILY_NONE) then ! all the other families
      select case(functl%type)
      case(XC_EXCHANGE)
        write(message(1), '(2x,a)') 'Exchange'
      case(XC_CORRELATION)
        write(message(1), '(2x,a)') 'Correlation'
      case(XC_EXCHANGE_CORRELATION)
        write(message(1), '(2x,a)') 'Exchange-correlation'
      end select

      call XC_F90(info_name)  (functl%info, s1)
      select case(functl%family)
        case (XC_FAMILY_LDA);      write(s2,'(a)') "LDA"
        case (XC_FAMILY_GGA);      write(s2,'(a)') "GGA"
        case (XC_FAMILY_HYB_GGA);  write(s2,'(a)') "Hybrid GGA"
        case (XC_FAMILY_MGGA);     write(s2,'(a)') "MGGA"
        case (XC_FAMILY_LCA);      write(s2,'(a)') "LCA"
      end select
      write(message(2), '(4x,4a)') trim(s1), ' (', trim(s2), ')'
      call messages_info(2, iunit)
      
      ii = 0
      call XC_F90(info_refs)(functl%info, ii, str, s1)
      do while(ii >= 0)
        write(message(1), '(4x,a,i1,2a)') '[', ii, '] ', trim(s1)
        call messages_info(1, iunit)
        call XC_F90(info_refs)(functl%info, ii, str, s1)
      end do
    end if

    POP_SUB(xc_functl_write_info)
  end subroutine xc_functl_write_info

end module xc_functl_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
