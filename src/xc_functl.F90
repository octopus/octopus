!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module xc_functl
  use global
  use lib_xc

  implicit none

  type xc_functl_type
    integer :: family              ! LDA, GGA, etc.
    integer :: id                  ! identifier
    
    integer :: ispin               ! UNPOLARIZED | SPIN_POLARIZED | SPINORS
    integer :: spin_channels       ! XC_UNPOLARIZED | XC_POLARIZED
    
    integer(POINTER_SIZE) :: conf  ! the pointer used to call the library
    integer(POINTER_SIZE) :: info  ! information about the functional
  end type xc_functl_type
  
contains

  ! -----------------------------------------------------------
  subroutine xc_functl_init(functl, ispin, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: ispin
    integer,              intent(in)  :: spin_channels

    functl%family = 0
    functl%id     = 0
    functl%spin_channels = spin_channels
    functl%ispin  = ispin
 
  end subroutine xc_functl_init


  ! -----------------------------------------------------------
  subroutine xc_functl_init_exchange(functl, ispin, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: ispin
    integer,              intent(in)  :: spin_channels
      
    integer :: rel, j
    FLOAT :: alpha
      
    ! initialize structure
    call xc_functl_init(functl, ispin, spin_channels)

    ! read input
    call loct_parse_int('XFunctional', XC_LDA_X, functl%id)

    ! initialize
    select case(functl%id)
    case(0)
        
    case(XC_LDA_X)
      functl%family = XC_FAMILY_LDA
      call loct_parse_int('LDAX', XC_NON_RELATIVISTIC, rel)
      call xc_lda_x_init(functl%conf, functl%info, &
         spin_channels, conf%dim, rel)
      
    case(XC_GGA_X_PBE, XC_GGA_XC_LB)
      functl%family = XC_FAMILY_GGA
      
      if(functl%id == XC_GGA_XC_LB) then
        call loct_parse_int  ("LB94_modified", 0, j)
        call loct_parse_float("LB94_threshold", CNST(1.0e-6), alpha)
        call xc_gga_lb_init(functl%conf, functl%info, &
           spin_channels, j, alpha)
      else
        call xc_gga_init(functl%conf, functl%info, functl%id, spin_channels)
      end if
      
    case(XC_OEP_X)
      functl%family = XC_FAMILY_OEP
      
    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
         "' is not a known exchange functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select
      
  end subroutine xc_functl_init_exchange
  

  ! -----------------------------------------------------------
  subroutine xc_functl_init_correlation(functl, ispin, spin_channels)
    type(xc_functl_type), intent(out) :: functl
    integer,              intent(in)  :: ispin
    integer,              intent(in)  :: spin_channels
    
    integer :: rel
    FLOAT :: alpha

    ! initialize structure
    call xc_functl_init(functl, ispin, spin_channels)

    ! read input
    call loct_parse_int('CFunctional', XC_LDA_C_PZ, functl%id)

    ! initialize
    select case(functl%id)
    case(0)
      
    case(XC_LDA_C_WIGNER, XC_LDA_C_RPA, XC_LDA_C_HL, XC_LDA_C_GL, XC_LDA_C_XALPHA, &
       XC_LDA_C_VWN, XC_LDA_C_PZ, XC_LDA_C_OB_PZ, XC_LDA_C_PW, XC_LDA_C_OB_PW,     &
       XC_LDA_C_LYP, XC_LDA_C_AMGB)
      
      if(functl%id==XC_LDA_C_AMGB.and.conf%dim.ne.2) then
        message(1) = 'Functional AMGB only allowed in 2D'
        call write_fatal(1)       
      end if
      
      functl%family = XC_FAMILY_LDA
      
      if(functl%id.ne.XC_LDA_C_XALPHA) then
        call xc_lda_init(functl%conf, functl%info, functl%id, spin_channels)
      else
        call loct_parse_int('LDAX', XC_NON_RELATIVISTIC, rel)
        ! WARNING: check what is the most convenient default for alpha
        call loct_parse_float('Xalpha', -M_ONE/M_THREE, alpha) 
        call xc_lda_c_xalpha_init(functl%conf, functl%info, &
           spin_channels, conf%dim, rel, alpha)
      end if
      
    case(XC_GGA_C_PBE)
      functl%family = XC_FAMILY_GGA
      call xc_gga_init(functl%conf, functl%info, functl%id, spin_channels)
      
    case default
      write(message(1), '(a,i3,a)') "'", functl%id, &
         "' is not a known correlation functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select
    
  end subroutine xc_functl_init_correlation

  
  ! -----------------------------------------------------------
  subroutine xc_functl_end(functl)
    type(xc_functl_type), intent(inout) :: functl
    
    select case(functl%family)
    case(XC_FAMILY_LDA); call xc_lda_end(functl%conf)
    case(XC_FAMILY_GGA); call xc_gga_end(functl%conf)
    end select

  end subroutine xc_functl_end


  ! -----------------------------------------------------------
  subroutine xc_functl_write_info(functl, iunit)
    type(xc_functl_type), intent(in) :: functl
    integer,              intent(in) :: iunit

    character(len=120) :: s1, s2
    integer(POINTER_SIZE) :: i
    
    if(functl%family==XC_FAMILY_LDA.or.functl%family==XC_FAMILY_GGA) then
      ! we hapilly call the xc library

      i = xc_info_kind(functl%info)
      select case(i)
      case(int(XC_EXCHANGE, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Exchange'
      case(int(XC_CORRELATION, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Correlation'
      case(int(XC_EXCHANGE_CORRELATION, POINTER_SIZE))
        write(iunit, '(2x,a)') 'Exchange-correlation'
      end select
      
      call xc_info_name  (functl%info, s1)
      call xc_info_family(functl%info, s2)
      write(iunit, '(4x,4a)') trim(s1), ' (', trim(s2), ')'
      
      i = 0;
      call xc_info_ref(functl%info, i, s1)
      do while(i>=0)
        write(iunit, '(4x,a,i1,2a)') '[', i, '] ', trim(s1)
        call xc_info_ref(functl%info, i, s1)
      end do
      
    else if(functl%family==XC_FAMILY_OEP) then 
      ! this is handled separately

      select case(functl%id)
      case(XC_OEP_X)
        write(iunit, '(2x,a)') 'Exchange'
      end select
      
      select case(functl%id)
      case(XC_OEP_X)
        s1 = 'Exact exchange'
      end select
      write(iunit, '(4x,2a)') trim(s1), ' (OEP)'
    end if
      
  end subroutine xc_functl_write_info


end module xc_functl
