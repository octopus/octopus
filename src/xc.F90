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

module xc
use global
use lib_oct_parser
use lib_basic_alg
use lib_adv_alg
use mesh
use poisson
use states
use lib_xc

implicit none

type xc_type
  logical :: nlcc ! repeated from system

  integer :: family(2)  ! 1: exchange, 2: correlation
  integer :: functl(2)

  integer(POINTER_SIZE) :: conf(2)
  integer(POINTER_SIZE) :: info(2)

  ! for the OEP
  integer  :: oep_level  ! 0 = Slater, 1 = KLI, 2 = full OEP
  FLOAT    :: oep_mixing ! how much of the function S(r) to add to vxc in every iteration
end type xc_type


type xc_oep_type
  integer           :: eigen_n
  integer,  pointer :: eigen_type(:), eigen_index(:)
  FLOAT          :: socc, sfact
  FLOAT, pointer :: vxc(:), lxc(:,:), uxc_bar(:)
end type xc_oep_type

FLOAT, parameter :: small     = CNST(1.0e-5)
FLOAT, parameter :: tiny      = CNST(1.0e-12)
FLOAT, parameter :: denom_eps = CNST(1.0e-20) ! added to denominators to avoid overflows...

contains

  subroutine xc_write_info(xcs, iunit)
    type(xc_type), intent(IN) :: xcs
    integer, intent(in) :: iunit
    
    integer :: i
    
#ifdef HAVE_MPI
    if(mpiv%node == 0) then
#endif
      
      if(xcs%family(1).ne.0) then
        write(iunit, '(2x,a)') 'Exchange'
        call get_info(xcs%info(1), xcs%family(1))
      end if
      
      if(xcs%family(2).ne.0) then
        write(iunit, '(2x,a)') 'Correlation'
        call get_info(xcs%info(2), xcs%family(2))
      end if
      
      if(any(xcs%family==XC_FAMILY_OEP)) then
        write(iunit, '(a)') 'The OEP equation will be handle at the level of:'
        select case(xcs%oep_level)
        case (0); write(iunit, '(a)') '      Slater approximation'
        case (1); write(iunit, '(a)') '      KLI approximation'
        case (2); write(iunit, '(a)') '      Full OEP'
        end select
      end if
      
#ifdef HAVE_MPI
    end if
#endif

  contains
    subroutine get_info(info, family)
      integer(POINTER_SIZE), intent(in) :: info
      integer, intent(in) :: family

      character(len=120) :: s1, s2
      integer :: i
      
      if(family==XC_FAMILY_LDA.or.family==XC_FAMILY_GGA) then
        call xc_info_name  (info, s1)
        call xc_info_family(info, s2)
        write(iunit, '(4x,4a)') trim(s1), ' (', trim(s2), ')'
        
        i = 0;
        call xc_info_ref(info, i, s1)
        do while(i>=0)
          write(iunit, '(4x,a,i1,2a)') '[', i, '] ', trim(s1)
          call xc_info_ref(info, i, s1)
        end do
      else if(family==XC_FAMILY_OEP) then
        write(iunit, *) "OEP"
      end if
      
    end subroutine get_info
  end subroutine xc_write_info
  

  subroutine xc_init(xcs, nlcc, spin_channels)
    type(xc_type), intent(out) :: xcs
    logical, intent(in) :: nlcc
    integer, intent(in) :: spin_channels
    
    integer :: func, i, j, rel
    FLOAT :: alpha
    
    call push_sub('xc_init')
    
    xcs%nlcc   = nlcc  ! make a copy of flag indicating non-local core corrections
    
    xcs%family = 0   ! initialize xc family and functional to zero
    xcs%functl = 0
    
    ! set the appropriate flags
    call loct_parse_int('XFunctional', XC_LDA_X, xcs%functl(1))
    select case(xcs%functl(1))
    case(0)
      
    case(XC_LDA_X)
      xcs%family(1) = XC_FAMILY_LDA
      call loct_parse_int('LDAX', XC_NON_RELATIVISTIC, rel)
      call xc_lda_x_init(xcs%conf(1), xcs%info(1), &
         spin_channels, conf%dim, rel)
      
    case(XC_GGA_X_PBE, XC_GGA_XC_LB)
      xcs%family(1) = XC_FAMILY_GGA
      
      if(xcs%functl(1) == XC_GGA_XC_LB) then
        call loct_parse_int  ("LB94_modified", 0, j)
        call loct_parse_float("LB94_threshold", CNST(1.0e-6), alpha)
        call xc_gga_lb_init(xcs%conf(1), xcs%info(1), &
           spin_channels, j, alpha)
      else
        call xc_gga_init(xcs%conf(1), xcs%info(1), xcs%functl(1), spin_channels)
      end if

    case(XC_OEP_X, XC_OEP_X_SIC)
      xcs%family(1) = XC_FAMILY_OEP

    case default
      write(message(1), '(a,i3,a)') "'", xcs%functl(1), &
       "' is not a known exchange functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select
    
    ! now the correlation
    call loct_parse_int('CFunctional', XC_LDA_C_PZ, xcs%functl(2))
    select case(xcs%functl(2))
    case(0)

    case(XC_LDA_C_WIGNER, XC_LDA_C_RPA, XC_LDA_C_HL, XC_LDA_C_GL, XC_LDA_C_XALPHA, &
       XC_LDA_C_VWN, XC_LDA_C_PZ, XC_LDA_C_OB_PZ, XC_LDA_C_PW, XC_LDA_C_OB_PW,     &
       XC_LDA_C_LYP, XC_LDA_C_AMGB)
      
      xcs%family(2) = XC_FAMILY_LDA
      
      if(xcs%functl(2).ne.XC_LDA_C_XALPHA) then
        call xc_lda_init(xcs%conf(2), xcs%info(2), xcs%functl(2), spin_channels)
      else
        call loct_parse_int('LDAX', XC_NON_RELATIVISTIC, rel)
        ! WARNING: check what is the most convenient default for alpha
        call loct_parse_float('Xalpha', -M_ONE/M_THREE, alpha) 
        call xc_lda_c_xalpha_init(xcs%conf(2), xcs%info(2), &
           spin_channels, conf%dim, rel, alpha)
      end if
      
    case(XC_GGA_C_PBE)
      xcs%family(2) = XC_FAMILY_GGA
      call xc_gga_init(xcs%conf(2), xcs%info(2), xcs%functl(2), spin_channels)
      
    case(XC_OEP_C_SIC)
      xcs%family(1) = XC_FAMILY_OEP

    case default
      write(message(1), '(a,i3,a)') "'", xcs%functl(2), &
         "' is not a known correlation functional!"
      message(2) = "Please check the manual for a list of possible values."
      call write_fatal(2)
    end select
    

    if(conf%dim.ne.2.and.any(xcs%functl==XC_LDA_C_AMGB)) then
      message(1) = 'Functional AMGB only allowed in 2D'
      call write_fatal(1)
    end if

    if(any(xcs%family==XC_FAMILY_OEP)) then
#if defined(HAVE_MPI)
      message(1) = "OEP is not allowed with the code parallelized on orbitals..."
      call write_fatal(1)
#endif

      call loct_parse_int("OEP_level", 1, xcs%oep_level)
      if(xcs%oep_level<0.or.xcs%oep_level>2) then
        message(1) = "OEP_level can only take the values:"
        message(2) = "0 (Slater), 1 (KLI), and 2 (full OEP)"
        call write_fatal(2)
      end if
      if(xcs%oep_level == 2) then
        call loct_parse_float("OEP_mixing", M_ONE, xcs%oep_mixing)
      end if
    end if
    
    call pop_sub()
  end subroutine xc_init

  subroutine xc_end(xcs)
    type(xc_type), intent(inout) :: xcs
    
    integer :: ixc
    
    do ixc = 1, 2
      select case(xcs%family(ixc))
      case(XC_FAMILY_LDA); call xc_lda_end(xcs%conf(ixc))
      case(XC_FAMILY_GGA); call xc_gga_end(xcs%conf(ixc))
      end select
    end do
    
  end subroutine xc_end

  ! A couple of auxiliary functions for oep
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_type), intent(out) :: oep
    integer,           intent(in)  :: nspin
    
    select case(nspin)
    case(1) ! we need to correct for the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case(2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default
      write(6,'(a,I2)') 'OEP: error cannot handle nspin=', nspin
    end select
  end subroutine xc_oep_SpinFactor
  
  subroutine xc_oep_AnalizeEigen(oep, st, is)
    type(xc_oep_type), intent(out) :: oep
    type(states_type), intent(in)  :: st
    integer,           intent(in)  :: is
    
    integer  :: i
    FLOAT :: max_eigen
    
    ! find out the top occupied state, to correct for the assymptotics
    ! of the potential
    max_eigen = CNST(-1e30)
    do i = 1, st%nst
      if((st%occ(i, is) .gt. small).and.(st%eigenval(i, is).gt.max_eigen)) then
        max_eigen = st%eigenval(i, is)
      end if
    end do
    
    oep%eigen_n = 1
    do i = 1, st%nst
      if(st%occ(i, is) .gt. small) then
        ! criterium for degeneracy
        if(abs(st%eigenval(i, is)-max_eigen).le.CNST(1e-3)) then
          oep%eigen_type(i) = 2
        else
          oep%eigen_type(i) = 1
          oep%eigen_index(oep%eigen_n) = i
          oep%eigen_n = oep%eigen_n +1
        end if
      else
        oep%eigen_type(i) = 0
      end if
    end do
    oep%eigen_n = oep%eigen_n - 1
    
  end subroutine xc_oep_AnalizeEigen
  
! include the xc potentials
#include "xc_pot.F90"
!#include "xc_MGGA.F90"

#include "undef.F90"
#include "real.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"

#include "undef.F90"

end module xc
