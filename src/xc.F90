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
use functions
use poisson
use states
use lib_xc
use xc_functl

implicit none

private
public :: xc_type, &
          xc_oep_type, &
          xc_init, &
          xc_end, &
          xc_write_info, &
          xc_oep_SpinFactor, &
          xc_oep_AnalizeEigen, &
          xc_get_vxc, &
          dxc_KLI_solve, zxc_KLI_solve, &
          doep_x, zoep_x, &
          doep_sic, zoep_sic

type xc_type
  logical :: nlcc                   ! repeated from system

  type(xc_functl_type) :: functl(2) ! 1: exchange, 2: correlation

  ! the SIC needs auxiliary LDAs
  integer              :: sic_correction
  type(xc_functl_type) :: sic_aux(2)

  ! the meta-GGA can be implemented in two ways
  integer :: mGGA_implementation  ! 1 => as a GGA like functional
                                  ! 2 => using the OEP method

  ! for the OEP
  integer  :: oep_level  ! 0 = no oep, 1 = Slater, 2 = KLI, 3 = CEDA, 4 = full OEP
  FLOAT    :: oep_mixing ! how much of the function S(r) to add to vxc in every iteration
end type xc_type

type xc_oep_type
  integer          :: eigen_n
  integer, pointer :: eigen_type(:), eigen_index(:)
  FLOAT            :: socc, sfact
  FLOAT,   pointer :: vxc(:), lxc(:,:), uxc_bar(:)
end type xc_oep_type

FLOAT, parameter :: small     = CNST(1.0e-5)
FLOAT, parameter :: tiny      = CNST(1.0e-12)
FLOAT, parameter :: denom_eps = CNST(1.0e-20) ! added to denominators to avoid overflows...

contains

  subroutine xc_write_info(xcs, iunit)
    type(xc_type), intent(in) :: xcs
    integer,       intent(in) :: iunit
    
    integer :: i
    
#ifdef HAVE_MPI
    if(mpiv%node == 0) then
#endif
      write(iunit,'(/,a)') &
         '**********************************************************************'
      write(iunit,'(a)') " Exchange and correlation:"
      do i = 1, 2
        call xc_functl_write_info(xcs%functl(i), iunit)
      end do
      
      if(xcs%sic_correction.ne.0) then
        write(iunit, '(2x,a)') 'Self-interaction corrections according to Perdew-Zunger'
      end if
        
      if(xcs%oep_level.ne.XC_OEP_NONE) then
        write(iunit, '(a)') 'The OEP equation will be handle at the level of:'
        select case(xcs%oep_level)
        case (XC_OEP_SLATER); write(iunit, '(a)') '    Slater approximation'
        case (XC_OEP_KLI);    write(iunit, '(a)') '    KLI approximation'
        case (XC_OEP_CEDA);   write(iunit, '(a)') '    CEDA approximation'
        case (XC_OEP_FULL);   write(iunit, '(a)') '    Full OEP'
        end select
      end if
      write(iunit,'(a,/)') &
         '**********************************************************************'
      
#ifdef HAVE_MPI
    end if
#endif
  end subroutine xc_write_info

  subroutine xc_init(xcs, nlcc, ispin, spin_channels)
    type(xc_type), intent(out) :: xcs
    logical,       intent(in)  :: nlcc
    integer,       intent(in)  :: ispin
    integer,       intent(in)  :: spin_channels
    
    integer :: func, i, j, rel
    integer(POINTER_SIZE) :: info_dummy
    FLOAT :: alpha
    
    call push_sub('xc_init')
    
    xcs%nlcc   = nlcc  ! make a copy of flag indicating non-local core corrections
    
    call xc_functl_init_exchange   (xcs%functl(1), ispin, spin_channels)
    call xc_functl_init_correlation(xcs%functl(2), ispin, spin_channels)
    
    if(any(xcs%functl(:)%family==XC_FAMILY_MGGA)) then
      call loct_parse_int("MGGAimplementation", 1, xcs%mGGA_implementation)
      if(xcs%mGGA_implementation.ne.1.and.xcs%mGGA_implementation.ne.2) then
        message(1) = 'MGGAimplementation can only assume the values:'
        message(2) = '  1 : GEA implementation'
        message(3) = '  2 : OEP implementation'
        call write_fatal(3)
      end if
    end if

    ! check for SIC
    xcs%sic_correction = 0
    if(any(xcs%functl(:)%family==XC_FAMILY_LDA) .or. &
       any(xcs%functl(:)%family==XC_FAMILY_GGA)) then

      call loct_parse_int("SICorrection", 0, xcs%sic_correction)

      ! we need some auxiliary functionals for the SIC
      if(xcs%sic_correction.ne.0) then
        call xc_functl_init_exchange   (xcs%sic_aux(1), SPIN_POLARIZED, XC_POLARIZED)
        call xc_functl_init_correlation(xcs%sic_aux(2), SPIN_POLARIZED, XC_POLARIZED)
      end if
    end if

    ! if OEP we need some extra variables
    if((xcs%sic_correction.ne.0).or.(any(xcs%functl(:)%family==XC_FAMILY_OEP))) then
#if defined(HAVE_MPI)
      message(1) = "OEP is not allowed with the code parallelized on orbitals..."
      call write_fatal(1)
#endif

      call loct_parse_int("OEP_level", XC_OEP_KLI, xcs%oep_level)
      if(xcs%oep_level<0.or.xcs%oep_level>XC_OEP_FULL) then
        message(1) = "OEP_level can only take the values:"
        message(2) = "0 (Slater), 1 (KLI), 2 (CEDA), or 3 (full OEP)"
        call write_fatal(2)
      end if
      if(xcs%oep_level == XC_OEP_FULL) then
        call loct_parse_float("OEP_mixing", M_ONE, xcs%oep_mixing)
      end if
      
    else
      xcs%oep_level = XC_OEP_NONE
    end if
    
    call pop_sub()

  end subroutine xc_init


  ! -----------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_type), intent(inout) :: xcs
    
    integer :: ixc
    
    call xc_functl_end(xcs%functl(1))
    call xc_functl_end(xcs%functl(2))
    
    if(xcs%sic_correction.ne.0) then
      call xc_functl_end(xcs%sic_aux(1))
      call xc_functl_end(xcs%sic_aux(2))
    end if

  end subroutine xc_end


  ! -----------------------------------------------------------
  ! A couple of auxiliary functions for oep
  ! -----------------------------------------------------------
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_type), intent(out) :: oep
    integer,           intent(in)  :: nspin
    
    select case(nspin)
    case(1) ! we need to correct or the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case(2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default
      write(6,'(a,I2)') 'OEP: error cannot handle nspin=', nspin
    end select

  end subroutine xc_oep_SpinFactor
  

  ! -----------------------------------------------------------
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
  

! -----------------------------------------------------------
#include "xc_pot.F90"

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
