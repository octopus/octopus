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
use lib_alg
use mesh
use mesh_function
use poisson
use states
use vxc

implicit none

!!! This parameters have to be updated when implementing
!!! new functionals
integer, parameter ::     &
     N_XC_FAMILIES = 4,   &
     N_X_FUNCTL    = 8,   &
     N_C_FUNCTL    = 6

!!! Families of xc functionals
integer, parameter ::     &
     XC_FAMILY_LDA  = ibset(0, 0),    &
     XC_FAMILY_GGA  = ibset(0, 1),    &
     XC_FAMILY_OEP  = ibset(0, 2),    &
     XC_FAMILY_MGGA = ibset(0, 3)

!!! Functionals
integer, parameter ::     &
     X_FUNC_LDA_REL   = ibset(0,  0), &
     X_FUNC_LDA_NREL  = ibset(0,  1), &
     X_FUNC_GGA_PBE   = ibset(0,  2), &
     X_FUNC_GGA_PBER  = ibset(0,  3), &
     X_FUNC_GGA_LB94  = ibset(0,  4), &
     X_FUNC_MGGA_PKZB = ibset(0,  5), &
     X_FUNC_OEP_X     = ibset(0,  6), &
     X_FUNC_OEP_SIC   = ibset(0,  7), &
     C_FUNC_LDA_RPA   = ibset(0,  8), &
     C_FUNC_LDA_PZ    = ibset(0,  9), &
     C_FUNC_LDA_PW92  = ibset(0, 10), &
     C_FUNC_GGA_PBE   = ibset(0, 11), &
     C_FUNC_OEP_SIC   = ibset(0, 12), &
     C_FUNC_MGGA_PKZB = ibset(0, 13)


character(len=4), parameter :: xc_name_families(N_XC_FAMILIES) =   &
     (/ 'LDA ', 'GGA ', 'OEP ', 'MGGA' /)

!!! Exchange functionals should come before correlation functionals
character(len=20), parameter :: xc_name_functionals(N_X_FUNCTL+N_C_FUNCTL) = (/ &
     'LDA  - relat.       ', &  ! X_FUNC_LDA_REL
     'LDA  - non-relat.   ', &  ! X_FUNC_LDA_NREL
     'GGA  - PBE          ', &  ! X_FUNC_GGA_PBE
     'GGA  - PBE - relat. ', &  ! X_FUNC_GGA_PBER
     'GGA  - LB94         ', &  ! X_FUNC_GGA_LB94
     'MGGA - PKZB         ', &  ! X_FUNC_MGGA_PKZB
     'OEP  - EXX          ', &  ! X_FUNC_OEP_X
     'OEP  - SIC (LDA)    ', &  ! X_FUNC_OEP_SIC
     'LDA  - RPA          ', &  ! C_FUNC_LDA_RPA
     'LDA  - PZ81         ', &  ! C_FUNC_LDA_PZ
     'LDA  - PW92         ', &  ! C_FUNC_LDA_PW92
     'GGA  - PBE          ', &  ! C_FUNC_GGA_PBE
     'OEP  - SIC (LDA)    ', &  ! C_FUNC_OEP_SIC
     'MGGA - PKZB         '  &  ! C_FUNC_MGGA_PKZB
     /)

type xc_type
  integer :: family, functl
  logical :: nlcc ! repeated from system

  ! For LB94 fine-tuning
  logical  :: lb94_modified
  FLOAT :: lb94_beta
  FLOAT :: lb94_threshold

  ! for the OEP
  integer  :: oep_level  ! 0 = Slater, 1 = KLI, 2 = full OEP
  FLOAT :: oep_mixing ! how much of the function S(r) to add to vxc in every iteration
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

    do i = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(btest(xcs%functl, i)) then
        if(i < N_X_FUNCTL) then
          write(iunit, '(6x,a,a)') 'Exchange    : ', xc_name_functionals(i+1)
        else
          write(iunit, '(6x,a,a)') 'Correlation : ', xc_name_functionals(i+1)
        end if
      end if
    end do

    if(iand(xcs%family, XC_FAMILY_OEP).ne.0) then
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
end subroutine xc_write_info

subroutine xc_init(xcs, nlcc)
  type(xc_type), intent(out) :: xcs
  logical, intent(in) :: nlcc

  character(len=50) :: func

  call push_sub('xc_init')

  xcs%nlcc   = nlcc  ! make a copy of falg indicating non-local core corrections
  xcs%family = 0     ! initialize xc family and functional to zero
  xcs%functl = 0

  ! set the appropriate flags
  call loct_parse_string('XFunctional', 'LDA', func)
  select case(func(1:4))
  case('ZER ')
  case('LDA ')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, X_FUNC_LDA_NREL)
  case('RLDA')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, X_FUNC_LDA_REL)
  case('PBE')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_PBE)
  case('RPBE')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_PBER)
  case('LB94')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, X_FUNC_GGA_LB94)
  case('PKZB')
    xcs%family = ior(xcs%family, XC_FAMILY_MGGA)
    xcs%functl = ior(xcs%functl, X_FUNC_MGGA_PKZB)
  case('EXX ')
    xcs%family = ior(xcs%family, XC_FAMILY_OEP)
    xcs%functl = ior(xcs%functl, X_FUNC_OEP_X)
  case('SIC ')
    xcs%family = ior(xcs%family, XC_FAMILY_OEP)
    xcs%functl = ior(xcs%functl, X_FUNC_OEP_SIC)
  case default
    write(message(1), '(a,a,a)') "'", trim(func), "' is not a known exchange functional!"
    message(2) = "Please check the manual for a list of possible values."
    call write_fatal(2)
  end select

  ! now the correlation
  call loct_parse_string('CFunctional', 'PZ81', func)
  select case(func(1:4))
  case('ZER ')
  case('PZ81')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, C_FUNC_LDA_PZ)
  case('PW92')
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, C_FUNC_LDA_PW92)
  case('PBE ')
    xcs%family = ior(xcs%family, XC_FAMILY_GGA)
    xcs%functl = ior(xcs%functl, C_FUNC_GGA_PBE)
  case('PKZB')
    xcs%family = ior(xcs%family, XC_FAMILY_MGGA)
    xcs%functl = ior(xcs%functl, C_FUNC_MGGA_PKZB)
  case('SIC ')
    xcs%family = ior(xcs%family, XC_FAMILY_OEP)
    xcs%functl = ior(xcs%functl, C_FUNC_OEP_SIC)
  case default
    write(message(1), '(a,a,a)') "'", trim(func), &
         "' is not a known correlation functional family!"
    message(2) = "(CFamily = ZER | LDA | GGA | OEP)"
    call write_fatal(2)
  end select

#if defined(HAVE_MPI) && defined(MPI_TD)
  if(btest(xcs%family, XC_FAMILY_OEP)) then
    message(1) = "OEP is not allowed with MPI_TD!"
    call write_fatal(1)
  end if
#endif

  ! the SIC potential has an LDA part, so let us add it
  if(iand(xcs%functl, X_FUNC_OEP_SIC).ne.0) then
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, X_FUNC_LDA_NREL)
  end if
  if(iand(xcs%functl, C_FUNC_OEP_SIC).ne.0) then
    xcs%family = ior(xcs%family, XC_FAMILY_LDA)
    xcs%functl = ior(xcs%functl, C_FUNC_LDA_PZ)
  end if

  ! Extra parameters to fine-tune LB94 potential.
  if(iand(xcs%functl, X_FUNC_GGA_LB94).ne.0) then
    call loct_parse_float("LB94_beta", CNST(0.05), xcs%lb94_beta)
    call loct_parse_float("LB94_threshold", CNST(1.0e-6), xcs%lb94_threshold)
    call loct_parse_logical("LB94_modified", .false., xcs%lb94_modified)
  end if

  if(iand(xcs%family, XC_FAMILY_OEP).ne.0) then
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
  return
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
#include "xc_LDA.F90"
#include "xc_GGA.F90"
!#include "xc_MGGA.F90"

#include "undef.F90"
#include "real.F90"
#include "xc_pot.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_pot.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"

#include "undef.F90"

end module xc
