!! Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module lib_xc
  implicit none

  integer, parameter :: &
     XC_UNPOLARIZED       =  1, &  ! Spin unpolarized
     XC_POLARIZED         =  2     ! Spin polarized

  integer, parameter :: &
     XC_NON_RELATIVISTIC  =  0, &  ! Functional includes or not realtivistic
     XC_RELATIVISTIC      =  1     ! corrections. Only available in some functionals.

  integer, parameter :: &
     XC_LDA_X             =  1, &  ! Exchange                  
     XC_LDA_C_WIGNER      =  2, &  ! Wigner parametrization    
     XC_LDA_C_RPA         =  3, &  ! Random Phase Approximation
     XC_LDA_C_HL          =  4, &  ! Hedin & Lundqvist         
     XC_LDA_C_GL          =  5, &  ! Gunnarson & Lundqvist     
     XC_LDA_C_XALPHA      =  6, &  ! Slaters Xalpha           
     XC_LDA_C_VWN         =  7, &  ! Vosko, Wilk, & Nussair    
     XC_LDA_C_PZ          =  8, &  ! Perdew & Zunger           
     XC_LDA_C_OB_PZ       =  9, &  ! Ortiz & Ballone (PZ)      
     XC_LDA_C_PW          = 10, &  ! Perdew & Wang             
     XC_LDA_C_OB_PW       = 11, &  ! Ortiz & Ballone (PW)      
     XC_LDA_C_LYP         = 12, &  ! Lee, Yang, & Parr LDA     
     XC_LDA_C_AMGB        = 13     ! Attacalite et al

  integer, parameter :: XC_LDA_N = XC_LDA_C_AMGB

  interface
    subroutine xc_lda_init(p, functional, nspin)
      integer(POINTER_SIZE), intent(out) :: p
      integer,               intent(in)  :: functional
      integer,               intent(in)  :: nspin
    end subroutine xc_lda_init

    subroutine xc_lda_end(p)
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine xc_lda_end

    subroutine xc_lda(p, rho, e, v)
      integer(POINTER_SIZE), intent(in)  :: p
      FLOAT,                 intent(in)  :: rho   ! rho(nspin) the density
      FLOAT,                 intent(out) :: e     ! the energy per unit particle
      FLOAT,                 intent(out) :: v     ! v(nspin) the potential
    end subroutine xc_lda

    subroutine xc_lda_x_init(p, nspin, dim, rel)
      integer(POINTER_SIZE), intent(out) :: p
      integer,               intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,               intent(in)  :: dim    ! 2 or 3 dimensions
      integer,               intent(in)  :: rel    ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_lda_x_init
    
    subroutine xc_lda_c_xalpha_init(p, nspin, dim, rel, alpha)
      integer(POINTER_SIZE), intent(out) :: p
      integer,               intent(in)  :: nspin  ! XC_UNPOLARIZED or XC_POLARIZED
      integer,               intent(in)  :: dim    ! 2 or 3 dimensions
      integer,               intent(in)  :: rel    ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
      FLOAT,                 intent(in)  :: alpha  ! Ec = alpha Ex
    end subroutine xc_lda_c_xalpha_init

  end interface
end module lib_xc
