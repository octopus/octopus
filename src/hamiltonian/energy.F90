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

module energy_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  public ::                          &
    energy_t,                        &
    energy_copy

  type energy_t
    ! Energies
    FLOAT :: total       !< Total energy E = Eii + Sum[Eigenvalues] - U + Ex + Ec - Int[n v_xc]
    FLOAT :: eigenvalues !< Sum[Eigenvalues]
    FLOAT :: exchange
    FLOAT :: correlation
    FLOAT :: xc_j
    FLOAT :: intnvxc     !< Int[n vxc]
    FLOAT :: hartree     !< Hartree      U = (1/2)*Int [n v_Hartree]
    FLOAT :: kinetic     !< Kinetic energy of the non-interacting (KS) system of electrons
    FLOAT :: extern      !< External     V = <Phi|V|Phi> = Int[n v] (if no non-local pseudos exist)
    FLOAT :: entropy
    FLOAT :: ts          !< TS
    FLOAT :: berry       !< Berry energy correction = -mu.E - <Vberry>
    FLOAT :: delta_xc    !< the XC derivative discontinuity

    !cmplxscl 
    FLOAT :: Imtotal       
    FLOAT :: Imeigenvalues 
    FLOAT :: Imexchange
    FLOAT :: Imcorrelation
    FLOAT :: Imxc_j
    FLOAT :: Imintnvxc    
    FLOAT :: Imhartree     
    FLOAT :: Imkinetic     
    FLOAT :: Imextern      
    FLOAT :: Imentropy
    FLOAT :: Imts          
    FLOAT :: Imberry       
    
  end type energy_t

contains

  subroutine energy_copy(ein, eout)
    type(energy_t), intent(in)  :: ein
    type(energy_t), intent(out) :: eout

    PUSH_SUB(energy_copy)

    eout%total        = ein%total
    eout%eigenvalues  = ein%eigenvalues
    eout%exchange     = ein%exchange
    eout%correlation  = ein%correlation
    eout%xc_j         = ein%xc_j
    eout%intnvxc      = ein%intnvxc
    eout%hartree      = ein%hartree
    eout%kinetic      = ein%kinetic
    eout%extern       = ein%extern
    eout%entropy      = ein%entropy
    eout%ts           = ein%ts
    eout%berry        = ein%berry
    eout%delta_xc     = ein%delta_xc

    eout%Imtotal = ein%Imtotal
    eout%Imeigenvalues = ein%Imeigenvalues
    eout%Imexchange = ein%Imexchange
    eout%Imcorrelation = ein%Imcorrelation
    eout%Imxc_j = ein%Imxc_j
    eout%Imintnvxc = ein%Imintnvxc
    eout%Imhartree = ein%Imhartree
    eout%Imkinetic = ein%Imkinetic
    eout%Imextern = ein%Imextern
    eout%Imentropy = ein%Imentropy
    eout%Imts = ein%Imts
    eout%Imberry = ein%Imberry

    
    POP_SUB(energy_copy)
  end subroutine energy_copy

end module energy_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
