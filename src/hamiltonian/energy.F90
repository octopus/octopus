!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012 M. Oliveira
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

module energy_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::         &
    energy_t,       &
    energy_nullify, &
    energy_copy

  type energy_t
    ! Energies
    FLOAT :: total       !< Total energy E = Eii + Sum[Eigenvalues] - U + Ex + Ec - Int[n v_xc] 
                         !!                - 1/2 Int[n^e v_pcm] + 1/2 Int[n^n v_pcm] - Int[n v_U]
    FLOAT :: eigenvalues !< Sum[Eigenvalues]
    FLOAT :: exchange
    FLOAT :: correlation
    FLOAT :: vdw
    FLOAT :: xc_j
    FLOAT :: intnvxc     !< Int[n vxc]
    FLOAT :: hartree     !< Hartree      U = (1/2)*Int [n v_Hartree]
    FLOAT :: int_ee_pcm  !< 1/2 [v_Hartree]*[q_pcm_e] dot product of vectors of dimension n_tesserae
    FLOAT :: int_en_pcm  !< 1/2 [v_Hartree]*[q_pcm_n] 
    FLOAT :: int_ne_pcm  !< 1/2 [v_n]*[q_pcm_e] 
    FLOAT :: int_nn_pcm  !< 1/2 [v_n]*[q_pcm_n]
    FLOAT :: int_e_ext_pcm  !< [v_Hartree]*[q_pcm_ext]
    FLOAT :: int_n_ext_pcm  !< [v_n]*[q_pcm_ext]
    FLOAT :: pcm_corr    !< Int[n (v_e_rs + v_n_rs)]
    FLOAT :: kinetic     !< Kinetic energy of the non-interacting (KS) system of electrons
    FLOAT :: extern      !< External     V = <Phi|V|Phi> = Int[n v] (if no non-local pseudos exist)
    FLOAT :: extern_local !< The local part of the external energy ( Int[n v] )
    FLOAT :: extern_non_local !< The part of the external energy coming from the non-local part of the pseudos
    FLOAT :: entropy
    FLOAT :: ts          !< TS
    FLOAT :: berry       !< Berry energy correction = -mu.E - <Vberry>
    FLOAT :: delta_xc    !< the XC derivative discontinuity
    FLOAT :: dft_u       !DFT+U contribution
    FLOAT :: int_dft_u !< Int[n v_U]

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
    FLOAT :: Imextern_local
    FLOAT :: Imextern_non_local
    FLOAT :: Imextern_dft_u
    FLOAT :: Imentropy
    FLOAT :: Imts
    FLOAT :: Imberry
    FLOAT :: Imint_dft_u
  end type energy_t

contains

  subroutine energy_nullify(this)
    type(energy_t), intent(out) :: this

    PUSH_SUB(energy_nullify)

    this%total        = M_ZERO
    this%eigenvalues  = M_ZERO
    this%exchange     = M_ZERO
    this%correlation  = M_ZERO
    this%vdw          = M_ZERO
    this%xc_j         = M_ZERO
    this%intnvxc      = M_ZERO
    this%hartree      = M_ZERO
    this%int_ee_pcm   = M_ZERO
    this%int_en_pcm   = M_ZERO
    this%int_ne_pcm   = M_ZERO
    this%int_nn_pcm   = M_ZERO
    this%int_e_ext_pcm   = M_ZERO
    this%int_n_ext_pcm   = M_ZERO
    this%pcm_corr     = M_ZERO
    this%kinetic      = M_ZERO
    this%extern       = M_ZERO
    this%extern_local = M_ZERO
    this%extern_non_local = M_ZERO
    this%entropy      = M_ZERO
    this%ts           = M_ZERO
    this%berry        = M_ZERO
    this%delta_xc     = M_ZERO
    this%dft_u        = M_ZERO
    this%int_dft_u    = M_ZERO

    this%Imtotal       = M_ZERO
    this%Imeigenvalues = M_ZERO
    this%Imexchange    = M_ZERO
    this%Imcorrelation = M_ZERO
    this%Imxc_j        = M_ZERO
    this%Imintnvxc     = M_ZERO
    this%Imhartree     = M_ZERO
    this%Imkinetic     = M_ZERO
    this%Imextern      = M_ZERO
    this%Imextern_local = M_ZERO
    this%Imextern_non_local = M_ZERO
    this%Imentropy     = M_ZERO
    this%Imts          = M_ZERO
    this%Imberry       = M_ZERO
    this%Imint_dft_u   = M_ZERO

    POP_SUB(energy_nullify)
  end subroutine energy_nullify

  subroutine energy_copy(ein, eout)
    type(energy_t), intent(in)  :: ein
    type(energy_t), intent(out) :: eout

    PUSH_SUB(energy_copy)

    eout%total        = ein%total
    eout%eigenvalues  = ein%eigenvalues
    eout%exchange     = ein%exchange
    eout%correlation  = ein%correlation
    eout%vdw          = ein%vdw
    eout%xc_j         = ein%xc_j
    eout%intnvxc      = ein%intnvxc
    eout%hartree      = ein%hartree
    eout%int_ee_pcm   = ein%int_ee_pcm
    eout%int_en_pcm   = ein%int_en_pcm
    eout%int_nn_pcm   = ein%int_nn_pcm
    eout%int_ne_pcm   = ein%int_ne_pcm
    eout%int_e_ext_pcm   = ein%int_e_ext_pcm
    eout%int_n_ext_pcm   = ein%int_n_ext_pcm
    eout%pcm_corr     = ein%pcm_corr
    eout%kinetic      = ein%kinetic
    eout%extern       = ein%extern
    eout%extern_local = ein%extern_local
    eout%extern_non_local = ein%extern_non_local
    eout%entropy      = ein%entropy
    eout%ts           = ein%ts
    eout%berry        = ein%berry
    eout%delta_xc     = ein%delta_xc
    eout%dft_u        = ein%dft_u
    eout%int_dft_u    = ein%int_dft_u

    eout%Imtotal = ein%Imtotal
    eout%Imeigenvalues = ein%Imeigenvalues
    eout%Imexchange = ein%Imexchange
    eout%Imcorrelation = ein%Imcorrelation
    eout%Imxc_j = ein%Imxc_j
    eout%Imintnvxc = ein%Imintnvxc
    eout%Imhartree = ein%Imhartree
    eout%Imkinetic = ein%Imkinetic
    eout%Imextern = ein%Imextern
    eout%Imextern_local = ein%Imextern_local
    eout%Imextern_non_local = ein%Imextern_non_local
    eout%Imextern_dft_u = ein%Imextern_dft_u
    eout%Imentropy = ein%Imentropy
    eout%Imts = ein%Imts
    eout%Imberry = ein%Imberry
    eout%Imint_dft_u = Ein%Imint_dft_u
    
    POP_SUB(energy_copy)
  end subroutine energy_copy

end module energy_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
