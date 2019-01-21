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

module propagator_base_oct_m
  use exponential_oct_m
  use potential_interpolation_oct_m
  use sparskit_oct_m
  
  implicit none

  private
  public ::                            &
    propagator_t

  integer, public, parameter ::        &
    PROP_ETRS                    = 2,  &
    PROP_AETRS                   = 3,  &
    PROP_EXPONENTIAL_MIDPOINT    = 4,  &
    PROP_CRANK_NICOLSON          = 5,  &
    PROP_CRANK_NICOLSON_SPARSKIT = 6,  &
    PROP_MAGNUS                  = 7,  &
    PROP_QOCT_TDDFT_PROPAGATOR   = 10, &
    PROP_CAETRS                  = 12, &
    PROP_RUNGE_KUTTA4            = 13, &
    PROP_RUNGE_KUTTA2            = 14, & 
    PROP_EXPLICIT_RUNGE_KUTTA4   = 15, &
    PROP_CFMAGNUS4               = 16

  type propagator_t
    integer             :: method           !< Which evolution method to use.
    type(exponential_t) :: te               !< How to apply the propagator \f$ e^{-i H \Delta t} \f$.
    !> Storage of the KS potential of previous iterations.
    type(potential_interpolation_t) :: vksold
    !> Auxiliary function to store the Magnus potentials.
    FLOAT, pointer      :: vmagnus(:, :, :) => null() 
    integer             :: scf_propagation_steps 
    logical             :: first
    type(sparskit_solver_t), pointer :: tdsk
    integer             :: tdsk_size
    FLOAT               :: scf_threshold
  end type propagator_t

end module propagator_base_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
