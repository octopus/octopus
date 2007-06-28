!! Copyright (C) 2007 the octopus team
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
!! $Id: response.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module pert_m
  use datasets_m
  use external_pot_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use projector_m
  use specie_pot_m
  use states_m
  use varinfo_m

  implicit none

  private
  public ::               &
     pert_t,         &
     pert_init,      &
     pert_end,       &
     pert_info,      &
     pert_setup_dir, &
     pert_setup_atom,&
     pert_type,      &
     dpert_apply,    &
     zpert_apply,    &
     dpert_apply_order_2,    &
     zpert_apply_order_2,    &
     dpert_expectation_value,&
     zpert_expectation_value,&
     dpert_expectation_density,&
     zpert_expectation_density


  integer, public, parameter :: &
     PERTURBATION_ELECTRIC = 1, &
     PERTURBATION_MAGNETIC = 2, &
     PERTURBATION_IONIC    = 3

  integer, public, parameter :: &
       GAUGE_GIPAW  = 1, &
       GAUGE_ICL    = 2

  type pert_t
     private
    integer :: pert_type
    integer :: dir
    integer :: dir2
    integer :: atom1, atom2
    integer :: gauge
  end type pert_t

  interface pert_init
    module procedure pert_init1
    module procedure pert_init2
  end interface

contains

  ! --------------------------------------------------------------------
  subroutine pert_init1(this, gr, geo)
    type(pert_t), intent(out)   :: this
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    integer :: ii

    !%Variable RespPerturbationType
    !%Type integer
    !%Section Linear Response
    !%Description
    !% Which perturbation to consider.
    !%Option electric 1
    !% Electric perturbation used to calculate, e.g., electric polarizabilities
    !% and hyperpolarizabilities
    !%Option magnetic 2
    !% Magnetic perturbation used to calculate magnetic susceptibilities
    !%Option ionic 3
    !% Displacements of the ions, used to calculate phonon frequencies and
    !% electron-phonon couplings
    !%End 
    call loct_parse_int(check_inp('RespPerturbationType'), PERTURBATION_ELECTRIC, ii)
    call pert_init2(this, ii, gr, geo)

  end subroutine pert_init1


  ! --------------------------------------------------------------------
  subroutine pert_init2(this, pert_type, gr, geo)
    type(pert_t), intent(out)   :: this
    integer,           intent(in)    :: pert_type
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    this%pert_type = pert_type
    if(.not.varinfo_valid_option('RespPerturbationType', this%pert_type)) call input_error('RespPerturbationType')

    this%dir = -1 
    this%dir2 = -1 
    this%atom1 = -1
    this%atom2 = -1

    if ( this%pert_type == PERTURBATION_MAGNETIC ) then

      !%Variable MagneticGaugeCorrection
      !%Type integer
      !%Default gipaw
      !%Section Linear Response
      !%Description
      !% How to describe the coupling of electrons to the magnetic
      !% field, that is how to handle gauge invariance. Default is
      !% gipaw.
      !%Option none 0
      !% No correction
      !%Option gipaw 1
      !% GIPAW correction: Pickard and Mauri, PRL 91 196401 (2003).
      !%End 
      
      call loct_parse_int(check_inp('MagneticGaugeCorrection'), GAUGE_GIPAW, this%gauge)

    end if

  end subroutine pert_init2

  ! --------------------------------------------------------------------
  subroutine pert_end(this)
    type(pert_t), intent(inout) :: this
    
  end subroutine pert_end

  ! --------------------------------------------------------------------
  subroutine pert_info(this, unit)
    type(pert_t), intent(in) :: this
    integer,           intent(in) :: unit

    call messages_print_var_option(stdout, 'RespPerturbationType', this%pert_type)
   
  end subroutine pert_info


  ! --------------------------------------------------------------------
  subroutine pert_setup_dir(this, dir, dir2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: dir
    integer, optional, intent(in)    :: dir2

    this%dir = dir
    if(present(dir2)) this%dir2 = dir2

  end subroutine pert_setup_dir
  
  subroutine pert_setup_atom(this, iatom, iatom2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer, optional, intent(in)    :: iatom2

    this%atom1 = iatom
    if(present(iatom2)) this%atom2 = iatom2

  end subroutine pert_setup_atom

  integer pure function pert_type(this)
    type(pert_t), intent(in) :: this

    pert_type = this%pert_type

  end function pert_type

#include "undef.F90"
#include "real.F90"
#include "pert_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "pert_inc.F90"

end module pert_m
