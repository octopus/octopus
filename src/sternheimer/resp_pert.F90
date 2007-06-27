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

module resp_pert_m
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
     resp_pert_t,         &
     resp_pert_init,      &
     resp_pert_end,       &
     resp_pert_info,      &
     resp_pert_setup_dir, &
     resp_pert_setup_atom,&
     resp_type,           &
     dresp_pert_apply,    &
     zresp_pert_apply,    &
     dresp_pert_apply_order_2,    &
     zresp_pert_apply_order_2,    &
     dresp_pert_expectation_value,&
     zresp_pert_expectation_value,&
     dresp_pert_expectation_density,&
     zresp_pert_expectation_density


  integer, public, parameter :: &
     RESP_PERTURBATION_ELECTRIC = 1, &
     RESP_PERTURBATION_MAGNETIC = 2, &
     RESP_PERTURBATION_IONIC    = 3

  integer, public, parameter :: &
       GAUGE_GIPAW  = 1, &
       GAUGE_ICL    = 2

  type resp_pert_t
     private
    integer :: resp_type
    integer :: dir
    integer :: dir2
    integer :: atom1, atom2
    integer :: gauge
  end type resp_pert_t

  interface resp_pert_init
    module procedure resp_pert_init1
    module procedure resp_pert_init2
  end interface

contains

  ! --------------------------------------------------------------------
  subroutine resp_pert_init1(this, gr, geo)
    type(resp_pert_t), intent(out)   :: this
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
    call loct_parse_int(check_inp('RespPerturbationType'), RESP_PERTURBATION_ELECTRIC, ii)
    call resp_pert_init2(this, ii, gr, geo)

  end subroutine resp_pert_init1


  ! --------------------------------------------------------------------
  subroutine resp_pert_init2(this, resp_type, gr, geo)
    type(resp_pert_t), intent(out)   :: this
    integer,           intent(in)    :: resp_type
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    this%resp_type = resp_type
    if(.not.varinfo_valid_option('RespPerturbationType', this%resp_type)) call input_error('RespPerturbationType')

    this%dir = -1 
    this%dir2 = -1 
    this%atom1 = -1
    this%atom2 = -1

    if ( this%resp_type == RESP_PERTURBATION_MAGNETIC ) then

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

  end subroutine resp_pert_init2

  ! --------------------------------------------------------------------
  subroutine resp_pert_end(this)
    type(resp_pert_t), intent(inout) :: this
    
  end subroutine resp_pert_end

  ! --------------------------------------------------------------------
  subroutine resp_pert_info(this, unit)
    type(resp_pert_t), intent(in) :: this
    integer,           intent(in) :: unit

    call messages_print_var_option(stdout, 'RespPerturbationType', this%resp_type)
   
  end subroutine resp_pert_info


  ! --------------------------------------------------------------------
  subroutine resp_pert_setup_dir(this, dir, dir2)
    type(resp_pert_t), intent(inout) :: this
    integer,           intent(in)    :: dir
    integer, optional, intent(in)    :: dir2

    this%dir = dir
    if(present(dir2)) this%dir2 = dir2

  end subroutine resp_pert_setup_dir
  
  subroutine resp_pert_setup_atom(this, iatom, iatom2)
    type(resp_pert_t), intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer, optional, intent(in)    :: iatom2

    this%atom1 = iatom
    if(present(iatom2)) this%atom2 = iatom2

  end subroutine resp_pert_setup_atom

  integer pure function resp_type(this)
    type(resp_pert_t), intent(in) :: this

    resp_type = this%resp_type

  end function resp_type

#include "undef.F90"
#include "real.F90"
#include "resp_pert_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "resp_pert_inc.F90"

end module resp_pert_m
