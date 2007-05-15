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
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use specie_pot_m
  use varinfo_m

  implicit none

  private
  public ::               &
     resp_pert_t,         &
     resp_pert_init,      &
     resp_pert_end,       &
     resp_pert_info,      &
     resp_pert_setup_dir, &
     dresp_pert_apply,    &
     zresp_pert_apply
  
  integer, public, parameter :: &
     RESP_PERTURBATION_ELECTRIC = 1, &
     RESP_PERTURBATION_MAGNETIC = 2, &
     RESP_PERTURBATION_DISPLACE = 3

  type resp_pert_t
    integer :: resp_type

    integer :: dir
    
    FLOAT, pointer :: dHdR  (:, :, :)    ! derivative of the Hamiltonian wrt atomic displacement
    FLOAT, pointer :: d2HdR2(:, :, :, :) ! second derivative
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
    !%Option displacements 3
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

    integer :: iatom

    this%resp_type = resp_type
    if(.not.varinfo_valid_option('RespPerturbationType', this%resp_type)) call input_error('RespPerturbationType')

    this%dir = -1

    if(this%resp_type == RESP_PERTURBATION_DISPLACE) then
      ! let us create the perturbations
      
      ALLOCATE(this%dHdR(NP_PART, NDIM, geo%natoms), NP_PART*NDIM*geo%natoms)
      ALLOCATE(this%d2HdR2(NP, NDIM, NDIM, geo%natoms), NP*NDIM*NDIM*geo%natoms)
 
      !calculate the derivative of the external potential with respect to the forces
      do iatom = 1, geo%natoms
        !first derivative
        call specie_get_glocal (geo%atom(iatom)%spec, gr, geo%atom(iatom)%x, this%dHdR  (:, :, iatom))
        call specie_get_g2local(geo%atom(iatom)%spec, gr, geo%atom(iatom)%x, this%d2HdR2(:, :, :, iatom))
      
        !the non-local part: tarea para la casa
      end  do
    end if
    
  end subroutine resp_pert_init2


  ! --------------------------------------------------------------------
  subroutine resp_pert_end(this)
    type(resp_pert_t), intent(inout) :: this

    if(this%resp_type == RESP_PERTURBATION_DISPLACE) then
      deallocate(this%dHdR);   nullify(this%dHdR)
      deallocate(this%d2HdR2); nullify(this%d2HdR2)
    end if

  end subroutine resp_pert_end

  ! --------------------------------------------------------------------
  subroutine resp_pert_info(this, unit)
    type(resp_pert_t), intent(in) :: this
    integer,           intent(in) :: unit

    call messages_print_var_option(stdout, 'RespPerturbationType', this%resp_type)
   
  end subroutine resp_pert_info


  ! --------------------------------------------------------------------
  subroutine resp_pert_setup_dir(this, dir)
    type(resp_pert_t), intent(inout) :: this
    integer,           intent(in)    :: dir

    this%dir = dir
  end subroutine resp_pert_setup_dir


#include "undef.F90"
#include "real.F90"
#include "resp_pert_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "resp_pert_inc.F90"

end module resp_pert_m
