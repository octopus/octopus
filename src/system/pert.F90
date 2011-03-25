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
  use derivatives_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_basic_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use parser_m
  use physics_op_m
  use poisson_m
  use profiling_m
  use projector_m
  use simul_box_m
  use species_m
  use states_m
  use states_dim_m
  use varinfo_m
  use vibrations_m

  implicit none

  private
  public ::                         &
     pert_t,                        &
     pert_init,                     &
     pert_end,                      &
     pert_info,                     &
     pert_setup_dir,                &
     pert_setup_atom,               &
     pert_setup_mixed_dir,          &
     pert_type,                     &
     dpert_apply,                   &
     zpert_apply,                   &
     dpert_apply_order_2,           &
     zpert_apply_order_2,           &
     dpert_expectation_value,       &
     zpert_expectation_value,       &
     dpert_expectation_density,     &
     zpert_expectation_density,     &
     dionic_pert_matrix_elements_2, &
     zionic_pert_matrix_elements_2

  integer, public, parameter ::  &
    PERTURBATION_ELECTRIC = 1,  &
    PERTURBATION_MAGNETIC = 2,  &
    PERTURBATION_IONIC    = 3,  &
    PERTURBATION_KDOTP    = 15, &   ! value chosen for compatibility with calculation mode kdotp
    PERTURBATION_NONE     = 0

  integer, public, parameter :: &
    GAUGE_GIPAW  = 1, &
    GAUGE_ICL    = 2

  type pert_ionic_t
    private
    !if pure_dir is .false. then the perturbation is a combination of
    !displacements of atoms
    logical :: pure_dir
    !if pure_dir is .true., these arrays are allocated
    FLOAT, pointer :: mix1(:,:) !mix1(natoms, ndim)
    FLOAT, pointer :: mix2(:,:)
    !if pure_dir is .false., atom, dir, atom2 and dir2 are used
  end type pert_ionic_t

  type pert_t
    private
    integer            :: pert_type
    integer            :: dir
    integer            :: dir2
    integer            :: atom1, atom2
    integer            :: gauge
    type(pert_ionic_t) :: ionic
    logical            :: use_nonlocalpps
  end type pert_t

  type(profile_t), save :: prof

contains

  ! --------------------------------------------------------------------
  subroutine pert_init(this, pert_type, gr, geo)
    type(pert_t),      intent(out)   :: this
    integer,           intent(in)    :: pert_type
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    PUSH_SUB(pert_init)
    
    this%pert_type = pert_type

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
      !% For magnetic linear response: how to handle gauge-invariance in the description
      !% of the coupling of electrons to the magnetic field.
      !%Option none 0
      !% No correction.
      !%Option gipaw 1
      !% GIPAW correction: C Pickard and F Mauri, <i>Phys. Rev. Lett.</i> <b>91</b>, 196401 (2003).
      !%Option icl 2
      !% ICL correction: S Ismail-Beigi, EK Chang, and SG Louie, <i>Phys. Rev. Lett.</i> <b>87</b>, 087402 (2001).
      !%End
      
      call parse_integer(datasets_check('MagneticGaugeCorrection'), GAUGE_GIPAW, this%gauge)
      if(.not.varinfo_valid_option('MagneticGaugeCorrection', this%gauge)) &
           call input_error('MagneticGaugeCorrection')

    end if

    if(this%pert_type == PERTURBATION_IONIC) then
      SAFE_ALLOCATE(this%ionic%mix1(1:geo%natoms, 1:gr%mesh%sb%dim))
      SAFE_ALLOCATE(this%ionic%mix2(1:geo%natoms, 1:gr%mesh%sb%dim))
    end if

    if(this%pert_type == PERTURBATION_KDOTP) then
      !%Variable KdotP_UseNonLocalPseudopotential
      !%Type logical
      !%Default true
      !%Section Linear Response::KdotP
      !%Description
      !% For testing purposes, set to false to ignore the term -i[r,V] in
      !% the k.p perturbation, which is due to non-local pseudopotentials.
      !%End

      call parse_logical(datasets_check('KdotP_UseNonLocalPseudopotential'), &
        .true., this%use_nonlocalpps)
    endif

    POP_SUB(pert_init)
  end subroutine pert_init

  ! --------------------------------------------------------------------
  subroutine pert_end(this)
    type(pert_t), intent(inout) :: this

    PUSH_SUB(pert_end)

    if(this%pert_type == PERTURBATION_IONIC) then
      SAFE_DEALLOCATE_P(this%ionic%mix1)
      SAFE_DEALLOCATE_P(this%ionic%mix2)
    end if

    POP_SUB(pert_end)
  end subroutine pert_end

  ! --------------------------------------------------------------------
  subroutine pert_info(this, unit)
    type(pert_t), intent(in) :: this
    integer,      intent(in) :: unit

    PUSH_SUB(pert_info)

    call messages_print_var_option(unit, 'RespPerturbationType', this%pert_type)
    
    if(this%pert_type == PERTURBATION_KDOTP) then
      if (.not. this%use_nonlocalpps) then
        write(message(1), '(a)') 'Ignoring non-local pseudopotential term.'
        call messages_info(1)
      endif
    endif
   
    POP_SUB(pert_info)
  end subroutine pert_info

  ! --------------------------------------------------------------------
  subroutine pert_setup_dir(this, dir, dir2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: dir
    integer, optional, intent(in)    :: dir2

    PUSH_SUB(pert_setup_dir)

    this%dir = dir
    if(present(dir2)) this%dir2 = dir2

    if(this%pert_type == PERTURBATION_IONIC) call pert_setup_ionic_pure(this)

    POP_SUB(pert_setup_dir)
  end subroutine pert_setup_dir

  ! --------------------------------------------------------------------
  subroutine pert_setup_ionic_pure(this)
    type(pert_t), intent(inout) :: this

    PUSH_SUB(pert_setup_ionic_pure)

    this%ionic%pure_dir = .true.

    this%ionic%mix1 = M_ZERO
    this%ionic%mix2 = M_ZERO

    if(this%dir  > 0 .and. this%atom1 > 0) this%ionic%mix1(this%atom1, this%dir ) = M_ONE
    if(this%dir2 > 0 .and. this%atom2 > 0) this%ionic%mix2(this%atom2, this%dir2) = M_ONE

    POP_SUB(pert_setup_ionic_pure)
  end subroutine pert_setup_ionic_pure

  ! --------------------------------------------------------------------
  subroutine pert_setup_atom(this, iatom, iatom2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer, optional, intent(in)    :: iatom2

    PUSH_SUB(pert_setup_atom)

    this%atom1 = iatom
    if(present(iatom2)) this%atom2 = iatom2

    if(this%pert_type == PERTURBATION_IONIC) call pert_setup_ionic_pure(this)

    POP_SUB(pert_setup_atom)
  end subroutine pert_setup_atom

  ! --------------------------------------------------------------------
  subroutine pert_setup_mixed_dir(this, iatom, idir, value, jatom, jdir, valuej)
    type(pert_t),      intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: idir
    FLOAT,             intent(in)    :: value
    integer, optional, intent(in)    :: jatom
    integer, optional, intent(in)    :: jdir
    FLOAT,   optional, intent(in)    :: valuej
    
    logical :: have_dir_2

    PUSH_SUB(pert_setup_mixed_dir)

    this%ionic%pure_dir = .false.
    
    this%ionic%mix1(iatom, idir) = value
    
    have_dir_2 = present(jatom) .and. present(jdir) .and. present(jatom)

    if(have_dir_2) then
      this%ionic%mix1(jatom, jdir) = valuej
    else
      ASSERT( .not. present(jatom) .and. .not. present(jdir) .and. .not. present(jatom) )
    end if

    POP_SUB(pert_setup_mixed_dir)
  end subroutine pert_setup_mixed_dir
    
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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
