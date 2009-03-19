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
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_basic_m
  use loct_parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use nl_operator_m
  use profiling_m
  use projector_m
  use simul_box_m
  use species_m
  use states_m
  use states_dim_m
  use varinfo_m
  use vibrations_m
  use mpi_m

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

  integer, public, parameter :: &
     PERTURBATION_ELECTRIC = 1, &
     PERTURBATION_MAGNETIC = 2, &
     PERTURBATION_IONIC    = 3, &
     PERTURBATION_KDOTP    = 15, &
     PERTURBATION_NONE     = 0
! kdotp value chosen for compatibility with calculation mode kdotp

  integer, public, parameter :: &
       GAUGE_GIPAW  = 1, &
       GAUGE_ICL    = 2

  type pert_ionic_t
    private
    !if pure_dir is .false. then the perturbation is a combination of
    !displacements of atoms
    logical :: pure_dir
    !in that case these arrays are allocated
    FLOAT, pointer :: mix1(:,:) !mix1(natoms, ndim)
    FLOAT, pointer :: mix2(:,:)
    !in the opposite case, atom, dir, atom2 and dir2 are used
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

  interface pert_init
    module procedure pert_init1
    module procedure pert_init2
  end interface

  type(profile_t), save :: prof

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
    !%Option kdotp 15
    !% Small change of k-point, as in k.p perturbation theory, for calculating
    !% effective masses, and dipoles and polarizabilities in periodic systems
    !%Option none 0
    !% Zero perturbation, for use in testing.
    !%End 
    call loct_parse_int(datasets_check('RespPerturbationType'), PERTURBATION_ELECTRIC, ii)

    call pert_init2(this, ii, gr, geo)

  end subroutine pert_init1


  ! --------------------------------------------------------------------
  subroutine pert_init2(this, pert_type, gr, geo)
    type(pert_t),      intent(out)   :: this
    integer,           intent(in)    :: pert_type
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    call push_sub('pert.pert_init2')
    
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
      !%Option icl 2
      !% ICL correction: Ismail-Beigi, Chang, and Louie, PRL 87, 087402 (2001).
      !%End
      
      call loct_parse_int(datasets_check('MagneticGaugeCorrection'), GAUGE_GIPAW, this%gauge)
      if(.not.varinfo_valid_option('MagneticGaugeCorrection', this%gauge)) &
           call input_error('MagneticGaugeCorrection')

    end if

    if(this%pert_type == PERTURBATION_IONIC) then
      ALLOCATE(this%ionic%mix1(geo%natoms, gr%mesh%sb%dim), geo%natoms*gr%mesh%sb%dim)
      ALLOCATE(this%ionic%mix2(geo%natoms, gr%mesh%sb%dim), geo%natoms*gr%mesh%sb%dim)
    end if

    if(this%pert_type == PERTURBATION_KDOTP) then
      !%Variable KdotP_UseNonLocalPseudopotential
      !%Type logical
      !%Default true
      !%Section Linear Response::KdotP
      !%Description
      !% For testing purposes, set to false to ignore the term -i[r,V] in
      !% the kdotp perturbation, which is due to non-local pseudopotentials.
      !%End

      call loct_parse_logical(datasets_check('KdotP_UseNonLocalPseudopotential'), &
        .true., this%use_nonlocalpps)
    endif

    call pop_sub()

  end subroutine pert_init2

  ! --------------------------------------------------------------------
  subroutine pert_end(this)
    type(pert_t), intent(inout) :: this

    if(this%pert_type == PERTURBATION_IONIC) then
      deallocate(this%ionic%mix1)
      deallocate(this%ionic%mix2)
    end if
  end subroutine pert_end

  ! --------------------------------------------------------------------
  subroutine pert_info(this, unit)
    type(pert_t), intent(in) :: this
    integer,      intent(in) :: unit

    call messages_print_var_option(unit, 'RespPerturbationType', this%pert_type)
    
    if(this%pert_type == PERTURBATION_KDOTP) then
      if (.not. this%use_nonlocalpps) then
        write(message(1), '(a)') 'Ignoring non-local pseudopotential term.'
        call write_info(1)
      endif
    endif
   
  end subroutine pert_info


  ! --------------------------------------------------------------------
  subroutine pert_setup_dir(this, dir, dir2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: dir
    integer, optional, intent(in)    :: dir2

    this%dir = dir
    if(present(dir2)) this%dir2 = dir2

    if(this%pert_type == PERTURBATION_IONIC) call pert_setup_ionic_pure(this)
  end subroutine pert_setup_dir

  subroutine pert_setup_ionic_pure(this)
    type(pert_t), intent(inout) :: this

    this%ionic%pure_dir = .true.

    this%ionic%mix1 = M_ZERO
    this%ionic%mix2 = M_ZERO

    if(this%dir  > 0 .and. this%atom1 > 0) this%ionic%mix1(this%atom1, this%dir ) = M_ONE
    if(this%dir2 > 0 .and. this%atom2 > 0) this%ionic%mix2(this%atom2, this%dir2) = M_ONE

  end subroutine pert_setup_ionic_pure

  subroutine pert_setup_atom(this, iatom, iatom2)
    type(pert_t), intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer, optional, intent(in)    :: iatom2

    this%atom1 = iatom
    if(present(iatom2)) this%atom2 = iatom2

    if(this%pert_type == PERTURBATION_IONIC) call pert_setup_ionic_pure(this)
  end subroutine pert_setup_atom

  subroutine pert_setup_mixed_dir(this, iatom, idir, value, jatom, jdir, valuej)
    type(pert_t),      intent(inout) :: this
    integer,           intent(in)    :: iatom
    integer,           intent(in)    :: idir
    FLOAT,             intent(in)    :: value
    integer, optional, intent(in)    :: jatom
    integer, optional, intent(in)    :: jdir
    FLOAT,   optional, intent(in)    :: valuej
    

    logical :: have_dir_2

    this%ionic%pure_dir = .false.
    
    this%ionic%mix1(iatom, idir) = value
    
    have_dir_2 = present(jatom) .and. present(jdir) .and. present(jatom)

    if(have_dir_2) then
      this%ionic%mix1(jatom, jdir) = valuej
    else
      ASSERT( .not. present(jatom) .and. .not. present(jdir) .and. .not. present(jatom) )
    end if

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
