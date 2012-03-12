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

module projector_m
  use atom_m
  use batch_m
  use comm_m
  use double_grid_m
  use geometry_m
  use global_m
  use grid_m
  use hgh_projector_m
  use io_m
  use kb_projector_m
  use kpoints_m
  use lalg_basic_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use multicomm_m
  use profiling_m
  use ps_m
  use rkb_projector_m
  use simul_box_m
  use species_m
  use states_m
  use states_dim_m
  use submesh_m
  use varinfo_m

  implicit none

  private
  public ::                       &
       projector_t,               &
       projector_null,            &
       projector_is_null,         &
       projector_is,              &
       projector_init,            &
       projector_init_phases,     &
       projector_build,           &
       projector_end,             &
       dproject_psi,              &
       zproject_psi,              &
       dproject_psi_batch,        &
       zproject_psi_batch,        &
       dproject_sphere,           &
       zproject_sphere,           &
       dprojector_matrix_element, &
       zprojector_matrix_element, &
       dprojector_commute_r,      &
       zprojector_commute_r

  integer, public, parameter ::  &
       M_NONE = 0,  &
       M_HGH  = 1,  &
       M_KB   = 2,  &
       M_RKB  = 3

  integer, parameter :: MAX_NPROJECTIONS = 24
  integer, parameter :: MAX_L = 5

  !> The projector data type is intended to hold the local and
  !! non-local parts of the pseudopotentials. The definition of the
  !! action of a projector (which is done through the X(project)
  !! subroutine) depends on the type of the projector. 
  !!
  !! There are four different types: 
  !! - local -> a local operator
  !! - HGH projector -> "normal"
  !! - Kleinman-Bylander projector (no spin-orbit) -> "relativistic"
  !! - Kleinman-Bylander projector (includes spin-orbit)

  type projector_t
    integer :: type = M_NONE
    integer :: nprojections
    integer :: lmax
    integer :: lloc
    integer :: nik
    integer :: reltype

    type(submesh_t)  :: sphere
    

    !> Only one of the following structures should be used at once
    !! The one to be used depends on the value of type variable
    type(hgh_projector_t), pointer :: hgh_p(:, :) => null()
    type(kb_projector_t),  pointer :: kb_p(:, :)  => null()
    type(rkb_projector_t), pointer :: rkb_p(:, :) => null()
    CMPLX,                 pointer :: phase(:, :) => null()
  end type projector_t

contains
  ! ---------------------------------------------------------
  subroutine projector_null(p)
    type(projector_t), intent(out) :: p

    PUSH_SUB(projector_null)

    p%type = M_NONE
    call submesh_null(p%sphere)

    POP_SUB(projector_null)

  end subroutine projector_null
  !---------------------------------------------------------


  !---------------------------------------------------------
  logical elemental function projector_is_null(p)
    type(projector_t), intent(in) :: p
    projector_is_null = (p%type == M_NONE)
  end function projector_is_null
  !---------------------------------------------------------

  logical elemental function projector_is(p, type)
    type(projector_t), intent(in) :: p
    integer,           intent(in) :: type
    projector_is = (p%type == type)
  end function projector_is

  !---------------------------------------------------------
  subroutine projector_init(p, mesh, atm, dim, reltype)
    type(projector_t), intent(inout) :: p
    type(mesh_t),      intent(in)    :: mesh
    type(atom_t),      intent(in)    :: atm
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: reltype

    type(ps_t), pointer :: ps
    
    PUSH_SUB(projector_init)

    call projector_null(p)
    ps => species_ps(atm%spec)

    nullify(p%phase)
    p%reltype = reltype
    p%lmax = ps%l_max

    if(p%lmax == 0) then
      p%type = M_NONE
      POP_SUB(projector_init)
      return
    end if

    p%lloc = ps%l_loc

    call submesh_init_sphere(p%sphere, mesh%sb, mesh, atm%x, ps%rc_max + mesh%spacing(1))

    select case (ps%kbc)
    case (1)
      p%type = M_KB
      if (reltype == 1) then
        write(message(1),'(a,a,a)') &
          "Spin-orbit coupling for species ", trim(species_label(atm%spec)), " is not available."
        call messages_warning(1)
      end if
    case (2)
      if (reltype == 0) then
        p%type = M_KB
      else
        p%type = M_RKB
      end if
    case (3)
      p%type = M_HGH
    end select
    
    select case(p%type)
    case(M_KB)
      p%nprojections = dim * ps%kbc
    case(M_RKB)
      p%nprojections = 4
    case(M_HGH)
      p%nprojections = 24
    end select

    POP_SUB(projector_init)
  end subroutine projector_init

  subroutine projector_init_phases(this, sb, std, vec_pot, vec_pot_var)
    type(projector_t),  intent(inout) :: this
    type(simul_box_t),  intent(in)    :: sb
    type(states_dim_t), intent(in)    :: std
    FLOAT, optional,    pointer       :: vec_pot(:)
    FLOAT, optional,    pointer       :: vec_pot_var(:, :)

    integer :: ns, iq, is, ikpoint
    FLOAT   :: kr, kpoint(1:MAX_DIM)
    integer :: ndim

    PUSH_SUB(projector_init_phases)

    ns = this%sphere%np
    ndim = sb%dim

    if(.not. associated(this%phase) .and. ns > 0) then
      SAFE_ALLOCATE(this%phase(1:ns, std%kpt%start:std%kpt%end))
    end if

    do iq = std%kpt%start, std%kpt%end
      ikpoint = states_dim_get_kpoint_index(std, iq)

      ! if this fails, it probably means that sb is not compatible with std
      ASSERT(ikpoint <= kpoints_number(sb%kpoints))
      
      kpoint = M_ZERO
      kpoint(1:ndim) = kpoints_get_point(sb%kpoints, ikpoint)
        
      do is = 1, ns
        ! this is only the correction to the global phase, that can
        ! appear if the sphere crossed the boundary of the cell.
        
        kr = sum(kpoint(1:ndim)*(this%sphere%x(is, 1:ndim) - this%sphere%mesh%x(this%sphere%map(is), 1:ndim)))

        if(present(vec_pot)) then
          if(associated(vec_pot)) kr = kr + &
            sum(vec_pot(1:ndim)*(this%sphere%x(is, 1:ndim)- this%sphere%mesh%x(this%sphere%map(is), 1:ndim)))
        end if

        if(present(vec_pot_var)) then
          if(associated(vec_pot_var)) kr = kr + sum(vec_pot_var(1:ndim, this%sphere%map(is))*this%sphere%x(is, 1:ndim))
        end if

        this%phase(is, iq) = exp(-M_zI*kr)
      end do

    end do

    POP_SUB(projector_init_phases)

  end subroutine projector_init_phases

  !---------------------------------------------------------
  subroutine projector_build(p, gr, a, so_strength)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(atom_t),      intent(in)    :: a
    FLOAT,             intent(in)    :: so_strength

    integer :: ll, mm

    PUSH_SUB(projector_build)

    select case (p%type)

    case (M_HGH)
      SAFE_ALLOCATE(p%hgh_p(0:p%lmax, -p%lmax:p%lmax))
      do ll = 0, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call hgh_projector_null(p%hgh_p(ll, mm))
          call hgh_projector_init(p%hgh_p(ll, mm), p%sphere, gr, a, ll, mm, so_strength)
        end do
      end do

    case (M_KB)
      SAFE_ALLOCATE(p%kb_p(0:p%lmax, -p%lmax:p%lmax))
      do ll = 0, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call kb_projector_null(p%kb_p(ll, mm))
          call kb_projector_init(p%kb_p(ll, mm), p%sphere, gr, a, ll, mm)
        end do
      end do

    case (M_RKB)
      SAFE_ALLOCATE(p%rkb_p(1:p%lmax, -p%lmax:p%lmax))
      do ll = 1, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call rkb_projector_null(p%rkb_p(ll, mm))
          call rkb_projector_init(p%rkb_p(ll, mm), p%sphere, a, ll, mm, so_strength)
        end do
      end do
      ! for rkb, l = 0 is a normal kb
      if(p%lloc /= 0) then
        SAFE_ALLOCATE(p%kb_p(1:1, 1:1))
        call kb_projector_null(p%kb_p(1, 1))
        call kb_projector_init(p%kb_p(1, 1), p%sphere, gr, a, 0, 0)
      end if

    end select

    POP_SUB(projector_build)
  end subroutine projector_build

  !---------------------------------------------------------
  subroutine projector_end(p)
    type(projector_t), intent(inout) :: p

    integer :: ll, mm

    PUSH_SUB(projector_end)

    call submesh_end(p%sphere)

    select case(p%type)
    case(M_HGH)
      do ll = 0, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call hgh_projector_end(p%hgh_p(ll, mm))
        end do
      end do
      SAFE_DEALLOCATE_P(p%hgh_p)

    case(M_KB)
      do ll = 0, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call kb_projector_end(p%kb_p(ll, mm))
        end do
      end do
      SAFE_DEALLOCATE_P(p%kb_p)

    case(M_RKB)
      do ll = 1, p%lmax
        if(ll == p%lloc) cycle
        do mm = -ll, ll
          call rkb_projector_end(p%rkb_p(ll, mm))
        end do
      end do
      SAFE_DEALLOCATE_P(p%rkb_p)
      if(p%lloc /= 0) then
        call kb_projector_end(p%kb_p(1, 1))
        SAFE_DEALLOCATE_P(p%kb_p)
      end if

    end select
    
    p%type = M_NONE

    SAFE_DEALLOCATE_P(p%phase)

    POP_SUB(projector_end)
  end subroutine projector_end

#include "undef.F90"
#include "real.F90"
#include "projector_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "projector_inc.F90"

end module projector_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
