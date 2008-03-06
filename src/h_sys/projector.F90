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
  use double_grid_m
  use global_m
  use grid_m
  use io_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use simul_box_m
  use ps_m
  use specie_m
  use specie_pot_m
  use geometry_m
  use states_m
  use submesh_m
  use mpi_m
  use mpi_debug_m
  use multicomm_m
  use varinfo_m
  use hgh_projector_m
  use kb_projector_m
  use rkb_projector_m

  implicit none

  private
  public ::                       &
       projector_t,               &
       projector_null,            &
       projector_init,            &
       projector_init_phases,     &
       projector_build,           &
#ifdef HAVE_MPI
       projector_broadcast,       &
#endif
       projector_end,             &
       dproject_psi,              &
       zproject_psi,              &
       dproject_sphere,           &
       zproject_sphere,           &
       dpsia_project_psib,        &
       zpsia_project_psib

  integer, public, parameter ::  &
       M_HGH = 1, &
       M_KB  = 2, &
       M_RKB = 3

  ! The projector data type is intended to hold the local and
  ! non-local parts of the pseudopotentials. The definition of the
  ! action of a projector (which is done through the X(project)
  ! subroutine) depends on the type of the projector. 

  ! There are four different types: 
  ! local -> a local operator
  ! HGH projector -> "normal"
  ! Kleinman-Bylander projector (no spin-orbit) -> "relativistc"
  ! Kleinman-Bylander projector (includes spin-orbit)

  type projector_t
    integer :: type
    integer :: iatom
    integer :: l
    integer :: lm

    type(submesh_t)  :: sphere

    integer          :: nik

    ! Only one of the following structures should be used at once
    ! The one to be used depends on the value of type variable
    type(hgh_projector_t)  :: hgh_p
    type(kb_projector_t)   :: kb_p
    type(rkb_projector_t)  :: rkb_p
    CMPLX, pointer         :: phase(:, :)
  end type projector_t

contains
  ! ---------------------------------------------------------
  subroutine projector_null(p)
    type(projector_t), intent(out) :: p

    p%type = 0
    nullify(p%phase)
    call submesh_null(p%sphere)

  end subroutine projector_null

  !---------------------------------------------------------
  subroutine projector_init(p, a, reltype, l, lm)
    type(projector_t), intent(inout) :: p
    type(atom_t),      intent(in)    :: a
    integer,           intent(in)    :: reltype
    integer,           intent(in)    :: l, lm

    call push_sub('projector.projector_init')

    nullify(p%phase)
    p%l = l
    p%lm = lm

    select case (a%spec%ps%kbc)
    case (1)
      p%type = M_KB
      if (reltype == 1) then
        write(message(1),'(a,a,a)') "Spin-orbit coupling for specie ", trim(a%spec%label), " is not available."
        call write_warning(1)
      end if
    case (2)
      if (l == 0 .or. reltype == 0) then
        p%type = M_KB
      else
        p%type = M_RKB
      end if
    case (3)
      p%type = M_HGH
    end select
    
    call pop_sub()
  end subroutine projector_init

  subroutine projector_init_phases(this, m, nkpoints, kpoints)
    type(projector_t), intent(inout) :: this
    type(mesh_t),      intent(in)    :: m
    integer,           intent(in)    :: nkpoints
    FLOAT,             intent(in)    :: kpoints(:, :)

    integer :: ns, ik, is
    FLOAT   :: kr

    ns = this%sphere%ns

    ASSERT(.not. associated(this%phase))
    ALLOCATE(this%phase(1:ns, 1:nkpoints), ns*nkpoints)

    do ik = 1, nkpoints
      do is = 1, ns
        kr = sum(kpoints(1:MAX_DIM, ik)*(this%sphere%x(is, 1:MAX_DIM)-m%x(this%sphere%jxyz(is), 1:MAX_DIM)))
        this%phase(is, ik) = exp(-M_zI*kr)
      end do
    end do

  end subroutine projector_init_phases

  !---------------------------------------------------------
  subroutine projector_build(p, gr, a)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(atom_t),      intent(in)    :: a

    integer :: l, lm

    call push_sub('projector.projector_build')

    l = p%l
    lm = p%lm

    select case (p%type)

    case (M_HGH)
      call hgh_projector_null(p%hgh_p)
      call hgh_projector_init(p%hgh_p, p%sphere, gr, a, l, lm)

    case (M_KB)
      call kb_projector_null(p%kb_p)
      call kb_projector_init(p%kb_p, p%sphere, gr, a, l, lm)

    case (M_RKB)
      call rkb_projector_null(p%rkb_p)
      call rkb_projector_init(p%rkb_p, p%sphere, gr, a, l, lm)

    end select

    call pop_sub()
  end subroutine projector_build

#ifdef HAVE_MPI
  !---------------------------------------------------------
  subroutine projector_broadcast(p, gr, mc, a, root)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    type(atom_t),      intent(in)    :: a
    integer,           intent(in)    :: root

    integer :: l, lm, rank

    call push_sub('projector.projector_broadcast')
    
    l = p%l
    lm = p%lm
    rank = mc%who_am_i(P_STRATEGY_STATES)

    select case (p%type)
    case (M_HGH)
      if ( rank /= root) then      
        call hgh_projector_null(p%hgh_p)
        !for the moment it is not copied but initialized
        call hgh_projector_init(p%hgh_p, p%sphere, gr, a, l, lm)
      end if
    case (M_KB)
      if ( rank /= root) then
        call kb_projector_null(p%kb_p)
      end if
      call kb_projector_broadcast(p%kb_p, p%sphere, gr, mc, a, l, lm, root)
    case (M_RKB)
      if ( rank /= root) then
        call rkb_projector_null(p%rkb_p)
        !for the moment it is not copied but initialized
        call rkb_projector_init(p%rkb_p, p%sphere, gr, a, l, lm)
      end if
    end select

    call pop_sub()
  end subroutine projector_broadcast
#endif

  !---------------------------------------------------------
  subroutine projector_end(p)
    type(projector_t), intent(inout) :: p

    call push_sub('projector.projector_end')

    call submesh_end(p%sphere)

    select case(p%type)
    case(M_HGH)
      call hgh_projector_end(p%hgh_p)
    case(M_KB)
      call kb_projector_end(p%kb_p)
    case(M_RKB)
      call rkb_projector_end(p%rkb_p)
    end select
    
    p%type = 0

    if(associated(p%phase)) deallocate(p%phase)

    call pop_sub()
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
