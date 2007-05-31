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
       local_t,                   &
       projector_t,               &
       projector_null,            &
       projector_init,            &
       projector_build,           &
#ifdef HAVE_MPI
       projector_broadcast,       &
#endif
       projector_end,             &
       projector_build_kb_sphere, &
       projector_copy_kb_sphere,  &
       dproject, zproject,        &
       dpsidprojectpsi,           &
       zpsidprojectpsi,           &
       dpsia_project_psib,        &
       zpsia_project_psib

  integer, public, parameter ::  &
       M_HGH = 1, &
       M_KB  = 2, &
       M_RKB = 3, &
       M_LOCAL = 4

  ! an oximoronic local projector type
  type local_t
     FLOAT, pointer :: v(:)
  end type local_t

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
    CMPLX,   pointer :: phases(:,:)

    ! Only one of the following structures should be used at once
    ! The one to be used depends on the value of type variable
    type(hgh_projector_t)  :: hgh_p
    type(kb_projector_t)   :: kb_p
    type(rkb_projector_t)  :: rkb_p
    type(local_t)          :: local_p
  end type projector_t

contains
  ! ---------------------------------------------------------
  subroutine projector_null(p)
    type(projector_t), intent(out) :: p

    p%type = 0
    nullify(p%phases)
    call submesh_null(p%sphere)

  end subroutine projector_null

  !---------------------------------------------------------
  subroutine projector_build_kb_sphere(p, sb, m, a, st)
    type(projector_t), intent(inout) :: p
    type(simul_box_t), intent(in)    :: sb
    type(mesh_t),      intent(in)    :: m
    type(atom_t),      intent(in)    :: a
    type(states_t),    intent(in)    :: st

    integer :: j, k
    FLOAT :: x(MAX_DIM)

    call push_sub('projector.projector_build_kb_sphere')

    call submesh_init_sphere(p%sphere, sb, m, a%x, a%spec%ps%rc_max)

    ! and here the phases for the periodic systems
    if (sb%periodic_dim /= 0) then
      p%nik = st%d%nik

      ALLOCATE(p%phases(p%sphere%ns, p%nik), p%sphere%ns*p%nik)

      do j = 1, p%sphere%ns
        x(:) = m%x(p%sphere%jxyz(j), :)
        do k = 1, p%nik
          p%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
        end do
      end do

    end if

    call pop_sub()
  end subroutine projector_build_kb_sphere

  !---------------------------------------------------------
  subroutine projector_copy_kb_sphere(pi, po)
    type(projector_t), intent(in)    :: pi
    type(projector_t), intent(inout) :: po

    call push_sub('projector.projector_copy_kb_sphere')

    call submesh_copy(pi%sphere, po%sphere)

    if (associated(pi%phases)) then
      po%nik = pi%nik
      ALLOCATE(po%phases(pi%sphere%ns, pi%nik), pi%sphere%ns*pi%nik)
      po%phases = pi%phases
    end if

    call pop_sub()
  end subroutine projector_copy_kb_sphere

  !---------------------------------------------------------
  subroutine projector_init(p, a, reltype, l, lm, force_type)
    type(projector_t), intent(inout) :: p
    type(atom_t),      intent(in)    :: a
    integer, optional, intent(in)    :: reltype
    integer, optional, intent(in)    :: l, lm
    integer, optional, intent(in)    :: force_type

    call push_sub('projector.projector_init')

    if(present(l)) p%l = l
    if(present(lm)) p%lm = lm

    if(present(force_type)) then 
      p%type = force_type
    else 
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
    end if

    call pop_sub()
  end subroutine projector_init

  !---------------------------------------------------------
  subroutine projector_build(p, gr, a, gen_grads)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(atom_t),      intent(in)    :: a
    logical,           intent(in)    :: gen_grads

    integer :: ns, l, lm

    call push_sub('projector.projector_build')

    l = p%l
    lm = p%lm

    select case (p%type)

    case(M_LOCAL)
      ns =  p%sphere%ns
      ALLOCATE(p%local_p%v(1:ns), ns)
      call double_grid_apply_local(gr%dgrid, a%spec, gr%m, p%sphere, a%x, p%local_p%v)

    case (M_HGH)
      call hgh_projector_null(p%hgh_p)
      call hgh_projector_init(p%hgh_p, p%sphere, gr, a, l, lm)

    case (M_KB)
      call kb_projector_null(p%kb_p)
      call kb_projector_init(p%kb_p, p%sphere, gr, a, l, lm, gen_grads)

    case (M_RKB)
      call rkb_projector_null(p%rkb_p)
      call rkb_projector_init(p%rkb_p, p%sphere, gr, a, l, lm)

    end select

    call pop_sub()
  end subroutine projector_build

#ifdef HAVE_MPI
  !---------------------------------------------------------
  subroutine projector_broadcast(p, gr, mc, a, gen_grads, root)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(multicomm_t), intent(in)    :: mc
    type(atom_t),      intent(in)    :: a
    logical,           intent(in)    :: gen_grads
    integer,           intent(in)    :: root

    integer :: ns, l, lm, mpi_err, rank

    call push_sub('projector.projector_broadcast')
    
    l = p%l
    lm = p%lm
    rank = mc%who_am_i(P_STRATEGY_STATES) - 1

    select case (p%type)
    case(M_LOCAL)
      ns =  p%sphere%ns
      if ( rank /= root) then
        ALLOCATE(p%local_p%v(1:ns), ns)
      end if
      call MPI_Bcast(p%local_p%v, ns, MPI_FLOAT, root, mc%group_comm(P_STRATEGY_STATES), mpi_err)

    case (M_HGH)
      if ( rank /= root) then      
        call hgh_projector_null(p%hgh_p)
        call hgh_projector_init(p%hgh_p, p%sphere, gr, a, l, lm)
      end if
    case (M_KB)
      if ( rank /= root) then
        call kb_projector_null(p%kb_p)
      end if
      call kb_projector_broadcast(p%kb_p, p%sphere, gr, mc, a, l, lm, gen_grads, root)
    case (M_RKB)
      if ( rank /= root) then
        call rkb_projector_null(p%rkb_p)
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
    if (associated(p%phases)) deallocate(p%phases)

    select case(p%type)
    case(M_LOCAL)
      deallocate(p%local_p%v)
    case(M_HGH)
      call hgh_projector_end(p%hgh_p)
    case(M_KB)
      call kb_projector_end(p%kb_p)
    case(M_RKB)
      call rkb_projector_end(p%rkb_p)
    end select
    
    p%type = 0

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
