!! Copyright (C) 2016 X. Andrade, N. Tancogne-Dejean
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

module scissor_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use global_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                  &
       scissor_t,            &
       scissor_init,         &
       dscissor_apply,       &
       zscissor_apply,       &
       scissor_commute_r,    &
       scissor_end

  type scissor_t
    private
    logical, public         :: apply = .false.
    FLOAT                   :: gap
    type(states_elec_t)     :: gs_st
  end type scissor_t
 
  interface scissor_commute_r
    module procedure dscissor_commute_r, zscissor_commute_r
  end interface

contains

 subroutine scissor_init(this, namespace, space, st, mesh, d, kpoints, gap, mc)
  type(scissor_t),          intent(inout) :: this
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(states_elec_t),      intent(in)    :: st
  type(mesh_t),             intent(in)    :: mesh
  type(kpoints_t),          intent(in)    :: kpoints
  type(states_elec_dim_t),  intent(in)    :: d
  FLOAT,                    intent(in)    :: gap
  type(multicomm_t),        intent(in)    :: mc

  CMPLX, allocatable   :: phase(:)
  type(restart_t) :: restart_gs
  integer :: ierr
  integer :: ist, ik, ip
  CMPLX, allocatable :: temp_state(:,:)
  FLOAT   :: kpoint(1:MAX_DIM)

  PUSH_SUB(scissor_init)

  ASSERT(.not. this%apply)
  ASSERT(.not. states_are_real(st))

  call messages_print_stress(stdout, "TDScissor", namespace=namespace)

  if(st%parallel_in_states) call messages_not_implemented("Scissor operator parallel in states", namespace=namespace)
  if(mesh%parallel_in_domains) call messages_not_implemented("Scissor operator parallel in domains", namespace=namespace)

  this%apply = .true.
  this%gap = gap

  write(message(1),'(a)')    'Start loading GS states.'
  call messages_info(1) 
  !We need to load GS states and to store them in this%gs_st
  call states_elec_copy(this%gs_st, st)
  
  call restart_init(restart_gs, namespace, RESTART_PROJ, RESTART_TYPE_LOAD, mc, ierr, mesh=mesh)
  if(ierr /= 0) then
     message(1) = "Unable to read states information."
     call messages_fatal(1, namespace=namespace)
  end if

  call states_elec_load(restart_gs, namespace, space, this%gs_st, mesh, kpoints, ierr, label = ': gs for TDScissor')
  if(ierr /= 0 .and. ierr /= (this%gs_st%st_end-this%gs_st%st_start+1)*this%gs_st%d%nik*this%gs_st%d%dim) then
    message(1) = "Unable to read wavefunctions for TDScissor."
    call messages_fatal(1, namespace=namespace)
  end if
  call restart_end(restart_gs)

  if (space%is_periodic() .and. .not. (kpoints_number(kpoints) == 1 .and. kpoints_point_is_gamma(kpoints, 1))) then

    write(message(1),'(a)')    'Adding the phase for GS states.'
    call messages_info(1)
  
    SAFE_ALLOCATE(temp_state(1:mesh%np_part, 1:this%gs_st%d%dim))
    SAFE_ALLOCATE(phase(1:mesh%np))
    ! We apply the phase to these states, as we need it for the projectors later
    do ik=this%gs_st%d%kpt%start, this%gs_st%d%kpt%end

      kpoint(1:space%dim) = kpoints%get_point(d%get_kpoint_index(ik))
      do ip = 1, mesh%np
        phase(ip) = exp(-M_zI * sum(mesh%x(ip, 1:space%dim) * kpoint(1:space%dim)))
      end do

      do ist = this%gs_st%st_start, this%gs_st%st_end
        call states_elec_get_state(this%gs_st, mesh, ist, ik, temp_state )
        call states_elec_set_phase(this%gs_st%d, temp_state, phase, mesh%np, .false.)
        call states_elec_set_state(this%gs_st, mesh, ist, ik,temp_state )
      end do
   
    end do

    SAFE_DEALLOCATE_A(temp_state) 
    SAFE_DEALLOCATE_A(phase)
  end if

  call messages_print_stress(stdout, namespace=namespace)

  POP_SUB(scissor_init)
 end subroutine scissor_init

 subroutine scissor_end(this)
   type(scissor_t), intent(inout) :: this

   PUSH_SUB(scissor_end)

   this%apply = .false.
   call states_elec_end(this%gs_st)

   POP_SUB(scissor_end)
 end subroutine scissor_end

#include "undef.F90"
#include "real.F90"
#include "scissor_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "scissor_inc.F90"
end module scissor_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

