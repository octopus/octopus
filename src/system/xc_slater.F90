!! Copyright (C) 2020 N. Tancogne-Dejean
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

module xc_slater_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use XC_F90(lib_m)
  use xc_functl_oct_m

  implicit none

  private
  public ::                     &
    xc_slater_t,                   &
    xc_slater_init,                &
    xc_slater_end,                 &
    dxc_slater_calc,               &
    zxc_slater_calc

  type xc_slater_t
    private
    FLOAT,   pointer,    public :: vxc(:,:)
  end type xc_slater_t

contains

  ! ---------------------------------------------------------
  subroutine xc_slater_init(slater, namespace, family, gr, st)
    type(xc_slater_t),   intent(out)   :: slater
    type(namespace_t),   intent(in)    :: namespace
    integer,             intent(in)    :: family
    type(grid_t),        intent(inout) :: gr
    type(states_elec_t), intent(in)    :: st

    PUSH_SUB(xc_slater_init)

    ! this routine is only prepared for finite systems. (Why not?)
    if(st%d%nik > st%d%ispin) &
      call messages_not_implemented("Slater for periodic systems", namespace=namespace)
    
    ! This variable will keep vxc across iterations
    SAFE_ALLOCATE(slater%vxc(1:gr%mesh%np,st%d%nspin))
    slater%vxc = M_ZERO

    POP_SUB(xc_slater_init)
  end subroutine xc_slater_init


  ! ---------------------------------------------------------
  subroutine xc_slater_end(slater)
    type(xc_slater_t), intent(inout) :: slater

    PUSH_SUB(xc_slater_end)

    SAFE_DEALLOCATE_P(slater%vxc)

    POP_SUB(xc_slater_end)
  end subroutine xc_slater_end

#include "undef.F90"
#include "real.F90"
#include "xc_slater_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_slater_inc.F90"


end module xc_slater_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
