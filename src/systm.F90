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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module system_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use mesh_m
  use simul_box_m
  use grid_m
  use geometry_m
  use states_m
  use v_ks_m
  use hamiltonian_m
  use output_m
  use multicomm_m
  use mpi_m
  use poisson_m
  use units_m

  implicit none

  private
  public ::               &
    system_t,             &
    system_init,          &
    system_end,           &
    system_h_setup

  type system_t
    type(geometry_t)            :: geo
    type(grid_t),     pointer :: gr    ! the mesh
    type(states_t),   pointer :: st    ! the states
    type(v_ks_t)              :: ks    ! the Kohn-Sham potentials
    type(output_t)            :: outp  ! the output
    type(multicomm_t)         :: mc    ! index and domain communicators
  end type system_t

contains

  !----------------------------------------------------------
  subroutine system_init(sys, parallel_mask)
    type(system_t), intent(out) :: sys
    integer, intent(in)         :: parallel_mask

    call push_sub('systm.system_init')

    ALLOCATE(sys%gr, 1)
    ALLOCATE(sys%st, 1)

    call geometry_init(sys%geo)
    call simul_box_init(sys%gr%sb, sys%geo)
    call states_init(sys%st, sys%gr, sys%geo)
    call grid_init_stage_1(sys%gr, sys%geo)

    call parallel_init()

    call grid_init_stage_2(sys%gr, sys%mc, sys%geo)
    call states_densities_init(sys%st, sys%gr, sys%geo)
    call output_init(sys%gr%sb, sys%outp)
    call poisson_init(sys%gr, sys%geo)
    call v_ks_init(sys%gr, sys%ks, sys%st%d)

    !print the mesh information if it is required
    call print_r()

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine parallel_init()
      integer :: index_dim, index_range(4)

      ! for the moment we need communicators for domains plus three
      ! index-dimensions: states, kpoints, and others
      index_dim = 4

      ! store the ranges for these two indices (serves as initial guess
      ! for parallelization strategy)
      index_range(1) = sys%gr%m%np_global     ! Number of points in mesh
      index_range(2) = sys%st%nst             ! Number of states
      index_range(3) = sys%st%d%nik           ! Number of kpoints
      index_range(4) = 100000                 ! Some large number

      ! create index and domain communicators
      call multicomm_init(sys%mc, parallel_mask, mpi_world%size, index_dim, &
         index_range, (/ 5000, 1, 1, 1 /))

    end subroutine parallel_init

    subroutine print_r()
      FLOAT :: u
      integer :: i, ierr
      character(len=80) :: fname
      
      
      if(iand(sys%outp%what, output_r).ne.0) then
        u = M_ONE/units_out%length%factor
        
        do i=1, sys%NDIM
          write(fname, '(a,i1)') 'r-', i
          call doutput_function(sys%outp%how, 'status/', fname, sys%gr%m, sys%gr%sb, sys%gr%m%x(:,i), u, ierr)
        end do
        
      end if
      
    end subroutine print_r
    
  end subroutine system_init


  !----------------------------------------------------------
  subroutine system_end(s)
    type(system_t), intent(inout) :: s

    call push_sub('systm.system_end')

#if defined(HAVE_MPI)
    call multicomm_end(s%mc)
#endif

    call v_ks_end(s%ks)

    if(associated(s%st)) then
      call states_end(s%st)
      deallocate(s%st); nullify(s%st)
    end if

!!!!NEW
    call geometry_end(s%geo)
!!!!ENDOFNEW

    call grid_end(s%gr)
    deallocate(s%gr);  nullify(s%gr)

    call pop_sub()
  end subroutine system_end


  !----------------------------------------------------------
  subroutine system_h_setup(sys, h)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h

    call push_sub('systm.hamiltonian_setup')

    call states_fermi(sys%st, sys%gr%m)
    call states_calc_dens(sys%st, sys%gr%m%np, sys%st%rho)

    call v_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
    call states_fermi(sys%st, sys%gr%m)                            ! occupations
    call hamiltonian_energy(h, sys%gr, sys%geo, sys%st, -1)            ! total energy

    call pop_sub()
  end subroutine system_h_setup

end module system_m
