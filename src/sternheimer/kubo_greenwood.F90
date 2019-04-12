!! Copyright (C) 2018 X. Andrade
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

module kubo_greenwood_oct_m
  use boundaries_oct_m
  use derivatives_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use restart_oct_m
  use smear_oct_m
  use states_oct_m
  use states_restart_oct_m
  use system_oct_m
  use unit_system_oct_m
  use parser_oct_m
  use profiling_oct_m
 
  implicit none

  private
  public ::                       &
       kubo_greenwood_run

contains
  
  subroutine kubo_greenwood_run(sys, hm)
    type(system_t),       intent(inout) :: sys
    type(hamiltonian_t),  intent(inout) :: hm

    type(restart_t) :: gs_restart
    integer :: ierr, nfreq, ifreq
    integer :: ist, jst, iqn, idim, idir, jdir
    CMPLX, allocatable :: psii(:, :), psij(:, :), gpsii(:, :, :), gpsij(:, :, :)
    CMPLX, allocatable :: tensor(:, :, :)
    CMPLX :: prod
    FLOAT :: eigi, eigj, occi, occj, df, width, dfreq, maxfreq
    type(mesh_t), pointer :: mesh
    
    PUSH_SUB(kubo_greewood_run)

    call messages_write('Info: Starting Kubo-Greenwood linear-response calculation.')
    call messages_info()
    
    call messages_experimental('Kubo Greenwood')

      !%Variable KuboGreenwoodwidth
      !%Type float
      !%Default 0.0
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% The width applied to the Kubo Greenwood conductivity.
      !% In units of energy. Cannot be negative.
      !%End

      call parse_variable('KuboGreenwoodWidth', M_ZERO, width, units_inp%energy)
      if(width < -M_EPSILON) then
        message(1) = "KuboGreenwoodWidth cannot be negative."
        call messages_fatal(1)
      end if

      !%Variable KuboGreenwoodMaxEnergy
      !%Type float
      !%Default 2.0
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% Maximum energy for Kubo Greenwood
      !% In units of energy. 
      !%End

      call parse_variable('KuboGreenwoodMaxEnergy', CNST(2.0),maxfreq , units_inp%energy)

      !%Variable KuboGreenwoodEnergySpacing
      !%Type float
      !%Default 0.01
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% Spacing of frequencies for Kubo Greenwood
      !% In units of energy. 
      !%End

      call parse_variable('KuboGreenwoodEnergySpacing', CNST(0.01), dfreq, units_inp%energy)

      
      call restart_init(gs_restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh = sys%gr%mesh, exact = .true.)
    if(ierr == 0) then
      call states_look_and_load(gs_restart, sys%st, sys%gr)
      call restart_end(gs_restart)
    else
      call messages_write("Cannot find occupied and unoccupied states, a previous ground state calculation is required.")
      call messages_fatal()
    end if

    ASSERT(.not. sys%st%parallel_in_states)    

    mesh => sys%gr%mesh
    
    SAFE_ALLOCATE(psii(1:mesh%np_part, 1:sys%st%d%dim))
    SAFE_ALLOCATE(psij(1:mesh%np_part, 1:sys%st%d%dim))
    SAFE_ALLOCATE(gpsii(1:mesh%np, 1:sys%gr%sb%dim, 1:sys%st%d%dim))
    SAFE_ALLOCATE(gpsij(1:mesh%np, 1:sys%gr%sb%dim, 1:sys%st%d%dim))
    nfreq = nint(maxfreq / (dfreq+1))
    
    SAFE_ALLOCATE(tensor(1:mesh%sb%dim, 1:mesh%sb%dim,1:nfreq))

    tensor = CNST(0.0)

    do iqn = sys%st%d%kpt%start, sys%st%d%kpt%end
    
       do ist = 1, sys%st%nst
        ! get the state i and calculate the gradient
        call states_get_state(sys%st, mesh, ist, iqn, psii)
        do idim = 1, sys%st%d%dim
          call boundaries_set(sys%gr%der%boundaries, psii(:, idim))
        end do
        if(associated(hm%hm_base%phase)) then 
          call states_set_phase(sys%st%d, psii, hm%hm_base%phase(:, iqn), mesh%np_part, .false.)
        end if
        do idim = 1, sys%st%d%dim
          call zderivatives_grad(sys%gr%der, psii(:, idim), gpsii(:, :, idim), set_bc = .false.)
        end do
        
        do jst = 1, sys%st%nst
          ! get the state j and calculate the gradient
          call states_get_state(sys%st, mesh, ist, iqn, psij)
          do idim = 1, sys%st%d%dim
            call boundaries_set(sys%gr%der%boundaries, psij(:, idim))
          end do
          if(associated(hm%hm_base%phase)) then 
            call states_set_phase(sys%st%d, psij, hm%hm_base%phase(:, iqn), mesh%np_part, .false.)
          end if
          do idim = 1, sys%st%d%dim
            call zderivatives_grad(sys%gr%der, psij(:, idim), gpsij(:, :, idim), set_bc = .false.)
          end do
          eigi = sys%st%eigenval(ist,iqn)
          eigj = sys%st%eigenval(jst,iqn)
          occi = sys%st%occ(ist,iqn)
          occj = sys%st%occ(jst,iqn)
          if (abs(eigi-eigj) < CNST(1.0e-10)) then
             df = smear_step_function_der(sys%st%smear,eigi)
          else
             df = (occj - occi)/(eigj - eigi)
          endif

          do idir = 1, mesh%sb%dim
             do jdir = 1, mesh%sb%dim
              prod = zmf_dotp(mesh, sys%st%d%dim, psii, gpsij(:, idir, :))*zmf_dotp(mesh, sys%st%d%dim, psij, gpsii(:, jdir, :))
              
              do ifreq = 1, nfreq
                 tensor(idir, jdir,ifreq) = tensor(idir, jdir,ifreq) - sys%st%d%kweights(iqn)*(CNST(2.0)/mesh%sb%rcell_volume)* &
                      real(prod,REAL_PRECISION) * (CNST(0.5)*width + M_ZI*(eigi-eigj - ifreq*dfreq))/ ((eigi-eigj - ifreq*dfreq)**2 + width/CNST(4.0))
              end do
            end do
          end do
          
        end do
      end do
      
    end do

    SAFE_DEALLOCATE_A(tensor)
    SAFE_DEALLOCATE_A(psii)
    SAFE_DEALLOCATE_A(psij)
    SAFE_DEALLOCATE_A(gpsii)
    SAFE_DEALLOCATE_A(gpsij)    
    
    POP_SUB(kubo_greewood_run)
  end subroutine kubo_greenwood_run

end module kubo_greenwood_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
