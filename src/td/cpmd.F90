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
!! $Id: td.F90 3694 2008-02-15 13:37:54Z marques $

#include "global.h"
  
module cpmd_m
  use global_m
  use io_m
  use datasets_m
  use io_function_m
  use loct_math_m
  use loct_parser_m
  use units_m
  use messages_m
  use mesh_m
  use external_pot_m
  use geometry_m
  use ground_state_m
  use hamiltonian_m
  use loct_m
  use profiling_m
  use scf_m
  use states_m
  use system_m
  use v_ks_m
  use grid_m

  implicit none

  private
  public ::               &
    cpmd_t,               &
    cpmd_init,            &
    cpmd_end,             &
    cpmd_propagate

  type cpmd_t
    private
    FLOAT          :: emass
    CMPLX, pointer :: psidot(:, :, :, :)
  end type cpmd_t

contains

  subroutine cpmd_init(this, gr, st)
    type(cpmd_t),        intent(out)   :: this
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    
    integer :: size

    !%Variable CPElectronicMass
    !%Type float
    !%Default 1.0
    !%Section Time Dependent::Propagation
    !%Description
    !% The fictious electronic mass used to propagate the electronic
    !% wavefunctions in the Carr-Parrinelo formalism.
    !%End
    
    call loct_parse_float(check_inp('CPElectronicMass'), CNST(1.0), this%emass)

    size = gr%m%np_part * st%d%dim * st%lnst * st%d%nik
    ALLOCATE(this%psidot(gr%m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), size)

  end subroutine cpmd_init
  
  subroutine cpmd_end(this)
    type(cpmd_t), intent(inout) :: this

    deallocate(this%psidot)

  end subroutine cpmd_end
  
  subroutine cpmd_propagate(this, gr, h, st, iter, dt)
    type(cpmd_t),        intent(inout) :: this
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: h
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer :: ik, ist, idim, ip

    CMPLX, allocatable :: hpsi(:, :)

    ALLOCATE(hpsi(1:gr%m%np, 1:st%d%dim), gr%m%np*st%d%dim)

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end

        call zHpsi(h, gr, st%zpsi(:, :, ist, ik), hpsi, ist, ik)
        
        do idim = 1, h%d%dim
          do ip = 1, gr%m%np
            ! the integration of the velocity is made in two steps

            if(iter > 1) then
              ! the remaining term from previous steps
              ! v(t) = v_unc(t) + 1/2*dt*a(t) 
              this%psidot(ip, idim, ist, ik) = this%psidot(ip, idim, ist, ik) - M_HALF*dt*hpsi(ip, idim)/this%emass
            end if

            ! x(t + dt) = x(t) + dt*v(t) + 1/2*dt^2*a(t)
            st%zpsi(ip, idim, ist, ik) = &
              st%zpsi(ip, idim, ist, ik) + dt*this%psidot(ip, idim, ist, ik) - M_HALF*dt*dt*hpsi(ip, idim)/this%emass

            ! the missing term is added in next iteration
            ! vunc(t + dt) = v(t) + 1/2*dt*a(t)
             this%psidot(ip, idim, ist, ik) = this%psidot(ip, idim, ist, ik) - M_HALF*dt*hpsi(ip, idim)/this%emass
          end do
        end do

      end do
      
      call zstates_gram_schmidt(st, st%st_end, gr%m, st%d%dim, st%zpsi(:, :, :, ik), start = st%st_start)

    end do


    deallocate(hpsi)

  end subroutine cpmd_propagate

end module cpmd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
