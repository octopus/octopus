!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module elf_oct_m
  use density_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use mesh_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use states_oct_m
  use states_dim_oct_m

  implicit none

  private
  public :: elf_init,               &
            elf_calc

  logical :: with_current_term = .true.

contains

  subroutine elf_init(parser)
    type(parser_t),    intent(in)    :: parser
    
    PUSH_SUB(elf_init)

    !%Variable ELFWithCurrentTerm
    !%Type logical
    !%Default true
    !%Section Output
    !%Description
    !% The ELF, when calculated for complex wavefunctions, should contain
    !% a term dependent on the current. This term is properly calculated by
    !% default; however, for research purposes it may be useful not to add it.
    !% If this feature proves to be useless, this option should go away.
    !%End
    call parse_variable(parser, 'ELFWithCurrentTerm', .true., with_current_term)

    POP_SUB(elf_init)
  end subroutine elf_init

  ! ---------------------------------------------------------
  !> (time-dependent) electron localization function, (TD)ELF.
  ! ---------------------------------------------------------
  subroutine elf_calc(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in)    :: gr
    !> elf(gr%mesh%np, 1) if st%d%ispin = 1, elf(gr%mesh%np, 3) otherwise.
    !! On output, it should contain the global ELF if st%d%ispin = 1,
    !! otherwise elf(:, 3) contains the global ELF, and 
    !! elf(:, 1) and elf(:, 2) the spin-resolved ELF.
    FLOAT,            intent(inout) :: elf(:,:) 
    FLOAT,  optional, intent(inout):: de(:,:)

    FLOAT :: factor, D0, dens
    integer :: ip, is, nelfs
    FLOAT, allocatable :: rho(:,:), grho(:,:,:), jj(:,:,:)
    FLOAT, allocatable :: kappa(:,:)
    FLOAT, parameter :: dmin = CNST(1e-10)

    PUSH_SUB(elf_calc)

    ASSERT(gr%sb%dim == 2 .or. gr%sb%dim == 3)

    ! We may or may not want the total ELF. 
    ! If we want it, the argument elf should have three components.
    nelfs = size(elf, 2)

    SAFE_ALLOCATE(  rho(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(kappa(1:gr%mesh%np, 1:st%d%nspin))
    rho = M_ZERO
    kappa = M_ZERO
    call density_calc(st, gr, rho)

    SAFE_ALLOCATE(grho(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))
    SAFE_ALLOCATE(  jj(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:st%d%nspin))

    call states_calc_quantities(gr%der, st, .false., kinetic_energy_density = kappa, &
                                paramagnetic_current = jj, density_gradient = grho)

    ! spin-dependent quantities
    if(st%d%ispin == UNPOLARIZED) then
      rho = rho/M_TWO
      kappa = kappa/M_TWO
      jj = jj/M_TWO
      grho = grho/M_TWO
    end if
    
    if(.not.with_current_term) jj = M_ZERO

    ! kappa will contain rho * D
    do_is: do is = 1, st%d%nspin
      do ip = 1, gr%mesh%np
        kappa(ip, is) = kappa(ip, is)*rho(ip, is)        &    ! + tau * rho
          - M_FOURTH*sum(grho(ip, 1:gr%mesh%sb%dim, is)**2)      &    ! - | nabla rho |^2 / 4
          - sum(jj(ip, 1:gr%mesh%sb%dim, is)**2)                      ! - j^2
      end do

      ! pass this information to the caller if requested
      if(present(de)) de(1:gr%mesh%np,is) = kappa(1:gr%mesh%np,is)

    end do do_is

    SAFE_DEALLOCATE_A(grho)
    SAFE_DEALLOCATE_A(jj)   ! these are no longer needed

    select case(gr%sb%dim)
      case(3); factor = M_THREE / M_FIVE * (CNST(6.0) * M_PI**2)**M_TWOTHIRD
      case(2); factor = M_TWO * M_Pi
    end select

    select case(st%d%ispin)
    case(UNPOLARIZED)
      do ip = 1, gr%mesh%np
        if(rho(ip, 1) >= dmin) then
          select case(gr%sb%dim)
            case(3); D0 = factor * rho(ip, 1)**(CNST(8.0) / M_THREE)
            case(2); D0 = factor * rho(ip, 1)**3
          end select
          elf(ip, 1) = D0 * D0 / (D0 * D0 + kappa(ip, 1)**2)
        else
          elf(ip, 1) = M_ZERO
        end if
      end do

    case(SPIN_POLARIZED, SPINORS)
      if(nelfs  ==  3) then
        do ip = 1, gr%mesh%np
          dens = rho(ip, 1) + rho(ip, 2)
          if( dens >= dmin ) then
            select case(gr%sb%dim)
              case(3); D0 = factor * dens ** (M_FIVE/M_THREE) * rho(ip, 1) * rho(ip, 2)
              case(2); D0 = factor * dens ** 2 * rho(ip, 1) * rho(ip, 2)
            end select
            elf(ip, 3) = D0 * D0 / (D0 * D0 + (kappa(ip, 1) * rho(ip, 2) + kappa(ip,2) * rho(ip, 1))**2)
          else
            elf(ip, 3) = M_ZERO
          end if
        end do
      end if
      do ip = 1, gr%mesh%np
        do is = 1, st%d%spin_channels
          if(rho(ip, is) >= dmin) then
            select case(gr%sb%dim)
              case(3); D0 = factor * rho(ip, is)**(CNST(8.0)/M_THREE)
              case(2); D0 = factor * rho(ip, is)**3
            end select
            elf(ip, is) = D0 * D0 / (D0 * D0 + kappa(ip,is)**2)
          else
            elf(ip, is) = M_ZERO
          end if
        end do
      end do
    end select

    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(kappa)

    POP_SUB(elf_calc)
  end subroutine elf_calc

end module elf_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
