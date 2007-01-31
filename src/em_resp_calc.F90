!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: em_resp.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module em_resp_calc_m
  use datasets_m
  use elf_m
  use functions_m
  use grid_m
  use global_m
  use lib_oct_parser_m
  use linear_response_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use poisson_m
  use states_m
  use system_m
  use xc_m

  implicit none

  private
  public ::                  &
    pol_props_t,             &
    pol_props_init,          &
    lr_calc_current,         &
    dlr_calc_elf,            &
    zlr_calc_elf,            &
    dlr_calc_polarizability, &
    zlr_calc_polarizability, &
    dlr_calc_beta,           &
    zlr_calc_beta

  type pol_props_t
    logical :: add_fxc
    logical :: add_hartree
    logical :: from_scratch
    logical :: orth_response
  end type pol_props_t

contains

  subroutine pol_props_init(props, ip_app)
    type(pol_props_t),  intent(out) :: props
    logical,            intent(in)  :: ip_app

    integer :: ham_var

    !%Variable PolOrthResponse
    !%Type logical
    !%Default false
    !%Section Linear Response::Polarizabilities
    !%Description
    !% Wheter variations should be orthogonalized or not against the
    !% occupied states.
    !%End

    call loct_parse_logical(check_inp('PolOrthResponse'), .true., props%orth_response)


    !%Variable PolHamiltonianVariation
    !%Type integer
    !%Default hartree+fxc
    !%Section Linear Response::Polarizabilities
    !%Description
    !% The terms are considered in the variation of the
    !% hamiltonian. V_ext is always considered. The default is to include
    !% the fxc and hartree terms. If you want to do RPA only include
    !% hartree.
    !%Option hartree 1 
    !% The variation of the hartree potential.
    !%Option fxc 2
    !% The exchange and correlation kernel, the variation of the
    !% exchange and correlation potential.
    !%End

    if(.not. ip_app) then 
      call loct_parse_int(check_inp('PolHamiltonianVariation'), 3, ham_var)    
      props%add_fxc = ((ham_var/2) == 1)
      props%add_hartree = (mod(ham_var, 2) == 1)
    else
      props%add_fxc = .false. 
      props%add_hartree = .false.
    end if
    
    message(1) = 'Hamiltonian variation: V_ext'
    if(props%add_hartree) write(message(1), '(2a)') trim(message(1)), ' + hartree'
    if(props%add_fxc)     write(message(1), '(2a)') trim(message(1)), ' + fxc'
    call write_info(1)

  end subroutine pol_props_init


  ! ---------------------------------------------------------
  subroutine lr_calc_current(st, gr, lr, lr_m)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    type(lr_t),       intent(inout) :: lr
    type(lr_t), optional, intent(inout) :: lr_m

    integer :: k, ist, ispin, idim, ndim, np

    CMPLX, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)

    call push_sub('em_resp.lr_calc_current')

    if(.not. associated(lr%dl_j)) ALLOCATE(lr%dl_j(gr%m%np, MAX_DIM, st%d%nspin), gr%m%np*MAX_DIM*st%d%nspin)

    np = NP
    ndim = NDIM

    ALLOCATE(   gpsi(1:np, 1:ndim), np*ndim)
    ALLOCATE(gdl_psi(1:np, 1:ndim), np*ndim)
    if(present(lr_m)) ALLOCATE(gdl_psi_m(1:np, 1:ndim), np*ndim)

    lr%dl_j = M_ZERO

    do ispin = 1, st%d%nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim

          call zf_gradient(gr%sb, gr%f_der, lr%zdl_psi(:, idim, ist, ispin), gdl_psi)
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, idim, ist, ispin), gpsi)

          if(present(lr_m)) then               

            call zf_gradient(gr%sb, gr%f_der, lr_m%zdl_psi(:, idim, ist, ispin), gdl_psi_m)

            do k = 1, NDIM 

              lr%dl_j(1:np,k,ispin) = lr%dl_j(1:np, k, ispin) + (           &
                + conjg(st%zpsi(1:np, idim, ist, ispin)) *      gdl_psi(1:np,k)   &
                -       st%zpsi(1:np, idim, ist, ispin) * conjg(gdl_psi_m(1:np,k))  &
                + conjg(lr_m%zdl_psi(1:np, idim, ist, ispin)) *     gpsi(1:np,k)   & 
                -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np,k))  &
                )/(M_TWO*M_zI)
            end do

          else 

            do k = 1, NDIM 

              lr%dl_j(1:np,k,ispin) = lr%dl_j(1:np, k, ispin) + (           &
                + conjg(st%zpsi(1:np, idim, ist, ispin)) *       gdl_psi(1:np,k)   &
                -       st%zpsi(1:np, idim, ist, ispin)  * conjg(gdl_psi(1:np,k))  &
                + conjg(lr%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np,k)   & 
                -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np,k))  &
                )/(M_TWO*M_zI)

            end do

          end if

        end do
      end do
    end do

    deallocate(gpsi)
    deallocate(gdl_psi)
    if(present(lr_m)) deallocate(gdl_psi_m)

    call pop_sub()

  end subroutine lr_calc_current

#include "undef.F90"
#include "real.F90"
#include "em_resp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "em_resp_calc_inc.F90"

end module em_resp_calc_m
