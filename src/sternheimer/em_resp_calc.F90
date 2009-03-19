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
!! $Id: em_resp_calc.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module em_resp_calc_m
  use datasets_m
  use derivatives_m
  use elf_m
  use grid_m
  use global_m
  use hamiltonian_m
  use loct_parser_m
  use linear_response_m
  use magnetic_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use poisson_m
  use pert_m
  use profiling_m
  use states_m
  use states_dim_m
  use sternheimer_m
  use system_m
  use xc_m

  implicit none

  private
  public ::                            &
     lr_calc_current,                  &
     dlr_calc_elf,                     &
     zlr_calc_elf,                     &
     dcalc_polarizability_finite,      &
     zcalc_polarizability_finite,      &
     zcalc_polarizability_periodic,    &
     dlr_calc_susceptibility,          &
     zlr_calc_susceptibility,          &
     dlr_calc_beta,                    &
     zlr_calc_beta,                    &
     freq2str,                         &
     em_wfs_tag,                       &
     em_rho_tag
  ! periodic version of polarizability in kdotp

  type(profile_t), save :: beta_prof

contains

  ! ---------------------------------------------------------
  subroutine lr_calc_current(st, gr, lr, lr_m)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    type(lr_t),       intent(inout) :: lr
    type(lr_t), optional, intent(inout) :: lr_m

    integer :: k, ist, ispin, idim, ndim, np

    CMPLX, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)

    call push_sub('em_resp_calc.lr_calc_current')

    if(.not. associated(lr%dl_j)) ALLOCATE(lr%dl_j(gr%mesh%np, MAX_DIM, st%d%nspin), gr%mesh%np*MAX_DIM*st%d%nspin)

    np = NP
    ndim = gr%mesh%sb%dim

    ALLOCATE(   gpsi(1:np, 1:ndim), np*ndim)
    ALLOCATE(gdl_psi(1:np, 1:ndim), np*ndim)
    if(present(lr_m)) ALLOCATE(gdl_psi_m(1:np, 1:ndim), np*ndim)

    lr%dl_j = M_ZERO

    do ispin = 1, st%d%nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim

          call zderivatives_grad(gr%der, lr%zdl_psi(:, idim, ist, ispin), gdl_psi)
          call zderivatives_grad(gr%der, st%zpsi(:, idim, ist, ispin), gpsi)

          if(present(lr_m)) then               

            call zderivatives_grad(gr%der, lr_m%zdl_psi(:, idim, ist, ispin), gdl_psi_m)

            do k = 1, gr%mesh%sb%dim 

              lr%dl_j(1:np,k,ispin) = lr%dl_j(1:np, k, ispin) + (           &
                + conjg(st%zpsi(1:np, idim, ist, ispin)) *      gdl_psi(1:np,k)   &
                -       st%zpsi(1:np, idim, ist, ispin) * conjg(gdl_psi_m(1:np,k))  &
                + conjg(lr_m%zdl_psi(1:np, idim, ist, ispin)) *     gpsi(1:np,k)   & 
                -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np,k))  &
                )/(M_TWO*M_zI)
            end do

          else 

            do k = 1, gr%mesh%sb%dim 

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


! ---------------------------------------------------------
subroutine zcalc_polarizability_periodic(sys, em_lr, kdotp_lr, nsigma, zpol, ndir)
  type(system_t),         intent(inout) :: sys
  type(lr_t),             intent(inout) :: em_lr(:,:)
  type(lr_t),             intent(inout) :: kdotp_lr(:)
  integer,                intent(in)    :: nsigma
  CMPLX,                  intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)
  integer, optional,      intent(in)    :: ndir

  integer :: dir1, dir2, ndir_, ist, ik, idim
  CMPLX :: term, subterm
  type(mesh_t), pointer :: m
  m => sys%gr%mesh

  call push_sub('em_resp_calc.zcalc_polarizability_periodic')

  ndir_ = sys%gr%mesh%sb%dim
  if(present(ndir)) ndir_ = ndir

  ! alpha_ij(w) = -e sum(m occ, k) [(<u_mk(0)|-id/dk_i)|u_mkj(1)(w)> + <u_mkj(1)(-w)|(-id/dk_i|u_mk(0)>)]
  ! Smearing is not implemented here yet?

  do dir1 = 1, ndir_
    do dir2 = 1, sys%gr%sb%dim

      zpol(dir1, dir2) = M_ZERO

      do ik = 1, sys%st%d%nik
        term = M_ZERO
        do ist = 1, sys%st%nst
          do idim = 1, sys%st%d%dim
            subterm = M_zI * zmf_dotp(m, kdotp_lr(dir1)%zdl_psi(1:m%np, idim, ist, ik), &
              em_lr(dir2, 1)%zdl_psi(1:m%np, idim, ist, ik))
            term = term + subterm

            if(nsigma == 1) then
              term = term + conjg(subterm)
            else
              term = term - M_zI * zmf_dotp(m, em_lr(dir2, 2)%zdl_psi(1:m%np, idim, ist, ik), & 
                kdotp_lr(dir1)%zdl_psi(1:m%np, idim, ist, ik))
            end if
          enddo
        enddo

        zpol(dir1, dir2) = zpol(dir1, dir2) + &
          term * sys%st%d%kweights(ik) * sys%st%smear%el_per_state
      enddo
    enddo
  enddo

  call pop_sub()

end subroutine zcalc_polarizability_periodic


! ---------------------------------------------------------
  character(len=12) function freq2str(w) result(str)
    FLOAT, intent(in) :: w

    call push_sub('em_resp_calc.freq2str')

    write(str, '(f11.4)') w
    str = trim(adjustl(str))

    call pop_sub()

  end function freq2str

  character(len=100) function em_rho_tag(w, dir) result(str)
    FLOAT, intent(in) :: w
    integer, intent(in) :: dir
    character(len=12) :: str_tmp

    !this function has to be consistent with oct_search_file_lr in liboct/oct_f.c

    call push_sub('em_resp_calc.em_rho_tag')

    str_tmp = freq2str(w)
    write(str, '(3a,i1)') 'rho_', trim(str_tmp), '_', dir

    call pop_sub()

  end function em_rho_tag
  
  character(len=100) function em_wfs_tag(dir, ifactor) result(str)
    integer, intent(in) :: dir, ifactor 

    call push_sub('em_resp_calc.em_wfs_tag')

    write(str, '(a,i1,a,i1)') "wfs_", dir, "_", ifactor

    call pop_sub()

  end function em_wfs_tag
  
#include "undef.F90"
#include "real.F90"
#include "em_resp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "em_resp_calc_inc.F90"

end module em_resp_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
