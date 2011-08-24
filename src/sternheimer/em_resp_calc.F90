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
  use density_m
  use derivatives_m
  use elf_m
  use geometry_m
  use grid_m
  use global_m
  use hamiltonian_m
  use linear_response_m
  use magnetic_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use pert_m
  use poisson_m
  use profiling_m
  use states_m
  use states_block_m
  use states_dim_m
  use sternheimer_m
  use system_m
  use utils_m
  use xc_m

  implicit none

  private
  public ::                            &
     lr_calc_current,                  &
     dlr_calc_elf,                     &
     zlr_calc_elf,                     &
     dcalc_polarizability_finite,      &
     zcalc_polarizability_finite,      &
     dcalc_polarizability_periodic,    &
     zcalc_polarizability_periodic,    &
     dlr_calc_susceptibility,          &
     zlr_calc_susceptibility,          &
     dlr_calc_beta,                    &
     zlr_calc_beta,                    &
     freq2str,                         &
     em_wfs_tag,                       &
     em_rho_tag,                       &
     dpost_orthogonalize,              &
     zpost_orthogonalize

  type(profile_t), save :: beta_prof

  type matrix_t
    FLOAT, pointer :: dmatrix(:, :)
    CMPLX, pointer :: zmatrix(:, :)
  end type matrix_t

contains

  ! ---------------------------------------------------------
  subroutine lr_calc_current(st, gr, lr, lr_m)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(lr_t),           intent(inout) :: lr
    type(lr_t), optional, intent(inout) :: lr_m

    integer :: idir, ist, ispin, idim, ndim, np

    CMPLX, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)

    PUSH_SUB(lr_calc_current)

    if(.not. associated(lr%dl_j)) then
      SAFE_ALLOCATE(lr%dl_j(1:gr%mesh%np, 1:MAX_DIM, 1:st%d%nspin))
    end if

    np = NP
    ndim = gr%mesh%sb%dim

    SAFE_ALLOCATE(   gpsi(1:np, 1:ndim))
    SAFE_ALLOCATE(gdl_psi(1:np, 1:ndim))
    if(present(lr_m)) then
      SAFE_ALLOCATE(gdl_psi_m(1:np, 1:ndim))
    end if

    lr%dl_j = M_ZERO

    do ispin = 1, st%d%nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim

          call zderivatives_grad(gr%der, lr%zdl_psi(:, idim, ist, ispin), gdl_psi)
          call zderivatives_grad(gr%der, st%zpsi(:, idim, ist, ispin), gpsi)

          if(present(lr_m)) then               

            call zderivatives_grad(gr%der, lr_m%zdl_psi(:, idim, ist, ispin), gdl_psi_m)

            do idir = 1, gr%mesh%sb%dim 

              lr%dl_j(1:np, idir, ispin) = lr%dl_j(1:np, idir, ispin) + (           &
                + conjg(st%zpsi(1:np, idim, ist, ispin)) *       gdl_psi  (1:np, idir)   &
                -       st%zpsi(1:np, idim, ist, ispin)  * conjg(gdl_psi_m(1:np, idir))  &
                + conjg(lr_m%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np, idir)   & 
                -       lr%zdl_psi  (1:np, idim, ist, ispin)  * conjg(gpsi(1:np, idir))  &
                ) / (M_TWO * M_zI)
            end do

          else 

            do idir = 1, gr%mesh%sb%dim 

              lr%dl_j(1:np, idir, ispin) = lr%dl_j(1:np, idir, ispin) + (           &
                + conjg(st%zpsi(1:np, idim, ist, ispin)) *       gdl_psi(1:np, idir)   &
                -       st%zpsi(1:np, idim, ist, ispin)  * conjg(gdl_psi(1:np, idir))  &
                + conjg(lr%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np, idir)   & 
                -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np, idir))  &
                ) / (M_TWO * M_zI)

            end do

          end if

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(gpsi)
    SAFE_DEALLOCATE_A(gdl_psi)
    if(present(lr_m)) then
      SAFE_DEALLOCATE_A(gdl_psi_m)
    end if

    POP_SUB(lr_calc_current)

  end subroutine lr_calc_current


! ---------------------------------------------------------
  character(len=12) function freq2str(freq) result(str)
    FLOAT, intent(in) :: freq

    PUSH_SUB(freq2str)

    ! some compilers (xlf) do not put a leading zero when the number
    ! is smaller than 1. We have to check and correct that behavior.

    write(str, '(f11.4)') freq
    str = adjustl(str)
    if(abs(freq) < M_ONE) then
      if(freq >= M_ZERO .and. str(1:1).ne.'0') str = "0"//trim(str)
      if(freq <  M_ZERO .and. str(2:2).ne.'0') str = "-0"//trim(str(2:))
    end if
    str = trim(adjustl(str))

    POP_SUB(freq2str)

  end function freq2str


! ---------------------------------------------------------
  character(len=100) function em_rho_tag(freq, dir) result(str)
    FLOAT,   intent(in) :: freq
    integer, intent(in) :: dir

    character(len=12) :: str_tmp

    !this function has to be consistent with oct_search_file_lr in liboct/oct_f.c

    PUSH_SUB(em_rho_tag)

    str_tmp = freq2str(freq)
    write(str, '(3a,i1)') 'rho_', trim(str_tmp), '_', dir

    POP_SUB(em_rho_tag)

  end function em_rho_tag
  

! ---------------------------------------------------------
  character(len=100) function em_wfs_tag(idir, ifactor, idir2) result(str)
    integer,           intent(in) :: idir
    integer,           intent(in) :: ifactor
    integer, optional, intent(in) :: idir2

    PUSH_SUB(em_wfs_tag)

    write(str, '(3a,i1)') "wfs_", index2axis(idir), "_f", ifactor
    if(present(idir2)) write(str, '(3a)') trim(str), "_", index2axis(idir2)

    POP_SUB(em_wfs_tag)

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
