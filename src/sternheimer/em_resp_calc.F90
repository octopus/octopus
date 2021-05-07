!! Copyright (C) 2004-2012 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca), David Strubbe
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

module em_resp_calc_oct_m
  use batch_oct_m
  use comm_oct_m
  use density_oct_m
  use derivatives_oct_m
  use elf_oct_m
  use grid_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use linear_response_oct_m
  use magnetic_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use pert_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use sternheimer_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use utils_oct_m
  use xc_oct_m

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
     dinhomog_B,                       &
     zinhomog_B,                       &
     dinhomog_KB_tot,                  &
     zinhomog_KB_tot,                  &
     dinhomog_KE_tot,                  &
     zinhomog_KE_tot,                  &
     dinhomog_K2_tot,                  &
     zinhomog_K2_tot,                  &
     dlr_calc_magnetization_periodic,  &
     zlr_calc_magnetization_periodic,  &
     dlr_calc_magneto_optics_finite,   &
     zlr_calc_magneto_optics_finite,   &
     dlr_calc_magneto_optics_periodic, & 
     zlr_calc_magneto_optics_periodic, &
     dlr_calc_susceptibility,          &
     zlr_calc_susceptibility,          &
     dlr_calc_susceptibility_periodic, &
     zlr_calc_susceptibility_periodic, &
     dlr_calc_beta,                    &
     zlr_calc_beta,                    &
     freq2str,                         &
     magn_dir,                         &
     em_wfs_tag,                       &
     em_rho_tag,                       &
     dpost_orthogonalize,              &
     zpost_orthogonalize,              &
     dem_resp_calc_eigenvalues,        &
     zem_resp_calc_eigenvalues

  type(profile_t), save :: beta_prof

  type matrix_t
    private
    FLOAT, allocatable :: dmatrix(:, :)
    CMPLX, allocatable :: zmatrix(:, :)
  end type matrix_t

contains

  ! ---------------------------------------------------------
  subroutine lr_calc_current(st, gr, lr, lr_m)
    type(states_elec_t),  intent(inout) :: st
    type(grid_t),         intent(in)    :: gr
    type(lr_t),           intent(inout) :: lr
    type(lr_t), optional, intent(inout) :: lr_m

    integer :: idir, ist, ispin, idim, ndim, np

    CMPLX, allocatable :: psi(:, :), gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)

    PUSH_SUB(lr_calc_current)

    if (.not. allocated(lr%dl_j)) then
      SAFE_ALLOCATE(lr%dl_j(1:gr%mesh%np, 1:gr%sb%dim, 1:st%d%nspin))
    end if

    np = gr%mesh%np
    ndim = gr%sb%dim

    SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:ndim))
    SAFE_ALLOCATE(gpsi(1:np, 1:ndim))
    SAFE_ALLOCATE(gdl_psi(1:np, 1:ndim))
    if(present(lr_m)) then
      SAFE_ALLOCATE(gdl_psi_m(1:np, 1:ndim))
    end if

    lr%dl_j = M_ZERO

    do ispin = 1, st%d%nspin
      do ist = 1, st%nst

        call states_elec_set_state(st, gr%mesh, ist, ispin, psi)
        
        do idim = 1, st%d%dim

          call zderivatives_grad(gr%der, lr%zdl_psi(:, idim, ist, ispin), gdl_psi)
          call zderivatives_grad(gr%der, psi(:, idim), gpsi)

          if(present(lr_m)) then

            call zderivatives_grad(gr%der, lr_m%zdl_psi(:, idim, ist, ispin), gdl_psi_m)

            do idir = 1, gr%sb%dim 

              lr%dl_j(1:np, idir, ispin) = lr%dl_j(1:np, idir, ispin) + (           &
                + conjg(psi(1:np, idim))*gdl_psi(1:np, idir)   &
                -       psi(1:np, idim)*conjg(gdl_psi_m(1:np, idir))  &
                + conjg(lr_m%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np, idir)   & 
                -       lr%zdl_psi  (1:np, idim, ist, ispin)  * conjg(gpsi(1:np, idir))  &
                ) / (M_TWO * M_zI)
            end do

          else 

            do idir = 1, gr%sb%dim 

              lr%dl_j(1:np, idir, ispin) = lr%dl_j(1:np, idir, ispin) + (           &
                + conjg(psi(1:np, idim))*gdl_psi(1:np, idir)   &
                -       psi(1:np, idim)*conjg(gdl_psi(1:np, idir))  &
                + conjg(lr%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np, idir)   & 
                -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np, idir))  &
                ) / (M_TWO * M_zI)

            end do

          end if

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
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
      if(freq >= M_ZERO .and. str(1:1) /= '0') str = "0"//trim(str)
      if(freq <  M_ZERO .and. str(2:2) /= '0') str = "-0"//trim(str(2:))
    end if
    str = trim(adjustl(str))

    POP_SUB(freq2str)

  end function freq2str


! ---------------------------------------------------------
  character(len=100) function em_rho_tag(freq, dir, dir2, ipert) result(str)
    FLOAT,             intent(in) :: freq
    integer,           intent(in) :: dir
    integer, optional, intent(in) :: dir2
    integer, optional, intent(in) :: ipert

    character(len=12) :: str_tmp

    !this function has to be consistent with oct_search_file_lr in liboct/oct_f.c

    PUSH_SUB(em_rho_tag)

    str_tmp = freq2str(freq)
    write(str, '(3a,i1)') 'rho_', trim(str_tmp), '_', dir
    if(present(dir2)) write(str, '(2a,i1)') trim(str), "_", dir2
    if(present(ipert)) write(str, '(3a)') trim(str), "_", index2pert(ipert) 

    POP_SUB(em_rho_tag)

  end function em_rho_tag
  

! ---------------------------------------------------------
  character(len=100) function em_wfs_tag(idir, ifactor, idir2, ipert) result(str)
    integer,           intent(in) :: idir
    integer,           intent(in) :: ifactor
    integer, optional, intent(in) :: idir2
    integer, optional, intent(in) :: ipert 

    PUSH_SUB(em_wfs_tag)

    write(str, '(3a,i1)') "wfs_", index2axis(idir), "_f", ifactor
    if(present(idir2)) write(str, '(3a)') trim(str), "_", index2axis(idir2)
    if(present(ipert)) write(str, '(3a)') trim(str), "_", index2pert(ipert) 

    POP_SUB(em_wfs_tag)

  end function em_wfs_tag

! ---------------------------------------------------------
! Provides indices of axes forming the right-hand system 
! to the given axis (to choose components of the position 
! and velocity, r_\alpha and V_\beta, for calculation of 
! the M_\gamma component of the magnetic dipole moment 
! M_\gamma = e_{\alpha \beta \gamma} r_\alpha V_\beta /2)
  integer pure function magn_dir(dir, ind) result(dir_out)
    integer, intent(in) :: dir, ind
 
    select case(dir)
      case(1)
        if(ind == 1) then
          dir_out = 2
        else
          dir_out = 3
        end if
      case(2)
        if(ind == 1) then
          dir_out = 3
        else
          dir_out = 1
        end if
      case(3)
        if(ind == 1) then
          dir_out = 1
        else
          dir_out = 2
        end if
      case default
        dir_out = 0
    end select

  end function magn_dir


! ---------------------------------------------------------

  character(len=2) pure function index2pert(ipert) result(ch)
    integer, intent(in) :: ipert
    
    select case(ipert)
      case(1)
        ch = 'B'
      case(2)
        ch = 'K2'
      case(3)
        ch = 'KB'
      case(4)
        ch = 'KE'
      case(5)
        ch = 'E'
    end select

  end function index2pert
  
#include "undef.F90"
#include "real.F90"
#include "em_resp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "em_resp_calc_inc.F90"

end module em_resp_calc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
