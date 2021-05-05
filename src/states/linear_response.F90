!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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

module linear_response_oct_m
  use comm_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use profiling_oct_m
  use restart_oct_m
  use smear_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m

  implicit none

  private

  public ::               &
       lr_t,              &
       lr_init,           &
       lr_allocate,       &
       lr_zero,           &
       lr_copy,           &
       dlr_orth_vector,   & 
       zlr_orth_vector,   & 
       dlr_build_dl_rho,  & 
       zlr_build_dl_rho,  &
       dlr_orth_response, &
       zlr_orth_response, &
       lr_alpha_j,        &
       lr_dealloc,        &
       lr_is_allocated,   &
       dlr_swap_sigma,    &
       zlr_swap_sigma,    &
       dlr_dump_rho,      &
       zlr_dump_rho,      &
       dlr_load_rho,      &
       zlr_load_rho



  type lr_t
    ! Components are public by default
    logical, private :: is_allocated, is_allocated_rho
     
    !> the real quantities
    FLOAT, allocatable :: ddl_rho(:,:)     !< response of the density
    FLOAT, allocatable :: ddl_psi(:,:,:,:) !< linear change of the real KS orbitals
    
    !> and the complex version
    CMPLX, allocatable :: zdl_rho(:,:)     !< response of the density
    CMPLX, allocatable :: zdl_psi(:,:,:,:) !< linear change of the complex KS orbitals

    !> other observables
    CMPLX, allocatable :: dl_j(:,:,:)     !< response of the current
    FLOAT, allocatable :: ddl_de(:,:)     !< unnormalized ELF
    FLOAT, allocatable :: ddl_elf(:,:)    !< normalized ELF
    CMPLX, allocatable :: zdl_de(:,:)     !< unnormalized ELF
    CMPLX, allocatable :: zdl_elf(:,:)    !< normalized ELF
    
  end type lr_t

contains

  ! ---------------------------------------------------------
  subroutine lr_init(lr)
    type(lr_t), intent(out) :: lr

    PUSH_SUB(lr_init)

    lr%is_allocated = .false.

    POP_SUB(lr_init)

  end subroutine lr_init



  ! ---------------------------------------------------------
  subroutine lr_allocate(lr, st, mesh, allocate_rho)
    type(lr_t),          intent(inout) :: lr
    type(states_elec_t), intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    logical, optional,      intent(in) :: allocate_rho

    PUSH_SUB(lr_allocate)

    lr%is_allocated_rho = .true.
    if(present(allocate_rho)) lr%is_allocated_rho = allocate_rho

    if (states_are_complex(st)) then
      SAFE_ALLOCATE(lr%zdl_psi(1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      if(lr%is_allocated_rho) SAFE_ALLOCATE(lr%zdl_rho(1:mesh%np, 1:st%d%nspin))
    else
      SAFE_ALLOCATE(lr%ddl_psi(1:mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
      if(lr%is_allocated_rho) SAFE_ALLOCATE(lr%ddl_rho(1:mesh%np, 1:st%d%nspin))
    end if

    lr%is_allocated = .true.
    
    call lr_zero(lr, st)

    POP_SUB(lr_allocate)

  end subroutine lr_allocate



  ! ---------------------------------------------------------
  subroutine lr_zero(lr, st)
    type(lr_t),          intent(inout) :: lr
    type(states_elec_t), intent(in)    :: st

    PUSH_SUB(lr_zero)

    ASSERT(lr%is_allocated)

    if (states_are_complex(st)) then
      lr%zdl_psi = M_ZERO
      if(lr%is_allocated_rho) lr%zdl_rho = M_ZERO
    else
      lr%ddl_psi = M_ZERO
      if(lr%is_allocated_rho) lr%ddl_rho = M_ZERO
    end if

    POP_SUB(lr_zero)

  end subroutine lr_zero


  ! ---------------------------------------------------------
  subroutine lr_dealloc(lr)
    type(lr_t), intent(inout) :: lr

    PUSH_SUB(lr_dealloc)

    SAFE_DEALLOCATE_A(lr%ddl_psi)
    SAFE_DEALLOCATE_A(lr%zdl_psi)

    SAFE_DEALLOCATE_A(lr%ddl_rho)
    SAFE_DEALLOCATE_A(lr%zdl_rho)

    SAFE_DEALLOCATE_A(lr%dl_j)
    SAFE_DEALLOCATE_A(lr%ddl_de)
    SAFE_DEALLOCATE_A(lr%ddl_elf)
    SAFE_DEALLOCATE_A(lr%zdl_de)
    SAFE_DEALLOCATE_A(lr%zdl_elf)

    POP_SUB(lr_dealloc)
  end subroutine lr_dealloc


  ! ---------------------------------------------------------
  subroutine lr_copy(st, mesh, src, dest)
    type(states_elec_t), intent(in) :: st
    type(mesh_t),        intent(in) :: mesh
    type(lr_t),          intent(in) :: src
    type(lr_t),       intent(inout) :: dest

    integer :: ik, ist

    PUSH_SUB(lr_copy)

    if(src%is_allocated_rho .and. dest%is_allocated_rho) then
      if(states_are_complex(st)) then
        call lalg_copy(mesh%np, st%d%nspin, src%zdl_rho, dest%zdl_rho)
      else
        call lalg_copy(mesh%np, st%d%nspin, src%ddl_rho, dest%ddl_rho)
      end if
    else
      if(dest%is_allocated_rho) then
        if(states_are_complex(st)) then
          dest%zdl_rho(:, :) = M_ZERO
        else
          dest%ddl_rho(:, :) = M_ZERO
        end if
      end if
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        if(states_are_complex(st)) then
          call lalg_copy(mesh%np_part, st%d%dim, src%zdl_psi(:, :, ist, ik), dest%zdl_psi(:, :, ist, ik))
        else 
          call lalg_copy(mesh%np_part, st%d%dim, src%ddl_psi(:, :, ist, ik), dest%ddl_psi(:, :, ist, ik))
        end if
      end do
    end do

    POP_SUB(lr_copy)

  end subroutine lr_copy



  ! ---------------------------------------------------------
  logical function lr_is_allocated(this) 
    type(lr_t), intent(in) :: this

    PUSH_SUB(lr_is_allocated)
    lr_is_allocated = this%is_allocated

    POP_SUB(lr_is_allocated)
  end function lr_is_allocated



  ! ---------------------------------------------------------
  FLOAT function lr_alpha_j(st, jst, ik) 
    type(states_elec_t), intent(in) :: st
    integer,             intent(in) :: jst
    integer,             intent(in) :: ik

    FLOAT :: dsmear

    PUSH_SUB(lr_alpha_j)

    if(st%smear%method == SMEAR_FIXED_OCC) then
      lr_alpha_j = st%occ(jst, ik) / st%smear%el_per_state
    else
      dsmear = max(CNST(1e-14), st%smear%dsmear)
      lr_alpha_j = max(st%smear%e_fermi + M_THREE*dsmear - st%eigenval(jst, ik), M_ZERO)
    end if

    POP_SUB(lr_alpha_j)
  end function lr_alpha_j

#include "undef.F90"
#include "real.F90"
#include "linear_response_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "linear_response_inc.F90"

end module linear_response_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
