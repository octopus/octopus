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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module linear_response_m
  use comm_m
  use datasets_m
  use global_m
  use grid_m
  use lalg_basic_m
  use loct_m
  use parser_m
  use mesh_m
  use messages_m
  use io_function_m
  use profiling_m
  use smear_m
  use states_m
  use states_calc_m
  use states_dim_m

  implicit none

  private

  public ::               &
       lr_t,              &
       lr_init,           &
       lr_allocate,       &
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
       zlr_swap_sigma


  type lr_t
    !the number of lr wfs
    integer :: nst
    logical :: is_allocated
     
    ! the real quantities
    FLOAT, pointer :: ddl_rho(:,:)     ! response of the density
    FLOAT, pointer :: ddl_psi(:,:,:,:) ! linear change of the real KS orbitals
    
    ! and the complex version
    CMPLX, pointer :: zdl_rho(:,:)     ! response of the density
    CMPLX, pointer :: zdl_psi(:,:,:,:) ! linear change of the complex KS orbitals

    !other observables
    CMPLX, pointer :: dl_j(:,:,:)     ! response of the current
    FLOAT, pointer :: ddl_de(:,:)     ! unnormalized ELF
    FLOAT, pointer :: ddl_elf(:,:)    ! normalized ELF
    CMPLX, pointer :: zdl_de(:,:)     ! unnormalized ELF
    CMPLX, pointer :: zdl_elf(:,:)    ! normalized ELF
    
  end type lr_t

contains

  ! ---------------------------------------------------------
  subroutine lr_init(lr)
    type(lr_t), intent(out) :: lr

    PUSH_SUB(lr_init)

    nullify(lr%ddl_rho, lr%ddl_psi)
    nullify(lr%zdl_rho, lr%zdl_psi)
    nullify(lr%dl_j, lr%ddl_de, lr%zdl_de, lr%ddl_elf, lr%zdl_elf)

    lr%is_allocated = .false.

    POP_SUB(lr_init)

  end subroutine lr_init


  ! ---------------------------------------------------------
  subroutine lr_allocate(lr, st, mesh)
    type(lr_t),     intent(inout) :: lr
    type(states_t), intent(in)    :: st
    type(mesh_t),   intent(in)    :: mesh

    PUSH_SUB(lr_allocate)

    if (states_are_complex(st)) then
      SAFE_ALLOCATE(lr%zdl_psi(1:mesh%np_part, 1:st%d%dim, 1:st%nst, 1:st%d%nik))
      SAFE_ALLOCATE(lr%zdl_rho(1:mesh%np, 1:st%d%nspin))
      
      lr%zdl_psi = M_ZERO
      lr%zdl_rho = M_ZERO
    else
      SAFE_ALLOCATE(lr%ddl_psi(1:mesh%np_part, 1:st%d%dim, 1:st%nst, 1:st%d%nik))
      SAFE_ALLOCATE(lr%ddl_rho(1:mesh%np, 1:st%d%nspin))

      lr%ddl_psi = M_ZERO
      lr%ddl_rho = M_ZERO
    end if

    lr%is_allocated = .true.
    
    POP_SUB(lr_allocate)

  end subroutine lr_allocate


  ! ---------------------------------------------------------
  subroutine lr_dealloc(lr)
    type(lr_t), intent(inout) :: lr

    PUSH_SUB(lr_dealloc)

    SAFE_DEALLOCATE_P(lr%ddl_psi)
    SAFE_DEALLOCATE_P(lr%zdl_psi)

    SAFE_DEALLOCATE_P(lr%ddl_rho)
    SAFE_DEALLOCATE_P(lr%zdl_rho)

    SAFE_DEALLOCATE_P(lr%dl_j)
    SAFE_DEALLOCATE_P(lr%ddl_de)
    SAFE_DEALLOCATE_P(lr%ddl_elf)
    SAFE_DEALLOCATE_P(lr%zdl_de)
    SAFE_DEALLOCATE_P(lr%zdl_elf)

    POP_SUB(lr_dealloc)

  end subroutine lr_dealloc


  ! ---------------------------------------------------------
  subroutine lr_copy(st, mesh, src, dest)
    type(states_t),    intent(in) :: st
    type(mesh_t),      intent(in) :: mesh
    type(lr_t),        intent(in) :: src
    type(lr_t),     intent(inout) :: dest

    integer :: ik, idim, ist

    PUSH_SUB(lr_copy)

    do ik = 1, st%d%nspin
      if(states_are_complex(st)) then
        call lalg_copy(mesh%np, src%zdl_rho(:, ik), dest%zdl_rho(:, ik))
      else
        call lalg_copy(mesh%np, src%ddl_rho(:, ik), dest%ddl_rho(:, ik))
      end if
    enddo

    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          if(states_are_complex(st)) then
            call lalg_copy(mesh%np_part, src%zdl_psi(:, idim, ist, ik), dest%zdl_psi(:, idim, ist, ik))
          else 
            call lalg_copy(mesh%np_part, src%ddl_psi(:, idim, ist, ik), dest%ddl_psi(:, idim, ist, ik))
          end if
        end do
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
    type(states_t), intent(in) :: st
    integer,        intent(in) :: jst
    integer,        intent(in) :: ik

    FLOAT :: dsmear

    PUSH_SUB(lr_alpha_j)

    if(st%smear%method == SMEAR_FIXED_OCC) then
      lr_alpha_j = st%occ(jst, ik) / st%smear%el_per_state
    else
      dsmear = max(CNST(1e-14), st%smear%dsmear)
      lr_alpha_j = max(st%smear%e_fermi + M_THREE*dsmear - st%eigenval(jst, ik), M_ZERO)
    endif

    POP_SUB(lr_alpha_j)
  end function lr_alpha_j

#include "undef.F90"
#include "real.F90"
#include "linear_response_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "linear_response_inc.F90"

end module linear_response_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
