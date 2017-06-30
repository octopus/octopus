!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$

#include "global.h"

module orbital_set_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use global_oct_m
  use kpoints_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use types_oct_m  
 
  implicit none

  private

  public ::                             &
       orbital_set_t,                   &
       orbital_set_nullify,             &
       orbital_set_init,                &
       orbital_set_end,                 &
       orbital_set_update_phase,        &
       dorbital_set_get_coefficients,   &
       zorbital_set_get_coefficients,   &
       dorbital_set_get_coeff_batch,    &
       zorbital_set_get_coeff_batch,    &
       dorbital_set_add_to_batch,       &
       zorbital_set_add_to_batch,       &
       dorbital_set_add_to_psi,         &
       zorbital_set_add_to_psi

  type orbital_set_t
    integer             :: nn, ll
    integer             :: norbs
    integer             :: iatom 
    type(submesh_t)     :: sphere             !> The submesh of the orbital
    CMPLX, pointer      :: phase(:,:)         !> Correction to the global phase 
                                              !> if the sphere cross the border of the box
    FLOAT               :: Ueff               !> The effective U of the simplified rotational invariant form
    FLOAT               :: Ubar, Jbar
    type(species_t), pointer :: spec          

    FLOAT, pointer      :: dS(:,:,:)             !> Overlap matrix for each orbital with similar orbital on other atomic sites    
    CMPLX, pointer      :: zS(:,:,:)

    FLOAT, pointer      :: dorb(:,:) !> The orbital, if real, on the submesh
    CMPLX, pointer      :: zorb(:,:) !> The orbital, if complex, on the submesh
    CMPLX, pointer      :: eorb(:,:,:) !> Orbitals with its phase factor

    type(poisson_t)  :: poisson               !> For computing the Coulomb integrals
  end type orbital_set_t

contains

 subroutine orbital_set_nullify(this)
  type(orbital_set_t),             intent(inout) :: this

  PUSH_SUB(orbital_set_nullify)

  nullify(this%phase)
  nullify(this%dS)
  nullify(this%zS)
  nullify(this%spec)
  nullify(this%dorb)
  nullify(this%zorb)
  nullify(this%eorb)

  POP_SUB(orbital_set_nullify)

 end subroutine orbital_set_nullify

 subroutine orbital_set_init(this)
  type(orbital_set_t),             intent(inout) :: this

  PUSH_SUB(orbital_set_init)

  POP_SUB(orbital_set_init)
 end subroutine orbital_set_init


 subroutine orbital_set_end(this)
   implicit none
   type(orbital_set_t), intent(inout) :: this

   integer :: iorb
  
   PUSH_SUB(orbital_set_end)  

   SAFE_DEALLOCATE_P(this%phase)
   SAFE_DEALLOCATE_P(this%dS)
   SAFE_DEALLOCATE_P(this%zS)
   SAFE_DEALLOCATE_P(this%dorb)
   SAFE_DEALLOCATE_P(this%zorb)
   SAFE_DEALLOCATE_P(this%eorb)
   nullify(this%spec)
   call submesh_end(this%sphere)
   
   POP_SUB(orbital_set_end)
 end subroutine orbital_set_end

  !> Build the phase correction to the global phase in case the orbital crosses the border of the simulaton box
  subroutine orbital_set_update_phase(os, sb, std, vec_pot, vec_pot_var)
    type(orbital_set_t),           intent(inout) :: os
    type(simul_box_t),             intent(in)    :: sb
    type(states_dim_t),            intent(in)    :: std
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (sb%dim)
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)

    integer :: ns, iq, is, ikpoint, im
    FLOAT   :: kr, kpoint(1:MAX_DIM)
    integer :: ndim

    PUSH_SUB(orbital_set_update_phase)

    ns = os%sphere%np
    ndim = sb%dim

    do iq = std%kpt%start, std%kpt%end
      ikpoint = states_dim_get_kpoint_index(std, iq)

      ! if this fails, it probably means that sb is not compatible with std
      ASSERT(ikpoint <= kpoints_number(sb%kpoints))

      kpoint = M_ZERO
      kpoint(1:ndim) = kpoints_get_point(sb%kpoints, ikpoint)

      do is = 1, ns
        ! this is only the correction to the global phase, that can
        ! appear if the sphere crossed the boundary of the cell.
        kr = sum(kpoint(1:ndim)*(os%sphere%x(is, 1:ndim) - os%sphere%mesh%x(os%sphere%map(is), 1:ndim)))

        if(present(vec_pot)) then
          if(allocated(vec_pot)) kr = kr + &
            sum(vec_pot(1:ndim)*(os%sphere%x(is, 1:ndim)- os%sphere%mesh%x(os%sphere%map(is), 1:ndim)))
        end if

        if(present(vec_pot_var)) then
          if(allocated(vec_pot_var)) kr = kr + sum(vec_pot_var(1:ndim, os%sphere%map(is))*os%sphere%x(is, 1:ndim))
        end if

        os%phase(is, iq) = exp(M_zI*kr)
      end do

      !We now compute the so-called Bloch sum of the localized orbitals
      do im = 1, os%norbs
        do is = 1, ns
          os%eorb(is,im,iq) = os%zorb(is,im)*os%phase(is, iq)
        end do
      end do
    end do

    POP_SUB(orbital_set_update_phase)
  end subroutine orbital_set_update_phase

#include "undef.F90"
#include "real.F90"
#include "orbital_set_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbital_set_inc.F90"

end module orbital_set_oct_m
