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

#include "global.h"

module linear_response
  use global
  use mesh
  use states
  use mix
  use hamiltonian
  use xc
  use functions

  implicit none

  type lr_type
    FLOAT   :: conv_abs_dens  ! convergence required
    FLOAT   :: abs_dens       ! convergence reached
    integer :: max_iter       ! maximum number of iterations
    integer :: iter           ! number of iterations used
    
    type(mix_type) :: mixer   ! can not live without it

    ! the real quantities
    FLOAT, pointer :: ddl_rho(:,:)     ! response of the density
    FLOAT, pointer :: ddl_psi(:,:,:,:) ! linear change of the real KS orbitals
    FLOAT, pointer :: ddl_Vhar(:)      ! linear change of the Hartree potential

    ! and the complex version
    CMPLX, pointer :: zdl_rho(:,:)     ! response of the density
    CMPLX, pointer :: zdl_psi(:,:,:,:) ! linear change of the real KS orbitals
    CMPLX, pointer :: zdl_Vhar(:)      ! linear change of the complex KS orbitals
    
    FLOAT, pointer :: dl_Vxc(:,:,:)    ! linear change of the xc potential (fxc)
  end type lr_type

contains

  ! ---------------------------------------------------------
  subroutine lr_init(lr)
    type(lr_type), intent(out) :: lr

    ! read some parameters from the input file
    call loct_parse_float("ConvAbsDens", CNST(1e-5), lr%conv_abs_dens)
    call loct_parse_int("MaximumIter", 50, lr%max_iter)    
  end subroutine lr_init


  ! ---------------------------------------------------------
  subroutine lr_build_fxc(m, st, xcs, fxc)
    type(mesh_type),   intent(in)  :: m
    type(states_type), intent(in)  :: st
    type(xc_type),     intent(in)  :: xcs
    FLOAT,             intent(out) :: fxc(:,:,:)

    FLOAT, allocatable :: rho(:, :)
    integer :: is

    call push_sub('lr_build_fxc')

    allocate(rho(m%np, st%d%nspin))
    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
      enddo
    else
      rho = st%rho
    endif
    fxc = M_ZERO
    call xc_get_fxc(xcs, m, rho, st%d%ispin, fxc)

    call pop_sub()
  end subroutine lr_build_fxc


  ! ---------------------------------------------------------
  ! orthogonalizes response of \alpha KS orbital to all occupied 
  ! \alpha KS orbitals  
  ! ---------------------------------------------------------
  subroutine lr_orth_response(m, st, lr)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(in)    :: st
    type(lr_type),     intent(inout) :: lr
    
    integer :: ist, ik
    call push_sub('lr_orth_response')

    do ik = 1, st%d%nspin
      do ist = 1, st%nst 
        if(st%occ(ist, ik) > M_ZERO) then  
          call X(lr_orth_vector) (m, st, lr%X(dl_psi)(:,:, ist, ik), ik) 
        endif
      end do
    end do
    
    call pop_sub()
  end subroutine lr_orth_response


#include "undef.F90"
#include "real.F90"
#include "linear_response_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "linear_response_inc.F90"


end module linear_response
