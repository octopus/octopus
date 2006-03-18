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
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use mesh_m
  use grid_m
  use states_m
  use mix_m
  use hamiltonian_m
  use xc_m
  use functions_m

  implicit none

  type lr_t
    FLOAT   :: conv_abs_dens  ! convergence required
    FLOAT   :: abs_dens       ! convergence reached
    integer :: max_iter       ! maximum number of iterations
    integer :: iter           ! number of iterations used

    type(mix_t) :: mixer   ! can not live without it

    ! the real quantities
    FLOAT, pointer :: ddl_rho(:,:)     ! response of the density
    FLOAT, pointer :: ddl_psi(:,:,:,:) ! linear change of the real KS orbitals
    FLOAT, pointer :: ddl_Vhar(:)      ! linear change of the Hartree potential

    ! and the complex version
    CMPLX, pointer :: zdl_rho(:,:)     ! response of the density
    CMPLX, pointer :: zdl_psi(:,:,:,:) ! linear change of the real KS orbitals
    CMPLX, pointer :: zdl_Vhar(:)      ! linear change of the complex KS orbitals

    FLOAT, pointer :: dl_Vxc(:,:,:)    ! linear change of the xc potential (fxc)
  end type lr_t

contains

  ! ---------------------------------------------------------
  subroutine lr_init(lr, prefix)
    type(lr_t),       intent(out) :: lr
    character(len=*), intent(in)  :: prefix

    ! read some parameters from the input file
    call loct_parse_float(check_inp(trim(prefix)//"ConvAbsDens"), &
        CNST(1e-5), lr%conv_abs_dens)
    call loct_parse_int  (check_inp(trim(prefix)//"MaximumIter"), 50, lr%max_iter)

    nullify(lr%ddl_rho, lr%ddl_psi, lr%ddl_Vhar, lr%dl_Vxc)
    nullify(lr%zdl_rho, lr%zdl_psi, lr%zdl_Vhar, lr%dl_Vxc)
  end subroutine lr_init


  ! ---------------------------------------------------------
  integer function dlr_alloc_psi(st, m, lr) result(r)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: m
    type(lr_t),  intent(inout) :: lr

    r =  1
    if(associated(lr%ddl_psi)) return ! do nothing

    ALLOCATE(lr%ddl_psi(m%np_part, st%d%dim, st%nst, st%d%nspin), m%np_part*st%d%dim*st%nst*st%d%nspin)

    if(associated(lr%zdl_psi)) then
      r = 2
      lr%ddl_psi = real(lr%zdl_psi, PRECISION)
      deallocate(lr%zdl_psi)
      nullify(lr%zdl_psi)
    else
      r = -1
      lr%ddl_psi = M_ZERO
    end if

  end function dlr_alloc_psi


  ! ---------------------------------------------------------
  integer function zlr_alloc_psi(st, m, lr) result(r)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: m
    type(lr_t),  intent(inout) :: lr

    r =  1
    if(associated(lr%zdl_psi)) return ! do nothing

    ALLOCATE(lr%zdl_psi(m%np_part, st%d%dim, st%nst, st%d%nspin), m%np_part*st%d%dim*st%nst*st%d%nspin)

    if(associated(lr%ddl_psi)) then
      r = 2
      lr%zdl_psi = lr%ddl_psi
      deallocate(lr%ddl_psi)
      nullify(lr%ddl_psi)
    else
      r = -1
      lr%zdl_psi = M_z0
    end if

  end function zlr_alloc_psi


  ! ---------------------------------------------------------
  subroutine lr_dealloc(lr)
    type(lr_t), intent(inout) :: lr

    if(associated(lr%ddl_rho)) then
      deallocate(lr%ddl_rho, lr%ddl_Vhar, lr%dl_Vxc)
      nullify   (lr%ddl_rho, lr%ddl_Vhar, lr%dl_Vxc)
    end if

    if(associated(lr%zdl_rho)) then
      deallocate(lr%zdl_rho, lr%zdl_Vhar, lr%dl_Vxc)
      nullify   (lr%zdl_rho, lr%zdl_Vhar, lr%dl_Vxc)
    end if

    if(associated(lr%ddl_psi)) then
      deallocate(lr%ddl_psi)
      nullify   (lr%ddl_psi)
    end if

    if(associated(lr%zdl_psi)) then
      deallocate(lr%zdl_psi)
      nullify   (lr%zdl_psi)
    end if

  end subroutine lr_dealloc


  ! ---------------------------------------------------------
  subroutine lr_build_fxc(m, st, xcs, fxc)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    type(xc_t),     intent(in)  :: xcs
    FLOAT,          intent(out) :: fxc(:,:,:)

    FLOAT, allocatable :: rho(:, :)
    integer :: is

    call push_sub('linear_response.lr_build_fxc')

    ALLOCATE(rho(m%np, st%d%nspin), m%np*st%d%nspin)
    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
      end do
    else
      rho = st%rho
    end if
    fxc = M_ZERO
    call xc_get_fxc(xcs, m, rho, st%d%ispin, fxc)

    call pop_sub()
  end subroutine lr_build_fxc



#include "undef.F90"
#include "real.F90"
#include "linear_response_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "linear_response_inc.F90"


end module linear_response_m
