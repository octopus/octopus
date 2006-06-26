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
  use output_m 
  use units_m
  use mesh_function_m
  use lib_oct_m
  use io_m
  use lib_basic_alg_m

  implicit none

  type lr_t
    FLOAT   :: conv_abs_dens  ! convergence required
    FLOAT   :: abs_dens       ! convergence reached
    integer :: max_iter       ! maximum number of iterations
    integer :: iter           ! number of iterations used
    integer :: solver         !the linear solver to use
    integer :: preconditioner !the linear solver to use
    logical :: ort_each_step  

    type(mix_t) :: mixer   ! can not live without it

    ! the real quantities
    FLOAT, pointer :: ddl_rho(:,:)     ! response of the density
    FLOAT, pointer :: ddl_psi(:,:,:,:) ! linear change of the real KS orbitals
    FLOAT, pointer :: ddl_Vhar(:)      ! linear change of the Hartree potential

    ! and the complex version
    CMPLX, pointer :: zdl_rho(:,:)     ! response of the density
    CMPLX, pointer :: zdl_psi(:,:,:,:) ! linear change of the complex KS orbitals
    CMPLX, pointer :: zdl_Vhar(:)      ! linear change of the Hartree potential

    FLOAT, pointer :: dl_Vxc(:,:,:)    ! linear change of the xc potential (fxc)
    
    !other observables
    CMPLX, pointer :: dl_j(:,:,:)     ! response of the current
    FLOAT, pointer :: ddl_de(:,:)     ! unnormalized elf
    FLOAT, pointer :: ddl_elf(:,:)    ! normalized elf
    CMPLX, pointer :: zdl_de(:,:)     ! unnormalized elf
    CMPLX, pointer :: zdl_elf(:,:)    ! normalized elf
    
  end type lr_t

  FLOAT, parameter :: lr_min_occ=CNST(1e-4) !the minimum value for a state to be considered occupied

  !solvers ( values must be < 100 )
  integer, parameter :: LR_CG = 5
  integer, parameter :: LR_BCG = 2
  integer, parameter :: LR_BICGSTAB = 3

  !preconditioners ( values must be multiples of 100 )
  integer, parameter :: LR_NONE = 0
  integer, parameter :: LR_DIAG = 100

contains

  ! ---------------------------------------------------------
  subroutine lr_init(lr, prefix, def_solver)
    type(lr_t),       intent(out) :: lr
    character(len=*), intent(in)  :: prefix
    integer, optional, intent(in) :: def_solver

    integer :: fsolver

    call push_sub('linear_response.lr_init')

    ! read some parameters from the input file
    call loct_parse_float(check_inp(trim(prefix)//"ConvAbsDens"), &
        CNST(1e-5), lr%conv_abs_dens)
    call loct_parse_int  (check_inp(trim(prefix)//"MaximumIter"), 1000, lr%max_iter)

    
    !%Variable LinearSolver
    !%Type integer
    !%Default cg
    !%Section Linear Response::Polarizabities
    !%Description
    !% To calculate response using density functional perturbation
    !% theory is necessary to solve the sterheimer equation, this is a self
    !% consistent linear equation where the operator is the shifted Kohn-Sham hamiltonian.
    !% This variable which method to use in order to solve this linear equation.
    !% An optional preconditioner can be added.
    !%Option cg 5
    !% Conjugated gradients. This is the faster solver but does not
    !& work when a imaginary part when an imaginary shift is added.
    !%Option bcg 2
    !% Biconjugated gradients. This solver is a generalization of the
    !% conjugated gradients for non-hermitian operators. It has some
    !% stability problems, so the bigstab should be prefered.
    !%Option bicgstab 3
    !% Biconjugated gradients stabilized. This is an improved version
    !% of bcg that is faster and more stable. This is the default when
    !% complex response is calculated.
    !%Option diag_prec 100
    !% Preconditioning using the diagonal of the operator. 
    !%End

    if(present(def_solver)) then
      call loct_parse_int  (check_inp(trim(prefix)//"LinearSolver"), def_solver, fsolver)
    else 
      call loct_parse_int  (check_inp(trim(prefix)//"LinearSolver"), LR_CG, fsolver)
    end if

    !the last 2 digits select the linear solver
    lr%solver = mod(fsolver, 100)

    !the next 2 digits select the preconditioner
    lr%preconditioner = (fsolver - lr%solver)

    call loct_parse_logical(check_inp(trim(prefix)//"OrtEachStep"), .false., lr%ort_each_step)

    nullify(lr%ddl_rho, lr%ddl_psi, lr%ddl_Vhar, lr%dl_Vxc)
    nullify(lr%zdl_rho, lr%zdl_psi, lr%zdl_Vhar, lr%dl_Vxc)
    
    nullify(lr%dl_j, lr%ddl_de, lr%zdl_de, lr%ddl_elf, lr%zdl_elf)
    
    call pop_sub()

  end subroutine lr_init

  ! ---------------------------------------------------------
  integer function dlr_alloc_psi(st, m, lr) result(r)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: m
    type(lr_t),  intent(inout) :: lr

    call push_sub('linear_response.dlr_alloc_psi')

    r = 1
    if(associated(lr%ddl_psi)) return ! do nothing

    ALLOCATE(lr%ddl_psi(m%np_part, st%d%dim, st%nst, st%d%nspin),
         m%np_part*st%d%dim*st%nst*st%d%nspin)

    if(associated(lr%zdl_psi)) then
      r = 2
      lr%ddl_psi = real(lr%zdl_psi, PRECISION)
      deallocate(lr%zdl_psi)
      nullify(lr%zdl_psi)
    else
      r = -1
      lr%ddl_psi = M_ZERO
    end if

    call pop_sub()

  end function dlr_alloc_psi


  ! ---------------------------------------------------------
  integer function zlr_alloc_psi(st, m, lr) result(r)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: m
    type(lr_t),  intent(inout) :: lr

    call push_sub('linear_response.zlr_alloc_psi')

    r = 1
    if(associated(lr%zdl_psi)) return ! do nothing

    ALLOCATE(lr%zdl_psi(m%np_part, st%d%dim, st%nst, st%d%nspin),
         m%np_part*st%d%dim*st%nst*st%d%nspin)

    if(associated(lr%ddl_psi)) then
      r = 2
      lr%zdl_psi = lr%ddl_psi
      deallocate(lr%ddl_psi)
      nullify(lr%ddl_psi)
    else
      r = -1
      lr%zdl_psi = M_z0
    end if

    call pop_sub()

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

    if(associated(lr%dl_j)) deallocate(lr%dl_j)
    if(associated(lr%ddl_de)) deallocate(lr%ddl_de)
    if(associated(lr%ddl_elf)) deallocate(lr%ddl_elf)
    if(associated(lr%zdl_de)) deallocate(lr%zdl_de)
    if(associated(lr%zdl_elf)) deallocate(lr%zdl_elf)

  end subroutine lr_dealloc


  ! ---------------------------------------------------------
  subroutine lr_build_fxc(m, st, xcs, fxc)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    type(xc_t),     intent(in)  :: xcs
    FLOAT,          intent(inout) :: fxc(:,:,:)

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
    deallocate(rho)

    call pop_sub()
  end subroutine lr_build_fxc


  ! ---------------------------------------------------------
  subroutine lr_alloc_fHxc(st, m, lr)
    type(states_t), intent(in)  :: st
    type(mesh_t),   intent(in)  :: m
    type(lr_t),     intent(inout) :: lr

    call push_sub('linear_response.lr_alloc_fHxc')

    lr%abs_dens = M_ZERO
    lr%iter     = 0

    ! allocate variables
    if (st%d%wfs_type == M_REAL) then
      ALLOCATE(lr%ddl_rho(m%np, st%d%nspin), m%np*st%d%nspin)
      ALLOCATE(lr%ddl_Vhar(m%np), m%np)
    else
      ALLOCATE(lr%zdl_rho(m%np, st%d%nspin), m%np*st%d%nspin)
      ALLOCATE(lr%zdl_Vhar(m%np), m%np)
    end if
    ALLOCATE(lr%dl_Vxc(m%np, st%d%nspin, st%d%nspin), m%np*st%d%nspin*st%d%nspin)

    call pop_sub()

  end subroutine lr_alloc_fHxc


#include "undef.F90"
#include "real.F90"

#include "linear_response_inc.F90"
#include "linear_response_out.F90" 
#include "linear_response_solvers.F90" 

#include "undef.F90"
#include "complex.F90"

#include "linear_response_inc.F90"
#include "linear_response_out.F90" 
#include "linear_response_solvers.F90" 

end module linear_response_m
