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
    FLOAT   :: conv_abs_psi   ! convergence required 
    FLOAT   :: abs_psi        ! convergence reached
    integer :: max_iter       ! maximum number of iterations
    integer :: max_scf_iter   ! maximum number of iterations
    integer :: iter           ! number of iterations used
    integer :: solver         ! the linear solver to use
    integer :: preconditioner ! the preconditioner solver to use
    integer :: ort_min_step   ! the step where to start orthogonalization

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
  integer, parameter :: LR_HX_FIXED = 11
  integer, parameter :: LR_HX = 10
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

    ! BEGIN INPUT PARSING

    !%Variable ConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Linear Response::SCF
    !%Description
    !% The tolerance in the variation of the density, to determine if
    !% the SCF for linear response is converged.
    !%End

    call loct_parse_float(check_inp(trim(prefix)//"ConvAbsDens"), &
        CNST(1e-5), lr%conv_abs_dens)
    
    !%Variable MaximumIter
    !%Type integer
    !%Default 200
    !%Section Linear Response::SCF
    !%Description
    !% The maximum number of SCF iterations to calculate response.
    !%End

    call loct_parse_int(check_inp('PolSCFIterations'), 200, lr%max_scf_iter)

    !%Variable LinearSolver
    !%Type integer
    !%Default cg
    !%Section Linear Response::Solver
    !%Description
    !% To calculate response using density functional perturbation
    !% theory is necessary to solve the sterheimer equation, this is a self
    !% consistent linear equation where the operator is the Kohn-Sham
    !% hamiltonian with a complex shift. This variable selects which
    !% method to use in order to solve this linear equation.
    !%Option hx 10
    !% This solver aproximates the solution of the non-hermitian
    !% equation by a power series in terms of the inverse of the
    !% hamiltonian. Each term implies solving a hermitian linear equation by
    !% conjugated gradients.
    !% Altough this might sound inefficient, the hermitian equations
    !% to solve are better conditioned than the original, so this can be more
    !% efficient than a non-hermitian solver.
    !% This version of the solvers apply the number of terms in the
    !% series as needed to be under the tolerance required.
    !%Option hx_fixed 11
    !% This is the same prevoius solver, but only two steps of the
    !% series are applied.
    !%Option cg 5
    !% Conjugated gradients. This is the fastest solver but does not
    !% work when an imaginary shift is added. This is the default.
    !%Option bcg 2
    !% Biconjugated gradients. This solver is a generalization of the
    !% conjugated gradients for non-hermitian operators. It has some
    !% stability problems, so the bigstab should be prefered.
    !%Option bicgstab 3
    !% Biconjugated gradients stabilized. This is an improved version
    !% of bcg that is faster and more stable. This is the default when
    !% complex polarizabilities are calculated.
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


    !%Variable LinearSolverMaxIter
    !%Type integer
    !%Default 1000
    !%Section Linear Response::Solver
    !%Description
    !% Maximum number of iterations the linear solver does, even if
    !% convergency is not achieved.
    !%End

    call loct_parse_int  (check_inp(trim(prefix)//"LinearSolverMaxIter"), 1000, lr%max_iter)

    !%Variable LinearSolverTol
    !%Type float
    !%Default 1e-6
    !%Section Linear Response::Solver
    !%Description
    !% This is the tolerance to determine that the linear solver has converged.
    !%End
    
    call loct_parse_float(check_inp(trim(prefix)//"LinearSolverTol"), &
        CNST(1e-6), lr%conv_abs_psi)
    
    !%Variable LinearSolverOrthogonalization
    !%Type integer
    !%Default 50
    !%Section Linear Response::Solver
    !%Description
    !% A good preconditioner to the Sternheimer equation is the
    !% projector onto the unoccupied state.
    !% The problem is that this operator is expensive to apply because it
    !% requires to orthogonalize with respect to all the occupied
    !% wavenfunctions.
    !% To overcome this, the operator is only applied is the linear solver
    !% has not converged after a certain number of steps. This variable
    !% controls that number. The default is to start to orthogonalize
    !% after 50 steps of linear solver. A value of 1 will activate
    !% orthogonalization always and 0 means that it will not be used
    !% at all.
    !%Option always 1 
    !% The ortogonalization preconditioner will be applied always.
    !%Option never 0 
    !% The ortogonalization will not be used.
    !%End

    call loct_parse_int(check_inp(trim(prefix)//"LinearSolverOrthogonalization"), 50, lr%ort_min_step)
    
    ! END INPUT PARSING

    nullify(lr%ddl_rho, lr%ddl_psi, lr%ddl_Vhar, lr%dl_Vxc)
    nullify(lr%zdl_rho, lr%zdl_psi, lr%zdl_Vhar, lr%dl_Vxc)
    
    nullify(lr%dl_j, lr%ddl_de, lr%zdl_de, lr%ddl_elf, lr%zdl_elf)
    
    !WRITE INFO
    
    write(message(1),'(a)') 'Linear Reponse'
    call messages_print_stress(stdout, trim(message(1)))
    
    
    ! solver 
    select case(lr%solver)
      case(LR_CG)
        message(1)='Linear Solver: Conjugated Gradients'

      case(LR_BCG)
        message(1)='Linear Solver: Biconjugated Gradients'

      case(LR_BICGSTAB)
        message(1)='Linear Solver: Biconjugated Gradients Stabilized'

      case(LR_HX_FIXED)
        message(1)='Linear Solver: Fixed Hermitian Expansion'

      case(LR_HX)
        message(1)='Linear Solver: Hermitian Expansion'

    end select

    call write_info(1)
    
    call messages_print_stress(stdout)

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
      rho(1:m%np, 1:st%d%nspin) = st%rho(1:m%np, 1:st%d%nspin)
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

  logical function precondition(lr) result(prec)
    type(lr_t),     intent(in) :: lr
    prec=(lr%preconditioner /= 0)
  end function precondition


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
