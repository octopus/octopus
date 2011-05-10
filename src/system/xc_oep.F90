!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module xc_oep_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use linear_response_m
  use linear_solver_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use states_dim_m
  use scf_tol_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)
  use xc_functl_m

  implicit none

  private
  public ::                     &
    xc_oep_t,                   &
    xc_oep_init,                &
    xc_oep_end,                 &
    xc_oep_messages_info,          &
    dxc_oep_calc,               &
    zxc_oep_calc

  ! the OEP levels
  integer, public, parameter :: &
    XC_OEP_NONE   = 1,          &
    XC_OEP_SLATER = 2,          &
    XC_OEP_KLI    = 3,          &
    XC_OEP_FULL   = 5

  type xc_oep_t
    integer       :: level      ! 0 = no oep, 1 = Slater, 2 = KLI, 4 = full OEP
    FLOAT         :: mixing     ! how much of the function S(r) to add to vxc in every iteration
    type(lr_t)    :: lr         ! to solve the equation H psi = b
    type(linear_solver_t) :: solver
    type(scf_tol_t) :: scftol
    integer          :: eigen_n
    integer, pointer :: eigen_type(:), eigen_index(:)
    FLOAT            :: socc, sfact
    FLOAT,   pointer :: vxc(:), uxc_bar(:)
    FLOAT,   pointer :: dlxc(:, :)
    CMPLX,   pointer :: zlxc(:, :)
  end type xc_oep_t

  FLOAT, parameter :: small     = CNST(1.0e-5)


contains

  ! ---------------------------------------------------------
  subroutine xc_oep_init(oep, family, gr, d)
    type(xc_oep_t),     intent(out)   :: oep
    integer,            intent(in)    :: family
    type(grid_t),       intent(inout) :: gr
    type(states_dim_t), intent(in)    :: d

    PUSH_SUB(xc_oep_init)

    if(iand(family, XC_FAMILY_OEP).eq.0) then
      oep%level = XC_OEP_NONE
    POP_SUB(xc_oep_init)
      return
    end if

#if defined(HAVE_MPI)
    if(oep%level == XC_OEP_FULL) then
      message(1) = "Full OEP is not allowed with the code parallel in states."
      call messages_fatal(1)
    end if
#endif

    !%Variable OEPLevel
    !%Type integer
    !%Default oep_kli
    !%Section Hamiltonian::XC
    !%Description
    !% At what level shall <tt>Octopus</tt> handle the optimized effective potential (OEP) equation.
    !%Option oep_none 1
    !% Do not solve OEP equation.
    !%Option oep_slater 2
    !% Slater approximation.
    !%Option oep_kli 3
    !% Krieger-Li-Iafrate (KLI) approximation
    !% (JB Krieger, Y Li, GJ Iafrate, <i>Phys. Rev. Lett. A</i> <b>146</b>, 256 (1990).
    !%Option oep_full 5
    !% (Experimental) Full solution of OEP equation using the Sternheimer approach.
    !%End
    call messages_obsolete_variable('OEP_Level', 'OEPLevel')
    call parse_integer(datasets_check('OEPLevel'), XC_OEP_KLI, oep%level)
    if(.not. varinfo_valid_option('OEPLevel', oep%level)) call input_error('OEP_level')

    if(oep%level .ne. XC_OEP_NONE) then
      if(oep%level == XC_OEP_FULL) then
        call messages_experimental("Full OEP")    
        !%Variable OEPMixing
        !%Type float
        !%Default 1.0
        !%Section Hamiltonian::XC
        !%Description
        !% The linear mixing factor used to solve the Sternheimer
        !% equation in the full OEP procedure. The default is 1.0.
        !%End
        call messages_obsolete_variable('OEP_Mixing', 'OEPMixing')
        call parse_float(datasets_check('OEPMixing'), M_ONE, oep%mixing)
      end if

      ! this routine is only prepared for finite systems, and ispin = 1, 2
      if(d%ispin > SPIN_POLARIZED .or. d%nik>d%ispin) then
        message(1) = "OEP only works for finite systems and collinear spin."
        call messages_fatal(1)
      end if
    
      ! obtain the spin factors
      call xc_oep_SpinFactor(oep, d%nspin)

      ! This variable will keep vxc across iterations
      SAFE_ALLOCATE(oep%vxc(1:gr%mesh%np))

      ! when performing full OEP, we need to solve a linear equation
      if(oep%level == XC_OEP_FULL) then 
        call scf_tol_init(oep%scftol, "OEP",def_maximumiter=10)
        call linear_solver_init(oep%solver, gr, "OEP")
        call lr_init(oep%lr)
      end if

      ! the linear equation has to be more converged if we are to attain the required precision
      !oep%lr%conv_abs_dens = oep%lr%conv_abs_dens / (oep%mixing)
    end if

    POP_SUB(xc_oep_init)
  end subroutine xc_oep_init


  ! ---------------------------------------------------------
  subroutine xc_oep_end(oep)
    type(xc_oep_t), intent(inout) :: oep

    PUSH_SUB(xc_oep_end)

    if(oep%level.ne.XC_OEP_NONE) then
      SAFE_DEALLOCATE_P(oep%vxc)

      if(oep%level == XC_OEP_FULL) then 
        call lr_dealloc(oep%lr)
        call linear_solver_end(oep%solver)
      end if
    end if

    POP_SUB(xc_oep_end)
  end subroutine xc_oep_end


  ! ---------------------------------------------------------
  subroutine xc_oep_messages_info(oep, iunit)
    type(xc_oep_t), intent(in) :: oep
    integer,        intent(in) :: iunit

    if(oep%level.eq.XC_OEP_NONE) return

    PUSH_SUB(xc_oep_messages_info)
    call messages_print_var_option(iunit, 'OEPLevel', oep%level)

    POP_SUB(xc_oep_messages_info)
  end subroutine xc_oep_messages_info


  ! ---------------------------------------------------------
  ! A couple of auxiliary functions for oep
  ! ---------------------------------------------------------
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_t), intent(inout) :: oep
    integer,        intent(in)    :: nspin

    PUSH_SUB(xc_oep_SpinFactor)

    select case(nspin)
    case(1) ! we need to correct or the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case(2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default ! cannot handle any other case
      ASSERT(.false.)
    end select

    POP_SUB(xc_oep_SpinFactor)
  end subroutine xc_oep_SpinFactor


  ! ---------------------------------------------------------
  subroutine xc_oep_AnalyzeEigen(oep, st, is)
    type(xc_oep_t), intent(inout) :: oep
    type(states_t), intent(in)    :: st
    integer,        intent(in)    :: is

    integer  :: ist
    FLOAT :: max_eigen
    FLOAT, allocatable :: eigenval(:), occ(:)

    PUSH_SUB(xc_oep_AnalyzeEigen)

    SAFE_ALLOCATE(eigenval(1:st%nst))
    SAFE_ALLOCATE     (occ(1:st%nst))
    eigenval = M_ZERO
    occ = M_ZERO

    do ist = st%st_start, st%st_end
      eigenval(ist) = st%eigenval(ist, is)
      occ(ist) = st%occ(ist, is)
    end do

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
      do ist = 1, st%nst
        call MPI_Bcast(eigenval(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
        call MPI_Bcast(occ(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    ! find out the top occupied state, to correct for the asymptotics
    ! of the potential
    max_eigen = CNST(-1e30)
    do ist = 1, st%nst
      if((occ(ist) .gt. small).and.(eigenval(ist).gt.max_eigen)) then
        max_eigen = eigenval(ist)
      end if
    end do

    oep%eigen_n = 1
    do ist = 1, st%nst
      if(occ(ist) .gt. small) then
        ! criterion for degeneracy
        if(abs(eigenval(ist)-max_eigen).le.CNST(1e-3)) then
          oep%eigen_type(ist) = 2
        else
          oep%eigen_type(ist) = 1
          oep%eigen_index(oep%eigen_n) = ist
          oep%eigen_n = oep%eigen_n + 1
        end if
      else
        oep%eigen_type(ist) = 0
      end if
    end do
    oep%eigen_n = oep%eigen_n - 1

    SAFE_DEALLOCATE_A(eigenval)
    SAFE_DEALLOCATE_A(occ)
    POP_SUB(xc_oep_AnalyzeEigen)
  end subroutine xc_oep_AnalyzeEigen


#include "undef.F90"
#include "real.F90"
#include "xc_kli_inc.F90"
#include "xc_oep_x_inc.F90"
#include "xc_oep_sic_inc.F90"
#include "xc_oep_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_kli_inc.F90"
#include "xc_oep_x_inc.F90"
#include "xc_oep_sic_inc.F90"
#include "xc_oep_inc.F90"

#include "undef.F90"

end module xc_oep_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
