!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module xc_OEP_m
  use global_m
  use messages_m
  use states_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use lib_xc_m
  use xc_m
  use mesh_m
  use functions_m
  use mesh_function_m
  use poisson_m
  use hamiltonian_m
  use linear_response_m
  use grid_m
  use mpi_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    xc_oep_t,                   &
    xc_oep_init,                &
    xc_oep_end,                 &
    xc_oep_write_info,          &
    dxc_oep_calc,               &
    zxc_oep_calc

  ! the OEP levels
  integer, public, parameter :: &
    XC_OEP_NONE   = 1,          &
    XC_OEP_SLATER = 2,          &
    XC_OEP_KLI    = 3,          &
    XC_OEP_CEDA   = 4,          &  ! not yet implemented
    XC_OEP_FULL   = 5              ! half implemented

  type xc_oep_t
    integer          :: level      ! 0 = no oep, 1 = Slater, 2 = KLI, 3 = CEDA, 4 = full OEP
    FLOAT            :: mixing     ! how much of the function S(r) to add to vxc in every iteration
    type(lr_t)    :: lr         ! to solve the equation H psi = b

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
  subroutine xc_oep_init(oep, family, m, d)
    type(xc_oep_t),     intent(out) :: oep
    integer,               intent(in)  :: family
    type(mesh_t),       intent(in)  :: m
    type(states_dim_t), intent(in)  :: d

    if(iand(family, XC_FAMILY_OEP).eq.0) then
      oep%level = XC_OEP_NONE
      return
    end if

#if defined(HAVE_MPI)
    if(oep%level == XC_OEP_FULL) then
      message(1) = "Full OEP is not allowed with the code parallelized on orbitals..."
      call write_fatal(1)
    end if
#endif

    !%Variable OEP_level
    !%Type integer
    !%Default oep_none
    !%Section Hamiltonian::XC
    !%Description
    !% At what level shall octopus handle the OEP equation
    !%Option oep_none 1
    !% Do not solve OEP equation
    !%Option oep_slater 2
    !% Slater approximation
    !%Option oep_kli 3
    !% KLI approximation
    !%Option oep_ceda 4
    !% CEDA (common energy denominator) approximation (not implemented)
    !%Option oep_full 5
    !% Full solution of OEP equation using the approach of S. Kuemmel (half implemented)
    !%End
    call loct_parse_int(check_inp('OEP_level'), XC_OEP_KLI, oep%level)
    if(.not.varinfo_valid_option('OEP_level', oep%level)) call input_error('OEP_level')

    if(oep%level == XC_OEP_FULL) then
      call loct_parse_float(check_inp('OEP_mixing'), M_ONE, oep%mixing)
    end if

    ! this routine is only prepared for finite systems, and ispin = 1, 2
    if(d%ispin > SPIN_POLARIZED .or. d%nik>d%ispin) then
      message(1) = "OEP only works for finite systems and collinear spin!"
      call write_fatal(1)
    end if

    ! obtain the spin factors
    call xc_oep_SpinFactor(oep, d%nspin)

    ! This variable will keep vxc across iterations
    ALLOCATE(oep%vxc(m%np), m%np)

    ! when performing full OEP, we need to solve a linear equation
    if(oep%level == XC_OEP_FULL) call lr_init(oep%lr, "OEP")

    ! the linear equation has to be more converged if we are to attain the required precision
    !oep%lr%conv_abs_dens = oep%lr%conv_abs_dens / (oep%mixing)

  end subroutine xc_oep_init


  ! ---------------------------------------------------------
  subroutine xc_oep_end(oep)
    type(xc_oep_t), intent(inout) :: oep

    if(oep%level.ne.XC_OEP_NONE) then
      deallocate(oep%vxc); nullify(oep%vxc)

      if(oep%level == XC_OEP_FULL) call lr_dealloc(oep%lr)
    end if

  end subroutine xc_oep_end


  ! ---------------------------------------------------------
  subroutine xc_oep_write_info(oep, iunit)
    type(xc_oep_t), intent(in) :: oep
    integer,           intent(in) :: iunit

    if(oep%level.eq.XC_OEP_NONE) return

    call messages_print_var_option(iunit, 'OEP_level', oep%level)

  end subroutine xc_oep_write_info


  ! ---------------------------------------------------------
  ! A couple of auxiliary functions for oep
  ! ---------------------------------------------------------
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_t), intent(inout) :: oep
    integer,           intent(in)    :: nspin

    select case(nspin)
    case(1) ! we need to correct or the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case(2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default ! can not handle any other case
      ASSERT(.false.)
    end select

  end subroutine xc_oep_SpinFactor


  ! ---------------------------------------------------------
  subroutine xc_oep_AnalizeEigen(oep, st, is)
    type(xc_oep_t), intent(inout) :: oep
    type(states_t), intent(in)    :: st
    integer,           intent(in)    :: is

    integer  :: i
    FLOAT :: max_eigen
    FLOAT, allocatable :: eigenval(:), occ(:)
#if defined(HAVE_MPI)
    integer  :: ierr
#endif

    ALLOCATE(eigenval(st%nst), st%nst)
    ALLOCATE     (occ(st%nst), st%nst)
    eigenval = M_ZERO; occ = M_ZERO

    do i = st%st_start, st%st_end
      eigenval(i) = st%eigenval(i, is)
      occ(i) = st%occ(i, is)
    end do

#if defined(HAVE_MPI)
    if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      call mpi_barrier(st%mpi_grp%comm, ierr)
      do i = 1, st%nst
        call mpi_bcast(eigenval(i), 1, R_MPITYPE, st%node(i), st%mpi_grp%comm, ierr)
        call mpi_bcast(occ(i), 1, R_MPITYPE, st%node(i), st%mpi_grp%comm, ierr)
      end do
    end if
#endif

    ! find out the top occupied state, to correct for the assymptotics
    ! of the potential
    max_eigen = CNST(-1e30)
    do i = 1, st%nst
      if((occ(i) .gt. small).and.(eigenval(i).gt.max_eigen)) then
        max_eigen = eigenval(i)
      end if
    end do

    oep%eigen_n = 1
    do i = 1, st%nst
      if(occ(i) .gt. small) then
        ! criterium for degeneracy
        if(abs(eigenval(i)-max_eigen).le.CNST(1e-3)) then
          oep%eigen_type(i) = 2
        else
          oep%eigen_type(i) = 1
          oep%eigen_index(oep%eigen_n) = i
          oep%eigen_n = oep%eigen_n + 1
        end if
      else
        oep%eigen_type(i) = 0
      end if
    end do
    oep%eigen_n = oep%eigen_n - 1

    deallocate(eigenval, occ)
  end subroutine xc_oep_AnalizeEigen


#include "undef.F90"
#include "real.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"
#include "xc_OEP_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_KLI.F90"
#include "xc_OEP_x.F90"
#include "xc_OEP_SIC.F90"
#include "xc_OEP_inc.F90"

#include "undef.F90"

end module xc_OEP_m
