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

module xc_OEP
  use global
  use messages
  use states
  use lib_oct_parser
  use lib_basic_alg
  use lib_adv_alg
  use lib_xc
  use xc
  use mesh
  use functions
  use mesh_function
  use poisson
  use hamiltonian
  use linear_response
  use grid

  implicit none

  private
  public ::                        &
     xc_oep_type,                  &
     xc_oep_init,                  &
     xc_oep_end,                   &
     xc_oep_write_info,            &
     dxc_oep_calc, zxc_oep_calc

  ! the OEP levels
  integer, public, parameter :: &
     XC_OEP_NONE   = 0, &
     XC_OEP_SLATER = 1, &
     XC_OEP_KLI    = 2, &
     XC_OEP_CEDA   = 3, &  ! not yet implemented
     XC_OEP_FULL   = 4     ! half implemented

  type xc_oep_type
    integer       :: level   ! 0 = no oep, 1 = Slater, 2 = KLI, 3 = CEDA, 4 = full OEP
    FLOAT         :: mixing  ! how much of the function S(r) to add to vxc in every iteration
    type(lr_type) :: lr      ! to solve the equation H psi = b

    integer          :: eigen_n
    integer, pointer :: eigen_type(:), eigen_index(:)
    FLOAT            :: socc, sfact
    FLOAT,   pointer :: vxc(:), uxc_bar(:)
    FLOAT,   pointer :: dlxc(:, :)
    CMPLX,   pointer :: zlxc(:, :)
  end type xc_oep_type

  FLOAT, parameter :: small     = CNST(1.0e-5)

contains


  ! -----------------------------------------------------------
  subroutine xc_oep_init(oep, family, m, d)
    type(xc_oep_type),     intent(out) :: oep
    integer,               intent(in)  :: family
    type(mesh_type),       intent(in)  :: m
    type(states_dim_type), intent(in)  :: d

    if(iand(family, XC_FAMILY_OEP).eq.0) then
      oep%level = XC_OEP_NONE
      return
    end if

#if defined(HAVE_MPI)
    if(oep%level == XC_OEP_FULL) then
      message(1) = "Full OEP is not allowed with the code parallelized on orbitals..."
      call write_fatal(1)
    endif
#endif

    call loct_parse_int(check_inp('OEP_level'), XC_OEP_KLI, oep%level)
    if(oep%level<0.or.oep%level>XC_OEP_FULL) then
      message(1) = "OEP_level can only take the values:"
      message(2) = "1 (Slater), 2 (KLI), 3 (CEDA), or 4 (full OEP)"
      call write_fatal(2)
    end if
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
    allocate(oep%vxc(m%np))

    ! when performing full OEP, we need to solve a linear equation
    if(oep%level == XC_OEP_FULL) call lr_init(oep%lr, "OEP")

    ! the linear equation has to be more converged if we are to attain the required precision
    !oep%lr%conv_abs_dens = oep%lr%conv_abs_dens / (oep%mixing)

  end subroutine xc_oep_init


  ! -----------------------------------------------------------
  subroutine xc_oep_end(oep)
    type(xc_oep_type), intent(inout) :: oep

    if(oep%level.ne.XC_OEP_NONE) then
      deallocate(oep%vxc); nullify(oep%vxc)

      if(oep%level == XC_OEP_FULL) call lr_dealloc(oep%lr)
    end if

  end subroutine xc_oep_end


  ! -----------------------------------------------------------
  subroutine xc_oep_write_info(oep, iunit)
    type(xc_oep_type), intent(in) :: oep
    integer,           intent(in) :: iunit

    if(oep%level.eq.XC_OEP_NONE) return

    write(iunit, '(2x,a)') 'The OEP equation will be handled at the level of:'
    select case(oep%level)
    case (XC_OEP_SLATER); write(iunit, '(a)') '    Slater approximation'
    case (XC_OEP_KLI);    write(iunit, '(a)') '    KLI approximation'
    case (XC_OEP_CEDA);   write(iunit, '(a)') '    CEDA approximation'
    case (XC_OEP_FULL);   write(iunit, '(a)') '    Full OEP'
    end select

  end subroutine xc_oep_write_info


  ! -----------------------------------------------------------
  ! A couple of auxiliary functions for oep
  ! -----------------------------------------------------------
  subroutine xc_oep_SpinFactor(oep, nspin)
    type(xc_oep_type), intent(inout) :: oep
    integer,           intent(in)    :: nspin

    select case(nspin)
    case(1) ! we need to correct or the spin occupancies
      oep%socc  = M_HALF
      oep%sfact = M_TWO
    case(2, 4)
      oep%socc  = M_ONE
      oep%sfact = M_ONE
    case default
      write(6,'(a,I2)') 'OEP: error cannot handle nspin=', nspin
    end select

  end subroutine xc_oep_SpinFactor


  ! -----------------------------------------------------------
  subroutine xc_oep_AnalizeEigen(oep, st, is)
    type(xc_oep_type), intent(inout) :: oep
    type(states_type), intent(in)    :: st
    integer,           intent(in)    :: is

    integer  :: i
    FLOAT :: max_eigen
    FLOAT, allocatable :: eigenval(:), occ(:)
#if defined(HAVE_MPI)
    integer  :: ierr
#endif

    allocate(eigenval(st%nst), occ(st%nst))
    eigenval = M_ZERO; occ = M_ZERO

    do i = st%st_start, st%st_end
       eigenval(i) = st%eigenval(i, is)
       occ(i) = st%occ(i, is)
    enddo

#if defined(HAVE_MPI)
    if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      call mpi_barrier(mpi_comm_world, ierr)
      do i = 1, st%nst
         call mpi_bcast(eigenval(i), 1, R_MPITYPE, st%node(i), MPI_COMM_WORLD, ierr)
         call mpi_bcast(occ(i), 1, R_MPITYPE, st%node(i), MPI_COMM_WORLD, ierr)
      enddo
    endif
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

end module xc_OEP
