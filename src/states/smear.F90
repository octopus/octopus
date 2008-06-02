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
!! $Id: states.F90 4223 2008-05-31 13:28:56Z acastro $

#include "global.h"

module smear_m
  use datasets_m
  use loct_parser_m
  use global_m
  use math_m
  use messages_m
  use units_m
  use varinfo_m

  implicit none
  
  private
  public ::                           &
    smear_t,                          &
    smear_init,                       &
    smear_copy,                       &
    smear_find_fermi_energy,          &
    smear_fill_occupations

  type smear_t
    integer :: method      ! which smearing function to take
    FLOAT   :: dsmear      ! the parameter defining this function
    FLOAT   :: e_fermi     ! the Fermi energy
    
    logical :: fixed_occ   ! Are occupations fixed, or are we allowed to change them
    FLOAT   :: ef_occ      ! Occupancy of the level at the Fermi energy
  end type smear_t

  integer, parameter ::               &
    SMEAR_SEMICONDUCTOR = 1,          &
    SMEAR_FERMI_DIRAC   = 2
  
contains

  !--------------------------------------------------
  subroutine smear_init(this, fixed_occ)
    type(smear_t), intent(out) :: this
    logical,       intent(in)  :: fixed_occ

    !%Variable SmearingFunction
    !%Type integer
    !%Default fermi_dirac
    !%Section States
    !%Description
    !% This is the function used to smear the electronic occupations
    !%Option fermi_dirac 1
    !% Simple Fermi-Dirac distribution. In this case, the <tt>Smearing</tt> has
    !% the meaning of an electronic temperature
    !%End
    call loct_parse_int(check_inp('SmearingFunction'), SMEAR_SEMICONDUCTOR, this%method)
    if(.not.varinfo_valid_option('SmearingFunction', this%method)) call input_error('SmearingFunction')
    call messages_print_var_option(stdout, 'SmearingFunction', this%method)

    !%Variable Smearing
    !%Type float
    !%Default 0.0
    !%Section States
    !%Description
    !% If <tt>Occupations</tt> is not set, <tt>Smearing</tt> is the
    !% smearing used in the <tt>SmearingFunction</tt> used to distribute the electrons
    !% among the existing states.
    !%End
    this%dsmear = M_ZERO
    if(this%method.ne.SMEAR_SEMICONDUCTOR) then
      call loct_parse_float(check_inp('Smearing'), CNST(0.1)/(M_TWO*P_Ry), this%dsmear)
      this%dsmear = this%dsmear * units_inp%energy%factor
    end if

    this%fixed_occ = fixed_occ
  end subroutine smear_init


  !--------------------------------------------------
  subroutine smear_copy(to, from)
    type(smear_t), intent(out) :: to
    type(smear_t), intent(in)  :: from

    to%method    = from%method
    to%dsmear    = from%dsmear
    to%e_fermi   = from%e_fermi
    to%fixed_occ = from%fixed_occ
    to%ef_occ    = from%ef_occ
  end subroutine smear_copy


  !--------------------------------------------------
  subroutine smear_find_fermi_energy(this, eigenvalues, occupations, &
    qtot, nik, nst, spin_channels, is_spinors, kweights)
    type(smear_t), intent(inout) :: this
    FLOAT,         intent(in)    :: eigenvalues(:,:), occupations(:,:)
    FLOAT,         intent(in)    :: qtot, kweights(:)
    integer,       intent(in)    :: nik, nst, spin_channels
    logical,       intent(in)    :: is_spinors

    integer, parameter :: nitmax = 200
    FLOAT, parameter   :: tol = CNST(1.0e-10)

    ! Local variables.
    integer            :: ist, ik, iter
    FLOAT              :: drange, xx, emin, emax, sumq, dsmear, el_per_state
    logical            :: conv

    FLOAT,   allocatable :: eigenval_list(:)
    integer, allocatable :: k_list(:), reorder(:)

    call push_sub('smear.smear_find_fermi_energy')

    ! Initializations
    emin = minval(eigenvalues)
    emax = maxval(eigenvalues)

    el_per_state = M_TWO
    if(is_spinors) el_per_state = M_ONE

    sumq = el_per_state*nst
    if (sumq < qtot) then ! not enough states
      message(1) = 'Not enough states'
      write(message(2),'(6x,a,f12.6,a,f12.6)')'(total charge = ', qtot, &
        ' max charge = ', sumq
      call write_fatal(2)
    end if

    conv = .true.
    if(this%fixed_occ) then ! Fermi energy last of occupied states
      ist_cycle: do ist = nst, 1, -1
        do ik = 1, nik
          if(occupations(ist, ik) > CNST(1e-5)) then
            this%e_fermi  = eigenvalues(ist, ik)
            this%ef_occ   = occupations(ist, ik)
            exit ist_cycle
          end if

        end do
      end do ist_cycle

    else if(this%method == SMEAR_SEMICONDUCTOR) then
      sumq = qtot

      ! first we sort the eigenvalues
      ALLOCATE(eigenval_list(nst*nik), nst*nik)
      ALLOCATE(k_list(nst*nik), nst*nik)
      ALLOCATE(reorder(nst*nik), nst*nik)

      iter = 1
      do ist = 1, nst
        do ik = 1, nik
          eigenval_list(iter) = eigenvalues(ist, ik)
          k_list(iter) = ik
          iter = iter + 1
        end do
      end do
      
      call sort(eigenval_list, reorder)

      do iter = 1, nst*nik
        xx = kweights(k_list(reorder(iter)))

        if(sumq <= xx*el_per_state/spin_channels) then
          this%e_fermi = eigenval_list(iter)
          this%ef_occ  = sumq / xx
          exit
        end if

        sumq = sumq - xx*el_per_state/spin_channels
      end do

      deallocate(eigenval_list)
      deallocate(k_list)
      deallocate(reorder)

    else ! bisection
      dsmear = max(CNST(1e-14), this%dsmear)
      drange = dsmear*sqrt(-log(tol*CNST(.01)))

      emin = emin - drange
      emax = emax + drange

      do iter = 1, nitmax
        this%e_fermi = M_HALF*(emin + emax)
        sumq         = M_ZERO

        do ik = 1, nik
          do ist = 1, nst
            xx   = (eigenvalues(ist, ik) - this%e_fermi)/dsmear
            sumq = sumq + kweights(ik)*(M_TWO/spin_channels) * &
              smear_step_function(this, xx)
          end do
        end do

        conv = (abs(sumq - qtot) <= tol)
        if(conv) exit

        if(sumq <= qtot ) emin = this%e_fermi
        if(sumq >= qtot ) emax = this%e_fermi

        this%ef_occ = (M_TWO/spin_channels) * smear_step_function(this, M_ZERO)
      end do

      if(.not.conv) then
        message(1) = 'Fermi: did not converge'
        call write_fatal(1)
      end if

    end if

    call pop_sub()
  end subroutine smear_find_fermi_energy


  ! ---------------------------------------------------------
  subroutine smear_fill_occupations(this, eigenvalues, occupations, nik, nst, spin_channels)
    type(smear_t), intent(in)  :: this
    FLOAT,         intent(in)  :: eigenvalues(:,:)
    FLOAT,         intent(out) :: occupations(:,:)
    integer,       intent(in)  :: nik, nst, spin_channels

    integer :: ik, ist
    FLOAT   :: dsmear, xx

    if(this%fixed_occ) then
      ! do nothing
    else if(this%method == SMEAR_SEMICONDUCTOR) then
      do ik = 1, nik
        do ist = 1, nst
          xx = eigenvalues(ist, ik) - this%e_fermi
          if(xx < M_ZERO) then
            occupations(ist, ik) = (M_TWO/spin_channels)
          else if(xx == M_ZERO) then
            occupations(ist, ik) = this%ef_occ
          else
            occupations(ist, ik) = M_ZERO
          end if
        end do
      end do

    else 
      dsmear = max(CNST(1e-14), this%dsmear)
      do ik = 1, nik
        do ist = 1, nst
          xx = (eigenvalues(ist, ik) - this%e_fermi)/dsmear
          occupations(ist, ik) = (M_TWO/spin_channels) * smear_step_function(this, xx)
        end do
      end do
    end if

  end subroutine smear_fill_occupations


  ! ---------------------------------------------------------
  FLOAT function smear_step_function(this, xx) result(stepf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    stepf = M_ZERO

    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)
      if(xx < M_ZERO) then
        stepf = M_ONE
      else if(xx == M_ZERO) then
        stepf = this%ef_occ
      end if

    case(SMEAR_FERMI_DIRAC)
      if (xx < CNST(-100.0)) then
        stepf = M_ONE
      else if(xx < CNST(100.0)) then
        stepf = M_ONE / ( M_ONE + exp(xx) )
      end if
    end select

  end function smear_step_function


end module smear_m
