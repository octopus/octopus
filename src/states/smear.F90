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

! Some pieces glued from PW/w(0,1)gauss.f90 from PWSCF

#include "global.h"

module smear_m
  use datasets_m
  use loct_math_m
  use parser_m
  use global_m
  use math_m
  use messages_m
  use profiling_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none
  
  private
  public ::                           &
    smear_t,                          &
    smear_init,                       &
    smear_copy,                       &
    smear_find_fermi_energy,          &
    smear_fill_occupations,           &
    smear_calc_entropy,               &
    smear_delta_function,             &
    smear_step_function,              &
    smear_entropy_function,           &
    smear_is_semiconducting

  type smear_t
    integer :: method       ! which smearing function to take
    FLOAT   :: dsmear       ! the parameter defining this function
    FLOAT   :: e_fermi      ! the Fermi energy
    
    FLOAT   :: el_per_state ! How many electrons can we put in each state
    logical :: fixed_occ    ! Are occupations fixed, or are we allowed to change them
    FLOAT   :: ef_occ       ! Occupancy of the level at the Fermi energy

    integer :: MP_n         ! order of Methfessel-Paxton smearing
  end type smear_t

  integer, parameter, public ::       &
    SMEAR_SEMICONDUCTOR     = 1,      &
    SMEAR_FERMI_DIRAC       = 2,      &
    SMEAR_COLD              = 3,      &
    SMEAR_METHFESSEL_PAXTON = 4,      &
    SMEAR_SPLINE            = 5
  
contains

  !--------------------------------------------------
  subroutine smear_init(this, ispin, fixed_occ)
    type(smear_t), intent(out) :: this
    integer,       intent(in)  :: ispin
    logical,       intent(in)  :: fixed_occ

    !%Variable SmearingFunction
    !%Type integer
    !%Default semiconducting
    !%Section States
    !%Description
    !% This is the function used to smear the electronic occupations
    !%Option semiconducting 1
    !% Semiconducting occupations, i.e., the lowest lying states are occupied
    !% until no more electrons are left 
    !%Option fermi_dirac 2
    !% Simple Fermi-Dirac distribution. In this case, the <tt>Smearing</tt> has
    !% the meaning of an electronic temperature
    !%Option cold_smearing 3
    !% N Marzari, D Vanderbilt, A De Vita, and MC Payne, Phys. Rev. Lett. 82, 3296 (1999).
    !%Option methfessel_paxton 4
    !% M Methfessel and AT Paxton, Phys. Rev. B 40, 3616 (1989).
    !% In this case, the variable <tt>SmearingMPOrder</tt> sets the order of the smearing.
    !%Option spline_smearing 5
    !% Nearly identical to Gaussian smearing
    !%End
    call parse_integer(datasets_check('SmearingFunction'), SMEAR_SEMICONDUCTOR, this%method)
    if(.not.varinfo_valid_option('SmearingFunction', this%method)) call input_error('SmearingFunction')
    call messages_print_var_option(stdout, 'SmearingFunction', this%method)

    !%Variable Smearing
    !%Type float
    !%Default 0.1 eV
    !%Section States
    !%Description
    !% If <tt>Occupations</tt> is not set, <tt>Smearing</tt> is the
    !% smearing used in the <tt>SmearingFunction</tt> used to distribute the electrons
    !% among the existing states.
    !%End
    this%dsmear = CNST(1e-14)
    if(this%method.ne.SMEAR_SEMICONDUCTOR) then
      call parse_float(datasets_check('Smearing'), CNST(0.1)/(M_TWO*P_Ry), this%dsmear)
      this%dsmear = units_to_atomic(units_inp%energy, this%dsmear)
    end if

    call messages_obsolete_variable("ElectronicTemperature", "Smearing")

    this%fixed_occ = fixed_occ
    this%el_per_state = M_ONE
    if(ispin == 1) & ! unpolarized
      this%el_per_state = M_TWO

    this%MP_n = 0
    if(this%method == SMEAR_METHFESSEL_PAXTON) then
      !%Variable SmearingMPOrder
      !%Type integer
      !%Default 1
      !%Section States
      !%Description
      !% Sets the order of the Methfessel-Paxton smearing function
      !%End
      call parse_integer(datasets_check('SmearingMPOrder'), 1, this%MP_n)
    end if

  end subroutine smear_init


  !--------------------------------------------------
  subroutine smear_copy(to, from)
    type(smear_t), intent(out) :: to
    type(smear_t), intent(in)  :: from

    to%method       = from%method
    to%dsmear       = from%dsmear
    to%e_fermi      = from%e_fermi
    to%el_per_state = from%el_per_state
    to%fixed_occ    = from%fixed_occ
    to%ef_occ       = from%ef_occ
    to%MP_n         = from%MP_n
  end subroutine smear_copy


  !--------------------------------------------------
  subroutine smear_find_fermi_energy(this, eigenvalues, occupations, &
    qtot, nik, nst, kweights)
    type(smear_t), intent(inout) :: this
    FLOAT,         intent(in)    :: eigenvalues(:,:), occupations(:,:)
    FLOAT,         intent(in)    :: qtot, kweights(:)
    integer,       intent(in)    :: nik, nst

    integer, parameter :: nitmax = 200
    FLOAT, parameter   :: tol = CNST(1.0e-10)

    ! Local variables.
    integer            :: ist, ik, iter
    FLOAT              :: drange, xx, emin, emax, sumq, dsmear
    logical            :: conv

    FLOAT,   allocatable :: eigenval_list(:)
    integer, allocatable :: k_list(:), reorder(:)

    call push_sub('smear.smear_find_fermi_energy')

    ! Initializations
    emin = minval(eigenvalues)
    emax = maxval(eigenvalues)

    sumq = this%el_per_state * nst* sum(kweights(:))

    if (sumq - qtot <= -CNST(1e-10)) then ! not enough states
      message(1) = 'Not enough states'
      write(message(2),'(6x,a,f12.6,a,f12.6)')'(total charge = ', qtot, &
        ' max charge = ', sumq
      call write_fatal(2)
    end if

    conv = .true.
    if(this%fixed_occ) then ! Fermi energy: last of occupied states
      ist_cycle: do ist = nst, 1, -1
        do ik = 1, nik
          if(occupations(ist, ik) > CNST(1e-5)) then
            this%e_fermi  = eigenvalues(ist, ik)
            this%ef_occ   = occupations(ist, ik) / this%el_per_state
            exit ist_cycle
          end if

        end do
      end do ist_cycle

    else if(this%method == SMEAR_SEMICONDUCTOR) then
      sumq = qtot
      ! first we sort the eigenvalues
      SAFE_ALLOCATE(eigenval_list(1:nst*nik))
      SAFE_ALLOCATE(       k_list(1:nst*nik))
      SAFE_ALLOCATE(      reorder(1:nst*nik))

      iter = 1
      do ist = 1, nst
        do ik = 1, nik
          eigenval_list(iter) = eigenvalues(ist, ik)
          k_list(iter) = ik
          iter = iter + 1
        end do
      end do
      
      call sort(eigenval_list, reorder)

      do iter = 1, nst * nik
        xx = kweights(k_list(reorder(iter)))

        this%e_fermi = eigenval_list(iter)
        this%ef_occ  = sumq / (xx * this%el_per_state)

        if(sumq - xx * this%el_per_state <= CNST(1e-10)) then
          exit
        end if

        sumq = sumq - xx * this%el_per_state
      end do

      SAFE_DEALLOCATE_A(eigenval_list)
      SAFE_DEALLOCATE_A(k_list)
      SAFE_DEALLOCATE_A(reorder)

    else ! bisection
      dsmear = max(CNST(1e-14), this%dsmear)
      drange = dsmear * sqrt(-log(tol * CNST(.01)))

      emin = emin - drange
      emax = emax + drange

      do iter = 1, nitmax
        this%e_fermi = M_HALF*(emin + emax)
        sumq         = M_ZERO

        do ik = 1, nik
          do ist = 1, nst
            xx   = (this%e_fermi - eigenvalues(ist, ik))/dsmear
            sumq = sumq + kweights(ik) * this%el_per_state * &
              smear_step_function(this, xx)
          end do
        end do

        conv = (abs(sumq - qtot) <= tol)
        if(conv) exit

        if(sumq <= qtot ) emin = this%e_fermi
        if(sumq >= qtot ) emax = this%e_fermi

        this%ef_occ = smear_step_function(this, M_ZERO)
      end do

      if(.not.conv) then
        message(1) = 'Fermi: did not converge'
        call write_fatal(1)
      end if

    end if

    call pop_sub()
  end subroutine smear_find_fermi_energy


  ! ---------------------------------------------------------
  subroutine smear_fill_occupations(this, eigenvalues, occupations, nik, nst)
    type(smear_t), intent(in)  :: this
    FLOAT,         intent(in)  :: eigenvalues(:,:)
    FLOAT,         intent(out) :: occupations(:,:)
    integer,       intent(in)  :: nik, nst

    integer :: ik, ist
    FLOAT   :: dsmear, xx

    if(this%fixed_occ) then
      ! do nothing
    else if(this%method == SMEAR_SEMICONDUCTOR) then
      do ik = 1, nik
        do ist = 1, nst
          xx = eigenvalues(ist, ik) - this%e_fermi
          if(xx < M_ZERO) then
            occupations(ist, ik) = this%el_per_state
          else if(xx == M_ZERO) then
            occupations(ist, ik) = this%ef_occ * this%el_per_state
          else
            occupations(ist, ik) = M_ZERO
          end if
        end do
      end do

    else 
      dsmear = max(CNST(1e-14), this%dsmear)
      do ik = 1, nik
        do ist = 1, nst
          xx = (this%e_fermi - eigenvalues(ist, ik))/dsmear
          occupations(ist, ik) = smear_step_function(this, xx) * this%el_per_state
        end do
      end do
    end if

  end subroutine smear_fill_occupations


  !--------------------------------------------------
  FLOAT function smear_calc_entropy(this, eigenvalues, &
    nik, nst, kweights) result(entropy)
    type(smear_t), intent(inout) :: this
    FLOAT,         intent(in)    :: eigenvalues(:,:)
    FLOAT,         intent(in)    :: kweights(:)
    integer,       intent(in)    :: nik, nst

    integer :: ist, ik
    FLOAT :: dsmear, xx

    dsmear = max(CNST(1e-14), this%dsmear)
    entropy = M_ZERO
    do ik = 1, nik
      do ist = 1, nst
        xx = (this%e_fermi - eigenvalues(ist, ik))/dsmear
        entropy = entropy + kweights(ik) * this%el_per_state *  &
          dsmear*smear_entropy_function(this, xx)
      end do
    end do

  end function smear_calc_entropy


  ! ---------------------------------------------------------
  FLOAT function smear_delta_function(this, xx) result(deltaf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, A
    integer :: ii, ni

    deltaf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)
      if(xx == M_ZERO) &
        deltaf = this%ef_occ

    case(SMEAR_FERMI_DIRAC)
      if (abs(xx) <= CNST(36.0)) &
        deltaf = M_ONE/(M_TWO + exp(-xx) + exp(xx))

    case(SMEAR_COLD)
      xp  = xx - M_ONE/sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      deltaf = exp(-arg)/sqrt(M_PI)*(M_TWO - sqrt(M_TWO)*xx)
      
    case(SMEAR_METHFESSEL_PAXTON)
      arg    = min(maxarg, xx**2)
      deltaf = exp(-arg)/sqrt(M_PI)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        hp = exp(-arg)
        ni = 0
        A = M_ONE/sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO*xx*hp - M_TWO*ni*hd
          ni = ni + 1
          A = -A/(M_FOUR*ii)
          hp = M_TWO*xx*hd - M_TWO*ni*hp
          ni = ni + 1
          deltaf = deltaf + A*hp
        end do
      end if

    case(SMEAR_SPLINE)
      xp     = abs(xx) + M_ONE/sqrt(M_TWO)
      deltaf = sqrt(M_E)*xp*exp(-xp*xp)

    end select
    
  end function smear_delta_function


  ! ---------------------------------------------------------
  FLOAT function smear_step_function(this, xx) result(stepf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, A
    integer :: ii, ni

    stepf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)
      if(xx > M_ZERO) then
        stepf = M_ONE
      else if(xx == M_ZERO) then
        stepf = this%ef_occ
      end if

    case(SMEAR_FERMI_DIRAC)
      if (xx > maxarg) then
        stepf = M_ONE
      else if(xx > -maxarg) then
        stepf = M_ONE/(M_ONE + exp(-xx))
      end if

    case(SMEAR_COLD)
      xp  = xx - M_ONE/sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      stepf = M_HALF*loct_erf(xp) + &
        M_ONE/sqrt(M_TWO*M_PI)*exp(-arg) + M_HALF
      
    case(SMEAR_METHFESSEL_PAXTON)
      stepf = M_HALF*loct_erfc(-xx)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        arg = min(maxarg, xx**2)
        hp = exp(-arg)
        ni = 0
        A = M_ONE/sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO*xx*hp - M_TWO*ni*hd
          ni = ni + 1
          A = -A/(M_FOUR*ii)
          stepf = stepf - A*hd
          hp = M_TWO*xx*hd - M_TWO*ni*hp
          ni = ni + 1
        end do
      end if

    case(SMEAR_SPLINE)
      if(xx <= M_ZERO) then
        xp = xx - M_ONE/sqrt(M_TWO)
        stepf = M_HALF*sqrt(M_E)*exp(-xp*xp)
      else
        xp = xx + M_ONE/sqrt(M_TWO)
        stepf = M_ONE - M_HALF*sqrt(M_E)*exp(-xp*xp)
      end if

    end select

  end function smear_step_function


  ! ---------------------------------------------------------
  ! This function is defined as \int_{-infty}^x y delta(y) dy
  FLOAT function smear_entropy_function(this, xx) result(entropyf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, hpm1, A
    integer :: ii, ni

    entropyf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)

    case(SMEAR_FERMI_DIRAC)
      if(abs(xx) <= 36.0) then
        xp = M_ONE/(M_ONE + exp(-xx))
        entropyf = xp*log(xp) + (M_ONE - xp)*log(M_ONE - xp)
      end if

    case(SMEAR_COLD)
      xp  = xx - M_ONE/sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      entropyf =  M_ONE/sqrt(M_TWO*M_PI)*xp*exp(-arg)
      
    case(SMEAR_METHFESSEL_PAXTON)
      arg = min(maxarg, xx**2)
      entropyf = -M_HALF*exp(-arg)/sqrt(M_PI)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        hp = exp(-arg)
        ni = 0
        A = M_ONE/sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO*xx*hp - M_TWO*ni*hd
          ni = ni + 1
          hpm1 = hp
          hp = M_TWO*xx*hd - M_TWO*ni*hp
          ni = ni + 1
          A = -A/(M_FOUR*ii)
          entropyf = entropyf - A*(M_HALF*hp + hpm1*ni)
        end do
      end if

    case(SMEAR_SPLINE)
      xp = abs(xx) + M_ONE/sqrt(M_TWO)
      entropyf = -sqrt(M_E)*(abs(xx)*exp(-xp*xp)/M_TWO + sqrt(M_PI)/M_FOUR*loct_erfc(xp))

    end select

  end function smear_entropy_function

  logical function smear_is_semiconducting(this) result(answer)
    type(smear_t), intent(in) :: this

    answer = this%method .eq. SMEAR_SEMICONDUCTOR

  end function smear_is_semiconducting

end module smear_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
