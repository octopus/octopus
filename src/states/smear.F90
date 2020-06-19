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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! Some pieces glued from PW/w(0,1)gauss.f90 from PWSCF

#include "global.h"

module smear_oct_m
  use global_oct_m
  use kpoints_oct_m
  use loct_math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use sort_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

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
    private
    integer, public :: method       !< which smearing function to take
    FLOAT,   public :: dsmear       !< the parameter defining this function
    FLOAT,   public :: e_fermi      !< the Fermi energy
    
    integer, public :: el_per_state !< How many electrons can we put in each state (1 or 2)
    FLOAT,   public :: ef_occ       !< Occupancy of the level at the Fermi energy
    logical, public :: integral_occs !< for fixed_occ, are they all integers?
    integer         :: MP_n         !< order of Methfessel-Paxton smearing
    integer         :: fermi_count  !< The number of occupied states at the fermi level
    integer         :: nik_factor   !< denominator, for treating k-weights as integers
    integer         :: nspins       !< = 2 if spin_polarized, else 1.
  end type smear_t

  integer, parameter, public ::       &
    SMEAR_SEMICONDUCTOR     = 1,      &
    SMEAR_FERMI_DIRAC       = 2,      &
    SMEAR_COLD              = 3,      &
    SMEAR_METHFESSEL_PAXTON = 4,      &
    SMEAR_SPLINE            = 5,      &
    SMEAR_FIXED_OCC         = 6
  
contains

  !--------------------------------------------------
  subroutine smear_init(this, namespace, ispin, fixed_occ, integral_occs, kpoints)
    type(smear_t),             intent(out) :: this
    type(namespace_t),         intent(in)  :: namespace
    integer,                   intent(in)  :: ispin
    logical,                   intent(in)  :: fixed_occ
    logical,                   intent(in)  :: integral_occs
    type(kpoints_t),           intent(in)  :: kpoints

    PUSH_SUB(smear_init)

    this%integral_occs = integral_occs

    !%Variable SmearingFunction
    !%Type integer
    !%Default semiconducting
    !%Section States
    !%Description
    !% This is the function used to smear the electronic occupations.
    !% It is ignored if the <tt>Occupations</tt> block is set.
    !%Option semiconducting 1
    !% Semiconducting occupations, <i>i.e.</i> the lowest lying states are occupied
    !% until no more electrons are left.
    !%Option fermi_dirac 2
    !% Simple Fermi-Dirac distribution. In this case, <tt>Smearing</tt> has
    !% the meaning of an electronic temperature. DN Mermin, <i>Phys. Rev.</i> <b>137</b>, A1441 (1965).
    !%Option cold_smearing 3
    !% N Marzari, D Vanderbilt, A De Vita, and MC Payne, <i>Phys. Rev. Lett.</i> <b>82</b>, 3296 (1999).
    !%Option methfessel_paxton 4
    !% M Methfessel and AT Paxton, <i>Phys. Rev. B</i> <b>40</b>, 3616 (1989).
    !% In this case, the variable <tt>SmearingMPOrder</tt> sets the order of the smearing.
    !% Occupations may be negative.
    !%Option spline_smearing 5
    !% Nearly identical to Gaussian smearing.
    !% JM Holender, MJ Gillan, MC Payne, and AD Simpson, <i>Phys. Rev. B</i> <b>52</b>, 967 (1995).
    !%End
    if(fixed_occ) then
      this%method = SMEAR_FIXED_OCC
    else
      call parse_variable(namespace, 'SmearingFunction', SMEAR_SEMICONDUCTOR, this%method)
      if(.not. varinfo_valid_option('SmearingFunction', this%method)) call messages_input_error(namespace, 'SmearingFunction')
      call messages_print_var_option(stdout, 'SmearingFunction', this%method)
    end if

    !%Variable Smearing
    !%Type float
    !%Default 0.1 eV
    !%Section States
    !%Description
    !% If <tt>Occupations</tt> is not set, <tt>Smearing</tt> is the
    !% smearing width used in the <tt>SmearingFunction</tt> to distribute the electrons
    !% among the existing states.
    !%End
    this%dsmear = CNST(1e-14)
    if(this%method /= SMEAR_SEMICONDUCTOR .and. this%method /= SMEAR_FIXED_OCC) then
      call parse_variable(namespace, 'Smearing', CNST(0.1) / (M_TWO * P_Ry), this%dsmear, units_inp%energy)
    end if

    call messages_obsolete_variable(namespace, 'ElectronicTemperature', 'Smearing')

    this%el_per_state = 1
    if(ispin == 1) & ! unpolarized
      this%el_per_state = 2

    if(ispin == 2) then
      this%nspins = 2
    else
      this%nspins = 1
    end if

    if(this%method == SMEAR_SEMICONDUCTOR) then
      this%nik_factor = kpoints_kweight_denominator(kpoints)

      if(this%nik_factor == 0) then
        message(1) = "k-point weights in KPoints or KPointsReduced blocks must be rational numbers for semiconducting smearing."
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    this%MP_n = 0
    if(this%method == SMEAR_METHFESSEL_PAXTON) then
      !%Variable SmearingMPOrder
      !%Type integer
      !%Default 1
      !%Section States
      !%Description
      !% Sets the order of the Methfessel-Paxton smearing function.
      !%End
      call parse_variable(namespace, 'SmearingMPOrder', 1, this%MP_n)
    end if

    POP_SUB(smear_init)
  end subroutine smear_init


  !--------------------------------------------------
  subroutine smear_copy(to, from)
    type(smear_t), intent(out) :: to
    type(smear_t), intent(in)  :: from

    PUSH_SUB(smear_copy)

    to%method       = from%method
    to%dsmear       = from%dsmear
    to%e_fermi      = from%e_fermi
    to%el_per_state = from%el_per_state
    to%ef_occ       = from%ef_occ
    to%MP_n         = from%MP_n
    to%fermi_count  = from%fermi_count
    to%nik_factor   = from%nik_factor

    POP_SUB(smear_copy)
  end subroutine smear_copy


  !--------------------------------------------------
  subroutine smear_find_fermi_energy(this, namespace, eigenvalues, occupations, &
    qtot, nik, nst, kweights)
    type(smear_t),     intent(inout) :: this
    type(namespace_t), intent(in)    :: namespace
    FLOAT,             intent(in)    :: eigenvalues(:,:), occupations(:,:)
    FLOAT,             intent(in)    :: qtot, kweights(:)
    integer,           intent(in)    :: nik, nst

    integer, parameter :: nitmax = 200
    FLOAT, parameter   :: tol = CNST(1.0e-10)
    integer            :: ist, ik, iter, maxq, weight, sumq_int
    FLOAT              :: drange, xx, emin, emax, sumq, dsmear, sumq_frac
    logical            :: conv
    FLOAT,   allocatable :: eigenval_list(:)
    integer, allocatable :: k_list(:), reorder(:)
    integer            :: fermi_count_up, fermi_count_down

    PUSH_SUB(smear_find_fermi_energy)

    maxq = this%el_per_state * nst * this%nspins
    if (maxq - qtot <= -tol) then ! not enough states
      message(1) = 'Not enough states'
      write(message(2),'(6x,a,f12.6,a,i10)')'(total charge = ', qtot, &
        ' max charge = ', maxq
      call messages_fatal(2, namespace=namespace)
    end if

    conv = .true.
    if(this%method == SMEAR_FIXED_OCC) then ! Fermi energy: last of occupied states
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
      ! first we sort the eigenvalues
      SAFE_ALLOCATE(eigenval_list(1:nst * nik))
      SAFE_ALLOCATE(       k_list(1:nst * nik))
      SAFE_ALLOCATE(      reorder(1:nst * nik))

      iter = 1
      do ist = 1, nst
        do ik = 1, nik
          eigenval_list(iter) = eigenvalues(ist, ik)
          k_list(iter) = ik
          reorder(iter) = iter
          iter = iter + 1
        end do
      end do
      
      call sort(eigenval_list, reorder)

      sumq_int = int(qtot) * this%nik_factor
      sumq_frac = qtot * this%nik_factor - sumq_int
      
      do iter = 1, nst * nik
        weight = int(kweights(k_list(reorder(iter))) * this%nik_factor + M_HALF)
        if(.not. weight > 0) cycle 
        this%e_fermi = eigenval_list(iter)
        this%ef_occ  = (sumq_int + sumq_frac) / (weight * this%el_per_state)

        if(sumq_int - weight * this%el_per_state <= 0) then
          ! count how many occupied states are at the fermi level,
          ! this is required later to fill the states
          this%fermi_count = 1
          fermi_count_up = 1
          fermi_count_down = 1
          do
            if(iter - fermi_count_down < 1) exit
            if(abs(this%e_fermi - eigenval_list(iter - fermi_count_down)) > CNST(1e-6)) exit
            fermi_count_down = fermi_count_down + 1
            this%ef_occ = this%ef_occ  &
                    + int(kweights(k_list(reorder(iter-fermi_count_down))) * this%nik_factor + M_HALF)
          end do
          do
            if(iter + fermi_count_up > nst*nik) exit
            if(abs(this%e_fermi - eigenval_list(iter + fermi_count_up)) > CNST(1e-6)) exit
            fermi_count_up = fermi_count_up + 1
          end do
          this%fermi_count = fermi_count_up + fermi_count_down - 1
          this%ef_occ = this%ef_occ / this%fermi_count
          exit
        end if

        sumq_int = sumq_int - weight * this%el_per_state
      end do
      ASSERT(this%ef_occ < M_ONE + M_EPSILON)

      SAFE_DEALLOCATE_A(eigenval_list)
      SAFE_DEALLOCATE_A(k_list)
      SAFE_DEALLOCATE_A(reorder)

    else ! bisection
      dsmear = max(CNST(1e-14), this%dsmear)
      drange = dsmear * sqrt(-log(tol * CNST(.01)))

      emin = minval(eigenvalues) - drange
      emax = maxval(eigenvalues) + drange

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
        message(1) = 'Fermi: did not converge.'
        call messages_fatal(1, namespace=namespace)
      end if

    end if

    POP_SUB(smear_find_fermi_energy)
  end subroutine smear_find_fermi_energy


  ! ---------------------------------------------------------
  subroutine smear_fill_occupations(this, eigenvalues, occupations, nik, nst)
    type(smear_t),   intent(in)    :: this
    FLOAT,           intent(in)    :: eigenvalues(:,:)
    FLOAT,           intent(inout) :: occupations(:,:)
    integer,         intent(in)    :: nik, nst

    integer :: ik, ist, ifermi
    FLOAT   :: dsmear, xx

    PUSH_SUB(smear_fill_occupations)

    if(this%method == SMEAR_FIXED_OCC) then
      ! do nothing
    else if(this%method == SMEAR_SEMICONDUCTOR) then
      ASSERT(this%fermi_count > 0 .and. this%fermi_count <= nik*nst)
      
      ifermi = 0
      do ik = 1, nik
        do ist = 1, nst
          xx = eigenvalues(ist, ik) - this%e_fermi
          if(xx < -CNST(1e-6)) then
            occupations(ist, ik) = this%el_per_state
          else if(abs(xx) <= CNST(1e-6) .and. ifermi < this%fermi_count) then
            occupations(ist, ik) = this%ef_occ * this%el_per_state
            ifermi = ifermi + 1
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

    POP_SUB(smear_fill_occupations)
  end subroutine smear_fill_occupations


  !--------------------------------------------------
  FLOAT function smear_calc_entropy(this, eigenvalues, &
    nik, nst, kweights, occ) result(entropy)
    type(smear_t),  intent(inout) :: this
    FLOAT,          intent(in)    :: eigenvalues(:,:)
    FLOAT,          intent(in)    :: kweights(:)
    integer,        intent(in)    :: nik, nst
    FLOAT,          intent(in)    :: occ(:, :) !< used if fixed_occ

    integer :: ist, ik
    FLOAT :: dsmear, xx, term, ff

    PUSH_SUB(smear_calc_entropy)

    entropy = M_ZERO

    if(this%method == SMEAR_FIXED_OCC .or. this%method == SMEAR_SEMICONDUCTOR) then
    ! Fermi-Dirac entropy, not quite the same as will be obtained with true smearing
    ! RM Wentzcovitch, JL Martins, and PB Allen, Phys. Rev. B 45, 11372 (1992) eqn (5)
    ! also N Marzari PhD thesis p 117, http://quasiamore.mit.edu/phd/Marzari_PhD.pdf
      do ik = 1, nik
        do ist = 1, nst
          ff = occ(ist, ik) / this%el_per_state
          if(ff > M_ZERO .and. ff  <  M_ONE) then
            term = ff * log(ff) + (1 - ff) * log (1 - ff)
          else ! we have semiconducting smearing, or perverse occupations as in Methfessel-Paxton
            term = M_ZERO      
          end if
          entropy = entropy - kweights(ik) * this%el_per_state * term
        end do
      end do
    else
      dsmear = max(CNST(1e-14), this%dsmear)
  
      do ik = 1, nik
        do ist = 1, nst
          if (eigenvalues(ist, ik) < HUGE(M_ONE)) then
            xx = (this%e_fermi - eigenvalues(ist, ik)) / dsmear
            entropy = entropy - kweights(ik) * this%el_per_state *  &
              smear_entropy_function(this, xx)
          end if
        end do
      end do
    end if

    POP_SUB(smear_calc_entropy)
  end function smear_calc_entropy


  ! ---------------------------------------------------------
  FLOAT function smear_delta_function(this, xx) result(deltaf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, aa
    integer :: ii, ni

    ! no PUSH_SUB, called too often

    ! smear_delta_function is not defined for SMEAR_FIXED_OCC
    ASSERT(this%method /= SMEAR_FIXED_OCC)
    
    deltaf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)
      if(abs(xx) <= M_EPSILON) &
        deltaf = this%ef_occ

    case(SMEAR_FERMI_DIRAC)
      if (abs(xx) <= CNST(36.0)) &
        deltaf = M_ONE / (M_TWO + exp(-xx) + exp(xx))

    case(SMEAR_COLD)
      xp  = xx - M_ONE / sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      deltaf = exp(-arg) / sqrt(M_PI) * (M_TWO - sqrt(M_TWO) * xx)
      
    case(SMEAR_METHFESSEL_PAXTON)
      arg    = min(maxarg, xx**2)
      deltaf = exp(-arg) / sqrt(M_PI)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        hp = exp(-arg)
        ni = 0
        aa = M_ONE / sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO * xx * hp - M_TWO * ni * hd
          ni = ni + 1
          aa = -aa / (M_FOUR * ii)
          hp = M_TWO * xx * hd - M_TWO * ni * hp
          ni = ni + 1
          deltaf = deltaf + aa * hp
        end do
      end if

    case(SMEAR_SPLINE)
      xp     = abs(xx) + M_ONE / sqrt(M_TWO)
      deltaf = sqrt(M_E) * xp * exp(-xp * xp)

    end select
    
  end function smear_delta_function


  ! ---------------------------------------------------------
  FLOAT function smear_step_function(this, xx) result(stepf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, aa
    integer :: ii, ni

    PUSH_SUB(smear_step_function)

    ! smear_step_function is not defined for SMEAR_FIXED_OCC
    ASSERT(this%method /= SMEAR_FIXED_OCC)

    stepf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)
      if(xx > M_ZERO) then
        stepf = M_ONE
      else if(abs(xx) <= M_EPSILON) then
        stepf = this%ef_occ
      end if

    case(SMEAR_FERMI_DIRAC)
      if (xx > maxarg) then
        stepf = M_ONE
      else if(xx > -maxarg) then
        stepf = M_ONE / (M_ONE + exp(-xx))
      end if

    case(SMEAR_COLD)
      xp  = xx - M_ONE / sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      stepf = M_HALF * loct_erf(xp) + &
        M_ONE / sqrt(M_TWO * M_PI) * exp(-arg) + M_HALF
      
    case(SMEAR_METHFESSEL_PAXTON)
      stepf = M_HALF * loct_erfc(-xx)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        arg = min(maxarg, xx**2)
        hp = exp(-arg)
        ni = 0
        aa = M_ONE / sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO * xx * hp - M_TWO * ni * hd
          ni = ni + 1
          aa = -aa / (M_FOUR * ii)
          stepf = stepf - aa * hd
          hp = M_TWO * xx * hd - M_TWO * ni * hp
          ni = ni + 1
        end do
      end if

    case(SMEAR_SPLINE)
      if(xx <= M_ZERO) then
        xp = xx - M_ONE / sqrt(M_TWO)
        stepf = M_HALF * sqrt(M_E) * exp(-xp * xp)
      else
        xp = xx + M_ONE / sqrt(M_TWO)
        stepf = M_ONE - M_HALF * sqrt(M_E) * exp(-xp * xp)
      end if

    end select

    POP_SUB(smear_step_function)
  end function smear_step_function


  ! ---------------------------------------------------------
  !> This function is defined as \int_{-infty}^x y delta(y) dy
  FLOAT function smear_entropy_function(this, xx) result(entropyf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) ::  xx

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: xp, arg, hd, hp, hpm1, aa
    integer :: ii, ni

    PUSH_SUB(smear_entropy_function)

    ! smear_entropy_function is not defined for SMEAR_FIXED_OCC
    ASSERT(this%method /= SMEAR_FIXED_OCC)

    entropyf = M_ZERO
    select case(this%method)
    case(SMEAR_SEMICONDUCTOR)

    case(SMEAR_FERMI_DIRAC)
      if(abs(xx) <= CNST(36.0)) then
        xp = M_ONE / (M_ONE + exp(-xx))
        entropyf = xp * log(xp) + (M_ONE - xp) * log(M_ONE - xp)
      end if

    case(SMEAR_COLD)
      xp  = xx - M_ONE / sqrt(M_TWO)
      arg = min(maxarg, xp**2)

      entropyf =  M_ONE / sqrt(M_TWO * M_PI) * xp * exp(-arg)
      
    case(SMEAR_METHFESSEL_PAXTON)
      arg = min(maxarg, xx**2)
      entropyf = -M_HALF * exp(-arg) / sqrt(M_PI)

      if(this%MP_n > 0) then ! recursion
        hd = M_ZERO
        hp = exp(-arg)
        ni = 0
        aa = M_ONE / sqrt(M_PI)
        do ii = 1, this%MP_n
          hd = M_TWO * xx * hp - M_TWO * ni * hd
          ni = ni + 1
          hpm1 = hp
          hp = M_TWO * xx * hd - M_TWO * ni * hp
          ni = ni + 1
          aa = -aa / (M_FOUR * ii)
          entropyf = entropyf - aa * (M_HALF * hp + hpm1 * ni)
        end do
      end if

    case(SMEAR_SPLINE)
      xp = abs(xx) + M_ONE / sqrt(M_TWO)
      entropyf = -sqrt(M_E) * (abs(xx) * exp(-xp * xp) / M_TWO + sqrt(M_PI) / M_FOUR * loct_erfc(xp))

    end select

    POP_SUB(smear_entropy_function)
  end function smear_entropy_function


  ! ---------------------------------------------------------
  logical function smear_is_semiconducting(this) result(answer)
    type(smear_t), intent(in) :: this

    PUSH_SUB(smear_is_semiconducting)

    answer = this%method  ==  SMEAR_SEMICONDUCTOR

    POP_SUB(smear_is_semiconducting)
  end function smear_is_semiconducting

end module smear_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
