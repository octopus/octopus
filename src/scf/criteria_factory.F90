!! Copyright (C) 2020 N. Tancogne-Dejean
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

#include "global.h"

module criteria_factory_oct_m
  use convergence_criteria_oct_m
  use density_criteria_oct_m
  use eigenval_criteria_oct_m
  use energy_criteria_oct_m
  use global_oct_m
  use namespace_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                               &
             criteria_factory_init

contains
  
  ! ---------------------------------------------------------
  subroutine criteria_factory_init(list, namespace, max_iter, check_conv)
    class(criteria_list_t), intent(inout) :: list
    type(namespace_t),      intent(in)    :: namespace
    integer,                intent(in)    :: max_iter
    logical,                intent(out)   :: check_conv

    FLOAT :: conv_abs_dens, conv_rel_dens, conv_abs_ev, conv_rel_ev
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev
    FLOAT :: conv_energy_diff
    FLOAT :: energy_diff
    class(convergence_criteria_t), pointer    :: crit, other
    type(criteria_iterator_t) :: iter

    PUSH_SUB(criteria_factory_init)
 
    !%Variable ConvEnergy
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Stop the SCF when the magnitude of change in energy during at
    !% one SCF iteration is smaller than this value.
    !%
    !%A zero value (the default) means do not use this criterion.
    !%
    !% If this criterion is used, the SCF loop will only stop once it is
    !% fulfilled for two consecutive iterations.
    !%End
    call parse_variable(namespace, 'ConvEnergy', M_ZERO, conv_energy_diff, unit = units_inp%energy)
    crit => energy_criteria_t(conv_energy_diff, M_ZERO)
    call list%add(crit)

    !%Variable ConvAbsDens
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the density: 
    !%
    !% <math>\varepsilon = \int {\rm d}^3r \left| \rho^{out}(\bf r) -\rho^{inp}(\bf r) \right|</math>.
    !%
    !% A zero value (the default) means do not use this criterion.
    !%
    !% If this criterion is used, the SCF loop will only stop once it is
    !% fulfilled for two consecutive iterations.
    !%End
    call parse_variable(namespace, 'ConvAbsDens', M_ZERO, conv_abs_dens)

    !%Variable ConvRelDens
    !%Type float
    !%Default 1e-6
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the density: 
    !%
    !% <math>\varepsilon = \frac{1}{N} \mathrm{ConvAbsDens}</math>.
    !% 
    !% <i>N</i> is the total number of electrons in the problem.  A
    !% zero value means do not use this criterion.
    !%
    !% If you reduce this value, you should also reduce
    !% <tt>EigensolverTolerance</tt> to a value of roughly 1/10 of
    !% <tt>ConvRelDens</tt> to avoid convergence problems.
    !%
    !% If this criterion is used, the SCF loop will only stop once it is
    !% fulfilled for two consecutive iterations.
    !%End
    call parse_variable(namespace, 'ConvRelDens', CNST(1e-6), conv_rel_dens)
    crit => density_criteria_t(conv_abs_dens, conv_rel_dens)
    call list%add(crit)

    !%Variable ConvAbsEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the sum of the eigenvalues:
    !%
    !% <math> \varepsilon = \left| \sum_{j=1}^{N_{occ}} \varepsilon_j^{out} -
    !% \sum_{j=1}^{N_{occ}} \varepsilon_j^{inp} \right| </math>
    !%
    !% A zero value (the default) means do not use this criterion.
    !%
    !% If this criterion is used, the SCF loop will only stop once it is
    !% fulfilled for two consecutive iterations.
    !%End
    call parse_variable(namespace, 'ConvAbsEv', M_ZERO, conv_abs_ev, unit = units_inp%energy)
    call list%add(crit)

    !%Variable ConvRelEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the sum of the eigenvalues:
    !%
    !% <math>\varepsilon = \frac{ \left| \sum_{j=1}^{N_{occ}} ( \varepsilon_j^{out} -  \varepsilon_j^{inp} ) \right|}
    !% {\left| \sum_{j=1}^{N_{occ}} \varepsilon_j^{out} \right|} </math>
    !%
    !%A zero value (the default) means do not use this criterion.
    !%
    !% If this criterion is used, the SCF loop will only stop once it is
    !% fulfilled for two consecutive iterations.
    !%End
    call parse_variable(namespace, 'ConvRelEv', M_ZERO, conv_rel_ev)
    crit => eigenval_criteria_t(conv_abs_ev, conv_rel_ev)
    call list%add(crit)

    call messages_obsolete_variable(namespace, 'ConvForce')
    call messages_obsolete_variable(namespace, 'ConvAbsForce')
    call messages_obsolete_variable(namespace, 'ConvRelForce')

    call iter%start(list)
    check_conv = iter%has_next()  !If the list is empty, this fails  
    do while (iter%has_next())
      other => iter%get_next()
      check_conv = check_conv .or. (other%tol_abs > M_ZERO) .or. (other%tol_rel > M_ZERO)
    end do

    POP_SUB(criteria_factory_init)
  end subroutine criteria_factory_init

end module criteria_factory_oct_m
  
