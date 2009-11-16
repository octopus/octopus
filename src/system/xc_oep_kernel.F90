
#include "global.h"

module xc_OEP_kernel_m
  use global_m
  use hamiltonian_m
  use linear_response_m
  use messages_m
  use poisson_m
  use profiling_m
  use states_m
  use system_m
  use xc_OEP_m

  implicit none

  private

  public :: &
       xc_oep_kernel_init, & 
       dxc_oep_kernel_calc, &
       zxc_oep_kernel_calc

contains


subroutine xc_oep_kernel_init(oep)
  type(xc_oep_t),     intent(in)   :: oep

  !%Variable OEP_kernel_level
  !%Type integer
  !%Default oep_kli
  !%Section Hamiltonian::XC
  !%Description
  !% At what level shall <tt>Octopus</tt> handle the OEP kernel.
  !%Option oep_kli 3
  !% Krieger-Li-Iafrate (KLI) approximation
  !% (JB Krieger, Y Li, GJ Iafrate, <i>Phys. Rev. Lett. A</i> <b>146</b>, 256 (1990).
  !%End

end subroutine xc_oep_kernel_init

#include "undef.F90"
#include "real.F90"
#include "xc_oep_kernel_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "xc_oep_kernel_inc.F90"

end module xc_oep_kernel_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
