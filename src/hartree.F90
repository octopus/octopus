#include "config_F90.h"

#if defined(ONE_D)
#  include "hartree1D.F90"
#elif defined(THREE_D)
#  include "hartree3D.F90"
#endif
