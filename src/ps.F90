#include "config_F90.h"

#if defined(ONE_D)
#  include "ps1D.F90"
#elif defined(THREE_D)
#  include "ps3D.F90"
#endif
