#include "config.h"

#if defined(ONE_D)
#  include "mesh1D.F90"
#elif defined(THREE_D)
#  include "mesh3D.F90"
#endif
