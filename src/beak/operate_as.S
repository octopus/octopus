#include <config.h>

#if OCT_ARCH == ia64
#include "operate_ia64.S"
#else
#error "There is no assembler version for this routine"	
#endif
