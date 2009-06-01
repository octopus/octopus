#include <config.h>

#if defined(OCT_ARCH_IA64)
#include "operate_ia64.S"
#else
#error "There is no assembler version for this routine"	
#endif
