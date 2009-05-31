#include <config.h>

#if defined(OCT_ARCH_IA64)
#include "operate_ia64.S"
#elif defined(OCT_ARCH_X86_64)
#include "operate_x86_64.S"
#else
#error "There is no assembler version for this routine"	
#endif
