#ifndef IN
#define IN in
#endif

!#define COMPLEX_WFNS

#define HAVE_FFTW
#define HAVE_LAPACK

#ifdef COMPLEX_WFNS
#include "complex.F90"
#else
#include "real.F90"
#endif
