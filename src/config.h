! src/config.h.  Generated automatically by configure.  
! Do we have FFTW
#define HAVE_FFTW 1

! Do we have GSL
#define HAVE_GSL 1

! Do we have LAPACK
#define HAVE_LAPACK 1
#define HAVE_BLAS 1

! the size of the pointers
#define POINTER_SIZE 4

! Dimensionality of the system
#define ONE_D 1
! #undef THREE_D 

! #undef BOUNDARIES_ZERO_DERIVATIVE 

! #undef COMPLEX_WFNS 

! what do you wish for dinner, dear?
#ifdef COMPLEX_WFNS
#include "complex.F90"
#else
#include "real.F90"
#endif
