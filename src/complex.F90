#define R_TCOMPLEX 1
#define R_FUNC(x) z ## x
#define R_TYPE complex(r8)
#define R_ABS(x) abs(x)
#define R_CONJ(x) conjg(x)
#define R_REAL(x) real(x, r8)
#define R_AIMAG(x) aimag(x, r8)
#define R_DOT zdotc
#define R_NRM2 dznrm2
#define REALORCOMPLEX(x) cmplx(x, 0.0_r8, r8)
