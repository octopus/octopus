#define R_TREAL 1
#define R_FUNC(x) d ## x
#define R_TYPE real(r8)
#define R_ABS(x) abs(x)
#define R_CONJ(x) (x)
#define R_REAL(x) (x)
#define R_AIMAG(x) (0._r8)
#define R_DOT ddot
#define R_NRM2 dnrm2
#define REALORCOMPLEX(x) real(x, r8)
