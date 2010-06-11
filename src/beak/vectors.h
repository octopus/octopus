/*
 Copyright (C) 2010 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: vectors.h 2146 2006-05-23 17:36:00Z xavier $
*/

#ifdef VEC_SSE2
#include <emmintrin.h>
#define VEC_SIZE 2
#define VEC_TYPE __m128d
#define VEC_LD(addr) _mm_load_pd(addr)
#define VEC_LDU(addr) _mm_loadu_pd(addr)
#define VEC_ST(addr, vec)  _mm_store_pd(addr, vec)
#define VEC_STU(addr, vec)  _mm_storeu_pd(addr, vec)
#define VEC_FMA(aa, bb, cc) _mm_add_pd(cc, _mm_mul_pd(aa, bb))
#define VEC_SCAL(aa) _mm_set1_pd(aa)
#define VEC_ZERO _mm_setzero_pd()
#endif

#ifdef VEC_BG
#define VEC_SIZE 2
#define VEC_TYPE double _Complex
#define VEC_LD(addr) __lfpd(addr)
#define VEC_LDU(addr) __cmplx((addr)[0], (addr)[1])
#define VEC_ST(addr, vec)  __stfpd(addr, vec)
#define VEC_STU(addr, vec)  (addr)[0] = __creal(vec); (addr)[1] = __cimag(vec)
#define VEC_FMA(aa, bb, cc) __fpmadd(cc, aa, bb)
#define VEC_SCAL(aa) __cmplx(aa, aa)
#define VEC_ZERO __cmplx(0.0, 0.0)
#endif
