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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#ifndef VECTORS_H
#define VECTORS_H

#ifdef HAVE_VEC

#ifdef HAVE_M512D
#include <immintrin.h>
#define VEC_SIZE 8
#define VEC_TYPE __m512d
#define VEC_LD(addr) _mm512_load_pd(addr)
#define VEC_LDU(addr) _mm512_loadu_pd(addr)
#define VEC_ST(addr, vec)  _mm512_stream_pd(addr, vec)
#define VEC_STU(addr, vec)  _mm512_storeu_pd(addr, vec)
#define VEC_FMA(aa, bb, cc) _mm512_fmadd_pd(aa, bb, cc)
#define VEC_SCAL(aa) _mm512_set1_pd(aa)
#define VEC_ZERO _mm512_setzero_pd()
#include <emmintrin.h>
#define FENCE _mm_mfence()

#define DEPTH 16

#elif defined(HAVE_M256D)
#include <immintrin.h>
#if defined(HAVE_FMA4) || defined(HAVE_FMA3)
#include <x86intrin.h>
#endif
#define VEC_SIZE 4
#define VEC_TYPE __m256d
#define VEC_LD(addr) _mm256_load_pd(addr)
#define VEC_LDU(addr) _mm256_loadu_pd(addr)
#define VEC_ST(addr, vec)  _mm256_stream_pd(addr, vec)
#define VEC_STU(addr, vec)  _mm256_storeu_pd(addr, vec)
#ifdef HAVE_FMA3
#define VEC_FMA(aa, bb, cc) _mm256_fmadd_pd(aa, bb, cc)
#elif defined(HAVE_FMA4)
#define VEC_FMA(aa, bb, cc) _mm256_macc_pd(aa, bb, cc)
#else
#define VEC_FMA(aa, bb, cc) _mm256_add_pd(cc, _mm256_mul_pd(aa, bb))
#endif
#define VEC_SCAL(aa) _mm256_set1_pd(aa)
#define VEC_ZERO _mm256_setzero_pd()
#include <emmintrin.h>
#define FENCE _mm_mfence()

#define DEPTH 16
#endif

#if !defined(HAVE_M256D) && defined(HAVE_M128D)
#include <emmintrin.h>
#if defined(HAVE_FMA4) || defined(HAVE_FMA3)
#include <x86intrin.h>
#endif
#define VEC_SIZE 2
#define VEC_TYPE __m128d
#define VEC_LD(addr) _mm_load_pd(addr)
#define VEC_LDU(addr) _mm_loadu_pd(addr)
#define VEC_ST(addr, vec)  _mm_stream_pd(addr, vec)
#define VEC_STU(addr, vec)  _mm_storeu_pd(addr, vec)
#ifdef HAVE_FMA3
#define VEC_FMA(aa, bb, cc) _mm_fmadd_pd(aa, bb, cc)
#elif defined(HAVE_FMA4)
#define VEC_FMA(aa, bb, cc) _mm_macc_pd(aa, bb, cc)
#else
#define VEC_FMA(aa, bb, cc) _mm_add_pd(cc, _mm_mul_pd(aa, bb))
#endif
#define VEC_SCAL(aa) _mm_set1_pd(aa)
#define VEC_ZERO _mm_setzero_pd()
#define FENCE _mm_mfence()

#define DEPTH 16
#endif

#ifdef HAVE_BLUE_GENE_Q
#define VEC_SIZE 4
#define VEC_TYPE vector4double
#define VEC_LD(addr) vec_ld(0, (double *) (addr))
#define VEC_LDU(addr) ((vector4double) {(addr)[0], (addr)[1], (addr)[2], (addr)[3]})
#define VEC_ST(addr, vec)  vec_st(vec, 0, (double *) (addr))
#define VEC_STU(addr, vec) (addr)[0] = vec_extract(vec, 0); (addr)[1] = vec_extract(vec, 1); (addr)[2] = vec_extract(vec, 2); (addr)[3] = vec_extract(vec, 3)
#define VEC_FMA(aa, bb, cc) vec_madd(aa, bb, cc)
#define VEC_SCAL(aa) ((vector4double) {aa, aa, aa, aa})
#define VEC_SCAL_LD(addr) vec_lds(0, (double*) (addr))
#define VEC_ZERO ((vector4double) {0.0, 0.0, 0.0, 0.0})

#define DEPTH 16
#endif

#if defined(HAVE_BLUE_GENE) && !defined(HAVE_BLUE_GENE_Q)
#define VEC_SIZE 2
#define VEC_TYPE double _Complex
#define VEC_LD(addr) __lfpd(addr)
#define VEC_LDU(addr) __cmplx((addr)[0], (addr)[1])
#define VEC_ST(addr, vec)  __stfpd(addr, vec)
#define VEC_STU(addr, vec)  (addr)[0] = __creal(vec); (addr)[1] = __cimag(vec)
#define VEC_FMA(aa, bb, cc) __fpmadd(cc, aa, bb)
#define VEC_SCAL(aa) __cmplx(aa, aa)
#define VEC_ZERO __cmplx(0.0, 0.0)

#define DEPTH 16
#endif

#endif

#ifndef VEC_SIZE

#define VEC_SIZE 1
#define VEC_TYPE double
#define VEC_LD(addr) (addr)[0]
#define VEC_LDU(addr) VEC_LD(addr)
#define VEC_ST(addr, vec) (addr)[0] = vec
#define VEC_STU(addr, vec) VEC_ST(addr, vec)
#define VEC_FMA(aa, bb, cc) aa*bb + cc
#define VEC_SCAL(aa) aa
#define VEC_ZERO 0.0

#define DEPTH 8
#endif

#define max1(x) (((x) > 0)?(x):1)

#endif
