#if 0
 Copyright (C) 2007 X. Andrade

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

 $Id: operate_ia64.S 2146 2006-05-23 17:36:00Z xavier $
#endif

#include <config.h>

#define ISIZE 4
#define ISHIFT 2
#define FSIZE 8
#define FSHIFT 3
#define LDI ld4

///////////////////////////////////////////////////////
//    DOPERATE
///////////////////////////////////////////////////////

.text
.align 32
.global FC_FUNC_(doperate_as, DOPERATE_AS)#
.proc FC_FUNC_(doperate_as, DOPERATE_AS)#
FC_FUNC_(doperate_as, DOPERATE_AS):

alloc r3 = ar.pfs, 8, 60, 0, 16
//save LC
mov r2 = ar.lc

// arguments

#define NN r65
#define ARG_WW r66
#define NR r67
#define ARG_RI r68
#define ARG_RM r69	
#define ARG_MX r70
#define ARG_FI r71
#define ARG_FO r72
#define ARG_WS r73

#define LL r80
#define II r81
#define WW r82
#define MAXII r83
#define RILL r84
#define WS r85
#define WS2 r86
#define MAXUNRII r87
#define NNM1 r88
#define FO_BASE r89
		
//copy or load input arguments
ld4 NN = [r32],4
;;
ld4 NR = [r32]
mov ARG_WW = r33
mov ARG_RI = r34
mov ARG_RM = r35
mov ARG_MX = r36	
mov ARG_FI = r37
mov FO_BASE = r38
mov ARG_WS = r39;;

cmp4.le p9,p0 = NR, r0;;
(p9) br.cond.spnt.few .dend_LL //if (ll == NR) do nothing

mov LL = r0
lfetch.nt1 [ARG_WW];;
adds NNM1 = -1, NN

.dloop_NR:
adds LL = 1, LL;;
cmp4.ne p2,p3 = LL, NR;;

LDI II = [ARG_RM],ISIZE            // II = ARG_RM[LL++]
shladd ARG_FO = II,FSHIFT, FO_BASE
	
LDI MAXII = [ARG_MX],ISIZE         // MAXII = ARG_MX[LL++]

//calculate WS[1:NN] = FI + II + ARG_RI[1:NN, LL]
mov RILL = ARG_RI
mov WS = ARG_WS
mov ar.ec = 5 + 1
mov pr.rot = 1 << 16;;
mov WW = ARG_WW	
mov ar.lc = NNM1;;

mov f33 = f0
mov f34 = f0
;; 

//load the addresses of the arrays in memory and perform the first iteration

.dloop_NN1:
(p16) LDI r38 = [RILL], 4
(p16) ldfd f45 = [WW], FSIZE
(p17) shladd r39 = r39, FSHIFT, ARG_FI
(p18) ldfd f56 = [r40], FSIZE
;;
(p16) sxt4 r38 = r38
(p17) shladd r39 = II, FSHIFT, r39
(p19) st8 [WS] = r41, 8
(p21) fma.d f37 = f50,f59, f39
br.ctop.sptk.few .dloop_NN1
;;

adds II = 1, II
fadd.d f8 = f38, f39
clrrrb;;
stfd.nta [ARG_FO] = f8, FSIZE
		
// FIRST, THE UNROLLED LOOP IN II
#define UNROLL 8
adds MAXUNRII = -UNROLL + 1, MAXII;;

//check whether we have to do at least one iteration
cmp.ge p7,p0 = II, MAXII        // p7 = II >= MAXII
cmp.ge p6,p0 = II, MAXUNRII;;   // p6 = II >= MAXUNRII
(p7) br.cond.dptk.few .dend_II   // if p7 goto next iteration
(p6) br.cond.dpnt.few .dloop_II  // if p6 goto rolled loop

.dloop_II_unrolled:
mov WS = ARG_WS;;
ld8 r33 = [WS],8 //load WS[0] in advance
mov WW = ARG_WW
mov WS2 = ARG_WS
adds II = UNROLL, II;;
adds r41 = FSIZE, r33
mov ar.lc = NNM1
mov ar.ec = 4 + 1
mov pr.rot = 1 << 16

//initialize accumulators
mov f8 = f0
mov f9 = f0
mov f10 = f0
mov f11 = f0
mov f12 = f0
mov f13 = f0
mov f14 = f0
mov f15 = f0;;

//inner loop in K
.dloop_NN2:

//cycle 1
(p16) ld8.s r32 = [WS],8 //load WS[K+1] 
(p16) ldfd f40 = [r33], FSIZE*2
(p20) fma.d f8 = f37,f44, f8

(p16) ldfd f33 = [WW], FSIZE
(p16) ldfd f45 = [r41], FSIZE*2
(p20) fma.d f9 = f37,f49, f9;;

//cycle 2
(p16) ldfd f50 = [r33], FSIZE*2
(p20) fma.d f10 = f37,f54, f10

(p16) adds r40 = FSIZE, r32
(p16) ldfd f55 = [r41], FSIZE*2
(p20) fma.d f11 = f37,f59, f11

 //cycle 3
(p17) st8 [WS2] = r34, 8 //store the incremented value of WS[K] for the next iteration
(p16) ldfd f60 = [r33], FSIZE*2
(p20) fma.d f12 = f37,f64, f12

(p16) ldfd f65 = [r41], FSIZE*2
(p20) fma.d f13 = f37,f69, f13

//cycle 4
(p16) ldfd f70 = [r33], FSIZE*2
(p20) fma.d f14 = f37,f74, f14

(p16) ldfd f75 = [r41], FSIZE*2
(p20) fma.d f15 = f37,f79, f15
br.ctop.dptk.few .dloop_NN2
;;
//create pointers for storing
adds r32 = FSIZE,   ARG_FO
adds r33 = FSIZE*2, ARG_FO
adds r34 = FSIZE*3, ARG_FO
adds r35 = FSIZE*4, ARG_FO
adds r36 = FSIZE*5, ARG_FO
adds r37 = FSIZE*6, ARG_FO
adds r38 = FSIZE*7, ARG_FO;;

//store
stfd.nta [ARG_FO] = f8, FSIZE*8
stfd.nta [r32]    = f9
stfd.nta [r33]    = f10
stfd.nta [r34]    = f11
stfd.nta [r35]    = f12
stfd.nta [r36]    = f13
stfd.nta [r37]    = f14
stfd.nta [r38]    = f15

//check the value of II
cmp4.lt p4,p0 = II, MAXUNRII
//reset register rotation
clrrrb;;

//iterate over II
(p4) br.cond.dptk.few .dloop_II_unrolled


// NOW WE HANDLE THE REST
.dloop_II:

//check the value of II
cmp4.lt p4,p5 = II, MAXII;;

//and go to end if we are already finished
(p5) br.cond.dpnt.few .dend_II

mov WW = ARG_WW
mov WS = ARG_WS
mov WS2 = ARG_WS
adds II = 1, II
mov ar.lc = NNM1
mov ar.ec = 7 + 1
mov pr.rot = 1 << 16
;;

//initialize accumulators
mov f33 = f0
mov f34 = f0
mov f35 = f0
mov f36 = f0

;;
//inner loop
.dloop_NN3:
(p16) ld8 r32 = [WS],8
(p16) ldfd f50 = [WW], FSIZE
(p23) fma.d f39 = f57,f67, f43
(p17) ldfd f61 = [r33], FSIZE
(p18) st8 [WS2] = r34, 8
br.ctop.sptk.few .dloop_NN3
;;
//reduce accumulators
fadd.d f40 = f40, f41
fadd.d f42 = f42, f43;;
fadd.d f8 = f40, f42;;

//store
stfd.nta [ARG_FO] = f8, FSIZE

//reset register rotation
clrrrb;;
//iterate over II
(p4) br.cond.sptk.few .dloop_II

.dend_II:

//increment ARG_RI in NN
shladd ARG_RI = NN, ISHIFT, ARG_RI
//Close the loop in LL
(p2) br.cond.sptk.few .dloop_NR

.dend_LL:

//restore LC
mov ar.lc = r2

br.ret.sptk.many b0
.endp FC_FUNC_(doperate_as, DOPERATE_AS)#




///////////////////////////////////////////////////////
//    ZOPERATE
///////////////////////////////////////////////////////

.text
.align 32
.global FC_FUNC_(zoperate_as, ZOPERATE_AS)#
.proc FC_FUNC_(zoperate_as, ZOPERATE_AS)#
FC_FUNC_(zoperate_as, ZOPERATE_AS):

alloc r3 = ar.pfs, 8, 60, 0, 16
//save LC
mov r2 = ar.lc

// arguments

#define NN r65
#define ARG_WW r66
#define NR r67
#define ARG_RI r68
#define ARG_RM r69
#define ARG_MX r70	
#define ARG_FI r71
#define ARG_FO r72
#define ARG_WS r73

#define LL r80
#define II r81
#define WW r82
#define MAXII r83
#define RILL r84
#define WS r85
#define WS2 r86
#define MAXUNRII r87
#define NNM1 r88
#define FO_BASE r89	

	
//copy or load input arguments
ld4 NN = [r32],4
ld4 NR = [r32]
mov ARG_WW = r33
mov ARG_RI = r34
mov ARG_RM = r35
mov ARG_MX = r36	
mov ARG_FI = r37
mov FO_BASE = r38
mov ARG_WS = r39;;

cmp4.le p9,p0 = NR,r0;;
(p9) br.cond.spnt.few .zend_LL //if (ll == 0) goto end

mov LL = r0

lfetch.nt1 [ARG_WW];;
adds NNM1 = -1, NN;;

.zloop_NR:
adds LL = 1, LL;;
cmp4.ne p2,p3 = LL, NR;;

LDI II = [ARG_RM],ISIZE                 // II = ARG_RM[LL++]
shladd ARG_FO = II,FSHIFT+1, FO_BASE
LDI MAXII = [ARG_MX],ISIZE              // MAXII = ARG_MX[LL++]

//calculate WS[1:NN] = FI + II + ARG_RI[1:NN, LL]
mov RILL = ARG_RI
mov WS = ARG_WS
mov ar.ec = 4 + 1
mov pr.rot = 1 << 16;;
mov ar.lc = NNM1;;

.zloop_NN1:
(p16) LDI r38 = [RILL], 4
(p17) sxt4 r39 = r39
(p18) shladd r40 = r40, FSHIFT+1, ARG_FI
(p19) shladd r41 = II, FSHIFT+1, r41
(p20) st8 [WS] = r42, 8
br.ctop.sptk.few .zloop_NN1
;;
clrrrb;;

// FIRST, THE UNROLLED LOOP IN II
#define UNROLL 8
adds MAXUNRII = -UNROLL + 1, MAXII;;

//check whether we have to do at least one iteration
cmp.ge p6,p0 = II, MAXUNRII;;   // if (II > MAXUNRII) 
(p6) br.cond.dpnt.few .zloop_II  // goto rolled loop

.zloop_II_unrolled:
mov WS = ARG_WS;;
ld8 r33 = [WS],8 //load WS[0] in advance
mov WS2 = ARG_WS
mov WW = ARG_WW
adds II = UNROLL, II
mov ar.ec = 2 + 1
mov pr.rot = 1 << 16;;
mov ar.lc = NNM1

//initialize accumulators to zero
mov f8 = f0
mov f9 = f0
mov f10 = f0
mov f11 = f0

mov f12 = f0
mov f13 = f0
mov f14 = f0
mov f15 = f0

mov f100 = f0
mov f102 = f0
mov f104 = f0
mov f106 = f0

mov f108 = f0
mov f110 = f0
mov f112 = f0
mov f114 = f0
;;

//inner loop in K
.zloop_NN2:

//cycle 1
(p16) ld8.s r32 = [WS],8 //load WS[K+1] 
(p16) ldfpd f35,f38 = [r33], FSIZE*2
(p18) fma.d f8 = f34,f37, f8

(p16) ldfd f32 = [WW], FSIZE
(p18) fma.d f9 = f34,f40, f9;;

//cycle 2
(p16) ldfpd f41,f44 = [r33], FSIZE*2
(p18) fma.d f10 = f34,f43, f10

(p17) st8 [WS2] = r34, 8 //store the incremented value of WS[K-1] for the next iteration
(p18) fma.d f11 = f34,f46, f11;;

//cycle 3
(p16) ldfpd f47,f50 = [r33], FSIZE*2
(p18) fma.d f12 = f34,f49, f12

(p18) fma.d f13 = f34,f52, f13;;

//cycle 4
(p16) ldfpd f53,f56 = [r33], FSIZE*2
(p18) fma.d f14 = f34,f55, f14

(p18) fma.d f15 = f34,f58, f15;;

 //cycle 5
(p16) ldfpd f59,f62 = [r33], FSIZE*2
(p18) fma.d f101 = f34,f61, f102

(p18) fma.d f103 = f34,f64, f104;;

 //cycle 6
(p16) ldfpd f65,f68= [r33], FSIZE*2
(p18) fma.d f105 = f34,f67, f106

(p18) fma.d f107 = f34,f70, f108;;

//cycle 7
(p16) ldfpd f71,f74 = [r33], FSIZE*2
(p18) fma.d f109 = f34,f73, f110

(p18) fma.d f111 = f34,f76, f112;;

//cycle 8
(p16) ldfpd f77,f80 = [r33], FSIZE*2
(p18) fma.d f113 = f34,f79, f114

(p18) fma.d f115 = f34,f82, f116
br.ctop.dptk.few .zloop_NN2
;;

//create pointers for storing
adds r32 = FSIZE,   ARG_FO
adds r33 = FSIZE*2, ARG_FO
adds r34 = FSIZE*3, ARG_FO

adds r35 = FSIZE*4, ARG_FO
adds r36 = FSIZE*5, ARG_FO
adds r37 = FSIZE*6, ARG_FO
adds r38 = FSIZE*7, ARG_FO

adds r39 = FSIZE*8, ARG_FO
adds r40 = FSIZE*9, ARG_FO
adds r41 = FSIZE*10, ARG_FO
adds r42 = FSIZE*11, ARG_FO

adds r43 = FSIZE*12, ARG_FO
adds r44 = FSIZE*13, ARG_FO
adds r45 = FSIZE*14, ARG_FO
adds r46 = FSIZE*15, ARG_FO;;

//store
stfd.nta [ARG_FO] = f8, FSIZE*16
stfd.nta [r32]    = f9
stfd.nta [r33]    = f10
stfd.nta [r34]    = f11

stfd.nta [r35]    = f12
stfd.nta [r36]    = f13
stfd.nta [r37]    = f14
stfd.nta [r38]    = f15

stfd.nta [r39]    = f102
stfd.nta [r40]    = f104
stfd.nta [r41]    = f106
stfd.nta [r42]    = f108

stfd.nta [r43]    = f110
stfd.nta [r44]    = f112
stfd.nta [r45]    = f114
stfd.nta [r46]    = f116

//check the value of II
cmp4.lt p4,p0 = II, MAXUNRII
//reset register rotation
clrrrb;;

//iterate over II
(p4) br.cond.dptk.few .zloop_II_unrolled


// NOW WE HANDLE THE REST
.zloop_II:

//check the value of II
cmp4.lt p4,p5 = II, MAXII;;

//and go to end if we are already finished
(p5) br.cond.dpnt.few .zend_II

mov WW = ARG_WW
mov WS = ARG_WS
mov WS2 = ARG_WS
adds II = 1, II
mov ar.lc = NNM1
mov ar.ec = 8 + 1
mov pr.rot = 1 << 16
;;

//initialize accumulators
mov f33 = f0
mov f34 = f0

mov f36 = f0
mov f37 = f0
;;
//inner loop
.zloop_NN3:
(p16) ld8 r32 = [WS],8
(p16) ldfd f50 = [WW], FSIZE
(p24) fma.d f40 = f58,f66, f42

(p17) ldfpd f59,f68 = [r33], FSIZE*2
(p24) fma.d f43 = f58,f75, f45
(p18) st8 [WS2] = r34, 8
br.ctop.sptk.few .zloop_NN3
;;

adds r32 = FSIZE, ARG_FO
//reduce accumulators
fadd.d f8 = f41, f42
fadd.d f9 = f44, f45;;

//store
stfd.nta [ARG_FO] = f8, FSIZE*2
stfd.nta [r32] = f9

//reset register rotation
clrrrb;;
//iterate over II
(p4) br.cond.sptk.few .zloop_II

.zend_II:

//increment ARG_RI in NN
shladd ARG_RI = NN, ISHIFT, ARG_RI
//Close the loop in LL
(p2) br.cond.sptk.few .zloop_NR

.zend_LL:

//restore LC
mov ar.lc = r2
;;

br.ret.sptk.many b0
.endp FC_FUNC_(zoperate_as, ZOPERATE_AS)#

// Local Variables:
// mode: C++
// End:
