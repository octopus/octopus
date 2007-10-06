/*
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
*/

#include <config.h>
#if defined(__ia64__) || defined(__ia64) || defined(_M_IA64)
.text
.align 32
.global FC_FUNC_(doperate_as, DOPERATE_AS)#
.proc FC_FUNC_(doperate_as, DOPERATE_AS)#
FC_FUNC_(doperate_as, DOPERATE_AS):
.prologue
alloc r3 = ar.pfs, 8, 60, 0, 8
.body
//save LC
mov r2 = ar.lc

#define ISIZE 4
#define ISHIFT 2
#define FSIZE 8
#define FSHIFT 3
#define LDI ld4
// arguments

#define NN r65
#define ARG_WW r66
#define NR r67
#define ARG_RI r68
#define ARG_RM r69
#define ARG_FI r70
#define ARG_FO r71
#define ARG_WS r72

//copy or load input arguments
ld4 NN = [r32]
mov ARG_WW = r33
ld4 NR = [r34]
mov ARG_RI = r35
mov ARG_RM = r36
mov ARG_FI = r37
mov ARG_FO = r38
mov ARG_WS = r39

#define LL r73
#define II r74
#define WW r75
#define MAXII r77
#define RILL r78
#define WS r79
#define WS2 r80

mov LL = r0
mov II = r0

lfetch [ARG_WS]
lfetch [ARG_RI]

;;
.loop_NR:
adds LL = 1, LL;;
cmp4.ne p2,p3 = LL, NR;;

// II = ARG_RM[LL]
// MAXII = ARG_RM[LL+1]
LDI II = [ARG_RM],4;;
LDI MAXII = [ARG_RM];;

//calculate WS[1:NN] = FI + II + ARG_RI[1:NN, LL]
mov RILL = ARG_RI
mov WS = ARG_WS
adds r90 = -1, NN;;
mov ar.lc = r90

.loop_NN1:
LDI r91 = [RILL], 4;;
sxt4 r91 = r91;;
shladd r91 = r91, FSHIFT, ARG_FI;;
shladd r91 = II, FSHIFT, r91;;
st8 [WS] = r91, 8
br.cloop.sptk.few .loop_NN1
.loop_II:
;;
adds r90 = -1, NN
mov f8 = f0
mov WW = ARG_WW
mov WS = ARG_WS
mov WS2 = ARG_WS
adds II = 1, II;;
mov ar.lc = r90
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
.loop_NN2:
.mmf
(p16) ld8 r32 = [WS],8
(p16) ldfd f50 = [WW], FSIZE
(p23) fma.d f39 = f57,f67, f43
.mmb
(p17) ldfd f61 = [r33], FSIZE
(p18) st8 [WS2] = r34, 8
br.ctop.sptk.few .loop_NN2
nop.i 0
;;
//reduce accumulators
fadd.d f40 = f40, f41
fadd.d f42 = f42, f43;;
fadd.d f8 = f40, f42;;
//store
stfd.nta [ARG_FO] = f8, FSIZE
//check the value of II
cmp4.ne p4,p0 = II, MAXII
//reset register rotation
clrrrb;;

//iterate over II
(p4) br.cond.sptk.few .loop_II

//increment ARG_RI in NN
shladd ARG_RI = NN, ISHIFT, ARG_RI
lfetch [ARG_RI]
//Close the loop in LL
(p2) br.cond.sptk.few .loop_NR

;;
//restore LC
mov   ar.lc = r2

br.ret.sptk.many b0
.endp FC_FUNC_(doperate_as, DOPERATE_AS)#

#endif

// Local Variables:
// mode: c
// End:
