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
alloc r3 = ar.pfs, 7, 89, 0, 0
.body
//save LC
mov r2 = ar.lc

// arguments
#define ARG_NN r32
#define ARG_WW r33
#define ARG_NR r34
#define ARG_OR r35
#define ARG_RM r36
#define ARG_FI r37
#define ARG_FO r38

#define NN r40
#define NR r41
#define LL r42
#define II r43
#define WW r44
#define OR r45
#define MAXII r46
#define ORLL r47
#define FO   r48
#define FIII r49

ld4 NN = [ARG_NN]
ld4 NR = [ARG_NR]
mov LL = r0
mov II = r0
mov OR = ARG_OR
mov FO = ARG_FO
;;
.loop_NR:
adds LL = 1, LL;;
cmp4.ne p2,p3 = LL, NR;;

// II = ARG_RM[LL]
// MAXII = ARG_RM[LL+1]
ld4 II = [ARG_RM],4;;
ld4 MAXII = [ARG_RM];;

.loop_II:

adds r69 = -1, NN
shladd FIII = II, 3, ARG_FI
mov f40 = f0
mov ORLL = OR
mov WW = ARG_WW
adds II = 1, II
mov ar.lc = r69
//inner loop
;;
.loop_NN:
ldfd f41 = [WW], 8
ld4 r68 = [ORLL], 4;;
sxt4 r68 = r68;;
shladd r63 = r68, 3, FIII;;
ldfd f42 = [r63];;
fma.d f40 = f41,f42, f40
br.cloop.sptk.few .loop_NN
;;

stfd.nta [FO] = f40, 8
cmp4.ne p4,p0 = II, MAXII;;
//iterate over II
(p4) br.cond.sptk.few .loop_II


shladd OR = NN, 2, OR
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
