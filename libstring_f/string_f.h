/*
 Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch

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

 $Id$
*/

/* --------------------- Fortran to C string compatibility ---------------------- */
#if defined(_CRAY)
#include <fortran.h>

#define STR_F_TYPE _fcd
#define TO_C_STR1(s) to_c_str(s)
#define TO_C_STR2(s) to_c_str(s)
#define TO_C_STR3(s) to_c_str(s)
#define TO_F_STR1(c, f) to_f_str(c, f)
#define TO_F_STR2(c, f) to_f_str(c, f)
#define TO_F_STR3(c, f) to_f_str(c, f)
#define STR_ARG1
#define STR_ARG2
#define STR_ARG3

char *to_c_str(STR_F_TYPE f);
void to_f_str(char *c, STR_F_TYPE f);

#else

#define STR_F_TYPE char *
#define TO_C_STR1(s) to_c_str(s, l1)
#define TO_C_STR2(s) to_c_str(s, l2)
#define TO_C_STR3(s) to_c_str(s, l3)
#define TO_F_STR1(c, f) to_f_str(c, f, l1)
#define TO_F_STR2(c, f) to_f_str(c, f, l2)
#define TO_F_STR3(c, f) to_f_str(c, f, l3)
#define STR_ARG1     , unsigned long l1
#define STR_ARG2     , unsigned long l1, unsigned long l2
#define STR_ARG3     , unsigned long l1, unsigned long l2, unsigned long l3

char *to_c_str(STR_F_TYPE f, unsigned long l);
void to_f_str(char *c, STR_F_TYPE f, unsigned long l);

#endif
