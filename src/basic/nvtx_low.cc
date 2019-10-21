/*
 Copyright (C) 2019 S. Ohlmann

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

/* This is a wrapper around the NVTX (NVIDIA Tools Extension) profiling functions */

#include <config.h>

#ifdef HAVE_NVTX
#include <nvToolsExt.h>
#endif

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include <fortran_types.h>

using namespace std;

extern "C" void FC_FUNC_(nvtx_range_push, NVTX_RANGE_PUSH)(STR_F_TYPE range_name STR_ARG1){
#ifdef HAVE_NVTX  
  char *range_name_c;
  TO_C_STR1(range_name, range_name_c);

  nvtxRangePushA(range_name_c);
  
  free(range_name_c);
#endif
}

extern "C" void FC_FUNC_(nvtx_range_pop, NVTX_RANGE_POP)(){
#ifdef HAVE_NVTX  
  nvtxRangePop();
#endif
}
