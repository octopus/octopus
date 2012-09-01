/*
 Copyright (C) 2012 X. Andrade

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

 $Id: cl_global.h 2146 2006-05-23 17:36:00Z xavier $
*/

#ifndef __CL_GLOBAL_H__
#define __CL_GLOBAL_H__

#ifdef EXT_KHR_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif EXT_AMD_FP64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
// We can use to printf from OpenCL, useful for debugging
//#pragma OPENCL EXTENSION cl_amd_printf:enable
#endif

#endif

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
