/*
 Copyright (C) 2010-2011 X. Andrade

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

 $Id: opencl.c 2146 2006-05-23 17:36:00Z xavier $
*/


#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>
#include <string_f.h>
#include <string.h>
#include <assert.h>
#include <string.h>

#include <string_f.h>

/* -----------------------------------------------------------------------*/

void FC_FUNC(flreleasecontext, FLRELEASECONTEXT)(cl_context * context, int * status){
  *status = (int) clReleaseContext(*context);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flcreatecontext_low, CLCREATECONTEXT_LOW)
     (const cl_platform_id * platform, const int * num_devices, const cl_device_id * devices, int * errcode_ret, cl_context * context){
  cl_int errcode_ret_cl;
  cl_context_properties context_properties[3];

  context_properties[0] = CL_CONTEXT_PLATFORM;
  context_properties[1] = (cl_context_properties) *platform;
  context_properties[2] = 0;
  
  *context = clCreateContext(context_properties, (cl_uint) *num_devices, devices, NULL, NULL, &errcode_ret_cl);
  *errcode_ret = (int) errcode_ret_cl;
  
}
