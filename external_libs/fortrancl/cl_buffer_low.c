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

void FC_FUNC_(flcreatebuffer_low, FLCREATEBUFFER_LOW)
     (cl_context * context, const int * flags, const size_t * size, int * errcode_ret, cl_mem * buffer){
  
  cl_int errcode_ret_cl;

  *buffer = clCreateBuffer(*context, (cl_mem_flags) *flags, *size, NULL, &errcode_ret_cl);
  *errcode_ret = (int) errcode_ret_cl;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC(flreleasememobject, FLRELEASEMEMOBJECT)(cl_mem * memobj, int * status){

  *status = (int)clReleaseMemObject(*memobj);
}

/* -----------------------------------------------------------------------*/
