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

/* clCreateCommandQueue */
void FC_FUNC_(flcreatecommandqueue_low, FLCREATECOMMANDQUEUE_LOW)
     (cl_context * context, cl_device_id * device, const int * properties, int * status, cl_command_queue * command_queue){
  cl_int status_cl;
  *command_queue = clCreateCommandQueue(*context, *device, (cl_command_queue_properties) *properties, &status_cl);
  *status = (int) status_cl;
}

/* -----------------------------------------------------------------------*/

/* clReleaseCommandQueue */
void FC_FUNC(flreleasecommandqueue, FLRELEASECOMMANDQUEUE)(cl_command_queue * command_queue, int * status){
  *status = clReleaseCommandQueue(*command_queue);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC(flfinish, FLFINISH)(cl_command_queue * cq, int * status){
  *status = clFinish(*cq);
}

/* -----------------------------------------------------------------------*/

/* clEnqueueNDRangeKernel*/
void FC_FUNC(flenqueuendrangekernel, FLENQUEUENDRANGEKERNEL)
     (cl_command_queue * command_queue, cl_kernel * kernel, const int * dim, 
      const size_t * globalsizes, const size_t * localsizes, int * status){

  *status = (int) clEnqueueNDRangeKernel(*command_queue, *kernel, *dim,
					 NULL,  globalsizes, localsizes, 0, NULL, NULL);

}

/* -----------------------------------------------------------------------*/

/* clEnqueueWriteBuffer */
void FC_FUNC(flenqueuewritebuffer, FLENQUEUEWRITEBUFFER)
     (cl_command_queue * command_queue, cl_mem * buffer, const int * blocking_write, 
      const cl_long * offset, const cl_long * cb, const void * ptr, int * status){

  *status = (int) clEnqueueWriteBuffer(*command_queue, *buffer, (cl_bool) * blocking_write, 
				       (size_t) *offset, (size_t) *cb, ptr, 0, NULL, NULL);

}

/* -----------------------------------------------------------------------*/

/* clEnqueueReadBuffer */
void FC_FUNC(flenqueuereadbuffer, FLENQUEUEREADBUFFER)
     (cl_command_queue * command_queue, cl_mem * buffer, const int * blocking_read, 
      const cl_long * offset, const cl_long * cb, void * ptr, int * status){

  *status = (int) clEnqueueReadBuffer(*command_queue, *buffer, (cl_bool) blocking_read, 
				      (size_t) *offset, (size_t) *cb, ptr, 0, NULL, NULL);
}
