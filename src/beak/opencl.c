/*
 Copyright (C) 2010 X. Andrade, N. Suberviola

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
#include <CL/cl.h>

typedef struct{
  int numerr;
  cl_context GPUContext;
  cl_device_id * GPUDevices;
  cl_command_queue GPUCommandQueue;
} opencl_t;

void FC_FUNC_(opencl_init,OPENCL_INIT)(opencl_t * this){
  size_t ParamDataBytes, leng;

  this->GPUContext = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU,NULL, NULL, &this->numerr);
  if (this->numerr != CL_SUCCESS) return;

  clGetContextInfo(this->GPUContext, CL_CONTEXT_DEVICES ,0 , NULL, &ParamDataBytes);
  this->GPUDevices = (cl_device_id*) malloc(ParamDataBytes);

  clGetContextInfo(this->GPUContext, CL_CONTEXT_DEVICES, ParamDataBytes, this->GPUDevices, NULL);

  this->GPUCommandQueue = clCreateCommandQueue(this->GPUContext, this->GPUDevices[0],CL_QUEUE_PROFILING_ENABLE ,&this->numerr);

}

void FC_FUNC_(opencl_end,OPENCL_END)(opencl_t * this){
  free(this->GPUDevices);
}
