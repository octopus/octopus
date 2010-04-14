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
#include <stdio.h>
#include <CL/cl.h>

typedef struct{
  int numerr;
  cl_context Context;
  cl_device_id * Devices;
  cl_command_queue CommandQueue;
} opencl_t;

void FC_FUNC_(opencl_init,OPENCL_INIT)(opencl_t ** thisptr){
  size_t ParamDataBytes;
  char device_string[2048];
  cl_uint dim;
  cl_ulong mem;
  cl_platform_id platform;
  cl_int status;
  cl_context_properties cps[3];
  opencl_t * this;

  this = (opencl_t *) malloc(sizeof(opencl_t));
  *thisptr = this;

  /* Just get the first platform */
  status = clGetPlatformIDs(1, &platform, NULL);

  if(status != CL_SUCCESS){
    printf("OpenCL initialization failed: %d\n", this->numerr);
    return;
  }
  
  cps[0] = CL_CONTEXT_PLATFORM;
  cps[1] = (cl_context_properties)platform;
  cps[2] = 0;
  
  this->Context = clCreateContextFromType(cps, CL_DEVICE_TYPE_ALL,NULL, NULL, &this->numerr);

  if (this->numerr != CL_SUCCESS){
    printf("OpenCL initialization failed: %d\n", this->numerr);
    return;
  };

  clGetContextInfo(this->Context, CL_CONTEXT_DEVICES ,0 , NULL, &ParamDataBytes);
  this->Devices = (cl_device_id*) malloc(ParamDataBytes);

  clGetContextInfo(this->Context, CL_CONTEXT_DEVICES, ParamDataBytes, this->Devices, NULL);

  /* print some info about the device */
  clGetDeviceInfo(this->Devices[0], CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
  printf("OpenCL device : %s\n", device_string);

  clGetDeviceInfo (this->Devices[0], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &dim, NULL);
  printf("Compute units : %d\n", dim);

  clGetDeviceInfo (this->Devices[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= (1024*1024); /* convert to megabytes */
  printf("Device memory : %d [Mb]\n", mem);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
  printf("Extensions    : %s\n", device_string);

  /* start command queue */
  this->CommandQueue = clCreateCommandQueue(this->Context, this->Devices[0], CL_QUEUE_PROFILING_ENABLE ,&this->numerr);

}

void FC_FUNC_(opencl_end,OPENCL_END)(opencl_t ** thisptr){
  opencl_t * this;

  this = *thisptr;

  free(this->Devices);
  free(this);
}
