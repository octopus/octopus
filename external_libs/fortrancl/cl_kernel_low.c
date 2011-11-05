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

void FC_FUNC_(f90_cl_create_kernel, F90_CL_CREATE_KERNEL)
     (cl_kernel * kernel, cl_program * program, STR_F_TYPE kernel_name_f, int * status STR_ARG1){
  char * kernel_name;

  TO_C_STR1(kernel_name_f, kernel_name);

  *kernel = clCreateKernel(*program, kernel_name, status);

  free(kernel_name);
}


/* -----------------------------------------------------------------------*/

void FC_FUNC(clreleasekernel, CLRELEASEKERNEL)(cl_kernel * kernel, int * status){
  *status = (int) clReleaseKernel(*kernel);
}

/* -----------------------------------------------------------------------*/

int FC_FUNC_(f90_cl_kernel_wgroup_size, F90_CL_KERNEL_WGROUP_SIZE)(cl_kernel * kernel, cl_device_id * device){
  size_t workgroup_size;
  clGetKernelWorkGroupInfo(*kernel, *device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
  return (int) workgroup_size;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_set_kernel_arg_buf, F90_CL_SET_KERNEL_ARG_BUF)
     (cl_kernel * kernel, const int * index, cl_mem * buffer, int * status){

  *status = clSetKernelArg(*kernel, *index, sizeof(cl_mem), buffer);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_set_kernel_arg_data, F90_CL_SET_KERNEL_ARG_DATA)
     (cl_kernel * kernel, const int * index, const int * sizeof_data, const void * data, int * status){
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *status = clSetKernelArg(*kernel, *index, *sizeof_data, data);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_set_kernel_arg_local, F90_CL_SET_KERNEL_ARG_LOCAL)
     (cl_kernel * kernel, const int * index, const int * size_of_local, int * status){
  
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *status = clSetKernelArg(*kernel, *index, *size_of_local, NULL);
}

/* -----------------------------------------------------------------------*/
