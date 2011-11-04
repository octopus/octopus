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

void FC_FUNC_(f90_cl_create_program_from_file, F90_CL_CREATE_PROGRAM_FROM_FILE)
     (cl_program * program, cl_context * context, STR_F_TYPE file_name_f STR_ARG1){
  FILE * source_file;
  size_t szSourceLength;
  char* cSourceString;
  char * file_name;
  cl_int status;

  TO_C_STR1(file_name_f, file_name);

  /* open the OpenCL source code file */
  source_file = fopen(file_name, "rb");
  if(source_file == 0){
    fprintf(stderr, "Error: Failed to open file %s\n", file_name);
    exit(1);
  } else {
    printf("Info: compiling OpenCL code %s\n", file_name);
  }

  /* get the length of the source code */
  fseek(source_file, 0, SEEK_END); 
  szSourceLength = ftell(source_file);
  fseek(source_file, 0, SEEK_SET); 
  
  /* allocate a buffer for the source code string and read it in */
  cSourceString = (char *) malloc((szSourceLength + 1)*sizeof(char));
  fread(cSourceString, szSourceLength, 1, source_file);
  fclose(source_file);
    
  cSourceString[szSourceLength] = '\0';

  *program = clCreateProgramWithSource(*context, 1, (const char**)&cSourceString, NULL, &status);

  if(status != CL_SUCCESS){
    fprintf(stderr, "Error: program creation %s failed.\n", file_name);
    exit(1);
  }

  free(file_name);

}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_build_program, F90_CL_BUILD_PROGRAM)
     (cl_program * program, cl_context * context, cl_device_id * device, STR_F_TYPE flags_f STR_ARG1){
  char * flags;
  cl_int status;
  size_t len;
  char buffer[5000];

  TO_C_STR1(flags_f, flags);

  status = clBuildProgram(*program, 0, NULL, flags, NULL, NULL);
  
  clGetProgramBuildInfo(*program, *device,
			CL_PROGRAM_BUILD_LOG, sizeof (buffer), buffer,
			&len);

  /* Print the compilation log */
  if(len > 2) printf("%s\n\n", buffer);

  if(status != CL_SUCCESS){
    fprintf(stderr, "Error: compilation failed.\n");
    exit(1);
  }

  free(flags);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_release_program, F90_CL_RELEASE_PROGRAM)
     (cl_program * program, int * status){

  *status = clReleaseProgram(*program);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_create_kernel, F90_CL_CREATE_KERNEL)
     (cl_kernel * kernel, cl_program * program, STR_F_TYPE kernel_name_f, int * status STR_ARG1){
  char * kernel_name;

  TO_C_STR1(kernel_name_f, kernel_name);

  *kernel = clCreateKernel(*program, kernel_name, status);

  free(kernel_name);
}


/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_release_kernel, F90_CL_RELEASE_KERNEL)(cl_kernel * kernel, int * status){
  *status = clReleaseKernel(*kernel);
}

/* -----------------------------------------------------------------------*/

int FC_FUNC_(f90_cl_kernel_wgroup_size, F90_CL_KERNEL_WGROUP_SIZE)(cl_kernel * kernel, cl_device_id * device){
  size_t workgroup_size;
  clGetKernelWorkGroupInfo(*kernel, *device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
  return (int) workgroup_size;
}

/* -----------------------------------------------------------------------*/

/* clEnqueueWriteBuffer */
void FC_FUNC(flenqueuewritebuffer, FLENQUEUEWRITEBUFFER)
     (cl_mem * buffer, cl_command_queue * cq, const size_t * size, const size_t * offset, const void * data, int * status){

  *status = clEnqueueWriteBuffer(*cq, *buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);

}

/* -----------------------------------------------------------------------*/

/* clEnqueueReadBuffer */
void FC_FUNC(flenqueuereadbuffer, FLENQUEUEREADBUFFER)
     (cl_mem * buffer, cl_command_queue * cq, const size_t * size, const size_t * offset, void * data, int * status){

  *status = clEnqueueReadBuffer(*cq, *buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);
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
