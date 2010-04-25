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
#include <string_f.h>
#include <string.h>

#include "opencl.h"

void FC_FUNC_(f90_opencl_env_init,F90_OPENCL_ENV_INIT)(opencl_env_t ** thisptr, STR_F_TYPE source_path_f STR_ARG1){
  size_t ParamDataBytes;
  char device_string[2048];
  cl_uint dim;
  cl_ulong mem;
  cl_platform_id platform;
  cl_int status;
  cl_context_properties cps[3];
  opencl_env_t * this;

  this = (opencl_env_t *) malloc(sizeof(opencl_env_t));
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
  clGetDeviceInfo(this->Devices[0], CL_DEVICE_VENDOR, sizeof(device_string), &device_string, NULL);
  printf("OpenCL device         : %s", device_string);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
  printf(" %s\n", device_string);

  clGetDeviceInfo (this->Devices[0], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &dim, NULL);
  printf("Compute units         : %d\n", dim);

  clGetDeviceInfo (this->Devices[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= (1024*1024); /* convert to megabytes */
  printf("Device memory         : %d [Mb]\n", mem);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
  printf("Extensions            : %s\n", device_string);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(this->max_workgroup_size), &this->max_workgroup_size, NULL);
  printf("Maximum workgroup size: %d\n", this->max_workgroup_size);

  /* start command queue */
  this->CommandQueue = clCreateCommandQueue(this->Context, this->Devices[0], CL_QUEUE_PROFILING_ENABLE ,&this->numerr);

  TO_C_STR1(source_path_f, this->source_path);
}

void FC_FUNC_(f90_opencl_env_end,F90_OPENCL_ENV_END)(opencl_env_t ** thisptr){
  opencl_env_t * this;

  this = *thisptr;

  free(this->source_path);
  free(this->Devices);
  free(this);
}

void opencl_kernel_init(opencl_kernel_t * this, opencl_env_t * env, const char * file_name, char * kernel_base_name){
  FILE * source_file;
  cl_int numerr;
  cl_program OpenCLProgram;
  int len;
  char * kernel_name;
  char * full_file_name;
  size_t szSourceLength;
  char* cSourceString;

  /* build the full path of the source file */
  full_file_name = (char *) malloc((strlen(env->source_path) + strlen(file_name) + 1)*sizeof(char));
  strcpy(full_file_name, env->source_path);
  strcat(full_file_name, file_name);

  /* open the OpenCL source code file */
  source_file = fopen(full_file_name, "rb");
  if(source_file == 0){
    fprintf(stderr, "Error: Failed to open file %s\n", full_file_name);
    exit(1);
  } else {
    printf("Info: compiling OpenCL code %s\n", full_file_name);
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

  OpenCLProgram = clCreateProgramWithSource(env->Context, 1, (const char**)&cSourceString, NULL, &numerr); 
  numerr = clBuildProgram(OpenCLProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
  
  if(numerr != CL_SUCCESS){
    size_t len;
    char buffer[2048];
    clGetProgramBuildInfo (OpenCLProgram, env->Devices[0],
			   CL_PROGRAM_BUILD_LOG, sizeof (buffer), buffer,
			   &len);    
    fprintf(stderr, "Error: compilation of file %s failed.\nCompilation log:\n%s\n", full_file_name, buffer);
    exit(1);
  }

  len = strlen(kernel_base_name);
  kernel_name = (char *) malloc((len + 2)*sizeof(char));

  /* Create a handle to the compiled OpenCL function (Kernel) */
  kernel_name[0] = 'd';
  kernel_name[1] = '\0';
  strcat(kernel_name, kernel_base_name);
  this->kernel_double = clCreateKernel(OpenCLProgram, kernel_name, &numerr);

  if(numerr != CL_SUCCESS){
    fprintf(stderr, "Error: creation of kernel '%s' failed with error %d.\n", kernel_name, numerr);
    exit(1);
  }

  kernel_name[0] = 'z';
  kernel_name[1] = '\0';
  strcat(kernel_name, kernel_base_name);
  this->kernel_complex = clCreateKernel(OpenCLProgram, kernel_name, &numerr);

  if(numerr != CL_SUCCESS){
    fprintf(stderr, "Error: creation of kernel  '%s' failed with error %d.\n", kernel_name, numerr);
    exit(1);
  }

  free(kernel_name);
  free(full_file_name);
}

void FC_FUNC_(f90_opencl_kernel_init, F90_OPENCL_KERNEL_INIT)
     (opencl_kernel_t ** this, opencl_env_t ** env, STR_F_TYPE file_name_f, STR_F_TYPE kernel_base_name_f STR_ARG2){
  char * file_name,  * kernel_base_name;

  *this = (opencl_kernel_t *) malloc(sizeof(opencl_kernel_t));
  TO_C_STR1(file_name_f, file_name);
  TO_C_STR2(kernel_base_name_f, kernel_base_name);
  opencl_kernel_init(*this, *env, file_name, kernel_base_name);
  free(file_name);
  free(kernel_base_name);
}

void opencl_kernel_end(opencl_kernel_t * this){
  clReleaseKernel(this->kernel_double);
  clReleaseKernel(this->kernel_complex);
}


void FC_FUNC_(f90_opencl_kernel_end, F90_OPENCL_KERNEL_END)(opencl_kernel_t ** this){
  opencl_kernel_end(*this);
  free(*this);
}

