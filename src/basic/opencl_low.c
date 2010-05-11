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

void FC_FUNC_(f90_cl_env_init,F90_CL_ENV_INIT)(opencl_env_t ** env, STR_F_TYPE source_path_f STR_ARG1){
  size_t ParamDataBytes;
  char device_string[2048];
  cl_uint dim;
  cl_ulong mem;
  cl_platform_id platform;
  cl_int status;
  cl_context_properties cps[3];
  opencl_env_t * this;
  size_t max_workgroup_size;

  this = (opencl_env_t *) malloc(sizeof(opencl_env_t));
  *env = this;

  /* Just get the first platform */
  status = clGetPlatformIDs(1, &platform, NULL);

  if(status != CL_SUCCESS){
    printf("OpenCL initialization failed. Error code: %d\n", status);
    exit(1);
  }
  
  cps[0] = CL_CONTEXT_PLATFORM;
  cps[1] = (cl_context_properties)platform;
  cps[2] = 0;
  
  this->Context = clCreateContextFromType(cps, CL_DEVICE_TYPE_ALL, NULL, NULL, &status);

  if (status != CL_SUCCESS){
    printf("OpenCL initialization failed. Error code: %d\n", status);
    exit(1);
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
  printf("Device memory         : %ld [Mb]\n", mem);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
  printf("Extensions            : %s\n", device_string);

  clGetDeviceInfo(this->Devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_workgroup_size), &max_workgroup_size, NULL);
  printf("Maximum workgroup size: %zd\n", max_workgroup_size);

  /* start command queue */
  this->CommandQueue = clCreateCommandQueue(this->Context, this->Devices[0], CL_QUEUE_PROFILING_ENABLE, &status);

  TO_C_STR1(source_path_f, this->source_path);

}

int FC_FUNC_(f90_cl_max_workgroup_size, F90_CL_MAX_WORKGROUP_SIZE)(opencl_env_t ** env){
  size_t max_workgroup_size;
  clGetDeviceInfo(env[0]->Devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_workgroup_size), &max_workgroup_size, NULL);
  return (int) max_workgroup_size;
}

void FC_FUNC_(f90_cl_env_end,F90_CL_ENV_END)(opencl_env_t ** env){
  opencl_env_t * this;

  this = *env;
  clReleaseCommandQueue(this->CommandQueue);
  clReleaseContext(this->Context);
  free(this->source_path);
  free(this->Devices);
  free(this);
}

void FC_FUNC_(f90_cl_build_program, F90_CL_BUILD_PROGRAM)
     (cl_program ** program, opencl_env_t ** env, STR_F_TYPE file_name_f STR_ARG1){
  FILE * source_file;
  char * full_file_name;
  size_t szSourceLength;
  char* cSourceString;
  char * file_name;
  cl_int status;

  *program = (cl_program *) malloc(sizeof(cl_program));

  TO_C_STR1(file_name_f, file_name);

  /* build the full path of the source file */
  full_file_name = (char *) malloc((strlen(env[0]->source_path) + strlen(file_name) + 1)*sizeof(char));
  strcpy(full_file_name, env[0]->source_path);
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

  **program = clCreateProgramWithSource(env[0]->Context, 1, (const char**)&cSourceString, NULL, &status);
  status = clBuildProgram(**program, 0, NULL, "-cl-mad-enable", NULL, NULL);
  
  if(status != CL_SUCCESS){
    size_t len;
    char buffer[2048];
    clGetProgramBuildInfo (**program, env[0]->Devices[0],
			   CL_PROGRAM_BUILD_LOG, sizeof (buffer), buffer,
			   &len);    
    fprintf(stderr, "Error: compilation of file %s failed.\nCompilation log:\n%s\n", full_file_name, buffer);
    exit(1);
  }

  free(file_name);

}

void FC_FUNC_(f90_cl_release_program, F90_CL_RELEASE_PROGRAM)
     (cl_program ** program, int * ierr){
  cl_int status;

  *ierr = clReleaseProgram(**program);
  if(*ierr /= CL_SUCCESS) free(*program);
}

void FC_FUNC_(f90_cl_create_kernel, F90_CL_CREATE_KERNEL)
     (cl_kernel ** kernel, cl_program ** program, STR_F_TYPE kernel_name_f, int * ierr STR_ARG1){
  char * kernel_name;
  cl_int status;

  TO_C_STR1(kernel_name_f, kernel_name);

  *kernel = (cl_kernel *) malloc(sizeof(cl_kernel));

  **kernel = clCreateKernel(**program, kernel_name, ierr);

  /*printf("kernel = %ld\n", *kernel);*/

  free(kernel_name);
}

void FC_FUNC_(f90_cl_release_kernel, F90_CL_RELEASE_KERNEL)(cl_kernel ** kernel, int * ierr){
  *ierr = clReleaseKernel(**kernel);
  free(*kernel);
}

int FC_FUNC_(f90_cl_kernel_wgroup_size, F90_CL_KERNEL_WGROUP_SIZE)(cl_kernel ** kernel, opencl_env_t ** env){
  size_t workgroup_size;
  clGetKernelWorkGroupInfo(**kernel, env[0]->Devices[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
  return (int) workgroup_size;
}

void FC_FUNC_(f90_cl_create_buffer, F90_CL_CREATE_BUFFER)
     (cl_mem ** buffer, opencl_env_t ** env, const int * flags, const size_t * size, int * ierr){

  *buffer = (cl_mem *) malloc(sizeof(cl_mem));
  
  /*
  printf("\nCreateBuffer\n");
  printf("queue=%ld buffer=%ld size=%d\n", env[0]->CommandQueue, **buffer, *size);
  */

  **buffer = clCreateBuffer(env[0]->Context, *flags, (size_t) *size, NULL, ierr);

}

void FC_FUNC_(f90_cl_release_buffer, F90_CL_RELEASE_BUFFER)(cl_mem ** buffer, int * ierr){

  *ierr = clReleaseMemObject(**buffer);
  free(*buffer);
}

void FC_FUNC_(f90_cl_write_buffer, F90_CL_WRITE_BUFFER)
     (cl_mem ** buffer, opencl_env_t ** env, const size_t * size, const size_t * offset, const void * data){
  cl_int ierr;

  /*  
      printf("\nWriteBuffer\n");
      printf("queue=%ld buffer=%ld data=%ld offest=%ld size=%d\n", env[0]->CommandQueue, **buffer, data, *offset, *size);
  */

  ierr = clEnqueueWriteBuffer(env[0]->CommandQueue, **buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);

  if(ierr != CL_SUCCESS){
    fprintf(stderr, "Error: buffer write failed. Error code %d\n", ierr);
    exit(1);
  }

}


void FC_FUNC_(f90_cl_read_buffer, F90_CL_READ_BUFFER)
     (cl_mem ** buffer, opencl_env_t ** env, const size_t * size, const size_t * offset, void * data){
  cl_int ierr;
  
  /*printf("queue=%ld buffer=%ld data=%ld offest=%ld size=%d\n", env[0]->CommandQueue, **buffer, data, *offset, *size);*/

  ierr = clEnqueueReadBuffer(env[0]->CommandQueue, **buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);

  if(ierr != CL_SUCCESS){
    fprintf(stderr, "Error: buffer read failed. Error code %d\n", ierr);
    exit(1);
  }

}


void FC_FUNC_(f90_cl_finish, F90_CL_FINISH)(opencl_env_t ** env){
  cl_int ierr;

  ierr = clFinish(env[0]->CommandQueue);

  if(ierr != CL_SUCCESS){
    fprintf(stderr, "Error: clFinish failed. Error code %d\n", ierr);
    exit(1);
  }

}

void FC_FUNC_(f90_cl_set_kernel_arg_buf, F90_CL_SET_KERNEL_ARG_BUF)
     (cl_kernel ** kernel, const int * index, cl_mem ** buffer, int * ierr){

  /*printf("index=%d\n", *index);*/

  *ierr = clSetKernelArg(**kernel, *index, sizeof(cl_mem), *buffer);
}

void FC_FUNC_(f90_cl_set_kernel_arg_data, F90_CL_SET_KERNEL_ARG_DATA)
     (cl_kernel ** kernel, const int * index, const int * sizeof_data, const void * data, int * ierr){
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *ierr = clSetKernelArg(**kernel, *index, *sizeof_data, data);
}

void FC_FUNC_(f90_cl_set_kernel_arg_local, F90_CL_SET_KERNEL_ARG_LOCAL)
     (cl_kernel ** kernel, const int * index, const int * size_of_local, int * ierr){
  
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *ierr = clSetKernelArg(**kernel, *index, *size_of_local, NULL);
}

void FC_FUNC_(f90_cl_kernel_run, F90_CL_KERNEL_RUN)
     (cl_kernel ** kernel, opencl_env_t ** env, const int * dim, const size_t * globalsizes, const size_t * localsizes, int * ierr){

  /*
  cl_uint numargs;

  clGetKernelInfo(**kernel, CL_KERNEL_NUM_ARGS, sizeof(numargs), &numargs, NULL);
  printf("Num args %d\n", numargs);
  */
  /*  printf("WG size  %d %d\n", globalsizes[0], localsizes[0]);
   */

  *ierr = clEnqueueNDRangeKernel(env[0]->CommandQueue, **kernel, *dim,
				NULL,  globalsizes, localsizes, 0, NULL, NULL);

}
