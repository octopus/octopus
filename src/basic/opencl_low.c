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
#include <assert.h>
#include <string.h>

#include <string_f.h>

#define MAX_PLATFORMS 4
#define MAX_DEVICES   16

int FC_FUNC(flgetplatformids, FLGETPLATFORMIDS)(const int * iplatform, cl_platform_id * platform){
  cl_platform_id all_platforms[MAX_PLATFORMS];
  cl_int status;
  cl_uint ret_platform;
  cl_uint ip;
  char info[2048];
  cl_device_id all_devices[MAX_DEVICES];
  cl_uint ret_devices;

  status = clGetPlatformIDs(MAX_PLATFORMS, all_platforms, &ret_platform);

  printf("Available OpenCL platforms: %d\n", ret_platform);

  for(ip = 0; ip < ret_platform; ip++){
    
    /*PLATFORM NAME*/
    status = clGetPlatformInfo(all_platforms[ip], CL_PLATFORM_NAME, sizeof(info), info, NULL);
    
    if (status != CL_SUCCESS){
      fprintf(stderr, "\nError: clGetPlatformInfo returned error code: %d\n", status);
      exit(1);
    }
    printf("%c Platform %d : %s", ((*iplatform == ip)?'*':' '), ip, info);

    /*PLATFORM VERSION*/
    status = clGetPlatformInfo(all_platforms[ip], CL_PLATFORM_VERSION, sizeof(info), info, NULL);
    
    if (status != CL_SUCCESS){
      fprintf(stderr, "\nError: clGetPlatformInfo returned error code: %d\n", status);
      exit(1);
    }

    printf(" %s\n", info);
  }

  printf("\n");

  *platform = all_platforms[*iplatform];

  return status;
}

int FC_FUNC_(f90_cl_get_number_of_devices, F90_CL_GET_NUMBER_OF_DEVICES)(const cl_platform_id * platform){
  cl_int status;
  cl_device_id all_devices[MAX_DEVICES];
  cl_uint ret_devices;

  status = clGetDeviceIDs(*platform, CL_DEVICE_TYPE_ALL, MAX_DEVICES, all_devices, &ret_devices);

  if (status != CL_SUCCESS){
    fprintf(stderr, "\nError: clGetDeviceIDs returned error code: %d\n", status);
    exit(1);
  }

  return ret_devices;
}

void FC_FUNC_(f90_cl_get_device_name, F90_CL_GET_DEVICE_NAME)
     (const cl_platform_id * platform, const int * idevice, STR_F_TYPE device_name STR_ARG1){
  char info[2048];
  cl_int status;
  cl_device_id all_devices[MAX_DEVICES];
  cl_uint ret_devices;

  status = clGetDeviceIDs(*platform, CL_DEVICE_TYPE_ALL, MAX_DEVICES, all_devices, &ret_devices);

  if (status != CL_SUCCESS){
    fprintf(stderr, "\nError: clGetDeviceIDs returned error code: %d\n", status);
    exit(1);
  }

  status = clGetDeviceInfo(all_devices[*idevice], CL_DEVICE_NAME, sizeof(info), info, NULL);

  if (status != CL_SUCCESS){
    fprintf(stderr, "\nError: clGetDeviceIDs returned error code: %d\n", status);
    exit(1);
  }

  TO_F_STR1(info, device_name);
}

void FC_FUNC_(f90_cl_init_context,F90_CL_INIT_CONTEXT)(const cl_platform_id * platform, cl_context * context){
  cl_int status;
  cl_context_properties cps[3];

  cps[0] = CL_CONTEXT_PLATFORM;
  cps[1] = (cl_context_properties) *platform;
  cps[2] = 0;
  
  *context = clCreateContextFromType(cps, CL_DEVICE_TYPE_ALL, NULL, NULL, &status);

  if (status != CL_SUCCESS){
    fprintf(stderr, "\nError: clCreateContextFromType returned error code: %d\n", status);
    exit(1);
  }
}

void FC_FUNC_(f90_cl_init_device,F90_CL_INIT_DEVICE)(const int * idevice, const cl_platform_id * platform, const cl_context * context, cl_device_id * device){
  size_t ParamDataBytes;
  char device_string[2048];
  cl_uint dim;
  cl_ulong mem;
  size_t max_workgroup_size;
  cl_device_id * Devices;
  cl_device_type device_type;
  cl_uint num_devices;

  if(*idevice >= 0){

    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0 , NULL, &ParamDataBytes);
    Devices = (cl_device_id*) malloc(ParamDataBytes);
    
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, ParamDataBytes, Devices, NULL);
    
    assert(sizeof(cl_device_id) == sizeof(void *));
    
    *device = Devices[*idevice];
    
    free(Devices);

  } else {

    device_type = CL_DEVICE_TYPE_DEFAULT;

    switch(*idevice){
      /* These values come from opencl.F90 */
    case -1:
      device_type = CL_DEVICE_TYPE_GPU;
      break;
    case -2:
      device_type = CL_DEVICE_TYPE_CPU;
      break;
    case -3:
      device_type = CL_DEVICE_TYPE_ACCELERATOR;
      break;
    case -4:
      device_type = CL_DEVICE_TYPE_DEFAULT;
      break;
    }

    clGetDeviceIDs(*platform, device_type, 1, device, &num_devices);

    /* Check if there were any devices of the type, otherwise just get the default.*/
    if(num_devices == 0) {
      clGetDeviceIDs(*platform, CL_DEVICE_TYPE_DEFAULT, 1, device, NULL);
    }
  }

  /* print some info about the device */
  clGetDeviceInfo(*device, CL_DEVICE_VENDOR, sizeof(device_string), &device_string, NULL);
  printf("device vendor           : %s\n", device_string);

  clGetDeviceInfo(*device, CL_DRIVER_VERSION, sizeof(device_string), &device_string, NULL);
  printf("driver version          : %s\n", device_string);

  clGetDeviceInfo(*device, CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
  printf("device name             : %s\n", device_string);

  clGetDeviceInfo (*device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &dim, NULL);
  printf("compute units           : %d\n", dim);

  clGetDeviceInfo (*device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &dim, NULL);
  printf("maximum clock frequency : %d MHz\n", dim);

  clGetDeviceInfo (*device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= (1024*1024); /* convert to megabytes */
  printf("device memory           : %ld Mb\n", mem);

  clGetDeviceInfo (*device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= 1024; /* convert to kilobytes */
  printf("local memory            : %ld Kb\n", mem);

  clGetDeviceInfo (*device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= 1024; /* convert to kilobytes */
  printf("constant memory         : %ld Kb\n", mem);


  clGetDeviceInfo (*device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(cl_ulong), &mem, NULL);
  mem /= (1024*1024); /* convert to megabytes */
  printf("device cache            : %ld Mb\n", mem);

  clGetDeviceInfo(*device, CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
  printf("cl_khr_fp64 extension   : ");
  if(strstr(device_string, "cl_khr_fp64")) printf("yes\n");
  else printf("no\n");

  printf("cl_amd_fp64 extension   : ");
  if(strstr(device_string, "cl_amd_fp64")) printf("yes\n");
  else printf("no\n");

  clGetDeviceInfo(*device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_workgroup_size), &max_workgroup_size, NULL);
  printf("maximum workgroup size  : %zd\n", max_workgroup_size);

  printf("\n");
}

/* clCreateCommandQueue */
void FC_FUNC(flcreatecommandqueue, FLCREATECOMMANDQUEUE)
     (cl_command_queue * command_queue, cl_context * context, cl_device_id * device, int * ierr){
  *command_queue = clCreateCommandQueue(*context, *device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, ierr);
}

/* clReleaseCommandQueue */
void FC_FUNC(flreleasecommandqueue, FLRELEASECOMMANDQUEUE)(cl_command_queue * command_queue, int * ierr){
  *ierr = clReleaseCommandQueue(*command_queue);
}

int FC_FUNC_(f90_cl_max_workgroup_size, F90_CL_MAX_WORKGROUP_SIZE)(cl_device_id * device){
  size_t max_workgroup_size;
  clGetDeviceInfo(*device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_workgroup_size), &max_workgroup_size, NULL);
  return (int) max_workgroup_size;
}

int FC_FUNC_(f90_cl_device_local_mem_size, F90_CL_DEVICE_LOCAL_MEM_SIZE)(cl_device_id * device){
  cl_ulong mem;

  clGetDeviceInfo (*device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem), &mem, NULL);

  return (int) mem;
}

int FC_FUNC_(f90_cl_device_max_constant_buffer_size, F90_CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE)(cl_device_id * device){
  cl_ulong mem;

  clGetDeviceInfo (*device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem), &mem, NULL);

  return (int) mem;
}

void FC_FUNC(flreleasecontext, FLRELEASECONTEXT)(cl_context * context){
  clReleaseContext(*context);
}

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

int FC_FUNC_(f90_cl_build_program, F90_CL_BUILD_PROGRAM)
     (cl_program * program, cl_context * context, cl_device_id * device, STR_F_TYPE flags_f STR_ARG1){
  char * flags;
  cl_int status;
  char device_string[2048];
  int ext_khr_fp64, ext_amd_fp64;
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

int FC_FUNC_(f90_cl_device_has_extension, F90_CL_DEVICE_HAS_EXTENSION)
     (cl_device_id * device, STR_F_TYPE extension_f STR_ARG1){
  char * extension;
  cl_int status;
  char device_string[2048];
  size_t len;
  char buffer[5000];
  int has_ext;

  TO_C_STR1(extension_f, extension);

  clGetDeviceInfo(*device, CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
  has_ext = strstr(device_string, extension) != NULL;
  free(extension);
  return has_ext;
}

void FC_FUNC_(f90_cl_release_program, F90_CL_RELEASE_PROGRAM)
     (cl_program * program, int * ierr){

  *ierr = clReleaseProgram(*program);
}

void FC_FUNC_(f90_cl_create_kernel, F90_CL_CREATE_KERNEL)
     (cl_kernel * kernel, cl_program * program, STR_F_TYPE kernel_name_f, int * ierr STR_ARG1){
  char * kernel_name;

  TO_C_STR1(kernel_name_f, kernel_name);

  *kernel = clCreateKernel(*program, kernel_name, ierr);

  free(kernel_name);
}

void FC_FUNC_(f90_cl_release_kernel, F90_CL_RELEASE_KERNEL)(cl_kernel * kernel, int * ierr){
  *ierr = clReleaseKernel(*kernel);
}

int FC_FUNC_(f90_cl_kernel_wgroup_size, F90_CL_KERNEL_WGROUP_SIZE)(cl_kernel * kernel, cl_device_id * device){
  size_t workgroup_size;
  clGetKernelWorkGroupInfo(*kernel, *device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
  return (int) workgroup_size;
}

void FC_FUNC_(f90_cl_create_buffer, F90_CL_CREATE_BUFFER)
     (cl_mem * buffer, cl_context * context, const int * flags, const size_t * size, int * ierr){

  *buffer = clCreateBuffer(*context, *flags, (size_t) *size, NULL, ierr);

}

void FC_FUNC_(f90_cl_release_buffer, F90_CL_RELEASE_BUFFER)(cl_mem * buffer, int * ierr){

  *ierr = clReleaseMemObject(*buffer);
}


/* clEnqueueWriteBuffer */
void FC_FUNC(flenqueuewritebuffer, FLENQUEUEWRITEBUFFER)
     (cl_mem * buffer, cl_command_queue * cq, const size_t * size, const size_t * offset, const void * data, int * ierr){

  *ierr = clEnqueueWriteBuffer(*cq, *buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);

}

/* clEnqueueReadBuffer */
void FC_FUNC(flenqueuereadbuffer, FLENQUEUEREADBUFFER)
     (cl_mem * buffer, cl_command_queue * cq, const size_t * size, const size_t * offset, void * data, int * ierr){

  *ierr = clEnqueueReadBuffer(*cq, *buffer, CL_TRUE, *offset, *size, data, 0, NULL, NULL);
}


void FC_FUNC(flfinish, FLFINISH)(cl_command_queue * cq, int * ierr){
  *ierr = clFinish(*cq);
}

void FC_FUNC_(f90_cl_set_kernel_arg_buf, F90_CL_SET_KERNEL_ARG_BUF)
     (cl_kernel * kernel, const int * index, cl_mem * buffer, int * ierr){

  *ierr = clSetKernelArg(*kernel, *index, sizeof(cl_mem), buffer);
}

void FC_FUNC_(f90_cl_set_kernel_arg_data, F90_CL_SET_KERNEL_ARG_DATA)
     (cl_kernel * kernel, const int * index, const int * sizeof_data, const void * data, int * ierr){
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *ierr = clSetKernelArg(*kernel, *index, *sizeof_data, data);
}

void FC_FUNC_(f90_cl_set_kernel_arg_local, F90_CL_SET_KERNEL_ARG_LOCAL)
     (cl_kernel * kernel, const int * index, const int * size_of_local, int * ierr){
  
  /* printf("kernel=%ld index=%d\n", *kernel, *index);*/

  *ierr = clSetKernelArg(*kernel, *index, *size_of_local, NULL);
}

/* clEnqueueNDRangeKernel*/
void FC_FUNC(flenqueuendrangekernel, FLENQUEUENDRANGEKERNEL)
     (cl_kernel * kernel, cl_command_queue * cq, const int * dim, const size_t * globalsizes, const size_t * localsizes, int * ierr){

  *ierr = clEnqueueNDRangeKernel(*cq, *kernel, *dim,
				NULL,  globalsizes, localsizes, 0, NULL, NULL);

}
