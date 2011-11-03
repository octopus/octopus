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

int FC_FUNC(flgetplatformids, FLGETPLATFORMIDS)(const int * iplatform, cl_platform_id * platform){
  cl_platform_id all_platforms[MAX_PLATFORMS];
  cl_int status;
  cl_uint ret_platform;
  cl_uint ip;
  char info[2048];

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

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceids_num, FLGETDEVICEIDS_NUM)
     (const cl_platform_id * platform, const int * device_type, int * num_devices, int * status){
  cl_uint unum_devices;

  *status = (int) clGetDeviceIDs(*platform, *device_type, 0, NULL, &unum_devices);
  *num_devices = (int) unum_devices;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceids_listall, FLGETDEVICEIDS_LISTALL)
     (const cl_platform_id * platform, const int * device_type, const int * num_entries, cl_device_id * devices, 
      int * num_devices, int * status){

  cl_uint unum_devices;

  *status = (int) clGetDeviceIDs(*platform, *device_type, (cl_uint) *num_entries, devices, &unum_devices);
  *num_devices = (int) unum_devices;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceids_getdev, FLGETDEVICEIDS_GETDEV)
     (const cl_device_id * alldevices, const int * idevice, cl_device_id * device){
  *device = alldevices[*idevice];
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceinfo_str, FLGETDEVICEINFO_STR)
     (const cl_device_id * device, const int * param_name, STR_F_TYPE param_value, int * status STR_ARG1){
  char info[2048];

  *status = (int) clGetDeviceInfo(*device, (cl_device_info) *param_name, sizeof(info), info, NULL);

  TO_F_STR1(info, param_value);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceinfo_int64, FLGETDEVICEINFO_INT64)
     (const cl_device_id * device, const int * param_name, cl_long * param_value, int * status){
  union { 
    cl_uint  val_uint;
    cl_bool  val_bool;
    cl_ulong val_ulong;
    size_t   val_size_t;
  } rval;
  
  *status = (int) clGetDeviceInfo(*device, (cl_device_info) *param_name, sizeof(rval), &rval, NULL);

  if(*status != CL_SUCCESS) return;
  
  switch(*param_name){
    /* return cl_uint*/
  case CL_DEVICE_ADDRESS_BITS:
  case CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE:
  case CL_DEVICE_MAX_CLOCK_FREQUENCY:
  case CL_DEVICE_MAX_COMPUTE_UNITS:
  case CL_DEVICE_MAX_CONSTANT_ARGS:
  case CL_DEVICE_MAX_READ_IMAGE_ARGS:
  case CL_DEVICE_MAX_SAMPLERS:
  case CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:
  case CL_DEVICE_MAX_WRITE_IMAGE_ARGS:
  case CL_DEVICE_MEM_BASE_ADDR_ALIGN:
  case CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_INT:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE:
  case CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE:
  case CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF:
  case CL_DEVICE_VENDOR_ID:
    *param_value = rval.val_uint;
    break;

    /* return cl_ulong */
  case CL_DEVICE_GLOBAL_MEM_CACHE_SIZE:
  case CL_DEVICE_GLOBAL_MEM_SIZE:
  case CL_DEVICE_LOCAL_MEM_SIZE:
  case CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE:
  case CL_DEVICE_MAX_MEM_ALLOC_SIZE:
    *param_value = rval.val_ulong;
    break;

    /* return size_t */
  case CL_DEVICE_IMAGE2D_MAX_HEIGHT:
  case CL_DEVICE_IMAGE2D_MAX_WIDTH:
  case CL_DEVICE_IMAGE3D_MAX_DEPTH:
  case CL_DEVICE_IMAGE3D_MAX_HEIGHT:
  case CL_DEVICE_IMAGE3D_MAX_WIDTH:
  case CL_DEVICE_MAX_PARAMETER_SIZE:
  case CL_DEVICE_MAX_WORK_GROUP_SIZE:
  case CL_DEVICE_PROFILING_TIMER_RESOLUTION:
   *param_value = rval.val_size_t;
    break;

    /* return cl_bool */
  case CL_DEVICE_AVAILABLE:
  case CL_DEVICE_COMPILER_AVAILABLE:
  case CL_DEVICE_ENDIAN_LITTLE:
  case CL_DEVICE_ERROR_CORRECTION_SUPPORT:
  case CL_DEVICE_HOST_UNIFIED_MEMORY:
  case CL_DEVICE_IMAGE_SUPPORT:
    *param_value = rval.val_bool;
    break;

  /* other */
  default:
    fprintf(stderr, "\nError: flGetDeviceInfo not implemented param_name.\n");
    exit(1);
    break;
  }

}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(flgetdeviceinfo_int, FLGETDEVICEINFO_INT)
     (const cl_device_id * device, const int * param_name, cl_int * param_value, int * status){
  cl_long param_value64;

  FC_FUNC_(flgetdeviceinfo_int64, FLGETDEVICEINFO_INT64)(device, param_name, &param_value64, status);
  
  *param_value = (cl_int) param_value64;
}

/* -----------------------------------------------------------------------*/

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

/* -----------------------------------------------------------------------*/

/* clCreateCommandQueue */
void FC_FUNC(flcreatecommandqueue, FLCREATECOMMANDQUEUE)
     (cl_command_queue * command_queue, cl_context * context, cl_device_id * device, int * status){
  *command_queue = clCreateCommandQueue(*context, *device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, status);
}

/* -----------------------------------------------------------------------*/

/* clReleaseCommandQueue */
void FC_FUNC(flreleasecommandqueue, FLRELEASECOMMANDQUEUE)(cl_command_queue * command_queue, int * status){
  *status = clReleaseCommandQueue(*command_queue);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC(flreleasecontext, FLRELEASECONTEXT)(cl_context * context, int * status){
  *status = (int) clReleaseContext(*context);
}

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

void FC_FUNC_(f90_cl_create_buffer, F90_CL_CREATE_BUFFER)
     (cl_mem * buffer, cl_context * context, const int * flags, const size_t * size, int * status){

  *buffer = clCreateBuffer(*context, *flags, (size_t) *size, NULL, status);

}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(f90_cl_release_buffer, F90_CL_RELEASE_BUFFER)(cl_mem * buffer, int * status){

  *status = clReleaseMemObject(*buffer);
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

void FC_FUNC(flfinish, FLFINISH)(cl_command_queue * cq, int * status){
  *status = clFinish(*cq);
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

/* clEnqueueNDRangeKernel*/
void FC_FUNC(flenqueuendrangekernel, FLENQUEUENDRANGEKERNEL)
     (cl_kernel * kernel, cl_command_queue * cq, const int * dim, const size_t * globalsizes, const size_t * localsizes, int * status){

  *status = clEnqueueNDRangeKernel(*cq, *kernel, *dim,
				   NULL,  globalsizes, localsizes, 0, NULL, NULL);

}
