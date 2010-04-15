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

 $Id: operate_opencl.c 2146 2006-05-23 17:36:00Z xavier $
*/


#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>

#include "string_f.h" /* fortran <-> c string compatibility issues */
#include "opencl.h"

typedef struct{
  cl_kernel kernel_double;
  cl_kernel kernel_complex;
} operate_opencl_t;

void FC_FUNC_(operate_opencl_init,OPERATE_OPENCL_INIT)(operate_opencl_t ** thisptr, opencl_t ** ocl, STR_F_TYPE sourcepath_f STR_ARG1){
  FILE * source_file;
  char* cSourceString;
  size_t szSourceLength;
  char * sourcepath;
  int numerr;
  cl_program OpenCLProgram;
  operate_opencl_t * this;

  this = (operate_opencl_t *) malloc(sizeof(operate_opencl_t));
  *thisptr = this;

  TO_C_STR1(sourcepath_f, sourcepath);

  // open the OpenCL source code file
  source_file = fopen(sourcepath, "rb");
  if(source_file == 0){       
    fprintf(stderr, "Error: Failed to open file %s\n", sourcepath);
    exit(1);
  } else {
    printf("Info: compiling OpenCL code %s\n", sourcepath);
  }

  /* get the length of the source code */
  fseek(source_file, 0, SEEK_END); 
  szSourceLength = ftell(source_file);
  fseek(source_file, 0, SEEK_SET); 
  
  /* allocate a buffer for the source code string and read it in */
  cSourceString = (char *)malloc(szSourceLength + 1);
  fread(cSourceString, szSourceLength, 1, source_file);
  fclose(source_file);
    
  cSourceString[szSourceLength] = '\0';

  OpenCLProgram = clCreateProgramWithSource(ocl[0]->Context,1, (const char**)&cSourceString, NULL, &numerr); 
  numerr = clBuildProgram(OpenCLProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
  
  if(numerr != CL_SUCCESS){
    fprintf(stderr, "Error: compilation of file %s failed.\n", sourcepath);
    exit(1);
  }

  /* Create a handle to the compiled OpenCL function (Kernel) */
  this->kernel_double = clCreateKernel(OpenCLProgram, "operatever7", NULL);

  free(sourcepath);
}

void FC_FUNC_(operate_opencl_end,OPERATE_OPENCL_END)(operate_opencl_t ** thisptr){
  free(*thisptr);
}
