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

void FC_FUNC_(clcreateprogramwithsource_low, CLCREATEPROGRAMWITHSOURCE_LOW)
     (cl_context * context, STR_F_TYPE string, int * retcode_err, cl_program * program STR_ARG1){
  char * string_c;
  cl_int retcode_err_cl;

  TO_C_STR1(string, string_c);

  *program = clCreateProgramWithSource(*context, 1, (const char**) &string_c, NULL, &retcode_err_cl);
  *retcode_err = (int) retcode_err_cl;

  free(string_c);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(clbuildprogram_nodevices,CLBUILDPROGRAM_NODEVICES)
     (cl_program * program, STR_F_TYPE options, int * retcode_err STR_ARG1){
  char * options_c;

  TO_C_STR1(options, options_c);

  *retcode_err = (int) clBuildProgram(*program, 0, NULL, options_c, NULL, NULL);

  free(options_c);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(clgetprogrambuildinfo_str,CLGETPROGRAMBUILDINFO_STR)
     (cl_program * program, cl_device_id * device, const int * param_name, 
      STR_F_TYPE param_value, int * retcode_err STR_ARG1){
  char param_value_c[2000];

  *retcode_err = (int) clGetProgramBuildInfo(*program, *device, (cl_program_build_info) *param_name,
					     sizeof(param_value_c), param_value_c, NULL);

  TO_F_STR1(param_value_c, param_value);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC(clreleaseprogram, CLRELEASEPROGRAM)
     (cl_program * program, int * status){

  *status = (int) clReleaseProgram(*program);
}

