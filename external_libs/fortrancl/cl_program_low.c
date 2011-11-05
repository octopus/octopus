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

void FC_FUNC(clreleaseprogram, CLRELEASEPROGRAM)
     (cl_program * program, int * status){

  *status = (int) clReleaseProgram(*program);
}

