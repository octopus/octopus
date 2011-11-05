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


void FC_FUNC(clgetplatformids_num, CLGETPLATFORMIDS_NUM)(int * num_platforms, int * status){
  cl_uint ret_platform;

  *status = (int) clGetPlatformIDs(0, NULL, &ret_platform);
  *num_platforms = (int) ret_platform;
}


/* -----------------------------------------------------------------------*/

void FC_FUNC_(clgetplatformids_listall, CLGETPLATFORMIDS_LISTALL)
     (const int * num_entries, cl_platform_id * platforms, int * num_platforms, int * status){

  cl_uint unum_platforms;

  *status = (int) clGetPlatformIDs((cl_uint) *num_entries, platforms, &unum_platforms);
  *num_platforms = (int) unum_platforms;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(clgetplatformids_getplat, CLGETPLATFORMIDS_GETPLAT)
     (const cl_platform_id * allplatforms, const int * iplatform, cl_platform_id * platform){
  *platform = allplatforms[*iplatform];
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(clgetplatforminfo, CLGETPLATFORMINFO)
     (const cl_platform_id * platform, const int * param_name, STR_F_TYPE param_value, int * status STR_ARG1){
  char info[2048];

  *status = (int) clGetPlatformInfo(*platform, (cl_platform_info) *param_name, sizeof(info), info, NULL);

  TO_F_STR1(info, param_value);
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(clgetdeviceids_setdev, CLGETDEVICEIDS_SETDEV)
     (cl_device_id * alldevices, const int * idevice, const cl_device_id * device){
  alldevices[*idevice] = *device;
}



