/*
 Copyright (C) 2016 X. Andrade

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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id$
*/

#ifndef __CUDA_COMPAT_H__
#define __CUDA_COMPAT_H__

#define __kernel extern "C" __global__
#define __global 
#define restrict __restrict__

__device__ __forceinline__ static int get_global_id(const int ii){
  switch(ii){
  case 0 :
    return threadIdx.x + blockDim.x*blockIdx.x;
  case 1 :
    return threadIdx.y + blockDim.y*blockIdx.y;
  case 2 :
    return threadIdx.z + blockDim.z*blockIdx.z;
  }
  return 0;
}

__device__ __forceinline__ static int get_global_size(const int ii){
  switch(ii){
  case 0 :
    return blockDim.x*gridDim.x;
  case 1 :
    return blockDim.y*gridDim.y;
  case 2 :
    return blockDim.z*gridDim.z;
  }
  return 0;
}

__device__ __forceinline__ static int get_local_id(const int ii){
  switch(ii){
  case 0 :
    return threadIdx.x;
  case 1 :
    return threadIdx.y;
  case 2 :
    return threadIdx.z;
  }
  return 0;
}

__device__ __forceinline__ static int get_local_size(const int ii){
  switch(ii){
  case 0 :
    return blockDim.x;
  case 1 :
    return blockDim.y;
  case 2 :
    return blockDim.z;
  }
  return 0;
}

__device__ __forceinline__ static double2 operator*(const double & a, const double2 & b){
  double2 c;
  c.x = a*b.x;
  c.y = a*b.y;
  return c;
}

__device__ __forceinline__ static double2 operator+(const double2 & a, const double2 & b){
  double2 c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return c;
}

__device__ __forceinline__ static double2 operator+=(double2 & a, const double2 & b){
  a.x += b.x;
  a.y += b.y;
  return a;
}

#endif
