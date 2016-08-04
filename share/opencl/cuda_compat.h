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
#define __constant       // perhaps there is a type of memory for this
#define __local __shared__
#define restrict __restrict__
#define barrier(x) __syncthreads()

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

#define double2 Double2

class __align__(16) double2{
 public:
  double x;
  double y;
  __forceinline__ double2(const double & a = 0, const double & b = 0):x(a), y(b){}
};

__forceinline__ static Double2 operator*(const double & a, const Double2 & b){
  Double2 c;
  c.x = a*b.x;
  c.y = a*b.y;
  return c;
}

__forceinline__ static Double2 operator*(const Double2 & a, const Double2 & b){
  Double2 c;
  c.x = a.x*b.x;
  c.y = a.y*b.y;
  return c;
}

__forceinline__ static Double2 operator+(const Double2 & a, const Double2 & b){
  Double2 c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return c;
}

__forceinline__ static Double2 operator+=(Double2 & a, const Double2 & b){
  a.x += b.x;
  a.y += b.y;
  return a;
}


__forceinline__ static Double2 operator-=(Double2 & a, const Double2 & b){
  a.x -= b.x;
  a.y -= b.y;
  return a;
}

__forceinline__ static Double2 operator/(const Double2 & a, const double & b){
  Double2 c;
  c.x = a.x/b;
  c.y = a.y/b;
  return c;
}

__forceinline__ static double sincos(const double & x, double * cosx){
  double sinx;
  sincos(x, &sinx, cosx);
  return sinx;
}

#endif
