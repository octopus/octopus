/*
 Copyright (C) 2019 S. Ohlmann

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

*/
#include <stdlib.h>
#include <config.h>
#include <stdio.h>
#include "vectors.h"

#ifdef HAVE_CUDA
#include <cuda_runtime.h>

inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG_ALLOC)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

#endif

void *my_allocate(int size_bytes) {
#ifndef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, unpinned.\n", (unsigned int)size_bytes);
#endif
  void *aligned;
  int status;
  // align on vector size to improve vectorization
  status = posix_memalign(&aligned, (size_t)sizeof(double)*VEC_SIZE, (unsigned int)size_bytes);
  if (status != 0) {
    printf("Error allocating aligned memory!\n");
    return NULL;
  }
  return aligned;
#else // HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, pinned.\n", (unsigned int)size_bytes);
#endif
  void *pinned;
  checkCuda(cudaMallocHost(&pinned, (unsigned int)size_bytes));
  return pinned;
#endif // HAVE_CUDA
}

void *dallocate_hardware_aware(int size) {
  return my_allocate(sizeof(double)*size);
}

void *zallocate_hardware_aware(int size) {
  return my_allocate(sizeof(double)*2*size);
}

void *sallocate_hardware_aware(int size) {
  return my_allocate(sizeof(float)*size);
}

void *callocate_hardware_aware(int size) {
  return my_allocate(sizeof(float)*2*size);
}

void deallocate_hardware_aware(void *array) {
#ifndef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Deallocating unpinned.\n");
#endif
  free(array);
#else
#ifdef DEBUG_ALLOC
  printf("Deallocating pinned.\n");
#endif
  checkCuda(cudaFreeHost(array));
#endif
}
