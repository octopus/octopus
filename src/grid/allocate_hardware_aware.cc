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
#include <iostream>

#ifdef HAVE_CUDA
#include <cuda.h>

#define CUDA_SAFE_CALL(x)                                         \
  do {                                                            \
    CUresult result = x;                                          \
    if (result != CUDA_SUCCESS) {                                 \
      const char *msg;                                            \
      cuGetErrorName(result, &msg);                               \
      std::cerr << "\nerror: " #x " failed with error "           \
                << msg << '\n';                                   \
      exit(1);                                                    \
    }                                                             \
  } while(0)

#endif

using namespace std;

void *allocate_aligned(int size_bytes) {
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
}

extern "C" void *dallocate_aligned(int size) {
  return allocate_aligned(sizeof(double)*size);
}

extern "C" void *zallocate_aligned(int size) {
  return allocate_aligned(sizeof(double)*2*size);
}

extern "C" void *sallocate_aligned(int size) {
  return allocate_aligned(sizeof(float)*size);
}

extern "C" void *callocate_aligned(int size) {
  return allocate_aligned(sizeof(float)*2*size);
}

extern "C" void deallocate_aligned(void *array) {
#ifdef DEBUG_ALLOC
  printf("Deallocating unpinned.\n");
#endif
  free(array);
}

void *allocate_pinned(int size_bytes) {
#ifdef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, pinned.\n", (unsigned int)size_bytes);
#endif
  void *pinned;
  CUDA_SAFE_CALL(cuMemAllocHost(&pinned, (unsigned int)size_bytes));
  return pinned;
#else
  printf("Error! Pinned memory requested, although CUDA not available. Returning aligned memory.");
  return allocate_aligned(size_bytes);
#endif
}

extern "C" void *dallocate_pinned(int size) {
  return allocate_pinned(sizeof(double)*size);
}

extern "C" void *zallocate_pinned(int size) {
  return allocate_pinned(sizeof(double)*2*size);
}

extern "C" void *sallocate_pinned(int size) {
  return allocate_pinned(sizeof(float)*size);
}

extern "C" void *callocate_pinned(int size) {
  return allocate_pinned(sizeof(float)*2*size);
}

extern "C" void deallocate_pinned(void *array) {
#ifdef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Deallocating pinned.\n");
#endif
  CUDA_SAFE_CALL(cuMemFreeHost(array));
#endif
}
