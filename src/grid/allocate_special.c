#include <stdlib.h>
#include <config.h>

#ifdef HAVE_CUDA
#include <cuda_runtime.h>

inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG_ALLOC)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
#endif
  return result;
}

#endif

void *my_allocate(int size_bytes) {
#ifndef HAVE_CUDA
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, unpinned.\n", (unsigned int)size_bytes);
#endif
  return malloc(size_bytes);
#else
#ifdef DEBUG_ALLOC
  printf("Allocating %d bytes, pinned.\n", (unsigned int)size_bytes);
#endif
  void *pinned;
  checkCuda(cudaMallocHost(&pinned, (unsigned int)size_bytes));
  return pinned;
#endif
}

void *dallocate_special(int size) {
  return my_allocate(sizeof(double)*size);
}

void *zallocate_special(int size) {
  return my_allocate(sizeof(double)*2*size);
}

void *sallocate_special(int size) {
  return my_allocate(sizeof(float)*size);
}

void *callocate_special(int size) {
  return my_allocate(sizeof(float)*2*size);
}

void deallocate_special(void *array) {
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
