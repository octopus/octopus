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

*/

#include <config.h>

#ifdef HAVE_CUDA
#include <cuda.h>
#include <nvrtc.h>
// all kernels and transfers are submitted to this non-blocking stream
// -> allows operations from libraries to overlap with this stream
CUstream * phStream;
int current_stream;
static int number_streams = 32;
#else
typedef int CUcontext;
typedef int CUdevice;
typedef int CUmodule;
typedef int CUfunction;
typedef int CUdeviceptr;
typedef int CUstream;
#endif

#include <stdlib.h> //we have to include this before cmath to workaround a bug in the PGI "compiler".
#include <cmath>

#include <iostream>

#include <fstream>

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include <vector>
#include <sstream>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <stdbool.h>

#include <fortran_types.h>

#define NVRTC_SAFE_CALL(x)                                        \
  do {                                                            \
    nvrtcResult result = x;                                       \
    if (result != NVRTC_SUCCESS) {                                \
      std::cerr << "\nerror: " #x " failed with error "           \
                << nvrtcGetErrorString(result) << '\n';           \
      exit(1);                                                    \
    }                                                             \
  } while(0)

#define CUDA_SAFE_CALL(x)                                         \
  do {                                                            \
    CUresult result = x;                                          \
    if (result != CUDA_SUCCESS) {                                 \
      const char *msg;                                            \
      cuGetErrorName(result, &msg);                               \
      std::cerr << "\nerror: " #x " failed with error "           \
                << msg << '\n';                                   \
      if (result == CUDA_ERROR_OUT_OF_MEMORY) {                   \
        std::cerr << "Octopus could not allocate enough memory on the GPU.\n"; \
        std::cerr << "Please use either more GPUs to distribute the memory or try StatesPack = no to keep the states mostly on the CPU.\n"; \
      }                                                           \
      exit(1);                                                    \
    }                                                             \
  } while(0)

using namespace std;

extern "C" void FC_FUNC_(cuda_init, CUDA_INIT)(CUcontext ** context, CUdevice ** device, CUstream ** stream, fint * device_number, fint * rank){

#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuInit(0));

  *context = new CUcontext;
  *device = new CUdevice;

  int ndevices;

  CUDA_SAFE_CALL(cuDeviceGetCount(&ndevices));
  
  if (ndevices == 0) {
    cerr << "Error: no CUDA devices available." << std::endl;
    exit(1);
  }
  
  *device_number = (*device_number + *rank) % ndevices;
  CUDA_SAFE_CALL(cuDeviceGet(*device, *device_number));

  CUDA_SAFE_CALL(cuCtxCreate(*context, 0, **device));

  CUDA_SAFE_CALL(cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_L1));

  phStream = new CUstream[number_streams];
  for(current_stream = 0; current_stream < number_streams; ++current_stream) {
    CUDA_SAFE_CALL(cuStreamCreate(&phStream[current_stream], CU_STREAM_NON_BLOCKING));
  }
  current_stream = 0;
  *stream = &phStream[current_stream];
#endif
}

extern "C" void FC_FUNC_(cuda_end, CUDA_END)(CUcontext ** context, CUdevice ** device){
#ifdef HAVE_CUDA

  CUDA_SAFE_CALL(cuStreamDestroy(phStream[current_stream]));
  CUDA_SAFE_CALL(cuCtxDestroy(**context));
  
  delete *context;
  delete *device;
#endif
}

extern "C" void FC_FUNC_(cuda_module_map_init, CUDA_MODULES_MAP_INIT)(map<string, CUmodule *> ** module_map){
  *module_map = new map<string, CUmodule *>;
}

extern "C" void FC_FUNC_(cuda_module_map_end, CUDA_MODULES_MAP_END)(map<string, CUmodule *> ** module_map){

  for(map<string, CUmodule *>::iterator map_it = (**module_map).begin(); map_it != (**module_map).end(); ++map_it){
    CUmodule * module = map_it->second;
#ifdef HAVE_CUDA
    CUDA_SAFE_CALL(cuModuleUnload(*module));
#endif
    delete module;
  }
  
  delete *module_map;
}

extern "C" void FC_FUNC_(cuda_build_program, CUDA_BUILD_PROGRAM)(map<string, CUmodule *> ** module_map, CUmodule ** module, CUdevice ** device,
								 STR_F_TYPE const fname, STR_F_TYPE const flags STR_ARG2){
#ifdef HAVE_CUDA
  char *fname_c;
  char *flags_c;
    
  TO_C_STR1(fname, fname_c);
  TO_C_STR2(flags, flags_c);

  string map_descriptor = string(fname_c) + string(flags_c);

  map<string, CUmodule *>::iterator map_it = (**module_map).find(map_descriptor);
  if(map_it != (**module_map).end()){
    *module = map_it->second;
    free(fname_c);
    return;
  }
  
  // read the source
  string source;

  source = "#include \"" + string(fname_c) + "\"\n";

  // cout << source << "|" << endl;

  nvrtcProgram prog;
  NVRTC_SAFE_CALL(nvrtcCreateProgram(&prog, source.c_str(), "kernel_include.c", 0, NULL, NULL));

  int major = 0, minor = 0;
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, **device));
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, **device));

  char compute_version[3];
  sprintf(compute_version, "%.1d%.1d", major, minor);

  string all_flags = "--gpu-architecture=compute_" + string(compute_version)
    + " --ftz=true --fmad=true -DCUDA -default-device " + string(flags_c);

  stringstream flags_stream(all_flags);
  istream_iterator<string> iter(flags_stream);
  istream_iterator<string> end;
  vector<string> tokens(iter, end);
  
  const char ** opts = new const char*[tokens.size()];
  for (unsigned ii = 0; ii < tokens.size(); ii++) opts[ii] = tokens[ii].c_str();



  nvrtcResult err = nvrtcCompileProgram(prog, tokens.size(), opts);

  free(flags_c);
  
  size_t logSize;
  NVRTC_SAFE_CALL(nvrtcGetProgramLogSize(prog, &logSize));
  char *log = new char[logSize];
  NVRTC_SAFE_CALL(nvrtcGetProgramLog(prog, log));

  if(logSize > 1){
    
    cout << "Cuda compilation messages" << endl;
    
    cout << "File    : " << fname_c << endl;
    
    cout << "Options : " << all_flags << endl;
    
    cout << log << endl;

  }

  if(NVRTC_SUCCESS != err){
    cerr << "Error in compiling" << endl;
    exit(1);
  }
  
  delete [] log;
  delete [] opts;

  // Obtain PTX from the program.
  size_t ptxSize;
  NVRTC_SAFE_CALL(nvrtcGetPTXSize(prog, &ptxSize));
  char *ptx = new char[ptxSize];
  NVRTC_SAFE_CALL(nvrtcGetPTX(prog, ptx));

  NVRTC_SAFE_CALL(nvrtcDestroyProgram(&prog));

  *module = new CUmodule;

  const int num_options = 2;
  CUjit_option options[num_options];
  void * option_values[num_options];

  unsigned log_size = 4096;
  char log_buffer[log_size];
  
  options[0] = CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES;
  option_values[0] = (void *) (long)log_size;

  options[1] = CU_JIT_ERROR_LOG_BUFFER;
  option_values[1] = (void *) log_buffer;
  
  CUresult result = cuModuleLoadDataEx(*module, ptx, num_options, options, option_values);

  if(result != CUDA_SUCCESS){
    std::cerr << log_buffer << std::endl;
    const char *msg;
    cuGetErrorName(result, &msg);
    std::cerr << "\nerror: cuModuleLoadDataEx failed with error " << msg << '\n';
    exit(1);
  }
      
  delete [] ptx;

  (**module_map)[map_descriptor] = *module;

  free(fname_c);
#endif
}

extern "C" void FC_FUNC_(cuda_create_kernel, CUDA_CREATE_KERNEL)(CUfunction ** kernel, CUmodule ** module, STR_F_TYPE kernel_name STR_ARG1){
#ifdef HAVE_CUDA  
  char *kernel_name_c;

  TO_C_STR1(kernel_name, kernel_name_c);
  
  *kernel = new CUfunction;

  CUDA_SAFE_CALL(cuModuleGetFunction(*kernel, **module, kernel_name_c));

  free(kernel_name_c);
#endif
}

extern "C" void FC_FUNC_(cuda_release_module, CUDA_RELEASE_MODULE)(CUmodule ** module){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuModuleUnload(**module));
  delete *module;
#endif
}

extern "C" void FC_FUNC_(cuda_release_kernel, CUDA_RELEASE_KERNEL)(CUfunction ** kernel){
#ifdef HAVE_CUDA
  delete *kernel;
#endif
}

extern "C" void FC_FUNC_(cuda_device_max_threads_per_block, CUDA_DEVICE_MAX_THREADS_PER_BLOCK)(CUdevice ** device, fint * max_threads){
#ifdef HAVE_CUDA
  int value;
  CUDA_SAFE_CALL(cuDeviceGetAttribute (&value, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, **device));
  *max_threads = value;
#endif
}

extern "C" void FC_FUNC_(cuda_device_total_memory, CUDA_DEVICE_TOTAL_MEMORY)(CUdevice ** device, fint8 * total_memory){
#ifdef HAVE_CUDA
  size_t mem;
  CUDA_SAFE_CALL(cuDeviceTotalMem(&mem, **device));
  *total_memory = mem;
#endif
}

extern "C" void FC_FUNC_(cuda_device_shared_memory, CUDA_DEVICE_SHARED_MEMORY)(CUdevice ** device, fint8 * shared_memory){
#ifdef HAVE_CUDA
  int mem;
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&mem, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, **device));
  *shared_memory = mem;
#endif
}

extern "C" void FC_FUNC_(cuda_mem_alloc, CUDA_MEM_ALLOC)(CUdeviceptr ** cuda_ptr, const fint8 * size){
#ifdef HAVE_CUDA
  *cuda_ptr = new CUdeviceptr;
  CUDA_SAFE_CALL(cuMemAlloc(*cuda_ptr, *size));
#endif
}

extern "C" void FC_FUNC_(cuda_mem_free, CUDA_MEM_FREE)(CUdeviceptr ** cuda_ptr){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuMemFree(**cuda_ptr));
  delete *cuda_ptr;
#endif
}

extern "C" void FC_FUNC_(cuda_memcpy_htod, CUDA_MEMCPY_HTOD)(CUdeviceptr ** cuda_ptr, const void * data, fint8 * size, fint8 * offset, bool * async){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuMemcpyHtoDAsync(**cuda_ptr + *offset, data, *size, phStream[current_stream]));
  if(!(*async)) CUDA_SAFE_CALL(cuStreamSynchronize(phStream[current_stream]));
#endif  
}

extern "C" void FC_FUNC_(cuda_memcpy_dtoh, CUDA_MEMCPY_DTOH)(CUdeviceptr ** cuda_ptr, void * data, fint8 * size, fint8 * offset, bool * async){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuMemcpyDtoHAsync(data, **cuda_ptr + *offset, *size, phStream[current_stream]));
  if(!(*async)) CUDA_SAFE_CALL(cuStreamSynchronize(phStream[current_stream]));
#endif  
}

extern "C" void FC_FUNC_(cuda_alloc_arg_array, CUDA_ALLOC_ARG_ARRAY)(vector<void *> ** arg_array){
  *arg_array = new vector<void *>;  
}

extern "C" void FC_FUNC_(cuda_free_arg_array, CUDA_FREE_ARG_ARRAY)(vector<void *> ** arg_array){

  for(unsigned ii = 0; ii < (**arg_array).size(); ii++) free((**arg_array)[ii]);
  delete *arg_array;

}

extern "C" void FC_FUNC_(cuda_kernel_set_arg_buffer, CUDA_KERNEL_SET_ARG_BUFFER)
  (vector<void *> ** arg_array, CUdeviceptr ** cuda_ptr, fint * arg_index){

  if(unsigned(*arg_index) >= (**arg_array).size()) (**arg_array).resize(*arg_index + 1, NULL);
  
  if((**arg_array)[*arg_index] == NULL) (**arg_array)[*arg_index] = malloc(sizeof(CUdeviceptr));
  
  memcpy((**arg_array)[*arg_index], *cuda_ptr, sizeof(CUdeviceptr));
}

extern "C" void FC_FUNC_(cuda_kernel_set_arg_value, CUDA_KERNEL_SET_ARG_VALUE)
  (vector<void *> ** arg_array, void * arg, fint * arg_index, fint * size){
  
  if(unsigned(*arg_index) >= (**arg_array).size()) (**arg_array).resize(*arg_index + 1, NULL);
  
  if((**arg_array)[*arg_index] == NULL) (**arg_array)[*arg_index] = malloc(*size);

  memcpy((**arg_array)[*arg_index], arg, *size);
  
}

extern "C" void FC_FUNC_(cuda_context_synchronize, CUDA_CONTEXT_SYNCHRONIZE)(){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuStreamSynchronize(phStream[current_stream]));
#endif
}

extern "C" void FC_FUNC_(cuda_synchronize_all_streams, CUDA_SYNCHRONIZE_ALL_STREAMS)(){
#ifdef HAVE_CUDA
  for(int i = 0; i < number_streams; ++i)
    CUDA_SAFE_CALL(cuStreamSynchronize(phStream[i]));
#endif
}

extern "C" void FC_FUNC_(cuda_launch_kernel, CUDA_LAUNCH_KERNEL)
  (CUfunction ** kernel, fint8 * griddim, fint8 * blockdim, fint8 * shared_mem, vector<void *> ** arg_array){
#ifdef HAVE_CUDA
  /*
  cout << "Kernel call" << endl;

  int nn;
  CUDA_SAFE_CALL(cuFuncGetAttribute(&nn, CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, **kernel));
  cout << "SIZE   " << nn << endl;
  CUDA_SAFE_CALL(cuFuncGetAttribute(&nn, CU_FUNC_ATTRIBUTE_PTX_VERSION, **kernel));
  cout << "PTX    " << nn << endl;
  CUDA_SAFE_CALL(cuFuncGetAttribute(&nn, CU_FUNC_ATTRIBUTE_BINARY_VERSION, **kernel));
  cout << "BINARY " << nn << endl;

  for(unsigned ii = 0; ii < (**arg_array).size(); ii++) cout << ii << " " << (**arg_array)[ii] << endl;

  cout << "GRID  " << griddim[0] << " " << griddim[1] << " " <<  griddim[2] << endl;
  cout << "BLOCK " << blockdim[0] << " " << blockdim[1] << " " <<  blockdim[2] << endl;
  */
  
  assert((**arg_array).size() > 0);
  for(unsigned ii = 0; ii < (**arg_array).size(); ii++) assert((**arg_array)[ii] != NULL);
  
  CUDA_SAFE_CALL(cuLaunchKernel(**kernel, griddim[0], griddim[1], griddim[2],
         blockdim[0], blockdim[1], blockdim[2], *shared_mem, phStream[current_stream], &(**arg_array)[0], NULL));

  // release the stored argument, this is not necessary in principle,
  // but it should help us to detect missing arguments.
  for(unsigned ii = 0; ii < (**arg_array).size(); ii++) free((**arg_array)[ii]);
  (**arg_array).resize(0);
#endif
}

extern "C" void FC_FUNC_(cuda_device_name, CUDA_DEVICE_NAME)(CUdevice ** device, STR_F_TYPE name STR_ARG1){
#ifdef HAVE_CUDA
  char devicename[200];
  CUDA_SAFE_CALL(cuDeviceGetName(devicename, sizeof(devicename), **device));
  TO_F_STR1(devicename, name);
  
#endif
  
}

extern "C" void FC_FUNC_(cuda_device_capability, CUDA_DEVICE_CAPABILITY)(CUdevice ** device, fint * major, fint * minor){
#ifdef HAVE_CUDA
  int cmajor = 0, cminor = 0;
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&cmajor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, **device));
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&cminor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, **device));
  *major = cmajor;
  *minor = cminor;
#endif
}

extern "C" void FC_FUNC_(cuda_driver_version, CUDA_DRIVER_VERSION)(fint * version){
#ifdef HAVE_CUDA
  int driverversion;
  CUDA_SAFE_CALL(cuDriverGetVersion(&driverversion));
  *version = driverversion;
#endif
}

extern "C" void FC_FUNC_(cuda_device_get_warpsize, CUDA_DEVICE_GET_WARPSIZE)(CUdevice ** device, fint * warpSize){
#ifdef HAVE_CUDA
  int cwarpSize=0;
  CUDA_SAFE_CALL(cuDeviceGetAttribute(&cwarpSize, CU_DEVICE_ATTRIBUTE_WARP_SIZE, **device));
  *warpSize = cwarpSize;
#endif
}

extern "C" void FC_FUNC_(cuda_deref, CUDA_DEREF)(CUdeviceptr ** cuda_ptr, void ** cuda_deref_ptr) {
  *cuda_deref_ptr = (void *) **cuda_ptr;
}

extern "C" void FC_FUNC_(cuda_set_stream, CUDA_SET_STREAM)(CUstream ** stream, fint * number) {
#ifdef HAVE_CUDA
  current_stream = (*number - 1) % number_streams;
  *stream = &phStream[current_stream];
#endif
}
