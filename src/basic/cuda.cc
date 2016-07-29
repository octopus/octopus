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

#include <config.h>

#ifdef HAVE_CUDA
#include <cuda.h>
#include <nvrtc.h>
#else
typedef int CUcontext;
typedef int CUdevice;
typedef int CUmodule;
typedef int CUfunction;
typedef int CUdeviceptr;
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
      exit(1);                                                    \
    }                                                             \
  } while(0)

using namespace std;

extern "C" void FC_FUNC_(cuda_init, CUDA_INIT)(CUcontext ** context, CUdevice ** device){

#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuInit(0));

  *context = new CUcontext;
  *device = new CUdevice;

  int ndevices;

  CUDA_SAFE_CALL(cuDeviceGetCount(&ndevices));
  
  cout << "Number of devices : " << ndevices << endl;

  if (ndevices == 0) {
    cerr << "Error: no CUDA devices" << std::endl;
    exit(1);
  }
  
  CUDA_SAFE_CALL(cuDeviceGet(*device, 0));

  char devicename[200];

  CUDA_SAFE_CALL(cuDeviceGetName(devicename, sizeof(devicename), **device));

  cout << "CUDA device name : " << devicename << endl;

  int major = 0, minor = 0;
  CUDA_SAFE_CALL(cuDeviceComputeCapability(&major, &minor, **device));
  cout << "CUDA compute capability : " <<  major << '.' <<  minor << endl;
  
  size_t mem;
  CUDA_SAFE_CALL(cuDeviceTotalMem(&mem, **device));
  cout << "CUDA device memory : " << mem/pow(1024.0, 3) << " GB" << endl;

  CUDA_SAFE_CALL(cuCtxCreate(*context, 0, **device));
#endif
}

extern "C" void FC_FUNC_(cuda_end, CUDA_END)(CUcontext ** context, CUdevice ** device){
#ifdef HAVE_CUDA

  CUDA_SAFE_CALL(cuCtxDestroy(**context));
  
  delete *context;
  delete *device;
#endif
}

extern "C" void FC_FUNC_(cuda_build_program, CUDA_BUILD_PROGRAM)(CUmodule ** module, STR_F_TYPE include_path, STR_F_TYPE const fname, STR_F_TYPE const flags STR_ARG3){
#ifdef HAVE_CUDA
  char *include_path_c;
  char *fname_c;
  char *flags_c;

  TO_C_STR1(include_path, include_path_c);
  TO_C_STR2(fname, fname_c);
  TO_C_STR3(flags, flags_c);
  
  // read the source
  
  string source;

  source = "#include \"" + string(fname_c) + "\"\n";

  // cout << source << "|" << endl;

  cout << fname_c << endl;

  nvrtcProgram prog;
  NVRTC_SAFE_CALL(nvrtcCreateProgram(&prog, source.c_str(), "kernel_include.c", 0, NULL, NULL));
  free(fname_c);

  string all_flags = "--gpu-architecture=compute_52 -DCUDA -default-device " + string("-I") + include_path_c + string(" ") + string(flags_c);

  stringstream flags_stream(all_flags);
  istream_iterator<string> iter(flags_stream);
  istream_iterator<string> end;
  vector<string> tokens(iter, end);
  
  const char ** opts = new const char*[tokens.size()];
  for (unsigned ii = 0; ii < tokens.size(); ii++) opts[ii] = tokens[ii].c_str();

  cout << "OPTS " << all_flags << endl;

  nvrtcResult err = nvrtcCompileProgram(prog, tokens.size(), opts);

  free(flags_c);
  
  size_t logSize;
  NVRTC_SAFE_CALL(nvrtcGetProgramLogSize(prog, &logSize));
  char *log = new char[logSize];
  NVRTC_SAFE_CALL(nvrtcGetProgramLog(prog, log));

  cout << "=====================================" << endl;

  cout << log << endl;

  if(NVRTC_SUCCESS != err){
    cerr << "Error in compiling" << endl;
    exit(1);
  }
  
  delete [] log;

  // Obtain PTX from the program.
  size_t ptxSize;
  NVRTC_SAFE_CALL(nvrtcGetPTXSize(prog, &ptxSize));
  char *ptx = new char[ptxSize];
  NVRTC_SAFE_CALL(nvrtcGetPTX(prog, ptx));

  NVRTC_SAFE_CALL(nvrtcDestroyProgram(&prog));

  *module = new CUmodule;

  CUDA_SAFE_CALL(cuModuleLoadDataEx(*module, ptx, 0, 0, 0));

  delete [] ptx;
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

extern "C" void FC_FUNC_(cuda_memcpy_htod, CUDA_MEMCPY_HTOD)(CUdeviceptr ** cuda_ptr, const void * data, fint8 * size, fint8 * offset){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuMemcpyHtoD(**cuda_ptr + *offset, data, *size));
#endif  
}

extern "C" void FC_FUNC_(cuda_memcpy_dtoh, CUDA_MEMCPY_DTOH)(CUdeviceptr ** cuda_ptr, void * data, fint8 * size, fint8 * offset){
#ifdef HAVE_CUDA
  CUDA_SAFE_CALL(cuMemcpyDtoH(data, **cuda_ptr + *offset, *size));
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
  CUDA_SAFE_CALL(cuCtxSynchronize());
#endif
}

extern "C" void FC_FUNC_(cuda_launch_kernel, CUDA_LAUNCH_KERNEL)
  (CUfunction ** kernel, fint8 * griddim, fint8 * blockdim, vector<void *> ** arg_array){
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
  				blockdim[0], blockdim[1], blockdim[2], 0, NULL, &(**arg_array)[0], NULL));

  // release the stored argument, this is not necessary in principle,
  // but it should help us to detect missing arguments.
  for(unsigned ii = 0; ii < (**arg_array).size(); ii++) free((**arg_array)[ii]);
  (**arg_array).resize(0);
 
#endif
}

